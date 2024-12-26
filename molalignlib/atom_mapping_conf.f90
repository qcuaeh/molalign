! MolAlignLib
! Copyright (C) 2022 José M. Vásquez

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module atom_mapping
use parameters
use globals
use random
use molecule
use strutils
use chemdata
use permutation
use rigid_body
use rotation
use alignment
use lap_driver
use lap_solvers
use adjacency
use partition
use bipartition
use partitioning
use bipartitioning
use biasing
use pruning
use printing
!use backtracking

implicit none

contains

subroutine map_atoms(mol1, mol2, eltypes, results)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_type), intent(in) :: eltypes
   type(registry_type), target, intent(inout) :: results

   ! Local variables
   integer :: natom1
   integer :: num_trials, num_steps
   integer, pointer :: lead_count
   real(rk) :: eigquat(4), totquat(4)
   integer, dimension(:), allocatable :: atomperm, auxperm
   real(rk), dimension(:,:), allocatable :: coords1, coords2
   type(bipartition_type) :: mnatypes
   type(metapartition_type) :: metatypes
   integer :: adjd
   real(rk) :: rmsd, dist
   integer :: h, i

   natom1 = size(mol1%atoms)
   coords1 = mol1%get_weighted_coords()
   coords2 = mol2%get_weighted_coords()
   allocate (atomperm(natom1))
   allocate (auxperm(natom1))

   ! Compute MNA types
   mnatypes = eltypes
   call compute_crossmnatypes(mol1, mol2, mnatypes)
!   call mnatypes%print_parts()

   ! Reflect atoms
   if (mirror_flag) then
      call reflect_coords(coords2)
   end if

   ! Translate atoms to their centroids
   call translate_coords(coords1, -centroid(coords1))
   call translate_coords(coords2, -centroid(coords2))

   ! Initialize random number generator
   call random_initialize()

   ! Optimize atom permutation

   num_trials = 0
   lead_count => results%records(1)%count

   do while (lead_count < max_count .and. num_trials < max_trials)

      num_trials = num_trials + 1

      ! Aply a random rotation to coords2
      call rotate_coords(coords2, randrotquat())

      ! Assign atoms with current orientation
      call assign_atoms_conf(mnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
      totquat = leasteigquat(atomperm, coords1, coords2)
      call rotate_coords(coords2, totquat)
      num_steps = 1

      do while (iter_flag)
         call assign_atoms_conf(mnatypes, mol1, mol2, coords1, coords2, auxperm, dist)
         if (all(auxperm == atomperm)) exit
         atomperm = auxperm
         eigquat = leasteigquat(atomperm, coords1, coords2)
         call rotate_coords(coords2, eigquat)
         totquat = quatmul(eigquat, totquat)
         num_steps = num_steps + 1
      end do

      ! Update results
      adjd = adjacencydiff(atomperm, mol1%adjmat, mol2%adjmat)
      rmsd = sqrt(sqdistsum(atomperm, coords1, coords2))
      call results%push(atomperm, num_steps, angle(totquat), adjd, rmsd)

   end do

   results%num_trials = num_trials

end subroutine

subroutine assign_atoms_conf( mnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_type), intent(in) :: mnatypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   integer, dimension(:), intent(out) :: atomperm
   real(rk), intent(out) :: dist
   ! Local variables
   type(bipartition_type) :: submnatypes
   type(metapartition_type) :: metatypes
   integer :: h, i

   submnatypes = mnatypes
   call collect_degenerated_mnatypes(mol1, submnatypes%first_partition(), metatypes)
   do while (metatypes%num_parts > 0)
      do i = 1, metatypes%num_parts
         h = random_element(metatypes%parts(i)%items)
         call minperm(submnatypes%parts(h), coords1, coords2, atomperm, dist)
         call split_crossmnatypes(h, atomperm, submnatypes)
      end do
      call compute_crossmnatypes(mol1, mol2, submnatypes)
      call collect_degenerated_mnatypes(mol1, submnatypes%first_partition(), metatypes)
   end do
   call assign_atoms(submnatypes, coords1, coords2, atomperm, dist)

end subroutine

!subroutine assign_atoms_conf( mnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
!   type(mol_type), intent(in) :: mol1, mol2
!   type(bipartition_type), intent(in) :: mnatypes
!   real(rk), dimension(:, :), intent(in) :: coords1, coords2
!   integer, dimension(:), intent(out) :: atomperm
!   real(rk), intent(out) :: dist
!   ! Local variables
!   type(bipartition_type) :: submnatypes
!   type(metapartition_type) :: metatypes
!   integer :: h, i
!
!   submnatypes = mnatypes
!   call assign_atoms_conf_rec(submnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
!   call assign_atoms(submnatypes, coords1, coords2, atomperm, dist)
!
!end subroutine

recursive subroutine assign_atoms_conf_rec( submnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_type), intent(inout) :: submnatypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   integer, dimension(:), intent(out) :: atomperm
   real(rk), intent(out) :: dist
   ! Local variables
   type(metapartition_type) :: metatypes
   integer :: h, i

   call collect_degenerated_mnatypes(mol1, submnatypes%first_partition(), metatypes)
!   call metatypes%print_parts()
!   call submnatypes%print_parts()

   do i = 1, metatypes%num_parts
!      write (stderr, *) 'loop:', i
      h = random_element(metatypes%parts(i)%items)
      call minperm(submnatypes%parts(h), coords1, coords2, atomperm, dist)
      call split_crossmnatypes(h, atomperm, submnatypes)
      call compute_crossmnatypes(mol1, mol2, submnatypes)
      call assign_atoms_conf_rec(submnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
   end do

end subroutine

end module
