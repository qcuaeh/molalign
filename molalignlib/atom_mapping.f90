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
use chemdata
use rotation
use rigid_body
use alignment
use lap_driver
use bipartition
use pruning
use printing
use registry

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
   real(rk), dimension(:, :), allocatable :: coords1, coords2
   type(boolmatrix_type), allocatable :: prunemask(:)
   real(rk) :: rmsd

   natom1 = size(mol1%atoms)
   coords1 = mol1%get_weighted_coords()
   coords2 = mol2%get_weighted_coords()
   allocate (auxperm(natom1))
   allocate (atomperm(natom1))

   ! Reflect atoms
   if (mirror_flag) then
      call reflect_coords(coords2)
   end if

   ! Translate atoms to their centroids
   call translate_coords(coords1, -centroid(coords1))
   call translate_coords(coords2, -centroid(coords2))

   ! Find unfeasible assignments
   call prune_procedure(eltypes, mol1%get_coords(), mol2%get_coords(), prunemask)

   ! Initialize counters
   num_trials = 0
   lead_count => results%records(1)%count

   ! Initialize random number generator
   call random_initialize()

   ! Loop for map searching
   do while (lead_count < max_count .and. num_trials < max_trials)

      num_trials = num_trials + 1

      ! Aply a random rotation to coords2
      call rotate_coords(coords2, randrotquat())

      ! Assign atoms with current orientation
      call assign_atoms_pruned(eltypes, coords1, coords2, prunemask, atomperm)
      totquat = leasteigquat(atomperm, coords1, coords2)
      call rotate_coords(coords2, totquat)
      num_steps = 1

      do while (iter_flag)
         call assign_atoms_pruned(eltypes, coords1, coords2, prunemask, auxperm)
         if (all(auxperm == atomperm)) exit
         atomperm = auxperm
         eigquat = leasteigquat(atomperm, coords1, coords2)
         call rotate_coords(coords2, eigquat)
         totquat = quatmul(eigquat, totquat)
         num_steps = num_steps + 1
      end do

      ! Update results
      rmsd = sqrt(leastsqdistsum(atomperm, coords1, coords2))
      call results%push(atomperm, num_steps, angle(totquat), 0, rmsd)

   end do

   results%num_trials = num_trials

end subroutine

end module
