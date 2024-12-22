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

   integer :: num_trials
   integer, allocatable :: atomperm(:)

   integer :: natom1
   real(rk), dimension(:, :), allocatable :: coords1, coords2

   type(bipartition_type) :: mnatypes, submnatypes
   type(metapartition_type) :: metatypes
   integer :: h, i
   integer :: adjd
   real(rk) :: rmsd, dist

   natom1 = size(mol1%atoms)
   coords1 = mol1%get_weighted_coords()
   coords2 = mol2%get_weighted_coords()
   allocate (atomperm(natom1))

   ! Compute mnatypes and metatypes

   mnatypes = eltypes
   call compute_crossmnatypes(mol1, mol2, mnatypes)
   call mnatypes%print_parts()

   num_trials = 0
   max_trials = 1000

   ! Loop for map searching

   do while (num_trials < max_trials)

      num_trials = num_trials + 1
      call rotate_coords(coords2, randrotquat())
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
!      call submnatypes%print_parts()
      call assign_atoms(submnatypes, coords1, coords2, atomperm, dist)

      ! Update results
      adjd = adjacencydiff(atomperm, mol1%adjmat, mol2%adjmat)
      rmsd = sqrt(leastsqdistsum(atomperm, coords1, coords2))
      call results%push(atomperm, 1, 0._rk, adjd, rmsd)

   end do

   results%num_trials = num_trials

end subroutine

end module
