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

module molalignlib
use stdio
use kinds
use molecule
use bounds
use random
use linalg
use sorting
use permutation
use rotation
use translation
use assorting
use adjacency
use alignment
use remapping
use assignment
use writemol
use biasing
!use reactivity

implicit none

contains

! Assign atoms0 and atoms1
subroutine molecule_remap( &
   mol1, &
   mol2, &
   eltypes, &
   nrec, &
   maplist, &
   countlist)

   type(molecule_type), intent(inout) :: mol1, mol2
   type(partition_type), intent(in) :: eltypes
   integer, intent(out) :: nrec
   integer, dimension(:, :), intent(inout) :: maplist
   integer, dimension(:), intent(inout) :: countlist

   real(rk) :: travec1(3), travec2(3)
   real(rk), allocatable, dimension(:, :) :: coords1, coords2

   ! Abort if molecules have different number of atoms

   if (size(mol1%atoms) /= size(mol2%atoms)) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Abort if molecules are not isomers

   if (any(sorted(mol1%atoms%elnum) /= sorted(mol2%atoms%elnum))) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Abort if there are conflicting atomic types

   if (any(sorted(eltypes%partition_map1) /= sorted(eltypes%partition_map2))) then
      write (stderr, '(a)') 'Error: There are conflicting atomic types'
      stop
   end if

   ! Backup coordinates

   coords1 = mol1%get_coords()
   coords2 = mol2%get_coords()

   ! Mirror coordinates

   if (mirror_flag) then
      call mol2%mirror_coords()
   end if

   ! Calculate centroids

   travec1 = -centroid(mol1)
   travec2 = -centroid(mol2)

   ! Center coordinates at the centroids

   call mol1%translate_coords(travec1)
   call mol2%translate_coords(travec2)

   ! Initialize random number generator

   call random_initialize()

   ! Optimize assignment to minimize the AdjD and RMSD

   call remap_atoms(mol1, mol2, eltypes, maplist, countlist, nrec)

   ! Remove bonds from reactive sites and reoptimize assignment

!   if (reac_flag) then
!      call remove_reactive_bonds(mol1, mol2, maplist(:, 1))
!      call remap_atoms(mol1, mol2, maplist, countlist, nrec)
!   end if

   ! Restore coordinates

   call mol1%set_coords(coords1)
   call mol2%set_coords(coords2)

end subroutine

! Align atoms
subroutine remapped_molecule_align( &
   mol1, &
   mol2, &
   mapping, &
   travec1, &
   travec2, &
   rotquat)

   type(molecule_type), intent(in) :: mol1, mol2
   integer, intent(in) :: mapping(:)
   real(rk), intent(out) :: travec1(3), travec2(3), rotquat(4)
   ! Local variables
   real(rk), allocatable :: weights(:)
   real(rk), allocatable, dimension(:, :) :: coords1, coords2

   weights = mol1%get_weights()
   coords1 = mol1%get_coords()
   coords2 = mol2%get_coords()

   ! Calculate centroids

   travec1 = -centroid(mol1)
   travec2 = -centroid(mol2)

   ! Calculate optimal rotation matrix

   rotquat = leastrotquat( &
      mol1%natom, &
      weights, &
      translated(mol1%natom, coords1, travec1), &
      translated(mol2%natom, coords2, travec2), &
      mapping &
   )

end subroutine

subroutine molecule_align( &
! Purpose: Align atoms0 and atoms1
   mol1, &
   mol2, &
   eltypes, &
   travec1, &
   travec2, &
   rotquat)

   type(molecule_type), intent(in) :: mol1, mol2
   type(partition_type), intent(in) :: eltypes
   real(rk), intent(out) :: travec1(3), travec2(3), rotquat(4)

   ! Abort if molecules have different number of atoms

   if (size(mol1%atoms) /= size(mol2%atoms)) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Abort if molecules are not isomers

!   if (any(sorted(elnums%partition_map1) /= sorted(elnums%partition_map2))) then
!      write (stderr, '(a)') '*Error: These molecules are not isomers'
!      stop
!   end if

   ! Abort if there are conflicting atomic types

   if (any(sorted(eltypes%partition_map1) /= sorted(eltypes%partition_map2))) then
      write (stderr, '(a)') 'Error: There are conflicting atomic types'
      stop
   end if

   ! Abort if atoms are not ordered

   if (any(mol1%atoms%elnum /= mol2%atoms%elnum)) then
!   if (any(elnums%partition_map1 /= elnums%partition_map2)) then
      write (stderr, '(a)') 'Error: The atoms are not in the same order'
      stop
   end if

   ! Abort if atomic types are not ordered

   if (any(eltypes%partition_map1 /= eltypes%partition_map2)) then
      write (stderr, '(a)') 'Error: Atomic types are not in the same order'
      stop
   end if

   ! Calculate centroids

   travec1 = -centroid(mol1)
   travec2 = -centroid(mol2)

   ! Calculate optimal rotation matrix

   rotquat = leastrotquat( &
      mol1%natom, &
      mol1%get_weights(), &
      translated(mol1%natom, mol1%get_coords(), travec1), &
      translated(mol2%natom, mol2%get_coords(), travec2), &
      identity_permutation(mol1%natom) &
   )

end subroutine

function get_rmsd(mol1, mol2, mapping) result(rmsd)
   type(molecule_type), intent(in) :: mol1, mol2
   integer :: mapping(:)
   real(rk) :: rmsd

   rmsd = sqrt(squaredist(mol1%natom, mol1%get_weights(), mol1%get_coords(), &
         mol2%get_coords(), mapping) / sum(mol1%get_weights()))

end function

function get_adjd(mol1, mol2, mapping) result(adjd)
   type(molecule_type), intent(in) :: mol1, mol2
   integer :: mapping(:)
   integer :: adjd

   adjd = adjacencydiff(mol1%natom, mol1%get_adjmatrix(), mol2%get_adjmatrix(), mapping)

end function

function centroid(mol)
! Purpose: Get the centroid coordinates
   type(molecule_type), intent(in) :: mol
   ! Local variables
   integer :: i
   real(rk) :: centroid(3)
   real(rk), allocatable :: coords(:, :)
   real(rk), allocatable :: weights(:)

   coords = mol%get_coords()
   weights = mol%get_weights()

! Calculate the coordinates of the center of mass

   centroid(:) = 0

   do i = 1, size(mol%atoms)
      centroid(:) = centroid(:) + weights(i)*coords(:, i)
   end do

   centroid(:) = centroid(:)/sum(weights)

end function

end module
