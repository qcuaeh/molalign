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
use bounds
use random
use linalg
use sorting
use molecule
use rotation
use translation
use permutation
use bipartition
use bipartitioning
use adjacency
use alignment
use assignment
!use reactivity

implicit none

contains

! Assign atoms0 and atoms1
subroutine molecule_remap( &
   mol1, &
   mol2, &
   nrec, &
   maplist, &
   countlist)

   type(mol_type), intent(inout) :: mol1, mol2
   integer, intent(out) :: nrec
   integer, dimension(:, :), intent(inout) :: maplist
   integer, dimension(:), intent(inout) :: countlist

   type(bipartition_type) :: eltypes
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

   ! Compute atomic types

   call compute_crosseltypes(mol1, mol2, eltypes)
!   call eltypes%print_parts()

   ! Abort if there are conflicting atomic types

   if (any(sorted(eltypes%indices1) /= sorted(eltypes%indices2))) then
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
   atomperm, &
   travec1, &
   travec2, &
   rotquat)

   type(mol_type), intent(in) :: mol1, mol2
   integer, intent(in) :: atomperm(:)
   real(rk), intent(out) :: travec1(3), travec2(3), rotquat(4)
   ! Local variables
   real(rk), allocatable, dimension(:, :) :: coords1, coords2

   coords1 = mol1%get_coords()
   coords2 = mol2%get_coords()

   ! Calculate centroids

   travec1 = -centroid(mol1)
   travec2 = -centroid(mol2)

   ! Calculate optimal rotation matrix

   rotquat = leastrotquat( &
      mol1%natom, &
      weights(mol1%atoms%elnum), &
      translated(mol1%natom, coords1, travec1), &
      translated(mol2%natom, coords2, travec2), &
      atomperm &
   )

end subroutine

subroutine molecule_align( &
! Purpose: Align atoms0 and atoms1
   mol1, &
   mol2, &
   travec1, &
   travec2, &
   rotquat)

   type(mol_type), intent(in) :: mol1, mol2
   real(rk), intent(out) :: travec1(3), travec2(3), rotquat(4)
   type(bipartition_type) :: eltypes

   ! Abort if molecules have different number of atoms

   if (size(mol1%atoms) /= size(mol2%atoms)) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Abort if molecules are not isomers

!   if (any(sorted(elnums%indices1) /= sorted(elnums%indices2))) then
!      write (stderr, '(a)') '*Error: These molecules are not isomers'
!      stop
!   end if

   ! Compute atomic types

   call compute_crosseltypes(mol1, mol2, eltypes)

   ! Abort if there are conflicting atomic types

   if (any(sorted(eltypes%indices1) /= sorted(eltypes%indices2))) then
      write (stderr, '(a)') 'Error: There are conflicting atomic types'
      stop
   end if

   ! Abort if atoms are not ordered

   if (any(mol1%atoms%elnum /= mol2%atoms%elnum)) then
!   if (any(elnums%indices1 /= elnums%indices2)) then
      write (stderr, '(a)') 'Error: The atoms are not in the same order'
      stop
   end if

   ! Abort if atomic types are not ordered

   if (any(eltypes%indices1 /= eltypes%indices2)) then
      write (stderr, '(a)') 'Error: Atomic types are not in the same order'
      stop
   end if

   ! Calculate centroids

   travec1 = -centroid(mol1)
   travec2 = -centroid(mol2)

   ! Calculate optimal rotation matrix

   rotquat = leastrotquat( &
      mol1%natom, &
      weights(mol1%atoms%elnum), &
      translated(mol1%natom, mol1%get_coords(), travec1), &
      translated(mol2%natom, mol2%get_coords(), travec2), &
      identity_perm(mol1%natom) &
   )

end subroutine

function get_rmsd(mol1, mol2, atomperm) result(rmsd)
   type(mol_type), intent(in) :: mol1, mol2
   integer :: atomperm(:)
   real(rk) :: rmsd

   rmsd = rmsdist(mol1%natom, weights(mol1%atoms%elnum), mol1%get_coords(), &
         mol2%get_coords(), atomperm)

end function

function get_adjd(mol1, mol2, atomperm) result(adjd)
   type(mol_type), intent(in) :: mol1, mol2
   integer :: atomperm(:)
   integer :: adjd

   adjd = adjdiff(mol1%natom, mol1%adjmat, mol2%adjmat, atomperm)

end function

function centroid(mol)
! Purpose: Get the centroid coordinates
   type(mol_type), intent(in) :: mol
   ! Local variables
   integer :: i
   real(rk) :: centroid(3)
   real(rk), allocatable :: coords(:, :)

   coords = mol%get_coords()

! Calculate the coordinates of the center of mass

   centroid(:) = 0

   do i = 1, size(mol%atoms)
      centroid(:) = centroid(:) + weights(mol%atoms(i)%elnum)*coords(:, i)
   end do

   centroid(:) = centroid(:)/sum(weights(mol%atoms%elnum))

end function

end module
