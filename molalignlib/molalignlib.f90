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
use rigid_body
use permutation
use bipartition
use bipartitioning
use adjacency
use alignment
use atom_mapping
!use reactivity

implicit none

contains

! Assign atoms0 and atoms1
subroutine molecule_remap( &
   mol1, &
   mol2, &
   nrec, &
   permlist, &
   countlist)

   type(mol_type), intent(inout) :: mol1, mol2
   integer, intent(out) :: nrec
   integer, dimension(:, :), intent(inout) :: permlist
   integer, dimension(:), intent(inout) :: countlist

   type(bipartition_type) :: eltypes
   real(rk) :: travec1(3), travec2(3)

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

   ! Initialize random number generator

   call random_initialize()

   ! Optimize assignment to minimize the AdjD and RMSD

   call map_atoms(mol1, mol2, eltypes, permlist, countlist, nrec)

   ! Remove bonds from reactive sites and reoptimize assignment

!   if (reac_flag) then
!      call remove_reactive_bonds(mol1, mol2, permlist(:, 1))
!      call map_atoms(mol1, mol2, permlist, countlist, nrec)
!   end if

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
   real(rk), dimension(:), allocatable :: weights1, weights2

   coords1 = mol1%get_coords()
   coords2 = mol2%get_coords()
   weights1 = element_weights(mol1%atoms%elnum)
   weights2 = element_weights(mol1%atoms%elnum)

   ! Calculate centroids

   travec1 = -centroid(coords1, weights1)
   travec2 = -centroid(coords2, weights2)


   ! Calculate optimal rotation matrix

   rotquat = leastrotquat( &
      size(mol1%atoms), &
      element_weights(mol1%atoms%elnum), &
      translated_coords(coords1, travec1), &
      translated_coords(coords2, travec2), &
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
   ! Local variables
   integer :: natom1
   real(rk), allocatable, dimension(:, :) :: coords1, coords2
   real(rk), dimension(:), allocatable :: weights1, weights2

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

   ! Abort if there are conflicting atomic types

   if (any(sorted(eltypes%indices1) /= sorted(eltypes%indices2))) then
      write (stderr, '(a)') 'Error: There are conflicting atomic types'
      stop
   end if

   ! Abort if atoms are not ordered

   if (any(mol1%atoms%elnum /= mol2%atoms%elnum)) then
!   if (any(atomic_numbers%indices1 /= atomic_numbers%indices2)) then
      write (stderr, '(a)') 'Error: The atoms are not in the same order'
      stop
   end if

   ! Abort if atomic types are not ordered

   if (any(eltypes%indices1 /= eltypes%indices2)) then
      write (stderr, '(a)') 'Error: Atomic types are not in the same order'
      stop
   end if

   natom1 = size(mol1%atoms)
   coords1 = mol1%get_coords()
   coords2 = mol2%get_coords()
   weights1 = element_weights(mol1%atoms%elnum)
   weights2 = element_weights(mol1%atoms%elnum)

   ! Calculate centroids

   travec1 = -centroid(coords1, weights1)
   travec2 = -centroid(coords2, weights2)

   ! Calculate optimal rotation matrix

   rotquat = leastrotquat( &
      natom1, &
      weights1, &
      translated_coords(coords1, travec1), &
      translated_coords(coords2, travec2), &
      identity_perm(natom1) &
   )

end subroutine

end module
