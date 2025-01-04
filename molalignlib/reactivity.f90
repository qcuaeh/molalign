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

module reactivity
use parameters
use sorting
use permutation
use alignment
use partition
use molecule
use tracking
use common_types
use bipartition

implicit none

contains

subroutine remove_reactive_bonds(mol1, mol2, molfrags1, molfrags2, mnatypes, atomperm)
!
! Remove mismatched and reactive bonds
!

   ! Arguments
   type(mol_type), intent(inout) :: mol1, mol2
   type(intlist_type), dimension(:) :: molfrags1, molfrags2
   type(bipartition_type), intent(in) :: mnatypes
   integer, dimension(:), intent(in) :: atomperm

   ! Local variables
!   integer, allocatable, dimension(:) :: adjidcs1, adjidcs2
   integer, allocatable, dimension(:) :: invatomperm
   integer :: iatom, jatom, katom
   integer :: j, k

   ! Initialization

   invatomperm = inverse_perm(atomperm)

   ! Remove mismatched bonds

   do iatom = 1, size(mol1%atoms)
      do j = 1, size(mol1%atoms(iatom)%adjlist)
         jatom = mol1%atoms(iatom)%adjlist(j)
         if (.not. mol2%adjmat(atomperm(iatom), atomperm(jatom))) then
!            write (stderr, *) 'remove mol1 bond:', iatom, jatom
            call mol1%remove_bond(iatom, jatom)
!            adjidcs1 = mnatypes%parts(mnatypes%idcs1(jatom))%items1
!            do k = 1, size(adjidcs1)
!               katom = adjidcs1(k)
!               call mol1%remove_bond(iatom, katom)
!               call mol2%remove_bond(atomperm(iatom), atomperm(katom))
!            end do
         end if
      end do
   end do

   do iatom = 1, size(mol2%atoms)
      do j = 1, size(mol2%atoms(iatom)%adjlist)
         jatom = mol2%atoms(iatom)%adjlist(j)
         if (.not. mol1%adjmat(invatomperm(iatom), invatomperm(jatom))) then
!            write (stderr, *) 'remove mol2 bond:', iatom, jatom
            call mol2%remove_bond(iatom, jatom)
!            adjidcs2 = mnatypes%parts(mnatypes%idcs2(jatom))%items2
!            do k = 1, size(adjidcs2)
!               katom = adjidcs2(k)
!               call mol2%remove_bond(iatom, katom)
!               call mol1%remove_bond(invatomperm(iatom), invatomperm(katom))
!            end do
         end if
      end do
   end do

   ! Dissociate water molecules
!
!   do iatom = 1, size(molfrags1)
!      if (all(sorted(mol1%atoms(molfrags1(iatom)%n)%elnum) == [1, 1, 8])) then
!         do j = 1, size(molfrags1(iatom)%n)
!            jatom = molfrags1(iatom)%n(j)
!            do k = 1, size(mol1%atoms(jatom)%adjlist)
!               katom = mol1%atoms(jatom)%adjlist(k)
!               call mol1%remove_bond(jatom, katom)
!            end do
!         end do
!      end if
!   end do
!
!   do iatom = 1, size(molfrags2)
!      if (all(sorted(mol2%atoms(molfrags2(iatom)%n)%elnum) == [1, 1, 8])) then
!         do j = 1, size(molfrags2(iatom)%n)
!            jatom = molfrags2(iatom)%n(j)
!            do k = 1, size(mol2%atoms(jatom)%adjlist)
!               katom = mol2%atoms(jatom)%adjlist(k)
!               call mol2%remove_bond(jatom, katom)
!            end do
!         end do
!      end if
!   end do

end subroutine

end module
