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
use kinds
use sorting
use permutation
use alignment
use partition
use molecule
use tracking

implicit none

contains

subroutine remove_reactive_bonds(mol1, mol2, mapping)

   type(molecule_type), intent(inout) :: mol1, mol2
   integer, dimension(:), intent(in) :: mapping

   integer :: i, j, k, j_, k_
   type(partition_type) :: mnatypes1, mnatypes2
   integer, allocatable, dimension(:) :: elnums1, elnums2
   integer, allocatable, dimension(:) :: mnatypemap1, mnatypemap2
   type(atomlist_type), allocatable, dimension(:) :: adjlists1, adjlists2
   type(atomlist_type), allocatable, dimension(:) :: molfragparts1, molfragparts2
   logical, allocatable, dimension(:, :) :: adjmat1, adjmat2
   integer, allocatable, dimension(:) :: mnatypepartidcs1, mnatypepartidcs2
   integer, allocatable, dimension(:) :: unmapping
   real(rk) :: rotquat(4)

   ! Align coordinates

   rotquat = leastrotquat(mol1%natom, mol1%get_weights(), mol1%get_coords(), mol2%get_coords(), mapping)
   call mol2%rotate_coords(rotquat)

!   write (stderr, *) sqrt(squaredist(mol1%natom, mol1%get_weights(), mol1%get_coords(), mol2%get_coords(), mapping) &
!         /sum(mol1%get_weights()))

   ! Initialization

   adjmat1 = mol1%get_adjmatrix()
   adjmat2 = mol2%get_adjmatrix()
   elnums1 = mol1%get_elnums()
   elnums2 = mol2%get_elnums()
   adjlists1 = mol1%get_adjlists()
   adjlists2 = mol2%get_adjlists()
   mnatypes1 = mol1%mnatypes
   mnatypes2 = mol2%mnatypes
   mnatypemap1 = mol1%mnatypes%partition_map
   mnatypemap2 = mol2%mnatypes%partition_map
   molfragparts1 = mol1%get_molfrags()
   molfragparts2 = mol2%get_molfrags()
   unmapping = inverse_permutation(mapping)

   ! Remove mismatched bonds

   do i = 1, size(mol1%atoms)
      do j_ = 1, size(adjlists1(i)%atomidcs)
         j = adjlists1(i)%atomidcs(j_)
         if (.not. adjmat2(mapping(i), mapping(j))) then
            mnatypepartidcs1 = mnatypes1%parts(mnatypemap1(j))%indices
            do k_ = 1, size(mnatypepartidcs1)
               k = mnatypepartidcs1(k_)
               call mol1%remove_bond(i, k)
               call mol2%remove_bond(mapping(i), mapping(k))
            end do
         end if
      end do
   end do

   do i = 1, size(mol2%atoms)
      do j_ = 1, size(adjlists2(i)%atomidcs)
         j = adjlists2(i)%atomidcs(j_)
         if (.not. adjmat1(unmapping(i), unmapping(j))) then
            mnatypepartidcs2 = mnatypes2%parts(mnatypemap2(j))%indices
            do k_ = 1, size(mnatypepartidcs2)
               k = mnatypepartidcs2(k_)
               call mol2%remove_bond(i, k)
               call mol1%remove_bond(unmapping(i), unmapping(k))
            end do
         end if
      end do
   end do

   ! Remove water bonds

   do i = 1, size(molfragparts1)
      if (all(sorted(elnums1(molfragparts1(i)%atomidcs)) == [1, 1, 8])) then
         do j_ = 1, size(molfragparts1(i)%atomidcs)
            j = molfragparts1(i)%atomidcs(j_)
            do k_ = 1, size(adjlists1(j)%atomidcs)
               k = adjlists1(j)%atomidcs(k_)
               call mol1%remove_bond(j, k)
            end do
         end do
      end if
   end do

   do i = 1, size(molfragparts2)
      if (all(sorted(elnums2(molfragparts2(i)%atomidcs)) == [1, 1, 8])) then
         do j_ = 1, size(molfragparts2(i)%atomidcs)
            j = molfragparts2(i)%atomidcs(j_)
            do k_ = 1, size(adjlists2(j)%atomidcs)
               k = adjlists2(j)%atomidcs(k_)
               call mol2%remove_bond(j, k)
            end do
         end do
      end if
   end do

   call find_molfrags(mol1)
   call find_molfrags(mol2)

end subroutine

end module
