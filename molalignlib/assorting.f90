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

module assorting
use stdio
use kinds
use hash_table
use partition
use permutation
use molecule
use flags
use bounds
use sorting
use chemdata

implicit none

contains

! Partition atoms by atomic number
subroutine compute_crosseltypes(mol1, mol2, eltypes)
   type(molecule_type), intent(in) :: mol1, mol2
   type(partition_type), intent(out) :: eltypes
   ! Local variables
   integer :: i, elnum, max_num_atoms
   type(partpointer_type), allocatable :: partlist(:)

   allocate (partlist(nelem))

   max_num_atoms = max(size(mol1%atoms), size(mol2%atoms))
   call eltypes%init(nelem, max_num_atoms)

   do i = 1, size(mol1%atoms)
      elnum = mol1%atoms(i)%elnum
      if (.not. associated(partlist(elnum)%ptr)) then
         partlist(elnum)%ptr => eltypes%new_part()
      end if
      call partlist(elnum)%ptr%subset1%add(i)
   end do

   do i = 1, size(mol2%atoms)
      elnum = mol2%atoms(i)%elnum
      if (.not. associated(partlist(elnum)%ptr)) then
         partlist(elnum)%ptr => eltypes%new_part()
      end if
      call partlist(elnum)%ptr%subset2%add(i)
   end do

end subroutine

! Level up MNA types
subroutine levelup_crossmnatypes(atoms1, atoms2, mnatypes, subtypes)
   type(atom_type), dimension(:), intent(in) :: atoms1, atoms2
   type(partition_type), intent(in) :: mnatypes
   type(partition_type), intent(out) :: subtypes
   ! Local variables
   integer :: h, i, index
   integer :: max_num_types, max_num_atoms
   type(dict_type) :: subtypedict
   type(partpointer_type), allocatable :: partlist(:)
   integer, allocatable :: neighborhood(:)

   max_num_types = size(atoms1) + size(atoms2)
   max_num_atoms = max(size(atoms1), size(atoms2))
   call subtypes%init(max_num_types, max_num_atoms)
   call subtypedict%init(maxcoord, mnatypes%partition_size, mnatypes%largest_part_size)
   allocate (partlist(subtypedict%num_slots))

   do h = 1, mnatypes%partition_size

      do i = 1, mnatypes%parts(h)%subset1%part_size
         index = mnatypes%parts(h)%subset1%indices(i)
         neighborhood = mnatypes%partition_map1(atoms1(index)%adjlist)
         if (.not. subtypedict%has_index(neighborhood)) then
            partlist(subtypedict%new_index(neighborhood))%ptr => subtypes%new_part()
         end if
         call partlist(subtypedict%get_index(neighborhood))%ptr%subset1%add(index)
      end do

      do i = 1, mnatypes%parts(h)%subset2%part_size
         index = mnatypes%parts(h)%subset2%indices(i)
         neighborhood = mnatypes%partition_map2(atoms2(index)%adjlist)
         if (.not. subtypedict%has_index(neighborhood)) then
            partlist(subtypedict%new_index(neighborhood))%ptr => subtypes%new_part()
         end if
         call partlist(subtypedict%get_index(neighborhood))%ptr%subset2%add(index)
      end do

      call subtypedict%reset()

   end do

end subroutine

! Partition atoms by atomic number
subroutine compute_eltypes(mol, eltypes)
   type(molecule_type), intent(in) :: mol
   type(partition_type), intent(out) :: eltypes
   ! Local variables
   integer :: i, elnum
   type(partpointer_type), allocatable :: partlist(:)

   allocate (partlist(nelem))

   call eltypes%init(nelem, size(mol%atoms))

   do i = 1, size(mol%atoms)
      elnum = mol%atoms(i)%elnum
      if (.not. associated(partlist(elnum)%ptr)) then
         partlist(elnum)%ptr => eltypes%new_part()
      end if
      call partlist(elnum)%ptr%subset1%add(i)
   end do

end subroutine

! Level up MNA types
subroutine levelup_mnatypes(mol, mnatypes, subtypes)
   type(molecule_type), intent(in) :: mol
   type(partition_type), intent(in) :: mnatypes
   type(partition_type), intent(out) :: subtypes
   ! Local variables
   integer :: h, i, index
   integer:: max_part_size, max_partition_size
   type(dict_type) :: subtypedict
   type(partpointer_type), allocatable :: partlist(:)
   integer, allocatable :: neighborhood(:)

   max_part_size = size(mol%atoms)
   max_partition_size = size(mol%atoms)
   call subtypes%init(max_partition_size, max_part_size)
   call subtypedict%init(maxcoord, mnatypes%partition_size, mnatypes%largest_part_size)
   allocate (partlist(subtypedict%num_slots))

   do h = 1, mnatypes%partition_size

      do i = 1, mnatypes%parts(h)%subset1%part_size
         index = mnatypes%parts(h)%subset1%indices(i)
         neighborhood = mnatypes%partition_map1(mol%atoms(index)%adjlist)
         if (.not. subtypedict%has_index(neighborhood)) then
            partlist(subtypedict%new_index(neighborhood))%ptr => subtypes%new_part()
         end if
         call partlist(subtypedict%get_index(neighborhood))%ptr%subset1%add(index)
      end do

      call subtypedict%reset()

   end do

end subroutine

! Iteratively compute MNA types
subroutine compute_mnatypes(mol, eltypes, mnatypes)
   type(molecule_type), intent(in) :: mol
   type(partition_type), intent(in) :: eltypes
   type(partition_type), intent(out) :: mnatypes
   ! Local variables
   type(partition_type) :: subtypes

   ! eltypes are our initial mnatypes
   mnatypes = eltypes

   do

      ! Compute next last_level MNA subtypes
      call levelup_mnatypes(mol, mnatypes, subtypes)

      ! Exit the loop if subtypes are unchanged
      if (subtypes == mnatypes) exit

      ! Update mnatypes
      mnatypes = subtypes

   end do

end subroutine

! Split MNA types
subroutine splitup_mnatypes(h0, mnatypes, subtypes)
   integer, intent(in) :: h0
   type(partition_type), intent(in) :: mnatypes
   type(partition_type), intent(out) :: subtypes
   ! Local variables
   integer :: h, i, index
   type(partpointer_type) :: newtype

   call subtypes%init(mnatypes%max_partition_size, mnatypes%max_part_size)

   do h = 1, mnatypes%partition_size

      if (h /= h0) then
         newtype%ptr => subtypes%new_part()
         do i = 1, mnatypes%parts(h)%subset1%part_size
            index = mnatypes%parts(h)%subset1%indices(i)
            call newtype%ptr%subset1%add(index)
         end do
      end if

      do i = 1, mnatypes%parts(h0)%subset1%part_size
         index = mnatypes%parts(h0)%subset1%indices(i)
         newtype%ptr => subtypes%new_part()
         call newtype%ptr%subset1%add(index)
      end do

   end do

end subroutine

end module
