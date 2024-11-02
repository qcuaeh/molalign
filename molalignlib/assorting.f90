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
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_type), intent(out) :: eltypes
   ! Local variables
   integer :: i, elnum, max_num_atoms
   type(bipartpointer_type), allocatable :: typelist(:)

   allocate (typelist(nelem))

   max_num_atoms = max(size(mol1%atoms), size(mol2%atoms))
   call eltypes%initialize(nelem, max_num_atoms)

   do i = 1, size(mol1%atoms)
      elnum = mol1%atoms(i)%elnum
      if (.not. associated(typelist(elnum)%ptr)) then
         typelist(elnum)%ptr => eltypes%new_part()
      end if
      call typelist(elnum)%ptr%subset1%add(i)
   end do

   do i = 1, size(mol2%atoms)
      elnum = mol2%atoms(i)%elnum
      if (.not. associated(typelist(elnum)%ptr)) then
         typelist(elnum)%ptr => eltypes%new_part()
      end if
      call typelist(elnum)%ptr%subset2%add(i)
   end do

end subroutine

! Level up MNA types
subroutine levelup_crossmnatypes(mol1, mol2, mnatypes, subtypes)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_type), intent(in) :: mnatypes
   type(bipartition_type), intent(out) :: subtypes
   ! Local variables
   integer :: h, i, index
   integer :: max_num_types, max_num_atoms
   type(dict_type) :: typedict
   type(bipartpointer_type), allocatable :: typelist(:)
   integer, allocatable :: neighborhood(:)

   max_num_types = size(mol1%atoms) + size(mol2%atoms)
   max_num_atoms = max(size(mol1%atoms), size(mol2%atoms))
   call subtypes%initialize(max_num_types, max_num_atoms)
   call typedict%init(maxcoord, mnatypes%partition_size, mnatypes%largest_part_size)
   allocate (typelist(typedict%num_slots))

   do h = 1, mnatypes%partition_size

      do i = 1, mnatypes%parts(h)%subset1%part_size
         index = mnatypes%parts(h)%subset1%indices(i)
         neighborhood = mnatypes%partition_map1(mol1%atoms(index)%adjlist)
         if (.not. typedict%has_index(neighborhood)) then
            typelist(typedict%new_index(neighborhood))%ptr => subtypes%new_part()
         end if
         call typelist(typedict%get_index(neighborhood))%ptr%subset1%add(index)
      end do

      do i = 1, mnatypes%parts(h)%subset2%part_size
         index = mnatypes%parts(h)%subset2%indices(i)
         neighborhood = mnatypes%partition_map2(mol2%atoms(index)%adjlist)
         if (.not. typedict%has_index(neighborhood)) then
            typelist(typedict%new_index(neighborhood))%ptr => subtypes%new_part()
         end if
         call typelist(typedict%get_index(neighborhood))%ptr%subset2%add(index)
      end do

      call typedict%reset()

   end do

end subroutine

! Partition atoms by atomic number
subroutine compute_eltypes(mol, eltypes)
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(out) :: eltypes
   ! Local variables
   integer :: i, elnum
   type(partpointer_type), allocatable :: typelist(:)

   allocate (typelist(nelem))

   call eltypes%initialize(nelem, size(mol%atoms))

   do i = 1, size(mol%atoms)
      elnum = mol%atoms(i)%elnum
      if (.not. associated(typelist(elnum)%ptr)) then
         typelist(elnum)%ptr => eltypes%new_part()
      end if
      call typelist(elnum)%ptr%add(i)
   end do

end subroutine

! Level up MNA types
subroutine levelup_mnatypes(mol, mnatypes, subtypes)
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(in) :: mnatypes
   type(partition_type), intent(out) :: subtypes
   ! Local variables
   integer :: h, i, index
   integer:: max_part_size, max_partition_size
   type(dict_type) :: typedict
   type(partpointer_type), allocatable :: typelist(:)
   integer, allocatable :: neighborhood(:)

   max_part_size = size(mol%atoms)
   max_partition_size = size(mol%atoms)
   call subtypes%initialize(max_partition_size, max_part_size)
   call typedict%init(maxcoord, mnatypes%partition_size, mnatypes%largest_part_size)
   allocate (typelist(typedict%num_slots))

   do h = 1, mnatypes%partition_size

      do i = 1, mnatypes%parts(h)%part_size
         index = mnatypes%parts(h)%indices(i)
         neighborhood = mnatypes%partition_map(mol%atoms(index)%adjlist)
         if (.not. typedict%has_index(neighborhood)) then
            typelist(typedict%new_index(neighborhood))%ptr => subtypes%new_part()
         end if
         call typelist(typedict%get_index(neighborhood))%ptr%add(index)
      end do

      call typedict%reset()

   end do

end subroutine

! Iteratively compute MNA types
subroutine compute_mnatypes(mol, eltypes, mnatypes)
   type(mol_type), intent(in) :: mol
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

   call subtypes%initialize(mnatypes%max_partition_size, mnatypes%max_part_size)

   do h = 1, mnatypes%partition_size

      if (h /= h0) then
         newtype%ptr => subtypes%new_part()
         do i = 1, mnatypes%parts(h)%part_size
            index = mnatypes%parts(h)%indices(i)
            call newtype%ptr%add(index)
         end do
      end if

      do i = 1, mnatypes%parts(h0)%part_size
         index = mnatypes%parts(h0)%indices(i)
         newtype%ptr => subtypes%new_part()
         call newtype%ptr%add(index)
      end do

   end do

end subroutine

end module
