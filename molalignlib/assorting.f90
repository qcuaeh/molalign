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
   call eltypes%init(nelem, max_num_atoms)

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
   call subtypes%init(max_num_types, max_num_atoms)
   call typedict%init(max_coord_num, mnatypes%largest_part_size)
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
   call eltypes%init(nelem, size(mol%atoms))

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
   type(dict_type) :: typedict
   type(partpointer_type), allocatable :: typelist(:)
   integer, allocatable :: neighborhood(:)

   call subtypes%init(mnatypes%max_partition_size, mnatypes%largest_part_size)
   call typedict%init(mnatypes%largest_part_size, max_coord_num)
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

   ! Level 0 mnatypes are eltypes
   mnatypes = eltypes

   do

!      write (stderr, *)
!      call mnatypes%print_parts()

      ! Compute next last_level MNA subtypes
      call levelup_mnatypes(mol, mnatypes, subtypes)

      ! Exit the loop if subtypes are unchanged
      if (subtypes == mnatypes) exit

      ! Update mnatypes
!      call mnatypes%update(subtypes)
      call mnatypes%reset()
      mnatypes = subtypes
      call subtypes%reset()

   end do

end subroutine

! Split MNA types
subroutine split_type(h, types, subtypes)
   integer, intent(in) :: h
   type(partition_type), intent(in) :: types
   type(partition_type), intent(inout) :: subtypes
   ! Local variables
   integer :: k, i
   type(part_type), pointer :: newtype

   call subtypes%init(types%max_partition_size, types%largest_part_size)

   do k = 1, h - 1
      call subtypes%add_part(types%parts(k))
   end do

   do i = 1, types%parts(h)%part_size
      newtype => subtypes%new_part()
      call newtype%add(types%parts(h)%indices(i))
   end do

   do k = h + 1, types%partition_size
      call subtypes%add_part(types%parts(k))
   end do

end subroutine

! Split MNA types
subroutine compute_submnatypes(mol, mnatypes)
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(in) :: mnatypes
   ! Local variables
   integer :: h
   type(partition_type) :: subtypes, submnatypes
!   type(partition_type), allocatable :: uniqsubmnatypelist

!   allocate (uniqsubmnatypelist(mnatypes%partition_size))

   do h = 1, mnatypes%partition_size
      call split_type(h, mnatypes, subtypes)
!      write (stderr, *)
!      call subtypes%print_parts()
      call compute_mnatypes(mol, subtypes, submnatypes)
!      if (.not. submnatypes .in. uniqsubmnatypelist(k)) then
!         uniqsubmnatypelist.append(submnatypes)
!      end if
      write (stderr, *)
      call submnatypes%print_parts()
      call submnatypes%reset()
      call subtypes%reset()
   end do

end subroutine

end module
