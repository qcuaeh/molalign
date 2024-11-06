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
use bipartition
use superpartition
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
subroutine levelup_crossmnatypes(mol1, mol2, mnatypes, submnatypes)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_type), intent(in) :: mnatypes
   type(bipartition_type), intent(out) :: submnatypes
   ! Local variables
   integer :: h, i, index
   integer :: max_num_types, max_num_atoms
   type(neighbordict_type) :: typedict
   type(bipartpointer_type), allocatable :: typelist(:)
   integer, allocatable :: neighborhood(:)

   max_num_types = size(mol1%atoms) + size(mol2%atoms)
   max_num_atoms = max(size(mol1%atoms), size(mol2%atoms))
   call submnatypes%init(max_num_types, max_num_atoms)
   call typedict%init(max_coord_num, mnatypes%largest_part_size)
   allocate (typelist(typedict%num_slots))

   do h = 1, mnatypes%partition_size

      do i = 1, mnatypes%parts(h)%subset1%part_size
         index = mnatypes%parts(h)%subset1%indices(i)
         neighborhood = mnatypes%partition_map1(mol1%atoms(index)%adjlist)
         if (.not. (neighborhood .in. typedict)) then
            typelist(typedict%new_index(neighborhood))%ptr => submnatypes%new_part()
         end if
         call typelist(typedict%get_index(neighborhood))%ptr%subset1%add(index)
      end do

      do i = 1, mnatypes%parts(h)%subset2%part_size
         index = mnatypes%parts(h)%subset2%indices(i)
         neighborhood = mnatypes%partition_map2(mol2%atoms(index)%adjlist)
         if (.not. (neighborhood .in. typedict)) then
            typelist(typedict%new_index(neighborhood))%ptr => submnatypes%new_part()
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
subroutine levelup_mnatypes(mol, mnatypes, submnatypes)
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(in) :: mnatypes
   type(partition_type), intent(out) :: submnatypes
   ! Local variables
   integer :: h, i, index
   type(neighbordict_type) :: typedict
   type(partpointer_type), allocatable :: typelist(:)
   integer, allocatable :: neighborhood(:)

   call submnatypes%init(mnatypes%max_partition_size, mnatypes%largest_part_size)
   call typedict%init(mnatypes%largest_part_size, max_coord_num)
   allocate (typelist(typedict%num_slots))

   do h = 1, mnatypes%partition_size

      do i = 1, mnatypes%parts(h)%part_size
         index = mnatypes%parts(h)%indices(i)
         neighborhood = mnatypes%partition_map(mol%atoms(index)%adjlist)
         if (.not. (neighborhood .in. typedict)) then
            typelist(typedict%new_index(neighborhood))%ptr => submnatypes%new_part()
         end if
         call typelist(typedict%get_index(neighborhood))%ptr%add(index)
      end do

      call typedict%reset()

   end do

end subroutine

! Iteratively compute MNA types
subroutine compute_mnatypes(mol, mnatypes)
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(inout) :: mnatypes
   ! Local variables
   type(partition_type) :: submnatypes

   do

!      write (stderr, *)
!      call mnatypes%print_parts()

      ! Compute next level MNA types
      call levelup_mnatypes(mol, mnatypes, submnatypes)

      ! Exit the loop if types did not change
      if (submnatypes == mnatypes) exit

      ! Update mnatypes
!      call mnatypes%update(submnatypes)
      call mnatypes%reset()
      mnatypes = submnatypes
      call submnatypes%reset()

   end do

end subroutine

! Split MNA types
subroutine split_mnatypes(h, mol, mnatypes, submnatypes)
   integer, intent(in) :: h
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(in) :: mnatypes
   type(partition_type), intent(out) :: submnatypes
   ! Local variables
   integer :: k, i
   type(part_type), pointer :: newtype

   call submnatypes%init(mnatypes%max_partition_size, mnatypes%largest_part_size)

   do k = 1, h - 1
      call submnatypes%add_part(mnatypes%parts(k))
   end do

   do i = 1, mnatypes%parts(h)%part_size
      newtype => submnatypes%new_part()
      call newtype%add(mnatypes%parts(h)%indices(i))
   end do

   do k = h + 1, mnatypes%partition_size
      call submnatypes%add_part(mnatypes%parts(k))
   end do

!   call submnatypes%print_parts()
   call compute_mnatypes(mol, submnatypes)
!   call submnatypes%print_parts()

end subroutine

! Split MNA types
subroutine compute_mnasupertypes(mol, mnatypes)
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(in) :: mnatypes
   ! Local variables
   integer :: h
   type(partition_type) :: submnatypes
   type(superpartition_type) :: mnasupertypes
   type(partitiondict_type) :: partitiondict
   type(partitionpointer_type), allocatable :: partitionlist(:)

   call mnasupertypes%init(mnatypes%partition_size, mnatypes%max_partition_size, mnatypes%largest_part_size)
   call partitiondict%init(mnatypes%partition_size)
   allocate (partitionlist(partitiondict%num_slots))

   do h = 1, mnatypes%partition_size
      call split_mnatypes(h, mol, mnatypes, submnatypes)
      if (.not. (submnatypes .in. partitiondict)) then
!         call submnatypes%print_parts()
         partitionlist(partitiondict%new_index(submnatypes))%ptr => mnasupertypes%new_partition()
      end if
      call partitionlist(partitiondict%get_index(submnatypes))%ptr%add_part(mnatypes%parts(h))
      call submnatypes%reset()
   end do

   do h = 1, mnasupertypes%superpartition_size
      call mnasupertypes%partitions(h)%print_parts()
   end do

end subroutine

end module
