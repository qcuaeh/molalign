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

module partitioning
use stdio
use kinds
use flags
use bounds
use sorting
use chemdata
use molecule
use neighbordict
use partition
use metapartition
use partitiondict
use permutation

implicit none

contains

! Partition atoms by atomic number
subroutine compute_eltypes(mol, eltypes)
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(out) :: eltypes
   ! Local variables
   integer :: i, elnum, atoms_size
   type(partpointer_type), allocatable :: typelist(:)

   atoms_size = size(mol%atoms)
   call eltypes%initialize(atoms_size)
   allocate (typelist(nelem))

   do i = 1, atoms_size
      elnum = mol%atoms(i)%elnum
      if (.not. associated(typelist(elnum)%ptr)) then
         typelist(elnum)%ptr => eltypes%new_part(atoms_size)
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
   integer :: h, i, iatom
   type(neighbordict_type) :: typedict
   type(partpointer_type), allocatable :: typelist(:)
   integer, allocatable :: neighborhood(:)

   call subtypes%initialize(mnatypes%num_items)
   call typedict%initialize(mnatypes%largest_part_size)
   allocate (typelist(typedict%num_slots))

   do h = 1, mnatypes%num_parts

      do i = 1, mnatypes%parts(h)%size
         iatom = mnatypes%parts(h)%list(i)
         neighborhood = mnatypes%indices(mol%atoms(iatom)%adjlist)
         if (.not. (neighborhood .in. typedict)) then
            typelist(typedict%new_index(neighborhood))%ptr => &
               subtypes%new_part(mnatypes%parts(h)%size)
         end if
         call typelist(typedict%get_index(neighborhood))%ptr%add(iatom)
      end do

      call typedict%reset()

   end do

end subroutine

! Iteratively compute MNA types
subroutine compute_mnatypes(mol, mnatypes)
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(inout) :: mnatypes
   ! Local variables
   type(partition_type) :: subtypes

   do

!      write (stderr, *)
!      call mnatypes%print_parts()

      ! Compute MNA upper level types
      call levelup_mnatypes(mol, mnatypes, subtypes)

      ! Exit loop if types did not change
      if (subtypes == mnatypes) then
         mnatypes = subtypes
         exit
      end if

      ! Update mnatypes
      mnatypes = subtypes

   end do

end subroutine

! Split MNA types
subroutine split_mnatypes(h, mol, mnatypes)
   integer, intent(in) :: h
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(inout) :: mnatypes
   ! Local variables
   integer :: k, i
   type(part_type), pointer :: newtype
   type(partition_type) :: subtypes

   call subtypes%initialize(mnatypes%num_items)

   do k = 1, h - 1
      call subtypes%add_part(mnatypes%parts(k))
   end do

   do i = 1, mnatypes%parts(h)%size
      newtype => subtypes%new_part(mnatypes%parts(h)%size)
      call newtype%add(mnatypes%parts(h)%list(i))
   end do

   do k = h + 1, mnatypes%num_parts
      call subtypes%add_part(mnatypes%parts(k))
   end do

   mnatypes = subtypes

end subroutine

! Split MNA types
subroutine collect_degenerated_mnatypes(mol, mnatypes, metatypes)
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(in) :: mnatypes
   type(metapartition_type), intent(out) :: metatypes
   ! Local variables
   integer :: h, k
   type(partition_type) :: subtypes
   type(partitiondict_type) :: partitiondict
   type(metapartition_type) :: mnametatypes
   type(metapartpointer_type), allocatable :: partlist(:)
   type(partition_type), allocatable :: subpartitionlist(:)
   logical, allocatable :: primality(:)

   call mnametatypes%initialize(mnatypes%num_parts)
   call partitiondict%initialize(mnatypes%num_parts)
   allocate (partlist(partitiondict%num_slots))
   allocate (subpartitionlist(mnatypes%num_parts))

   do h = 1, mnatypes%num_parts
      if (mnatypes%parts(h)%size > 1) then
         subtypes = mnatypes
         call split_mnatypes(h, mol, subtypes)
         call compute_mnatypes(mol, subtypes)
         if (.not. (subtypes .in. partitiondict)) then
!            call subtypes%print_parts()
            partlist(partitiondict%new_index(subtypes))%ptr => &
               mnametatypes%new_part(mnatypes%num_parts)
            subpartitionlist(mnametatypes%num_parts) = subtypes
         end if
         call partlist(partitiondict%get_index(subtypes))%ptr%add(h)
      end if
   end do

   allocate (primality(mnametatypes%num_parts))
   primality(:) = .true.

   do h = 1, mnametatypes%num_parts
      if (primality(h)) then
         do k = h + 1, mnametatypes%num_parts
            if (primality(k)) then
               if (subpartitionlist(h) < subpartitionlist(k)) then
                  primality(k) = .false.
               else if (subpartitionlist(k) < subpartitionlist(h)) then
                  primality(h) = .false.
               end if
            end if
         end do
      end if
   end do

   call metatypes%initialize(count(primality))

   do h = 1, mnametatypes%num_parts
      if (primality(h)) then
         call metatypes%add_part(mnametatypes%parts(h))
      end if
   end do

end subroutine

end module
