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
use parameters
use sorting
use chemdata
use molecule
use tupledict
use partition
use metapartition
use partitiondict
use permutation

implicit none

contains

! Partition atoms by atomic number and label
subroutine compute_eltypes(mol, eltypes)
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(out) :: eltypes
   ! Local variables
   integer :: eltype(2)
   type(tupledict_type) :: typedict
   type(partpointer_type), allocatable :: typelist(:)
   integer :: i, num_atoms

   num_atoms = size(mol%atoms)
   call eltypes%initialize(num_atoms)
   call typedict%initialize(num_atoms, 'ordered')
   allocate (typelist(typedict%num_slots))

   do i = 1, num_atoms
      eltype(1) = mol%atoms(i)%elnum
      eltype(2) = mol%atoms(i)%label
      if (.not. (eltype .in. typedict)) then
         typelist(typedict%new_index(eltype))%ptr => &
            eltypes%new_part(num_atoms)
      end if
      call typelist(typedict%get_index(eltype))%ptr%add(i)
   end do

end subroutine

! Level up MNA types
subroutine levelup_mnatypes(mol, mnatypes, subtypes)
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(in) :: mnatypes
   type(partition_type), intent(out) :: subtypes
   ! Local variables
   integer :: h, i, iatom
   type(tupledict_type) :: typedict
   type(partpointer_type), allocatable :: typelist(:)
   integer, allocatable :: neighborhood(:)

   call subtypes%initialize(mnatypes%num_items)
   call typedict%initialize(mnatypes%largest_part_size, 'unordered')
   allocate (typelist(typedict%num_slots))

   do h = 1, mnatypes%num_parts

      do i = 1, mnatypes%parts(h)%part_size
         iatom = mnatypes%parts(h)%items(i)
         neighborhood = mnatypes%idcs(mol%atoms(iatom)%adjlist)
         if (.not. (neighborhood .in. typedict)) then
            typelist(typedict%new_index(neighborhood))%ptr => &
               subtypes%new_part(mnatypes%parts(h)%part_size)
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

end module
