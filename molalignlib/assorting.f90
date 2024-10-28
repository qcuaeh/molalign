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
subroutine compute_eltypes(mol1, mol2, eltypes)
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

!   call eltypes%print_parts()

end subroutine

! Level up MNA types
subroutine levelup_mnatypes(atoms1, atoms2, mnatypes, subtypes)
   type(atom_type), dimension(:), intent(in) :: atoms1, atoms2
   type(partition_type), intent(in) :: mnatypes
   type(partition_type), intent(out) :: subtypes
   ! Local variables
   integer :: h, i, iatom
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
         iatom = mnatypes%parts(h)%subset1%indices(i)
         neighborhood = mnatypes%partition_map1(atoms1(iatom)%adjlist)
         if (.not. subtypedict%has_index(neighborhood)) then
            partlist(subtypedict%new_index(neighborhood))%ptr => subtypes%new_part()
         end if
         call partlist(subtypedict%get_index(neighborhood))%ptr%subset1%add(iatom)
      end do

      do i = 1, mnatypes%parts(h)%subset2%part_size
         iatom = mnatypes%parts(h)%subset2%indices(i)
         neighborhood = mnatypes%partition_map2(atoms2(iatom)%adjlist)
         if (.not. subtypedict%has_index(neighborhood)) then
            partlist(subtypedict%new_index(neighborhood))%ptr => subtypes%new_part()
         end if
         call partlist(subtypedict%get_index(neighborhood))%ptr%subset2%add(iatom)
      end do

      call subtypedict%reset()

   end do

end subroutine

end module
