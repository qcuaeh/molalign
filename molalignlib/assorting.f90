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
   type(partpointer_type), allocatable :: typelist(:)

   allocate (typelist(nelem))

   max_num_atoms = max(size(mol1%atoms), size(mol2%atoms))
   call eltypes%init(nelem, max_num_atoms)

   do i = 1, size(mol1%atoms)
      elnum = mol1%atoms(i)%elnum
      if (.not. associated(typelist(elnum)%ptr)) then
         typelist(elnum)%ptr => eltypes%get_new_part()
      end if
      call typelist(elnum)%ptr%subset1%add(i)
   end do

   do i = 1, size(mol2%atoms)
      elnum = mol2%atoms(i)%elnum
      if (.not. associated(typelist(elnum)%ptr)) then
         typelist(elnum)%ptr => eltypes%get_new_part()
      end if
      call typelist(elnum)%ptr%subset2%add(i)
   end do

!   call eltypes%print_parts()

end subroutine

!! Level up MNA types
!subroutine levelup_mnatypes(adjlists1, adjlists2, mnatypes, subtypes)
!   type(list_type), dimension(:), intent(in) :: adjlists1, adjlists2
!   type(partition_type), intent(in) :: mnatypes
!   type(partition_type), intent(out) :: subtypes
!   ! Local variables
!   integer :: h, i
!   integer :: max_num_types
!   type(dict_type) :: subtypedict
!   type(partpointer_type), allocatable :: subtypelist(:)
!   integer, allocatable :: neighborhood(:)
!
!   max_num_types = size(adjlists1) + size(adjlists2)
!   max_num_atoms = max(size(adjlists1), size(adjlists2))
!   call subtypes%init(max_num_types, max_num_atoms)
!
!   call subtypedict%init(maxcoord, mnatypes%partition_size, mnatypes%largest_subset_size)
!   allocate (subtypelist(subtypedict%num_slots))
!
!   do h = 1, mnatypes%partition_size
!      do i = 1, mnatypes%parts(h)%subset1%part_size
!         iatom = mnatypes%parts(h)%subset1%indices(i)
!         neighborhood = mnatypes%partition_map1%s(adjlists1(iatom)%n)
!         if (.not. subtypedict%has_index(neighborhood)) then
!            subtypelist(subtypedict%get_new_index(neighborhood))%ptr => subtypes%get_new_part()
!         end if
!         call subtypelist(subtypedict%get_index(neighborhood))%ptr%subset1%add(iatom)
!      end do
!      do i = 1, mnatypes%parts(h)%subset2%part_size
!         iatom = mnatypes%parts(h)%subset2%indices(i)
!         neighborhood = mnatypes%partition_map2%s(adjlists2(iatom)%n)
!         if (.not. subtypedict%has_index(neighborhood)) then
!            subtypelist(subtypedict%get_new_index(neighborhood))%ptr => subtypes%get_new_part()
!         end if
!         call subtypelist(subtypedict%get_index(neighborhood))%ptr%subset2%add(iatom)
!      end do
!      call subtypedict%reset()
!   end do
!
!end subroutine
!
!! Iteratively level up MNA types
!subroutine compute_mnatypes(mol1, mol2, mnatypes)
!   type(molecule_type), intent(in) :: mol1, mol2
!   type(partition_type), intent(out) :: mnatypes
!   ! Local variables
!   type(partition_type) :: subtypes
!   type(intlist_type), dimension(:), allocatable :: adjlists1, adjlists2
!
!   do i = 1, size(mol1%atoms)
!      adjlists1(i)%n = mol1%atoms(i)%adjlist
!   end do
!
!   do i = 1, size(mol2%atoms)
!      adjlists2(i)%n = mol2%atoms(i)%adjlist
!   end do
!
!   mnatypes = mol%eltypes
!
!   do
!      ! Compute next level MNA subtypes
!      call levelup_mnatypes(mol%atoms, mnatypes, subtypes)
!      ! Exit the loop if subtypes are unchanged
!      if (subtypes == mnatypes) exit
!      mnatypes = subtypes
!   end do
!
!   call mnatypes%print_parts()
!
!end subroutine

! Compute next level cross MNA types between mol1 and mol2
subroutine compute_crossmnatypes(adjlists1, adjlists2, nintype, intypes1, intypes2, &
      ntype, types1, types2)
   type(atomlist_type), dimension(:), intent(in) :: adjlists1, adjlists2
   integer, intent(in) :: nintype
   integer, dimension(:), intent(in) :: intypes1, intypes2
   integer, intent(out) :: ntype
   integer, dimension(:), intent(out) :: types1, types2
   ! Local variables
   integer :: h, i, j
   logical :: untyped(size(adjlists1))
   integer :: archeatom(size(adjlists1))

   ntype = 0
   untyped(:) = .true.

   do i = 1, size(adjlists1)
      if (untyped(i)) then
         ntype = ntype + 1
         types1(i) = ntype
         archeatom(ntype) = i
         do j = i + 1, size(adjlists1)
            if (untyped(j)) then
               if (intypes1(j) == intypes1(i)) then
                  if (same_adjacency(nintype, intypes1, adjlists1(i)%atomidcs, intypes1, adjlists1(j)%atomidcs)) then
                     types1(j) = ntype
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

   untyped(:) = .true.

   do h = 1, ntype
      do i = 1, size(adjlists1)
         if (untyped(i)) then
            if (intypes2(i) == intypes1(archeatom(h))) then
               if (same_adjacency(nintype, intypes1, adjlists1(archeatom(h))%atomidcs, intypes2, adjlists2(i)%atomidcs)) then
                  types2(i) = h
                  untyped(i) = .false.
               end if
            end if
         end if
      end do
   end do

   do i = 1, size(adjlists1)
      if (untyped(i)) then
         ntype = ntype + 1
         types2(i) = ntype
         do j = i + 1, size(adjlists1)
            if (untyped(j)) then
               if (intypes2(j) == intypes2(i)) then
                  if (same_adjacency(nintype, intypes2, adjlists2(i)%atomidcs, intypes2, adjlists2(j)%atomidcs)) then
                     types2(j) = ntype
                     untyped(j) = .false.
                  end if
               end if
            end if
         end do
      end if
   end do

end subroutine

! Test if adjacent atoms to a pair of atoms in mol1 and mol2 are the same
function same_adjacency(neltype, atomtype0, adjlist1, atomtype1, adjlist2) result(sameadj)
   integer, intent(in) :: neltype
   integer, dimension(:), intent(in) :: adjlist1, adjlist2
   integer, dimension(:) :: atomtype0, atomtype1
   ! Result variable
   logical :: sameadj
   ! Local variables
   integer :: i1, i2
   integer, dimension(neltype) :: n1, n2

   sameadj = .true.

   if (size(adjlist1) /= size(adjlist2)) then
      sameadj = .false.
      return
   end if

!  Find if adjacent atoms are the same

   n1(:) = 0
   n2(:) = 0

   do i1 = 1, size(adjlist1)
      n1(atomtype0(adjlist1(i1))) = n1(atomtype0(adjlist1(i1))) + 1
   end do

   do i2 = 1, size(adjlist2)
      n2(atomtype1(adjlist2(i2))) = n2(atomtype1(adjlist2(i2))) + 1
   end do

   if (any(n1 /= n2)) then
      sameadj = .false.
      return
   end if

end function

end module
