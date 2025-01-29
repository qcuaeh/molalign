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
use partition_tree
use metapartition
use partitiondict
use permutation

implicit none

type :: tree_node_pointer
   type(tree_node), pointer :: ptr => null()
end type

type :: neighborhood_item
   type(tree_node_pointer), allocatable :: neighborhood(:)
   type(tree_node), pointer :: node
end type

type :: neighborhood_table
   integer :: num_items
   type(neighborhood_item), allocatable :: items(:)
end type

type :: atomtype_item
   integer :: elnum
   integer :: label
   type(tree_node), pointer :: node
end type

type :: atomtype_table
   integer :: num_items
   type(atomtype_item), allocatable :: items(:)
end type

interface find_item
   module procedure atomtypetable_find_item
   module procedure neighborhoodtable_find_item
end interface

interface append_item
   module procedure atomtypetable_append_item
   module procedure neighborhoodtable_append_item
end interface

interface operator (==)
   module procedure neighborhood_equality
end interface

contains

function neighborhood_equality( array1, array2) result(equality)
   type(tree_node_pointer), dimension(:), intent(in) :: array1, array2
   logical :: equality
   integer :: i, j, matches
   
   ! Check if sizes are equal
   if (size(array1) /= size(array2)) then
      equality = .false.
      return
   end if

   ! Early exit: Check if any value appears more times in one array
   ! than it does in the other
   do i = 1, size(array1)
       matches = 0
       do j = 1, size(array1)
           if (associated(array1(i)%ptr, array2(j)%ptr)) matches = matches + 1
           if (associated(array1(i)%ptr, array1(j)%ptr)) matches = matches - 1
       end do
       if (matches /= 0) then
           equality = .false.
           return
       end if
   end do

   equality = .true.

end function

subroutine atomtypetable_append_item(atomtypetable, elnum, label, node)
   type(atomtype_table), intent(inout) :: atomtypetable
   integer, intent(in) :: elnum
   integer, intent(in) :: label
   type(tree_node), pointer, intent(in) :: node

   atomtypetable%num_items = atomtypetable%num_items + 1
   atomtypetable%items(atomtypetable%num_items)%elnum = elnum
   atomtypetable%items(atomtypetable%num_items)%label = label
   atomtypetable%items(atomtypetable%num_items)%node => node

end subroutine

function atomtypetable_find_item(atomtypetable, elnum, label) result(node)
   type(atomtype_table), intent(in) :: atomtypetable
   integer, intent(in) :: elnum
   integer, intent(in) :: label
   type(tree_node), pointer :: node
   integer :: i

   do i = 1, atomtypetable%num_items
      if (atomtypetable%items(i)%elnum == elnum .and. &
          atomtypetable%items(i)%label == label) then
         node => atomtypetable%items(i)%node
         return
      end if
   end do

   node => null()

end function

subroutine neighborhoodtable_append_item(neighborhoodtable, neighborhood, node)
   type(neighborhood_table), intent(inout) :: neighborhoodtable
   type(tree_node_pointer), intent(in) :: neighborhood(:)
   type(tree_node), pointer, intent(in) :: node

   neighborhoodtable%num_items = neighborhoodtable%num_items + 1
   neighborhoodtable%items(neighborhoodtable%num_items)%neighborhood = neighborhood
   neighborhoodtable%items(neighborhoodtable%num_items)%node => node

end subroutine

function neighborhoodtable_find_item(neighborhoodtable, neighborhood) result(node)
   type(neighborhood_table), intent(in) :: neighborhoodtable
   type(tree_node_pointer), intent(in) :: neighborhood(:)
   type(tree_node), pointer :: node
   integer :: i

   do i = 1, neighborhoodtable%num_items
      if (neighborhoodtable%items(i)%neighborhood == neighborhood) then
         node => neighborhoodtable%items(i)%node
         return
      end if
   end do

   node => null()

end function

! Partition atoms by atomic number and label
subroutine compute_eltypes(mol, eltypes, typemap, num_leaves)
   type(mol_type), intent(in) :: mol
   type(tree_node), pointer, intent(out) :: eltypes
   type(tree_node_pointer), intent(out) :: typemap(:)
   integer, intent(out) :: num_leaves
   ! Local variables
   type(tree_node), pointer :: inode
   type(atomtype_table) :: atomtypetable
   integer :: i, elnum, label

   eltypes => new_tree()
   allocate (atomtypetable%items(size(mol%atoms)))
   atomtypetable%num_items = 0
   num_leaves = 0

   do i = 1, size(mol%atoms)
      call add_item(eltypes, i)
      elnum = mol%atoms(i)%elnum
      label = mol%atoms(i)%label
      inode => find_item(atomtypetable, elnum, label)
      if (.not. associated(inode)) then
         inode => new_leaf(eltypes)
         call append_item(atomtypetable, elnum, label, inode)
         num_leaves = num_leaves + 1
      end if
      call add_item(inode, i)
      typemap(i)%ptr => inode
   end do

end subroutine

! Compute Next Level MNA Types
subroutine nextlevel_mnatypes(mol, mnatypes, typemap, num_leaves)
   type(mol_type), intent(in) :: mol
   type(tree_node), pointer, intent(inout) :: mnatypes
   type(tree_node_pointer), intent(inout) :: typemap(:)
   integer, intent(out) :: num_leaves
   ! Local variables
   type(tree_node_pointer), allocatable :: curr_typemap(:)
   type(tree_node_pointer), allocatable :: neighborhood(:)
   type(neighborhood_table) :: neighborhoodtable

   allocate (curr_typemap(size(mol%atoms)))
   allocate (neighborhoodtable%items(size(mol%atoms)))

   num_leaves = 0
   curr_typemap = typemap
   call traverse_tree(mnatypes)
   
   contains

   recursive subroutine traverse_tree(inode)
      type(tree_node), intent(inout) :: inode
      type(tree_node), pointer :: jnode
      type(item_node), pointer :: item

      ! Check if node is a branch or a leaf
      if (associated(inode%first_child)) then

         ! Not a leaf, process its children
         jnode => inode%first_child
         call traverse_tree(jnode)
         do while (associated(jnode%next_sibling))
            jnode => jnode%next_sibling
            call traverse_tree(jnode)
         end do

      else

         ! A leaf, process its items
         item => inode%first_item
         neighborhoodtable%num_items = 0
         do while (associated(item))
            ! Process each item
            neighborhood = curr_typemap(mol%atoms(item%idx)%adjlist)
            jnode => find_item(neighborhoodtable, neighborhood)
            if (.not. associated(jnode)) then
               jnode => new_leaf(inode)
               call append_item(neighborhoodtable, neighborhood, jnode)
               num_leaves = num_leaves + 1
            end if
            call add_item(jnode, item%idx)
            typemap(item%idx)%ptr => jnode
            item => item%next_item
         end do

      end if
      
   end subroutine

end subroutine

! Iteratively compute MNA types
subroutine compute_mnatypes_tree(mol, mnatypes)
   type(mol_type), intent(in) :: mol
   type(tree_node), pointer, intent(inout) :: mnatypes
   ! Local variables
   integer :: num_leaves, prev_num_leaves
   type(tree_node_pointer), allocatable :: typemap(:)

   allocate (typemap(size(mol%atoms)))
   call compute_eltypes(mol, mnatypes, typemap, num_leaves)

   do

!      call mnatypes%print_tree()
      prev_num_leaves = num_leaves

      ! Compute MNA upper level types
      call nextlevel_mnatypes(mol, mnatypes, typemap, num_leaves)

      ! Exit loop if types did not change
      if (num_leaves == prev_num_leaves) exit

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
