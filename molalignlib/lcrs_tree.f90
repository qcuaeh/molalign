module lcrs_tree
use iso_fortran_env, only: stdout => error_unit
implicit none
private

type, public :: item_node
   integer :: idx
   type(item_node), pointer :: next_item => null()
end type

type, public :: tree_node
   ! First item list
   type(item_node), pointer :: first_item1 => null()
   type(item_node), pointer :: last_item1 => null()
   integer :: item_count1 = 0

   ! Second item list
   type(item_node), pointer :: first_item2 => null()
   type(item_node), pointer :: last_item2 => null()
   integer :: item_count2 = 0

   ! Tree structure pointers
   type(tree_node), pointer :: first_child => null()
   type(tree_node), pointer :: last_child => null()
   type(tree_node), pointer :: next_sibling => null()
end type

type, public :: treenode_ptr
   type(tree_node), pointer :: ptr => null()
end type

type, public :: tree_type
   type(tree_node), pointer :: tree
   type(treenode_ptr), allocatable :: itemdir1(:)
   type(treenode_ptr), allocatable :: itemdir2(:)
end type

! Make types and procedures public
public :: new_tree
public :: new_child
public :: delete_tree
public :: print_tree
public :: print_items
public :: prune_node
public :: flatten_tree
public :: update_itemdir
public :: add_item1
public :: add_new_item1
public :: add_item2
public :: add_new_item2

contains

! First list procedures
subroutine add_new_item1(node, val)
   type(tree_node), intent(inout) :: node
   integer, intent(in) :: val
   type(item_node), pointer :: new_item

   allocate(new_item)
   new_item%idx = val
   new_item%next_item => null()

   if (associated(node%first_item1)) then
      node%last_item1%next_item => new_item
   else
      node%first_item1 => new_item
   end if

   node%last_item1 => new_item
   node%item_count1 = node%item_count1 + 1
end subroutine

subroutine add_item1(node, item)
   type(tree_node), intent(inout) :: node
   type(item_node), pointer, intent(in) :: item

   item%next_item => null()

   if (associated(node%first_item1)) then
      node%last_item1%next_item => item
   else
      node%first_item1 => item
   end if

   node%last_item1 => item
   node%item_count1 = node%item_count1 + 1
end subroutine

! Second list procedures
subroutine add_new_item2(node, val)
   type(tree_node), intent(inout) :: node
   integer, intent(in) :: val
   type(item_node), pointer :: new_item

   allocate(new_item)
   new_item%idx = val
   new_item%next_item => null()

   if (associated(node%first_item2)) then
      node%last_item2%next_item => new_item
   else
      node%first_item2 => new_item
   end if

   node%last_item2 => new_item
   node%item_count2 = node%item_count2 + 1
end subroutine

subroutine add_item2(node, item)
   type(tree_node), intent(inout) :: node
   type(item_node), pointer, intent(in) :: item

   item%next_item => null()

   if (associated(node%first_item2)) then
      node%last_item2%next_item => item
   else
      node%first_item2 => item
   end if

   node%last_item2 => item
   node%item_count2 = node%item_count2 + 1
end subroutine

function new_tree() result(root)
   type(tree_node), pointer :: root
   allocate(root)
   root%first_item1 => null()
   root%last_item1 => null()
   root%first_item2 => null()
   root%last_item2 => null()
   root%first_child => null()
   root%last_child => null()
   root%next_sibling => null()
   root%item_count1 = 0
   root%item_count2 = 0
end function

function new_child(parent) result(child)
   type(tree_node), target, intent(inout) :: parent
   type(tree_node), pointer :: child

   allocate(child)
   child%first_item1 => null()
   child%last_item1 => null()
   child%first_item2 => null()
   child%last_item2 => null()
   child%first_child => null()
   child%last_child => null()
   child%next_sibling => null()
   child%item_count1 = 0
   child%item_count2 = 0

   if (associated(parent%first_child)) then
      parent%last_child%next_sibling => child
   else
      parent%first_child => child
   end if

   parent%last_child => child
end function

recursive subroutine delete_tree(node)
   type(tree_node), pointer :: node
   type(item_node), pointer :: curr_item, next_item
   type(tree_node), pointer :: curr_child, next_child

   if (.not. associated(node)) return

   ! Delete first list items
   curr_item => node%first_item1
   do while (associated(curr_item))
      next_item => curr_item%next_item
      deallocate(curr_item)
      curr_item => next_item
   end do

   ! Delete second list items
   curr_item => node%first_item2
   do while (associated(curr_item))
      next_item => curr_item%next_item
      deallocate(curr_item)
      curr_item => next_item
   end do

   ! Delete children
   curr_child => node%first_child
   do while (associated(curr_child))
      next_child => curr_child%next_sibling
      call delete_tree(curr_child)
      curr_child => next_child
   end do

   deallocate(node)
   node => null()
end subroutine

subroutine prune_node(node)
   type(tree_node), pointer :: node
   type(tree_node), pointer :: curr_child, next_child

   if (.not. associated(node)) return

   curr_child => node%first_child
   do while (associated(curr_child))
      next_child => curr_child%next_sibling
      call delete_tree(curr_child)
      curr_child => next_child
   end do

   node%first_child => null()
   node%last_child => null()
end subroutine

subroutine print_items(node)
   type(tree_node), intent(in) :: node
   type(item_node), pointer :: item

   if (.not. associated(node%first_item1) .and. &
       .not. associated(node%first_item2)) then
      write(stdout, '(A)', advance='no') '()'
      return
   end if

   write(stdout, '(A)', advance='no') '('

   ! Print first list
   item => node%first_item1
   do while (associated(item))
      write(stdout, '(1X,I0)', advance='no') item%idx
      item => item%next_item
   end do

   write(stdout, '(1X,A)', advance='no') '|'

   ! Print second list
   item => node%first_item2
   do while (associated(item))
      write(stdout, '(1X,I0)', advance='no') item%idx
      item => item%next_item
   end do

   write(stdout, '(1X,A)', advance='no') ')'
end subroutine

subroutine print_tree(node)
   type(tree_node), intent(in) :: node
   write(stdout, *)
   call traverse_tree(node, 1)

contains
   recursive subroutine traverse_tree(node, indent)
      type(tree_node), intent(in) :: node
      integer, intent(in) :: indent
      type(tree_node), pointer :: child

      if (associated(node%first_child)) then
         child => node%first_child
         call print_items(node)
         write(stdout, '(A)', advance='no') '---'
         call traverse_tree(child, indent + 1)

         child => child%next_sibling
         do
            if (.not. associated(child)) return
            write(stdout, '(A,A)', advance='no') repeat('     ', indent)
            call traverse_tree(child, indent + 1)
            child => child%next_sibling
         end do
      end if

      ! Leaf node
      call print_items(node)
      write (stdout, *)
   end subroutine
end subroutine

subroutine flatten_tree(tree)
   type(tree_type), target, intent(inout) :: tree
   type(tree_node), pointer :: flat_root, curr_node, key_node
   type(item_node), pointer :: curr_item

   ! Create temporary root for flattened tree
   flat_root => new_tree()

   ! Collect leaves and sort them by total item count
   call collect_and_sort_leaves(tree%tree, flat_root)

   ! Update the original tree with the flattened and sorted structure
   tree%tree%first_child => flat_root%first_child
   tree%tree%last_child => flat_root%last_child
   flat_root%first_child => null()
   flat_root%last_child => null()

   ! Clean up temporary root
   deallocate(flat_root)

   ! Update itemdir to point to new structure
   call update_itemdir(tree%itemdir1, tree%itemdir2, tree%tree)

contains
   recursive subroutine collect_and_sort_leaves(node, flat_parent)
      type(tree_node), target, intent(in) :: node
      type(tree_node), target, intent(inout) :: flat_parent
      type(tree_node), pointer :: child, new_leaf
      type(tree_node), pointer :: curr, prev
      type(item_node), pointer :: src_item

      if (.not. associated(node%first_child)) then
         ! Leaf node - create new leaf
         new_leaf => new_tree()

         ! Copy items from first list in original order
         src_item => node%first_item1
         do while (associated(src_item))
            call add_new_item1(new_leaf, src_item%idx)
            src_item => src_item%next_item
         end do

         ! Copy items from second list in original order
         src_item => node%first_item2
         do while (associated(src_item))
            call add_new_item2(new_leaf, src_item%idx)
            src_item => src_item%next_item
         end do

         ! Insert into sorted position using insertion sort based on total items
         if (.not. associated(flat_parent%first_child)) then
            ! First child
            flat_parent%first_child => new_leaf
            flat_parent%last_child => new_leaf
         else
            ! Find insertion point
            prev => null()
            curr => flat_parent%first_child

            ! Traverse until we find the right position - sort only by first list count
            do while (associated(curr))
               if (new_leaf%item_count1 < curr%item_count1) exit
               prev => curr
               curr => curr%next_sibling
            end do

            ! Insert node
            if (.not. associated(prev)) then
               ! Insert at beginning
               new_leaf%next_sibling => flat_parent%first_child
               flat_parent%first_child => new_leaf
            else if (.not. associated(curr)) then
               ! Insert at end
               prev%next_sibling => new_leaf
               flat_parent%last_child => new_leaf
            else
               ! Insert in middle
               new_leaf%next_sibling => prev%next_sibling
               prev%next_sibling => new_leaf
            end if
         end if
      else
         ! Non-leaf node - process children
         child => node%first_child
         do while (associated(child))
            call collect_and_sort_leaves(child, flat_parent)
            child => child%next_sibling
         end do
      end if
   end subroutine
end subroutine

subroutine update_itemdir(itemdir1, itemdir2, root)
   type(treenode_ptr), intent(inout) :: itemdir1(:)
   type(treenode_ptr), intent(inout) :: itemdir2(:)
   type(tree_node), target, intent(in) :: root

   call traverse_leaves(root)

contains
   recursive subroutine traverse_leaves(node)
      type(tree_node), target, intent(in) :: node
      type(tree_node), pointer :: child
      type(item_node), pointer :: curr_item

      if (associated(node%first_child)) then
         child => node%first_child
         do
            call traverse_leaves(child)
            if (.not. associated(child%next_sibling)) return
            child => child%next_sibling
         end do
      end if

      ! Process first list items
      curr_item => node%first_item1
      do while (associated(curr_item))
         itemdir1(curr_item%idx)%ptr => node
         curr_item => curr_item%next_item
      end do

      ! Process second list items
      curr_item => node%first_item2
      do while (associated(curr_item))
         itemdir2(curr_item%idx)%ptr => node
         curr_item => curr_item%next_item
      end do
   end subroutine
end subroutine

end module
