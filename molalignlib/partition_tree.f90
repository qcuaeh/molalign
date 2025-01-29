module partition_tree
use iso_fortran_env, only: stdout => error_unit
implicit none
private

type, public :: item_node
   integer :: idx
   type(item_node), pointer :: next_item => null()
end type item_node

type, public :: tree_node
   type(item_node), pointer :: first_item => null()
   type(tree_node), pointer :: first_child => null()
   type(tree_node), pointer :: next_sibling => null()
   type(tree_node), pointer :: prev_node => null()
end type tree_node

public :: new_tree
public :: delete_tree
public :: delete_branch
public :: new_leaf
public :: new_child
public :: new_sibling
public :: add_item
public :: print_tree

contains

function new_tree() result(root)
   type(tree_node), pointer :: root
   allocate(root)
   root%first_item => null()
   root%first_child => null()
   root%next_sibling => null()
   root%prev_node => null()
end function new_tree

function new_child(parent) result(child)
   type(tree_node), target, intent(inout) :: parent
   type(tree_node), pointer :: child
   
   allocate(child)
   child%first_item => null()
   child%first_child => null()
   child%next_sibling => null()
   child%prev_node => parent
   
   if (.not. associated(parent%first_child)) then
      parent%first_child => child
   else
      child%next_sibling => parent%first_child
      parent%first_child%prev_node => child
      parent%first_child => child
   end if
end function new_child

function new_sibling(node) result(sibling)
   type(tree_node), target, intent(inout) :: node
   type(tree_node), pointer :: sibling
   type(tree_node), pointer :: current
   
   allocate(sibling)
   sibling%first_item => null()
   sibling%first_child => null()
   sibling%next_sibling => null()
   sibling%prev_node => node
   
   current => node
   do while (associated(current%next_sibling))
      current => current%next_sibling
   end do
   current%next_sibling => sibling
end function new_sibling

function new_leaf(node) result(new_node)
   type(tree_node), target, intent(inout) :: node
   type(tree_node), pointer :: new_node
   
   if (.not. associated(node%first_child)) then
      new_node => new_child(node)
   else
      new_node => new_sibling(node%first_child)
   end if
end function new_leaf

subroutine delete_branch(node)
   type(tree_node), pointer :: node
   type(item_node), pointer :: curr_item, next_item
   type(tree_node), pointer :: curr_child, next_child
   type(tree_node), pointer :: next_sibling
   
   if (.not. associated(node)) return
   
   next_sibling => node%next_sibling
   
   curr_item => node%first_item
   do while (associated(curr_item))
      next_item => curr_item%next_item
      deallocate(curr_item)
      curr_item => next_item
   end do
   
   curr_child => node%first_child
   do while (associated(curr_child))
      next_child => curr_child%next_sibling
      call delete_tree(curr_child)
      curr_child => next_child
   end do
   
   if (associated(node%prev_node)) then
      if (associated(node%prev_node%first_child, node)) then
         node%prev_node%first_child => next_sibling
      else
         node%prev_node%next_sibling => next_sibling
      end if
   end if
   
   if (associated(next_sibling)) then
      next_sibling%prev_node => node%prev_node
   end if
   
   deallocate(node)
   node => next_sibling
end subroutine delete_branch

subroutine add_item(node, val)
   type(tree_node), intent(inout) :: node
   integer, intent(in) :: val
   type(item_node), pointer :: new_item
   type(item_node), pointer :: current
   
   allocate(new_item)
   new_item%idx = val
   new_item%next_item => null()
   
   if (.not. associated(node%first_item)) then
      node%first_item => new_item
   else
      current => node%first_item
      do while (associated(current%next_item))
         current => current%next_item
      end do
      current%next_item => new_item
   end if
end subroutine add_item

subroutine print_items(node)
   type(tree_node), intent(in) :: node
   type(item_node), pointer :: item

   write(stdout, '(A)', advance='no') '['
   item => node%first_item
   do while (associated(item))
      write(stdout, '(1X,I0)', advance='no') item%idx
      item => item%next_item
   end do
   write(stdout, '(1X,A)') ']'
end subroutine print_items

subroutine print_tree(node)
   type(tree_node), intent(in) :: node

   write (stdout, *)
   write(stdout, '(A,A)', advance='no') repeat('    ', 0), '*---'
   call traverse_tree(node, 1)

   contains
   recursive subroutine traverse_tree(node, indent)
      type(tree_node), intent(in) :: node
      integer, intent(in) :: indent
      type(tree_node), pointer :: child

      if (.not. associated(node%first_child)) then
         ! Leaf node - print items
         call print_items(node)
      else if (.not. associated(node%first_child%next_sibling)) then
         ! Single child - continue chain
         write(stdout, '(A)', advance='no') '*---'
         call traverse_tree(node%first_child, indent)
      else
         ! Multiple children - print items and process children
         call print_items(node)
         child => node%first_child
         do while (associated(child))
            write(stdout, '(A,A)', advance='no') repeat('    ', indent), '*---'
            call traverse_tree(child, indent + 1)
            child => child%next_sibling
         end do
      end if
   end subroutine traverse_tree
end subroutine print_tree

recursive subroutine delete_tree(node)
   type(tree_node), pointer :: node
   type(item_node), pointer :: curr_item, next_item
   type(tree_node), pointer :: curr_child, next_child
   
   if (.not. associated(node)) return
   
   curr_item => node%first_item
   do while (associated(curr_item))
      next_item => curr_item%next_item
      deallocate(curr_item)
      curr_item => next_item
   end do
   
   curr_child => node%first_child
   do while (associated(curr_child))
      next_child => curr_child%next_sibling
      call delete_tree(curr_child)
      curr_child => next_child
   end do
   
   deallocate(node)
   node => null()
end subroutine delete_tree

end module
