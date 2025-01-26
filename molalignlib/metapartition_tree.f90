module metapartition_tree
use parameters
use partition
implicit none
private

type :: item_node
   integer :: value
   type(item_node), pointer :: next => null()
end type item_node

type, public :: tree_node
   type(item_node), pointer :: first_item => null()
   type(tree_node), pointer :: first_child => null()
   type(tree_node), pointer :: next_sibling => null()
   type(tree_node), pointer :: prev_node => null()
   type(partition_type) :: subtypes
contains
   procedure :: new_leaf
   procedure :: add_item
   procedure :: print_tree
end type tree_node

type, public :: nodepointer_type
   type(tree_node), pointer :: ptr => null()
end type

public :: new_tree
public :: delete_tree
public :: delete_branch

contains

function new_tree() result(root)
   type(tree_node), pointer :: root
   allocate(root)
end function new_tree

function new_leaf(this, subtypes) result(new_node)
   class(tree_node), target, intent(inout) :: this
   class(partition_type), target, intent(inout) :: subtypes
   type(tree_node), pointer :: new_node
   type(tree_node), pointer :: current
   
   allocate(new_node)
   new_node%first_item => null()
   new_node%first_child => null()
   new_node%next_sibling => null()
   new_node%prev_node => null()
   new_node%subtypes = subtypes
   
   if (.not. associated(this%first_child)) then
      this%first_child => new_node
      new_node%prev_node => this
   else
      current => this%first_child
      do while (associated(current%next_sibling))
         current => current%next_sibling
      end do
      current%next_sibling => new_node
      new_node%prev_node => current
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
      next_item => curr_item%next
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

subroutine add_item(this, val)
   class(tree_node), intent(inout) :: this
   integer, intent(in) :: val
   type(item_node), pointer :: new_item
   type(item_node), pointer :: current
   
   allocate(new_item)
   new_item%value = val
   new_item%next => null()
   
   if (.not. associated(this%first_item)) then
      this%first_item => new_item
   else
      current => this%first_item
      do while (associated(current%next))
         current => current%next
      end do
      current%next => new_item
   end if
end subroutine add_item

recursive subroutine print_tree(this, indent)
   class(tree_node), intent(in) :: this
   integer, intent(in), optional :: indent
   type(tree_node), pointer :: child
   type(item_node), pointer :: item
   integer :: i, current_indent
   character(len=4) :: indent_str
   
   current_indent = 0
   if (present(indent)) current_indent = indent
   
   indent_str = '    '
   do i = 1, current_indent
      write(stderr, '(A)', advance='no') indent_str
   end do
   
   write(stderr, '(A)', advance='no') '+-['
   if (associated(this%first_item)) then
      item => this%first_item
      do while (associated(item))
         write(stderr, '(I0)', advance='no') item%value
         item => item%next
         if (associated(item)) write(stderr, '(A)', advance='no') ','
      end do
   end if
   write(stderr, '(A)') ']'
   
   if (associated(this%first_child)) then
      child => this%first_child
      do while (associated(child))
         call child%print_tree(current_indent + 1)
         child => child%next_sibling
      end do
   end if
end subroutine print_tree

recursive subroutine delete_tree(node)
   type(tree_node), pointer :: node
   type(item_node), pointer :: curr_item, next_item
   type(tree_node), pointer :: curr_child, next_child
   
   if (.not. associated(node)) return
   
   curr_item => node%first_item
   do while (associated(curr_item))
      next_item => curr_item%next
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
