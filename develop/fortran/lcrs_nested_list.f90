module lcrs_nested_list
   implicit none
   private

   ! Item node: part of a linked list holding values at each tree node
   type :: item_node
      integer :: value
      type(item_node), pointer :: next => null()
   end type item_node

   ! Tree node: LCRS tree with item lists at each node
   type, public :: tree_node
      type(item_node), pointer :: first_item => null()     ! Head of item list
      type(tree_node), pointer :: first_child => null()    ! Leftmost child
      type(tree_node), pointer :: next_sibling => null()   ! Next node at same level
      integer :: index = 0                                 ! Node's position among siblings
   contains
      procedure :: new_branch  ! Add child or sibling
      procedure :: add_item    ! Add value to item list
      procedure :: print => print_node
   end type tree_node

   public :: new_branch
   public :: add_item
   public :: print_node
   public :: delete_node

contains
   ! Creates new node: first creation becomes child (index 1), 
   ! subsequent ones become siblings (incrementing indices)
   function new_branch(this) result(new_node)
      class(tree_node), intent(inout) :: this
      type(tree_node), pointer :: new_node
      type(tree_node), pointer :: current
      
      allocate(new_node)
      new_node%first_item => null()
      new_node%first_child => null()
      new_node%next_sibling => null()
      
      if (.not. associated(this%first_child)) then
         this%first_child => new_node        ! First becomes child
         new_node%index = 1                  ! First child gets index 1
      else
         current => this%first_child
         do while (associated(current%next_sibling))  ! Find last sibling
            current => current%next_sibling
         end do
         new_node%index = current%index + 1  ! Increment index for new sibling
         current%next_sibling => new_node    ! Add as last sibling
      end if
   end function new_branch

   ! Adds value to end of node's item list
   subroutine add_item(this, val)
      class(tree_node), intent(inout) :: this
      integer, intent(in) :: val
      type(item_node), pointer :: new_item
      type(item_node), pointer :: current
      
      allocate(new_item)
      new_item%value = val
      new_item%next => null()
      
      if (.not. associated(this%first_item)) then
         this%first_item => new_item         ! First item in list
      else
         current => this%first_item
         do while (associated(current%next))  ! Find last item
            current => current%next
         end do
         current%next => new_item            ! Add to end
      end if
   end subroutine add_item

   ! Prints items lists with index trace and braces for specified level
   ! Format: [index_trace]: {item1,item2,...}
   ! level = -1 prints all levels (default)
   recursive subroutine print_node(this, level, index_trace, current_level)
      class(tree_node), intent(in) :: this
      integer, intent(in) :: level         ! Level to print
      character(len=*), intent(in), optional :: index_trace
      integer, intent(in), optional :: current_level
      character(len=100) :: current_trace
      type(tree_node), pointer :: child
      type(item_node), pointer :: item
      logical :: need_comma, has_items
      integer :: curr_level
      
      ! Set current level
      curr_level = 0
      if (present(current_level)) curr_level = current_level
      
      ! Build current index trace
      if (present(index_trace)) then
         write(current_trace, '(A,".",I0)') trim(index_trace), this%index
      else
         write(current_trace, '(I0)') this%index
      end if

      ! Print items if at target level
      if (level == curr_level) then
         write(*, '(A,": ")', advance='no') trim(current_trace)
         write(*, '(A)', advance='no') '{'
         has_items = associated(this%first_item)
         if (has_items) then
            need_comma = .false.
            item => this%first_item
            do while (associated(item))
               if (need_comma) write(*, '(A)', advance='no') ','
               write(*, '(I0)', advance='no') item%value
               need_comma = .true.
               item => item%next
            end do
         end if
         write(*, '(A)') '}'
      end if
      
      ! Recursively print children
      if (associated(this%first_child)) then
         child => this%first_child
         do while (associated(child))
            call child%print(level, trim(current_trace), curr_level + 1)
            child => child%next_sibling
         end do
      end if
   end subroutine print_node

   ! Recursively deletes node and all its subtrees
   recursive subroutine delete_node(node)
      type(tree_node), pointer :: node       ! Will be deallocated and nullified
      type(item_node), pointer :: curr_item, next_item
      type(tree_node), pointer :: curr_child, next_child
      
      if (.not. associated(node)) return
      
      ! Delete item list
      curr_item => node%first_item
      do while (associated(curr_item))
         next_item => curr_item%next         ! Save next before dealloc
         deallocate(curr_item)
         curr_item => next_item
      end do
      
      ! Delete child subtrees
      curr_child => node%first_child
      do while (associated(curr_child))
         next_child => curr_child%next_sibling  ! Save next before recursive delete
         call delete_node(curr_child)           ! Recursively delete subtree
         curr_child => next_child
      end do
      
      deallocate(node)                       ! Delete node itself
      node => null()
   end subroutine delete_node

end module lcrs_nested_list

program test_lcrs
   use lcrs_nested_list
   implicit none
   
   type(tree_node), pointer :: root, node1, node2
   integer :: i
   
   ! Create a single tree with test cases sharing branches where possible
   allocate(root)
   root%index = 0  ! Root node at level 0
   
   ! Level 0: Root list
   call root%add_item(1)
   call root%add_item(2)
   call root%add_item(3)
   call root%add_item(4)
   call root%add_item(5)
   
   ! Level 1: Two sublists
   node1 => root%new_branch()     ! First branch
   call node1%add_item(1)
   call node1%add_item(2)
   node2 => node1%new_branch()    ! Child of first branch with {1,2}
   call node2%add_item(1)
   call node2%add_item(2)
   
   node1 => root%new_branch()     ! Second branch
   call node1%add_item(3)
   call node1%add_item(4)
   call node1%add_item(5)
   node2 => node1%new_branch()    ! First child of second branch with {3,4}
   call node2%add_item(3)
   call node2%add_item(4)
   
   node2 => node1%new_branch()    ! Second child of second branch with {5}
   call node2%add_item(5)
   
   ! Print each level
   do i = 0, 2
      write(*, '(A,I0)') "Level ", i
      call root%print(i)
      write(*, *)
   end do
   
   call delete_node(root)
   
end program test_lcrs
