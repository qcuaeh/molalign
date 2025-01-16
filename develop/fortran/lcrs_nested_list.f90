module lcrs_nested_list
   implicit none
   private

   ! ID manager type to handle unique IDs for each tree
   type :: id_manager
      integer :: next_id = 0    ! Starts at 0 so first ID will be 1
   contains
      procedure :: get_next_id
   end type id_manager

   ! Item node: part of a linked list holding values at each tree node
   type :: item_node
      integer :: value
      type(item_node), pointer :: next => null()
   end type item_node

   ! Tree node: LCRS tree with item lists at each node
   type, public :: tree_node
      private
      type(item_node), pointer :: first_item => null()     ! Head of item list
      type(tree_node), pointer :: first_child => null()    ! Leftmost child
      type(tree_node), pointer :: next_sibling => null()   ! Next node at same level
      type(id_manager), pointer :: id_mgr => null()        ! Pointer to tree's ID manager
      integer :: node_id = 0                              ! Node identifier
   contains
      procedure :: new_branch      ! Add child or sibling
      procedure :: add_item        ! Add value to item list
      procedure :: print => print_node
   end type tree_node

   public :: new_tree         ! Constructor for new trees
   public :: delete_node

contains
   ! ID manager methods
   function get_next_id(this) result(id)
      class(id_manager), intent(inout) :: this
      integer :: id
      this%next_id = this%next_id + 1
      id = this%next_id
   end function get_next_id

   ! Constructor for new trees
   function new_tree() result(root)
      type(tree_node), pointer :: root
      
      allocate(root)
      allocate(root%id_mgr)
      root%id_mgr%next_id = 0     ! Start at 0 so first ID will be 1
      root%node_id = root%id_mgr%get_next_id()  ! Root gets ID 1
   end function new_tree

   ! Creates new node as child or sibling
   function new_branch(this) result(new_node)
      class(tree_node), intent(inout) :: this
      type(tree_node), pointer :: new_node
      type(tree_node), pointer :: current
      
      allocate(new_node)
      new_node%first_item => null()
      new_node%first_child => null()
      new_node%next_sibling => null()
      new_node%id_mgr => this%id_mgr           ! Share ID manager with parent
      new_node%node_id = this%id_mgr%get_next_id()
      
      if (.not. associated(this%first_child)) then
         this%first_child => new_node        ! First becomes child
      else
         current => this%first_child
         do while (associated(current%next_sibling))  ! Find last sibling
            current => current%next_sibling
         end do
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

   ! Print node with just ID
   recursive subroutine print_node(this, level, current_level)
      class(tree_node), intent(in) :: this
      integer, intent(in) :: level
      integer, intent(in), optional :: current_level
      type(tree_node), pointer :: child
      type(item_node), pointer :: item
      logical :: need_comma
      integer :: curr_level
      
      curr_level = 0
      if (present(current_level)) curr_level = current_level

      if (level == curr_level) then
         write(*, '("Node ",I0,": {")', advance='no') this%node_id
         if (associated(this%first_item)) then
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
      
      if (associated(this%first_child)) then
         child => this%first_child
         do while (associated(child))
            call child%print(level, curr_level + 1)
            child => child%next_sibling
         end do
      end if
   end subroutine print_node

   ! Recursively deletes node and all its subtrees
   recursive subroutine delete_node(node)
      type(tree_node), pointer :: node
      type(item_node), pointer :: curr_item, next_item
      type(tree_node), pointer :: curr_child, next_child
      
      if (.not. associated(node)) return
      
      ! Delete item list
      curr_item => node%first_item
      do while (associated(curr_item))
         next_item => curr_item%next
         deallocate(curr_item)
         curr_item => next_item
      end do
      
      ! Delete child subtrees
      curr_child => node%first_child
      do while (associated(curr_child))
         next_child => curr_child%next_sibling
         call delete_node(curr_child)
         curr_child => next_child
      end do
      
      ! Delete ID manager if this is a root node (has id 1)
      if (node%node_id == 1) then
         deallocate(node%id_mgr)
         node%id_mgr => null()
      end if
      
      deallocate(node)
      node => null()
   end subroutine delete_node

end module lcrs_nested_list

! Test program demonstrating multiple trees
program test_lcrs
   use lcrs_nested_list
   implicit none
   
   type(tree_node), pointer :: root1, root2, node1, node2, node3
   integer :: i
   
   ! Create first tree (original example)
   root1 => new_tree()
   
   ! Level 0: Root list
   call root1%add_item(1)
   call root1%add_item(2)
   call root1%add_item(3)
   call root1%add_item(4)
   call root1%add_item(5)
   
   ! Level 1: Two sublists
   node1 => root1%new_branch()     ! First branch
   call node1%add_item(1)
   call node1%add_item(2)
   node2 => node1%new_branch()     ! Child of first branch with {1,2}
   call node2%add_item(1)
   call node2%add_item(2)
   
   node1 => root1%new_branch()     ! Second branch
   call node1%add_item(3)
   call node1%add_item(4)
   call node1%add_item(5)
   node2 => node1%new_branch()     ! First child of second branch with {3,4}
   call node2%add_item(3)
   call node2%add_item(4)
   
   node2 => node1%new_branch()     ! Second child of second branch with {5}
   call node2%add_item(5)
   
   ! Create second tree (three levels deep)
   root2 => new_tree()
   
   ! Level 0: Root
   call root2%add_item(10)
   call root2%add_item(20)
   
   ! Level 1: First child branch
   node1 => root2%new_branch()
   call node1%add_item(100)
   call node1%add_item(200)
   
   ! Level 2: Children of first branch
   node2 => node1%new_branch()
   call node2%add_item(1000)
   
   node2 => node1%new_branch()
   call node2%add_item(2000)
   
   ! Level 1: Second child branch
   node1 => root2%new_branch()
   call node1%add_item(300)
   call node1%add_item(400)
   
   ! Level 2: Child of second branch
   node2 => node1%new_branch()
   call node2%add_item(3000)
   
   ! Level 3: Child of previous node
   node3 => node2%new_branch()
   call node3%add_item(30000)
   
   ! Print first tree
   write(*, '(A)') "First Tree (Original Example):"
   do i = 0, 2
      write(*, '(A,I0)') "Level ", i
      call root1%print(i)
      write(*, *)
   end do
   
   ! Print second tree
   write(*, '(A)') "Second Tree (Three Levels):"
   do i = 0, 3
      write(*, '(A,I0)') "Level ", i
      call root2%print(i)
      write(*, *)
   end do
   
   ! Clean up
   call delete_node(root1)
   call delete_node(root2)
   
end program test_lcrs
