module partition
use stdio
use kinds
use part
implicit none

type, public :: partition_type
   integer :: partition_size
   integer :: max_partition_size
   integer :: max_part_size
   integer, pointer :: largest_part_size
   integer, pointer :: partition_map(:)
   type(part_type), pointer :: parts(:) => null()
contains
   procedure :: init => partition_init
   procedure :: reset => partition_reset
   procedure :: new_part => partition_new_part
   procedure :: add_part => partition_add_part
   procedure :: print_parts => partition_print_parts
end type

interface assignment (=)
   module procedure partition_assignment
end interface

interface operator (==)
   module procedure partition_equality
end interface

contains

subroutine partition_assignment(left, right)
   class(partition_type), intent(inout) :: left
   type(partition_type), intent(in) :: right
   ! Local variables
   integer :: h

   call left%init(right%max_partition_size, right%max_part_size)

   do h = 1, right%partition_size
      call left%add_part(right%parts(h))
   end do

end subroutine

subroutine partition_init(self, max_partition_size, max_part_size)
   class(partition_type), intent(inout) :: self
   integer, intent(in) :: max_partition_size, max_part_size

   if (associated(self%parts)) then
      error stop 'Memory leak'
   end if

   allocate (self%parts(max_partition_size))
   allocate (self%partition_map(max_partition_size))
   allocate (self%largest_part_size)

   self%partition_size = 0
   self%largest_part_size = 0
   self%max_partition_size = max_partition_size
   self%max_part_size = max_part_size

end subroutine

subroutine partition_reset(self)
   class(partition_type), intent(inout) :: self
   integer :: h

   do h = 1, self%partition_size
      deallocate (self%parts(h)%indices)
   end do

   deallocate (self%parts)
   deallocate (self%partition_map)
   deallocate (self%largest_part_size)

end subroutine

function partition_new_part(self) result(part)
   class(partition_type), intent(inout) :: self
   ! Result variable
   type(part_type), pointer :: part

   ! Increase partition size
   self%partition_size = self%partition_size + 1

   ! Initialize new part
   call self%parts(self%partition_size)%init( &
      self%partition_size, & ! index
      self%max_part_size, &
      self%largest_part_size, &
      self%partition_map)

   ! Return pointer to new part
   part => self%parts(self%partition_size)

end function

subroutine partition_add_part(self, part)
   class(partition_type), intent(inout) :: self
   class(part_type), intent(in) :: part
   integer :: i

   ! Increase partition size
   self%partition_size = self%partition_size + 1

   ! Initialize new part
   call self%parts(self%partition_size)%init( &
      self%partition_size, & ! index
      self%max_part_size, &
      self%largest_part_size, &
      self%partition_map)

   ! Update new part size
   self%parts(self%partition_size)%part_size = part%part_size

   ! Update new part indices
   self%parts(self%partition_size)%indices(:part%part_size) = part%indices(:part%part_size)

   ! Add indices to part map
   do i = 1, part%part_size
      self%partition_map(part%indices(i)) = self%partition_size
   end do

   ! Update max part size
   if (part%part_size > self%largest_part_size) then
      self%largest_part_size = part%part_size
   end if

end subroutine

subroutine partition_print_parts(self)
   class(partition_type), intent(in) :: self
   integer :: h
   character(:), allocatable :: fmtstr

   write (stderr, *)
   do h = 1, self%partition_size
      fmtstr = "('{'" // repeat(',1x,i3', self%parts(h)%part_size) // ",' }')"
      write (stderr, fmtstr) self%parts(h)%indices(:self%parts(h)%part_size)
   end do

end subroutine

function partition_equality(left, right) result(equality)
   type(partition_type), intent(in) :: left, right
   ! Result variable
   logical :: equality
   ! Local variables
   integer :: h, i

   if (left%partition_size /= right%partition_size) then
      equality = .false.
      return
   end if

   do h = 1, left%partition_size

      if (left%parts(h)%part_size /= right%parts(h)%part_size) then
         equality = .false.
         return
      end if

      do i = 1, left%parts(h)%part_size
         if (left%parts(h)%indices(i) /= right%parts(h)%indices(i)) then
            equality = .false.
            return
         end if
      end do

   end do

   equality = .true.

end function

end module
