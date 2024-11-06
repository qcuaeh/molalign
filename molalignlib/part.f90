module part
use stdio
use kinds
implicit none

type, public :: part_type
   integer :: index
   integer :: part_size
   integer, allocatable :: indices(:)
   integer, pointer :: largest_part_size
   integer, pointer :: partition_map(:)
contains
   procedure :: add => part_add
   procedure :: init => part_init
end type

type, public :: partpointer_type
   type(part_type), pointer :: ptr => null()
end type

contains

subroutine part_init(self, part_index, max_part_size, largest_part_size, partition_map)
   class(part_type), intent(inout) :: self
   integer, intent(in) :: part_index, max_part_size
   integer, pointer, intent(in) :: largest_part_size
   integer, pointer, intent(in) :: partition_map(:)

   ! Initialize new part size
   self%part_size = 0

   ! Set part index to partition current index
   self%index = part_index

   ! Allocate new part indices
   allocate (self%indices(max_part_size))

   ! Point largest part size pointer to partition largest part size
   self%largest_part_size => largest_part_size

   ! Point map pointer to partition map
   self%partition_map => partition_map

end subroutine

subroutine part_add(self, element)
   class(part_type), intent(inout) :: self
   integer, intent(in) :: element

   ! Increase part size
   self%part_size = self%part_size + 1

   ! Add element to part
   self%indices(self%part_size) = element

   ! Add index to part map
   self%partition_map(element) = self%index

   ! Update max part size
   if (self%part_size > self%largest_part_size) then
      self%largest_part_size = self%part_size
   end if

end subroutine

end module
