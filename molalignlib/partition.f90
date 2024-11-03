module partition
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
end type

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

type, public :: bipart_type
   integer :: index
   type(part_type) :: subset1
   type(part_type) :: subset2
end type

type, public :: bipartition_type
   integer :: partition_size
   integer :: max_partition_size
   integer :: max_part_size
   integer, pointer :: largest_part_size
   integer, pointer :: partition_map1(:)
   integer, pointer :: partition_map2(:)
   type(bipart_type), pointer :: parts(:)
contains
   procedure :: init => bipartition_init
   procedure :: new_part => bipartition_new_part
   procedure :: print_parts => bipartition_print_parts
end type

type, public :: partpointer_type
   type(part_type), pointer :: ptr => null()
end type

type, public :: bipartpointer_type
   type(bipart_type), pointer :: ptr => null()
end type

interface assignment (=)
   module procedure partition_assignment
end interface

interface operator (==)
   module procedure partition_equality
   module procedure bipartition_equality
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

   ! Initialize new part size
   self%parts(self%partition_size)%part_size = 0

   ! Allocate new part indices
   allocate (self%parts(self%partition_size)%indices(self%max_part_size))

   ! Set new part index to partition current num part
   self%parts(self%partition_size)%index = self%partition_size

   ! Point new part max part size pointer to partition max part size
   self%parts(self%partition_size)%largest_part_size => self%largest_part_size

   ! Point new part map pointer to partition maps
   self%parts(self%partition_size)%partition_map => self%partition_map

   ! Return pointer to new part
   part => self%parts(self%partition_size)

end function

subroutine partition_add_part(self, part)
   class(partition_type), intent(inout) :: self
   class(part_type), intent(in) :: part
   integer :: i

   ! Increase partition size
   self%partition_size = self%partition_size + 1

   ! Initialize new part size
   self%parts(self%partition_size)%part_size = part%part_size

   ! Allocate new part indices
   allocate (self%parts(self%partition_size)%indices(self%max_part_size))

   ! Populate new part indices
   self%parts(self%partition_size)%indices(:part%part_size) = part%indices(:part%part_size)

   ! Set new part index to partition current num part
   self%parts(self%partition_size)%index = self%partition_size

   ! Point new part max part size pointer to partition max part size
   self%parts(self%partition_size)%largest_part_size => self%largest_part_size

   ! Point new part map pointer to partition maps
   self%parts(self%partition_size)%partition_map => self%partition_map

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

   do h = 1, self%partition_size
      fmtstr = "('{'" // repeat(',1x,i3', self%parts(h)%part_size) // ",' }')"
      write (stderr, fmtstr) self%parts(h)%indices(:self%parts(h)%part_size)
   end do

end subroutine

function partition_equality(self, other) result(equality)
   type(partition_type), intent(in) :: self, other
   ! Result variable
   logical :: equality
   ! Local variables
   integer :: h, i

   if (self%partition_size /= other%partition_size) then
      equality = .false.
      return
   end if

   do h = 1, self%partition_size

      if (self%parts(h)%part_size /= other%parts(h)%part_size) then
         equality = .false.
         return
      end if

      do i = 1, self%parts(h)%part_size
         if (self%parts(h)%indices(i) /= other%parts(h)%indices(i)) then
            equality = .false.
            return
         end if
      end do

   end do

   equality = .true.

end function

subroutine bipartition_init(self, max_partition_size, max_part_size)
   class(bipartition_type), intent(inout) :: self
   integer, intent(in) :: max_partition_size, max_part_size

   allocate (self%parts(max_partition_size))
   allocate (self%partition_map1(max_part_size))
   allocate (self%partition_map2(max_part_size))
   allocate (self%largest_part_size)

   self%partition_size = 0
   self%largest_part_size = 0
   self%max_partition_size = max_partition_size
   self%max_part_size = max_part_size

end subroutine

function bipartition_new_part(self) result(part)
   class(bipartition_type), intent(inout) :: self
   ! Result variable
   type(bipart_type), pointer :: part

   ! Increase partition size
   self%partition_size = self%partition_size + 1

   ! Initialize new part size
   self%parts(self%partition_size)%subset1%part_size = 0
   self%parts(self%partition_size)%subset2%part_size = 0

   ! Allocate new part indices
   allocate (self%parts(self%partition_size)%subset1%indices(self%max_part_size))
   allocate (self%parts(self%partition_size)%subset2%indices(self%max_part_size))

   ! Set new part index to partition current num part
   self%parts(self%partition_size)%index = self%partition_size
   self%parts(self%partition_size)%subset1%index = self%partition_size
   self%parts(self%partition_size)%subset2%index = self%partition_size

   ! Point new part max part size pointer to partition max part size
   self%parts(self%partition_size)%subset1%largest_part_size => self%largest_part_size
   self%parts(self%partition_size)%subset2%largest_part_size => self%largest_part_size

   ! Point new part map pointer to partition maps
   self%parts(self%partition_size)%subset1%partition_map => self%partition_map1
   self%parts(self%partition_size)%subset2%partition_map => self%partition_map2

   ! Return pointer to new part
   part => self%parts(self%partition_size)

end function

subroutine bipartition_print_parts(self)
   class(bipartition_type), intent(in) :: self
   integer :: h
   character(:), allocatable :: fmtstr

   do h = 1, self%partition_size
      fmtstr = "('{'" // repeat(',1x,i3', self%parts(h)%subset1%part_size) &
            // ",' }   {'" // repeat(',1x,i3', self%parts(h)%subset2%part_size) // ",' }')"
      write (stderr, fmtstr) &
            self%parts(h)%subset1%indices(:self%parts(h)%subset1%part_size), &
            self%parts(h)%subset2%indices(:self%parts(h)%subset2%part_size)
   end do

end subroutine

function bipartition_equality(self, other) result(equality)
   type(bipartition_type), intent(in) :: self, other
   ! Result variable
   logical :: equality
   ! Local variables
   integer :: h

   if (self%partition_size /= other%partition_size) then
      equality = .false.
      return
   end if

   do h = 1, self%partition_size
      if (self%parts(h)%subset1%part_size /= other%parts(h)%subset1%part_size) then
         equality = .false.
         return
      end if
      if (self%parts(h)%subset2%part_size /= other%parts(h)%subset2%part_size) then
         equality = .false.
         return
      end if
   end do

   equality = .true.

end function

end module
