module bipartition
use stdio
use kinds
use part
implicit none

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

type, public :: bipartpointer_type
   type(bipart_type), pointer :: ptr => null()
end type

interface operator (==)
   module procedure bipartition_equality
end interface

contains

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

   write (stderr, *)
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
