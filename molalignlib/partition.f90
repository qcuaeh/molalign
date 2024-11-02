module partition
use stdio
use kinds
use sorting

implicit none

type, public :: intlist_type
   integer, allocatable :: n(:)
end type

type, public :: reallist_type
   real(rk), allocatable :: x(:)
end type

type, public :: nested_intlist_type
   type(intlist_type), allocatable :: s(:)
end type

type, public :: nested_reallist_type
   type(reallist_type), allocatable :: s(:)
end type

type :: boolmatrix_type
   logical, allocatable :: b(:, :)
end type

type :: intmatrix_type
   integer, allocatable :: n(:, :)
end type

type :: realmatrix_type
   real(rk), allocatable :: x(:, :)
end type

type, public :: atomlist_type
   integer, allocatable :: atomidcs(:)
end type

type, public :: partpointer_type
   type(part_type), pointer :: ptr => null()
end type

type, public :: subset_type
   integer :: index
   integer :: part_size
   integer, pointer :: indices(:)
   integer, pointer :: allocation(:)
   integer, pointer :: largest_part_size
   integer, pointer :: partition_map(:)
contains
   procedure :: add => subset_add
end type

type, public :: part_type
   integer :: index
   type(subset_type) :: subset1
   type(subset_type) :: subset2
end type

type, public :: partition_type
   integer :: partition_size
   integer :: max_partition_size
   integer :: max_part_size
   integer, pointer :: largest_part_size
   integer, pointer :: partition_map1(:)
   integer, pointer :: partition_map2(:)
   type(part_type), pointer :: parts(:)
contains
   procedure :: init => partition_init
   procedure :: new_part => partition_new_part
   procedure :: print_parts => partition_print_parts
end type

interface operator (==)
   module procedure partition_equality
end interface

contains

function partition_equality(self, other) result(equality)
   type(partition_type), intent(in) :: self, other
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

subroutine partition_init(self, max_partition_size, max_part_size)
   class(partition_type), intent(inout) :: self
   integer, intent(in) :: max_partition_size, max_part_size

   allocate (self%largest_part_size)
   allocate (self%parts(max_partition_size))
   allocate (self%partition_map1(max_part_size))
   allocate (self%partition_map2(max_part_size))

   self%partition_size = 0
   self%largest_part_size = 0
   self%max_partition_size = max_partition_size
   self%max_part_size = max_part_size

end subroutine

function partition_new_part(self) result(part)
   class(partition_type), intent(inout) :: self
   ! Result variable
   type(part_type), pointer :: part

   ! Increase partition size
   self%partition_size = self%partition_size + 1

   ! Initialize new part size
   self%parts(self%partition_size)%subset1%part_size = 0
   self%parts(self%partition_size)%subset2%part_size = 0

   ! Create new part allocation
   allocate (self%parts(self%partition_size)%subset1%allocation(self%max_part_size))
   allocate (self%parts(self%partition_size)%subset2%allocation(self%max_part_size))

   ! Point new part indices to part allocation with size 0
   self%parts(self%partition_size)%subset1%indices => self%parts(self%partition_size)%subset1%allocation(:0)
   self%parts(self%partition_size)%subset2%indices => self%parts(self%partition_size)%subset2%allocation(:0)

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

subroutine subset_add(self, element)
   class(subset_type), intent(inout) :: self
   integer, intent(in) :: element

   ! Increase part size
   self%part_size = self%part_size + 1

   ! Update max part size
   if (self%part_size > self%largest_part_size) then
      self%largest_part_size = self%part_size
   end if

   ! Update part pointer size
   self%indices => self%allocation(:self%part_size)

   ! Add element to part
   self%indices(self%part_size) = element

   ! Add index to part map
   self%partition_map(element) = self%index

end subroutine

subroutine partition_print_parts(self)
   class(partition_type), intent(in) :: self
   integer :: h
   character(:), allocatable :: fmtstr

   do h = 1, self%partition_size
      fmtstr = "('{'" // repeat(',1x,i3', self%parts(h)%subset1%part_size) &
            // ",' }   {'" // repeat(',1x,i3', self%parts(h)%subset2%part_size) // ",' }')"
      write (stderr, fmtstr) self%parts(h)%subset1%indices, self%parts(h)%subset2%indices
   end do

end subroutine

end module
