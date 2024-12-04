module multisetdict
use stdio
implicit none
private

public multisetdict_type
public operator (.in.)

real, parameter :: MAX_LOAD_FACTOR = 0.5

type, public :: multiset_type
   integer, allocatable :: items(:)
end type

type, public :: multisetdict_type
   type(multiset_type), allocatable :: multisets(:)
   logical, allocatable :: occupied(:)
   integer :: num_slots
   integer :: num_occupied
contains
   procedure :: get_index
   procedure :: new_index
   procedure :: initialize
   procedure :: reset
end type

interface operator (.in.)
   module procedure in_dict
end interface

contains

subroutine initialize( this, min_dict_size)
   class(multisetdict_type), intent(inout) :: this
   integer, intent(in) :: min_dict_size
   
   this%num_slots = int(min_dict_size / MAX_LOAD_FACTOR)
   
   allocate (this%occupied(this%num_slots))
   allocate (this%multisets(this%num_slots))
   
   this%occupied = .false.
   this%num_occupied = 0

end subroutine

subroutine reset(this)
   class(multisetdict_type), intent(inout) :: this

   this%occupied = .false.
   this%num_occupied = 0

end subroutine

function new_index( this, multiset) result(index)
   class(multisetdict_type), intent(inout) :: this
   integer, intent(in) :: multiset(:)
   integer :: index

   index = modulo(hash(multiset), this%num_slots) + 1

   do while (this%occupied(index))
      if (multiset_equality(this%multisets(index)%items, multiset)) then
         error stop "Key already in hash table"
      end if
      index = modulo(index, this%num_slots) + 1
   end do

   if (this%num_occupied >= this%num_slots*MAX_LOAD_FACTOR) then
      error stop "Hash table too full"
   end if

   this%multisets(index)%items = multiset
   this%occupied(index) = .true.
   this%num_occupied = this%num_occupied + 1

end function

function get_index( this, multiset) result(index)
   class(multisetdict_type), intent(in) :: this
   integer, intent(in) :: multiset(:)
   integer :: index

   index = modulo(hash(multiset), this%num_slots) + 1

   do while (this%occupied(index))
      if (multiset_equality(this%multisets(index)%items, multiset)) then
         return
      end if
      index = modulo(index, this%num_slots) + 1
   end do

   error stop "Key not found"

end function

function in_dict( multiset, dict)
   integer, intent(in) :: multiset(:)
   class(multisetdict_type), intent(in) :: dict
   logical :: in_dict
   integer :: index

   index = modulo(hash(multiset), dict%num_slots) + 1

   do while (dict%occupied(index))
      if (multiset_equality(dict%multisets(index)%items, multiset)) then
         in_dict = .true.
         return
      end if
      index = modulo(index, dict%num_slots) + 1
   end do

   in_dict = .false.

end function

integer function hash( multiset)
   integer, parameter :: HASH_CONSTANT = 5381
   integer, intent(in) :: multiset(:)
   integer :: i

   hash = 1
   do i = 1, size(multiset)
      hash = iand(hash * (HASH_CONSTANT + 2 * multiset(i)), 2**16 - 1)
   end do
   hash = hash / 2

end function

function multiset_equality( multiset1, multiset2) result(equality)
   integer, intent(in) :: multiset1(:), multiset2(:)
   logical :: equality
   integer :: i, j, matches
   
   ! Check if sizes are equal
   if (size(multiset1) /= size(multiset2)) then
      equality = .false.
      return
   end if

   ! Quick check for exact match first
   if (all(multiset1 == multiset2)) then
       equality = .true.
       return
   end if
   
   ! Early exit: Check if any value appears more times in one array
   ! than it does in the other
   do i = 1, size(multiset1)
       matches = 0
       do j = 1, size(multiset1)
           if (multiset1(i) == multiset2(j)) matches = matches + 1
           if (multiset1(i) == multiset1(j)) matches = matches - 1
       end do
       if (matches /= 0) then
           equality = .false.
           return
       end if
   end do

   equality = .true.

end function

end module
