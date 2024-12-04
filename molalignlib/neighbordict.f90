module neighbordict
use stdio
use bounds
implicit none
private

public neighbordict_type
public operator (.in.)

real, parameter :: MAX_LOAD_FACTOR = 0.5

type, public :: neighbordict_type
   integer, allocatable :: keys(:,:)
   integer, allocatable :: key_lengths(:)
   logical, allocatable :: occupied(:)
   integer :: num_slots
   integer :: num_occupied
contains
   procedure :: get_index
   procedure :: new_index
   procedure :: initialize
   procedure :: reset
end type neighbordict_type

interface operator (.in.)
   module procedure key_in_dict
end interface

contains

subroutine initialize( this, min_dict_size)
   class(neighbordict_type), intent(inout) :: this
   integer, intent(in) :: min_dict_size
   
   this%num_slots = int(min_dict_size / MAX_LOAD_FACTOR)
   
   allocate(this%keys(max_coord_num, this%num_slots))
   allocate(this%key_lengths(this%num_slots))
   allocate(this%occupied(this%num_slots))
   
   this%occupied = .false.
   this%num_occupied = 0

end subroutine

subroutine reset(this)
   class(neighbordict_type), intent(inout) :: this

   this%occupied = .false.
   this%num_occupied = 0

end subroutine

function new_index( this, key) result(index)
   class(neighbordict_type), intent(inout) :: this
   integer, intent(in) :: key(:)
   integer :: index

   index = modulo(hash(key), this%num_slots) + 1

   do while (this%occupied(index))
      if (same_keys(this%keys(:, index), this%key_lengths(index), key, size(key))) then
         error stop "Key already in hash table"
      end if
      index = modulo(index, this%num_slots) + 1
   end do

   if (this%num_occupied >= this%num_slots*MAX_LOAD_FACTOR) then
      error stop "Hash table too full"
   end if

   this%keys(:size(key), index) = key
   this%key_lengths(index) = size(key)
   this%occupied(index) = .true.
   this%num_occupied = this%num_occupied + 1

end function

function get_index( this, key) result(index)
   class(neighbordict_type), intent(in) :: this
   integer, intent(in) :: key(:)
   integer :: index

   index = modulo(hash(key), this%num_slots) + 1

   do while (this%occupied(index))
      if (same_keys(this%keys(:, index), this%key_lengths(index), key, size(key))) then
         return
      end if
      index = modulo(index, this%num_slots) + 1
   end do

   error stop "Key not found"

end function

function key_in_dict( key, dict)
   integer, intent(in) :: key(:)
   class(neighbordict_type), intent(in) :: dict
   logical :: key_in_dict
   integer :: index

   index = modulo(hash(key), dict%num_slots) + 1

   do while (dict%occupied(index))
      if (same_keys(dict%keys(:, index), dict%key_lengths(index), key, size(key))) then
         key_in_dict = .true.
         return
      end if
      index = modulo(index, dict%num_slots) + 1
   end do

   key_in_dict = .false.

end function

integer function hash(key)
   integer, parameter :: HASH_CONSTANT = 5381
   integer, intent(in) :: key(:)
   integer :: i

   hash = 1
   do i = 1, size(key)
      hash = iand(hash * (HASH_CONSTANT + 2 * key(i)), 2**16 - 1)
   end do
   hash = hash / 2

end function

function same_keys(key1, key1_len, key2, key2_len)
   integer, intent(in) :: key1(:), key2(:)
   integer, intent(in) :: key1_len, key2_len
   logical :: same_keys
   integer :: i, j, matches
   
   ! Check if sizes are equal
   if (key1_len /= key2_len) then
      same_keys = .false.
      return
   end if

   ! Quick check for exact match first
   if (all(key1(1:key1_len) == key2(1:key1_len))) then
       same_keys = .true.
       return
   end if
   
   ! Early exit: Check if any value appears more times in one array
   ! than it does in the other
   do i = 1, key1_len
       matches = 0
       do j = 1, key1_len
           if (key1(i) == key2(j)) matches = matches + 1
           if (key1(i) == key1(j)) matches = matches - 1
       end do
       if (matches /= 0) then
           same_keys = .false.
           return
       end if
   end do

   same_keys = .true.

end function

end module
