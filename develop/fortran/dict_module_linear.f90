module dict_mod
   use, intrinsic :: iso_fortran_env, only: int32, int64
   implicit none
   private

   public :: dict

   integer, parameter :: MAX_KEY_LEN = 16
   real, parameter :: MAX_LOAD_FACTOR = 0.5

   type, public :: dict
      integer, allocatable :: keys(:,:)
      integer, allocatable :: key_lengths(:)
      integer, allocatable :: values(:)
      logical, allocatable :: occupied(:)
      integer :: num_slots
      integer :: num_occupied
   contains
      procedure :: add => dict_add
      procedure :: get => dict_get
      procedure :: has => dict_has
      procedure :: init => dict_init
      procedure :: reset => dict_reset
   end type dict

contains

   subroutine dict_init(this, min_dict_size)
      class(dict), intent(inout) :: this
      integer, intent(in) :: min_dict_size
      
      this%num_slots = int(min_dict_size / MAX_LOAD_FACTOR)
      
      allocate(this%keys(MAX_KEY_LEN, this%num_slots))
      allocate(this%key_lengths(this%num_slots))
      allocate(this%values(this%num_slots))
      allocate(this%occupied(this%num_slots))
      
      this%occupied = .false.
      this%num_occupied = 0

   end subroutine dict_init

   subroutine dict_reset(this)
      class(dict), intent(inout) :: this
      this%occupied = .false.
      this%num_occupied = 0

   end subroutine dict_reset

   subroutine dict_add(this, key, value)
      class(dict), intent(inout) :: this
      integer, intent(in) :: key(:)
      integer, intent(in) :: value
      integer :: index

      index = modulo(hash(key), this%num_slots) + 1

      do while (this%occupied(index))
         if (same_keys(this%keys(:, index), this%key_lengths(index), key, size(key))) then
            this%values(index) = value
            return
         end if
         index = modulo(index, this%num_slots) + 1
      end do

      if (this%num_occupied >= this%num_slots*MAX_LOAD_FACTOR) then
         error stop "Dictionary is too full"
      end if

      this%keys(:size(key), index) = key
      this%key_lengths(index) = size(key)
      this%values(index) = value
      this%occupied(index) = .true.
      this%num_occupied = this%num_occupied + 1

   end subroutine dict_add

   function dict_get(this, key) result(value)
      class(dict), intent(in) :: this
      integer, intent(in) :: key(:)
      integer :: value
      integer :: index

      index = modulo(hash(key), this%num_slots) + 1

      do while (this%occupied(index))
         if (same_keys(this%keys(:, index), this%key_lengths(index), key, size(key))) then
            value = this%values(index)
            return
         end if
         index = modulo(index, this%num_slots) + 1
      end do

      error stop "Key not found"

   end function dict_get

   function dict_has(this, key) result(exists)
      class(dict), intent(in) :: this
      integer, intent(in) :: key(:)
      logical :: exists
      integer :: index

      index = modulo(hash(key), this%num_slots) + 1

      exists = .false.
      do while (this%occupied(index))
         if (same_keys(this%keys(:, index), this%key_lengths(index), key, size(key))) then
            exists = .true.
            return
         end if
         index = modulo(index, this%num_slots) + 1
      end do

   end function dict_has

   function hash(key)
      integer, parameter :: HASH_CONSTANT = 5381
      integer, intent(in) :: key(:)
      integer :: hash
      integer :: i

      hash = 1
      do i = 1, size(key)
         hash = iand(hash * (HASH_CONSTANT + 2 * key(i)), 2**16 - 1)
      end do
      hash = hash / 2

   end function hash

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

   end function same_keys

end module dict_mod
