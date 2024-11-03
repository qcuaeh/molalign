module dict_mod
   use, intrinsic :: iso_fortran_env, only: int32, int64
   implicit none
   private

   public :: dict

   real, parameter :: MAX_LOAD_FACTOR = 0.5

   type, public :: dict
      integer, allocatable :: keys(:,:)
      integer, allocatable :: key_lengths(:)
      integer, allocatable :: values(:)
      logical, allocatable :: occupied(:)
      integer :: num_slots
      integer :: num_occupied
      integer :: max_key_len
   contains
      procedure :: add => dict_add
      procedure :: get => dict_get
      procedure :: has => dict_has
      procedure :: init => dict_init
      procedure :: reset => dict_reset
   end type dict

contains

   subroutine dict_init(this, min_dict_size, max_key_len)
      class(dict), intent(inout) :: this
      integer, intent(in) :: min_dict_size
      integer, intent(in) :: max_key_len
      
      this%num_slots = int(min_dict_size / MAX_LOAD_FACTOR)
      this%max_key_len = max_key_len
      
      allocate(this%keys(max_key_len, this%num_slots))
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
      integer :: hash, index
      integer, intent(in) :: value

      hash = compute_hash(key)
      index = modulo(hash, this%num_slots) + 1

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
      integer :: hash, index
      integer :: value

      hash = compute_hash(key)
      index = modulo(hash, this%num_slots) + 1

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
      integer :: hash, index

      hash = compute_hash(key)
      index = modulo(hash, this%num_slots) + 1

      exists = .false.
      do while (this%occupied(index))
         if (same_keys(this%keys(:, index), this%key_lengths(index), key, size(key))) then
            exists = .true.
            return
         end if
         index = modulo(index, this%num_slots) + 1
      end do

   end function dict_has

   function compute_hash(key) result(hash)
      integer, parameter :: HASH_CONSTANT = 1779033703
      integer, intent(in) :: key(:)
      integer :: hash
      integer :: i

      hash = 1
      do i = 1, size(key)
         hash = hash * (HASH_CONSTANT + 2 * key(i))
      end do
      hash = hash / 2

   end function compute_hash

   function same_keys(key1, key1_len, key2, key2_len)
      integer, intent(in) :: key1(:), key2(:)
      integer, intent(in) :: key1_len, key2_len
      logical :: same_keys
      
      ! Check if sizes are equal
      if (key1_len /= key2_len) then
         same_keys = .false.
         return
      end if

      ! Handle special cases based on size
      select case(key1_len)

      case(1)

         same_keys = (key1(1) == key2(1))

      case(2)

         same_keys = (key1(1) == key2(1) .and. key1(2) == key2(2)) .or. &
                     (key1(1) == key2(2) .and. key1(2) == key2(1))

      case(3)

      block
         integer :: k1_1, k1_2, k1_3
         integer :: k2_1, k2_2, k2_3
         
         ! Cache values
         k1_1 = key1(1); k1_2 = key1(2); k1_3 = key1(3)
         k2_1 = key2(1); k2_2 = key2(2); k2_3 = key2(3)
         
         ! Check if it's already in order
         if (k1_1 == k2_1) then
            if (k1_2 == k2_2) then
               same_keys = (k1_3 == k2_3)
               return
            else if (k1_2 == k2_3) then
               same_keys = (k1_3 == k2_2)
               return
            end if
         ! Check rotations
         else if (k1_1 == k2_2) then
            if (k1_2 == k2_1) then
               same_keys = (k1_3 == k2_3)
               return
            else if (k1_2 == k2_3) then
               same_keys = (k1_3 == k2_1)
               return
            end if
         else if (k1_1 == k2_3) then
            if (k1_2 == k2_1) then
               same_keys = (k1_3 == k2_2)
               return
            else if (k1_2 == k2_2) then
               same_keys = (k1_3 == k2_1)
               return
            end if
         end if
         same_keys = .false.
      end block

      case(4)

      block
         integer :: k1_1, k1_2, k1_3, k1_4
         integer :: k2_1, k2_2, k2_3, k2_4
         
         ! Cache values
         k1_1 = key1(1); k1_2 = key1(2); k1_3 = key1(3); k1_4 = key1(4)
         k2_1 = key2(1); k2_2 = key2(2); k2_3 = key2(3); k2_4 = key2(4)
         
         ! Check if already in order
         if (k1_1 == k2_1) then
            if (k1_2 == k2_2) then
               if (k1_3 == k2_3) then
                  same_keys = (k1_4 == k2_4)
                  return
               else if (k1_3 == k2_4) then
                  same_keys = (k1_4 == k2_3)
                  return
               end if
            else if (k1_2 == k2_3) then
               if (k1_3 == k2_2) then
                  same_keys = (k1_4 == k2_4)
                  return
               else if (k1_3 == k2_4) then
                  same_keys = (k1_4 == k2_2)
                  return
               end if
            else if (k1_2 == k2_4) then
               if (k1_3 == k2_2) then
                  same_keys = (k1_4 == k2_3)
                  return
               else if (k1_3 == k2_3) then
                  same_keys = (k1_4 == k2_2)
                  return
               end if
            end if
         ! k1_1 matches k2_2
         else if (k1_1 == k2_2) then
            if (k1_2 == k2_1) then
               if (k1_3 == k2_3) then
                  same_keys = (k1_4 == k2_4)
                  return
               else if (k1_3 == k2_4) then
                  same_keys = (k1_4 == k2_3)
                  return
               end if
            else if (k1_2 == k2_3) then
               if (k1_3 == k2_1) then
                  same_keys = (k1_4 == k2_4)
                  return
               else if (k1_3 == k2_4) then
                  same_keys = (k1_4 == k2_1)
                  return
               end if
            else if (k1_2 == k2_4) then
               if (k1_3 == k2_1) then
                  same_keys = (k1_4 == k2_3)
                  return
               else if (k1_3 == k2_3) then
                  same_keys = (k1_4 == k2_1)
                  return
               end if
            end if
         ! k1_1 matches k2_3
         else if (k1_1 == k2_3) then
            if (k1_2 == k2_1) then
               if (k1_3 == k2_2) then
                  same_keys = (k1_4 == k2_4)
                  return
               else if (k1_3 == k2_4) then
                  same_keys = (k1_4 == k2_2)
                  return
               end if
            else if (k1_2 == k2_2) then
               if (k1_3 == k2_1) then
                  same_keys = (k1_4 == k2_4)
                  return
               else if (k1_3 == k2_4) then
                  same_keys = (k1_4 == k2_1)
                  return
               end if
            else if (k1_2 == k2_4) then
               if (k1_3 == k2_1) then
                  same_keys = (k1_4 == k2_2)
                  return
               else if (k1_3 == k2_2) then
                  same_keys = (k1_4 == k2_1)
                  return
               end if
            end if
         ! k1_1 matches k2_4
         else if (k1_1 == k2_4) then
            if (k1_2 == k2_1) then
               if (k1_3 == k2_2) then
                  same_keys = (k1_4 == k2_3)
                  return
               else if (k1_3 == k2_3) then
                  same_keys = (k1_4 == k2_2)
                  return
               end if
            else if (k1_2 == k2_2) then
               if (k1_3 == k2_1) then
                  same_keys = (k1_4 == k2_3)
                  return
               else if (k1_3 == k2_3) then
                  same_keys = (k1_4 == k2_1)
                  return
               end if
            else if (k1_2 == k2_3) then
               if (k1_3 == k2_1) then
                  same_keys = (k1_4 == k2_2)
                  return
               else if (k1_3 == k2_2) then
                  same_keys = (k1_4 == k2_1)
                  return
               end if
            end if
         end if
         same_keys = .false.
      end block

      case default

      block
         integer :: i
         same_keys = .true.
         do i = 1, key1_len
            if (count(key1(:key1_len) == key1(i)) /= count(key2(:key2_len) == key1(i))) then
               same_keys = .false.
               return
            end if
         end do
      end block

      end select

   end function same_keys

end module dict_mod
