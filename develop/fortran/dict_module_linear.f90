module dict_mod
   use, intrinsic :: iso_fortran_env, only: int32, int64
   implicit none
   private

   public :: dict

   real, parameter :: MAX_LOAD_FACTOR = 0.7

   type, public :: dict
      integer, allocatable :: keys(:,:)
      integer, allocatable :: key_lengths(:)
      integer, allocatable :: values(:)
      logical, allocatable :: occupied(:)
      integer :: num_slots
      integer :: num_occupied
      integer :: max_key_len
      integer :: max_key_val
      integer :: max_occupied
   contains
      procedure :: add => dict_add
      procedure :: get => dict_get
      procedure :: has => dict_has
      procedure :: init => dict_init
      procedure :: reset => dict_reset
   end type dict

contains

   subroutine dict_init(this, max_key_len, max_key_val, max_occupied)
      class(dict), intent(inout) :: this
      integer, intent(in) :: max_occupied
      integer, intent(in) :: max_key_len
      integer, intent(in) :: max_key_val
      
      this%num_slots = int(max_occupied / MAX_LOAD_FACTOR)
      this%max_occupied = max_occupied
      this%max_key_len = max_key_len
      this%max_key_val = max_key_val
      
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
      integer, intent(in) :: value
      integer :: hash, index

      hash = compute_hash(key)
      index = modulo(hash, this%num_slots) + 1

      do while (this%occupied(index))
         if (same_keys(this%keys(:, index), this%key_lengths(index), key, size(key), this%max_key_val)) then
            this%values(index) = value
            return
         end if
         index = modulo(index, this%num_slots) + 1
      end do

      if (this%num_occupied >= this%max_occupied) then
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
      integer :: hash, index

      hash = compute_hash(key)
      index = modulo(hash, this%num_slots) + 1

      do while (this%occupied(index))
         if (same_keys(this%keys(:, index), this%key_lengths(index), key, size(key), this%max_key_val)) then
            value = this%values(index)
            return
         end if
         index = modulo(index, this%num_slots) + 1
         if (index == modulo(hash, this%num_slots) + 1) then
            error stop "Key not found"
         end if
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

      do while (this%occupied(index))
         if (same_keys(this%keys(:, index), this%key_lengths(index), key, size(key), this%max_key_val)) then
            exists = .true.
            return
         end if
         index = modulo(index, this%num_slots) + 1
         if (index == modulo(hash, this%num_slots) + 1) then
            exists = .false.
            return
         end if
      end do

      exists = .false.

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

   function same_keys(key1, key1_len, key2, key2_len, max_key_val)
      integer, intent(in) :: key1(:), key2(:)
      integer, intent(in) :: key1_len, key2_len
      integer, intent(in) :: max_key_val
      logical :: same_keys
      
      ! Check if sizes are equal
      if (key1_len /= key2_len) then
         same_keys = .false.
         return
      end if

      ! Handle special cases based on size
      select case(key1_len)
         case(0)
            same_keys = .true.
            return

         case(1)
            same_keys = (key1(1) == key2(1))
            return

         case(2)
            same_keys = (key1(1) == key2(1) .and. key1(2) == key2(2)) .or. &
                        (key1(1) == key2(2) .and. key1(2) == key2(1))
            return

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
            return
            
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
            return
            
         case default
            block
               integer :: count1(max_key_val), count2(max_key_val)
               integer :: i
               
               ! Initialize count arrays to zero
               count1 = 0
               count2 = 0
               
               ! Count occurrences in both arrays
               do i = 1, key1_len
                  count1(key1(i)) = count1(key1(i)) + 1
                  count2(key2(i)) = count2(key2(i)) + 1
               end do
               
               ! Compare count arrays
               same_keys = all(count1 == count2)
            end block
      end select

   end function same_keys

end module dict_mod
