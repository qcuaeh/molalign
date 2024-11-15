program test_dict_part2b
   use dict_mod
   use, intrinsic :: iso_fortran_env, only: int64
   implicit none

   integer, parameter :: KEY_LENGTH = 6
   integer, parameter :: MAX_KEY_VALUE = 1000
   integer, parameter :: PROGRESS_INTERVAL = 100000
   integer :: NUM_KEYS = 100000 ! Will be set from command line or default
   integer, parameter :: NUM_REPETITIONS = 10000

   type(dict) :: d
   integer :: key(KEY_LENGTH)
   integer, allocatable :: keys(:,:)
   integer :: i, j, value
   logical :: has_key
   integer(int64) :: start_time, end_time, count_rate, count_max
   real :: insertion_time, retrieval_time, lookup_time
   integer :: total_operations, current_operation
   character(len=32) :: arg

   ! Get number of keys from command line or use default
   if (command_argument_count() > 0) then
      call get_command_argument(1, arg)
      read(arg, *) NUM_KEYS
   end if

   print '(A,I0)', "Number of keys to be generated: ", NUM_KEYS

   ! Create dictionary
   call d%init(NUM_KEYS)
   print '(A,I0)', "Actual dictionary size: ", d%num_slots

   ! Allocate arrays for keys
   allocate(keys(KEY_LENGTH, NUM_KEYS))

   total_operations = 3 * NUM_REPETITIONS * NUM_KEYS
   current_operation = 0

   insertion_time = 0.0
   retrieval_time = 0.0
   lookup_time = 0.0

   ! Print initial progress
   write(*, '(A)', advance='no') "Progress:   0%"
   flush(6)

   do j = 1, NUM_REPETITIONS
      call d%reset()  ! Reset the dictionary for each repetition

      ! Generate all keys for this repetition
      do i = 1, NUM_KEYS
         call generate_random_key(keys(:,i), MAX_KEY_VALUE)
      end do

      ! Addition loop
      call system_clock(start_time, count_rate, count_max)
      do i = 1, NUM_KEYS
         value = i
         call d%add(keys(:, i), value)

         current_operation = current_operation + 1
         if (mod(current_operation, PROGRESS_INTERVAL) == 0) then
            write(*, '(A,I3,A)', advance='no') char(13)//"Progress: ", &
               int(real(current_operation) / real(total_operations) * 100), "%"
            flush(6)
         end if
      end do
      call system_clock(end_time)
      insertion_time = insertion_time + real(end_time - start_time) / real(count_rate)

      ! Retrieval loop
      call system_clock(start_time, count_rate, count_max)
      do i = 1, NUM_KEYS
         value = d%get(keys(:, i))

         current_operation = current_operation + 1
         if (mod(current_operation, PROGRESS_INTERVAL) == 0) then
            write(*, '(A,I3,A)', advance='no') char(13)//"Progress: ", &
               int(real(current_operation) / real(total_operations) * 100), "%"
            flush(6)
         end if
      end do
      call system_clock(end_time)
      retrieval_time = retrieval_time + real(end_time - start_time) / real(count_rate)

      ! Lookup loop
      call system_clock(start_time, count_rate, count_max)
      do i = 1, NUM_KEYS
         call generate_random_key(key, MAX_KEY_VALUE)
         has_key = d%has(key)

         current_operation = current_operation + 1
         if (mod(current_operation, PROGRESS_INTERVAL) == 0) then
            write(*, '(A,I3,A)', advance='no') char(13)//"Progress: ", &
               int(real(current_operation) / real(total_operations) * 100), "%"
            flush(6)
         end if
      end do
      call system_clock(end_time)
      lookup_time = lookup_time + real(end_time - start_time) / real(count_rate)

   end do

   ! Print final progress
   write(*, '(A)') char(13)//"Progress: 100%"

   print '(A,I0)', "Repetitions: ", NUM_REPETITIONS
   print '(A,F10.3,A)', "Total insertion time: ", insertion_time, " seconds"
   print '(A,F10.3,A)', "Total retrieval time: ", retrieval_time, " seconds"
   print '(A,F10.3,A)', "Total lookup time: ", lookup_time, " seconds"

contains

   subroutine generate_random_key(key, max_value)
      integer, intent(out) :: key(:)        ! Array to store the key
      integer, intent(in) :: max_value      ! Maximum allowed key value
      integer :: i
      real :: r

      ! Generate random values for the key
      do i = 1, size(key)
         call random_number(r)
         key(i) = int(r * max_value) + 1
      end do
!      print *, key
   end subroutine generate_random_key

end program test_dict_part2b
