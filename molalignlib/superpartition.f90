module superpartition
use stdio
use kinds
use partition
implicit none

real, parameter :: MAX_LOAD_FACTOR = 0.5

type, public :: superpartition_type
   integer :: superpartition_size
   integer :: max_superpartition_size
   integer :: max_partition_size
   integer :: max_part_size
   type(partition_type), pointer :: partitions(:)
contains
   procedure :: init => superpartition_init
   procedure :: new_partition => superpartition_new_partition
end type

type, public :: partitiondict_type
   type(partition_type), allocatable :: partitions(:)
   logical, allocatable :: occupied(:)
   integer :: num_slots
   integer :: num_occupied
contains
   procedure :: get_index
   procedure :: new_index
   procedure :: init => dict_init
   procedure :: reset => dict_reset
end type

type, public :: partitionpointer_type
   type(partition_type), pointer :: ptr
end type

interface operator (.in.)
   module procedure partition_in_dict
end interface

contains

subroutine superpartition_init(self, max_superpartition_size, max_partition_size, max_part_size)
   class(superpartition_type), intent(inout) :: self
   integer, intent(in) :: max_superpartition_size, max_partition_size, max_part_size

   self%superpartition_size = 0
   self%max_superpartition_size = max_superpartition_size
   self%max_partition_size = max_partition_size
   self%max_part_size = max_part_size
   allocate (self%partitions(max_superpartition_size))

end subroutine

function superpartition_new_partition(self) result(partition)
   class(superpartition_type), intent(inout) :: self
   ! Result variable
   type(partition_type), pointer :: partition

   self%superpartition_size = self%superpartition_size + 1
   call self%partitions(self%superpartition_size)%init(self%max_partition_size, self%max_part_size)

   ! Return pointer to new partition
   partition => self%partitions(self%superpartition_size)

end function

subroutine dict_init( this, min_dict_size)
   class(partitiondict_type), intent(inout) :: this
   integer, intent(in) :: min_dict_size
   
   this%num_slots = int(min_dict_size / MAX_LOAD_FACTOR)
   
   allocate(this%partitions(this%num_slots))
   allocate(this%occupied(this%num_slots))
   
   this%occupied = .false.
   this%num_occupied = 0

end subroutine dict_init

subroutine dict_reset(this)
   class(partitiondict_type), intent(inout) :: this

   this%occupied = .false.
   this%num_occupied = 0

end subroutine dict_reset

function new_index( this, partition) result(index)
   class(partitiondict_type), intent(inout) :: this
   type(partition_type), intent(in) :: partition
   integer :: index

   index = compute_hash(partition, this%num_slots) + 1

   do while (this%occupied(index))
      if (this%partitions(index) == partition) then
         error stop "Key already in hash table"
         return
      end if
      index = modulo(index, this%num_slots) + 1
   end do

   if (this%num_occupied >= this%num_slots*MAX_LOAD_FACTOR) then
      error stop "Hash table too full"
   end if

   this%partitions(index) = partition
   this%occupied(index) = .true.
   this%num_occupied = this%num_occupied + 1

end function new_index

function get_index( this, partition) result(index)
   class(partitiondict_type), intent(in) :: this
   type(partition_type), intent(in) :: partition
   integer :: index

   index = compute_hash(partition, this%num_slots) + 1

   do while (this%occupied(index))
      if (this%partitions(index) == partition) then
         return
      end if
      index = modulo(index, this%num_slots) + 1
   end do

   error stop "Key not found"

end function get_index

function partition_in_dict( partition, dict)
   class(partitiondict_type), intent(in) :: dict
   type(partition_type), intent(in) :: partition
   logical :: partition_in_dict
   integer :: index

   index = compute_hash(partition, dict%num_slots) + 1

   partition_in_dict = .false.
   do while (dict%occupied(index))
      if (dict%partitions(index) == partition) then
         partition_in_dict = .true.
         return
      end if
      index = modulo(index, dict%num_slots) + 1
   end do

end function partition_in_dict

function compute_hash(partition, num_slots) result(hash)
    type(partition_type), intent(in) :: partition
    integer, intent(in) :: num_slots  ! Number of slots in hash table
    integer :: hash
    integer :: h, i
    
    ! DJB2 initial value
    hash = 5381
    
    ! Hash the partition structure
    do h = 1, partition%partition_size
        ! Hash the part size first
        hash = modulo(hash * 33 + partition%parts(h)%part_size, num_slots)
        
        ! Hash each index in the part
        do i = 1, partition%parts(h)%part_size
            hash = modulo(hash * 33 + partition%parts(h)%indices(i), num_slots)
        end do
    end do
    
end function compute_hash

end module
