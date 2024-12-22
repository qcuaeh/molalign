! MolAlignLib
! Copyright (C) 2022 José M. Vásquez

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module registry
use parameters
use alignment

implicit none

type :: record_type
   real(rk) :: msd
   real(rk) :: steps
   real(rk) :: rotang
   integer :: count
   integer, allocatable :: atomperm(:)
end type

type :: registry_type
   logical :: overflow
   integer :: num_trials
   integer :: num_records
   integer :: total_steps
   type(record_type), allocatable :: records(:)
contains
   procedure :: initialize => registry_initialize
   procedure :: push => registry_push
end type

contains

subroutine registry_initialize(self, max_records)
   class(registry_type), intent(inout) :: self
   integer, intent(in) :: max_records
   integer :: i

   if (max_records < 1) then
      error stop 'max_records < 1'
   end if

   self%num_records = 0
   self%total_steps = 0
   self%overflow = .false.

   allocate (self%records(max_records))
   self%records%count = 0
   self%records%msd = huge(self%records(1)%msd)

end subroutine

subroutine registry_push(self, atomperm, steps, rotang, coords1, coords2)
   class(registry_type), intent(inout) :: self
   integer, intent(in) :: atomperm(:)
   integer, intent(in) :: steps
   real(rk), intent(in) :: rotang
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   ! Local variables
   integer :: i, j
   real(rk) :: msd

   self%total_steps = self%total_steps + steps

   do i = 1, self%num_records
      if (allocated(self%records(i)%atomperm)) then
         if (all(atomperm == self%records(i)%atomperm)) then
            self%records(i)%count = self%records(i)%count + 1
            self%records(i)%steps = self%records(i)%steps + (steps - self%records(i)%steps) / self%records(i)%count
            self%records(i)%rotang = self%records(i)%rotang + (rotang - self%records(i)%rotang) / self%records(i)%count
            return
         end if
      end if
   end do

   msd = squaredist(atomperm, coords1, coords2)
   do i = 1, size(self%records)
      if (msd < self%records(i)%msd) then
         do j = size(self%records), i + 1, -1
            self%records(j) = self%records(j - 1)
         end do
         self%records(i)%atomperm = atomperm
         self%records(i)%count = 1
         self%records(i)%msd = msd
         self%records(i)%steps = steps
         self%records(i)%rotang = rotang
         exit
      end if
   end do

   if (.not. self%overflow) then
      if (self%num_records < size(self%records)) then
         self%num_records = self%num_records + 1
      else
         self%overflow = .true.
      end if
   end if

end subroutine

end module
