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

implicit none

type :: record_type
   integer :: count
   integer :: adjd
   real(rk) :: rmsd
   real(rk) :: steps
   real(rk) :: rotang
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
   procedure :: push_adjd => registry_push_adjd
   procedure :: push_rmsd => registry_push_rmsd
end type

contains

subroutine registry_initialize(self, max_records)
   class(registry_type), intent(inout) :: self
   integer, intent(in) :: max_records

   if (max_records < 1) then
      error stop 'max_records < 1'
   end if

   self%num_records = 0
   self%total_steps = 0
   self%overflow = .false.

   allocate (self%records(max_records))

   self%records%count = 0
   self%records%rmsd = huge(self%records(1)%rmsd)
   self%records%adjd = huge(self%records(1)%adjd)

end subroutine

subroutine registry_push_rmsd(self, atomperm, steps, rotang, rmsd)
   class(registry_type), target, intent(inout) :: self
   integer, intent(in) :: atomperm(:)
   integer, intent(in) :: steps
   real(rk), intent(in) :: rotang, rmsd
   ! Local variables
   type(record_type), pointer :: record
   integer :: i, j

   self%total_steps = self%total_steps + steps

   do i = 1, self%num_records
      record => self%records(i)
      if (allocated(record%atomperm)) then
         if (all(atomperm == record%atomperm)) then
            record%count = record%count + 1
            record%steps = record%steps + (steps - record%steps) / record%count
            record%rotang = record%rotang + (rotang - record%rotang) / record%count
            return
         end if
      end if
   end do

   do i = 1, size(self%records)
      record => self%records(i)
      if (rmsd < record%rmsd) then
         do j = size(self%records), i + 1, -1
            self%records(j) = self%records(j - 1)
         end do
         record%atomperm = atomperm
         record%count = 1
         record%rmsd = rmsd
         record%steps = steps
         record%rotang = rotang
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

subroutine registry_push_adjd(self, atomperm, steps, rotang, adjd)
   class(registry_type), target, intent(inout) :: self
   integer, intent(in) :: atomperm(:)
   integer, intent(in) :: steps, adjd
   real(rk), intent(in) :: rotang
   ! Local variables
   type(record_type), pointer :: record
   integer :: i, j

   self%total_steps = self%total_steps + steps

   do i = 1, self%num_records
      record => self%records(i)
      if (allocated(record%atomperm)) then
         if (adjd == record%adjd) then
            record%count = record%count + 1
            record%steps = record%steps + (steps - record%steps) / record%count
            record%rotang = record%rotang + (rotang - record%rotang) / record%count
            return
         end if
      end if
   end do

   do i = 1, size(self%records)
      record => self%records(i)
      if (adjd < record%adjd) then
         do j = size(self%records), i + 1, -1
            self%records(j) = self%records(j - 1)
         end do
         record%atomperm = atomperm
         record%count = 1
         record%adjd = adjd
         record%steps = steps
         record%rotang = rotang
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
