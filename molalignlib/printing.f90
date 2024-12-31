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

module printing
use parameters
use registry

implicit none

contains

subroutine print_stats_rmsd(results)
   type(registry_type), intent(in) :: results
   ! Parameters
   character(*), parameter :: line = repeat('-', 42)
   ! Local variables
   type(record_type) :: record
   integer :: i

   write (stdout, '(2x,a,4x,a,4x,a,4x,a,7x,a)') '#', 'Count', 'Steps', 'Rot-θ', 'RMSD'
   write (stdout, '(a)') line
   do i = 1, results%num_records
      record = results%records(i)
      write (stdout, '(i3,4x,i4,4x,f5.1,5x,f5.1,4x,f8.4)') &
         i, record%count, record%steps, record%rotang, record%rmsd
   end do
   write (stdout, '(a)') line
   flush(stdout)

end subroutine

subroutine print_stats_adjd(results)
   type(registry_type), intent(in) :: results
   ! Parameters
   character(*), parameter :: line = repeat('-', 49)
   ! Local variables
   type(record_type) :: record
   integer :: i

   write (stdout, '(2x,a,4x,a,4x,a,4x,a,4x,a,6x,a)') '#', 'Count', 'Steps', 'Rot-θ', 'Δadj', 'RMSD'
   write (stdout, '(a)') line
   do i = 1, results%num_records
      record = results%records(i)
      write (stdout, '(i3,4x,i4,4x,f5.1,5x,f5.1,3x,i4,4x,f8.4)') &
         i, record%count, record%steps, record%rotang, record%adjd, record%rmsd
   end do
   write (stdout, '(a)') line
   flush(stdout)

end subroutine

subroutine print_final_stats(results)
   type(registry_type), intent(in) :: results

   write (stdout, '(a,1x,i0)') 'Random trials =', results%num_trials
   write (stdout, '(a,1x,i0)') 'Minimization steps =', results%total_steps
   if (results%overflow) then
      write (stdout, '(a,1x,i0)') 'Visited local minima >', results%num_records
   else
      write (stdout, '(a,1x,i0)') 'Visited local minima =', results%num_records
   end if
   flush(stdout)

end subroutine

end module
