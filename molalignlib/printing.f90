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

interface print_stats
   module procedure print_stats_rmsd
   module procedure print_stats_adjd
end interface

contains

subroutine print_stats_rmsd(results)
   type(registry_type), intent(in) :: results
   character(*), parameter :: line = repeat('-', 42)
   ! Local variables
   integer :: i
   type(record_type) :: record

   write (stdout, '(2x,a,3x,a,4x,a,5x,a,8x,a)') '#', 'Count', 'Steps', 'Rot-θ', 'RMSD'
   write (stdout, '(a)') line
   do i = 1, results%num_records
      record = results%records(i)
      write (stdout, '(i3,4x,i4,4x,f5.1,5x,f5.1,4x,f8.4)') &
         i, record%count, record%steps, record%rotang, sqrt(record%msd)
   end do
   write (stdout, '(a)') line
   write (stdout, '(a,1x,i0)') 'Random trials =', results%num_trials
   write (stdout, '(a,1x,i0)') 'Minimization steps =', results%total_steps
   if (results%overflow) then
      write (stdout, '(a,1x,i0)') 'Visited local minima >', results%num_records
   else
      write (stdout, '(a,1x,i0)') 'Visited local minima =', results%num_records
   end if
   flush(stdout)

end subroutine

subroutine print_stats_adjd(nrec, matches, avgsteps, avgrot, adjd, rmsd)
   integer, intent(in) :: nrec
   integer, dimension(:), intent(in) :: matches, adjd
   real(rk), dimension(:), intent(in) :: avgsteps, avgrot, rmsd
   character(*), parameter :: line = repeat('-', 49)
   integer :: irec

   write (stdout, '(2x,a,3x,a,4x,a,5x,a,6x,a,6x,a)') '#', 'Count', 'Steps', 'Rot-θ', 'Δadj', 'RMSD'
   write (stdout, '(a)') line
   do irec = 1, nrec
      write (stdout, '(i3,4x,i4,4x,f5.1,5x,f5.1,4x,i4,4x,f8.4)') &
         irec, matches(irec), avgsteps(irec), 90./asin(1.)*avgrot(irec), adjd(irec), rmsd(irec)
   end do
   write (stdout, '(a)') line
   flush(stdout)

end subroutine

subroutine print_final_stats(overflow, max_records, nrec, trials, nstep)
   logical, intent(in) :: overflow
   integer, intent(in) :: max_records, nrec, trials, nstep

   write (stdout, '(a,1x,i0)') 'Random trials =', trials
   write (stdout, '(a,1x,i0)') 'Minimization steps =', nstep
   if (overflow) then
      write (stdout, '(a,1x,i0)') 'Visited local minima >', max_records
   else
      write (stdout, '(a,1x,i0)') 'Visited local minima =', nrec
   end if
   flush(stdout)

end subroutine

end module
