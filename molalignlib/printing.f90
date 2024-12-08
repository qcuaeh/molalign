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
use stdio
use kinds
use flags

implicit none

character(*), parameter :: line1 = repeat('-', 41)
character(*), parameter :: line2 = repeat('-', 48)

interface print_stats
   module procedure print_stats_rmsd
   module procedure print_stats_adjd
end interface

contains

subroutine print_stats_rmsd(nrec, matches, avgsteps, avgrot, rmsd)
   integer, intent(in) :: nrec
   integer, dimension(:), intent(in) :: matches
   real(rk), dimension(:), intent(in) :: avgsteps, avgrot, rmsd
   integer :: irec
   write (stderr, '(a,3x,a,4x,a,5x,a,8x,a)') 'Map', 'Count', 'Steps', 'Rot.', 'RMSD'
   write (stderr, '(a)') line1
   do irec = 1, nrec
      write (stderr, '(i3,4x,i4,4x,f5.1,5x,f5.1,3x,f8.4)') &
         irec, matches(irec), avgsteps(irec), 90./asin(1.)*avgrot(irec), rmsd(irec)
   end do
   write (stderr, '(a)') line1
   flush(stderr)
end subroutine

subroutine print_stats_adjd(nrec, matches, avgsteps, avgrot, adjd, rmsd)
   integer, intent(in) :: nrec
   integer, dimension(:), intent(in) :: matches, adjd
   real(rk), dimension(:), intent(in) :: avgsteps, avgrot, rmsd
   integer :: irec
   write (stderr, '(a,3x,a,4x,a,5x,a,6x,a,5x,a)') 'Map', 'Count', 'Steps', 'Rot.', 'Δadj', 'RMSD'
   write (stderr, '(a)') line2
   do irec = 1, nrec
      write (stderr, '(i3,4x,i4,4x,f5.1,5x,f5.1,3x,i4,3x,f8.4)') &
         irec, matches(irec), avgsteps(irec), 90./asin(1.)*avgrot(irec), adjd(irec), rmsd(irec)
   end do
   write (stderr, '(a)') line2
   flush(stderr)
end subroutine

subroutine print_final_stats(overflow, maxrec, nrec, ntrial, nstep)
   logical, intent(in) :: overflow
   integer, intent(in) :: maxrec, nrec, ntrial, nstep
   write (stderr, '(a,1x,i0)') 'Random trials =', ntrial
   write (stderr, '(a,1x,i0)') 'Minimization steps =', nstep
   if (overflow) then
      write (stderr, '(a,1x,i0)') 'Visited local minima >', maxrec
   else
      write (stderr, '(a,1x,i0)') 'Visited local minima =', nrec
   end if
   flush(stderr)
end subroutine

end module
