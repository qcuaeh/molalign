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

module chemutils
use stdio
use kinds
use chemdata
use strutils

implicit none

private
public elsym2num
public split_tag

contains

function elsym2num(symbol) result(z)
   character(*), intent(in) :: symbol
   integer :: z

   do z = 1, num_elems
      if (uppercase(symbol) == uppercase(element_symbols(z))) then
         return
      end if
   end do

   select case (uppercase(symbol))
   case ('LJ')
      z = 1001
   case default
      z = -1
   end select

end function

subroutine split_tag( atom_tag, symbol, label)
   character(*), intent(in) :: atom_tag
   character(:), allocatable, intent(out) :: symbol, label
   integer :: m, n

   n = len_trim(atom_tag)
   m = verify(uppercase(trim(atom_tag)), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')

   if (m == 0) then
      label = '0'
      symbol = atom_tag
   else if (verify(atom_tag(m:n), '1234567890') == 0) then
      read (atom_tag(m:n), *) label
      symbol = atom_tag(1:m-1)
   else
      write (stderr, '(a,1x,a)') 'Invalid atomic symbol/label:', atom_tag
      stop
   end if

end subroutine

end module
