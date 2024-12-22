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

module strutils
use parameters

implicit none

private
public int
public str
public lowercase
public uppercase

interface int
   module procedure str_int
end interface

interface str
   module procedure int_str
   module procedure real_str
end interface

contains

function str_int(x) result(n)
   character(*), intent(in) :: x
   integer :: n
   read (x, *) n
end function

function int_str(n) result(str)
   integer, intent(in) :: n
   character(:), allocatable :: str
   integer :: l
   l = floor(log10(real(max(abs(n), 1)))) + 1
   if (n < 0) l = l + 1
   allocate (character(l) :: str)
   write (str, '(i0)') n
end function

function real_str(x, n) result(str)
   real(rk), intent(in) :: x
   integer, intent(in) :: n
   integer :: l
   character(32) format
   character(:), allocatable :: str
   l = floor(log10(max(abs(x), 1.0_rk))) + n + 2
   if (x < 0) l = l + 1
   allocate (character(l) :: str)
   write (format, '(a,i0,a,i0,a)') '(f', l, '.', n, ')'
   write (str, format) x
end function

function lowercase(str)
   character(*), intent(in) :: str
   character(26), parameter :: lowerchars = 'abcdefghijklmnopqrstuvwxyz'
   character(26), parameter :: upperchars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   character(len(str)) :: lowercase
   integer :: i, j
   lowercase = str
   do j = 1, len(str)
      i = index(upperchars, str(j:j))
      if (i > 0) lowercase(j:j) = lowerchars(i:i)
   end do
end function

function uppercase(str)
   character(*), intent(in) :: str
   character(26), parameter :: lowerchars = 'abcdefghijklmnopqrstuvwxyz'
   character(26), parameter :: upperchars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   character(len(str)) :: uppercase
   integer :: i, j
   uppercase = str
   do j = 1, len(str)
      i = index(lowerchars, str(j:j))
      if (i > 0) uppercase(j:j) = upperchars(i:i)
   end do
end function

end module
