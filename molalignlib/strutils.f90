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
use kinds
use stdio

implicit none

private
public lowercase
public uppercase
public intstr
public realstr
public get_extension

contains

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

function intstr(x) result(str)
   integer, intent(in) :: x
   character(:), allocatable :: str
   integer :: l
   l = floor(log10(real(max(abs(x), 1)))) + 1
   if (x < 0) l = l + 1
   allocate (character(l) :: str)
   write (str, '(i0)') x
end function

function realstr(x, n) result(str)
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

function get_extension(filepath) result(extension)
   character(*), intent(in) :: filepath
   character(:), allocatable :: filename, basename, extension
   integer :: pos
   pos = index(filepath, '/', back=.true.)
   filename = filepath(pos+1:)
   if (len(filename) == 0) then
      write (stderr, '(a,1x,a)') 'Missing file name'
      stop
   end if
   pos = index(filename, '.', back=.true.)
   if (pos /= 0) then
      basename = filename(:pos-1)
      extension = filename(pos+1:)
   else
      basename = filename
      extension = ''
   end if
   if (len(basename) == 0) then
      write (stderr, '(a,1x,a)') 'Missing base name of', filename
      stop
   end if
   if (len(extension) == 0) then
      write (stderr, '(a,1x,a)') 'Missing extension of', filename
      stop
   end if
end function

end module
