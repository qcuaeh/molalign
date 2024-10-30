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

module permutation
use kinds

implicit none

contains

! Get an identity permutation
function identityperm(n)
   integer, intent(in) :: n
   integer :: identityperm(n)
   integer :: i

   do i = 1, n
      identityperm(i) = i
   end do

end function

! Get the inverse permutation
function inverseperm(perm)
   integer, intent(in) :: perm(:)
   integer :: inverseperm(size(perm))
   integer :: i

   do i = 1, size(perm)
      inverseperm(perm(i)) = i
   end do

end function

! Check if array is a valid permutation
logical function isperm(arr)
   implicit none
   integer, intent(in) :: arr(:)
   integer :: i, N
   logical :: seen(size(arr))
   
   N = size(arr)
   seen = .false.
   isperm = .true.
   
   do i = 1, N
      if (arr(i) < 1 .or. arr(i) > N .or. seen(arr(i))) then
         isperm = .false.
         return
      end if
      seen(arr(i)) = .true.
   end do

end function

end module
