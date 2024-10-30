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

module alignment
use kinds
use jacobi
use rotation

implicit none

private
public squaredist
public leastsquaredist
public leastrotquat

contains

real(wp) function squaredist(natom, weights, coords1, coords2, atomperm)
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: atomperm
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords1, coords2

   squaredist = sum(weights*sum((coords1 - coords2(:, atomperm))**2, dim=1))

end function

real(wp) function leastsquaredist(natom, weights, coords1, coords2, atomperm)
! Purpose: Calculate least square distance from eigenvalues
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: atomperm
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords1, coords2
   real(wp) :: residuals(4, 4)

   call kearsley(natom, weights, coords1, coords2, atomperm, residuals)
   ! eigenvalue can be negative due to numerical errors
   leastsquaredist = max(leasteigval(residuals), 0._wp)

end function

!real(wp) function leastsquaredist(natom, weights, coords1, coords2, atomperm)
!! Purpose: Calculate least square distance from aligned coordinates
!   integer, intent(in) :: natom
!   integer, dimension(:), intent(in) :: atomperm
!   real(wp), dimension(:), intent(in) :: weights
!   real(wp), dimension(:, :), intent(in) :: coords1, coords2
!   real(wp) :: residuals(4, 4), eigval(4)
!
!   call kearsley(natom, weights, coords1, coords2, atomperm, residuals)
!   call syevec4(residuals, eigval)
!   leastsquaredist = squaredist(natom, weights, coords1, rotated(natom, coords2, residuals(:, 1)), atomperm)
!
!end function

function leastrotquat(natom, weights, coords1, coords2, atomperm)
! Purpose: Calculate rotation quaternion which minimzes the square distance
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: atomperm
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords1, coords2
   real(wp) :: leastrotquat(4), residuals(4, 4)

   call kearsley(natom, weights, coords1, coords2, atomperm, residuals)
   leastrotquat = leasteigvec(residuals)

end function

subroutine kearsley(natom, weights, coords1, coords2, atomperm, residuals)
! Purpose: Find the best orientation by least squares minimization
! Reference: Acta Cryst. (1989). A45, 208-210
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: atomperm
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords1, coords2

   integer :: i
   real(wp) :: residuals(4, 4), p(3, natom), q(3, natom)

   residuals = 0

   do i = 1, natom
      p(:, i) = coords1(:, i) + coords2(:, atomperm(i))
      q(:, i) = coords1(:, i) - coords2(:, atomperm(i))
   end do

   ! Calculate upper matrix elements

   do i = 1, natom
      residuals(1, 1) = residuals(1, 1) + weights(i)*(q(1, i)**2 + q(2, i)**2 + q(3, i)**2)
      residuals(1, 2) = residuals(1, 2) + weights(i)*(p(2, i)*q(3, i) - q(2, i)*p(3, i))
      residuals(1, 3) = residuals(1, 3) + weights(i)*(q(1, i)*p(3, i) - p(1, i)*q(3, i))
      residuals(1, 4) = residuals(1, 4) + weights(i)*(p(1, i)*q(2, i) - q(1, i)*p(2, i))
      residuals(2, 2) = residuals(2, 2) + weights(i)*(p(2, i)**2 + p(3, i)**2 + q(1, i)**2)
      residuals(2, 3) = residuals(2, 3) + weights(i)*(q(1, i)*q(2, i) - p(1, i)*p(2, i))
      residuals(2, 4) = residuals(2, 4) + weights(i)*(q(1, i)*q(3, i) - p(1, i)*p(3, i))
      residuals(3, 3) = residuals(3, 3) + weights(i)*(p(1, i)**2 + p(3, i)**2 + q(2, i)**2)
      residuals(3, 4) = residuals(3, 4) + weights(i)*(q(2, i)*q(3, i) - p(2, i)*p(3, i))
      residuals(4, 4) = residuals(4, 4) + weights(i)*(p(1, i)**2 + p(2, i)**2 + q(3, i)**2)
   end do

   ! Symmetrize matrix

   residuals(2, 1) = residuals(1, 2)
   residuals(3, 1) = residuals(1, 3)
   residuals(4, 1) = residuals(1, 4)
   residuals(3, 2) = residuals(2, 3)
   residuals(4, 2) = residuals(2, 4)
   residuals(4, 3) = residuals(3, 4)

end subroutine

end module
