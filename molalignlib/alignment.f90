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
use parameters
use eigen
use rotation

implicit none

private
public sqdistsum
public leastsqdistsum
public leasteigquat

contains

real(rk) function sqdistsum(atomperm, coords1, coords2)
   integer, dimension(:), intent(in) :: atomperm
   real(rk), dimension(:,:), intent(in) :: coords1, coords2

   sqdistsum = sum(sum((coords1 - coords2(:, atomperm))**2, dim=1))

end function

real(rk) function leastsqdistsum(atomperm, coords1, coords2)
! Purpose: Calculate least square distance from eigenvalues
   integer, dimension(:), intent(in) :: atomperm
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   real(rk) :: residuals(4, 4)

   call kearsley(atomperm, coords1, coords2, residuals)
   ! eigenvalue can be negative due to numerical errors
   leastsqdistsum = max(leasteigval(residuals), 0._rk)

end function

function leasteigquat(atomperm, coords1, coords2)
! Purpose: Calculate rotation quaternion which minimzes the square distance
   integer, dimension(:), intent(in) :: atomperm
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   real(rk) :: leasteigquat(4), residuals(4, 4)

   call kearsley(atomperm, coords1, coords2, residuals)
   leasteigquat = leasteigvec(residuals)

end function

subroutine kearsley(atomperm, coords1, coords2, residuals)
! Purpose: Find the best orientation by least squares minimization
! Reference: Acta Cryst. (1989). A45, 208-210
   integer, dimension(:), intent(in) :: atomperm
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   ! Local variables
   integer :: i, num_atoms
   real(rk) :: residuals(4, 4)
   real(rk), dimension(:,:), allocatable :: p, q

   num_atoms = size(atomperm)

   allocate (p(3, num_atoms))
   allocate (q(3, num_atoms))

   do i = 1, num_atoms
      p(:, i) = coords1(:, i) + coords2(:, atomperm(i))
      q(:, i) = coords1(:, i) - coords2(:, atomperm(i))
   end do

   ! Calculate upper matrix elements

   residuals = 0

   do i = 1, num_atoms
      residuals(1, 1) = residuals(1, 1) + (q(1, i)**2 + q(2, i)**2 + q(3, i)**2)
      residuals(1, 2) = residuals(1, 2) + (p(2, i)*q(3, i) - q(2, i)*p(3, i))
      residuals(1, 3) = residuals(1, 3) + (q(1, i)*p(3, i) - p(1, i)*q(3, i))
      residuals(1, 4) = residuals(1, 4) + (p(1, i)*q(2, i) - q(1, i)*p(2, i))
      residuals(2, 2) = residuals(2, 2) + (p(2, i)**2 + p(3, i)**2 + q(1, i)**2)
      residuals(2, 3) = residuals(2, 3) + (q(1, i)*q(2, i) - p(1, i)*p(2, i))
      residuals(2, 4) = residuals(2, 4) + (q(1, i)*q(3, i) - p(1, i)*p(3, i))
      residuals(3, 3) = residuals(3, 3) + (p(1, i)**2 + p(3, i)**2 + q(2, i)**2)
      residuals(3, 4) = residuals(3, 4) + (q(2, i)*q(3, i) - p(2, i)*p(3, i))
      residuals(4, 4) = residuals(4, 4) + (p(1, i)**2 + p(2, i)**2 + q(3, i)**2)
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
