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

module lap_driver
use parameters
use globals
use common_types
use bipartition
use permutation
use lap_solvers

implicit none

private
public assign_atoms
public assign_atoms_biased
public assign_atoms_nearest
public assign_atoms_pruned

contains

! Find best correspondence between points sets with fixed orientation
subroutine assign_atoms_nearest( eltypes, coords1, coords2, pruned, mnadists, atomperm)
   type(bipartition_type), intent(in) :: eltypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   type(boolmatrix_type), dimension(:), intent(in) :: pruned
   type(realmatrix_type), dimension(:), intent(in) :: mnadists
   integer, dimension(:), intent(out) :: atomperm
   ! Local variables
   integer :: h, part_size1
   integer, allocatable :: auxperm(:)
   integer, allocatable, dimension(:) :: atomidcs1, atomidcs2
   real(rk) :: dist

   allocate (auxperm(eltypes%largest_part_size))

   ! Fill distance matrix for each block

   do h = 1, eltypes%num_parts
      part_size1 = eltypes%parts(h)%size1
      atomidcs1 = eltypes%parts(h)%items1
      atomidcs2 = eltypes%parts(h)%items2
      call minperm_nearest(part_size1, atomidcs1, atomidcs2, coords1, coords2, auxperm, dist)
      atomperm(atomidcs1) = atomidcs2(auxperm(:part_size1))
   end do

end subroutine

! Find best correspondence between points sets with fixed orientation
subroutine assign_atoms_pruned( eltypes, coords1, coords2, pruned, atomperm)
   type(bipartition_type), intent(in) :: eltypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   type(boolmatrix_type), dimension(:), intent(in) :: pruned
   integer, dimension(:), intent(out) :: atomperm
   ! Local variables
   integer :: h, part_size1
   integer, allocatable :: auxperm(:)
   integer, allocatable, dimension(:) :: atomidcs1, atomidcs2
   real(rk) :: dist

   allocate (auxperm(eltypes%largest_part_size))

   ! Optimize atomperm for each block
   do h = 1, eltypes%num_parts
      part_size1 = eltypes%parts(h)%size1
      atomidcs1 = eltypes%parts(h)%items1
      atomidcs2 = eltypes%parts(h)%items2
      call minperm_pruned(part_size1, atomidcs1, atomidcs2, coords1, coords2, pruned(h)%b, auxperm, dist)
      atomperm(atomidcs1) = atomidcs2(auxperm(:part_size1))
   end do

end subroutine

! Find best correspondence between points sets with fixed orientation
subroutine assign_atoms( eltypes, coords1, coords2, atomperm, dist)
   type(bipartition_type), intent(in) :: eltypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   integer, dimension(:), intent(out) :: atomperm
   real(rk), intent(out) :: dist
   ! Local variables
   integer :: h
   integer, allocatable :: auxperm(:)

   allocate (auxperm(eltypes%largest_part_size))

   ! Optimize atomperm for each block
   do h = 1, eltypes%num_parts
      call minperm(eltypes%parts(h), coords1, coords2, auxperm, dist)
      atomperm(eltypes%parts(h)%items1) = eltypes%parts(h)%items2(auxperm(:eltypes%parts(h)%size1))
   end do

end subroutine

! Find best correspondence between points sets with fixed orientation
subroutine assign_atoms_biased( eltypes, coords1, coords2, pruned, mnadists, atomperm)
   type(bipartition_type), intent(in) :: eltypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   type(boolmatrix_type), dimension(:), intent(in) :: pruned
   type(realmatrix_type), dimension(:), intent(in) :: mnadists
   integer, dimension(:), intent(out) :: atomperm
   ! Local variables
   integer :: h
   integer, allocatable :: auxperm(:)
   real(rk) :: dist

   allocate (auxperm(eltypes%largest_part_size))

   ! Optimize atomperm for each block
   do h = 1, eltypes%num_parts
      call minperm_biased(eltypes%parts(h), coords1, coords2, mnadists(h)%x, auxperm, dist)
      atomperm(eltypes%parts(h)%items1) = eltypes%parts(h)%items2(auxperm(:eltypes%parts(h)%size1))
   end do

end subroutine

end module
