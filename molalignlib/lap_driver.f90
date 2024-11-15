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
use stdio
use types
use kinds
use bounds
use bipartition
use permutation
use lap_solvers

implicit none

private
public assign_atoms_proc
public assign_atoms_generic
public assign_atoms
public assign_atoms_biased
public assign_atoms_nearest
public assign_atoms_pruned

abstract interface
   subroutine assign_atoms_proc( eltypes, coords1, coords2, prunemask, mnadists, atomperm)
      use kinds
      use types
      use bipartition
      type(bipartition_type), intent(in) :: eltypes
      real(rk), dimension(:, :), intent(in) :: coords1, coords2
      type(boolmatrix_type), dimension(:), intent(in) :: prunemask
      type(realmatrix_type), dimension(:), intent(in) :: mnadists
      integer, dimension(:), intent(out) :: atomperm
   end subroutine
end interface

procedure(assign_atoms_proc), pointer :: assign_atoms_generic

contains

! Find best correspondence between points sets with fixed orientation
subroutine assign_atoms_nearest( eltypes, coords1, coords2, prunemask, mnadists, atomperm)
   type(bipartition_type), intent(in) :: eltypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   type(boolmatrix_type), dimension(:), intent(in) :: prunemask
   type(realmatrix_type), dimension(:), intent(in) :: mnadists
   integer, dimension(:), intent(out) :: atomperm
   ! Local variables
   integer :: h, part_size1
   integer, allocatable :: auxmap(:)
   integer, allocatable, dimension(:) :: atomidcs1, atomidcs2
   real(rk) :: dist

   allocate (auxmap(eltypes%largest_part_size))

   ! Fill distance matrix for each block

   do h = 1, eltypes%num_parts
      part_size1 = eltypes%parts(h)%size1
      atomidcs1 = eltypes%parts(h)%list1
      atomidcs2 = eltypes%parts(h)%list2
      call minperm_nearest(part_size1, atomidcs1, atomidcs2, coords1, coords2, auxmap, dist)
      atomperm(atomidcs1) = atomidcs2(auxmap(:part_size1))
   end do

end subroutine

! Find best correspondence between points sets with fixed orientation
subroutine assign_atoms_pruned( eltypes, coords1, coords2, prunemask, mnadists, atomperm)
   type(bipartition_type), intent(in) :: eltypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   type(boolmatrix_type), dimension(:), intent(in) :: prunemask
   type(realmatrix_type), dimension(:), intent(in) :: mnadists
   integer, dimension(:), intent(out) :: atomperm
   ! Local variables
   integer :: h, part_size1
   integer, allocatable :: auxmap(:)
   integer, allocatable, dimension(:) :: atomidcs1, atomidcs2
   real(rk) :: dist

   allocate (auxmap(eltypes%largest_part_size))

   ! Optimize atomperm for each block
   do h = 1, eltypes%num_parts
      part_size1 = eltypes%parts(h)%size1
      atomidcs1 = eltypes%parts(h)%list1
      atomidcs2 = eltypes%parts(h)%list2
      call minperm_pruned(part_size1, atomidcs1, atomidcs2, coords1, coords2, prunemask(h)%b, auxmap, dist)
      atomperm(atomidcs1) = atomidcs2(auxmap(:part_size1))
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
   integer, allocatable :: auxmap(:)

   allocate (auxmap(eltypes%largest_part_size))

   ! Optimize atomperm for each block
   do h = 1, eltypes%num_parts
      call minperm(eltypes%parts(h), coords1, coords2, auxmap, dist)
      atomperm(eltypes%parts(h)%list1) = eltypes%parts(h)%list2(auxmap(:eltypes%parts(h)%size1))
   end do

end subroutine

! Find best correspondence between points sets with fixed orientation
subroutine assign_atoms_biased( eltypes, coords1, coords2, prunemask, mnadists, atomperm)
   type(bipartition_type), intent(in) :: eltypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   type(boolmatrix_type), dimension(:), intent(in) :: prunemask
   type(realmatrix_type), dimension(:), intent(in) :: mnadists
   integer, dimension(:), intent(out) :: atomperm
   ! Local variables
   integer :: h
   integer, allocatable :: auxmap(:)
   real(rk) :: dist

   allocate (auxmap(eltypes%largest_part_size))

   ! Optimize atomperm for each block
   do h = 1, eltypes%num_parts
      call minperm_biased(eltypes%parts(h), coords1, coords2, mnadists(h)%x, auxmap, dist)
      atomperm(eltypes%parts(h)%list1) = eltypes%parts(h)%list2(auxmap(:eltypes%parts(h)%size1))
   end do

end subroutine

end module
