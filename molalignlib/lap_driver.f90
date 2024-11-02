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
use partition
use permutation
use lap_solvers

implicit none

private
public assign_atoms
public assign_atoms_nearest
public assign_atoms_pruned
public assign_atoms_biased
public assign_proc

abstract interface
   subroutine assign_proc( eltypes, coords1, coords2, prunemask, mnadists, atomperm)
      use kinds
      use types
      use partition
      type(bipartition_type), intent(in) :: eltypes
      real(rk), dimension(:, :), intent(in) :: coords1, coords2
      type(boolmatrix_type), dimension(:), intent(in) :: prunemask
      type(realmatrix_type), dimension(:), intent(in) :: mnadists
      integer, dimension(:), intent(out) :: atomperm
   end subroutine
end interface

procedure(assign_proc), pointer :: assign_atoms

contains

! Find best correspondence between points sets with fixed orientation
subroutine assign_atoms_nearest( eltypes, coords1, coords2, prunemask, mnadists, atomperm)
   type(bipartition_type), intent(in) :: eltypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   type(boolmatrix_type), dimension(:), intent(in) :: prunemask
   type(realmatrix_type), dimension(:), intent(in) :: mnadists
   integer, dimension(:), intent(out) :: atomperm
   ! Local variables
   integer :: h, subset1_size
   integer, allocatable :: auxmap(:)
   integer, allocatable, dimension(:) :: atomidcs1, atomidcs2
   real(rk) :: dist

   allocate (auxmap(eltypes%largest_part_size))

   ! Fill distance matrix for each block

   do h = 1, eltypes%partition_size
      subset1_size = eltypes%parts(h)%subset1%part_size
      atomidcs1 = eltypes%parts(h)%subset1%indices
      atomidcs2 = eltypes%parts(h)%subset2%indices
      call minperm_nearest(subset1_size, atomidcs1, atomidcs2, coords1, coords2, auxmap, dist)
      atomperm(atomidcs1) = atomidcs2(auxmap(:subset1_size))
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
   integer :: h, subset1_size
   integer, allocatable :: auxmap(:)
   integer, allocatable, dimension(:) :: atomidcs1, atomidcs2
   real(rk) :: dist

   allocate (auxmap(eltypes%largest_part_size))

   ! Optimize atomperm for each block
   do h = 1, eltypes%partition_size
      subset1_size = eltypes%parts(h)%subset1%part_size
      atomidcs1 = eltypes%parts(h)%subset1%indices
      atomidcs2 = eltypes%parts(h)%subset2%indices
      call minperm_pruned(subset1_size, atomidcs1, atomidcs2, coords1, coords2, prunemask(h)%b, auxmap, dist)
      atomperm(atomidcs1) = atomidcs2(auxmap(:subset1_size))
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
   integer :: h, subset1_size
   integer, allocatable :: auxmap(:)
   integer, allocatable, dimension(:) :: atomidcs1, atomidcs2
   real(rk) :: dist

   allocate (auxmap(eltypes%largest_part_size))

   ! Optimize atomperm for each block
   do h = 1, eltypes%partition_size
      subset1_size = eltypes%parts(h)%subset1%part_size
      atomidcs1 = eltypes%parts(h)%subset1%indices
      atomidcs2 = eltypes%parts(h)%subset2%indices
      call minperm_biased(subset1_size, atomidcs1, atomidcs2, coords1, coords2, mnadists(h)%x, auxmap, dist)
      atomperm(atomidcs1) = atomidcs2(auxmap(:subset1_size))
   end do

end subroutine

end module
