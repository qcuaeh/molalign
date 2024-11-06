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

module pruning
use kinds
use types
use flags
use bounds
use sorting
use partition
use bipartition

implicit none

real(rk) :: prune_tol
procedure(prune_proc), pointer :: prune_procedure

abstract interface
   subroutine prune_proc( eltypes, coords1, coords2, prunemask)
      use kinds
      use types
      use partition
      use bipartition
      type(bipartition_type), intent(in) :: eltypes
      real(rk), dimension(:, :), intent(in) :: coords1, coords2
      type(boolmatrix_type), dimension(:), allocatable, intent(out) :: prunemask
   end subroutine
end interface

contains

subroutine prune_none( eltypes, coords1, coords2, prunemask)
   type(bipartition_type), intent(in) :: eltypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   type(boolmatrix_type), dimension(:), allocatable, intent(out) :: prunemask
   ! Local variables
   integer :: h, i, j
   integer :: subset1_size, subset2_size

   allocate (prunemask(eltypes%partition_size))
   do h = 1, eltypes%partition_size
      subset1_size = eltypes%parts(h)%subset1%part_size
      subset2_size = eltypes%parts(h)%subset2%part_size
      allocate (prunemask(h)%b(subset1_size, subset2_size))
      do i = 1, subset1_size
         do j = 1, subset2_size
            prunemask(h)%b(j, i) = .true.
         end do
      end do
   end do

end subroutine

subroutine prune_rd( eltypes, coords1, coords2, prunemask)
   type(bipartition_type), intent(in) :: eltypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   type(boolmatrix_type), dimension(:), allocatable, intent(out) :: prunemask
   ! Local variables
   integer :: h, i, j, k
   integer :: iatom, jatom
   integer :: subset1_size, subset2_size
   type(nested_reallist_type), allocatable, dimension(:) :: dists1, dists2
   logical :: pruned

   allocate (dists1(size(coords1, dim=2)))
   allocate (dists2(size(coords2, dim=2)))
   do i = 1, size(coords1, dim=2)
      allocate (dists1(i)%s(eltypes%partition_size))
      allocate (dists2(i)%s(eltypes%partition_size))
      do h = 1, eltypes%partition_size
         allocate (dists1(i)%s(h)%x(eltypes%parts(h)%subset1%part_size))
         allocate (dists2(i)%s(h)%x(eltypes%parts(h)%subset2%part_size))
      end do
   end do

   do i = 1, size(coords1, dim=2)
      do h = 1, eltypes%partition_size
         do j = 1, eltypes%parts(h)%subset1%part_size
            jatom = eltypes%parts(h)%subset1%indices(j)
            dists1(i)%s(h)%x(j) = sqrt(sum((coords1(:, jatom) - coords1(:, i))**2))
         end do
         call sort(dists1(i)%s(h)%x)
      end do
   end do

   do i = 1, size(coords2, dim=2)
      do h = 1, eltypes%partition_size
         do j = 1, eltypes%parts(h)%subset2%part_size
            jatom = eltypes%parts(h)%subset2%indices(j)
            dists2(i)%s(h)%x(j) = sqrt(sum((coords2(:, jatom) - coords2(:, i))**2))
         end do
         call sort(dists2(i)%s(h)%x)
      end do
   end do

   allocate (prunemask(eltypes%partition_size))
   do h = 1, eltypes%partition_size
      subset1_size = eltypes%parts(h)%subset1%part_size
      subset2_size = eltypes%parts(h)%subset2%part_size
      allocate (prunemask(h)%b(subset1_size, subset2_size))
      do i = 1, subset1_size
         iatom = eltypes%parts(h)%subset1%indices(i)
         do j = 1, subset2_size
            jatom = eltypes%parts(h)%subset2%indices(j)
            pruned = .false.
            do k = 1, eltypes%partition_size
               if (any(abs(dists2(jatom)%s(k)%x - dists1(iatom)%s(k)%x) > prune_tol)) then
                  pruned = .true.
                  exit
               end if
            end do
            prunemask(h)%b(j, i) = .not. pruned
         end do
      end do
   end do

end subroutine

end module
