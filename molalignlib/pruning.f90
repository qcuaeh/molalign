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
use parameters
use common_types
use globals
use sorting
use bipartition

implicit none

real(rk) :: prune_tol
procedure(prune_proc), pointer :: prune_procedure

abstract interface
   subroutine prune_proc( eltypes, coords1, coords2, pruned)
      use parameters
      use common_types
      use bipartition
      type(bipartition_type), intent(in) :: eltypes
      real(rk), dimension(:, :), intent(in) :: coords1, coords2
      type(boolmatrix_type), dimension(:), allocatable, intent(out) :: pruned
   end subroutine
end interface

contains

subroutine prune_none( eltypes, coords1, coords2, pruned)
   type(bipartition_type), intent(in) :: eltypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   type(boolmatrix_type), dimension(:), allocatable, intent(out) :: pruned
   ! Local variables
   integer :: h, i, j
   integer :: part_size1, part_size2

   allocate (pruned(eltypes%num_parts))
   do h = 1, eltypes%num_parts
      part_size1 = eltypes%parts(h)%part_size1
      part_size2 = eltypes%parts(h)%part_size2
      allocate (pruned(h)%b(part_size1, part_size2))
      do i = 1, part_size1
         do j = 1, part_size2
            pruned(h)%b(j, i) = .true.
         end do
      end do
   end do

end subroutine

subroutine prune_rd( eltypes, coords1, coords2, pruned)
   type(bipartition_type), intent(in) :: eltypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   type(boolmatrix_type), dimension(:), allocatable, intent(out) :: pruned
   ! Local variables
   integer :: h, i, j, k
   integer :: iatom, jatom
   integer :: part_size1, part_size2
   type(nested_reallist_type), allocatable, dimension(:) :: dists1, dists2

   allocate (dists1(size(coords1, dim=2)))
   allocate (dists2(size(coords2, dim=2)))
   do i = 1, size(coords1, dim=2)
      allocate (dists1(i)%s(eltypes%num_parts))
      allocate (dists2(i)%s(eltypes%num_parts))
      do h = 1, eltypes%num_parts
         allocate (dists1(i)%s(h)%x(eltypes%parts(h)%part_size1))
         allocate (dists2(i)%s(h)%x(eltypes%parts(h)%part_size2))
      end do
   end do

   do i = 1, size(coords1, dim=2)
      do h = 1, eltypes%num_parts
         do j = 1, eltypes%parts(h)%part_size1
            jatom = eltypes%parts(h)%items1(j)
            dists1(i)%s(h)%x(j) = sqrt(sum((coords1(:, jatom) - coords1(:, i))**2))
         end do
         call sort(dists1(i)%s(h)%x)
      end do
   end do

   do i = 1, size(coords2, dim=2)
      do h = 1, eltypes%num_parts
         do j = 1, eltypes%parts(h)%part_size2
            jatom = eltypes%parts(h)%items2(j)
            dists2(i)%s(h)%x(j) = sqrt(sum((coords2(:, jatom) - coords2(:, i))**2))
         end do
         call sort(dists2(i)%s(h)%x)
      end do
   end do

   allocate (pruned(eltypes%num_parts))
   do h = 1, eltypes%num_parts
      part_size1 = eltypes%parts(h)%part_size1
      part_size2 = eltypes%parts(h)%part_size2
      allocate (pruned(h)%b(part_size1, part_size2))
      pruned(h)%b = .false.
      do i = 1, part_size1
         iatom = eltypes%parts(h)%items1(i)
         do j = 1, part_size2
            jatom = eltypes%parts(h)%items2(j)
            do k = 1, eltypes%num_parts
               if (any(abs(dists2(jatom)%s(k)%x - dists1(iatom)%s(k)%x) > prune_tol)) then
                  pruned(h)%b(j, i) = .true.
                  exit
               end if
            end do
         end do
      end do
   end do

end subroutine

end module
