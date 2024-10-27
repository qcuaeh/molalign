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

module biasing
use kinds
use flags
use bounds
use sorting
use partition
use assorting

implicit none

abstract interface
   subroutine f_bias( eltypes, adjlists1, adjlists2, mnadists)
      use kinds
      use partition
      type(partition_type), intent(in) :: eltypes
      type(atomlist_type), allocatable, dimension(:), intent(in) :: adjlists1, adjlists2
      type(rematrix_type), dimension(:), allocatable, intent(out) :: mnadists
   end subroutine
end interface

real(rk) :: bias_scale
procedure(f_bias), pointer :: bias_procedure

contains

subroutine bias_none( eltypes, adjlists1, adjlists2, mnadists)
   type(partition_type), intent(in) :: eltypes
   type(atomlist_type), allocatable, dimension(:), intent(in) :: adjlists1, adjlists2
   type(rematrix_type), dimension(:), allocatable, intent(out) :: mnadists
   ! Local variables
   integer :: h, i, j
   integer :: iatom, jatom
   integer :: subset1_size, subset2_size

   allocate (mnadists(eltypes%partition_size))
   do h = 1, eltypes%partition_size
      subset1_size = eltypes%parts(h)%subset1%part_size
      subset2_size = eltypes%parts(h)%subset2%part_size
      allocate (mnadists(h)%x(subset1_size, subset2_size))
      do i = 1, subset1_size
         do j = 1, subset2_size
            mnadists(h)%x(j, i) = 0
         end do
      end do
   end do

end subroutine

subroutine bias_mna( eltypes, adjlists1, adjlists2, mnadists)
   type(partition_type), intent(in) :: eltypes
   type(atomlist_type), allocatable, dimension(:), intent(in) :: adjlists1, adjlists2
   type(rematrix_type), dimension(:), allocatable, intent(out) :: mnadists
   ! Local variables
   integer :: h, i, j
   integer :: maxequiv
   integer :: iatom, jatom
   integer :: subset1_size, subset2_size
   integer, allocatable :: equivmat(:, :)

   ! Calculate MNA equivalence matrix
   call compute_equivmat(eltypes, adjlists1, adjlists2, equivmat)

   maxequiv = 0
   do h = 1, eltypes%partition_size
      do i = 1, eltypes%parts(h)%subset1%part_size
         iatom = eltypes%parts(h)%subset1%indices(i)
         do j = 1, eltypes%parts(h)%subset2%part_size
            jatom = eltypes%parts(h)%subset2%indices(j)
            maxequiv = max(maxequiv, equivmat(jatom, iatom))
         end do
      end do
   end do

   allocate (mnadists(eltypes%partition_size))
   do h = 1, eltypes%partition_size
      subset1_size = eltypes%parts(h)%subset1%part_size
      subset2_size = eltypes%parts(h)%subset2%part_size
      allocate (mnadists(h)%x(subset1_size, subset2_size))
      do i = 1, subset1_size
         iatom = eltypes%parts(h)%subset1%indices(i)
         do j = 1, subset2_size
            jatom = eltypes%parts(h)%subset2%indices(j)
            mnadists(h)%x(j, i) = bias_scale**2*(maxequiv - equivmat(jatom, iatom))
         end do
      end do
   end do

end subroutine

! Calculate the maximum common MNA level for all atom cross assignments
subroutine compute_equivmat( eltypes, adjlists1, adjlists2, equivmat)
   type(partition_type), intent(in) :: eltypes
   type(atomlist_type), allocatable, dimension(:), intent(in) :: adjlists1, adjlists2
   integer, allocatable, dimension(:, :), intent(out) :: equivmat
   ! Local variables
   integer :: h, i, j
   integer :: iatom, jatom
   integer :: subset1_size, subset2_size
   integer :: ntype, nintype, level
   integer, allocatable, dimension(:) :: types1, types2
   integer, allocatable, dimension(:) :: intypes1, intypes2

   allocate (equivmat(eltypes%total_num_elems1, eltypes%total_num_elems2))
   allocate (types1(eltypes%total_num_elems1), types2(eltypes%total_num_elems2))
   allocate (intypes1(eltypes%total_num_elems1), intypes2(eltypes%total_num_elems2))

   level = 1
   nintype = eltypes%partition_size
   intypes1 = eltypes%partition_map1
   intypes2 = eltypes%partition_map2

   do

      do h = 1, eltypes%partition_size
         do i = 1, eltypes%parts(h)%subset1%part_size
            iatom = eltypes%parts(h)%subset1%indices(i)
            do j = 1, eltypes%parts(h)%subset2%part_size
               jatom = eltypes%parts(h)%subset2%indices(j)
               if (intypes1(iatom) == intypes2(jatom)) then
                  equivmat(jatom, iatom) = level
               end if
            end do
         end do
      end do

      call compute_crossmnatypes(adjlists1, adjlists2, nintype, intypes1, intypes2, ntype, types1, types2)

      if (all(types1 == intypes1) .and. all(types2 == intypes2)) exit

      nintype = ntype
      intypes1 = types1
      intypes2 = types2
      level = level + 1

   end do

end subroutine

end module
