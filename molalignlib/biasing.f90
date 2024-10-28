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
use assorting
use partition
use molecule

implicit none

abstract interface
   subroutine f_bias( atoms1, atoms2, eltypes, mnadists)
      use kinds
      use partition
      use molecule
      type(atom_type), dimension(:), intent(in) :: atoms1, atoms2
      type(partition_type), intent(in) :: eltypes
      type(realmatrix_type), dimension(:), allocatable, intent(out) :: mnadists
   end subroutine
end interface

real(rk) :: bias_scale
procedure(f_bias), pointer :: bias_procedure

contains

subroutine bias_none( atoms1, atoms2, eltypes, mnadists)
   type(atom_type), dimension(:), intent(in) :: atoms1, atoms2
   type(partition_type), intent(in) :: eltypes
   type(realmatrix_type), dimension(:), allocatable, intent(out) :: mnadists
   ! Local variables
   integer :: h, i, j
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

subroutine bias_mna( atoms1, atoms2, eltypes, mnadists)
   type(atom_type), dimension(:), intent(in) :: atoms1, atoms2
   type(partition_type), intent(in) :: eltypes
   type(realmatrix_type), dimension(:), allocatable, intent(out) :: mnadists
   ! Local variables
   integer :: h, i, j, last_level
   integer :: subset1_size, subset2_size
   integer, allocatable :: equivmat(:, :)
   type(intmatrix_type), dimension(:), allocatable :: last_common_levels

   ! Calculate MNA equivalence matrix
   call compute_equivmat(atoms1, atoms2, eltypes, last_level, last_common_levels)

   allocate (mnadists(eltypes%partition_size))
   do h = 1, eltypes%partition_size
      subset1_size = eltypes%parts(h)%subset1%part_size
      subset2_size = eltypes%parts(h)%subset2%part_size
      allocate (mnadists(h)%x(subset1_size, subset2_size))
      do i = 1, subset1_size
         do j = 1, subset2_size
            mnadists(h)%x(j, i) = bias_scale**2*(last_level - last_common_levels(h)%n(j, i))
         end do
      end do
   end do

end subroutine

! Iteratively compute MNA types
subroutine compute_equivmat(atoms1, atoms2, eltypes, last_level, last_common_levels)
   type(atom_type), dimension(:), intent(in) :: atoms1, atoms2
   type(partition_type), intent(in) :: eltypes
   integer, intent(out) :: last_level
   type(intmatrix_type), dimension(:), allocatable, intent(out) :: last_common_levels
   ! Local variables
   integer :: h, i, j, iatom, jatom
   integer :: subset1_size, subset2_size
   type(partition_type) :: mnatypes, subtypes

   allocate (last_common_levels(eltypes%partition_size))
   do h = 1, eltypes%partition_size
      subset1_size = eltypes%parts(h)%subset1%part_size
      subset2_size = eltypes%parts(h)%subset2%part_size
      allocate (last_common_levels(h)%n(subset1_size, subset2_size))
   end do

   ! eltypes are our initial mnatypes
   mnatypes = eltypes
   last_level = 1

   do

      do h = 1, eltypes%partition_size
         do i = 1, eltypes%parts(h)%subset1%part_size
            iatom = eltypes%parts(h)%subset1%indices(i)
            do j = 1, eltypes%parts(h)%subset2%part_size
               jatom = eltypes%parts(h)%subset2%indices(j)
               if (mnatypes%partition_map1(iatom) == mnatypes%partition_map2(jatom)) then
                  last_common_levels(h)%n(j, i) = last_level
               end if
            end do
         end do
      end do

      ! Compute next last_level MNA subtypes
      call levelup_mnatypes(atoms1, atoms2, mnatypes, subtypes)

      ! Exit the loop if subtypes are unchanged
      if (subtypes == mnatypes) exit

      ! Update mnatypes
      mnatypes = subtypes
      last_level = last_level + 1

   end do

!   call mnatypes%print_parts()

end subroutine

end module
