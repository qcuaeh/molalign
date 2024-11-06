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
use stdio
use kinds
use types
use flags
use bounds
use sorting
use strutils
use assorting
use partition
use bipartition
use molecule

implicit none

abstract interface
   subroutine bias_proc( mol1, mol2, eltypes, mnadists)
      use kinds
      use types
      use partition
      use bipartition
      use molecule
      type(mol_type), intent(in) :: mol1, mol2
      type(bipartition_type), intent(in) :: eltypes
      type(realmatrix_type), dimension(:), allocatable, intent(out) :: mnadists
   end subroutine
end interface

real(rk) :: bias_scale
procedure(bias_proc), pointer :: bias_procedure

contains

subroutine bias_none( mol1, mol2, eltypes, mnadists)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_type), intent(in) :: eltypes
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

subroutine bias_mna( mol1, mol2, eltypes, mnadists)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_type), intent(in) :: eltypes
   type(realmatrix_type), dimension(:), allocatable, intent(out) :: mnadists
   ! Local variables
   integer :: h, i, j
   integer :: subset1_size, subset2_size
   type(intmatrix_type), dimension(:), allocatable :: mnadiffs

   ! Calculate MNA equivalence matrix
   call compute_equivmat(mol1, mol2, eltypes, mnadiffs)

   allocate (mnadists(eltypes%partition_size))
   do h = 1, eltypes%partition_size
      subset1_size = eltypes%parts(h)%subset1%part_size
      subset2_size = eltypes%parts(h)%subset2%part_size
      allocate (mnadists(h)%x(subset1_size, subset2_size))
      do i = 1, subset1_size
         do j = 1, subset2_size
            mnadists(h)%x(j, i) = mnadiffs(h)%n(j, i)*bias_scale**2
         end do
      end do
   end do

end subroutine

! Iteratively compute MNA types
subroutine compute_equivmat(mol1, mol2, eltypes, mnadiffs)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_type), intent(in) :: eltypes
   type(intmatrix_type), dimension(:), allocatable, intent(out) :: mnadiffs
   ! Local variables
   integer :: h, i, j, level
   integer :: iatom, jatom
   integer :: subset1_size, subset2_size
   type(bipartition_type) :: mnatypes, submnatypes

   allocate (mnadiffs(eltypes%partition_size))
   do h = 1, eltypes%partition_size
      subset1_size = eltypes%parts(h)%subset1%part_size
      subset2_size = eltypes%parts(h)%subset2%part_size
      allocate (mnadiffs(h)%n(subset1_size, subset2_size))
      mnadiffs(h)%n = 0
   end do

   ! eltypes are our initial mnatypes
   mnatypes = eltypes

   level = 0
   do

!      write (stderr, *)
!      write (stderr, '(a)') repeat('-- level '//intstr(level)//' --', 6)
!      call mnatypes%print_parts()
      level = level + 1

      do h = 1, eltypes%partition_size
         do i = 1, eltypes%parts(h)%subset1%part_size
            iatom = eltypes%parts(h)%subset1%indices(i)
            do j = 1, eltypes%parts(h)%subset2%part_size
               jatom = eltypes%parts(h)%subset2%indices(j)
               if (mnatypes%partition_map1(iatom) /= mnatypes%partition_map2(jatom)) then
                  mnadiffs(h)%n(j, i) = mnadiffs(h)%n(j, i) + 1
               end if
            end do
         end do
      end do

      ! Compute next level MNA types
      call levelup_crossmnatypes(mol1, mol2, mnatypes, submnatypes)

      ! Exit the loop if types did not change
      if (submnatypes == mnatypes) exit

      ! Update mnatypes
      mnatypes = submnatypes

   end do

!block
!   use lap_solvers
!   integer, allocatable :: perm(:)
!   real :: dist
!   do h = 1, eltypes%partition_size
!      write (stderr, *)
!      write (stderr, '(a)') repeat('-- block '//intstr(h)//' --', 6)
!      do i = 1, eltypes%parts(h)%subset1%part_size
!         write (stderr, '(*(1x,i2))') mnadiffs(h)%n(:, i)
!      end do
!      allocate(perm(eltypes%parts(h)%subset1%part_size))
!      call minperm(eltypes%parts(h)%subset1%part_size, mnadiffs(h)%n, perm, dist)
!      write (stderr, *)
!      do i = 1, eltypes%parts(h)%subset1%part_size
!         write (stderr, '(*(i3,1x,i3,3x,i3,1x,i3,3x,i2))') &
!               i, &
!               perm(i), &
!               eltypes%parts(h)%subset1%indices(i), &
!               eltypes%parts(h)%subset2%indices(perm(i)), &
!               mnadiffs(h)%n(perm(i), i)
!      end do
!      deallocate(perm)
!   end do
!end block

end subroutine

end module
