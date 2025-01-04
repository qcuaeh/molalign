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
use parameters
use globals
use common_types
use sorting
use strutils
use bipartition
use bipartitioning
use molecule

implicit none

abstract interface
   subroutine bias_proc( mol1, mol2, eltypes, mnadiffs)
      use parameters
      use common_types
      use bipartition
      use molecule
      type(mol_type), intent(in) :: mol1, mol2
      type(bipartition_type), intent(in) :: eltypes
      type(intmatrix_type), dimension(:), allocatable, intent(out) :: mnadiffs
   end subroutine
end interface

real(rk) :: bias_scale
procedure(bias_proc), pointer :: bias_procedure

contains

subroutine bias_none( mol1, mol2, eltypes, mnadiffs)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_type), intent(in) :: eltypes
   type(intmatrix_type), dimension(:), allocatable, intent(out) :: mnadiffs
   ! Local variables
   integer :: h, i, j
   integer :: part_size1, part_size2

   allocate (mnadiffs(eltypes%num_parts))
   do h = 1, eltypes%num_parts
      part_size1 = eltypes%parts(h)%part_size1
      part_size2 = eltypes%parts(h)%part_size2
      allocate (mnadiffs(h)%n(part_size1, part_size2))
      do i = 1, part_size1
         do j = 1, part_size2
            mnadiffs(h)%n(j, i) = 0
         end do
      end do
   end do

end subroutine

! Iteratively compute MNA types
subroutine bias_mna( mol1, mol2, eltypes, mnadiffs)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_type), intent(in) :: eltypes
   type(intmatrix_type), dimension(:), allocatable, intent(out) :: mnadiffs
   ! Local variables
   integer :: h, i, j, level
   integer :: iatom, jatom
   integer :: part_size1, part_size2
   type(bipartition_type) :: mnatypes, subtypes

   allocate (mnadiffs(eltypes%num_parts))
   do h = 1, eltypes%num_parts
      part_size1 = eltypes%parts(h)%part_size1
      part_size2 = eltypes%parts(h)%part_size2
      allocate (mnadiffs(h)%n(part_size1, part_size2))
      mnadiffs(h)%n = 0
   end do

   ! eltypes are our initial mnatypes
   mnatypes = eltypes

   level = 0
   do

!      write (stderr, *)
!      write (stderr, '(a)') repeat('-- level '//str(level)//' --', 6)
!      call mnatypes%print_parts()
      level = level + 1

      do h = 1, eltypes%num_parts
         do j = 1, eltypes%parts(h)%part_size2
            jatom = eltypes%parts(h)%items2(j)
            do i = 1, eltypes%parts(h)%part_size1
               iatom = eltypes%parts(h)%items1(i)
               if (mnatypes%idcs1(iatom) /= mnatypes%idcs2(jatom)) then
                  mnadiffs(h)%n(i, j) = mnadiffs(h)%n(i, j) + 1
               end if
            end do
         end do
      end do

      ! Compute next level MNA types
      call levelup_crossmnatypes(mol1, mol2, mnatypes, subtypes)

      ! Exit the loop if types did not change
      if (subtypes == mnatypes) exit

      ! Update mnatypes
      mnatypes = subtypes

   end do

end subroutine

end module
