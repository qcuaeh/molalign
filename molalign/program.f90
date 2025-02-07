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

program molalign

   use kinds
   use flags
   use bounds
   use biasing
   use strutils
   use argparse
   use chemutils
   use readmol
   use writemol
   use translation
   use rotation
   use alignment
   use permutation
   use library

   implicit none

   integer :: i
   integer :: nrec, error
   integer :: natom0, natom1
   integer :: read_unit0, read_unit1, write_unit
   integer, allocatable, dimension(:) :: znums0, znums1
   integer, allocatable, dimension(:) :: types0, types1
   integer, allocatable, dimension(:) :: countlist
   integer, allocatable, dimension(:, :) :: permlist
   character(:), allocatable :: arg
   character(:), allocatable :: pathin1, pathin2, pathout
   character(:), allocatable :: fmtin, fmtin0, fmtin1, fmtout
   character(:), allocatable :: title0, title1
   character(ll) :: posargs(2)
   character(wl), allocatable, dimension(:) :: labels0, labels1
   real(wp) :: rmsd, minrmsd, travec(3), rotmat(3, 3)
   real(wp), allocatable, dimension(:) :: weights0, weights1
   real(wp), allocatable, dimension(:, :) :: coords0, coords1, aligned1
   logical :: remap_flag, mirror_flag, stdin_flag, stdout_flag

   procedure(f_realint), pointer :: weight_function

   ! Set default options

   iter_flag = .false.
   remap_flag = .false.
   trial_flag = .false.
   stdin_flag = .false.
   stdout_flag = .false.
   test_flag = .false.
   stats_flag = .false.
   mirror_flag = .false.
   live_flag = .false.
   bias_flag = .false.

   maxrec = 1
   maxcount = 10
   maxcoord = 32
   bias_tol = 0.35
   bias_scale = 1.e3
   pathout = 'aligned.xyz'

   weight_function => unity

   ! Get user options

   call initarg()

   do while (getarg(arg))

      select case (arg)
      case ('-live')
         live_flag = .true.
      case ('-stats')
         stats_flag = .true.
      case ('-test')
         test_flag = .true.
      case ('-remap', '-sort')
         remap_flag = .true.
      case ('-fast')
         iter_flag = .true.
         bias_flag = .true.
      case ('-mass')
         weight_function => stdmass
      case ('-mirror')
         mirror_flag = .true.
      case ('-count')
         call readoptarg(arg, maxcount)
      case ('-trials')
         trial_flag = .false.
         call readoptarg(arg, maxtrials)
      case ('-tol')
         call readoptarg(arg, bias_tol)
!      case ('-scale')
!         call readoptarg(arg, bias_scale)
      case ('-rec')
         call readoptarg(arg, maxrec)
      case ('-out')
         call readoptarg(arg, pathout)
      case ('-stdin')
         stdin_flag = .true.
         call readoptarg(arg, fmtin)
      case ('-stdout')
         stdout_flag = .true.
         call readoptarg(arg, fmtout)
      case default
         call readposarg(arg, posargs)
      end select

   end do

   if (stdin_flag) then
      fmtin0 = fmtin
      fmtin1 = fmtin
      read_unit0 = input_unit
      read_unit1 = input_unit
   else
      select case (ipos)
      case (0)
         write (error_unit, '(a)') 'Error: Missing arguments'
         stop
      case (1)
         write (error_unit, '(a)') 'Error: Too few arguments'
         stop
      case (2)
         pathin1 = trim(posargs(1))
         pathin2 = trim(posargs(2))
         call open2read(pathin1, read_unit0, fmtin0)
         call open2read(pathin2, read_unit1, fmtin1)
      case default
         write (error_unit, '(a)') 'Error: Too many arguments'
         stop
      end select
   end if

   ! Read coordinates

   call readfile(read_unit0, fmtin0, title0, natom0, labels0, coords0)
   call readfile(read_unit1, fmtin1, title1, natom1, labels1, coords1)

   if (mirror_flag) then
      coords1(1, :) = -coords1(1, :)
   end if

   ! Allocate arrays

   allocate(znums0(natom0), znums1(natom1))
   allocate(types0(natom0), types1(natom1))
   allocate(weights0(natom0), weights1(natom1))
   allocate(aligned1(3, natom1))
   allocate(permlist(natom0, maxrec))
   allocate(countlist(maxrec))

   ! Get atomic numbers, types and weights

   do i = 1, natom0
      call readlabel(labels0(i), znums0(i), types0(i))
      weights0(i) = weight_function(znums0(i))
   end do

   do i = 1, natom1
      call readlabel(labels1(i), znums1(i), types1(i))
      weights1(i) = weight_function(znums1(i))
   end do

   if (stdout_flag) then
      write_unit = output_unit
   else
      fmtout = baseext(pathout)
      if (len(fmtout) > 0) then
         call open2write(pathout, write_unit)
      else
         write (error_unit, '(a)') 'Error: Output file must have an extension'
         stop
      end if
   end if

   ! Sort atoms to minimize MSD

   if (remap_flag) then

      call assign_atoms( &
         natom0, &
         znums0, &
         types0, &
         weights0, &
         coords0, &
         natom1, &
         znums1, &
         types1, &
         weights1, &
         coords1, &
         permlist, &
         countlist, &
         nrec, &
         error)

      if (error /= 0) stop

      call writefile(write_unit, fmtout, 'Reference', natom0, znums0, coords0)

      minrmsd = huge(rmsd)

      do i = 1, nrec

         call align_atoms( &
            natom0, &
            znums0, &
            types0, &
            weights0, &
            coords0, &
            natom1, &
            znums1(permlist(:, i)), &
            types1(permlist(:, i)), &
            weights1(permlist(:, i)), &
            coords1(:, permlist(:, i)), &
            travec, &
            rotmat, &
            error)

         aligned1 = translated(natom1, rotated(natom1, coords1, rotmat), travec)
         rmsd = sqrt(sum(weights0*sum((aligned1(:, permlist(:, i)) - coords0)**2, dim=1))/sum(weights0))
         minrmsd = min(minrmsd, rmsd)

         call writefile(write_unit, fmtout, 'RMSD '//realstr(rmsd, 4), natom1, &
            znums1(permlist(:, i)), aligned1(:, permlist(:, i)))

      end do

      if (.not. stats_flag) write (output_unit, '(a)') realstr(minrmsd, 4)

   else

      ! Align atoms to minimize RMSD

      call align_atoms( &
         natom0, &
         znums0, &
         types0, &
         weights0, &
         coords0, &
         natom1, &
         znums1, &
         types1, &
         weights1, &
         coords1, &
         travec, &
         rotmat, &
         error)

      if (error /= 0) stop

      aligned1 = translated(natom1, rotated(natom1, coords1, rotmat), travec)
      rmsd = sqrt(sum(weights0*sum((aligned1 - coords0)**2, dim=1))/sum(weights0))
      write (output_unit, '(a)') realstr(rmsd, 4)
      call writefile(write_unit, fmtout, 'Reference', natom0, znums0, coords0)
      call writefile(write_unit, fmtout, 'RMSD '//realstr(rmsd, 4), natom1, znums1, aligned1)

   end if

end program
