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
use molecule
use printing
use rotation
use translation
use strutils
use chemutils
use alignment
use adjacency
use permutation
use fileio
use argparse
use molalignlib
use biasing
use pruning

implicit none

integer :: i, nrec
integer :: read_unit1, read_unit2, write_unit
integer, allocatable :: countlist(:)
integer, allocatable :: maplist(:, :)
integer, allocatable :: atomperm(:)
character(:), allocatable :: arg, optarg
character(:), allocatable :: fmtin1, fmtin2, fmtout
character(:), allocatable :: optfmtin, optfmtout
character(:), allocatable, target :: pathin1, pathin2
character(:), allocatable :: pathout
logical :: fmtin_flag, fmtout_flag
logical :: remap_flag, pipe_flag, nrec_flag
real(rk) :: travec1(3), travec2(3), rotquat(4)
integer :: adjd, minadjd
real(rk) :: rmsd, minrmsd
type(p_char) :: posargs(2)
type(molecule_type) :: mol1, mol2, auxmol1
type(partition_type) :: eltypes

! Set default options

iter_flag = .false.
bond_flag = .false.
back_flag = .false.
test_flag = .false.
reac_flag = .false.
print_flag = .true.
stats_flag = .false.
trial_flag = .false.
mirror_flag = .false.
remap_flag = .false.
pipe_flag = .false.
fmtin_flag = .false.
fmtout_flag = .false.
nrec_flag = .false.

maxrec = 1
maxcount = 10
maxcoord = 16
maxlevel = 16

prune_tol = 0.5
bias_scale = 1.e3

weights => ones
print_stats => print_stats_dist
assign_atoms => assign_atoms_pruned
bias_procedure => bias_none
prune_procedure => prune_none

posargs(1)%var => pathin1
posargs(2)%var => pathin2
pathout = 'aligned.xyz'

! Get user options

call init_args()

do while (get_arg(arg))

   select case (arg)
   case ('-stats')
      stats_flag = .true.
   case ('-test')
      test_flag = .true.
   case ('-remap')
      remap_flag = .true.
   case ('-near')
      iter_flag = .false.
      bias_procedure => bias_none
      prune_procedure => prune_none
      assign_atoms => assign_atoms_nearest
   case ('-bias')
      iter_flag = .true.
      assign_atoms => assign_atoms_biased
      call read_optarg(arg, optarg)
      select case (optarg)
      case ('mna')
         bias_procedure => bias_mna
      case default
         write (stderr, '(a,1x,a)') 'Error: Unknown -bias option:', optarg
         stop
      end select
   case ('-prune')
      iter_flag = .true.
      assign_atoms => assign_atoms_pruned
      call read_optarg(arg, optarg)
      select case (optarg)
      case ('rd')
         prune_procedure => prune_rd
      case default
         write (stderr, '(a,1x,a)') 'Error: Unknown -prune option:', optarg
         stop
      end select
   case ('-bond')
      bond_flag = .true.
      print_stats => print_stats_diff
   case ('-back')
      back_flag = .true.
   case ('-reac')
      reac_flag = .true.
   case ('-mass')
      weights => stdmasses
   case ('-mirror')
      mirror_flag = .true.
   case ('-count')
      call read_optarg(arg, maxcount)
   case ('-trials')
      trial_flag = .true.
      call read_optarg(arg, maxtrials)
   case ('-tol')
      call read_optarg(arg, prune_tol)
!      case ('-scale')
!         call read_optarg(arg, bias_scale)
   case ('-N')
      nrec_flag = .true.
      call read_optarg(arg, maxrec)
   case ('-out')
      call read_optarg(arg, pathout)
   case ('-fmtin')
      fmtin_flag = .true.
      call read_optarg(arg, optfmtin)
   case ('-fmtout')
      fmtout_flag = .true.
      call read_optarg(arg, optfmtout)
   case ('-pipe')
      pipe_flag = .true.
   case default
      call read_posarg(arg, posargs)
   end select

end do

if (pipe_flag) then
   read_unit1 = stdin
   read_unit2 = stdin
   fmtin1 = 'xyz'
   fmtin2 = 'xyz'
else
   select case (ipos)
   case (0)
      write (stderr, '(a)') 'Error: Missing arguments'
      stop
   case (1)
      write (stderr, '(a)') 'Error: Too few arguments'
      stop
   case (2)
      call open2read(pathin1, read_unit1, fmtin1)
      call open2read(pathin2, read_unit2, fmtin2)
   case default
      write (stderr, '(a)') 'Error: Too many arguments'
      stop
   end select
end if

if (fmtin_flag) then
   fmtin1 = optfmtin
   fmtin2 = optfmtin
end if

! Read coordinates

call readfile(read_unit1, fmtin1, mol1)
call readfile(read_unit2, fmtin2, mol2)

! Allocate arrays

allocate (maplist(mol1%natom, maxrec))
allocate (countlist(maxrec))

if (pipe_flag) then
   write_unit = stdout
   fmtout = 'xyz'
else
   call open2write(pathout, write_unit, fmtout)
end if

if (fmtout_flag) then
   fmtout = optfmtout
end if

! Compute types

call compute_eltypes(mol1, mol2, eltypes)

if (remap_flag) then

   ! Remap atoms to minimize the MSD

   call molecule_remap( &
      mol1, &
      mol2, &
      eltypes, &
      nrec, &
      maplist, &
      countlist)

   minrmsd = huge(rmsd)
   minadjd = huge(adjd)

   if (.not. nrec_flag) then
      call writefile(write_unit, fmtout, mol1)
   end if

   do i = 1, nrec

      call remapped_molecule_align( &
         mol1, &
         mol2, &
         maplist(:, i), &
         travec1, &
         travec2, &
         rotquat)

      auxmol1 = mol2
      call auxmol1%permutate_atoms(maplist(:, i))
      call auxmol1%translate_coords(travec2)
      call auxmol1%rotate_coords(rotquat)
      call auxmol1%translate_coords(-travec1)

!      atomperm = maplist(:, i)
      atomperm = identity_perm(size(mol1%atoms))

      rmsd = get_rmsd(mol1, auxmol1, atomperm)
      adjd = get_adjd(mol1, auxmol1, atomperm)

      minrmsd = min(minrmsd, rmsd)
      minadjd = min(minadjd, adjd)

      auxmol1%title = 'Map='//intstr(i)//' RMSD='//realstr(rmsd, 4)
      call writefile(write_unit, fmtout, auxmol1)

   end do

else

   ! Align atoms to minimize RMSD

   call molecule_align( &
      mol1, &
      mol2, &
      eltypes, &
      travec1, &
      travec2, &
      rotquat)

   call mol2%translate_coords(travec2)
   call mol2%rotate_coords(rotquat)
   call mol2%translate_coords(-travec1)

   atomperm = identity_perm(size(mol1%atoms))

   minrmsd = get_rmsd(mol1, mol2, atomperm)
   minadjd = get_adjd(mol1, mol2, atomperm)

   mol2%title = 'RMSD='//realstr(rmsd, 4)
   call writefile(write_unit, fmtout, mol2)

end if

if (print_flag) then
   if (bond_flag) then
      write (stderr, "(a,1x,i0)") realstr(minrmsd, 4), minadjd
   else
      write (stderr, "(a)") realstr(minrmsd, 4)
   end if
end if

end program
