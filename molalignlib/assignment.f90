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

module assignment
use kinds
use flags
use bounds
use molecule
use strutils
use chemdata
use permutation
use translation
use rotation
use alignment
use lap_driver
use adjacency
use assorting
use biasing
use pruning
use printing
!use backtracking

implicit none

contains

subroutine remap_atoms(mol1, mol2, eltypes, maplist, countlist, nrec)
   type(molecule_type), intent(in) :: mol1, mol2
   type(partition_type), intent(in) :: eltypes
   integer, intent(out) :: maplist(:, :)
   integer, intent(out) :: countlist(:)
   integer, intent(out) :: nrec

   ! Local variables

   logical :: visited, overflow
   integer :: irec, krec, ntrial, nstep, steps
   integer :: adjd, recadjd(maxrec)
   integer, dimension(mol1%natom) :: atomperm, newmapping
   real(rk) :: rmsd, totalrot
   real(rk) :: workcoords(3, mol2%natom)
   real(rk), dimension(4) :: rotquat, prodquat
   real(rk), dimension(maxrec) :: recrmsd, avgsteps, avgtotalrot, avgrealrot
   type(boolmatrix_type), dimension(:), allocatable :: prunemask
   type(realmatrix_type), dimension(:), allocatable :: mnadists

   integer :: natom1
   type(atom_type), dimension(:), allocatable :: atoms1, atoms2
   logical, dimension(:, :), allocatable :: adjmat1, adjmat2
   real(rk), dimension(:, :), allocatable :: coords1, coords2
   real(rk), dimension(:), allocatable :: weights1

   natom1 = mol1%natom
   atoms1 = mol1%atoms
   atoms2 = mol2%atoms
   coords1 = mol1%get_coords()
   coords2 = mol2%get_coords()
   adjmat1 = mol1%get_adjmat()
   adjmat2 = mol2%get_adjmat()
   weights1 = weights(mol1%atoms%elnum)

   ! Calculate prune matrix

   call prune_procedure(eltypes, coords1, coords2, prunemask)

   ! Calculate bias matrix

   call bias_procedure(atoms1, atoms2, eltypes, mnadists)

   ! Initialize loop variables

   nrec = 0
   nstep = 0
   ntrial = 0
   countlist(1) = 0
   overflow = .false.

   ! Loop for map searching

   do while (countlist(1) < maxcount .and. (.not. trial_flag .or. ntrial < maxtrials))

      ntrial = ntrial + 1

      ! Work with a copy of coords2

      workcoords = coords2

      ! Aply a random rotation to workcoords

      call rotate(natom1, workcoords, randrotquat())

      ! Minimize the euclidean distance

      call assign_atoms(eltypes, coords1, workcoords, prunemask, mnadists, atomperm)
      rotquat = leastrotquat(natom1, weights1, coords1, workcoords, atomperm)
      prodquat = rotquat
      totalrot = rotangle(rotquat)
      call rotate(natom1, workcoords, rotquat)
!      print *, sqrt(leastsquaredist(natom1, weights1, coords1, coords2, atomperm)/sum(weights1))
!      stop

      steps = 1

      do while (iter_flag)
         call assign_atoms(eltypes, coords1, workcoords, prunemask, mnadists, newmapping)
         if (all(newmapping == atomperm)) exit
         rotquat = leastrotquat(natom1, weights1, coords1, workcoords, newmapping)
         prodquat = quatmul(rotquat, prodquat)
         call rotate(natom1, workcoords, rotquat)
         totalrot = totalrot + rotangle(rotquat)
         atomperm = newmapping
         steps = steps + 1
      end do

      nstep = nstep + steps

!      if (back_flag) then
!         call minadjdiff(mol1, mol2, atomperm)
!         call eqvatomperm(mol1, mol2, workcoords, atomperm)
!      end if

      adjd = adjacencydiff(natom1, adjmat1, adjmat2, atomperm)
      rmsd = sqrt(leastsquaredist(natom1, weights1, coords1, coords2, atomperm)/sum(weights1))

      ! Check for new best maplist

      visited = .false.

      do irec = 1, nrec
         if (all(atomperm == maplist(:, irec))) then
            countlist(irec) = countlist(irec) + 1
            avgsteps(irec) = avgsteps(irec) + (steps - avgsteps(irec))/countlist(irec)
            avgrealrot(irec) = avgrealrot(irec) + (rotangle(prodquat) - avgrealrot(irec))/countlist(irec)
            avgtotalrot(irec) = avgtotalrot(irec) + (totalrot - avgtotalrot(irec))/countlist(irec)
            visited = .true.
            exit
         end if
      end do

      if (.not. visited) then
         krec = nrec + 1
         do irec = nrec, 1, -1
            if (adjd < recadjd(irec) .or. (adjd == recadjd(irec) .and. rmsd < recrmsd(irec))) then
               krec = irec
            else
               exit
            end if
         end do
         if (nrec < maxrec) then
            nrec = nrec + 1
         else
            overflow = .true.
         end if
         if (krec <= maxrec) then
            do irec = nrec, krec + 1, -1
               countlist(irec) = countlist(irec - 1)
               recrmsd(irec) = recrmsd(irec - 1)
               recadjd(irec) = recadjd(irec - 1)
               avgsteps(irec) = avgsteps(irec - 1)
               avgrealrot(irec) = avgrealrot(irec - 1)
               avgtotalrot(irec) = avgtotalrot(irec - 1)
               maplist(:, irec) = maplist(:, irec - 1)
            end do
            countlist(krec) = 1
            recrmsd(krec) = rmsd
            recadjd(krec) = adjd
            avgsteps(krec) = steps
            avgrealrot(krec) = rotangle(prodquat)
            avgtotalrot(krec) = totalrot
            maplist(:, krec) = atomperm
         end if
      end if

   end do

   ! Reorder back to default atom ordering

!   do irec = 1, nrec
!      maplist(:, irec) = mol2%atomorder(maplist(mol1%atomordermap, irec))
!   end do

   ! Print stats if requested

   if (stats_flag) then
      print_flag = .false.
      call print_stats(nrec, countlist, avgsteps, avgtotalrot, avgrealrot, recadjd, recrmsd)
      call print_final_stats(overflow, maxrec, nrec, ntrial, nstep)
   end if

end subroutine

end module
