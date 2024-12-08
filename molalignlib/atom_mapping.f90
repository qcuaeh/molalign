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

module atom_mapping
use kinds
use flags
use bounds
use molecule
use strutils
use chemdata
use permutation
use rigid_body
use rotation
use alignment
use lap_driver
use adjacency
use bipartition
use bipartitioning
use biasing
use pruning
use printing
!use backtracking

implicit none

contains

subroutine map_atoms(mol1, mol2, eltypes, permlist, countlist, nrec)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_type), intent(in) :: eltypes
   integer, intent(out) :: permlist(:, :)
   integer, intent(out) :: countlist(:)
   integer, intent(out) :: nrec

   ! Local variables

   logical :: visited, overflow
   integer :: irec, krec, ntrial, nstep, steps
   integer, dimension(mol1%natom) :: atomperm, newmapping
   real(rk) :: rmsd, totalrot
   real(rk) :: workcoords(3, mol2%natom)
   real(rk), dimension(4) :: rotquat, prodquat
   real(rk), dimension(maxrec) :: rmsdlist, avgsteps, avgtotalrot, avgrealrot
   type(boolmatrix_type), dimension(:), allocatable :: prunemask
   type(realmatrix_type), dimension(:), allocatable :: mnadists

   integer :: natom1
   real(rk), dimension(:, :), allocatable :: coords1, coords2
   real(rk), dimension(:), allocatable :: weights1, weights2
   real(rk) :: travec1(3), travec2(3)

   natom1 = size(mol1%atoms)
   coords1 = mol1%get_coords()
   coords2 = mol2%get_coords()
   weights1 = element_weights(mol1%atoms%elnum)
   weights2 = element_weights(mol1%atoms%elnum)

   ! Mirror coordinates

   if (mirror_flag) then
      call reflect_coords(coords2)
   end if

   ! Calculate centroids

   travec1 = -centroid(coords1, weights1)
   travec2 = -centroid(coords2, weights2)

   ! Center coordinates at the centroids

   call translate_coords(coords1, travec1)
   call translate_coords(coords2, travec2)

   ! Calculate prune matrix

   call prune_procedure(eltypes, coords1, coords2, prunemask)

   ! Calculate bias matrix

   call bias_procedure(mol1, mol2, eltypes, mnadists)

   ! Initialize loop variables

   nrec = 0
   nstep = 0
   ntrial = 0
   countlist(1) = 0
   overflow = .false.

   ! Loop for map searching

   do while (countlist(1) < maxcount .and. ntrial < maxtrials)

      ntrial = ntrial + 1

      ! Work with a copy of coords2

      workcoords = coords2

      ! Aply a random rotation to workcoords

      call rotate_coords(workcoords, randrotquat())

      ! Minimize the euclidean distance

      call assign_atoms_generic(eltypes, coords1, workcoords, prunemask, mnadists, atomperm)
      rotquat = leastrotquat(natom1, weights1, coords1, workcoords, atomperm)
      prodquat = rotquat
      totalrot = rotangle(rotquat)
      call rotate_coords(workcoords, rotquat)
!      print *, sqrt(leastsquaredist(natom1, weights1, coords1, coords2, atomperm))
!      stop

      steps = 1

      do while (iter_flag)
         call assign_atoms_generic(eltypes, coords1, workcoords, prunemask, mnadists, newmapping)
         if (all(newmapping == atomperm)) exit
         rotquat = leastrotquat(natom1, weights1, coords1, workcoords, newmapping)
         prodquat = quatmul(rotquat, prodquat)
         call rotate_coords(workcoords, rotquat)
         totalrot = totalrot + rotangle(rotquat)
         atomperm = newmapping
         steps = steps + 1
      end do

      nstep = nstep + steps

      rmsd = sqrt(leastsquaredist(natom1, weights1, coords1, coords2, atomperm))

      ! Check for new best permlist

      visited = .false.

      do irec = 1, nrec
         if (all(atomperm == permlist(:, irec))) then
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
            if (rmsd < rmsdlist(irec)) then
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
               rmsdlist(irec) = rmsdlist(irec - 1)
               avgsteps(irec) = avgsteps(irec - 1)
               avgrealrot(irec) = avgrealrot(irec - 1)
               avgtotalrot(irec) = avgtotalrot(irec - 1)
               permlist(:, irec) = permlist(:, irec - 1)
            end do
            countlist(krec) = 1
            rmsdlist(krec) = rmsd
            avgsteps(krec) = steps
            avgrealrot(krec) = rotangle(prodquat)
            avgtotalrot(krec) = totalrot
            permlist(:, krec) = atomperm
         end if
      end if

   end do

   ! Print stats if requested

   if (stats_flag) then
      print_flag = .false.
      call print_stats(nrec, countlist, avgsteps, avgrealrot, rmsdlist)
      call print_final_stats(overflow, maxrec, nrec, ntrial, nstep)
   end if

end subroutine

end module
