! template for module backtracking in molalign program
module backtracking
use stdio
use sorting
use molecule
use partition
use permutation
use adjacency
use alignment
use random

implicit none
private
public minadjdiff
public eqvatomperm

logical, parameter :: printInfo = .false.

contains

! Find best correspondence between points of graphs
subroutine minadjdiff (mol1, mol2, atomperm)
   type(mol_type), intent(in) :: mol1, mol2
   integer, dimension(:), intent(inout) :: atomperm

   ! Local variables

   integer :: eltypemap1(mol1%natom)
   integer :: h, i, moldiff
   integer :: ntrack, track(mol1%natom)
   integer :: unmapping(mol1%natom)
   integer, dimension(mol1%natom) :: mnatypemap1, mnatypemap2
   logical :: tracked(mol1%natom)
   real(rk) moldist

   integer :: natom
   integer :: neltype1
   integer :: nmnatype1, nmnatype2

   type(bipartition_type) :: eltypes1
   type(bipartition_type) :: mnatypes1, mnatypes2

   integer, allocatable, dimension(:) :: eltypepartsizes1
   integer, allocatable, dimension(:) :: mnatypepartsizes1, mnatypepartsizes2
   integer, allocatable, dimension(:) :: fragroots1, fragroots2

   logical, dimension(:, :), allocatable :: adjmat1, adjmat2
   real(rk), dimension(:, :), allocatable :: coords1, coords2
   real(rk), dimension(:), allocatable :: weights

   integer :: nadjs1(mol1%natom)
   integer :: nadjs2(mol2%natom)
   integer :: adjlists1(maxcoord, mol1%natom)
   integer :: adjlists2(maxcoord, mol2%natom)
   integer :: nadjmnatypes1(mol1%natom)
   integer :: nadjmnatypes2(mol2%natom)
   integer :: adjmnatypepartlens1(maxcoord, mol1%natom)
   integer :: adjmnatypepartlens2(maxcoord, mol2%natom)

   natom = mol1%get_natom()

   eltypes1 = mol1%gather_eltypes()
   eltypemap1 = eltypes1%partition_map

   neltype1 = eltypes1%size
   allocate (eltypepartsizes1(neltype1))
   do h = 1, neltype1
      eltypepartsizes1(h) = eltypes1%parts(h)%size
   end do

   mnatypes1 = mol1%gather_mnatypes()
   mnatypes2 = mol2%gather_mnatypes()
   mnatypemap1 = mnatypes1%partition_map
   mnatypemap2 = mnatypes2%partition_map

   nmnatype1 = size(mnatypes1%parts)
   nmnatype2 = size(mnatypes2%parts)
   allocate (mnatypepartsizes1(nmnatype1))
   allocate (mnatypepartsizes2(nmnatype2))
   do h = 1, nmnatype1
      mnatypepartsizes1(h) = mnatypes1%parts(h)%size
   end do
   do h = 1, nmnatype2
      mnatypepartsizes2(h) = mnatypes2%parts(h)%size
   end do

!   adjlists1 = mol1%gather_adjlists()
!   adjlists2 = mol2%gather_adjlists()
   nadjs1 = mol1%gather_nadjs()
   nadjs2 = mol2%gather_nadjs()
   adjlists1 = mol1%gather_adjlists()
   adjlists2 = mol2%gather_adjlists()

!   adjpartitions0 = mol1%gather_adjpartitions()
!   adjpartitions1 = mol2%gather_adjpartitions()
   nadjmnatypes1 = mol1%gather_nadjmnatypes()
   nadjmnatypes2 = mol2%gather_nadjmnatypes()
   adjmnatypepartlens1 = mol1%gather_adjmnatypepartlens()
   adjmnatypepartlens2 = mol2%gather_adjmnatypepartlens()

   coords1 = mol1%gather_coords()
   coords2 = mol2%gather_coords()
   weights = mol1%gather_weights()
   adjmat1 = mol1%gather_adjmatrix()
   adjmat2 = mol2%gather_adjmatrix()
   fragroots1 = mol1%gather_molfragroots()
   fragroots2 = mol2%gather_molfragroots()

   !  initialization

   ntrack = 0
   tracked(:) = .false.
   unmapping = inverse_perm(atomperm)
   moldiff = adjacencydiff (natom, adjmat1, adjmat2, atomperm)
   moldist = squaredist (natom, weights, coords1, coords2, atomperm)

   if ( printInfo ) then
      print '(a,1x,i0)', "moldiff:", moldiff
      print '(a,1x,f0.4)', "moldist:", moldist
   end if

   do i = 1, size(fragroots1)
      call recursive_backtrack (fragroots1(i), atomperm, unmapping, tracked, moldiff, moldist, ntrack, track)
!        print *, ntrack
   end do

   if ( printInfo ) then
      print '(a,1x,i0)', "countFrag:", size(fragroots1)
      print '(a,1x,i0,1x,i0)', "moldiff:", adjacencydiff (natom, adjmat1, adjmat2, atomperm), moldiff
      print '(a,1x,f0.4,1x,f0.4)', "moldist:", squaredist (natom, weights, coords1, coords2, atomperm), moldist
   end if

!    if (adjacencydiff (natom, adjmat1, adjmat2, atomperm) /= moldiff) then
!        print '(a,x,i0,x,i0)', "moldiff:", adjacencydiff (natom, adjmat1, adjmat2, atomperm), moldiff
!    end if

   contains    

! classify the atoms connected to node as matches or unmatched
   subroutine nodematch (node, atomperm, tracked, nmatch, matches, &
                   nmismatch0, mismatches0, nmismatch1, mismatches1)
      integer, intent(in) :: node, atomperm(:)
      logical, intent(in) :: tracked(:)
      integer, intent(out) :: nmatch, matches(natom)
      integer, intent(out) :: nmismatch0, mismatches0(natom)
      integer, intent(out) :: nmismatch1, mismatches1(natom)
      integer :: i, adj0(natom), adj1(natom)

      adj0(:nadjs1(node)) = adjlists1(:nadjs1(node), node)
      adj1(:nadjs2(atomperm(node))) = adjlists2(:nadjs2(atomperm(node)), atomperm(node))
      nmatch = 0
      nmismatch0 = 0
      nmismatch1 = 0

      do i = 1, nadjs1(node)
         if (findInd(atomperm(adj0(i)), nadjs2(atomperm(node)), adj1) /= 0) then
            call addInd(adj0(i), nmatch, matches)
         else
            call addInd(adj0(i), nmismatch0, mismatches0)
         end if
      end do

      do i = 1, nadjs2(atomperm(node))
         if (findInd(adj1(i), nmatch, atomperm(matches(:nmatch))) == 0) then
            call addInd(adj1(i), nmismatch1, mismatches1)
         end if
      end do
   end subroutine nodematch

! backtracks structure to find assignments that minimize moldiff
   recursive subroutine recursive_backtrack (node, atomperm, unmapping, tracked, moldiff, moldist, &
                                  ntrack, track)
      integer, intent(in) :: node
      integer, intent(inout) :: atomperm(natom), unmapping(natom)
      logical, intent(inout) :: tracked(natom)
      integer, intent(inout) :: moldiff, ntrack, track(natom)
      real(rk), intent(inout) :: moldist

      logical, parameter :: printInfo = .false.

      integer :: nmatch, matches(natom)
      integer :: nmismatch0, mismatches0(natom)
      integer :: nmismatch1, mismatches1(natom)
      integer :: mapping_branch(natom), moldiff_branch, unmapping_branch(natom)
      integer :: ntrack_branch, track_branch(natom)
      logical :: tracked_branch(natom)
      integer :: i, j
      logical :: matched0(natom), matched1(natom)
      real(rk) :: moldist_branch, mol0local(3,natom), mol1local(3,natom)

      ! reserve node node as tracked
      ntrack = ntrack + 1
      track(ntrack) = node
      tracked(node) = .true.

      ! classify neighbor atoms as matches or mismatched for coords1/coords2
      call nodematch (node, atomperm, tracked, nmatch, matches, &
                  nmismatch0, mismatches0, nmismatch1, mismatches1)

      if (printInfo) then
         print *, "node:", node
         print *, "matches:", matches(:nmatch)
         print *, "mismatches0:", mismatches0(:nmismatch0)
         print *, "mismatches1:", mismatches1(:nmismatch1)
      end if

      ! shuffle indices
      call shuffle (matches(:nmatch))
      call shuffle (mismatches0(:nmismatch0))
      call shuffle (mismatches1(:nmismatch1))

      ! run over matches neighbors
      do i = 1, nmatch
         if (.not. tracked(matches(i))) then
            call recursive_backtrack(matches(i), atomperm, unmapping, tracked, moldiff, &
                            moldist, ntrack, track)
         end if
      end do

      matched0(:nmismatch0) = .false.
      matched1(:nmismatch1) = .false.

      ! run over mismatched neighbors
      do i = 1, nmismatch0
         if (.not. tracked(mismatches0(i))) then
            do j = 1, nmismatch1
               if (.not. matched1(j)) then
                  if (eltypemap1(mismatches0(i)) == eltypemap1(mismatches1(j))) then

                     ntrack_branch = ntrack
                     track_branch(:) = track(:)
                     tracked_branch(:) = tracked(:)
                     mapping_branch(:) = atomperm(:)
                     unmapping_branch(:) = unmapping(:)

                     ! Apply swap to atomperm branch 
                     mapping_branch(mismatches0(i)) = mismatches1(j)
                     mapping_branch(unmapping(mismatches1(j))) = atomperm(mismatches0(i))

                     ! Apply swap to unmapping branch 
                     unmapping_branch(mismatches1(j)) = mismatches0(i)
                     unmapping_branch(atomperm(mismatches0(i))) = unmapping(mismatches1(j))

                     ! Update ssd with swap
                     moldist_branch = moldist + weights(mismatches0(i))*( &
                        - sum((coords2(:, atomperm(mismatches0(i))) - coords1(:, mismatches0(i)))**2) &
                        - sum((coords2(:, mismatches1(j)) - coords1(:, unmapping(mismatches1(j))))**2) &
                        + sum((coords2(:, mismatches1(j)) - coords1(:, mismatches0(i)))**2) &
                        + sum((coords2(:, atomperm(mismatches0(i))) - coords1(:, unmapping(mismatches1(j))))**2))

                     ! Update adjd with swap
                     moldiff_branch = moldiff + adjacencydelta(nadjs1, adjlists1, adjmat2, &
                                   atomperm, mismatches0(i), unmapping(mismatches1(j)))

                     ! backtrack swapped index
                     call recursive_backtrack(mismatches0(i), mapping_branch, unmapping_branch, &
                        tracked_branch, moldiff_branch, moldist_branch, ntrack_branch, track_branch)

                     if ( &
                        moldiff_branch < moldiff &
                        .and. ( &
                           mnatypemap1(mismatches0(i)) == mnatypemap1(unmapping(mismatches1(j))) &
                           .and. mnatypemap2(atomperm(mismatches0(i))) == mnatypemap2(mismatches1(j)) &
                        ) &
                     ) then
                        ntrack = ntrack_branch
                        track(:) = track_branch(:)
                        tracked(:) = tracked_branch(:)
                        atomperm(:) = mapping_branch(:)
                        unmapping(:) = unmapping_branch(:)
                        moldiff = moldiff_branch
                        moldist = moldist_branch
                        matched0(i) = .true.
                        matched1(j) = .true.
                        exit   ! exits inner do loop
                     end if
                  end if
               end if
            end do
         end if
      end do

      ! run over non matches neighbors
      do i = 1, nmismatch0
         if (.not. matched0(i)) then
            if (.not. tracked(mismatches0(i))) then
               call recursive_backtrack(mismatches0(i), atomperm, unmapping, &
                 tracked, moldiff, moldist, ntrack, track)
            end if
         end if
      end do

   end subroutine

end subroutine

! Find best correspondence between points of graphs
subroutine eqvatomperm (mol1, mol2, workcoords, atomperm)
   type(mol_type), intent(in) :: mol1, mol2
   real(rk), intent(in) :: workcoords(:, :)
   integer, intent(inout) :: atomperm(:)

   ! Local variables

   integer :: unmap(mol1%natom), track(mol1%natom)
   integer, dimension(mol1%natom) :: eqvos, eqvidx
   integer :: permcount, h, i, n, diff, fragcount, ntrack
   logical :: tracked(mol1%natom), held(mol1%natom)
   real(rk) :: dist

   integer :: natom
   integer :: neltype1
   integer :: nmnatype1, nmnatype2

   type(bipartition_type) :: eltypes1
   type(bipartition_type) :: mnatypes1, mnatypes2

   integer, allocatable, dimension(:) :: eltypepartsizes1
   integer, allocatable, dimension(:) :: mnatypepartsizes1, mnatypepartsizes2
   integer, allocatable, dimension(:) :: fragroots1, fragroots2

   logical, dimension(:, :), allocatable :: adjmat1, adjmat2
   real(rk), dimension(:, :), allocatable :: coords1, coords2
   real(rk), dimension(:), allocatable :: weights

   integer :: nadjs1(mol1%natom)
   integer :: nadjs2(mol2%natom)
   integer :: adjlists1(maxcoord, mol1%natom)
   integer :: adjlists2(maxcoord, mol2%natom)
   integer :: nadjmnatypes1(mol1%natom)
   integer :: nadjmnatypes2(mol2%natom)
   integer :: adjmnatypepartlens1(maxcoord, mol1%natom)
   integer :: adjmnatypepartlens2(maxcoord, mol2%natom)

   natom = mol1%get_natom()

   eltypes1 = mol1%gather_eltypes()

   neltype1 = eltypes1%size
   allocate (eltypepartsizes1(neltype1))
   do h = 1, neltype1
      eltypepartsizes1(h) = eltypes1%parts(h)%size
   end do

   mnatypes1 = mol1%gather_mnatypes()
   mnatypes2 = mol2%gather_mnatypes()

   nmnatype1 = size(mnatypes1%parts)
   nmnatype2 = size(mnatypes2%parts)
   allocate (mnatypepartsizes1(nmnatype1))
   allocate (mnatypepartsizes2(nmnatype2))
   do h = 1, nmnatype1
      mnatypepartsizes1(h) = mnatypes1%parts(h)%size
   end do
   do h = 1, nmnatype2
      mnatypepartsizes2(h) = mnatypes2%parts(h)%size
   end do

!   adjlists1 = mol1%gather_adjlists()
!   adjlists2 = mol2%gather_adjlists()
   nadjs1 = mol1%gather_nadjs()
   nadjs2 = mol2%gather_nadjs()
   adjlists1 = mol1%gather_adjlists()
   adjlists2 = mol2%gather_adjlists()

!   adjpartitions0 = mol1%gather_adjpartitions()
!   adjpartitions1 = mol2%gather_adjpartitions()
   nadjmnatypes1 = mol1%gather_nadjmnatypes()
   nadjmnatypes2 = mol2%gather_nadjmnatypes()
   adjmnatypepartlens1 = mol1%gather_adjmnatypepartlens()
   adjmnatypepartlens2 = mol2%gather_adjmnatypepartlens()

   coords1 = mol1%gather_coords()
   coords2 = mol2%gather_coords()
   weights = mol1%gather_weights()
   adjmat1 = mol1%gather_adjmatrix()
   adjmat2 = mol2%gather_adjmatrix()
   fragroots1 = mol1%gather_molfragroots()
   fragroots2 = mol2%gather_molfragroots()

   ! set equivalence group offsets

   eqvos(1) = 0
   do h = 1, nmnatype1 - 1
      eqvos(h+1) = eqvos(h) + mnatypepartsizes1(h)
   end do

   ! set atoms equivalence indices

   do h = 1, nmnatype1
      eqvidx(eqvos(h)+1:eqvos(h)+mnatypepartsizes1(h)) = h
   end do

   ! initialization

   permcount = 1
   fragcount = 0
   tracked(:) = .false.
   held(:) = .false.
   ntrack = 0

!    i = 1
!    do while ( i <= natom )
!        if ( tracked(i) ) then
!            i = i + 1
!        else
!            call recursive_permut (i, atomperm, tracked, held, 0, 0)
!            i = 1   ! restart loop to find unconnected, untracked atoms
!            fragcount = fragcount + 1
!        end if
!    end do

   do i = 1, size(fragroots1)
      call recursive_permut (fragroots1(i), atomperm, tracked, held, ntrack, track)
!        print *, ntrack
   end do

!   print '(a,i0)', "Natoms: ", natom
!   print '(a,f8.4)', "dist: ", sqrt(leastsquaredist (natom, weights, coords1, workcoords, atomperm)/sum(weights))
!   print '(a,i0)', "permcount: ", permcount
!   print '(a,i0)', "fragcount: ", fragcount

   contains    

   recursive subroutine recursive_remap (nodea, nodeb, mapping_ref, atomperm, held)
      integer, intent(in) :: nodea, nodeb, mapping_ref(natom)
      integer, intent(inout) :: atomperm(natom)
      logical, dimension(natom), intent(inout) :: held
      logical, dimension(natom) :: locked_c
      integer :: meqvnei, equiva(natom), equivb(natom)
      integer :: h, i, offset, first, last

      first = eqvos(eqvidx(nodea)) + 1
      last = eqvos(eqvidx(nodea)) + mnatypepartsizes1(eqvidx(nodea))
      held(first:last) = .true.
      offset = 0
      do h = 1, nadjmnatypes1(nodea)
         ! find not tracked adjlists1 in group
         meqvnei = 0
         do i = 1, adjmnatypepartlens1(h, nodea)
            if (.not. held(adjlists1(offset+i, nodea))) then 
               meqvnei = meqvnei + 1
               equiva(meqvnei) = adjlists1(offset+meqvnei, nodea)
               equivb(meqvnei) = adjlists1(offset+meqvnei, nodeb)
            end if
         end do
         locked_c(:) = held(:)
         do i = 1, meqvnei
            if ( equiva(i) /= equivb(i) ) then
               atomperm(equiva(i)) = mapping_ref(equivb(i))
               call recursive_remap (equiva(i), equivb(i), mapping_ref, &
                  atomperm, locked_c)
            end if
         end do
         offset = offset + adjmnatypepartlens1(h, nodea)
      end do
   end subroutine recursive_remap

   recursive subroutine recursive_permut (node, atomperm, tracked, held, ntrack, track)
      integer, intent(in) :: node
      integer, intent(inout) :: ntrack
      logical, dimension(natom), intent(inout) :: tracked, held
      integer, dimension(natom), intent(inout) :: atomperm, track

      logical :: locked_c(natom)
      integer :: meqvnei, moldiff, track4ind(4), track4ind_c(4)
      integer, dimension(natom) :: mapping_p, mapping_min, equiv, perm, perm_min
      real(rk) :: moldist, moldist_p, moldist_min, dihed0(maxcoord), dihed1(maxcoord)
      logical :: more, calcd, printInfo = .false.
      integer :: h, i, j, offset, first, last, rank
      character(len=80) :: strfmt

      if ( printInfo ) then   ! print debugging info
         moldist = sqrt(leastsquaredist(natom, weights, coords1, workcoords, atomperm)/sum(weights))
         moldiff = adjacencydiff(natom, adjmat1, adjmat2, atomperm)
         write (strfmt, '(a,i0,a)') '(',1,'(2x),i0,a,i0,f8.4)'
         print strfmt, node,": ",moldiff,moldist
      end if

      ! reserves node as tracked
      tracked(node) = .true.
      ntrack = ntrack + 1
      track(ntrack) = node

      ! reserves atoms with the same type as n
      first = eqvos(eqvidx(node)) + 1
      last = eqvos(eqvidx(node)) + mnatypepartsizes1(eqvidx(node))
      held(first:last) = .true.

      offset = 0
      ! run over groups of atoms with equivalent type
      do h = 1, nadjmnatypes1(node)
         ! find not tracked adjlists1 in group
         meqvnei = 0
         do i = 1, adjmnatypepartlens1(h, node)
            if (.not. tracked(adjlists1(offset+i, node)) &
               .and. .not. held(adjlists1(offset+i, node))) then
               meqvnei = meqvnei + 1
               equiv(meqvnei) = adjlists1(offset+meqvnei, node)
            end if
         end do
         ! check permutations
         if ( meqvnei > 0 ) then
            ! save initial state before permutations
            mapping_min(:) = atomperm(:)
            ! run over all permutarions for the meqvnei equivalent atoms
            more = .false.
            call perm1_next3 (meqvnei, perm, more, rank)
            perm_min(:) = perm(:)
            ! initialize state to test new permutation
            mapping_min(:) = atomperm(:)
            moldist_min = 0.0
            ! apply permutation
            do i = 1, meqvnei
               mapping_min(equiv(i)) = atomperm(equiv(perm_min(i)))
               moldist_min = moldist_min &
                        + sum((workcoords(:, mapping_min(equiv(i))) - coords1(:, equiv(i)))**2)
            end do
            do while ( more )
               call perm1_next3 (meqvnei, perm, more, rank)
               permcount = permcount + 1
               ! initialize state to test new permutation
               mapping_p(:) = atomperm(:)
               moldist_p = 0.0
               ! apply permutation
               do i = 1, meqvnei
                  mapping_p(equiv(i)) = atomperm(equiv(perm(i)))
                  moldist_p = moldist_p &
                         + sum((workcoords(:, mapping_p(equiv(i))) - coords1(:, equiv(i)))**2) 
               end do
               ! save min dist permut
               if ( moldist_p < moldist_min ) then
                  mapping_min(:) = mapping_p(:)
                  moldist_min = moldist_p
                  perm_min(:) = perm(:)
               end if
            end do
            ! check chirality
!                if ( printInfo ) then 
!
!                    if ( meqvnei >= 3 ) then
!                        track4ind(1) = equiv(1)
!                        track4ind(2) = node
!                        track4ind(3) = track(ntrack)
!
!                        do j = 1, meqvnei
!                            track4ind(4) = equiv(j)
!                            call calc_dihedral (0, track4ind, coords1, calcd, dihed0(j))
!                            track4ind_c(:) = mapping_min(track4ind(:))
!                            call calc_dihedral (0, track4ind_c, workcoords, calcd, dihed1(j))
!                        end do
!                        print '(a,3f10.3)', "coords1: ", dihed0(:meqvnei)
!                        print '(a,3f10.3)', "workcoords: ", dihed1(:meqvnei)
!                    end if
!                end if
            ! remap connectivity of permuted branch path
            locked_c(:) = held(:)
            do i = 1, meqvnei
               if ( equiv(i) /= equiv(perm_min(i)) ) then
                  call recursive_remap (equiv(i), equiv(perm_min(i)), &
                               atomperm, mapping_min, locked_c)
               end if
            end do
            ! update atomperm
            atomperm(:) = mapping_min(:)
            do i = 1, meqvnei
               call recursive_permut (equiv(i), atomperm, tracked, &
                                locked_c, ntrack, track)
            end do
         end if
         offset = offset + adjmnatypepartlens1(h, node)
      end do
   end subroutine

end subroutine

! Adds a new element to a list of indices in growing order, no repeats
subroutine addInd(ind, nInd, vecInd)
   integer, intent(in) :: ind
   integer, intent(inout) :: nInd, vecInd(:)
   integer :: p, e
   p = 1
   do while (p <= nInd)   ! find position for growing order
      if (ind == vecInd(p)) return
      if (ind < vecInd(p)) then
         exit
      else
         p = p + 1
      end if
   end do
   nInd = nInd + 1
   e = nInd
   do while (e > p)
      vecInd(e) = vecInd(e - 1)   ! shifts one position from the end
      e = e - 1
   end do
   vecInd(p) = ind
end subroutine addInd

! Returns the position in vecInd where the element ind is found or 0 otherwise
function findInd(ind, nInd, vecInd) result(pos)
   integer, intent(in) :: ind, nInd, vecInd(nInd)
   integer :: pos
   integer :: i
   pos = 0
   do i = 1, nInd
      if (ind == vecInd(i)) then
         pos = i
         return
      end if
   end do
end function findInd

end module
