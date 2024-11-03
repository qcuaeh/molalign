module adjacency
use stdio
use kinds
use flags
use bounds
use sorting
use chemdata

implicit none

contains

subroutine adjmat2list(natom, adjmat, nadjs, adjlists)
   integer, intent(in) :: natom
   logical, dimension(:, :), intent(in) :: adjmat
   integer, dimension(:), intent(out) :: nadjs
   integer, dimension(:, :), intent(out) :: adjlists
   integer :: i, j

   nadjs(:) = 0

   do i = 1, natom
      do j = i + 1, natom
         if (adjmat(i, j)) then
            nadjs(i) = nadjs(i) + 1
            nadjs(j) = nadjs(j) + 1
            if (nadjs(i) > max_coord_num .or. nadjs(j) > max_coord_num) then
               write (stderr, '("Maximum coordination number exceeded!")')
               stop
            end if
            adjlists(nadjs(i), i) = j
            adjlists(nadjs(j), j) = i
         end if
      end do
   end do

end subroutine

subroutine adjmat2bonds(natom, adjmat, nbond, bonds)
   integer, intent(in) :: natom
   logical, dimension(:, :), intent(in) :: adjmat
   integer, intent(out) :: nbond
   integer, dimension(:, :), intent(out) :: bonds
   integer :: i, j

   nbond = 0

   do i = 1, natom
      do j = i + 1, natom
         if (adjmat(i, j)) then
            nbond = nbond + 1
            bonds(1, nbond) = i
            bonds(2, nbond) = j
         end if
      end do
   end do

end subroutine

function adjacencydiff(natom, adjmat1, adjmat2, atomperm) result(diff)
! Purpose: Check if two graphs are equal.
! Return the number of differences between graphs.
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: atomperm
   logical, dimension(:, :), intent(in) :: adjmat1, adjmat2
   integer :: diff

   integer :: i, j

   diff = 0

! Check differences element by element

   do i = 1, natom
      do j = i + 1, natom
         if (adjmat1(i, j) .neqv. adjmat2(atomperm(i), atomperm(j))) then
            diff = diff + 1
!                print *, i, j, adjmat1(i, j), adjmat2(atomperm(i), atomperm(j))
         end if
      end do
   end do

end function

function adjacencydelta(nadjs1, adjlists1, adjmat2, atomperm, k, l) result(delta)
   integer, intent(in) :: k, l
   integer, dimension(:), intent(in) :: atomperm, nadjs1
   integer, dimension(:, :), intent(in) :: adjlists1
   logical, dimension(:, :), intent(in) :: adjmat2
   integer :: i, nkk, nkl, nll, nlk, delta

   nkk = 0
   nkl = 0

   do i = 1, nadjs1(k)
      if (adjlists1(i, k) /= l) then
         if (adjmat2(atomperm(k), atomperm(adjlists1(i, k)))) nkk = nkk + 1
         if (adjmat2(atomperm(l), atomperm(adjlists1(i, k)))) nkl = nkl + 1
      end if
   end do

   nll = 0
   nlk = 0

   do i = 1, nadjs1(l)
      if (adjlists1(i, l) /= k) then
         if (adjmat2(atomperm(l), atomperm(adjlists1(i, l)))) nll = nll + 1
         if (adjmat2(atomperm(k), atomperm(adjlists1(i, l)))) nlk = nlk + 1
      end if
   end do

!        dkk = nadjs1(k) + nadjs2(atomperm(k)) - 2*nkk
!        dll = nadjs1(l) + nadjs2(atomperm(l)) - 2*nll
!        dkl = nadjs1(k) + nadjs2(atomperm(l)) - 2*nkl
!        dlk = nadjs1(l) + nadjs2(atomperm(k)) - 2*nlk
!        delta = dkl + dlk - dkk - dll

   ! Notice that dkl + dlk - dkk - dll == 2*(nkk + nll - nkl - nlk)
   delta = 2*(nkk + nll - nkl - nlk)

end function

end module
