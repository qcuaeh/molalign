module adjacency
use parameters
use sorting
use chemdata

implicit none

contains

function adjacencydiff(atomperm, adjmat1, adjmat2) result(diff)
! Purpose: Check if two graphs are equal.
! Return the number of differences between graphs.
   integer, dimension(:), intent(in) :: atomperm
   logical, dimension(:, :), intent(in) :: adjmat1, adjmat2
   integer :: diff

   integer :: i, j

   diff = 0

! Check differences element by element

   do i = 1, size(atomperm)
      do j = i + 1, size(atomperm)
         if (adjmat1(i, j) .neqv. adjmat2(atomperm(i), atomperm(j))) then
!            print *, i, j, adjmat1(i, j), adjmat2(atomperm(i), atomperm(j))
            diff = diff + 1
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
