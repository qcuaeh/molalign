module biasing
use globals
use sorting
implicit none

contains

subroutine setadjbias(natom, nblock, blocksize, blocktype, mol0, ncoord0, adjlist0, &
mol1, ncoord1, adjlist1, weights, bias)
    integer, intent(in) :: natom, nblock
    integer, dimension(:), intent(in) :: ncoord0, ncoord1
    integer, dimension(:), intent(in) :: blocktype, blocksize
    integer, dimension(:, :), intent(in) :: adjlist0, adjlist1
    real, dimension(:), intent(in) :: weights
    real, dimension(:, :), intent(in) :: mol0, mol1
    real, dimension(:, :), intent(out) :: bias

    integer h, i, j, offset
    real d0(natom, natom), d1(natom, natom)

    do i = 1, natom
        offset = 0
        do h = 1, nblock
            do j = offset + 1, offset + blocksize(h)
                d0(j, i) = sqrt(sum((mol0(:, j) - mol0(:, i))**2))
            end do
            call sort(d0(:, i), offset + 1, offset + blocksize(h))
            offset = offset + blocksize(h)
        end do
    end do

    do i = 1, natom
        offset = 0
        do h = 1, nblock
            do j = offset + 1, offset + blocksize(h)
                d1(j, i) = sqrt(sum((mol1(:, j) - mol1(:, i))**2))
            end do
            call sort(d1(:, i), offset + 1, offset + blocksize(h))
            offset = offset + blocksize(h)
        end do
    end do

    offset = 0
    bias(:, :) = 0.

    if (biased) then
        do h = 1, nblock
            do i = offset + 1, offset + blocksize(h)
                do j = offset + 1, offset + blocksize(h)
                    if (any(abs(d1(:, j) - d0(:, i)) > tolerance)) then
                        bias(i, j) = 1.
                    end if
                end do
            end do
            offset = offset + blocksize(h)
        end do
    end if

end subroutine

end module
