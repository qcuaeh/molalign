module assignment
use globals
use hungarian
implicit none
contains

subroutine assignatoms(natom, weights, mol0, mol1, nblock, blocksize, bias, atomap)
! Purpose: Find best correspondence between points shapes

! nblock: Number of block atoms
! blocksize: Number of atoms in each block
! atomap: Map between correspondent points in the adjmat
! offset: First element of current block

    integer, intent(in) :: natom, nblock
    integer, dimension(:), intent(in) :: blocksize
    real, dimension(:), intent(in) :: weights
    real, dimension(:, :), intent(in) :: mol0
    real, dimension(:, :), intent(in) :: bias
    real, dimension(:, :), intent(inout) :: mol1
    integer, dimension(:), intent(out) :: atomap

    integer h, i, j, offset
    integer, dimension(natom) :: blockmap
    real, dimension(natom, natom) :: costs
    real blocksum

! Fill distance matrix for each block

    offset = 0

    do h = 1, nblock

        do i = offset + 1, offset + blocksize(h)
            do j = offset + 1, offset + blocksize(h)
                costs(i - offset, j - offset) = weights(i)*(sum((mol1(:, j) - mol0(:, i))**2) + scaling**2*bias(i, j))
            end do
        end do

! Find correspondence between points in the current block

        call assndx(1, costs, blocksize(h), blocksize(h), blockmap, blocksum)
        atomap(offset+1:offset+blocksize(h)) = blockmap(:blocksize(h)) + offset

        offset = offset + blocksize(h)

    end do

end subroutine

function biasingdist(natom, weights, bias, mapping)
    integer, intent(in) :: natom
    real, dimension(:), intent(in) :: weights
    integer, dimension(:), intent(in) :: mapping
    real, dimension(:, :), intent(in) :: bias
    real biasingdist
    integer i

    biasingdist = 0.
    do i = 1, natom
        biasingdist = biasingdist + weights(i)*scaling**2*bias(i, mapping(i))
    end do

end function

end module
