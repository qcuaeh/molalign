program ralign

use iso_fortran_env, only: input_unit
use iso_fortran_env, only: output_unit

use options
use chemdata
use strutils
use optparse
use chemutils
use readwrite
use translation
use rotation
use alignment
!use library

implicit none

integer i, n, u
integer natom0, natom1
integer nrecord, maxrecord
real travec(3), rotmat(3, 3)
character(title_len) title0, title1
integer, allocatable :: mapcount(:)
integer, allocatable :: maplist(:, :)
integer, dimension(:), allocatable :: znums0, znums1, types0, types1
real, dimension(:, :), allocatable :: coords0, coords1
real, dimension(:), allocatable :: weights0, weights1
character(label_len), dimension(:), allocatable :: labels0, labels1
character(arg_len) arg, path, weighting
real, dimension(nelem) :: property
integer file_unit(2)
logical stdin

procedure(writeabstractfile), pointer :: writefile => null()

! Set default options

live = .false.
stdin = .false.
biased = .false.
iterative = .false.
bounded = .false.
converged = .false.
testing = .false.
maxrecord = 1
lenscale = 1000.0
weighting = 'none'
formatout = 'xyz'

! Get user options

call initarg()

do while (getarg(arg))

    select case (arg)
    case ('-live')
        live = .true.
    case ('-test')
        testing = .true.
    case ('-iter')
        iterative = .true.
    case ('-bias')
        biased = .true.
        call readoptarg(arg, tolerance)
    case ('-trial')
        bounded = .true.
        call readoptarg(arg, maxtrial)
    case ('-count')
        converged = .true.
        call readoptarg(arg, mincount)
    case ('-scale')
        call readoptarg(arg, lenscale)
    case ('-weight')
        call readoptarg(arg, weighting)
    case ('-print')
        call readoptarg(arg, maxrecord)
    case ('-out')
        call readoptarg(arg, formatout)
    case ('-stdin')
        stdin = .true.
    case default
        call openfile(arg, file_unit)
    end select

end do

if (stdin) then
    file_unit(1) = input_unit
    file_unit(2) = input_unit
else
    if (ifile == 0) then
        write (error_unit, '(a)') 'Error: No file paths were specified'
        stop
    else if (ifile == 1) then
        file_unit(2) = file_unit(1)
    end if
end if

! Read coordinates

call readxyzfile(file_unit(1), natom0, title0, labels0, coords0)
call readxyzfile(file_unit(2), natom1, title1, labels1, coords1)

! Allocate arrays

allocate(znums0(natom0), znums1(natom1))
allocate(types0(natom0), types1(natom1))
allocate(weights0(natom0), weights1(natom1))
allocate(maplist(natom0, maxrecord), mapcount(maxrecord))

! Get atomic numbers and types

do i = 1, natom0
    call getznum(labels0(i), znums0(i), types0(i))
end do

do i = 1, natom1
    call getznum(labels1(i), znums1(i), types1(i))
end do

! Select output format

select case (formatout)
case ('xyz')
    writefile => writexyzfile
case ('mol2')
    writefile => writemol2file
case default
    write (error_unit, '(a,x,a)') 'Invalid format:', trim(formatout)
    stop
end select

! Select weighting property

select case (weighting)
case ('none')
    property = [(1.0, i=1, nelem)]
case ('mass')
    property = stdmatom
case default
    write (error_unit, '(a,x,a)') 'Invalid weighting option:', trim(weighting)
    stop
end select

! Get normalized weights

weights0 = property(znums0)/sum(property(znums0))
weights1 = property(znums1)/sum(property(znums1))

! Superpose atoms

if (bounded .or. converged) then

    call remap(natom0, natom1, znums0, znums1, types0, types1, weights0, weights1, &
        coords0, coords1, maxrecord, nrecord, maplist, mapcount)

    ! Write aligned coordinates

    do i = 1, nrecord

        call align( &
            natom0, natom1, &
            znums0, znums1(maplist(:, i)), &
            types0, types1(maplist(:, i)), &
            weights0, weights1(maplist(:, i)), &
            coords0, coords1(:, maplist(:, i)), &
            travec, rotmat &
        )

        open(newunit=u, file='aligned_'//str(i)//'.'//trim(formatout), action='write', status='replace')
        call writefile(u, natom0, title0, znums0, coords0)
        call writefile(u, natom1, title1, znums1(maplist(:, i)), &
            translated(natom1, rotated(natom1, coords1(:, maplist(:, i)), rotmat), travec))
        close(u)

    end do

else

    call align(natom0, natom1, znums0, znums1, types0, types1, weights0, weights1, &
        coords0, coords1, travec, rotmat)

    call rotate(natom1, coords1, rotmat)
    call translate(natom1, coords1, travec)

    write (output_unit, '(a,x,f0.4,x,a)') 'RMSD:', &
        squaredist(natom0, weights0, coords0, coords1, identitymap(natom0)), &
        '(only alignment performed)'

    open(newunit=u, file='aligned.'//trim(formatout), action='write', status='replace')
    call writefile(u, natom0, title0, znums0, coords0)
    call writefile(u, natom1, title1, znums1, coords1) 
    close(u)

end if

end program
