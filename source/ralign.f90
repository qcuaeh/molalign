program ralign

use options
use strutils
use messages
use decoding
use rotation
use readwrite
!use superposition

implicit none

integer natom0, natom1
integer nrecord, maxrecord
character(512) title0, title1
real(wp), dimension(3) :: center0, center1
integer, dimension(:, :), allocatable :: atomaplist
real(wp), dimension(:, :, :), allocatable :: rotmatlist
real(wp), dimension(:, :), allocatable :: atoms0, atoms1
character(32), dimension(:), allocatable :: labels0, labels1

integer i
character(32) arg

! Set default options

live = .false.
remap = .false.
biased = .false.
iterative = .false.
trialing = .false.
matching = .false.
testing = .false.
maxrecord = 9
scale = 1000.0_wp
weighter = 'none'
outformat = 'xyz'

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
    case ('-trials')
        remap = .true.
        trialing = .true.
        call readoptarg(arg, maxtrial)
    case ('-matches')
        remap = .true.
        matching = .true.
        call readoptarg(arg, maxmatch)
    case ('-scale')
        call readoptarg(arg, scale)
    case ('-weight')
        call readoptarg(arg, weighter)
    case ('-records')
        call readoptarg(arg, maxrecord)
    case ('-out')
        call readoptarg(arg, outformat)
    case default
        call error('Invalid option: '//trim(arg))
    end select

end do

! Read coordinates

call readxyzfile(natom0, title0, labels0, atoms0)
call readxyzfile(natom1, title1, labels1, atoms1)

! Check number of atoms

if (natom0 /= natom1) then
    call error('The molecules do not have the same number of atoms!')
end if

! Allocate records

allocate(rotmatlist(3, 3, maxrecord))
allocate(atomaplist(natom0, maxrecord))

! Superpose atoms

call superpose(natom0, atoms0, atoms1, labels0, labels1, center0, center1, maxrecord, nrecord, atomaplist, rotmatlist)

! Write aligned coordinates

do i = 1, nrecord
    open (file_unit, file='aligned_'//str(i)//'.'//trim(outformat), action='write', status='replace')
    call writexyzfile(file_unit, natom0, [(i, i=1, natom0)], title0, labels0, atoms0)
    call writexyzfile(file_unit, natom1, atomaplist(:, i), title1, labels1, rotated(natom1, rotmatlist(:, :, i), atoms1))
    close (file_unit)
end do

end program
