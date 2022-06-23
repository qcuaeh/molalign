program main

    use globals
    use decoding
    use readwrite
    use utilities
    use rotation

    implicit none

    integer mapcount
    integer, dimension(:, :), allocatable :: atomaplist
    real, dimension(:, :), allocatable :: rotquatlist
    integer, dimension(:, :), allocatable :: bonds0, bonds1
    real, dimension(:, :), allocatable :: mol0, mol1
    character(32), dimension(:), allocatable :: label0, label1
    character(512) title0, title1

    integer i
    character(32) arg

    ! Set default options 

    biased = .false.
    iterative = .false.
    testable = .false.
    live = .false.

    maxcoord = 0
    maxcount = -1
    maxrecord = 9

    scaling = 1000.
    tolerance = 0.1

    weighting = 'none'

    stop_test => trial_stop_test

    ! Get user options 

    call initarg()

    do while (getarg(arg))

        select case (arg)
        case ('-iter')
            iterative = .true.
        case ('-converge')
            stop_test => match_stop_test
        case ('-bias')
            biased = .true.
        case ('-test')
            testable = .true.
        case ('-live')
            live = .true.
        case ('-scale')
            call readoptarg(arg, scaling)
        case ('-tol')
            call readoptarg(arg, tolerance)
        case ('-weight')
            call readoptarg(arg, weighting)
        case ('-out')
            call readoptarg(arg, output_format)
        case ('-n')
            call readoptarg(arg, maxrecord)
        case default
            call readarg(arg, maxcount)
        end select

    end do

    ! Read coordinates

    call readmol('xyz', title0, label0, mol0, bonds0)
    call readmol('xyz', title1, label1, mol1, bonds1)

    ! Align atoms

    call ralign(mol0, mol1, label0, label1, bonds0, bonds1, mapcount, atomaplist, rotquatlist)

    ! Write aligned coordinates

    do i = 1, mapcount
        open (file_unit, file='aligned_'//str(i)//'.'//trim(output_format), action='write', status='replace')
        call writemol(file_unit, [(i, i=1, size(label0))], title0, label0, mol0, bonds0)
        call writemol(file_unit, atomaplist(:, i), title1, label1, &
            rotated(size(mol1, dim=1), rotquatlist(:, i), mol1), &
            bonds1)
        close (file_unit)
    end do

contains

    logical function trial_stop_test(trials, matches) result(stop_test)
        integer, intent(in) :: trials, matches
        stop_test = trials < maxcount
    end function

    logical function match_stop_test(trials, matches) result(stop_test)
        integer, intent(in) :: trials, matches
        stop_test = matches < maxcount
    end function

end program
