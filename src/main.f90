program loess3d

    use iso_fortran_env, only: RK => real64

    use OMP_LIB
    
    use ioFunc
    use mathFunc
    
    implicit none (type, external)

    ! Physical parameters
    real(RK) :: f, w
    real(RK), dimension(:), allocatable :: x, y, z, O, Oout, Wout
    integer  :: totL, npoints, n, l, m
    
    real(RK), parameter :: frac = 0.1
    integer, parameter :: degree = 1

    ! Computational parameters
    integer  :: j, d, Nth

    integer :: num_args
    character(len=4096), dimension(2) :: args

    real(RK) :: systemtime
    integer  :: iTimes1, iTimes2, rate

    num_args = command_argument_count()
    if (num_args .gt. 2) STOP "Only two args are contemplated, the path of the data file and the path for the output!"
    do j = 1, 2
        call get_command_argument(j,args(j))
        if (args(j) .eq. '') then
            if (j.eq.1) then
                args(j) = 'data.dat'
            else if (j .eq. 2) then
                args(j) = 'output.dat'
            end if
        end if
    end do

    call readParams(args(1), totL, n, l, m, Nth)
    allocate(x(totL), y(totL), z(totL), O(totL))
    call readData(args(1), totL, x, y, z, O)

    npoints = int(ceiling(frac*totL))
    d = 2*(degree**2 + 1)

    allocate(Oout(totL), Wout(totL))

    CALL system_clock(count_rate=rate)
    call SYSTEM_CLOCK(iTimes1)

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j, f, w) NUM_THREADS(Nth)
    do j = 1, totL

        call compute_loess(j, totL, npoints, d, x, y, z, O, f, w)
        Oout(j) = f
        Wout(j) = w

    end do
    !$OMP END PARALLEL DO

    call saveOutput(O, Oout, Wout, x, y, z, n, l, m, totL, args(2))
    
    deallocate(x, y, z, O, Oout, Wout)

    call SYSTEM_CLOCK(iTimes2)
    systemtime = real(iTimes2-iTimes1)/real(rate)
    write(*, *) "Total system runtime:", systemtime, " seconds."
    write(*, *) "System runtime for iteration:", systemtime/totL, " seconds."
    
    stop
end program loess3d