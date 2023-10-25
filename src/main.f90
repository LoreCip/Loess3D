program loess3d

    use hdf5
    use iso_fortran_env, only: RK => real64

    use OMP_LIB
    
    use ioFunc
    use ioH5
    use mathFunc
    
    implicit none (type, external)

    ! Physical parameters
    real(RK) :: f, w
    real(RK), dimension(:), allocatable :: x, y, z, O, Oout, Wout
    integer  :: totL, npoints, n, l, m
    
    real(RK), parameter :: frac = 0.1
    integer, parameter :: degree = 1

    ! Computational parameters
    integer(hid_t) :: file_id
    integer  :: j, d, Nth, status

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
                args(j) = 'output.h5'
            end if
        end if
    end do

    call readParams(args(1), totL, n, l, m, Nth)
    allocate(x(totL), y(totL), z(totL), O(totL))
    call readData(args(1), totL, x, y, z, O)

    npoints = int(ceiling(frac*totL))
    d = 2*(degree**2 + 1)

    allocate(Oout(totL), Wout(totL))

    call create_hdf5_file(args(2), n, m, l, totL, x, y, z, O, file_id)

    CALL system_clock(count_rate=rate)
    call SYSTEM_CLOCK(iTimes1)

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j, f, w) NUM_THREADS(Nth)
    do j = 1, totL

        call compute_loess(j, totL, npoints, d, x, y, z, O, f, w)
        Oout(j) = f
        Wout(j) = w

    end do
    !$OMP END PARALLEL DO

    call SYSTEM_CLOCK(iTimes2)
    systemtime = real(iTimes2-iTimes1)/real(rate)
    write(*, *) "Total system runtime:", systemtime, " seconds."
    write(*, *) "System runtime for iteration:", systemtime/totL, " seconds."
    
    call write_real_1d_array(Oout, "f_out", file_id, status)
    call close_hdf5_file(file_id)

    deallocate(x, y, z, O, Oout, Wout)

    stop
end program loess3d