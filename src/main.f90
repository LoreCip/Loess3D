program loess3d

    use iso_fortran_env, only: RK => real64
    
    use OMP_LIB

    use ioFunc
    use ioH5
    use mathFunc
    
    implicit none (type, external)

    ! Physical parameters
    real(RK) :: f, w
    real(RK), dimension(:), allocatable :: x, y, z, O, LO, LW
    real(RK), dimension(:,:,:), allocatable :: Oout, Wout, Oin, Xout, Yout, Zout
    integer  :: totL, npoints, n, l, m
    
    real(RK), parameter :: frac = 0.1
    integer, parameter :: degree = 1

    ! Computational parameters
    integer(hid_t) :: file_id
    integer  :: ii, i, j, k, d, Nth, status

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
    allocate(x(totL), y(totL), z(totL), O(totL), LO(totL), LW(totL))
    call readData(args(1), totL, x, y, z, O)

    npoints = int(ceiling(frac*totL))
    d = 2*(degree**2 + 1)

    CALL system_clock(count_rate=rate)
    call SYSTEM_CLOCK(iTimes1)

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j, f, w) NUM_THREADS(Nth)
    do ii = 1, totL

        call compute_loess(ii, totL, npoints, d, x, y, z, O, f, w)
        LO(ii) = f
        LW(ii) = w

    end do
    !$OMP END PARALLEL DO

    call SYSTEM_CLOCK(iTimes2)
    systemtime = real(iTimes2-iTimes1)/real(rate)
    write(*, *) "Total system runtime:", systemtime, " seconds."
    write(*, *) "System runtime for iteration:", systemtime/totL, " seconds."
    
    allocate(Oout(n,m,l), Wout(n,m,l), Oin(n,m,l), Xout(n,m,l), Yout(n,m,l), Zout(n,m,l))

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii, i, j, k) NUM_THREADS(Nth)
    do ii = 1, totL
        k = mod(ii - 1, l) + 1
        j = mod((ii - 1 - k + 1) / l, m) + 1
        i = mod(( (ii-k)/l - (j-1)) / m , n) + 1

        Xout(i,j,k) = x(ii)
        Yout(i,j,k) = y(ii)
        Zout(i,j,k) = z(ii)
        Oin(i,j,k)  = O(ii)
        Oout(i,j,k) = LO(ii)
        Wout(i,j,k) = LW(ii)
    end do
    !$OMP END PARALLEL DO

    deallocate(x, y, z, O, LO)

!!!! OUTPUT
    call create_hdf5file(args(2), file_id, status)
    call write_to_hdf5(n, "n", file_id, status)
    call write_to_hdf5(m, "m", file_id, status)
    call write_to_hdf5(l, "l", file_id, status)
    call write_to_hdf5(Xout, "Yq",      file_id, status)
    call write_to_hdf5(Yout, "logtemp", file_id, status)
    call write_to_hdf5(Zout, "lognb",   file_id, status)
    call write_to_hdf5(Oin, "f_in", file_id, status)
    call write_to_hdf5(Oout, "f_out", file_id, status)
    call close_hdf5file(file_id, status)
!!! END OUTPUT

    deallocate(Oout, Wout, Xout, Yout, Zout, Oin)

    stop
end program loess3d