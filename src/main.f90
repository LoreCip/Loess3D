program loess3d

    use iso_fortran_env, only: RK => real64
    
    use OMP_LIB

    use TimerModule
    use utils
    use ioH5
    use mathFunc
    
    implicit none (type, external)

    ! Physical parameters
    real(RK) :: f, w, frac
    real(RK), dimension(:), allocatable :: xtmp, ytmp, ztmp, x, y, z, O, LO, LW
    real(RK), dimension(:,:,:), allocatable :: Oout, Wout, Xin, Yin, Zin, Oin
    integer  :: totL, npoints, n, l, m
    
    ! Computational parameters
    type(TimerClass) :: timer
    logical :: got
    integer(hid_t) :: file_id
    integer  :: ii, j, d, Nth, status, degree

    integer :: num_args
    character(len=4096), dimension(2) :: args

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

    call open_hdf5file(args(1), file_id, status)
    
    call read_from_hdf5(n, 'n', file_id, got, status)
    call read_from_hdf5(m, 'm', file_id, got, status)
    call read_from_hdf5(l, 'l', file_id, got, status)
    call read_from_hdf5(Nth, 'Nth', file_id, got, status)
    call read_from_hdf5(degree, 'degree', file_id, got, status)
    call read_from_hdf5(frac, 'frac', file_id, got, status)
    
    totL = n*m*l
    allocate(xtmp(n), ytmp(m), ztmp(l), Oin(n,m,l))
             
    call read_from_hdf5(xtmp, 'Yq', file_id, got, status)
    call read_from_hdf5(ytmp, 'logtemp', file_id, got, status)
    call read_from_hdf5(ztmp, 'lognb', file_id, got, status)
    call read_from_hdf5(Oin, 'LogEntropy', file_id, got, status)
    call close_hdf5file(file_id, status)

    allocate(Xin(n,m,l), Yin(n,m,l), Zin(n,m,l))
    
    do ii = 1, n
        Xin(ii, :, :) = xtmp(ii)
    end do
    do ii = 1, m
        Yin(:, ii, :) = ytmp(ii)
    end do
    do ii = 1, l
        Zin(:, :, ii) = ztmp(ii)
    end do

    deallocate(xtmp, ytmp, ztmp)
    allocate(x(totL), y(totL), z(totL), O(totL), LO(totL), LW(totL))

    call Cflatten(Xin, x)
    call Cflatten(Yin, y)
    call Cflatten(Zin, z)
    call Cflatten(Oin, O)

    deallocate(Xin, Yin, Zin, Oin)

    npoints = int(ceiling(frac*real(totL)))
    d = 2*(degree**2 + 1)

    call timer%initialize()
    call timer%start()

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii, f, w) NUM_THREADS(Nth)
    do ii = 1, totL

        call compute_loess(ii, totL, npoints, d, x, y, z, O, f, w)
        LO(ii) = f
        LW(ii) = w

    end do
    !$OMP END PARALLEL DO

    call timer%stop()    
    call timer%printTime()
    
    allocate(Oout(n,m,l), Wout(n,m,l))

    call Cpack_3d(LO, Oout, n, m, l, Nth)
    call Cpack_3d(LW, Wout, n, m, l, Nth)

    deallocate(x, y, z, O, LO)

!!!! OUTPUT

    call open_hdf5file(args(1), file_id, status)
    call write_to_hdf5(Oout, "S_LogEntropy", file_id, status)
    call write_to_hdf5(Wout, "W_LogEntropy", file_id, status)
    call close_hdf5file(file_id, status)

!!! END OUTPUT

    deallocate(Oout, Wout)

    stop
end program loess3d

