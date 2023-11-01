program loess3d

    use iso_fortran_env, only: RK => real64
    
#ifdef _OPENMP
    use omp_lib
#endif
#ifdef USE_MPI
    use mpi_f08
#endif

    use TimerModule
    use utils
    use ioH5
    use mathFunc
    
    implicit none (type, external)

    ! Physical parameters
    real(RK) :: f, w, frac
    real(RK), dimension(:), allocatable :: x, y, z, O, LO, LW
    real(RK), dimension(:,:,:), allocatable :: Oout, Wout, Xin, Yin, Zin, Oin
    integer  :: totL, npoints, n, l, m
    
    ! Computational parameters
    type(TimerClass) :: timer
    logical :: got
    integer(hid_t) :: file_id
    integer  :: ii, d, Nth, status, degree

    ! MPI parameters
    real(RK), dimension(:), allocatable :: mpi_LO, mpi_LW
    integer :: provided, mpi_len, mpi_rank, mpi_size

    integer :: num_args
    character(len=4096), dimension(1) :: args

    !! MPI STARTS

#ifdef USE_MPI
    call mpi_init_thread(MPI_THREAD_FUNNELED, provided, status)
    call mpi_comm_size(MPI_COMM_WORLD, mpi_size, status)
    call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, status)
#else
    mpi_rank = 0
    mpi_size = 1
#endif

    if (mpi_rank.eq.0) then
        num_args = command_argument_count()
        if (num_args .gt. 1) STOP "Only one arg is expected, the path of the data file."
        call get_command_argument(1,args(1))
        if (args(1) .eq. '') then
            args(1) = 'data.h5'
        end if

        call open_hdf5file(args(1), file_id, status)

        call read_from_hdf5(n, 'n', file_id, got, status)
        call read_from_hdf5(m, 'm', file_id, got, status)
        call read_from_hdf5(l, 'l', file_id, got, status)
        call read_from_hdf5(Nth, 'Nth', file_id, got, status)
        call read_from_hdf5(degree, 'degree', file_id, got, status)
        call read_from_hdf5(frac, 'frac', file_id, got, status)

        totL = n*m*l
        allocate(x(n), y(m), z(l), Oin(n,m,l))

        call read_from_hdf5(x, 'Yq', file_id, got, status)
        call read_from_hdf5(y, 'logtemp', file_id, got, status)
        call read_from_hdf5(z, 'lognb', file_id, got, status)
        call read_from_hdf5(Oin, 'LogEntropy', file_id, got, status)
        call close_hdf5file(file_id, status)

        allocate(Xin(n,m,l), Yin(n,m,l), Zin(n,m,l))

        do ii = 1, n
            Xin(ii, :, :) = x(ii)
        end do
        do ii = 1, m
            Yin(:, ii, :) = y(ii)
        end do
        do ii = 1, l
            Zin(:, :, ii) = z(ii)
        end do

        deallocate(x, y, z)
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
    end if

#ifdef USE_MPI
    call mpi_barrier(MPI_COMM_WORLD)

    ! Scalar broadcast
    call mpi_bcast(totL,    1, MPI_INTEGER, 0, MPI_COMM_WORLD)
    call mpi_bcast(d,       1, MPI_INTEGER, 0, MPI_COMM_WORLD)
    call mpi_bcast(npoints, 1, MPI_INTEGER, 0, MPI_COMM_WORLD)
    call mpi_bcast(Nth,     1, MPI_INTEGER, 0, MPI_COMM_WORLD)

    ! Array broadcast
    if (mpi_rank.ne.0) allocate(x(totL), y(totL), z(totL), O(totL), LO(totL), LW(totL))
    call mpi_bcast(x, totL, MPI_DOUBLE, 0, MPI_COMM_WORLD)
    call mpi_bcast(y, totL, MPI_DOUBLE, 0, MPI_COMM_WORLD)
    call mpi_bcast(z, totL, MPI_DOUBLE, 0, MPI_COMM_WORLD)
    call mpi_bcast(O, totL, MPI_DOUBLE, 0, MPI_COMM_WORLD)
#endif

    mpi_len = int(real(totL) / real(mpi_size))

#ifdef USE_MPI
    allocate(mpi_LO(mpi_len), mpi_LW(mpi_len))
#endif

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii, f, w) NUM_THREADS(Nth)
    do ii = 1, mpi_len

        call compute_loess(mpi_rank*mpi_len + ii, totL, npoints, d, x, y, z, O, f, w)
#ifdef USE_MPI
        mpi_LO(ii) = f
        mpi_LW(ii) = w
#else
        LO(ii) = f
        LW(ii) = w
#endif
    
    end do
    !$OMP END PARALLEL DO
    
#ifdef USE_MPI
    call mpi_gather(mpi_LO, mpi_len, MPI_DOUBLE, LO, mpi_len, MPI_DOUBLE, 0, MPI_COMM_WORLD)
    call mpi_gather(mpi_LW, mpi_len, MPI_DOUBLE, LW, mpi_len, MPI_DOUBLE, 0, MPI_COMM_WORLD)

    if (mpi_rank.ne.0) deallocate(x, y, z, O, mpi_LO, mpi_LW)
    call mpi_barrier(MPI_COMM_WORLD)
#endif

    if (mpi_rank.eq.0) then
        call timer%stop()    
        call timer%printTime()
    
        allocate(Oout(n,m,l), Wout(n,m,l))

        call Cpack_3d(LO, Oout, n, m, l, Nth)
        call Cpack_3d(LW, Wout, n, m, l, Nth)

        deallocate(x, y, z, O, LO, LW)

        !!! OUTPUT
        call open_hdf5file(args(1), file_id, status)
        call write_to_hdf5(Oout, "S_LogEntropy", file_id, status)
        call write_to_hdf5(Wout, "W_LogEntropy", file_id, status)
        call close_hdf5file(file_id, status)
        !!! END OUTPUT

        deallocate(Oout, Wout)
    end if

#ifdef USE_MPI
    call mpi_finalize(status)
#endif
    !! MPI ENDS

    stop
end program loess3d

