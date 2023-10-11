program loess3d

    use iso_fortran_env, only: RK => real64

    use OMP_LIB
    
    use ioFunc
    use mathFunc
    
    implicit none (type, external)

    ! Physical parameters

    real(RK), dimension(:), allocatable :: x, y, z, O, Oout, Wout, Ofit
    integer  :: totL, npoints, n, l, m
    
    real(RK), parameter :: frac = 0.1
    integer, parameter :: degree = 1

    ! Computational parameters

    real(RK), dimension(:), allocatable :: dist, dist_weights, aerr, uu, biweights, tot_weights, coeff
    real(RK) :: xj, yj, zj, mad 
    logical, dimension(:), allocatable :: bad, bad_old
    integer, dimension(:), allocatable :: inds, Tinds
    integer  :: j, p, d

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

    call readParams(args(1), totL, n, l, m)
    allocate(x(totL), y(totL), z(totL), O(totL))
    call readData(args(1), totL, x, y, z, O)

    npoints = int(ceiling(frac*totL))
    d = 2*(degree**2 + 1)

    allocate(Oout(totL), Wout(totL), dist(totL), Tinds(totL))
    allocate(dist_weights(npoints), Ofit(npoints), aerr(npoints),   \
             uu(npoints), biweights(npoints), tot_weights(npoints), \
             bad(npoints), bad_old(npoints), inds(npoints))
    allocate(coeff(d))

    CALL system_clock(count_rate=rate)
    call SYSTEM_CLOCK(iTimes1)

    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(totL, d, npoints, x, y, z, O, Oout, Wout)
    do j = 1, size(x)

        xj = x(j)
        yj = y(j)
        zj = z(j)

        dist = sqrt( (x(:) - xj)**2_RK + (y(:) - yj)**2_RK + (z(:) - zj)**2_RK )
        call merge_argsort(totL, dist, Tinds)
        inds = Tinds(:npoints)  
        dist_weights = (1_RK - (dist(inds) / dist(inds(npoints)))**3_RK )**3_RK
        
        call comp_Ofit(npoints, x(inds), y(inds), z(inds), O(inds), d, dist_weights, Ofit, coeff)

        bad(:) = .false.
        do p = 1, 10
            aerr(:) = abs(Ofit(:) - O(inds))
            call Median(npoints, aerr, mad)
            uu = ( aerr(:) / (6_RK*mad) )**2_RK
            
            uu(:) = max(0.0_RK, min(uu(:), 1.0_RK))

            biweights(:) = (1 - uu(:))**2
            tot_weights(:) = dist_weights(:)*biweights(:)

            call comp_Ofit(npoints, x(inds), y(inds), z(inds), O(inds), d, tot_weights, Ofit, coeff)

            bad_old(:) = bad(:)
            bad(:) = (biweights(:).lt.0.34_RK)

            if (all(bad.eqv.bad_old)) then
                exit
            end if
        end do

        Oout(j) = Ofit(1)
        Wout(j) = biweights(1)
    end do
    !$OMP END PARALLEL DO

    call saveOutput(O, Oout, Wout, x, y, z, n, l, m, totL, args(2))
    
    deallocate(dist_weights, Ofit, aerr, uu, biweights,    \
               tot_weights, coeff, bad, bad_old, dist, inds)
    deallocate(x, y, z, O, Oout, Wout, Tinds)

    call SYSTEM_CLOCK(iTimes2)
    systemtime = real(iTimes2-iTimes1)/real(rate)
    write(*, *) "Total system runtime:", systemtime, " seconds."
    write(*, *) "System runtime for iteration:", systemtime/totL, " seconds."
    
    stop
end program loess3d