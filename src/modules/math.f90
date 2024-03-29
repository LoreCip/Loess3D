module mathFunc
    use sort_interface
    use omp_lib
    use iso_fortran_env, only: RK => real64

    implicit none (type, external)
    
    public :: compute_loess
    private :: comp_coeff, comp_Ofit, comp_a

contains

    subroutine compute_loess(j, totL, npoints, d, x, y, z, O, fit, w)

        integer, intent(in) :: j, totL, npoints, d
        real(RK), dimension(totL), intent(in) :: x, y, z, O
        real(RK), intent(out) :: fit, w

        real(RK), dimension(totL)    :: dist
        real(RK), dimension(npoints) :: dist_weights, aerr, uu, biweights, tot_weights, Ofit
        real(RK) :: xj, yj, zj, mad
        logical, dimension(npoints) :: bad, bad_old
        integer, dimension(npoints) :: inds
        integer  :: p

        xj = x(j)
        yj = y(j)
        zj = z(j)

        dist = sqrt( (x(:) - xj)**2_RK + (y(:) - yj)**2_RK + (z(:) - zj)**2_RK )
        call pargsort(dist, inds, npoints)
        dist_weights = (1_RK - (dist(inds) / dist(inds(npoints)))**3_RK )**3_RK

        call comp_Ofit(npoints, x(inds), y(inds), z(inds), O(inds), d, dist_weights, Ofit)

        bad(:) = .false.
        do p = 1, 10
            aerr(:) = abs(Ofit(:) - O(inds))
            mad = median(aerr)
            uu = ( aerr(:) / (6_RK*mad) )**2_RK

            uu(:) = max(0.0_RK, min(uu(:), 1.0_RK))

            biweights(:) = (1_RK - uu(:))**2_RK
            tot_weights(:) = dist_weights(:)*biweights(:)

            call comp_Ofit(npoints, x(inds), y(inds), z(inds), O(inds), d, tot_weights, Ofit)

            bad_old(:) = bad(:)
            bad(:) = (biweights(:).lt.0.34_RK)

            if (all(bad.eqv.bad_old)) then
                exit
            end if
        end do

        fit = Ofit(1)
        w = biweights(1)

    end subroutine compute_loess

    subroutine comp_Ofit(npoints, x, y, z, O, d, weights, Ofit)
        
        integer, intent(in) :: npoints, d
        real(RK), dimension(npoints), intent(in) :: x, y, z, O, weights
        real(RK), dimension(npoints), intent(out) :: Ofit
        
        real(RK), dimension(d) :: coeff
        real(RK), dimension(npoints, d) :: a

        call comp_a(npoints, x, y, z, d, a)
        call comp_coeff(npoints, d, O, a, weights, coeff)
        Ofit = matmul(a, coeff)
        
    end subroutine comp_Ofit

    subroutine comp_a(npoints, x, y, z, d, a)

        integer, intent(in) :: npoints, d
        real(RK), dimension(npoints), intent(in) :: x, y, z
        real(RK), dimension(npoints, d), intent(out) :: a
        
        a(:, 1) = 1
        a(:, 2) = x
        a(:, 3) = y
        a(:, 4) = z
        if (d.eq.10) then
            a(:, 5) = x*y
            a(:, 6) = x*z
            a(:, 7) = y*z
            a(:, 8) = x**2
            a(:, 9) = y**2
            a(:, 10) = z**2
        end if

    end subroutine comp_a

    subroutine comp_coeff(npoints, d, O, a, weights, coeff)

        integer, intent(in) :: npoints, d
        real(RK), dimension(npoints), intent(in) :: O, weights
        real(RK), dimension(npoints, d), intent(in) :: a
        real(RK), dimension(d), intent(out) :: coeff

        real(RK), dimension(npoints) :: sqw
        real(RK), dimension(npoints, d) :: LHS

        integer :: i

        ! LAPACK
        real(RK), dimension(npoints) :: b
        real(RK), dimension(2*d) :: work
        integer :: lwork, info
        external :: dgels

        sqw = sqrt(weights)

        lwork = 132
        do i = 1, npoints
            LHS(i,:) = a(i,:)*sqw(i)
        end do
        b(:) = O(:)*sqw(:)
        
        call dgels('N', npoints, d, 1, LHS, npoints, b, npoints, work, lwork, info)
        coeff = b(:d)

    end subroutine comp_coeff

end module mathFunc