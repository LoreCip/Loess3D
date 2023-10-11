module mathFunc

    use iso_fortran_env, only: RK => real64
    implicit none (type, external)
    
contains

    subroutine comp_Ofit(npoints, x, y, z, O, d, weights, Ofit, coeff)
        
        integer, intent(in) :: npoints, d
        real(RK), dimension(npoints), intent(in) :: x, y, z, O, weights
        real(RK), dimension(npoints), intent(out) :: Ofit
        real(RK), dimension(d), intent(out) :: coeff

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

    subroutine Median(npoints, arr, med)
        
        integer, intent(in)   :: npoints
        real(RK), dimension(npoints), intent(in) :: arr
        real(RK), intent(out) :: med

        integer, dimension(npoints) :: inds
        integer :: n

        call merge_argsort(npoints, arr, inds)
        n = npoints / 2        
        if (mod(npoints,2) == 0) then           ! compute the median
            med = (arr(inds(n)) + arr(inds(n+1))) / 2.0_RK
        else
            med = arr(inds(n+1))
        end if

    end subroutine Median

    subroutine merge_argsort(length, array, indices)

        integer, intent(in) :: length
        real(RK), dimension(length), intent(in) :: array
        integer, dimension(length), intent(out) :: indices
      
        integer, dimension(length) :: il

        integer :: stepsize
        integer :: i, j, left, k, ksize
        
        do i = 1, length
            indices(i) = i
        end do
      
        stepsize = 1
        do while (stepsize.lt.length)
            do left= 1, length-stepsize, stepsize*2
                i = left
                j = left + stepsize
                ksize = min(stepsize*2, length-left+1)
                k=1
          
                do while ( (i.lt.left+stepsize) .and. (j.lt.left+ksize) )
                    if ( array(indices(i)).gt.array(indices(j)) ) then
                        il(k)=indices(i)
                        i=i+1
                        k=k+1
                    else
                        il(k)=indices(j)
                        j=j+1
                        k=k+1
                    endif
                enddo
          
                if ( i.lt.left+stepsize ) then
                    ! fill up remaining from left
                    il(k:ksize) = indices(i:left+stepsize-1)
                else
                    ! fill up remaining from right
                    il(k:ksize) = indices(j:left+ksize-1)
                endif
                indices(left:left+ksize-1) = il(1:ksize)
            end do
            stepsize=stepsize*2
        end do

        indices = indices(length:1:-1)

    end subroutine    

end module mathFunc