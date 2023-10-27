module mathFunc
    use mrgrnk_mod
    use omp_lib
    use iso_fortran_env, only: RK => real64
    implicit none (type, external)
    
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
        integer, dimension(totL)    :: Tinds
        integer  :: p

        xj = x(j)
        yj = y(j)
        zj = z(j)

        dist = sqrt( (x(:) - xj)**2_RK + (y(:) - yj)**2_RK + (z(:) - zj)**2_RK )
        call mrgrnk(dist, Tinds)
        inds = Tinds(:npoints)
        dist_weights = (1_RK - (dist(inds) / dist(inds(npoints)))**3_RK )**3_RK

        call comp_Ofit(npoints, x(inds), y(inds), z(inds), O(inds), d, dist_weights, Ofit)

        bad(:) = .false.
        do p = 1, 10
            aerr(:) = abs(Ofit(:) - O(inds))
            call Median(npoints, aerr, mad)
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

    subroutine Median(npoints, arr, med)
        
        integer, intent(in)   :: npoints
        real(RK), dimension(npoints), intent(in) :: arr
        real(RK), intent(out) :: med

        integer, dimension(npoints) :: inds
        integer :: n

        call mrgrnk(arr, inds)
        n = npoints / 2        
        if (mod(npoints,2) == 0) then           ! compute the median
            med = (arr(inds(n)) + arr(inds(n+1))) / 2.0_RK
        else
            med = arr(inds(n+1))
        end if

    end subroutine Median  

end module mathFunc


Module mrgrnk_mod
    use iso_fortran_env, only: RK => real64
    public :: mrgrnk
    private :: R_mrgrnk, I_mrgrnk, D_mrgrnk
    interface mrgrnk
      module procedure D_mrgrnk, R_mrgrnk, I_mrgrnk
    end interface mrgrnk
    contains
    
    Subroutine D_mrgrnk (XDONT, IRNGT)
    ! __________________________________________________________
    !   MRGRNK = Merge-sort ranking of an array
    !   For performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.
    ! __________________________________________________________
    ! __________________________________________________________
          Real (RK), Dimension (:), Intent (In) :: XDONT
          Integer, Dimension (:), Intent (Out) :: IRNGT
    ! __________________________________________________________
          Real (RK) :: XVALA, XVALB
    !
          Integer, Dimension (SIZE(IRNGT)) :: JWRKT
          Integer :: LMTNA, LMTNC, IRNG1, IRNG2
          Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
    !
          NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
          Select Case (NVAL)
          Case (:0)
             Return
          Case (1)
             IRNGT (1) = 1
             Return
          Case Default
             Continue
          End Select
    !
    !  Fill-in the index array, creating ordered couples
    !
          Do IIND = 2, NVAL, 2
             If (XDONT(IIND-1) <= XDONT(IIND)) Then
                IRNGT (IIND-1) = IIND - 1
                IRNGT (IIND) = IIND
             Else
                IRNGT (IIND-1) = IIND
                IRNGT (IIND) = IIND - 1
             End If
          End Do
          If (Modulo(NVAL, 2) /= 0) Then
             IRNGT (NVAL) = NVAL
          End If
    !
    !  We will now have ordered subsets A - B - A - B - ...
    !  and merge A and B couples into     C   -   C   - ...
    !
          LMTNA = 2
          LMTNC = 4
    !
    !  First iteration. The length of the ordered subsets goes from 2 to 4
    !
          Do
             If (NVAL <= 2) Exit
    !
    !   Loop on merges of A and B into C
    !
             Do IWRKD = 0, NVAL - 1, 4
                If ((IWRKD+4) > NVAL) Then
                   If ((IWRKD+2) >= NVAL) Exit
    !
    !   1 2 3
    !
                   If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
    !
    !   1 3 2
    !
                   If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                      IRNG2 = IRNGT (IWRKD+2)
                      IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                      IRNGT (IWRKD+3) = IRNG2
    !
    !   3 1 2
    !
                   Else
                      IRNG1 = IRNGT (IWRKD+1)
                      IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                      IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                      IRNGT (IWRKD+2) = IRNG1
                   End If
                   Exit
                End If
    !
    !   1 2 3 4
    !
                If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
    !
    !   1 3 x x
    !
                If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                   IRNG2 = IRNGT (IWRKD+2)
                   IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                   If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
    !   1 3 2 4
                      IRNGT (IWRKD+3) = IRNG2
                   Else
    !   1 3 4 2
                      IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                      IRNGT (IWRKD+4) = IRNG2
                   End If
    !
    !   3 x x x
    !
                Else
                   IRNG1 = IRNGT (IWRKD+1)
                   IRNG2 = IRNGT (IWRKD+2)
                   IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                   If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                      IRNGT (IWRKD+2) = IRNG1
                      If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
    !   3 1 2 4
                         IRNGT (IWRKD+3) = IRNG2
                      Else
    !   3 1 4 2
                         IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                         IRNGT (IWRKD+4) = IRNG2
                      End If
                   Else
    !   3 4 1 2
                      IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                      IRNGT (IWRKD+3) = IRNG1
                      IRNGT (IWRKD+4) = IRNG2
                   End If
                End If
             End Do
    !
    !  The Cs become As and Bs
    !
             LMTNA = 4
             Exit
          End Do
    !
    !  Iteration loop. Each time, the length of the ordered subsets
    !  is doubled.
    !
          Do
             If (LMTNA >= NVAL) Exit
             IWRKF = 0
             LMTNC = 2 * LMTNC
    !
    !   Loop on merges of A and B into C
    !
             Do
                IWRK = IWRKF
                IWRKD = IWRKF + 1
                JINDA = IWRKF + LMTNA
                IWRKF = IWRKF + LMTNC
                If (IWRKF >= NVAL) Then
                   If (JINDA >= NVAL) Exit
                   IWRKF = NVAL
                End If
                IINDA = 1
                IINDB = JINDA + 1
    !
    !   Shortcut for the case when the max of A is smaller
    !   than the min of B. This line may be activated when the
    !   initial set is already close to sorted.
    !
             IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
    !
    !  One steps in the C subset, that we build in the final rank array
    !
    !  Make a copy of the rank array for the merge iteration
    !
                JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
    !
                XVALA = XDONT (JWRKT(IINDA))
                XVALB = XDONT (IRNGT(IINDB))
    !
                Do
                   IWRK = IWRK + 1
    !
    !  We still have unprocessed values in both A and B
    !
                   If (XVALA > XVALB) Then
                      IRNGT (IWRK) = IRNGT (IINDB)
                      IINDB = IINDB + 1
                      If (IINDB > IWRKF) Then
    !  Only A still with unprocessed values
                         IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                         Exit
                      End If
                      XVALB = XDONT (IRNGT(IINDB))
                   Else
                      IRNGT (IWRK) = JWRKT (IINDA)
                      IINDA = IINDA + 1
                      If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                      XVALA = XDONT (JWRKT(IINDA))
                   End If
    !
                End Do
             End Do
    !
    !  The Cs become As and Bs
    !
             LMTNA = 2 * LMTNA
          End Do
    !
          Return
    !
    End Subroutine D_mrgrnk
    
    Subroutine R_mrgrnk (XDONT, IRNGT)
    ! __________________________________________________________
    !   MRGRNK = Merge-sort ranking of an array
    !   For performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.
    ! __________________________________________________________
    ! _________________________________________________________
          Real, Dimension (:), Intent (In) :: XDONT
          Integer, Dimension (:), Intent (Out) :: IRNGT
    ! __________________________________________________________
          Real :: XVALA, XVALB
    !
          Integer, Dimension (SIZE(IRNGT)) :: JWRKT
          Integer :: LMTNA, LMTNC, IRNG1, IRNG2
          Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
    !
          NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
          Select Case (NVAL)
          Case (:0)
             Return
          Case (1)
             IRNGT (1) = 1
             Return
          Case Default
             Continue
          End Select
    !
    !  Fill-in the index array, creating ordered couples
    !
          Do IIND = 2, NVAL, 2
             If (XDONT(IIND-1) <= XDONT(IIND)) Then
                IRNGT (IIND-1) = IIND - 1
                IRNGT (IIND) = IIND
             Else
                IRNGT (IIND-1) = IIND
                IRNGT (IIND) = IIND - 1
             End If
          End Do
          If (Modulo(NVAL, 2) /= 0) Then
             IRNGT (NVAL) = NVAL
          End If
    !
    !  We will now have ordered subsets A - B - A - B - ...
    !  and merge A and B couples into     C   -   C   - ...
    !
          LMTNA = 2
          LMTNC = 4
    !
    !  First iteration. The length of the ordered subsets goes from 2 to 4
    !
          Do
             If (NVAL <= 2) Exit
    !
    !   Loop on merges of A and B into C
    !
             Do IWRKD = 0, NVAL - 1, 4
                If ((IWRKD+4) > NVAL) Then
                   If ((IWRKD+2) >= NVAL) Exit
    !
    !   1 2 3
    !
                   If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
    !
    !   1 3 2
    !
                   If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                      IRNG2 = IRNGT (IWRKD+2)
                      IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                      IRNGT (IWRKD+3) = IRNG2
    !
    !   3 1 2
    !
                   Else
                      IRNG1 = IRNGT (IWRKD+1)
                      IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                      IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                      IRNGT (IWRKD+2) = IRNG1
                   End If
                   Exit
                End If
    !
    !   1 2 3 4
    !
                If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
    !
    !   1 3 x x
    !
                If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                   IRNG2 = IRNGT (IWRKD+2)
                   IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                   If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
    !   1 3 2 4
                      IRNGT (IWRKD+3) = IRNG2
                   Else
    !   1 3 4 2
                      IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                      IRNGT (IWRKD+4) = IRNG2
                   End If
    !
    !   3 x x x
    !
                Else
                   IRNG1 = IRNGT (IWRKD+1)
                   IRNG2 = IRNGT (IWRKD+2)
                   IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                   If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                      IRNGT (IWRKD+2) = IRNG1
                      If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
    !   3 1 2 4
                         IRNGT (IWRKD+3) = IRNG2
                      Else
    !   3 1 4 2
                         IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                         IRNGT (IWRKD+4) = IRNG2
                      End If
                   Else
    !   3 4 1 2
                      IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                      IRNGT (IWRKD+3) = IRNG1
                      IRNGT (IWRKD+4) = IRNG2
                   End If
                End If
             End Do
    !
    !  The Cs become As and Bs
    !
             LMTNA = 4
             Exit
          End Do
    !
    !  Iteration loop. Each time, the length of the ordered subsets
    !  is doubled.
    !
          Do
             If (LMTNA >= NVAL) Exit
             IWRKF = 0
             LMTNC = 2 * LMTNC
    !
    !   Loop on merges of A and B into C
    !
             Do
                IWRK = IWRKF
                IWRKD = IWRKF + 1
                JINDA = IWRKF + LMTNA
                IWRKF = IWRKF + LMTNC
                If (IWRKF >= NVAL) Then
                   If (JINDA >= NVAL) Exit
                   IWRKF = NVAL
                End If
                IINDA = 1
                IINDB = JINDA + 1
    !
    !   Shortcut for the case when the max of A is smaller
    !   than the min of B. This line may be activated when the
    !   initial set is already close to sorted.
    !
             IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
    !
    !  One steps in the C subset, that we build in the final rank array
    !
    !  Make a copy of the rank array for the merge iteration
    !
                JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
    !
                XVALA = XDONT (JWRKT(IINDA))
                XVALB = XDONT (IRNGT(IINDB))
    !
                Do
                   IWRK = IWRK + 1
    !
    !  We still have unprocessed values in both A and B
    !
                   If (XVALA > XVALB) Then
                      IRNGT (IWRK) = IRNGT (IINDB)
                      IINDB = IINDB + 1
                      If (IINDB > IWRKF) Then
    !  Only A still with unprocessed values
                         IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                         Exit
                      End If
                      XVALB = XDONT (IRNGT(IINDB))
                   Else
                      IRNGT (IWRK) = JWRKT (IINDA)
                      IINDA = IINDA + 1
                      If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                      XVALA = XDONT (JWRKT(IINDA))
                   End If
    !
                End Do
             End Do
    !
    !  The Cs become As and Bs
    !
             LMTNA = 2 * LMTNA
          End Do
    !
          Return
    !
    End Subroutine R_mrgrnk

    Subroutine I_mrgrnk (XDONT, IRNGT)
    ! __________________________________________________________
    !   MRGRNK = Merge-sort ranking of an array
    !   For performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.
    ! __________________________________________________________
    ! __________________________________________________________
          Integer, Dimension (:), Intent (In)  :: XDONT
          Integer, Dimension (:), Intent (Out) :: IRNGT
    ! __________________________________________________________
          Integer :: XVALA, XVALB
    !
          Integer, Dimension (SIZE(IRNGT)) :: JWRKT
          Integer :: LMTNA, LMTNC, IRNG1, IRNG2
          Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
    !
          NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
          Select Case (NVAL)
          Case (:0)
             Return
          Case (1)
             IRNGT (1) = 1
             Return
          Case Default
             Continue
          End Select
    !
    !  Fill-in the index array, creating ordered couples
    !
          Do IIND = 2, NVAL, 2
             If (XDONT(IIND-1) <= XDONT(IIND)) Then
                IRNGT (IIND-1) = IIND - 1
                IRNGT (IIND) = IIND
             Else
                IRNGT (IIND-1) = IIND
                IRNGT (IIND) = IIND - 1
             End If
          End Do
          If (Modulo(NVAL, 2) /= 0) Then
             IRNGT (NVAL) = NVAL
          End If
    !
    !  We will now have ordered subsets A - B - A - B - ...
    !  and merge A and B couples into     C   -   C   - ...
    !
          LMTNA = 2
          LMTNC = 4
    !
    !  First iteration. The length of the ordered subsets goes from 2 to 4
    !
          Do
             If (NVAL <= 2) Exit
    !
    !   Loop on merges of A and B into C
    !
             Do IWRKD = 0, NVAL - 1, 4
                If ((IWRKD+4) > NVAL) Then
                   If ((IWRKD+2) >= NVAL) Exit
    !
    !   1 2 3
    !
                   If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
    !
    !   1 3 2
    !
                   If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                      IRNG2 = IRNGT (IWRKD+2)
                      IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                      IRNGT (IWRKD+3) = IRNG2
    !
    !   3 1 2
    !
                   Else
                      IRNG1 = IRNGT (IWRKD+1)
                      IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                      IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                      IRNGT (IWRKD+2) = IRNG1
                   End If
                   Exit
                End If
    !
    !   1 2 3 4
    !
                If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
    !
    !   1 3 x x
    !
                If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                   IRNG2 = IRNGT (IWRKD+2)
                   IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                   If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
    !   1 3 2 4
                      IRNGT (IWRKD+3) = IRNG2
                   Else
    !   1 3 4 2
                      IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                      IRNGT (IWRKD+4) = IRNG2
                   End If
    !
    !   3 x x x
    !
                Else
                   IRNG1 = IRNGT (IWRKD+1)
                   IRNG2 = IRNGT (IWRKD+2)
                   IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                   If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                      IRNGT (IWRKD+2) = IRNG1
                      If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
    !   3 1 2 4
                         IRNGT (IWRKD+3) = IRNG2
                      Else
    !   3 1 4 2
                         IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                         IRNGT (IWRKD+4) = IRNG2
                      End If
                   Else
    !   3 4 1 2
                      IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                      IRNGT (IWRKD+3) = IRNG1
                      IRNGT (IWRKD+4) = IRNG2
                   End If
                End If
             End Do
    !
    !  The Cs become As and Bs
    !
             LMTNA = 4
             Exit
          End Do
    !
    !  Iteration loop. Each time, the length of the ordered subsets
    !  is doubled.
    !
          Do
             If (LMTNA >= NVAL) Exit
             IWRKF = 0
             LMTNC = 2 * LMTNC
    !
    !   Loop on merges of A and B into C
    !
             Do
                IWRK = IWRKF
                IWRKD = IWRKF + 1
                JINDA = IWRKF + LMTNA
                IWRKF = IWRKF + LMTNC
                If (IWRKF >= NVAL) Then
                   If (JINDA >= NVAL) Exit
                   IWRKF = NVAL
                End If
                IINDA = 1
                IINDB = JINDA + 1
    !
    !   Shortcut for the case when the max of A is smaller
    !   than the min of B. This line may be activated when the
    !   initial set is already close to sorted.
    !
             IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
    !
    !  One steps in the C subset, that we build in the final rank array
    !
    !  Make a copy of the rank array for the merge iteration
    !
                JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
    !
                XVALA = XDONT (JWRKT(IINDA))
                XVALB = XDONT (IRNGT(IINDB))
    !
                Do
                   IWRK = IWRK + 1
    !
    !  We still have unprocessed values in both A and B
    !
                   If (XVALA > XVALB) Then
                      IRNGT (IWRK) = IRNGT (IINDB)
                      IINDB = IINDB + 1
                      If (IINDB > IWRKF) Then
    !  Only A still with unprocessed values
                         IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                         Exit
                      End If
                      XVALB = XDONT (IRNGT(IINDB))
                   Else
                      IRNGT (IWRK) = JWRKT (IINDA)
                      IINDA = IINDA + 1
                      If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                      XVALA = XDONT (JWRKT(IINDA))
                   End If
    !
                End Do
             End Do
    !
    !  The Cs become As and Bs
    !
             LMTNA = 2 * LMTNA
          End Do
    !
          Return
    !
    End Subroutine I_mrgrnk
end module 
