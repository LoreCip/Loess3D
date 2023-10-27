module sort_interface

    use iso_fortran_env, only: RK => real64
    
    public :: argsort, &
              sort

    private :: R_mrgrnk, I_mrgrnk, D_mrgrnk, &
               sort_real8, sort_real4, sort_integer
    
    interface argsort
      module procedure D_mrgrnk, &
                       R_mrgrnk, &
                       I_mrgrnk
    end interface argsort

    interface sort
        module procedure sort_real8, &
                         sort_real4, &
                         sort_integer
    end interface sort
    
contains

!------------------------------------------------------------------------
!               ARGSORTING ARRAYS                                       |
!------------------------------------------------------------------------
!--------------------------------------------------------------
! The mrgrnk subroutine applies a merge-sort ranking algorithm to an input array. It initially pairs adjacent
! elements based on their values and iteratively merges these pairs, progressively increasing the ordered subset
! size. This process continues until the entire array is sorted. The subroutine efficiently handles ordered
! subsets and minimizes unnecessary comparisons. It is O(n log n) in time and O(n) in memory.
!---------------------------------------------------------------  

    subroutine d_mrgrnk(xdont, irngt)
        real(RK), dimension(:), intent(in) :: xdont
        integer, dimension(:), intent(out) :: irngt
        real(RK) :: xvala, xvalb
        integer, dimension(size(irngt)) :: jwrkt
        integer :: lmtna, lmtnc, irng1, irng2
        integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
    
        nval = min(size(xdont), size(irngt))
        select case(nval)
        case(:0)
            return
        case(1)
            irngt(1) = 1
            return
        case default
            continue
        end select
    
        do iind = 2, nval, 2
            if (xdont(iind-1) <= xdont(iind)) then
                irngt(iind-1) = iind - 1
                irngt(iind) = iind
            else
                irngt(iind-1) = iind
                irngt(iind) = iind - 1
            end if
        end do
        if (modulo(nval, 2) /= 0) then
            irngt(nval) = nval
        end if
    
        lmtna = 2
        lmtnc = 4
    
        do
            if (nval <= 2) exit
    
            do iwrkd = 0, nval - 1, 4
                if ((iwrkd+4) > nval) then
                    if ((iwrkd+2) >= nval) exit
    
                    if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) exit
    
                    if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                        irng2 = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irng2
                    else
                        irng1 = irngt(iwrkd+1)
                        irngt(iwrkd+1) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irng1
                    end if
                    exit
                end if
    
                if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) cycle
    
                if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+2) = irngt(iwrkd+3)
                    if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                        irngt(iwrkd+3) = irng2
                    else
                        irngt(iwrkd+3) = irngt(iwrkd+4)
                        irngt(iwrkd+4) = irng2
                    end if
                else
                    irng1 = irngt(iwrkd+1)
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+1) = irngt(iwrkd+3)
                    if (xdont(irng1) <= xdont(irngt(iwrkd+4))) then
                        irngt(iwrkd+2) = irng1
                        if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                            irngt(iwrkd+3) = irng2
                        else
                            irngt(iwrkd+3) = irngt(iwrkd+4)
                            irngt(iwrkd+4) = irng2
                        end if
                    else
                        irngt(iwrkd+2) = irngt(iwrkd+4)
                        irngt(iwrkd+3) = irng1
                        irngt(iwrkd+4) = irng2
                    end if
                end if
            end do
    
            lmtna = 4
            exit
        end do
    
        do
            if (lmtna >= nval) exit
            iwrkf = 0
            lmtnc = 2 * lmtnc
    
            do
                iwrk = iwrkf
                iwrkd = iwrkf + 1
                jinda = iwrkf + lmtna
                iwrkf = iwrkf + lmtnc
                if (iwrkf >= nval) then
                    if (jinda >= nval) exit
                    iwrkf = nval
                end if
                iinda = 1
                iindb = jinda + 1
    
                if (xdont(irngt(jinda)) <= xdont(irngt(iindb))) cycle
    
                jwrkt(1:lmtna) = irngt(iwrkd:jinda)
                xvala = xdont(jwrkt(iinda))
                xvalb = xdont(irngt(iindb))
    
                do
                    iwrk = iwrk + 1
    
                    if (xvala > xvalb) then
                        irngt(iwrk) = irngt(iindb)
                        iindb = iindb + 1
                        if (iindb > iwrkf) then
                            irngt(iwrk+1:iwrkf) = jwrkt(iinda:lmtna)
                            exit
                        end if
                        xvalb = xdont(irngt(iindb))
                    else
                        irngt(iwrk) = jwrkt(iinda)
                        iinda = iinda + 1
                        if (iinda > lmtna) exit
                        xvala = xdont(jwrkt(iinda))
                    end if
                end do
            end do
    
            lmtna = 2 * lmtna
        end do
    
        return
    end subroutine d_mrgrnk
        

    subroutine r_mrgrnk(xdont, irngt)
        real, dimension(:), intent(in) :: xdont
        integer, dimension(:), intent(out) :: irngt
        real :: xvala, xvalb
        integer, dimension(size(irngt)) :: jwrkt
        integer :: lmtna, lmtnc, irng1, irng2
        integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
    
        nval = min(size(xdont), size(irngt))
        select case(nval)
        case(:0)
            return
        case(1)
            irngt(1) = 1
            return
        case default
            continue
        end select
    
        do iind = 2, nval, 2
            if (xdont(iind-1) <= xdont(iind)) then
                irngt(iind-1) = iind - 1
                irngt(iind) = iind
            else
                irngt(iind-1) = iind
                irngt(iind) = iind - 1
            end if
        end do
        if (modulo(nval, 2) /= 0) then
            irngt(nval) = nval
        end if
    
        lmtna = 2
        lmtnc = 4
    
        do
            if (nval <= 2) exit
    
            do iwrkd = 0, nval - 1, 4
                if ((iwrkd+4) > nval) then
                    if ((iwrkd+2) >= nval) exit
    
                    if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) exit
    
                    if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                        irng2 = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irng2
                    else
                        irng1 = irngt(iwrkd+1)
                        irngt(iwrkd+1) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irng1
                    end if
                    exit
                end if
    
                if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) cycle
    
                if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+2) = irngt(iwrkd+3)
                    if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                        irngt(iwrkd+3) = irng2
                    else
                        irngt(iwrkd+3) = irngt(iwrkd+4)
                        irngt(iwrkd+4) = irng2
                    end if
                else
                    irng1 = irngt(iwrkd+1)
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+1) = irngt(iwrkd+3)
                    if (xdont(irng1) <= xdont(irngt(iwrkd+4))) then
                        irngt(iwrkd+2) = irng1
                        if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                            irngt(iwrkd+3) = irng2
                        else
                            irngt(iwrkd+3) = irngt(iwrkd+4)
                            irngt(iwrkd+4) = irng2
                        end if
                    else
                        irngt(iwrkd+2) = irngt(iwrkd+4)
                        irngt(iwrkd+3) = irng1
                        irngt(iwrkd+4) = irng2
                    end if
                end if
            end do
    
            lmtna = 4
            exit
        end do
    
        do
            if (lmtna >= nval) exit
            iwrkf = 0
            lmtnc = 2 * lmtnc
    
            do
                iwrk = iwrkf
                iwrkd = iwrkf + 1
                jinda = iwrkf + lmtna
                iwrkf = iwrkf + lmtnc
                if (iwrkf >= nval) then
                    if (jinda >= nval) exit
                    iwrkf = nval
                end if
                iinda = 1
                iindb = jinda + 1
    
                if (xdont(irngt(jinda)) <= xdont(irngt(iindb))) cycle
    
                jwrkt(1:lmtna) = irngt(iwrkd:jinda)
                xvala = xdont(jwrkt(iinda))
                xvalb = xdont(irngt(iindb))
    
                do
                    iwrk = iwrk + 1
    
                    if (xvala > xvalb) then
                        irngt(iwrk) = irngt(iindb)
                        iindb = iindb + 1
                        if (iindb > iwrkf) then
                            irngt(iwrk+1:iwrkf) = jwrkt(iinda:lmtna)
                            exit
                        end if
                        xvalb = xdont(irngt(iindb))
                    else
                        irngt(iwrk) = jwrkt(iinda)
                        iinda = iinda + 1
                        if (iinda > lmtna) exit
                        xvala = xdont(jwrkt(iinda))
                    end if
                end do
            end do
    
            lmtna = 2 * lmtna
        end do
    
        return
    end subroutine r_mrgrnk
    

    subroutine i_mrgrnk(xdont, irngt)
        integer, dimension(:), intent(in) :: xdont
        integer, dimension(:), intent(out) :: irngt
        integer :: xvala, xvalb
        integer, dimension(size(irngt)) :: jwrkt
        integer :: lmtna, lmtnc, irng1, irng2
        integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
    
        nval = min(size(xdont), size(irngt))
        select case(nval)
        case(:0)
            return
        case(1)
            irngt(1) = 1
            return
        case default
            continue
        end select
    
        do iind = 2, nval, 2
            if (xdont(iind-1) <= xdont(iind)) then
                irngt(iind-1) = iind - 1
                irngt(iind) = iind
            else
                irngt(iind-1) = iind
                irngt(iind) = iind - 1
            end if
        end do
        if (modulo(nval, 2) /= 0) then
            irngt(nval) = nval
        end if
    
        lmtna = 2
        lmtnc = 4
    
        do
            if (nval <= 2) exit
    
            do iwrkd = 0, nval - 1, 4
                if ((iwrkd+4) > nval) then
                    if ((iwrkd+2) >= nval) exit
    
                    if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) exit
    
                    if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                        irng2 = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irng2
                    else
                        irng1 = irngt(iwrkd+1)
                        irngt(iwrkd+1) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irng1
                    end if
                    exit
                end if
    
                if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) cycle
    
                if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+2) = irngt(iwrkd+3)
                    if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                        irngt(iwrkd+3) = irng2
                    else
                        irngt(iwrkd+3) = irngt(iwrkd+4)
                        irngt(iwrkd+4) = irng2
                    end if
                else
                    irng1 = irngt(iwrkd+1)
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+1) = irngt(iwrkd+3)
                    if (xdont(irng1) <= xdont(irngt(iwrkd+4))) then
                        irngt(iwrkd+2) = irng1
                        if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                            irngt(iwrkd+3) = irng2
                        else
                            irngt(iwrkd+3) = irngt(iwrkd+4)
                            irngt(iwrkd+4) = irng2
                        end if
                    else
                        irngt(iwrkd+2) = irngt(iwrkd+4)
                        irngt(iwrkd+3) = irng1
                        irngt(iwrkd+4) = irng2
                    end if
                end if
            end do
    
            lmtna = 4
            exit
        end do
    
        do
            if (lmtna >= nval) exit
            iwrkf = 0
            lmtnc = 2 * lmtnc
    
            do
                iwrk = iwrkf
                iwrkd = iwrkf + 1
                jinda = iwrkf + lmtna
                iwrkf = iwrkf + lmtnc
                if (iwrkf >= nval) then
                    if (jinda >= nval) exit
                    iwrkf = nval
                end if
                iinda = 1
                iindb = jinda + 1
    
                if (xdont(irngt(jinda)) <= xdont(irngt(iindb))) cycle
    
                jwrkt(1:lmtna) = irngt(iwrkd:jinda)
                xvala = xdont(jwrkt(iinda))
                xvalb = xdont(irngt(iindb))
    
                do
                    iwrk = iwrk + 1
    
                    if (xvala > xvalb) then
                        irngt(iwrk) = irngt(iindb)
                        iindb = iindb + 1
                        if (iindb > iwrkf) then
                            irngt(iwrk+1:iwrkf) = jwrkt(iinda:lmtna)
                            exit
                        end if
                        xvalb = xdont(irngt(iindb))
                    else
                        irngt(iwrk) = jwrkt(iinda)
                        iinda = iinda + 1
                        if (iinda > lmtna) exit
                        xvala = xdont(jwrkt(iinda))
                    end if
                end do
            end do
    
            lmtna = 2 * lmtna
        end do
    
        return
    end subroutine i_mrgrnk
    
!------------------------------------------------------------------------
!              SORTING ARRAY                                            |
!------------------------------------------------------------------------

    subroutine sort_real8(xin, xout)
        real(RK), dimension(:), intent(in)  :: xin
        real(RK), dimension(:), intent(out) :: xout
        integer, dimension(size(xin)) :: inds
        
        call D_mrgrnk(xin, inds)
        xout = xin(inds)
    end subroutine sort_real8
    subroutine sort_real4(xin, xout)
        real, dimension(:), intent(in)  :: xin
        real, dimension(:), intent(out) :: xout
        integer, dimension(size(xin)) :: inds
        
        call R_mrgrnk(xin, inds)
        xout = xin(inds)
    end subroutine sort_real4
    subroutine sort_integer(xin, xout)
        integer, dimension(:), intent(in)  :: xin
        integer, dimension(:), intent(out) :: xout
        integer, dimension(size(xin)) :: inds
        
        call I_mrgrnk(xin, inds)
        xout = xin(inds)
    end subroutine sort_integer

end module sort_interface
