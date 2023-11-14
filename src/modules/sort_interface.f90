module sort_interface

    use iso_fortran_env, only: RK => real64
    
    public :: argsort, sort, &
              pargsort, psort

    private :: R_mrgrnk, I_mrgrnk, D_mrgrnk, &
               sort_real8, sort_real4, sort_integer, &
               R_rnkpar, I_rnkpar, D_rnkpar, &
               psort_real8, psort_real4, psort_integer

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

    interface pargsort
        module procedure d_rnkpar, &
                         r_rnkpar, &
                         i_rnkpar
    end interface pargsort

    interface psort
        module procedure psort_real8, &
                         psort_real4, &
                         psort_integer
    end interface psort

    interface median
         module procedure d_valmed, r_valmed, i_valmed
    end interface median
    
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

    
 !------------------------------------------------------------------------
 !           PARTIAL ARGSORTING ARRAY                                    |
 !------------------------------------------------------------------------


    subroutine d_rnkpar (xdont, irngt, nord)
        real (RK), dimension (:), intent (in) :: xdont
        integer, dimension (:), intent (out) :: irngt
        integer, intent (in) :: nord
        real (RK) :: xpiv, xpiv0, xwrk, xwrk1, xmin, xmax
        integer, dimension (size(xdont)) :: ilowt, ihigt
        integer :: ndon, jhig, jlow, ihig, iwrk, iwrk1, iwrk2, iwrk3
        integer :: ideb, jdeb, imil, ifin, nwrk, icrs, idcr, ilow
        integer :: jlm2, jlm1, jhm2, jhm1
        ndon = size (xdont)
        if (ndon < 2) then
           if (nord >= 1) irngt (1) = 1
           return
        end if

        if (xdont(2) < xdont(1)) then
           ilowt (1) = 2
           ihigt (1) = 1
        else
           ilowt (1) = 1
           ihigt (1) = 2
        end if

        if (ndon < 3) then
           if (nord >= 1) irngt (1) = ilowt (1)
           if (nord >= 2) irngt (2) = ihigt (1)
           return
        end if

        if (xdont(3) < xdont(ihigt(1))) then
           ihigt (2) = ihigt (1)
           if (xdont(3) < xdont(ilowt(1))) then
              ihigt (1) = ilowt (1)
              ilowt (1) = 3
           else
              ihigt (1) = 3
           end if

        else
           ihigt (2) = 3
        end if

        if (ndon < 4) then
           if (nord >= 1) irngt (1) = ilowt (1)
           if (nord >= 2) irngt (2) = ihigt (1)
           if (nord >= 3) irngt (3) = ihigt (2)
           return
        end if

        if (xdont(ndon) < xdont(ihigt(1))) then
           ihigt (3) = ihigt (2)
           ihigt (2) = ihigt (1)
           if (xdont(ndon) < xdont(ilowt(1))) then
              ihigt (1) = ilowt (1)
              ilowt (1) = ndon
           else
              ihigt (1) = ndon
           end if

        else
           ihigt (3) = ndon
        end if

        if (ndon < 5) then
           if (nord >= 1) irngt (1) = ilowt (1)
           if (nord >= 2) irngt (2) = ihigt (1)
           if (nord >= 3) irngt (3) = ihigt (2)
           if (nord >= 4) irngt (4) = ihigt (3)
           return
        end if

        jdeb = 0
        ideb = jdeb + 1
        jlow = ideb
        jhig = 3
        xpiv = xdont (ilowt(ideb)) + real(2*nord)/real(ndon+nord) * &
                                     (xdont(ihigt(3))-xdont(ilowt(ideb)))
        if (xpiv >= xdont(ihigt(1))) then
           xpiv = xdont (ilowt(ideb)) + real(2*nord)/real(ndon+nord) * &
                                        (xdont(ihigt(2))-xdont(ilowt(ideb)))
           if (xpiv >= xdont(ihigt(1))) &
               xpiv = xdont (ilowt(ideb)) + real (2*nord) / real (ndon+nord) * &
                                            (xdont(ihigt(1))-xdont(ilowt(ideb)))
        end if

        xpiv0 = xpiv
        if (xdont(ndon) > xpiv) then
           icrs = 3
           do
              icrs = icrs + 1
              if (xdont(icrs) > xpiv) then
                 if (icrs >= ndon) exit
                 jhig = jhig + 1
                 ihigt (jhig) = icrs
              else
                 jlow = jlow + 1
                 ilowt (jlow) = icrs
                 if (jlow >= nord) exit
              end if

           end do

           if (icrs < ndon-1) then
              do
                 icrs = icrs + 1
                 if (xdont(icrs) <= xpiv) then
                    jlow = jlow + 1
                    ilowt (jlow) = icrs
                 else if (icrs >= ndon) then
                    exit
                 end if

              end do

           end if

        else
           do icrs = 4, ndon - 1
              if (xdont(icrs) > xpiv) then
                 jhig = jhig + 1
                 ihigt (jhig) = icrs
              else
                 jlow = jlow + 1
                 ilowt (jlow) = icrs
                 if (jlow >= nord) exit
              end if

           end do

           if (icrs < ndon-1) then
              do
                 icrs = icrs + 1
                 if (xdont(icrs) <= xpiv) then
                    if (icrs >= ndon) exit
                    jlow = jlow + 1
                    ilowt (jlow) = icrs
                 end if

              end do

           end if

        end if

        jlm2 = 0
        jlm1 = 0
        jhm2 = 0
        jhm1 = 0
        do
           if (jlow == nord) exit
           if (jlm2 == jlow .and. jhm2 == jhig) then
             if (nord > jlow) then
                  xmin = xdont (ihigt(1))
                  ihig = 1
                  do icrs = 2, jhig
                     if (xdont(ihigt(icrs)) < xmin) then
                        xmin = xdont (ihigt(icrs))
                        ihig = icrs
                     end if

                  end do

                  jlow = jlow + 1
                  ilowt (jlow) = ihigt (ihig)
                  ihigt (ihig) = ihigt (jhig)
                  jhig = jhig - 1
               else
                  ilow = ilowt (jlow)
                  xmax = xdont (ilow)
                  do icrs = 1, jlow
                     if (xdont(ilowt(icrs)) > xmax) then
                        iwrk = ilowt (icrs)
                        xmax = xdont (iwrk)
                        ilowt (icrs) = ilow
                        ilow = iwrk
                     end if

                  end do

                  jlow = jlow - 1
               end if

           end if

           jlm2 = jlm1
           jlm1 = jlow
           jhm2 = jhm1
           jhm1 = jhig
          select case (nord-jlow)
           case (2:)
              select case (jhig)
              case (2)
                 if (xdont(ihigt(1)) <= xdont(ihigt(2))) then
                    jlow = jlow + 1
                    ilowt (jlow) = ihigt (1)
                    jlow = jlow + 1
                    ilowt (jlow) = ihigt (2)
                 else
                    jlow = jlow + 1
                    ilowt (jlow) = ihigt (2)
                    jlow = jlow + 1
                    ilowt (jlow) = ihigt (1)
                 end if

                 exit
              case (3)
                 iwrk1 = ihigt (1)
                 iwrk2 = ihigt (2)
                 iwrk3 = ihigt (3)
                 if (xdont(iwrk2) < xdont(iwrk1)) then
                    ihigt (1) = iwrk2
                    ihigt (2) = iwrk1
                    iwrk2 = iwrk1
                 end if

                 if (xdont(iwrk2) > xdont(iwrk3)) then
                    ihigt (3) = iwrk2
                    ihigt (2) = iwrk3
                    iwrk2 = iwrk3
                    if (xdont(iwrk2) < xdont(ihigt(1))) then
                       ihigt (2) = ihigt (1)
                       ihigt (1) = iwrk2
                    end if

                 end if

                 jhig = 0
                 do icrs = jlow + 1, nord
                    jhig = jhig + 1
                    ilowt (icrs) = ihigt (jhig)
                 end do

                 jlow = nord
                 exit
              case (4:)
                 xpiv0 = xpiv
                 ifin = jhig
                 iwrk1 = ihigt (1)
                 iwrk2 = ihigt (2)
                 iwrk3 = ihigt (ifin)
                 if (xdont(iwrk2) < xdont(iwrk1)) then
                    ihigt (1) = iwrk2
                    ihigt (2) = iwrk1
                    iwrk2 = iwrk1
                 end if

                 if (xdont(iwrk2) > xdont(iwrk3)) then
                    ihigt (ifin) = iwrk2
                    ihigt (2) = iwrk3
                    iwrk2 = iwrk3
                    if (xdont(iwrk2) < xdont(ihigt(1))) then
                       ihigt (2) = ihigt (1)
                       ihigt (1) = iwrk2
                    end if

                 end if

                 jdeb = jlow
                 nwrk = nord - jlow
                 iwrk1 = ihigt (1)
                 jlow = jlow + 1
                 ilowt (jlow) = iwrk1
                 xpiv = xdont (iwrk1) + real (nwrk) / real (nord+nwrk) * &
                                        (xdont(ihigt(ifin))-xdont(iwrk1))
                 jhig = 0
                 do icrs = 2, ifin
                    if (xdont(ihigt(icrs)) <= xpiv) then
                       jlow = jlow + 1
                       ilowt (jlow) = ihigt (icrs)
                       if (jlow >= nord) exit
                    else
                       jhig = jhig + 1
                       ihigt (jhig) = ihigt (icrs)
                    end if

                 end do

                 do icrs = icrs + 1, ifin
                    if (xdont(ihigt(icrs)) <= xpiv) then
                       jlow = jlow + 1
                       ilowt (jlow) = ihigt (icrs)
                    end if

                 end do

             end select

           case (1)
              xmin = xdont (ihigt(1))
              ihig = 1
              do icrs = 2, jhig
                 if (xdont(ihigt(icrs)) < xmin) then
                    xmin = xdont (ihigt(icrs))
                    ihig = icrs
                 end if

              end do

              jlow = jlow + 1
              ilowt (jlow) = ihigt (ihig)
              exit
           case (0)
              exit
           case (-5:-1)
              irngt (1) = ilowt (1)
              do icrs = 2, nord
                 iwrk = ilowt (icrs)
                 xwrk = xdont (iwrk)
                 do idcr = icrs - 1, 1, - 1
                    if (xwrk < xdont(irngt(idcr))) then
                       irngt (idcr+1) = irngt (idcr)
                    else
                       exit
                    end if

                 end do

                 irngt (idcr+1) = iwrk
              end do

              xwrk1 = xdont (irngt(nord))
              do icrs = nord + 1, jlow
                 if (xdont(ilowt (icrs)) < xwrk1) then
                    xwrk = xdont (ilowt (icrs))
                    do idcr = nord - 1, 1, - 1
                       if (xwrk >= xdont(irngt(idcr))) exit
                       irngt (idcr+1) = irngt (idcr)
                    end do

                    irngt (idcr+1) = ilowt (icrs)
                    xwrk1 = xdont (irngt(nord))
                 end if

              end do

              return
           case (:-6)
              ideb = jdeb + 1
              imil = (jlow+ideb) / 2
              ifin = jlow
              if (xdont(ilowt(imil)) < xdont(ilowt(ideb))) then
                 iwrk = ilowt (ideb)
                 ilowt (ideb) = ilowt (imil)
                 ilowt (imil) = iwrk
              end if

              if (xdont(ilowt(imil)) > xdont(ilowt(ifin))) then
                 iwrk = ilowt (ifin)
                 ilowt (ifin) = ilowt (imil)
                 ilowt (imil) = iwrk
                 if (xdont(ilowt(imil)) < xdont(ilowt(ideb))) then
                    iwrk = ilowt (ideb)
                    ilowt (ideb) = ilowt (imil)
                    ilowt (imil) = iwrk
                 end if

              end if

              if (ifin <= 3) exit
              xpiv = xdont (ilowt(1)) + real(nord)/real(jlow+nord) * &
                                        (xdont(ilowt(ifin))-xdont(ilowt(1)))
              if (jdeb > 0) then
                 if (xpiv <= xpiv0) &
                     xpiv = xpiv0 + real(2*nord-jdeb)/real (jlow+nord) * &
                                    (xdont(ilowt(ifin))-xpiv0)
              else
                 ideb = 1
              end if

              jhig = 0
              jlow = jdeb
              if (xdont(ilowt(ifin)) > xpiv) then
                 icrs = jdeb
                 do
                   icrs = icrs + 1
                    if (xdont(ilowt(icrs)) > xpiv) then
                       jhig = jhig + 1
                       ihigt (jhig) = ilowt (icrs)
                       if (icrs >= ifin) exit
                    else
                       jlow = jlow + 1
                       ilowt (jlow) = ilowt (icrs)
                       if (jlow >= nord) exit
                    end if

                 end do

                 if (icrs < ifin) then
                    do
                       icrs = icrs + 1
                       if (xdont(ilowt(icrs)) <= xpiv) then
                          jlow = jlow + 1
                          ilowt (jlow) = ilowt (icrs)
                       else
                          if (icrs >= ifin) exit
                       end if

                    end do

                 end if

             else
                 do icrs = ideb, ifin
                    if (xdont(ilowt(icrs)) > xpiv) then
                       jhig = jhig + 1
                       ihigt (jhig) = ilowt (icrs)
                    else
                       jlow = jlow + 1
                       ilowt (jlow) = ilowt (icrs)
                       if (jlow >= nord) exit
                    end if

                 end do

                 do icrs = icrs + 1, ifin
                    if (xdont(ilowt(icrs)) <= xpiv) then
                       jlow = jlow + 1
                       ilowt (jlow) = ilowt (icrs)
                    end if

                 end do

              end if

           end select

        end do

        irngt (1) = ilowt (1)
        do icrs = 2, nord
           iwrk = ilowt (icrs)
           xwrk = xdont (iwrk)
           do idcr = icrs - 1, 1, - 1
              if (xwrk < xdont(irngt(idcr))) then
                 irngt (idcr+1) = irngt (idcr)
              else
                 exit
              end if

           end do

           irngt (idcr+1) = iwrk
        end do

       return
    end subroutine d_rnkpar

  
    subroutine r_rnkpar (xdont, irngt, nord)
        real, dimension (:), intent (in) :: xdont
        integer, dimension (:), intent (out) :: irngt
        integer, intent (in) :: nord
        real    :: xpiv, xpiv0, xwrk, xwrk1, xmin, xmax
        integer, dimension (size(xdont)) :: ilowt, ihigt
        integer :: ndon, jhig, jlow, ihig, iwrk, iwrk1, iwrk2, iwrk3
        integer :: ideb, jdeb, imil, ifin, nwrk, icrs, idcr, ilow
        integer :: jlm2, jlm1, jhm2, jhm1
        ndon = size (xdont)
        if (ndon < 2) then
           if (nord >= 1) irngt (1) = 1
           return
        end if

        if (xdont(2) < xdont(1)) then
           ilowt (1) = 2
           ihigt (1) = 1
        else
           ilowt (1) = 1
           ihigt (1) = 2
        end if

        if (ndon < 3) then
           if (nord >= 1) irngt (1) = ilowt (1)
           if (nord >= 2) irngt (2) = ihigt (1)
           return
        end if

        if (xdont(3) < xdont(ihigt(1))) then
           ihigt (2) = ihigt (1)
           if (xdont(3) < xdont(ilowt(1))) then
              ihigt (1) = ilowt (1)
              ilowt (1) = 3
           else
              ihigt (1) = 3
           end if

        else
           ihigt (2) = 3
        end if

        if (ndon < 4) then
           if (nord >= 1) irngt (1) = ilowt (1)
           if (nord >= 2) irngt (2) = ihigt (1)
           if (nord >= 3) irngt (3) = ihigt (2)
           return
        end if

        if (xdont(ndon) < xdont(ihigt(1))) then
           ihigt (3) = ihigt (2)
           ihigt (2) = ihigt (1)
           if (xdont(ndon) < xdont(ilowt(1))) then
              ihigt (1) = ilowt (1)
              ilowt (1) = ndon
           else
              ihigt (1) = ndon
           end if

        else
           ihigt (3) = ndon
        end if

        if (ndon < 5) then
           if (nord >= 1) irngt (1) = ilowt (1)
           if (nord >= 2) irngt (2) = ihigt (1)
           if (nord >= 3) irngt (3) = ihigt (2)
           if (nord >= 4) irngt (4) = ihigt (3)
           return
        end if

        jdeb = 0
        ideb = jdeb + 1
        jlow = ideb
        jhig = 3
        xpiv = xdont (ilowt(ideb)) + real(2*nord)/real(ndon+nord) * &
                                     (xdont(ihigt(3))-xdont(ilowt(ideb)))
        if (xpiv >= xdont(ihigt(1))) then
           xpiv = xdont (ilowt(ideb)) + real(2*nord)/real(ndon+nord) * &
                                        (xdont(ihigt(2))-xdont(ilowt(ideb)))
           if (xpiv >= xdont(ihigt(1))) &
               xpiv = xdont (ilowt(ideb)) + real (2*nord) / real (ndon+nord) * &
                                            (xdont(ihigt(1))-xdont(ilowt(ideb)))
        end if

        xpiv0 = xpiv
        if (xdont(ndon) > xpiv) then
           icrs = 3
           do
              icrs = icrs + 1
              if (xdont(icrs) > xpiv) then
                 if (icrs >= ndon) exit
                 jhig = jhig + 1
                 ihigt (jhig) = icrs
              else
                 jlow = jlow + 1
                 ilowt (jlow) = icrs
                 if (jlow >= nord) exit
              end if

           end do

           if (icrs < ndon-1) then
              do
                 icrs = icrs + 1
                 if (xdont(icrs) <= xpiv) then
                    jlow = jlow + 1
                    ilowt (jlow) = icrs
                 else if (icrs >= ndon) then
                    exit
                 end if

              end do

           end if

        else
           do icrs = 4, ndon - 1
              if (xdont(icrs) > xpiv) then
                 jhig = jhig + 1
                 ihigt (jhig) = icrs
              else
                 jlow = jlow + 1
                 ilowt (jlow) = icrs
                 if (jlow >= nord) exit
              end if

           end do

           if (icrs < ndon-1) then
              do
                 icrs = icrs + 1
                 if (xdont(icrs) <= xpiv) then
                    if (icrs >= ndon) exit
                    jlow = jlow + 1
                    ilowt (jlow) = icrs
                 end if

              end do

           end if

        end if

        jlm2 = 0
        jlm1 = 0
        jhm2 = 0
        jhm1 = 0
        do
           if (jlow == nord) exit
           if (jlm2 == jlow .and. jhm2 == jhig) then
             if (nord > jlow) then
                  xmin = xdont (ihigt(1))
                  ihig = 1
                  do icrs = 2, jhig
                     if (xdont(ihigt(icrs)) < xmin) then
                        xmin = xdont (ihigt(icrs))
                        ihig = icrs
                     end if

                  end do

                  jlow = jlow + 1
                  ilowt (jlow) = ihigt (ihig)
                  ihigt (ihig) = ihigt (jhig)
                  jhig = jhig - 1
               else
                  ilow = ilowt (jlow)
                  xmax = xdont (ilow)
                  do icrs = 1, jlow
                     if (xdont(ilowt(icrs)) > xmax) then
                        iwrk = ilowt (icrs)
                        xmax = xdont (iwrk)
                        ilowt (icrs) = ilow
                        ilow = iwrk
                     end if

                  end do

                  jlow = jlow - 1
               end if

           end if

           jlm2 = jlm1
           jlm1 = jlow
           jhm2 = jhm1
           jhm1 = jhig
          select case (nord-jlow)
           case (2:)
              select case (jhig)
              case (2)
                 if (xdont(ihigt(1)) <= xdont(ihigt(2))) then
                    jlow = jlow + 1
                    ilowt (jlow) = ihigt (1)
                    jlow = jlow + 1
                    ilowt (jlow) = ihigt (2)
                 else
                    jlow = jlow + 1
                    ilowt (jlow) = ihigt (2)
                    jlow = jlow + 1
                    ilowt (jlow) = ihigt (1)
                 end if

                 exit
              case (3)
                 iwrk1 = ihigt (1)
                 iwrk2 = ihigt (2)
                 iwrk3 = ihigt (3)
                 if (xdont(iwrk2) < xdont(iwrk1)) then
                    ihigt (1) = iwrk2
                    ihigt (2) = iwrk1
                    iwrk2 = iwrk1
                 end if

                 if (xdont(iwrk2) > xdont(iwrk3)) then
                    ihigt (3) = iwrk2
                    ihigt (2) = iwrk3
                    iwrk2 = iwrk3
                    if (xdont(iwrk2) < xdont(ihigt(1))) then
                       ihigt (2) = ihigt (1)
                       ihigt (1) = iwrk2
                    end if

                 end if

                 jhig = 0
                 do icrs = jlow + 1, nord
                    jhig = jhig + 1
                    ilowt (icrs) = ihigt (jhig)
                 end do

                 jlow = nord
                 exit
              case (4:)
                 xpiv0 = xpiv
                 ifin = jhig
                 iwrk1 = ihigt (1)
                 iwrk2 = ihigt (2)
                 iwrk3 = ihigt (ifin)
                 if (xdont(iwrk2) < xdont(iwrk1)) then
                    ihigt (1) = iwrk2
                    ihigt (2) = iwrk1
                    iwrk2 = iwrk1
                 end if

                 if (xdont(iwrk2) > xdont(iwrk3)) then
                    ihigt (ifin) = iwrk2
                    ihigt (2) = iwrk3
                    iwrk2 = iwrk3
                    if (xdont(iwrk2) < xdont(ihigt(1))) then
                       ihigt (2) = ihigt (1)
                       ihigt (1) = iwrk2
                    end if

                 end if

                 jdeb = jlow
                 nwrk = nord - jlow
                 iwrk1 = ihigt (1)
                 jlow = jlow + 1
                 ilowt (jlow) = iwrk1
                 xpiv = xdont (iwrk1) + real (nwrk) / real (nord+nwrk) * &
                                        (xdont(ihigt(ifin))-xdont(iwrk1))
                 jhig = 0
                 do icrs = 2, ifin
                    if (xdont(ihigt(icrs)) <= xpiv) then
                       jlow = jlow + 1
                       ilowt (jlow) = ihigt (icrs)
                       if (jlow >= nord) exit
                    else
                       jhig = jhig + 1
                       ihigt (jhig) = ihigt (icrs)
                    end if

                 end do

                 do icrs = icrs + 1, ifin
                    if (xdont(ihigt(icrs)) <= xpiv) then
                       jlow = jlow + 1
                       ilowt (jlow) = ihigt (icrs)
                    end if
                 end do
             end select

           case (1)
              xmin = xdont (ihigt(1))
              ihig = 1
              do icrs = 2, jhig
                 if (xdont(ihigt(icrs)) < xmin) then
                    xmin = xdont (ihigt(icrs))
                    ihig = icrs
                 end if
              end do
              jlow = jlow + 1
              ilowt (jlow) = ihigt (ihig)
              exit

           case (0)
              exit

           case (-5:-1)
              irngt (1) = ilowt (1)
              do icrs = 2, nord
                 iwrk = ilowt (icrs)
                 xwrk = xdont (iwrk)
                 do idcr = icrs - 1, 1, - 1
                    if (xwrk < xdont(irngt(idcr))) then
                       irngt (idcr+1) = irngt (idcr)
                    else
                       exit
                    end if
                 end do

                 irngt (idcr+1) = iwrk
              end do

              xwrk1 = xdont (irngt(nord))
              do icrs = nord + 1, jlow
                    if (xdont(ilowt (icrs)) < xwrk1) then
                        xwrk = xdont (ilowt (icrs))
                        do idcr = nord - 1, 1, - 1
                           if (xwrk >= xdont(irngt(idcr))) exit
                           irngt (idcr+1) = irngt (idcr)
                        end do 
                        irngt (idcr+1) = ilowt (icrs)
                        xwrk1 = xdont (irngt(nord))
                    end if
              end do
              return

           case (:-6)
              ideb = jdeb + 1
              imil = (jlow+ideb) / 2
              ifin = jlow
              if (xdont(ilowt(imil)) < xdont(ilowt(ideb))) then
                 iwrk = ilowt (ideb)
                 ilowt (ideb) = ilowt (imil)
                 ilowt (imil) = iwrk
              end if

              if (xdont(ilowt(imil)) > xdont(ilowt(ifin))) then
                 iwrk = ilowt (ifin)
                 ilowt (ifin) = ilowt (imil)
                 ilowt (imil) = iwrk
                 if (xdont(ilowt(imil)) < xdont(ilowt(ideb))) then
                    iwrk = ilowt (ideb)
                    ilowt (ideb) = ilowt (imil)
                    ilowt (imil) = iwrk
                 end if
              end if

              if (ifin <= 3) exit
              xpiv = xdont (ilowt(1)) + real(nord)/real(jlow+nord) * &
                                        (xdont(ilowt(ifin))-xdont(ilowt(1)))
              if (jdeb > 0) then
                 if (xpiv <= xpiv0) &
                     xpiv = xpiv0 + real(2*nord-jdeb)/real (jlow+nord) * &
                                    (xdont(ilowt(ifin))-xpiv0)
              else
                 ideb = 1
              end if

              jhig = 0
              jlow = jdeb
              if (xdont(ilowt(ifin)) > xpiv) then
                 icrs = jdeb
                 do
                   icrs = icrs + 1
                    if (xdont(ilowt(icrs)) > xpiv) then
                       jhig = jhig + 1
                       ihigt (jhig) = ilowt (icrs)
                       if (icrs >= ifin) exit
                    else
                       jlow = jlow + 1
                       ilowt (jlow) = ilowt (icrs)
                       if (jlow >= nord) exit
                    end if
                 end do

                 if (icrs < ifin) then
                    do
                       icrs = icrs + 1
                       if (xdont(ilowt(icrs)) <= xpiv) then
                          jlow = jlow + 1
                          ilowt (jlow) = ilowt (icrs)
                       else
                          if (icrs >= ifin) exit
                       end if
                    end do
                 end if

             else

                 do icrs = ideb, ifin
                    if (xdont(ilowt(icrs)) > xpiv) then
                       jhig = jhig + 1
                       ihigt (jhig) = ilowt (icrs)
                    else
                       jlow = jlow + 1
                       ilowt (jlow) = ilowt (icrs)
                       if (jlow >= nord) exit
                    end if

                 end do

                 do icrs = icrs + 1, ifin
                    if (xdont(ilowt(icrs)) <= xpiv) then
                       jlow = jlow + 1
                       ilowt (jlow) = ilowt (icrs)
                    end if
                 end do

              end if
           end select
        end do

        irngt (1) = ilowt (1)
        do icrs = 2, nord
           iwrk = ilowt (icrs)
           xwrk = xdont (iwrk)
           do idcr = icrs - 1, 1, - 1
              if (xwrk < xdont(irngt(idcr))) then
                 irngt (idcr+1) = irngt (idcr)
              else
                 exit
              end if
           end do
           irngt (idcr+1) = iwrk
        end do

       return
    end subroutine r_rnkpar

    subroutine i_rnkpar (xdont, irngt, nord)
        integer, dimension (:), intent (in)  :: xdont
        integer, dimension (:), intent (out) :: irngt
        integer, intent (in) :: nord
        integer :: xpiv, xpiv0, xwrk, xwrk1, xmin, xmax
        integer, dimension (size(xdont)) :: ilowt, ihigt
        integer :: ndon, jhig, jlow, ihig, iwrk, iwrk1, iwrk2, iwrk3
        integer :: ideb, jdeb, imil, ifin, nwrk, icrs, idcr, ilow
        integer :: jlm2, jlm1, jhm2, jhm1
        
        ndon = size (xdont)
        if (ndon < 2) then
           if (nord >= 1) irngt (1) = 1
           return
        end if

        if (xdont(2) < xdont(1)) then
           ilowt (1) = 2
           ihigt (1) = 1
        else
           ilowt (1) = 1
           ihigt (1) = 2
        end if

        if (ndon < 3) then
           if (nord >= 1) irngt (1) = ilowt (1)
           if (nord >= 2) irngt (2) = ihigt (1)
           return
        end if

        if (xdont(3) < xdont(ihigt(1))) then
           ihigt (2) = ihigt (1)
           if (xdont(3) < xdont(ilowt(1))) then
              ihigt (1) = ilowt (1)
              ilowt (1) = 3
           else
              ihigt (1) = 3
           end if

        else
           ihigt (2) = 3
        end if

        if (ndon < 4) then
           if (nord >= 1) irngt (1) = ilowt (1)
           if (nord >= 2) irngt (2) = ihigt (1)
           if (nord >= 3) irngt (3) = ihigt (2)
           return
        end if

        if (xdont(ndon) < xdont(ihigt(1))) then
           ihigt (3) = ihigt (2)
           ihigt (2) = ihigt (1)
           if (xdont(ndon) < xdont(ilowt(1))) then
              ihigt (1) = ilowt (1)
              ilowt (1) = ndon
           else
              ihigt (1) = ndon
           end if

        else
           ihigt (3) = ndon
        end if

        if (ndon < 5) then
           if (nord >= 1) irngt (1) = ilowt (1)
           if (nord >= 2) irngt (2) = ihigt (1)
           if (nord >= 3) irngt (3) = ihigt (2)
           if (nord >= 4) irngt (4) = ihigt (3)
           return
        end if

        jdeb = 0
        ideb = jdeb + 1
        jlow = ideb
        jhig = 3
        xpiv = xdont (ilowt(ideb)) + real(2*nord)/real(ndon+nord) * &
                                     (xdont(ihigt(3))-xdont(ilowt(ideb)))
        if (xpiv >= xdont(ihigt(1))) then
           xpiv = xdont (ilowt(ideb)) + real(2*nord)/real(ndon+nord) * &
                                        (xdont(ihigt(2))-xdont(ilowt(ideb)))
           if (xpiv >= xdont(ihigt(1))) &
               xpiv = xdont (ilowt(ideb)) + real (2*nord) / real(ndon+nord) * &
                                            (xdont(ihigt(1))-xdont(ilowt(ideb)))
        end if

        xpiv0 = xpiv
        if (xdont(ndon) > xpiv) then
           icrs = 3
           do
              icrs = icrs + 1
              if (xdont(icrs) > xpiv) then
                 if (icrs >= ndon) exit
                 jhig = jhig + 1
                 ihigt (jhig) = icrs
              else
                 jlow = jlow + 1
                 ilowt (jlow) = icrs
                 if (jlow >= nord) exit
              end if

           end do

           if (icrs < ndon-1) then
              do
                 icrs = icrs + 1
                 if (xdont(icrs) <= xpiv) then
                    jlow = jlow + 1
                    ilowt (jlow) = icrs
                 else if (icrs >= ndon) then
                    exit
                 end if
              end do
           end if

        else
           do icrs = 4, ndon - 1
              if (xdont(icrs) > xpiv) then
                 jhig = jhig + 1
                 ihigt (jhig) = icrs
              else
                 jlow = jlow + 1
                 ilowt (jlow) = icrs
                 if (jlow >= nord) exit
              end if
           end do

           if (icrs < ndon-1) then
              do
                 icrs = icrs + 1
                 if (xdont(icrs) <= xpiv) then
                    if (icrs >= ndon) exit
                    jlow = jlow + 1
                    ilowt (jlow) = icrs
                 end if
              end do
           end if
        end if

        jlm2 = 0
        jlm1 = 0
        jhm2 = 0
        jhm1 = 0
        do
           if (jlow == nord) exit
           if (jlm2 == jlow .and. jhm2 == jhig) then
             if (nord > jlow) then
                  xmin = xdont (ihigt(1))
                  ihig = 1
                  do icrs = 2, jhig
                     if (xdont(ihigt(icrs)) < xmin) then
                        xmin = xdont (ihigt(icrs))
                        ihig = icrs
                     end if
                  end do

                  jlow = jlow + 1
                  ilowt (jlow) = ihigt (ihig)
                  ihigt (ihig) = ihigt (jhig)
                  jhig = jhig - 1
               else
                  ilow = ilowt (jlow)
                  xmax = xdont (ilow)
                  do icrs = 1, jlow
                     if (xdont(ilowt(icrs)) > xmax) then
                        iwrk = ilowt (icrs)
                        xmax = xdont (iwrk)
                        ilowt (icrs) = ilow
                        ilow = iwrk
                     end if
                  end do
                  jlow = jlow - 1
               end if

           end if

           jlm2 = jlm1
           jlm1 = jlow
           jhm2 = jhm1
           jhm1 = jhig
           select case (nord-jlow)
           case (2:)

              select case (jhig)
              case (2)
                 if (xdont(ihigt(1)) <= xdont(ihigt(2))) then
                    jlow = jlow + 1
                    ilowt (jlow) = ihigt (1)
                    jlow = jlow + 1
                    ilowt (jlow) = ihigt (2)
                 else
                    jlow = jlow + 1
                    ilowt (jlow) = ihigt (2)
                    jlow = jlow + 1
                    ilowt (jlow) = ihigt (1)
                 end if
                 exit

              case (3)

                 iwrk1 = ihigt (1)
                 iwrk2 = ihigt (2)
                 iwrk3 = ihigt (3)
                 if (xdont(iwrk2) < xdont(iwrk1)) then
                    ihigt (1) = iwrk2
                    ihigt (2) = iwrk1
                    iwrk2 = iwrk1
                 end if

                 if (xdont(iwrk2) > xdont(iwrk3)) then
                    ihigt (3) = iwrk2
                    ihigt (2) = iwrk3
                    iwrk2 = iwrk3
                    if (xdont(iwrk2) < xdont(ihigt(1))) then
                       ihigt (2) = ihigt (1)
                       ihigt (1) = iwrk2
                    end if
                 end if

                 jhig = 0
                 do icrs = jlow + 1, nord
                    jhig = jhig + 1
                    ilowt (icrs) = ihigt (jhig)
                 end do
                 jlow = nord
                 exit

              case (4:)
                 xpiv0 = xpiv
                 ifin = jhig
                 iwrk1 = ihigt (1)
                 iwrk2 = ihigt (2)
                 iwrk3 = ihigt (ifin)
                 if (xdont(iwrk2) < xdont(iwrk1)) then
                    ihigt (1) = iwrk2
                    ihigt (2) = iwrk1
                    iwrk2 = iwrk1
                 end if

                 if (xdont(iwrk2) > xdont(iwrk3)) then
                    ihigt (ifin) = iwrk2
                    ihigt (2) = iwrk3
                    iwrk2 = iwrk3
                    if (xdont(iwrk2) < xdont(ihigt(1))) then
                       ihigt (2) = ihigt (1)
                       ihigt (1) = iwrk2
                    end if
                 end if

                 jdeb = jlow
                 nwrk = nord - jlow
                 iwrk1 = ihigt (1)
                 jlow = jlow + 1
                 ilowt (jlow) = iwrk1
                 xpiv = xdont (iwrk1) + real (nwrk) / real (nord+nwrk) * (xdont(ihigt(ifin))-xdont(iwrk1))
                 jhig = 0
                 do icrs = 2, ifin
                    if (xdont(ihigt(icrs)) <= xpiv) then
                       jlow = jlow + 1
                       ilowt (jlow) = ihigt (icrs)
                       if (jlow >= nord) exit
                    else
                       jhig = jhig + 1
                       ihigt (jhig) = ihigt (icrs)
                    end if
                 end do

                 do icrs = icrs + 1, ifin
                    if (xdont(ihigt(icrs)) <= xpiv) then
                       jlow = jlow + 1
                       ilowt (jlow) = ihigt (icrs)
                    end if
                 end do
             end select

           case (1)
              xmin = xdont (ihigt(1))
              ihig = 1
              do icrs = 2, jhig
                 if (xdont(ihigt(icrs)) < xmin) then
                    xmin = xdont (ihigt(icrs))
                    ihig = icrs
                 end if
              end do
              jlow = jlow + 1
              ilowt (jlow) = ihigt (ihig)
              exit

           case (0)
              exit

           case (-5:-1)
              irngt (1) = ilowt (1)
              do icrs = 2, nord
                 iwrk = ilowt (icrs)
                 xwrk = xdont (iwrk)
                 do idcr = icrs - 1, 1, - 1
                    if (xwrk < xdont(irngt(idcr))) then
                       irngt (idcr+1) = irngt (idcr)
                    else
                       exit
                    end if
                 end do
                 irngt (idcr+1) = iwrk
              end do

              xwrk1 = xdont (irngt(nord))
              do icrs = nord + 1, jlow
                 if (xdont(ilowt (icrs)) < xwrk1) then
                    xwrk = xdont (ilowt (icrs))
                    do idcr = nord - 1, 1, - 1
                       if (xwrk >= xdont(irngt(idcr))) exit
                       irngt (idcr+1) = irngt (idcr)
                    end do

                    irngt (idcr+1) = ilowt (icrs)
                    xwrk1 = xdont (irngt(nord))
                 end if
              end do
              return

           case (:-6)
              ideb = jdeb + 1
              imil = (jlow+ideb) / 2
              ifin = jlow
              if (xdont(ilowt(imil)) < xdont(ilowt(ideb))) then
                 iwrk = ilowt (ideb)
                 ilowt (ideb) = ilowt (imil)
                 ilowt (imil) = iwrk
              end if

              if (xdont(ilowt(imil)) > xdont(ilowt(ifin))) then
                 iwrk = ilowt (ifin)
                 ilowt (ifin) = ilowt (imil)
                 ilowt (imil) = iwrk
                 if (xdont(ilowt(imil)) < xdont(ilowt(ideb))) then
                    iwrk = ilowt (ideb)
                    ilowt (ideb) = ilowt (imil)
                    ilowt (imil) = iwrk
                 end if

              end if

              if (ifin <= 3) exit
              xpiv = xdont (ilowt(1)) + real(nord)/real(jlow+nord) * &
                                        (xdont(ilowt(ifin))-xdont(ilowt(1)))
              if (jdeb > 0) then
                 if (xpiv <= xpiv0) &
                     xpiv = xpiv0 + real(2*nord-jdeb)/real (jlow+nord) * &
                                    (xdont(ilowt(ifin))-xpiv0)
              else
                 ideb = 1
              end if

              jhig = 0
              jlow = jdeb
              if (xdont(ilowt(ifin)) > xpiv) then
                 icrs = jdeb
                 do
                   icrs = icrs + 1
                    if (xdont(ilowt(icrs)) > xpiv) then
                       jhig = jhig + 1
                       ihigt (jhig) = ilowt (icrs)
                       if (icrs >= ifin) exit
                    else
                       jlow = jlow + 1
                       ilowt (jlow) = ilowt (icrs)
                       if (jlow >= nord) exit
                    end if
                 end do

                 if (icrs < ifin) then
                    do
                       icrs = icrs + 1
                       if (xdont(ilowt(icrs)) <= xpiv) then
                          jlow = jlow + 1
                          ilowt (jlow) = ilowt (icrs)
                       else
                          if (icrs >= ifin) exit
                       end if
                    end do
                 end if

             else

                 do icrs = ideb, ifin
                    if (xdont(ilowt(icrs)) > xpiv) then
                       jhig = jhig + 1
                       ihigt (jhig) = ilowt (icrs)
                    else
                       jlow = jlow + 1
                       ilowt (jlow) = ilowt (icrs)
                       if (jlow >= nord) exit
                    end if
                 end do

                 do icrs = icrs + 1, ifin
                    if (xdont(ilowt(icrs)) <= xpiv) then
                       jlow = jlow + 1
                       ilowt (jlow) = ilowt (icrs)
                    end if
                 end do

              end if
           end select
        end do

        irngt (1) = ilowt (1)
        do icrs = 2, nord
           iwrk = ilowt (icrs)
           xwrk = xdont (iwrk)
           do idcr = icrs - 1, 1, - 1
              if (xwrk < xdont(irngt(idcr))) then
                 irngt (idcr+1) = irngt (idcr)
              else
                 exit
              end if
           end do
           irngt (idcr+1) = iwrk
        end do

       return
    end subroutine i_rnkpar


!------------------------------------------------------------------------
!              PARTIAL SORTING ARRAY                                    |
!------------------------------------------------------------------------

    subroutine psort_real8(xin, points, xout)
        integer, intent(in) :: points
        real(RK), dimension(:), intent(in)  :: xin
        real(RK), dimension(:), intent(out) :: xout
        integer, dimension(size(xin)) :: inds
        
        call D_rnkpar(xin, inds, points)
        xout = xin(inds)
    end subroutine psort_real8
    subroutine psort_real4(xin, points, xout)
        integer, intent(in) :: points
        real, dimension(:), intent(in)  :: xin
        real, dimension(:), intent(out) :: xout
        integer, dimension(size(xin)) :: inds
        
        call R_rnkpar(xin, inds, points)
        xout = xin(inds)
    end subroutine psort_real4
    subroutine psort_integer(xin, points, xout)
        integer, intent(in) :: points
        integer, dimension(:), intent(in)  :: xin
        integer, dimension(:), intent(out) :: xout
        integer, dimension(size(xin)) :: inds
        
        call I_rnkpar(xin, inds, points)
        xout = xin(inds)
    end subroutine psort_integer

 !------------------------------------------------------------------------
 !               MEDIAN ARRAYS                                       |
 !------------------------------------------------------------------------
 !--------------------------------------------------------------
 !  Finds the median of XDONT using the recursive procedure
 !  described in Knuth, The Art of Computer Programming,
 !  vol. 3, 5.3.3 - This procedure is linear in time, and
 !  does not require to be able to interpolate in the
 !  set as the one used in INDNTH. It also has better worst
 !  case behavior than INDNTH, but is about 30% slower in
 !  average for random uniformly distributed values.
 !---------------------------------------------------------------  

    Recursive Function D_valmed (XDONT) Result (res_med)

! __________________________________________________________
! __________________________________________________________
      Real (RK), Dimension (:), Intent (In) :: XDONT
      Real (RK) :: res_med
! __________________________________________________________
      Real (RK), Parameter :: XHUGE = HUGE (XDONT)
      Real (RK), Dimension (SIZE(XDONT)+6) :: XWRKT
      Real (RK) :: XWRK, XWRK1, XMED7
!
      Integer, Dimension ((SIZE(XDONT)+6)/7) :: ISTRT, IENDT, IMEDT
      Integer :: NDON, NTRI, NMED, NORD, NEQU, NLEQ, IMED, IDON, IDON1
      Integer :: IDEB, IWRK, IDCR, ICRS, ICRS1, ICRS2, IMED1
!
      NDON = SIZE (XDONT)
      NMED = (NDON+1) / 2
!      write(unit=*,fmt=*) NMED, NDON
!
!  If the number of values is small, then use insertion sort
!
      If (NDON < 35) Then
!
!  Bring minimum to first location to save test in decreasing loop
!
         IDCR = NDON
         If (XDONT (1) < XDONT (NDON)) Then
            XWRK = XDONT (1)
            XWRKT (IDCR) = XDONT (IDCR)
         Else
            XWRK = XDONT (IDCR)
            XWRKT (IDCR) = XDONT (1)
         Endif
         Do IWRK = 1, NDON - 2
            IDCR = IDCR - 1
            XWRK1 = XDONT (IDCR)
            If (XWRK1 < XWRK) Then
                XWRKT (IDCR) = XWRK
                XWRK = XWRK1
            Else
                XWRKT (IDCR) = XWRK1
            Endif
         End Do
         XWRKT (1) = XWRK
!
! Sort the first half, until we have NMED sorted values
!
         Do ICRS = 3, NMED
            XWRK = XWRKT (ICRS)
               IDCR = ICRS - 1
               Do
                  If (XWRK >= XWRKT(IDCR)) Exit
                  XWRKT (IDCR+1) = XWRKT (IDCR)
                  IDCR = IDCR - 1
               End Do
            XWRKT (IDCR+1) = XWRK
         End Do
!
!  Insert any value less than the current median in the first half
!
         Do ICRS = NMED+1, NDON
            XWRK = XWRKT (ICRS)
            If (XWRK < XWRKT (NMED)) Then
               IDCR = NMED - 1
               Do
                  If (XWRK >= XWRKT(IDCR)) Exit
                  XWRKT (IDCR+1) = XWRKT (IDCR)
                  IDCR = IDCR - 1
               End Do
               XWRKT (IDCR+1) = XWRK
            End If
         End Do
         res_med = XWRKT (NMED)
         Return
      End If
!
!  Make sorted subsets of 7 elements
!  This is done by a variant of insertion sort where a first
!  pass is used to bring the smallest element to the first position
!  decreasing disorder at the same time, so that we may remove
!  remove the loop test in the insertion loop.
!
      DO IDEB = 1, NDON-6, 7
         IDCR = IDEB + 6
         If (XDONT (IDEB) < XDONT (IDCR)) Then
            XWRK = XDONT (IDEB)
            XWRKT (IDCR) = XDONT (IDCR)
         Else
            XWRK = XDONT (IDCR)
            XWRKT (IDCR) = XDONT (IDEB)
         Endif
         Do IWRK = 1, 5
            IDCR = IDCR - 1
            XWRK1 = XDONT (IDCR)
            If (XWRK1 < XWRK) Then
                XWRKT (IDCR) = XWRK
                XWRK = XWRK1
            Else
                XWRKT (IDCR) = XWRK1
            Endif
         End Do
         XWRKT (IDEB) = XWRK
         Do ICRS = IDEB+2, IDEB+6
            XWRK = XWRKT (ICRS)
            If (XWRK < XWRKT(ICRS-1)) Then
               XWRKT (ICRS) = XWRKT (ICRS-1)
               IDCR = ICRS - 1
               XWRK1 = XWRKT (IDCR-1)
               Do
                  If (XWRK >= XWRK1) Exit
                  XWRKT (IDCR) = XWRK1
                  IDCR = IDCR - 1
                  XWRK1 = XWRKT (IDCR-1)
               End Do
               XWRKT (IDCR) = XWRK
            EndIf
         End Do
      End Do
!
!  Add-up alternatively + and - HUGE values to make the number of data
!  an exact multiple of 7.
!
      IDEB = 7 * (NDON/7)
      NTRI = NDON
      If (IDEB < NDON) Then
!
         XWRK1 = XHUGE
         Do ICRS = IDEB+1, IDEB+7
            If (ICRS <= NDON) Then
               XWRKT (ICRS) = XDONT (ICRS)
            Else
               If (XWRK1 /= XHUGE) NMED = NMED + 1
               XWRKT (ICRS) = XWRK1
               XWRK1 = - XWRK1
            Endif
         End Do
!
         Do ICRS = IDEB+2, IDEB+7
            XWRK = XWRKT (ICRS)
            Do IDCR = ICRS - 1, IDEB+1, - 1
               If (XWRK >= XWRKT(IDCR)) Exit
               XWRKT (IDCR+1) = XWRKT (IDCR)
            End Do
            XWRKT (IDCR+1) = XWRK
         End Do
!
         NTRI = IDEB+7
      End If
!
!  Make the set of the indices of median values of each sorted subset
!
         IDON1 = 0
         Do IDON = 1, NTRI, 7
            IDON1 = IDON1 + 1
            IMEDT (IDON1) = IDON + 3
         End Do
!
!  Find XMED7, the median of the medians
!
         XMED7 = D_valmed (XWRKT (IMEDT))
!
!  Count how many values are not higher than (and how many equal to) XMED7
!  This number is at least 4 * 1/2 * (N/7) : 4 values in each of the
!  subsets where the median is lower than the median of medians. For similar
!  reasons, we also have at least 2N/7 values not lower than XMED7. At the
!  same time, we find in each subset the index of the last value < XMED7,
!  and that of the first > XMED7. These indices will be used to restrict the
!  search for the median as the Kth element in the subset (> or <) where
!  we know it to be.
!
         IDON1 = 1
         NLEQ = 0
         NEQU = 0
         Do IDON = 1, NTRI, 7
            IMED = IDON+3
            If (XWRKT (IMED) > XMED7) Then
                  IMED = IMED - 2
                  If (XWRKT (IMED) > XMED7) Then
                     IMED = IMED - 1
                  Else If (XWRKT (IMED) < XMED7) Then
                     IMED = IMED + 1
                  Endif
            Else If (XWRKT (IMED) < XMED7) Then
                  IMED = IMED + 2
                  If (XWRKT (IMED) > XMED7) Then
                     IMED = IMED - 1
                  Else If (XWRKT (IMED) < XMED7) Then
                     IMED = IMED + 1
                  Endif
            Endif
            If (XWRKT (IMED) > XMED7) Then
               NLEQ = NLEQ + IMED - IDON
               IENDT (IDON1) = IMED - 1
               ISTRT (IDON1) = IMED
            Else If (XWRKT (IMED) < XMED7) Then
               NLEQ = NLEQ + IMED - IDON + 1
               IENDT (IDON1) = IMED
               ISTRT (IDON1) = IMED + 1
            Else                    !       If (XWRKT (IMED) == XMED7)
               NLEQ = NLEQ + IMED - IDON + 1
               NEQU = NEQU + 1
               IENDT (IDON1) = IMED - 1
               Do IMED1 = IMED - 1, IDON, -1
                  If (XWRKT (IMED1) == XMED7) Then
                     NEQU = NEQU + 1
                     IENDT (IDON1) = IMED1 - 1
                  Else
                     Exit
                  End If
               End Do
               ISTRT (IDON1) = IMED + 1
               Do IMED1 = IMED + 1, IDON + 6
                  If (XWRKT (IMED1) == XMED7) Then
                     NEQU = NEQU + 1
                     NLEQ = NLEQ + 1
                     ISTRT (IDON1) = IMED1 + 1
                  Else
                     Exit
                  End If
               End Do
            Endif
            IDON1 = IDON1 + 1
         End Do
!
!  Carry out a partial insertion sort to find the Kth smallest of the
!  large values, or the Kth largest of the small values, according to
!  what is needed.
!
        If (NLEQ - NEQU + 1 <= NMED) Then
            If (NLEQ < NMED) Then   !      Not enough low values
                XWRK1 = XHUGE
                NORD = NMED - NLEQ
                IDON1 = 0
                ICRS1 = 1
                ICRS2 = 0
                IDCR = 0
               Do IDON = 1, NTRI, 7
                   IDON1 = IDON1 + 1
                   If (ICRS2 < NORD) Then
                      Do ICRS = ISTRT (IDON1), IDON + 6
                         If (XWRKT(ICRS) < XWRK1) Then
                            XWRK = XWRKT (ICRS)
                            Do IDCR = ICRS1 - 1, 1, - 1
                               If (XWRK >= XWRKT(IDCR)) Exit
                               XWRKT (IDCR+1) = XWRKT (IDCR)
                            End Do
                            XWRKT (IDCR+1) = XWRK
                            XWRK1 = XWRKT(ICRS1)
                         Else
                           If (ICRS2 < NORD) Then
                              XWRKT (ICRS1) = XWRKT (ICRS)
                              XWRK1 = XWRKT(ICRS1)
                           Endif
                         End If
                         ICRS1 = MIN (NORD, ICRS1 + 1)
                         ICRS2 = MIN (NORD, ICRS2 + 1)
                      End Do
                   Else
                      Do ICRS = ISTRT (IDON1), IDON + 6
                         If (XWRKT(ICRS) >= XWRK1) Exit
                         XWRK = XWRKT (ICRS)
                         Do IDCR = ICRS1 - 1, 1, - 1
                               If (XWRK >= XWRKT(IDCR)) Exit
                               XWRKT (IDCR+1) = XWRKT (IDCR)
                         End Do
                         XWRKT (IDCR+1) = XWRK
                         XWRK1 = XWRKT(ICRS1)
                      End Do
                   End If
                End Do
                res_med = XWRK1
                Return
            Else
                res_med = XMED7
                Return
            End If
         Else                       !      If (NLEQ > NMED)
!                                          Not enough high values
                XWRK1 = -XHUGE
                NORD = NLEQ - NEQU - NMED + 1
                IDON1 = 0
                ICRS1 = 1
                ICRS2 = 0
                Do IDON = 1, NTRI, 7
                   IDON1 = IDON1 + 1
                   If (ICRS2 < NORD) Then
!
                      Do ICRS = IDON, IENDT (IDON1)
                         If (XWRKT(ICRS) > XWRK1) Then
                            XWRK = XWRKT (ICRS)
                            IDCR = ICRS1 - 1
                            Do IDCR = ICRS1 - 1, 1, - 1
                               If (XWRK <= XWRKT(IDCR)) Exit
                               XWRKT (IDCR+1) = XWRKT (IDCR)
                            End Do
                            XWRKT (IDCR+1) = XWRK
                            XWRK1 = XWRKT(ICRS1)
                         Else
                            If (ICRS2 < NORD) Then
                               XWRKT (ICRS1) = XWRKT (ICRS)
                               XWRK1 = XWRKT (ICRS1)
                            End If
                         End If
                         ICRS1 = MIN (NORD, ICRS1 + 1)
                         ICRS2 = MIN (NORD, ICRS2 + 1)
                      End Do
                   Else
                      Do ICRS = IENDT (IDON1), IDON, -1
                         If (XWRKT(ICRS) > XWRK1) Then
                            XWRK = XWRKT (ICRS)
                            IDCR = ICRS1 - 1
                            Do IDCR = ICRS1 - 1, 1, - 1
                               If (XWRK <= XWRKT(IDCR)) Exit
                               XWRKT (IDCR+1) = XWRKT (IDCR)
                            End Do
                            XWRKT (IDCR+1) = XWRK
                            XWRK1 = XWRKT(ICRS1)
                         Else
                            Exit
                         End If
                      End Do
                   Endif
                End Do
!
                res_med = XWRK1
                Return
         End If
!
End Function D_valmed

Recursive Function R_valmed (XDONT) Result (res_med)
!  Finds the median of XDONT using the recursive procedure
!  described in Knuth, The Art of Computer Programming,
!  vol. 3, 5.3.3 - This procedure is linear in time, and
!  does not require to be able to interpolate in the
!  set as the one used in INDNTH. It also has better worst
!  case behavior than INDNTH, but is about 30% slower in
!  average for random uniformly distributed values.
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (In) :: XDONT
      Real :: res_med
! __________________________________________________________
      Real, Parameter :: XHUGE = HUGE (XDONT)
      Real, Dimension (SIZE(XDONT)+6) :: XWRKT
      Real :: XWRK, XWRK1, XMED7
!
      Integer, Dimension ((SIZE(XDONT)+6)/7) :: ISTRT, IENDT, IMEDT
      Integer :: NDON, NTRI, NMED, NORD, NEQU, NLEQ, IMED, IDON, IDON1
      Integer :: IDEB, IWRK, IDCR, ICRS, ICRS1, ICRS2, IMED1
!
      NDON = SIZE (XDONT)
      NMED = (NDON+1) / 2
!      write(unit=*,fmt=*) NMED, NDON
!
!  If the number of values is small, then use insertion sort
!
      If (NDON < 35) Then
!
!  Bring minimum to first location to save test in decreasing loop
!
         IDCR = NDON
         If (XDONT (1) < XDONT (NDON)) Then
            XWRK = XDONT (1)
            XWRKT (IDCR) = XDONT (IDCR)
         Else
            XWRK = XDONT (IDCR)
            XWRKT (IDCR) = XDONT (1)
         Endif
         Do IWRK = 1, NDON - 2
            IDCR = IDCR - 1
            XWRK1 = XDONT (IDCR)
            If (XWRK1 < XWRK) Then
                XWRKT (IDCR) = XWRK
                XWRK = XWRK1
            Else
                XWRKT (IDCR) = XWRK1
            Endif
         End Do
         XWRKT (1) = XWRK
!
! Sort the first half, until we have NMED sorted values
!
         Do ICRS = 3, NMED
            XWRK = XWRKT (ICRS)
               IDCR = ICRS - 1
               Do
                  If (XWRK >= XWRKT(IDCR)) Exit
                  XWRKT (IDCR+1) = XWRKT (IDCR)
                  IDCR = IDCR - 1
               End Do
            XWRKT (IDCR+1) = XWRK
         End Do
!
!  Insert any value less than the current median in the first half
!
         Do ICRS = NMED+1, NDON
            XWRK = XWRKT (ICRS)
            If (XWRK < XWRKT (NMED)) Then
               IDCR = NMED - 1
               Do
                  If (XWRK >= XWRKT(IDCR)) Exit
                  XWRKT (IDCR+1) = XWRKT (IDCR)
                  IDCR = IDCR - 1
               End Do
               XWRKT (IDCR+1) = XWRK
            End If
         End Do
         res_med = XWRKT (NMED)
         Return
      End If
!
!  Make sorted subsets of 7 elements
!  This is done by a variant of insertion sort where a first
!  pass is used to bring the smallest element to the first position
!  decreasing disorder at the same time, so that we may remove
!  remove the loop test in the insertion loop.
!
      DO IDEB = 1, NDON-6, 7
         IDCR = IDEB + 6
         If (XDONT (IDEB) < XDONT (IDCR)) Then
            XWRK = XDONT (IDEB)
            XWRKT (IDCR) = XDONT (IDCR)
         Else
            XWRK = XDONT (IDCR)
            XWRKT (IDCR) = XDONT (IDEB)
         Endif
         Do IWRK = 1, 5
            IDCR = IDCR - 1
            XWRK1 = XDONT (IDCR)
            If (XWRK1 < XWRK) Then
                XWRKT (IDCR) = XWRK
                XWRK = XWRK1
            Else
                XWRKT (IDCR) = XWRK1
            Endif
         End Do
         XWRKT (IDEB) = XWRK
         Do ICRS = IDEB+2, IDEB+6
            XWRK = XWRKT (ICRS)
            If (XWRK < XWRKT(ICRS-1)) Then
               XWRKT (ICRS) = XWRKT (ICRS-1)
               IDCR = ICRS - 1
               XWRK1 = XWRKT (IDCR-1)
               Do
                  If (XWRK >= XWRK1) Exit
                  XWRKT (IDCR) = XWRK1
                  IDCR = IDCR - 1
                  XWRK1 = XWRKT (IDCR-1)
               End Do
               XWRKT (IDCR) = XWRK
            EndIf
         End Do
      End Do
!
!  Add-up alternatively + and - HUGE values to make the number of data
!  an exact multiple of 7.
!
      IDEB = 7 * (NDON/7)
      NTRI = NDON
      If (IDEB < NDON) Then
!
         XWRK1 = XHUGE
         Do ICRS = IDEB+1, IDEB+7
            If (ICRS <= NDON) Then
               XWRKT (ICRS) = XDONT (ICRS)
            Else
               If (XWRK1 /= XHUGE) NMED = NMED + 1
               XWRKT (ICRS) = XWRK1
               XWRK1 = - XWRK1
            Endif
         End Do
!
         Do ICRS = IDEB+2, IDEB+7
            XWRK = XWRKT (ICRS)
            Do IDCR = ICRS - 1, IDEB+1, - 1
               If (XWRK >= XWRKT(IDCR)) Exit
               XWRKT (IDCR+1) = XWRKT (IDCR)
            End Do
            XWRKT (IDCR+1) = XWRK
         End Do
!
         NTRI = IDEB+7
      End If
!
!  Make the set of the indices of median values of each sorted subset
!
         IDON1 = 0
         Do IDON = 1, NTRI, 7
            IDON1 = IDON1 + 1
            IMEDT (IDON1) = IDON + 3
         End Do
!
!  Find XMED7, the median of the medians
!
         XMED7 = R_valmed (XWRKT (IMEDT))
!
!  Count how many values are not higher than (and how many equal to) XMED7
!  This number is at least 4 * 1/2 * (N/7) : 4 values in each of the
!  subsets where the median is lower than the median of medians. For similar
!  reasons, we also have at least 2N/7 values not lower than XMED7. At the
!  same time, we find in each subset the index of the last value < XMED7,
!  and that of the first > XMED7. These indices will be used to restrict the
!  search for the median as the Kth element in the subset (> or <) where
!  we know it to be.
!
         IDON1 = 1
         NLEQ = 0
         NEQU = 0
         Do IDON = 1, NTRI, 7
            IMED = IDON+3
            If (XWRKT (IMED) > XMED7) Then
                  IMED = IMED - 2
                  If (XWRKT (IMED) > XMED7) Then
                     IMED = IMED - 1
                  Else If (XWRKT (IMED) < XMED7) Then
                     IMED = IMED + 1
                  Endif
            Else If (XWRKT (IMED) < XMED7) Then
                  IMED = IMED + 2
                  If (XWRKT (IMED) > XMED7) Then
                     IMED = IMED - 1
                  Else If (XWRKT (IMED) < XMED7) Then
                     IMED = IMED + 1
                  Endif
            Endif
            If (XWRKT (IMED) > XMED7) Then
               NLEQ = NLEQ + IMED - IDON
               IENDT (IDON1) = IMED - 1
               ISTRT (IDON1) = IMED
            Else If (XWRKT (IMED) < XMED7) Then
               NLEQ = NLEQ + IMED - IDON + 1
               IENDT (IDON1) = IMED
               ISTRT (IDON1) = IMED + 1
            Else                    !       If (XWRKT (IMED) == XMED7)
               NLEQ = NLEQ + IMED - IDON + 1
               NEQU = NEQU + 1
               IENDT (IDON1) = IMED - 1
               Do IMED1 = IMED - 1, IDON, -1
                  If (XWRKT (IMED1) == XMED7) Then
                     NEQU = NEQU + 1
                     IENDT (IDON1) = IMED1 - 1
                  Else
                     Exit
                  End If
               End Do
               ISTRT (IDON1) = IMED + 1
               Do IMED1 = IMED + 1, IDON + 6
                  If (XWRKT (IMED1) == XMED7) Then
                     NEQU = NEQU + 1
                     NLEQ = NLEQ + 1
                     ISTRT (IDON1) = IMED1 + 1
                  Else
                     Exit
                  End If
               End Do
            Endif
            IDON1 = IDON1 + 1
         End Do
!
!  Carry out a partial insertion sort to find the Kth smallest of the
!  large values, or the Kth largest of the small values, according to
!  what is needed.
!
        If (NLEQ - NEQU + 1 <= NMED) Then
            If (NLEQ < NMED) Then   !      Not enough low values
                XWRK1 = XHUGE
                NORD = NMED - NLEQ
                IDON1 = 0
                ICRS1 = 1
                ICRS2 = 0
                IDCR = 0
               Do IDON = 1, NTRI, 7
                   IDON1 = IDON1 + 1
                   If (ICRS2 < NORD) Then
                      Do ICRS = ISTRT (IDON1), IDON + 6
                         If (XWRKT(ICRS) < XWRK1) Then
                            XWRK = XWRKT (ICRS)
                            Do IDCR = ICRS1 - 1, 1, - 1
                               If (XWRK >= XWRKT(IDCR)) Exit
                               XWRKT (IDCR+1) = XWRKT (IDCR)
                            End Do
                            XWRKT (IDCR+1) = XWRK
                            XWRK1 = XWRKT(ICRS1)
                         Else
                           If (ICRS2 < NORD) Then
                              XWRKT (ICRS1) = XWRKT (ICRS)
                              XWRK1 = XWRKT(ICRS1)
                           Endif
                         End If
                         ICRS1 = MIN (NORD, ICRS1 + 1)
                         ICRS2 = MIN (NORD, ICRS2 + 1)
                      End Do
                   Else
                      Do ICRS = ISTRT (IDON1), IDON + 6
                         If (XWRKT(ICRS) >= XWRK1) Exit
                         XWRK = XWRKT (ICRS)
                         Do IDCR = ICRS1 - 1, 1, - 1
                               If (XWRK >= XWRKT(IDCR)) Exit
                               XWRKT (IDCR+1) = XWRKT (IDCR)
                         End Do
                         XWRKT (IDCR+1) = XWRK
                         XWRK1 = XWRKT(ICRS1)
                      End Do
                   End If
                End Do
                res_med = XWRK1
                Return
            Else
                res_med = XMED7
                Return
            End If
         Else                       !      If (NLEQ > NMED)
!                                          Not enough high values
                XWRK1 = -XHUGE
                NORD = NLEQ - NEQU - NMED + 1
                IDON1 = 0
                ICRS1 = 1
                ICRS2 = 0
                Do IDON = 1, NTRI, 7
                   IDON1 = IDON1 + 1
                   If (ICRS2 < NORD) Then
!
                      Do ICRS = IDON, IENDT (IDON1)
                         If (XWRKT(ICRS) > XWRK1) Then
                            XWRK = XWRKT (ICRS)
                            IDCR = ICRS1 - 1
                            Do IDCR = ICRS1 - 1, 1, - 1
                               If (XWRK <= XWRKT(IDCR)) Exit
                               XWRKT (IDCR+1) = XWRKT (IDCR)
                            End Do
                            XWRKT (IDCR+1) = XWRK
                            XWRK1 = XWRKT(ICRS1)
                         Else
                            If (ICRS2 < NORD) Then
                               XWRKT (ICRS1) = XWRKT (ICRS)
                               XWRK1 = XWRKT (ICRS1)
                            End If
                         End If
                         ICRS1 = MIN (NORD, ICRS1 + 1)
                         ICRS2 = MIN (NORD, ICRS2 + 1)
                      End Do
                   Else
                      Do ICRS = IENDT (IDON1), IDON, -1
                         If (XWRKT(ICRS) > XWRK1) Then
                            XWRK = XWRKT (ICRS)
                            IDCR = ICRS1 - 1
                            Do IDCR = ICRS1 - 1, 1, - 1
                               If (XWRK <= XWRKT(IDCR)) Exit
                               XWRKT (IDCR+1) = XWRKT (IDCR)
                            End Do
                            XWRKT (IDCR+1) = XWRK
                            XWRK1 = XWRKT(ICRS1)
                         Else
                            Exit
                         End If
                      End Do
                   Endif
                End Do
!
                res_med = XWRK1
                Return
         End If
!
End Function R_valmed
Recursive Function I_valmed (XDONT) Result (res_med)
!  Finds the median of XDONT using the recursive procedure
!  described in Knuth, The Art of Computer Programming,
!  vol. 3, 5.3.3 - This procedure is linear in time, and
!  does not require to be able to interpolate in the
!  set as the one used in INDNTH. It also has better worst
!  case behavior than INDNTH, but is about 30% slower in
!  average for random uniformly distributed values.
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (In)  :: XDONT
      Integer :: res_med
! __________________________________________________________
      Integer, Parameter :: XHUGE = HUGE (XDONT)
      Integer, Dimension (SIZE(XDONT)+6) :: XWRKT
      Integer :: XWRK, XWRK1, XMED7
!
      Integer, Dimension ((SIZE(XDONT)+6)/7) :: ISTRT, IENDT, IMEDT
      Integer :: NDON, NTRI, NMED, NORD, NEQU, NLEQ, IMED, IDON, IDON1
      Integer :: IDEB, IWRK, IDCR, ICRS, ICRS1, ICRS2, IMED1
!
      NDON = SIZE (XDONT)
      NMED = (NDON+1) / 2
!      write(unit=*,fmt=*) NMED, NDON
!
!  If the number of values is small, then use insertion sort
!
      If (NDON < 35) Then
!
!  Bring minimum to first location to save test in decreasing loop
!
         IDCR = NDON
         If (XDONT (1) < XDONT (NDON)) Then
            XWRK = XDONT (1)
            XWRKT (IDCR) = XDONT (IDCR)
         Else
            XWRK = XDONT (IDCR)
            XWRKT (IDCR) = XDONT (1)
         Endif
         Do IWRK = 1, NDON - 2
            IDCR = IDCR - 1
            XWRK1 = XDONT (IDCR)
            If (XWRK1 < XWRK) Then
                XWRKT (IDCR) = XWRK
                XWRK = XWRK1
            Else
                XWRKT (IDCR) = XWRK1
            Endif
         End Do
         XWRKT (1) = XWRK
!
! Sort the first half, until we have NMED sorted values
!
         Do ICRS = 3, NMED
            XWRK = XWRKT (ICRS)
               IDCR = ICRS - 1
               Do
                  If (XWRK >= XWRKT(IDCR)) Exit
                  XWRKT (IDCR+1) = XWRKT (IDCR)
                  IDCR = IDCR - 1
               End Do
            XWRKT (IDCR+1) = XWRK
         End Do
!
!  Insert any value less than the current median in the first half
!
         Do ICRS = NMED+1, NDON
            XWRK = XWRKT (ICRS)
            If (XWRK < XWRKT (NMED)) Then
               IDCR = NMED - 1
               Do
                  If (XWRK >= XWRKT(IDCR)) Exit
                  XWRKT (IDCR+1) = XWRKT (IDCR)
                  IDCR = IDCR - 1
               End Do
               XWRKT (IDCR+1) = XWRK
            End If
         End Do
         res_med = XWRKT (NMED)
         Return
      End If
!
!  Make sorted subsets of 7 elements
!  This is done by a variant of insertion sort where a first
!  pass is used to bring the smallest element to the first position
!  decreasing disorder at the same time, so that we may remove
!  remove the loop test in the insertion loop.
!
      DO IDEB = 1, NDON-6, 7
         IDCR = IDEB + 6
         If (XDONT (IDEB) < XDONT (IDCR)) Then
            XWRK = XDONT (IDEB)
            XWRKT (IDCR) = XDONT (IDCR)
         Else
            XWRK = XDONT (IDCR)
            XWRKT (IDCR) = XDONT (IDEB)
         Endif
         Do IWRK = 1, 5
            IDCR = IDCR - 1
            XWRK1 = XDONT (IDCR)
            If (XWRK1 < XWRK) Then
                XWRKT (IDCR) = XWRK
                XWRK = XWRK1
            Else
                XWRKT (IDCR) = XWRK1
            Endif
         End Do
         XWRKT (IDEB) = XWRK
         Do ICRS = IDEB+2, IDEB+6
            XWRK = XWRKT (ICRS)
            If (XWRK < XWRKT(ICRS-1)) Then
               XWRKT (ICRS) = XWRKT (ICRS-1)
               IDCR = ICRS - 1
               XWRK1 = XWRKT (IDCR-1)
               Do
                  If (XWRK >= XWRK1) Exit
                  XWRKT (IDCR) = XWRK1
                  IDCR = IDCR - 1
                  XWRK1 = XWRKT (IDCR-1)
               End Do
               XWRKT (IDCR) = XWRK
            EndIf
         End Do
      End Do
!
!  Add-up alternatively + and - HUGE values to make the number of data
!  an exact multiple of 7.
!
      IDEB = 7 * (NDON/7)
      NTRI = NDON
      If (IDEB < NDON) Then
!
         XWRK1 = XHUGE
         Do ICRS = IDEB+1, IDEB+7
            If (ICRS <= NDON) Then
               XWRKT (ICRS) = XDONT (ICRS)
            Else
               If (XWRK1 /= XHUGE) NMED = NMED + 1
               XWRKT (ICRS) = XWRK1
               XWRK1 = - XWRK1
            Endif
         End Do
!
         Do ICRS = IDEB+2, IDEB+7
            XWRK = XWRKT (ICRS)
            Do IDCR = ICRS - 1, IDEB+1, - 1
               If (XWRK >= XWRKT(IDCR)) Exit
               XWRKT (IDCR+1) = XWRKT (IDCR)
            End Do
            XWRKT (IDCR+1) = XWRK
         End Do
!
         NTRI = IDEB+7
      End If
!
!  Make the set of the indices of median values of each sorted subset
!
         IDON1 = 0
         Do IDON = 1, NTRI, 7
            IDON1 = IDON1 + 1
            IMEDT (IDON1) = IDON + 3
         End Do
!
!  Find XMED7, the median of the medians
!
         XMED7 = I_valmed (XWRKT (IMEDT))
!
!  Count how many values are not higher than (and how many equal to) XMED7
!  This number is at least 4 * 1/2 * (N/7) : 4 values in each of the
!  subsets where the median is lower than the median of medians. For similar
!  reasons, we also have at least 2N/7 values not lower than XMED7. At the
!  same time, we find in each subset the index of the last value < XMED7,
!  and that of the first > XMED7. These indices will be used to restrict the
!  search for the median as the Kth element in the subset (> or <) where
!  we know it to be.
!
         IDON1 = 1
         NLEQ = 0
         NEQU = 0
         Do IDON = 1, NTRI, 7
            IMED = IDON+3
            If (XWRKT (IMED) > XMED7) Then
                  IMED = IMED - 2
                  If (XWRKT (IMED) > XMED7) Then
                     IMED = IMED - 1
                  Else If (XWRKT (IMED) < XMED7) Then
                     IMED = IMED + 1
                  Endif
            Else If (XWRKT (IMED) < XMED7) Then
                  IMED = IMED + 2
                  If (XWRKT (IMED) > XMED7) Then
                     IMED = IMED - 1
                  Else If (XWRKT (IMED) < XMED7) Then
                     IMED = IMED + 1
                  Endif
            Endif
            If (XWRKT (IMED) > XMED7) Then
               NLEQ = NLEQ + IMED - IDON
               IENDT (IDON1) = IMED - 1
               ISTRT (IDON1) = IMED
            Else If (XWRKT (IMED) < XMED7) Then
               NLEQ = NLEQ + IMED - IDON + 1
               IENDT (IDON1) = IMED
               ISTRT (IDON1) = IMED + 1
            Else                    !       If (XWRKT (IMED) == XMED7)
               NLEQ = NLEQ + IMED - IDON + 1
               NEQU = NEQU + 1
               IENDT (IDON1) = IMED - 1
               Do IMED1 = IMED - 1, IDON, -1
                  If (XWRKT (IMED1) == XMED7) Then
                     NEQU = NEQU + 1
                     IENDT (IDON1) = IMED1 - 1
                  Else
                     Exit
                  End If
               End Do
               ISTRT (IDON1) = IMED + 1
               Do IMED1 = IMED + 1, IDON + 6
                  If (XWRKT (IMED1) == XMED7) Then
                     NEQU = NEQU + 1
                     NLEQ = NLEQ + 1
                     ISTRT (IDON1) = IMED1 + 1
                  Else
                     Exit
                  End If
               End Do
            Endif
            IDON1 = IDON1 + 1
         End Do
!
!  Carry out a partial insertion sort to find the Kth smallest of the
!  large values, or the Kth largest of the small values, according to
!  what is needed.
!
        If (NLEQ - NEQU + 1 <= NMED) Then
            If (NLEQ < NMED) Then   !      Not enough low values
                XWRK1 = XHUGE
                NORD = NMED - NLEQ
                IDON1 = 0
                ICRS1 = 1
                ICRS2 = 0
                IDCR = 0
               Do IDON = 1, NTRI, 7
                   IDON1 = IDON1 + 1
                   If (ICRS2 < NORD) Then
                      Do ICRS = ISTRT (IDON1), IDON + 6
                         If (XWRKT(ICRS) < XWRK1) Then
                            XWRK = XWRKT (ICRS)
                            Do IDCR = ICRS1 - 1, 1, - 1
                               If (XWRK >= XWRKT(IDCR)) Exit
                               XWRKT (IDCR+1) = XWRKT (IDCR)
                            End Do
                            XWRKT (IDCR+1) = XWRK
                            XWRK1 = XWRKT(ICRS1)
                         Else
                           If (ICRS2 < NORD) Then
                              XWRKT (ICRS1) = XWRKT (ICRS)
                              XWRK1 = XWRKT(ICRS1)
                           Endif
                         End If
                         ICRS1 = MIN (NORD, ICRS1 + 1)
                         ICRS2 = MIN (NORD, ICRS2 + 1)
                      End Do
                   Else
                      Do ICRS = ISTRT (IDON1), IDON + 6
                         If (XWRKT(ICRS) >= XWRK1) Exit
                         XWRK = XWRKT (ICRS)
                         Do IDCR = ICRS1 - 1, 1, - 1
                               If (XWRK >= XWRKT(IDCR)) Exit
                               XWRKT (IDCR+1) = XWRKT (IDCR)
                         End Do
                         XWRKT (IDCR+1) = XWRK
                         XWRK1 = XWRKT(ICRS1)
                      End Do
                   End If
                End Do
                res_med = XWRK1
                Return
            Else
                res_med = XMED7
                Return
            End If
         Else                       !      If (NLEQ > NMED)
!                                          Not enough high values
                XWRK1 = -XHUGE
                NORD = NLEQ - NEQU - NMED + 1
                IDON1 = 0
                ICRS1 = 1
                ICRS2 = 0
                Do IDON = 1, NTRI, 7
                   IDON1 = IDON1 + 1
                   If (ICRS2 < NORD) Then
!
                      Do ICRS = IDON, IENDT (IDON1)
                         If (XWRKT(ICRS) > XWRK1) Then
                            XWRK = XWRKT (ICRS)
                            IDCR = ICRS1 - 1
                            Do IDCR = ICRS1 - 1, 1, - 1
                               If (XWRK <= XWRKT(IDCR)) Exit
                               XWRKT (IDCR+1) = XWRKT (IDCR)
                            End Do
                            XWRKT (IDCR+1) = XWRK
                            XWRK1 = XWRKT(ICRS1)
                         Else
                            If (ICRS2 < NORD) Then
                               XWRKT (ICRS1) = XWRKT (ICRS)
                               XWRK1 = XWRKT (ICRS1)
                            End If
                         End If
                         ICRS1 = MIN (NORD, ICRS1 + 1)
                         ICRS2 = MIN (NORD, ICRS2 + 1)
                      End Do
                   Else
                      Do ICRS = IENDT (IDON1), IDON, -1
                         If (XWRKT(ICRS) > XWRK1) Then
                            XWRK = XWRKT (ICRS)
                            IDCR = ICRS1 - 1
                            Do IDCR = ICRS1 - 1, 1, - 1
                               If (XWRK <= XWRKT(IDCR)) Exit
                               XWRKT (IDCR+1) = XWRKT (IDCR)
                            End Do
                            XWRKT (IDCR+1) = XWRK
                            XWRK1 = XWRKT(ICRS1)
                         Else
                            Exit
                         End If
                      End Do
                   Endif
                End Do
!
                res_med = XWRK1
                Return
         End If
!
End Function I_valmed
end module sort_interface
