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

end module sort_interface
