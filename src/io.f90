module ioFunc

    use iso_fortran_env, only: RK => real64
    implicit none
    
contains

    subroutine readParams(fpath, totL, n, l, m, Nth)
        character(4096), intent(in) :: fpath
        integer, intent(out) :: totL, n, l, m, Nth

        integer :: nu, ios

        open(newunit=nu, file = fpath, status='old', iostat=ios)
        if ( ios /= 0 ) stop "Error opening data file."

        ! Read Nth
        read(nu, *, iostat=ios) Nth
        if (ios /= 0) STOP "Error while reading Nth from data file."

        ! Read n, m, l
        read(nu, *, iostat=ios) n, m, l
        if (ios /= 0) STOP "Error while reading n,m,l from data file."

        close(nu) 
        totL = n*m*l
    end subroutine readParams

    subroutine readData(fpath, totL, x, y, z, O)
        
        character(4096), intent(in) :: fpath
        integer, intent(in) :: totL
        real(RK), dimension(totL), intent(out) :: x, y, z, O
    
        integer :: nu, ios, i, nd, md, ld
        
        open(newunit=nu, file = fpath, status='old', iostat=ios)
        if ( ios /= 0 ) stop "Error opening data file."

        ! Skip Nth, n, m, l
        read(nu, *, iostat=ios) nd
        if (ios /= 0) STOP "Error skipping reading Nth from data file."
        read(nu, *, iostat=ios) nd, md, ld
        if (ios /= 0) STOP "Error skipping reading n,m,l from data file."

        ! Read x, y, z, O
        read(nu, *, iostat=ios) (x(i), i = 1, totL)
        if (ios /= 0) STOP "Error while reading x from data file."
        read(nu, *, iostat=ios) (y(i), i = 1, totL)
        if (ios /= 0) STOP "Error while reading y from data file."
        read(nu, *, iostat=ios) (z(i), i = 1, totL)
        if (ios /= 0) STOP "Error while reading z from data file."
        read(nu, *, iostat=ios) (O(i), i = 1, totL)
        if (ios /= 0) STOP "Error while reading O from data file."
        
        close(nu)

        return
    end subroutine readData


    subroutine saveOutput(Oin, Oout, Wout, x, y, z, n, l, m, totL, fpath)
        
        integer, intent(in) :: totL, n, m, l
        real(RK), dimension(totL), intent(in) :: Oin, Oout, Wout, x, y, z
        character(len=4096), intent(in) :: fpath

        integer :: nu, ios, i
        
        open(newunit=nu, file = fpath, status='new', iostat=ios)
        if ( ios /= 0 ) stop "Error opening output file."

        write(nu, *, iostat=ios) n, m, l
        if ( ios /= 0 ) stop "Error writing output file."
        
        write(nu, *, iostat=ios) (x(i), i = 1, totL)
        if ( ios /= 0 ) stop "Error writing output file."
        write(nu, *, iostat=ios) (y(i), i = 1, totL)
        if ( ios /= 0 ) stop "Error writing output file."
        write(nu, *, iostat=ios) (z(i), i = 1, totL)
        if ( ios /= 0 ) stop "Error writing output file."

        write(nu, *, iostat=ios) (Oin(i), i = 1, totL)
        if ( ios /= 0 ) stop "Error writing output file."

        write(nu, *, iostat=ios) (Oout(i), i = 1, totL)
        if ( ios /= 0 ) stop "Error writing output file."

        write(nu, *, iostat=ios) (Wout(i), i = 1, totL)
        if ( ios /= 0 ) stop "Error opening output file."

        close(nu)

    end subroutine saveOutput
    
end module ioFunc