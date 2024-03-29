module utils

    use iso_fortran_env, only: RK => real64
    implicit none

    public :: Cflatten, Cpack_3d, checkInputArgs

contains

    subroutine checkInputArgs(args)

        character(len=*), dimension(2), intent(out) :: args
        integer :: num_args, ii


        num_args = command_argument_count()
        if (num_args.ne.size(args)) then
            write(*,'(A, I1, A)') "I expect ", size(args), " input arguments: the H5 file with the data and the name of the variable to smooth."
            STOP
        end if

        do ii = 1, num_args
            call get_command_argument(ii,args(ii))
        end do
        
    end subroutine checkInputArgs

    subroutine Cflatten(xi, xo)

        real(RK), intent(in) :: xi(:, :, :)
        real(RK), intent(inout) :: xo(:)

        integer, dimension(3) :: sh
        integer :: i1, j1, k1, l1

        sh = shape(xi)
        do k1 = 1, sh(3)
            do j1 = 1, sh(2)
                do i1 = 1, sh(1)
                    l1 = k1 + sh(3) * (j1-1) + sh(3)*sh(2)*(i1-1)
                    xo(l1) = xi(i1, j1, k1)
                end do
            end do
        end do

    end subroutine Cflatten

    subroutine Cpack_3d(xi, xo, d1, d2, d3)
        real(RK), intent(in) :: xi(:)
        real(RK), intent(inout) :: xo(:, :, :)
        integer, intent(in) :: d1, d2, d3

        integer :: ii, i, j, k

        !$OMP PARALLEL DO SIMD DEFAULT(SHARED) PRIVATE(ii, i,j,k) NUM_THREADS(4)
        do ii = 1, size(xi)
            k = mod(ii - 1, d3) + 1
            j = mod((ii - k) / d3, d2) + 1
            i = mod(( (ii-k)/d3 - (j-1) ) / d2 , d1) + 1

            xo(i,j,k) = xi(ii)
        end do
        !$OMP END PARALLEL DO SIMD

    end subroutine Cpack_3d
    
end module utils


module TimerModule
    use iso_fortran_env, only: RK => real64
    implicit none
  
    type :: TimerClass
      private
      integer :: rate
      integer :: startTime
      real(RK) :: system_time, secondsPerHour, secondsPerMinute
      logical :: alreadyStarted
    contains
      procedure :: start, stop
      procedure :: printTime, initializeTimer
    end type TimerClass
  
    
contains
  
    subroutine initializeTimer(this)
        class(TimerClass), intent(inout) :: this
        this%secondsPerHour = 3600
        this%secondsPerMinute = 60
        this%alreadyStarted = .false.
        call system_clock(count_rate=this%rate)
    end subroutine initializeTimer

    subroutine start(this)
      class(TimerClass), intent(inout) :: this
    
      if (this%alreadyStarted.eqv..true.) then
        write(*, '(A)') 'Timer already started.'
        return
      end if

      this%alreadyStarted = .true.
      call system_clock(this%startTime)
    end subroutine start
  
    subroutine stop(this)
      class(TimerClass), intent(inout) :: this  
      integer :: endTime

      call system_clock(endTime)

      if (this%alreadyStarted.eqv..false.) then
        write(*, '(A)') 'Timer not started.'
        return
      end if

      this%system_time = real(endTime - this%startTime) / real(this%rate)
      this%alreadyStarted = .false.

    end subroutine stop

    subroutine printTime(this)
        class(TimerClass), intent(inout) :: this
        integer :: hours, minutes, seconds

        hours = int(this%system_time / this%secondsPerHour)
        minutes = int(mod(this%system_time, this%secondsPerHour) / this%secondsPerMinute)
        seconds = int(mod(this%system_time, this%secondsPerMinute))
        write(*, '(A, I0.3, A, I0.2, A, I0.2)') "Total system runtime ", hours, ':', minutes, ':', seconds
    end subroutine printTime
  
end module TimerModule
