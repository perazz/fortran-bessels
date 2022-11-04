program bessels_test
    use bessels
    implicit none

    integer :: nfailed=0,npassed=0

    call add_test(test_cuberoot())
    call add_test(test_bessel_j0())
    call add_test(test_bessel_j1())
    call add_test(test_bessel_y0())
    call add_test(test_bessel_j0_cputime())
    call add_test(test_bessel_j1_cputime())
    call add_test(test_bessel_y0_cputime())

    print 1, npassed+nfailed,npassed,nfailed
    if (nfailed>0) then
        stop -1
    else
        stop 0
    endif

    1 format('[bessels] ',i0,' test completed: ',i0,' passed, ',i0,' failed.')

    contains

    subroutine add_test(success)
        logical, intent(in) :: success
        if (success) then
            npassed = npassed+1
        else
            nfailed = nfailed+1
        end if
    end subroutine add_test

    ! Test bessel j0 cpu time
    logical function test_bessel_j0_cputime() result(success)

        integer, parameter :: nsize = 100000
        integer, parameter :: ntest = 100
        real(BK), parameter :: xmin = -1e+6_BK
        real(BK), parameter :: xmax =  1e+6_BK
        real(BK), allocatable :: x(:),intrin(:),packge(:),z(:)
        integer :: i
        real(BK) :: time,timep,c_start,c_end
        allocate(x(nsize),intrin(nsize),packge(nsize),z(ntest))

        call random_number(x)
        x    = xmin*(ONE-x) + xmax*x

        time = ZERO
        do i=1,ntest
            call cpu_time(c_start)
            intrin = bessel_j0(x)
            call cpu_time(c_end)
            z(i) = sum(intrin)
            time = time+c_end-c_start
        end do
        print "('[bessel_j0] INTRINSIC time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*time/(nsize*ntest),sum(z)

        timep = ZERO
        do i=1,ntest
            call cpu_time(c_start)
            packge = besselj0(x)
            call cpu_time(c_end)
            z(i) = sum(packge)
            timep = timep+c_end-c_start
        end do
        print "('[bessel_j0] PACKAGE   time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*timep/(nsize*ntest),sum(z)

        success = timep<time

    end function test_bessel_j0_cputime

    ! Test bessel j0 function
    logical function test_bessel_j0() result(success)

      integer, parameter :: NTEST = 100000

      real(BK), parameter :: xmin = -1e+6_BK
      real(BK), parameter :: xmax =  1e+6_BK
      real(BK), parameter :: RTOL =  1e-6_BK
      real(BK), parameter :: ATOL =  1e-10_BK
      real(BK) :: x(NTEST),fun(NTEST),intr(NTEST),err(NTEST)
      integer :: i

      ! Randoms in range
      call random_number(x)
      x    = xmin*(ONE-x) + xmax*x
      fun  = besselj0(x)
      intr = bessel_j0(x)
      err  = abs(fun-intr)*rewt(intr,RTOL,ATOL)

      success = all(err<one)

      if (.not.success) then
         do i=1,NTEST
            if (err(i)>=one) &
            print *, '[bessel_j0] x=',x(i),' package=',fun(i),' intrinsic=',intr(i),' relerr=',err(i)
         end do
      end if

    end function test_bessel_j0

    ! Test bessel j0 function
    logical function test_bessel_j1() result(success)

      integer, parameter :: NTEST = 100000

      real(BK), parameter :: xmin = -1e+6_BK
      real(BK), parameter :: xmax =  1e+6_BK
      real(BK), parameter :: RTOL =  1e-6_BK
      real(BK), parameter :: ATOL =  1e-10_BK
      real(BK) :: x(NTEST),fun(NTEST),intr(NTEST),err(NTEST)
      integer :: i

      ! Randoms in range
      call random_number(x)
      x    = xmin*(ONE-x) + xmax*x
      fun  = besselj1(x)
      intr = bessel_j1(x)
      err  = abs(fun-intr)*rewt(intr,RTOL,ATOL)

      success = all(err<one)

      if (.not.success) then
         do i=1,NTEST
            if (err(i)>=one) &
            print *, '[bessel_j1] x=',x(i),' package=',fun(i),' intrinsic=',intr(i),' relerr=',err(i)
         end do
      end if

    end function test_bessel_j1

    ! Test bessel j0 cpu time
    logical function test_bessel_j1_cputime() result(success)

        integer, parameter :: nsize = 100000
        integer, parameter :: ntest = 100
        real(BK), parameter :: xmin = -1e+6_BK
        real(BK), parameter :: xmax =  1e+6_BK
        real(BK), allocatable :: x(:),intrin(:),packge(:),z(:)
        integer :: i
        real(BK) :: time,timep,c_start,c_end
        allocate(x(nsize),intrin(nsize),packge(nsize),z(ntest))

        call random_number(x)
        x    = xmin*(ONE-x) + xmax*x

        time = ZERO
        do i=1,ntest
            call cpu_time(c_start)
            intrin = bessel_j1(x)
            call cpu_time(c_end)
            z(i) = sum(intrin)
            time = time+c_end-c_start
        end do
        print "('[bessel_j1] INTRINSIC time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*time/(nsize*ntest),sum(z)

        timep = ZERO
        do i=1,ntest
            call cpu_time(c_start)
            packge = besselj1(x)
            call cpu_time(c_end)
            z(i) = sum(packge)
            timep = timep+c_end-c_start
        end do
        print "('[bessel_j1] PACKAGE   time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*timep/(nsize*ntest),sum(z)

        success = timep<time

    end function test_bessel_j1_cputime


    ! Test bessel j0 function
    logical function test_bessel_y0() result(success)

      integer, parameter :: NTEST = 100000

      real(BK), parameter :: xmin =    0.0_BK
      real(BK), parameter :: xmax =  1e+12_BK
      real(BK), parameter :: RTOL =  1e-6_BK
      real(BK), parameter :: ATOL =  1e-10_BK
      real(BK) :: x(NTEST),fun(NTEST),intr(NTEST),err(NTEST)
      integer :: i

      ! Randoms in range
      call random_number(x)
      x    = xmin*(ONE-x) + xmax*x
      fun  = bessely0(x)
      intr = bessel_y0(x)
      err  = abs(fun-intr)*rewt(intr,RTOL,ATOL)

      success = all(err<one)

      if (.not.success) then
         do i=1,NTEST
            if (err(i)>=one) &
            print *, '[bessel_y0] x=',x(i),' package=',fun(i),' intrinsic=',intr(i),' relerr=',err(i)
         end do
      end if

    end function test_bessel_y0

    ! Test bessel j0 cpu time
    logical function test_bessel_y0_cputime() result(success)

        integer, parameter :: nsize = 100000
        integer, parameter :: ntest = 100
        real(BK), parameter :: xmin =    0.0_BK
        real(BK), parameter :: xmax =  1e+12_BK
        real(BK), allocatable :: x(:),intrin(:),packge(:),z(:)
        integer :: i
        real(BK) :: time,timep,c_start,c_end
        allocate(x(nsize),intrin(nsize),packge(nsize),z(ntest))

        call random_number(x)
        x    = xmin*(ONE-x) + xmax*x

        time = ZERO
        do i=1,ntest
            call cpu_time(c_start)
            intrin = bessel_y0(x)
            call cpu_time(c_end)
            z(i) = sum(intrin)
            time = time+c_end-c_start
        end do
        print "('[bessel_y0] INTRINSIC time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*time/(nsize*ntest),sum(z)

        timep = ZERO
        do i=1,ntest
            call cpu_time(c_start)
            packge = bessely0(x)
            call cpu_time(c_end)
            z(i) = sum(packge)
            timep = timep+c_end-c_start
        end do
        print "('[bessel_y0] PACKAGE   time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*timep/(nsize*ntest),sum(z)

        success = timep<time

    end function test_bessel_y0_cputime

    ! ode-like inverse error weight
    elemental real(BK) function rewt(x,RTOL,ATOL)
       real(BK), intent(in) :: x,RTOL,ATOL
       rewt = ONE/(RTOL*abs(x)+ATOL)
    end function rewt

    ! Test approximated cube root
    logical function test_cuberoot() result(success)

      integer, parameter :: NTEST = 1000

      real(BK), parameter :: xmin = -1e+6_BK
      real(BK), parameter :: xmax =  1e+6_BK
      real(BK), parameter :: RTOL =  1e-6_BK
      real(BK), parameter :: ATOL =  1e-20_BK
      real(BK) :: x(NTEST),fun(NTEST),intr(NTEST)
      integer :: i

      ! Randoms in range
      call random_number(x)
      x = xmin*(ONE-x) + xmax*x

      fun  = cbrt(x)
      intr = sign(abs(x)**THIRD,x)

      success = all(abs(fun-intr)<=RTOL*abs(intr)+ATOL)

      if (.not.success) then
         do i=1,NTEST
            print *, 'x=',x(i),' cbrt=',fun(i),' intrinsic=',intr(i),' relerr=',abs(fun(i)-intr(i))
         end do
      end if

    end function test_cuberoot

end program bessels_test
