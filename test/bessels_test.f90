program bessels_test
    use bessels
    implicit none

    integer :: nfailed=0,npassed=0

    call add_test(test_cuberoot())
    call add_test(test_bessel_j0())
    call add_test(test_bessel_j1())
    call add_test(test_bessel_y0())
    call add_test(test_bessel_y1())
    call add_test(test_bessel_k0())
    call add_test(test_bessel_k1())
    call add_test(test_bessel_i0())
    call add_test(test_bessel_i1())
    call add_test(test_bessel_j0_cputime())
    call add_test(test_bessel_j1_cputime())
    call add_test(test_bessel_y0_cputime())
    call add_test(test_bessel_y1_cputime())
    call add_test(test_bessel_k0_cputime())
    call add_test(test_bessel_k1_cputime())
    call add_test(test_bessel_i0_cputime())
    call add_test(test_bessel_i1_cputime())

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
        real(BK), parameter :: xmin = -1e+3_BK
        real(BK), parameter :: xmax =  1e+3_BK
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

      real(BK), parameter :: xmin = -1e+3_BK
      real(BK), parameter :: xmax =  1e+3_BK
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

      real(BK), parameter :: xmin = -1e+3_BK
      real(BK), parameter :: xmax =  1e+3_BK
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
        real(BK), parameter :: xmin = -1e+3_BK
        real(BK), parameter :: xmax =  1e+3_BK
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

      real(BK), parameter :: xmin =   0.0_BK
      real(BK), parameter :: xmax =  1e+3_BK
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
        real(BK), parameter :: xmax =  1e+3_BK
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

    ! Test bessel j0 function
    logical function test_bessel_y1() result(success)

      integer, parameter :: NTEST = 100000

      real(BK), parameter :: xmin =   0.0_BK
      real(BK), parameter :: xmax =  1e+3_BK
      real(BK), parameter :: RTOL =  1e-6_BK
      real(BK), parameter :: ATOL =  1e-10_BK
      real(BK) :: x(NTEST),fun(NTEST),intr(NTEST),err(NTEST)
      integer :: i

      ! Randoms in range
      call random_number(x)
      x    = xmin*(ONE-x) + xmax*x
      fun  = bessely1(x)
      intr = bessel_y1(x)
      err  = abs(fun-intr)*rewt(intr,RTOL,ATOL)

      success = all(err<one)

      if (.not.success) then
         do i=1,NTEST
            if (err(i)>=one) &
            print *, '[bessel_y1] x=',x(i),' package=',fun(i),' intrinsic=',intr(i),' relerr=',err(i)
         end do
      end if

    end function test_bessel_y1

    ! Test bessel j0 cpu time
    logical function test_bessel_y1_cputime() result(success)

        integer, parameter :: nsize = 100000
        integer, parameter :: ntest = 100
        real(BK), parameter :: xmin =    0.0_BK
        real(BK), parameter :: xmax =  1e+3_BK
        real(BK), allocatable :: x(:),intrin(:),packge(:),z(:)
        integer :: i
        real(BK) :: time,timep,c_start,c_end
        allocate(x(nsize),intrin(nsize),packge(nsize),z(ntest))

        call random_number(x)
        x    = xmin*(ONE-x) + xmax*x

        time = ZERO
        do i=1,ntest
            call cpu_time(c_start)
            intrin = bessel_y1(x)
            call cpu_time(c_end)
            z(i) = sum(intrin)
            time = time+c_end-c_start
        end do
        print "('[bessel_y1] INTRINSIC time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*time/(nsize*ntest),sum(z)

        timep = ZERO
        do i=1,ntest
            call cpu_time(c_start)
            packge = bessely1(x)
            call cpu_time(c_end)
            z(i) = sum(packge)
            timep = timep+c_end-c_start
        end do
        print "('[bessel_y1] PACKAGE   time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*timep/(nsize*ntest),sum(z)

        success = timep<time

    end function test_bessel_y1_cputime

    ! Test bessel k0 function
    logical function test_bessel_k0() result(success)
      use rjk, only: RKBESL

      integer, parameter :: NTEST = 100000

      real(BK), parameter :: xmin =   0.0_BK
      real(BK), parameter :: xmax =  1e+6_BK
      real(BK), parameter :: RTOL =  1e-6_BK
      real(BK), parameter :: ATOL =  1e-10_BK
      real(BK) :: x(NTEST),fun(NTEST),intr(NTEST),err(NTEST)
      integer :: i,ierr

      ! Randoms in range
      call random_number(x)
      x    = xmin*(ONE-x) + xmax*x
      fun  = besselk0(x)

      do i=1,NTEST
         CALL RKBESL(X=x(i), ALPHA=zero, NB=1, IZE=1, BK=intr(i), NCALC=ierr)
      end do

      err  = abs(fun-intr)*rewt(intr,RTOL,ATOL)

      success = all(err<one)

      if (.not.success) then
         do i=1,NTEST
            if (err(i)>=one) &
            print *, '[bessel_k0] x=',x(i),' package=',fun(i),' intrinsic=',intr(i),' relerr=',err(i)
         end do
      end if

    end function test_bessel_k0

    ! Test bessel j0 cpu time
    logical function test_bessel_k0_cputime() result(success)
        use rjk

        integer, parameter :: nsize = 100000
        integer, parameter :: ntest = 100
        real(BK), parameter :: xmin =    0.0_BK
        real(BK), parameter :: xmax =  1e+2_BK
        real(BK), allocatable :: x(:),intrin(:),packge(:),z(:)
        integer :: i,j,ierr
        real(BK) :: time,timep,c_start,c_end
        allocate(x(nsize),intrin(nsize),packge(nsize),z(ntest))

        call random_number(x)
        x    = xmin*(ONE-x) + xmax*x

        time = ZERO
        do i=1,ntest
            call cpu_time(c_start)

            do j=1,nsize
               CALL RKBESL(X=x(j), ALPHA=zero, NB=1, IZE=1, BK=intrin(j), NCALC=ierr)
               if (ierr/=1) then
                  print *, 'RKBESL error: x=',x(i),' ierr=',ierr
                  stop 'RKBESL error'
               endif
            end do

            call cpu_time(c_end)
            z(i) = sum(intrin)
            time = time+c_end-c_start
        end do
        print "('[bessel_k0] NETLIB    time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*time/(nsize*ntest),sum(z)

        timep = ZERO
        do i=1,ntest
            call cpu_time(c_start)
            packge = besselk0(x)
            call cpu_time(c_end)
            z(i) = sum(packge)
            timep = timep+c_end-c_start
        end do
        print "('[bessel_k0] PACKAGE   time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*timep/(nsize*ntest),sum(z)

        success = timep<time

    end function test_bessel_k0_cputime

    ! Test bessel k0 function
    logical function test_bessel_k1() result(success)
      use rjk, only: RKBESL

      integer, parameter :: NTEST = 100000

      real(BK), parameter :: xmin =   0.0_BK

      ! Limit the max x range to RKBESL validity
      real(BK), parameter :: xmax =  10.0_BK
      real(BK), parameter :: RTOL =  1e-6_BK
      real(BK), parameter :: ATOL =  1e-10_BK
      real(BK) :: x(NTEST),fun(NTEST),intr(0:NTEST),err(NTEST),this(2)
      integer :: i,ierr

      ! Randoms in range
      call random_number(x)
      x    = xmin*(ONE-x) + xmax*x
      fun  = besselk1(x)

      do i=1,NTEST
         CALL RKBESL(X=x(i), ALPHA=ZERO, NB=2, IZE=1, BK=this, NCALC=ierr)
         intr(i) = this(2)
         if (ierr/=2) then
            print *, 'RKBESL error: x=',x(i),' ierr=',ierr
            stop 'RKBESL error'
         endif
      end do

      err  = abs(fun-intr(1:))*rewt(intr(1:),RTOL,ATOL)

      success = all(err<one)

      if (.not.success) then
         do i=1,NTEST
            if (err(i)>=one) &
            print *, '[bessel_k1] x=',x(i),' package=',fun(i),' intrinsic=',intr(i),' relerr=',err(i)
         end do
      end if

    end function test_bessel_k1

    ! Test bessel j0 cpu time
    logical function test_bessel_k1_cputime() result(success)
        use rjk

        integer, parameter :: nsize = 100000
        integer, parameter :: ntest = 100
        real(BK), parameter :: xmin =    0.0_BK
        ! Limit the max x range to RKBESL validity
        real(BK), parameter :: xmax =  600.0_BK
        real(BK), allocatable :: x(:),intrin(:),packge(:),z(:)
        integer :: i,j,ierr
        real(BK) :: time,timep,c_start,c_end,this(2)
        allocate(x(nsize),intrin(nsize+1),packge(nsize),z(ntest))

        call random_number(x)
        x    = xmin*(ONE-x) + xmax*x

        time = ZERO
        do i=1,ntest
            call cpu_time(c_start)

            do j=1,nsize
               CALL RKBESL(X=x(j), ALPHA=ZERO, NB=2, IZE=1, BK=intrin(j:j+1), NCALC=ierr)
            end do

            call cpu_time(c_end)
            z(i) = sum(intrin(1:nsize))
            time = time+c_end-c_start
        end do
        print "('[bessel_k1] NETLIB    time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*time/(nsize*ntest),sum(z)

        timep = ZERO
        do i=1,ntest
            call cpu_time(c_start)
            packge = besselk1(x)
            call cpu_time(c_end)
            z(i) = sum(packge)
            timep = timep+c_end-c_start
        end do
        print "('[bessel_k1] PACKAGE   time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*timep/(nsize*ntest),sum(z)

        success = timep<time

    end function test_bessel_k1_cputime

    ! Test bessel i0 function
    logical function test_bessel_i0() result(success)
      use rji, only: RIBESL

      integer, parameter :: NTEST = 2000
      real(BK), parameter :: xmin =   0.0_BK

      ! Limit the max x range to RKBESL validity
      real(BK), parameter :: xmax =  500.0_BK
      real(BK), parameter :: RTOL =  1e-6_BK
      real(BK), parameter :: ATOL =  1e-10_BK
      real(BK) :: x(NTEST),fun(NTEST),intr(NTEST),err(NTEST),this(2)
      integer :: i,ierr

      ! Randoms in range
      call random_number(x)
      x    = xmin*(ONE-x) + xmax*x
      fun  = besseli0(x)

      do i=1,NTEST
         CALL RIBESL(X=x(i), ALPHA=ZERO, NB=1, IZE=1, B=this, NCALC=ierr)
         intr(i) = this(1)
         if (ierr/=1) then
            print *, 'RIBESL error: x=',x(i),' ierr=',ierr
            stop 'RIBESL error'
         endif
      end do

      err  = abs(fun-intr)*rewt(intr,RTOL,ATOL)

      success = all(err<one)

      if (.not.success) then
         do i=1,NTEST
            if (err(i)>=one) &
            print *, '[bessel_i0] x=',x(i),' package=',fun(i),' intrinsic=',intr(i),' relerr=',err(i)
         end do
      end if

    end function test_bessel_i0

    ! Test bessel j0 cpu time
    logical function test_bessel_i0_cputime() result(success)
        use rji

        integer, parameter :: nsize = 100000
        integer, parameter :: ntest = 100
        real(BK), parameter :: xmin =    0.0_BK
        ! Limit the max x range to RKBESL validity
        real(BK), parameter :: xmax =  600.0_BK
        real(BK), allocatable :: x(:),intrin(:),packge(:),z(:)
        integer :: i,j,ierr
        real(BK) :: time,timep,c_start,c_end,this(2)
        allocate(x(nsize),intrin(nsize+1),packge(nsize),z(ntest))

        call random_number(x)
        x    = xmin*(ONE-x) + xmax*x

        time = ZERO
        do i=1,ntest
            call cpu_time(c_start)

            do j=1,nsize
               CALL RIBESL(X=x(j), ALPHA=ZERO, NB=1, IZE=1, B=intrin(j), NCALC=ierr)
            end do

            call cpu_time(c_end)
            z(i) = sum(intrin(1:nsize))
            time = time+c_end-c_start
        end do
        print "('[bessel_i0] NETLIB    time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*time/(nsize*ntest),sum(z)

        timep = ZERO
        do i=1,ntest
            call cpu_time(c_start)
            packge = besseli0(x)
            call cpu_time(c_end)
            z(i) = sum(packge)
            timep = timep+c_end-c_start
        end do
        print "('[bessel_i0] PACKAGE   time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*timep/(nsize*ntest),sum(z)

        success = timep<time

    end function test_bessel_i0_cputime

    ! Test bessel k0 function
    logical function test_bessel_i1() result(success)
      use rji

      integer, parameter :: NTEST = 2000

      real(BK), parameter :: xmin =   0.0_BK

      ! Limit the max x range to RKBESL validity
      real(BK), parameter :: xmax =  50.0_BK
      real(BK), parameter :: RTOL =  1e-6_BK
      real(BK), parameter :: ATOL =  1e-10_BK
      real(BK) :: x(NTEST),fun(NTEST),intr(0:NTEST),err(NTEST),this(2)
      integer :: i,ierr

      ! Randoms in range
      call random_number(x)
      x    = xmin*(ONE-x) + xmax*x
      fun  = besseli1(x)

      do i=1,NTEST
         CALL RIBESL(X=x(i), ALPHA=ZERO, NB=2, IZE=1, B=this, NCALC=ierr)
         intr(i) = this(2)
         if (ierr/=2) then
            print *, 'RIBESL error: x=',x(i),' ierr=',ierr
            stop 'RIBESL error'
         endif
      end do

      err  = abs(fun-intr(1:))*rewt(intr(1:),RTOL,ATOL)

      success = all(err<one)

      if (.not.success) then
         do i=1,NTEST
            if (err(i)>=one) &
            print *, '[bessel_i1] x=',x(i),' package=',fun(i),' intrinsic=',intr(i),' relerr=',err(i)
         end do
      end if

    end function test_bessel_i1

    ! Test bessel j0 cpu time
    logical function test_bessel_i1_cputime() result(success)
        use rji

        integer, parameter :: nsize = 100000
        integer, parameter :: ntest = 100
        real(BK), parameter :: xmin =    0.0_BK
        real(BK), parameter :: xmax =  100.0_BK
        real(BK), allocatable :: x(:),intrin(:),packge(:),z(:)
        integer :: i,j,ierr
        real(BK) :: time,timep,c_start,c_end,this(2)
        allocate(x(nsize),intrin(nsize+1),packge(nsize),z(ntest))

        call random_number(x)
        x    = xmin*(ONE-x) + xmax*x

        time = ZERO
        do i=1,ntest
            call cpu_time(c_start)

            do j=1,nsize
               CALL RIBESL(X=x(j), ALPHA=ZERO, NB=2, IZE=1, B=intrin(j:j+1), NCALC=ierr)

               if (ierr/=2) then
                  print *, 'RIBESL error: x=',x(i),' ierr=',ierr
                  stop 'RIBESL error'
               endif

            end do

            call cpu_time(c_end)
            z(i) = sum(intrin(1:nsize))
            time = time+c_end-c_start
        end do
        print "('[bessel_i1] NETLIB    time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*time/(nsize*ntest),sum(z)

        timep = ZERO
        do i=1,ntest
            call cpu_time(c_start)
            packge = besseli1(x)
            call cpu_time(c_end)
            z(i) = sum(packge)
            timep = timep+c_end-c_start
        end do
        print "('[bessel_i1] PACKAGE   time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*timep/(nsize*ntest),sum(z)

        success = timep<time

    end function test_bessel_i1_cputime

    ! ode-like inverse error weight
    elemental real(BK) function rewt(x,RTOL,ATOL)
       real(BK), intent(in) :: x,RTOL,ATOL
       rewt = ONE/(RTOL*abs(x)+ATOL)
    end function rewt

    ! Test approximated cube root
    logical function test_cuberoot() result(success)

      integer, parameter :: NTEST = 1000

      real(BK), parameter :: xmin = -1e+3_BK
      real(BK), parameter :: xmax =  1e+3_BK
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
