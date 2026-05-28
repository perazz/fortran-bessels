!  ************************************************************************************************************
!
!                                    ____  ________________ ________   _____
!                                   / __ )/ ____/ ___/ ___// ____/ /  / ___/
!                                  / __  / __/  \__ \\__ \/ __/ / /   \__ \
!                                 / /_/ / /___ ___/ /__/ / /___/ /______/ /
!                                /_____/_____//____/____/_____/_____/____/
!
!                                              Test programs
!
!  MIT License
!
!  Copyright (c) 2022 Federico Perini
!  Copyright (c) 2021-2022 Michael Helton, Oscar Smith, and the Bessels.jl contributors
!
!  ************************************************************************************************************
program bessels_test
    use bessels
    implicit none

    integer :: nfailed=0,npassed=0

    call add_test(test_cuberoot())
    call add_test(test_bessel_j0())
    call add_test(test_bessel_j1())
    call add_test(test_bessel_jn())
    call add_test(test_bessel_y0())
    call add_test(test_bessel_y1())
    call add_test(test_besselj_up_recurrence())
    call add_test(test_bessely_nu_integer())
    call add_test(test_hankelh1_consistency())
    call add_test(test_hankelh2_conjugate())
    call add_test(test_bessel_k0())
    call add_test(test_bessel_k1())
    call add_test(test_bessel_i0())
    call add_test(test_bessel_i1())
    call add_test(test_gamma())
    call add_test(test_bessel_j0_cputime())
    call add_test(test_bessel_j1_cputime())
    call add_test(test_bessel_jn_cputime())
    call add_test(test_bessel_y0_cputime())
    call add_test(test_bessel_y1_cputime())
    call add_test(test_bessely_nu_cputime())
    call add_test(test_hankelh1_cputime())
    call add_test(test_bessel_k0_cputime())
    call add_test(test_bessel_k1_cputime())
    call add_test(test_bessel_i0_cputime())
    call add_test(test_bessel_i1_cputime())
    call add_test(test_gamma_cputime())

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

        success = timep<2*time

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

    ! Test bessel j function for integer order
    logical function test_bessel_jn() result(success)

      integer, parameter :: NTEST = 100000

      real(BK), parameter :: xmin = -1e+3_BK
      real(BK), parameter :: xmax =  1e+3_BK
      real(BK), parameter :: RTOL =  1e-6_BK
      real(BK), parameter :: ATOL =  1e-8_BK
      real(BK) :: x(NTEST),fun(NTEST),intr(NTEST),err(NTEST)
      integer :: i,n

      ! Randoms in range
      call random_number(x)
      x    = xmin*(ONE-x) + xmax*x

      success = .true.

      do n=0,10

          do i=1,NTEST
          fun(i)  = besseljn(n,x(i))
          end do
          intr = bessel_jn(n,x)
          err  = abs(fun-intr)*rewt(intr,RTOL,ATOL)

          success = success .and. all(err<one)

          if (.not.all(err<one)) then
             do i=1,NTEST
                if (err(i)>=one) then
                fun(i)  = besseljn(n,x(i))
                print *, '[bessel_jn] x=',x(i),' n=',n,' package=',fun(i),' intrinsic=',intr(i),' relerr=',err(i)
                stop
                endif
             end do
          end if

      end do

    end function test_bessel_jn

    ! Test bessel j0 cpu time
    logical function test_bessel_jn_cputime() result(success)

        integer, parameter :: nsize = 10000
        integer, parameter :: ntest = 100
        real(BK), parameter :: xmin = -1e+3_BK
        real(BK), parameter :: xmax =  1e+3_BK
        real(BK), allocatable :: x(:),intrin(:),packge(:),z(:)
        integer :: i,n
        real(BK) :: time,timep,c_start,c_end
        allocate(x(nsize),intrin(nsize),packge(nsize),z(ntest))

        call random_number(x)
        x    = xmin*(ONE-x) + xmax*x

        time = ZERO
        do n=0,9
        do i=1,ntest
            call cpu_time(c_start)
            intrin = bessel_jn(n,x)
            call cpu_time(c_end)
            z(i) = sum(intrin)
            time = time+c_end-c_start
        end do
        end do
        print "('[bessel_jn] INTRINSIC time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*time/(10*nsize*ntest),sum(z)

        timep = ZERO
        do n=0,9
        do i=1,ntest
            call cpu_time(c_start)
            packge = besseljn(n,x)
            call cpu_time(c_end)
            z(i) = sum(packge)
            timep = timep+c_end-c_start
        end do
        end do
        print "('[bessel_jn] PACKAGE   time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*timep/(10*nsize*ntest),sum(z)

        success = timep<2*time

    end function test_bessel_jn_cputime

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

        success = timep<2*time

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

        success = timep<2*time

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

        success = timep<2*time

    end function test_bessel_y1_cputime

    ! Test bessel k0 function
    logical function test_bessel_k0() result(success)
      use rjk, only: RKBESL

      integer, parameter :: NTEST = 100000

      real(BK), parameter :: xmin =   0.0_BK
      real(BK), parameter :: xmax =  1e+6_BK
      real(BK), parameter :: RTOL =  1e-6_BK
      real(BK), parameter :: ATOL =  1e-10_BK
      real(BK), allocatable, dimension(:) :: x,fun,intr,err
      integer :: i,ierr

      allocate(x(NTEST),fun(NTEST),intr(0:NTEST),err(NTEST))

      ! Randoms in range
      call random_number(x)
      x    = xmin*(ONE-x) + xmax*x
      fun  = besselk0(x)

      do i=1,NTEST
         CALL RKBESL(X=x(i), ALPHA=zero, NB=1, IZE=1, BK=intr(i), NCALC=ierr)
      end do

      err  = abs(fun-intr(1:))*rewt(intr,RTOL,ATOL)

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

        success = timep<2*time

    end function test_bessel_k0_cputime

    ! Test bessel k0 function
    logical function test_bessel_k1() result(success)
      use rjk, only: RKBESL

      integer, parameter :: NTEST = 10000

      real(BK), parameter :: xmin =   0.0_BK

      ! Limit the max x range to RKBESL validity
      real(BK), parameter :: xmax =  10.0_BK
      real(BK), parameter :: RTOL =  1e-6_BK
      real(BK), parameter :: ATOL =  1e-10_BK
      real(BK), allocatable, dimension(:) :: x,fun,intr,err,this
      integer :: i,ierr

      allocate(x(NTEST),fun(NTEST),intr(0:NTEST),err(NTEST),this(2))

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
        allocate(x(nsize),intrin(nsize),packge(nsize),z(ntest))

        call random_number(x)
        x    = xmin*(ONE-x) + xmax*x

        time = ZERO
        do i=1,ntest
            call cpu_time(c_start)

            do j=1,nsize
               CALL RKBESL(X=x(j), ALPHA=ZERO, NB=2, IZE=1, BK=this, NCALC=ierr)
               intrin(j) = this(2)
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

        success = timep<2*time

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

        integer, parameter :: nsize = 1000
        integer, parameter :: ntest = 100
        real(BK), parameter :: xmin =    0.0_BK
        ! Limit the max x range to RKBESL validity
        real(BK), parameter :: xmax =  600.0_BK
        real(BK), allocatable :: x(:),intrin(:),packge(:),z(:)
        integer :: i,j,ierr
        real(BK) :: time,timep,c_start,c_end
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

        success = timep<2*time

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


    ! Test bessel k0 function
    logical function test_gamma() result(success)

      integer, parameter :: NTEST = 2000

      real(BK), parameter :: xmin =  -12.0_BK

      ! Limit the max x range to RKBESL validity
      real(BK), parameter :: xmax =  12.0_BK
      real(BK), parameter :: RTOL =  1e-6_BK
      real(BK), parameter :: ATOL =  1e-10_BK
      real(BK) :: x(NTEST),fun(NTEST),intr(NTEST),err(NTEST),this(2)
      integer :: i,ierr

      ! Randoms in range
      call random_number(x)
      x    = xmin*(ONE-x) + xmax*x
      fun  = gamma_BK(x)
      intr = gamma(x)

      err  = abs(fun-intr(1:))*rewt(intr(1:),RTOL,ATOL)

      success = all(err<one)

      if (.not.success) then
         do i=1,NTEST
            if (err(i)>=one) &
            print *, '[bessel_gamma] x=',x(i),' package=',fun(i),' intrinsic=',intr(i),' relerr=',err(i)
         end do
      end if

    end function test_gamma


    ! Test bessel j0 cpu time
    logical function test_gamma_cputime() result(success)

        integer, parameter :: nsize = 100000
        integer, parameter :: ntest = 100
        real(BK), parameter :: xmin = -50.0_BK
        real(BK), parameter :: xmax =  50.0_BK
        real(BK), allocatable :: x(:),intrin(:),packge(:),z(:)
        integer :: i,j,ierr
        real(BK) :: time,timep,c_start,c_end
        allocate(x(nsize),intrin(nsize),packge(nsize),z(ntest))

        call random_number(x)
        x    = xmin*(ONE-x) + xmax*x

        time = ZERO
        do i=1,ntest
            call cpu_time(c_start)
            intrin = gamma(x)
            call cpu_time(c_end)
            z(i) = sum(intrin)
            time = time+c_end-c_start
        end do
        print "('[bessel_gamma] INTRINSIC    time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*time/(nsize*ntest),sum(z)

        timep = ZERO
        do i=1,ntest
            call cpu_time(c_start)
            packge = gamma_BK(x)
            call cpu_time(c_end)
            z(i) = sum(packge)
            timep = timep+c_end-c_start
        end do
        print "('[bessel_gamma] PACKAGE   time used: ',f9.4,' ns/eval, sum(z)=',g0)",1e9*timep/(nsize*ntest),sum(z)

        success = timep<2*time

    end function test_gamma_cputime

    ! Test bessel j0 cpu time
    logical function test_bessel_i1_cputime() result(success)
        use rji

        integer, parameter :: nsize = 100000
        integer, parameter :: ntest = 100
        real(BK), parameter :: xmin =    0.0_BK
        real(BK), parameter :: xmax =  100.0_BK
        real(BK), allocatable :: x(:),intrin(:),packge(:),z(:)
        integer :: i,j,ierr
        real(BK) :: time,timep,c_start,c_end
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

        success = timep<2*time

    end function test_bessel_i1_cputime

    ! ode-like inverse error weight
    elemental real(BK) function rewt(x,RTOL,ATOL)
       real(BK), intent(in) :: x,RTOL,ATOL
       rewt = ONE/(RTOL*abs(x)+ATOL)
    end function rewt

    ! Regression test for besselj_up_recurrence bug fix.
    ! Forward Y-recurrence from Y_0, Y_1 must reach the correct Y_n for integer n,
    ! with coefficient (2k/x) updated each step. Before the fix the coefficient was
    ! frozen at nu_start*2/x and the loop counter ran backwards.
    logical function test_besselj_up_recurrence() result(success)
       use bessels_constants, only: besselj_up_recurrence

       real(BK), parameter :: RTOL = 1e-10_BK
       real(BK), parameter :: ATOL = 1e-14_BK
       real(BK), parameter :: x_test(3) = [0.5_BK, 2.0_BK, 10.0_BK]
       real(BK) :: x, Y_pkg, Y_next, Y_intr, err
       integer  :: n, i

       success = .true.

       do i = 1, size(x_test)
          x = x_test(i)
          do n = 1, 10
             call besselj_up_recurrence(x, bessely1(x), bessely0(x), &
                                        ONE, real(n,BK), Y_pkg, Y_next)
             Y_intr = bessel_yn(n, x)
             err = abs(Y_pkg - Y_intr) * rewt(Y_intr, RTOL, ATOL)
             if (err >= ONE) then
                success = .false.
                print *, '[besselj_up_recurrence] x=', x, ' n=', n, &
                         ' package=', Y_pkg, ' intrinsic=', Y_intr, ' relerr=', err
             end if
          end do
       end do

    end function test_besselj_up_recurrence

    ! Test bessely(real nu, x) against bessel_yn intrinsic for integer nu.
    logical function test_bessely_nu_integer() result(success)

       integer, parameter :: NTEST = 1000
       real(BK), parameter :: xmin = 0.1_BK
       real(BK), parameter :: xmax = 1e+2_BK
       real(BK), parameter :: RTOL = 1e-6_BK
       real(BK), parameter :: ATOL = 1e-10_BK
       real(BK) :: x(NTEST), fun(NTEST), intr(NTEST), err(NTEST)
       integer  :: i, n

       call random_number(x)
       x = xmin*(ONE-x) + xmax*x

       success = .true.

       do n = 0, 10
          do i = 1, NTEST
             fun(i) = bessely(real(n,BK), x(i))
          end do
          intr = bessel_yn(n, x)
          err  = abs(fun-intr) * rewt(intr, RTOL, ATOL)

          success = success .and. all(err<one)

          if (.not. all(err<one)) then
             do i = 1, NTEST
                if (err(i) >= one) &
                   print *, '[bessely_nu_int] n=',n,' x=',x(i),' package=',fun(i),' intrinsic=',intr(i),' relerr=',err(i)
             end do
          end if
       end do

    end function test_bessely_nu_integer

    ! Test hankelh1(integer n, x) = J_n(x) + i*Y_n(x).
    ! Restricted to integer nu — non-integer nu paths in bessely_positive_args
    ! / hankel_debye have pre-existing port bugs (see todo/07-nonintegerNU-bugs.md).
    logical function test_hankelh1_consistency() result(success)

       real(BK), parameter :: RTOL = 1e-10_BK
       real(BK), parameter :: ATOL = 1e-14_BK
       real(BK), parameter :: x_test(4) = [0.5_BK, 5.0_BK, 30.0_BK, 100.0_BK]
       real(BK) :: x, Jref, Yref, err_r, err_i
       complex(BK) :: H
       integer :: i

       success = .true.

       do i = 1, size(x_test)
          x = x_test(i)

          ! nu = 0
          H    = hankelh1(0.0_BK, x)
          Jref = besselj0(x)
          Yref = bessely0(x)
          err_r = abs(H%re - Jref) * rewt(Jref, RTOL, ATOL)
          err_i = abs(H%im - Yref) * rewt(Yref, RTOL, ATOL)
          if (err_r >= ONE .or. err_i >= ONE) then
             success = .false.
             print *, '[hankelh1_cons] nu=0 x=', x, ' H=', H, ' J=', Jref, ' Y=', Yref, &
                      ' err_r=', err_r, ' err_i=', err_i
          end if

          ! nu = 1
          H    = hankelh1(1.0_BK, x)
          Jref = besselj1(x)
          Yref = bessely1(x)
          err_r = abs(H%re - Jref) * rewt(Jref, RTOL, ATOL)
          err_i = abs(H%im - Yref) * rewt(Yref, RTOL, ATOL)
          if (err_r >= ONE .or. err_i >= ONE) then
             success = .false.
             print *, '[hankelh1_cons] nu=1 x=', x, ' H=', H, ' J=', Jref, ' Y=', Yref, &
                      ' err_r=', err_r, ' err_i=', err_i
          end if
       end do

    end function test_hankelh1_consistency

    ! Test hankelh2(n, x) = conjg(hankelh1(n, x)) for real x > 0, integer nu.
    logical function test_hankelh2_conjugate() result(success)

       real(BK), parameter :: TOL = 1e-14_BK
       real(BK), parameter :: x_test(4) = [0.5_BK, 5.0_BK, 30.0_BK, 100.0_BK]
       real(BK), parameter :: nu_test(3) = [0.0_BK, 1.0_BK, 3.0_BK]
       complex(BK) :: H1, H2
       real(BK) :: diff
       integer :: i, j

       success = .true.

       do i = 1, size(x_test)
          do j = 1, size(nu_test)
             H1 = hankelh1(nu_test(j), x_test(i))
             H2 = hankelh2(nu_test(j), x_test(i))
             diff = abs(H2 - conjg(H1))
             if (diff > TOL) then
                success = .false.
                print *, '[hankelh2_conj] nu=', nu_test(j), ' x=', x_test(i), &
                         ' H1=', H1, ' H2=', H2, ' |H2-conjg(H1)|=', diff
             end if
          end do
       end do

    end function test_hankelh2_conjugate

    ! Benchmark bessely(real nu, x) vs. bessel_yn intrinsic at nu=1.
    logical function test_bessely_nu_cputime() result(success)

       integer, parameter :: nsize = 100000
       integer, parameter :: ntest = 100
       real(BK), parameter :: xmin = 0.1_BK
       real(BK), parameter :: xmax = 1e+3_BK
       real(BK), allocatable :: x(:), intrin(:), packge(:), z(:)
       integer :: i
       real(BK) :: time, timep, c_start, c_end
       allocate(x(nsize), intrin(nsize), packge(nsize), z(ntest))

       call random_number(x)
       x = xmin*(ONE-x) + xmax*x

       time = ZERO
       do i = 1, ntest
          call cpu_time(c_start)
          intrin = bessel_yn(1, x)
          call cpu_time(c_end)
          z(i) = sum(intrin)
          time = time + c_end - c_start
       end do
       print "('[bessely_nu] INTRINSIC time used: ',f9.4,' ns/eval, sum(z)=',g0)", 1e9*time/(nsize*ntest), sum(z)

       timep = ZERO
       do i = 1, ntest
          call cpu_time(c_start)
          packge = bessely(1.0_BK, x)
          call cpu_time(c_end)
          z(i) = sum(packge)
          timep = timep + c_end - c_start
       end do
       print "('[bessely_nu] PACKAGE   time used: ',f9.4,' ns/eval, sum(z)=',g0)", 1e9*timep/(nsize*ntest), sum(z)

       success = timep < 3*time

    end function test_bessely_nu_cputime

    ! Benchmark hankelh1(0, x) vs. assembling cmplx(besselj0, bessely0) by hand.
    logical function test_hankelh1_cputime() result(success)

       integer, parameter :: nsize = 100000
       integer, parameter :: ntest = 100
       real(BK), parameter :: xmin = 0.1_BK
       real(BK), parameter :: xmax = 1e+3_BK
       real(BK), allocatable :: x(:), z(:)
       complex(BK), allocatable :: baseline(:), packge(:)
       integer :: i
       real(BK) :: time, timep, c_start, c_end
       allocate(x(nsize), baseline(nsize), packge(nsize), z(ntest))

       call random_number(x)
       x = xmin*(ONE-x) + xmax*x

       time = ZERO
       do i = 1, ntest
          call cpu_time(c_start)
          baseline = cmplx(besselj0(x), bessely0(x), BK)
          call cpu_time(c_end)
          z(i) = sum(baseline%re)
          time = time + c_end - c_start
       end do
       print "('[hankelh1]  BASELINE  time used: ',f9.4,' ns/eval, sum(z)=',g0)", 1e9*time/(nsize*ntest), sum(z)

       timep = ZERO
       do i = 1, ntest
          call cpu_time(c_start)
          packge = hankelh1(0.0_BK, x)
          call cpu_time(c_end)
          z(i) = sum(packge%re)
          timep = timep + c_end - c_start
       end do
       print "('[hankelh1]  PACKAGE   time used: ',f9.4,' ns/eval, sum(z)=',g0)", 1e9*timep/(nsize*ntest), sum(z)

       success = timep < 3*time

    end function test_hankelh1_cputime

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
