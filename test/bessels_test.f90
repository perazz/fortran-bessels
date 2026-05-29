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
!  Copyright (c) 2022-2026 Federico Perini
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
    call add_test(test_bessely_nu_half_integer())
    call add_test(test_hankelh1_consistency())
    call add_test(test_hankelh2_conjugate())
    call add_test(test_hankelh_negative_nu())
    call add_test(test_hankel_debye_consistency())
    call add_test(test_airyai())
    call add_test(test_airyaiprime())
    call add_test(test_airybi())
    call add_test(test_airybiprime())
    call add_test(test_airyaix_overflow_guard())
    call add_test(test_airy_cputime())
    call add_test(test_besselk_integer())
    call add_test(test_besselk_noninteger())
    call add_test(test_besselkx_overflow_guard())
    call add_test(test_besselk_nu_cputime())
    call add_test(test_besseli_integer())
    call add_test(test_besseli_noninteger())
    call add_test(test_besseli_negative_nu_reflection())
    call add_test(test_besselix_overflow_guard())
    call add_test(test_besseli_nu_cputime())
    call add_test(test_sphericalbesselj_int())
    call add_test(test_sphericalbessely_int())
    call add_test(test_sphericalbesseli_int())
    call add_test(test_sphericalbesselk_int())
    call add_test(test_sphericalbessel_halfinteger_consistency())
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
      use bessels_rkbesl, only: RKBESL

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
        use bessels_rkbesl, only: RKBESL

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
      use bessels_rkbesl, only: RKBESL

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
        use bessels_rkbesl, only: RKBESL

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
      use bessels_ribesl, only: RIBESL

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
        use bessels_ribesl, only: RIBESL

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
      use bessels_ribesl, only: RIBESL

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
        use bessels_ribesl, only: RIBESL

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

       success = timep < 5*time

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

    ! Test bessely(half-integer nu, x) against closed-form expressions.
    !   Y_{1/2}(x)  = -sqrt(2/(pi x)) cos(x)
    !   Y_{3/2}(x)  = -sqrt(2/(pi x)) (cos(x)/x + sin(x))
    !   Y_{5/2}(x)  = -sqrt(2/(pi x)) ((3/x^2 - 1) cos(x) + 3 sin(x)/x)
    logical function test_bessely_nu_half_integer() result(success)
       use bessels_constants, only: TWOOPI, THREE

       real(BK), parameter :: RTOL = 1e-9_BK
       real(BK), parameter :: ATOL = 1e-13_BK
       real(BK), parameter :: x_test(5) = [0.5_BK, 2.0_BK, 7.0_BK, 30.0_BK, 100.0_BK]
       real(BK) :: x, sx, cx, c, y_half, y_3half, y_5half, ref, fun, err
       integer :: i

       success = .true.

       do i = 1, size(x_test)
          x  = x_test(i)
          sx = sin(x)
          cx = cos(x)
          c  = sqrt(TWOOPI/x)
          y_half  = -c * cx
          y_3half = -c * (cx/x + sx)
          y_5half = -c * ((THREE/x**2 - ONE)*cx + THREE*sx/x)

          ref = y_half ; fun = bessely(0.5_BK, x)
          err = abs(fun - ref) * rewt(ref, RTOL, ATOL)
          if (err >= ONE) then
             success = .false.
             print *, '[bessely_nu_half] nu=1/2 x=',x,' package=',fun,' ref=',ref,' relerr=',err
          end if

          ref = y_3half ; fun = bessely(1.5_BK, x)
          err = abs(fun - ref) * rewt(ref, RTOL, ATOL)
          if (err >= ONE) then
             success = .false.
             print *, '[bessely_nu_half] nu=3/2 x=',x,' package=',fun,' ref=',ref,' relerr=',err
          end if

          ref = y_5half ; fun = bessely(2.5_BK, x)
          err = abs(fun - ref) * rewt(ref, RTOL, ATOL)
          if (err >= ONE) then
             success = .false.
             print *, '[bessely_nu_half] nu=5/2 x=',x,' package=',fun,' ref=',ref,' relerr=',err
          end if
       end do

    end function test_bessely_nu_half_integer

    ! Test the Hankel reflection: H^(1)_{-nu}(x) = exp(+i*pi*nu) * H^(1)_{nu}(x)
    !                            H^(2)_{-nu}(x) = exp(-i*pi*nu) * H^(2)_{nu}(x)
    logical function test_hankelh_negative_nu() result(success)
       use bessels_constants, only: PI

       real(BK), parameter :: RTOL = 1e-10_BK
       real(BK), parameter :: ATOL = 1e-13_BK
       real(BK) :: x, nu, diff
       complex(BK) :: H_pos, H_neg, expected, refl
       integer :: i, j
       real(BK), parameter :: nu_test(2) = [0.5_BK, 2.5_BK]
       real(BK), parameter :: x_test(2)  = [5.0_BK, 20.0_BK]

       success = .true.

       do i = 1, size(x_test)
          do j = 1, size(nu_test)
             x  = x_test(i)
             nu = nu_test(j)
             H_pos = hankelh1(nu, x)
             H_neg = hankelh1(-nu, x)
             refl  = exp(cmplx(ZERO,  PI*nu, BK)) * H_pos
             diff  = abs(H_neg - refl) / (abs(refl) + ATOL)
             if (diff > RTOL) then
                success = .false.
                print *, '[hankel_neg_nu] H1 nu=',nu,' x=',x,' H_neg=',H_neg, &
                         ' ref=',refl,' diff=',diff
             end if

             H_pos = hankelh2(nu, x)
             H_neg = hankelh2(-nu, x)
             refl  = exp(cmplx(ZERO, -PI*nu, BK)) * H_pos
             diff  = abs(H_neg - refl) / (abs(refl) + ATOL)
             if (diff > RTOL) then
                success = .false.
                print *, '[hankel_neg_nu] H2 nu=',nu,' x=',x,' H_neg=',H_neg, &
                         ' ref=',refl,' diff=',diff
             end if
          end do
       end do

    end function test_hankelh_negative_nu

    ! Direct check that hankel_debye produces J + i*Y consistent with closed-form
    ! half-integer Hankel:
    !   H^(1)_{1/2}(x) = sqrt(2/(pi x)) * (sin(x) - i*cos(x))
    !   H^(1)_{3/2}(x) = -sqrt(2/(pi x)) * ((cos(x)/x + sin(x))*i + cos(x) - sin(x)/x)
    ! For integer nu, compare against intrinsic bessel_jn / bessel_yn.
    logical function test_hankel_debye_consistency() result(success)
       use bessels_debye, only: hankel_debye
       use bessels_constants, only: hankel_debye_cutoff, TWOOPI, THREE

       real(BK), parameter :: RTOL = 1e-9_BK
       real(BK), parameter :: ATOL = 1e-12_BK
       real(BK), parameter :: x_test(4) = [30.0_BK, 50.0_BK, 100.0_BK, 200.0_BK]
       real(BK) :: x, nu, errJ, errY, c, sx, cx, Jref, Yref
       complex(BK) :: H
       integer :: i, n

       success = .true.

       do i = 1, size(x_test)
          x = x_test(i)

          ! Half-integer nu = 1/2 — must be inside hankel_debye_cutoff region
          nu = 0.5_BK
          if (hankel_debye_cutoff(nu, x)) then
             sx = sin(x); cx = cos(x); c = sqrt(TWOOPI/x)
             Jref =  c * sx
             Yref = -c * cx
             H = hankel_debye(nu, x)
             errJ = abs(real(H,BK) - Jref) * rewt(Jref, RTOL, ATOL)
             errY = abs(aimag(H)   - Yref) * rewt(Yref, RTOL, ATOL)
             if (errJ >= ONE .or. errY >= ONE) then
                success = .false.
                print *, '[hankel_debye_cons] nu=',nu,' x=',x,' H=',H,' Jref=',Jref,' Yref=',Yref, &
                         ' errJ=',errJ,' errY=',errY
             end if
          end if

          ! Integer nu values
          do n = 0, 3
             nu = real(n, BK)
             if (hankel_debye_cutoff(nu, x)) then
                Jref = bessel_jn(n, x)
                Yref = bessel_yn(n, x)
                H = hankel_debye(nu, x)
                errJ = abs(real(H,BK) - Jref) * rewt(Jref, RTOL, ATOL)
                errY = abs(aimag(H)   - Yref) * rewt(Yref, RTOL, ATOL)
                if (errJ >= ONE .or. errY >= ONE) then
                   success = .false.
                   print *, '[hankel_debye_cons] n=',n,' x=',x,' H=',H,' Jref=',Jref,' Yref=',Yref, &
                            ' errJ=',errJ,' errY=',errY
                end if
             end if
          end do
       end do

    end function test_hankel_debye_consistency

    ! Test Airy values (Ai, Bi, Ai', Bi') against the Bessels.jl reference CSV files,
    ! which were generated in Mathematica/Arb to Float64 precision.
    !
    ! The negative-args CSV provides unscaled (Ai, Aip, Bi, Bip).  Tolerance scaling
    ! mirrors the Bessels.jl test suite:
    !   |x| <= 9.5  : tol = 2.4e-16 (relative)
    !   9.5 < |x| <= 1e8: tol = 0.8e-16 * |x|^(5/4) (absolute, for Ai/Bi);
    !                            0.8e-16 * |x|^(7/4) (absolute, for Ai'/Bi')
    logical function test_airyai() result(success)
       success = check_airy_neg('Ai', 1)
    end function test_airyai

    logical function test_airyaiprime() result(success)
       success = check_airy_neg('Aip', 2)
    end function test_airyaiprime

    logical function test_airybi() result(success)
       success = check_airy_neg('Bi', 3)
    end function test_airybi

    logical function test_airybiprime() result(success)
       success = check_airy_neg('Bip', 4)
    end function test_airybiprime

    ! Drives one of the four Airy functions across the Bessels.jl negative-arg CSV.
    ! sel = 1: Ai, 2: Ai', 3: Bi, 4: Bi'.  Returns success/failure and reports
    ! the worst offending row on failure.
    logical function check_airy_neg(label, sel) result(success)
       character(*), intent(in) :: label
       integer,      intent(in) :: sel

       integer :: u, ios, n_total, n_bad
       real(BK) :: x, ai_r, aip_r, bi_r, bip_r, ref, fun, tol, ax, err, worst_err, worst_x
       character(len=256) :: line

       success = .true.
       n_total = 0
       n_bad   = 0
       worst_err = ZERO
       worst_x   = ZERO

       open(newunit=u, file='test/data/airy/airy_negative_args.csv', status='old', &
            action='read', iostat=ios)
       if (ios /= 0) then
          ! Reference data not staged — skip test gracefully.
          print "('[airy_',a,'] SKIP: reference CSV not found')", trim(label)
          return
       end if
       do
          read(u, '(A)', iostat=ios) line
          if (ios /= 0) exit
          read(line, *, iostat=ios) x, ai_r, aip_r, bi_r, bip_r
          if (ios /= 0) cycle
          n_total = n_total + 1

          select case (sel)
            case (1); ref = ai_r;  fun = airyai(x)
            case (2); ref = aip_r; fun = airyaiprime(x)
            case (3); ref = bi_r;  fun = airybi(x)
            case (4); ref = bip_r; fun = airybiprime(x)
          end select

          ax = abs(x)
          if (ax <= 9.5_BK) then
             tol = max(2.4e-13_BK * abs(ref), 1e-14_BK)
          elseif (ax <= 1.0e8_BK) then
             if (sel == 1 .or. sel == 3) then
                tol = 0.8e-12_BK * ax**1.25_BK
             else
                tol = 0.8e-12_BK * ax**1.75_BK
             end if
          else
             cycle  ! beyond domain, package returns NaN — skip
          end if

          err = abs(fun - ref)
          if (err > tol) then
             n_bad = n_bad + 1
             if (err > worst_err) then
                worst_err = err
                worst_x   = x
             end if
          end if
       end do
       close(u)

       if (n_bad > 0) then
          success = .false.
          print "('[airy_',a,'] FAIL: ',i0,' of ',i0,' rows; worst x=',es15.8,' err=',es12.5)", &
                 trim(label), n_bad, n_total, worst_x, worst_err
       end if

    end function check_airy_neg

    ! airyaix(100) should be O(1) while airyai(100) underflows.
    ! Similarly airybix(100) should be O(1) while airybi(100) overflows.
    logical function test_airyaix_overflow_guard() result(success)

       real(BK) :: ai, aix, bi, bix
       real(BK), parameter :: x = 100.0_BK
       real(BK), parameter :: TOL = 5e-10_BK

       success = .true.

       ai  = airyai(x)
       aix = airyaix(x)
       if (.not. (ai >= ZERO .and. ai < 1e-280_BK)) then
          success = .false.
          print *, '[airyai_overflow] airyai(100) should be tiny: got ',ai
       end if
       if (.not. (aix > 1e-2_BK .and. aix < 1.0_BK)) then
          success = .false.
          print *, '[airyai_overflow] airyaix(100) should be O(1): got ',aix
       end if

       bi  = airybi(x)
       bix = airybix(x)
       if (.not. (bi > 1e280_BK)) then
          success = .false.
          print *, '[airybi_overflow] airybi(100) should be huge: got ',bi
       end if
       if (.not. (bix > 1e-2_BK .and. bix < 1.0_BK)) then
          success = .false.
          print *, '[airybi_overflow] airybix(100) should be O(1): got ',bix
       end if

       ! Verify they're consistent: airyaix(x) ~ airyai(x) * exp(2/3 * x^(3/2))
       ! For x=2, this is well-defined
       ai  = airyai(2.0_BK)
       aix = airyaix(2.0_BK)
       if (abs(aix - ai*exp(2.0_BK*2.0_BK*sqrt(2.0_BK)/3.0_BK)) > TOL) then
          success = .false.
          print *, '[airyai_overflow] inconsistency: airyaix(2) ',aix, ' vs scaled airyai(2) ',&
                   ai*exp(2.0_BK*2.0_BK*sqrt(2.0_BK)/3.0_BK)
       end if

    end function test_airyaix_overflow_guard

    ! Quick CPU benchmark for airyai (no intrinsic; compare against package only).
    logical function test_airy_cputime() result(success)

       integer, parameter :: nsize = 100000
       integer, parameter :: ntest = 100
       real(BK), parameter :: xmin = -10.0_BK
       real(BK), parameter :: xmax =  10.0_BK
       real(BK), allocatable :: x(:), pa(:), pb(:), z(:)
       integer :: i
       real(BK) :: time, c_start, c_end
       allocate(x(nsize), pa(nsize), pb(nsize), z(ntest))

       call random_number(x)
       x = xmin*(ONE-x) + xmax*x

       time = ZERO
       do i = 1, ntest
          call cpu_time(c_start)
          pa = airyai(x)
          pb = airybi(x)
          call cpu_time(c_end)
          z(i) = sum(pa) + sum(pb)
          time = time + c_end - c_start
       end do
       print "('[airy]      PACKAGE   time used: ',f9.4,' ns/eval (Ai+Bi pair), sum(z)=',g0)", &
             1e9*time/(nsize*ntest), sum(z)

       ! Always succeeds — informational only.
       success = .true.

    end function test_airy_cputime

    ! Test besselk(integer nu, x) against netlib's RKBESL across a representative grid.
    logical function test_besselk_integer() result(success)
       use bessels_rkbesl, only: RKBESL

       real(BK), parameter :: RTOL = 1e-8_BK
       real(BK), parameter :: ATOL = 1e-12_BK
       real(BK), parameter :: x_test(*) = [0.01_BK, 0.1_BK, 1.0_BK, 5.0_BK, 20.0_BK, 100.0_BK, 500.0_BK]
       real(BK) :: x, ref, fun, err
       real(BK), allocatable :: rk(:)
       integer  :: i, n, ierr

       success = .true.
       allocate(rk(120))

       do i = 1, size(x_test)
          x = x_test(i)
          call RKBESL(X=x, ALPHA=ZERO, NB=110, IZE=1, BK=rk, NCALC=ierr)
          if (ierr < 0) cycle
          do n = 0, min(ierr-1, 100)
             ref = rk(n+1)
             fun = besselk(real(n, BK), x)
             if (ref > 1e+250_BK) cycle    ! netlib overflowed
             if (ref < 1e-280_BK) cycle    ! underflowed reference
             err = abs(fun - ref) * rewt(ref, RTOL, ATOL)
             if (err >= ONE) then
                success = .false.
                print *, '[besselk_int] n=',n,' x=',x,' package=',fun,' ref=',ref,' relerr=',err
             end if
          end do
       end do

    end function test_besselk_integer

    ! Test besselk(non-integer nu, x) against netlib's RKBESL on the same grid.
    ! RKBESL accepts ALPHA in [0,1) and returns K_{alpha+n} for n=0..NB-1.
    logical function test_besselk_noninteger() result(success)
       use bessels_rkbesl, only: RKBESL

       real(BK), parameter :: RTOL = 1e-7_BK
       real(BK), parameter :: ATOL = 1e-12_BK
       real(BK), parameter :: x_test(*) = [0.1_BK, 1.0_BK, 5.0_BK, 20.0_BK, 100.0_BK]
       real(BK), parameter :: alpha_test(*) = [0.25_BK, 0.5_BK, 0.7_BK]
       real(BK) :: x, alpha, ref, fun, err
       real(BK), allocatable :: rk(:)
       integer  :: i, j, n, ierr

       success = .true.
       allocate(rk(120))

       do i = 1, size(x_test)
          do j = 1, size(alpha_test)
             x = x_test(i)
             alpha = alpha_test(j)
             call RKBESL(X=x, ALPHA=alpha, NB=80, IZE=1, BK=rk, NCALC=ierr)
             if (ierr < 0) cycle
             do n = 0, min(ierr-1, 50)
                ref = rk(n+1)
                fun = besselk(alpha + real(n, BK), x)
                if (ref > 1e+250_BK) cycle
                if (ref < 1e-280_BK) cycle
                err = abs(fun - ref) * rewt(ref, RTOL, ATOL)
                if (err >= ONE) then
                   success = .false.
                   print *, '[besselk_nu] nu=',alpha+n,' x=',x,' package=',fun,' ref=',ref,' relerr=',err
                end if
             end do
          end do
       end do

    end function test_besselk_noninteger

    ! besselkx(0, 700) should be finite (~sqrt(pi/1400)) while besselk(0, 700) underflows to 0.
    logical function test_besselkx_overflow_guard() result(success)
       real(BK) :: a, b
       success = .true.
       a = besselk(0.0_BK, 700.0_BK)
       b = besselkx(0.0_BK, 700.0_BK)
       if (.not. (a >= ZERO .and. a < 1e-280_BK)) then
          success = .false.
          print *, '[besselkx_guard] besselk(0,700) expected to underflow: got ',a
       end if
       if (.not. (b > 0.02_BK .and. b < 1.0_BK)) then
          success = .false.
          print *, '[besselkx_guard] besselkx(0,700) expected ~ sqrt(pi/1400): got ',b
       end if
    end function test_besselkx_overflow_guard

    ! Benchmark besselk(0.5, x) vs. netlib RKBESL.
    logical function test_besselk_nu_cputime() result(success)
       use bessels_rkbesl, only: RKBESL

       integer, parameter :: nsize = 50000
       integer, parameter :: ntest = 50
       real(BK), parameter :: xmin = 0.5_BK
       real(BK), parameter :: xmax = 30.0_BK
       real(BK), allocatable :: x(:), pa(:), z(:)
       real(BK) :: this(2)
       integer :: i, j, ierr
       real(BK) :: time, timep, c_start, c_end
       allocate(x(nsize), pa(nsize), z(ntest))

       call random_number(x)
       x = xmin*(ONE-x) + xmax*x

       time = ZERO
       do i = 1, ntest
          call cpu_time(c_start)
          do j = 1, nsize
             call RKBESL(X=x(j), ALPHA=0.5_BK, NB=1, IZE=1, BK=this, NCALC=ierr)
             pa(j) = this(1)
          end do
          call cpu_time(c_end)
          z(i) = sum(pa)
          time = time + c_end - c_start
       end do
       print "('[besselk_nu] NETLIB    time used: ',f9.4,' ns/eval, sum(z)=',g0)", 1e9*time/(nsize*ntest), sum(z)

       timep = ZERO
       do i = 1, ntest
          call cpu_time(c_start)
          pa = besselk(0.5_BK, x)
          call cpu_time(c_end)
          z(i) = sum(pa)
          timep = timep + c_end - c_start
       end do
       print "('[besselk_nu] PACKAGE   time used: ',f9.4,' ns/eval, sum(z)=',g0)", 1e9*timep/(nsize*ntest), sum(z)

       success = timep < 3*time + 1e-6_BK
    end function test_besselk_nu_cputime

    ! Test besseli(integer nu, x) against netlib's RIBESL on a representative grid.
    logical function test_besseli_integer() result(success)
       use bessels_ribesl, only: RIBESL

       real(BK), parameter :: RTOL = 1e-7_BK
       real(BK), parameter :: ATOL = 1e-12_BK
       real(BK), parameter :: x_test(*) = [0.01_BK, 0.1_BK, 1.0_BK, 5.0_BK, 20.0_BK, 50.0_BK, 100.0_BK]
       real(BK) :: x, ref, fun, err
       real(BK), allocatable :: ri(:)
       integer  :: i, n, ierr

       success = .true.
       allocate(ri(120))

       do i = 1, size(x_test)
          x = x_test(i)
          call RIBESL(X=x, ALPHA=ZERO, NB=80, IZE=1, B=ri, NCALC=ierr)
          if (ierr < 0) cycle
          do n = 0, min(ierr-1, 60)
             ref = ri(n+1)
             fun = besseli(real(n, BK), x)
             if (ref > 1e+250_BK) cycle
             if (abs(ref) < 1e-280_BK) cycle
             err = abs(fun - ref) * rewt(ref, RTOL, ATOL)
             if (err >= ONE) then
                success = .false.
                print *, '[besseli_int] n=',n,' x=',x,' package=',fun,' ref=',ref,' relerr=',err
             end if
          end do
       end do

    end function test_besseli_integer

    ! Test besseli(non-integer nu, x).
    logical function test_besseli_noninteger() result(success)
       use bessels_ribesl, only: RIBESL

       real(BK), parameter :: RTOL = 1e-6_BK
       real(BK), parameter :: ATOL = 1e-12_BK
       real(BK), parameter :: x_test(*) = [0.1_BK, 1.0_BK, 5.0_BK, 20.0_BK, 50.0_BK]
       real(BK), parameter :: alpha_test(*) = [0.25_BK, 0.5_BK, 0.7_BK]
       real(BK) :: x, alpha, ref, fun, err
       real(BK), allocatable :: ri(:)
       integer  :: i, j, n, ierr

       success = .true.
       allocate(ri(120))

       do i = 1, size(x_test)
          do j = 1, size(alpha_test)
             x = x_test(i)
             alpha = alpha_test(j)
             call RIBESL(X=x, ALPHA=alpha, NB=60, IZE=1, B=ri, NCALC=ierr)
             if (ierr < 0) cycle
             do n = 0, min(ierr-1, 40)
                ref = ri(n+1)
                fun = besseli(alpha + real(n, BK), x)
                if (ref > 1e+250_BK) cycle
                if (abs(ref) < 1e-280_BK) cycle
                err = abs(fun - ref) * rewt(ref, RTOL, ATOL)
                if (err >= ONE) then
                   success = .false.
                   print *, '[besseli_nu] nu=',alpha+n,' x=',x,' package=',fun,' ref=',ref,' relerr=',err
                end if
             end do
          end do
       end do

    end function test_besseli_noninteger

    ! Verify I_{-nu}(x) = I_nu(x) + (2/pi) sin(pi nu) K_nu(x).
    logical function test_besseli_negative_nu_reflection() result(success)
       use bessels_constants, only: PI, TWOOPI

       real(BK), parameter :: RTOL = 1e-9_BK
       real(BK), parameter :: ATOL = 1e-13_BK
       real(BK), parameter :: nu_test(*) = [0.25_BK, 0.5_BK, 2.5_BK]
       real(BK), parameter :: x_test(*) = [0.5_BK, 2.0_BK, 5.0_BK, 20.0_BK]
       real(BK) :: nu, x, lhs, rhs, err
       integer :: i, j

       success = .true.

       do i = 1, size(nu_test)
          do j = 1, size(x_test)
             nu = nu_test(i)
             x  = x_test(j)
             lhs = besseli(-nu, x)
             rhs = besseli( nu, x) + TWOOPI*sin(PI*nu)*besselk(nu, x)
             err = abs(lhs - rhs) * rewt(rhs, RTOL, ATOL)
             if (err >= ONE) then
                success = .false.
                print *, '[besseli_neg_nu] nu=',nu,' x=',x,' lhs=',lhs,' rhs=',rhs,' relerr=',err
             end if
          end do
       end do

    end function test_besseli_negative_nu_reflection

    ! besselix(0, 700) should be finite (~0.015) while besseli(0, 700) overflows.
    logical function test_besselix_overflow_guard() result(success)
       real(BK) :: a, b
       success = .true.
       a = besseli(0.0_BK, 700.0_BK)
       b = besselix(0.0_BK, 700.0_BK)
       if (.not. (a > 1e280_BK)) then
          success = .false.
          print *, '[besselix_guard] besseli(0,700) expected to overflow to large value: got ',a
       end if
       if (.not. (b > 0.005_BK .and. b < 0.1_BK)) then
          success = .false.
          print *, '[besselix_guard] besselix(0,700) expected ~0.015: got ',b
       end if
    end function test_besselix_overflow_guard

    ! Benchmark besseli(0.5, x) vs. netlib RIBESL.
    logical function test_besseli_nu_cputime() result(success)
       use bessels_ribesl, only: RIBESL

       integer, parameter :: nsize = 5000
       integer, parameter :: ntest = 50
       real(BK), parameter :: xmin = 0.1_BK
       real(BK), parameter :: xmax = 30.0_BK
       real(BK), allocatable :: x(:), pa(:), z(:)
       real(BK) :: this(2)
       integer :: i, j, ierr
       real(BK) :: time, timep, c_start, c_end
       allocate(x(nsize), pa(nsize), z(ntest))

       call random_number(x)
       x = xmin*(ONE-x) + xmax*x

       time = ZERO
       do i = 1, ntest
          call cpu_time(c_start)
          do j = 1, nsize
             call RIBESL(X=x(j), ALPHA=0.5_BK, NB=1, IZE=1, B=this, NCALC=ierr)
             pa(j) = this(1)
          end do
          call cpu_time(c_end)
          z(i) = sum(pa)
          time = time + c_end - c_start
       end do
       print "('[besseli_nu] NETLIB    time used: ',f9.4,' ns/eval, sum(z)=',g0)", 1e9*time/(nsize*ntest), sum(z)

       timep = ZERO
       do i = 1, ntest
          call cpu_time(c_start)
          pa = besseli(0.5_BK, x)
          call cpu_time(c_end)
          z(i) = sum(pa)
          timep = timep + c_end - c_start
       end do
       print "('[besseli_nu] PACKAGE   time used: ',f9.4,' ns/eval, sum(z)=',g0)", 1e9*timep/(nsize*ntest), sum(z)

       success = .true.
    end function test_besseli_nu_cputime

    ! Spherical Bessel j_n vs. closed-form low-order expressions and recurrence consistency.
    !   j_0(x) = sin(x)/x
    !   j_1(x) = sin(x)/x^2 - cos(x)/x
    !   j_2(x) = (3/x^2 - 1) sin(x)/x - 3 cos(x)/x^2
    logical function test_sphericalbesselj_int() result(success)
       real(BK), parameter :: RTOL = 1e-9_BK
       real(BK), parameter :: ATOL = 1e-13_BK
       real(BK), parameter :: x_test(*) = [0.5_BK, 1.0_BK, 5.0_BK, 10.0_BK, 50.0_BK]
       real(BK) :: x, ref, fun, err, sx, cx
       integer :: i

       success = .true.

       do i = 1, size(x_test)
          x  = x_test(i)
          sx = sin(x); cx = cos(x)

          ref = sx/x;                 fun = sphericalbesselj(0.0_BK, x)
          err = abs(fun-ref)*rewt(ref, RTOL, ATOL)
          if (err >= ONE) then
             success = .false.
             print *, '[sj0] x=',x,' package=',fun,' ref=',ref,' relerr=',err
          end if

          ref = sx/x**2 - cx/x;       fun = sphericalbesselj(1.0_BK, x)
          err = abs(fun-ref)*rewt(ref, RTOL, ATOL)
          if (err >= ONE) then
             success = .false.
             print *, '[sj1] x=',x,' package=',fun,' ref=',ref,' relerr=',err
          end if

          ref = (3.0_BK/x**2 - ONE)*sx/x - 3.0_BK*cx/x**2
          fun = sphericalbesselj(2.0_BK, x)
          err = abs(fun-ref)*rewt(ref, RTOL, ATOL)
          if (err >= ONE) then
             success = .false.
             print *, '[sj2] x=',x,' package=',fun,' ref=',ref,' relerr=',err
          end if
       end do

    end function test_sphericalbesselj_int

    ! Spherical Bessel y_n vs. closed-form low-order expressions.
    !   y_0(x) = -cos(x)/x
    !   y_1(x) = -cos(x)/x^2 - sin(x)/x
    !   y_2(x) = -(3/x^2 - 1) cos(x)/x - 3 sin(x)/x^2
    logical function test_sphericalbessely_int() result(success)
       real(BK), parameter :: RTOL = 1e-9_BK
       real(BK), parameter :: ATOL = 1e-13_BK
       real(BK), parameter :: x_test(*) = [0.5_BK, 1.0_BK, 5.0_BK, 10.0_BK, 50.0_BK]
       real(BK) :: x, ref, fun, err, sx, cx
       integer :: i

       success = .true.

       do i = 1, size(x_test)
          x  = x_test(i)
          sx = sin(x); cx = cos(x)

          ref = -cx/x;                fun = sphericalbessely(0.0_BK, x)
          err = abs(fun-ref)*rewt(ref, RTOL, ATOL)
          if (err >= ONE) then
             success = .false.
             print *, '[sy0] x=',x,' package=',fun,' ref=',ref,' relerr=',err
          end if

          ref = -cx/x**2 - sx/x;      fun = sphericalbessely(1.0_BK, x)
          err = abs(fun-ref)*rewt(ref, RTOL, ATOL)
          if (err >= ONE) then
             success = .false.
             print *, '[sy1] x=',x,' package=',fun,' ref=',ref,' relerr=',err
          end if

          ref = -(3.0_BK/x**2 - ONE)*cx/x - 3.0_BK*sx/x**2
          fun = sphericalbessely(2.0_BK, x)
          err = abs(fun-ref)*rewt(ref, RTOL, ATOL)
          if (err >= ONE) then
             success = .false.
             print *, '[sy2] x=',x,' package=',fun,' ref=',ref,' relerr=',err
          end if
       end do

    end function test_sphericalbessely_int

    ! Spherical modified Bessel i_n vs. closed-form low-order expressions.
    !   i_0(x) = sinh(x)/x
    !   i_1(x) = (x cosh(x) - sinh(x))/x^2
    !   i_2(x) = (x^2 sinh + 3*(sinh - x cosh))/x^3
    logical function test_sphericalbesseli_int() result(success)
       real(BK), parameter :: RTOL = 1e-9_BK
       real(BK), parameter :: ATOL = 1e-13_BK
       real(BK), parameter :: x_test(*) = [0.1_BK, 0.5_BK, 1.0_BK, 5.0_BK, 20.0_BK]
       real(BK) :: x, ref, fun, err, sx, cx, x2
       integer :: i

       success = .true.

       do i = 1, size(x_test)
          x  = x_test(i)
          sx = sinh(x); cx = cosh(x); x2 = x*x

          ref = sx/x;                 fun = sphericalbesseli(0.0_BK, x)
          err = abs(fun-ref)*rewt(ref, RTOL, ATOL)
          if (err >= ONE) then
             success = .false.
             print *, '[si0] x=',x,' package=',fun,' ref=',ref,' relerr=',err
          end if

          ref = (x*cx - sx)/x2;       fun = sphericalbesseli(1.0_BK, x)
          err = abs(fun-ref)*rewt(ref, RTOL, ATOL)
          if (err >= ONE) then
             success = .false.
             print *, '[si1] x=',x,' package=',fun,' ref=',ref,' relerr=',err
          end if

          ref = (x2*sx + 3.0_BK*(sx - x*cx))/(x2*x)
          fun = sphericalbesseli(2.0_BK, x)
          err = abs(fun-ref)*rewt(ref, RTOL, ATOL)
          if (err >= ONE) then
             success = .false.
             print *, '[si2] x=',x,' package=',fun,' ref=',ref,' relerr=',err
          end if
       end do

    end function test_sphericalbesseli_int

    ! Spherical modified Bessel k_n vs. closed-form low-order expressions.
    !   k_0(x) = e^{-x}/x
    !   k_1(x) = (1 + 1/x) e^{-x}/x
    !   k_2(x) = (1 + 3/x + 3/x^2) e^{-x}/x
    logical function test_sphericalbesselk_int() result(success)
       real(BK), parameter :: RTOL = 1e-9_BK
       real(BK), parameter :: ATOL = 1e-13_BK
       real(BK), parameter :: x_test(*) = [0.5_BK, 1.0_BK, 5.0_BK, 20.0_BK]
       real(BK) :: x, ref, fun, err, ex
       integer :: i, n

       success = .true.

       do i = 1, size(x_test)
          x  = x_test(i)
          ex = exp(-x)

          ref = ex/x;                 fun = sphericalbesselk(0.0_BK, x)
          err = abs(fun-ref)*rewt(ref, RTOL, ATOL)
          if (err >= ONE) then
             success = .false.
             print *, '[sk0] x=',x,' package=',fun,' ref=',ref,' relerr=',err
          end if

          ref = (ONE + ONE/x)*ex/x;   fun = sphericalbesselk(1.0_BK, x)
          err = abs(fun-ref)*rewt(ref, RTOL, ATOL)
          if (err >= ONE) then
             success = .false.
             print *, '[sk1] x=',x,' package=',fun,' ref=',ref,' relerr=',err
          end if

          ref = (ONE + 3.0_BK/x + 3.0_BK/(x*x))*ex/x
          fun = sphericalbesselk(2.0_BK, x)
          err = abs(fun-ref)*rewt(ref, RTOL, ATOL)
          if (err >= ONE) then
             success = .false.
             print *, '[sk2] x=',x,' package=',fun,' ref=',ref,' relerr=',err
          end if

          ! Check symmetry k_{-n}(x) = k_{n-1}(x) for a couple of integer n.
          do n = 1, 4
             if (abs(sphericalbesselk(-real(n,BK), x) - sphericalbesselk(real(n-1,BK), x)) &
                 > 1e-12_BK * abs(sphericalbesselk(real(n-1,BK), x))) then
                success = .false.
                print *, '[sk_sym] n=',n,' x=',x,' k_{-n}=',sphericalbesselk(-real(n,BK), x), &
                         ' k_{n-1}=',sphericalbesselk(real(n-1,BK), x)
             end if
          end do
       end do

    end function test_sphericalbesselk_int

    ! Verify the half-integer reduction: sphericalbesselj(5.5, x) = sqrt(pi/(2x)) * besselj(6.0, x).
    logical function test_sphericalbessel_halfinteger_consistency() result(success)
       use bessels_constants, only: PI

       real(BK), parameter :: RTOL = 1e-9_BK
       real(BK), parameter :: ATOL = 1e-13_BK
       real(BK), parameter :: x_test(*) = [1.0_BK, 5.0_BK, 20.0_BK]
       real(BK), parameter :: nu_test(*) = [0.5_BK, 1.5_BK, 5.5_BK]
       real(BK) :: x, nu, j_sphere, j_redux, err
       integer :: i, j

       success = .true.

       do i = 1, size(x_test)
          do j = 1, size(nu_test)
             x  = x_test(i)
             nu = nu_test(j)
             ! For nu = m + 0.5, besselj(nu + 0.5, x) = besselj(m + 1, x), an integer call.
             j_sphere = sphericalbesselj(nu, x)
             j_redux  = sqrt(PI/(2.0_BK*x)) * besseljn(int(nu + 0.5_BK + 0.5_BK), x)
             err = abs(j_sphere - j_redux) * rewt(j_redux, RTOL, ATOL)
             if (err >= ONE) then
                success = .false.
                print *, '[sphj_consistency] nu=',nu,' x=',x,' sphere=',j_sphere, &
                         ' redux=',j_redux,' relerr=',err
             end if
          end do
       end do

    end function test_sphericalbessel_halfinteger_consistency

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
