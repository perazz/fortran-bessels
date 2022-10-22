program bessels_test
    use bessels
    implicit none


    integer :: nfailed,npassed

    nfailed = 0
    npassed = 0

    call add_test(test_cuberoot())

    print "('[bessels] ',i0,' test completed: ',i0,' passed, ',i0,' failed.')", npassed+nfailed,npassed,nfailed
    if (nfailed>0) then
        stop -1
    else
        stop 0
    endif

    contains

    subroutine add_test(success)
        logical, intent(in) :: success
        if (success) then
            npassed = npassed+1
        else
            nfailed = nfailed+1
        end if
    end subroutine add_test

    ! Test approximated cube root
    logical function test_cuberoot() result(success)

      integer, parameter :: NTEST = 1000

      real(BKIND), parameter :: xmin = -1e+6_BKIND
      real(BKIND), parameter :: xmax =  1e+6_BKIND
      real(BKIND), parameter :: RTOL =  1e-6_BKIND
      real(BKIND), parameter :: ATOL =  1e-20_BKIND
      real(BKIND) :: x(NTEST),fun(NTEST),intr(NTEST)
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
