!  ************************************************************************************************************
!
!                                    ____  ________________ ________   _____
!                                   / __ )/ ____/ ___/ ___// ____/ /  / ___/
!                                  / __  / __/  \__ \\__ \/ __/ / /   \__ \
!                                 / /_/ / /___ ___/ /__/ / /___/ /______/ /
!                                /_____/_____//____/____/_____/_____/____/
!
!            A Modern Fortran port of the Bessels.jl library for fast Bessel function calculations
!
!  MIT License
!
!  Copyright (c) 2022 Federico Perini
!  Copyright (c) 2021-2022 Michael Helton, Oscar Smith, and the Bessels.jl contributors
!
!  ************************************************************************************************************
module bessels
    use bessels_constants
    use bessels_debye

    implicit none
    private

    ! Todo: make one module per real precision
    public :: BK,BSIZE

    public :: besseli0,besseli1
    public :: besselj0,besselj1,besseljn
    public :: bessely0,bessely1
    public :: besselk0,besselk1


    public :: cbrt
    public :: ZERO,ONE,THIRD


    contains

    ! Calculation of besselj0 is done in three branches using polynomial approximations
    !
    ! Branch 1: x <= pi/2
    !           besselj0 is calculated using a 9 term, even minimax polynomial
    !
    ! Branch 2: pi/2 < x < 26.0
    !           besselj0 is calculated by one of 16 different degree 13 minimax polynomials
    !           Each polynomial is an expansion around either a root or extrema of the besselj0.
    !           This ensures accuracy near the roots. Method taken from [2]
    !
    ! Branch 3: x >= 26.0
    !           besselj0 = sqrt(2/(pi*x))*beta(x)*(cos(x - pi/4 - alpha(x))
    !           See modified expansions given in [2]. Exact coefficients are used.

    elemental real(BK) function besselj0(x)
       real(BK), intent(in) :: x

       real(BK) :: ax,xinv,x2,xn,p,q,a,b,r
       intrinsic :: sqrt
       integer :: n

       real(BK), parameter :: ppoly(*) = [ONE, -1.0_BK/16.0_BK, 53.0_BK/512.0_BK, -4447.0_BK/8192.0_BK, 3066403.0_BK/524288.0_BK, &
                                          -896631415.0_BK/8388608.0_BK, 796754802993.0_BK/268435456.0_BK, &
                                          -500528959023471.0_BK/4294967296.0_BK]
       real(BK), parameter :: qpoly(*) = [-1.0_BK/8.0_BK, 25.0_BK/384.0_BK, -1073.0_BK/5120.0_BK, 375733.0_BK/229376.0_BK, &
                                          -55384775.0_BK/2359296.0_BK, 24713030909.0_BK/46137344.0_BK, &
                                          -7780757249041.0_BK/436207616.0_BK]

       ax = abs(x)

       if (ax <= PIO2) then

          besselj0 = evalpoly(size(J0_POLY_PIO2), ax**2, J0_POLY_PIO2)

       elseif (ax < 26.0_BK) then

           n        = int(TWOOPI*ax) ! 1 < n < 16
           r        = ax - sum(J0_ROOTS(:,n))
           besselj0 = evalpoly(size(J0_POLYS,1), r, J0_POLYS(:,n))

       elseif (ax<huge(ax)) then

           xinv = ONE/ax
           x2   = xinv**2

           ! Cut to 5-th order when we know we'll have enough accuracy
           if (ax < 120.0_BK) then

               p = evalpoly8(x2, ppoly)
               q = evalpoly7(x2, qpoly)

           else
               p = evalpoly4(x2, ppoly(1:4))
               q = evalpoly4(x2, qpoly(1:4))

           endif

           a  = SQ2OPI * sqrt(xinv) * p
           xn = muladd(xinv, q, -PIO4)

           ! the following computes b = cos(x + xn) more accurately see src/misc.jl
           b = cos(ax+xn)

           besselj0 = a*b

       else

           besselj0 = ZERO

       end if

    end function besselj0

    !
    ! Calculation of besselj1 is done in a similar way as besselj0.
    ! For details on similarities:
    ! Harrison, John. "Fast and accurate Bessel function computation." 2009 19th IEEE Symposium on Computer
    !                 Arithmetic. IEEE, 2009.

    elemental real(BK) function besselj1(x)
       real(BK), intent(in) :: x

       real(BK) :: ax,xinv,x2,xn,p,q,a,b,r,s
       intrinsic :: sqrt
       integer :: n

       real(BK), parameter :: ppoly(8) = [ONE, 3.0_BK/16.0_BK, &
                                          -99.0_BK/512.0_BK, 6597.0_BK/8192.0_BK, &
                                          -4057965.0_BK/524288.0_BK, 1113686901.0_BK/8388608.0_BK, &
                                          -951148335159.0_BK/268435456.0_BK, 581513783771781.0_BK/4294967296.0_BK]

       real(BK), parameter :: qpoly(7) = [3.0_BK/8.0_BK, -21.0_BK/128.0_BK, 1899.0_BK/5120.0_BK, &
                                          -543483.0_BK/229376.0_BK, 8027901.0_BK/262144.0_BK, &
                                          -30413055339.0_BK/46137344.0_BK, 9228545313147.0_BK/436207616.0_BK]

       ax = abs(x)
       s  = sign(ONE,x)

       if (ax<=PIO2) then

           besselj1 = ax*s*evalpoly(size(J1_POLY_PIO2), ax**2, J1_POLY_PIO2)

       elseif (ax <= 26.0_BK) then

           n        = int(TWOOPI*ax) ! 1 < n < 16
           r        = ax - sum(J1_ROOTS(:,n))
           besselj1 = s*evalpoly(size(J1_POLYS,1), r, J1_POLYS(:,n))

       elseif (ax<huge(ax)) then

           xinv = ONE/ax
           x2   = xinv**2

           ! Cut to 5-th order when we know we'll have enough accuracy
           if (ax < 120.0_BK) then

               p = evalpoly8(x2, ppoly)
               q = evalpoly7(x2, qpoly)

           else

               p = evalpoly4(x2, ppoly(1:4))
               q = evalpoly4(x2, qpoly(1:4))

           endif

           a  = SQ2OPI * sqrt(xinv) * p
           xn = muladd(xinv, q, -3*PIO4)

           ! the following computes b = cos(x + xn) more accurately see src/misc.jl
           b = cos(ax+xn)

           besselj1 = s*a*b

       else

           besselj1 = ZERO

       end if


    end function besselj1


    elemental real(BK) function besseljn(nu, x)
        integer, intent(in) :: nu
        real(BK), intent(in) :: x

        real(BK) :: ax,cpi,spi,Jnu,Ynu,ranu
        integer :: sgn,anu

        anu  = abs(nu)
        ranu = real(anu,BK)
        ax   = abs(x)

        ! +1 if even; this could be faster!
        sgn = merge(1,-1,mod(anu,2)==0)

        Jnu = besselj_positive_args(ranu, ax)
        if (nu >= ZERO) then

            besseljn = merge(Jnu,Jnu*sgn,x>=ZERO)
        else
            if (x >= ZERO) then

                besseljn = Jnu *sgn

            else

                Ynu = bessely_positive_args(ranu, ax)
                spi = sin(PI*anu)
                cpi = cos(PI*anu)
                besseljn = (cpi*Jnu - spi*Ynu) * sgn

            endif
        endif
    end function besseljn

    ! Recurrence J_{nu}(x)

    ! At this point we must fill the region when x Å v with recurrence
    ! Backward recurrence is always stable and forward recurrence is stable when x > nu
    ! However, we only use backward recurrence by shifting the order up and using `besseljy_debye` to
    ! generate start values. Both `besseljy_debye` and `hankel_debye` get more accurate for large orders,
    ! however `besseljy_debye` is slightly more efficient (no complex variables) and we need no branches
    ! if only consider one direction. On the other hand, shifting the order down avoids any concern about
    ! underflow for large orders.
    ! Shifting the order too high while keeping x fixed could result in numerical underflow
    ! Therefore we need to shift up only until the `besseljy_debye` is accurate and need to test that no
    ! underflow occurs. Shifting the order up decreases the value substantially for high orders and results
    ! in a stable forward recurrence as the values rapidly increase
    elemental real(BK) function besselj_recurrence(nu, x)
        real(BK), intent(in) :: nu, x

        real(BK) :: nu_shift,v,jnu,jnup1,dummy,jnu2(2)

        ! shift order up to where expansions are valid see U_polynomials.jl
        nu_shift = ceiling(besseljy_debye_fit(x)) - floor(nu)

        v = nu + nu_shift

        ! compute jnu and jnup1 then use downard recurrence starting from order v down to nu
        call besseljy_debye(v    , x, jnu  ,dummy);
        call besseljy_debye(v+ONE, x, jnup1,dummy);

        jnu2 = besselj_down_recurrence(x,jnu,jnup1,v,nu)
        besselj_recurrence = jnu2(1)

    end function besselj_recurrence

    ! Bessel function of the first kind of order nu, ``J_{nu}(x)``.
    ! nu and x must be real and nu and x must be positive.
    ! No checks on arguments are performed and should only be called if certain nu, x >= 0.
    elemental real(BK) function besselj_positive_args(nu, x)
        real(BK), intent(in) :: nu, x

        real(BK) :: besselj,bessely

        if (nu==ZERO) then

            besselj_positive_args = besselj0(x)

        elseif (nu==ONE) then

            besselj_positive_args = besselj1(x)

        elseif (besseljy_debye_cutoff(nu, x)) then

            ! x < ~nu branch see src/U_polynomials.jl
            call besseljy_debye(nu, x, besselj,bessely)
            besselj_positive_args = besselj

        elseif (besseljy_large_argument_cutoff(nu, x)) then

            ! large argument branch see src/asymptotics.jl
            call besseljy_large_argument(nu, x, besselj, bessely)
            besselj_positive_args = besselj

        elseif (hankel_debye_cutoff(nu, x)) then

            ! x > ~nu branch see src/U_polynomials.jl on computing Hankel function
            besselj_positive_args = real(hankel_debye(nu, x), BK)

        elseif (besselj_series_cutoff(nu, x)) then

            ! use power series for small x and for when nu > x
            besselj_positive_args = besselj_power_series(nu, x)

        else

            ! shift nu up and use downward recurrence
            besselj_positive_args = besselj_recurrence(nu, x)

        end if

    end function besselj_positive_args

    ! Fallback for Y_{nu}(x)
    elemental real(BK) function bessely_fallback(nu, x)
        real(BK), intent(in) :: nu, x

        complex(BK) :: cheb
        real(BK) :: nu_shift,v2

        ! for x in (6, 19) we use Chebyshev approximation and forward recurrence
        if (besseljy_chebyshev_cutoff(x)) then

            cheb = bessely_chebyshev(nu, x)

        else

            ! at this point x > 19.0 (for Float64) and fairly close to nu
            ! shift nu down and use the debye expansion for Hankel function (valid x > nu) then
            ! use forward recurrence
            nu_shift = ceiling(nu) - floor(hankel_debye_fit(x)) + FOUR
            v2 = max(nu - nu_shift, ONE + nu-int(nu))
            call besselj_up_recurrence(x, aimag(hankel_debye(v2, x)), aimag(hankel_debye(v2 - 1, x)), v2, nu, cheb%re, cheb%im)

        endif

        bessely_fallback = cheb%re

    end function bessely_fallback

     ! Chebyshev approximation for Y_{nu}(x)
     ! Computes ``Y_{nu}(x)`` for medium arguments x in (6, 19) for any positive order using a
     ! Chebyshev approximation. Forward recurrence is used to fill orders starting at low orders nu in (0, 2).
     elemental complex(BK) function bessely_chebyshev(nu, x) result(cheb)
        real(BK), intent(in) :: nu, x

        real(BK) :: nu_floor
        complex(BK) :: Y

        nu_floor = nu - int(nu)

        Y = bessely_chebyshev_low_orders(nu_floor, x)

        call besselj_up_recurrence(x, Y%im, Y%re, nu_floor + ONE, nu, cheb%re, cheb%im)

     end function bessely_chebyshev

     ! only implemented for Float64 so far
     elemental logical function besseljy_chebyshev_cutoff(x)
        real(BK), intent(in) :: x
        besseljy_chebyshev_cutoff = x <= 19.0_BK .and. x >= 6.0_BK
     end function besseljy_chebyshev_cutoff

     ! compute bessely for x in (6, 19) and nu in (0, 2) using chebyshev approximation with a (16, 28) grid
     ! optimized to return both (nu, nu + 1) in around the same time, therefore nu must be in (0, 1)
     ! no checks are performed on arguments
     elemental complex(BK) function bessely_chebyshev_low_orders(nu, x) result(cheb)
        real(BK), intent(in) :: nu, x

        real(BK), parameter :: TWOO13 = TWO/13.0_BK
        real(BK) :: x1,nu1,nu2,a(size(bessely_cheb_weights,2))
        integer :: i
        ! need to rescale inputs according to
        !x0 = (x - lb) * 2 / (ub - lb) - 1
        x1  = (x - SIX)*TWOO13 - ONE
        nu1 = nu - ONE
        nu2 = nu
        forall(i=1:size(bessely_cheb_weights,2)) a(i) = clenshaw_chebyshev(x1, bessely_cheb_weights(:,i))

        cheb%re = clenshaw_chebyshev(nu1, a)
        cheb%im = clenshaw_chebyshev(nu2, a)

     end function bessely_chebyshev_low_orders

     ! use the Clenshaw algorithm to recursively evaluate a linear combination of Chebyshev polynomials
     pure real(BK) function clenshaw_chebyshev(x, c) result(cheb)
        real(BK), intent(in) :: x, c(:)

        real(BK) :: x2,c0,c1,a,b
        integer  :: lc,i

        lc = size(c)
        x2 = 2*x

        c0 = c(lc-1)
        c1 = c(lc)
        do i=lc-2,1,-1
           a = c(i) - c1
           b = c0 + c1 * x2
           c0 = a
           c1 = b
        end do

        cheb = c0 + c1 * x
     end function clenshaw_chebyshev

    ! Bessel function of the second kind of order nu, ``Y_{nu}(x)``.
    ! nu and x must be real and nu and x must be positive.
    ! No checks on arguments are performed and should only be called if certain nu, x >= 0.
    elemental real(BK) function bessely_positive_args(nu, x)
       real(BK), intent(in) :: nu, x

       real(BK) :: dummy,YJ(2)

       if (x==ZERO) then

          bessely_positive_args = ieee_value(x,ieee_negative_inf)

       elseif (isinteger(nu) .and. nu<250.0_BK) then

          ! use forward recurrence if nu is an integer up until it becomes inefficient
          call besselj_up_recurrence(x, bessely1(x), bessely0(x), ONE, nu, bessely_positive_args, dummy)

       elseif (besseljy_debye_cutoff(nu, x)) then

          ! x < ~nu branch see src/U_polynomials.jl
          call besseljy_debye(nu, x, dummy, bessely_positive_args)

       elseif (besseljy_large_argument_cutoff(nu, x)) then

          ! large argument branch see src/asymptotics.jl
          call besseljy_large_argument(nu, x, dummy, bessely_positive_args)

       elseif (hankel_debye_cutoff(nu, x)) then

          ! x > ~nu branch see src/U_polynomials.jl on computing Hankel function
          bessely_positive_args = aimag(hankel_debye(nu, x))

       elseif (bessely_series_cutoff(nu, x)) then

          ! use power series for small x and for when nu > x
          YJ = bessely_power_series(nu, x)
          bessely_positive_args = YJ(1)

       else

          ! shift nu down and use upward recurrence starting from either chebyshev approx or hankel expansion
          bessely_positive_args = bessely_fallback(nu,x)

       endif

    end function bessely_positive_args


    ! Bessel functions of the second kind of order zero
    ! Calculation of bessely0 is done in three branches using polynomial approximations
    !
    ! Branch 1: x <= 5.0
    !           bessely0 = R(x^2) + 2*log(x)*besselj0(x) / pi
    !           where r1 and r2 are zeros of J0 and P3 and Q8 are a 3 and 8 degree polynomial respectively
    !           Polynomial coefficients are from [1] which is based on [2]. For tiny arugments the power
    !           series  expansion is used.
    !
    ! Branch 2: 5.0 < x < 25.0
    !           bessely0 = sqrt(2/(pi*x))*(sin(x - pi/4)*R7(x) - cos(x - pi/4)*R8(x))
    !           Hankel's asymptotic expansion is used where R7 and R8 are rational functions (Pn(x)/Qn(x))
    !           of degree 7 and 8 respectively See section 4 of [3] for more details and [1] for coefficients
    !           of polynomials
    !
    ! Branch 3: x >= 25.0
    !           bessely0 = sqrt(2/(pi*x))*beta(x)*(sin(x - pi/4 - alpha(x))
    !           See modified expansions given in [3]. Exact coefficients are used.
    !
    ! [1] https://github.com/deepmind/torch-cephes
    ! [2] Cephes Math Library Release 2.8:  June, 2000 by Stephen L. Moshier
    ! [3] Harrison, John. "Fast and accurate Bessel function computation.", 2009 19th IEEE Symposium on
    !     Computer Arithmetic. IEEE, 2009.
    !
    elemental real(BK) function bessely0(x)
       real(BK), intent(in) :: x

       real(BK) :: z,w,p,q,xn,xinv,x2,a,b
       intrinsic :: sqrt

       real(BK), parameter :: ppoly(8) = [ONE, -1.0_BK/16.0_BK, 53.0_BK/512.0_BK, -4447.0_BK/8192.0_BK, &
                                          3066403.0_BK/524288.0_BK, -896631415.0_BK/8388608.0_BK, &
                                          796754802993.0_BK/268435456.0_BK, -500528959023471.0_BK/4294967296.0_BK]

       real(BK), parameter :: qpoly(7) = [-1.0_BK/8.0_BK, 25.0_BK/384.0_BK, -1073.0_BK/5120.0_BK, &
                                          375733.0_BK/229376.0_BK, -55384775.0_BK/2359296.0_BK, &
                                          24713030909.0_BK/46137344.0_BK, -7780757249041.0_BK/436207616.0_BK]

       ! Handle
       if (x<ZERO) then

          bessely0 = ieee_value(bessely0,ieee_quiet_nan)

       elseif (x==ZERO) then

          bessely0 = -huge(bessely0)

       elseif (x<=FIVE) then

          z = x**2

          ! TODO: replace the two polynomials with a single one
          w = evalpoly(size(YP_Y0),z,YP_Y0) / evalpoly(size(YQ_Y0),z,YQ_Y0)
          bessely0 = w + TWOOPI * log(x) * besselj0(x)

       elseif (x<25.0_BK) then

          w = FIVE / x
          z = w**2

          ! TODO: replace the two polynomials with a single one
          p = evalpoly(size(PP_Y0), z, PP_Y0) / evalpoly(size(PQ_Y0), z, PQ_Y0)
          q = evalpoly(size(QP_Y0), z, QP_y0) / evalpoly(size(QQ_Y0), z, QQ_Y0)
          xn = x - PIO4
          p = p * sin(xn) + w * q * cos(xn)
          bessely0 = p * SQ2OPI / sqrt(x)

       elseif (x<huge(bessely0)) then

          xinv = ONE/x
          x2 = xinv**2
          if (x < 120.0_BK) then
              p = evalpoly(8, x2, ppoly)
              q = evalpoly(7, x2, qpoly)
          else
              p = evalpoly(4, x2, ppoly)
              q = evalpoly(4, x2, qpoly)
          endif

          a  = SQ2OPI * sqrt(xinv) * p
          xn = muladd(xinv, q, -PIO4)
          b = sin(x+xn)
          bessely0 = a*b

       else

          bessely0 = ZERO

       endif

    end function bessely0


    ! Bessel functions of the second kind of order one, ``Y_1(x)``
    ! Calculation of bessely1 is done similar to bessely0
    !
    elemental real(BK) function bessely1(x)
       real(BK), intent(in) :: x

       real(BK) :: z,w,p,q,xn,xinv,x2,a,b
       intrinsic :: sqrt

       real(BK), parameter :: ppoly(8) = [ONE, 3.0_BK/16.0_BK, -99.0_BK/512.0_BK, 6597.0_BK/8192.0_BK, &
                                          -4057965.0_BK/524288.0_BK, 1113686901.0_BK/8388608.0_BK, &
                                          -951148335159.0_BK/268435456.0_BK, 581513783771781.0_BK/4294967296.0_BK]

       real(BK), parameter :: qpoly(7) = [3.0_BK/8.0_BK, -21.0_BK/128.0_BK, 1899.0_BK/5120.0_BK, &
                                          -543483.0_BK/229376.0_BK, 8027901.0_BK/262144.0_BK, &
                                          -30413055339.0_BK/46137344.0_BK, 9228545313147.0_BK/436207616.0_BK]

       ! Handle
       if (x<ZERO) then

          bessely1 = ieee_value(bessely1,ieee_quiet_nan)

       elseif (x==ZERO) then

          bessely1 = -huge(bessely1)

       elseif (x<=FIVE) then

          z = x**2

          ! TODO: replace the two polynomials with a single one
          w = x*evalpoly(size(YP_Y1),z,YP_Y1) / evalpoly(size(YQ_Y1),z,YQ_Y1)
          bessely1 = w + TWOOPI * (besselj1(x)*log(x) - ONE/x)


       elseif (x<25.0_BK) then

          w = FIVE / x
          z = w**2

          ! TODO: replace the two polynomials with a single one
          p = evalpoly(size(PP_J1), z, PP_J1) / evalpoly(size(PQ_j1), z, PQ_J1)
          q = evalpoly(size(QP_J1), z, QP_J1) / evalpoly(size(QQ_j1), z, QQ_j1)
          xn = x - THPIO4
          p = p * sin(xn) + w * q * cos(xn)
          bessely1 = p * SQ2OPI / sqrt(x)

       elseif (x<huge(bessely1)) then

          xinv = ONE/x
          x2 = xinv**2
          if (x < 130.0_BK) then
              p = evalpoly(8, x2, ppoly)
              q = evalpoly(7, x2, qpoly)
          else
              p = evalpoly(4, x2, ppoly)
              q = evalpoly(4, x2, qpoly)
          endif

          a  = SQ2OPI * sqrt(xinv) * p
          xn = muladd(xinv, q, -3*PIO4)
          b = sin(x+xn)
          bessely1 = a*b

       else

          bessely1 = ZERO

       endif

    end function bessely1

    ! Modified Bessel functions of the second kind of order zero besselk0
    !
    ! Calculation of besselk0 is done in two branches using polynomial approximations [1]
    !
    ! Branch 1: x < 1.0
    !              besselk0(x) + log(x)*besseli0(x) = P7(x^2)
    !              besseli0(x) = [x/2]^2 * P6([x/2]^2) + 1
    ! Branch 2: x >= 1.0
    !              sqrt(x) * exp(x) * besselk0(x) = P22(1/x) / Q2(1/x)
    !    where P7, P6, P22, and Q2 are 7, 6, 22, and 2 degree polynomials respectively.
    !
    ! The polynomial coefficients are taken from boost math library [3].
    !
    ! [1] "Rational Approximations for the Modified Bessel Function of the Second Kind
    !     - K0(x) for Computations with Double Precision" by Pavel Holoborodko
    ! [2] "Rational Approximations for the Modified Bessel Function of the Second Kind
    !     - K1(x) for Computations with Double Precision" by Pavel Holoborodko
    ! [3] https://github.com/boostorg/math/tree/develop/include/boost/math/special_functions/detail

    elemental real(BK) function besselk0(x)
       real(BK), intent(in) :: x

       real(BK) :: a,s,x2,rx
       intrinsic :: sqrt

       if (x<=ZERO) then

           besselk0 = ieee_value(besselk0,ieee_quiet_nan)

       elseif (x<=ONE) then

           x2 = x**2
           a  = FOURTH*x2
           s  = evalpoly(size(P1_K0),a,P1_K0)/evalpoly(size(Q1_K0),a,Q1_K0) + Y_K0
           besselk0 = evalpoly(size(P2_K0), x2, P2_K0) - (ONE+s*a)*log(x)

       else

           s  = exp(-HALF*x)
           rx = ONE/x
           a  = (evalpoly(size(P3_K0),rx,P3_K0)/evalpoly(size(Q3_K0),rx,Q3_K0) + ONE) * s / sqrt(x)
           besselk0 = a*s

       end if

    end function besselk0


    ! Modified Bessel function of the second kind of order one, ``K_1(x)``.
    !
    ! Calculation of besselk1 is done in two branches using polynomial approximations [2]
    !
    ! Branch 1: x < 1.0
    !           besselk1(x) - log(x)*besseli1(x) - 1/x = x*P8(x^2)
    !                         besseli1(x) = [x/2]^2 * (1 + 0.5 * (*x/2)^2 + (x/2)^4 * P5([x/2]^2))
    ! Branch 2: x >= 1.0
    !           sqrt(x) * exp(x) * besselk1(x) = P22(1/x) / Q2(1/x)
    !           where P8, P5, P22, and Q2 are 8, 5, 22, and 2 degree polynomials respectively.
    elemental real(BK) function besselk1(x)
        real(BK), intent(in) :: x

        real(BK) :: z,a,pq,rx,s

        if (x<=ZERO) then

            besselk1 = ieee_value(besselk1,ieee_quiet_nan)

        elseif (x<=ONE) then

            z  = x**2
            rx = ONE/x
            a  = FOURTH*z
            pq = evalpoly(size(P1_k1), a, P1_k1)/evalpoly(size(Q1_k1), a, Q1_k1) + Y_K1
            pq = pq*a**2 + (HALF*a + ONE)
            a  = HALF*pq*x
            pq = x*evalpoly(size(P2_k1), z, P2_k1)/evalpoly(size(Q2_k1), z, Q2_k1) + rx
            besselk1 = a*log(x)+pq

        else

            s  = exp(-HALF*x)
            rx = ONE/x
            a  = evalpoly(size(P3_k1), rx, P3_k1)/evalpoly(size(Q3_k1), rx, Q3_k1) + Y2_K1
            besselk1 = a * s**2/sqrt(x)

        end if

    end function besselk1

    ! Modified Bessel function of the first kind of order zero, ``I_0(x)``.
    !
    !  Calculation of besseli0 is done in two branches using polynomial approximations [1]
    !
    !  Branch 1: x < 7.75
    !            besseli0 = [x/2]^2 P16([x/2]^2)
    !  Branch 2: x >= 7.75
    !            sqrt(x) * exp(-x) * besseli0(x) = P22(1/x)
    !            where P16 and P22 are a 16 and 22 degree polynomial respectively.
    !
    !  Remez.jl is then used to approximate the polynomial coefficients of
    !  P22(y) = sqrt(1/y) * exp(-inv(y)) * besseli0(inv(y))
    !  N,D,E,X = ratfn_minimax(g, (1/1e6, 1/7.75), 21, 0)
    !
    elemental real(BK) function besseli0(x)
       real(BK), intent(in) :: x
       real(BK) :: ax,a,s
        ax = abs(x)

        if (ax < 7.75_BK) then
            a = FOURTH*ax**2
            besseli0 = a*evalpoly(size(besseli0_small_coefs), a, besseli0_small_coefs) + ONE
        else
            a = exp(HALF*ax)
            s = a * evalpoly(size(besseli0_med_coefs), one/ax, besseli0_med_coefs) / sqrt(ax)
            besseli0 = a * s
        endif

    end function besseli0

    ! Modified Bessel function of the first kind of order one, ``I_1(x)``.
    !
    ! Calculation of besseli1 is done in two branches using polynomial approximations [2]
    !
    !  Branch 1: x < 7.75
    !            besseli1 = x / 2 * (1 + 1/2 * (x/2)^2 + (x/2)^4 * P13([x/2]^2)
    !  Branch 2: x >= 7.75
    !            sqrt(x) * exp(-x) * besseli1(x) = P22(1/x)
    !            where P13 and P22 are a 16 and 22 degree polynomial respectively.
    !
    !  Remez.jl is then used to approximate the polynomial coefficients of
    !  P13(y) = (besseli1(2 * sqrt(y)) / sqrt(y) - 1 - y/2) / y^2
    !  N,D,E,X = ratfn_minimax(g, (1/1e6, 1/7.75), 21, 0)
    !
    !  A third branch is used for scaled functions for large values
    !
    elemental real(BK) function besseli1(x)
       real(BK), intent(in) :: x

       real(BK) :: z,inner(3),s,a

       z = abs(x)
       if (z < 7.75_BK) then

           a = FOURTH*z**2
           inner = [ONE,HALF,evalpoly(size(besseli1_small_coefs),a,besseli1_small_coefs)]
           z = HALF*z*evalpoly(3,a,inner)

       else

           a = exp(HALF*z)
           s = a*evalpoly(size(besseli1_med_coefs), one/z, besseli1_med_coefs)/sqrt(z)
           z = a*s

       end if

       besseli1 = sign(ONE,x)*z

   end function besseli1


end module bessels
