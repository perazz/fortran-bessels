module bessels_rkbesl
!  Modernized free-form Fortran wrapper around the netlib RKBESL reference routine.
!  Test-only: used by bessels_test to compare the modern bessels library against a
!  known-good reference. Algorithm and coefficients preserved verbatim from the
!  original (W. J. Cody & L. Stoltz, Argonne National Laboratory, May 30 1989).

   use iso_fortran_env, only: real64

   implicit none
   private

   public :: rkbesl

contains

   subroutine rkbesl(x, alpha, nb, ize, bk, ncalc)
      !-------------------------------------------------------------------
      !  This routine calculates modified Bessel functions of the second
      !  kind, K SUB(N+ALPHA) (X), for non-negative argument X, and
      !  non-negative order N+ALPHA, with or without exponential scaling.
      !
      !  Explanation of variables in the calling sequence:
      !
      !  X     - Working precision non-negative real argument for which
      !          K's or exponentially scaled K's (K*EXP(X)) are to be
      !          calculated. If K's are to be calculated, X must not be
      !          greater than XMAX.
      !  ALPHA - Working precision fractional part of order for which
      !          K's or exponentially scaled K's (K*EXP(X)) are to be
      !          calculated.  0 .LE. ALPHA .LT. 1.0.
      !  NB    - Integer number of functions to be calculated, NB .GT. 0.
      !  IZE   - Integer type.  IZE = 1 for unscaled K's, 2 for
      !          exponentially scaled K's.
      !  BK    - Working precision output vector of length NB.
      !  NCALC - Integer output variable indicating possible errors.
      !          Before using BK, the user should check that NCALC = NB.
      !
      !  Machine-dependent constants (IEEE double precision values used
      !  below): EPS=2.22D-16, SQXMIN=1.49D-154, XINF=1.79D+308,
      !  XMIN=2.23D-308, XMAX=705.342.
      !
      !  Acknowledgement
      !
      !  This program is based on a program written by J. B. Campbell
      !  that computes values of the Bessel functions K of real argument
      !  and real order. Modifications include the addition of non-scaled
      !  functions, parameterization of machine dependencies, and the use
      !  of more accurate approximations for SINH and SIN.
      !
      !  References:
      !    "On Temme's Algorithm for the Modified Bessel Functions of the
      !     Third Kind," Campbell, J. B., TOMS 6(4), Dec. 1980,
      !     pp. 581-586.
      !    "A FORTRAN IV Subroutine for the Modified Bessel Functions of
      !     the Third Kind of Real Order and Real Argument," Campbell,
      !     J. B., Report NRC/ERB-925, National Research Council, Canada.
      !-------------------------------------------------------------------
      real(real64), intent(in)  :: x
      real(real64), intent(in)  :: alpha
      integer,      intent(in)  :: nb
      integer,      intent(in)  :: ize
      real(real64), intent(out) :: bk(nb)
      integer,      intent(out) :: ncalc

      integer      :: i, iend, itemp, j, k, m, mplus1
      real(real64) :: blpha, bk1, bk2, c, dm, d1, d2, d3, enu, ex
      real(real64) :: f0, f1, f2, p0, q0, ratio, twonu, twox, t1, t2
      real(real64) :: wminf, x2by4

      !-------------------------------------------------------------------
      !  Mathematical constants
      !    A = log(2) - Euler's constant
      !    D = sqrt(2/pi)
      !-------------------------------------------------------------------
      real(real64), parameter :: ZERO  = 0.0_real64
      real(real64), parameter :: HALF  = 0.5_real64
      real(real64), parameter :: ONE   = 1.0_real64
      real(real64), parameter :: TWO   = 2.0_real64
      real(real64), parameter :: FOUR  = 4.0_real64
      real(real64), parameter :: TINYX = 1.0e-10_real64
      real(real64), parameter :: A     = 0.11593151565841244881_real64
      real(real64), parameter :: D     = 0.797884560802865364_real64

      !-------------------------------------------------------------------
      !  Machine-dependent parameters (IEEE double precision)
      !-------------------------------------------------------------------
      real(real64), parameter :: EPS    = 2.22e-16_real64
      real(real64), parameter :: SQXMIN = 1.49e-154_real64
      real(real64), parameter :: XINF   = 1.79e+308_real64
      real(real64), parameter :: XMIN   = 2.23e-308_real64
      real(real64), parameter :: XMAX   = 705.342_real64

      !-------------------------------------------------------------------
      !  P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA + Euler's constant
      !         Coefficients converted from hex to decimal and modified
      !         by W. J. Cody, 2/26/82
      !  R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2*ALPHA)
      !  T    - Approximation for SINH(Y)/Y
      !-------------------------------------------------------------------
      real(real64), parameter :: P(8) = [ &
         0.805629875690432845e+00_real64,  0.204045500205365151e+02_real64, &
         0.157705605106676174e+03_real64,  0.536671116469207504e+03_real64, &
         0.900382759291288778e+03_real64,  0.730923886650660393e+03_real64, &
         0.229299301509425145e+03_real64,  0.822467033424113231e+00_real64 ]
      real(real64), parameter :: Q(7) = [ &
         0.294601986247850434e+02_real64,  0.277577868510221208e+03_real64, &
         0.120670325591027438e+04_real64,  0.276291444159791519e+04_real64, &
         0.344374050506564618e+04_real64,  0.221063190113378647e+04_real64, &
         0.572267338359892221e+03_real64 ]
      real(real64), parameter :: R(5) = [ &
        -0.48672575865218401848e+00_real64, 0.13079485869097804016e+02_real64, &
        -0.10196490580880537526e+03_real64, 0.34765409106507813131e+03_real64, &
         0.34958981245219347820e-03_real64 ]
      real(real64), parameter :: S(4) = [ &
        -0.25579105509976461286e+02_real64, 0.21257260432226544008e+03_real64, &
        -0.61069018684944109624e+03_real64, 0.42269668805777760407e+03_real64 ]
      real(real64), parameter :: T(6) = [ &
         0.16125990452916363814e-09_real64, 0.25051878502858255354e-07_real64, &
         0.27557319615147964774e-05_real64, 0.19841269840928373686e-03_real64, &
         0.83333333333334751799e-02_real64, 0.16666666666666666446e+00_real64 ]
      real(real64), parameter :: ESTM(6) = [ &
         5.20583e+01_real64, 5.7607e+00_real64, 2.7782e+00_real64, &
         1.44303e+01_real64, 1.853004e+02_real64, 9.3715e+00_real64 ]
      real(real64), parameter :: ESTF(7) = [ &
         4.18341e+01_real64, 7.1075e+00_real64, 6.4306e+00_real64, &
         4.25110e+01_real64, 1.35633e+00_real64, 8.45096e+01_real64, &
         2.0e+01_real64 ]

      !-------------------------------------------------------------------
      ex = x
      enu = alpha
      ncalc = min(nb, 0) - 2
      if ((nb > 0) .and. ((enu >= ZERO) .and. (enu < ONE)) .and. &
          ((ize >= 1) .and. (ize <= 2)) .and. &
          ((ize /= 1) .or. (ex <= XMAX)) .and. (ex > ZERO)) then

         k = 0
         if (enu < SQXMIN) enu = ZERO
         if (enu > HALF) then
            k = 1
            enu = enu - ONE
         end if
         twonu = enu + enu
         iend = nb + k - 1
         c = enu*enu
         d3 = -c
         if (ex <= ONE) then
            !---------------------------------------------------------------------
            !  Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA
            !                 Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA
            !---------------------------------------------------------------------
            d1 = ZERO
            d2 = P(1)
            t1 = ONE
            t2 = Q(1)
            do i = 2, 7, 2
               d1 = c*d1 + P(i)
               d2 = c*d2 + P(i+1)
               t1 = c*t1 + Q(i)
               t2 = c*t2 + Q(i+1)
            end do
            d1 = enu*d1
            t1 = enu*t1
            f1 = log(ex)
            f0 = A + enu*(P(8) - enu*(d1+d2)/(t1+t2)) - f1
            q0 = exp(-enu*(A - enu*(P(8) + enu*(d1-d2)/(t1-t2)) - f1))
            f1 = enu*f0
            p0 = exp(f1)
            !---------------------------------------------------------------------
            !  Calculation of F0
            !---------------------------------------------------------------------
            d1 = R(5)
            t1 = ONE
            do i = 1, 4
               d1 = c*d1 + R(i)
               t1 = c*t1 + S(i)
            end do
            if (abs(f1) <= HALF) then
               f1 = f1*f1
               d2 = ZERO
               do i = 1, 6
                  d2 = f1*d2 + T(i)
               end do
               d2 = f0 + f0*f1*d2
            else
               d2 = sinh(f1) / enu
            end if
            f0 = d2 - enu*d1/(t1*p0)
            if (ex <= TINYX) then
               !--------------------------------------------------------------------
               !  X .LE. 1.0E-10
               !  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X)
               !--------------------------------------------------------------------
               bk(1) = f0 + ex*f0
               if (ize == 1) bk(1) = bk(1) - ex*bk(1)
               ratio = p0 / f0
               c = ex*XINF
               if (k /= 0) then
                  !--------------------------------------------------------------------
                  !  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X),
                  !  ALPHA .GE. 1/2
                  !--------------------------------------------------------------------
                  ncalc = -1
                  if (bk(1) >= c/ratio) goto 500
                  bk(1) = ratio*bk(1)/ex
                  twonu = twonu + TWO
                  ratio = twonu
               end if
               ncalc = 1
               if (nb == 1) goto 500
               !--------------------------------------------------------------------
               !  Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X), L = 1, 2, ... , NB-1
               !--------------------------------------------------------------------
               ncalc = -1
               do i = 2, nb
                  if (ratio >= c) goto 500
                  bk(i) = ratio/ex
                  twonu = twonu + TWO
                  ratio = twonu
               end do
               ncalc = 1
               goto 420
            else
               !--------------------------------------------------------------------
               !  1.0E-10 .LT. X .LE. 1.0
               !--------------------------------------------------------------------
               c = ONE
               x2by4 = ex*ex / FOUR
               p0 = HALF*p0
               q0 = HALF*q0
               d1 = -ONE
               d2 = ZERO
               bk1 = ZERO
               bk2 = ZERO
               f1 = f0
               f2 = p0
               do
                  d1 = d1 + TWO
                  d2 = d2 + ONE
                  d3 = d1 + d3
                  c = x2by4*c / d2
                  f0 = (d2*f0 + p0 + q0) / d3
                  p0 = p0 / (d2 - enu)
                  q0 = q0 / (d2 + enu)
                  t1 = c*f0
                  t2 = c*(p0 - d2*f0)
                  bk1 = bk1 + t1
                  bk2 = bk2 + t2
                  if ((abs(t1/(f1+bk1)) <= EPS) .and. &
                      (abs(t2/(f2+bk2)) <= EPS)) exit
               end do
               bk1 = f1 + bk1
               bk2 = TWO*(f2 + bk2) / ex
               if (ize == 2) then
                  d1 = exp(ex)
                  bk1 = bk1*d1
                  bk2 = bk2*d1
               end if
               wminf = ESTF(1)*ex + ESTF(2)
            end if
         else if (EPS*ex > ONE) then
            !--------------------------------------------------------------------
            !  X .GT. ONE/EPS
            !--------------------------------------------------------------------
            ncalc = nb
            bk1 = ONE / (D*sqrt(ex))
            do i = 1, nb
               bk(i) = bk1
            end do
            goto 500
         else
            !--------------------------------------------------------------------
            !  X .GT. 1.0
            !--------------------------------------------------------------------
            twox = ex + ex
            blpha = ZERO
            ratio = ZERO
            if (ex <= FOUR) then
               !--------------------------------------------------------------------
               !  Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 .LE. X .LE. 4.0
               !--------------------------------------------------------------------
               d2 = aint(ESTM(1)/ex + ESTM(2))
               m = int(d2)
               d1 = d2 + d2
               d2 = d2 - HALF
               d2 = d2*d2
               do i = 2, m
                  d1 = d1 - TWO
                  d2 = d2 - d1
                  ratio = (d3 + d2) / (twox + d1 - ratio)
               end do
               !--------------------------------------------------------------------
               !  Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward
               !    recurrence and K(ALPHA,X) from the wronskian
               !--------------------------------------------------------------------
               d2 = aint(ESTM(3)*ex + ESTM(4))
               m = int(d2)
               c = abs(enu)
               d3 = c + c
               d1 = d3 - ONE
               f1 = XMIN
               f0 = (TWO*(c + d2)/ex + HALF*ex/(c + d2 + ONE))*XMIN
               do i = 3, m
                  d2 = d2 - ONE
                  f2 = (d3 + d2 + d2)*f0
                  blpha = (ONE + d1/d2)*(f2 + blpha)
                  f2 = f2/ex + f1
                  f1 = f0
                  f0 = f2
               end do
               f1 = (d3 + TWO)*f0/ex + f1
               d1 = ZERO
               t1 = ONE
               do i = 1, 7
                  d1 = c*d1 + P(i)
                  t1 = c*t1 + Q(i)
               end do
               p0 = exp(c*(A + c*(P(8) - c*d1/t1) - log(ex))) / ex
               f2 = (c + HALF - ratio)*f1 / ex
               bk1 = p0 + (d3*f0 - f2 + f0 + blpha) / (f2 + f1 + f0) * p0
               if (ize == 1) bk1 = bk1*exp(-ex)
               wminf = ESTF(3)*ex + ESTF(4)
            else
               !--------------------------------------------------------------------
               !  Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by backward
               !  recurrence, for  X .GT. 4.0
               !--------------------------------------------------------------------
               dm = aint(ESTM(5)/ex + ESTM(6))
               m = int(dm)
               d2 = dm - HALF
               d2 = d2*d2
               d1 = dm + dm
               do i = 2, m
                  dm = dm - ONE
                  d1 = d1 - TWO
                  d2 = d2 - d1
                  ratio = (d3 + d2) / (twox + d1 - ratio)
                  blpha = (ratio + ratio*blpha) / dm
               end do
               bk1 = ONE / ((D + D*blpha)*sqrt(ex))
               if (ize == 1) bk1 = bk1*exp(-ex)
               wminf = ESTF(5)*(ex - abs(ex - ESTF(7))) + ESTF(6)
            end if
            !--------------------------------------------------------------------
            !  Calculation of K(ALPHA+1,X) from K(ALPHA,X) and
            !    K(ALPHA+1,X)/K(ALPHA,X)
            !--------------------------------------------------------------------
            bk2 = bk1 + bk1*(enu + HALF - ratio) / ex
         end if
         !--------------------------------------------------------------------
         !  Calculation of 'NCALC', K(ALPHA+I,X), I = 0, 1, ... , NCALC-1,
         !  K(ALPHA+I,X)/K(ALPHA+I-1,X), I = NCALC, NCALC+1, ... , NB-1
         !--------------------------------------------------------------------
         ncalc = nb
         bk(1) = bk1
         if (iend == 0) goto 500
         j = 2 - k
         if (j > 0) bk(j) = bk2
         if (iend == 1) goto 500
         m = min(int(wminf - enu), iend)
         itemp = 1
         do i = 2, m
            t1 = bk1
            bk1 = bk2
            twonu = twonu + TWO
            if (ex < ONE) then
               if (bk1 >= (XINF/twonu)*ex) goto 195
               goto 187
            else
               if (bk1/ex >= XINF/twonu) goto 195
            end if
187         bk2 = twonu/ex*bk1 + t1
            itemp = i
            j = j + 1
            if (j > 0) bk(j) = bk2
         end do
195      m = itemp
         if (m == iend) goto 500
         ratio = bk2 / bk1
         mplus1 = m + 1
         ncalc = -1
         do i = mplus1, iend
            twonu = twonu + TWO
            ratio = twonu/ex + ONE/ratio
            j = j + 1
            if (j > 1) then
               bk(j) = ratio
            else
               if (bk2 >= XINF/ratio) goto 500
               bk2 = ratio*bk2
            end if
         end do
         ncalc = max(mplus1 - k, 1)
         if (ncalc == 1) bk(1) = bk2
         if (nb == 1) goto 500
420      j = ncalc + 1
         do i = j, nb
            if (bk(ncalc) >= XINF/bk(i)) goto 500
            bk(i) = bk(ncalc)*bk(i)
            ncalc = i
         end do
      end if
500   return

   end subroutine rkbesl

end module bessels_rkbesl
