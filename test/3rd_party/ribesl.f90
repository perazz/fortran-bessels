module bessels_ribesl
!  Modernized free-form Fortran wrapper around the netlib RIBESL reference routine.
!  Test-only: used by bessels_test to compare the modern bessels library against a
!  known-good reference. Algorithm and coefficients preserved verbatim from the
!  original (W. J. Cody & L. Stoltz, Argonne National Laboratory, May 30 1989).

   use iso_fortran_env, only: real64

   implicit none
   private

   public :: ribesl

contains

   subroutine ribesl(x, alpha, nb, ize, b, ncalc)
      !-------------------------------------------------------------------
      !  This routine calculates Bessel functions I SUB(N+ALPHA) (X)
      !  for non-negative argument X, and non-negative order N+ALPHA,
      !  with or without exponential scaling.
      !
      !  Explanation of variables in the calling sequence:
      !
      !  X     - Working precision non-negative real argument for which
      !          I's or exponentially scaled I's (I*EXP(-X)) are to be
      !          calculated. If I's are to be calculated, X must be less
      !          than EXPARG.
      !  ALPHA - Working precision fractional part of order for which
      !          I's or exponentially scaled I's (I*EXP(-X)) are to be
      !          calculated.  0 .LE. ALPHA .LT. 1.0.
      !  NB    - Integer number of functions to be calculated, NB .GT. 0.
      !          The first function is of order ALPHA, the last is of
      !          order (NB - 1 + ALPHA).
      !  IZE   - Integer type.  IZE = 1 for unscaled I's, 2 for
      !          exponentially scaled I's.
      !  B     - Working precision output vector of length NB.
      !  NCALC - Integer output variable indicating possible errors.
      !          Before using B, the user should check that NCALC = NB.
      !
      !  Machine-dependent constants (IEEE double precision values used
      !  below): NSIG=16, ENTEN=1.0D+308, ENSIG=1.0D+16, RTNSIG=1.0D-4,
      !  ENMTEN=8.9D-308, XLARGE=1.0D+4, EXPARG=709.
      !
      !  Acknowledgement
      !
      !  This program is based on a program written by David J. Sookne
      !  that computes values of the Bessel functions J or I of real
      !  argument and integer order.
      !
      !  References:
      !    "A Note on Backward Recurrence Algorithms," Olver, F. W. J.,
      !     and Sookne, D. J., Math. Comp. 26, 1972, pp 941-947.
      !    "Bessel Functions of Real Argument and Integer Order," Sookne,
      !     D. J., NBS Jour. of Res. B. 77B, 1973, pp 125-132.
      !    "ALGORITHM 597, Sequence of Modified Bessel Functions of the
      !     First Kind," Cody, W. J., Trans. Math. Soft., 1983,
      !     pp. 242-245.
      !-------------------------------------------------------------------
      real(real64), intent(in)  :: x
      real(real64), intent(in)  :: alpha
      integer,      intent(in)  :: nb
      integer,      intent(in)  :: ize
      real(real64), intent(out) :: b(nb)
      integer,      intent(out) :: ncalc

      integer      :: k, l, magx, n, nbmx, nend, nstart
      real(real64) :: em, empal, emp2al, en, halfx
      real(real64) :: p, plast, pold, psave, psavel
      real(real64) :: sum, tempa, tempb, tempc, test, tover

      !-------------------------------------------------------------------
      !  Mathematical constants
      !-------------------------------------------------------------------
      real(real64), parameter :: ZERO  = 0.0_real64
      real(real64), parameter :: HALF  = 0.5_real64
      real(real64), parameter :: ONE   = 1.0_real64
      real(real64), parameter :: TWO   = 2.0_real64
      real(real64), parameter :: CONST = 1.585_real64

      !-------------------------------------------------------------------
      !  Machine-dependent parameters (IEEE double precision)
      !-------------------------------------------------------------------
      integer,      parameter :: NSIG   = 16
      real(real64), parameter :: XLARGE = 1.0e4_real64
      real(real64), parameter :: EXPARG = 709.0_real64
      real(real64), parameter :: ENTEN  = 1.0e308_real64
      real(real64), parameter :: ENSIG  = 1.0e16_real64
      real(real64), parameter :: RTNSIG = 1.0e-4_real64
      real(real64), parameter :: ENMTEN = 8.9e-308_real64

      !-------------------------------------------------------------------
      !  Check for X, NB, or IZE out of range.
      !-------------------------------------------------------------------
      if ((nb > 0) .and. (x >= ZERO) .and. &
          (alpha >= ZERO) .and. (alpha < ONE) .and. &
          (((ize == 1) .and. (x <= EXPARG)) .or. &
           ((ize == 2) .and. (x <= XLARGE)))) then

         !-------------------------------------------------------------------
         !  Use 2-term ascending series for small X
         !-------------------------------------------------------------------
         ncalc = nb
         magx = int(x)
         if (x >= RTNSIG) then
            !-------------------------------------------------------------------
            !  Initialize the forward sweep, the P-sequence of Olver
            !-------------------------------------------------------------------
            nbmx = nb - magx
            n = magx + 1
            en = real(n + n, real64) + (alpha + alpha)
            plast = ONE
            p = en / x
            !-------------------------------------------------------------------
            !  Calculate general significance test
            !-------------------------------------------------------------------
            test = ENSIG + ENSIG
            if (2*magx > 5*NSIG) then
               test = sqrt(test*p)
            else
               test = test / CONST**magx
            end if
            if (nbmx >= 3) then
               !-------------------------------------------------------------------
               !  Calculate P-sequence until N = NB-1.  Check for possible overflow.
               !-------------------------------------------------------------------
               tover = ENTEN / ENSIG
               nstart = magx + 2
               nend = nb - 1
               do k = nstart, nend
                  n = k
                  en = en + TWO
                  pold = plast
                  plast = p
                  p = en * plast / x + pold
                  if (p > tover) then
                     !-------------------------------------------------------------------
                     !  To avoid overflow, divide P-sequence by TOVER.  Calculate
                     !  P-sequence until ABS(P) .GT. 1.
                     !-------------------------------------------------------------------
                     tover = ENTEN
                     p = p / tover
                     plast = plast / tover
                     psave = p
                     psavel = plast
                     nstart = n + 1
                     do
                        n = n + 1
                        en = en + TWO
                        pold = plast
                        plast = p
                        p = en * plast / x + pold
                        if (p > ONE) exit
                     end do
                     tempb = en / x
                     !-------------------------------------------------------------------
                     !  Calculate backward test, and find NCALC, the highest N
                     !  such that the test is passed.
                     !-------------------------------------------------------------------
                     test = pold * plast / ENSIG
                     test = test * (HALF - HALF / (tempb*tempb))
                     p = plast * tover
                     n = n - 1
                     en = en - TWO
                     nend = min(nb, n)
                     ncalc = nend + 1
                     do l = nstart, nend
                        pold = psavel
                        psavel = psave
                        psave = en * psavel / x + pold
                        if (psave*psavel > test) then
                           ncalc = l
                           exit
                        end if
                     end do
                     ncalc = ncalc - 1
                     goto 120
                  end if
               end do
               n = nend
               en = real(n + n, real64) + (alpha + alpha)
               !-------------------------------------------------------------------
               !  Calculate special significance test for NBMX .GT. 2.
               !-------------------------------------------------------------------
               test = max(test, sqrt(plast*ENSIG) * sqrt(p + p))
            end if
            !-------------------------------------------------------------------
            !  Calculate P-sequence until significance test passed.
            !-------------------------------------------------------------------
            do
               n = n + 1
               en = en + TWO
               pold = plast
               plast = p
               p = en * plast / x + pold
               if (p >= test) exit
            end do
            !-------------------------------------------------------------------
            !  Initialize the backward recursion and the normalization sum.
            !-------------------------------------------------------------------
120         n = n + 1
            en = en + TWO
            tempb = ZERO
            tempa = ONE / p
            em = real(n, real64) - ONE
            empal = em + alpha
            emp2al = (em - ONE) + (alpha + alpha)
            sum = tempa * empal * emp2al / em
            nend = n - nb
            if (nend < 0) then
               !-------------------------------------------------------------------
               !  N .LT. NB, so store B(N) and set higher orders to zero.
               !-------------------------------------------------------------------
               b(n) = tempa
               nend = -nend
               do l = 1, nend
                  b(n + l) = ZERO
               end do
            else
               if (nend > 0) then
                  !-------------------------------------------------------------------
                  !  Recur backward via difference equation, calculating (but
                  !  not storing) B(N), until N = NB.
                  !-------------------------------------------------------------------
                  do l = 1, nend
                     n = n - 1
                     en = en - TWO
                     tempc = tempb
                     tempb = tempa
                     tempa = (en*tempb) / x + tempc
                     em = em - ONE
                     emp2al = emp2al - ONE
                     if (n == 1) goto 150
                     if (n == 2) emp2al = ONE
                     empal = empal - ONE
                     sum = (sum + tempa*empal) * emp2al / em
                  end do
               end if
               !-------------------------------------------------------------------
               !  Store B(NB)
               !-------------------------------------------------------------------
150            b(n) = tempa
               if (nb <= 1) then
                  sum = (sum + sum) + tempa
                  goto 230
               end if
               !-------------------------------------------------------------------
               !  Calculate and Store B(NB-1)
               !-------------------------------------------------------------------
               n = n - 1
               en = en - TWO
               b(n) = (en*tempa) / x + tempb
               if (n == 1) goto 220
               em = em - ONE
               emp2al = emp2al - ONE
               if (n == 2) emp2al = ONE
               empal = empal - ONE
               sum = (sum + b(n)*empal) * emp2al / em
            end if
            nend = n - 2
            if (nend > 0) then
               !-------------------------------------------------------------------
               !  Calculate via difference equation and store B(N), until N = 2.
               !-------------------------------------------------------------------
               do l = 1, nend
                  n = n - 1
                  en = en - TWO
                  b(n) = (en*b(n+1)) / x + b(n+2)
                  em = em - ONE
                  emp2al = emp2al - ONE
                  if (n == 2) emp2al = ONE
                  empal = empal - ONE
                  sum = (sum + b(n)*empal) * emp2al / em
               end do
            end if
            !-------------------------------------------------------------------
            !  Calculate B(1)
            !-------------------------------------------------------------------
            b(1) = TWO*empal*b(2) / x + b(3)
220         sum = (sum + sum) + b(1)
            !-------------------------------------------------------------------
            !  Normalize.  Divide all B(N) by sum.
            !-------------------------------------------------------------------
230         if (alpha /= ZERO) sum = sum * gamma(ONE + alpha) * (x*HALF)**(-alpha)
            if (ize == 1) sum = sum * exp(-x)
            tempa = ENMTEN
            if (sum > ONE) tempa = tempa * sum
            do n = 1, nb
               if (b(n) < tempa) b(n) = ZERO
               b(n) = b(n) / sum
            end do
            return
         else
            !-------------------------------------------------------------------
            !  Two-term ascending series for small X.
            !-------------------------------------------------------------------
            tempa = ONE
            empal = ONE + alpha
            halfx = ZERO
            if (x > ENMTEN) halfx = HALF * x
            if (alpha /= ZERO) tempa = halfx**alpha / gamma(empal)
            if (ize == 2) tempa = tempa * exp(-x)
            tempb = ZERO
            if ((x + ONE) > ONE) tempb = halfx * halfx
            b(1) = tempa + tempa*tempb / empal
            if ((x /= ZERO) .and. (b(1) == ZERO)) ncalc = 0
            if (nb > 1) then
               if (x == ZERO) then
                  do n = 2, nb
                     b(n) = ZERO
                  end do
               else
                  !-------------------------------------------------------------------
                  !  Calculate higher-order functions.
                  !-------------------------------------------------------------------
                  tempc = halfx
                  tover = (ENMTEN + ENMTEN) / x
                  if (tempb /= ZERO) tover = ENMTEN / tempb
                  do n = 2, nb
                     tempa = tempa / empal
                     empal = empal + ONE
                     tempa = tempa * tempc
                     if (tempa <= tover*empal) tempa = ZERO
                     b(n) = tempa + tempa*tempb / empal
                     if ((b(n) == ZERO) .and. (ncalc > n)) ncalc = n - 1
                  end do
               end if
            end if
         end if
      else
         ncalc = min(nb, 0) - 1
      end if

   end subroutine ribesl

end module bessels_ribesl
