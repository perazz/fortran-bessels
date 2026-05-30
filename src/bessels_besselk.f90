!  ************************************************************************************************************
!
!                                    ____  ________________ ________   _____
!                                   / __ )/ ____/ ___/ ___// ____/ /  / ___/
!                                  / __  / __/  \__ \\__ \/ __/ / /   \__ \
!                                 / /_/ / /___ ___/ /__/ / /___/ /______/ /
!                                /_____/_____//____/____/_____/_____/____/
!
!                              Modified Bessel function K_nu(x), variable order
!
!  MIT License
!  Copyright (c) 2022-2026 Federico Perini
!  Copyright (c) 2021-2022 Michael Helton, Oscar Smith, and the Bessels.jl contributors
!
!  Port of Bessels.jl/src/BesselFunctions/besselk.jl.  Three branches:
!    - Large argument asymptotic series in 1/x (NIST DLMF 10.40.E2),
!    - Debye uniform asymptotic for large order (DLMF 10.41.4) via the U_k polynomials,
!    - Power series / Temme series for small x or near-integer nu, with forward recurrence
!      bringing the seed values up to the target order.
!
!  Reflection: K_{-nu}(x) = K_{nu}(x), absorbed by abs(nu) at the top of the dispatcher.
!  Domain:    x <= 0 returns ieee_quiet_nan (matches besselk0/besselk1).
!  ************************************************************************************************************
module bessels_besselk
    use bessels_constants
    use bessels_debye,  only: Uk_poly10
    use bessels_gamma,  only: gamma_BK
    implicit none
    private

    public :: besselk, besselkx
    public :: besselk_up_recurrence

    contains

    ! "near integer": fractional part within 1e-5 of an integer.
    elemental logical function is_nearint(v) result(r)
        real(BK), intent(in) :: v
        r = abs(v - real(nint(v), BK)) < 1.0e-5_BK
    end function is_nearint

    ! ----------------------------------------------------------------------------------------------
    ! Public: K_{nu}(x)
    ! ----------------------------------------------------------------------------------------------
    elemental real(BK) function besselk(nu, x) result(y)
        real(BK), intent(in) :: nu, x
        real(BK) :: anu

        if (x /= x) then
            y = x
            return
        end if
        if (x <= ZERO) then
            y = ieee_value(y, ieee_quiet_nan)
            return
        end if

        anu = abs(nu)

        if (besselk_large_args_cutoff(anu, x)) then
            y = besselk_large_args(anu, x)
        elseif (besselk_debye_cutoff(anu, x)) then
            y = besselk_debye(anu, x)
        elseif (isinteger(anu)) then
            y = besselk_integer(anu, x, scaled=.false.)
        else
            y = besselk_noninteger(anu, x, scaled=.false.)
        end if
    end function besselk

    ! Public: e^x K_{nu}(x).
    elemental real(BK) function besselkx(nu, x) result(y)
        real(BK), intent(in) :: nu, x
        real(BK) :: anu

        if (x /= x) then
            y = x
            return
        end if
        if (x <= ZERO) then
            y = ieee_value(y, ieee_quiet_nan)
            return
        end if

        anu = abs(nu)

        if (besselk_large_args_cutoff(anu, x)) then
            y = besselkx_large_args(anu, x)
        elseif (besselk_debye_cutoff(anu, x)) then
            y = besselk_debye_scaled(anu, x)
        elseif (isinteger(anu)) then
            y = besselk_integer(anu, x, scaled=.true.)
        else
            y = besselk_noninteger(anu, x, scaled=.true.)
        end if
    end function besselkx

    ! ----------------------------------------------------------------------------------------------
    ! Cutoffs
    ! ----------------------------------------------------------------------------------------------
    elemental logical function besselk_large_args_cutoff(nu, x)
        real(BK), intent(in) :: nu, x
        besselk_large_args_cutoff = x > nu*nu/36.0_BK + 18.0_BK
    end function besselk_large_args_cutoff

    elemental logical function besselk_debye_cutoff(nu, x)
        real(BK), intent(in) :: nu, x
        besselk_debye_cutoff = (nu > 25.0_BK) .or. (x > 35.0_BK)
    end function besselk_debye_cutoff

    ! ----------------------------------------------------------------------------------------------
    ! Branches
    ! ----------------------------------------------------------------------------------------------

    ! NIST 10.40.E2: K_nu(x) = sqrt(pi/(2x)) exp(-x) sum_{k=0}^inf (4nu^2 - (2k-1)^2 ... )/(k! (8x)^k).
    elemental real(BK) function besselk_large_args(nu, x) result(y)
        real(BK), intent(in) :: nu, x
        y = exp(-HALF*x) * besselkx_large_args(nu, x) * exp(-HALF*x)
    end function besselk_large_args

    elemental real(BK) function besselkx_large_args(nu, x) result(y)
        real(BK), intent(in) :: nu, x
        integer, parameter :: MaxIter = 75
        real(BK) :: t, s, invx, fnu2
        integer :: i
        invx  = ONE/(8.0_BK*x)
        fnu2  = FOUR*nu*nu
        t = ONE
        s = ONE
        do i = 1, MaxIter
            t = t * invx * ((fnu2 - real((2*i - 1)**2, BK)) / real(i, BK))
            s = s + t
            if (abs(t) < epsilon(ONE) * abs(s)) exit
        end do
        y = s * sqrt(PIO2 / x)
    end function besselkx_large_args

    ! Debye uniform asymptotic for K_nu(x) — DLMF 10.41.4.
    elemental real(BK) function besselk_debye(nu, x) result(y)
        real(BK), intent(in) :: nu, x
        real(BK) :: coef
        coef = besselk_debye_prefactor(nu, x, scaled=.false.)
        y    = coef * besselk_debye_polynomial(nu, x)
    end function besselk_debye

    elemental real(BK) function besselk_debye_scaled(nu, x) result(y)
        real(BK), intent(in) :: nu, x
        real(BK) :: coef
        coef = besselk_debye_prefactor(nu, x, scaled=.true.)
        y    = coef * besselk_debye_polynomial(nu, x)
    end function besselk_debye_scaled

    elemental real(BK) function besselk_debye_prefactor(nu, x, scaled) result(coef)
        real(BK), intent(in) :: nu, x
        logical, intent(in)  :: scaled
        real(BK) :: z, zs, n
        z  = x/nu
        zs = sqrt(ONE + z*z)
        n  = zs + log(z) - log(ONE + zs)
        if (scaled) then
            coef = SQPIO2 * sqrt(ONE/nu) * exp(-nu*n + x) / sqrt(zs)
        else
            coef = SQPIO2 * sqrt(ONE/nu) * exp(-nu*n)     / sqrt(zs)
        end if
    end function besselk_debye_prefactor

    elemental real(BK) function besselk_debye_polynomial(nu, x) result(uk)
        real(BK), intent(in) :: nu, x
        real(BK) :: p, p2, mx, mn, Uk_In, Uk_Kn
        p  = ONE/sqrt(ONE + (x/nu)**2)
        mx = max(nu, x)
        mn = min(nu, x)
        p2 = nu*nu / (mx*mx + mn*mn)
        call Uk_poly10(p, nu, p2, Uk_In, Uk_Kn)
        uk = Uk_Kn
    end function besselk_debye_polynomial

    ! ----------------------------------------------------------------------------------------------
    ! Integer-nu fast path: K_0, K_1 from the public scalar routines + forward recurrence.
    ! For nu > 50 we route to the Debye expansion instead.
    ! ----------------------------------------------------------------------------------------------
    elemental real(BK) function besselk_integer(nu, x, scaled) result(y)
        real(BK), intent(in) :: nu, x
        logical,  intent(in) :: scaled
        real(BK) :: k0v, k1v
        integer :: n

        n = nint(nu)
        if (n == 0) then
            if (scaled) then
                y = besselk0x_local(x)
            else
                y = besselk0_local(x)
            end if
            return
        end if
        if (n == 1) then
            if (scaled) then
                y = besselk1x_local(x)
            else
                y = besselk1_local(x)
            end if
            return
        end if

        if (nu > 50.0_BK) then
            if (scaled) then
                y = besselk_debye_scaled(nu, x)
            else
                y = besselk_debye(nu, x)
            end if
            return
        end if

        if (scaled) then
            k0v = besselk0x_local(x)
            k1v = besselk1x_local(x)
        else
            k0v = besselk0_local(x)
            k1v = besselk1_local(x)
        end if
        ! forward recurrence from K_0, K_1 to K_n
        y = besselk_forward_int(x, k0v, k1v, n)
    end function besselk_integer

    ! ----------------------------------------------------------------------------------------------
    ! Non-integer-nu path:
    !   * x > 1.5  : Levin-accelerated asymptotic series for K_{vf}, K_{vf+1} + forward recurrence
    !   * x <= 1.5 + near integer: Temme series + forward recurrence
    !   * x <= 1.5 + general: power series
    ! vf is the fractional part of nu, taken in [0, 1).
    ! ----------------------------------------------------------------------------------------------
    elemental real(BK) function besselk_noninteger(nu, x, scaled) result(y)
        real(BK), intent(in) :: nu, x
        logical,  intent(in) :: scaled
        real(BK) :: vf, kv, kvp1, k_out, k_out_s, vt
        real(BK), parameter :: LEVIN_MIN = 1.5_BK

        vf = nu - real(int(nu), BK)   ! fractional part in [0, 1)

        if (x > LEVIN_MIN) then
            ! Levin acceleration on the asymptotic series of K_{vf}*e^x and K_{vf+1}*e^x.
            kv   = besselkx_levin(vf,       x)
            kvp1 = besselkx_levin(vf + ONE, x)
            k_out_s = besselk_forward(x, kv, kvp1, vf, nu)
            if (scaled) then
                y = k_out_s
            else
                y = k_out_s * exp(-x)
            end if
        else if (is_nearint(vf)) then
            ! near integer; shift to [-0.5, 0.5]
            vt = vf
            if (vt > HALF) vt = vt - ONE
            call besselk_temme_series(vt, x, kv, kvp1)
            k_out = besselk_forward(x, kv, kvp1, vt, nu)
            if (scaled) then
                y = k_out * exp(x)
            else
                y = k_out
            end if
        else
            k_out = besselk_power_series(nu, x)
            if (scaled) then
                y = k_out * exp(x)
            else
                y = k_out
            end if
        end if
    end function besselk_noninteger

    ! Levin sequence transform applied to the asymptotic series
    !   K_v(x) e^x = sqrt(pi/(2x)) * sum_{k>=0} a_k(v) / (8x)^k * k! ...
    ! The partial sums s_i and the next term w_i are fed into a u-Levin transform of length N=16.
    ! See Bessels.jl src/BesselFunctions/besselk.jl + Math/Math.jl.
    pure function besselkx_levin(v, x) result(y)
        real(BK), intent(in) :: v, x
        real(BK) :: y
        integer, parameter :: N = 16
        real(BK) :: s, t, fnu2, scale
        real(BK) :: s_arr(N), w_arr(N), a_num(N), a_den(N)
        integer :: i, k

        fnu2 = FOUR*v*v
        s = ZERO
        t = ONE
        do i = 1, N
            s = s + t
            s_arr(i) = s
            w_arr(i) = t
            t = t * (fnu2 - real((2*i - 1)**2, BK)) / (8.0_BK * x * real(i, BK))
        end do

        do i = 1, N
            if (w_arr(i) == ZERO) then
                y = s_arr(i) * sqrt(PIO2/x)
                return
            end if
            a_num(i) = s_arr(i) / w_arr(i)
            a_den(i) = ONE      / w_arr(i)
        end do

        do k = 1, N - 1
            do i = 1, N - k
                scale = -(real(i + k, BK) * real(i + k - 1, BK)) / &
                         (real(i + 2*k - 1, BK) * real(i + 2*k - 2, BK))
                a_num(i) = a_num(i)*scale + a_num(i+1)
                a_den(i) = a_den(i)*scale + a_den(i+1)
            end do
        end do
        y = (a_num(1)/a_den(1)) * sqrt(PIO2/x)
    end function besselkx_levin

    ! Forward recurrence K_{n+1} = (2n/x) K_n + K_{n-1}.
    ! Given (kv, kvp1) at orders (v_start, v_start+1), recur to nu.
    ! Returns K at order nu.
    elemental real(BK) function besselk_forward(x, kv, kvp1, v_start, nu) result(out)
        real(BK), intent(in) :: x, kv, kvp1, v_start, nu
        real(BK) :: k_prev, k_curr, k_next, nu_curr
        integer :: nstep, i

        nstep = nint(nu - v_start)
        if (nstep <= 0) then
            out = kv
            return
        end if

        k_prev = kv
        k_curr = kvp1
        nu_curr = v_start + ONE

        do i = 1, nstep - 1
            ! advance from K_{nu_curr-1}, K_{nu_curr} to K_{nu_curr}, K_{nu_curr+1}
            k_next = (TWO*nu_curr/x) * k_curr + k_prev
            k_prev = k_curr
            k_curr = k_next
            nu_curr = nu_curr + ONE
        end do
        out = k_curr
    end function besselk_forward

    ! Integer-order forward recurrence specialised to integer n.
    elemental real(BK) function besselk_forward_int(x, k0v, k1v, n) result(out)
        real(BK), intent(in) :: x, k0v, k1v
        integer,  intent(in) :: n
        real(BK) :: k_prev, k_curr, k_next
        integer :: i

        if (n == 0) then; out = k0v; return; end if
        if (n == 1) then; out = k1v; return; end if

        k_prev = k0v
        k_curr = k1v
        do i = 1, n-1
            k_next = (TWO*real(i, BK)/x) * k_curr + k_prev
            k_prev = k_curr
            k_curr = k_next
        end do
        out = k_curr
    end function besselk_forward_int

    ! ----------------------------------------------------------------------------------------------
    ! Power series — only valid when nu is non-integer (sin(pi nu) appears in denominator).
    ! ----------------------------------------------------------------------------------------------
    elemental real(BK) function besselk_power_series(nu, x) result(y)
        real(BK), intent(in) :: nu, x

        integer, parameter :: MaxIter = 200
        real(BK) :: gam, ngam, sinpv, xx, s1, s2, t1, t2, xpv, k_real, num1, num2
        integer :: k

        gam = gamma_BK(nu)
        sinpv = sin(-PI*abs(nu))
        ngam = PI / (sinpv * gam * nu)

        xx = x*x

        s1 = ZERO; s2 = ZERO
        t1 = ONE;  t2 = ONE
        do k = 1, MaxIter
            s1 = s1 + t1
            s2 = s2 + t2
            k_real = real(k, BK)
            num1 = FOUR*k_real*(k_real - nu)
            num2 = FOUR*k_real*(k_real + nu)
            t1 = t1 * xx / num1
            t2 = t2 * xx / num2
            if (abs(t1) < epsilon(ONE) * abs(s1) .and. &
                abs(t2) < epsilon(ONE) * abs(s2)) exit
        end do

        xpv = (HALF*x)**nu
        y = (gam * s1 + xpv*xpv * ngam * s2) / (TWO*xpv)
    end function besselk_power_series

    ! ----------------------------------------------------------------------------------------------
    ! Temme series: computes K_v(x) and K_{v+1}(x) for |v| <= 0.5 (or near zero).
    ! Uses the f_0 local expansion to handle the cancellation as v -> integer.
    ! ----------------------------------------------------------------------------------------------
    elemental subroutine besselk_temme_series(v, x, K_v, K_vp1)
        real(BK), intent(in) :: v, x
        real(BK), intent(out) :: K_v, K_vp1

        integer, parameter :: MaxIter = 500
        ! gamma(1+v) ≈ Taylor around v=0 with given coefficients
        real(BK), parameter :: GAMMA_1PV_TAYLOR(4) = &
            [ONE, -0.5772156649015329_BK, 0.9890559953279725_BK, -0.23263776388631713_BK]
        real(BK) :: z, zz, fk, zv, pk, qk, ck, out_v, out_vp1, term_v, term_vp1
        real(BK) :: k_real
        integer :: k

        z  = HALF*x
        zz = z*z
        fk = f0_local_expansion_v0(v, x)
        zv = z**v

        pk = evalpoly(4,  v, GAMMA_1PV_TAYLOR) / (TWO*zv)
        qk = evalpoly(4, -v, GAMMA_1PV_TAYLOR) * zv * HALF
        ck = ONE
        out_v   = ZERO
        out_vp1 = ZERO

        do k = 1, MaxIter
            term_v   = ck * fk
            term_vp1 = ck * (pk - real(k-1, BK)*fk)
            out_v   = out_v   + term_v
            out_vp1 = out_vp1 + term_vp1
            if (abs(term_v)   < epsilon(ONE)*abs(out_v)   .and. &
                abs(term_vp1) < epsilon(ONE)*abs(out_vp1)) exit
            k_real = real(k, BK)
            fk = (k_real*fk + pk + qk) / (k_real*k_real - v*v)
            pk = pk / (k_real - v)
            qk = qk / (k_real + v)
            ck = ck * zz / k_real
        end do
        K_v   = out_v
        K_vp1 = out_vp1 / z
    end subroutine besselk_temme_series

    elemental real(BK) function f0_local_expansion_v0(v, x) result(f0)
        real(BK), intent(in) :: v, x
        real(BK), parameter :: SP_COEF(4) = &
            [ONE, 1.6449340668482264_BK, 1.8940656589944918_BK, 1.9711021825948702_BK]
        real(BK), parameter :: G1_COEF(3) = &
            [-0.5772156649015329_BK, 0.04200263503409518_BK, 0.042197734555544306_BK]
        real(BK), parameter :: G2_COEF(3) = &
            [ONE, -0.6558780715202539_BK, 0.16653861138229145_BK]
        real(BK), parameter :: SH_COEF(5) = &
            [ONE, 0.16666666666666666_BK, 0.008333333333333333_BK, &
             0.0001984126984126984_BK, 2.7557319223985893e-6_BK]
        real(BK) :: l2dx, mu, vv, sp, g1, g2, sh
        l2dx = log(TWO) - log(x)
        mu   = v*l2dx
        vv   = v*v
        sp = evalpoly(4, vv,    SP_COEF)
        g1 = evalpoly(3, vv,    G1_COEF)
        g2 = evalpoly(3, vv,    G2_COEF)
        sh = evalpoly(5, mu*mu, SH_COEF)
        f0 = sp * (g1 * cosh(mu) + g2 * sh * l2dx)
    end function f0_local_expansion_v0

    ! ----------------------------------------------------------------------------------------------
    ! Local implementations of besselk0/k1 and their scaled variants — duplicating bessels.f90 so
    ! we avoid a circular module dependency.  Coefficients are identical and cross-checked by tests.
    ! ----------------------------------------------------------------------------------------------
    elemental real(BK) function besselk0_local(x) result(y)
        real(BK), intent(in) :: x
        real(BK) :: a, s, x2, rx
        if (x<=ZERO) then
            y = ieee_value(y, ieee_quiet_nan); return
        elseif (x<=ONE) then
            x2 = x*x
            a  = FOURTH*x2
            s  = evalpoly(size(P1_K0),a,P1_K0)/evalpoly(size(Q1_K0),a,Q1_K0) + Y_K0
            y  = evalpoly(size(P2_K0), x2, P2_K0) - (ONE+s*a)*log(x)
        else
            s  = exp(-HALF*x)
            rx = ONE/x
            a  = (evalpoly(size(P3_K0),rx,P3_K0)/evalpoly(size(Q3_K0),rx,Q3_K0) + ONE) * s / sqrt(x)
            y  = a*s
        end if
    end function besselk0_local

    elemental real(BK) function besselk0x_local(x) result(y)
        real(BK), intent(in) :: x
        real(BK) :: a, s, x2, rx
        if (x<=ZERO) then
            y = ieee_value(y, ieee_quiet_nan); return
        elseif (x<=ONE) then
            x2 = x*x
            a  = FOURTH*x2
            s  = evalpoly(size(P1_K0),a,P1_K0)/evalpoly(size(Q1_K0),a,Q1_K0) + Y_K0
            y  = (evalpoly(size(P2_K0), x2, P2_K0) - (ONE+s*a)*log(x)) * exp(x)
        else
            rx = ONE/x
            y  = (evalpoly(size(P3_K0),rx,P3_K0)/evalpoly(size(Q3_K0),rx,Q3_K0) + ONE) / sqrt(x)
        end if
    end function besselk0x_local

    elemental real(BK) function besselk1_local(x) result(y)
        real(BK), intent(in) :: x
        real(BK) :: z, a, pq, rx, s
        if (x<=ZERO) then
            y = ieee_value(y, ieee_quiet_nan); return
        elseif (x<=ONE) then
            z  = x*x
            rx = ONE/x
            a  = FOURTH*z
            pq = evalpoly(size(P1_k1), a, P1_k1)/evalpoly(size(Q1_k1), a, Q1_k1) + Y_K1
            pq = pq*a*a + (HALF*a + ONE)
            a  = HALF*pq*x
            pq = x*evalpoly(size(P2_k1), z, P2_k1)/evalpoly(size(Q2_k1), z, Q2_k1) + rx
            y = a*log(x)+pq
        else
            s  = exp(-HALF*x)
            rx = ONE/x
            a  = evalpoly(size(P3_k1), rx, P3_k1)/evalpoly(size(Q3_k1), rx, Q3_k1) + Y2_K1
            y = a * s*s/sqrt(x)
        end if
    end function besselk1_local

    elemental real(BK) function besselk1x_local(x) result(y)
        real(BK), intent(in) :: x
        real(BK) :: z, a, pq, rx
        if (x<=ZERO) then
            y = ieee_value(y, ieee_quiet_nan); return
        elseif (x<=ONE) then
            z  = x*x
            rx = ONE/x
            a  = FOURTH*z
            pq = evalpoly(size(P1_k1), a, P1_k1)/evalpoly(size(Q1_k1), a, Q1_k1) + Y_K1
            pq = pq*a*a + (HALF*a + ONE)
            a  = HALF*pq*x
            pq = x*evalpoly(size(P2_k1), z, P2_k1)/evalpoly(size(Q2_k1), z, Q2_k1) + rx
            y = (a*log(x)+pq) * exp(x)
        else
            rx = ONE/x
            y  = (evalpoly(size(P3_k1), rx, P3_k1)/evalpoly(size(Q3_k1), rx, Q3_k1) + Y2_K1) / sqrt(x)
        end if
    end function besselk1x_local

    ! Public alias used internally and exported for the besseli reflection path.
    elemental subroutine besselk_up_recurrence(x, kvp1, kv, v_start, v_end, k_out, k_outp1)
        real(BK), intent(in)  :: x, kvp1, kv, v_start, v_end
        real(BK), intent(out) :: k_out, k_outp1
        real(BK) :: k_prev, k_curr, k_next, nu_curr
        integer :: nstep, i

        nstep = nint(v_end - v_start)
        k_prev = kv
        k_curr = kvp1
        nu_curr = v_start + ONE

        if (nstep <= 0) then
            k_out   = kv
            k_outp1 = kvp1
            return
        end if

        do i = 1, nstep - 1
            k_next = (TWO*nu_curr/x) * k_curr + k_prev
            k_prev = k_curr
            k_curr = k_next
            nu_curr = nu_curr + ONE
        end do
        k_out   = k_curr
        ! one extra step to expose K at v_end+1
        k_next = (TWO*nu_curr/x) * k_curr + k_prev
        k_outp1 = k_next
    end subroutine besselk_up_recurrence

end module bessels_besselk
