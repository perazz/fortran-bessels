!  ************************************************************************************************************
!
!                                    ____  ________________ ________   _____
!                                   / __ )/ ____/ ___/ ___// ____/ /  / ___/
!                                  / __  / __/  \__ \\__ \/ __/ / /   \__ \
!                                 / /_/ / /___ ___/ /__/ / /___/ /______/ /
!                                /_____/_____//____/____/_____/_____/____/
!
!                              Modified Bessel function I_nu(x), variable order
!
!  MIT License
!  Copyright (c) 2022 Federico Perini
!  Copyright (c) 2021-2022 Michael Helton, Oscar Smith, and the Bessels.jl contributors
!
!  Port of Bessels.jl/src/BesselFunctions/besseli.jl.
!    - Large argument asymptotic series for x > nu^2/2 + 19 (DLMF 10.40.E1)
!    - Debye uniform asymptotic for nu > 25 or x > 35 (DLMF 10.41) — shares the U_k polynomial
!      machinery with besselk (Uk_poly10 returns Uk_In as the first split component)
!    - Power series otherwise (DLMF 10.25.E2)
!
!  Negative arguments:
!    integer nu:     I_n(-x) = (-1)^n I_n(x)
!    real    nu:     I_v(-x) is complex-valued — return NaN
!  Negative orders:
!    I_{-nu}(x) = I_nu(x) + (2/pi) sin(pi nu) K_nu(x)   (uses bessels_besselk for K_v)
!  ************************************************************************************************************
module bessels_besseli
    use bessels_constants
    use bessels_debye,    only: Uk_poly10
    use bessels_gamma,    only: gamma_BK
    use bessels_besselk,  only: besselk
    implicit none
    private

    public :: besseli, besselix

    contains

    elemental real(BK) function besseli(nu, x) result(y)
        real(BK), intent(in) :: nu, x
        real(BK) :: anu, ax
        integer  :: n

        if (x /= x) then; y = x; return; end if
        if (.not. is_finite(x)) then; y = x; return; end if

        anu = abs(nu)
        ax  = abs(x)

        if (nu >= ZERO) then
            if (x >= ZERO) then
                y = besseli_positive(anu, ax)
            else
                if (isinteger(anu)) then
                    n = nint(anu)
                    if (mod(n, 2) == 0) then
                        y =  besseli_positive(anu, ax)
                    else
                        y = -besseli_positive(anu, ax)
                    end if
                else
                    y = ieee_value(y, ieee_quiet_nan)
                end if
            end if
        else
            if (x >= ZERO) then
                y = besseli_positive(anu, ax) + TWOOPI*sin(PI*anu)*besselk(anu, ax)
            else
                if (isinteger(anu)) then
                    ! I_{-n}(-x) = (-1)^n I_n(x) * (-1)^n = I_n(x)
                    y = besseli_positive(anu, ax)
                else
                    y = ieee_value(y, ieee_quiet_nan)
                end if
            end if
        end if
    end function besseli

    elemental real(BK) function besselix(nu, x) result(y)
        real(BK), intent(in) :: nu, x
        real(BK) :: anu, ax
        integer  :: n

        if (x /= x) then; y = x; return; end if
        if (.not. is_finite(x)) then; y = x; return; end if

        anu = abs(nu)
        ax  = abs(x)

        if (nu >= ZERO) then
            if (x >= ZERO) then
                y = besselix_positive(anu, ax)
            else
                if (isinteger(anu)) then
                    n = nint(anu)
                    ! I_n(-x) e^{-(-x)} = (-1)^n I_n(x) e^x ; we want e^{-x} I (with x being |-x|... reflect convention)
                    ! Caller passes negative x; besselix(nu, -|x|) = e^{-(-|x|)} I_nu(-|x|) = e^{|x|} I_nu(-|x|).
                    ! For integer n: I_n(-|x|) = (-1)^n I_n(|x|).  So result = (-1)^n * e^{|x|} I_n(|x|).
                    ! This rapidly diverges; users should pass |x|.  We honor the relation though.
                    if (mod(n, 2) == 0) then
                        y =  exp(TWO*ax) * besselix_positive(anu, ax)
                    else
                        y = -exp(TWO*ax) * besselix_positive(anu, ax)
                    end if
                else
                    y = ieee_value(y, ieee_quiet_nan)
                end if
            end if
        else
            if (x >= ZERO) then
                y = besselix_positive(anu, ax) + TWOOPI*sin(PI*anu)*besselk(anu, ax)*exp(-ax)
            else
                if (isinteger(anu)) then
                    y = exp(TWO*ax) * besselix_positive(anu, ax)
                else
                    y = ieee_value(y, ieee_quiet_nan)
                end if
            end if
        end if
    end function besselix

    elemental logical function is_finite(x) result(r)
        real(BK), intent(in) :: x
        r = .not. (abs(x) > huge(x))
    end function is_finite

    ! Positive-(nu, x) dispatcher
    elemental real(BK) function besseli_positive(nu, x) result(y)
        real(BK), intent(in) :: nu, x
        if (besseli_large_args_cutoff(nu, x)) then
            y = besseli_large_args(nu, x)
        elseif (besseli_debye_cutoff(nu, x)) then
            y = besseli_debye(nu, x)
        else
            y = besseli_power_series(nu, x)
        end if
    end function besseli_positive

    elemental real(BK) function besselix_positive(nu, x) result(y)
        real(BK), intent(in) :: nu, x
        if (besseli_large_args_cutoff(nu, x)) then
            y = besselix_large_args(nu, x)
        elseif (besseli_debye_cutoff(nu, x)) then
            y = besseli_debye_scaled(nu, x)
        else
            y = besseli_power_series(nu, x) * exp(-x)
        end if
    end function besselix_positive

    ! ----------------------------------------------------------------------------------------------
    ! Cutoffs (Float64).
    ! ----------------------------------------------------------------------------------------------
    elemental logical function besseli_large_args_cutoff(nu, x)
        real(BK), intent(in) :: nu, x
        besseli_large_args_cutoff = x > HALF*nu*nu + 19.0_BK
    end function besseli_large_args_cutoff

    elemental logical function besseli_debye_cutoff(nu, x)
        real(BK), intent(in) :: nu, x
        besseli_debye_cutoff = (nu > 25.0_BK) .or. (x > 35.0_BK)
    end function besseli_debye_cutoff

    ! ----------------------------------------------------------------------------------------------
    ! Large argument asymptotic series.  Signs alternate (vs. K's same sign).
    ! ----------------------------------------------------------------------------------------------
    elemental real(BK) function besseli_large_args(nu, x) result(y)
        real(BK), intent(in) :: nu, x
        y = exp(HALF*x) * besselix_large_args(nu, x) * exp(HALF*x)
    end function besseli_large_args

    elemental real(BK) function besselix_large_args(nu, x) result(y)
        real(BK), intent(in) :: nu, x
        integer, parameter :: MaxIter = 1000
        real(BK) :: t, s, invx, fnu2
        integer :: i
        invx = ONE/(8.0_BK*x)
        fnu2 = FOUR*nu*nu
        t = ONE
        s = ONE
        do i = 1, MaxIter
            t = -t * invx * ((fnu2 - real((2*i - 1)**2, BK)) / real(i, BK))
            s = s + t
            if (abs(t) < epsilon(ONE) * abs(s)) exit
        end do
        y = s / sqrt(TWO*PI*x)
    end function besselix_large_args

    ! ----------------------------------------------------------------------------------------------
    ! Debye uniform asymptotic (DLMF 10.41).  Same η as the K case; sign in the exponent flipped.
    ! ----------------------------------------------------------------------------------------------
    elemental real(BK) function besseli_debye(nu, x) result(y)
        real(BK), intent(in) :: nu, x
        y = besseli_debye_prefactor(nu, x, scaled=.false.) * besseli_debye_polynomial(nu, x)
    end function besseli_debye

    elemental real(BK) function besseli_debye_scaled(nu, x) result(y)
        real(BK), intent(in) :: nu, x
        y = besseli_debye_prefactor(nu, x, scaled=.true.) * besseli_debye_polynomial(nu, x)
    end function besseli_debye_scaled

    elemental real(BK) function besseli_debye_prefactor(nu, x, scaled) result(coef)
        real(BK), intent(in) :: nu, x
        logical,  intent(in) :: scaled
        real(BK) :: z, zs, n
        z  = x/nu
        zs = sqrt(ONE + z*z)
        n  = zs + log(z) - log(ONE + zs)
        if (scaled) then
            coef = SQ1O2PI * sqrt(ONE/nu) * exp(nu*n - x) / sqrt(zs)
        else
            coef = SQ1O2PI * sqrt(ONE/nu) * exp(nu*n)     / sqrt(zs)
        end if
    end function besseli_debye_prefactor

    elemental real(BK) function besseli_debye_polynomial(nu, x) result(uk)
        real(BK), intent(in) :: nu, x
        real(BK) :: p, p2, mx, mn, Uk_In, Uk_Kn
        p  = ONE/sqrt(ONE + (x/nu)**2)
        mx = max(nu, x)
        mn = min(nu, x)
        p2 = nu*nu / (mx*mx + mn*mn)
        call Uk_poly10(p, nu, p2, Uk_In, Uk_Kn)
        ! Uk_In is the first split output (matches Julia: Uk_poly10(p, v, p2)[1]).
        uk = Uk_In
    end function besseli_debye_polynomial

    ! ----------------------------------------------------------------------------------------------
    ! Power series I_nu(x) = (x/2)^nu / Gamma(nu+1) * sum_{m=0..} ((x/2)^2)^m / (m! * Pochhammer(nu+1, m))
    ! ----------------------------------------------------------------------------------------------
    elemental real(BK) function besseli_power_series(nu, x) result(y)
        real(BK), intent(in) :: nu, x
        integer, parameter :: MaxIter = 3000
        real(BK) :: s, t, xx
        integer :: i
        s = ZERO
        t = ONE
        xx = HALF*HALF*x*x   ! (x/2)^2
        do i = 0, MaxIter
            s = s + t
            if (abs(t) < epsilon(ONE) * abs(s)) exit
            t = t * xx / (real(i+1, BK)*(nu + real(i+1, BK)))
        end do
        y = s * (HALF*x)**nu / gamma_BK(nu + ONE)
    end function besseli_power_series

end module bessels_besseli
