!  ************************************************************************************************************
!
!                                    ____  ________________ ________   _____
!                                   / __ )/ ____/ ___/ ___// ____/ /  / ___/
!                                  / __  / __/  \__ \\__ \/ __/ / /   \__ \
!                                 / /_/ / /___ ___/ /__/ / /___/ /______/ /
!                                /_____/_____//____/____/_____/_____/____/
!
!                                         Constants and parameters
!
!  MIT License
!
!  Copyright (c) 2022 Federico Perini
!  Copyright (c) 2021-2022 Michael Helton, Oscar Smith, and the Bessels.jl contributors
!
!  ************************************************************************************************************
module bessels_gamma
    use iso_fortran_env, only: real32,real64
    use ieee_arithmetic, only: ieee_quiet_nan,ieee_value
    use bessels_constants
    implicit none
    public

    contains

    ! real64 version adapted from Cephes Mathematical Library (MIT license
    ! https://en.smath.com/view/CephesMathLibrary/license) by Stephen L. Moshier
    elemental real(BK) function gamma_BK(x)
       real(BK), intent(in) :: x

       real(BK) :: ax,s,v,w,res,p,q,z

       real(BK), parameter :: ppoly(8) = [1.000000000000000000009e0_BK,  8.378004301573126728826e-1_BK, &
                                          3.629515436640239168939e-1_BK, 1.113062816019361559013e-1_BK, &
                                          2.385363243461108252554e-2_BK, 4.092666828394035500949e-3_BK, &
                                          4.542931960608009155600e-4_BK, 4.212760487471622013093e-5_BK]

       real(BK), parameter :: qpoly(9) = [ 9.999999999999999999908e-1_BK, 4.150160950588455434583e-1_BK, &
                                          -2.243510905670329164562e-1_BK, -4.633887671244534213831e-2_BK, &
                                           2.773706565840072979165e-2_BK, -7.955933682494738320586e-4_BK, &
                                          -1.237799246653152231188e-3_BK, 2.346584059160635244282e-4_BK, &
                                          -1.397148517476170440917e-5_BK]

       real(BK), parameter :: coefs(10) = [ONE                           ,  8.333333333333331800504e-2_BK, &
                                            3.472222222230075327854e-3_BK, -2.681327161876304418288e-3_BK, &
                                           -2.294719747873185405699e-4_BK,  7.840334842744753003862e-4_BK, &
                                            6.989332260623193171870e-5_BK, -5.950237554056330156018e-4_BK, &
                                           -2.363848809501759061727e-5_BK,  7.147391378143610789273e-4_BK]

       if (x<ZERO) then
          ax = -x
          s  = ax*sin(PI*x)
       else
          ax = x
       endif

       if (ax<=11.5_BK) then

            z = ONE
            do while (ax >= THREE)
                ax = ax-ONE
                z  = z*ax
            end do
            do while (ax < TWO)
                z  = z/ax
                ax = ax+ONE
            end do

            ax = ax-TWO

            ! This polynomial is for (x+2)
            p = evalpoly(size(ppoly), ax, ppoly)
            q = evalpoly(size(qpoly), ax, qpoly)

            res = z*p/q

            if (x<ZERO) then
                gamma_BK = PI/(s*res)
            else
                gamma_BK = res
            end if

       elseif (ax<huge(ax)) then

            w = evalpoly(size(coefs), ONE/ax, coefs)

            ! avoid overflow
            v = ax**muladd(HALF,ax,-FOURTH)

            res = SQ2PI*v*(v/exp(ax))*w

            if (x<ZERO) then
                gamma_BK = PI/(s*res)
            else
                gamma_BK = res
            end if

       else

           gamma_BK = huge(ax)

       end if

   end function gamma_BK
!
!function gamma(_x::Float32)
!    x = Float64(_x)
!    if _x < 0
!        s = sinpi(x)
!        s == 0 && throw(DomainError(_x, "NaN result for non-NaN input."))
!        x = 1 - x
!    end
!    if x < 5
!        z = 1.0
!        while x > 1
!            x -= 1
!            z *= x
!        end
!        num = evalpoly(x, (1.0, 0.41702538904450015, 0.24081703455575904, 0.04071509011391178, 0.015839573267537377))
!        den = x*evalpoly(x, (1.0, 0.9942411061082665, -0.17434932941689474, -0.13577921102050783, 0.03028452206514555))
!        res = z * num / den
!    else
!        x -= 1
!        w = evalpoly(inv(x), (2.506628299028453, 0.20888413086840676, 0.008736513049552962, -0.007022997182153692, 0.0006787969600290756))
!        res = @fastmath sqrt(x) * exp(log(x*1/ℯ) * x) * w
!    end
!    return Float32(_x < 0 ? π / (s * res) : res)
!end
!
!function gamma(_x::Float16)
!    x = Float32(_x)
!    if _x < 0
!        s = sinpi(x)
!        s == 0 && throw(DomainError(_x, "NaN result for non-NaN input."))
!        x = 1 - x
!    end
!    x > 14 && return Float16(ifelse(_x > 0, Inf32, 0f0))
!    z = 1f0
!    while x > 1
!        x -= 1
!        z *= x
!    end
!    num = evalpoly(x, (1.0f0, 0.4170254f0, 0.24081704f0, 0.04071509f0, 0.015839573f0))
!    den = x*evalpoly(x, (1.0f0, 0.9942411f0, -0.17434932f0, -0.13577922f0, 0.030284522f0))
!    return Float16(_x < 0 ? Float32(π)*den / (s*z*num) : z * num / den)
!end
!
!function gamma(n::Integer)
!    n < 0 && throw(DomainError(n, "`n` must not be negative."))
!    n == 0 && return Inf*one(n)
!    n > 20 && return gamma(float(n))
!    @inbounds return Float64(factorial(n-1))
!end
end module bessels_gamma
