!  ************************************************************************************************************
!
!                                    ____  ________________ ________   _____
!                                   / __ )/ ____/ ___/ ___// ____/ /  / ___/
!                                  / __  / __/  \__ \\__ \/ __/ / /   \__ \
!                                 / /_/ / /___ ___/ /__/ / /___/ /______/ /
!                                /_____/_____//____/____/_____/_____/____/
!
!                                              Airy functions
!
!  MIT License
!
!  Copyright (c) 2022 Federico Perini
!  Copyright (c) 2021-2022 Michael Helton, Oscar Smith, and the Bessels.jl contributors
!
!  Routines for computing Airy functions Ai, Bi, Ai', Bi' and their exponentially-scaled variants
!  for real arguments.  Port of Bessels.jl/src/AiryFunctions/airy.jl which combines:
!    * minimax polynomial / rational approximations (Cephes-style) for |x| small,
!    * tabulated asymptotic expansions for large positive x,
!    * piecewise Taylor series around the Airy zeros for moderate negative x,
!    * an Euler-formula asymptotic that stays in real arithmetic for large negative x.
!  See Jentschura & Lötstedt, Comput. Phys. Commun. 183 (2012) 506-519.
!  ************************************************************************************************************
module bessels_airy
    use bessels_constants
    implicit none
    private

    public :: airyai, airyaiprime, airybi, airybiprime
    public :: airyaix, airyaiprimex, airybix, airybiprimex

    real(BK), parameter :: NEG_LIMIT = -1.0e8_BK   ! beyond which we return NaN (loss of precision)

    contains

    ! Public Ai(x).  NaN beyond the precision-loss limit; NaN for NaN.
    elemental real(BK) function airyai(x) result(y)
        real(BK), intent(in) :: x
        real(BK) :: dummy
        if (x /= x) then
            y = x
        elseif (x >= ZERO) then
            call airyai_aiprime_pos(x, y, dummy)
        elseif (x >= NEG_LIMIT) then
            y = airyai_neg(x)
        else
            y = ieee_value(y, ieee_quiet_nan)
        end if
    end function airyai

    elemental real(BK) function airyaiprime(x) result(y)
        real(BK), intent(in) :: x
        real(BK) :: dummy
        if (x /= x) then
            y = x
        elseif (x >= ZERO) then
            call airyai_aiprime_pos(x, dummy, y)
        elseif (x >= NEG_LIMIT) then
            y = airyaiprime_neg(x)
        else
            y = ieee_value(y, ieee_quiet_nan)
        end if
    end function airyaiprime

    elemental real(BK) function airybi(x) result(y)
        real(BK), intent(in) :: x
        real(BK) :: dummy
        if (x /= x) then
            y = x
        elseif (x >= ZERO) then
            call airybi_biprime_pos(x, y, dummy)
        elseif (x >= NEG_LIMIT) then
            y = airybi_neg(x)
        else
            y = ieee_value(y, ieee_quiet_nan)
        end if
    end function airybi

    elemental real(BK) function airybiprime(x) result(y)
        real(BK), intent(in) :: x
        real(BK) :: dummy
        if (x /= x) then
            y = x
        elseif (x >= ZERO) then
            call airybi_biprime_pos(x, dummy, y)
        elseif (x >= NEG_LIMIT) then
            y = airybiprime_neg(x)
        else
            y = ieee_value(y, ieee_quiet_nan)
        end if
    end function airybiprime

    ! Scaled variants.  For Ai the scaling factor is exp(2/3 * x^(3/2)) (real only for x >= 0);
    ! for Bi the unscaled value overflows for moderate positive x, so the scaled variant exp(-2/3 * x^(3/2)) * Bi
    ! is the practically useful one there.  Negative-x Bi scaled stays bounded so we evaluate it directly.
    elemental real(BK) function airyaix(x) result(y)
        real(BK), intent(in) :: x
        real(BK) :: dummy
        if (x /= x) then
            y = x
        elseif (x >= ZERO) then
            call airyaix_aiprimex_pos(x, y, dummy)
        else
            ! Ai exponentially grows for x<0, no real scaled form; return NaN.
            y = ieee_value(y, ieee_quiet_nan)
        end if
    end function airyaix

    elemental real(BK) function airyaiprimex(x) result(y)
        real(BK), intent(in) :: x
        real(BK) :: dummy
        if (x /= x) then
            y = x
        elseif (x >= ZERO) then
            call airyaix_aiprimex_pos(x, dummy, y)
        else
            y = ieee_value(y, ieee_quiet_nan)
        end if
    end function airyaiprimex

    elemental real(BK) function airybix(x) result(y)
        real(BK), intent(in) :: x
        real(BK) :: dummy
        if (x /= x) then
            y = x
        elseif (x >= ZERO) then
            call airybix_biprimex_pos(x, y, dummy)
        elseif (x >= NEG_LIMIT) then
            ! For x < 0 Bi is bounded; the scaled form equals Bi.
            y = airybi_neg(x)
        else
            y = ieee_value(y, ieee_quiet_nan)
        end if
    end function airybix

    elemental real(BK) function airybiprimex(x) result(y)
        real(BK), intent(in) :: x
        real(BK) :: dummy
        if (x /= x) then
            y = x
        elseif (x >= ZERO) then
            call airybix_biprimex_pos(x, dummy, y)
        elseif (x >= NEG_LIMIT) then
            y = airybiprime_neg(x)
        else
            y = ieee_value(y, ieee_quiet_nan)
        end if
    end function airybiprimex

    ! -----------------------------------------------------------------------------------------------
    ! Positive-x: paired (Ai, Ai')
    ! -----------------------------------------------------------------------------------------------
    elemental subroutine airyai_aiprime_pos(x, ai, aiprime)
        real(BK), intent(in)  :: x
        real(BK), intent(out) :: ai, aiprime
        real(BK) :: c, aix, aiprimex_
        if (x >= 2.06_BK) then
            if (x > huge(x)) then
                ai = ZERO
                aiprime = ZERO
                return
            end if
            c = exp(-x*sqrt(x)/THREE)
            call airyaix_large_pos(x, aix, aiprimex_)
            ai = (c*aix)*c
            aiprime = (c*aiprimex_)*c
        else
            call airyai_small_pos(x, ai, aiprime)
        end if
    end subroutine airyai_aiprime_pos

    elemental subroutine airyaix_aiprimex_pos(x, aix, aiprimex_)
        real(BK), intent(in)  :: x
        real(BK), intent(out) :: aix, aiprimex_
        real(BK) :: c, ai, aiprime
        if (x >= 2.06_BK) then
            call airyaix_large_pos(x, aix, aiprimex_)
        else
            c = exp(TWO*x*sqrt(x)/THREE)
            call airyai_small_pos(x, ai, aiprime)
            aix = c*ai
            aiprimex_ = c*aiprime
        end if
    end subroutine airyaix_aiprimex_pos

    elemental subroutine airyai_small_pos(x, ai, aiprime)
        real(BK), intent(in)  :: x
        real(BK), intent(out) :: ai, aiprime

        real(BK) :: p1, q1, p2, q2

        real(BK), parameter :: ai_low(*)  = [ 0.3550280538878172_BK,  -0.2588194037928022_BK,  -9.261034519371688e-13_BK, &
                                              0.0591713423857866_BK,  -0.021568286307597485_BK, 5.279319501583404e-8_BK, &
                                              0.001971826010638854_BK,-0.0005109862484427906_BK]
        real(BK), parameter :: aip_low(*) = [-0.2588194037928068_BK,   7.026072286612935e-16_BK, 0.1775140269437213_BK, &
                                             -0.08627313457837005_BK, -9.890864876107356e-10_BK, 0.011834297039414418_BK, &
                                             -0.0035951989961820843_BK, 4.799599033246355e-6_BK, 0.0002209248351055058_BK]
        ! 0.1 <= x < 2.06 — Cephes-style rational
        real(BK), parameter :: ai_p(*) = [ 0.35502805388781744_BK,    -0.10482341989661095_BK,   -0.0486577708334799_BK, &
                                           0.027567409888166392_BK,   -0.003289115420423328_BK,  -0.0007305591670017397_BK, &
                                           0.00028980078696335705_BK, -4.0935293322551856e-5_BK,  2.843352302971319e-6_BK, &
                                          -8.16006902427902e-8_BK]
        real(BK), parameter :: ai_q(*) = [ 1.0_BK,                     0.43375722625250424_BK,    0.1791605343840418_BK, &
                                           0.04159189705098047_BK,     0.009514630923507036_BK,   0.0013695863161550792_BK, &
                                           0.00021134924004087647_BK,  1.643834172337603e-5_BK,   1.820647551730948e-6_BK]
        real(BK), parameter :: aip_p(*) = [-0.25881940379280843_BK,   -0.1281771862490461_BK,     0.1215372422831532_BK, &
                                           -0.012741106612258057_BK,  -0.007817713517073376_BK,   0.0024841727869283543_BK, &
                                           -0.00022800908016912754_BK,-1.0310082365339424e-5_BK,  3.234478554177193e-6_BK, &
                                           -1.6817588355348419e-7_BK]
        real(BK), parameter :: aip_q(*) = [ 1.0_BK,                    0.4952379318195294_BK,     0.2162773881682495_BK, &
                                            0.05555863041440847_BK,    0.013462108119986573_BK,   0.0021389588233556925_BK, &
                                            0.0003499286802943533_BK,  3.0276702049000568e-5_BK,  3.5914680052925263e-6_BK]

        if (x < 0.1_BK) then
            ai      = evalpoly(size(ai_low),  x, ai_low)
            aiprime = evalpoly(size(aip_low), x, aip_low)
        else
            p1 = evalpoly(size(ai_p),  x, ai_p)
            q1 = evalpoly(size(ai_q),  x, ai_q)
            p2 = evalpoly(size(aip_p), x, aip_p)
            q2 = evalpoly(size(aip_q), x, aip_q)
            ai      = p1/q1
            aiprime = p2/q2
        end if
    end subroutine airyai_small_pos

    elemental subroutine airyaix_large_pos(x, aix, aiprimex_)
        real(BK), intent(in)  :: x
        real(BK), intent(out) :: aix, aiprimex_

        real(BK) :: invx3, xsqr, xsqrx, zinv, a, c, p, q, p1, q1

        ! "Asymptotic" branch for x > 1000 — small invx3 polynomials
        real(BK), parameter :: P_asy(*)  = [ 1.5707963267948966_BK, 0.13124057851910487_BK,  0.4584353787485384_BK]
        real(BK), parameter :: Q_asy(*)  = [ 0.1636246173744684_BK, 0.20141783231057064_BK,  1.3848568733028765_BK]
        real(BK), parameter :: P1_asy(*) = [-1.5707963267948966_BK, 0.15510250188621483_BK,  0.4982993247266722_BK]
        real(BK), parameter :: Q1_asy(*) = [ 0.22907446432425577_BK,0.22511404787652015_BK,  1.4803642438754887_BK]

        ! Rational approximations (Cephes) for 2.06 <= x <= 1000
        real(BK), parameter :: AI_P(*) = [ 9.99999999999999995305e-1_BK, 1.40264691163389668864e1_BK,  &
                                           7.05360906840444183113e1_BK,  1.59756391350164413639e2_BK,  &
                                           1.68089224934630576269e2_BK,  7.62796053615234516538e1_BK,  &
                                           1.20075952739645805542e1_BK,  3.46538101525629032477e-1_BK]
        real(BK), parameter :: AI_Q(*) = [ 1.00000000000000000470e0_BK,  1.40959135607834029598e1_BK,  &
                                           7.14778400825575695274e1_BK,  1.64234692871529701831e2_BK,  &
                                           1.77318088145400459522e2_BK,  8.45138970141474626562e1_BK,  &
                                           1.47562562584847203173e1_BK,  5.67594532638770212846e-1_BK]
        real(BK), parameter :: AIP_P(*) = [1.00000000000000000550e0_BK,  1.39470856980481566958e1_BK,  &
                                           6.99778599330103016170e1_BK,  1.59317847137141783523e2_BK,  &
                                           1.71184781360976385540e2_BK,  8.20584123476060982430e1_BK,  &
                                           1.47454670787755323881e1_BK,  6.13759184814035759225e-1_BK]
        real(BK), parameter :: AIP_Q(*) = [9.99999999999999994502e-1_BK, 1.38498634758259442477e1_BK,  &
                                           6.86752304592780337944e1_BK,  1.53206427475809220834e2_BK,  &
                                           1.58778084372838313640e2_BK,  7.11727352147859965283e1_BK,  &
                                           1.11810297306158156705e1_BK,  3.34203677749736953049e-1_BK]

        if (x > 1000.0_BK) then
            invx3 = ONE/(x*x*x)
            p  = evalpoly(size(P_asy),  invx3, P_asy)
            q  = evalpoly(size(Q_asy),  invx3, Q_asy)
            p1 = evalpoly(size(P1_asy), invx3, P1_asy)
            q1 = evalpoly(size(Q1_asy), invx3, Q1_asy)
            xsqr = sqrt(x)
            xsqrx = ONE/(x*xsqr)
            a = muladd(xsqrx, -q, p)
            c = muladd(xsqrx, -q1, p1)
            xsqr = sqrt(xsqr)
            aix = a/(PIPOW3O2*xsqr)
            aiprimex_ = c*xsqr/PIPOW3O2
        else
            xsqr  = sqrt(x)
            zinv  = THREE/(TWO*x*xsqr)
            xsqr  = sqrt(xsqr)
            p  = evalpoly(size(AI_P),  zinv, AI_P)
            q  = evalpoly(size(AI_Q),  zinv, AI_Q)
            p1 = evalpoly(size(AIP_P), zinv, AIP_P)
            q1 = evalpoly(size(AIP_Q), zinv, AIP_Q)
            aix       = (ONEOSQPI*HALF)*p/(q*xsqr)
            aiprimex_ = -(ONEOSQPI*HALF)*xsqr*p1/q1
        end if
    end subroutine airyaix_large_pos

    ! -----------------------------------------------------------------------------------------------
    ! Negative-x: Ai
    ! -----------------------------------------------------------------------------------------------
    elemental real(BK) function airyai_neg(x) result(y)
        real(BK), intent(in) :: x
        real(BK) :: ai, aip
        if (x < -10.0_BK) then
            call airyai_large_neg(x, ai, aip)
            y = ai
        elseif (x > -1.0_BK) then
            y = evalpoly_taylor_ai_neg_minus_half(x + 0.5_BK)
        elseif (x > -3.2_BK) then
            y = evalpoly_taylor_ai_neg_a1(x + 2.338107410459767_BK)
        elseif (x > -4.75_BK) then
            y = evalpoly_taylor_ai_neg_a2(x + 4.08794944413097_BK)
        elseif (x > -6.1_BK) then
            y = evalpoly_taylor_ai_neg_a3(x + 5.520559828095551_BK)
        elseif (x > -7.2_BK) then
            y = evalpoly_taylor_ai_neg_a4(x + 6.786708090071759_BK)
        elseif (x > -8.5_BK) then
            y = evalpoly_taylor_ai_neg_a5(x + 7.944133587120853_BK)
        else
            y = evalpoly_taylor_ai_neg_a6(x + 9.02265085334098_BK)
        end if
    end function airyai_neg

    elemental real(BK) function airyaiprime_neg(x) result(y)
        real(BK), intent(in) :: x
        real(BK) :: ai, aip
        if (x < -10.0_BK) then
            call airyai_large_neg(x, ai, aip)
            y = aip
        elseif (x > -2.1_BK) then
            y = evalpoly_taylor_aip_neg_b1(x + 1.018792971647471_BK)
        elseif (x > -4.03_BK) then
            y = evalpoly_taylor_aip_neg_b2(x + 3.2481975821798366_BK)
        elseif (x > -5.5_BK) then
            y = evalpoly_taylor_aip_neg_b3(x + 4.820099211178736_BK)
        elseif (x > -6.8_BK) then
            y = evalpoly_taylor_aip_neg_b4(x + 6.163307355639486_BK)
        elseif (x > -7.9_BK) then
            y = evalpoly_taylor_aip_neg_b5(x + 7.37217725504777_BK)
        elseif (x > -9.0_BK) then
            y = evalpoly_taylor_aip_neg_b6(x + 8.488486734019721_BK)
        else
            y = evalpoly_taylor_aip_neg_b7(x + 9.535449052433547123541757173_BK)
        end if
    end function airyaiprime_neg

    elemental subroutine airyai_large_neg(x, ai, aip)
        real(BK), intent(in)  :: x
        real(BK), intent(out) :: ai, aip

        real(BK) :: invx3, xabs, xsqr, xsqrx, z, spc, smc, a, c, arg, p, q, p1, q1

        real(BK), parameter :: P_l(*)  = [ 1.5707963267948966_BK, 0.13124057851910487_BK,  0.4584353787485384_BK, &
                                           5.217255928936184_BK,  123.97197893818594_BK,   5038.313653002081_BK, &
                                           312467.7049060495_BK,  2.746439545069411e7_BK,  3.2482560591146026e9_BK, &
                                           4.97462635569055e11_BK,9.57732308323407e13_BK]
        real(BK), parameter :: Q_l(*)  = [ 0.1636246173744684_BK, 0.20141783231057064_BK,  1.3848568733028765_BK, &
                                           23.555289417250567_BK, 745.2667344964557_BK,    37835.063701047824_BK, &
                                           2.8147130917899106e6_BK, 2.8856687720069575e8_BK, 3.8998976239149216e10_BK, &
                                           6.718472897263214e12_BK, 1.4370735281142392e15_BK]
        real(BK), parameter :: P1_l(*) = [-1.5707963267948966_BK, 0.15510250188621483_BK,  0.4982993247266722_BK, &
                                           5.515384839161109_BK,  129.24738229725767_BK,   5209.103946324185_BK, &
                                           321269.61208650155_BK, 2.812618811215662e7_BK,  3.3166403972012258e9_BK, &
                                           5.0676100258903735e11_BK, 9.738286496397669e13_BK]
        real(BK), parameter :: Q1_l(*) = [ 0.22907446432425577_BK,0.22511404787652015_BK,  1.4803642438754887_BK, &
                                           24.70432792540913_BK,  773.390007496322_BK,     38999.21950723391_BK, &
                                           2.8878225227454924e6_BK, 2.950515261265541e8_BK, 3.97712331943799e10_BK, &
                                           6.837383921993536e12_BK, 1.460066704564067e15_BK]

        invx3 = ONE/(x*x*x)
        p  = evalpoly(size(P_l),  invx3, P_l)
        q  = evalpoly(size(Q_l),  invx3, Q_l)
        p1 = evalpoly(size(P1_l), invx3, P1_l)
        q1 = evalpoly(size(Q1_l), invx3, Q1_l)

        xabs  = -x
        xsqr  = sqrt(xabs)
        xsqrx = xabs*xsqr
        z     = -TWO*xsqrx/THREE

        arg = modulo(z, TWO*PI) + PIO4
        spc = sin(arg)
        smc = cos(arg)
        ! spc = sin(arg), smc = cos(arg).  Original convention: spc = sin(z)+cos(z)
        ! is encoded as sin(z+pi/4)*sqrt(2); the sincos call returns those rotated values.
        ! Convert: sin(z)+cos(z) = sqrt(2)*sin(z+pi/4),  sin(z)-cos(z) = -sqrt(2)*cos(z+pi/4).
        ! Below we use the same combinations as the Julia code which absorbed the sqrt(2).
        a = -p*smc + q/xsqrx*spc
        c =  p1*spc + q1/xsqrx*smc

        xsqr = sqrt(xsqr)
        ai  = -TWO/PIPOW3O2*a/xsqr
        aip =  TWO/PIPOW3O2*c*xsqr
    end subroutine airyai_large_neg

    ! -----------------------------------------------------------------------------------------------
    ! Positive-x: paired (Bi, Bi')
    ! -----------------------------------------------------------------------------------------------
    elemental subroutine airybi_biprime_pos(x, bi, biprime)
        real(BK), intent(in)  :: x
        real(BK), intent(out) :: bi, biprime
        real(BK) :: c, bix, biprimex_

        real(BK), parameter :: bi_low(*) = [ 0.6149266274460008_BK,  0.4482883573538208_BK,   3.5122876076055943e-13_BK, &
                                             0.10248777123224333_BK, 0.03735736322749984_BK, -9.112176579053885e-10_BK, &
                                             0.003416263825552168_BK,0.0008894434900642912_BK,4.656589540929903e-8_BK, &
                                             4.7356203781096e-5_BK,  1.0019617843212616e-5_BK,-1.5501028624002673e-7_BK, &
                                             4.936240473816792e-7_BK,-2.4890774184202023e-8_BK, 4.351871351672967e-8_BK, &
                                            -1.397416816153122e-8_BK, 4.207858977853583e-9_BK, -6.320382000910312e-10_BK, &
                                             5.7683713281742576e-11_BK]
        real(BK), parameter :: bip_low(*) = [0.4482883573538264_BK,  -1.3251373159736958e-14_BK, 0.30746331372383734_BK, &
                                             0.1494294524303411_BK,   2.758926834393856e-10_BK, 0.020497552037528_BK, &
                                             0.006226238906796048_BK,-4.342796143830116e-8_BK,  0.0004271490047546154_BK, &
                                             9.859632186626943e-5_BK, 3.5043952855799877e-7_BK, 3.911502040592847e-6_BK, &
                                             1.1754828402935787e-6_BK,-2.3404028248712297e-7_BK,1.4233554924731257e-7_BK, &
                                            -3.8239033234612784e-8_BK,1.0757888521998427e-8_BK,-1.6237703662454916e-9_BK, &
                                             1.5114681466422037e-10_BK]

        if (x <= TWO) then
            bi      = evalpoly(size(bi_low),  x, bi_low)
            biprime = evalpoly(size(bip_low), x, bip_low)
        else
            if (x > huge(x)) then
                bi = x
                biprime = x
                return
            end if
            c = exp(x*sqrt(x)/THREE)
            if (x <= 6.0_BK) then
                bix = bix_cheb_6(x)
                biprimex_ = biprimex_cheb_6(x)
            elseif (x <= 10.0_BK) then
                bix = bix_cheb_10(x)
                biprimex_ = biprimex_cheb_10(x)
            else
                call airybix_large_pos(x, bix, biprimex_)
            end if
            bi      = (c*bix)*c
            biprime = (c*biprimex_)*c
        end if
    end subroutine airybi_biprime_pos

    elemental subroutine airybix_biprimex_pos(x, bix, biprimex_)
        real(BK), intent(in)  :: x
        real(BK), intent(out) :: bix, biprimex_
        real(BK) :: c, bi, biprime

        ! Reuse the small-x polynomials and downscale
        real(BK), parameter :: bi_low(*) = [ 0.6149266274460008_BK,  0.4482883573538208_BK,   3.5122876076055943e-13_BK, &
                                             0.10248777123224333_BK, 0.03735736322749984_BK, -9.112176579053885e-10_BK, &
                                             0.003416263825552168_BK,0.0008894434900642912_BK,4.656589540929903e-8_BK, &
                                             4.7356203781096e-5_BK,  1.0019617843212616e-5_BK,-1.5501028624002673e-7_BK, &
                                             4.936240473816792e-7_BK,-2.4890774184202023e-8_BK, 4.351871351672967e-8_BK, &
                                            -1.397416816153122e-8_BK, 4.207858977853583e-9_BK, -6.320382000910312e-10_BK, &
                                             5.7683713281742576e-11_BK]
        real(BK), parameter :: bip_low(*) = [0.4482883573538264_BK,  -1.3251373159736958e-14_BK, 0.30746331372383734_BK, &
                                             0.1494294524303411_BK,   2.758926834393856e-10_BK, 0.020497552037528_BK, &
                                             0.006226238906796048_BK,-4.342796143830116e-8_BK,  0.0004271490047546154_BK, &
                                             9.859632186626943e-5_BK, 3.5043952855799877e-7_BK, 3.911502040592847e-6_BK, &
                                             1.1754828402935787e-6_BK,-2.3404028248712297e-7_BK,1.4233554924731257e-7_BK, &
                                            -3.8239033234612784e-8_BK,1.0757888521998427e-8_BK,-1.6237703662454916e-9_BK, &
                                             1.5114681466422037e-10_BK]
        if (x <= TWO) then
            c = exp(-x*sqrt(x)/THREE)
            bi = evalpoly(size(bi_low),  x, bi_low)
            biprime = evalpoly(size(bip_low), x, bip_low)
            bix = (c*bi)*c
            biprimex_ = (c*biprime)*c
        else
            if (x <= 6.0_BK) then
                bix = bix_cheb_6(x)
                biprimex_ = biprimex_cheb_6(x)
            elseif (x <= 10.0_BK) then
                bix = bix_cheb_10(x)
                biprimex_ = biprimex_cheb_10(x)
            else
                call airybix_large_pos(x, bix, biprimex_)
            end if
        end if
    end subroutine airybix_biprimex_pos

    pure real(BK) function bix_cheb_6(x) result(y)
        real(BK), intent(in) :: x
        real(BK), parameter :: c(*) = [ 0.41739858494452137_BK, -0.06495162894671667_BK,  0.01346635223349219_BK, &
                                       -0.0034498974696635805_BK, 0.0009181061897552018_BK,-0.0002204042930330767_BK, &
                                        3.7976551651085086e-5_BK,-1.3803071542564344e-8_BK,-3.398553997511632e-6_BK, &
                                        1.6746478140542785e-6_BK,-5.176623166786395e-7_BK,  1.0972618211669806e-7_BK, &
                                       -1.2473514618542336e-8_BK,-1.4217711195129223e-9_BK, 1.146316775379731e-9_BK, &
                                       -3.175280659616182e-10_BK, 5.034723214458251e-11_BK,-2.0319170143404405e-12_BK, &
                                       -1.4129596145607766e-12_BK, 4.802366525142762e-13_BK,-8.22183850663261e-14_BK, &
                                        5.914052924936658e-15_BK, 1.0066800053200003e-15_BK,-3.986773142283109e-16_BK, &
                                        6.067332401977342e-17_BK]
        y = clenshaw_chebyshev(muladd(x, HALF, -TWO), c)
    end function bix_cheb_6

    pure real(BK) function biprimex_cheb_6(x) result(y)
        real(BK), intent(in) :: x
        real(BK), parameter :: c(*) = [ 0.7658021754601785_BK,    0.1218178967795153_BK,    -0.016874295932807537_BK, &
                                        0.0038022995437998084_BK,-0.0009201686158669025_BK,  0.00018122863592354458_BK, &
                                       -1.089109598429424e-5_BK, -1.2116775872566089e-5_BK,  7.326775924782161e-6_BK, &
                                       -2.5316832797383597e-6_BK, 5.814914990161355e-7_BK,  -6.775973510473317e-8_BK, &
                                       -1.0915897862521443e-8_BK, 8.299026220548351e-9_BK,  -2.4078260629342257e-9_BK, &
                                        4.037468537238015e-10_BK,-1.7933493682633233e-11_BK,-1.216289673737319e-11_BK, &
                                        4.331771122129753e-12_BK,-7.552310555170731e-13_BK,  4.643024726419001e-14_BK, &
                                        1.4365062751453502e-14_BK,-5.2995707246535766e-15_BK,8.910837004634583e-16_BK, &
                                       -6.308304600047433e-17_BK]
        y = clenshaw_chebyshev(muladd(x, HALF, -TWO), c)
    end function biprimex_cheb_6

    pure real(BK) function bix_cheb_10(x) result(y)
        real(BK), intent(in) :: x
        real(BK), parameter :: c(*) = [ 0.3388998275567636_BK,    -0.022239535106198405_BK,   0.001847776356735787_BK, &
                                       -0.00018642620881156231_BK, 2.0614467307054072e-5_BK, -2.413394583210988e-6_BK, &
                                        2.9461272726896696e-7_BK, -3.722819620752169e-8_BK,   4.854723628008316e-9_BK, &
                                       -6.531662239526833e-10_BK,  9.08219420289064e-11_BK,  -1.3082391956275184e-11_BK, &
                                        1.953106497996366e-12_BK, -3.007137783175837e-13_BK,  4.709451540463225e-14_BK, &
                                       -7.327307135432242e-15_BK,  1.0856338884553918e-15_BK,-1.4783310001069625e-16_BK]
        y = clenshaw_chebyshev(muladd(x, HALF, -FOUR), c)
    end function bix_cheb_10

    pure real(BK) function biprimex_cheb_10(x) result(y)
        real(BK), intent(in) :: x
        real(BK), parameter :: c(*) = [ 0.9393816481647195_BK,     0.0621197721009847_BK,    -0.0031821373579424077_BK, &
                                        0.0002596583618312277_BK, -2.5490889161776127e-5_BK,  2.7833396581896174e-6_BK, &
                                       -3.270630123629566e-7_BK,   4.070105276526434e-8_BK,  -5.322537351397486e-9_BK, &
                                        7.293189148846469e-10_BK, -1.0468150599102346e-10_BK, 1.5727589039481186e-11_BK, &
                                       -2.460184674017499e-12_BK,  3.9506685986702333e-13_BK,-6.352712208186561e-14_BK, &
                                        9.868492209739352e-15_BK, -1.3994310387891848e-15_BK, 1.6674714479000505e-16_BK]
        y = clenshaw_chebyshev(muladd(x, HALF, -FOUR), c)
    end function biprimex_cheb_10

    elemental subroutine airybix_large_pos(x, bix, biprimex_)
        real(BK), intent(in)  :: x
        real(BK), intent(out) :: bix, biprimex_

        real(BK) :: invx3, xsqr, xsqrx, p, p1, q, q1

        real(BK), parameter :: P_l(*)  = [ 1.5707963267948966_BK, 0.13124057851910487_BK,  0.4584353787485384_BK, &
                                           5.217255928936184_BK,  123.97197893818594_BK,   5038.313653002081_BK, &
                                           312467.7049060495_BK,  2.746439545069411e7_BK,  3.2482560591146026e9_BK, &
                                           4.97462635569055e11_BK,9.57732308323407e13_BK,  2.2640712393216476e16_BK]
        real(BK), parameter :: P1_l(*) = [ 0.1636246173744684_BK, 0.20141783231057064_BK,  1.3848568733028765_BK, &
                                           23.555289417250567_BK, 745.2667344964557_BK,    37835.063701047824_BK, &
                                           2.8147130917899106e6_BK, 2.8856687720069575e8_BK, 3.8998976239149216e10_BK, &
                                           6.718472897263214e12_BK, 1.4370735281142392e1_BK, 3.7367429394637446e17_BK]
        real(BK), parameter :: Q_l(*)  = [-1.5707963267948966_BK, 0.15510250188621483_BK,  0.4982993247266722_BK, &
                                           5.515384839161109_BK,  129.24738229725767_BK,   5209.103946324185_BK, &
                                           321269.61208650155_BK, 2.812618811215662e7_BK,  3.3166403972012258e9_BK, &
                                           5.0676100258903735e11_BK, 9.738286496397669e13_BK, 2.298637212441062e16_BK]
        real(BK), parameter :: Q1_l(*) = [ 0.22907446432425577_BK,0.22511404787652015_BK,  1.4803642438754887_BK, &
                                           24.70432792540913_BK,  773.390007496322_BK,     38999.21950723391_BK, &
                                           2.8878225227454924e6_BK, 2.950515261265541e8_BK, 3.97712331943799e10_BK, &
                                           6.837383921993536e12_BK, 1.460066704564067e15_BK, 3.7912939312807334e17_BK]

        invx3 = ONE/(x*x*x)
        p  = evalpoly(size(P_l),  invx3, P_l)
        p1 = evalpoly(size(P1_l), invx3, P1_l)
        q  = evalpoly(size(Q_l),  invx3, Q_l)
        q1 = evalpoly(size(Q1_l), invx3, Q1_l)

        xsqr  = sqrt(x)
        xsqrx = ONE/(x*xsqr)
        xsqr  = sqrt(xsqr)
        bix       =  TWO/PIPOW3O2 * muladd(xsqrx, p1, p) / xsqr
        biprimex_ = -TWO/PIPOW3O2 * muladd(xsqrx, q1, q) * xsqr
    end subroutine airybix_large_pos

    ! -----------------------------------------------------------------------------------------------
    ! Negative-x: Bi and Bi'
    ! -----------------------------------------------------------------------------------------------
    elemental real(BK) function airybi_neg(x) result(y)
        real(BK), intent(in) :: x
        real(BK) :: bi, bip, x3
        real(BK), parameter :: bi_pow_a(*) = [ &
                  0.6149266274460007_BK,    0.10248777124100013_BK,   0.0034162590413666706_BK, &
                  4.7448042241203764e-5_BK, 3.5945486546366485e-7_BK, 1.7116898355412612e-9_BK, &
                  5.5937576324877815e-12_BK,1.3318470553542337e-14_BK,2.412766404627235e-17_BK, &
                  3.4369891803806764e-20_BK,3.9505622762996283e-23_BK,3.7410627616473754e-26_BK, &
                  2.9690974298788695e-29_BK,2.0034395613217741e-32_BK,1.163437608200798e-35_BK, &
                  5.875947516165647e-39_BK, 2.604586664967042e-42_BK, 1.0214065352811929e-45_BK, &
                  3.568855818592568e-49_BK, 1.1180625998097016e-52_BK,3.1583689260161066e-56_BK, &
                  8.08594195088609e-60_BK,  1.884834953586501e-63_BK]
        real(BK), parameter :: bi_pow_b(*) = [ &
                  0.4482883573538264_BK,    0.03735736311281886_BK,   0.0008894610264956872_BK, &
                  9.882900294396524e-6_BK,  6.335192496408029e-8_BK,  2.6396635401700117e-10_BK, &
                  7.718314444941556e-13_BK, 1.6706308322384318e-15_BK,2.7843847203973864e-18_BK, &
                  3.683048571954215e-21_BK, 3.960267281671199e-24_BK, 3.52964998366417e-27_BK, &
                  2.6498873751232507e-30_BK,1.698645753284135e-33_BK, 9.405568955061657e-37_BK, &
                  4.5437531183872734e-40_BK,1.9318678224435686e-43_BK,7.284569466227635e-47_BK, &
                  2.4527169919958367e-50_BK,7.418986666654073e-54_BK, 2.0270455373371784e-57_BK, &
                  5.0273946858560976e-61_BK,1.136905175453663e-64_BK]

        if (x < -10.0_BK) then
            call airybi_large_neg(x, bi, bip)
            y = bi
        elseif (x > -0.2_BK) then
            x3 = x*x*x
            bi  = evalpoly(size(bi_pow_a), x3, bi_pow_a)
            bip = evalpoly(size(bi_pow_b), x3, bi_pow_b)
            y = muladd(x, bip, bi)
        elseif (x > -2.2_BK) then
            y = evalpoly_taylor_bi_neg_c1(x + 1.173713222709128_BK)
        elseif (x > -4.05_BK) then
            y = evalpoly_taylor_bi_neg_c2(x + 3.271093302836353_BK)
        elseif (x > -5.5_BK) then
            y = evalpoly_taylor_bi_neg_c3(x + 4.830737841662016_BK)
        elseif (x > -6.7_BK) then
            y = evalpoly_taylor_bi_neg_c4(x + 6.169852128310251_BK)
        elseif (x > -7.9_BK) then
            y = evalpoly_taylor_bi_neg_c5(x + 7.376762079367763_BK)
        elseif (x > -9.0_BK) then
            y = evalpoly_taylor_bi_neg_c6(x + 8.491948846509388_BK)
        else
            y = evalpoly_taylor_bi_neg_c7(x + 9.538194379346239_BK)
        end if
    end function airybi_neg

    elemental real(BK) function airybiprime_neg(x) result(y)
        real(BK), intent(in) :: x
        real(BK) :: bi, bip, x3, x2
        real(BK), parameter :: bip_pow_a(*) = [ &
                  0.4482883573538264_BK,    0.14942945245127545_BK,   0.0062262271854698105_BK, &
                  9.882900294396525e-5_BK,  8.235750245330437e-7_BK,  4.223461664272019e-9_BK, &
                  1.4664797445388956e-11_BK,3.67538783092455e-14_BK,  6.960961800993466e-17_BK, &
                  1.0312536001471802e-19_BK,1.2276828573180715e-22_BK,1.2000809944458178e-25_BK, &
                  9.804583287956028e-29_BK, 6.79458301313654e-32_BK,  4.0443946506765124e-35_BK, &
                  2.0901264344581458e-38_BK,9.466152329973487e-42_BK, 3.78797612243837e-45_BK, &
                  1.3489943455977101e-48_BK,4.3030122666593625e-52_BK,1.236497777775679e-55_BK, &
                  3.2175325989479025e-59_BK,7.617264675539542e-63_BK]
        real(BK), parameter :: bip_pow_b(*) = [ &
                  0.30746331372300034_BK,   0.020497554248200024_BK,  0.0004270323801708338_BK, &
                  4.313458385563978e-6_BK,  2.567534753311892e-8_BK,  1.0068763738478007e-10_BK, &
                  2.796878816243891e-13_BK, 5.790639371105364e-16_BK, 9.279870787027826e-19_BK, &
                  1.1851686828898885e-21_BK,1.2345507113436338e-24_BK,1.068875074756393e-27_BK, &
                  7.813414289154919e-31_BK, 4.8864379544433515e-34_BK,2.6441763822745408e-37_BK, &
                  1.2502015991841802e-40_BK,5.209173329934083e-44_BK, 1.9271821420399865e-47_BK, &
                  6.372956818915299e-51_BK, 1.895021355609664e-54_BK, 5.0941434290582364e-58_BK, &
                  1.2439910693670906e-61_BK,2.771816108215443e-65_BK]

        if (x < -10.0_BK) then
            call airybi_large_neg(x, bi, bip)
            y = bip
        elseif (x > -1.5_BK) then
            x2 = x*x
            x3 = x*x2
            bi  = evalpoly(size(bip_pow_a), x3, bip_pow_a)
            bip = evalpoly(size(bip_pow_b), x3, bip_pow_b)
            y = muladd(x2, bip, bi)
        elseif (x > -3.2_BK) then
            y = evalpoly_taylor_bip_neg_d1(x + 2.294439682614123_BK)
        elseif (x > -4.8_BK) then
            y = evalpoly_taylor_bip_neg_d2(x + 4.073155089071828_BK)
        elseif (x > -6.1_BK) then
            y = evalpoly_taylor_bip_neg_d3(x + 5.5123957296635995_BK)
        elseif (x > -7.4_BK) then
            y = evalpoly_taylor_bip_neg_d4(x + 6.781294445990305_BK)
        elseif (x > -8.5_BK) then
            y = evalpoly_taylor_bip_neg_d5(x + 7.940178689168579_BK)
        else
            y = evalpoly_taylor_bip_neg_d6(x + 9.01958335879424_BK)
        end if
    end function airybiprime_neg

    elemental subroutine airybi_large_neg(x, bi, bip)
        real(BK), intent(in)  :: x
        real(BK), intent(out) :: bi, bip

        real(BK) :: invx3, xabs, xsqr, xsqrx, z, spc, smc, b, d, arg, p, p1, q, q1

        real(BK), parameter :: P_l(*)  = [ 1.5707963267948966_BK, 0.13124057851910487_BK,  0.4584353787485384_BK, &
                                           5.217255928936184_BK,  123.97197893818594_BK,   5038.313653002081_BK, &
                                           312467.7049060495_BK,  2.746439545069411e7_BK,  3.2482560591146026e9_BK, &
                                           4.97462635569055e11_BK,9.57732308323407e13_BK,  2.2640712393216476e16_BK]
        real(BK), parameter :: P1_l(*) = [ 0.1636246173744684_BK, 0.20141783231057064_BK,  1.3848568733028765_BK, &
                                           23.555289417250567_BK, 745.2667344964557_BK,    37835.063701047824_BK, &
                                           2.8147130917899106e6_BK, 2.8856687720069575e8_BK, 3.8998976239149216e10_BK, &
                                           6.718472897263214e12_BK, 1.4370735281142392e1_BK, 3.7367429394637446e17_BK]
        real(BK), parameter :: Q_l(*)  = [-1.5707963267948966_BK, 0.15510250188621483_BK,  0.4982993247266722_BK, &
                                           5.515384839161109_BK,  129.24738229725767_BK,   5209.103946324185_BK, &
                                           321269.61208650155_BK, 2.812618811215662e7_BK,  3.3166403972012258e9_BK, &
                                           5.0676100258903735e11_BK, 9.738286496397669e13_BK, 2.298637212441062e16_BK]
        real(BK), parameter :: Q1_l(*) = [ 0.22907446432425577_BK,0.22511404787652015_BK,  1.4803642438754887_BK, &
                                           24.70432792540913_BK,  773.390007496322_BK,     38999.21950723391_BK, &
                                           2.8878225227454924e6_BK, 2.950515261265541e8_BK, 3.97712331943799e10_BK, &
                                           6.837383921993536e12_BK, 1.460066704564067e15_BK, 3.7912939312807334e17_BK]

        invx3 = ONE/(x*x*x)
        p  = evalpoly(size(P_l),  invx3, P_l)
        p1 = evalpoly(size(P1_l), invx3, P1_l)
        q  = evalpoly(size(Q_l),  invx3, Q_l)
        q1 = evalpoly(size(Q1_l), invx3, Q1_l)

        xabs  = -x
        xsqr  = sqrt(xabs)
        xsqrx = xabs*xsqr
        z     = TWO*xsqrx/THREE

        arg = modulo(z, TWO*PI) + PIO4
        spc = sin(arg)
        smc = cos(arg)
        b = p1*spc/xsqrx + p*smc
        d = q*spc - q1*smc/xsqrx

        xsqr = sqrt(xsqr)
        bi  =  TWO/PIPOW3O2 * b / xsqr
        bip = -TWO/PIPOW3O2 * d * xsqr
    end subroutine airybi_large_neg

    ! -----------------------------------------------------------------------------------------------
    ! Taylor series for Ai/Ai' around the first several zeros
    ! Argument here is (x - root); the polynomial in x_shift gives the value directly.
    ! -----------------------------------------------------------------------------------------------
    pure real(BK) function evalpoly_taylor_ai_neg_minus_half(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
            0.4757280916105396_BK, -0.20408167033954738_BK, -0.1189320229026349_BK,    0.09629482113005221_BK, &
           -0.012051304907352494_BK,-0.00835397167338305_BK,  0.0034106824527909488_BK,-0.00018748378739668977_BK, &
           -0.00017963058749604507_BK, 4.867256036790685e-5_BK,-1.0852054849851912e-6_BK,-1.85424425163635e-6_BK, &
            3.728421447757534e-7_BK, -1.0133548664552323e-9_BK,-1.1212446835297949e-8_BK, 1.777851534328481e-9_BK, &
            1.9136952296640593e-11_BK,-4.449034045022864e-11_BK, 5.778702804510329e-12_BK, 1.2100035825074537e-13_BK, &
           -1.2468339961179948e-13_BK, 1.3614768155678468e-14_BK, 3.9684428150788985e-16_BK,-2.598632088728038e-16_BK, &
            2.430497466471834e-17_BK,  8.779598099071528e-19_BK,-4.184856864694815e-19_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_ai_neg_minus_half

    pure real(BK) function evalpoly_taylor_ai_neg_a1(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
            2.743319340666283e-17_BK, 0.7012108227206914_BK,    -3.207087639834719e-17_BK,-0.2732510368163064_BK, &
            0.058434235226724286_BK,  0.03194451370480103_BK,    -0.01366255184081532_BK, -0.0003870349759524901_BK, &
            0.0011408754894571796_BK, -0.00017718920132546787_BK,-3.3939160136269284e-5_BK, 1.4137844310270033e-5_BK, &
           -7.41180299288454e-7_BK,   -4.2945486337203914e-7_BK,  8.72022168160613e-8_BK,  1.252053805808083e-9_BK, &
           -2.6389292196591304e-9_BK,  3.0983375196473186e-10_BK, 2.4255404477032373e-11_BK,-9.8343678688258e-12_BK, &
            6.661105552981143e-13_BK,  1.1249812587695572e-13_BK,-2.4657588515917296e-14_BK,7.965965484631788e-16_BK, &
            3.0824314548929273e-16_BK,-4.420019468170955e-17_BK,  1.167553319557466e-19_BK, 5.863076185446756e-19_BK, &
           -5.882695924413488e-20_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_ai_neg_a1

    pure real(BK) function evalpoly_taylor_ai_neg_a2(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
           -2.720348378642871e-16_BK,-0.803111369654864_BK,       5.560323321157856e-16_BK, 0.5471797795259772_BK, &
           -0.06692594747123885_BK,  -0.11184216377764623_BK,     0.027358988976298886_BK,  0.009292361042237976_BK, &
           -0.003994362992058797_BK, -0.00014760712751364074_BK,  0.00028467905572535616_BK,-3.082684106536036e-5_BK, &
           -9.934550872135151e-6_BK,  2.6326770738641654e-6_BK,   5.3764289286128445e-8_BK,-9.855619834673579e-8_BK, &
            1.0053714072345167e-8_BK, 1.6788861968137084e-9_BK,  -4.5638978170010107e-10_BK, 9.328982974624003e-12_BK, &
            9.327854082162345e-12_BK,-1.1774433153941012e-12_BK, -6.234375094261117e-14_BK,  2.7947001637990885e-14_BK, &
           -1.6713500242449534e-15_BK,-2.9431613458960554e-16_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_ai_neg_a2

    pure real(BK) function evalpoly_taylor_ai_neg_a3(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
            2.313678943005095e-16_BK, 0.8652040258941519_BK,    -6.386401513932251e-16_BK,-0.7960684314096329_BK, &
            0.07210033549117963_BK,   0.21973717014275287_BK,    -0.0398034215704817_BK,  -0.027165996636626166_BK, &
            0.007847756076526893_BK,  0.0015301123354424383_BK, -0.0007832222619265922_BK,-5.448369227182677e-6_BK, &
            4.434801281139783e-5_BK,  -4.82784752334856e-6_BK,   -1.3751331165365516e-6_BK, 3.3809730430936423e-7_BK, &
            1.1515237992029015e-8_BK,-1.1917718796669944e-8_BK,   8.97146222351664e-10_BK,  2.2604595796334607e-10_BK, &
           -4.439596892555847e-11_BK,-8.351286011527614e-13_BK,   1.0197761050717805e-12_BK,-7.862765122280811e-14_BK, &
           -1.1711709421130035e-14_BK, 2.423074596316538e-15_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_ai_neg_a3

    pure real(BK) function evalpoly_taylor_ai_neg_a4(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
           -8.710477837103708e-17_BK,-0.9108507370496018_BK,      2.9557735202731247e-16_BK, 1.0302796776637262_BK, &
           -0.07590422808746698_BK,  -0.3496103711718467_BK,      0.05151398388318634_BK,   0.05468569776946418_BK, &
           -0.012486084684708815_BK, -0.004439192827500768_BK,    0.0015491678856944288_BK, 0.00016037655628253796_BK, &
           -0.00011327987159259477_BK, 2.9534552161146124e-6_BK,  5.10535152341918e-6_BK,  -6.348767142926879e-7_BK, &
           -1.3206281362722718e-7_BK,  3.461056785480714e-8_BK,   8.542321939505024e-10_BK,-1.0729667675129367e-9_BK, &
            7.582406135085136e-11_BK,  1.9371772465258795e-11_BK,-3.436282550010513e-12_BK,-1.0997332719001945e-13_BK, &
            7.734206349128181e-14_BK, -4.483209467796499e-15_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_ai_neg_a4

    pure real(BK) function evalpoly_taylor_ai_neg_a5(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
           -3.222967925030853e-17_BK, 0.9473357094415678_BK,     1.28018438717254e-16_BK, -1.25429357127562_BK, &
            0.0789446424534639_BK,    0.49821378438402086_BK,    -0.06271467856378098_BK,  -0.09235552894376721_BK, &
            0.017793349442286454_BK,  0.00931902751214878_BK,   -0.0025967585986179237_BK, -0.0005112568183297996_BK, &
            0.00022687897509904778_BK,9.389319637951164e-6_BK,  -1.2712163212229125e-5_BK,  7.25185550490412e-7_BK, &
            4.5990184323826115e-7_BK,-6.791593419402475e-8_BK,  -9.56972591227924e-9_BK,    2.922324845525447e-9_BK, &
            2.1334860119088188e-11_BK,-7.805967826213263e-11_BK,  5.9585235209456126e-12_BK, 1.2676904585002835e-12_BK, &
           -2.2716482806979397e-13_BK,-6.8536313808318734e-15_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_ai_neg_a5

    pure real(BK) function evalpoly_taylor_ai_neg_a6(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
            2.1834671977219237e-16_BK,-0.9779228085694986_BK,   -9.850331087383878e-16_BK, 1.4705760105401993_BK, &
           -0.0814935673807908_BK,    -0.6634246948201652_BK,    0.07352880052700973_BK,   0.14057990051109173_BK, &
           -0.02369373910072015_BK,   -0.016595479983083087_BK,  0.00393733595363381_BK,   0.0011458316593658873_BK, &
           -0.0003948536938259645_BK, -4.103271183031362e-5_BK,  2.587065207093167e-5_BK, -1.1728505435852704e-7_BK, &
           -1.1435607200608036e-6_BK,  9.900321384824919e-8_BK,  3.3335503439036955e-8_BK,-5.955649567170201e-9_BK, &
           -5.309773544803389e-10_BK,  2.0731250021063097e-10_BK,-2.5212690187520218e-12_BK,-4.7460190937036376e-12_BK, &
            4.1677722875756535e-13_BK, 6.716734034504283e-14_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_ai_neg_a6

    pure real(BK) function evalpoly_taylor_aip_neg_b1(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
           -1.1246873724687218e-17_BK,-0.5457232363649821_BK,    0.26782832800784995_BK,   0.09266316627889255_BK, &
           -0.09095387272749701_BK,    0.013134992740413161_BK,  0.006949737470916942_BK, -0.002917297275227855_BK, &
            0.00014721097333898997_BK, 0.0001515927648587197_BK,-3.813263266242759e-5_BK,  8.296458876147209e-8_BK, &
            1.55758560220824e-6_BK,   -2.672035857158013e-7_BK, -8.22515911965634e-9_BK,   9.283928640130678e-9_BK, &
           -1.1579575909261867e-9_BK, -6.7029059958855e-11_BK,   3.609115436353033e-11_BK,-3.38533335716781e-12_BK, &
           -2.829534500722224e-13_BK,  9.866581579686296e-14_BK,-7.069976274435534e-15_BK,-7.844811532823979e-16_BK, &
            1.9991570243463872e-16_BK,-1.0963571103138449e-17_BK,-1.5705239700998077e-18_BK, 3.120825352080613e-19_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_aip_neg_b1

    pure real(BK) function evalpoly_taylor_aip_neg_b2(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
           -3.315151654629452e-17_BK, 1.3610450626413026_BK,    -0.20950773901628192_BK,  -0.7368238802848806_BK, &
            0.22684084377355043_BK,   0.10570029472060026_BK,   -0.055261791021366045_BK, -0.0016934864099735572_BK, &
            0.005407468330176707_BK, -0.0008007715851877276_BK, -0.00021632997518904743_BK,7.826693155932615e-5_BK, &
           -1.3497442526091854e-6_BK,-3.142453810419922e-6_BK,   4.899638015045217e-7_BK,  4.168447683589849e-8_BK, &
           -2.0660058467837985e-8_BK, 1.4236346987697045e-9_BK,  3.6404479934989394e-10_BK,-7.748421874625937e-11_BK, &
            8.427265218804001e-13_BK, 1.511640723657024e-12_BK, -1.8202548037745256e-13_BK,-7.958994836861e-15_BK, &
            3.934069548729346e-15_BK,-2.7347873985422944e-16_BK,-3.2414238181600296e-17_BK, 7.093654290955294e-18_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_aip_neg_b2

    pure real(BK) function evalpoly_taylor_aip_neg_b3(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
            1.9002007513487082e-16_BK,-1.8335969193618502_BK,   0.19020323431407618_BK,   1.473019844105969_BK, &
           -0.30559948656030816_BK,  -0.34232487381035187_BK,   0.11047648830794765_BK,   0.030555249293994464_BK, &
           -0.016640833099863316_BK, -0.00029195123942834367_BK,0.001273168021708355_BK, -0.00015529618742302247_BK, &
           -4.892380712504251e-5_BK,  1.370162563867412e-5_BK,  3.713190536062472e-7_BK, -5.653827230320543e-7_BK, &
            5.3710493778794144e-8_BK, 1.1475273803238557e-8_BK, -2.8091799082246527e-9_BK, 4.555486497839958e-12_BK, &
            6.750872335321682e-11_BK,-7.0928318552774136e-12_BK,-6.939729077051773e-13_BK, 2.073351337864048e-13_BK, &
           -7.373577967817808e-15_BK,-2.872535929393802e-15_BK,  3.869468846775618e-16_BK, 8.799696696329233e-18_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_aip_neg_b3

    pure real(BK) function evalpoly_taylor_aip_neg_b4(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
            3.3524788679156663e-16_BK, 2.205896662123774_BK,    -0.17895397185614714_BK, -2.2659365205746744_BK, &
            0.3676494436872962_BK,     0.6863528964430919_BK,   -0.1699452390431007_BK,  -0.09021486959017666_BK, &
            0.0330030319232237_BK,     0.005024983874975287_BK, -0.0033877728633306464_BK, 5.181377628766877e-5_BK, &
            0.00020005581568389722_BK,-2.573779771192161e-5_BK,  -6.466340051454468e-6_BK, 1.7813080120881245e-6_BK, &
            5.115762237080756e-8_BK,  -6.572124324183966e-8_BK,  5.154704616302637e-9_BK,  1.3427693455239035e-9_BK, &
           -2.661643481572335e-10_BK, -6.78546505935419e-12_BK,  6.6025118995813095e-12_BK,-4.684148897527137e-13_BK, &
           -8.657102435817311e-14_BK,  1.6294270949749456e-14_BK, 7.020254520958461e-17_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_aip_neg_b4

    pure real(BK) function evalpoly_taylor_aip_neg_b5(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
            1.0028785856694853e-16_BK,-2.523505448425921_BK,    0.17115062220581154_BK,  3.1006215783124498_BK, &
           -0.42058424140431994_BK,  -1.1315065523268806_BK,    0.23254661837343366_BK,  0.18659442331706105_BK, &
           -0.05418689050759992_BK,  -0.01541443888390889_BK,   0.006771045421706508_BK, 0.00048573017836110537_BK, &
           -0.000506615377447726_BK,  2.439555149034597e-5_BK,  2.3412450856604934e-5_BK,-3.454448194852507e-6_BK, &
           -6.10261028364471e-7_BK,   1.854414901637662e-7_BK,  2.707847274247799e-9_BK, -5.88674348492048e-9_BK, &
            4.6258175023029366e-10_BK,1.1011543328016952e-10_BK,-2.076042276366767e-11_BK,-6.466028459201257e-13_BK, &
            4.858155738241608e-13_BK,-2.8160298411012274e-14_BK,-6.546251007171164e-15_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_aip_neg_b5

    pure real(BK) function evalpoly_taylor_aip_neg_b6(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
            1.9865763636596917e-15_BK, 2.8052430870313785_BK,   -0.16523811457399187_BK,-3.968711454994397_BK, &
            0.4675405145052357_BK,     1.6734018525386696_BK,   -0.29765335912458146_BK,-0.3248476382988628_BK, &
            0.07998087056358093_BK,    0.033573466676480326_BK, -0.011604112798228065_BK,-0.001782911750352928_BK, &
            0.0010260012946360286_BK,  1.586660639280247e-5_BK, -5.84653082865604e-5_BK, 4.620195206034647e-6_BK, &
            2.1386747050487544e-6_BK, -3.734612544877612e-7_BK, -4.32848184829349e-8_BK, 1.5890643286120207e-8_BK, &
           -7.049077495383936e-11_BK, -4.296440036986022e-10_BK, 3.7410249910157256e-11_BK,7.061620449058568e-12_BK, &
           -1.389003055566532e-12_BK, -3.4842814979282226e-14_BK,2.945598253096211e-14_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_aip_neg_b6

    pure real(BK) function evalpoly_taylor_aip_neg_b7(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
           -1.0626912677608902e-15_BK,-3.0610916737763576_BK,    0.16051114409736328_BK,  4.864813950020498_BK, &
           -0.510181945629397_BK,     -2.3087085355595987_BK,    0.36486104625153865_BK,  0.5095798638807557_BK, &
           -0.11022512873634545_BK,   -0.061695677462806485_BK,  0.01804803829143557_BK,  0.004234760217630525_BK, &
           -0.001817889354798794_BK,  -0.00013263827730311434_BK,0.00012045081443999045_BK,-3.2998164849963075e-6_BK, &
           -5.377771018267827e-6_BK,   5.880371349671046e-7_BK,  1.5612224361123518e-7_BK,-3.304476207996624e-8_BK, &
           -2.2841837768922437e-9_BK,  1.1415139263443008e-9_BK, -2.7957320677518033e-11_BK,-2.624071604789366e-11_BK, &
            2.644903133984302e-12_BK,  3.6840692435819356e-13_BK,-8.085295122966928e-14_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_aip_neg_b7

    ! Bi Taylor series
    pure real(BK) function evalpoly_taylor_bi_neg_c1(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
           -6.5985241890579e-17_BK,   0.6019578879762396_BK,    3.8723875455316416e-17_BK,-0.11775432210529542_BK, &
            0.05016315733135329_BK,   0.006910490244306752_BK, -0.00588771610526477_BK,    0.001001243418004953_BK, &
            0.0002468032230109554_BK,-9.809567700177314e-5_BK,  7.906302352775167e-6_BK,   3.2903583290776384e-6_BK, &
           -8.134500652724367e-7_BK,  2.5925418426189934e-8_BK, 2.3324810036821433e-8_BK, -4.0184717699215364e-9_BK, &
           -6.046748046691523e-12_BK, 1.030931745891974e-10_BK,-1.3109067391453412e-11_BK,-3.7148704746518476e-13_BK, &
            3.117880534840298e-13_BK,-3.017392412328179e-14_BK,-1.5961835682297216e-15_BK, 6.86173097250094e-16_BK, &
           -5.126895355684697e-17_BK,-4.002590009232389e-18_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_bi_neg_c1

    pure real(BK) function evalpoly_taylor_bi_neg_c2(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
            1.4912020792246877e-16_BK,-0.7603101414928011_BK,  -2.43893056726376e-16_BK,  0.4145075686526103_BK, &
           -0.06335917845773335_BK,   -0.06779464658972667_BK,  0.020725378432630507_BK,  0.0037715103802399754_BK, &
           -0.0024212373782045234_BK,  0.00011650577897569078_BK,0.00012990670836290934_BK,-2.5475805923214892e-5_BK, &
           -2.3365998844222516e-6_BK,  1.366925939119848e-6_BK, -9.79811521421597e-8_BK, -3.241877223500228e-8_BK, &
            7.030964290401067e-9_BK,   2.9645869860049095e-11_BK,-1.811036354176953e-10_BK,2.02748242223541e-11_BK, &
            1.6369809441952335e-12_BK,-5.891058979705288e-13_BK,3.229460350416958e-14_BK,  7.043480833415071e-15_BK, &
           -1.2585952159615458e-15_BK, 1.5424534202215543e-17_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_bi_neg_c2

    pure real(BK) function evalpoly_taylor_bi_neg_c3(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
            1.9733757157816445e-16_BK, 0.8369910126192611_BK,   -4.766430373021629e-16_BK,-0.6738806929651456_BK, &
            0.06974925105160529_BK,    0.16276704821360755_BK,  -0.033694034648257314_BK,-0.017060373526892724_BK, &
            0.005813108864771701_BK,   0.0006766688519582255_BK,-0.0005015775388627367_BK, 2.3129991230391897e-5_BK, &
            2.3482260980834454e-5_BK, -3.931490146004496e-6_BK, -4.961904148214753e-7_BK,  2.0225837715837952e-7_BK, &
           -6.393851385235093e-9_BK,  -5.4163515115941455e-9_BK, 7.61913062088228e-10_BK,  5.781030065954972e-11_BK, &
           -2.3939352033369288e-11_BK, 1.1491587025101967e-12_BK, 3.754437974718945e-13_BK,-5.828189024438837e-14_BK, &
           -1.2038258285408947e-15_BK, 1.0949805502650894e-15_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_bi_neg_c3

    pure real(BK) function evalpoly_taylor_bi_neg_c4(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
           -2.4652322944876206e-16_BK,-0.8894799014265397_BK,    7.605059359461805e-16_BK, 0.914659910484288_BK, &
           -0.07412332511887869_BK,   -0.28216581976907734_BK,   0.045732995524214476_BK,  0.03968566805999452_BK, &
           -0.010077350706038485_BK,  -0.002765579278043845_BK,  0.0011317936862374988_BK, 6.350785898737327e-5_BK, &
           -7.385287092336652e-5_BK,   4.7433306878194e-6_BK,    2.8525777573391674e-6_BK,-4.910405707756589e-7_BK, &
           -5.3569384414458937e-8_BK,  2.162582892656635e-8_BK, -5.245927788687565e-10_BK,-5.467764650394962e-10_BK, &
            6.54276021046472e-11_BK,   6.7831837087265114e-12_BK,-2.057262110682013e-12_BK, 4.6593600123797446e-14_BK, &
            3.5282946957989075e-14_BK,-3.907896222619103e-15_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_bi_neg_c4

    pure real(BK) function evalpoly_taylor_bi_neg_c5(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
            3.9229436221051404e-16_BK, 0.9299836385680267_BK,  -1.446931087552141e-15_BK,-1.143378006570179_BK, &
            0.07749863654733644_BK,    0.42172137606250004_BK, -0.05716890032850917_BK,  -0.07222475282022156_BK, &
            0.015061477716517889_BK,   0.006605776631491822_BK,-0.00203699656109637_BK,  -0.0003060705894858609_BK, &
            0.00016388042135750047_BK, 1.415470237388889e-6_BK,-8.32405201816858e-6_BK,   7.306611151707449e-7_BK, &
            2.6175008963208153e-7_BK, -5.0418989799370194e-8_BK,-3.922245164373553e-9_BK, 1.8528625194856889e-9_BK, &
           -5.65408431706205e-11_BK,  -4.188183603904142e-11_BK,4.913313565632472e-12_BK, 4.988369500720917e-13_BK, &
           -1.415329370181925e-13_BK,  2.055853447588659e-15_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_bi_neg_c5

    pure real(BK) function evalpoly_taylor_bi_neg_c6(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
           -2.5434712503223577e-16_BK,-0.9632344301904238_BK,   1.0799513875152367e-15_BK, 1.3632895847289495_BK, &
           -0.08026953584920274_BK,   -0.5788492708248631_BK,   0.06816447923644768_BK,    0.11512592528178081_BK, &
           -0.020673188243745148_BK,  -0.01263165262727079_BK,  0.0032297953593547093_BK,  0.0007872196337628573_BK, &
           -0.0003034765879060938_BK, -2.214893270209827e-5_BK, 1.848531480425613e-5_BK,  -5.494713542623637e-7_BK, &
           -7.463553351311001e-7_BK,   8.511543175468515e-8_BK, 1.891679729735395e-8_BK,  -4.2957638235100794e-9_BK, &
           -1.987501137742938e-10_BK,  1.3189572367449318e-10_BK,-5.644995723162025e-12_BK,-2.606327771326931e-12_BK, &
            3.257839467288938e-13_BK,  2.747967731363848e-14_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_bi_neg_c6

    pure real(BK) function evalpoly_taylor_bi_neg_c7(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
            1.8020634719987473e-16_BK, 0.9915863705176604_BK,   -8.59421583992181e-16_BK,-1.576323924317981_BK, &
            0.08263219754313905_BK,    0.7517641997479385_BK,   -0.07881619621589926_BK, -0.16875811588300885_BK, &
            0.02684872141956927_BK,    0.02126154883567329_BK,  -0.004720515995775305_BK,-0.0015995278561913206_BK, &
            0.000502172333137813_BK,   6.753904879331083e-5_BK, -3.5106292206536236e-5_BK,-6.763249640425508e-7_BK, &
            1.676623699155756e-6_BK,  -1.0535063689666924e-7_BK,-5.447152845623752e-8_BK,  7.840580561006179e-9_BK, &
            1.0900247096289847e-9_BK, -3.0775359498580484e-10_BK,-5.5330887396656266e-12_BK,7.95541169870788e-12_BK, &
           -4.619165200519967e-13_BK, -1.3568891981611074e-13_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_bi_neg_c7

    ! Bi' Taylor series
    pure real(BK) function evalpoly_taylor_bip_neg_d1(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
            1.3050830085951108e-16_BK, 1.0438424472052532_BK,   -0.22747219181982883_BK, -0.39917225554412844_BK, &
            0.17397374120087558_BK,    0.030629020377963098_BK, -0.029937919165809637_BK, 0.0032974297534822746_BK, &
            0.0018647251224384134_BK, -0.0005802849783048269_BK,-6.3210093935846245e-6_BK,3.093950621960744e-5_BK, &
           -4.725835616276457e-6_BK,  -4.992594707816475e-7_BK, 2.4374144602236595e-7_BK,-1.8780193767765568e-8_BK, &
           -4.559045448346791e-9_BK,   1.1142680058276989e-9_BK, -3.1024513355303083e-11_BK,-2.1590191000016282e-11_BK, &
            3.2825148908037815e-12_BK, 4.0190496443783935e-14_BK,-6.537063414339622e-14_BK,6.613854651768023e-15_BK, &
            3.4783752581066423e-16_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_bip_neg_d1

    pure real(BK) function evalpoly_taylor_bip_neg_d2(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
           -3.5706417820884903e-16_BK,-1.615099007771363_BK,    0.19826141804723305_BK,  1.0964247904764644_BK, &
           -0.26918316796189407_BK,   -0.21007798288587065_BK,  0.08223185928573487_BK,  0.012682390560458464_BK, &
           -0.010357751717268665_BK,   0.0005878045014194119_BK,0.0006272935388783455_BK,-0.000126389381050211_BK, &
           -1.4458173626233763e-5_BK,  7.68669086448864e-6_BK, -4.287442087180428e-7_BK,-2.2323535638602655e-7_BK, &
            4.1592007781675354e-8_BK,  1.661562764553658e-9_BK,-1.3287524977803896e-9_BK, 1.0897892878497395e-10_BK, &
            1.885812317231568e-11_BK, -4.387083138454629e-12_BK, 8.141951267221141e-14_BK, 7.435849614474649e-14_BK, &
           -8.909656360099548e-15_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_bip_neg_d2

    pure real(BK) function evalpoly_taylor_bip_neg_d3(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
            4.2827700216989146e-17_BK, 2.0283916344286097_BK,   -0.18398458074347976_BK, -1.8635495639516062_BK, &
            0.3380652724047683_BK,     0.5013654942009328_BK,   -0.13976621729637048_BK,-0.05614396862733763_BK, &
            0.02420309123708635_BK,    0.00207992990472809_BK,  -0.002184210905396979_BK, 0.0001402447897449079_BK, &
            0.00010854664971815054_BK,-2.0229873802934234e-5_BK,-2.4528584220804398e-6_BK,1.0876736284147158e-6_BK, &
           -3.397391040630207e-8_BK,  -3.166202125524556e-8_BK,  4.388662970046761e-9_BK,4.0514975156612263e-10_BK, &
           -1.516133406608186e-10_BK,  5.6816654076685864e-12_BK,2.729783486500342e-12_BK,-3.7579567525699357e-13_BK, &
           -1.6499499832774852e-14_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_bip_neg_d3

    pure real(BK) function evalpoly_taylor_bip_neg_d4(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
           -7.721317833829156e-16_BK,-2.370056419850032_BK,     0.17474955841590484_BK,  2.6786750727687814_BK, &
           -0.39500940330834017_BK,  -0.8965942491179204_BK,    0.20090063045765894_BK,  0.13347757899714746_BK, &
           -0.04300701702659608_BK,  -0.00938263950825763_BK,   0.00490895024525397_BK,  0.0001440078792007367_BK, &
           -0.0003303783370671734_BK, 2.80683244741357e-5_BK,   1.3167040208955681e-5_BK,-2.600626800799143e-6_BK, &
           -2.4673488268918906e-7_BK, 1.1647227562465254e-7_BK, -3.5620394975059554e-9_BK,-3.07333770503881e-9_BK, &
            3.871005164858819e-10_BK, 4.0694506307856404e-11_BK,-1.2666768772243378e-11_BK,2.5607202530829446e-13_BK, &
            2.3268359659136705e-13_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_bip_neg_d4

    pure real(BK) function evalpoly_taylor_bip_neg_d5(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
            8.51891254034672e-16_BK,   2.6681083909107306_BK,  -0.16801312006683342_BK,-3.5308762309835418_BK, &
            0.44468473181845736_BK,    1.3905885354962473_BK,  -0.26481571732376624_BK,-0.2501880899230309_BK, &
            0.06651852511898834_BK,    0.02338738341397691_BK, -0.008995895297201397_BK,-0.0010162775547694785_BK, &
            0.0007360237718500377_BK,-1.1181272935399847e-5_BK,-3.816004186244733e-5_BK, 4.197249000708146e-6_BK, &
            1.2125733995252362e-6_BK, -2.7217261662439e-7_BK,  -1.6890432688596642e-8_BK, 1.0073100751030303e-8_BK, &
           -4.03105957534568e-10_BK,  -2.3276576844590606e-10_BK,2.982140502613569e-11_BK, 2.817984811976686e-12_BK, &
           -8.698067306477899e-13_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_bip_neg_d5

    pure real(BK) function evalpoly_taylor_bip_neg_d6(t) result(y)
        real(BK), intent(in) :: t
        real(BK), parameter :: c(*) = [ &
            1.4408236312418457e-15_BK,-2.935962201853191_BK,    0.1627548682163053_BK,    4.413525969647323_BK, &
           -0.4893270336421936_BK,   -1.9795579449241065_BK,    0.3310144477235477_BK,    0.411133225177042_BK, &
           -0.09455529773724823_BK,  -0.04624928014006227_BK,   0.014615269647565455_BK,  0.0028371617779749344_BK, &
           -0.0013840746291783367_BK,-6.183388868085215e-5_BK,  8.548004625717443e-5_BK, -4.442028547574055e-6_BK, &
           -3.4885208239390457e-6_BK, 4.82514569024414e-7_BK,   8.740310143336716e-8_BK, -2.3525755287654844e-8_BK, &
           -7.34259540066286e-10_BK,  7.242756561007131e-10_BK,-3.9132744422986586e-11_BK,-1.4430610543508328e-11_BK, &
            2.0111563898805826e-12_BK, 1.4887321103652866e-13_BK,-5.1033351829108886e-14_BK]
        y = evalpoly(size(c), t, c)
    end function evalpoly_taylor_bip_neg_d6

end module bessels_airy
