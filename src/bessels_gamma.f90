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
    private

    public :: gamma_BK

    interface gamma_BK
        module procedure gamma_BK
        module procedure gamma_integer
    end interface

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

            if (x<ZERO) then
                gamma_BK = PI*q/(s*z*p)
            else
                gamma_BK = z*p/q
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

   ! real64 version adapted from Cephes Mathematical Library (MIT license
   ! https://en.smath.com/view/CephesMathLibrary/license) by Stephen L. Moshier
   elemental real(BK) function gamma_integer(x) result(gam)
      integer, intent(in) :: x

      integer :: i
      integer, parameter  :: SMALL = 16
      real(BK), parameter :: rsmall   (SMALL) = [(real(i,BK),i=1,SMALL)]
      real(BK), parameter :: factorial(SMALL) = [(product(rsmall(1:i)),i=1,SMALL)]

      if (x<0) then

         gam = ieee_value(gam,ieee_quiet_nan)

      elseif (x==0) then

         gam = huge(gam)

      elseif (x<SMALL) then

         gam = factorial(x+1)

      else

         gam = gamma_BK(real(x,BK))

      end if

   end function gamma_integer

end module bessels_gamma
