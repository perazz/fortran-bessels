!  ************************************************************************************************************
!                                    ____  ________________ ________   _____
!                                   / __ )/ ____/ ___/ ___// ____/ /  / ___/
!                                  / __  / __/  \__ \\__ \/ __/ / /   \__ \
!                                 / /_/ / /___ ___/ /__/ / /___/ /______/ /
!                                /_____/_____//____/____/_____/_____/____/
!
!                          Documentation helper: dump "x,f(x)" CSV for one function
!
!  Prints a dense grid of x,f(x) pairs to stdout so doc/generate_plots.py can compare the
!  library's output against an mpmath arbitrary-precision reference and plot the relative error.
!  Test-time tooling only — not part of the library.
!
!  Usage:  fpm run --example dump_values -- <name> [order]
!          <name> in {besselj0,besselj1,besseljn,bessely0,bessely1,
!                     besseli0,besseli1,besselk0,besselk1,gamma}
!
!  MIT License
!  Copyright (c) 2022 Federico Perini
!  Copyright (c) 2021-2022 Michael Helton, Oscar Smith, and the Bessels.jl contributors
!  ************************************************************************************************************
program dump_values
    use bessels
    implicit none

    integer, parameter :: NPTS = 500
    character(len=64) :: name, arg
    integer  :: i, nu, ios
    real(BK) :: x0, x1, x, dx, f

    if (command_argument_count() < 1) then
        write(*,'(a)') '# usage: dump_values <name> [order]'
        stop 1
    end if

    call get_command_argument(1, name)

    ! optional integer order (for besseljn); default 5
    nu = 5
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read(arg, *, iostat=ios) nu
    end if

    ! function-appropriate domain
    select case (trim(name))
       case ('besselj0','besselj1','besseljn'); x0 =  0.0_BK; x1 = 30.0_BK
       case ('bessely0','bessely1');            x0 =  0.1_BK; x1 = 30.0_BK
       case ('besseli0','besseli1');            x0 =  0.0_BK; x1 = 15.0_BK
       case ('besselk0','besselk1');            x0 =  0.1_BK; x1 = 15.0_BK
       case ('gamma');                          x0 =  0.1_BK; x1 =  8.0_BK
       case default
          write(*,'(a)') '# unknown function: '//trim(name)
          stop 1
    end select

    dx = (x1 - x0) / real(NPTS - 1, BK)

    do i = 1, NPTS
        x = x0 + real(i - 1, BK)*dx
        select case (trim(name))
           case ('besselj0'); f = besselj0(x)
           case ('besselj1'); f = besselj1(x)
           case ('besseljn'); f = besseljn(nu, x)
           case ('bessely0'); f = bessely0(x)
           case ('bessely1'); f = bessely1(x)
           case ('besseli0'); f = besseli0(x)
           case ('besseli1'); f = besseli1(x)
           case ('besselk0'); f = besselk0(x)
           case ('besselk1'); f = besselk1(x)
           case ('gamma');    f = gamma_BK(x)
        end select
        write(*,'(es24.16e3,",",es24.16e3)') x, f
    end do

end program dump_values
