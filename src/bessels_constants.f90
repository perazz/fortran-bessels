module bessels_constants
    use iso_fortran_env
    implicit none
    public

    integer, parameter :: BKIND = real64


    real(BKIND), parameter :: ZERO     = 0.00_BKIND
    real(BKIND), parameter :: ONE      = 1.00_BKIND
    real(BKIND), parameter :: TWO      = 2.00_BKIND
    real(BKIND), parameter :: FOUR     = 4.00_BKIND
    real(BKIND), parameter :: HALF     = 0.50_BKIND
    real(BKIND), parameter :: FOURTH   = 0.25_BKIND
    real(BKIND), parameter :: TWOTHD   = 2.0_BKIND/3.0_BKIND
    real(BKIND), parameter :: THIRD    = 1.0_BKIND/3.0_BKIND
    real(BKIND), parameter :: SIXTH    = 1.0_BKIND/6.0_BKIND


    real(BKIND), parameter :: PI       = acos(-1.0_BKIND)
    real(BKIND), parameter :: ONEOSQPI = ONE/SQRT(PI)
    real(BKIND), parameter :: TWOOPI   = TWO/PI
    real(BKIND), parameter :: PIO2     = PI*HALF
    real(BKIND), parameter :: PIO4     = PI*FOURTH
    real(BKIND), parameter :: SQPIO2   = 1.253314137315500251207882642405522626503493370304969158314961788171146827303924_BKIND
    real(BKIND), parameter :: SQ1O2PI  = 0.3989422804014326779399460599343818684758586311649346576659258296706579258993008_BKIND
    real(BKIND), parameter :: SQ2OPI   = 0.7978845608028653558798921198687637369517172623298693153318516593413158517986017_BKIND
    real(BKIND), parameter :: SQ2O2    = 0.707106781186547524400844362104849039284835937688474036588339868995366239231051_BKIND

    real(BKIND), parameter :: THPIO4   = 2.35619449019234492885_BKIND
    real(BKIND), parameter :: SQ2PI    = 2.5066282746310007_BKIND


    real(BKIND), parameter :: GAMMA_TWO_THIRDS  = gamma(TWOTHD)
    real(BKIND), parameter :: GAMMA_ONE_THIRD   = gamma(THIRD)
    real(BKIND), parameter :: GAMMA_ONE_SIXTH   = gamma(SIXTH)
    real(BKIND), parameter :: GAMMA_FIVE_SIXTHS = gamma(5.0_BKIND*SIXTH)

    contains

      ! cubic root function in double precision
      elemental real(BKIND) function cbrt(x)
         real(BKIND), intent(in) :: x
         real(BKIND), parameter :: c(0:23) = &
              [real(BKIND) ::  1.5319394088521e-3, -1.8843445653409e-2, &
                               1.0170534986000e-1, -3.1702448761286e-1, &
                               6.3520892642253e-1, -8.8106985991189e-1, &
                               1.0517503764540d0, 4.2674123235580e-1, &
                               1.5079083659190e-5, -3.7095709111375e-4, &
                               4.0043972242353e-3, -2.4964114079723e-2, &
                               1.0003913718511e-1, -2.7751961573273e-1, &
                               6.6256121926465e-1, 5.3766026150315e-1, &
                               1.4842542902609e-7, -7.3027601203435e-6, &
                               1.5766326109233e-4, -1.9658008013138e-3, &
                               1.5755176844105e-2, -8.7413201405100e-2, &
                               4.1738741349777e-1, 6.7740948115980e-1 ]


         real(BKIND), parameter :: TWO_POW_3   = 8.0_BKIND
         real(BKIND), parameter :: TWO_POW_16  = 65536.0_BKIND
         real(BKIND), parameter :: TWO_POW_48  = 281474976710656.0_BKIND
         real(BKIND), parameter :: TWO_POW_M3  = ONE/TWO_POW_3
         real(BKIND), parameter :: TWO_POW_M16 = ONE/TWO_POW_16
         real(BKIND), parameter :: TWO_POW_M48 = ONE/TWO_POW_48

         integer :: k
         real(BKIND) :: w,y,u

         if (x==ZERO) then
            cbrt = ZERO
            return
         end if

         w = abs(x)
         y = sign(HALF,x)

         if (w > TWO_POW_3) then
             do while (w > TWO_POW_48)
                 w = w * TWO_POW_M48
                 y = y * TWO_POW_16
             end do
             do while (w > TWO_POW_3)
                 w = w * TWO_POW_M3
                 y = y * TWO
             end do
         else if (w < ONE) then
             do while (w < TWO_POW_M48)
                 w = w * TWO_POW_48
                 y = y * TWO_POW_M16
             end do
             do while (w < ONE)
                 w = w * TWO_POW_3
                 y = y * HALF
             end do
         end if

         k = merge(0,merge(8,16,w<FOUR),w<TWO)

         u = ((((((c(k)   *w + c(k+1))*w +  c(k+2))*w+ c(k+3))*w+ &
                   c(k+4))*w + c(k+5))*w +  c(k+6))*w+ c(k+7)

         cbrt = y*(u+3*u*w/(w+2*u*u*u))
      end function cbrt


end module bessels_constants
