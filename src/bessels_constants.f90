module bessels_constants
    use iso_fortran_env
    implicit none
    private

    integer, parameter :: BKIND = real64

    real(BKIND), parameter :: ONE      = 1.00_BKIND
    real(BKIND), parameter :: TWO      = 2.00_BKIND
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



end module bessels_constants
