SUBROUTINE subtract_poly ( locwvl, npts, ll_rad, lu_rad )

  USE OMSAO_precision_module
  USE OMSAO_indices_module,     ONLY : max_rs_idx, solar_idx
  USE OMSAO_parameters_module,  ONLY : max_spec_pts,  elsunc_np, elsunc_nw
  USE OMSAO_variables_module,   ONLY : poly_x, poly_y, poly_w, database, poly_order
  USE bounded_nonlin_LS,        ONLY : elsunc

  IMPLICIT NONE

  EXTERNAL poly_specfit

  INTEGER,                          INTENT (IN) :: npts, ll_rad, lu_rad
  REAL (KIND=dp), DIMENSION (npts), INTENT (IN) :: locwvl

  INTEGER                          :: i, nlower, nupper, nfitted
  REAL (KIND=dp)                   :: locavg, chisq
  REAL (KIND=dp), DIMENSION (npts) :: x, ptemp, sig

  ! ================
  ! ELSUNC variables
  ! ================
  INTEGER                                      :: exval
  REAL (KIND=dp), DIMENSION (max_spec_pts)     :: fitres
  REAL (KIND=dp)                               :: rms, chisq2
  INTEGER                                      :: elbnd
  INTEGER,        DIMENSION (elsunc_np)        :: p
  REAL (KIND=dp), DIMENSION (elsunc_nw)        :: w
  REAL (KIND=dp), DIMENSION (poly_order)        :: blow, bupp
  REAL (KIND=dp), DIMENSION (max_spec_pts)           :: f
  REAL (KIND=dp), DIMENSION (max_spec_pts,poly_order) :: dfda
  REAL (KIND=dp), DIMENSION (poly_order)        :: par


  ! ======================
  ! Assign fitting weights
  ! ======================
  sig = 1.0_dp

  ! Find limits for polynomial fitting, with ~1 nm overlap
  nlower = MINVAL ( MINLOC ( locwvl(1:npts), MASK=(locwvl(1:npts) >= locwvl(ll_rad)-1.0_dp) ) )
  nupper = MAXVAL ( MAXLOC ( locwvl(1:npts), MASK=(locwvl(1:npts) <= locwvl(lu_rad)+1.0_dp) ) )
  nfitted = nupper - nlower + 1

  ! Find average position over fitted region
  locavg = SUM ( locwvl(1+nlower-1:nfitted+nlower-1) ) / REAL ( nfitted, KIND=dp )

  ! Load temporary position file: re-define positions in order to fit
  ! about mean position
  DO i = 1, nfitted
     ptemp(i) = locwvl(i+nlower-1) - locavg
  END DO

  !     Load and fit database spectra nos. 2-11
  DO i = 1, max_rs_idx
     IF ( i /= solar_idx ) THEN
        ! ===============================================================
        ! ELBND: 0 = unconstrained
        !        1 = all variables have same lower bound
        !        else: lower and upper bounds must be supplied by the use
        ! ===============================================================  
        elbnd = 0  ;  exval = 0
        p   = -1 ; p(1) = 0  ;  p(3) = 5  ; w = -1.0
        blow(1:poly_order) = 0.0_dp  ;  bupp(1:poly_order) = 0.0_dp
  
        poly_x(1:nfitted) = ptemp(1:nfitted)
        poly_y(1:nfitted) = database(i, 1+nlower-1:nfitted+nlower-1)
        poly_w(1:nfitted) = sig(1:nfitted)

        par = 0.0_dp ; f = 0.0_dp ; dfda = 0.0_dp
        CALL elsunc ( &
             par, poly_order, nfitted, nfitted, poly_specfit, elbnd, blow(1:poly_order), &
             bupp(1:poly_order), p, w, exval, f(1:nfitted), dfda(1:nfitted,1:poly_order) )
        chisq = SUM  ( f(1:nfitted)**2 ) ! This gives the same CHI**2 as the NR routines
        x(1:npts) = locwvl(1:npts) - locavg

        poly_x(1:npts) = x(1:npts)
        exval = 3
        CALL poly_specfit ( &
             par(1:poly_order), poly_order, poly_y(1:npts), npts, exval, dfda, 0 )
        database(i,1:npts) = database(i,1:npts) - poly_y(1:npts)
     END IF
  END DO

  RETURN
END SUBROUTINE subtract_poly

SUBROUTINE subtract_poly_meas ( locwvl, npoints, locspec, ll_rad, lu_rad )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : max_spec_pts, elsunc_np, elsunc_nw
  USE OMSAO_variables_module,  ONLY : poly_x, poly_y, poly_w, poly_order
  USE bounded_nonlin_LS,       ONLY : elsunc

  IMPLICIT NONE

  EXTERNAL poly_specfit

  INTEGER,                             INTENT (IN) :: npoints, ll_rad, lu_rad
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN) :: locwvl

  REAL (KIND=dp), DIMENSION (npoints), INTENT (INOUT) :: locspec

  REAL (KIND=dp), DIMENSION (npoints)              :: tmp, ptmp, sig
  REAL (KIND=dp), DIMENSION (poly_order)           :: par, polfunc, r
  REAL (KIND=dp), DIMENSION (poly_order,poly_order):: covar

  INTEGER                             :: i, j, nlower, nupper, nfitted
  REAL (KIND=dp)                      :: locavg, chisq
  REAL (KIND=dp), DIMENSION (npoints) :: x

  ! ================
  ! ELSUNC variables
  ! ================
  INTEGER                                       :: exval
  REAL (KIND=dp), DIMENSION (npoints)           :: fitres
  REAL (KIND=dp)                                :: rms, chisq2
  INTEGER                                       :: elbnd
  INTEGER,        DIMENSION (11)                :: p
  REAL (KIND=dp), DIMENSION (6)                 :: w
  REAL (KIND=dp), DIMENSION (poly_order)         :: blow, bupp
  REAL (KIND=dp), DIMENSION (npoints)           :: f, fitspec
  REAL (KIND=dp), DIMENSION (npoints,poly_order) :: dfda


  ! ======================
  ! Assign fitting weights
  ! ======================
  sig = 1.0_dp

  !     Find limits for polynomial fitting, with ~1 nm overlap
  nlower = MINVAL(MINLOC( locwvl(1:npoints), MASK=(locwvl(1:npoints) >= locwvl(ll_rad)-1.0) ))
  nupper = MAXVAL(MAXLOC( locwvl(1:npoints), MASK=(locwvl(1:npoints) <= locwvl(lu_rad)+1.0) ))
  nfitted = nupper - nlower + 1

  !     Find average position over fitted region
  locavg = SUM ( locwvl(1+nlower-1:nfitted+nlower-1) ) / REAL ( nfitted, KIND=dp )

  !     Load temporary position file: re-define positions in order to fit
  !     about mean position
  DO i = 1, nfitted
     ptmp(i) = locwvl(i+nlower-1) - locavg
  END DO

  !     Load and fit spectrum
  tmp(1:nfitted) = locspec(1+nlower-1:nfitted+nlower-1)

  ! ===============================================================
  ! ELBND: 0 = unconstrained
  !        1 = all variables have same lower bound
  !        else: lower and upper bounds must be supplied by the use
  ! ===============================================================  
  elbnd = 0  ;  exval = 0
  p   = -1 ; p(1) = 0  ;  p(3) = 5  ; w = -1.0
  blow(1:poly_order) = 0.0  ;  bupp(1:poly_order) = 0.0
  
  poly_x(1:nfitted) = ptmp(1:nfitted)
  poly_y(1:nfitted) =  tmp(1:nfitted)
  poly_w(1:nfitted) =  sig(1:nfitted)

  r = 0.0 ; f = 0.0 ; dfda = 0.0

  CALL elsunc ( &
       r, poly_order, nfitted, nfitted, poly_specfit, elbnd, blow(1:poly_order), &
       bupp(1:poly_order), p, w, exval, f(1:nfitted), dfda(1:nfitted,1:poly_order) )
  chisq2 = SUM  ( f(1:nfitted)**2 ) ! This gives the same CHI**2 as the NR routines

  ! Re-load spec with high-pass filtered data, over whole  spectral region
  x(1:npoints) = locwvl(1:npoints) - locavg

  poly_x(1:npoints) = x(1:npoints)
  exval = 3
  CALL poly_specfit ( &
       r(1:poly_order), poly_order, poly_y(1:npoints), npoints, exval, dfda, 0 )
  locspec(1:npoints) = locspec(1:npoints) - poly_y(1:npoints)

  RETURN
END SUBROUTINE subtract_poly_meas


SUBROUTINE poly_fit (locwvl, npoints, locspec, ll_rad, lu_rad, r )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : max_spec_pts, elsunc_np, elsunc_nw
  USE OMSAO_variables_module,  ONLY : poly_x, poly_y, poly_w, poly_order
  USE bounded_nonlin_LS,       ONLY : elsunc

  IMPLICIT NONE

  EXTERNAL poly_specfit

  INTEGER,                             INTENT (IN) :: npoints, ll_rad, lu_rad
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN) :: locwvl

  REAL (KIND=dp), DIMENSION (npoints),    INTENT (INOUT) :: locspec
  REAL (KIND=dp), DIMENSION (poly_order), INTENT (OUT)   :: r

  REAL (KIND=dp), DIMENSION (npoints)              :: tmp, ptmp, sig
  REAL (KIND=dp), DIMENSION (poly_order)           :: par, polfunc
  REAL (KIND=dp), DIMENSION (poly_order,poly_order):: covar

  INTEGER                             :: i, j, nlower, nupper, nfitted
  REAL (KIND=dp)                      :: locavg, chisq
  REAL (KIND=dp), DIMENSION (npoints) :: x

  ! ================
  ! ELSUNC variables
  ! ================
  INTEGER                                       :: exval
  REAL (KIND=dp), DIMENSION (npoints)           :: fitres
  REAL (KIND=dp)                                :: rms, chisq2
  INTEGER                                       :: elbnd
  INTEGER,        DIMENSION (11)                :: p
  REAL (KIND=dp), DIMENSION (6)                 :: w
  REAL (KIND=dp), DIMENSION (poly_order)         :: blow, bupp
  REAL (KIND=dp), DIMENSION (npoints)           :: f, fitspec
  REAL (KIND=dp), DIMENSION (npoints,poly_order) :: dfda

  ! ======================
  ! Assign fitting weights
  ! ======================
  sig = 1.0_dp

  !     Find limits for polynomial fitting, with ~1 nm overlap
  nlower = MINVAL(MINLOC( locwvl(1:npoints), MASK=(locwvl(1:npoints) >= locwvl(ll_rad)-1.0) ))
  nupper = MAXVAL(MAXLOC( locwvl(1:npoints), MASK=(locwvl(1:npoints) <= locwvl(lu_rad)+1.0) ))
  nfitted = nupper - nlower + 1

  !     Find average position over fitted region
  locavg = SUM ( locwvl(1+nlower-1:nfitted+nlower-1) ) / REAL ( nfitted, KIND=dp )

  !     Load temporary position file: re-define positions in order to fit
  !     about mean position
  DO i = 1, nfitted
     ptmp(i) = locwvl(i+nlower-1) - locavg
  END DO

  !     Load and fit spectrum
  tmp(1:nfitted) = locspec(1+nlower-1:nfitted+nlower-1)

  ! ===============================================================
  ! ELBND: 0 = unconstrained
  !        1 = all variables have same lower bound
  !        else: lower and upper bounds must be supplied by the use
  ! ===============================================================  
  elbnd = 0  ;  exval = 0
  p   = -1 ; p(1) = 0  ;  p(3) = 5  ; w = -1.0
  blow(1:poly_order) = 0.0  ;  bupp(1:poly_order) = 0.0
  
  poly_x(1:nfitted) = ptmp(1:nfitted)
  poly_y(1:nfitted) =  tmp(1:nfitted)
  poly_w(1:nfitted) =  sig(1:nfitted)

  r = 0.0 ; f = 0.0 ; dfda = 0.0

  CALL elsunc ( &
       r, poly_order, nfitted, nfitted, poly_specfit, elbnd, blow(1:poly_order), &
       bupp(1:poly_order), p, w, exval, f(1:nfitted), dfda(1:nfitted,1:poly_order) )
  chisq2 = SUM  ( f(1:nfitted)**2 ) ! This gives the same CHI**2 as the NR routines

  ! Re-load spec with high-pass filtered data, over whole  spectral region
  x(1:npoints) = locwvl(1:npoints) - locavg

  poly_x(1:npoints) = x(1:npoints)
  exval = 3
  CALL poly_specfit ( &
       r(1:poly_order), poly_order, poly_y(1:npoints), npoints, exval, dfda, 0 )
  locspec(1:npoints) = poly_y(1:npoints)

  RETURN
END SUBROUTINE poly_fit
