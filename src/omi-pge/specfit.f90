SUBROUTINE specfit ( &
     nfitvar, fitvar, n_fitvar_rad, nspecpts, lowbnd, uppbnd, max_itnum, &
     rms, chisq, covar, fitspec, fitres, stderr, exval, fitfunc )

  USE OMSAO_precision_module
  USE OMSAO_indices_module,     ONLY: elsunc_userdef
  USE OMSAO_parameters_module,  ONLY: max_spec_pts, elsunc_np, elsunc_nw
  USE OMSAO_variables_module,   ONLY: tol,  epsrel,  epsabs,  epsx, fitwavs, &
       currspec, fitweights
  USE bounded_nonlin_LS,        ONLY: elsunc

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                              INTENT (IN) :: nfitvar, nspecpts, &
       n_fitvar_rad, max_itnum
  REAL (KIND=dp), DIMENSION (nfitvar),  INTENT (IN) :: lowbnd, uppbnd

  ! ==================
  ! Modified variables
  ! ==================
  REAL (KIND=dp), DIMENSION (nfitvar),  INTENT (INOUT) :: fitvar

  ! ================
  ! Output variables
  ! ================
  INTEGER,                                     INTENT (OUT) :: exval
  REAL (KIND=dp), DIMENSION (nfitvar, nfitvar),INTENT (OUT) :: covar
  REAL (KIND=dp), DIMENSION (nspecpts),        INTENT (OUT) :: fitres, fitspec
  REAL (KIND=dp),                              INTENT (OUT) :: rms, chisq
  REAL (KIND=dp), DIMENSION (nfitvar),         INTENT (OUT) :: stderr

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                                      :: elbnd, i, j
  INTEGER,        DIMENSION (elsunc_np)        :: p ! elsunc_np 11
  REAL (KIND=dp), DIMENSION (elsunc_nw)        :: w ! elsunc_nw 6
  REAL (KIND=dp), DIMENSION (nfitvar)          :: blow, bupp
  REAL (KIND=dp), DIMENSION (nspecpts)         :: f
  REAL (KIND=dp), DIMENSION (nspecpts,nfitvar) :: dfda
  REAL (KIND=dp), DIMENSION (nfitvar, nfitvar) :: correl

  EXTERNAL fitfunc

  ! ===============================================================
  ! ELBND: 0 = unconstrained
  !        1 = all variables have same lower bound
  !        else: lower and upper bounds must be supplied by the use
  ! ===============================================================  
  elbnd = elsunc_userdef !(2)

  exval = 0
  
  p   = -1    ;  p(1)   = 0  ;  p(3) = max_itnum
  w   = -1.0  ;  w(1:4) = (/ tol,  epsrel,  epsabs,  epsx /)
 
  blow(1:nfitvar) = lowbnd(1:nfitvar)
  bupp(1:nfitvar) = uppbnd(1:nfitvar)
 
  CALL elsunc ( fitvar(1:nfitvar), nfitvar, nspecpts, nspecpts, fitfunc,    &
       elbnd, blow(1:nfitvar),  bupp(1:nfitvar), p, w, exval, f(1:nspecpts),&
       dfda(1:nspecpts,1:nfitvar) )

  ! ------------------------------------------------------------------
  ! Call to ELSUNC fitting function to obtain complete fitted spectrum
  ! ------------------------------------------------------------------
 
  CALL fitfunc ( fitvar, nfitvar, fitspec, nspecpts, 3, &
       dfda(1:nspecpts,1:nfitvar), 0) 

  ! ---------------------------------------------------------------
  ! Compute fitting residual
  ! FITRES is the negative of the returned function F = Model-Data.
  ! ---------------------------------------------------------------
  fitres(1:nspecpts) = -f(1:nspecpts)

  ! Covariance matrix.
  ! ------------------
  covar(1:nfitvar,1:nfitvar) = dfda(1:nfitvar,1:nfitvar)

  ! Fitting RMS and CHI**2
  ! ----------------------
  rms  = SQRT(SUM(fitres(1:nspecpts)**2 ) / nspecpts)
  ! This gives the same CHI**2 as the NR routines
  chisq = SUM  (fitres(1:nspecpts)**2 ) 

  ! compute standard deviation for each variable
  DO i = 1, nfitvar
     stderr(i) = rms * SQRT(covar(i, i) * nspecpts / (nspecpts - nfitvar))
  END DO

  !compulte the correlation matrix
  !WRITE(*, *) 'Correlation matrix = '
  !DO i = 1, nfitvar
  !   correl(i, i) = 1.0
  !   DO j = 1, i - 1
  !      correl(i, j) = covar(i, j)/SQRT (covar(i, i) * covar(j, j))       
  !   END DO
  !   WRITE(*, '(i3, 2d12.4, 3x, 40d10.2)') i, fitvar(i), stderr(i), correl(i, 1:i)
  !END DO
  
  RETURN
END SUBROUTINE specfit
