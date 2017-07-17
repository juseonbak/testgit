SUBROUTINE specfit_func_sol ( fitvar, nfitvar, ymod, npoints, ctrl, dyda, mdy )

  !
  ! Calculates the Solar spectrum and its derivatives for ELSUNC
  !
  ! NOTE: the variable DYDA as required here is the transpose of that 
  !       rquired for the Numerical Recipes
  !

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : max_spec_pts
  USE OMSAO_variables_module, ONLY : &
       yn_smooth, yn_doas, fitwavs, fitweights, currspec, sol_wav_avg

  IMPLICIT NONE

  INTEGER,                                     INTENT (IN)    :: nfitvar, npoints, mdy
  INTEGER,                                     INTENT (INOUT) :: ctrl
  REAL (KIND=dp), DIMENSION (nfitvar),         INTENT (IN)    :: fitvar
  REAL (KIND=dp), DIMENSION (npoints),         INTENT (INOUT) :: ymod
  REAL (KIND=dp), DIMENSION (npoints,nfitvar), INTENT (INOUT) :: dyda

  INTEGER                             :: i
  REAL (KIND=dp), DIMENSION (nfitvar) :: vartmp
  REAL (KIND=dp), DIMENSION (npoints) :: dyplus, dyminus, locwvl

  locwvl(1:npoints) = fitwavs(1:npoints)

  SELECT CASE ( ABS ( ctrl ) )
  CASE ( 1 )
     ! -----------------------------------------------------------------------
     ! Calculate the weighted difference between fitted and measured spectrum.
     ! -----------------------------------------------------------------------
     CALL spectrum_solar ( &
          npoints, nfitvar, yn_smooth, sol_wav_avg, locwvl(1:npoints), ymod(1:npoints), &
          fitvar(1:nfitvar), yn_doas )
     ymod(1:npoints) = ( ymod(1:npoints) - currspec(1:npoints) ) / fitweights(1:npoints)
!   write(*,'(a,100e15.7)') 'ymod1',ymod(1:10) 
  CASE ( 2 )
     ! ---------------------------------------------------------------------
     ! The following sets up ELSUNC for numerical computation of the fitting
     ! function derivative. It is faster and more flexible than the original
     ! "manual" (AUTODIFF) scheme, and gives better fitting uncertainties.
     ! ---------------------------------------------------------------------
     ctrl = 0; RETURN

  CASE ( 3 )
     ! Calculate the spectrum, without weighting
     CALL spectrum_solar ( &
          npoints, nfitvar, yn_smooth, sol_wav_avg, locwvl(1:npoints), ymod(1:npoints), &
          fitvar(1:nfitvar), yn_doas )
!   write(*,'(a,100e15.7)') 'ymod3',ymod(1:10) 
  CASE DEFAULT
     WRITE (*,'(A,I4)') &
          "ERROR in function ELSUNC_SPECFIT_FUNC. Don't know how to handle CTRL = ", ctrl
  END SELECT

  RETURN
END SUBROUTINE specfit_func_sol


SUBROUTINE specfit_func ( vars, npars, ymod, npoints, ctrl, dyda, mdy )

  !
  !     Calculates the spectrum and its derivatives for ELSUNC
  !
  ! NOTE: the variable DYDA as required here is the transpose of that 
  !       rquired for the Numerical Recipes
  !

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : max_spec_pts
  USE OMSAO_variables_module, ONLY : &
       database, yn_doas, yn_smooth, rad_wav_avg, fitwavs, fitweights, currspec

  IMPLICIT NONE

  INTEGER,                                   INTENT (IN)    :: npars, npoints, mdy
  INTEGER,                                   INTENT (INOUT) :: ctrl
  REAL (KIND=dp), DIMENSION (npars),         INTENT (IN)    :: vars
  REAL (KIND=dp), DIMENSION (npoints),       INTENT (INOUT) :: ymod
  REAL (KIND=dp), DIMENSION (npoints,npars), INTENT (INOUT) :: dyda

  INTEGER                                   :: i
  REAL (KIND=dp), DIMENSION (npars)         :: vartmp
  REAL (KIND=dp), DIMENSION (npoints)       :: dyplus, dyminus, locwvl

  locwvl(1:npoints) = fitwavs(1:npoints)

  SELECT CASE ( ABS ( ctrl ) )
  CASE ( 1 )
     ! Calculate the weighted difference between fitted and measured spectrum.
     CALL spectrum_earthshine ( &
          npoints, npars, yn_smooth, rad_wav_avg, locwvl(1:npoints), &
          ymod(1:npoints), vars(1:npars),database, yn_doas )
     ymod(1:npoints) = ( ymod(1:npoints) - currspec(1:npoints) )&
          / fitweights(1:npoints)

  CASE ( 2 )
     ! ---------------------------------------------------------------------
     ! The following sets up ELSUNC for numerical computation of the fitting
     ! function derivative. It is faster and more flexible than the original
     ! "manual" (AUTODIFF) scheme, and gives better fitting uncertainties.
     ! ---------------------------------------------------------------------
     ctrl = 0; RETURN

  CASE ( 3 )
     ! Calculate the spectrum, without weighting
     CALL spectrum_earthshine ( &
          npoints, npars, yn_smooth, rad_wav_avg, locwvl(1:npoints), ymod(1:npoints), &
          vars(1:npars), database, yn_doas )

  CASE DEFAULT
     WRITE (*,'(A,I4)') &
          "ERROR in function ELSUNC_SPECFIT_FUNC. Don't know how to handle CTRL = ", ctrl
  END SELECT

  RETURN
END SUBROUTINE specfit_func


SUBROUTINE cubic_func ( x, afunc, ma )

  ! ***************************************************
  !
  !   Computes a third-order polynomial. Used in LFIT
  !
  ! ***************************************************

  USE OMSAO_precision_module

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,        INTENT (IN) :: ma
  REAL (KIND=dp), INTENT (IN) :: x

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=dp), DIMENSION (ma), INTENT (OUT) :: afunc

  
  ! ===============
  ! Local variables
  ! ===============
  INTEGER :: i


  afunc(1) = 1.0_dp
  DO i = 2, ma
     afunc(i) = afunc(i-1) * x
  END DO

  RETURN
END SUBROUTINE cubic_func

SUBROUTINE cubic_specfit ( a, na, y, m, ctrl, dyda, mdy )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : max_spec_pts
  USE OMSAO_variables_module, ONLY  : cubic_x, cubic_y, cubic_w

  IMPLICIT NONE

  ! Input parameters
  ! ================
  INTEGER,                         INTENT (IN)  :: na, m, mdy
  REAL (KIND=dp), DIMENSION (na),  INTENT (IN)  :: a

  ! Modified parameters
  ! ===================
  INTEGER, INTENT (INOUT) :: ctrl

  ! Output parameters
  ! =================
  REAL (KIND=dp), DIMENSION (m),    INTENT (OUT)  :: y
  REAL (KIND=dp), DIMENSION (m,na), INTENT (OUT)  :: dyda

  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION (m) :: x, y0

  x  = cubic_x(1:m)
  y0 = a(1) + a(2)*x + a(3)*x*x + a(4)*x*x*x


  SELECT CASE ( ABS(ctrl) )
  CASE ( 1 )
     y  = ( y0 - cubic_y(1:m) ) / cubic_w(1:m)
  CASE ( 2 )
     dyda = 0.0_dp
     dyda(1:m,1) = 1.0_dp
     dyda(1:m,2) = x(1:m)
     dyda(1:m,3) = x(1:m)*x(1:m)
     dyda(1:m,4) = x(1:m)*x(1:m)*x(1:m)
  CASE ( 3 )
     ! This CASE is included to get the complete fitted spectrum
     y  = y0
  CASE DEFAULT
     WRITE (*, '(A,I3)') "Don't know how to handle CTRL = ", ctrl
  END SELECT

  RETURN
END SUBROUTINE cubic_specfit


SUBROUTINE poly_specfit ( a, na, y, m, ctrl, dyda, mdy )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : max_spec_pts
  USE OMSAO_variables_module, ONLY  : poly_x, poly_y, poly_w

  IMPLICIT NONE

  ! Input parameters
  ! ================
  INTEGER,                         INTENT (IN)  :: na, m, mdy
  REAL (KIND=dp), DIMENSION (na),  INTENT (IN)  :: a

  ! Modified parameters
  ! ===================
  INTEGER, INTENT (INOUT) :: ctrl

  ! Output parameters
  ! =================
  REAL (KIND=dp), DIMENSION (m),    INTENT (OUT)  :: y
  REAL (KIND=dp), DIMENSION (m,na), INTENT (OUT)  :: dyda

  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION (m) :: x, y0
  INTEGER                       :: i

  x  = poly_x(1:m)

  y0 = 0.0_dp
  DO i = 1, na
     y0 = y0 + a(i) * (x ** (i-1))
  END DO
 
  SELECT CASE ( ABS(ctrl) )
  CASE ( 1 )
     y  = ( y0 - poly_y(1:m) ) / poly_w(1:m)
  CASE ( 2 )
     dyda = 0.0_dp
     DO i = 1, na
        dyda(1:m, i) = x(1:m) ** (i-1)
     END DO
  CASE ( 3 )
     ! This CASE is included to get the complete fitted spectrum
     y  = y0
  CASE DEFAULT
     WRITE (*, '(A,I3)') "Don't know how to handle CTRL = ", ctrl
  END SELECT

  RETURN
END SUBROUTINE poly_specfit

