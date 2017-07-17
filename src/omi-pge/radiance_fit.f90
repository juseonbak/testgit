SUBROUTINE radiance_fit ( fitcol, rms, dfitcol, exval )

  ! ***************************************************************
  !
  !   Perform solar wavelength calibration and slit width fitting
  !
  ! ***************************************************************

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: solar_idx, n_max_fitpars, wvl_idx, spc_idx, sig_idx
  USE OMSAO_parameters_module,   ONLY: missing_value_dp
  USE OMSAO_variables_module,    ONLY: &
       n_mol_fit, fitcol_idx, n_fincol_idx, fincol_idx, pm_one, database, refspec_norm, &
       renorm, weight_rad, yn_doas, yn_smooth, &
       rad_wav_avg, fitvar_rad, fitvar_rad_init, fitvar_rad_saved, &
       n_fitvar_rad, chisq, sol_zen_eff, &
       sza_atm, vza_atm, n_rad_wvl,  &
       lo_radbnd, up_radbnd, fitweights, currspec, fitwavs, rad_winwav_idx,    &
       mask_fitvar_rad, max_itnum_rad, curr_rad_spec, fitvar_rad_str

  IMPLICIT NONE

  ! *******************************************************************
  ! CAREFUL: Assumes that radiance and solar wavelength arrays have the
  ! same number of points. That must not be the case if we read in a
  ! general EL1 file. Examine and adjust! (tpk, note to himself)
  ! *******************************************************************

  ! ===============
  ! Input variables
  ! ===============

  ! ==================
  ! Modified variables
  ! ==================
  INTEGER,            INTENT (OUT) :: exval
  REAL     (KIND=dp), INTENT (OUT) :: fitcol, rms, dfitcol

  ! ===============
  ! Local variables
  ! ===============
  INTEGER        :: i, j, idx, n_fitwav_rad, j1, j2, k1, k2, l, ll_rad, lu_rad
  REAL (KIND=dp) :: asum, ssum, remult, sigsum
  REAL (KIND=dp), DIMENSION (n_rad_wvl) :: fitres, fitspec, tmp, spec, sig, pos
  REAL (KIND=dp), DIMENSION (n_max_fitpars) :: fitvar, lobnd, upbnd
  REAL (KIND=dp), DIMENSION (n_fitvar_rad, n_fitvar_rad) :: covar_matrix

  INTEGER, DIMENSION(n_fitvar_rad) :: iw
  REAL (KIND=dp), DIMENSION(n_fitvar_rad*(n_rad_wvl+5)+n_rad_wvl) :: wa
  INTEGER        :: iopt, nprint, info, lwa
  REAL (KIND=dp) :: tol

  EXTERNAL specfit_func


  ! ============================================================
  ! Assign LL_RAD, LU_RAD, and SIG for each earthshine radiance
  ! ============================================================
  ll_rad = rad_winwav_idx(2)  ;  lu_rad = rad_winwav_idx(3)

  pos (1:n_rad_wvl) = curr_rad_spec(wvl_idx,1:n_rad_wvl)
  spec(1:n_rad_wvl) = curr_rad_spec(spc_idx,1:n_rad_wvl)
  sig (1:n_rad_wvl) = curr_rad_spec(sig_idx,1:n_rad_wvl)

  !   renormalize, if requested (normally don't do for DOAS since we
  !   are dividing by the solar spectrum just below, and
  !   logarithmizing)
  IF (renorm) THEN
     remult = 0.0_dp  ; sigsum = 0.0_dp
     IF ( weight_rad ) THEN
        remult = SUM ( spec(1:n_rad_wvl) / (sig(1:n_rad_wvl)*sig(1:n_rad_wvl) ) )
        sigsum = SUM (             1.0_dp / (sig(1:n_rad_wvl)*sig(1:n_rad_wvl) ) )
        remult = remult / sigsum
     ELSE
        remult = SUM ( spec(1:n_rad_wvl) ) / REAL ( n_rad_wvl, KIND=dp )
     END IF
     spec (1:n_rad_wvl) = spec (1:n_rad_wvl) / remult
  END IF

  !   High pass filtering for DOAS. First, take log (rad/irrad), then
  !   filter by subtracting a cubic, then re-add the log (irradiance).
  !   This way we are fitting the log (rad) with the proper filtering
  !   having been done. The spectrum subroutine will start by
  !   subtracting the irradiance spectrum, where the proper shifting
  !   and squeezing can take place. In this high-pass filtering, we
  !   are ignoring the small extra wavelength calibration change, as
  !   for Ring effect, above.
  IF ( yn_doas ) THEN
     spec(1:n_rad_wvl) = LOG ( spec(1:n_rad_wvl) / database(solar_idx,1:n_rad_wvl) )
     CALL subtract_cubic_meas (pos, n_rad_wvl, spec, ll_rad, lu_rad)
     spec(1:n_rad_wvl) = spec(1:n_rad_wvl) + LOG ( database(solar_idx,1:n_rad_wvl) )
  END IF

  !   Apply smoothing (1/16,1/4,3/8,1/4,1/16); 2/98 uhe/ife recommendation
  IF ( yn_smooth ) THEN
     tmp(1:n_rad_wvl) = spec(1:n_rad_wvl)
     spec (3:n_rad_wvl-2) = 0.375_dp * tmp (3:n_rad_wvl-2) +  &
          0.25_dp   * (tmp (4:n_rad_wvl-1) + tmp (2:n_rad_wvl-3)) +  &
          0.0625_dp * (tmp (5:n_rad_wvl) + tmp (1:n_rad_wvl-4))
  END IF
  IF ( weight_rad ) THEN
     asum = SUM ( pos(1:n_rad_wvl) / ( sig(1:n_rad_wvl)*sig(1:n_rad_wvl) ) )
     ssum = SUM (         1.0_dp   / ( sig(1:n_rad_wvl)*sig(1:n_rad_wvl) ) )
     rad_wav_avg = asum / ssum
  ELSE
     rad_wav_avg = (pos (n_rad_wvl) + pos (1)) / 2.0_dp
  END IF


  exval = 0
  !     BOAS or DOAS fitting now done here        
     
  ! ELSUNC fitting
  fitwavs   (1:n_rad_wvl) = pos (1:n_rad_wvl)
  fitweights(1:n_rad_wvl) = sig (1:n_rad_wvl)
  currspec  (1:n_rad_wvl) = spec(1:n_rad_wvl)
  !fitvar_rad(1:n_max_fitpars) = fitvar_rad_init(1:n_max_fitpars)
  fitvar_rad(1:n_max_fitpars) = fitvar_rad_saved(1:n_max_fitpars)

  ! -----------------------------------------------------------------
  ! Create a condensed array of fitting variables that only contains
  ! the ones that are varied. This considerably reduces the execution
  ! time of the fitting routine.
  ! -----------------------------------------------------------------
  fitvar = 0.0_dp ; lobnd = 0.0_dp ; upbnd = 0.0_dp
  DO i = 1, n_fitvar_rad
     idx       = mask_fitvar_rad(i)
     fitvar(i) = fitvar_rad(idx)
     lobnd(i)  = lo_radbnd(idx)
     upbnd(i)  = up_radbnd(idx)
  END DO

  iopt   = 1
  tol    = 1.0E-3_dp
  nprint = 0
  info   = 0
  lwa = n_fitvar_rad*(n_rad_wvl+5) + n_rad_wvl
  
  CALL specfit ( &
       n_fitvar_rad, fitvar(1:n_fitvar_rad), n_fitvar_rad, n_rad_wvl, &
       lobnd(1:n_fitvar_rad), upbnd(1:n_fitvar_rad), max_itnum_rad, &
       rms, chisq, covar_matrix(1:n_fitvar_rad, 1:n_fitvar_rad),    &
       fitspec(1:n_rad_wvl), fitres(1:n_rad_wvl), exval, specfit_func )

  !WRITE (99,'(I6)')             n_rad_wvl
  !WRITE (99,'(200(1F8.3:))')    fitwavs(1:n_rad_wvl)
  !WRITE (99,'(200(1P1E13.5:))') spec   (1:n_rad_wvl)
  !WRITE (99,'(200(1P1E13.5:))') sig    (1:n_rad_wvl)
  !WRITE (99,'(200(1P1E13.5:))') fitspec(1:n_rad_wvl)
  !WRITE (99,'(200(1P1E13.5:))') fitres (1:n_rad_wvl)

  fitvar_rad_saved(1:n_max_fitpars) = fitvar_rad(1:n_max_fitpars)
  !fitvar_rad_saved(1:n_max_fitpars) = fitvar_rad_init(1:n_max_fitpars)

  ! --------------------------------------------------------------------
  ! At this point we have to make sure we only report the fitted columns
  ! if EXVAL, the exit variable from the fitting routine, is >0. Because
  ! all negative values mean trouble and are likely to have produced
  ! negative uncertainties. In any case, the values for the fitted 
  ! columns should not be trusted.
  ! --------------------------------------------------------------------
  SELECT CASE ( exval )
  CASE ( 1: )
     ! =====================================================================
     ! Compute the actual number of radiance wavelengths used in the fitting
     ! =====================================================================
     n_fitwav_rad = &
          INT ( SUM ( 1.0_dp / (fitweights(1:n_rad_wvl)*fitweights(1:n_rad_wvl)) ) )

     ! --------------------------------------------------------------------------
     ! Assign total column. We have done the preliminary work with FITCOL_IDX and
     ! PM_ONE, so we don't have to "IF DOAS" here.
     !
     ! Reminder about the (admittedly confusing) variable names:
     !
     !   FITCOL:       Fitted column
     !   DFITCOL:      Uncertainty of fitted column
     !   N_FINCOL_IDX: Number of fitting variables that make up the final FITCOL
     !                 (e.g., O3 at more than one temperature).
     !   FINCOL_IDX:   Array of dimension (2,N_FINCOL_IDX*MXS_IDX), where MXS_IDX
     !                 is the maximum fitting sub-index (e.g., AD1_IDX, LBE_IDX).
     !                 FINCOL_IDX(1,*) carries the relative indices of the varied
     !                 final column variables in the array that was passed to
     !                 the fitting routine. FINCOL_IDX(2,*) contains the index
     !                 for the associated reference spectrum; this we need for
     !                 access to the normalization factor. See subroutine
     !                 READ_CONTROL_FILE for assignment of these indices.
     ! --------------------------------------------------------------------------
     fitcol = 0.0_dp  ;  dfitcol = 0.0_dp
     DO i = 1, n_fincol_idx
        ! --------------------------------------------------
        ! First add the contribution of the diagonal element
        ! --------------------------------------------------
        j1 = fincol_idx(1, i) ; k1 = fincol_idx(2,i)
        fitcol  = fitcol  + pm_one * fitvar(j1) / refspec_norm(k1)
        dfitcol = dfitcol + covar_matrix(j1,j1) / (refspec_norm(k1)*refspec_norm(k1))
        !WRITE (*,'(a,2I3,1p3E15.5)') '(a)', i, j1, fitvar(j1), covar_matrix(j1,j1), refspec_norm(k1)

        ! --------------------------------------------------------------------
        ! Then add the contributions from off-diagonal elements (correlations)
        ! --------------------------------------------------------------------
        DO l = i+1, n_fincol_idx
           j2 = fincol_idx(1,l) ; k2 = fincol_idx(2,l)

           dfitcol = dfitcol + 2.0_dp * covar_matrix(j2,l) / (refspec_norm(k1)*refspec_norm(k2))
        END DO
     END DO
  
     !WRITE (*,'(a,2I3,1p2E15.5 //)') '(c)', n_rad_wvl, n_fitvar_rad, dfitcol, rms
     dfitcol = rms * SQRT ( dfitcol *  &
          REAL ( n_rad_wvl, KIND=dp ) / REAL ( n_fitwav_rad - n_fitvar_rad, KIND=dp ) )

     !DO j1 = 1, n_fitvar_rad
     !   DO j2 = 1, n_fitvar_rad
     !      WRITE (*,'(a,2I3,1p2E15.5)') 'covar:', j1, j2, covar_matrix(j1,j2), fitvar(j1)
     !   END DO
     !END DO

     !WRITE (*,*)
     !DO j1 = 1, n_fitvar_rad
     !   !IF ( fitvar(j1) /= 0.0_dp ) WRITE (*,'(I3,1pE12.4 $)') j1, fitvar(j1)
     !   idx       = mask_fitvar_rad(j1)
     !   WRITE (*,'(1X,A,1X,1pE12.4 $)') TRIM(ADJUSTL(fitvar_rad_str(idx))), fitvar(j1)
     !END DO
     !WRITE (*,*)

     ! ======================================================
     ! Anything but EXVAL > 0 means that trouble has occurred, 
     ! and the fit most likely is screwed.
     ! ======================================================
  CASE ( :0 )
     fitcol = missing_value_dp ; dfitcol = missing_value_dp
  CASE DEFAULT
     ! -------------------------------------------------------------------
     ! We should never reach here, because the above CASE statements cover
     ! all possible values of EXVAL. But better safe than sorry.
     fitcol = missing_value_dp ; dfitcol = missing_value_dp
     ! -------------------------------------------------------------------
  END SELECT

  ! --------------------------------------------------------------------------
  ! Write fitting results to file. Should be moved to MAIN program eventually,
  ! since it's here that we will have to meddle with HDF5 output.
  ! --------------------------------------------------------------------------
       
  RETURN
END SUBROUTINE radiance_fit
