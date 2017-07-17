! ********************************************************************
! Note: make sure currspec, fitwavs, fitweights, slitcal are set properly
! ********************************************************************

SUBROUTINE cal_fit_one (n_fit_pts, n_fitvar_sol, wrt_to_screen, &
     wrt_to_file, slitcal, slit_unit, avgwav, varstd, solfit_exval)

  USE OMSAO_precision_module
  USE OMSAO_indices_module,     ONLY : max_calfit_idx, hwe_idx, asy_idx, &
       shi_idx, squ_idx, wvl_idx, spc_idx, sig_idx, sin_idx, bl0_idx,    &
       bl1_idx, bl2_idx, bl3_idx, sc0_idx, sc1_idx, sc2_idx, sc3_idx,    &
       vgl_idx, vgr_idx, hwl_idx, hwr_idx, solar_idx,hwn_idx, ops_idx
  USE OMSAO_variables_module,   ONLY: weight_sun, fitvar_sol, &
       fitvar_sol_saved,  fitvar_sol_init, chisq,  sol_wav_avg, which_slit, &
       mask_fitvar_sol, fitwavs, fitweights, currspec, lo_sunbnd, up_sunbnd,&
       max_itnum_sol, n_irrad_wvl, n_refspec_pts, refspec_orig_data, &
       refspec_norm, debug_boreas, here_stop
  USE bounded_nonlin_ls, ONLY: niter     
  IMPLICIT NONE

  ! ================
  ! Input variables
  ! ================
  INTEGER, INTENT(IN) :: slit_unit, n_fit_pts, n_fitvar_sol
  LOGICAL, INTENT(IN) :: wrt_to_screen, wrt_to_file, slitcal

  ! ================
  ! Output variables
  ! ================
  INTEGER,            INTENT (OUT)              :: solfit_exval
  REAL (KIND = dp),   INTENT (OUT)              :: avgwav  
  REAL (KIND=dp), DIMENSION (max_calfit_idx, 2) :: varstd

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp)  :: asum, ssum, rms, sum_sig2
  REAL (KIND=dp), DIMENSION (n_fit_pts)                 :: fitres, fitspec
  REAL (KIND=dp), DIMENSION (n_fitvar_sol, n_fitvar_sol):: covar
  REAL (KIND=dp)  :: hw1e, e_asym, vgl, vgr, hwl, hwr, sin, shi, squ, ops,pwr
  REAL (KIND=dp), DIMENSION (n_fitvar_sol) :: fitvar, lobnd, upbnd, stderr
  INTEGER         :: i, ref_fidx, ref_lidx, ref_all_pts, n_ref_pts

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=11), PARAMETER :: modulename = 'cal_fit_one'

  EXTERNAL specfit_func_sol
    
  ! -------------------------------------------------------
  ! Calculate SOL_WAV_AVG of measured solar spectra here,
  ! for use in calculated spectra.
  ! -------------------------------------------------------
  IF ( weight_sun ) THEN
     asum = SUM ( fitwavs(1:n_fit_pts) / (fitweights(1:n_fit_pts) &
          * fitweights(1:n_fit_pts)))
     ssum = SUM (1.0/ (fitweights(1:n_fit_pts)*fitweights(1:n_fit_pts)))
     sol_wav_avg = asum / ssum
  ELSE
     sol_wav_avg = ( fitwavs(1) + fitwavs(n_fit_pts)) / 2.0
  END IF
  !print * , n_fitvar_sol
  ! --------------------------------------------------------------
  ! ELSUNC FIT: Calculate and iterate on the irradiance spectrum.
  ! --------------------------------------------------------------     
  ! initialize the sin variables for faster convergence
  ref_all_pts = n_refspec_pts(solar_idx)   
  ref_fidx=MINVAL(MAXLOC(refspec_orig_data(solar_idx, 1:ref_all_pts, wvl_idx),&
       MASK=(refspec_orig_data(solar_idx, 1:ref_all_pts, wvl_idx) <= fitwavs(1))))
  
  ref_lidx=MINVAL(MINLOC(refspec_orig_data(solar_idx, 1:ref_all_pts, wvl_idx),&
       MASK=(refspec_orig_data(solar_idx,1:ref_all_pts,wvl_idx)>=fitwavs(n_fit_pts))))
 
  IF (ref_fidx > ref_all_pts .OR. ref_fidx <= 0 .OR. ref_lidx > ref_all_pts &
       .OR. ref_lidx <= 0 .OR. ref_fidx > ref_lidx) THEN
     WRITE(*, *) ref_fidx, ref_lidx, n_fit_pts
     WRITE(*, *) fitwavs(1), fitwavs(n_fit_pts)
     WRITE(*, *) refspec_orig_data(solar_idx, ref_fidx, wvl_idx), refspec_orig_data(solar_idx, ref_lidx, wvl_idx)
     WRITE(*, *) modulename,' : Solar spectra not cover fitting window!!!'; STOP
  END IF
     
  fitvar_sol(sin_idx) = SUM(currspec(1:n_fit_pts)) * (ref_lidx-ref_fidx + 1.0) / &
       SUM(refspec_orig_data(solar_idx, ref_fidx:ref_lidx, spc_idx)) /n_fit_pts

  fitvar = 0.0; lobnd = 0.0; upbnd = 0.0
  fitvar(1:n_fitvar_sol) = fitvar_sol(mask_fitvar_sol(1:n_fitvar_sol))
  lobnd(1:n_fitvar_sol) = lo_sunbnd(mask_fitvar_sol(1:n_fitvar_sol)) 
  upbnd(1:n_fitvar_sol) = up_sunbnd(mask_fitvar_sol(1:n_fitvar_sol))  
  fitspec(1:n_fit_pts) = 0.0
  CALL specfit ( &
       n_fitvar_sol, fitvar(1:n_fitvar_sol), n_fitvar_sol, n_fit_pts, &
       lobnd(1:n_fitvar_sol), upbnd(1:n_fitvar_sol), max_itnum_sol,   &
       rms, chisq, covar(1:n_fitvar_sol,1:n_fitvar_sol), fitspec(1:n_fit_pts), &
       fitres(1:n_fit_pts), stderr, solfit_exval, specfit_func_sol)
  
 !DO i = 1, n_fit_pts
 !   WRITE(90, '(f10.4, 3D14.6)') fitwavs(i), currspec(i), fitspec(i), fitres(i)
 !ENDDO

  fitvar_sol_saved = fitvar_sol
  !standard variables variables and standard error
  varstd(1:max_calfit_idx, 1:2) = 0.0
  avgwav = ( fitwavs(1) + fitwavs(n_fit_pts)) / 2.0  ! wincal_wav
  
  DO i = 1, max_calfit_idx
     varstd(i, 1) = fitvar_sol(i)
  ENDDO
  DO i = 1, n_fitvar_sol
     varstd(mask_fitvar_sol(i), 2) = stderr(i)
  END DO
  
  IF (wrt_to_screen .OR. (wrt_to_file .AND. solfit_exval > 0)) THEN
     IF (which_slit == 3) THEN
        hw1e = fitvar_sol(hwe_idx); e_asym = 0.0
     ELSE IF (which_slit == 5) THEN
        hw1e = 0.0; e_asym = 0.0; ops=fitvar_sol(ops_idx)
     ELSE IF (which_slit == 2) THEN
        vgl  = fitvar_sol(vgl_idx);  vgr    = fitvar_sol(vgr_idx)
        hwl  = fitvar_sol(hwl_idx);  hwr    = fitvar_sol(hwr_idx)
     ELSE
        hw1e = fitvar_sol(hwe_idx); e_asym =  fitvar_sol(asy_idx); pwr = fitvar_sol (hwn_idx)
     ENDIF
     shi = fitvar_sol(shi_idx); squ = fitvar_sol(squ_idx); sin = fitvar_sol(sin_idx)
  ENDIF
  
  IF (wrt_to_screen) THEN
     IF ( which_slit /= 2) THEN
        WRITE(*, '(5(A, 1pd14.6), A, I6,A,I6)') 'hwle = ',  hw1e,  &
             ' e_asym =', e_asym, ' rms = ', rms,'hwn =',pwr,' ops=',ops, ' exval = ', solfit_exval, 'Niter= ', niter
     ELSE
        WRITE(*, '(5(A,1pd14.6),A,I6)') 'vgl = ',  vgl, ' vgr = ',  vgr,  &
             ' hwl = ', hwl, ' hwr = ', hwr,' rms = ',rms,' exval = ', solfit_exval
     END IF
  
     WRITE(*, '(3(A, 1pd14.6))') 'shi = ', shi, ' squ =', squ, ' sin = ', sin
     WRITE(*, '(4(A, 1pd14.6))') 'b10 = ', fitvar_sol(bl0_idx), ' b11 = ', &
          fitvar_sol(bl1_idx),' b12 = ', fitvar_sol(bl2_idx), &
          ' b13 = ', fitvar_sol(bl3_idx)
     WRITE(*, '(4(A, 1pd14.6))') 'sc0 = ', fitvar_sol(sc0_idx), ' sc1 = ', &
          fitvar_sol(sc1_idx), ' sc2 = ', fitvar_sol(sc2_idx), &
          ' sc3 = ', fitvar_sol(sc3_idx)
  END IF
  IF (wrt_to_file ) WRITE(slit_unit, '(i6)') solfit_exval
  IF (wrt_to_file .AND. solfit_exval > 0) THEN
     IF ( slitcal) THEN   ! for solar slit width calibration
        IF (which_slit /= 2) THEN
           WRITE(slit_unit, '(f8.3,1p8d11.3)') avgwav, shi, squ, sin, rms, & 
                                               hw1e, e_asym, pwr,ops
        ELSE
           WRITE(slit_unit, '(f8.3,1p8d11.3)') avgwav, shi, squ, sin, rms, & 
                                                  vgl, vgr, hwl, hwr
        END IF
     ELSE                ! for solar/rad wavelength calibration
        WRITE(slit_unit, '(f8.3,1p4d11.3, I6, 1pd11.3)') avgwav, shi, squ, sin, rms, & 
                                                         solfit_exval, varstd(shi_idx, 2)
     END IF
     DO i = 1, n_fit_pts
        write(slit_unit,'(f8.3, 1p3d11.3)') fitwavs(i), currspec(i), fitspec(i), fitres(i)
     ENDDO
  END IF
  currspec = fitspec
  RETURN

END SUBROUTINE cal_fit_one
