  ! *********************** Modification History ********
  ! xliu: 
  ! 1. Add call to asym_gauss_vary if vary slit width 
  !    option is set
  ! 2. Add variables hw1earr, slitwav, e_asymarr, n_slit_pts
  !    slitwav,  from OMSAO_variables_module
  ! 3. Add arguments vgl, vgr, hwl, hwr
  ! 4. Use voigt_gauss to replace asym_gauss
  ! *****************************************************
SUBROUTINE undersample (n_gome_pts, curr_wvl, phase, pge_error_status )

  !  Convolves input spectrum with Gaussian slit function of specified
  !  HW1e, and samples at a particular input phase to give the OMI
  !  undersampling spectrum. This version calculates both phases of the
  !  undersampling spectrum, phase1 - i.e., underspec (1, i) - being the
  !  more common in OMI spectra.

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: solar_idx, us1_idx, us2_idx, wvl_idx, &
       spc_idx, hwe_idx, hwr_idx, hwl_idx, asy_idx, ops_idx, hwn_idx
  USE OMSAO_parameters_module, ONLY: max_spec_pts, maxchlen
  USE OMSAO_variables_module,  ONLY: n_refspec_pts, refspec_orig_data, &
       database, yn_varyslit, which_slit, database_shiwf, &
      nslit, nslit_rad, nslit_sol,slitwav, slitwav_sol, slitwav_rad, slitfit, &
      solslitfit, radslitfit, slit_rad, numwin, nradpix,solwinfit, &
      sring_fidx, sring_lidx, nsol_ring, sol_spec_ring, have_undersampling, refnhextra, &
      instrument_idx, omps_idx, omi_idx
  USE ozprof_data_module,     ONLY: nsl
  USE OMSAO_slitfunction_module
  !USE OMPS_slit_module
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,                                INTENT (IN) :: n_gome_pts
  REAL (KIND=dp),                         INTENT (IN) :: phase
  REAL (KIND=dp), DIMENSION (n_gome_pts), INTENT (IN) :: curr_wvl

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER, INTENT (OUT) :: pge_error_status

  ! ---------------
  ! Local variables
  ! ---------------
  REAL (KIND=dp), DIMENSION (2,n_gome_pts+4)  :: underspec
  REAL (KIND=dp), DIMENSION (max_spec_pts)    :: locwvl, locspec, specmod, specmod1
  REAL (KIND=dp), DIMENSION (n_gome_pts + 4)  :: tmpwav, over, under, resample, &
       resample1, subwav, tmpspec
  INTEGER :: npts, i, j, errstat, iwin, fidx, lidx, npoints

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=11), PARAMETER :: modulename = 'undersample'

  errstat = pge_errstat_ok

  IF (have_undersampling .OR. sring_fidx > 0 .OR. sring_lidx < nsol_ring) THEN
  
     ! Use derived slit from irradiance
     IF (slit_rad) THEN
        nslit = nslit_sol; slitwav = slitwav_sol; slitfit = solslitfit
     ENDIF
     
     ! ==================================================
     ! Assign solar reference spectrum to local variables
     ! ==================================================
     npts = n_refspec_pts(solar_idx)
     locwvl (1:npts) = refspec_orig_data(solar_idx,1:npts,wvl_idx)
     locspec(1:npts) = refspec_orig_data(solar_idx,1:npts,spc_idx)
     
     IF (.NOT. yn_varyslit ) THEN
        IF (which_slit == 0) THEN
           solwinfit(1:numwin, hwe_idx,1) = solwinfit(1:numwin, hwe_idx,1) + 0.01
           IF (nsl > 0) CALL gauss_multi (locwvl, locspec, specmod1, npts)
           solwinfit(1:numwin, hwe_idx,1) = solwinfit(1:numwin, hwe_idx,1) - 0.01
           CALL gauss_multi (locwvl, locspec, specmod, npts)
        ELSE IF (which_slit == 1) THEN
           solwinfit(1:numwin, hwe_idx:asy_idx,1) = solwinfit(1:numwin, hwe_idx:asy_idx,1) + 0.01
           IF (nsl > 0) CALL asym_gauss_multi (locwvl, locspec, specmod1, npts)
           solwinfit(1:numwin, hwe_idx:asy_idx,1) = solwinfit(1:numwin, hwe_idx:asy_idx,1) - 0.01
           CALL asym_gauss_multi (locwvl, locspec, specmod, npts)
        ELSE IF (which_slit == 2) THEN           
           solwinfit(1:numwin, hwl_idx:hwr_idx,1) = solwinfit(1:numwin, hwl_idx:hwr_idx,1) + 0.01
           IF (nsl > 0) CALL asym_voigt_multi (locwvl, locspec, specmod1, npts)
           solwinfit(1:numwin, hwl_idx:hwr_idx,1) = solwinfit(1:numwin, hwl_idx:hwr_idx,1) - 0.01
           CALL asym_voigt_multi (locwvl, locspec, specmod, npts)
        ELSE  IF (which_slit == 3) THEN
           solwinfit(1:numwin, hwe_idx,1) = solwinfit(1:numwin, hwe_idx,1) + 0.01
           IF (nsl > 0) CALL triangle_multi (locwvl, locspec, specmod1, npts)
           solwinfit(1:numwin, hwe_idx,1) = solwinfit(1:numwin, hwe_idx,1) - 0.01
           CALL triangle_multi (locwvl, locspec, specmod, npts)
        ELSE IF (which_slit == 4) THEN 
           solwinfit(1:numwin, hwe_idx,1) = solwinfit(1:numwin, hwe_idx,1) + 0.01
           solwinfit(1:numwin, hwn_idx,1) = solwinfit(1:numwin, hwn_idx,1) + 0.01
           IF (nsl > 0) CALL super_gauss_multi (locwvl, locspec, specmod1, npts)
           solwinfit(1:numwin, hwe_idx,1) = solwinfit(1:numwin, hwe_idx,1) - 0.01
           solwinfit(1:numwin, hwn_idx,1) = solwinfit(1:numwin, hwn_idx,1) - 0.01
           CALL super_gauss_multi (locwvl, locspec, specmod, npts)
        ELSE IF (which_slit == 5) THEN 
          SELECT CASE( instrument_idx)
          CASE (omi_idx)
            CALL omislit_multi (locwvl, locspec, specmod, npts)
          CASE (omps_idx)
           !solwinfit(1:numwin, ops_idx,1) = solwinfit(1:numwin, ops_idx,1) + 0.01
           !IF (nsl > 0 ) call ompsslit_multi (locwvl, locspec, specmod1, npts)
           ! CALL ompsslit_multi (locwvl, locspec, specmod, npts)
           !solwinfit(1:numwin, ops_idx,1) = solwinfit(1:numwin, ops_idx,1) - 0.01
          END SELECT
        ENDIF
     ELSE 
        IF (which_slit == 0) THEN
           slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, hwe_idx, 1) + 0.01
           IF (nsl > 0) CALL gauss_vary (locwvl, locspec, specmod1, npts)
           slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, hwe_idx, 1) - 0.01
           CALL gauss_vary (locwvl, locspec, specmod, npts)     
        ELSE IF (which_slit == 1) THEN
           slitfit(1:nslit, hwe_idx:asy_idx, 1) = slitfit(1:nslit, hwe_idx:asy_idx, 1) + 0.01
           IF (nsl > 0) CALL asym_gauss_vary (locwvl, locspec, specmod1, npts)
           slitfit(1:nslit, hwe_idx:asy_idx, 1) = slitfit(1:nslit, hwe_idx:asy_idx, 1) - 0.01
           CALL asym_gauss_vary (locwvl, locspec, specmod, npts)
        ELSE IF (which_slit == 2) THEN
           slitfit(1:nslit, hwl_idx:hwr_idx, 1) = slitfit(1:nslit, hwl_idx:hwr_idx, 1) + 0.01
           IF (nsl > 0) CALL asym_voigt_vary (locwvl, locspec, specmod1, npts)
           slitfit(1:nslit, hwl_idx:hwr_idx, 1) = slitfit(1:nslit, hwl_idx:hwr_idx, 1) - 0.01
           CALL asym_voigt_vary (locwvl, locspec, specmod, npts)
        ELSE IF (which_slit == 3) THEN
           slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, hwe_idx, 1) + 0.01
           IF (nsl > 0) CALL triangle_vary (locwvl, locspec, specmod1, npts)
           slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, hwe_idx, 1) - 0.01
           CALL triangle_vary (locwvl, locspec, specmod, npts)
        ELSE IF (which_slit == 4) THEN 
           slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, hwe_idx, 1) + 0.01
           slitfit(1:nslit, hwn_idx, 1) = slitfit(1:nslit, hwn_idx, 1) + 0.01
           IF (nsl > 0) CALL super_gauss_vary (locwvl, locspec, specmod1, npts)
           slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, hwe_idx, 1) - 0.01
           slitfit(1:nslit, hwn_idx, 1) = slitfit(1:nslit, hwn_idx, 1) - 0.01
           CALL super_gauss_vary (locwvl, locspec, specmod, npts)     
        ELSE IF (which_slit == 5) THEN 
          SELECT CASE( instrument_idx)
          CASE (omi_idx)
            CALL omislit_vary (locwvl, locspec, specmod, npts)
          CASE (omps_idx)
            !slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, ops_idx, 1) + 0.01
            !IF (nsl > 0) CALL triangle_vary (locwvl, locspec, specmod1, npts)
            !slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, ops_idx, 1) - 0.01
            !CALL ompsslit_vary (locwvl, locspec, specmod, npts)
          END SELECT
        ENDIF
     ENDIF
     
     ! Append Ring Source Spectrum
     IF (sring_fidx > 0 .OR. sring_lidx < nsol_ring) THEN
        CALL append_solring(nsol_ring, sring_fidx, sring_lidx, sol_spec_ring(wvl_idx, 1:nsol_ring), &
             sol_spec_ring(spc_idx, 1:nsol_ring), npts, locwvl(1:npts), specmod(1:npts), pge_error_status)
        IF (pge_error_status > pge_errstat_warning) RETURN
     ENDIF
  ENDIF
 
  IF (.NOT. have_undersampling) RETURN
  
  ! do it separately for each window
  lidx = 0
  DO iwin = 1, numwin
     npoints = nradpix(iwin) + refnhextra*2 + 4 ! reference has 4 more wavelengths and add 4 extra pixels
     fidx =  lidx + 1            ! fidx:lidx refers to position in curr_wvl
     lidx =  fidx + npoints - 5

     !IF (iwin == 2) THEN
     !   database(us1_idx, fidx:lidx) = 0.0;  database(us2_idx, fidx:lidx) = 0.0; CYCLE
     !ENDIF

     subwav(3:npoints-2) = curr_wvl(fidx:lidx)

     ! Add extra wavelengths
     subwav(2) = 2 * subwav(3) - subwav(4); subwav(1) = 2 * subwav(2) - subwav(3)
     subwav(npoints-1) = 2 * subwav(npoints-2) - subwav(npoints-3)
     subwav(npoints)   = 2 * subwav(npoints-1) - subwav(npoints-2)
 
     ! Phase1 calculation: Calculate spline derivatives for KPNO data
     !                     Calculate solar spectrum at OMI positions 
     CALL interpolation (npts, locwvl(1:npts), specmod(1:npts), npoints, &
          subwav(1:npoints), resample(1:npoints),  errstat)

     IF ( errstat > pge_errstat_warning ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; 
        pge_error_status = pge_errstat_error; RETURN
     END IF
     
     ! weighting function for slit width
     IF (nsl > 0) THEN
        CALL interpolation (npts, locwvl(1:npts), specmod1(1:npts), npoints, &
             subwav(1:npoints),  resample1(1:npoints), errstat)
        
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; 
           pge_error_status = pge_errstat_error; RETURN
        END IF
        database_shiwf(1,fidx:lidx)=(resample1(3:npoints-2)- resample(3:npoints-2)) / 0.01
     ENDIF
     
     ! Calculate solar spectrum at OMI + phase positions, original and resampled.
     ! ------------------------------------------------------------------------------
     ! The original ("modified K.C.) scheme to compute the UNDERSPEC wavelength array
     ! ------------------------------------------------------------------------------
     ! ( assumes ABS(PHASE) < 1.0 )
     ! ----------------------------
     tmpwav(1:npoints-1) = (1.0-phase) * subwav(1:npoints-1) + phase * subwav(2:npoints)

     tmpwav(1) = subwav(1);  tmpwav(npoints)   = subwav(npoints)
     IF ( tmpwav(2) <= tmpwav(1) ) tmpwav(2) = (tmpwav(1)+tmpwav(3))/2.0
    
     CALL interpolation (  npts, locwvl(1:npts), specmod(1:npts), npoints, &
          tmpwav(1:npoints), over(1:npoints), errstat )

     IF ( errstat /= pge_errstat_ok ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0)
        pge_error_status = pge_errstat_error; RETURN
     END IF
     
     CALL interpolation (npoints, subwav(1:npoints), resample(1:npoints),&
          npoints, tmpwav(1:npoints), under(1:npoints), errstat )
    
     IF ( errstat /= pge_errstat_ok ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0)
        pge_error_status = pge_errstat_error; RETURN
     END IF
   
     underspec(1,1:npoints) = over(1:npoints) - under(1:npoints)
     print * , underspec(1,1:npoints)
     ! --------------------------------------------------------------
     ! Phase2 calculation: Calculate solar spectrum at OMI positions, 
     ! original and resampled.
     ! -------------------------------------------------------------- 
     tmpspec(1:npoints)   = resample(1:npoints)
     resample (1:npoints) = over(1:npoints)
     over(1:npoints)      = tmpspec(1:npoints)
     
     CALL interpolation (npoints, tmpwav(1:npoints), &
          resample(1:npoints), npoints, subwav(1:npoints), &
          under(1:npoints), errstat )
     IF ( errstat /= pge_errstat_ok ) THEN    
        errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0)
        pge_error_status = pge_errstat_error; RETURN
     END IF
  
     ! ====================================
     ! Compute final undersampling spectrum
     ! ====================================
     underspec(2,1:npoints) = over(1:npoints) - under(1:npoints)
  
     ! check for missing F and thus spikes
     !dwavmax = subwav(2) - subwav(1)
     !DO i = 2, npoints
     !   IF (subwav(i) - subwav(i-1) > dwavmax * 1.2) THEN
     !      underspec(1, i-1:i) = 0.0
     !      underspec(2, i-1:i) = 0.0
     !   ENDIF
     !ENDDO
   
     refspec_orig_data(us1_idx, fidx:lidx, wvl_idx) = tmpwav   (   3:npoints-2)
     refspec_orig_data(us1_idx, fidx:lidx, spc_idx) = underspec(1, 3:npoints-2)
     refspec_orig_data(us2_idx, fidx:lidx, wvl_idx) = subwav   (   3:npoints-2)
     refspec_orig_data(us2_idx, fidx:lidx, spc_idx) = underspec(2, 3:npoints-2)
     
     database(us1_idx, fidx:lidx) = underspec(1, 3:npoints-2)
     database(us2_idx, fidx:lidx) = underspec(2, 3:npoints-2)    

  ENDDO ! end window loop

  n_refspec_pts (us1_idx)  = n_gome_pts
  n_refspec_pts (us2_idx)  = n_gome_pts

  !WRITE(90, *) n_gome_pts, nradpix(1:numwin) + 4
  !DO i = 1, n_gome_pts
  !   WRITE(91, *) curr_wvl(i), database(us1_idx, i),  database(us2_idx, i)
  !ENDDO
  !STOP

  ! Use derived slit from radiance later on if it is available
  IF (slit_rad) THEN
     nslit = nslit_rad; slitwav = slitwav_rad; slitfit = radslitfit
  ENDIF

  RETURN
END SUBROUTINE undersample



