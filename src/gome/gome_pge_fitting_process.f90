  ! *********************** Modification History **************
  ! xliu: 
  ! 1. Call get_reg_matrix to setup regularization matrix 
  !    for ozone profile retrieval
  ! 2. Calculate average fitting statistic for ozone
  !    profile retrieval
  ! 3. Call solar_fit_vary if yn_varyslit is set
  ! 4. ***** to be changed *****
  !    For radiance spectra where no calibration is performed
  !    and variable slit width is used, then it is better to 
  !    update wavelength position according to shiarr, sswav, ls

  !    and squarr
  ! **********************************************************

SUBROUTINE gome_pge_fitting_process (l1funit, l2funit, l1_inputs_fname_sol,  &
     l1_inputs_fname_rad, l2_output_fname, l2_hdf_flag, pge_error_status )

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: max_rs_idx, refspec_strings,           &
       wvl_idx, spc_idx, sig_idx, gome_data_input_lun, gome_data_output_lun
  USE OMSAO_parameters_module, ONLY: maxchlen, max_fit_pts
  USE OMSAO_variables_module,  ONLY: refspec_fname,       &
       linenum_lim, pixnum_lim, radwavcal_freq, curr_rad_spec,         &
       curr_sol_spec, n_irrad_wvl, sza_atm, vza_atm, fitvar_rad, fitvar_sol, &
       lo_sunbnd, up_sunbnd,  phase, n_rad_wvl, slit_rad, slit_redo,         &
       wavcal_redo, amf, amfgeo, sol_zen_eff, yn_varyslit,      &
       wavcal_sol, do_bandavg, npix_fitted, n_radwvl_sav, radwvl_sav,        &
       npix_fitting, numwin, band_selectors, nradpix, nradpix_sav, the_lat, &
       the_lon, the_lons, the_lats, refidx, refwvl, ozabs_convl, so2crs_convl,&
       instrument_idx, gome_idx, scia_idx, gome2_idx
  USE OMSAO_gome_data_module, ONLY: &
       lm_gome_solspec, lm_gome_eshine, lm_gome_groundpix, n_gome_max_pts,   &
       lm_gome_chan_1, lm_gome_chan_2, lm_gome_chan_3, lm_gome_chan_4,       &
       lm_gome_band_1, lm_gome_band_2a, lm_gome_band_2b, lm_gome_band_3,     &
       lm_gome_band_4, lat_idx, lon_idx, gome_spec_missing,                  &
       gome_data_input_fname, gome_data_output_fname, n_gome_data_dim,       &
       n_gome_ang, n_gome_geo, gome_radspec, gome_geoloc, gome_curpix,       &
       gome_curscan, n_gome_radpts, gome_npix, gome_stpix, gome_endpix,      &
       gome_angles_wrtn
  USE ozprof_data_module, ONLY : ozprof_flag, lcurve_write, ozwrtint,        &
       lcurve_fname, ozwrtint_fname, lcurve_unit, ozwrtint_unit,             &
       div_rad, div_sun, algorithm_name, algorithm_version, fullorb, do_ch2reso, &
       nsaa_spike, radcalwrt, do_simu, calunit
  USE OMSAO_errstat_module

  ! GOME-2
  USE BOREAS_gome2_data_module, ONLY: &
       n_gome2_scans, gome2_curr_nxtrack, gome2_utc_year, cloud, &
       gome2_cloud_press, gome2_cloud_frac, n_xtrack_per_record
  USE file_handling, ONLY: open_file, fileptr
  USE gome_mdr_1b_earthshine_module
  USE mphr_module


  IMPLICIT NONE

  ! -------------------------
  ! Name of module/subroutine
  ! -------------------------
  CHARACTER (LEN=24), PARAMETER :: modulename = 'gome_pge_fitting_process'

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,           INTENT (IN) :: l1funit, l2funit
  CHARACTER (LEN=*), INTENT (IN) :: l1_inputs_fname_sol,l1_inputs_fname_rad, &
       l2_output_fname

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER, INTENT (INOUT) :: pge_error_status, l2_hdf_flag

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER :: funit, j, i, ical, exval, initval, npix, nn_rad_wvl, fidx, lidx
  LOGICAL :: error
  CHARACTER (LEN=maxchlen) :: lm_gome_chan, lm_gome_band
  CHARACTER (LEN=4)        :: endpix_str
  INTEGER                  :: bndc, numwav_sav

  ! ------------------------
  ! GOME-2 related variables
  ! ------------------------
  TYPE(fileptr) :: fptr
  INTEGER       :: status
  INTEGER       :: gome2_chan, gome2_band
  INTEGER       :: iscan 
  INTEGER       :: ixtrack
  INTEGER       :: n_total_pixels
  INTEGER       :: n_gome2_xtrack

  ! ======================================
  ! Variable for I/O status of file access
  ! ======================================
  INTEGER :: file_read_stat
  REAL (KIND=dp), DIMENSION(3)    :: fitcol
  REAL (KIND=dp), DIMENSION(3, 2) :: dfitcol
  REAL (KIND=dp) :: fitcol_avg, rms_avg, dfitcol_avg, drel_fitcol_avg, &
       rms, drel_fitcol, tmpwav, tmpcorr
  INTEGER        :: n_fitted_col, ifitpix, solfit_exval

  ! -------------------------
  ! OMI L1B related variables
  ! -------------------------
  INTEGER :: errstat, OMI_SMF_setmsg
  error = .FALSE.
 

  
  ! -------------------------------------------------------
  ! Determine if selected bands are valid
  ! -------------------------------------------------------
  DO i = 1, numwin	
     bndc = band_selectors(i)
     IF (bndc /= 1 .AND. bndc /= 2 .AND. bndc /= 3 .AND. bndc /= 4) THEN
        WRITE (*, '(A, I2)') 'ERROR...unknown GOME band selected: ', bndc
        pge_error_status = pge_errstat_error ; RETURN   
     ENDIF
  ENDDO

  ! -------------------------------------------------------------------
  ! Open GOME solar data file, read GOME solar spectrum, and close file
  ! -------------------------------------------------------------------

  IF ( (instrument_idx == gome_idx) .OR. (instrument_idx == scia_idx) ) THEN
     OPEN (UNIT=l1funit, FILE=TRIM(ADJUSTL(l1_inputs_fname_sol)), &
          STATUS='OLD', IOSTAT=errstat )
     IF ( errstat /= pge_errstat_ok ) THEN
        pge_error_status = pge_errstat_error; RETURN
     END IF
     IF (instrument_idx == gome_idx) THEN
        CALL gome_read_el1data_sol (l1funit, file_read_stat, error )
     ELSE
        CALL scia_read_el1data_sol (l1funit, file_read_stat, error )
     END IF  
     CLOSE ( l1funit )
  ELSE IF (instrument_idx == gome2_idx) THEN
     fptr%filename = TRIM(ADJUSTL(l1_inputs_fname_sol))
     CALL open_file(fptr, status)
     CALL read_metadata(fptr, status)
     gome2_utc_year = sensing_start_year  !Read from metadata
     CALL gome2_read_l1b_sol(fptr, file_read_stat, error) 
  END IF


  ! -----------------------------------------------------------
  ! Perform solar wavelength calibration and slit width fitting
  ! -----------------------------------------------------------
  CALL adj_solar_data ( error)


  PRINT *, 'Performing solar wavelength calibration'
  IF (yn_varyslit) THEN
     CALL solar_fit_vary (calunit, error )
     IF (wavcal_sol) CALL solar_wavcal_vary(calunit, error)
  ELSE 
     CALL solar_fit ( error )  !xliu
  END IF
  IF (error) RETURN



  ifitpix = 0; npix_fitting=0     !number of pixels (failure + success)
  n_fitted_col = 0; npix_fitted=0 !number of successfully fitted pixels 
   
  fitcol_avg = 0.0; rms_avg = 0.0; dfitcol_avg = 0.0; drel_fitcol_avg = 0.0



  ! ---------------------------------------------------------------
  ! Open GOME radiance data file; reading is done inside pixel loop
  ! ---------------------------------------------------------------
  IF ( (instrument_idx == gome_idx) .OR. (instrument_idx == scia_idx) ) THEN
     OPEN (UNIT=l1funit, FILE=TRIM(ADJUSTL(l1_inputs_fname_rad)), &
          STATUS='OLD', IOSTAT=errstat)
     IF ( errstat /= pge_errstat_ok ) THEN
        WRITE(*, *) modulename, ': Cannot open radiance file!!!'
        pge_error_status = pge_errstat_error; RETURN
     END IF
  ELSE IF (instrument_idx == gome2_idx) THEN
     fptr%filename = TRIM(ADJUSTL(l1_inputs_fname_rad))
     CALL open_file(fptr, status)
     CALL read_metadata(fptr, status)
     n_gome2_scans = n_earthview_records(fptr)
     CALL which_gome2_bands(fptr)
     CALL gome2_count_pixels( n_total_pixels, n_gome2_scans, fptr )
     iscan = linenum_lim(1)
     ixtrack = 0
  END IF

  ! ---------------------------------------------
  ! Open L2 output file for slant column fitting
  ! ---------------------------------------------
  IF ( (l2_hdf_flag == 0) .OR. (instrument_idx /= gome2_idx)) THEN
     IF (.NOT. ozprof_flag) THEN
        OPEN (UNIT=l2funit, FILE=TRIM(ADJUSTL(l2_output_fname)//'.out'), &
             STATUS='UNKNOWN', IOSTAT=errstat)
        IF ( errstat /= pge_errstat_ok ) THEN
           WRITE(*, *) modulename, ': Cannot open output file!!!'
           pge_error_status = pge_errstat_error; RETURN
        END IF
     ENDIF

     IF (ozprof_flag) THEN
        IF (lcurve_write) OPEN(UNIT=lcurve_unit, &
             file=TRIM(ADJUSTL(lcurve_fname)), status='unknown')
        IF (ozwrtint)     OPEN(UNIT=ozwrtint_unit, &
             file=TRIM(ADJUSTL(ozwrtint_fname)), status='unknown')
     END IF
  ELSE  
     ! If writing to HDF file, initialize output variables
     IF (.NOT. ozprof_flag) THEN
        !CALL gome2_init_output(n_gome2_scans, n_total_pixels, l2_output_fname)
     ELSE
        IF (l2_hdf_flag >= 1) CALL gome2_init_output_profile(n_gome2_scans, &
             n_total_pixels, l2_output_fname, l2_hdf_flag)
     ENDIF
  ENDIF



  ifitpix = 0
  file_read_stat = file_read_ok
  gome_groundpix: DO WHILE ( file_read_stat /= file_read_eof )
 
     ! Need to convolve high-resolution ozone absorption cross section (for this position)
     ! Once the xsection is convolved, it will be set to false in ROUTINE getabs_crs
     ozabs_convl = .TRUE.; so2crs_convl = .TRUE.

     ! ------------------
     ! Read radiance data
     ! ------------------
     SELECT CASE (instrument_idx)
     CASE (gome_idx)
        IF (fullorb .AND. ozprof_flag .AND. do_ch2reso) THEN
           CALL gome_read_el1data_allradch1 (l1funit, file_read_stat, error)
        ELSE IF (fullorb .AND. ozprof_flag) THEN
           CALL gome_read_el1data_radavg (l1funit, file_read_stat, error )  
        ELSE     
           CALL gome_read_el1data_rad (l1funit, file_read_stat, error )
        ENDIF
     CASE (scia_idx)
        CALL scia_read_el1data_rad (l1funit, file_read_stat, error)  
     CASE (gome2_idx)
        ixtrack = ixtrack + 1
        gome_curpix = ifitpix + 1
        gome2_curr_nxtrack = n_gome2_xtrack(iscan, fptr)
        CALL gome2_read_l1b_rad (fptr, iscan, ixtrack, file_read_stat, error)
        IF (file_read_stat == file_read_eof) EXIT gome_groundpix
     END SELECT


     ! -----------------------------------------------------------------------
     ! Cases of "missing" data. This includes any restrictions imposed on the 
     ! kind of GOME pixel (not) to process.
     ! -----------------------------------------------------------------------
     SELECT CASE (instrument_idx)
     CASE (gome2_idx)  !Linenum_lim is by scan # (usually 32 pixels/scan)
        ! Band GOME-2 angles show up as -2147.48
        IF ( (file_read_stat == file_read_ok)                           .AND. &
             ( (n_gome_radpts .LE. 0 )                                  .OR.  &
             (SUM(gome_radspec(n_gome_data_dim,1:n_gome_radpts)) > 0.0) .OR.  &
             (iscan < linenum_lim(1))          )                        .OR.  &
             ANY ( gome_geoloc <= -2000.0_sp )                          .OR.  &
             ANY ( gome_angles_wrtn <=  0.0_sp ) )  &
             file_read_stat = file_read_missing
     CASE DEFAULT
        IF ( (file_read_stat == file_read_ok)                           .AND. &
             ( (n_gome_radpts == 0)                                     .OR.  &
             (SUM(gome_radspec(n_gome_data_dim,1:n_gome_radpts)) > 0.0) .OR.  &
             (gome_curpix < linenum_lim(1))          ) )                    &
             file_read_stat = file_read_missing
     END SELECT


     ! ---------------------------------------------------------------------------------
     ! End Of File can also mean having read past the last GOME pixel number to process.
     ! ---------------------------------------------------------------------------------
     SELECT CASE (instrument_idx)
     CASE (gome2_idx)
        IF ( iscan > linenum_lim(2) ) file_read_stat = file_read_eof
     CASE DEFAULT
        IF ( gome_curpix > linenum_lim(2) ) file_read_stat = file_read_eof
     END SELECT

     SELECT CASE ( file_read_stat )
     CASE ( file_read_eof, file_read_failed ) 
        EXIT gome_groundpix
     CASE ( file_read_ok )
      

        !===================================================
        ! Open level 2 out file for ozone profile retrieval
        !===================================================
        IF ( (l2_hdf_flag == 0) .OR. (instrument_idx /= gome2_idx)) THEN
           IF (ozprof_flag) THEN          
              IF (lcurve_write) WRITE(lcurve_unit, '(A8,2I4,2F10.2)') &
                   'Pixels: ', gome_stpix, gome_endpix, gome_geoloc(lat_idx, 5), gome_geoloc(lon_idx, 5)
              IF (ozwrtint)     WRITE(ozwrtint_unit, '(A20,I4,A13,I4,A7)') &
                   '======Start Pixel = ', gome_stpix, ' End Pixel = ', gome_endpix, ' ======'
           ENDIF
        ENDIF


        ! ---------------------
        ! Convert angles to TOA
        ! ---------------------
        CALL adj_earthshine_data (error)

        ifitpix = ifitpix + 1


        ! South Atlantic Anomaly region: detect spikes        
        IF ((the_lat > -50.0 .AND. the_lat < 5.0 .AND. the_lon > -90.0 .AND. the_lon < 0.0) & 
             .OR. (the_lat > -35.0 .AND. the_lat < -10.0 .AND. the_lon > 0.0 .AND. the_lon < 35.0)) THEN
             CALL ROUGH_SPIKE_DETECT(n_rad_wvl, curr_rad_spec(1, 1:n_rad_wvl), &
             curr_rad_spec(2, 1:n_rad_wvl), curr_sol_spec(2, 1:n_rad_wvl), nsaa_spike)
        ELSE
           nsaa_spike = 0
        ENDIF

        ! Obtain measurement error in term sun-normalized radiance
        IF (ozprof_flag) THEN
           CALL adj_rad_sig (curr_rad_spec(wvl_idx:sig_idx, 1:n_rad_wvl), &
                curr_sol_spec(wvl_idx:sig_idx,1:n_rad_wvl), n_rad_wvl)
           
           !CALL reduce_resolution(curr_rad_spec(wvl_idx:sig_idx, 1:n_rad_wvl), &
           !curr_sol_spec(wvl_idx:sig_idx,1:n_rad_wvl), n_rad_wvl, nn_rad_wvl)
           !n_rad_wvl = nn_rad_wvl  
           !nradpix(1) = nn_rad_wvl
           !n_irrad_wvl = n_rad_wvl
        ENDIF

           
        ! save the original grids for later obtain ozone cross section
        n_radwvl_sav = n_rad_wvl; nradpix_sav = nradpix
        radwvl_sav(1:n_rad_wvl) = curr_rad_spec(wvl_idx, 1:n_rad_wvl)   
      

        ! average and subsampling for selected bands
        ! update nradpix
        IF (ozprof_flag) THEN
           IF (do_bandavg) THEN 
              !WRITE(90, *) n_rad_wvl
              !DO i = 1, n_rad_wvl
              !   WRITE(90, '(F10.4, D14.6)') curr_rad_spec(wvl_idx:spc_idx, i)
              !ENDDO
              CALL avg_band_radspec (curr_rad_spec(wvl_idx:sig_idx,1:n_rad_wvl), &
                   n_rad_wvl, pge_error_status )
              !WRITE(91, *) n_rad_wvl
              !DO i = 1, n_rad_wvl
              !   WRITE(91, '(F10.4, D14.6)') curr_rad_spec(wvl_idx:spc_idx, i)
              !ENDDO
              !IF ( pge_error_status >= pge_errstat_error ) RETURN
           ENDIF            
        ENDIF

!!$        fidx = 1
!!$        DO i = 1, numwin
!!$           lidx = fidx + nradpix(i) - 1
!!$           refidx(fidx:lidx) = (/(j, j = fidx + 2 + (i - 1) * 4, lidx + 2 + (i - 1) * 4)/)
!!$           fidx = lidx  + 1
!!$        ENDDO



        ! ----------------------------------------
        ! Write geolocation if writing to HDF file
        ! ----------------------------------------
        IF ( (instrument_idx == gome2_idx) .AND. (l2_hdf_flag /= 0) ) THEN
           CALL gome2_geoloc_output_profile(iscan)
        ENDIF

        ! --------------------------
        ! Fit radiance spectra
        ! --------------------------           
        PRINT *, 'Fitting GOME pixel number ', gome_curpix, 'Scan', iscan, 'Xtrack', ixtrack
        
        
        IF (ozprof_flag) THEN
           CALL specfit_ozprof (initval, fitcol, dfitcol, rms, exval)
        ELSE
           CALL radiance_fit (fitcol(1), dfitcol(1, 1), rms, exval) 
        ENDIF

        ! Write ascii file or write to HDF

        IF ( ( exval > 0 .AND. fitcol(1) > 0.0 .AND. dfitcol(1, 1) >= 0.0 ) .OR. &
             (radcalwrt .AND. do_simu) )THEN

           IF ( (l2_hdf_flag == 0) .OR. (instrument_idx /= gome2_idx)) THEN
              IF (ozprof_flag) THEN
                 WRITE(endpix_str, '(I5.5)') gome_endpix
                 IF (gome_npix == 1) THEN 
                    OPEN (UNIT=l2funit, FILE=TRIM(ADJUSTL(l2_output_fname)) // '_' // endpix_str &
                         // '.out', STATUS='UNKNOWN', IOSTAT=errstat)
                 ELSE
                    OPEN (UNIT=l2funit, FILE=TRIM(ADJUSTL(l2_output_fname)) // '_' // endpix_str &
                         // '_0avg.out', STATUS='UNKNOWN', IOSTAT=errstat)
                 ENDIF
                 IF ( errstat /= pge_errstat_ok ) THEN
                    WRITE(*, *) modulename, ': Cannot open output file!!!'
                    pge_error_status = pge_errstat_error; RETURN
                 END IF
                 WRITE(l2funit, '(A10,A2,A11)') algorithm_name, ', ', algorithm_version 
              ENDIF

              CALL gome_write_intermed (l2funit, gome_curpix, gome_curscan, fitcol,&
                   dfitcol,  rms, amf, amfgeo, sol_zen_eff, n_gome_ang, &
                   sza_atm(1:n_gome_ang), vza_atm(1:n_gome_ang), n_gome_geo, &
                   gome_geoloc(lat_idx,1:n_gome_geo), &
                   gome_geoloc(lon_idx,1:n_gome_geo), exval)
              IF (ozprof_flag) CLOSE ( l2funit )

           ELSE  ! Write to HDF file

              IF (.NOT. ozprof_flag) THEN
                 CALL gome2_column_output(fitcol, dfitcol, rms, amfgeo, amf)
              ELSE
                 CALL gome2_output_profile(gome_curpix, fitcol, dfitcol, &
                      amfgeo, exval, rms, l2_hdf_flag)
              ENDIF

           ENDIF

           
           ! ----------------------------------------------------------------------
           ! Some general statistics on the average fitted column and uncertainty.
           ! Again, we make sure that only "good" fits are included in the average.
           ! ----------------------------------------------------------------------
           
           fitcol_avg      = fitcol_avg + fitcol(1)
           rms_avg         = rms_avg + rms
           dfitcol_avg     = dfitcol_avg + dfitcol(1, 1)
           drel_fitcol_avg = drel_fitcol_avg + dfitcol(1, 1) / fitcol(1)
           n_fitted_col    = n_fitted_col + 1; npix_fitted = n_fitted_col
        END IF


        
     CASE DEFAULT
        ! -----------------------
        ! Nothing to be done here
        ! -----------------------
     END SELECT

     ! If come to end of cross-track read in GOME-2, reset counter for next scan
     IF ( (instrument_idx .EQ. gome2_idx) .AND. &
          ( (ixtrack .GE. gome2_curr_nxtrack) .OR. &
          (ixtrack .GE. pixnum_lim(2) ) ) ) THEN
        ixtrack = 0
        iscan = iscan + 1
     END IF

     npix_fitting = ifitpix

  END DO gome_groundpix

  ! ------------
  ! Final output
  ! ------------
  IF ( n_fitted_col == 0 ) n_fitted_col = 1
  
  CALL write_final(fitcol_avg, rms_avg, dfitcol_avg, drel_fitcol_avg, &
       n_fitted_col)
  WRITE(*, '(1x,A,I5,A,I5)') 'Number of pixels = ', ifitpix, &
       '      Number of fitted pixels = ', npix_fitting


  ! -----------------------------------------
  ! Close L1 radiance file and L2 output file
  ! -----------------------------------------
  CLOSE ( l1funit )
  IF ( (l2_hdf_flag == 0) .OR. (instrument_idx /= gome2_idx) ) THEN
     IF (.NOT. ozprof_flag) CLOSE ( l2funit )
     IF(ozprof_flag) THEN
        IF (lcurve_write) CLOSE (lcurve_unit)
        IF (ozwrtint)     CLOSE (ozwrtint_unit)
     END IF
  ELSE
     CALL gome2_close_output
  ENDIF


  RETURN
END SUBROUTINE gome_pge_fitting_process
