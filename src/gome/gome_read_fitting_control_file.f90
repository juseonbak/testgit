  ! *********************** Modification History **********************
  ! Xiong Liu; July, 2003  (!xliu)
  ! 1. Read whether to do ozone profile retrieval and set a flag 
  ! 2. Read additional fitting control variables if ozprof_flag is set
  ! 3. Read an option about variable slit width, add yn_varyslit variable 
  !    n_slit_itnerval, slit_fname, slit_redo, wavcal_redo, wavcal_fname, 
  !    in USE OMSAO_variables_module
  ! 4. Add 1 nm more for winwav_min, winwav_max to avoid interpolation
  !    out of bounds for solar spectrum calibration
  ! 5. Read option use_meas_sig
  ! 6. Read option use_pixel_bin
  ! *******************************************************************

SUBROUTINE gome_read_fitting_control_file (fit_ctrl_unit, fit_ctrl_file, &
     instrument_idx, l1_inputs_fname_sol, l1_inputs_fname_rad, l2_output_fname, &
     pge_error_status)     

  ! ***********************************************************
  !
  !   Read fitting control parameters from input control file
  !
  ! ***********************************************************

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: &
       max_rs_idx, calfit_strings, max_calfit_idx, radfit_strings, &
       mns_idx, mxs_idx, ad1_idx, bro_idx, lbe_idx, n_max_fitpars, refspec_strings, &
       genline_str, socline_str, racline_str, rafline_str, rspline_str,  &
       iofline_str, molline_str, eoi3str, solar_idx, amf_idx, us1_idx, us2_idx
  USE OMSAO_parameters_module,   ONLY: &
       maxchlen, forever, n_sol_winwav, n_rad_winwav, max_mol_fit, &
       vb_lev_omidebug, vb_lev_develop, maxwin
  USE OMSAO_variables_module,    ONLY: &
       fitcol_idx, fincol_idx, n_fincol_idx, n_mol_fit, max_itnum_sol, max_itnum_rad,  &
       yn_smooth, yn_doas, weight_sun, fitvar_sol_saved, fitvar_sol_init,              &
       fitvar_rad_init, fitvar_rad_saved, yn_varyslit, which_slit, n_slit_step,         &
       slit_fname, slit_redo, wavcal_redo, wavcal_fname, swavcal_fname, use_meas_sig,  &
       smooth_slit, slit_fit_pts, wavcal_fit_pts, n_wavcal_step,                       &
       wavcal_sol, mask_fitvar_rad, mask_fitvar_sol, n_fitvar_sol, renorm, weight_rad, &
       szamax, zatmos, slit_rad, rslit_fname,  n_fitvar_rad, lo_sunbnd, up_sunbnd,     &
       lo_sunbnd_init,up_sunbnd_init, lo_radbnd, up_radbnd, n_refspec, fitctrl_fname,  &
       refspec_fname, fitpar_idxname, radwavcal_freq, tol,  epsrel,  epsabs,  epsx,    &
       pm_one, phase, pixnum_lim,  static_input_fnames, &
       fitvar_rad_str, verb_thresh_lev, winwav_min, winwav_max, have_amftable,         &
       have_undersampling, winpix, sol_identifier, rad_identifier, numwin,             &
       band_selectors, do_bandavg, n_band_avg, n_band_samp, outdir, atmdbdir, refdbdir,&
       gome_idx, which_instrument, max_instrument_idx
  USE OMSAO_gome_data_module, ONLY: lm_gome_solspec, &
       lm_gome_eshine, n_gome_max_pts, n_gome_data_dim, gome_spec_missing, gome_orbc
  USE OMSAO_errstat_module

  !xliu
  USE ozprof_data_module,     ONLY: ozprof_str, ozprof_flag, ozprof_input_fname,      &
       fullorb, do_ch2reso, b1ab_change, colprof, l1l2inp_unit

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,           INTENT (IN) :: fit_ctrl_unit
  CHARACTER (LEN=*), INTENT (IN) :: fit_ctrl_file

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER,           INTENT (INOUT) :: instrument_idx
  INTEGER,           INTENT (OUT)   :: pge_error_status
  CHARACTER (LEN=*), INTENT (OUT)   :: l1_inputs_fname_sol, l1_inputs_fname_rad, &
       l2_output_fname

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER                  :: i, j, k, file_read_stat, sidx, ridx, idx, nidx
  CHARACTER (LEN=maxchlen) :: tmpchar, l1l2_files
  CHARACTER (LEN=3)        :: idxchar
  CHARACTER (LEN=5)        :: idxchar1
  LOGICAL                  :: yn_eoi, rw_l1l2_here
  REAL      (KIND=dp)      :: vartmp, lotmp, uptmp

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=30), PARAMETER :: modulename = 'gome_read_fitting_control_file'

  ! ========================
  ! Error handling variables
  ! ========================
  INTEGER :: errstat, version, ios

  ! -----------------------
  ! GOME specific additions
  ! -----------------------
  CHARACTER (LEN=25), PARAMETER :: lm_instrument = 'Satellite instrument name'
  CHARACTER (LEN=19), PARAMETER :: lm_l1inputs   = 'Level 1 input files'
  CHARACTER (LEN=19), PARAMETER :: lm_l2output   = 'Level 2 output file'
  CHARACTER (LEN=21), PARAMETER :: lm_amftable   = 'Air mass factor table'
  CHARACTER (LEN=30), PARAMETER :: lm_atmdb      = 'Atmospheric database directory'
  CHARACTER (LEN=27), PARAMETER :: lm_refdb      = 'Reference spectra directory'
  CHARACTER (LEN=30), PARAMETER :: lm_bandselect = 'GOME radiance bands to be used'

  ! =================================
  ! External OMI and Toolkit routines
  ! =================================
  INTEGER :: OMI_SMF_setmsg

  ! -----------------------------------------------------------
  ! Initialize array with reference spectrum names to Zero_Spec
  ! -----------------------------------------------------------
  DO j = solar_idx, max_rs_idx
     refspec_fname(j) = 'OMSAO_Zero_Spec.dat'
  END DO

  ! -------------------------
  ! Open fitting control file
  ! -------------------------

  OPEN ( UNIT=fit_ctrl_unit, FILE=TRIM(ADJUSTL(fit_ctrl_file)), &
       STATUS='OLD', IOSTAT=errstat)
           print * , fit_ctrl_file
  IF ( errstat /= pge_errstat_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
          TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF

  ! -----------------------------------------------
  ! Position cursor to read instrument name
  ! -----------------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, lm_instrument, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_instrument, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF
  READ (fit_ctrl_unit, '(A)') tmpchar
  CALL string2index ( which_instrument, max_instrument_idx, tmpchar, instrument_idx )

  ! --------------------------
  ! For now let's do GOME only
  ! --------------------------
  IF ( instrument_idx /= gome_idx ) THEN
     OPEN ( UNIT=fit_ctrl_unit )  ;  RETURN
  END IF

  ! -------------------------------------------
  ! Position cursor to read Level 1 input files
  ! -------------------------------------------
  REWIND ( fit_ctrl_unit )
  rw_l1l2_here = .FALSE.             ! read l1 (write l2) fname from another file
  CALL skip_to_filemark ( fit_ctrl_unit, lm_l1inputs, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_l2output, modulename, 0)
     pge_error_status = pge_errstat_error ; RETURN
  ELSE
     READ (fit_ctrl_unit, *)        rw_l1l2_here, fullorb, do_ch2reso
     READ (fit_ctrl_unit, '(A)')    l1l2_files
     IF (rw_l1l2_here) THEN
        READ (fit_ctrl_unit, '(A)') l1_inputs_fname_sol
        READ (fit_ctrl_unit, '(A)') l1_inputs_fname_rad
     END IF
  END IF

  ! -------------------------------------------
  ! Position cursor to read Level 2 output file
  ! -------------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, lm_l2output, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_w_read_fitctrl_file, lm_l2output, modulename, 0)
     pge_error_status = pge_errstat_warning
     l2_output_fname = '../out/L2-default-output.dat'
  ELSE
     IF (rw_l1l2_here) READ (fit_ctrl_unit, '(A)') l2_output_fname
  END IF

  IF (.NOT. rw_l1l2_here) THEN
     OPEN(UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(l1l2_files)), STATUS='OLD', IOSTAT=errstat)
     IF ( errstat /= pge_errstat_ok ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
             TRIM(ADJUSTL(l1l2_files)), modulename, 0)
        pge_error_status = pge_errstat_error; RETURN
     ELSE
        READ(l1l2inp_unit, '(A)') l1_inputs_fname_sol
        READ(l1l2inp_unit, '(A)') l1_inputs_fname_rad
        READ(l1l2inp_unit, '(A)') l2_output_fname  
        !READ(l1l2inp_unit, *) colprof  ! read collocated sage/ozonesond profiles
        CLOSE(UNIT=l1l2inp_unit) 
     END IF
  END IF
  
  i = INDEX(l1_inputs_fname_sol, 'lv1_') + 11
  gome_orbc = l1_inputs_fname_rad(i-2:i)

  !xliu, 09/23/05 Add direcotry, remove hard code directory
  ! ----------------------------------------------------------
  ! Position cursor to read database directory
  ! ----------------------------------------------------------  
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, lm_atmdb, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_atmdb, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  ELSE
     READ (fit_ctrl_unit, '(A)') atmdbdir
  ENDIF

  REWIND ( fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, lm_refdb, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_refdb, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  ELSE
     READ (fit_ctrl_unit, '(A)') refdbdir
  ENDIF

  ! --------------------------------------
  ! Position cursor to read AMF table file
  ! --------------------------------------
  REWIND ( fit_ctrl_unit )
  have_amftable = .FALSE.
  CALL skip_to_filemark ( fit_ctrl_unit, lm_amftable, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_amftable, modulename, 0)
     pge_error_status = pge_errstat_warning
  ELSE
     READ (fit_ctrl_unit, *) have_amftable
     IF ( have_amftable ) THEN
        READ (fit_ctrl_unit, '(A)') static_input_fnames(amf_idx)
        static_input_fnames(amf_idx) = TRIM(ADJUSTL(refdbdir)) // static_input_fnames(amf_idx)
     ENDIF
  END IF

  !xliu: add the following block 
  ! ----------------------------------------------------------
  ! Position cursor to read whether to retrieve ozone profile
  ! ----------------------------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, ozprof_str, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     ozprof_flag = .FALSE.
  ELSE 
     READ (fit_ctrl_unit, *) ozprof_flag
     IF (ozprof_flag)  THEN
        READ (fit_ctrl_unit, '(A)') ozprof_input_fname
     END IF
  END IF
  
  ! -----------------------------------------------
  ! Position cursor to read molecule name(s) to fit
  ! -----------------------------------------------
  IF (.NOT. ozprof_flag) THEN  !xliu
     REWIND ( fit_ctrl_unit )
     CALL skip_to_filemark ( fit_ctrl_unit, molline_str, tmpchar, file_read_stat )
     IF ( file_read_stat /= file_read_ok ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, molline_str, modulename, 0)
        pge_error_status = pge_errstat_error; RETURN
     END IF
     READ (fit_ctrl_unit, '(A)') tmpchar
     CALL get_mols_for_fitting ( tmpchar, n_mol_fit, fitcol_idx, errstat )
     IF ( errstat /= pge_errstat_ok ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_get_molfitname, '', modulename, 0)
        pge_error_status = pge_errstat_error; RETURN
     END IF
  END IF  !xliu
     
  ! -----------------------------------------------------
  ! Position cursor to read GOME channel used for fitting
  ! -----------------------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, lm_bandselect, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_bandselect, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF
  READ (fit_ctrl_unit, *) numwin, do_bandavg, winwav_min, winwav_max
  IF (numwin > maxwin .OR. numwin < 1) THEN
	   WRITE(*, *) 'Number of windows exceeds maxwin or less than 1!!!'
		pge_error_status = pge_errstat_error; RETURN
  ENDIF	
  DO i = 1, numwin
     READ(fit_ctrl_unit, *) band_selectors(i), winpix(i, 1), winpix(i, 2), &
          n_band_avg(i), n_band_samp(i)
  END DO

  IF (do_bandavg) THEN
     IF (ANY(n_band_avg(1:numwin) < 1)) THEN
        WRITE(*, *) 'Number of points for averaging must >= 1!!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF

     IF (ANY(n_band_samp(1:numwin) < 1)) THEN
        WRITE(*, *) 'Number of points for sampling must >= 1!!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF

  ENDIF
  
  ! ------------------------------------------------
  ! Position cursor to read general input parameters
  ! ------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, genline_str, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, genline_str, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF

  READ (fit_ctrl_unit, *) which_slit
  READ (fit_ctrl_unit, *) yn_varyslit, wavcal_sol, smooth_slit, slit_rad
  READ (fit_ctrl_unit, *) slit_fit_pts, n_slit_step, slit_redo
  READ (fit_ctrl_unit, *) wavcal_fit_pts, n_wavcal_step, wavcal_redo
  READ (fit_ctrl_unit, *) yn_smooth
  READ (fit_ctrl_unit, *) yn_doas
  READ (fit_ctrl_unit, *) use_meas_sig 
  READ (fit_ctrl_unit, *) pixnum_lim
  READ (fit_ctrl_unit, *) tol
  READ (fit_ctrl_unit, *) epsrel
  READ (fit_ctrl_unit, *) epsabs
  READ (fit_ctrl_unit, *) epsx

  ! generate identifer for irradiance and radiance spectrum
  i = INDEX(l1_inputs_fname_sol, 'lv1') + 4; sol_identifier = l1_inputs_fname_sol(i:i+7)
  i = INDEX(l1_inputs_fname_rad, 'lv1') + 4; rad_identifier = l1_inputs_fname_rad(i:i+7)
  j = INDEX(l2_output_fname, 'lv2')        ; outdir = l2_output_fname(1:j-1)

  slit_fname    = TRIM(ADJUSTL(outdir)) // 'slit_'    // sol_identifier // '.dat'
  swavcal_fname = TRIM(ADJUSTL(outdir)) // 'swavcal_' // sol_identifier // '.dat'
  rslit_fname   = TRIM(ADJUSTL(outdir)) // 'rslit_'   // rad_identifier // '.dat'
  wavcal_fname  = TRIM(ADJUSTL(outdir)) // 'wavcal_'  // rad_identifier // '.dat'

  ! To determine whether b1a/b1b boundary is changed from 307 nm to 282 nm
  ! according to the date (Does not work after July 1 2005)
  IF (LGT(rad_identifier(1:5), '80606') .OR. LLT(rad_identifier(1:5), '50630')) THEN
     b1ab_change = .TRUE.
  ELSE
     b1ab_change = .FALSE.
  ENDIF
  
  !WRITE(*, '(A100)') l2_output_fname
  !WRITE(*, '(A100)') ch1_out_file
  !WRITE (*, '(A100)') slit_fname  
  !WRITE (*, '(A100)') swavcal_fname   
  !WRITE (*, '(A100)') rslit_fname  
  !WRITE (*, '(A100)') wavcal_fname
 
  ! ------------------------------------------------
  ! Check for consistency of pixel limits to process
  ! ------------------------------------------------
  IF ( ALL ( pixnum_lim < 0 ) )        pixnum_lim(1:2) = (/ 0, forever /)
  IF ( pixnum_lim(1) > pixnum_lim(2) ) pixnum_lim(1) = pixnum_lim(2)

  ! ----------------------------------------------------------
  ! Position cursor to read solar calibration input parameters
  ! ----------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, socline_str, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, socline_str, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF

  ! ---------------------------------------------------
  ! First thing to read is the Solar Reference Spectrum
  ! ---------------------------------------------------
  READ (fit_ctrl_unit, '(A)') refspec_fname(solar_idx)  
  refspec_fname(solar_idx) = TRIM(ADJUSTL(refdbdir)) // refspec_fname(solar_idx)  
  READ (fit_ctrl_unit, *) weight_sun
  READ (fit_ctrl_unit, *) max_itnum_sol

  n_fitvar_sol = 0;  fitvar_sol_init = 0.0
  solpars: DO i = 1, max_calfit_idx

     READ (fit_ctrl_unit, *) idxchar, vartmp, lotmp, uptmp
     ! ---------------------------------------------------------
     ! Check for consitency of bounds and adjust where necessary
     ! ---------------------------------------------------------
     IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
        lotmp = vartmp ; uptmp = vartmp
     END IF
     IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
        uptmp = vartmp ; lotmp = vartmp
     END IF

     IF ( idxchar == eoi3str ) EXIT solpars

     CALL string2index ( calfit_strings, max_calfit_idx, idxchar, sidx )
     IF ( sidx > 0 ) THEN
        fitvar_sol_init(sidx) = vartmp
        lo_sunbnd(sidx) = lotmp ; up_sunbnd(sidx) = uptmp
        IF ( lotmp < uptmp ) THEN
           n_fitvar_sol = n_fitvar_sol + 1
           mask_fitvar_sol(n_fitvar_sol) = i
        ENDIF
     END IF
  END DO solpars
  fitvar_sol_saved = fitvar_sol_init
  lo_sunbnd_init = lo_sunbnd; up_sunbnd_init = up_sunbnd

  ! -------------------------------------------------------------
  ! Position cursor to read radiance calibration input parameters
  ! -------------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, racline_str, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, racline_str, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF
  READ (fit_ctrl_unit, *) renorm
  READ (fit_ctrl_unit, *) weight_rad
  READ (fit_ctrl_unit, *) max_itnum_rad
  READ (fit_ctrl_unit, *) radwavcal_freq
  READ (fit_ctrl_unit, *) szamax
  READ (fit_ctrl_unit, *) zatmos
  READ (fit_ctrl_unit, *) phase

  fitvar_rad_init = 0.0
  fitvar_rad_str = '      '
  radpars: DO i = 1, max_calfit_idx

     READ (fit_ctrl_unit, *) idxchar, vartmp, lotmp, uptmp
     ! ---------------------------------------------------------
     ! Check for consitency of bounds and adjust where necessary
     ! ---------------------------------------------------------
     IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
        lotmp = vartmp ; uptmp = vartmp
     END IF
     IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
        uptmp = vartmp ; lotmp = vartmp
     END IF

     IF ( idxchar == eoi3str ) EXIT radpars
     CALL string2index ( calfit_strings, max_calfit_idx, idxchar, sidx )
     IF ( sidx > 0 ) THEN
        fitvar_rad_init(sidx) = vartmp
        fitvar_rad_str (sidx) = TRIM(ADJUSTL(idxchar))
        lo_radbnd(sidx) = lotmp ; up_radbnd(sidx) = uptmp
     END IF
  END DO radpars

  ! ---------------------------------------------------------
  ! Position cursor to read radiance fitting input parameters
  ! ---------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, rafline_str, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, rafline_str, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF

  ! By default we set the undersampling spectrum to FALSE. Only if we select
  ! it to be included in the fitting does it become trues. This way we save
  ! computation time for cases where we don't include the undersampling.
  have_undersampling = .FALSE.

  ! -------------------------------------------------------------
  ! Now keep reading spectrum blocks until EOF. This, obviously, 
  ! has to be the last READ action performed from the input file.
  ! -------------------------------------------------------------
  getpars: DO j = 1, max_rs_idx

     ! Read the spectrum identification string (SIS)
     READ (UNIT=fit_ctrl_unit, FMT='(A)', IOSTAT=errstat) tmpchar
     IF ( errstat /= file_read_ok ) THEN
        errstat = OMI_SMF_setmsg ( &
             omsao_e_read_fitctrl_file, 'radiance fitting parameters', modulename, 0)
        pge_error_status = pge_errstat_error; RETURN
     END IF
     CALL check_for_endofinput ( TRIM(ADJUSTL(tmpchar)), yn_eoi )
     IF ( yn_eoi ) EXIT getpars

     ! Convert SIS to index
     CALL string2index ( refspec_strings, max_rs_idx, tmpchar, ridx )

     ! GOME specific: the first line in the "spectrum parameter block" is
     ! the name of the corresponding reference spectrum
     READ (UNIT=fit_ctrl_unit, FMT='(A)', IOSTAT=errstat) refspec_fname(ridx)
     refspec_fname(ridx) = TRIM(ADJUSTL(refdbdir)) // refspec_fname(ridx)  

     ! Read the block of fitting parameters for current reference spectrum
     DO k = 1, mxs_idx

        READ (fit_ctrl_unit, *) idxchar, vartmp, lotmp, uptmp
        ! ---------------------------------------------------------
        ! Check for consitency of bounds and adjust where necessary
        ! ---------------------------------------------------------
        IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
           lotmp = vartmp ; uptmp = vartmp
        END IF
        IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
           uptmp = vartmp ; lotmp = vartmp
        END IF

        CALL string2index ( radfit_strings, mxs_idx, idxchar, sidx )
        IF ( sidx > 0 ) THEN
           i = max_calfit_idx + (ridx-1)*mxs_idx + sidx
           fitvar_rad_init(i) = vartmp
           fitvar_rad_str (i) = TRIM(ADJUSTL(tmpchar))
           lo_radbnd (i) = lotmp ; up_radbnd (i) = uptmp
           IF ( (ridx == us1_idx .OR. ridx == us2_idx) .AND. &
                ANY ( (/ vartmp,lotmp,uptmp /) /= 0.0 ) ) have_undersampling = .TRUE.
        END IF
     END DO

     ! read shift parameter
     READ (fit_ctrl_unit, *) idxchar1, vartmp, lotmp, uptmp
     IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
        lotmp = vartmp ; uptmp = vartmp
     END IF
     IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
        uptmp = vartmp ; lotmp = vartmp
     END IF

     i =  max_calfit_idx + (ridx-1)*mxs_idx + 1
     IF  (ALL(lo_radbnd(i:i+2) - up_radbnd(i:i+2) >= 0.0)) THEN
        vartmp = 0.0; uptmp = 0.0; lotmp = 0.0
     END IF

     i =  max_calfit_idx + max_rs_idx * mxs_idx + ridx
     fitvar_rad_init(i) = vartmp
     fitvar_rad_str(i) = TRIM(ADJUSTL(idxchar1))

     lo_radbnd (i) = lotmp ; up_radbnd (i) = uptmp
        
  END DO getpars

  ! -----------------------------------------------------
  ! For safety, copy REFSPEC_FNAME to STATIC_INPUT_FNAMES
  ! (in the OMI branch this is the other way around since
  !  there we get file names from the PCF)
  ! -----------------------------------------------------
  DO j = solar_idx, max_rs_idx
     static_input_fnames(j) = refspec_fname(j)
  END DO

  ! -------------------------------------------------------------
  ! Find the indices of those variables that are actually varied
  ! during the fitting, and save those in MASK_FITVAR_RAD. Save
  ! the number of varied parameters in N_FITVAR_RAD.
  !
  ! In addition, we need to determine the fitting indices that
  ! will make up the final fitted column of the molecule(s) in
  ! question. For this we have to jump through a double loop:
  ! Since we are compressing the fitting parameter array to 
  ! include only the varied parameters, the final covariance 
  ! matrix, which is crucial for determining the uncertainties,
  ! only knows the compressed indices. Therefore we have to have
  ! an "index of an index" type array, that remembers the index
  ! position of indices of the final molecule(s), AS THEY APPEAR
  ! IN THE COMPRESSED FITTING PARAMETER LIST.
  !
  ! For the latter task it is easier to split the loops into
  ! calibration parameters (unrelated to the final fitted column)
  ! and reference spectra parameters.
  ! -------------------------------------------------------------
  n_fitvar_rad = 0 ; mask_fitvar_rad = 0
  ! --------------------------------
  ! First the calibration parameters
  ! --------------------------------
  DO i = 1, max_calfit_idx
     IF ( lo_radbnd(i) < up_radbnd(i) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = i
     END IF
  END DO

  ! ------------------------------------
  ! Now the reference spectra parameters
  ! ------------------------------------
  n_fincol_idx = 0
  DO i = 1, max_rs_idx
     idx = max_calfit_idx + (i-1) * mxs_idx
     DO j = mns_idx, mxs_idx        
        ! ----------------------------------------------
        ! Assign only entries that are varied in the fit
        ! ----------------------------------------------
        IF ( lo_radbnd(idx+j) < up_radbnd(idx+j) ) THEN
           n_fitvar_rad = n_fitvar_rad + 1
           mask_fitvar_rad(n_fitvar_rad) = idx+j

           ! -------------------------------------------------
           ! And here the loop over the final column molecules.
           ! We have to match the FITCOL_IDX with the current
           ! molecule index, and then remember the position of
           ! the fitting index in the MASK_FITVAR_RAD array.
           ! The second index remembers the reference spectrum
           ! that is associated with this molecule, so that we
           ! can easily access its normalization factor.
           ! -------------------------------------------------
           IF (.NOT. ozprof_flag) THEN  !xliu
              getfincol: DO k = 1, n_mol_fit
                 IF ( fitcol_idx(k) == i ) THEN
                    n_fincol_idx = n_fincol_idx + 1
                    fincol_idx (1,n_fincol_idx) = n_fitvar_rad
                    fincol_idx (2,n_fincol_idx) = i
                    EXIT getfincol
                 END IF
              END DO getfincol
           ENDIF   !xliu

        END IF
     END DO
  END DO

  DO i = 1, max_rs_idx
     j = max_calfit_idx + max_rs_idx * mxs_idx + i
     IF ( lo_radbnd(j) < up_radbnd(j) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = j
     END IF
  END DO

  ! -----------------------------------------------
  ! Close fitting control file, report SUCCESS read
  ! -----------------------------------------------
  CLOSE ( UNIT=fit_ctrl_unit )

  !xliu: add the following block
  ! ------------------------------------------------------------------------------
  ! Read fitting conrol parameters from input file for ozone profile variables
  IF (ozprof_flag) THEN 
     CALL read_ozprof_input ( &
          fit_ctrl_unit, ozprof_input_fname, pge_error_status )
     IF ( pge_error_status >= pge_errstat_error ) RETURN 
  END IF
  ! ------------------------------------------------------------------------------

  fitvar_rad_saved = fitvar_rad_init

  IF ( yn_doas ) THEN
     pm_one     = -1.D0
  ELSE
     pm_one     = 1.D0
  END IF

  RETURN
END SUBROUTINE gome_read_fitting_control_file
