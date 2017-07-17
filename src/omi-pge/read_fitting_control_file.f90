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

SUBROUTINE read_fitting_control_file (fit_ctrl_unit, fit_ctrl_file,  &
     instrument_idx, l1_inputs_fname_sol, l1_inputs_fname_rad, l2_output_fname, &
     l2_cld_fname, l2_hdf_flag, pge_error_status)     

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
       molline_str, eoi3str, solar_idx, amf_idx, us1_idx, us2_idx, &
       shift_offset, comm_idx, com1_idx, comfidx, cm1fidx, comvidx, cm1vidx
  USE OMSAO_parameters_module,   ONLY: maxchlen, maxwin, max_fit_pts,max_mol_fit 
  USE OMSAO_variables_module,    ONLY: use_backup, use_solcomp, avg_solcomp, avgsol_allorb,use_meas_sig, &
       fitcol_idx, fincol_idx, n_fincol_idx, n_mol_fit, &
       weight_sun, max_itnum_sol,weight_rad,max_itnum_rad,renorm, szamax, zatmos,  &
       which_slit, slit_trunc_limit, yn_varyslit, wavcal, wavcal_sol,smooth_slit,slit_rad, &
       slit_fit_pts,n_slit_step,slit_redo,wavcal_fit_pts, n_wavcal_step, wavcal_redo,   &
       yn_smooth, yn_doas,pm_one,  tol,  epsrel,  epsabs,  epsx,radwavcal_freq, phase, &
       n_fitvar_sol, fitvar_sol_init, fitvar_sol_saved, mask_fitvar_rad,fitvar_rad_str,      &
       n_fitvar_rad, fitvar_rad_init, fitvar_rad_saved, mask_fitvar_sol,fitvar_rad_unit,rmask_fitvar_rad,&  
       lo_sunbnd, up_sunbnd, lo_sunbnd_init,up_sunbnd_init, lo_radbnd, up_radbnd, &   
       static_input_fnames,refspec_fname, wavcal_fname, swavcal_fname,slit_fname,rslit_fname,&  
       l2_swathname,l1b_rad_filename,outdir, atmdbdir, refdbdir, &
       reduce_resolution,reduce_slit, redsampr, redlam, redfixwav,use_redfixwav,redfixwav_fname,nredfixwav, &
       rm_mgline,numwin, do_bandavg,wcal_bef_coadd,band_selectors,  n_band_avg, n_band_samp,& 
       winlim,winwav_min, winwav_max,linenum_lim,  pixnum_lim,coadd_uv2, &   
       have_amftable,have_undersampling, sol_identifier, rad_identifier,   &
       scnwrt, database_indices, radnhtrunc, refnhextra, &
       omi_idx, scia_idx, gome2_idx, gome_idx, which_instrument, max_instrument_idx, write_is_slit
  USE OMSAO_gome_data_module, ONLY:lm_gome_solspec,lm_gome_eshine, &
                                   n_gome_max_pts, n_gome_data_dim, gome_spec_missing, gome_orbc       
  USE OMSAO_errstat_module
  USE OMSAO_omidata_module, ONLY: orbc, orbnum, orbcsol, orbnumsol, mswath, nswath,    &
       upper_wvls, lower_wvls, nxtrack_max, ntimes_max, ncoadd, do_xbin, do_ybin,      &
       nxbin, nybin, ncoadd, nlines_max, omiraddate, retlbnd, retubnd, omisol_version, &
       omi_redslw
  USE ozprof_data_module,   ONLY: ozprof_str, ozprof_flag, ozprof_input_fname,         &
       fullorb, do_ch2reso, l1l2inp_unit, nos, nsh, nsl, do_simu, radcalwrt

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
  INTEGER,           INTENT (OUT)   :: pge_error_status, l2_hdf_flag
  CHARACTER (LEN=*), INTENT (OUT)   :: l1_inputs_fname_sol, l1_inputs_fname_rad, &
       l2_output_fname, l2_cld_fname

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER                  :: i, j, k, file_read_stat, sidx, ridx, idx, nidx, cldorb, ntsh
  INTEGER, DIMENSION(2)    :: pixlim, linelim
  CHARACTER (LEN=maxchlen) :: tmpchar, l1l2_files, slitdir, slitc
  CHARACTER (LEN=3)        :: idxchar, xbinchar, ybinchar
  CHARACTER (LEN=5)        :: idxchar1, cldorbc
  CHARACTER (LEN=4)        :: slinechar, elinechar
  CHARACTER (LEN=2)        :: sxchar, exchar
  LOGICAL                  :: yn_eoi, rw_l1l2_here, select_lonlat
  REAL      (KIND=dp)      :: vartmp, lotmp, uptmp, slat, elat, slon, elon

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=30), PARAMETER :: modulename = 'read_fitting_control_file'

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
  CHARACTER (LEN=16), PARAMETER :: lm_l2hdf      = 'Write HDF output'
  CHARACTER (LEN=21), PARAMETER :: lm_amftable   = 'Air mass factor table'
  CHARACTER (LEN=30), PARAMETER :: lm_atmdb      = 'Atmospheric database directory'
  CHARACTER (LEN=27), PARAMETER :: lm_refdb      = 'Reference spectra directory'
  CHARACTER (LEN=30), PARAMETER :: lm_bandselect = 'OMI radiance bands to be used'
  CHARACTER (LEN=27), PARAMETER :: lm_reduceres  = 'Reduce spectral resolution'
  CHARACTER (LEN=22), PARAMETER :: lm_wrtfit     = 'Write Fitting Results?'
  

  ! =================================
  ! External OMI and Toolkit routines
  ! =================================
  INTEGER :: OMI_SMF_setmsg

  pge_error_status = pge_errstat_ok
  write_is_slit = .false.
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

  ! -------------------------------------------
  ! Position cursor to read Level 1 input files
  ! -------------------------------------------
  select_lonlat = .FALSE.
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
        READ (fit_ctrl_unit, '(A)') l2_cld_fname
        pixlim = -9999; linelim = -9999
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
        READ(l1l2inp_unit, '(A)') l2_cld_fname
        READ(l1l2inp_unit, '(A)') l2_output_fname
        READ(l1l2inp_unit, *)     linelim, pixlim 
        READ(l1l2inp_unit, *, IOSTAT=errstat)  select_lonlat
        IF ( errstat /= pge_errstat_ok ) select_lonlat = .FALSE.
        IF (select_lonlat) THEN
           READ(l1l2inp_unit, *)   slat, elat, slon, elon
           IF (slat >= elat .OR. slon >= elon) THEN
              WRITE(*, *) 'Incorrect lat/lon range!!!'
              pge_error_status = pge_errstat_error; RETURN
           ENDIF
           pixlim(1)  = -5;  pixlim(2) = -5
           linelim(1) = -5; linelim(2) = -5
        ENDIF
        CLOSE(UNIT=l1l2inp_unit) 
     END IF
  END IF
  l1b_rad_filename = l1_inputs_fname_rad

  IF (instrument_idx == omi_idx) THEN
     ! obtain orbit number from irradiance file
     i = INDEX(l1_inputs_fname_sol, '-o') + 2
     orbcsol = l1_inputs_fname_sol(i : i + 5)
     READ (orbcsol, *) orbnumsol
     READ (l1_inputs_fname_sol(i+7 : i + 9), *) omisol_version 

     ! obtain orbit number from radiance file
     i = INDEX(l1_inputs_fname_rad, '-o') + 2
     orbc = l1_inputs_fname_rad(i : i + 5)
     READ (orbc, *) orbnum

     ! check cloud file
     i = INDEX(l2_cld_fname, '-o') + 2
     cldorbc = l2_cld_fname(i : i + 5)
     READ (cldorbc, *) cldorb

     IF (cldorb /= orbnum) THEN
        WRITE(*, *) 'Inconsistent orbit number between radiance and cloud file!!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
    
     ! generate identifer for irradiance and radiance spectrum
     i = INDEX(l2_output_fname, 'OMIO3PROF')   ; outdir = l2_output_fname(1:i-1)
     rad_identifier = 'o' // orbc

     slit_fname    = TRIM(ADJUSTL(outdir)) // 'slit_o'    // orbcsol
     swavcal_fname = TRIM(ADJUSTL(outdir)) // 'swavcal_o' // orbcsol
     rslit_fname   = TRIM(ADJUSTL(outdir)) // 'rslit_o'   // orbc
     wavcal_fname  = TRIM(ADJUSTL(outdir)) // 'wavcal_o'  // orbc
  ELSE
     IF (instrument_idx == gome_idx) THEN
        i = INDEX(l1_inputs_fname_sol, 'lv1_') + 11
        orbc = l1_inputs_fname_rad(i-2:i)
        orbcsol = orbc
        i = INDEX(l1_inputs_fname_sol, 'lv1') + 4; sol_identifier = l1_inputs_fname_sol(i:i+7)
        i = INDEX(l1_inputs_fname_rad, 'lv1') + 4; rad_identifier = l1_inputs_fname_rad(i:i+7)
        j = INDEX(l2_output_fname, 'lv2')        ; outdir = l2_output_fname(1:j-1)
        !ELSEIF (instrument_idx == scia_idx) THEN
        !   i = INDEX(l1_inputs_fname_rad, 'Ch1orb') + 28
        !   orbc = l1_inputs_fname_rad(i-2:i) 
        !   i = INDEX(l1_inputs_fname_sol, 'Ch1orb') + 20
        !   sol_identifier = l1_inputs_fname_sol(i:i+4) // l1_inputs_fname_sol(i+6:i+8) 
        !   rad_identifier = sol_identifier
        !   j = INDEX(l2_output_fname, 'lv2')        ; outdir = l2_output_fname(1:j-1)
        !   sciaorb_identifier =  l1_inputs_fname_sol(i-3:i+11)
     ELSEIF (instrument_idx == gome2_idx) THEN
        orbc = '0'
        i = INDEX(l1_inputs_fname_sol, 'GOME_xxx_1B') + 18; sol_identifier = l1_inputs_fname_sol(i:i+7)
        i = INDEX(l1_inputs_fname_rad, 'GOME_xxx_1B') + 18; rad_identifier = l1_inputs_fname_rad(i:i+7)
        j = INDEX(l2_output_fname, 'lv2')            ; outdir = l2_output_fname(1:j-1)
     ENDIF
     slit_fname    = TRIM(ADJUSTL(outdir)) // 'slit_'    // sol_identifier // '.dat'
     swavcal_fname = TRIM(ADJUSTL(outdir)) // 'swavcal_' // sol_identifier // '.dat'
     rslit_fname   = TRIM(ADJUSTL(outdir)) // 'rslit_'   // rad_identifier // '.dat'
     wavcal_fname  = TRIM(ADJUSTL(outdir)) // 'wavcal_'  // rad_identifier // '.dat'
  ENDIF
  ! ----------------------------------------------
  ! Position cursor to read HDF output flags (CRN)
  ! ----------------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, lm_l2hdf, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_w_read_fitctrl_file, lm_l2hdf, modulename, 0)
     pge_error_status = pge_errstat_warning
     l2_hdf_flag = 0  ! Write ASCII
  ELSE
     READ (fit_ctrl_unit, '(I)') l2_hdf_flag
     
     ! Check output data format
     IF (instrument_idx == gome2_idx .AND. (l2_hdf_flag < 1 .OR. l2_hdf_flag > 2)) THEN
        WRITE(*, *) 'GOME-2 data can only be written in HDF'
        pge_error_status = pge_errstat_error; RETURN
     ELSE IF ((instrument_idx == gome_idx .OR. instrument_idx == scia_idx) .AND. l2_hdf_flag  > 0) THEN
        WRITE(*, *) 'GOME-1/SCIA data can only be written in ASCII'
        pge_error_status = pge_errstat_error; RETURN
     ELSE IF (instrument_idx == omi_idx .AND. (l2_hdf_flag /= 0 .AND. l2_hdf_flag /= 3)) THEN
        WRITE(*, *) 'OMI data can only be written in ASCII or HDF-EOS5'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF     
  END IF




 ! ----------------------------------------------
  ! Position cursor to write Fitting results to file JBAK
  ! ----------------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, lm_wrtfit, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_w_read_fitctrl_file, lm_wrtfit, modulename,0)
     pge_error_status = pge_errstat_warning
  ELSE
     READ (fit_ctrl_unit,     *) write_is_slit
     READ (fit_ctrl_unit, '(A)') slitdir
  END IF

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
     WRITE(*, *) 'This algorithm is only for ozone profile retrieval!!!'
     pge_error_status = pge_errstat_error; RETURN
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
     
  !xliu, 01/03/2007, read options to degrade spectral resolution
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, lm_reduceres, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_bandselect, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF
  READ (fit_ctrl_unit, *) reduce_resolution
  READ (fit_ctrl_unit, *) reduce_slit
  READ (fit_ctrl_unit, *) omi_redslw(1:mswath)
  IF ( reduce_slit == 1 ) omi_redslw(1:mswath) = omi_redslw(1:mswath) / 1.66511  ! convert from FWHM to hw1e 
  READ (fit_ctrl_unit, *) use_redfixwav
  IF (.NOT. reduce_resolution) use_redfixwav = .FALSE.
  READ (fit_ctrl_unit, '(A)') redfixwav_fname
  IF (use_redfixwav) THEN
     OPEN(UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(redfixwav_fname)), STATUS='OLD', IOSTAT=errstat)
     IF ( errstat /= pge_errstat_ok ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
             TRIM(ADJUSTL(redfixwav_fname)), modulename, 0)
        pge_error_status = pge_errstat_error; RETURN
     ELSE
        READ(l1l2inp_unit, *) nredfixwav
        IF (nredfixwav > max_fit_pts) THEN
           pge_error_status = pge_errstat_error; RETURN
        ENDIF
        READ(l1l2inp_unit, *) (redfixwav(i), i = 1, nredfixwav)
        CLOSE(UNIT=l1l2inp_unit) 
     END IF
  ENDIF
  READ (fit_ctrl_unit, *) redsampr 
  READ (fit_ctrl_unit, *) redlam

  ! -----------------------------------------------------
  ! Position cursor to read OMI channels used for fitting
  ! -----------------------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, lm_bandselect, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_bandselect, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF
  READ (fit_ctrl_unit, *) do_xbin, nxbin
  READ (fit_ctrl_unit, *) do_ybin, nybin
  READ (fit_ctrl_unit, *) rm_mgline
  READ (fit_ctrl_unit, *) numwin, do_bandavg, wcal_bef_coadd !, winwav_min, winwav_max
  IF (reduce_resolution) THEN
     do_bandavg = .FALSE.; rm_mgline = .FALSE.
  ENDIF
  IF (numwin > maxwin .OR. numwin < 1) THEN
     WRITE(*, *) 'Number of windows exceeds maxwin or less than 1!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  
  retlbnd = 1000.0; retubnd = 0.0
  DO i = 1, numwin
     READ(fit_ctrl_unit, *) band_selectors(i), winlim(i, 1), winlim(i, 2), &
          n_band_avg(i), n_band_samp(i)
     
     IF ((band_selectors(i) < 0) .OR. (band_selectors(i) > mswath)) THEN
        WRITE(*, *) 'No such bands exist !!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
     
     IF (numwin > 1) THEN
        IF (winlim(i, 1) < lower_wvls(band_selectors(i)) .OR. winlim(i, 2) > upper_wvls(band_selectors(i))) THEN
           WRITE(*, *) 'Specified fitting windows does not make sense!!!'
           pge_error_status = pge_errstat_error; RETURN
        ENDIF
     ELSE
        ! Allow 2 extra nm for radiance calibration
        IF (winlim(i, 1) < lower_wvls(band_selectors(i)) - 1.0 .OR. winlim(i, 2) > upper_wvls(band_selectors(i)) + 1.0) THEN
           WRITE(*, *) 'Specified fitting windows does not make sense!!!'
           pge_error_status = pge_errstat_error; RETURN
        ENDIF
     ENDIF

     IF (i > 1) THEN
        IF  (band_selectors(i) < band_selectors(i-1) .OR. winlim(i, 1) < winlim(i-1, 2))  THEN
           WRITE(*, *) 'Incorrect band selection (must be in increasing wavelength) !!!'
           pge_error_status = pge_errstat_error; RETURN
        ENDIF
     ENDIF

     IF (winlim(i, 1) < retlbnd(band_selectors(i))) retlbnd(band_selectors(i)) = winlim(i, 1)
     IF (winlim(i, 2) > retubnd(band_selectors(i))) retubnd(band_selectors(i)) = winlim(i, 2)
  END DO 
  IF (MAXVAL(n_band_avg(1:numwin)) == 1) do_bandavg = .FALSE.

  DO i = 1, mswath
     IF (retlbnd(i) == 1000.0) retlbnd(i) = lower_wvls(i)
     IF (retubnd(i) == 0.0)    retubnd(i) = upper_wvls(i)
  ENDDO
  
  IF (ANY(band_selectors(1:numwin) == 1) .AND. ANY(band_selectors(1:numwin) == 2)) THEN
     coadd_uv2 = .TRUE.; nswath = 2
  ELSE
     coadd_uv2 = .FALSE.; nswath = 1
  ENDIF
  IF (nswath == 1 .AND. band_selectors(1) == 1) THEN
     nswath = 2 ! have to read measurements around 370 nm
     coadd_uv2 = .TRUE.
  ENDIF

  IF (coadd_uv2) THEN 
     ncoadd = 2 
  ELSE
     ncoadd = 1
  ENDIF
  
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
  winwav_min = winlim(1, 1) - 10.0
  winwav_max = winlim(numwin, 2) + 10.0
    
  ! ------------------------------------------------
  ! Position cursor to read general input parameters
  ! ------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, genline_str, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, genline_str, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF

  READ (fit_ctrl_unit, *) use_backup, use_solcomp, avg_solcomp, avgsol_allorb
  READ (fit_ctrl_unit, *) which_slit
  IF (reduce_resolution .AND. which_slit /= 0 .AND. which_slit /= 3) THEN
     WRITE(*, *) 'Have to use consistent slit function!!!'
     IF (reduce_slit == 1) which_slit = 0
     IF (reduce_slit == 2) which_slit = 3
  ENDIF

  READ (fit_ctrl_unit, *) slit_trunc_limit
  READ (fit_ctrl_unit, *) yn_varyslit, wavcal, wavcal_sol, smooth_slit, slit_rad
  IF (reduce_resolution .AND. yn_varyslit) THEN
     yn_varyslit = .FALSE.
     WRITE(*, *) 'Could not use variable slit function when to reduce resolution!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
 ! IF (.NOT. wavcal) wavcal_sol = .FALSE.
  READ (fit_ctrl_unit, *) slit_fit_pts, n_slit_step, slit_redo
  READ (fit_ctrl_unit, *) wavcal_fit_pts, n_wavcal_step, wavcal_redo
  READ (fit_ctrl_unit, *) yn_smooth
  READ (fit_ctrl_unit, *) yn_doas
  READ (fit_ctrl_unit, *) use_meas_sig 
  READ (fit_ctrl_unit, *) linenum_lim
  READ (fit_ctrl_unit, *) pixnum_lim
  READ (fit_ctrl_unit, *) tol
  READ (fit_ctrl_unit, *) epsrel
  READ (fit_ctrl_unit, *) epsabs
  READ (fit_ctrl_unit, *) epsx


  !====================
  ! slit_names
  !=====================
   select case (which_slit)
    case (0) ; slitc = 'sga'
    case (1) ; slitc = 'aga'
    case (2) ; slitc = 'voi'
    case (3) ; slitc = 'tra'
    case (4) ; slitc = 'spg'
    case (5) ; slitc = 'ins'
  end select
  write(tmpchar, '(f5.1, "_",f5.1)') winlim(1,1), winlim(numwin,2)
  slitc='_'//TRIM(ADJUSTL(slitc))//'_'//TRIM(ADJUSTL(tmpchar))//'.dat'
  slit_fname    = TRIM(ADJUSTL(slitdir))//'slit_o'//orbcsol//TRIM(ADJUSTL(slitc))
  swavcal_fname = TRIM(ADJUSTL(slitdir))//'wavcal_o'//orbcsol//TRIM(ADJUSTL(slitc))
  rslit_fname   = TRIM(ADJUSTL(slitdir))//'rslit_o'//orbc//TRIM(ADJUSTL(slitc))
  wavcal_fname  = TRIM(ADJUSTL(slitdir))//'rwavcal_o'//orbc//TRIM(ADJUSTL(slitc))
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
        !print * , idxchar, sidx
     IF ( sidx > 0 ) THEN
        fitvar_sol_init(sidx) = vartmp
        lo_sunbnd(sidx) = lotmp ; up_sunbnd(sidx) = uptmp
        IF ( lotmp < uptmp ) THEN
           n_fitvar_sol = n_fitvar_sol + 1
           mask_fitvar_sol(n_fitvar_sol) = i
           ! print *, 'fitsol:', n_fitvar_sol,i, idxchar, vartmp
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


  ! Read date from radiance file (used for correcting sun-earth distance when using backupirradiance)
  i = INDEX(l1_inputs_fname_rad, '-o') -14
  omiraddate = l1_inputs_fname_rad(i : i + 8)

  ! ------------------------------------------------
  ! Check for consistency of pixel limits to process
  ! ------------------------------------------------

  IF (select_lonlat) THEN
     CALL find_scan_line_range(slat, elat, slon, elon, linelim(1), &
          linelim(2), pixlim(1), pixlim(2), pge_error_status )
     IF (pixlim(1) < 0 .OR. linelim(1) < 0 .OR. pge_error_status >= pge_errstat_error) THEN
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
     ! pixlim is based on UV1, if both channels are selected
     IF (coadd_uv2) THEN
        pixlim(1) = pixlim(1) * ncoadd - 1; pixlim(2) = pixlim(2) * ncoadd
     ENDIF
  ENDIF

 
  IF (linelim(1) /= -9999) linenum_lim =  linelim
  IF (pixlim(1)  /= -9999) pixnum_lim  =  pixlim

  IF ( ALL ( linenum_lim < 0 ) )         linenum_lim(1:2) = (/ 1, ntimes_max /)
  IF ( linenum_lim(1) > linenum_lim(2) ) linenum_lim([1, 2]) = linenum_lim([2, 1])   
  IF ( linenum_lim(1) < 1 )              linenum_lim(1) = 1
  IF ( linenum_lim(2) > ntimes_max )     linenum_lim(2) = ntimes_max

  IF ( ALL ( pixnum_lim < 0 ) )          pixnum_lim(1:2) = (/ 1, nxtrack_max /)
  IF ( pixnum_lim(1) > pixnum_lim(2) )   pixnum_lim([1, 2]) = pixnum_lim([2, 1])   
  IF ( pixnum_lim(1) < 1 )               pixnum_lim(1) = 1
  IF ( pixnum_lim(2) > nxtrack_max )     pixnum_lim(2) = nxtrack_max
  

  ! check for selected across track position (must start from odd positions)
  IF (coadd_uv2)  THEN
     i = pixnum_lim(2)-pixnum_lim(1) + 1
     IF ( MOD(pixnum_lim(1), ncoadd) /= 1 .OR. MOD(i, ncoadd) /= 0 ) THEN
        WRITE(*, '(A,2I4)') 'Incorrect across track positions to be coadded: ', pixnum_lim(1:2)
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
     
     pixnum_lim = CEILING(1.0 * pixnum_lim / ncoadd)  !nint 
  ENDIF

  WRITE(slinechar, '(I4.4)') linenum_lim(1); WRITE(elinechar, '(I4.4)') linenum_lim(2)
  WRITE(sxchar, '(I2.2)')    pixnum_lim(1) ; WRITE(exchar, '(I2.2)')    pixnum_lim(2)
   
  ! must divide and must start from odd coadded positions
  IF (do_xbin .AND. nxbin > 1) THEN
     i = pixnum_lim(2)-pixnum_lim(1) + 1
     IF ( MOD (i, nxbin) /= 0 .OR. MOD(pixnum_lim(1), nxbin) /= 1 ) THEN
        WRITE(*, '(A,2I4)') 'Incorrect across track binning option: ', pixnum_lim(1:2)
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
  ELSE
     nxbin = 1
  ENDIF

  ! Could start from any positions, adjust line positions if necessary
  IF (do_ybin .AND. nybin > 1)  THEN
     i = linenum_lim(2)-linenum_lim(1) 
     IF (MOD(i, nybin) /= 0) THEN
        linenum_lim(2) = NINT(1.0 * i / nybin) * nybin + linenum_lim(1) - 1
        IF (linenum_lim(2) > ntimes_max) linenum_lim(2) = linenum_lim(2) - nybin
     ENDIF
  ELSE
     nybin = 1
  ENDIF

  IF (nxbin == 1) do_xbin = .FALSE.
  IF (nybin == 1) do_ybin = .FALSE.

  IF ( .NOT. (linenum_lim(1) == 1 .AND. linenum_lim(2) == ntimes_max)) THEN
     l2_output_fname = TRIM(ADJUSTL(l2_output_fname)) // '_L' // slinechar // '-' // elinechar
  ENDIF
  
  IF ( .NOT. (pixnum_lim(1) == 1 .AND. (pixnum_lim(2) == nxtrack_max .OR. &
       (coadd_uv2 .AND. pixnum_lim(2) == nxtrack_max / ncoadd)))) THEN
     l2_output_fname = TRIM(ADJUSTL(l2_output_fname)) // '_X' // sxchar // '-' // exchar
  ENDIF

  IF (do_xbin) THEN
     WRITE(xbinchar, '(A2,I1)') 'BX', nxbin
     l2_output_fname = TRIM(ADJUSTL(l2_output_fname)) // '-' // xbinchar 
  ENDIF
  IF (do_ybin) THEN
     WRITE(ybinchar, '(A2,I1)') 'BY', nybin
     l2_output_fname = TRIM(ADJUSTL(l2_output_fname)) // '-' // ybinchar 
  ENDIF
     
  IF ( instrument_idx /= gome2_idx) THEN
     IF (instrument_idx == omi_idx .AND. l2_hdf_flag == 3) THEN
        l2_output_fname = TRIM(ADJUSTL(l2_output_fname)) // '.he5'
        l2_swathname = 'OzoneProfile'
     ELSE      
        l2_output_fname = TRIM(ADJUSTL(l2_output_fname)) // '.out'
     ENDIF
  ELSE 
     l2_output_fname = TRIM(ADJUSTL(l2_output_fname))
  ENDIF
           
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
  comvidx = 0; cm1vidx = 0; comfidx = 0; cm1fidx = 0
  fitvar_rad_unit = 'NoUnits'
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

           IF ( ridx == comm_idx .AND. lotmp < uptmp ) THEN
              comvidx = i
           ENDIF

           IF ( ridx == com1_idx .AND. lotmp < uptmp ) THEN
              cm1vidx = i
           ENDIF
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

     i =  shift_offset + ridx
     fitvar_rad_init(i) = vartmp
     fitvar_rad_str(i) = TRIM(ADJUSTL(idxchar1))
     fitvar_rad_unit(i) = 'nm'
     
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
  rmask_fitvar_rad = 0; database_indices = 0
  ! --------------------------------
  ! First the calibration parameters
  ! --------------------------------
  DO i = 1, max_calfit_idx
     IF ( lo_radbnd(i) < up_radbnd(i) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = i
        rmask_fitvar_rad(i) = n_fitvar_rad      
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
           rmask_fitvar_rad(idx+j) = n_fitvar_rad
           database_indices(n_fitvar_rad) = i

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

           IF (idx + j == comvidx) comfidx = n_fitvar_rad
           IF (idx + j == cm1vidx) cm1fidx = n_fitvar_rad  
        !print *,'fitref' , n_fitvar_rad,idx+j,fitvar_rad_str(idx+j),fitvar_rad_init(idx+j), fitvar_rad_unit(idx+j)

        END IF
     END DO
  END DO

  ! For shift
  ntsh = 0
  DO i = 1, max_rs_idx
     j = shift_offset + i
     IF ( lo_radbnd(j) < up_radbnd(j) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = j
        rmask_fitvar_rad(j) = n_fitvar_rad
        ntsh = ntsh + 1
        !print *,'shift' , n_fitvar_rad,idx+j,fitvar_rad_str(idx+j),fitvar_rad_init(idx+j), fitvar_rad_unit(idx+j)
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

  ! refnhextra must >= 1 and radnhtrunc > refnhextra, if interpolation is performed
  ! radnhtrunc should be 
  IF ( (ntsh == 0 .AND. nsh == 0 .AND. nos == 0 .AND. nsl == 0) &
       .OR. (do_simu .AND. .NOT. radcalwrt)) THEN
     radnhtrunc = 3; refnhextra = 2
  !ELSE IF (reduce_resolution .AND. use_redfixwav) THEN
  !   radnhtrunc = 2; refnhextra = 1
  ELSE
     radnhtrunc = 2; refnhextra = 1
  ENDIF
     

  ! ------------------------------------------------------------------------------
  fitvar_rad_saved = fitvar_rad_init

  IF ( yn_doas ) THEN
     pm_one     = -1.D0
  ELSE
     pm_one     = 1.D0
  END IF

  RETURN
END SUBROUTINE read_fitting_control_file


SUBROUTINE get_mols_for_fitting ( tmpchar, n_mol_fit, fitcol_idx, errstat )

  USE OMSAO_indices_module,    ONLY: refspec_strings, max_rs_idx
  USE OMSAO_parameters_module, ONLY: max_mol_fit
  USE OMSAO_errstat_module,    ONLY: pge_errstat_error

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  CHARACTER (LEN=*), INTENT (IN) :: tmpchar

  ! ================
  ! Output variables
  ! ================
  INTEGER,                          INTENT (OUT) :: n_mol_fit, errstat
  INTEGER, DIMENSION (max_mol_fit), INTENT (OUT) :: fitcol_idx

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                      :: i, ncl, sstart, sidx
  LOGICAL                      :: yn_eoc
  CHARACTER (LEN=LEN(tmpchar)) :: tmpsub

  ! ----------------------------
  ! Initialize output quantities
  ! ----------------------------
  n_mol_fit = 0  ;  fitcol_idx = 0

  ! -----------------------------------------------
  ! Get names and indices for main molecules to fit
  ! -----------------------------------------------
  sstart = 0 ; ncl = 0 ; yn_eoc = .FALSE.
  getmolnames: DO i = 1, max_mol_fit
     ! ---------------------------------------------------------
     ! Extract index string, find index, then extract file name.
     ! ---------------------------------------------------------
     CALL get_substring ( tmpchar, sstart, tmpsub, ncl, yn_eoc )
     IF ( ncl > 0 ) THEN
        CALL string2index ( refspec_strings, max_rs_idx, tmpsub, sidx )
        IF ( sidx > 0 ) THEN
           n_mol_fit = n_mol_fit + 1
           fitcol_idx(n_mol_fit) = sidx
        END IF
     END IF
     IF ( yn_eoc ) EXIT getmolnames
  END DO getmolnames
  IF ( n_mol_fit == 0 .OR. ALL(fitcol_idx == 0) ) errstat = pge_errstat_error

  RETURN
END SUBROUTINE get_mols_for_fitting


