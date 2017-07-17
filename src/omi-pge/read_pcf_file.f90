SUBROUTINE read_pcf_file ( pge_error_status )

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: &
       max_rs_idx, amf_idx, refspec_strings, verbosity_lun, pge_mol_id, &
       pge_l1b_radiance_lun, &
       pge_l1b_irradiance_lun, icf_idx, pge_static_input_luns, pge_l2_output_luns, &
       orbitnumber_lun, granule_s_lun, granule_e_lun
  USE OMSAO_parameters_module,     ONLY: zerospec_string, vb_lev_omidebug, vb_lev_develop
  USE OMSAO_variables_module,      ONLY: &
       pge_idx, pge_name, verb_thresh_char, verb_thresh_lev, pcf_orbit_number, &
       l1b_rad_filename, l1b_irrad_filename, l2_filename, static_input_fnames, have_amftable
  USE OMSAO_metadata_module,       ONLY: pcf_granule_s_time,  pcf_granule_e_time
  USE OMSAO_errstat_module


  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=13), PARAMETER :: modulename = 'read_pcf_file'


  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER, INTENT (OUT) :: pge_error_status

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER :: i, strlen
  CHARACTER (LEN=PGSd_PC_VALUE_LENGTH_MAX) :: tmpchar

  ! ------------------------
  ! Error handling variables
  ! ------------------------
  INTEGER :: errstat, version  

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: &
       pgs_pc_getnumberoffiles, pgs_pc_getreference, pgs_pc_getconfigdata, &
       pgs_smf_teststatuslevel, OMI_SMF_setmsg

  ! --------------------------------------------------------------------
  ! Get the verbosity threshold. The setting of this variable determines
  ! how "chatty" in terms of screen output the PGE will be.
  ! --------------------------------------------------------------------
  errstat = PGS_PC_GetConfigData ( verbosity_lun, verb_thresh_char )
  IF ( PGS_SMF_TestStatusLevel(errstat) /= pgs_smf_mask_lev_s ) THEN
     errstat = OMI_SMF_setmsg (omsao_w_getlun, 'VERBOSITY_LUN', modulename, 0)
     verb_thresh_char = '1'
     pge_error_status = pge_errstat_warning
  ENDIF
  READ (verb_thresh_char, '(I1)') verb_thresh_lev  ! Convert character to integer

  ! ----------------------------------------------------------------------
  ! Get the orbit number. This will be checked against L1B metadata later.
  ! ----------------------------------------------------------------------
  errstat = PGS_PC_GetConfigData ( orbitnumber_lun, tmpchar )
  IF ( PGS_SMF_TestStatusLevel(errstat) /= pgs_smf_mask_lev_s ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_getlun, 'ORBITNUMBER_LUN', modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  READ (tmpchar, '(I5)') pcf_orbit_number  ! Convert character to integer

  ! -----------------------------------------------------------------
  ! Get Granule Start and End time. This is provided by the
  ! Scheduler and has to be checked against the L1B Metadata
  ! fields for RangeBeginningTime and RangeEndingTime.
  ! -----------------------------------------------------------------
  errstat = PGS_PC_GetConfigData ( granule_s_lun, pcf_granule_s_time )
  errstat = PGS_PC_GetConfigData ( granule_e_lun, pcf_granule_e_time )  
  IF ( PGS_SMF_TestStatusLevel(errstat) /= pgs_smf_mask_lev_s ) THEN
     errstat = OMI_SMF_setmsg ( OMSAO_W_GETLUN, 'PGE_MOL_ID', modulename, 0 )
     pge_error_status = pge_errstat_warning
  END IF

  ! -------------------------------------------------------------------------
  ! Get the SAO PGE Name string. This string is converted to the PGE
  ! index number (10, 11, 12), which in turn is used to access LUNs and other
  ! elements that are PGE specific. All of this has to be done here, because
  ! we require the PGE index to identify LUNs and other PGE specific items.
  ! -------------------------------------------------------------------------
  errstat = PGS_PC_GetConfigData ( pge_mol_id, tmpchar )
  IF ( PGS_SMF_TestStatusLevel(errstat) /= pgs_smf_mask_lev_s ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_getlun, 'PGE_MOL_ID', modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

  ! ---------------------------------------------------------------
  ! Translate molecule string into index (required for LUN look-up)
  ! ---------------------------------------------------------------
  CALL get_pge_ident ( TRIM(ADJUSTL(tmpchar)), pge_name, pge_idx, errstat )
  IF ( errstat /= pge_errstat_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_get_molindex, '', modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  ELSE
     IF ( verb_thresh_lev >= vb_lev_develop ) &
          errstat = OMI_SMF_setmsg (omsao_s_get_molindex, TRIM(ADJUSTL(pge_name)), modulename, 0)
  END IF

  ! --------------------------------------------------------------
  ! Read names of L1B input file from PCF
  ! --------------------------------------------------------------
  version = 1
  errstat = PGS_PC_GetReference (pge_l1b_radiance_lun, version, l1b_rad_filename)
  IF ( PGS_SMF_TestStatusLevel(errstat) /= pgs_smf_mask_lev_s .OR. &
       LEN(TRIM(ADJUSTL(l1b_rad_filename))) == 0 ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_getlun, 'PGE_L1B_RADIANCE_LUN', modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

  version = 1
  errstat = PGS_PC_GetReference (pge_l1b_irradiance_lun, version, l1b_irrad_filename)
  IF ( PGS_SMF_TestStatusLevel(errstat) /= pgs_smf_mask_lev_s .OR. &
       LEN(TRIM(ADJUSTL(l1b_irrad_filename))) == 0 ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_getlun, 'PGE_L1B_IRRADIANCE_LUN', modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

  ! ------------------------------------------------------------------
  ! Read names of static input files from PCF. The first entry is the 
  ! algorithm control file (initial fitting values, etc.), the other
  ! entries are reference spectra.
  ! ------------------------------------------------------------------
  static_input_fnames = zerospec_string
  DO i = icf_idx, max_rs_idx
     version = 1
     errstat = PGS_PC_GetReference (pge_static_input_luns(i,pge_idx), version, tmpchar)
     tmpchar = TRIM(ADJUSTL(tmpchar)) ; strlen = LEN(TRIM(ADJUSTL(tmpchar)))
     IF ( PGS_SMF_TestStatusLevel(errstat) /= pgs_smf_mask_lev_s .OR. strlen == 0 ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_getlun, 'PGE_STATIC_INPUT_LUNS(i,pge_idx)', modulename, 0)
        pge_error_status = pge_errstat_error; RETURN
     ELSE IF ( INDEX ( TRIM(ADJUSTL(tmpchar)), zerospec_string ) == 0 ) THEN
        static_input_fnames(i) = TRIM(ADJUSTL(tmpchar))
     END IF
  END DO

  ! -----------------------------------------------------------------------
  ! Read name of AMF table file. Store as last entry of static input files.
  ! Remember that a missing AMF table is not a fatal problem, since in that
  ! case the slant columns will be written out.
  ! -----------------------------------------------------------------------
  version = 1
  errstat = PGS_PC_GetReference (pge_static_input_luns(amf_idx,pge_idx), version, tmpchar)
  tmpchar = TRIM(ADJUSTL(tmpchar)) ; strlen = LEN(TRIM(ADJUSTL(tmpchar)))
  IF ( ( PGS_SMF_TestStatusLevel(errstat) /= pgs_smf_mask_lev_s         ) .OR. &
       ( strlen == 0                                                    ) .OR. &
       ( INDEX ( TRIM(ADJUSTL(tmpchar)), zerospec_string ) /= 0 )       ) THEN
     IF ( verb_thresh_lev >= vb_lev_omidebug ) &
          errstat = OMI_SMF_setmsg (omsao_w_getlun, 'PGE_STATIC_INPUT_LUNS(amf_idx,pge_idx)', &
          modulename, 0)
     pge_error_status = pge_errstat_warning
     have_amftable = .FALSE.
  ELSE
     static_input_fnames(i) = TRIM(ADJUSTL(tmpchar))
     have_amftable = .TRUE.
  END IF

  ! ---------------------------
  ! Read name of L2 output file
  ! ---------------------------
  version = 1
  errstat = PGS_PC_GetReference (pge_l2_output_luns(pge_idx), version, l2_filename)
  IF ( PGS_SMF_TestStatusLevel(errstat) /= pgs_smf_mask_lev_s .OR. &
       LEN(TRIM(ADJUSTL(l2_filename))) == 0 ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_getlun, 'PGE_L2_OUTPUT_LUNS(pge_idx)', modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

  RETURN
END SUBROUTINE read_pcf_file
