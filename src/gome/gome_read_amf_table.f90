SUBROUTINE gome_read_amf_table ( amfunit )

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: amf_idx, amftab_ang_idx, amftab_amf_idx, &
       n_amftab_ang_max, n_amftab_dim_max, pge_static_input_luns
  USE OMSAO_parameters_module,       ONLY: lm_start_of_table, vb_lev_omidebug, vb_lev_develop
  USE OMSAO_variables_module,        ONLY: &
       n_amftab_dim, n_amftab_ang, amf_table, have_amftable, amf_esza_min, amf_esza_max, &
       static_input_fnames, verb_thresh_lev
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER, INTENT (IN) :: amfunit


  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=19), PARAMETER :: modulename = 'gome_read_amf_table'

  ! ========================
  ! Error handling variables
  ! ========================
  INTEGER :: errstat, version, file_read_stat

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                  :: i, j, funit, omi_lun
  CHARACTER (LEN=maxchlen) :: lastline

  ! =================================
  ! External OMI and Toolkit routines
  ! =================================
  INTEGER :: OMI_SMF_setmsg

  ! ----------------------------------------------------------------------------------
  ! First check whether HAVE_AMFTABLE=.TRUE. from READ_PCF_FILE. HAVE_AMFTABLE=.FALSE.
  ! means that we failed to get the AMF table file from the PCF and thus can't 
  ! determine the AMFs. HAVE_AMFTABLE will also be set .FALSE. if any error occurs in
  ! this routine, in which case no slant-to-vertical conversion of the column amount 
  ! will be performed.
  ! ----------------------------------------------------------------------------------
  IF ( .NOT. have_amftable ) RETURN

  ! -----------------------
  ! Open AMF table file
  ! -----------------------
  OPEN ( UNIT=amfunit, FILE=TRIM(ADJUSTL(static_input_fnames(amf_idx))), &
       STATUS='OLD', IOSTAT=errstat )
  IF ( errstat /= pge_errstat_ok ) THEN
     errstat = OMI_SMF_setmsg ( omsao_w_open_amftable_file, &
          TRIM(ADJUSTL(static_input_fnames(amf_idx))), modulename, verb_thresh_lev )
     have_amftable = .FALSE.  ;  CLOSE ( amfunit ) ; RETURN
  END IF

  ! ----------------------------
  ! Initialize output quantities
  ! ----------------------------
  n_amftab_dim = 0      ; n_amftab_ang = 0
  amf_esza_min = 0.0_dp ; amf_esza_max = 0.0_dp

  ! --------------------------------------
  ! Skip comments header to start of table
  ! --------------------------------------
  CALL skip_to_filemark ( funit, lm_start_of_table, lastline, file_read_stat )
  IF ( file_read_stat /= file_read_ok .AND. &
       verb_thresh_lev >= vb_lev_omidebug ) THEN
     errstat = OMI_SMF_setmsg (omsao_w_read_amftable_file, '', modulename, 0)
     have_amftable = .FALSE.  ;  CLOSE ( amfunit ) ; RETURN
  END IF

  ! -----------------------------------------------------------------
  ! Read table dimensions and MIN/MAX of effective solar zenith angle
  ! -----------------------------------------------------------------
  READ (UNIT=funit, FMT=*, IOSTAT=errstat) n_amftab_ang, n_amftab_dim, amf_esza_min, amf_esza_max

  ! ------------------------------------------
  ! Check whether we can proceed to read table
  ! ------------------------------------------
  IF ( (errstat      /= file_read_ok)     .OR. &   ! READ successful
       (n_amftab_dim <= 0)                .OR. &   ! Table dimensions >  0
       (n_amftab_dim >  n_amftab_dim_max) .OR. &   ! Table dimensions <= maximum table dim
       (n_amftab_ang <= 0)                .OR. &   ! Table dimensions >  0
       (n_amftab_ang >  n_amftab_ang_max) .OR. &   ! Table dimensions <= maximum table dim
       (amf_esza_min <  0.0_dp)           .OR. &   ! Minimum ESZA >= 0
       (amf_esza_max <= 0.0_dp)           .OR. &   ! Maximum ESZA >= 0
       (amf_esza_max <  amf_esza_min)          &   ! Maximum ESZA >= Miminum ESZA
       ) THEN
     IF ( verb_thresh_lev >= vb_lev_omidebug ) &
          errstat = OMI_SMF_setmsg (omsao_w_read_amftable_file, '', modulename, 0)
     have_amftable = .FALSE.  ;  CLOSE ( amfunit ) ; RETURN
  ELSE
     DO i = 1, n_amftab_ang
        READ (UNIT=funit, FMT=*, IOSTAT=errstat) amf_table(i,1:n_amftab_dim)
        IF ( i /= n_amftab_ang .AND. errstat /= 0 .AND. &
             verb_thresh_lev >= vb_lev_omidebug ) THEN
           errstat = OMI_SMF_setmsg (omsao_w_read_amftable_file, '', modulename, 0)
           have_amftable = .FALSE. ; CLOSE ( amfunit ) ; RETURN
        END IF
     END DO
  END IF

  ! ------------------------------------
  ! Report successful reading of spectra
  ! ------------------------------------
  IF ( verb_thresh_lev >= vb_lev_develop ) &
       errstat = OMI_SMF_setmsg (omsao_s_read_amftable_file, '', modulename, 0)

  CLOSE ( amfunit )
  RETURN
END SUBROUTINE gome_read_amf_table
