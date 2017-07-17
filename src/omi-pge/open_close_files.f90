SUBROUTINE open_close_files ( omi_lun, fname, openclose, act, omi_funit )

  ! ****************************************************
  ! OPEN or CLOSE files. Report either a warning or and 
  ! error and stop if action ACT fails or is unknown.
  ! ****************************************************

  USE OMSAO_parameters_module, ONLY: vb_lev_omidebug, vb_lev_develop
  USE OMSAO_variables_module,  ONLY: verb_thresh_lev
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,           INTENT (IN) :: omi_lun
  CHARACTER (LEN=*), INTENT (IN) :: fname, openclose, act

  ! ==================
  ! Modified variables
  ! ==================
  INTEGER, INTENT (INOUT) :: omi_funit

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=16), PARAMETER :: modulename = 'open_close_files'

  ! ========================
  ! Error handling variables
  ! ========================
  INTEGER :: errstat, version

  ! =================================
  ! External OMI and Toolkit routines
  ! =================================
  INTEGER :: OMI_SMF_setmsg, pgs_smf_teststatuslevel, pgs_io_gen_openf, pgs_io_gen_closef

  version = 1
  SELECT CASE ( openclose )
  CASE ( 'OPEN' )
     SELECT CASE ( act )
     CASE ( 'READ' )
        errstat = PGS_IO_GEN_OPENF ( omi_lun, PGSd_IO_Gen_RSeqFrm, 0, omi_funit, version )
        IF ( PGS_SMF_TESTSTATUSLEVEL(errstat) /= pgs_smf_mask_lev_s ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_open_rancil_file, &
                TRIM(ADJUSTL(fname)), modulename, 0) ; STOP 1
        ELSE
           IF ( verb_thresh_lev >= vb_lev_develop ) &
                errstat = OMI_SMF_setmsg (omsao_s_open_rancil_file, &
                TRIM(ADJUSTL(fname)), modulename, 0)
        END IF
     CASE ( 'WRITE' )
        errstat = PGS_IO_GEN_OPENF ( omi_lun, PGSd_IO_Gen_WSeqFrm, 0, omi_funit, version )
        IF ( PGS_SMF_TESTSTATUSLEVEL(errstat) /= pgs_smf_mask_lev_s ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_open_wancil_file, &
                TRIM(ADJUSTL(fname)), modulename, 0) ; STOP 1
        ELSE
           IF ( verb_thresh_lev >= vb_lev_develop ) &
                errstat = OMI_SMF_setmsg (omsao_s_open_wancil_file, &
                TRIM(ADJUSTL(fname)), modulename, 0)
        END IF
     CASE DEFAULT
        errstat = OMI_SMF_setmsg (omsao_e_open_uancil_file, &
             TRIM(ADJUSTL(act))//'|'//TRIM(ADJUSTL(fname)), modulename, 0) ; STOP
     END SELECT
  CASE ( 'CLOSE' )
     errstat = PGS_IO_GEN_CLOSEF ( omi_funit )
     IF ( PGS_SMF_TESTSTATUSLEVEL(errstat) /= pgs_smf_mask_lev_s .AND. &
          verb_thresh_lev >= vb_lev_omidebug ) THEN
        errstat = OMI_SMF_setmsg (omsao_w_close_ancil_file, &
             TRIM(ADJUSTL(fname)), modulename, 0)
     ELSE
        IF ( verb_thresh_lev >= vb_lev_develop ) &
             errstat = OMI_SMF_setmsg (omsao_s_close_ancil_file, &
             TRIM(ADJUSTL(fname)), modulename, 0)
     END IF
  CASE DEFAULT
     errstat = OMI_SMF_setmsg (omsao_e_uaction_ancil_file, &
          TRIM(ADJUSTL(act))//'|'//TRIM(ADJUSTL(fname)), modulename, 0); STOP 1
  END SELECT


  RETURN
END SUBROUTINE open_close_files
