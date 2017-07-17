PROGRAM BOREAS_main

  ! ************************************************************************
  !
  ! This is the main program of the SAO Product Generation Executives (PGEs)
  ! for  the Ozone Monitoring Istrument (OMI).
  !
  ! Authors: Thomas P. Kurosu, Kelly Chance
  !          Smithsonian Astrophysical Observatory
  !          60 Garden Street (MS 50)
  !          Cambridge, MA 02138 (USA)
  !
  !          EMail: tkurosu@cfa.harvard.edu
  !                 kchance@cfa.harvard.edu
  !
  ! ************************************************************************ 

  USE OMSAO_precision_module
  USE OMSAO_indices_module,   ONLY: icf_idx
  USE OMSAO_variables_module, ONLY: pge_idx, static_input_fnames
  USE OMSAO_gome_data_module, ONLY: gome_idx, omi_idx, instrument_idx
  USE OMSAO_errstat_module


  IMPLICIT NONE

  ! -------------------------
  ! Name of module/subroutine
  ! -------------------------
  CHARACTER (LEN=11), PARAMETER :: modulename = 'BOREAS_main'

  ! -------------------------
  ! OMI L1B related variables
  ! -------------------------
  INTEGER :: errstat, OMI_SMF_setmsg

  ! ------------------------------------------------------------------
  ! The general PGE error status variable. This is a relatively late 
  ! addition and is not used consistently in all routines yet. However,
  ! it is envisaged that it will be used ubiquitously througout the
  ! PGE once the PGE developer gets around to implementing it as such.
  ! ------------------------------------------------------------------
  INTEGER :: pge_error_status = pge_errstat_ok

  ! ----------------------------------------------------------------------------
  ! Some variables/parameters that are specific for the GOME way of doing things
  ! ----------------------------------------------------------------------------
  INTEGER, PARAMETER       :: fcunit = 11, specunit = 12, amfunit = 13
  INTEGER, PARAMETER       :: l1funit = 21, l2funit = 22, calunit = 23
  CHARACTER (LEN=maxchlen) :: l1_inputs_fname_sol, l1_inputs_fname_rad, l2_output_fname

  ! *****************
  ! >>> Start Of Code
  ! *****************

  instrument_idx = 0

  ! ------------------------------------------------
  ! Assign fitting input control file
  static_input_fnames(icf_idx) = 'INP/BOREAS.inp'
  ! ------------------------------------------------

  ! --------------------------------
  ! Make PGE write STD/IO unbuffered
  CALL unbufferSTDout()
  ! --------------------------------

  ! ------------------------------------------------------------------------------
  ! Read fitting conrol parameters from input file
  CALL gome_read_fitting_control_file (fcunit, static_input_fnames(icf_idx), &
       instrument_idx, l1_inputs_fname_sol, l1_inputs_fname_rad, l2_output_fname, &
       pge_error_status )

  ! ------------------------------------------------------------------------------
  IF ( pge_error_status >= pge_errstat_error ) GOTO 666

  ! -----------------------------------------------
  ! Different instruments require different actions
  ! -----------------------------------------------
  SELECT CASE ( instrument_idx )

  CASE ( gome_idx)     ! GOME fitting branch

     ! -------------------------------------------------------
     ! Read reference spectra
     CALL gome_read_reference_spectra ( specunit, pge_error_status )
     ! -------------------------------------------------------
     IF ( pge_error_status >= pge_errstat_error ) GOTO 666


     ! ----------------------------------
     ! Read AMF table
!     CALL gome_read_amf_table ( amfunit )
!    ! ----------------------------------
!     IF ( pge_error_status >= pge_errstat_error ) GOTO 666

     ! ----------------------------------------------------------------
     ! GOME fitting
     CALL gome_pge_fitting_process (l1funit, l2funit, calunit, l1_inputs_fname_sol,&
          l1_inputs_fname_rad, l2_output_fname, pge_error_status )
     ! ----------------------------------------------------------------

  CASE ( omi_idx )     ! OMI fitting branch


     ! -------------------------------------
     ! Read PCF file
!     CALL read_pcf_file ( pge_error_status )
!     ! -------------------------------------
!     IF ( pge_error_status >= pge_errstat_error ) GOTO 666

!     ! ----------------------------------------------------------
!     ! Read fitting conrol parameters from input file
!     !CALL read_fitting_control_file ( pge_idx, pge_error_status )
!     ! ----------------------------------------------------------
!     IF ( pge_error_status >= pge_errstat_error ) GOTO 666

!     ! -----------------------------
!     ! Read AMF table
!     CALL read_amf_table ( pge_idx )
!     ! -----------------------------
!     IF ( pge_error_status >= pge_errstat_error ) GOTO 666

     ! ---------------------------------------------------------
     ! Where all the work is done
     CALL omi_pge_fitting_process  ( pge_idx, pge_error_status )
     ! ---------------------------------------------------------
     IF ( pge_error_status >= pge_errstat_error ) GOTO 666
  END SELECT
  
  ! ------------------------------------
  ! Write END_OF_RUN message to log file
  ! ------------------------------------
666 SELECT CASE ( pge_error_status )
  CASE ( pge_errstat_ok )
     ! ----------------------------------------------------------------
     ! PGE execution completed successfully. All is well.
     ! ----------------------------------------------------------------
     errstat = OMI_SMF_setmsg ( OMSAO_S_ENDOFRUN, '', modulename, 0 )
!     STOP 0
  CASE ( pge_errstat_warning )
     ! ----------------------------------------------------------------
     ! PGE execution raised non-terminal warnings. Nothing serious, we
     ! hope, so execution completed but with a non-zero exit status.
     ! ----------------------------------------------------------------
     errstat = OMI_SMF_setmsg ( OMSAO_W_ENDOFRUN, '', modulename, 0 )
!     STOP 1
  CASE ( pge_errstat_error )
     ! ----------------------------------------------------------------
     ! PGE execution encountered an error that lead to termination.
     ! ----------------------------------------------------------------
     errstat = OMI_SMF_setmsg ( OMSAO_E_ENDOFRUN, '', modulename, 0 )
!     STOP 1
  CASE DEFAULT
     ! ----------------------------------------------------------------
     ! If we ever reach here, then PGE_ERRSTAT has been set to a funny
     ! value. This should never happen, but we buffer this case anyway.
     ! ----------------------------------------------------------------
     errstat = OMI_SMF_setmsg ( OMSAO_U_ENDOFRUN, '', modulename, 0 )
!     STOP
  END SELECT


!  STOP

END PROGRAM BOREAS_main
