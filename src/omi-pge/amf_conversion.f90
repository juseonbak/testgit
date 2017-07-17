SUBROUTINE amf_conversion ( n_arr_pix, sza, los, sol_zen_eff, amf, amfgeo, have_amf )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: deg2rad, rad2deg
  USE OMSAO_variables_module,  ONLY: szamax, amf_esza_min, amf_esza_max, have_amftable

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=14), PARAMETER :: modulename = 'amf_conversion'

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,                               INTENT (IN) :: n_arr_pix
  REAL (KIND=dp), DIMENSION (n_arr_pix), INTENT (IN) :: sza, los

  ! ----------------
  ! Output variables
  ! ----------------
  REAL (KIND=dp), INTENT (OUT) :: sol_zen_eff, amf, amfgeo
  LOGICAL,        INTENT (OUT) :: have_amf

  ! ------------------
  ! External functions
  ! ------------------
  REAL (KIND=dp) :: esza2amf

  have_amf = .TRUE.
  
  IF ( ALL (ABS(sza(1:n_arr_pix)) <= szamax) .AND. &       ! Solar zenith angles within bounds
       ALL (ABS(los(1:n_arr_pix)) <= szamax)       ) THEN  ! LineOfSight zenith angles within bounds

     ! -------------------------------------------------------------------------------
     ! Compute the geometric AMF as the average over all geometric AMFs for that pixel
     ! -------------------------------------------------------------------------------
     amfgeo = ( SUM(1.0_dp/COS(deg2rad*sza)) + SUM(1.0_dp/COS(deg2rad*los)) ) &
          / REAL(n_arr_pix, KIND=dp)

     ! ------------------------------------------------
     ! Compute effective solar zenith angle from AMFGEO
     ! ------------------------------------------------
     sol_zen_eff = ACOS ( 1.0_dp/ (amfgeo - 1.0_dp) ) * rad2deg

     ! -------------------------------------------------------------------------
     ! Table look-up of AMF if ESZA is within table bounds, otherwise AMF=AMFGEO
     ! -------------------------------------------------------------------------
     IF ( have_amftable                .AND. &
          sol_zen_eff <= amf_esza_max  .AND. &
          sol_zen_eff >= amf_esza_min           ) THEN
        amf = esza2amf ( sol_zen_eff )
     ELSE
        amf = amfgeo
     END IF

  ELSE
     amfgeo = 1.0_dp ; amf = 1.0_dp ; have_amf = .FALSE.
  END IF

  RETURN
END SUBROUTINE amf_conversion

FUNCTION esza2amf ( sol_zen_eff )  RESULT ( amf )

  ! --------------------------------------------------------------------------
  ! Function to convert an Effective Solar Zenith Angle to and Air Mass Factor
  ! --------------------------------------------------------------------------

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: amftab_ang_idx, amftab_amf_idx
  USE OMSAO_parameters_module, ONLY: vb_lev_omidebug
  USE OMSAO_variables_module,  ONLY: n_amftab_ang, amf_table, verb_thresh_lev
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  REAL (KIND=dp), INTENT(IN) :: sol_zen_eff

  ! ---------------
  ! Return variable
  ! ---------------
  REAL (KIND=dp) :: amf

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER :: errstat

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=8), PARAMETER :: modulename = 'esza2amf'

  CALL interpolation ( &
       n_amftab_ang, amf_table(1:n_amftab_ang,amftab_ang_idx), &
       amf_table(1:n_amftab_ang,amftab_amf_idx), 1, sol_zen_eff, amf, errstat )
  ! -----------------------------------------------------------------------------------------
  ! Set AMF to Zero if ERRSTAT is not 0, i.e., if the interpolation has encountered an error.
  ! -----------------------------------------------------------------------------------------
  IF ( errstat /= pge_errstat_ok .AND. &
       verb_thresh_lev >= vb_lev_omidebug ) THEN
     errstat = OMI_SMF_setmsg (omsao_w_interpol, modulename, '', 0) ; amf = 0.0_dp
  END IF

  RETURN
END FUNCTION esza2amf
