MODULE OMSAO_metadata_module

  ! ==============================
  ! Module for Metadata parameters
  ! ==============================

  USE OMSAO_precision_module,  ONLY: i4, r8
  USE OMSAO_indices_module,    ONLY: n_sao_pge, pge_l1b_radiance_lun, pge_l1b_irradiance_lun
  USE OMSAO_parameters_module, ONLY: maxchlen

  IMPLICIT NONE

  INCLUDE "PGS_MET.f"
  INCLUDE "PGS_PC.f"

  ! ==============================================================
  ! The following are lists of MetaData Objects from the MCF that
  ! are either obtained from the PCF or have to be set in the PGE.
  ! ==============================================================
  !
  ! PCF:  ReprocessingActual                string
  !       OrbitNumber                       integer
  !       PGEVersion                        string
  !
  ! PGE:  LocalGranuleID                    string
  !       AutomaticQualityFlagExplanation   string 
  !       AutomaticQualityFlag              string
  !       OperationalQualityFlagExplanation string
  !       OperationalQualityFlag            string
  !       QAPercentMissingData              integer
  !       QAPercentOutofBoundsData          integer
  !       QAStats                           integer
  !       ParameterName                     string
  !       EquatorCrossingDate               date
  !       EquatorCrossingTime               time
  !       EquatorCrossingLongitude          double
  !       InputPointer                      string
  !       RangeEndingDate                   date
  !       RangeEndingTime                   time
  !       RangeBeginningDate                date
  !       RangeBeginningTime                time
  !       AdditionalAttributeName           string
  !       InformationContent                string
  !       (?) AssociatedSensorShortName
  !       (?) AssociatedPlatformShortName
  !       (?) AssociatedInstrumentShortName
  ! ==============================================================

  ! -----------------------------------------------------
  ! Lists of MetaData Object field names, with dimensions
  ! --------------------------------------------------------------------------
  ! NOTE: F90 intrinsically, the character fields are truncated to the length
  !       of the first (or shortest?) entry. Therefore blanks should be filled
  !       in to bring all fields to the same length.
  ! --------------------------------------------------------------------------
  INTEGER, PARAMETER :: n_mdata_l1b_str = 14, n_mdata_l1b_int = 4, n_mdata_l1b_dbl = 1
  INTEGER, PARAMETER :: n_mdata_mcf_str =  4, n_mdata_mcf_int = 1, n_arcmdata_mcf_str = 2

  ! STRING Objects in L1B file
  CHARACTER (LEN=PGSd_MET_NAME_L), DIMENSION (n_mdata_l1b_str), PARAMETER :: &
       mdata_l1b_objects_str = (/             &
       "LocalGranuleID                   ",   &
       "AutomaticQualityFlag.1           ",   &
       "AutomaticQualityFlagExplanation.1",   &
       "EquatorCrossingDate.1            ",   &
       "EquatorCrossingTime.1            ",   &
       "RangeBeginningDate               ",   &
       "RangeEndingDate                  ",   &
       "RangeBeginningTime               ",   &
       "RangeEndingTime                  ",   &
       "InputPointer                     ",   &
       "ShortName                        ",   &
       "AssociatedSensorShortName.1      ",   &
       "AssociatedPlatformShortName.1    ",   &
       "AssociatedInstrumentShortName.1  "     /)
!       "ParameterName.1                  ",   &

  ! INTEGER Objects in L1B file
  CHARACTER (LEN=PGSd_MET_NAME_L), DIMENSION (n_mdata_l1b_int), PARAMETER:: &
       mdata_l1b_objects_int = (/             &
       "QAPercentMissingData.1    ",          &
       "QAPercentOutofBoundsData.1",          &
       "OrbitNumber.1             ",          &
       "VersionID                 "           /)

  ! DOUBLE Objects in L1B file
  CHARACTER (LEN=PGSd_MET_NAME_L), DIMENSION (n_mdata_l1b_dbl), PARAMETER:: &
       mdata_l1b_objects_dbl = (/ "EquatorCrossingLongitude.1" /)


  ! STRING Objects in MCF file (Inventory Metadata)
  CHARACTER (LEN=PGSd_MET_NAME_L), DIMENSION (n_mdata_mcf_str), PARAMETER :: &
       mdata_mcf_objects_str = (/           &
       "ShortName                      ",   &
       "AssociatedSensorShortName.1    ",   &
       "AssociatedPlatformShortName.1  ",   &
       "AssociatedInstrumentShortName.1"     /)

  ! STRING Objects in MCF file (Archive Metadata)
  CHARACTER (LEN=PGSd_MET_NAME_L), DIMENSION (n_arcmdata_mcf_str), PARAMETER :: &
       arcmdata_mcf_objects_str = (/ &
       "LongName              ", &
       "ESDTDescriptorRevision"   /)

  ! INTEGER Objects in MCF file
  CHARACTER (LEN=PGSd_MET_NAME_L), DIMENSION (n_mdata_mcf_int), PARAMETER:: &
       mdata_mcf_objects_int = (/ "VersionID" /)


  ! --------------------------------------------------------------------
  ! A number of PSA fields to be copied from the L1b file (after K.Yang)
  ! --------------------------------------------------------------------
  CHARACTER(LEN=1000) ::  psa_2bcopied_from_l1b = &
       "NrMeasurements"               // &
       "NrZoom"                       // &
       "NrSpatialZoom"                // &
       "NrSpectralZoom"               // &
       "ExpeditedData"                // &
       "SouthAtlanticAnomalyCrossing" // &
       "SpacecraftManeuverFlag"       // &
       "SolarEclipse"                 // &
       "MasterClockPeriods"           // &
       "ExposureTimes"                // &
       "PathNr"                       // &
       "StartBlockNr"                 // &
       "EndBlockNr"

  ! ------------------------------------
  ! Parameter values set by the SAO PGEs
  ! ------------------------------------
  INTEGER (KIND=I4), PARAMETER :: n_l2_output_paras = 1
  CHARACTER (LEN=80), DIMENSION (n_l2_output_paras) :: &
       pge_parameter_names = (/ "ColumnAmount"   /)

!  INTEGER (KIND=I4), PARAMETER :: n_l2_output_paras = 11, n_l2_op_mol = 4
!  CHARACTER (LEN=80), DIMENSION (n_l2_output_paras) :: pge_parameter_names = (/ &
!       "ColumnAmount        ", &  ! Molecule name to be appended to this field
!       "ColumnUncertainty   ", &  ! Molecule name to be appended to this field
!       "FitConversionFlag   ", &  ! Molecule name to be appended to this field
!       "AirMassFactor       ", &  ! Molecule name to be appended to this field
!       "AirMassFactorGeo    ", &
!       "OMISlitWidthFit     ", &
!       "Latitude            ", &
!       "Longitude           ", &
!       "SolarZenithAngle    ", &
!       "ViewingZenithAngle  ", &
!       "RelativeAzimuthAngle"   /)

  ! -------------------------------------------------------------
  ! Arrays and variables to hold Metatdata values read from files
  ! -------------------------------------------------------------
  CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L), DIMENSION(n_mdata_l1b_str)    :: mdata_l1b_fieldval_str
  CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L), DIMENSION(n_mdata_mcf_str)    :: mdata_mcf_fieldval_str
  CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L), DIMENSION(n_arcmdata_mcf_str) :: arcmdata_mcf_fieldval_str
  INTEGER   (KIND=i4),                       DIMENSION(n_mdata_l1b_int)    :: mdata_l1b_fieldval_int
  INTEGER   (KIND=i4),                       DIMENSION(n_mdata_mcf_int)    :: mdata_mcf_fieldval_int
  REAL      (KIND=r8),                       DIMENSION(n_mdata_l1b_dbl)    :: mdata_l1b_fieldval_dbl
  INTEGER   (KIND=i4)                                                      :: mcf_versionid

  ! -----------------------------------------------------------------------
  ! Some fixed strings. They will be used as short hand for special task,
  ! like extracting the date in the required format from the Metadata field
  ! 'RageBeginningDate', or are recurring strings, like 'CoreMetadata.xxx'
  ! that we don't want to spell out all the time they appear.
  ! -----------------------------------------------------------------------
  CHARACTER (LEN=12), PARAMETER :: cmd_str  = "CoreMetadata"
  CHARACTER (LEN=18), PARAMETER :: rbd_str  = "RangeBeginningDate"
  CHARACTER (LEN=18), PARAMETER :: rbt_str  = "RangeBeginningTime"
  CHARACTER (LEN=27), PARAMETER :: apsn_str = "AssociatedPlatformShortName"
  CHARACTER (LEN=29), PARAMETER :: aisn_str = "AssociatedInstrumentShortName"
  CHARACTER (LEN= 9), PARAMETER :: vid_str  = "VersionID"
  CHARACTER (LEN= 9), PARAMETER :: sn_str   = "ShortName"
  CHARACTER (LEN=14), PARAMETER :: lgid_str = "LocalGranuleID"
  CHARACTER (LEN=13), PARAMETER :: par_str  = "ParameterName"

  ! -------------------------------------------------------
  ! Parameters and variables for L2 Metadata initialization
  ! -------------------------------------------------------
  INTEGER (KIND=i4),                              PARAMETER :: n_lun_inp_point = 2
  INTEGER (KIND=i4), DIMENSION (n_lun_inp_point), PARAMETER :: &
       lun_input_pointer = (/ pge_l1b_radiance_lun, pge_l1b_irradiance_lun /)

  CHARACTER (LEN=PGSd_MET_GROUP_NAME_L), DIMENSION (PGSd_MET_NUM_OF_GROUPS) :: mcf_groups
  CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX)                                     :: l2_local_granule_id

  ! -------------------------------------------------
  ! STRING variables for various purposes:
  !  * ShortName from MCF
  !  * Granule Start/End Time from PCF
  !    (will be compared against RangeBeginningTime
  !     and RangeEndTime)
  ! -------------------------------------------------
  CHARACTER (LEN=maxchlen) :: mcf_shortname
  CHARACTER (LEN=maxchlen) :: pcf_granule_s_time, pcf_granule_e_time


  ! --------------------------------------------------
  ! Some variables that will have to be set by the PGE
  ! --------------------------------------------------
  INTEGER (KIND=i4) :: QAPercentMissingData, QAPercentOutofBoundsData, NumberOfInputSamples

END MODULE OMSAO_metadata_module
