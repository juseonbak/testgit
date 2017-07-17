MODULE BOREAS_gome2_data_module

  ! Module to define GOME-2 data.  
  ! Based on a version by Trevor Beck (NOAA).
  ! C Nowlan, cnowlan@cfa.harvard.edu

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxchlen
  ! The following modules are defined in library created from g2_l1b_reader routines.
  ! They contain routines for reading GOME-2 level 1B data from binary files.
  USE gome_mdr_1b_earthshine_module
  USE gome_viadr_smr_module
  IMPLICIT NONE


  ! =================
  ! Define CLOUD type
  ! =================
  TYPE (cloud_struct)     :: CLOUD	

  TYPE (polss_struct)     :: POL_SS
  TYPE (pcd_basic_struct) :: PCD_BASIC

  ! =================
  ! PMD Parameters
  ! =================
  ! GOME-2 PMD numbers are different from GOME-1
  INTEGER, PARAMETER :: n_gome2_pmd = 4 
  INTEGER, PARAMETER :: n_gome2_pmd_pix = 64 ! FIXME, not certain of values for G2.

  ! =============================
  ! Scan parameters and variables
  ! =============================
  ! max cross-track for bands 1-6, not including POL bands 6 through 10.
  INTEGER, PARAMETER :: max_xtrack = 32 
  INTEGER :: gome2_curr_ixtrack ! current cross track index
  INTEGER :: gome2_curr_nxtrack ! current number xtracks per scanline.
  INTEGER, DIMENSION(max_xtrack) :: xtrack_indices ! for current scan, 1..4, 1..8, 1..16, or 1..32
  INTEGER :: i
  INTEGER, PARAMETER, DIMENSION(max_xtrack) :: xtrack_index = (/(i,i=1,max_xtrack)/)
  INTEGER :: n_gome2_scans
  
  ! ===============
  ! Main data types
  ! ===============
  TYPE(gome_mdr_1b_earthshine_struct) :: g2_earthview
  TYPE(gome_mdr_1b_earthshine_struct) ::g2_earthview_nextscan !For indexing issue in R_O data
  TYPE(gome_viadr_smr_struct) :: g2_smr
  INTEGER :: n_gome2_ev

  
  REAL(KIND=sp) :: integration_time
  REAL(KIND=dp) :: gome2_cloud_press, gome2_cloud_frac

  REAL(KIND=dp) :: gome2_surf_elevation

  INTEGER(KIND=4) :: gome2_utc_day
  INTEGER(KIND=4) :: gome2_utc_year
  INTEGER(KIND=8) :: gome2_utc_millisec
  INTEGER :: gome2_band  ! current earthview band
  INTEGER :: gome2_smr_band !current SMR band.
  INTEGER :: gome2_lastscan=0
  INTEGER , PARAMETER :: gome2_max_negrad = 5 !if more rads than this are neg. or zero. then BADPIX
  REAL(KIND=dp) :: gome2_relazm

  INTEGER, PARAMETER :: maxscan=1000  ! max number of along track scans, typically 600 scans per orbit
  INTEGER, DIMENSION(maxscan) :: scan_start_index
  INTEGER, DIMENSION(maxscan) :: nxtrack_per_scan

  INTEGER :: glint_flg


END MODULE BOREAS_gome2_data_module
