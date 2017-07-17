
MODULE OMSAO_omidata_module

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxchlen, n_rad_winwav, maxwin, max_fit_pts, max_ring_pts, mrefl
  USE OMSAO_indices_module,    ONLY: n_max_fitpars, max_rs_idx, max_calfit_idx, o3_t1_idx, o3_t3_idx, spc_idx
  IMPLICIT NONE

  ! maximum number of channels (UV-1 and UV-2) required for OMI ozone profile retrievals
  INTEGER, PARAMETER :: mswath  = 2  
  INTEGER            :: nswath           ! actually number of swath used
  LOGICAL            :: do_xbin, do_ybin ! binning across the track and along the track
  INTEGER            :: nxbin, nybin , ncoadd    
  LOGICAL            :: zoom_mode
  INTEGER, PARAMETER :: zoom_p1 = 16, zoom_p2 = 45

  ! ------------------------------------------------------------
  ! Boundary wavelengths (approximate) for UV-2 and VIS channels
  ! ------------------------------------------------------------
  REAL (KIND=r4), DIMENSION(mswath),  PARAMETER :: upper_wvls = (/310.0, 387.0/), &
                                                   lower_wvls = (/260.0, 310.0/)
  REAL (KIND=dp), DIMENSION(mswath)             :: reduce_ubnd, reduce_lbnd, retubnd, retlbnd

  ! ---------------------------------------
  ! Minimum OMI spectral resolution (in nm)
  ! ---------------------------------------
  REAL (KIND=r8), PARAMETER :: omi_min_specres = 0.5_r8

  CHARACTER (LEN=5)         :: orbc, orbcsol
  CHARACTER (LEN=9)         :: omiraddate
  INTEGER (KIND=i4)         :: orbnum, orbnumsol, omisol_version

  ! ---------------------------------
  ! Maximum OMI data/swath dimensions
  ! ---------------------------------
  INTEGER (KIND=i4), PARAMETER :: &
       ntimes_max     = 2000, nxtrack_max  =  60, nwavel_max  =  max_fit_pts, &
       nwavelcoef_max =    5, nlines_max   = 100

  ! --------------------------------------------------
  ! Parameters defined by the NISE snow cover approach
  ! --------------------------------------------------
  INTEGER (KIND=i2), PARAMETER :: &
       NISE_snowfree =   0, NISE_allsnow = 100, NISE_permice = 101, NISE_drysnow = 103, &
       NISE_ocean    = 104, NISE_suspect = 125, NISE_error   = 127

  ! --------------------------------------------------------------------
  ! A diagnostic array that shows how the AMF was computed, with values
  ! that indicate missing cloud products, glint, and geometric or no AMF
  ! --------------------------------------------------------------------
  INTEGER (KIND=i2), PARAMETER :: &
       omi_cfr_addmiss = 1000, omi_ctp_addmiss = 2000, omi_glint_add = 10000, &
       omi_geo_amf = -1, omi_oobview_amf = -2
  INTEGER (KIND=i2), DIMENSION (nxtrack_max,0:nlines_max-1) :: amf_diagnostic

  ! retrieval status for a block
  INTEGER, DIMENSION (nxtrack_max, 0:nlines_max-1)                       :: omi_exitval, omi_initval
  REAL (KIND=dp), DIMENSION (nxtrack_max, 0:nlines_max-1, n_max_fitpars) :: omi_fitvar

  ! -----------------------
  ! Arrays for OMI L1b data
  ! -----------------------
  REAL    (KIND=r4), DIMENSION (0:nlines_max-1)                        :: omi_auraalt, omi_auralon, omi_auralat
  INTEGER (KIND=i2)                                                    :: omi_mflg
  REAL    (KIND=r8), DIMENSION (0:nlines_max-1)                        :: omi_time
  INTEGER (KIND=i2), DIMENSION (nxtrack_max, 0:nlines_max-1)           :: omi_radpix_errstat
  INTEGER (KIND=i2), DIMENSION (nxtrack_max)                           :: omi_solpix_errstat
  INTEGER (KIND=i2), DIMENSION (0:nlines_max-1)                        :: omi_radiance_errstat
  INTEGER (KIND=i2), DIMENSION (0:nlines_max-1)                        :: omi_saa_flag
  INTEGER (KIND=i2), DIMENSION (nxtrack_max,0:nlines_max-1)            ::  &
       omi_height, omi_geoflg, land_water_flg, glint_flg, snow_ice_flg
  INTEGER (KIND=i1), DIMENSION (nxtrack_max,0:nlines_max-1)            ::  omi_xtrackqflg, &
  rowanomaly_flg, waveshift_flg, blockage_flg, straysun_flg, strayearth_flg
  REAL    (KIND=r4), DIMENSION (nxtrack_max, 0:nlines_max-1)           ::  &
       omi_latitude, omi_longitude, omi_szenith, omi_sazimuth, omi_vzenith, omi_vazimuth, omi_eaza, omi_esca
  REAL    (KIND=r4), DIMENSION (nwavel_max,nxtrack_max,0:nlines_max-1) ::  &
       omi_radiance_spec, omi_radiance_prec, omi_radiance_wavl
  REAL (KIND=r8), DIMENSION(nxtrack_max)                               :: omi_solnorm
  REAL (KIND=r8), DIMENSION(nxtrack_max, 0:nlines_max-1)               :: omi_radnorm
  INTEGER (KIND=i2), DIMENSION (nwavel_max,nxtrack_max,0:nlines_max-1) :: omi_radiance_qflg, radwind
  INTEGER (KIND=i2), DIMENSION (nwavel_max,nxtrack_max) :: omi_irradiance_qflg, irradwind
  REAL    (KIND=r4), DIMENSION (nwavel_max,nxtrack_max) :: omi_irradiance_prec, &
       omi_irradiance_wavl, omi_irradiance_spec, omi_irrad_stray, omi_rad_stray

  ! radwind: the indices of selected irradiances in the L1b data
  ! irradwind: the indices of selected radiances in the selected irradiances

  REAL    (KIND=r4), DIMENSION (spc_idx, max_ring_pts, nxtrack_max)   :: omi_solspec_ring
  INTEGER, DIMENSION (nxtrack_max)              :: solring_lin, solring_uin, omi_nsolring, omi_solring_ndiv
  REAL    (KIND=r4), DIMENSION (spc_idx, mrefl, nxtrack_max)              :: omi_solspecr
  REAL    (KIND=r4), DIMENSION (spc_idx, mrefl, nxtrack_max, 0:nlines_max-1) :: omi_specr  

  ! ----------------------------------------
  ! Arrays for fitting and/or derived output
  ! ----------------------------------------
  !REAL    (KIND=r8), DIMENSION (nxtrack_max,0:nlines_max-1) :: omi_column_amount, omi_column_uncert
  !REAL    (KIND=r8), DIMENSION (nxtrack_max,0:nlines_max-1) :: omi_amf_col, omi_amf_geo, omi_fit_rms
  !REAL    (KIND=r8), DIMENSION (nxtrack_max,0:nlines_max-1) :: omi_ezenith, omi_radfit_chisq
  !REAL    (KIND=r4), DIMENSION (nxtrack_max,0:nlines_max-1) :: omi_razimuth
  !INTEGER (KIND=i2), DIMENSION (nxtrack_max,0:nlines_max-1) :: omi_fitconv_flag
  !INTEGER (KIND=i4), DIMENSION (nxtrack_max,0:nlines_max-1) :: temp_int
  !INTEGER (KIND=i2), DIMENSION (nxtrack_max,0:nlines_max-1) :: omi_itnum_flag, omi_main_qa_flag

  ! ----------------------------------------------------------------------------
  ! Correlations with main output product. Due to a bug in the HDF-EOS5 routines
  ! (non-TLCF implementation), STRING fields cannot be written to file directly.
  ! A work-around solution is to convert the CHARACTERs to INTEGERs. Thus the
  ! need for the additional array CORRELATION_NAMES_INT.
  ! ----------------------------------------------------------------------------
  !CHARACTER (LEN=maxchlen), DIMENSION (n_max_fitpars) :: correlation_names
  !CHARACTER (LEN=n_max_fitpars*maxchlen)              :: correlation_names_concat

  ! --------------------------------------------------------
  ! Ozone is a special case: We can have up to 3 temperatues
  ! --------------------------------------------------------
  !REAL (KIND=r8), DIMENSION (o3_t1_idx:o3_t3_idx, nxtrack_max,0:nlines_max-1) :: &
  !     omi_o3_amount, omi_o3_uncert

  ! ---------------------------------
  ! Dimensions for measurement swaths
  ! ---------------------------------
  ! nfxtrack: is the number of pixels after coadding UV-2
  INTEGER  :: ntimes, nxtrack, nwavel, ntimes_loop, nfxtrack, offset_line
  INTEGER, DIMENSION (maxwin, nxtrack_max)                :: omi_nsolpix
  INTEGER, DIMENSION (maxwin, nxtrack_max,0:nlines_max-1) :: omi_nradpix
  INTEGER, DIMENSION (nxtrack_max)                        :: omi_nwav_irrad
  INTEGER, DIMENSION (nxtrack_max,0:nlines_max-1)         :: omi_nwav_rad
  REAL (KIND=dp), DIMENSION(maxwin, nxtrack_max, 2)       :: omisol_winpix
  REAL (KIND=dp), DIMENSION(nxtrack_max, 2)               :: omisolr_winpix

  ! ---------------------------------------
  ! Swath attributes for measurement swaths
  ! ---------------------------------------
  INTEGER (KIND=i4), DIMENSION (nxtrack_max)                         :: n_omi_database_wvl
  INTEGER (KIND=i2), DIMENSION (nxtrack_max)                         :: &
       omi_solcal_itnum, omi_radcal_itnum, omi_radref_itnum,            &
       omi_solcal_xflag, omi_radcal_xflag, omi_radref_xflag
  REAL    (KIND=r8), DIMENSION (max_calfit_idx, nxtrack_max)         :: &
       omi_solcal_pars,  omi_radcal_pars,  omi_radref_pars
  REAL    (KIND=r8), DIMENSION (max_rs_idx, nwavel_max, nxtrack_max) :: omi_database
  REAL    (KIND=r8), DIMENSION (            nwavel_max, nxtrack_max) :: omi_database_wvl
  REAL    (KIND=r8), DIMENSION (nxtrack_max)                         :: omi_sol_wav_avg
  REAL    (KIND=r8), DIMENSION (nxtrack_max)                         :: &
       omi_solcal_chisq, omi_radcal_chisq, omi_radref_chisq, &
       omi_radref_col,   omi_radref_dcol,  omi_radref_rms
  REAL    (KIND=r8), DIMENSION (2,nxtrack_max,0:nlines_max-1)        :: omi_wavwin_rad, omi_fitwin_rad
  REAL    (KIND=r8), DIMENSION (2,nxtrack_max)                       :: omi_wavwin_sol, omi_fitwin_sol

  INTEGER (KIND=i4), DIMENSION (nxtrack_max)                         :: omi_solfit_xflag
  REAL    (KIND=dp), DIMENSION (max_calfit_idx, nxtrack_max)         :: omi_solfit_pars
  INTEGER (KIND=i4), DIMENSION (nxtrack_max, n_rad_winwav)           :: omi_rad_winwav_idx

  ! ---------------------------------
  ! OMI swath names (Uv-1 and UV-2)
  ! --------------------------------  
  CHARACTER (LEN=maxchlen), DIMENSION(mswath) :: omi_radiance_swathname, omi_irradiance_swathname
  INTEGER, DIMENSION(mswath)                  :: omichs

  ! ------------------------------
  ! Distance between Earth and Sun
  ! ------------------------------
  REAL (KIND=r4) :: EarthSunDistance

  ! ---------------------------
  ! OMI L2 output data QA flags
  ! ----------------------------------------------------------------
  ! Per Centages of output column data that are used to classify the
  ! scientific data quality:
  !         Good data >= QA_PERCENT_PASS   : "Passed"
  !         Good data >= QA_PERCENT_SUSPECT: "Suspect"
  !         Good data <  QA_PERCENT_SUSPECT: "Failed"
  ! ----------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: qa_percent_passed = 80, qa_percent_suspect = 50

  ! calibration variables for all the cross track pixels
  REAL (KIND=dp), DIMENSION (max_fit_pts, nxtrack_max) ::  omi_slitwav_sol=0.0, omi_sswav_sol=0.0, &
     omi_slitwav_rad=0.0, omi_sswav_rad=0.0 
  REAL (KIND=dp), DIMENSION (max_fit_pts, max_calfit_idx, 2, nxtrack_max)  :: omi_solslitfit=0.0,  &
       omi_solwavfit=0.0, omi_radslitfit=0.0,  omi_radwavfit=0.0
  INTEGER, DIMENSION (nxtrack_max) :: omi_nslit_sol, omi_nwavcal_sol, omi_nslit_rad, omi_nwavcal_rad
  REAL (KIND=dp), DIMENSION(maxwin,max_calfit_idx,2, nxtrack_max) :: omi_solwinfit, omi_radwinfit
  REAL (KIND=dp), DIMENSION(maxwin)                               :: omi_redslw
  REAL (KIND=dp), DIMENSION(maxwin, nxtrack_max)                  :: omi_wincal_wav
  REAL (KIND=dp), DIMENSION(max_fit_pts, 3,max_fit_pts ,nxtrack_max) ::omi_cali_sol

!  ! ------------------------------------------------------
!  ! Finally some variables that will be initialized in the
!  ! course of the processing.
!  ! ------------------------------------------------------
!  INTEGER (KIND=i4) :: &
!       n_omi_radwvl, n_omi_radrefwvl, n_omi_irradwvl,  &
!       ntimes_irrad, nxtrack_irrad,   nwavelcoef_irrad, &
!       ntimes_rad,   nxtrack_rad,     nwavelcoef_rad,   &
!       ntimes_smapix_irrad, ntimes_smpix_rad, nclenfit
!
  ! --------------------------------
  ! Current cross-track pixel number
  ! --------------------------------
  INTEGER (KIND=i4) :: curr_xtrack_pixnum

  LOGICAL, DIMENSION (nxtrack_max) ::  omi_cross_track_skippix = .FALSE.


END MODULE OMSAO_omidata_module
