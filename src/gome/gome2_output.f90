MODULE g2_h5_output

  USE h5_util !, only : HID_T
  !USE OMSAO_he5_module
  USE OMSAO_indices_module, only : refspec_strings

  IMPLICIT NONE

  INTEGER :: Meta_gid, Geo_gid, Col_gid
  INTEGER(HID_T) :: fid       ! File1 identifier 
  integer, parameter :: dims0(MAXRANK)=(/-1,-1,-1,-1/)


  TYPE (DSh5_T) :: ds_fnamesol  = DSH5_T("/Metadata/L1B_Filename_Solar"     ,-1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_fnamerad  = DSH5_T("/Metadata/L1B_Filename_EarthShine",-1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_npix      = DSH5_T("/Metadata/nGround_Pix",       -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_nscan     = DSH5_T("/Metadata/nScan",             -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_iscan     = DSH5_T("/Metadata/Scan_Index",        -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_Xtrack_dim= DSH5_T("/Metadata/Xtrack_per_Scan",   -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_Xtrack_idx= DSH5_T("/Metadata/Xtrack_Indices",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_int_time  = DSH5_T("/Metadata/IntegrationTime",   -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_loncornA=   DSH5_T("/Geolocation/Lon_CornerA",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_latcornA=   DSH5_T("/Geolocation/Lat_CornerA",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_loncornB=   DSH5_T("/Geolocation/Lon_CornerB",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_latcornB=   DSH5_T("/Geolocation/Lat_CornerB",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_loncornC=   DSH5_T("/Geolocation/Lon_CornerC",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_latcornC=   DSH5_T("/Geolocation/Lat_CornerC",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_loncornD=   DSH5_T("/Geolocation/Lon_CornerD",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_latcornD=   DSH5_T("/Geolocation/Lat_CornerD",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_lonc=       DSH5_T("/Geolocation/Lon_Centre",     -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_latc=       DSH5_T("/Geolocation/Lat_Centre",     -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_szac=       DSH5_T("/Geolocation/SZA_Centre",     -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_vzac=       DSH5_T("/Geolocation/VZA_Centre",     -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_razc=       DSH5_T("/Geolocation/RAZA_Centre",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_vazc=       DSH5_T("/Geolocation/VAZA_Centre",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_sazc=       DSH5_T("/Geolocation/SAZA_Centre",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_utc_year=   DSH5_T("/Geolocation/UTC_Year",       -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_utc_doy=    DSH5_T("/Geolocation/UTC_DOY",        -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_utc_sec=    DSH5_T("/Geolocation/UTC_Sec",        -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_scd=        DSH5_T("/Column/SCD",                 -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_vcd=        DSH5_T("/Column/VCD",                 -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_dscd=       DSH5_T("/Column/DSCD",                -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_dvcd=       DSH5_T("/Column/DVCD",                -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rms=        DSH5_T("/Column/RMS",                 -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_amf=        DSH5_T("/Column/Amf",                 -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_amfgeo=     DSH5_T("/Column/AmfGeo",              -1,-1,-1,1,1,dims0 )

  TYPE (DSh5_T) :: ds_f_sunglint= DSH5_T("/Column/SunGlintFlag",        -1,-1,-1,1,1,dims0 )
  !TYPE (DSh5_T) :: ds_f_saa=      DSH5_T("/Column/SAAFlag",        -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_fitmode=    DSH5_T("/Column/CloudFitMode",        -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_cldp=       DSH5_T("/Column/CloudTopPressure",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_cldf=       DSH5_T("/Column/CloudFraction",       -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_cldR=       DSH5_T("/Column/CloudReflectivity",   -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_cldchisq=   DSH5_T("/Column/CloudChiSquare",      -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_Errcldp=    DSH5_T("/Column/ErrCloudTopPressure", -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_Errcldf=    DSH5_T("/Column/ErrCloudFraction",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_surfR758=   DSH5_T("/Column/SurfReflect_758nm",   -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_surfR772=   DSH5_T("/Column/SurfReflect_772nm",   -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_surfP=      DSH5_T("/Column/SurfacePressure",     -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_qpolss=     DSH5_T("/Column/SingleScatPolQ",      -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_upolss=     DSH5_T("/Column/SingleScatPolU",      -1,-1,-1,1,1,dims0 )

  TYPE (DSh5_T) :: ds_sname=      DSH5_T("/Metadata/SpeciesName",       -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_ring=    DSH5_T("/Metadata/Refspec_RRaman",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_o3_t1=   DSH5_T("/Metadata/Refspec_O3_T1",     -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_o3_t2=   DSH5_T("/Metadata/Refspec_O3_T2",     -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_o3_t3=   DSH5_T("/Metadata/Refspec_O3_T3",     -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_no2_t1=  DSH5_T("/Metadata/Refspec_NO2_T1",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_no2_t2=  DSH5_T("/Metadata/Refspec_NO2_T2",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_o4=      DSH5_T("/Metadata/Refspec_O2-O2",     -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_so2=     DSH5_T("/Metadata/Refspec_SO2",       -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_so2_v=   DSH5_T("/Metadata/Refspec_SO2_Vert",  -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_bro_t1=  DSH5_T("/Metadata/Refspec_BrO_T1",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_bro_t2=  DSH5_T("/Metadata/Refspec_BrO_T2",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_oclo=    DSH5_T("/Metadata/Refspec_OClO",      -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_hcho=    DSH5_T("/Metadata/Refspec_HCHO",      -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_h2o=     DSH5_T("/Metadata/Refspec_H2O",       -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_chocho=  DSH5_T("/Metadata/Refspec_CHOCHO",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_io=      DSH5_T("/Metadata/Refspec_IO",        -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_vram=    DSH5_T("/Metadata/Refspec_VRaman",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_rs_common=  DSH5_T("/Metadata/Refspec_CommonMode",-1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_llim=       DSH5_T("/Metadata/Fit_Window_Lower",  -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_ulim=       DSH5_T("/Metadata/Fit_Window_Upper",  -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_algdesc =   DSH5_T("/Metadata/Algorithm",         -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_proccentre =DSH5_T("/Metadata/ProcessingCentre",  -1,-1,-1,1,1,dims0 )

  ! CRN: HDF output parameters for SAOPROF retrievals
  TYPE (DSh5_T) :: ds_sza_atm=   DSH5_T("/Geolocation/SZA_Effective",   -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_vza_atm=   DSH5_T("/Geolocation/VZA_Effective",   -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_aza_atm=   DSH5_T("/Geolocation/AZA_Effective",   -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_sca_atm=   DSH5_T("/Geolocation/SCA_Effective",   -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_exval=     DSH5_T("/Column/ExitStatus",           -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_numiter=   DSH5_T("/Column/nIterations",          -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_saaflag=   DSH5_T("/Column/SaoprofSaaFlag",       -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_nsaa_spike=DSH5_T("/Column/SaoprofNsaaSpike",     -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_avgres=    DSH5_T("/Column/AvgResidual",          -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3dfs=     DSH5_T("/Column/O3DFS",                -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3info=    DSH5_T("/Column/O3Info",               -1,-1,-1,1,1,dims0 )

  TYPE (DSh5_T) :: ds_whichalb=  DSH5_T("/Column/WhichAlbedo",          -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_nalb=      DSH5_T("/Column/nAlbedo",              -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_albinit=   DSH5_T("/Column/AlbedoInit",           -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_alb=       DSH5_T("/Column/Albedo",               -1,-1,-1,1,1,dims0 )

  TYPE (DSh5_T) :: ds_newcldf=   DSH5_T("/Column/saoprof_CldFraction",  -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_cldod=     DSH5_T("/Column/saoprof_CldOD",        -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_newcldp=   DSH5_T("/Column/saoprof_CldPressure",  -1,-1,-1,1,1,dims0 )

  TYPE (DSh5_T) :: ds_nlay=      DSH5_T("/Column/nLayers",              -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_ntp=       DSH5_T("/Column/nTp",                  -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_p=         DSH5_T("/Column/Pressure",             -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_z=         DSH5_T("/Column/Zkm",                  -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_T=         DSH5_T("/Column/Temperature",          -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3ap=      DSH5_T("/Column/O3_Ap",                -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3apstd=   DSH5_T("/Column/O3_ApStd",             -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3prof=    DSH5_T("/Column/O3_Prof",              -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3profstd= DSH5_T("/Column/O3_ProfStd",           -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3profnstd=DSH5_T("/Column/O3_ProfNstd",          -1,-1,-1,1,1,dims0 )

  TYPE (DSh5_T) :: ds_o3totap=   DSH5_T("/Column/O3_TotalAp",           -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3tot=     DSH5_T("/Column/O3_Total",             -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3dtot1=   DSH5_T("/Column/O3_TotalE1",           -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3dtot2=   DSH5_T("/Column/O3_TotalE2",           -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3stratap= DSH5_T("/Column/O3_StratAp",           -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3strat=   DSH5_T("/Column/O3_Strat",             -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3dstrat1= DSH5_T("/Column/O3_StratE1",           -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3dstrat2= DSH5_T("/Column/O3_StratE2",           -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3tropap=  DSH5_T("/Column/O3_TropAp",            -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3trop=    DSH5_T("/Column/O3_Trop",              -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3dtrop1=  DSH5_T("/Column/O3_TropE1",            -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_o3dtrop2=  DSH5_T("/Column/O3_TropE2",            -1,-1,-1,1,1,dims0 )

  TYPE (DSh5_T) :: ds_ai=        DSH5_T("/Column/Aerosol_Index",        -1,-1,-1,1,1,dims0 )

  TYPE (DSh5_T) :: ds_fitvarstr= DSH5_T("/Metadata/FitVariableStr",     -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_fitvar=    DSH5_T("/Column/FitVariableRet",       -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_fitvarstd= DSH5_T("/Column/FitVariableRetStd",    -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_fitvarnstd=DSH5_T("/Column/FitVariableRetNStd",   -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_fitvarap=  DSH5_T("/Column/FitVariableAp",        -1,-1,-1,1,1,dims0 )

  TYPE (DSh5_T) :: ds_wvlrad    = DSH5_T("/Column/WvlRad",              -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_obsrad    = DSH5_T("/Column/ObsRad",              -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_simrad    = DSH5_T("/Column/SimRad",              -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_resrad    = DSH5_T("/Column/ResRad",              -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_radspc    = DSH5_T("/Column/Radiance",              -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_stokes    = DSH5_T("/Column/StokesFraction",      -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_wvlsol    = DSH5_T("/Column/WvlSol",              -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_obssol    = DSH5_T("/Column/ObsSol",              -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_fitwgt    = DSH5_T("/Column/FitWeight",           -1,-1,-1,1,1,dims0 )

  TYPE (DSh5_T) :: ds_divsun    = DSH5_T("/Column/DivSun",              -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_divrad    = DSH5_T("/Column/DivRad",              -1,-1,-1,1,1,dims0 )

  TYPE (DSh5_T) :: ds_so2prof   = DSH5_T("/Column/SO2_Prof",            -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_so2avgk   = DSH5_T("/Column/SO2_AvgK",            -1,-1,-1,1,1,dims0 )

  TYPE (DSh5_T) :: ds_so2amf    = DSH5_T("/Column/SO2_amf",             -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_so2vdfs    = DSH5_T("/Column/so2v_dfs",             -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_so2zdfs    = DSH5_T("/Column/so2z_dfs",             -1,-1,-1,1,1,dims0 )
  TYPE (DSh5_T) :: ds_so2offset  = DSH5_T("/Column/so2_offset",         -1,-1,-1,1,1,dims0 )

END MODULE g2_h5_output


SUBROUTINE gome2_init_output( Nscans, N_pixels, outfile)
  USE h5_util
  USE g2_h5_output
  USE HDF5 !, ONLY: H5F_ACC_TRUNC_F, H5T_NATIVE_INTEGER, H5T_NATIVE_CHARACTER, & 
  !H5T_NATIVE_REAL, h5fcreate_f, h5gcreate_f,  h5open_f
  USE OMSAO_variables_module,   ONLY: &
       refspec_fname, sol_winwav_lim, rad_winwav_lim, pixnum_lim, radwavcal_freq, &
       sza_atm, vza_atm, aza_atm, fitvar_rad, fitvar_sol, lo_sunbnd, &
       up_sunbnd, phase, verb_thresh_lev, amf, amfgeo, sol_zen_eff, fitcol_idx, &
       numwin, winlim
  USE OMSAO_indices_module, only:  min_rs_idx, max_rs_idx, ring_idx, &
       o3_t1_idx, o3_t2_idx, o3_t3_idx, no2_t1_idx, no2_t2_idx, o2o2_idx, & 
       so2_idx, so2v_idx, bro_idx, bro2_idx, oclo_idx, hcho_idx, h2o_idx, glyox_idx, io_idx,  &
       vraman_idx, comm_idx 
  !USE OMSAO_variables_module,   ONLY: fitcol_idx, n_mol_fit

  IMPLICIT NONE

  INTEGER,            INTENT(IN) :: Nscans, N_pixels
  CHARACTER(LEN=255), INTENT(IN) :: outfile

  INTEGER                        :: npixels, nscan, nchar
  INTEGER,            PARAMETER  :: one = 1
  CHARACTER(LEN=127), PARAMETER  :: geogroup  = "Geolocation"
  CHARACTER(LEN=127), PARAMETER  :: colgroup  = "Column"
  CHARACTER(LEN=127), PARAMETER  :: metagroup = "Metadata"
  CHARACTER(LEN=127)             :: colname, proc_centre_string, alg_info_string
  INTEGER                        :: i 
  INTEGER                        :: error 

  npixels=N_pixels
  nscan=Nscans
  proc_centre_string = "Harvard/Smithsonian, cnowlan@cfa.harvard.edu"
  alg_info_string    = "SAOPROF, Harvard/Smithsonian SAO, Dec. 2008"
  !
  ! Initialize FORTRAN interface.
  !
  CALL h5open_f(error)
  CALL h5fcreate_f(TRIM(outfile), H5F_ACC_TRUNC_F, fid, error)
  call check_msg(error, "h5fcreate_f failed A subr.gome2_output.f90 "//TRIM(outfile))
  CALL h5gcreate_f(fid, geogroup,  geo_gid,   error)
  call check_msg(error, "h5gcreate_f failed A subr.gome2_output.f90 "//TRIM(geogroup))
  CALL h5gcreate_f(fid, metagroup, meta_gid, error)
  call check_msg(error, "h5gcreate_f failed A  subr.gome2_output.f90 "//TRIM(metagroup))
  CALL h5gcreate_f(fid, colgroup,  col_gid,   error)
  call check_msg(error, "h5gcreate_f failed A subr.gome2_output.f90 "//TRIM(colgroup))


  colname = TRIM(ADJUSTL(refspec_strings(fitcol_idx(1))))


  call createDS( fid, H5T_NATIVE_INTEGER, one,      ds_nscan      )
  call createDS( fid, H5T_NATIVE_INTEGER, one,      ds_npix       )
  call createDS( fid, H5T_NATIVE_INTEGER, npixels,  ds_iscan      )
  call createDS( fid, H5T_NATIVE_INTEGER, nscan,    ds_Xtrack_dim )
  call createDS( fid, H5T_NATIVE_INTEGER, npixels,  ds_Xtrack_idx )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_int_time   )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_loncornA   )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_latcornA   )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_loncornB   )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_latcornB   )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_loncornC   )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_latcornC   )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_loncornD   )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_latcornD   )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_lonc       )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_latc       )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_szac       )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_vzac       )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_razc       )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_sazc       )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_vazc       )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_scd        )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_vcd        )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_dscd       )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_dvcd       )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_rms        )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_amf        )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_amfgeo     )
  call createDS( fid, H5T_NATIVE_INTEGER, nscan,    ds_utc_doy    )
  call createDS( fid, H5T_NATIVE_INTEGER, nscan,    ds_utc_year   )
  call createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_utc_sec    )
  call createDS( fid, H5T_NATIVE_REAL,    one,      ds_llim       )
  call createDS( fid, H5T_NATIVE_REAL,    one,      ds_ulim       )

  call createDS(  fid, H5T_NATIVE_INTEGER, npixels,  ds_f_sunglint)
  !call createDS(  fid, H5T_NATIVE_INTEGER, npixels,  ds_f_saa     )
  call createDS(  fid, H5T_NATIVE_INTEGER, npixels,  ds_fitmode   )
  call createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_cldp      )
  call createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_cldf      )
  call createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_cldR      )
  call createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_cldchisq  )
  call createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_Errcldp   )
  call createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_Errcldf   )
  call createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_surfR758  )
  call createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_surfR772  )
  call createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_surfP     )
  call createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_qpolss    )
  call createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_upolss    )


  nchar= len(trim(adjustl(colname)))
  call write_scalar_ik4( ds_nscan,   Nscans)
  call write_scalar_ik4( ds_npix,    N_pixels)
  call write_scalar_r4(  ds_llim,    real(winlim(numwin,1) ,4 ))
  call write_scalar_r4(  ds_ulim,    real(winlim(numwin,2) ,4 ))


  call write_string( fid, ds_sname,      meta_gid,  trim(adjustl(colname))             ,      file_bname=.FALSE.    )
  call write_string( fid, ds_rs_ring,    meta_gid,  trim(adjustl(refspec_fname(ring_idx))),   file_bname=.TRUE.     )
  call write_string( fid, ds_rs_o3_t1,   meta_gid,  trim(adjustl(refspec_fname(o3_t1_idx))),  file_bname=.TRUE.     )
  call write_string( fid, ds_rs_o3_t2,   meta_gid,  trim(adjustl(refspec_fname(o3_t2_idx))),  file_bname=.TRUE.     )
  call write_string( fid, ds_rs_o3_t3,   meta_gid,  trim(adjustl(refspec_fname(o3_t3_idx))),  file_bname=.TRUE.     )
  call write_string( fid, ds_rs_no2_t1,  meta_gid,  trim(adjustl(refspec_fname(no2_t1_idx))), file_bname=.TRUE.     )
  call write_string( fid, ds_rs_no2_t2,  meta_gid,  trim(adjustl(refspec_fname(no2_t2_idx))), file_bname=.TRUE.     )
  call write_string( fid, ds_rs_o4,      meta_gid,  trim(adjustl(refspec_fname(o2o2_idx)))  , file_bname=.TRUE.     )
  call write_string( fid, ds_rs_so2,     meta_gid,  trim(adjustl(refspec_fname(so2_idx)))   , file_bname=.TRUE.     )
  call write_string( fid, ds_rs_bro_t1,  meta_gid,  trim(adjustl(refspec_fname(bro_idx)))   , file_bname=.TRUE.     )
  call write_string( fid, ds_rs_oclo,    meta_gid,  trim(adjustl(refspec_fname(oclo_idx)))  , file_bname=.TRUE.     )
  call write_string( fid, ds_rs_hcho,    meta_gid,  trim(adjustl(refspec_fname(hcho_idx)))  , file_bname=.TRUE.     )
  call write_string( fid, ds_rs_h2o,     meta_gid,  trim(adjustl(refspec_fname(h2o_idx)))   , file_bname=.TRUE.     )
  call write_string( fid, ds_rs_chocho,  meta_gid,  trim(adjustl(refspec_fname(glyox_idx))) , file_bname=.TRUE.     )
  call write_string( fid, ds_rs_io,      meta_gid,  trim(adjustl(refspec_fname(io_idx)))    , file_bname=.TRUE.     )
  call write_string( fid, ds_rs_vram,    meta_gid,  trim(adjustl(refspec_fname(vraman_idx))), file_bname=.TRUE.     )
  call write_string( fid, ds_rs_common,  meta_gid,  trim(adjustl(refspec_fname(comm_idx)))  , file_bname=.TRUE.     )
  call write_string( fid, ds_algdesc  ,  meta_gid,  alg_info_string                                                 )
  call write_string( fid, ds_proccentre, meta_gid,  proc_centre_string                                              )


END SUBROUTINE gome2_init_output


SUBROUTINE gome2_init_output_profile(Nscans, N_pixels, outfile, l2_hdf_flag)

  USE h5_util
  USE g2_h5_output
  USE HDF5
  USE OMSAO_variables_module,   ONLY: &
       refspec_fname, rad_winwav_lim,  n_fitvar_rad, curr_sol_spec, curr_rad_spec, &
       numwin, winlim
  USE OMSAO_gome_data_module,   ONLY: n_gome_solpts
  USE OMSAO_indices_module,     ONLY: &
       min_rs_idx, max_rs_idx, ring_idx, o3_t1_idx, o3_t2_idx, &
       o3_t3_idx, no2_t1_idx, no2_t2_idx, o2o2_idx, so2_idx,   &
       so2v_idx, bro_idx, bro2_idx, oclo_idx, hcho_idx, comm_idx 
  USE ozprof_data_module,       ONLY: nalb, nlay, div_sun, &
     atmwrt, ozwrtcorr, ozwrtcovar, ozwrtcontri, ozwrtres, ozwrtavgk, ozwrtsnr 

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER, INTENT(IN) :: Nscans, N_pixels, l2_hdf_flag
  CHARACTER(LEN=255), INTENT(IN) :: outfile

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER                        :: npixels, nscan, nchar, nlayers, nvar, nalbtmp, nwvl
  INTEGER,            PARAMETER  :: one = 1
  CHARACTER(LEN=127), PARAMETER  :: geogroup  = "Geolocation"
  CHARACTER(LEN=127), PARAMETER  :: colgroup  = "Column"
  CHARACTER(LEN=127), PARAMETER  :: metagroup = "Metadata"
  CHARACTER(LEN=127) :: colname, proc_centre_string, alg_info_string
  INTEGER :: i 
  INTEGER :: error 
  REAL    (KIND=4), DIMENSION (1) :: valr

  npixels = N_pixels
  nlayers = nlay
  nscan = Nscans
  nalbtmp = nalb
  nwvl = n_gome_solpts
  nvar = n_fitvar_rad
  proc_centre_string = "Harvard/Smitsonian, cnowlan@cfa.harvard.edu"
  alg_info_string    = "BOREAS, Harvard/Smithsonian SAO, Sep. 2008"



  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  CALL h5fcreate_f(TRIM(outfile), H5F_ACC_TRUNC_F, fid, error)
  CALL check_msg(error, "h5fcreate_f failed A subr.gome2_output.f90 "//TRIM(outfile))
  CALL h5gcreate_f(fid, geogroup,  geo_gid,   error)
  CALL check_msg(error, "h5gcreate_f failed A subr.gome2_output.f90 "//TRIM(geogroup))
  CALL h5gcreate_f(fid, metagroup, meta_gid, error)
  CALL check_msg(error, "h5gcreate_f failed A  subr.gome2_output.f90 "//TRIM(metagroup))
  CALL h5gcreate_f(fid, colgroup,  col_gid,   error)
  CALL check_msg(error, "h5gcreate_f failed A subr.gome2_output.f90 "//TRIM(colgroup))


  CALL createDS( fid, H5T_NATIVE_INTEGER, one,      ds_nscan      )
  CALL createDS( fid, H5T_NATIVE_INTEGER, one,      ds_npix       )
  CALL createDS( fid, H5T_NATIVE_INTEGER, npixels,  ds_iscan      )
  CALL createDS( fid, H5T_NATIVE_INTEGER, nscan,    ds_Xtrack_dim )
  CALL createDS( fid, H5T_NATIVE_INTEGER, npixels,  ds_Xtrack_idx )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_int_time   )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_loncornA   )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_latcornA   )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_loncornB   )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_latcornB   )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_loncornC   )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_latcornC   )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_loncornD   )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_latcornD   )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_lonc       )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_latc       )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_szac       )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_vzac       )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_razc       )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_sazc       )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_vazc       )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_rms        )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_amfgeo     )
  CALL createDS( fid, H5T_NATIVE_INTEGER, nscan,    ds_utc_doy    )
  CALL createDS( fid, H5T_NATIVE_INTEGER, nscan,    ds_utc_year   )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels,  ds_utc_sec    )
  CALL createDS( fid, H5T_NATIVE_REAL,    one,      ds_llim       )
  CALL createDS( fid, H5T_NATIVE_REAL,    one,      ds_ulim       )

  CALL createDS(  fid, H5T_NATIVE_INTEGER, npixels,  ds_f_sunglint)
  CALL createDS(  fid, H5T_NATIVE_INTEGER, npixels,  ds_fitmode   )
  CALL createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_cldp      )
  CALL createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_cldf      )
  CALL createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_cldR      )
  CALL createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_cldchisq  )
  CALL createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_Errcldp   )
  CALL createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_Errcldf   )
  CALL createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_surfR758  )
  CALL createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_surfR772  )
  CALL createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_surfP     )
  CALL createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_qpolss    )
  CALL createDS(  fid, H5T_NATIVE_REAL,    npixels,  ds_upolss    )

  CALL write_scalar_ik4( ds_nscan,   Nscans)
  CALL write_scalar_ik4( ds_npix,    N_pixels)
  CALL write_scalar_r4(  ds_llim,    real(winlim(1,1) ,4 ))
  CALL write_scalar_r4(  ds_ulim,    real(winlim(numwin,2) ,4 ))
  CALL write_string( fid, ds_sname,      meta_gid,  TRIM(ADJUSTL(colname))             ,      file_bname=.FALSE.    )
  CALL write_string( fid, ds_rs_ring,    meta_gid,  TRIM(ADJUSTL(refspec_fname(ring_idx))),   file_bname=.TRUE.     )
  CALL write_string( fid, ds_rs_o3_t1,   meta_gid,  TRIM(ADJUSTL(refspec_fname(o3_t1_idx))),  file_bname=.TRUE.     )
  CALL write_string( fid, ds_rs_o3_t2,   meta_gid,  TRIM(ADJUSTL(refspec_fname(o3_t2_idx))),  file_bname=.TRUE.     )
  CALL write_string( fid, ds_rs_o3_t3,   meta_gid,  TRIM(ADJUSTL(refspec_fname(o3_t3_idx))),  file_bname=.TRUE.     )
  CALL write_string( fid, ds_rs_no2_t1,  meta_gid,  TRIM(ADJUSTL(refspec_fname(no2_t1_idx))), file_bname=.TRUE.     )
  CALL write_string( fid, ds_rs_no2_t2,  meta_gid,  TRIM(ADJUSTL(refspec_fname(no2_t2_idx))), file_bname=.TRUE.     )
  CALL write_string( fid, ds_rs_o4,      meta_gid,  TRIM(ADJUSTL(refspec_fname(o2o2_idx)))  , file_bname=.TRUE.     )
  CALL write_string( fid, ds_rs_so2,     meta_gid,  TRIM(ADJUSTL(refspec_fname(so2_idx)))   , file_bname=.TRUE.     )
  CALL write_string( fid, ds_rs_so2_v,   meta_gid,  TRIM(ADJUSTL(refspec_fname(so2v_idx)))  , file_bname=.TRUE.     )
  CALL write_string( fid, ds_rs_bro_t1,  meta_gid,  TRIM(ADJUSTL(refspec_fname(bro_idx)))   , file_bname=.TRUE.     )
  CALL write_string( fid, ds_rs_bro_t2,  meta_gid,  TRIM(ADJUSTL(refspec_fname(bro2_idx)))  , file_bname=.TRUE.     )
  CALL write_string( fid, ds_rs_oclo,    meta_gid,  TRIM(ADJUSTL(refspec_fname(oclo_idx)))  , file_bname=.TRUE.     )
  CALL write_string( fid, ds_rs_hcho,    meta_gid,  TRIM(ADJUSTL(refspec_fname(hcho_idx)))  , file_bname=.TRUE.     )
  CALL write_string( fid, ds_rs_common,  meta_gid,  TRIM(ADJUSTL(refspec_fname(comm_idx)))  , file_bname=.TRUE.     )
  CALL write_string( fid, ds_algdesc  ,  meta_gid,  TRIM(ADJUSTL(alg_info_string))          , file_bname=.TRUE.     )                              
  CALL write_string( fid, ds_proccentre, meta_gid,  TRIM(ADJUSTL(proc_centre_string))       , file_bname=.TRUE.     )     


  ! Variables associated with SAOPROF profile retrieval
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels, ds_sza_atm)
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels, ds_vza_atm)
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels, ds_aza_atm)
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels, ds_sca_atm)
  CALL createDS( fid, H5T_NATIVE_INTEGER, npixels, ds_exval  )
  CALL createDS( fid, H5T_NATIVE_INTEGER, npixels, ds_numiter  )
  CALL createDS( fid, H5T_NATIVE_INTEGER, npixels, ds_saaflag  )
  CALL createDS( fid, H5T_NATIVE_INTEGER, npixels, ds_nsaa_spike  )
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels, ds_avgres)
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels, ds_o3dfs)
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels, ds_o3info)
  !CALL createDS( fid, H5T_NATIVE_INTEGER, npixels, ds_whichalb)
  CALL createDS( fid, H5T_NATIVE_INTEGER, npixels, ds_nalb)
  CALL createDS_ndims( fid, H5T_NATIVE_REAL, 2, (/nalbtmp,npixels/), ds_albinit)
  CALL createDS_ndims( fid, H5T_NATIVE_REAL, 2, (/nalbtmp,npixels/), ds_alb)
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels, ds_newcldf)
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels, ds_cldod)
  CALL createDS( fid, H5T_NATIVE_REAL,    npixels, ds_newcldp)
  CALL createDS( fid, H5T_NATIVE_INTEGER, npixels, ds_nlay)
  CALL createDS( fid, H5T_NATIVE_INTEGER, npixels, ds_ntp)
  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_o3totap)
  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_o3tot)
  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_o3dtot1)
  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_o3dtot2)
  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_o3stratap)
  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_o3strat)
  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_o3dstrat1)
  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_o3dstrat2)
  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_o3tropap)
  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_o3trop)
  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_o3dtrop1)
  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_o3dtrop2)
  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_ai)
  CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nvar,npixels/), ds_fitvar)
  CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nvar,npixels/), ds_fitvarstd)
  CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nvar,npixels/), ds_fitvarnstd)
  CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nvar,npixels/), ds_fitvarap)
  !CALL createDS_ndims( fid, H5T_NATIVE_CHAR,    2, [n_fitvar_rad,6], ds_fitvarstr )

  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_so2amf)
  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_so2vdfs)
  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_so2zdfs)
  CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_so2offset)

  IF (l2_hdf_flag == 2) THEN
    IF (atmwrt) THEN
     CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nlayers,npixels/), ds_p)
     CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nlayers,npixels/), ds_z)
     CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nlayers,npixels/), ds_T)
    ENDIF
     !CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nlayers,npixels/), ds_o3ap)
     !CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nlayers,npixels/), ds_o3apstd)
     !CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nlayers,npixels/), ds_o3prof)
     !CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nlayers,npixels/), ds_o3profstd)
     !CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nlayers,npixels/), ds_o3profnstd)
    IF (ozwrtres) THEN
     CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nwvl,npixels/), ds_wvlrad)
     CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nwvl,npixels/), ds_obsrad)
     CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nwvl,npixels/), ds_simrad)
     CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nwvl,npixels/), ds_radspc)
     !CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nwvl,npixels/), ds_stokes)
     CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nwvl,one/),     ds_wvlsol)
     CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nwvl,one/),     ds_obssol)
     CALL writeN_r4( ds_wvlsol, one, nwvl, REAL(curr_sol_spec(1,1:n_gome_solpts),4) )
     CALL writeN_r4( ds_obssol, one, nwvl, REAL(curr_sol_spec(2,1:n_gome_solpts),4) )

     CALL createDS( fid, H5T_NATIVE_REAL, one, ds_divsun)
     CALL createDS( fid, H5T_NATIVE_REAL, npixels, ds_divrad)
     valr =  REAL(div_sun,    4) ; CALL writeone_r4(  ds_divsun,    one,  valr )
    ENDIF 
    IF (ozwrtavgk) THEN  
     CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nlayers,npixels/), ds_so2prof)
     CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nlayers,npixels/), ds_so2avgk)
    ENDIF
    IF (ozwrtsnr) THEN
      CALL createDS_ndims( fid, H5T_NATIVE_REAL,    2, (/nwvl,npixels/), ds_fitwgt)
    ENDIF
 ENDIF


END SUBROUTINE gome2_init_output_profile



SUBROUTINE gome2_geoloc_output
  USE BOREAS_gome2_data_module
  USE OMSAO_gome_data_module
  USE h5_util
  USE g2_h5_output

  IMPLICIT NONE

  INTEGER           :: error
  INTEGER           :: ipix, iscan ! Elements coordinates to write to
  INTEGER           :: year_rad, month_rad,jday_rad, day_rad ! required in call to time conversion routine
  INTEGER           :: hour_rad, min_rad,sec_rad ! required in call to time conversion routine
  REAL(KIND=sp)     :: seconds
  INTEGER(KIND=4)   :: mon, day, jday, hour, min, sec
  CHARACTER(LEN=25) :: datestring

  INTEGER (KIND=4), DIMENSION (1) :: vali
  REAL    (KIND=4), DIMENSION (1) :: valr

  ipix = gome_curpix
  iscan = gome_curscan

  ! epoch_time_to_utc is a function defined in g2_l1b_reader library.
  CALL epoch_time_to_utc(gome2_utc_day, gome2_utc_millisec, gome2_utc_year,mon,& 
       day, jday, hour, min, sec, datestring)
  seconds = sec + 0.001 *  MOD(REAL(gome2_utc_millisec,4),1000.0)  + 60.0 * min + 3600.0*hour

  valr = REAL(seconds,                      4)  ; CALL writeone_r4(  ds_utc_sec,   ipix, valr )
  valr = REAL(gome_geoloc(lon_idx,1),       4)  ; CALL writeone_r4(  ds_loncornA,  ipix, valr )
  valr = REAL(gome_geoloc(lat_idx,1),       4)  ; CALL writeone_r4(  ds_latcornA,  ipix, valr )
  valr = REAL(gome_geoloc(lon_idx,2),       4)  ; CALL writeone_r4(  ds_loncornB,  ipix, valr )
  valr = REAL(gome_geoloc(lat_idx,2),       4)  ; CALL writeone_r4(  ds_latcornB,  ipix, valr )
  valr = REAL(gome_geoloc(lon_idx,3),       4)  ; CALL writeone_r4(  ds_loncornC,  ipix, valr )
  valr = REAL(gome_geoloc(lat_idx,3),       4)  ; CALL writeone_r4(  ds_latcornC,  ipix, valr )
  valr = REAL(gome_geoloc(lon_idx,4),       4)  ; CALL writeone_r4(  ds_loncornD,  ipix, valr )
  valr = REAL(gome_geoloc(lat_idx,4),       4)  ; CALL writeone_r4(  ds_latcornD,  ipix, valr )
  valr = REAL(gome_geoloc(lon_idx,5),       4)  ; CALL writeone_r4(  ds_lonc,      ipix, valr )
  valr = REAL(gome_geoloc(lat_idx,5),       4)  ; CALL writeone_r4(  ds_latc,      ipix, valr )
  valr = REAL(gome_angles_wrtn(zen0_idx,2), 4)  ; CALL writeone_r4(  ds_szac,      ipix, valr )
  valr = REAL(gome_angles_wrts(zen_idx,2),  4)  ; CALL writeone_r4(  ds_vzac,      ipix, valr )
  valr = REAL(gome_angles_wrtn(azm_idx,2),  4)  ; CALL writeone_r4(  ds_vazc,      ipix, valr )
  valr = REAL(gome_angles_wrtn(azm0_idx,2), 4)  ; CALL writeone_r4(  ds_sazc,      ipix, valr )
  valr = REAL(gome2_relazm,                 4)  ; CALL writeone_r4(  ds_razc,      ipix, valr )
  valr = REAL(integration_time,             4)  ; CALL writeone_r4(  ds_int_time,  ipix, valr )
  vali = INT(gome2_curr_ixtrack,            4)  ; CALL writeone_ik4( ds_Xtrack_idx,ipix, vali )
  vali = INT(iscan,                         4)  ; CALL writeone_ik4( ds_iscan,     ipix, vali )
  
  IF (gome2_curr_ixtrack .EQ. 1 ) THEN 
     vali = INT(gome2_utc_year,                4) ; CALL writeone_ik4( ds_utc_year,   iscan, vali )
     vali = INT(jday,                          4) ; CALL writeone_ik4( ds_utc_doy,    iscan, vali )
  ENDIF

  RETURN

END SUBROUTINE gome2_geoloc_output


SUBROUTINE gome2_geoloc_output_profile(curr_iscan)

  ! ******************
  ! Same as gome2_geoloc_output, but contains the_sza_atm, the_vza_atm, the_aza_atm etc 
  ! from Xiong Liu's SAOPROF code.
  ! ******************

  USE BOREAS_gome2_data_module
  USE OMSAO_gome_data_module
  USE OMSAO_variables_module, ONLY: &
       the_sza_atm, the_vza_atm, the_aza_atm, the_sca_atm
  USE h5_util
  USE g2_h5_output

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: curr_iscan

  INTEGER           :: error
  INTEGER           :: ipix,iscan ! Elements coordinates to write to
  INTEGER           :: year_rad, month_rad,jday_rad, day_rad ! required in call to time conversion routine
  INTEGER           :: hour_rad, min_rad,sec_rad ! required in call to time conversion routine
  REAL(KIND=sp)     :: seconds
  INTEGER(KIND=4)   :: mon, day, jday, hour, min, sec
  CHARACTER(LEN=25) :: datestring

  INTEGER (KIND=4), DIMENSION (1) :: vali
  REAL    (KIND=4), DIMENSION (1) :: valr

  ipix = gome_curpix
  iscan = curr_iscan


  ! epoch_time_to_utc is a function defined in g2_l1b_reader library.
  CALL epoch_time_to_utc(gome2_utc_day, gome2_utc_millisec, gome2_utc_year,mon,& 
       day, jday, hour, min, sec, datestring)
  seconds = sec + 0.001 *  MOD(REAL(gome2_utc_millisec,4),1000.0)  + 60.0 * min + 3600.0*hour

  valr = REAL(seconds,                      4)  ; CALL writeone_r4(  ds_utc_sec,   ipix, valr )
  valr = REAL(gome_geoloc(lon_idx,1),       4)  ; CALL writeone_r4(  ds_loncornA,  ipix, valr )
  valr = REAL(gome_geoloc(lat_idx,1),       4)  ; CALL writeone_r4(  ds_latcornA,  ipix, valr )
  valr = REAL(gome_geoloc(lon_idx,2),       4)  ; CALL writeone_r4(  ds_loncornB,  ipix, valr )
  valr = REAL(gome_geoloc(lat_idx,2),       4)  ; CALL writeone_r4(  ds_latcornB,  ipix, valr )
  valr = REAL(gome_geoloc(lon_idx,3),       4)  ; CALL writeone_r4(  ds_loncornC,  ipix, valr )
  valr = REAL(gome_geoloc(lat_idx,3),       4)  ; CALL writeone_r4(  ds_latcornC,  ipix, valr )
  valr = REAL(gome_geoloc(lon_idx,4),       4)  ; CALL writeone_r4(  ds_loncornD,  ipix, valr )
  valr = REAL(gome_geoloc(lat_idx,4),       4)  ; CALL writeone_r4(  ds_latcornD,  ipix, valr )
  valr = REAL(gome_geoloc(lon_idx,5),       4)  ; CALL writeone_r4(  ds_lonc,      ipix, valr )
  valr = REAL(gome_geoloc(lat_idx,5),       4)  ; CALL writeone_r4(  ds_latc,      ipix, valr )
  valr = REAL(gome_angles_wrtn(zen0_idx,2), 4)  ; CALL writeone_r4(  ds_szac,      ipix, valr )
  valr = REAL(gome_angles_wrts(zen_idx,2),  4)  ; CALL writeone_r4(  ds_vzac,      ipix, valr )
  valr = REAL(gome_angles_wrtn(azm_idx,2),  4)  ; CALL writeone_r4(  ds_vazc,      ipix, valr )
  valr = REAL(gome_angles_wrtn(azm0_idx,2), 4)  ; CALL writeone_r4(  ds_sazc,      ipix, valr )
  valr = REAL(integration_time,             4)  ; CALL writeone_r4(  ds_int_time,  ipix, valr )
  vali = INT(gome2_curr_ixtrack,            4)  ; CALL writeone_ik4( ds_Xtrack_idx,ipix, vali )
  vali = INT(iscan,                         4)  ; CALL writeone_ik4( ds_iscan,     ipix, vali )


  IF (gome2_curr_ixtrack .EQ. 1 ) THEN 
     vali = INT(gome2_utc_year,                4) ; CALL writeone_ik4( ds_utc_year,   iscan, vali )
     vali = INT(jday,                          4) ; CALL writeone_ik4( ds_utc_doy,    iscan, vali )
  ENDIF

  ! SAOPROF extra variables
  valr = REAL(the_sza_atm,                  4)  ; CALL writeone_r4(  ds_sza_atm,      ipix, valr )
  valr = REAL(the_vza_atm,                  4)  ; CALL writeone_r4(  ds_vza_atm,      ipix, valr )
  valr = REAL(the_aza_atm,                  4)  ; CALL writeone_r4(  ds_aza_atm,      ipix, valr )
  valr = REAL(the_sca_atm,                  4)  ; CALL writeone_r4(  ds_sca_atm,      ipix, valr )

  ! Variables from L1B file
  vali = INT(pcd_basic%F_SUNGLINT,                        4) ; CALL writeone_ik4( ds_f_sunglint,ipix, vali  )
  vali = INT(cloud%FIT_MODE(gome2_curr_ixtrack),          4) ; CALL writeone_ik4( ds_fitmode,   ipix, vali  )
  valr = REAL(cloud%FIT_1(gome2_curr_ixtrack),            4) ; CALL writeone_r4(  ds_cldp,      ipix, valr  )
  valr = REAL(cloud%FIT_2(gome2_curr_ixtrack),            4) ; CALL writeone_r4(  ds_cldf,      ipix, valr  )
  valr = REAL(cloud%CLOUD_ALBEDO(gome2_curr_ixtrack),     4) ; CALL writeone_r4(  ds_cldR,      ipix, valr  )
  valr = REAL(cloud%GOOD_FIT(gome2_curr_ixtrack),         4) ; CALL writeone_r4(  ds_cldchisq,  ipix, valr  )
  valr = REAL(cloud%E_FIT_1(gome2_curr_ixtrack),          4) ; CALL writeone_r4(  ds_Errcldp,   ipix, valr  )
  valr = REAL(cloud%E_FIT_2(gome2_curr_ixtrack),          4) ; CALL writeone_r4(  ds_Errcldf,   ipix, valr  )
  valr = REAL(cloud%SURFACE_ALBEDO(gome2_curr_ixtrack,1), 4) ; CALL writeone_r4(  ds_surfR758,  ipix, valr  )
  valr = REAL(cloud%SURFACE_ALBEDO(gome2_curr_ixtrack,2), 4) ; CALL writeone_r4(  ds_surfR772,  ipix, valr  )
  valr = REAL(cloud%SURFACE_PRESSURE(gome2_curr_ixtrack), 4) ; CALL writeone_r4(  ds_surfP,     ipix, valr  )
  valr = REAL(pol_ss%Q_POL_SS,                            4) ; CALL writeone_r4(  ds_qpolss,    ipix, valr  )
  valr = REAL(pol_ss%U_POL_SS,                            4) ; CALL writeone_r4(  ds_upolss,    ipix, valr  )

  valr =   REAL(gome_radspec(4,1), 4) ; CALL writeone_r4(  ds_amfgeo,      ipix,valr )

  RETURN

END SUBROUTINE gome2_geoloc_output_profile




SUBROUTINE gome2_column_output( fitcol, dfitcol, rms, amfgeo, amf )

  USE BOREAS_gome2_data_module
  USE OMSAO_gome_data_module, ONLY: gome_curpix
  USE h5_util, only : writeone_r4
  USE g2_h5_output

  IMPLICIT NONE

  INTEGER :: error
  INTEGER :: ipix

  INTEGER (KIND=4), DIMENSION (1) :: vali
  REAL    (KIND=4), DIMENSION (1) :: valr

  ! ---------------
  ! Input variables
  ! ---------------
  REAL (KIND=dp) :: fitcol, dfitcol, rms, amfgeo, amf

  ipix = gome_curpix

  vali = int((/pcd_basic%F_SUNGLINT/),                   KIND=4) ; call writeone_ik4( ds_f_sunglint,ipix,  vali  )
  vali = int((/cloud%FIT_MODE(gome2_curr_ixtrack)/),          4) ; call writeone_ik4( ds_fitmode,   ipix,  vali  )
  valr = real((/cloud%FIT_1(gome2_curr_ixtrack)/),            4) ; call writeone_r4(  ds_cldp,      ipix,  valr  )
  valr = real((/cloud%FIT_2(gome2_curr_ixtrack)/),            4) ; call writeone_r4(  ds_cldf,      ipix,  valr  )
  valr = real((/cloud%CLOUD_ALBEDO(gome2_curr_ixtrack)/),     4) ; call writeone_r4(  ds_cldR,      ipix,  valr  )
  valr = real((/cloud%GOOD_FIT(gome2_curr_ixtrack)/),         4) ; call writeone_r4(  ds_cldchisq,  ipix,  valr  )
  valr = real((/cloud%E_FIT_1(gome2_curr_ixtrack)/),          4) ; call writeone_r4(  ds_Errcldp,   ipix,  valr  )
  valr = real((/cloud%E_FIT_2(gome2_curr_ixtrack)/),          4) ; call writeone_r4(  ds_Errcldf,   ipix,  valr  )
  valr = real((/cloud%SURFACE_ALBEDO(gome2_curr_ixtrack,1)/), 4) ; call writeone_r4(  ds_surfR758,  ipix,  valr  )
  valr = real((/cloud%SURFACE_ALBEDO(gome2_curr_ixtrack,2)/), 4) ; call writeone_r4(  ds_surfR772,  ipix,  valr  )
  valr = real((/cloud%SURFACE_PRESSURE(gome2_curr_ixtrack)/), 4) ; call writeone_r4(  ds_surfP,     ipix,  valr  )
  valr = real((/pol_ss%Q_POL_SS/),                            4) ; call writeone_r4(  ds_qpolss,    ipix,  valr  )
  valr = real((/pol_ss%U_POL_SS/),                            4) ; call writeone_r4(  ds_upolss,    ipix,  valr  )
  valr = real((/fitcol/),                                     4) ; call writeone_r4(  ds_scd,       ipix,  valr  )
  valr = real((/dfitcol/),                                    4) ; call writeone_r4(  ds_dscd,      ipix,  valr  )
  valr = real((/fitcol/amf/),                                 4) ; call writeone_r4(  ds_vcd,       ipix,  valr  )
  valr = real((/dfitcol/amf/),                                4) ; call writeone_r4(  ds_dvcd,      ipix,  valr  )
  valr = real((/rms/),                                        4) ; call writeone_r4(  ds_rms,       ipix,  valr  )
  valr = real((/amfgeo/),                                     4) ; call writeone_r4(  ds_amfgeo,    ipix,  valr  )
  valr = real((/amf/),                                        4) ; call writeone_r4(  ds_amf,       ipix,  valr  )

  RETURN


END SUBROUTINE gome2_column_output



!******************************************************************************
SUBROUTINE gome2_output_profile(ipix, fitcol, dfitcol, amfgeo, exval, rms, l2_hdf_flag)
!******************************************************************************

  USE BOREAS_gome2_data_module
  USE h5_util, ONLY: writeone_r4, writeone_ik4
  USE g2_h5_output
  USE OMSAO_indices_module, ONLY: so2_idx, so2v_idx
  USE OMSAO_variables_module,  ONLY : fitwavs, fitres_rad, fitspec_rad, &
       n_rad_wvl, n_fitvar_rad, fitvar_rad, mask_fitvar_rad, &
       fitvar_rad_std, fitvar_rad_nstd, fitvar_rad_apriori, curr_rad_spec, currspec, &
       fitweights
  USE ozprof_data_module,      ONLY :  atmosprof, eff_alb, &
       eff_alb_init, nalb, nlay, nsaa_spike, ntp, num_iter, &
       ozdfs, ozinfo, ozprof, ozprof_ap, ozprof_apstd, ozprof_nstd, &
       ozprof_std, saa_flag, the_cfrac, the_cod, the_ctp, which_alb, mgasprof, &
       div_rad, tracegas, avg_kernel, the_ai, fgasidxs, gasidxs, ngas, &
       trace_prof, trace_avgk, trace_contri, &
       atmwrt, ozwrtcorr, ozwrtcovar, ozwrtcontri, ozwrtres, ozwrtavgk, ozwrtsnr, so2zfind, &
       vcd_so2offsetcorr
  USE OMSAO_gome_data_module, ONLY: gome_radspec

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,                         INTENT(IN) :: ipix, exval, l2_hdf_flag
  REAL (KIND=dp), DIMENSION(3),    INTENT(IN) :: fitcol
  REAL (KIND=dp), DIMENSION(3,2),  INTENT(IN) :: dfitcol
  REAL (KIND=dp),                  INTENT(IN) :: rms, amfgeo


  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER                              :: igas, error
  REAL (KIND=dp)                       :: avgres
  REAL (KIND=8), DIMENSION(n_rad_wvl)  :: simrad 

  INTEGER (KIND=4), DIMENSION (1)            :: vali
  REAL    (KIND=4), DIMENSION (1)            :: valr
  REAL    (KIND=4), DIMENSION (nalb)         :: valr_alb
  REAL    (KIND=4), DIMENSION (n_fitvar_rad) :: valr_nfitvar
  REAL    (KIND=4), DIMENSION (nlay)         :: valr_nlay
  REAL    (KIND=4), DIMENSION (nlay+2)       :: valr_nlay2
  REAL    (KIND=4), DIMENSION (n_rad_wvl)    :: valr_nwvl

  avgres = SQRT(SUM((ABS(fitres_rad(1:n_rad_wvl)) / &
       fitspec_rad(1:n_rad_wvl))**2.0)/n_rad_wvl)*100.0

  valr = REAL(rms,                                        4)  ; CALL writeone_r4(  ds_rms,       ipix, valr  )
  valr = REAL(amfgeo,                                     4)  ; CALL writeone_r4(  ds_amfgeo,    ipix, valr  )

  ! SAOPROF output
  vali = INT (exval,                      4)  ; CALL writeone_ik4( ds_exval,     ipix,  vali )
  vali = INT (num_iter,                   4)  ; CALL writeone_ik4( ds_numiter,   ipix,  vali )


  !CALL writeone_ik4( ds_saaflag,   ipix,  INT ((/saa_flag/),   4)  )
  vali = INT (nsaa_spike,                 4) ;  CALL writeone_ik4( ds_nsaa_spike,ipix, vali)
  valr = REAL(avgres,                     4) ;  CALL writeone_r4(  ds_avgres,    ipix, valr ) 
  valr = REAL(ozdfs,                      4) ;  CALL writeone_r4(  ds_o3dfs,     ipix, valr ) 
  valr = REAL(ozinfo,                     4) ;  CALL writeone_r4(  ds_o3info,    ipix, valr )
  !valr = REAL(which_alb,                  4) ;  CALL writeone_r4(  ds_whichalb,  ipix, valr ) 
  valr = REAL(nalb,                       4) ;  CALL writeone_r4(  ds_nalb,      ipix, valr ) 
  valr_alb = REAL(eff_alb_init(1:nalb),   4) ;  CALL writeN_r4(    ds_albinit,   ipix,  nalb, valr_alb ) 
  valr_alb = REAL(eff_alb(1:nalb),        4) ;  CALL writeN_r4(    ds_alb,       ipix,  nalb, valr_alb )
  valr = REAL(the_cfrac,                  4) ;  CALL writeone_r4(  ds_newcldf,   ipix, valr ) 
  valr = REAL(the_cod,                    4) ;  CALL writeone_r4(  ds_cldod,     ipix, valr ) 
  valr = REAL(the_ctp,                    4) ;  CALL writeone_r4(  ds_newcldp,   ipix, valr ) 
  valr = REAL(nlay,                       4) ;  CALL writeone_r4(  ds_nlay,      ipix, valr ) 
  valr = REAL(ntp,                        4) ;  CALL writeone_r4(  ds_ntp,       ipix, valr ) 

  valr = REAL(SUM(ozprof_ap(1:nlay)),     4) ;  CALL writeone_r4(  ds_o3totap,   ipix, valr )
  valr = REAL(fitcol(1),                  4) ;  CALL writeone_r4(  ds_o3tot,     ipix, valr )
  valr = REAL(dfitcol(1,1),               4) ;  CALL writeone_r4(  ds_o3dtot1,   ipix, valr )
  valr = REAL(dfitcol(1,2),               4) ;  CALL writeone_r4(  ds_o3dtot2,   ipix, valr )
  valr = REAL(SUM(ozprof_ap(1:ntp)),      4) ;  CALL writeone_r4(  ds_o3stratap, ipix, valr )
  valr = REAL(fitcol(2),                  4) ;  CALL writeone_r4(  ds_o3strat,   ipix, valr )
  valr = REAL(dfitcol(2,1),               4) ;  CALL writeone_r4(  ds_o3dstrat1, ipix, valr )
  valr = REAL(dfitcol(2,2),               4) ;  CALL writeone_r4(  ds_o3dstrat2, ipix, valr )
  valr = REAL(SUM(ozprof_ap(ntp+1:nlay)), 4) ;  CALL writeone_r4(  ds_o3tropap,  ipix, valr )
  valr = REAL(fitcol(3),                  4) ;  CALL writeone_r4(  ds_o3trop,    ipix, valr )
  valr = REAL(dfitcol(3,1),               4) ;  CALL writeone_r4(  ds_o3dtrop1,  ipix, valr )
  valr = REAL(dfitcol(3,2),               4) ;  CALL writeone_r4(  ds_o3dtrop2,  ipix, valr )

  valr = REAL(the_ai,                     4) ;  CALL writeone_r4(  ds_ai,        ipix, valr )

  valr_nfitvar = REAL(fitvar_rad(mask_fitvar_rad(1:n_fitvar_rad)),        4) ; CALL writeN_r4(ds_fitvar,    ipix, n_fitvar_rad, valr_nfitvar)
  valr_nfitvar = REAL(fitvar_rad_std(mask_fitvar_rad(1:n_fitvar_rad)),    4) ; CALL writeN_r4(ds_fitvarstd, ipix, n_fitvar_rad, valr_nfitvar)
  valr_nfitvar = REAL(fitvar_rad_nstd(mask_fitvar_rad(1:n_fitvar_rad)),   4) ; CALL writeN_r4(ds_fitvarnstd,ipix, n_fitvar_rad, valr_nfitvar)
  valr_nfitvar = REAL(fitvar_rad_apriori(mask_fitvar_rad(1:n_fitvar_rad)),4) ; CALL writeN_r4(ds_fitvarap,  ipix, n_fitvar_rad, valr_nfitvar)
 
  DO igas = 1, ngas
     IF (fgasidxs(igas) > 0 .AND. (gasidxs(igas) == so2_idx .OR. gasidxs(igas) == so2v_idx)) THEN
        valr = REAL(tracegas(igas,7),     4) ; CALL writeone_r4(  ds_so2amf,     ipix, valr )
        valr = REAL(avg_kernel(4,4),      4) ; CALL writeone_r4(  ds_so2vdfs,    ipix, valr )
     ENDIF
  ENDDO

  IF (so2zfind > 0) THEN
	valr = REAL(avg_kernel(so2zfind,so2zfind), 4) ; CALL writeone_r4( ds_so2zdfs, ipix, valr)
  ENDIF

  valr = REAL(vcd_so2offsetcorr, 4) ;  CALL writeone_r4( ds_so2offset, ipix, valr )


  IF (l2_hdf_flag == 2) THEN
     IF (atmwrt) THEN
     valr_nlay = REAL(atmosprof(1,1:nlay),   4) ; CALL writeN_r4(    ds_p,         ipix, nlay, valr_nlay )
     valr_nlay = REAL(atmosprof(2,1:nlay),   4) ; CALL writeN_r4(    ds_z,         ipix, nlay, valr_nlay )
     valr_nlay = REAL(atmosprof(3,1:nlay),   4) ; CALL writeN_r4(    ds_t,         ipix, nlay, valr_nlay )
     ENDIF
     !valr_nlay = REAL(ozprof_ap(1:nlay),     4) ; CALL writeN_r4(    ds_o3ap,      ipix, nlay, valr_nlay )
     !valr_nlay = REAL(ozprof_apstd(1:nlay),  4) ; CALL writeN_r4(    ds_o3apstd,   ipix, nlay, valr_nlay )
     !valr_nlay = REAL(ozprof(1:nlay),        4) ; CALL writeN_r4(    ds_o3prof,    ipix, nlay, valr_nlay )
     !valr_nlay = REAL(ozprof_std(1:nlay),    4) ; CALL writeN_r4(    ds_o3profstd, ipix, nlay, valr_nlay )
     !valr_nlay = REAL(ozprof_nstd(1:nlay),   4) ; CALL writeN_r4(    ds_o3profnstd,ipix, nlay, valr_nlay )     
     IF (ozwrtres) THEN
     simrad = fitspec_rad(1:n_rad_wvl) - fitres_rad(1:n_rad_wvl)
     valr_nwvl = REAL(fitwavs(1:n_rad_wvl),        4) ; CALL writeN_r4(    ds_wvlrad,    ipix, n_rad_wvl, valr_nwvl )
     valr_nwvl = REAL(simrad(1:n_rad_wvl),         4) ; CALL writeN_r4(    ds_simrad,    ipix, n_rad_wvl, valr_nwvl )
     valr_nwvl = REAL(fitspec_rad(1:n_rad_wvl),    4) ; CALL writeN_r4(    ds_obsrad,    ipix, n_rad_wvl, valr_nwvl )
     valr_nwvl = REAL(curr_rad_spec(2,1:n_rad_wvl),4) ; CALL writeN_r4(    ds_radspc,    ipix, n_rad_wvl, valr_nwvl )
     !!CALL writeN_r4(    ds_stokes,    ipix, n_rad_wvl, REAL(gome_radspec(4,1:n_rad_wvl), 4) )
     valr = REAL(div_rad,    4) ; CALL writeone_r4(  ds_divrad,    ipix, valr ) 
     ENDIF
     IF (ozwrtavgk) THEN
     DO igas = 1, ngas
        IF (fgasidxs(igas) > 0 .AND. (gasidxs(igas) == so2_idx .OR. gasidxs(igas) == so2v_idx)) THEN	
            valr_nlay = REAL(trace_prof(igas,1:nlay),   4) ; CALL writeN_r4( ds_so2prof, ipix, nlay, valr_nlay )
            valr_nlay = REAL(trace_avgk(igas,1:nlay),   4) ; CALL writeN_r4( ds_so2avgk, ipix, nlay, valr_nlay )
	ENDIF
      ENDDO  
     ENDIF
     IF (ozwrtsnr) THEN
        valr_nwvl = REAL(fitweights(1:n_rad_wvl), 4) ; CALL writeN_r4( ds_fitwgt, ipix, n_rad_wvl, valr_nwvl)
     ENDIF
  ENDIF

  RETURN

END SUBROUTINE gome2_output_profile



SUBROUTINE gome2_postprocessing_output( fitcol, dfitcol, amfgeo, amf, curpix )

  ! Caroline Nowlan, 10-Oct-2008

  USE OMSAO_precision_module
  USE h5_util, only : writeone_r4
  USE g2_h5_output

  IMPLICIT NONE


  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER(KIND=i4), INTENT(IN) :: curpix
  REAL(KIND=r4),    INTENT(IN) :: fitcol, dfitcol, amf, amfgeo

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER :: ipix


  ipix = curpix

  CALL writeone_r4(  ds_vcd,      ipix,  real((/fitcol/amf/),   4)  )
  CALL writeone_r4(  ds_dvcd,     ipix,  real((/dfitcol/amf/),  4)  )
  CALL writeone_r4(  ds_amfgeo,   ipix,  real((/amfgeo/),       4)  )
  CALL writeone_r4(  ds_amf,      ipix,  real((/amf/),          4)  )

END SUBROUTINE gome2_postprocessing_output







SUBROUTINE gome2_close_output

  USE h5_util
  USE g2_h5_output
  USE HDF5, ONLY:  h5gclose_f, h5fclose_f, h5close_f

  IMPLICIT NONE

  INTEGER :: error

  CALL h5gclose_f(geo_gid,  error)
  CALL check_msg(error, "h5gclose geo failed A subr.gome2_close_output" )
  CALL h5gclose_f(meta_gid, error)
  CALL check_msg(error, "h5gclose meta failed A subr.gome2_close_output" )
  CALL h5gclose_f(col_gid,  error)
  CALL check_msg(error, "h5gclose col failed A subr.gome2_close_output" )
  CALL h5fclose_f(fid, error)
  CALL h5close_f(error)

  RETURN

END SUBROUTINE gome2_close_output



SUBROUTINE gome2_count_pixels( total_pixels , nscans , fptr )

  USE control_input, ONLY: bands
  USE BOREAS_gome2_data_module, ONLY :  &
       nxtrack_per_scan, scan_start_index, n_xtrack_per_record
  USE file_handling, ONLY: open_file, fileptr
  USE gome_mdr_1b_earthshine_module

  IMPLICIT NONE

  TYPE(FILEPTR),INTENT(INOUT) :: fptr
  INTEGER, INTENT(INOUT) :: total_pixels
  INTEGER, INTENT(IN) :: nscans

  INTEGER :: irecord, record_size
  INTEGER :: curr_nxtrack
  INTEGER :: n_gome2_xtrack
  
  total_pixels = 0
  DO irecord = 1, nscans
     curr_nxtrack = n_gome2_xtrack(irecord, fptr)
     total_pixels = total_pixels + curr_nxtrack
     nxtrack_per_scan(irecord) = curr_nxtrack
  ENDDO

!!$  total_pixels = 0
!!$  DO irecord = 1, nscans
!!$     gome2_curr_nXtrack = n_xtrack_per_record(irecord,iband,fptr)
!!$     total_pixels = total_pixels + gome2_curr_nXtrack
!!$     nxtrack_per_scan(irecord) = gome2_curr_nXtrack
!!$  ENDDO

  scan_start_index(1) = 1
  DO irecord = 2, nscans
     scan_start_index(irecord) = scan_start_index(irecord-1) + nxtrack_per_scan(irecord-1)
  ENDDO

END SUBROUTINE gome2_count_pixels


