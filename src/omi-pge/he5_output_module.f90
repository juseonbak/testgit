! Subroutines to be defined
! 1. InitSwath: create file, setup swath (dimensions)
! 2. Write all the geolocation fields (once)
! 3. Write data fields (one pixel)
! 4. Write missing/fill data fields (one pixel)
! 5. Write Attributes
! 6. Close File/Swath

MODULE he5_output_module
  USE OMSAO_parameters_module, ONLY: maxwin, max_fit_pts, maxlay
  USE OMSAO_HE5_module
  USE he5_l2writer_class
  USE OMSAO_omidata_module,    ONLY: mTimes=>ntimes, mXtrack=>nfxtrack,                      &
       nxbin, nybin, nlines_max, offset_line, omi_nwav_irrad, mswath, omi_redslw, orbnum
  USE OMSAO_variables_module,  ONLY: l2_filename,  numwin,        &
       radnhtrunc, l2_swathname, fitvar_rad, mask_fitvar_rad, n_fitvar_rad,                  &
       fitvar_rad_std, n_rad_wvl, fitvar_rad_apriori, fitspec_rad, fitres_rad, fitwavs,      &
       fitvar_rad_str, fitvar_rad_nstd, simspec_rad, clmspec_rad, actspec_rad, fitweights,   &
       database, refidx, nradpix, fothvarpos, fitvar_rad_aperror, reduce_resolution,         &
       fitvar_rad_unit, which_slit, rm_mgline, winlim, n_band_avg, n_band_samp, wavcal,      &
       use_backup, yn_varyslit, refspec_fname, do_bandavg, use_redfixwav, redfixwav_fname,   &
       redsampr, redlam, wcal_bef_coadd, szamax, zatmos, radwavcal_freq, max_itnum_rad,      &
       n_wavcal_step, n_slit_step, slit_fit_pts, wavcal_fit_pts, smooth_slit, wavcal_sol,    &
       slit_trunc_limit, l1b_irrad_filename, l1b_rad_filename, l2_cld_filename,              &
       TAI93At0ZOfGranule, TAI93StartOfGranule, GranuleYear, GranuleMonth, GranuleDay,       &
       GranuleJDay
  USE ozprof_data_module,      ONLY: l2funit, l2swathunit, nGas=>nfgas, nTgas=>ngas,         &
       nLayer=>nlay, ozfit_start_index, ozfit_end_index, aerosol, do_lambcld, saa_flag,      &
       ozwrtavgk, ozwrtcorr, ozwrtcovar, ozwrtcontri, ozwrtres, ozwrtwf, ozwrtsnr,           &
       wrtring, ozwrtvar, gaswrt, ozprof, ozprof_std, ozprof_init, ozprof_ap, ozprof_apstd,  &
       eff_alb, eff_alb_init, nlay, ozdfs, ozinfo, num_iter, covar, contri, avg_kernel,      &
       use_lograd, use_oe, nalb, atmosprof, ntp, nlay_fit, start_layer, end_layer,           &
       the_ctp, the_cfrac, the_cod, lambcld_refl, do_lambcld, which_cld, which_alb,          &
       ozprof_nstd, strataod, stratsca, tropaod, tropsca, aerwavs, maxawin, actawin,         &
       the_cld_flg, fgasidxs, fgaspos, tracegas, saa_flag, nsaa_spike, nsfc, toms_fwhm,      &
       weight_function, do_simu, glintprob, thealbidx, use_effcrs, biascorr, degcorr,        &
       biasfname, degfname, use_logstate, smooth_ozbc, which_clima, which_aperr, norm_tropo3,&
       which_alb, min_serr, min_terr, strat_aerosol, nmom, hres_samprate, pos_alb, pst0,     &
       ozabs_fname, polcorr, use_reg_presgrid, presgrid_fname, use_tropopause, which_toz,    &
       adjust_trop_layer, radc_samprate, radc_lambnd, radc_nsegsr, lambcld_initalb,          &
       loose_aperr, ndiv, VlidortNstream, fit_atanring, use_large_so2_aperr, use_flns,       &
       ring_on_line, ring_convol, which_toz, norm_tropo3, which_aerosol, scale_aod,          &
       scacld_initcod, fixed_ptrop, ozcrs_alb_fname, ntp0, useasy, algorithm_name,           &
       algorithm_version, scaled_aod
  USE he5_l2_fs
  USE OMSAO_pixelcorner_module,ONLY: omi_alllat, omi_alllon, omi_allclat, omi_allclon,       &
       omi_allsza, omi_allvza, omi_allaza, omi_alltime, omi_allSecondsInDay, omi_allSpcftAlt,&
       omi_allSpcftLon, omi_allSpcftLat, omi_allGeoFlg, omi_allMflg, omi_allNSPC
  USE OMSAO_indices_module,    ONLY: n_max_fitpars, ring_idx, solar_idx
  
  USE OMSAO_precision_module 
  USE OMSAO_errstat_module 
  IMPLICIT NONE

  !INTEGER (KIND=i4), EXTERNAL           :: OMI_SMF_setmsg
  INTEGER (KIND=i4), PARAMETER, PRIVATE :: mDim = 12
  INTEGER (KIND=i4), PRIVATE            :: nTimes, nXtrack, XOffset, YOffset, nTimesp1, &
       nXtrackp1, nLayerp1, nOth, mWavel, nAer, nDim, l2_fileid, l2_swid, nFitvar, nCh
  INTEGER (KIND=i4), DIMENSION(mDim)    :: dims
  CHARACTER (LEN=256)                   :: dimList
  
CONTAINS
 
  SUBROUTINE He5_L2WrtInit (spix, epix, sline, eline, pge_error_status )

    INTEGER, INTENT(IN)  :: spix, epix, sline, eline     
    INTEGER, INTENT(OUT) :: pge_error_status  
    INTEGER (KIND=i4)    :: ierr, status
    CHARACTER (LEN = 14) :: modulename = 'he5_output_init'
    
    pge_error_status = pge_errstat_ok

    ! Obtain the dimensions of the output
    nFitvar = n_fitvar_rad; nCh = numwin
    nTimes  = eline - sline + 1;  nTimesp1  = nTimes + 1
    nXtrack = epix  - spix  + 1;  nXtrackp1 = nXtrack + 1
    XOffset = spix;               YOffset   = offset_line
    nLayerp1 = nLayer + 1;        nAer = nCh + 2
    nOth = nFitvar - nGas - (ozfit_end_index - ozfit_start_index + 1)
    mWavel = MAXVAL(omi_nwav_irrad(spix:epix)) - nCh * 2 * radnhtrunc
        
    !! Create the L2 output file
    status = L2_createFile(l2funit, l2_filename )
   
    IF ( status /= omi_s_success ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_creatFile fails.", modulename, 0 )
       pge_error_status = pge_errstat_error; RETURN
    ENDIF
    
    !! Setup the swath in the L2 output file
    nDim = 10
    dimList = "nXtrack,nXtrackp1,nTimes,nTimesp1,nLayer,nLayerp1,nFitvar,nOth,mWavel,nCh"
    dims(1:10) = (/nXtrack,nXtrackp1,nTimes,nTimesp1,nLayer,nLayerp1,nFitvar,nOth,mWavel,nCh/)

    IF (nGas > 0) THEN
       nDim = nDim + 1
       dimList = TRIM(ADJUSTL(dimList)) // ",nGas"; dims(nDim) = nGas
    ENDIF
    IF (aerosol) THEN
       nDim = nDim + 1
       dimList = TRIM(ADJUSTL(dimList)) // ",nAer"; dims(nDim) = nAer
    ENDIF

    status = L2_setupSwath(l2_filename, dimList, dims, l2swathunit, &
         TRIM(ADJUSTL(l2_swathname)), l2_fileid, l2_swid )
    IF ( status /= omi_s_success ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_setupSwath fails.", modulename, 0 )
       pge_error_status = pge_errstat_error; RETURN
    ENDIF
  
    ! Define Geolocation and Data Fields
    CALL DefineFields (nGas, nOth, ozwrtavgk, ozwrtcorr, ozwrtcovar, ozwrtcontri, &
         ozwrtres, ozwrtwf, ozwrtsnr, wrtring, gaswrt, ozwrtvar, do_lambcld, &
         aerosol, saa_flag, reduce_resolution )
    
    RETURN
  END SUBROUTINE He5_L2WrtInit

  SUBROUTINE He5_L2SetGeoFields (spix, epix, pge_error_status )
    INTEGER, INTENT(IN)    :: spix, epix
    INTEGER, INTENT(OUT)   :: pge_error_status  
    TYPE (L2_generic_type) :: GeoBlk, GeoBlk1
    INTEGER (KIND=i4)      :: ierr, status, ig, Ls, Le, bsize, sLine
    INTEGER (KIND=i1)      :: I1
    CHARACTER (LEN = 18)   :: modulename = 'He5_L2SetGeoFields'
    
    pge_error_status = pge_errstat_ok

    !! Define the geofield in the L2 output file
    status = L2_defSWgeofields( l2_swid, gf_o3prof )
    IF( status /= OMI_S_SUCCESS ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_defSWgeofields failed.", modulename, 0 )
       pge_error_status = pge_errstat_error; RETURN
    END IF

    !! Allocate the memory for storing the L2 output information for a block.
    status = L2_newBlockW( GeoBlk, l2_filename, l2_swathname, l2_fileid, &
         l2_swid, gf_o3prof, nTimes )
    IF (status /= OMI_S_SUCCESS ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_newBlockW for Geo failed.", modulename, 0 )
       pge_error_status = pge_errstat_error; RETURN
    END IF
          

    ! Transfer data to the block
    DO ig = 1, GeoBlk%nFields
       bsize = GeoBlk%BlkSize(ig)
       Ls    = GeoBlk%accuBlkSize(ig-1) 
       Le    = Ls + bsize
       Ls    = Ls + 1       !! fortran index scheme

       IF ( ig == 1 ) THEN
          GeoBlk%data( Ls:Le ) = TRANSFER( omi_allgeoflg(spix:epix, 0:nTimes-1), I1, bsize)
       ELSE IF( ig == 2 ) THEN
          GeoBlk%data( Ls:Le ) = TRANSFER( REAL(omi_alllat(spix:epix, 0:nTimes-1), KIND=4), I1, bsize)
       ELSE IF( ig == 3 ) THEN
          GeoBlk%data( Ls:Le ) = TRANSFER( REAL(omi_alllon(spix:epix, 0:nTimes-1), KIND=4), I1, bsize)
       ELSE IF( ig == 4 ) THEN
          GeoBlk%data( Ls:Le ) = TRANSFER( REAL(omi_allsza(spix:epix, 0:nTimes-1), KIND=4), I1, bsize)
       ELSE IF( ig == 5 ) THEN
          GeoBlk%data( Ls:Le ) = TRANSFER( REAL(omi_allvza(spix:epix, 0:nTimes-1), KIND=4), I1, bsize)
       ELSE IF( ig == 6 ) THEN
          GeoBlk%data( Ls:Le ) = TRANSFER( REAL(omi_allaza(spix:epix, 0:nTimes-1), KIND=4), I1, bsize)
       ELSE IF( ig == 7 ) THEN
          GeoBlk%data( Ls:Le ) = TRANSFER( omi_alltime(0:nTimes-1), I1, bsize )
       ELSE IF( ig == 8 ) THEN
          GeoBlk%data( Ls:Le ) = TRANSFER( omi_allSecondsInDay(0:nTimes-1), I1, bsize )
       ELSE IF( ig == 9 ) THEN
          GeoBlk%data( Ls:Le ) = TRANSFER( omi_allSpcftLat(0:nTimes-1), I1, bsize )
       ELSE IF( ig == 10 ) THEN
          GeoBlk%data( Ls:Le ) = TRANSFER( omi_allSpcftLon(0:nTimes-1), I1, bsize )
       ELSE IF( ig == 11 ) THEN
          GeoBlk%data( Ls:Le ) = TRANSFER( omi_allSpcftAlt(0:nTimes-1), I1, bsize )
       ENDIF
    ENDDO

    sLine = 0
    status = L2_writeBlock (GeoBlk, sLine, nTimes) 
    IF (status /= OMI_S_SUCCESS ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_writeBlockW for Geo failed.", modulename, 0 )
       pge_error_status = pge_errstat_error; RETURN
    END IF

    ! Deallocate the block
    CALL L2_disposeBlockW( GeoBlk  )

    !! Define the geofield in the L2 output file
    status = L2_defSWgeofields( l2_swid, gf1_o3prof )
    IF( status /= OMI_S_SUCCESS ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_defSWgeofields failed.", modulename, 0 )
       pge_error_status = pge_errstat_error; RETURN
    END IF

    !! Allocate the memory for storing the L2 output information for a block.
    status = L2_newBlockW( GeoBlk1, l2_filename, l2_swathname, l2_fileid, &
         l2_swid, gf1_o3prof, nTimes + 1 )
    IF (status /= OMI_S_SUCCESS ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_newBlockW for Geo failed.", modulename, 0 )
       pge_error_status = pge_errstat_error; RETURN
    END IF
   
    ! Transfer data to the block
    DO ig = 1, GeoBlk1%nFields
       bsize = GeoBlk1%BlkSize(ig)
       Ls    = GeoBlk1%accuBlkSize(ig-1) 
       Le    = Ls + bsize
       Ls    = Ls + 1       !! fortran index scheme
 
       IF ( ig == 1 ) THEN
          GeoBlk1%data( Ls:Le ) = TRANSFER( REAL(omi_allclat(spix-1:epix, 0:nTimes), KIND=4), I1, bsize)
       ELSE IF( ig == 2 ) THEN
          GeoBlk1%data( Ls:Le ) = TRANSFER( REAL(omi_allclon(spix-1:epix, 0:nTimes), KIND=4), I1, bsize)
       ENDIF
    ENDDO
    
    sLine = 0
    status = L2_writeBlock (GeoBlk1, sLine, nTimesp1) 
    IF (status /= OMI_S_SUCCESS ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_writeBlockW for Geo failed.", modulename, 0 )
       pge_error_status = pge_errstat_error; RETURN
    END IF

    ! Deallocate the block
    CALL L2_disposeBlockW( GeoBlk1  )
        
    RETURN

  END SUBROUTINE He5_L2SetGeoFields

  ! currpix: curr x-track position (spatially-coadded)
  ! spix:    first position to be processed
  ! epix:    last position to be processed
  ! iloop:   curr scan line (0 based) within a block
  ! iline:   curr scan line (among all lines to be processed)
  ! nloop:   number of scan lines within a block
  SUBROUTINE He5_L2SetDataFields (currpix, spix, epix, iloop, iline, nloop, &
       exval, fitcol, dfitcol, pge_error_status )
    INTEGER, INTENT(IN)    :: currpix, spix, epix, iloop, iline, nloop, exval
    INTEGER, INTENT(OUT)   :: pge_error_status  
    REAL(KIND=dp), DIMENSION(3), INTENT (IN)    :: fitcol
    REAL(KIND=dp), DIMENSION(3, 2), INTENT (IN) :: dfitcol
    TYPE (L2_generic_type) :: DataBlk, WavBlk, L1bBlk
    INTEGER (KIND=i4)      :: ierr, status, ix, Ls, Le, bsize, id, errstat
    INTEGER (KIND=1)       :: I1
    CHARACTER (LEN = 19)   :: modulename = 'He5_L2SetDataFields'   
    LOGICAL, SAVE          :: first = .TRUE.
    INTEGER (KIND=i4), SAVE:: sWrtLine = 0
       
    pge_error_status = pge_errstat_ok

    IF ( first ) THEN

       ! Write common wavelength grid
       IF (reduce_resolution) THEN
          status = L2_defSWdatafields( l2_swid, df_wav )
          IF( status /= OMI_S_SUCCESS ) THEN
             ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_defSWdatafields failed.", modulename, 0 )
             pge_error_status = pge_errstat_error; RETURN
          END IF
          
          !! Allocate the memory for storing the L2 output information for a block.
          status = L2_newBlockW( WavBlk, l2_filename, l2_swathname, l2_fileid, &
               l2_swid, df_wav, mWavel)
          IF (status /= OMI_S_SUCCESS ) THEN
             ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_newBlockW for WavBlk failed.", modulename, 0 )
             pge_error_status = pge_errstat_error; RETURN
          ENDIF
          
          WavBlk%data(:) = TRANSFER( REAL(fitwavs(1:n_rad_wvl), KIND=4), &
               I1, n_rad_wvl * WavBlk%elmSize(1))

          status = L2_writeBlock ( WavBlk, sWrtLine, n_rad_wvl) 
          IF (status /= OMI_S_SUCCESS ) THEN
             ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_writeBlockW for WavBlk failed.", modulename, 0 )
             pge_error_status = pge_errstat_error; RETURN
          ENDIF

          CALL L2_disposeBlockW( WavBlk  )
       ENDIF

       ! Write MeasurementQualityFlags and NumberSmallPixelColumns
       status = L2_defSWdatafields( l2_swid, df_l1b )
       IF( status /= OMI_S_SUCCESS ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_defSWdatafields failed.", modulename, 0 )
          pge_error_status = pge_errstat_error; RETURN
       ENDIF
       
       !! Allocate the memory for storing the L2 output information for a block.
       status = L2_newBlockW( L1bBlk, l2_filename, l2_swathname, l2_fileid, &
            l2_swid, df_l1b, nTimes)
       IF (status /= OMI_S_SUCCESS ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_newBlockW for L1B fields failed.", modulename, 0 )
          pge_error_status = pge_errstat_error; RETURN
       ENDIF

       ! Transfer data to the block
       DO id = 1, L1bBlk%nFields
          bsize = L1bBlk%BlkSize(id)
          Ls    = L1bBlk%accuBlkSize(id-1) 
          Le    = Ls + bsize
          Ls    = Ls + 1       !! fortran index scheme
          
          IF ( id == 1 ) THEN
             L1bBlk%data( Ls:Le ) = TRANSFER( INT(omi_allMflg(0:nTimes-1), KIND=1), I1, bsize)
          ELSE IF( id == 2 ) THEN
             L1bBlk%data( Ls:Le ) = TRANSFER( omi_allNSPC(0:nTimes-1), I1, bsize)
          ENDIF
       ENDDO
    
       status = L2_writeBlock (L1bBlk, sWrtLine, nTimes) 
       IF (status /= OMI_S_SUCCESS ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_writeBlockW for L1B fields failed.", modulename, 0 )
          pge_error_status = pge_errstat_error; RETURN
       END IF

       ! Deallocate the block
       CALL L2_disposeBlockW( L1bBlk  )
       
       !! Define the data field in the L2 output file
       status = L2_defSWdatafields( l2_swid, df_o3prof(1:nDataF) )
       IF( status /= OMI_S_SUCCESS ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_defSWdatafields failed.", modulename, 0 )
          pge_error_status = pge_errstat_error; RETURN
       END IF
     
       !! Allocate the memory for storing the L2 output information for a block.
       status = L2_newBlockW( DataBlk, l2_filename, l2_swathname, l2_fileid, &
            l2_swid, df_o3prof(1:nDataF), nloop)
       IF (status /= OMI_S_SUCCESS ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_newBlockW for Data failed.", modulename, 0 )
          pge_error_status = pge_errstat_error; RETURN
       END IF
       
       first = .FALSE.
    ENDIF

    ix = currpix - spix + 1 
    IF (exval >= 0) THEN  ! Retrieval finished.
       CALL He5_L2SetDataPix(DataBlk, ix, iloop, exval, fitcol, dfitcol)       
    ELSE                  ! Retrieval not finished. Filling in them as missing values.
       CALL He5_L2FillDataPix(DataBlk, ix, iloop, exval)
    ENDIF
    
    ! Write a block
    IF (currpix == epix .AND. iloop == nloop - 1) THEN
       status = L2_writeBlock (DataBlk, sWrtLine, nloop) 
       IF (status /= OMI_S_SUCCESS ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_FAILURE, "L2_writeBlockW for Data failed.", modulename, 0 )
          pge_error_status = pge_errstat_error; RETURN
       ENDIF
       sWrtLine = sWrtLine + nloop
    ENDIF

    ! End of processing, deallocate the block
    IF ( currpix == epix .AND. iline == nTimes - 1 ) THEN 
       CALL L2_disposeBlockW( DataBlk  )

       CALL He5_L2WrtClose ( errstat )
       IF ( errstat /= pge_errstat_ok ) THEN
          WRITE(*, *) modulename, ' : Cannot close HE5 output file!!!'
          pge_error_status = pge_errstat_error; RETURN
       END IF

       ! Write swath attributes
       CALL he5_writeAttribute( errstat )

       IF ( errstat /= pge_errstat_ok ) THEN
          WRITE(*, *) modulename, ' : Cannot finish writing the attributes!!!'
          pge_error_status = pge_errstat_error; RETURN
       END IF
    ENDIF

    RETURN

  END SUBROUTINE He5_L2SetDataFields

  ! ix: start X-track position within a block (1 based)
  ! iy: start along-track position within a block (0-based)
  SUBROUTINE He5_L2SetDataPix (DataBlk, ix, iy, exval, fitcol, dfitcol)

    IMPLICIT NONE
    TYPE (L2_generic_type), INTENT (INOUT)       :: DataBlk
    INTEGER, INTENT(IN)                          :: ix, iy, exval  
    REAL (KIND=dp), DIMENSION(3), INTENT (IN)    :: fitcol
    REAL (KIND=dp), DIMENSION(3, 2), INTENT (IN) :: dfitcol

    REAL (KIND=dp), DIMENSION (maxwin)           :: allrms, allavgres
    REAL (KIND=dp), DIMENSION (5, n_max_fitpars) :: tempvar
    REAL (KIND=dp), DIMENSION (n_max_fitpars, n_max_fitpars) :: correl
    REAL (KIND=dp), DIMENSION (max_fit_pts)      :: tempring

    REAL (KIND=dp)                               :: avgres
    INTEGER (KIND=i4)      :: id, Ls, Le, bsize, LL, fidx, lidx, i, fid, j
    INTEGER (KIND=i1)      :: I1
    CHARACTER (LEN = 16)   :: modulename = 'He5_L2SetDataPix'  

    fidx = 1
    DO i = 1, nCh
       lidx = fidx + nradpix(i) - 1
       allrms(i) = SQRT(SUM((fitres_rad(fidx:lidx)/fitweights(fidx:lidx))**2) / nradpix(i))
       fidx = lidx + 1
    ENDDO
    simspec_rad(1:n_rad_wvl) = fitspec_rad(1:n_rad_wvl) - fitres_rad(1:n_rad_wvl)
    IF (use_lograd) THEN
       fitspec_rad(1:n_rad_wvl) = EXP(fitspec_rad(1:n_rad_wvl))
       simspec_rad(1:n_rad_wvl) = EXP(simspec_rad(1:n_rad_wvl))
       fitres_rad(1:n_rad_wvl) = fitspec_rad(1:n_rad_wvl) -  simspec_rad(1:n_rad_wvl)
    END IF
    avgres = SQRT(SUM((ABS(fitres_rad(1:n_rad_wvl)) / &
         fitspec_rad(1:n_rad_wvl))**2.0)/n_rad_wvl)*100.0
    
    fidx = 1
    DO i = 1, nCh
       lidx = fidx + nradpix(i) - 1
       allavgres(i) = SQRT(SUM((ABS(fitres_rad(fidx:lidx)) / &
            fitspec_rad(fidx:lidx))**2.0)/nradpix(i))*100.0
       fidx = lidx + 1
    ENDDO   
    
    LL    = iy * nXtrack + (ix - 1)     !! 0 based index
    ! First write necessary Data fields
    DO id = 1, 31     
       bsize = DataBlk%pixSize(id)
       Ls    = DataBlk%accuBlkSize(id-1) + LL * bsize
       Le    = Ls + bsize
       Ls    = Ls + 1       !! fortran index scheme

       IF (id == 1) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(atmosprof(1, 0:nlay), KIND=4), I1, bsize )
       ELSE IF (id == 2) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(atmosprof(2, 0:nlay), KIND=4), I1, bsize )
       ELSE IF (id == 3) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(atmosprof(3, 0:nlay), KIND=4), I1, bsize )
       ELSE IF (id == 4) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(ozprof_ap(1:nlay), KIND=4), I1, bsize )
       ELSE IF (id == 5) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(ozprof_apstd(1:nlay), KIND=4), I1, bsize )
       ELSE IF (id == 6) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(ozprof(1:nlay), KIND=4), I1, bsize )
       ELSE IF (id == 7) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(ozprof_nstd(1:nlay), KIND=4), I1, bsize )
       ELSE IF (id == 8) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(ozprof_std(1:nlay), KIND=4), I1, bsize )
       ELSE IF (id == 9) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(fitcol(1:1), KIND=4), I1, bsize )
       ELSE IF (id == 10) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(dfitcol(1:1, 2:2), KIND=4), I1, bsize )
       ELSE IF (id == 11) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(dfitcol(1:1, 1:1), KIND=4), I1, bsize )
       ELSE IF (id == 12) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(fitcol(2:2), KIND=4), I1, bsize )
       ELSE IF (id == 13) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(dfitcol(2:2, 2:2), KIND=4), I1, bsize )
       ELSE IF (id == 14) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(dfitcol(2:2, 1:1), KIND=4), I1, bsize )
       ELSE IF (id == 15) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(fitcol(3:3), KIND=4), I1, bsize )
       ELSE IF (id == 16) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(dfitcol(3:3, 2:2), KIND=4), I1, bsize )
       ELSE IF (id == 17) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(dfitcol(3:3, 1:1), KIND=4), I1, bsize )
       ELSE IF (id == 18) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( INT((/ntp/), KIND=1), I1, bsize )
       ELSE IF (id == 19) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( INT((/nsfc/), KIND=1), I1, bsize )
       ELSE IF (id == 20) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( INT((/exval/), KIND=1), I1, bsize )
       ELSE IF (id == 21) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( INT((/num_iter/), KIND=1), I1, bsize )
       ELSE IF (id == 22) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL((/ozinfo/), KIND=4), I1, bsize )
       ELSE IF (id == 23) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL((/the_cfrac/), KIND=4), I1, bsize )
       ELSE IF (id == 24) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL((/the_ctp/), KIND=4), I1, bsize )
       ELSE IF (id == 25) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( INT((/the_cld_flg/), KIND=1), I1, bsize )
       ELSE IF (id == 26) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL((/eff_alb(thealbidx)/), KIND=4), I1, bsize )
       ELSE IF (id == 27) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL((/glintprob/), KIND=4), I1, bsize )
       ELSE IF (id == 28) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(allrms(1:nCh), KIND=4), I1, bsize )
       ELSE IF (id == 29) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(allavgres(1:nCh), KIND=4), I1, bsize )
       ELSE IF (id == 30) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( INT((/n_rad_wvl/), KIND=2), I1, bsize )
       ELSE IF (id == 31) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( INT(nradpix(1:nCh), KIND=2), I1, bsize )
       ENDIF
    ENDDO

    ! Write trace gas if selected
    fid = 31
    IF ( nGas > 0 .AND. (gaswrt .OR. ozwrtvar) ) THEN
       tempvar(1, 1:nGas) = tracegas(fgaspos(1:nGas), 2)
       tempvar(2, 1:nGas) = tracegas(fgaspos(1:nGas), 3)
       tempvar(3, 1:nGas) = tracegas(fgaspos(1:nGas), 4)
       tempvar(4, 1:nGas) = tracegas(fgaspos(1:nGas), 6)
       tempvar(5, 1:nGas) = tracegas(fgaspos(1:nGas), 5)

       DO id = fid + 1, fid + 5 
          bsize = DataBlk%pixSize(id)
          Ls    = DataBlk%accuBlkSize(id-1) + LL * bsize
          Le    = Ls + bsize
          Ls    = Ls + 1       
         
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(tempvar(id-fid, 1:nGas), KIND=4), I1, bsize )
       ENDDO  

       fid = fid + 5
    ENDIF

    ! Write other variables if selected
    IF (nOth > 0 .AND. ozwrtvar) THEN
       tempvar(1, 1:nOth) = fitvar_rad_apriori(mask_fitvar_rad(fothvarpos(1:nOth)))
       tempvar(2, 1:nOth) = fitvar_rad_aperror(mask_fitvar_rad(fothvarpos(1:nOth)))       
       tempvar(3, 1:nOth) = fitvar_rad(mask_fitvar_rad(fothvarpos(1:nOth)))
       tempvar(4, 1:nOth) = fitvar_rad_nstd(mask_fitvar_rad(fothvarpos(1:nOth)))
       tempvar(5, 1:nOth) = fitvar_rad_std(mask_fitvar_rad(fothvarpos(1:nOth)))

       DO id = fid + 1, fid + 5
          bsize = DataBlk%pixSize(id)
          Ls    = DataBlk%accuBlkSize(id-1) + LL * bsize
          Le    = Ls + bsize
          Ls    = Ls + 1       
          
          DataBlk%data( Ls:Le ) = TRANSFER( REAL(tempvar(id-fid, 1:nOth), KIND=4), I1, bsize )
       ENDDO  

       fid = fid + 5
    ENDIF

    IF (ozwrtavgk) THEN
       bsize = DataBlk%pixSize(fid + 1)
       Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
       Le    = Ls + bsize
       Ls    = Ls + 1       

       i = ozfit_start_index; j = ozfit_end_index
       DataBlk%data( Ls:Le ) = TRANSFER( REAL(avg_kernel(i:j, i:j), KIND=4), I1, bsize )
       fid = fid + 1
    ENDIF

    IF (ozwrtcorr) THEN
       DO i = 1, nFitvar
          correl(i, i) = 1.0
          DO j = 1, i - 1
             correl(i, j) = covar(i, j) / SQRT(covar(i, i) * covar(j, j))
             correl(j, i) = correl(i, j)
          ENDDO
       ENDDO

   
       bsize = DataBlk%pixSize(fid + 1)
       Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
       Le    = Ls + bsize
       Ls    = Ls + 1       

       DataBlk%data( Ls:Le ) = TRANSFER( REAL(correl(1:nFitvar, 1:nFitvar), KIND=4), I1, bsize )
       fid = fid + 1
    ENDIF

    IF (ozwrtcovar) THEN  
       bsize = DataBlk%pixSize(fid + 1)
       Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
       Le    = Ls + bsize
       Ls    = Ls + 1       

       i = ozfit_start_index; j = ozfit_end_index
       DataBlk%data( Ls:Le ) = TRANSFER( REAL(covar(i:j, i:j), KIND=4), I1, bsize )
       fid = fid + 1
    ENDIF

    IF (ozwrtcontri) THEN  
       bsize = DataBlk%pixSize(fid + 1)
       Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
       Le    = Ls + bsize
       Ls    = Ls + 1       

       IF (mWavel > n_rad_wvl) THEN
          contri(1:nFitvar, n_rad_wvl+1:mWavel) = 0.D0
       ENDIF

       DataBlk%data( Ls:Le ) = TRANSFER( REAL(contri(1:nFitvar, 1:mWavel), KIND=4), I1, bsize )
       fid = fid + 1
    ENDIF

    IF (ozwrtwf) THEN  
       bsize = DataBlk%pixSize(fid + 1)
       Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
       Le    = Ls + bsize
       Ls    = Ls + 1       

       IF (mWavel > n_rad_wvl) THEN
          weight_function(n_rad_wvl+1:mWavel, 1:nFitvar) = 0.D0
       ENDIF

       !WRITE(*, *) mWavel, n_rad_wvl
       !DO j = 1, mWavel
       !   WRITE(90, '(40D16.7)') fitwavs(j), weight_function(j, 1:nFitvar)
       !ENDDO

       DataBlk%data( Ls:Le ) = TRANSFER( REAL(weight_function(1:mWavel, 1:nFitvar), KIND=4), I1, bsize )
       fid = fid + 1
    ENDIF

     IF ( (.NOT. reduce_resolution) .AND. &
          (ozwrtcontri .OR. ozwrtres .OR. ozwrtwf .OR. ozwrtsnr .OR. wrtring)) THEN
        IF (mWavel > n_rad_wvl) THEN
           fitwavs(n_rad_wvl + 1 : mWavel) = 0.D0
        ENDIF
  
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1       
        
        DataBlk%data( Ls:Le ) = TRANSFER( REAL(fitwavs(1:mWavel), KIND=4), I1, bsize )
        fid = fid + 1
     ENDIF
    
     IF (ozwrtres) THEN
        IF (mWavel > n_rad_wvl) THEN
           fitspec_rad(n_rad_wvl + 1 : mWavel) = 0.D0
           simspec_rad(n_rad_wvl + 1 : mWavel) = 0.D0
        ENDIF
   
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1       
        
        DataBlk%data( Ls:Le ) = TRANSFER( REAL(simspec_rad(1:mWavel), KIND=4), I1, bsize )
        fid = fid + 1
 
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1       
        
        DataBlk%data( Ls:Le ) = TRANSFER( REAL(fitspec_rad(1:mWavel), KIND=4), I1, bsize )
        fid = fid + 1
     ENDIF

     IF (ozwrtsnr) THEN
        IF (mWavel > n_rad_wvl) THEN
           fitweights(n_rad_wvl + 1 : mWavel) = 0.D0
        ENDIF
 
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1       
        
        DataBlk%data( Ls:Le ) = TRANSFER( REAL(fitweights(1:mWavel), KIND=4), I1, bsize )
        fid = fid + 1
     ENDIF

     IF (wrtring) THEN
        tempring(1:n_rad_wvl) = database(ring_idx, refidx(1:n_rad_wvl))
        IF (mWavel > n_rad_wvl) THEN
           tempring(n_rad_wvl + 1 : mWavel) = 0.D0
        ENDIF
   
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1        
        
        DataBlk%data( Ls:Le ) = TRANSFER( REAL(tempring(1:mWavel), KIND=4), I1, bsize )
        fid = fid + 1
     ENDIF

     IF (.NOT. do_lambcld) THEN  
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1 
        
        DataBlk%data( Ls:Le ) = TRANSFER( REAL((/the_cod/), KIND=4), I1, bsize )
        fid = fid + 1
     ENDIF

     IF (.NOT. saa_flag) THEN   
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1 
        
        DataBlk%data( Ls:Le ) = TRANSFER( INT((/nsaa_spike/), KIND=2), I1, bsize )
        fid = fid + 1
     ENDIF

     IF (aerosol) THEN  
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1 
        DataBlk%data( Ls:Le ) = TRANSFER( REAL(aerwavs(1:nAer), KIND=4), I1, bsize )
        fid = fid + 1
  
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1
        DataBlk%data( Ls:Le ) = TRANSFER( REAL(tropaod(1:nAer), KIND=4), I1, bsize )
        fid = fid + 1
   
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1 
        DataBlk%data( Ls:Le ) = TRANSFER( REAL(tropsca(1:nAer), KIND=4), I1, bsize )
        fid = fid + 1        
     ENDIF

    RETURN
    
  END SUBROUTINE He5_L2SetDataPix
  
  SUBROUTINE He5_L2FillDataPix (DataBlk, ix, iy, exval)
    
    IMPLICIT NONE
    TYPE (L2_generic_type), INTENT (INOUT) :: DataBlk
    INTEGER, INTENT(IN)                    :: ix, iy, exval
    INTEGER (KIND=i4)      :: id, Ls, Le, bsize, LL, fid, i, j
    INTEGER (KIND=i1)      :: I1
    CHARACTER (LEN = 17)   :: modulename = 'He5_L2FillDataPix' 

    REAL (KIND=r4), DIMENSION(0:maxlay)      :: tmp1D_layer
    REAL (KIND=r4), DIMENSION(n_max_fitpars) :: tmp1D_fitvar
    REAL (KIND=r4), DIMENSION(max_fit_pts)   :: tmp1D_fitpts
    REAL (KIND=r4), DIMENSION(maxwin)        :: tmp1D_numwin
    INTEGER (KIND=2), DIMENSION(maxwin)      :: tmp1D_num
    REAL (KIND=r4), DIMENSION(maxwin+2)      :: tmp1D_aer
    REAL (KIND=r4), DIMENSION(n_max_fitpars, n_max_fitpars) :: tmp2D_fitvar
    REAL (KIND=r4), DIMENSION(n_max_fitpars, max_fit_pts)   :: tmp2D_contri
    REAL (KIND=r4), DIMENSION(max_fit_pts, n_max_fitpars)   :: tmp2D_wf
    
    tmp1D_layer(0:nlay)     = fill_float32
    tmp1D_numwin(1:nCh)     = fill_float32
    tmp1D_num(1:nCh)        = fill_uint16
    tmp1D_fitvar(1:nFitvar) = fill_float32
    tmp1D_fitpts(1:mWavel)  = fill_float32

    LL    = iy * nXtrack + (ix - 1)     !! 0 based index
    ! First write necessary Data fields
    DO id = 1, 31    
       bsize = DataBlk%pixSize(id)
       Ls    = DataBlk%accuBlkSize(id-1) + LL * bsize
       Le    = Ls + bsize
       Ls    = Ls + 1       !! fortran index scheme

       IF (id == 1) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(0:nlay), I1, bsize )
       ELSE IF (id == 2) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(0:nlay), I1, bsize )
       ELSE IF (id == 3) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(0:nlay), I1, bsize )
       ELSE IF (id == 4) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1:nlay), I1, bsize )
       ELSE IF (id == 5) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1:nlay), I1, bsize )
       ELSE IF (id == 6) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1:nlay), I1, bsize )
       ELSE IF (id == 7) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1:nlay), I1, bsize )
       ELSE IF (id == 8) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1:nlay), I1, bsize )
       ELSE IF (id == 9) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1), I1, bsize )
       ELSE IF (id == 10) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1), I1, bsize )
       ELSE IF (id == 11) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1), I1, bsize )
       ELSE IF (id == 12) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1), I1, bsize )
       ELSE IF (id == 13) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1), I1, bsize )
       ELSE IF (id == 14) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1), I1, bsize )
       ELSE IF (id == 15) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1), I1, bsize )
       ELSE IF (id == 16) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1), I1, bsize )
       ELSE IF (id == 17) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1), I1, bsize )
       ELSE IF (id == 18) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( fill_uint8, I1, bsize )
       ELSE IF (id == 19) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( fill_uint8, I1, bsize )
       ELSE IF (id == 20) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( INT((/exval/), KIND=1), I1, bsize )
       ELSE IF (id == 21) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( fill_uint8, I1, bsize )
       ELSE IF (id == 22) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1), I1, bsize )
       ELSE IF (id == 23) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1), I1, bsize )
       ELSE IF (id == 24) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1), I1, bsize )
       ELSE IF (id == 25) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( fill_uint8, I1, bsize )
       ELSE IF (id == 26) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1), I1, bsize )
       ELSE IF (id == 27) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1), I1, bsize )
       ELSE IF (id == 28) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_numwin(1:nCh), I1, bsize )
       ELSE IF (id == 29) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_numwin(1:nCh), I1, bsize )
       ELSE IF (id == 30) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( fill_uint16, I1, bsize )
       ELSE IF (id == 31) THEN
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_num(1:nCh), I1, bsize )
       ENDIF
    ENDDO

    ! Write trace gas if selected
    fid = 31
    IF ( nGas > 0 .AND. (gaswrt .OR. ozwrtvar) ) THEN
       DO id = fid + 1, fid + 5 
          bsize = DataBlk%pixSize(id)
          Ls    = DataBlk%accuBlkSize(id-1) + LL * bsize
          Le    = Ls + bsize
          Ls    = Ls + 1       
         
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_fitvar(1:nGas), I1, bsize )
       ENDDO  

       fid = fid + 5
    ENDIF

    ! Write other variables if selected
    IF (nOth > 0 .AND. ozwrtvar) THEN
       DO id = fid + 1, fid + 5
          bsize = DataBlk%pixSize(id)
          Ls    = DataBlk%accuBlkSize(id-1) + LL * bsize
          Le    = Ls + bsize
          Ls    = Ls + 1       
          
          DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_fitvar(1:nOth), I1, bsize )
       ENDDO  

       fid = fid + 5
    ENDIF

    IF (ozwrtavgk) THEN
       bsize = DataBlk%pixSize(fid + 1)
       Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
       Le    = Ls + bsize
       Ls    = Ls + 1       

       i = ozfit_start_index; j = ozfit_end_index
       tmp2D_fitvar(i:j, i:j) = fill_float32
       DataBlk%data( Ls:Le ) = TRANSFER( tmp2D_fitvar(i:j, i:j), I1, bsize )
       fid = fid + 1
    ENDIF

    IF (ozwrtcorr) THEN
       tmp2D_fitvar(1:nFitvar, 1:nFitvar) = fill_float32   
       bsize = DataBlk%pixSize(fid + 1)
       Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
       Le    = Ls + bsize
       Ls    = Ls + 1       

       DataBlk%data( Ls:Le ) = TRANSFER( tmp2D_fitvar(1:nFitvar, 1:nFitvar), I1, bsize )
       fid = fid + 1
    ENDIF

    IF (ozwrtcovar) THEN  
       bsize = DataBlk%pixSize(fid + 1)
       Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
       Le    = Ls + bsize
       Ls    = Ls + 1       

       i = ozfit_start_index; j = ozfit_end_index
       tmp2D_fitvar(i:j, i:j) = fill_float32
       DataBlk%data( Ls:Le ) = TRANSFER( tmp2D_fitvar(i:j, i:j), I1, bsize )
       fid = fid + 1
    ENDIF

    IF (ozwrtcontri) THEN  
       bsize = DataBlk%pixSize(fid + 1)
       Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
       Le    = Ls + bsize
       Ls    = Ls + 1       

       tmp2D_contri(1:nFitvar, 1:mWavel) = fill_float32
       DataBlk%data( Ls:Le ) = TRANSFER( tmp2D_contri(1:nFitvar, 1:mWavel), I1, bsize )
       fid = fid + 1
    ENDIF

    IF (ozwrtwf) THEN  
       bsize = DataBlk%pixSize(fid + 1)
       Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
       Le    = Ls + bsize
       Ls    = Ls + 1       
       
       tmp2D_wf(1:mWavel, 1:nFitvar) = fill_float32

       DataBlk%data( Ls:Le ) = TRANSFER( tmp2D_wf(1:mWavel, 1:nFitvar), I1, bsize )
       fid = fid + 1
    ENDIF

     IF ( (.NOT. reduce_resolution) .AND. &
          (ozwrtcontri .OR. ozwrtres .OR. ozwrtwf .OR. ozwrtsnr .OR. wrtring)) THEN  
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1       
        
        DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_fitpts(1:mWavel), I1, bsize )
        fid = fid + 1
     ENDIF
    
     IF (ozwrtres) THEN  
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1       
        
        DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_fitpts(1:mWavel), I1, bsize )
        fid = fid + 1
 
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1       
        
        DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_fitpts(1:mWavel), I1, bsize )
        fid = fid + 1
     ENDIF

     IF (ozwrtsnr) THEN 
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1       
        
        DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_fitpts(1:mWavel), I1, bsize )
        fid = fid + 1
     ENDIF

     IF (wrtring) THEN 
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1        
        
        DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_fitpts(1:mWavel), I1, bsize )
        fid = fid + 1
     ENDIF

     IF (.NOT. do_lambcld) THEN  
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1 
        
        DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_layer(1), I1, bsize )
        fid = fid + 1
     ENDIF

     IF (.NOT. saa_flag) THEN   
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1 
        
        DataBlk%data( Ls:Le ) = TRANSFER( fill_uint8, I1, bsize )
        fid = fid + 1
     ENDIF

     IF (aerosol) THEN  
        tmp1D_aer(1:nAer) = fill_float32
        
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1 
        DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_aer(1:nAer), I1, bsize )
        fid = fid + 1
        
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1
        DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_aer(1:nAer), I1, bsize )
        fid = fid + 1
        
        bsize = DataBlk%pixSize(fid + 1)
        Ls    = DataBlk%accuBlkSize(fid) + LL * bsize
        Le    = Ls + bsize
        Ls    = Ls + 1 
        DataBlk%data( Ls:Le ) = TRANSFER( tmp1D_aer(1:nAer), I1, bsize )
        fid = fid + 1        
     ENDIF
     
     RETURN    
  END SUBROUTINE He5_L2FillDataPix
  
  SUBROUTINE He5_writeAttribute( errstat )

    INTEGER (KIND=4), INTENT(OUT)                  :: errstat
    INTEGER (KIND=4)                               :: ierr, status, count, i
    CHARACTER (LEN = PGS_SMF_MAX_MSG_SIZE  )       :: msg
    CHARACTER (LEN = 4)                            :: tempc
    CHARACTER (LEN = 8), DIMENSION(n_max_fitpars)  :: varnames
    CHARACTER (LEN = 4), DIMENSION(n_max_fitpars)  :: varnames_nNum
    CHARACTER (LEN = 20), DIMENSION(n_max_fitpars) :: units
    CHARACTER (LEN = 256)                          :: tempc1
    CHARACTER (LEN = 1000)                         :: varnamec, gasnamec, othnamec, othunitc
    CHARACTER (LEN = 18)                           :: modulename='He5_writeAttribute'

    errstat = pge_errstat_ok

    l2_fileid = he5_swopen( l2_filename, HE5F_ACC_RDWR )
    IF( l2_fileid == -1 ) THEN
       WRITE( msg,'(A)' ) "he5_swopen:"// TRIM(l2_filename) // " failed."
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, msg, "He5_writeSwathAttribute", 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    ! Write global attributes (algorithm options)
    tempc1 = algorithm_name
    count=LEN(TRIM(tempc1))
    status = he5_ehwrglatt( l2_fileid, "AlgorithmName", HE5T_NATIVE_CHAR, count, tempc1 )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute AlgorithmName failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc1 = algorithm_version
    count=LEN(TRIM(tempc1))
    status = he5_ehwrglatt( l2_fileid, "AlgorithmVersion", HE5T_NATIVE_CHAR, count, tempc1 )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute AlgorithmVersion failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "OrbitNumber", HE5T_NATIVE_UINT16, 1, INT(orbnum, KIND=2))
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute OrbitNumber failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    
    IF (.NOT. use_backup) THEN
       tempc1 = l1b_irrad_filename; count=LEN(TRIM(tempc1))
       status = he5_ehwrglatt( l2_fileid, "L1_SolarIrradianceFilename", &
            HE5T_NATIVE_CHAR, count, tempc1 )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, &
               "Write Global Attribute L1_SolarIrradianceFilename failed.", modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF
    ENDIF
    
    tempc1 = l1b_rad_filename; count=LEN(TRIM(tempc1))
    status = he5_ehwrglatt( l2_fileid, "L1_RadianceFilename", HE5T_NATIVE_CHAR, count, tempc1)
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute L1_RadianceFilename failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    
    tempc1 = l2_cld_filename; count=LEN(TRIM(tempc1))
    status = he5_ehwrglatt( l2_fileid, "L2_CloudFilename", HE5T_NATIVE_CHAR, count, tempc1 )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute L2_CloudFilename failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "GranuleYear", HE5T_NATIVE_INT16, 1, INT(GranuleYear, KIND=2) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute GranuleYear failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "GranuleMonth", HE5T_NATIVE_UINT8, 1, INT(GranuleMonth, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute GranuleMonth failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "GranuleDay", HE5T_NATIVE_UINT8, 1, INT(GranuleDay, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute GranuleMonth failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "GranuleJulDay", HE5T_NATIVE_INT16, 1, INT(GranuleJDay, KIND=2) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute GranuleJulDay failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "TAI93At0ZOfGranule", &
         HE5T_NATIVE_DOUBLE, 1, TAI93At0ZOfGranule )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute TAI93At0ZOfGranule failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "TAI93StartOfGranule", &
         HE5T_NATIVE_DOUBLE, 1, TAI93StartOfGranule )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute TAI93StartOfGranule failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    
    tempc = '.T.'; IF (.NOT. reduce_resolution) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "reduce_resolution", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute reduce_resolution failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    
    WRITE(tempc1, '(I5, A90)') which_slit, &
         ' 0: Symmetric Gaussian 1: Gaussian 2: Voigt 3: Triagnular 4: OMI slit'
    count = LEN(TRIM(tempc1))
    status = he5_ehwrglatt( l2_fileid, "which_slit", HE5T_NATIVE_CHAR, count, tempc1 )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute which_slit failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    
    IF (reduce_resolution) THEN
       status = he5_ehwrglatt( l2_fileid, "AdditionalSlitWidths", HE5T_NATIVE_FLOAT, &
            mswath, omi_redslw(1:mswath) )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, &
               "Write Global Attribute AdditionalSlitWidths failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF
       
       tempc = '.T.'; IF (.NOT. use_redfixwav) tempc = '.F.'
       count=LEN(TRIM(tempc))
       status = he5_ehwrglatt( l2_fileid, "UsedFixedWavelengths", HE5T_NATIVE_CHAR, count, tempc )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute UseFixedWavelengths failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF

       IF (use_redfixwav) THEN
          tempc1 = redfixwav_fname
          count=LEN(TRIM(tempc1))
          status = he5_ehwrglatt( l2_fileid, "FixedWavelengthFilename", HE5T_NATIVE_CHAR, count, tempc1 )
          IF( status == -1 ) THEN
             ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute FixedWavelengthFilename failed.", &
                  modulename, 0 )
             errstat = pge_errstat_error; RETURN
          ENDIF
       ELSE
          status = he5_ehwrglatt( l2_fileid, "ReducedSampRate", &
               HE5T_NATIVE_FLOAT, 1, REAL(redsampr, KIND=4) )
          IF( status == -1 ) THEN
             ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute ReducedSampRate failed.", &
                  modulename, 0 )
             errstat = pge_errstat_error; RETURN
          ENDIF
          
          status = he5_ehwrglatt( l2_fileid, "ReducedFineGridDw", &
               HE5T_NATIVE_FLOAT, 1, REAL(redlam, KIND=4) )
          IF( status == -1 ) THEN
             ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute ReducedFineGridDw failed.", &
                  modulename, 0 )
             errstat = pge_errstat_error; RETURN
          ENDIF
       ENDIF
    ENDIF
    
    tempc = '.T.'; IF (.NOT. rm_mgline) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "RemoveMgLine", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute RemoveMgLine failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "NumberOfChannels", &
         HE5T_NATIVE_INT8, 1, INT(nCh, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute NumberOfChannels failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF    

    status = he5_ehwrglatt( l2_fileid, "StartWavelengths", &
         HE5T_NATIVE_FLOAT, nCh, REAL(winlim(1:nCh, 1), KIND=4) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute StartWavelengths failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF 

    status = he5_ehwrglatt( l2_fileid, "EndWavelengths", &
         HE5T_NATIVE_FLOAT, nCh, REAL(winlim(1:nCh, 2), KIND=4) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute EndWavelengths failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF 

    tempc = '.T.'; IF (.NOT. do_bandavg) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "SpectralCoadding", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute SpectralCoadding failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. wcal_bef_coadd) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "WaveCalibrationBefCoadd", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute WaveCalibrationBefCoaddfailed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    IF (do_bandavg) THEN
       status = he5_ehwrglatt( l2_fileid, "NumberOfSpectralAverage", &
            HE5T_NATIVE_INT8, nCh, INT(n_band_avg(1:nCh), KIND=1) )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, &
               "Write Global Attribute NumberOfSpectralAverage failed.", modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF

       status = he5_ehwrglatt( l2_fileid, "NumberOfSpectralSampling", &
            HE5T_NATIVE_INT8, nCh, INT(n_band_samp(1:nCh), KIND=1) )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, &
               "Write Global Attribute NumberOfSpectralSampling failed.", modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF
    ENDIF

    tempc = '.T.'; IF (.NOT. use_backup) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "UseMeanSoalrIrradiance", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute UseMeanSolarIrradiance failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "SlitTruncationLimit", &
         HE5T_NATIVE_FLOAT, 1, REAL(slit_trunc_limit, KIND=4) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute SlitTruncationLimit failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. yn_varyslit) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "VariableSlitWidth", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute VarySlitWidth failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. wavcal) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "PreWavelengthCalibration", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute PreWavelengthCalibration failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. wavcal_sol) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "SolarWavelengthCalibration", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute SolarWavelengthCalibration failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    IF (yn_varyslit) THEN
       tempc = '.T.'; IF (.NOT. smooth_slit) tempc = '.F.'
       count=LEN(TRIM(tempc))
       status = he5_ehwrglatt( l2_fileid, "SmoothSlit", HE5T_NATIVE_CHAR, count, tempc )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute SmoothSlit failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF

       status = he5_ehwrglatt( l2_fileid, "NumberSlitFitPoints", &
            HE5T_NATIVE_UINT8, 1, INT(slit_fit_pts, KIND=1) )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute NumberSlitFitPoints failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF

       status = he5_ehwrglatt( l2_fileid, "NumberSlitStepPoints", &
            HE5T_NATIVE_UINT8, 1, n_slit_step)
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute NumberSlitStepPoints failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF

       status = he5_ehwrglatt( l2_fileid, "NumberWavcalPoints", &
            HE5T_NATIVE_UINT8, 1, INT(wavcal_fit_pts, KIND=1) )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute NumberWavcalPoints failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF
       
       status = he5_ehwrglatt( l2_fileid, "NumberWavcalStepPoints", &
            HE5T_NATIVE_UINT8, 1, n_wavcal_step)
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute NumberWavcalStepPoints failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF
    ENDIF

    tempc1 = refspec_fname(solar_idx)
    count=LEN(TRIM(tempc1))
    status = he5_ehwrglatt( l2_fileid, "SolarReferenceFilename", HE5T_NATIVE_CHAR, count, tempc1 )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute SolarReferenceFilename failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "MaxNumberIterations", &
         HE5T_NATIVE_UINT8, 1, INT(max_itnum_rad, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute MaxNumberIterations failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "RadWavcalFreq", &
         HE5T_NATIVE_INT16, 1, INT(radwavcal_freq, KIND=2) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute RadWavcalFreq failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "MaxSZAProcessed", &
         HE5T_NATIVE_FLOAT, 1, REAL(szamax, KIND=4) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute MaxSZAProcessed failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "ViewGeometryAltitude", &
         HE5T_NATIVE_FLOAT, 1, REAL(zatmos, KIND=4) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute ViewGeometryAltitude failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. do_simu) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "DoRadianceSimulation", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute DoRadianceSimulation failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. do_simu) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "DoRadianceCalibration", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute DoRadianceCalibration failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    
    tempc = '.T.'; IF (.NOT. use_oe) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "UseOptimalEstimation", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute UseOptimalEstimation failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. biascorr) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "UseSoftCalibration", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute UseSoftCalibration failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    IF (biascorr) THEN
       tempc1 = biasfname
       count=LEN(TRIM(tempc1))
       status = he5_ehwrglatt( l2_fileid, "SoftCalibrationFilename", HE5T_NATIVE_CHAR, count, tempc1 )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute SoftCalibrationFilename failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF
    ENDIF

    tempc = '.T.'; IF (.NOT. degcorr) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "UseDegradationCorrection", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute UseDegradationCorrection failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    IF (degcorr) THEN
       tempc1 = degfname
       count=LEN(TRIM(tempc1))
       status = he5_ehwrglatt( l2_fileid, "DegradationCorrectionFilename", HE5T_NATIVE_CHAR, count, tempc1 )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute DegradationCorrectionFilename failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF
    ENDIF

    tempc = '.T.'; IF (.NOT. use_lograd) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "UseLogRad", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute UseLogRad failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. use_logstate) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "UseLogState", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute UseLogState failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. smooth_ozbc) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "SmoothOzoneBelowClouds", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute SmoothOzoneBelowClouds failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. use_large_so2_aperr) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "UseLargeSO2Aperr", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute UseLargeSO2Aperr failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. use_flns) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "UseFloorNoise", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute UseFloorNoise failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. ring_on_line) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "RingOnlineCalculation", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute RingOnlineCalculation failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. ring_convol) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "RingConvolution", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute RingConvolution failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. fit_atanring) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "FitAtanRing", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute FitAtanRing failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "WhichOzoneCLima", &
         HE5T_NATIVE_INT8, 1, INT(which_clima, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute WhichOzoneClima failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF 

    status = he5_ehwrglatt( l2_fileid, "WhichTOZ", &
         HE5T_NATIVE_INT8, 1, INT(which_toz, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute WhichTOZ failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF 

    tempc = '.T.'; IF (.NOT. norm_tropo3) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "NormalizeTropO3", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute NormalizeTropO3 failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "WhichOzoneAprioriError", &
         HE5T_NATIVE_INT8, 1, INT(which_aperr, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute WhichOzoneAprioriError failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. loose_aperr) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "LooseAprioriError", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute LooseAprioriError failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    IF (loose_aperr) THEN
       status = he5_ehwrglatt( l2_fileid, "MinStratError", &
            HE5T_NATIVE_FLOAT, 1, REAL(min_serr, KIND=4) )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute MinStratError failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF

       status = he5_ehwrglatt( l2_fileid, "MinTropoError", &
            HE5T_NATIVE_FLOAT, 1, REAL(min_terr, KIND=4) )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute MinTropoError failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "WhichAlbedo", &
         HE5T_NATIVE_INT8, 1, INT(which_alb, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute WhichAlbedo failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF 

    status = he5_ehwrglatt( l2_fileid, "WhichCloud", &
         HE5T_NATIVE_INT8, 1, INT(which_cld, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute WhichCloud failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF 

    tempc = '.T.'; IF (.NOT. aerosol) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "UseAerosol", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute UseAerosol failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    IF (aerosol) THEN
       tempc = '.T.'; IF (.NOT. strat_aerosol) tempc = '.F.'
       count=LEN(TRIM(tempc))
       status = he5_ehwrglatt( l2_fileid, "UseStratAerosol", HE5T_NATIVE_CHAR, count, tempc )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute UseStratAerosol failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF   

       status = he5_ehwrglatt( l2_fileid, "WhichAerosol", HE5T_NATIVE_UINT8, 1, INT(which_aerosol, KIND=1) )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute WhichAerosol failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF    

       tempc = '.T.'; IF (.NOT. scale_aod) tempc = '.F.'
       count=LEN(TRIM(tempc))
       status = he5_ehwrglatt( l2_fileid, "ScaleAOD", HE5T_NATIVE_CHAR, count, tempc )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute ScaleAOD failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF  
       
       IF (scale_aod) THEN
          status = he5_ehwrglatt( l2_fileid, "AerosolAOD", &
               HE5T_NATIVE_FLOAT, 1, REAL(scaled_aod, KIND=4) )
          IF( status == -1 ) THEN
             ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute AerosolAOD failed.", &
                  modulename, 0 )
             errstat = pge_errstat_error; RETURN
          ENDIF
       ENDIF
    ENDIF

    tempc = '.T.'; IF (.NOT. do_lambcld) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "DoLambertianCloud", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute DoLambertianCloud failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    IF (do_lambcld) THEN
       status = he5_ehwrglatt( l2_fileid, "LambertianCloudAlbedo", &
            HE5T_NATIVE_FLOAT, 1, REAL(lambcld_initalb, KIND=4) )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute LambertianCloudAlbedo failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF
    ELSE
       status = he5_ehwrglatt( l2_fileid, "ScatterCloudCOD", &
            HE5T_NATIVE_FLOAT, 1, REAL(scacld_initcod, KIND=4) )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute ScatterCloudCOD failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF
    ENDIF

    tempc = '.T.'; IF (.NOT. useasy) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "UseAsymmetryFactor", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute UseAsymmetryFactor failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "NumberOfPhaseMoments", HE5T_NATIVE_UINT8, 1, nmom )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute NumberOfPhaseMoments failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. use_effcrs) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "UseEffectiveCrossSection", &
         HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute UseEffectiveCrossSection failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    IF (.NOT. use_effcrs) THEN
       status = he5_ehwrglatt( l2_fileid, "HighResolutionSampRate", &
            HE5T_NATIVE_FLOAT, 1, REAL(hres_samprate, KIND=4) )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute HighResolutionSampRate failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF
       status = he5_ehwrglatt( l2_fileid, "NumberRadcalRegions", &
            HE5T_NATIVE_UINT8, 1, INT(radc_nsegsr, KIND=1) )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute NumberRadcalRegions failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF

       status = he5_ehwrglatt( l2_fileid, "RadcalRegions", &
            HE5T_NATIVE_FLOAT, radc_nsegsr, REAL(radc_lambnd(1:radc_nsegsr), KIND=4) )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute RadcalRegions failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF

       status = he5_ehwrglatt( l2_fileid, "RadcalSampRate", &
            HE5T_NATIVE_FLOAT, radc_nsegsr, REAL(radc_samprate(1:radc_nsegsr), KIND=4) )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute RadcalSampRate failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "WhichPolarizationCorr", &
         HE5T_NATIVE_UINT8, 1, INT(polcorr, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute WhichPolarizationCorr failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "NumberOfHalfStreams", &
         HE5T_NATIVE_UINT8, 1, INT(VlidortNstream, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute NumberOfHalfStreams failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc1 = ozabs_fname
    count=LEN(TRIM(tempc1))
    status = he5_ehwrglatt( l2_fileid, "OzoneXsectionFilename", HE5T_NATIVE_CHAR, count, tempc1 )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute OzoneXsectionFilename failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "RadcalNDiv", &
         HE5T_NATIVE_UINT8, 1, INT(ndiv, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute RadcalNDiv failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc = '.T.'; IF (.NOT. use_reg_presgrid) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "UseRegPresGrid", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute UseRegPresGrid failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    IF (.NOT. use_reg_presgrid) THEN
       tempc1 = presgrid_fname
       count=LEN(TRIM(tempc1))
       status = he5_ehwrglatt( l2_fileid, "PresGridFilename", HE5T_NATIVE_CHAR, count, tempc1 )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute PresGridFilename failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF
    ENDIF

    tempc = '.T.'; IF (.NOT. use_tropopause) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "UseTropopause", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute UseTropopause failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    IF (.NOT. use_tropopause) THEN
       tempc = '.T.'; IF (.NOT. fixed_ptrop) tempc = '.F.'
       count=LEN(TRIM(tempc))
       status = he5_ehwrglatt( l2_fileid, "FixedPresTropopause", HE5T_NATIVE_CHAR, count, tempc )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute FixedPresTropopause failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF
       
       IF (fixed_ptrop) THEN
          status = he5_ehwrglatt( l2_fileid, "FixedTropopause (mb)", HE5T_NATIVE_FLOAT, 1, REAL(pst0, KIND=4) )
          IF( status == -1 ) THEN
             ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute FixedTropopause failed.", &
                  modulename, 0 )
             errstat = pge_errstat_error; RETURN
          ENDIF
       ELSE
          status = he5_ehwrglatt( l2_fileid, "TheTropopauseLayer", HE5T_NATIVE_UINT8, 1, INT(ntp0, KIND=1) )
          IF( status == -1 ) THEN
             ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute THeTropopauseLayer failed.", &
                  modulename, 0 )
             errstat = pge_errstat_error; RETURN
          ENDIF
       ENDIF
    ENDIF

    tempc = '.T.'; IF (.NOT. adjust_trop_layer) tempc = '.F.'
    count=LEN(TRIM(tempc))
    status = he5_ehwrglatt( l2_fileid, "RedivideTroposphericLayers", HE5T_NATIVE_CHAR, count, tempc )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute RedivideTroposphericLayers failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
        
    status = he5_ehwrglatt( l2_fileid, "InitCldWavelength", HE5T_NATIVE_FLOAT, 1, REAL(pos_alb, KIND=4) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute InitCldWavelength failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_ehwrglatt( l2_fileid, "InitCldWaveFWHM", HE5T_NATIVE_FLOAT, 1, REAL(toms_fwhm, KIND=4) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute InitCldWaveFWHM failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    tempc1 = ozcrs_alb_fname
    count=LEN(TRIM(tempc1))
    status = he5_ehwrglatt( l2_fileid, "InitCldWaveXSecFilename", HE5T_NATIVE_CHAR, count, tempc1 )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Global Attribute InitCldWaveXSecFilename failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    !! attach the swath
    l2_swid = he5_swattach( l2_fileid, TRIM(ADJUSTL(l2_swathname)) )
    IF( l2_swid == -1 ) THEN
       WRITE( msg,'(A)' ) "he5_swattach:"// TRIM(l2_swathname) // &
            " failed in file " // TRIM(l2_filename )
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, msg, modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    
    ! Write Swath attributes   
    count  = 1
    status = he5_swwrattr( l2_swid, "nXtrack", HE5T_NATIVE_UINT8, count, INT(nXtrack, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute nXtrack failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_swwrattr( l2_swid, "nXtrackp1", HE5T_NATIVE_UINT8, count, INT(nXtrackp1, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute nXtrackp1 failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    status = he5_swwrattr( l2_swid, "nTimes", HE5T_NATIVE_UINT16, &
         count, INT(nTimes, KIND=2) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute nTimes failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    status = he5_swwrattr( l2_swid, "nTimesp1", HE5T_NATIVE_UINT16, count, INT(nTimesp1, KIND=2) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute nTimesp1 failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    status = he5_swwrattr( l2_swid, "nLayer", HE5T_NATIVE_UINT8, count, INT(nLayer, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute nLayer failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    status = he5_swwrattr( l2_swid, "nLayerp1", HE5T_NATIVE_UINT8, count, INT(nLayerp1, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute nLayerp1 failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

    status = he5_swwrattr( l2_swid, "nFitvar", HE5T_NATIVE_INT, count, INT(nFitvar) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute nFitvar failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    status = he5_swwrattr( l2_swid, "nOth", HE5T_NATIVE_INT, count, INT(nOth) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute nOth failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    status = he5_swwrattr( l2_swid, "mWavel", HE5T_NATIVE_UINT16, count, INT(mWavel, KIND=2) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute mWavel failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    status = he5_swwrattr( l2_swid, "nCh", HE5T_NATIVE_UINT8, count, INT(nCh, KIND=2))
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute nCh failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    IF (ngas > 0 .AND. (gaswrt .OR. ozwrtvar)) THEN
       status = he5_swwrattr( l2_swid, "nGas", HE5T_NATIVE_INT, count, INT(nGas) )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute nGas failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF
    ENDIF
    IF (aerosol) THEN
       status = he5_swwrattr( l2_swid, "nAer", HE5T_NATIVE_UINT8, count, INT(nAer, KIND=1) )
       IF( status == -1 ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute nAer failed.", &
               modulename, 0 )
          errstat = pge_errstat_error; RETURN
       ENDIF
    ENDIF
    status = he5_swwrattr( l2_swid, "mXtrack", HE5T_NATIVE_UINT8, count, INT(mXtrack, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute mXtrack failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    status = he5_swwrattr( l2_swid, "mTimes", HE5T_NATIVE_UINT16, count, INT(mTimes, KIND=2) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute mTimes failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    status = he5_swwrattr( l2_swid, "XOffset", HE5T_NATIVE_UINT8, count, INT(XOffset,KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute XOffset failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    status = he5_swwrattr( l2_swid, "YOffset", HE5T_NATIVE_UINT16, count, INT(YOffset, KIND=2) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute YOffset failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    status = he5_swwrattr( l2_swid, "nXbin", HE5T_NATIVE_UINT8, count, INT(nxbin, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute nXbin failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF
    status = he5_swwrattr( l2_swid, "nYbin", HE5T_NATIVE_UINT8, count, INT(nYbin, KIND=1) )
    IF( status == -1 ) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute nYbin failed.", &
            modulename, 0 )
       errstat = pge_errstat_error; RETURN
    ENDIF

!    !! Write variable names, gas names, other name and units
!    DO i = 1, nFitvar
!       WRITE(varnames(i), '(I2,A2,A4)') i, '. ', TRIM(fitvar_rad_str(mask_fitvar_rad(i)))
!       WRITE(varnames_nNum(i), '(A4)') TRIM(fitvar_rad_str(mask_fitvar_rad(i)))
!       WRITE(units(i), '(A20)' ) TRIM(fitvar_rad_unit(mask_fitvar_rad(i)))
!    ENDDO
!    varnamec = varnames(1)
!    DO i = 2, nFitvar
!       varnamec = TRIM(ADJUSTL(varnamec)) // ' ' // TRIM(ADJUSTL(varnames(i)))      
!    ENDDO
!    count=LEN(TRIM(varnamec))
!    status = he5_swwrattr( l2_swid, "varnames", HE5T_NATIVE_CHAR, count, TRIM(varnamec) )
!    IF( status == -1 ) THEN
!       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute varnames failed.", &
!            modulename, 0 )
!       errstat = pge_errstat_error; RETURN
!    ENDIF
!
!    IF (ngas > 0 .AND. (gaswrt .OR. ozwrtvar)) THEN
!       gasnamec = '1. ' // TRIM(ADJUSTL(varnames_nNum(fgasidxs(fgaspos(1)))))
!       DO i = 2, nGas
!       WRITE(tempc, '(I2,A2)') i, '. '
!          gasnamec = TRIM(ADJUSTL(gasnamec)) // ' ' // tempc // &
!               TRIM(ADJUSTL(varnames_nNum(fgasidxs(fgaspos(i)))))      
!       ENDDO
!       count=LEN(TRIM(gasnamec))
!       status = he5_swwrattr( l2_swid, "gasnames", HE5T_NATIVE_CHAR, count, TRIM(gasnamec) )
!       IF( status == -1 ) THEN
!          ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute varnames failed.", &
!               modulename, 0 )
!          errstat = pge_errstat_error; RETURN
!       ENDIF
!    ENDIF
!
!    othnamec = '1. ' // TRIM(ADJUSTL(varnames_nNum(fothvarpos(1))))
!    othunitc = '1. ' // TRIM(ADJUSTL(units(fothvarpos(1))))
!    DO i = 2, nOth
!       WRITE(tempc, '(I2,A2)') i, '. '
!       othnamec = TRIM(ADJUSTL(othnamec)) // ' ' // tempc // TRIM(ADJUSTL(varnames_nNum(fothvarpos(i)))) 
!       othunitc = TRIM(ADJUSTL(othunitc)) // ' ' // tempc // TRIM(ADJUSTL(units(fothvarpos(i))))       
!    ENDDO
!    count=LEN(TRIM(othnamec))
!    status = he5_swwrattr( l2_swid, "othnames", HE5T_NATIVE_CHAR, count, TRIM(othnamec) )
!    IF( status == -1 ) THEN
!       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute varnames failed.", &
!            modulename, 0 )
!       errstat = pge_errstat_error; RETURN
!    ENDIF
!    count=LEN(TRIM(othunitc))
!    status = he5_swwrattr( l2_swid, "othunits", HE5T_NATIVE_CHAR, count, TRIM(othunitc) )
!    IF( status == -1 ) THEN
!       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "Write Swath Attribute varnames failed.", &
!            modulename, 0 )
!       errstat = pge_errstat_error; RETURN
!    ENDIF
!    print *, TRIM(varnamec), TRIM(gasnamec), TRIM(othnamec), TRIM(othunitc)
    
    CALL He5_L2WrtClose (errstat)

    RETURN
  END SUBROUTINE He5_writeAttribute

  SUBROUTINE He5_L2WrtClose (pge_error_status )
    INTEGER, INTENT(OUT)   :: pge_error_status  
    INTEGER (KIND=i4)      :: ierr, status
    CHARACTER (LEN = 16)   :: modulename = 'he5_output_close'

    pge_error_status = pge_errstat_ok

    status = HE5_SWdetach( l2_swid )
    IF (status /= 0) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "HE5_SWdetach:"// TRIM(L2_filename) &
            // " failed.", modulename, 0 )
       pge_error_status = pge_errstat_error
    ENDIF
    
    status = HE5_SWclose ( l2_fileid )   
    IF (status /= 0) THEN
       ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "HE5_SWclose:"// TRIM(L2_filename) &
            // " failed.", modulename, 0 )
       pge_error_status = pge_errstat_error
    ENDIF

  END SUBROUTINE He5_L2WrtClose

END MODULE he5_output_module
