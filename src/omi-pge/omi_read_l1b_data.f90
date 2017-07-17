SUBROUTINE omi_read_radiance_paras (pge_error_status )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxchlen, vb_lev_omidebug
  USE OMSAO_variables_module,  ONLY: l1b_rad_filename, verb_thresh_lev, coadd_uv2
  USE OMSAO_omidata_module,    ONLY: omi_radiance_swathname, nxtrack_max, ntimes_max, &
       nxtrack, nfxtrack, ntimes, ncoadd, nswath, zoom_mode
  USE OMSAO_errstat_module
  USE hdfeos4_parameters
  USE L1B_Reader_class
  IMPLICIT NONE

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER, INTENT (OUT) :: pge_error_status

  ! ---------------
  ! Local variables
  ! ---------------
  TYPE (l1b_block_type) :: omi_data_block
  INTEGER (KIND=i4)     :: errstat, iline
  REAL (KIND=r4), DIMENSION (1:nxtrack_max) :: tmp_vza  

  ! Exteranl functions
  INTEGER               :: estat

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=23), PARAMETER :: modulename = 'omi_read_radiance_paras'

  ! --------------------------
  ! Initialize OUTPUT variable
  ! --------------------------
  pge_error_status = pge_errstat_ok
  errstat          = omi_s_success
  ntimes = 0;      nxtrack = 0
  zoom_mode = .FALSE.
  
  ! Determine if UV2 data are observed by zoom mode
  IF (nswath == 2 .OR. omi_radiance_swathname(1) == 'Earth UV-2 Swath')  THEN

     ! -----------------------------------------------------------------------
     ! Open data block called 'omi_data_block' with default size of 100 lines
     ! -----------------------------------------------------------------------
     errstat = L1Br_open ( omi_data_block, l1b_rad_filename, omi_radiance_swathname(nswath))
     IF ( errstat /= omi_s_success ) THEN
        estat = OMI_SMF_setmsg (omsao_e_open_l1b_file, 'L1Br_open failed.', modulename, 0); STOP 1
     END IF
     
     ! -------------------------------
     ! Get dimensions of current Swath
     ! -------------------------------
     errstat = L1Br_getSWdims ( omi_data_block, NumTimes_k=ntimes, nXtrack_k=nfxtrack)
     IF ( errstat /= omi_s_success ) THEN
        estat = OMI_SMF_setmsg (omsao_e_read_l1b_file, 'L1Br_getSWdims failed.', modulename, 0); STOP 1
     END IF

     iline = ntimes / 2
     errstat = L1Br_getGEOline ( omi_data_block, iline,    &
          ViewingZenithAngle_k    = tmp_vza (1:nfxtrack))

     IF (tmp_vza(nfxtrack) == tmp_vza(nfxtrack-1) .AND. tmp_vza(nfxtrack-1) &
          == tmp_vza(nfxtrack-2) ) zoom_mode = .TRUE.

     ! --------------------------
     ! Close data block structure
     ! --------------------------
     errstat = L1Br_CLOSE ( omi_data_block )
     IF ( errstat /= omi_s_success .AND. verb_thresh_lev >= vb_lev_omidebug ) THEN
        estat = OMI_SMF_setmsg ( omsao_w_clos_l1b_file, 'L1Br_CLOSE failed.', modulename, 0); STOP 1    
     END IF
  ENDIF

  IF (nswath ==2 .OR. omi_radiance_swathname(1) == 'Earth UV-1 Swath') THEN
     ! -----------------------------------------------------------------------
     ! Open data block called 'omi_data_block' with default size of 100 lines
     ! -----------------------------------------------------------------------
     errstat = L1Br_open ( omi_data_block, l1b_rad_filename, omi_radiance_swathname(1))
     IF ( errstat /= omi_s_success ) THEN
        estat = OMI_SMF_setmsg (omsao_e_open_l1b_file, 'L1Br_open failed.', modulename, 0); STOP 1
     END IF
     
     ! -------------------------------
     ! Get dimensions of current Swath
     ! -------------------------------
     errstat = L1Br_getSWdims ( omi_data_block, NumTimes_k=ntimes, nXtrack_k=nfxtrack)
     IF ( errstat /= omi_s_success ) THEN
        estat = OMI_SMF_setmsg (omsao_e_read_l1b_file, 'L1Br_getSWdims failed.', modulename, 0); STOP 1
     END IF
     
     ! --------------------------
     ! Close data block structure
     ! --------------------------
     errstat = L1Br_CLOSE ( omi_data_block )
     IF ( errstat /= omi_s_success .AND. verb_thresh_lev >= vb_lev_omidebug ) THEN
        estat = OMI_SMF_setmsg ( omsao_w_clos_l1b_file, 'L1Br_CLOSE failed.', modulename, 0); STOP 1    
     END IF
  ENDIF
  
  ! Note that nfxtrack means the number of pixels for UV-1 if it is selected
  nxtrack = nfxtrack * ncoadd

  IF (ntimes > ntimes_max) THEN
     pge_error_status = pge_errstat_error
     WRITE(*, '(A,I5)') 'Need to increase ntimes_max >= ', ntimes 
  ENDIF
  
  IF (nxtrack > nxtrack_max) THEN
     pge_error_status = pge_errstat_error
     WRITE(*, '(A)') 'Need to increase nxtrack_max!!!'
  ENDIF
  RETURN
END SUBROUTINE omi_read_radiance_paras

! Find the scan line/x track position based on inut lat/lon range
! Avoid including descending orbits
SUBROUTINE find_scan_line_range ( slat, elat, slon, elon, sline, eline, spix, epix, pge_error_status )

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: l1b_rad_filename, szamax
  USE OMSAO_omidata_module,   ONLY: nswath, omi_radiance_swathname, nxtrack_max, &
       ntimes_max, ntimes, nxtrack, nfxtrack
  USE L1B_Reader_class
  USE OMSAO_errstat_module
  IMPLICIT NONE

  ! -----------------------
  ! Input/Output variables
  ! -----------------------
  REAL (KIND=dp), INTENT(IN) :: slat, elat, slon, elon
  INTEGER,      INTENT (OUT) :: pge_error_status, sline, eline, spix, epix

  REAL (KIND=r4), DIMENSION (1:nxtrack_max, 0:ntimes_max) :: lons, lats, szas
  REAL (KIND=r4), DIMENSION (1:ntimes_max)                :: latdf
  LOGICAL, DIMENSION(1:nxtrack_max, 0:ntimes_max)         :: any_pixels
  REAL (KIND=dp)                                          :: tmplon, dlat
  INTEGER, DIMENSION (1:nxtrack_max, 2)                   :: okline
  INTEGER, DIMENSION (1:ntimes_max)                       :: lines
  INTEGER                                                 :: iline, errstat, ix

  ! Exteranl functions
  INTEGER                    :: estat

  TYPE (L1B_block_type)      :: omi_data_block   

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=20), PARAMETER :: modulename = 'find_scan_line_range'

  ! --------------------------
  ! Initialize OUTPUT variable
  ! --------------------------
  pge_error_status = pge_errstat_ok


  CALL omi_set_fitting_parameters ( pge_error_status )
  CALL omi_read_radiance_paras (pge_error_status )

  IF ( pge_error_status >= pge_errstat_error ) RETURN
   
  ! -----------------------------------------------------------
  ! Open data block (UV-1, if both are selected) 
  ! called 'omi_data_block' with nTimes lines
  ! -----------------------------------------------------------
  errstat = L1Br_open ( omi_data_block, l1b_rad_filename, omi_radiance_swathname(1))
  IF( errstat /= omi_s_success ) THEN
     estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, "L1Br_open failed.", modulename, 0 )
     pge_error_status = pge_errstat_error; RETURN
  END IF

  ! ---------------------------------
  ! Read all Latitudes and Longitudes
  ! ---------------------------------
  sline = -5; eline = -5; spix =-5; epix =-5
  DO iline = 0, ntimes - 1    
     errstat = L1Br_getGEOline ( omi_data_block, iline,  &
          Latitude_k              = lats (1:nfxtrack, iline), & 
          Longitude_k             = lons (1:nfxtrack, iline), &
          SolarZenithANgle_k      = szas (1:nfxtrack, iline))  
     
     IF( errstat /= omi_s_success ) THEN
        estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
             "L1Br_getGEOline failed.", modulename, 0 )
        pge_error_status = pge_errstat_error; RETURN
     END IF
  ENDDO

  lines(1:ntimes-1) = (/(iline, iline=1, ntimes-1)/)
  DO ix = 1, nfxtrack 
     latdf(1:ntimes-1) = lats(ix, 1:ntimes-1) - lats(ix, 0:ntimes-2)
     dlat  = latdf(ntimes / 2)
     
     okline(ix, 1) = MINVAL(MINLOC( lines(1:ntimes-1), MASK = (latdf(1:ntimes-1) > dlat/2.)))
     okline(ix, 2) = MINVAL(MAXLOC( lines(1:ntimes-1), MASK = (latdf(1:ntimes-1) > dlat/2.)))
     !print *, ix, okline(ix, 1), okline(ix, 2), ntimes, dlat
  ENDDO

  DO iline = 0, ntimes - 1      
     IF (ANY(lats(1:nfxtrack, iline) >= slat) .AND. ANY(lats(1:nfxtrack, iline) <= elat)) THEN
        any_pixels(1:nfxtrack, iline) = .FALSE.
        DO ix = 1, nfxtrack
           tmplon = lons(ix, iline)
           IF (elon > 180 .AND. tmplon < 0) tmplon = tmplon + 360.
           IF (tmplon >= slon .AND. tmplon <= elon .AND. szas(ix, iline) <= szamax &
                .AND. iline >= okline(ix, 1) .AND. iline <= okline(ix, 2)) any_pixels(ix, iline) = .TRUE.
          !IF (tmplon >= slon .AND. tmplon <= elon .AND. szas(ix, iline) <= szamax) any_pixels(ix, iline) = .TRUE.
        ENDDO
        IF ( ANY(any_pixels(1:nfxtrack, iline)) ) THEN
           IF (sline < 0) sline = iline
           !WRITE(*, '(I5,60L2)') iline, any_pixels(1:nfxtrack, iline)
           !WRITE(*, '(I5,60F6.1)') iline, szas(1:nfxtrack, iline)
           eline = iline
        ENDIF
     ENDIF
  ENDDO

  ! --------------------------
  ! Close data block structure
  ! --------------------------
  errstat = L1Br_close ( omi_data_block )
  IF( errstat /= omi_s_success ) THEN
     estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
          "L1Br_close failed.", modulename, 0 )
     pge_error_status = pge_errstat_error; RETURN
  END IF

  IF (sline >=0 .AND. eline >= sline) THEN
     DO ix = 1, nfxtrack
        IF (ANY(any_pixels(ix, sline:eline))) THEN
           IF (spix < 0) spix = ix
           epix = ix
        ENDIF
     ENDDO
  ENDIF

  !WRITE(*, '(4F8.2)') slat, elat, slon, elon
  !WRITE(*, '(4I8)')  sline, eline, spix, epix
  !STOP

  RETURN
END SUBROUTINE find_scan_line_range


SUBROUTINE omi_read_irradiance_data (lun, nxcoadd, first_pix, last_pix, pge_error_status ) 

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx, maxwin
  USE OMSAO_parameters_module, ONLY: maxchlen, maxwin, max_ring_pts, mrefl, vb_lev_omidebug
  USE OMSAO_variables_module,  ONLY: verb_thresh_lev, l1b_irrad_filename, wcal_bef_coadd, currpix, &
       numwin, coadd_uv2, band_selectors, winpix, winlim, scnwrt, l2_filename, refdbdir, use_backup, &
       reduce_resolution, redlam, redsampr, reduce_slit, rm_mgline, dwavmax, use_redfixwav, &
       which_slit
  USE OMSAO_omidata_module,    ONLY: orbnum, nswath, mswath, omi_nsolpix, nwavel_max,  &
       omi_irradiance_swathname, omi_irradiance_spec, omi_irradiance_qflg, omi_irradiance_prec, &
       omi_irradiance_wavl, omi_nwav_irrad, ntimes, nxtrack, nfxtrack, ncoadd, nxtrack_max, &
       omi_nsolring, omi_solspec_ring, solring_lin, solring_uin, omi_solspecr, omisol_winpix, &
       omisolr_winpix, omi_solnorm, omi_solpix_errstat, omi_solring_ndiv, irradwind, nxbin, &
       omiraddate, reduce_ubnd, reduce_lbnd, retlbnd, retubnd, orbnumsol, omichs, &
       omisol_version, zoom_p1, omi_irrad_stray, omi_rad_stray, omi_redslw
  USE ozprof_data_module,     ONLY: pos_alb, toms_fwhm, nrefl
  USE hdfeos4_parameters
  USE L1B_Reader_class
  USE OMSAO_errstat_module
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER, INTENT (IN)        :: nxcoadd, first_pix, last_pix, lun

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER , INTENT (OUT)      :: pge_error_status

  ! ---------------
  ! Local variables
  ! ---------------
  TYPE (L1B_block_type)                   :: omi_data_block
  INTEGER, PARAMETER                      :: nbits = 16
  INTEGER (KIND=i2), DIMENSION(0:nbits-1) :: mflgbits
  INTEGER (KIND=i2), DIMENSION(nxcoadd, nwavel_max, 0:nbits-1) :: flgbits
  INTEGER (KIND=i2), DIMENSION(nwavel_max)                     :: flgmsks
  INTEGER (KIND=i2)                     :: mflg
  INTEGER (KIND=i4), DIMENSION(mswath)  :: nwls
  INTEGER (KIND=i4)                     :: nx, nw, nwc, nt
  INTEGER (KIND=i4), DIMENSION(maxwin)  :: nwbin
  !REAL      (KIND=r8)                  :: tim
  !REAL      (KIND=r4)                  :: lat, lon, sazmin, sazmax, selmin, selmax, salt

  INTEGER   (KIND=i2), DIMENSION (mswath) :: spos, epos
  INTEGER  :: nwavel, is, ix, i, j, iix, nomi, fidx, lidx, ch, idum,  &
       iw, ic, idx, noff1, noff2, nring, irefl, nbin, theyear, themon, &
       theday, thedoy, np, npos, nsolbin
  INTEGER, DIMENSION (nwavel_max)                        :: idxs
  INTEGER (KIND=i2), DIMENSION (nwavel_max, nxtrack_max) :: tmpqflg
  INTEGER (KIND=i1)                                      :: tmpNinteg

  INTEGER (KIND=i4)                             :: errstat, iline
  REAL (KIND = dp)                              :: wcenter, normsc, tmpsampr, retswav, retewav, fdum
  REAL (KIND = dp), DIMENSION (maxwin, nxcoadd) :: wshis, wsqus                       
  REAL (KIND = dp), DIMENSION (nxcoadd * 2, sig_idx, nwavel_max) :: omispec
  REAL (KIND = dp), DIMENSION (nxcoadd * 2, 2, nwavel_max)       :: strayspec 
  REAL (KIND = dp), DIMENSION (sig_idx, nwavel_max, nxtrack_max) :: tmpspec
  REAL (KIND = dp), DIMENSION (spc_idx, max_ring_pts)            :: omirsol 
  LOGICAL                                 :: error

  CHARACTER (LEN=maxchlen)                :: bkfname, straylight_fname
  CHARACTER (LEN=3)                       :: opfc

  ! Exteranl functions
  INTEGER                                 :: estat
    
  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=24), PARAMETER :: modulename = 'omi_read_irradiance_data' 

  !nt = 0; nx = 0 ; nw = 0 ; nwc = 0 
  j = 1; nwavel = 0; iline = 0
  pge_error_status = pge_errstat_ok; errstat = omi_s_success
  omi_nsolpix = 0; omi_nwav_irrad = 0

  ! For zoom-in global products (once every 32 days), solar irradiance (not in radiance)
  ! in UV1 are provided at 60 Xtrack positions. Under this condition, solar irradiance needs
  ! to be rebinned additionally to the normal.
  ! In UV-2, also provided at 60 across-track position, but corresponding to positions 16-45
  ! in normal mode after coadd every two spectra
  nsolbin = 1  

  ! ----------------------------
  ! Initialize irradiance arrays
  ! ----------------------------  
  omi_irradiance_spec (1:nwavel_max,1:nxtrack) = 0.0
  omi_irradiance_prec (1:nwavel_max,1:nxtrack) = 0.0
  omi_irradiance_qflg (1:nwavel_max,1:nxtrack) = 0
  omi_irradiance_wavl (1:nwavel_max,1:nxtrack) = 0.0
  omi_solpix_errstat(1:nxtrack) = pge_errstat_ok
 
  IF (.NOT. use_backup) THEN 
      DO is = 1, nswath
         ch = omichs(is)
         
         ! ------------------------------------------------------
         ! Open data block structure with default size of 1 lines
         ! ------------------------------------------------------
         errstat = L1Br_open ( omi_data_block, l1b_irrad_filename, &
              TRIM(ADJUSTL(omi_irradiance_swathname(is))) )
         !errstat = L1Br_open ( omi_data_block, l1b_irrad_filename, 'Sun Volume VIS Swath ' )
         IF( errstat /= omi_s_success ) THEN
            estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, "L1Br_open failed.", modulename, 0); STOP 1
         END IF
         ! ----------------------------------
         ! Obtain irradiance swath dimensions
         ! ----------------------------------
         errstat = L1Br_getSWdims ( omi_data_block, NumTimes_k=nt, nXtrack_k=nx, nWavel_k=nw, nWavelCoef_k=nwc)
         IF( errstat /=  omi_s_success ) THEN
            estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, "IL1Br_getSWdims failed.", modulename, 0); STOP 1
         END IF

         ! For zoom-in global products, solar irradiance needs to be binned by a factor of 2 additionally
         IF (ch == 1 .AND. nx == nfxtrack * 2) nsolbin = 2
            
         ! ----------------------------------------------------------
         ! Obtain time, geolocation, and angular information on block
         ! ----------------------------------------------------------
         !errstat = L1Br_getGEOline ( omi_data_block, iline, Time_k= tim, SpacecraftLatitude_k = lat, &
         !     SpacecraftLongitude_k = lon, SolarElevationMinimum_k = selmin, &
         !     SolarElevationMaximum_k = selmax, SolarAzimuthMinimum_k = sazmin, &
         !     SolarAzimuthMaximum_k = sazmax,  SpacecraftAltitude_k = salt)
         !IF( errstat /= omi_s_success ) THEN
         !   estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, "L1Br_getGEOline failed.", modulename, 0); STOP 1
         !END IF
    
         errstat = L1Br_getDATA ( omi_data_block, iline, MeasurementQualityFlags_k = mflg)
         IF( errstat /= omi_s_success ) THEN
            estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, "L1Br_getDATA failed.", modulename, 0); STOP 1
         END IF
    
         CALL convert_2bytes_to_16bits ( nbits, 1, mflg, mflgbits(0:nbits-1))

         IF (mflgbits(0) == 1 .OR. mflgbits(1) == 1 .OR. mflgbits(3) == 1 .OR. mflgbits(12) == 1) THEN
            WRITE(*, *) 'All irradiances could not be used, use backup irradiances: ', &
                 TRIM(ADJUSTL(omi_irradiance_swathname(is)))
            omi_solpix_errstat(1:nxtrack) = pge_errstat_error
            pge_error_status = pge_errstat_error; RETURN
         ELSE IF (ANY(mflgbits == 1)) THEN
            WRITE(*, *) 'Warning set on all irradiances: ', TRIM(ADJUSTL(omi_irradiance_swathname(is)))
            omi_solpix_errstat(1:nxtrack) = pge_errstat_error  
            pge_error_status = pge_errstat_error; RETURN
         ENDIF
        
         errstat = L1Br_getSIGline ( omi_data_block, iline,     &
              Signal_k            = omi_irradiance_spec(j:, :), &
              SignalPrecision_k   = omi_irradiance_prec(j:, :), &
              PixelQualityFlags_k = omi_irradiance_qflg(j:, :), &
              Wavelength_k        = omi_irradiance_wavl(j:, :), &
              NumberSmallPixelColumns_k = tmpNinteg,            &
              Nwl_k  = nwls(ch) )
         IF( errstat /= omi_s_success ) THEN
            estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, "L1Br_getSIGline failed.", modulename, 0); STOP 1
         END IF

         nwavel = nwavel + nwls(ch); spos(ch) = j; j = nwavel + 1; epos(ch) = nwavel
     
         ! --------------------------
         ! Close data block structure
         ! --------------------------
         errstat = L1Br_CLOSE ( omi_data_block )
         IF( errstat /= omi_s_success .AND. verb_thresh_lev >= vb_lev_omidebug ) THEN
            estat = OMI_SMF_setmsg ( omsao_w_clos_l1b_file, "L1Br_CLOSE failed.", modulename, 0); STOP 1
         END IF 
         ! Need to sort the data in increasing wavelength
         IF (omi_irradiance_wavl(spos(ch), 1) > omi_irradiance_wavl(epos(ch), 1)) THEN
            idxs(spos(ch):epos(ch)) = (/ (i, i = epos(ch), spos(ch), -1) /)
            omi_irradiance_wavl(spos(ch):epos(ch), :) = omi_irradiance_wavl(idxs(spos(ch):epos(ch)), :)
            omi_irradiance_spec(spos(ch):epos(ch), :) = omi_irradiance_spec(idxs(spos(ch):epos(ch)), :)
            omi_irradiance_prec(spos(ch):epos(ch), :) = omi_irradiance_prec(idxs(spos(ch):epos(ch)), :)
            omi_irradiance_qflg(spos(ch):epos(ch), :) = omi_irradiance_qflg(idxs(spos(ch):epos(ch)), :)     
         ENDIF

         !IF ( ch == 2 ) CALL corruv2wav(nwls(ch), nx, omi_irradiance_wavl(spos(ch):epos(ch), 1:nx)) 
       
         !OPEN(unit=90, FILE='/data/dumbo/xliu/OMIHCLD/OMIL1BBIRR-o05168_vis.dat', STATUS='old')
         !WRITE(90, *) nx, nw
         !DO ix = 1, nx 
         !   WRITE(90, *) ix
         !   DO iw = spos(ch), epos(ch)
         !      !CALL convert_2bytes_to_16bits ( nbits, 1, omi_irradiance_qflg(iw, ix), mflgbits(0:nbits-1))
         !      WRITE(90, '(F10.4,D14.6,1X)') omi_irradiance_wavl(iw, ix), omi_irradiance_spec(iw, ix) !, &
         !      !mflgbits(0:nbits-1)
         !   ENDDO
         !ENDDO
      ENDDO
  ELSE
     ! Determine sun-earth distance correction
     READ(omiraddate, '(I4,1X,2I2)') theyear, themon, theday
     CALL GET_DOY(theyear, themon,  theday, thedoy)
     IF (thedoy == 366) thedoy = 365

     OPEN (UNIT=lun, FILE= ADJUSTL(TRIM(refdbdir)) // 'solar-distance.dat', &
          STATUS='UNKNOWN', IOSTAT=errstat)
     IF ( errstat /= pge_errstat_ok ) THEN
        WRITE(*, '(2A)') modulename, ': Cannot open Sun-Earth Distance datafile!!!'
        pge_error_status = pge_errstat_error; RETURN
     END IF
     DO i = 1, 12
        READ(LUN, *)
     ENDDO
     DO i = 1, thedoy
        READ(LUN, *) normsc, normsc
     ENDDO
     CLOSE(LUN)
     
     normsc = 1.0 / normsc ** 2  ! solar energy is inversely proportional to square distance
    
     ! Determine backup filename
     IF (omisol_version == 2) THEN
        IF (orbnum >= 6551) THEN
           opfc = 'aft'
        ELSE
           opfc = 'bef'
        ENDIF
        bkfname = ADJUSTL(TRIM(refdbdir)) // 'OMI/omisol_'  // opfc // '_6551_avg_backup.dat'
     ELSE
        bkfname = ADJUSTL(TRIM(refdbdir)) // 'OMI/omisol_v003_avg_nshi_backup.dat'
     ENDIF
     OPEN (UNIT=lun, FILE=TRIM(ADJUSTL(bkfname)), STATUS='UNKNOWN', IOSTAT=errstat)
     IF ( errstat /= pge_errstat_ok ) THEN
        WRITE(*, '(2A)') modulename, ': Cannot open solar backup file!!!'
        pge_error_status = pge_errstat_error; RETURN
     END IF

    
     ! Open the file and read all the data
     nwavel = 0
     DO is = 1, 2
        READ(lun, *) nx, nwls(is)
        spos(is) = nwavel + 1; epos(is) = nwavel + nwls(is)
        DO i = 1, nx
           READ(lun, *) 
           DO j = 1, nwls(is)
              READ(lun, *) omi_irradiance_wavl(nwavel + j, i), omi_irradiance_spec(nwavel + j, i), &
                   omi_irradiance_prec(nwavel + j, i), idum, idum
              IF (idum > 0) omi_irradiance_prec(nwavel + j, i) = omi_irradiance_prec(nwavel + j, i) &
                   / SQRT( REAL(idum, KIND=dp) )
           ENDDO
        ENDDO
        omi_irradiance_spec(spos(is):epos(is), 1:nx) = omi_irradiance_spec(spos(is):epos(is), 1:nx) * normsc
        omi_irradiance_prec(spos(is):epos(is), 1:nx) = omi_irradiance_prec(spos(is):epos(is), 1:nx) * normsc
        nwavel = epos(is)
        !IF ( is == 2 ) CALL corruv2wav(nwls(is), nx, omi_irradiance_wavl(spos(is):epos(is), 1:nx)) 
     ENDDO

     IF (nswath == 1) THEN
        IF (band_selectors(1) == 1) THEN
           nwavel = nwls(1)         
        ELSE
           nwavel = nwls(2)
           omi_irradiance_wavl(1:nwavel, :) = omi_irradiance_wavl(spos(2):epos(2), :)
           omi_irradiance_spec(1:nwavel, :) = omi_irradiance_spec(spos(2):epos(2), :) 
           omi_irradiance_prec(1:nwavel, :) = omi_irradiance_prec(spos(2):epos(2), :)
           spos(2) = spos(2) - nwls(1); epos(2) = epos(2) - nwls(1) 
        ENDIF
     ENDIF
     CLOSE(LUN)
  ENDIF
          
!  ! xliu: Feb/19/2008, read straylight spectra
!  straylight_fname = ADJUSTL(TRIM(refdbdir)) // 'OMI/omi_irrad_sl_v3.dat' 
!  OPEN (UNIT=lun, FILE=TRIM(ADJUSTL(straylight_fname)), STATUS='UNKNOWN', IOSTAT=errstat)
!  IF ( errstat /= pge_errstat_ok ) THEN
!     WRITE(*, '(2A)') modulename, ': Cannot open irradiance straylight file!!!'
!     pge_error_status = pge_errstat_error; RETURN
!  END IF
!  
!  ! Open the file and read all the data
!  nwavel = 0
!  DO is = 1, 2
!     READ(lun, *) nx, nwls(is)
!     spos(is) = nwavel + 1; epos(is) = nwavel + nwls(is)
!     DO i = 1, nx
!        READ(lun, *) 
!        DO j = 1, nwls(is)
!           READ(lun, *) fdum, omi_irrad_stray(nwavel + j, i)
!        ENDDO
!     ENDDO
!     nwavel = epos(is)
!  ENDDO
!
!  IF (nswath == 1) THEN
!     IF (band_selectors(1) == 1) THEN
!        nwavel = nwls(1)         
!     ELSE
!        nwavel = nwls(2)
!        omi_irrad_stray(1:nwavel, :) = omi_irrad_stray(spos(2):epos(2), :) 
!        spos(2) = spos(2) - nwls(1); epos(2) = epos(2) - nwls(1) 
!     ENDIF
!  ENDIF
!  CLOSE(LUN)
!
!  straylight_fname = ADJUSTL(TRIM(refdbdir)) // 'OMI/omi_rad_sl_v3.dat' 
!  OPEN (UNIT=lun, FILE=TRIM(ADJUSTL(straylight_fname)), STATUS='UNKNOWN', IOSTAT=errstat)
!  IF ( errstat /= pge_errstat_ok ) THEN
!     WRITE(*, '(2A)') modulename, ': Cannot open radiance straylight file!!!'
!     pge_error_status = pge_errstat_error; RETURN
!  END IF
!  
!  ! Open the file and read all the data
!  nwavel = 0
!  DO is = 1, 2
!     READ(lun, *) nx, nwls(is)
!     spos(is) = nwavel + 1; epos(is) = nwavel + nwls(is)
!     DO i = 1, nx
!        READ(lun, *) 
!        DO j = 1, nwls(is)
!           READ(lun, *) fdum, omi_rad_stray(nwavel + j, i)
!        ENDDO
!     ENDDO
!     nwavel = epos(is)
!  ENDDO
!
!  IF (nswath == 1) THEN
!     IF (band_selectors(1) == 1) THEN
!        nwavel = nwls(1)         
!     ELSE
!        nwavel = nwls(2)
!        omi_irrad_stray(1:nwavel, :) = omi_rad_stray(spos(2):epos(2), :) 
!        spos(2) = spos(2) - nwls(1); epos(2) = epos(2) - nwls(1) 
!     ENDIF
!  ENDIF
!  CLOSE(LUN) 

  ! Do not coadd wavelengths with a gap (e.g., filter Mg absorption lines), need to determine
  ! delta-lamda in UV-1
  ! Note in OMI delta-lamda varies with wavelength (largest for the first two pixels in each channel)
  dwavmax = ( omi_irradiance_wavl(2, 1) - omi_irradiance_wavl(1, 1) ) * 1.1
  
  ! Degrade spectral resolution if necessary
  IF (reduce_resolution) THEN
     nwavel = 0; j = 1
     DO is = 1, nswath
        ch = omichs(is)
        IF (coadd_uv2 .AND. is == 1) THEN
           npos = nfxtrack * nsolbin
        ELSE
           npos = nxtrack
        ENDIF
        !IF (nswath == 2) THEN
        !   retswav = retlbnd(ch); retewav = retubnd(ch)
        !ELSE
        !   retswav = retlbnd(band_selectors(1)); retewav = retubnd(band_selectors(1))
        !ENDIF
        retswav = retlbnd(ch); retewav = retubnd(ch)

        tmpsampr = redsampr ; IF (is == 1 .AND. band_selectors(1) == 1) tmpsampr = redsampr / 3.0
        np = nwls(ch)

        tmpspec(wvl_idx, 1:np, 1:npos) = omi_irradiance_wavl(spos(ch):epos(ch), 1:npos)
        tmpspec(spc_idx, 1:np, 1:npos) = omi_irradiance_spec(spos(ch):epos(ch), 1:npos)
        tmpspec(sig_idx, 1:np, 1:npos) = omi_irradiance_prec(spos(ch):epos(ch), 1:npos)
        tmpqflg(   1:np, 1:npos) = omi_irradiance_qflg(spos(ch):epos(ch), 1:npos)       

        DO ix = 1, npos
           CALL convert_2bytes_to_16bits ( nbits, np, tmpqflg(1:np, ix), flgbits(1, 1:np, 0:nbits-1))
           tmpqflg(1:np, ix) = flgbits(1, 1:np, 0) &   ! Missing
                             + flgbits(1, 1:np, 1) &   ! Bad
                             + flgbits(1, 1:np, 2)     ! Processing error
        ENDDO
        
        CALL reduce_irrad_resolution (tmpspec(:, 1:np, 1:npos), tmpqflg(1:np, 1:npos), np, npos, &
            reduce_slit, omi_redslw(is), tmpsampr, redlam, retswav, retewav, reduce_lbnd(ch), &
             reduce_ubnd(ch), nwls(ch), pge_error_status)

        IF (pge_error_status == pge_errstat_error) RETURN

        nwavel = nwavel + nwls(ch); spos(ch) = j; j = nwavel + 1; epos(ch) = nwavel
        omi_irradiance_wavl(spos(ch):epos(ch), 1:npos) = tmpspec(wvl_idx, 1:nwls(ch), 1:npos)
        omi_irradiance_spec(spos(ch):epos(ch), 1:npos) = tmpspec(spc_idx, 1:nwls(ch), 1:npos)
        omi_irradiance_prec(spos(ch):epos(ch), 1:npos) = tmpspec(sig_idx, 1:nwls(ch), 1:npos)  
        omi_irradiance_qflg(spos(ch):epos(ch), 1:npos) = 0   ! All data are good  (pre filtered) 
     ENDDO
  ENDIF

  IF (nwavel > nwavel_max) THEN
     WRITE(*, '(A)') "Need to increase nwavel_max!!!"
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

  ! Determine number of wavelengths to be read for deteriming cloud fraction
  fidx = MAXVAL ( MINLOC ( omi_irradiance_wavl(1:nwavel, 1), MASK = &
       (omi_irradiance_wavl(1:nwavel, 1) > pos_alb - toms_fwhm * 1.4) ))
  lidx = MAXVAL ( MAXLOC ( omi_irradiance_wavl(1:nwavel, 1), MASK = &
       (omi_irradiance_wavl(1:nwavel, 1) < pos_alb + toms_fwhm * 1.4) ))
  IF (fidx <1 .OR. lidx > nwavel) THEN
     WRITE(*, '(2A)') modulename, ': Need to change pos_alb/toms_fwhm!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  nrefl = lidx - fidx + 1 
  IF (nrefl > mrefl ) THEN
     WRITE(*, '(2A)') modulename, ': Need to increase mrefl!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
            
  ! Determine number of binning for different fitting windows
  DO iw = 1, numwin
     ch = band_selectors(iw)
     IF (ch == 1 .OR. .NOT. coadd_uv2) THEN
        nwbin(iw) = nxbin 
        IF (ch == 1) nwbin(iw) = nxbin 
     ELSE
        nwbin(iw) = nxbin * ncoadd 
     ENDIF
     nwbin(iw) = nwbin(iw) * nsolbin

  ENDDO

  ! Subset and coadd irradiance spectrum
  DO ix = first_pix, last_pix
     currpix = ix          ! Global variable, to be used in wavlength calibraiton with omi slit before coadding
      
     ! Indices for UV-1 are from 1 to nfxtrack
     ! Indices for UV-2 are from 1 to nxtrack (nfxtrack * 2)
     ! UV2 indices corresponding to UV-1 pixel ix are ix * 2 -1 & ix * 2 (or iix+1, iix+2), respectively
     ! If additional across track coadding (e.g., nxbin) is performed, then for a particular ix
     ! UV1: nbin = nxbin; UV-2: nbin = nxbin * 2
     ! Coadded original across track pixels are: iix + 1 : iix + nbin (iix = (ix-1) * nbin)

     ! Get quality flag bits, coadd flags if necessary to avoid coadding inconsistent # of pixels 
     flgmsks = 0
     DO is = 1, nswath
        ch = omichs(is)

        !Do not use nwbin
        IF (is == 1) THEN
           nbin = nxbin 
        ELSE
           nbin = nxbin * ncoadd
        ENDIF
        nbin = nbin * nsolbin  ! Zoom mode UV1 solar irradiance got 60 positions
        iix = (ix - 1) * nbin

        ! Shift the position by 15 
        IF (ch == 2 .AND. nsolbin == 2) THEN
           iix = iix - (zoom_p1 - 1) * nsolbin
        ENDIF
        
        IF (.NOT. reduce_resolution) THEN
           ! properly align cross track positions to be coadded (should be within one pixel)
       
           IF (nbin / nsolbin > 2) CALL prespec_align(nwls(ch), nbin, &
                omi_irradiance_wavl(spos(ch):epos(ch), iix+1:iix+nbin), &
                omi_irradiance_spec(spos(ch):epos(ch), iix+1:iix+nbin), &
!                omi_irrad_stray(spos(ch):epos(ch), iix+1:iix+nbin), &
!                omi_rad_stray(spos(ch):epos(ch), iix+1:iix+nbin), &
                omi_irradiance_prec(spos(ch):epos(ch), iix+1:iix+nbin), &       
                omi_irradiance_qflg(spos(ch):epos(ch), iix+1:iix+nbin))

           DO ic = 1, nbin
              
              CALL convert_2bytes_to_16bits ( nbits, nwls(ch), omi_irradiance_qflg(spos(ch):epos(ch), iix + ic ), &
                   flgbits(ic, spos(ch):epos(ch), 0:nbits-1))
              flgmsks(spos(ch):epos(ch)) = flgmsks(spos(ch):epos(ch)) &
                   + flgbits(ic, spos(ch):epos(ch), 0)                &   ! Missing
                   + flgbits(ic, spos(ch):epos(ch), 1)                &   ! Bad 
                   + flgbits(ic, spos(ch):epos(ch), 2)                    ! Processing error
           ENDDO

           !DO i = spos(is), epos(is)
           !   WRITE(90, '(F10.4, D14.6, 16I2)') omi_irradiance_wavl(i, iix+1), &
           !        omi_irradiance_spec(i, iix+1), flgbits(1, i, 0:nbits-1)
           !ENDDO

       ELSE
          ! Already aligned because of using common wavelength scale
          DO ic = 1, nbin
             flgmsks(spos(ch):epos(ch)) = flgmsks(spos(ch):epos(ch)) + &
                  omi_irradiance_qflg(spos(ch):epos(ch), iix + ic )
          ENDDO
       ENDIF       
    ENDDO
   
     ! Subset valid data
     nomi = 0; omispec = 0.0
     ! strayspec = 0.0
     DO iw = 1, numwin
        ch = band_selectors(iw)
        nbin = nwbin(iw);  iix = (ix - 1) * nbin

        ! Shift the position by 15 
        IF (ch == 2 .AND. nsolbin == 2) THEN
           iix = iix - (zoom_p1 - 1) * nsolbin
        ENDIF
                        
        winpix(iw, 1) = MINVAL ( MINLOC ( omi_irradiance_wavl(spos(ch):epos(ch),iix + 1), &
             MASK = omi_irradiance_wavl(spos(ch):epos(ch),iix + 1) >= winlim(iw, 1)) )
        winpix(iw, 2) = MAXVAL ( MAXLOC ( omi_irradiance_wavl(spos(ch):epos(ch),iix + 1), &
             MASK = omi_irradiance_wavl(spos(ch):epos(ch),iix + 1) <= winlim(iw, 2)) )

        omi_nsolpix  (iw, ix) = nomi

        fidx = winpix(iw, 1) + spos(ch) - 1
        lidx = winpix(iw, 2) + spos(ch) - 1
        omisol_winpix(iw, ix, 1:2) = 0       

        IF (rm_mgline .AND. winlim(iw, 1) < 286.0 .AND. winlim(iw, 2) > 286.0) THEN
           DO i = fidx, lidx
              IF (ALL(omi_irradiance_spec(i, iix+1:iix+nbin)   > 0.0)    .AND. &
                   ALL(omi_irradiance_spec(i, iix+1:iix+nbin)  < 4.0E14) .AND. flgmsks(i) == 0 .AND. &
                   !(ALL(omi_irradiance_wavl(i, iix+1:iix+nbin) < 273.8)   .OR. &
                   !ALL(omi_irradiance_wavl(i, iix+1:iix+nbin)  > 275.2))  .AND. &
                   (ALL(omi_irradiance_wavl(i, iix+1:iix+nbin) < 278.8)   .OR. &
                   ALL(omi_irradiance_wavl(i, iix+1:iix+nbin)  > 281.0))  .AND. &
                   (ALL(omi_irradiance_wavl(i, iix+1:iix+nbin) < 284.7)   .OR. &
                   ALL(omi_irradiance_wavl(i, iix+1:iix+nbin)  > 285.7))) THEN
                 !(ALL(omi_irradiance_wavl(i, iix+1:iix+nbin) < 278.0)   .OR. &
                 !ALL(omi_irradiance_wavl(i, iix+1:iix+nbin)  > 282.0))  .AND. &
                 !(ALL(omi_irradiance_wavl(i, iix+1:iix+nbin) < 284.0)   .OR. &
                 !ALL(omi_irradiance_wavl(i, iix+1:iix+nbin)  > 286.0))) THEN
                 nomi = nomi + 1
                 omispec(1:nbin, wvl_idx, nomi) = omi_irradiance_wavl(i, iix+1:iix+nbin)
                 omispec(1:nbin, spc_idx, nomi) = omi_irradiance_spec(i, iix+1:iix+nbin)
                 omispec(1:nbin, sig_idx, nomi) = omi_irradiance_prec(i, iix+1:iix+nbin)
!                 strayspec(1:nbin, 1, nomi)     = omi_irrad_stray(i, iix+1:iix+nbin)
!                 strayspec(1:nbin, 2, nomi)     = omi_rad_stray(i, iix+1:iix+nbin)
                 IF (omisol_winpix(iw, ix, 1) == 0) omisol_winpix(iw, ix, 1) = i
                 omisol_winpix(iw, ix, 2) = i
                 irradwind(nomi, ix) = i
              ENDIF
           ENDDO
        ELSE
           DO i = fidx, lidx
              IF (ALL(omi_irradiance_spec(i, iix+1:iix+nbin)  > 0.0)   .AND. &
                  ALL(omi_irradiance_spec(i, iix+1:iix+nbin)  < 4.0E14) .AND. flgmsks(i) == 0 ) THEN
                 nomi = nomi + 1
                 omispec(1:nbin, wvl_idx, nomi) = omi_irradiance_wavl(i, iix+1:iix+nbin)
                 omispec(1:nbin, spc_idx, nomi) = omi_irradiance_spec(i, iix+1:iix+nbin)
                 omispec(1:nbin, sig_idx, nomi) = omi_irradiance_prec(i, iix+1:iix+nbin)
!                 strayspec(1:nbin, 1, nomi)     = omi_irrad_stray(i, iix+1:iix+nbin)
!                 strayspec(1:nbin, 2, nomi)     = omi_rad_stray(i, iix+1:iix+nbin)
                 IF (omisol_winpix(iw, ix, 1) == 0) omisol_winpix(iw, ix, 1) = i
                 omisol_winpix(iw, ix, 2) = i
                 irradwind(nomi, ix) = i
              ENDIF
           ENDDO
        ENDIF
           
        omi_nsolpix(iw, ix) = nomi - omi_nsolpix(iw, ix)
    
        !WRITE(*, '(2I5, 2F8.3, 2I5, 2F8.3, 6I5)') ix, iw, omi_irradiance_wavl(spos(ch), iix+1), &
        !     omi_irradiance_wavl(epos(ch), iix+1), spos(ch), epos(ch), winlim(iw, 1), winlim(iw, 2), &
        !     winpix(iw, 1), winpix(iw, 2), fidx, lidx, lidx - fidx + 1, omi_nsolpix(iw, ix)
     ENDDO

     omi_nwav_irrad(ix) = nomi

     ! Perform coadding when necessary
     fidx = 1
     DO iw = 1, numwin       
        ch = band_selectors(iw);  nbin = nwbin(iw)
        lidx = fidx + omi_nsolpix(iw, ix) - 1 
        IF (nbin > 1) THEN
!           CALL stray_coadd(omi_nsolpix(iw, ix), nbin, strayspec(1:nbin, 1:2, fidx:lidx))
           CALL solwavcal_coadd(wcal_bef_coadd, omi_nsolpix(iw, ix), nbin, &
                omispec(1:nbin, :, fidx:lidx), wshis(iw, 1:nbin), wsqus(iw, 1:nbin), error)
           IF (error) THEN
              WRITE(*, '(A)') 'No solar wavelength calibration before coadding!!!'
              omi_solpix_errstat(ix) = pge_errstat_warning
           ENDIF
        ENDIF
        fidx = lidx + 1    
     ENDDO
                 
     ! Subset solar spectrum for Ring effect
     IF (.NOT. reduce_resolution .OR. (reduce_resolution .AND. .NOT. use_redfixwav)) THEN
        omirsol = 0.0
        ch = band_selectors(1)
    
        IF (band_selectors(1) == 1) THEN
           noff1 = 12
        ELSE
           noff1 = 25
        ENDIF
       
        nring = omi_nsolpix(1, ix) + noff1
        omirsol(1:spc_idx, noff1+1 : nring) = omispec(1, 1:spc_idx, 1:omi_nsolpix(1, ix)) 
        omi_solring_ndiv(ix) = 0
        
        ! add extra spectra before first window (uncoadded)
        ! if unavailable, needed to ammened with solar reference spectrum
        noff1 = noff1 + 1 ;     nbin = nwbin(1); iix = (ix - 1) * nbin
        
        IF (ch == 2 .AND. nsolbin == 2) THEN  ! Shift the position by 15 
           iix = iix - (zoom_p1 - 1) * nsolbin
        ENDIF
        
        DO i = winpix(1, 1) - 1, 1, -1
           IF (ALL(omi_irradiance_spec(i, iix+1:iix+nbin) > 0.0) .AND. &
                ALL(omi_irradiance_spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0) THEN
              noff1 = noff1 - 1
              omirsol(wvl_idx, noff1) = SUM(omi_irradiance_wavl(i, iix+1:iix+nbin)) / nbin
              omirsol(spc_idx, noff1) = SUM(omi_irradiance_spec(i, iix+1:iix+nbin)) / nbin
              IF (noff1 == 1) EXIT
           ENDIF
        ENDDO
       
        DO iw = 2, numwin
           ch = band_selectors(iw)
           nbin = nwbin(iw - 1) ; iix = (ix - 1) * nbin
           IF (ch == 2 .AND. nsolbin == 2) THEN  ! Shift the position by 15 
              iix = iix - (zoom_p1 - 1) * nsolbin
           ENDIF
            
           IF (ch == band_selectors(iw - 1) ) THEN   
              DO i = winpix(iw-1, 2)+1, winpix(iw, 1)-1 
                 IF (ALL(omi_irradiance_spec(i, iix+1:iix+nbin) > 0.0) .AND. &
                      ALL(omi_irradiance_spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0) THEN
                    nring = nring + 1

                    omirsol(wvl_idx, nring) = SUM(omi_irradiance_wavl(i, iix+1:iix+nbin)) / nbin
                    omirsol(spc_idx, nring) = SUM(omi_irradiance_spec(i, iix+1:iix+nbin)) / nbin
                   
                 ENDIF
              ENDDO
           ELSE  ! first channel 1 and second channel 2
              wcenter = (winlim(iw-1, 2) + winlim(iw, 1)) / 2.0 
              idx = MAXVAL ( MAXLOC ( omi_irradiance_wavl(spos(ch-1):epos(ch-1), iix+1), &
                   MASK = omi_irradiance_wavl(spos(ch-1):epos(ch-1), iix+1) < wcenter ) ) + spos(ch-1) - 1
               
              DO i = winpix(iw-1, 2)+1, idx 

                 IF (ALL(omi_irradiance_spec(i, iix+1:iix+nbin) > 0.0) .AND. &
                      ALL(omi_irradiance_spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0 .AND. &
                      ALL(omi_irradiance_wavl(i, iix+1:iix+nbin) > omirsol(wvl_idx, nring)) ) THEN
                    nring = nring + 1

                    omirsol(wvl_idx, nring) = SUM(omi_irradiance_wavl(i, iix+1:iix+nbin)) / nbin
                    omirsol(spc_idx, nring) = SUM(omi_irradiance_spec(i, iix+1:iix+nbin)) / nbin                   
                 ENDIF
              ENDDO
              
              omi_solring_ndiv(ix) = nring               ! 1:nring is from the same channel
              nbin = nwbin(iw) ; iix = (ix - 1) * nbin
              IF (ch == 2 .AND. nsolbin == 2) THEN  ! Shift the position by 15 
                 iix = iix - (zoom_p1 - 1) * nsolbin
              ENDIF
              
              idx = MAXVAL ( MINLOC ( omi_irradiance_wavl(spos(ch):epos(ch), iix+1),   &
                   MASK = omi_irradiance_wavl(spos(ch):epos(ch), iix+1) > wcenter ) ) + spos(ch) - 1

              DO i = idx, winpix(iw, 1) - 2 + spos(ch)                 
                 IF (ALL(omi_irradiance_spec(i, iix+1:iix+nbin) > 0.0) .AND. &
                      ALL(omi_irradiance_spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0 .AND. &
                      ALL(omi_irradiance_wavl(i, iix+1:iix+nbin) > omirsol(wvl_idx, nring)) ) THEN
                    nring = nring + 1

                    omirsol(wvl_idx, nring) = SUM(omi_irradiance_wavl(i, iix+1:iix+nbin)) / nbin
                    omirsol(spc_idx, nring) = SUM(omi_irradiance_spec(i, iix+1:iix+nbin)) / nbin                   
                 ENDIF
              ENDDO
           ENDIF
           
           idx = SUM(omi_nsolpix(1:iw-1, ix))
           omirsol(wvl_idx:spc_idx, nring+1:nring+omi_nsolpix(iw, ix)) = omispec(1, wvl_idx:spc_idx, &
                idx+1:idx+omi_nsolpix(iw, ix))
           nring = nring + omi_nsolpix(iw, ix)
        ENDDO

        ! Add extra spectra after fitting window
        noff2 = nring
        IF (ch == 2) THEN
           nring = nring + 25
        ELSE
           nring = nring + 12
        ENDIF
        
        DO i = spos(ch) + winpix(numwin, 2), epos(ch)
           IF (ch == 1 .OR. .NOT. coadd_uv2) THEN
              nbin = nxbin
           ELSE
              nbin = nxbin * ncoadd
           ENDIF
           iix = (ix - 1) * nbin
           
           IF (ALL(omi_irradiance_spec(i, iix+1:iix+nbin) > 0.0) .AND. &
                ALL(omi_irradiance_spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0) THEN
              noff2 = noff2 + 1

              omirsol(wvl_idx, noff2) = SUM(omi_irradiance_wavl(i, iix+1:iix+nbin)) / nbin
              omirsol(spc_idx, noff2) = SUM(omi_irradiance_spec(i, iix+1:iix+nbin)) / nbin
              
           ENDIF
           
           IF (noff2 == nring) EXIT
        ENDDO 
    
        omi_nsolring(ix) = nring; solring_lin(ix) = noff1; solring_uin(ix) = noff2
     ELSE
        ! Need to convole the saved solar spectra with additional slit width
        nring = omi_nsolring(ix)
        omirsol(1, 1:nring) = omi_solspec_ring(1, 1:nring, ix)
        omirsol(2, 1:nring) = omi_solspec_ring(2, 1:nring, ix)

        IF (which_slit == 0) THEN
           CALL gauss_uneven(omirsol(1, 1:nring), omirsol(2, 1:nring), nring, &
                nswath, omi_redslw(omichs(1:nswath)), retlbnd(omichs(1:nswath)), retubnd(omichs(1:nswath)))
        ELSE IF (which_slit == 3) THEN
           CALL triangle_uneven(omirsol(1, 1:nring), omirsol(2, 1:nring), nring, &
                nswath, omi_redslw(omichs(1:nswath)), retlbnd(omichs(1:nswath)), retubnd(omichs(1:nswath)))
        ELSE
           WRITE(*, *) 'This type of slit convolution is not implemented!!!'
           pge_error_status = pge_errstat_error
        ENDIF
     ENDIF

     ! Get data for surface albedo & cloud fraction at 370.2 nm +/- 15 pixels
     irefl = 0; omisolr_winpix(ix, 1:2) = 0
     IF (.NOT. coadd_uv2) THEN
        nbin = nxbin
     ELSE
        nbin = nxbin * ncoadd
     ENDIF
     nbin = nbin * nsolbin
     iix = (ix - 1) * nbin
     IF ( nsolbin == 2 ) THEN  ! Shift the position by 15 
        iix = iix - (zoom_p1 - 1) * nsolbin
     ENDIF

     idx = MAXVAL ( MINLOC ( omi_irradiance_wavl(1:nwavel, iix+1), MASK = &
          (omi_irradiance_wavl(1:nwavel, iix+1) > pos_alb - toms_fwhm * 1.4) ))
     DO i  = idx, nwavel
        IF (ALL(omi_irradiance_spec(i, iix+1:iix+nbin) > 0.0) .AND. &
             ALL(omi_irradiance_spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0) THEN
           irefl = irefl + 1
           omi_solspecr(wvl_idx, irefl, ix) = SUM(omi_irradiance_wavl(i, iix+1:iix+nbin)) / nbin
           omi_solspecr(spc_idx, irefl, ix) = SUM(omi_irradiance_spec(i, iix+1:iix+nbin)) / nbin
          
           IF (omisolr_winpix(ix, 1) == 0) omisolr_winpix(ix, 1) = i
           omisolr_winpix(ix, 2) = i
        ENDIF
        IF (irefl == nrefl) EXIT
     ENDDO
        
     IF (irefl /= nrefl) THEN
        WRITE(*, *) 'Could not get enough irradiance points for cloud fraction!!!'
        omi_solpix_errstat(ix) = pge_errstat_error; CYCLE
     ENDIF

     IF (scnwrt) THEN
       WRITE(*, *) 'End Of Reading Irradiance Spectrum: ', ix
       DO i = 1, numwin
          WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', i, winlim(i,1), winlim(i,2), omi_nsolpix(i, ix)
          IF (omi_nsolpix(i, ix) < 4) THEN
             WRITE(*, '(A,f8.3,A3,f8.3)') ' Not enough points (>=4)  in window: ', winlim(i,1), ' - ', winlim(i,2)
             pge_error_status = pge_errstat_error
          ENDIF
       ENDDO
     ENDIF

     omi_solnorm(ix) = SUM ( omispec(1, spc_idx, 1:nomi) ) / nomi
     
     IF ( omi_solnorm(ix) <= 0.0 ) THEN 
        omi_solpix_errstat(ix) = pge_errstat_error; CYCLE
     ENDIF
     
     omi_irradiance_wavl(1:nomi, ix)  = omispec(1, wvl_idx, 1:nomi)
     omi_irradiance_spec(1:nomi, ix)  = omispec(1, spc_idx, 1:nomi) / omi_solnorm(ix)
     omi_irradiance_prec(1:nomi, ix)  = omispec(1, sig_idx, 1:nomi) / omi_solnorm(ix) 
!     omi_irrad_stray(1:nomi, ix)      = strayspec(1, 1,     1:nomi) / omi_solnorm(ix)
!     omi_rad_stray  (1:nomi, ix)      = strayspec(1, 2,     1:nomi) / omi_solnorm(ix)
     
     omi_solspec_ring(1, 1:nring, ix) = omirsol(1, 1:nring)
     omi_solspec_ring(2, 1:nring, ix) = omirsol(2, 1:nring) / omi_solnorm(ix) 
     
     !DO i = 1, nomi
     !   WRITE(90, '(F10.4, 2D16.7)') omi_irradiance_wavl(i, ix), &
     !        omi_irradiance_spec(i, ix) * omi_solnorm(ix) , omi_irradiance_prec(i, ix) * omi_solnorm(ix)
     !ENDDO
     !STOP

     ! Back up solar irradiance from OMTO3
     ! IF (orbnumsol == 99999) omi_irradiance_prec(1:nomi, ix)  = 1.0
  ENDDO

  RETURN
END SUBROUTINE omi_read_irradiance_data


SUBROUTINE omi_read_radiance_lines ( iline, ny, first_line, last_line, first_pix, last_pix, nxcoadd, pge_error_status )

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx
  USE OMSAO_parameters_module, ONLY: maxchlen, maxwin, mrefl, vb_lev_omidebug, missing_value_sp
  USE OMSAO_variables_module,  ONLY: verb_thresh_lev, l1b_rad_filename, wcal_bef_coadd,   &
       numwin, coadd_uv2, band_selectors, winpix, winlim, szamax, scnwrt, currpix, &
       reduce_resolution, redlam, redsampr, reduce_slit
  USE OMSAO_omidata_module,    ONLY: nswath, mswath, omi_nradpix, nwavel_max, nxtrack_max, &
       omi_radiance_swathname, omi_radiance_spec, omi_radiance_qflg, omi_radiance_prec,   &
       omi_radiance_wavl, omi_nwav_rad, ntimes, nxtrack, nfxtrack, ncoadd,      &
       omi_specr, omi_height, omi_geoflg, omi_latitude, omi_longitude, omi_szenith,     &
       omi_sazimuth, omi_vzenith, omi_vazimuth, omi_mflg, omi_time, omi_auraalt, omi_eaza, omi_esca,  &
       omi_auralat, omi_auralon, land_water_flg, glint_flg, snow_ice_flg, omi_radiance_errstat, &
       omi_radpix_errstat, omi_nsolpix, omi_nwav_irrad, omisol_winpix, omisolr_winpix, omichs, &
       omi_saa_flag, omi_radnorm, irradwind, radwind, nybin, nxbin, reduce_lbnd, reduce_ubnd, &
       retlbnd, retubnd, zoom_mode, zoom_p1, zoom_p2, omi_redslw, omi_xtrackqflg, rowanomaly_flg, &
       waveshift_flg, blockage_flg, straysun_flg, strayearth_flg
  USE OMSAO_pixelcorner_module, ONLY: omi_alllon, omi_alllat, omi_allsza, omi_allvza, &
       omi_allaza, omi_allsca, omi_allXTrackQFlg, omi_allGeoFlg
  USE ozprof_data_module,       ONLY:  nrefl,use_flns
  USE OMSAO_errstat_module
  USE hdfeos4_parameters
  USE L1B_Reader_class
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,  INTENT (IN) :: iline, ny, first_pix, last_pix, first_line, last_line, nxcoadd

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER, INTENT (OUT) :: pge_error_status


  ! ---------------
  ! Local variables
  ! ---------------
  TYPE (l1b_block_type) :: omi_data_block
  INTEGER (KIND=i4)     :: blockline, errstat
  INTEGER, PARAMETER    :: nbits = 16
  INTEGER (KIND=i2), DIMENSION(0:nbits-1)                      :: mflgbits , tmp_mflgbits
  INTEGER (KIND=i2), DIMENSION(nxcoadd, nwavel_max, 0:nbits-1) :: flgbits
  INTEGER (KIND=i2), DIMENSION(nwavel_max)                     :: flgmsks
  INTEGER (KIND=i4), DIMENSION(maxwin)                         :: nwbin

  REAL (KIND=r8)                                         :: tmp_time
  REAL (KIND=r4)                                         :: tmp_alt, tmp_lat, tmp_lon, tmp_ExposureTime
  INTEGER (KIND=i2), DIMENSION(nxtrack_max)              :: tmp_height
  REAL (KIND=r4), DIMENSION (nwavel_max, nxtrack_max)    :: tmp_rspec, tmp_rprec, tmp_rwavl
  INTEGER (KIND=I2), DIMENSION (nwavel_max, nxtrack_max) :: tmp_rqflg
  REAL (KIND= dp)                                        :: tmpNinteg
  
  INTEGER   (KIND=4), DIMENSION (mswath)    :: nwls
  INTEGER   (KIND=i2), DIMENSION(mswath)    :: spos, epos
  INTEGER, DIMENSION (nwavel_max)           :: idxs
  INTEGER                                   :: nwavel, is, iloop, nwl, i, j, ix, iix, ii, &
       nomi, fidx, lidx, ch, iw, ic, idx, irefl, nx, nbin, fpix, lpix, np, npos
  LOGICAL                                   :: error
  LOGICAL, DIMENSION (maxwin, nxtrack_max)  :: wavcals 
  REAL (KIND = dp)                          :: tmpsampr, retswav, retewav
  REAL (KIND = dp), DIMENSION (maxwin, nxtrack_max, nxcoadd) :: wshis, wsqus  
  REAL (KIND = dp), DIMENSION (nxcoadd, sig_idx, nwavel_max) :: omispec
  REAL (KIND = dp), DIMENSION (sig_idx, nwavel_max, nxtrack_max) :: tmpspec

  LOGICAL                                   :: correct_merr = .TRUE.


  ! Exteranl functions
  INTEGER                                   :: estat
  
  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=23), PARAMETER :: modulename = 'omi_read_radiance_lines'
 
  pge_error_status = pge_errstat_ok;  errstat = omi_s_success

  ! Initialize all local data arrays
  nx = nfxtrack / nxbin
  
  omi_latitude  (1:nx, 0:ny-1) = omi_alllat (1:nx, iline:iline+ny-1)
  omi_longitude (1:nx, 0:ny-1) = omi_alllon (1:nx, iline:iline+ny-1)
  omi_szenith   (1:nx, 0:ny-1) = omi_allsza (1:nx, iline:iline+ny-1)
  omi_vzenith   (1:nx, 0:ny-1) = omi_allvza (1:nx, iline:iline+ny-1)
  omi_eaza      (1:nx, 0:ny-1) = omi_allaza (1:nx, iline:iline+ny-1)
  omi_esca      (1:nx, 0:ny-1) = omi_allsca (1:nx, iline:iline+ny-1)
  omi_xtrackqflg(1:nx, 0:ny-1) = omi_allXTrackQFlg(1:nx, iline:iline+ny-1)
  omi_geoflg(1:nx, 0:ny-1)     = omi_allGeoFlg(1:nx, iline:iline+ny-1)
  ! Derive Row Anomaly Related Flags
  DO i = 0, ny - 1
     CALL convert_xtrackqfag_info ( nx, omi_xtrackqflg(1:nx, i), rowanomaly_flg(1:nx, i), &
          waveshift_flg(1:nx, i), blockage_flg(1:nx, i), straysun_flg(1:nx, i), strayearth_flg(1:nx, i) )    
  ENDDO
  !print *, nx, ny
  !DO ix = 1, nx
  !   WRITE(*, '(I5, 6F10.4)') ix, omi_longitude(ix, 0), omi_latitude(ix, 0), omi_szenith(ix, 0), &
  !        omi_vzenith(ix, 0), omi_eaza(ix, 0), omi_esca(ix, 0)
  !ENDDO
  !STOP
  IF  (.NOT. use_flns  ) then 
  correct_merr = .FALSE.
  print * , 'correct_merr set to "F" because of not using floor noise'
  ENDIF
  omi_time      (0:ny-1)       = 0.0
  omi_auralat   (0:ny-1)       = 0.0
  omi_auralon   (0:ny-1)       = 0.0
  omi_auraalt   (0:ny-1)       = 0.0
  omi_height    (1:nx, 0:ny-1) = 0.0
  omi_geoflg    (1:nx, 0:ny-1) = 0
  omi_radiance_spec (1:nwavel_max, 1:nxtrack, 0:ny-1) = 0.0
  omi_radiance_prec (1:nwavel_max, 1:nxtrack, 0:ny-1) = 0.0
  omi_radiance_qflg (1:nwavel_max, 1:nxtrack, 0:ny-1) = 0
  omi_radiance_wavl (1:nwavel_max, 1:nxtrack, 0:ny-1) = 0.0
  tmp_rspec (1:nwavel_max, 1:nxtrack) = 0.0
  tmp_rprec (1:nwavel_max, 1:nxtrack) = 0.0
  tmp_rqflg (1:nwavel_max, 1:nxtrack) = 0.0
  tmp_rwavl (1:nwavel_max, 1:nxtrack) = 0.0

  omi_radiance_errstat(0:ny-1)          = pge_errstat_ok
  omi_radpix_errstat(1:nxtrack, 0:ny-1) = pge_errstat_ok
  omi_nradpix = 0;  omi_nwav_rad = 0

  j = 1; nwavel = 0;  wavcals = .TRUE.

  DO is = 1, nswath

     ch = omichs(is)

     ! Open data block called 'omi_data_block' with default size of 100 lines
     errstat = L1Br_OPEN ( omi_data_block, l1b_rad_filename, omi_radiance_swathname(is), ny )
     IF ( errstat /= omi_s_success ) THEN
        estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, 'L1Br_OPEN failed.', modulename, 0 ) ; STOP 1
     END IF
     
     DO iloop = 0, ny - 1
        ! The current scan line number we are reading
        mflgbits = 0
        
        DO i = 0, nybin - 1
           blockline = first_line + iloop * nybin + i

           ! Get geolocation fields for UV-1 (if both UV-1 and UV-2 are selected)  
           IF ((coadd_uv2 .AND. is == 1) .OR. .NOT. coadd_uv2) THEN    
              errstat = L1Br_getGEOline ( omi_data_block, blockline,          &
                   Time_k                    = tmp_time,                      &
                   SpacecraftAltitude_k      = tmp_alt,                       &
                   SpacecraftLongitude_k     = tmp_lon,                       &
                   SpacecraftLatitude_k      = tmp_lat,                       &
                   TerrainHeight_k           = tmp_height(1:nfxtrack),        &
                   GroundPixelQualityFlags_k = omi_geoflg(1:nfxtrack, iloop))
              
              IF ( errstat /= omi_s_success ) THEN
                 estat = OMI_SMF_setmsg (omsao_e_read_l1b_file, 'L1Br_getGEOline failed.', modulename, 0); STOP 1
              END IF

              omi_time(iloop)    = tmp_time !+ omi_time(iloop) 
              omi_auraalt(iloop) = omi_auraalt(iloop) + tmp_alt
              omi_auralon(iloop) = omi_auralon(iloop) + tmp_lon
              omi_auralat(iloop) = omi_auralat(iloop) + tmp_lat
              omi_height(1:nfxtrack, iloop) = omi_height(1:nfxtrack, iloop) + tmp_height(1:nfxtrack)
                         
              ! Could not be properly coadded, just use the last line of these binned lines
              IF (i == nybin - 1) CALL convert_gpqualflag_info (nfxtrack, omi_geoflg(1:nfxtrack, iloop), &
                   land_water_flg(1:nfxtrack, iloop), glint_flg(1:nfxtrack, iloop), snow_ice_flg(1:nfxtrack, iloop))
           ENDIF

           errstat = L1Br_getDATA ( omi_data_block, blockline, MeasurementQualityFlags_k = omi_mflg)
           IF( errstat /= omi_s_success ) THEN
              estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, "L1Br_getDATA failed.", modulename, 0 ); STOP 1
           END IF
         
           CALL convert_2bytes_to_16bits ( nbits, 1, omi_mflg, tmp_mflgbits(0:nbits-1))
           mflgbits = mflgbits + tmp_mflgbits(0:nbits-1)
           
           ! Get radiances associated with wavelength range
           errstat = L1Br_getSIGline ( omi_data_block, blockline,   &
                Signal_k            = tmp_rspec(j:, :),  &
                SignalPrecision_k   = tmp_rprec(j:, :),  &
                PixelQualityFlags_k = tmp_rqflg(j:, :),  &
                Wavelength_k        = tmp_rwavl(j:, :),  &
                Nwl_k  = nwls(ch) )
           IF( errstat /= omi_s_success ) THEN
              estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, "L1Br_getSIGline failed.", modulename, 0); STOP 1
           END IF

           IF (correct_merr) THEN
              errstat = L1Br_getDATA ( omi_data_block, blockline,   &
                   ExposureTime_k      = tmp_ExposureTime)
              IF( errstat /= omi_s_success ) THEN
                 estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, "L1Br_getDATA failed.", modulename, 0); STOP 1
              END IF
              tmpNinteg = 2.d0 / tmp_ExposureTime
!              PRINT * , 'tmpNinteg', tmpNinteg, tmp_exposureTime
           ENDIF
           
           nwl = j + nwls(ch) - 1

           !print *, j, nwl, nwls(ch), ny, iloop
           !print *, tmp_rwavl(j, 6), tmp_rwavl(nwl, 6)
           !print *, omi_radiance_wavl(j, 6, iloop), omi_radiance_wavl(nwl, 6, iloop)

           omi_radiance_spec(j:nwl, :, iloop) = omi_radiance_spec(j:nwl, :, iloop) + tmp_rspec(j:nwl, :)
           IF (correct_merr) THEN
              omi_radiance_prec(j:nwl, :, iloop) = omi_radiance_prec(j:nwl, :, iloop) + tmp_rprec(j:nwl, :) / SQRT( tmpNinteg ) 
           ELSE
              omi_radiance_prec(j:nwl, :, iloop) = omi_radiance_prec(j:nwl, :, iloop) + tmp_rprec(j:nwl, :)
           ENDIF
           omi_radiance_wavl(j:nwl, :, iloop) = omi_radiance_wavl(j:nwl, :, iloop) + tmp_rwavl(j:nwl, :)

           DO ix = 1, nxtrack
              CALL coadd_2bytes_qflgs(nbits, nwls(ch), omi_radiance_qflg(j:nwl, ix, iloop), tmp_rqflg(j:nwl, ix)) 
           ENDDO
           
           IF( errstat /= omi_s_success ) THEN
              estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, "L1Br_getSIGline failed.", modulename, 0); STOP 1
           END IF
        ENDDO

        IF ((coadd_uv2 .AND. is == 1) .OR. .NOT. coadd_uv2) THEN  
           !omi_time(iloop)    = omi_time(iloop)    / nybin
           omi_auraalt(iloop) = omi_auraalt(iloop) / nybin
           omi_auralon(iloop) = omi_auralon(iloop) / nybin
           omi_auralat(iloop) = omi_auralat(iloop) / nybin
           omi_height(1:nfxtrack, iloop) = NINT (1.0 * omi_height(1:nfxtrack, iloop) / nybin)
        ENDIF
        
        IF (mflgbits(0) >= 1 .OR. mflgbits(1) >= 1 .OR. mflgbits(3) >= 1 .OR. mflgbits(12) >= 1) THEN
           WRITE(*, *) 'All radiances could not be used: line ', blockline, ' Swath ', is
           omi_radiance_errstat(iloop) = pge_errstat_error      
        ELSE IF (ANY(mflgbits >= 1)) THEN
           !WRITE(*, *) 'Warning set on all radiances: line',  blockline, ' Swath ', is
           IF (omi_radiance_errstat(iloop) /= pge_errstat_error) &
                omi_radiance_errstat(iloop) = pge_errstat_warning      
       
           ! Over SAA region
           IF (mflgbits(10) >= 1) omi_saa_flag(iloop) = 1
        ENDIF
        
        omi_radiance_spec(j:nwl, :, iloop) = omi_radiance_spec(j:nwl, :, iloop) / nybin
        IF ( correct_merr ) then !jbak
            omi_radiance_prec(j:nwl, :, iloop) = omi_radiance_prec(j:nwl, :, iloop) / nybin / SQRT(1.0D0 * nybin)
        ELSE
            omi_radiance_prec(j:nwl, :, iloop) = omi_radiance_prec(j:nwl, :, iloop) / nybin !/ SQRT(1.0D0 * nybin)
        ENDIF
        omi_radiance_wavl(j:nwl, :, iloop) = omi_radiance_wavl(j:nwl, :, iloop) / nybin 
     END DO     ! end iloop
     
     nwavel = nwavel + nwls(ch); spos(ch) = j; j = nwavel + 1; epos(ch) = nwavel
     
     ! Close data block structure
     errstat = L1Br_CLOSE ( omi_data_block )
     IF ( errstat /= omi_s_success .AND. verb_thresh_lev >= vb_lev_omidebug ) THEN
        estat = OMI_SMF_setmsg ( omsao_w_clos_l1b_file, 'L1Br_CLOSE failed.', modulename, 0 )
     END IF

     ! Sort data in wavelength increasing order   
     IF (omi_radiance_wavl(spos(ch), 1, 0) > omi_radiance_wavl(epos(ch), 1, 0)) THEN
        idxs(spos(ch):epos(ch)) = (/ (i, i = epos(ch), spos(ch), -1) /)
        omi_radiance_wavl(spos(ch):epos(ch), :, :) = omi_radiance_wavl(idxs(spos(ch):epos(ch)), :, :)
        omi_radiance_spec(spos(ch):epos(ch), :, :) = omi_radiance_spec(idxs(spos(ch):epos(ch)), :, :)
        omi_radiance_prec(spos(ch):epos(ch), :, :) = omi_radiance_prec(idxs(spos(ch):epos(ch)), :, :)
        omi_radiance_qflg(spos(ch):epos(ch), :, :) = omi_radiance_qflg(idxs(spos(ch):epos(ch)), :, :)     
     ENDIF

     !WRITE(*, '(10F10.4)') omi_radiance_wavl(150, 1:60, 0)

     ! Aligh zoom mode data
     IF (zoom_mode .AND. ch == 2) THEN
        omi_radiance_wavl(spos(ch):epos(ch), zoom_p1:zoom_p2, :) = &
             omi_radiance_wavl(spos(ch):epos(ch), 1:nxtrack/2, :)
        omi_radiance_spec(spos(ch):epos(ch), zoom_p1:zoom_p2, :) = &
             omi_radiance_spec(spos(ch):epos(ch), 1:nxtrack/2, :)
        omi_radiance_prec(spos(ch):epos(ch), zoom_p1:zoom_p2, :) = &
             omi_radiance_prec(spos(ch):epos(ch), 1:nxtrack/2, :)
        omi_radiance_qflg(spos(ch):epos(ch), zoom_p1:zoom_p2, :) = &
             omi_radiance_qflg(spos(ch):epos(ch), 1:nxtrack/2, :)  
     ENDIF

     !IF ( ch == 2 ) THEN
     !   DO iloop = 0, ny - 1 
     !      CALL corruv2wav(nwls(ch), nxtrack, omi_radiance_wavl(spos(ch):epos(ch), 1:nxtrack, iloop)) 
     !   ENDDO
     !ENDIF

     !print *, ny
     !IF (is == 1) THEN 
     !   nx = 30
     !ELSE
     !   nx = 60
     !ENDIF
     !WRITE(90, *) nx, nwls(ch)
     !DO ix = 1, nx
     !   WRITE(90, *) ix
     !   DO iw = spos(ch), epos(ch)
     !      WRITE(90, '(F10.5,2D14.6,I6)') omi_radiance_wavl(iw, ix, 0), omi_radiance_spec(iw, ix, 0),&
     !           omi_radiance_prec(iw, ix, 0), omi_radiance_qflg(iw, ix, 0) !, mflgbits(0:nbits-1)
     !   ENDDO
     !ENDDO

  ENDDO ! end swath loop

  IF (zoom_mode .AND. nswath == 1) THEN
     omi_height(zoom_p1:zoom_p2, 0:ny-1) = omi_height(1:nxtrack/2, 0:ny-1)
     glint_flg (zoom_p1:zoom_p2, 0:ny-1) = glint_flg  (1:nxtrack/2, 0:ny-1)
     land_water_flg(zoom_p1:zoom_p2, 0:ny-1) = land_water_flg(1:nxtrack/2, 0:ny-1)
     snow_ice_flg(zoom_p1:zoom_p2, 0:ny-1) = snow_ice_flg(1:nxtrack/2, 0:ny-1)     
  ENDIF

  ! Degrade spectral resolution if necessary
  IF (reduce_resolution) THEN
     nwavel = 0; j = 1
     DO is = 1, nswath
        ch = omichs(is)
        
        IF (coadd_uv2 .AND. is == 1) THEN
           npos = nfxtrack
        ELSE
           npos = nxtrack
        ENDIF
        fpix = first_pix; lpix = last_pix
        IF (is == 2) THEN
           fpix = first_pix * 2 -1; lpix = last_pix * 2
        ENDIF
        !IF (nswath == 2) THEN
        !   retswav = retlbnd(ch); retewav = retubnd(ch)
        !ELSE
        !   retswav = retlbnd(band_selectors(1)); retewav = retubnd(band_selectors(1))
        !ENDIF
        retswav = retlbnd(ch); retewav = retubnd(ch)
        tmpsampr = redsampr; IF (is == 1 .AND. band_selectors(1) == 1) tmpsampr = redsampr / 3.0
        np = nwls(ch)

        DO iloop = 0, ny - 1
           tmpspec(wvl_idx, 1:np, fpix:lpix) = omi_radiance_wavl(spos(ch):epos(ch), fpix:lpix, iloop)
           tmpspec(spc_idx, 1:np, fpix:lpix) = omi_radiance_spec(spos(ch):epos(ch), fpix:lpix, iloop)
           tmpspec(sig_idx, 1:np, fpix:lpix) = omi_radiance_prec(spos(ch):epos(ch), fpix:lpix, iloop)
           tmp_rqflg(       1:np, fpix:lpix) = omi_radiance_qflg(spos(ch):epos(ch), fpix:lpix, iloop)       
           
           DO ix = fpix, lpix
              CALL convert_2bytes_to_16bits ( nbits, np, tmp_rqflg(1:np, ix), flgbits(1, 1:np, 0:nbits-1))
              tmp_rqflg(1:np, ix) = flgbits(1, 1:np, 0) &   ! Missing
                   + flgbits(1, 1:np, 1)                &   ! Bad
                   + flgbits(1, 1:np, 2)                &   ! Processing error
                   + flgbits(1, 1:np, 4)                &   ! RTS_Pixel_Warning Flag
                   + flgbits(1, 1:np, 5)                &   ! Saturation Possibility Flag
                   + flgbits(1, 1:np, 7)                    ! Dark Current Warning Flag
           ENDDO
           
           CALL reduce_rad_resolution (tmpspec(:, 1:np, fpix:lpix), tmp_rqflg(1:np, fpix:lpix),   &
                np, lpix-fpix+1, reduce_slit, omi_redslw(is), tmpsampr, redlam, retswav, retewav, reduce_lbnd(ch), &
                reduce_ubnd(ch), nwls(ch), pge_error_status)
        
           IF (pge_error_status == pge_errstat_error) RETURN

           nwavel = j + nwls(ch) - 1
           omi_radiance_wavl(j:nwavel, fpix:lpix, iloop) = tmpspec(wvl_idx, 1:nwls(ch), fpix:lpix)
           omi_radiance_spec(j:nwavel, fpix:lpix, iloop) = tmpspec(spc_idx, 1:nwls(ch), fpix:lpix)
           omi_radiance_prec(j:nwavel, fpix:lpix, iloop) = tmpspec(sig_idx, 1:nwls(ch), fpix:lpix)  
           omi_radiance_qflg(j:nwavel, fpix:lpix, iloop) = 0   ! All data are good  (pre filtered)  
        ENDDO
        spos(ch) = j; epos(ch) = nwavel; j = nwavel + 1
     ENDDO
  ENDIF

  IF (nwavel > nwavel_max) THEN
     WRITE(*, *) "Need to increase nwavel_max!!!"
     pge_error_status = pge_errstat_error; RETURN
  ENDIF  

  ! Resample the water/land, sea glint, and snow/oce flags if xbin, surface altitude
  IF (nxbin > 1) THEN
     DO ix = 1, nfxtrack / nxbin
        iix = (ix - 1) * nxbin + NINT(nxbin / 2.0)
        land_water_flg(ix, 0:ny-1) = land_water_flg(iix, 0:ny-1)
        snow_ice_flg(ix, 0:ny-1) = snow_ice_flg(iix, 0:ny-1) 
        
        iix = (ix - 1) * nxbin
        DO i = 0, ny - 1
           omi_height(ix, i) = NINT(SUM(1.0 * omi_height(iix+1:iix+nxbin, i)) / nxbin)
           glint_flg(ix, i)  = NINT(SUM(1.0 * glint_flg(iix+1:iix+nxbin, i))  / nxbin)
        ENDDO   
     ENDDO
  ENDIF

  ! Determine number of binning for different fitting windows
  DO iw = 1, numwin
     ch = band_selectors(iw)
     IF (ch == 1 .OR. .NOT. coadd_uv2) THEN
        nwbin(iw) = nxbin
     ELSE
        nwbin(iw) = nxbin * ncoadd
     ENDIF
  ENDDO

  ! Subset and coadd radiance spectrum
  DO iloop = 0, ny - 1 

     IF (omi_radiance_errstat(iloop) == pge_errstat_error) THEN
        omi_radpix_errstat(first_pix:last_pix, iloop) = pge_errstat_error; CYCLE
     ENDIF
     
     DO ix = first_pix, last_pix
        currpix = ix

        IF (omi_szenith(ix, iloop) > szamax .OR. omi_szenith(ix, iloop) < 0 ) THEN
           omi_radpix_errstat(ix, iloop) = pge_errstat_error; CYCLE
        ENDIF

        IF (rowanomaly_flg(ix, iloop) == 1) THEN
           print * , 'rowanomaly_flg', ix, iloop
           omi_radpix_errstat(ix, iloop) = pge_errstat_error; CYCLE
        ENDIF

        ! Get quality flag bits
        ! Coadd uv-2 flags if necessary to avoid coadding inconsistent # of pixels 
        flgmsks = 0
        DO is = 1, nswath
           ch = omichs(is)
           IF (is == 1) THEN
              nbin = nxbin
           ELSE
              nbin = nxbin * ncoadd
           ENDIF
           iix = (ix - 1) * nbin 

           ! properly align cross track positions to be coadded (should be within one pixel)
           IF (.NOT. reduce_resolution) THEN
              IF (nbin > 2) CALL prespec_align(nwls(ch), nbin, omi_radiance_wavl(spos(ch):epos(ch),&
                   iix+1:iix+nbin, iloop), omi_radiance_spec(spos(ch):epos(ch), iix+1:iix+nbin, iloop), &
                   omi_radiance_prec(spos(ch):epos(ch), iix+1:iix+nbin, iloop), &       
                   omi_radiance_qflg(spos(ch):epos(ch), iix+1:iix+nbin, iloop))
              
              DO ic = 1, nbin
                 CALL convert_2bytes_to_16bits ( nbits, nwls(ch), omi_radiance_qflg(spos(ch):epos(ch), &
                      iix + ic, iloop), flgbits(ic, spos(ch):epos(ch), 0:nbits-1))
                 flgmsks(spos(ch):epos(ch)) = flgmsks(spos(ch):epos(ch)) &
                      + flgbits(ic, spos(ch):epos(ch), 0)                &   ! Missing
                      + flgbits(ic, spos(ch):epos(ch), 1)                &   ! Bad 
                      + flgbits(ic, spos(ch):epos(ch), 2)                &   ! Processing error
                      + flgbits(ic, spos(ch):epos(ch), 4)                &   ! RTS_Pixel_Warning Flag
                      + flgbits(ic, spos(ch):epos(ch), 5)                &   ! Saturation Possibility Flag
                      + flgbits(ic, spos(ch):epos(ch), 7)                     ! Dark Current Warning Flag
              ENDDO

           !DO i = spos(ch), epos(ch)
           !   WRITE(91, '(F10.4, D14.6, 16I2)') omi_radiance_wavl(i, iix+1, iloop), &
           !        omi_radiance_spec(i, iix+1, iloop), flgbits(1, i, 0:nbits-1)
           !ENDDO

           ELSE
              ! Already aligned because of using common wavelength scale
              DO ic = 1, nbin
                 flgmsks(spos(ch):epos(ch)) = flgmsks(spos(ch):epos(ch)) + &
                      omi_radiance_qflg(spos(ch):epos(ch), iix + ic, iloop)
              ENDDO
           ENDIF
        ENDDO

        ! Subset valid data
        nomi = 0; omispec = 0.0; fidx = 1
        DO iw = 1, numwin
           ch = band_selectors(iw)                  
           omi_nradpix(iw, ix, iloop) = nomi
           nbin = nwbin(iw);  iix = (ix - 1) * nbin

           !fidx = omisol_winpix(iw, ix, 1) 
           !lidx = omisol_winpix(iw, ix, 2) 
           lidx = fidx + omi_nsolpix(iw, ix) - 1
           
           DO ii = fidx, lidx
              i = irradwind(ii, ix)
              IF (ALL(omi_radiance_spec(i, iix+1:iix+nbin, iloop) > 0.0) .AND. &
                   ALL(omi_radiance_spec(i, iix+1:iix+nbin, iloop) < 4.0E14) .AND. flgmsks(i) == 0 ) THEN
                 nomi = nomi + 1
                 omispec(1:nbin, wvl_idx, nomi) = omi_radiance_wavl(i, iix+1:iix+nbin, iloop)
                 omispec(1:nbin, spc_idx, nomi) = omi_radiance_spec(i, iix+1:iix+nbin, iloop)
                 omispec(1:nbin, sig_idx, nomi) = omi_radiance_prec(i, iix+1:iix+nbin, iloop)
                 radwind(nomi, ix, iloop) = ii
              ENDIF
           ENDDO
           fidx = lidx + 1
           omi_nradpix(iw, ix, iloop) = nomi - omi_nradpix(iw, ix, iloop)

           ! If the # of wavelengths is <= 75% of the # of irradiances, stop processing this pixel
           ! 90 % to 80 %
           
         !  WRITE(*, '(A,3I5,F9.2)') 'rad/irrad waves: ', iw, &
         !           omi_nradpix(iw, ix, iloop), omi_nsolpix(iw, ix),omi_nradpix(iw, ix, iloop)*100./omi_nsolpix(iw, ix)
           IF (omi_nradpix(iw, ix, iloop) <= omi_nsolpix(iw, ix) * 0.70 ) THEN
              WRITE(*, '(A,5I5,F9.2)') 'Too fewer radiance points: ', ix, iloop, iw, &
                    omi_nradpix(iw, ix, iloop), omi_nsolpix(iw, ix),omi_nradpix(iw, ix, iloop)*100./omi_nsolpix(iw, ix)
              omi_radpix_errstat(ix, iloop) = pge_errstat_error; EXIT
           ENDIF

           !WRITE(*, '(2I5, 2F8.3, 2I5, 2F8.3, 6I5)') ix, iw, omi_radiance_wavl(spos(ch), iix+1, iloop), &
           !     omi_radiance_wavl(epos(ch), iix+1, iloop), spos(ch), epos(ch),  &
           !     fidx, lidx, lidx - fidx + 1, omi_nradpix(iw, ix, iloop)
        ENDDO
          
        IF (omi_radpix_errstat(ix, iloop) == pge_errstat_error) CYCLE   ! This pixel will not be processed.     
        omi_nwav_rad(ix, iloop) = nomi
    
        ! Perform coadding if UV-2 is selected with UV-1
        fidx = 1
        DO iw = 1, numwin
           ch = band_selectors(iw);  nbin = nwbin(iw)
           lidx = fidx + omi_nradpix(iw, ix, iloop) - 1 

           IF (nbin > 1) THEN
              CALL radwavcal_coadd(wcal_bef_coadd, wavcals(iw, ix), iw, ix, omi_nradpix(iw, ix, iloop), nbin, &
                   omispec(1:nbin, :, fidx:lidx), wshis(iw, ix, 1:nbin), wsqus(iw, ix, 1:nbin), error)
              wavcals(iw, ix) = .FALSE.
              IF (error) THEN
                 WRITE(*, '(A)') 'No radiance wavelength calibration before coadding!!!'
                 pge_error_status = pge_errstat_warning
              ENDIF
           ENDIF
           fidx = lidx + 1
        ENDDO

        ! Get data for surface albedo & cloud fraction at 370.2 nm +/- 20 pixels
        irefl = 0; fidx = omisolr_winpix(ix, 1)
        IF (.NOT. coadd_uv2) THEN
           nbin = nxbin
        ELSE
           nbin = nxbin * ncoadd
        ENDIF
        iix = (ix - 1) * nbin
        
        DO i  = fidx, nwavel
           IF ( ALL(omi_radiance_spec(i, iix+1:iix+nbin, iloop) > 0.0) .AND. &
                ALL(omi_radiance_spec(i, iix+1:iix+nbin, iloop) < 4.0E14) .AND. flgmsks(i) == 0) THEN
              irefl = irefl + 1
              omi_specr(wvl_idx, irefl, ix, iloop) = SUM(omi_radiance_wavl(i, iix+1:iix+nbin, iloop)) / nbin
              omi_specr(spc_idx, irefl, ix, iloop) = SUM(omi_radiance_spec(i, iix+1:iix+nbin, iloop)) / nbin
           ENDIF
           IF (irefl == nrefl) EXIT
        ENDDO

        IF (irefl /= nrefl) THEN
           WRITE(*, '(A, 2I5, F9.2)') 'Number of rad/sol points (cloud fraction) do not match: ', &
                ix, iloop, omi_szenith(ix, iloop)
           omi_radpix_errstat(ix, iloop) = pge_errstat_error  
        ENDIF
   
        !IF (scnwrt) THEN
        !   WRITE(*, *) 'End Of Reading Radiance Spectrum: ', ix, iloop + iline
        !   DO i = 1, numwin
        !      WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', i, winlim(i, 1), winlim(i, 2), omi_nradpix(i, ix, iloop)
        !   
        !      IF (omi_nradpix(i, ix, iloop) <= 20) THEN
        !         WRITE(*, '(A,f8.3,A3,f8.3)') ' Not enough points in window: ', winlim(i, 1), ' - ', winlim(i, 2)
        !         pge_error_status = pge_errstat_error
        !        
        !      ENDIF
        !   ENDDO
        !   print *, omispec(1, 1, 138:139)
        !ENDIF

        omi_radnorm(ix, iloop) = SUM ( omispec(1, spc_idx, 1:nomi) ) / nomi
        !omi_radnorm(ix, iloop) = 1.0E11 
        IF ( omi_radnorm(ix, iloop) <= 0.0 ) THEN 
           pge_error_status = pge_errstat_error; RETURN
        ENDIF
        omi_radiance_wavl(1:nomi, ix, iloop) = omispec(1, wvl_idx, 1:nomi)
        omi_radiance_spec(1:nomi, ix, iloop) = omispec(1, spc_idx, 1:nomi) / omi_radnorm(ix, iloop)
        omi_radiance_prec(1:nomi, ix, iloop) = omispec(1, sig_idx, 1:nomi) / omi_radnorm(ix, iloop)       
     ENDDO

  ENDDO
  
  !DO ix = 1, nfxtrack
  !   WRITE(90, '(10I5)') ix, omi_nsolpix(1:numwin, ix), omi_nwav_irrad(ix)
  !   DO i = 0, nloop - 1
  !      WRITE(90, '(I5, F10.3, 3I5)') i, omi_szenith(ix, i), omi_nradpix(1:numwin, ix, i), omi_nwav_rad(ix, i)
  !   ENDDO
  !ENDDO

  RETURN
END SUBROUTINE omi_read_radiance_lines

! Replace Solar Composite with original OMI solar irradiance
SUBROUTINE replace_solar_irradiance (lun, nxcoadd, first_pix, last_pix, pge_error_status )

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx
  USE OMSAO_parameters_module, ONLY: maxchlen, maxwin, mrefl
  USE OMSAO_variables_module,  ONLY: avg_solcomp, avgsol_allorb, numwin, coadd_uv2, band_selectors, refdbdir
  USE OMSAO_omidata_module,    ONLY: nswath, mswath, orbnum, nxtrack_max, nxtrack, nfxtrack, nwavel_max,   &
       omi_irradiance_wavl, omi_irradiance_spec, ncoadd, omi_nsolpix, omi_solnorm, omi_solspecr, nxbin, &
       omi_nwav_irrad, omi_solspec_ring, solring_lin, solring_uin, omi_nsolring, omi_solring_ndiv
  USE ozprof_data_module,      ONLY: nrefl
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,  INTENT (IN) :: first_pix, last_pix, lun, nxcoadd

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER, INTENT (OUT) :: pge_error_status

  ! ----------------
  ! Loacal variables
  ! ----------------
  INTEGER, PARAMETER :: ntype = 2, norbtype = 2, mscpt = 7000
  CHARACTER(LEN=3), DIMENSION(ntype) :: comp_types = (/'med', 'pc0'/)  ! average and 1st principal compoment
  CHARACTER(LEN=3), DIMENSION(mswath):: channels   = (/'uv1', 'uv2'/)  ! average and 1st principal compoment
  CHARACTER(LEN=4), DIMENSION(ntype) :: orb_types  = (/'comp','1st7'/) ! all orbits and first 7 days after OPF change

  ! For the solar composite data
  INTEGER,        DIMENSION(mswath)                      :: nscpts   
  REAL (KIND=dp), DIMENSION(mswath)                      :: snorms                
  REAL (KIND=dp), DIMENSION(mswath, nxtrack_max, mscpt)  :: solcomp
  REAL (KIND=dp), DIMENSION(mswath, mscpt)               :: solcomp_wvl

  INTEGER :: i, j, ch, nx, errstat, ix, iix, iw, fidx, lidx, sidx, eidx, npts, npts1, nbin
  INTEGER (KIND=i4), DIMENSION(maxwin)  :: nwbin
  REAL (KIND=dp)                        :: swav, ewav
  CHARACTER(LEN=maxchlen)               :: scfname
  CHARACTER(LEN=3)                      :: chc, typec
  CHARACTER(LEN=4)                      :: orbtypec, opfc
  REAL (KIND=dp), DIMENSION(nwavel_max) :: tmpspec, tmpwvl

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg
  
  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=24), PARAMETER :: modulename = 'replace_solar_irradiance'

  pge_error_status = pge_errstat_ok

  ! Selectively read the solar composite data
  DO i = 1, nswath
     IF (nswath == 1) THEN
        ch = band_selectors(1)
     ELSE
        ch = i
     ENDIF

     chc = channels(ch)
     IF (avg_solcomp) THEN 
        typec = comp_types(1)
     ELSE
        typec = comp_types(2)
     ENDIF

     IF (orbnum >= 6551) THEN
        opfc = 'ge25'
     ELSE
        opfc = 'lt25'
     ENDIF

     IF (avgsol_allorb) THEN 
        orbtypec = orb_types(1)
     ELSE
        orbtypec = orb_types(2)
     ENDIF

     scfname = ADJUSTL(TRIM(refdbdir)) // 'OMI/SolarComposite/' // chc // '_' // typec // '_' &
          // opfc // '_' // orbtypec // '.dat'
     !print *, i, ch, ADJUSTL(TRIM(scfname))
     
     OPEN (UNIT=lun, FILE=TRIM(ADJUSTL(scfname)), STATUS='UNKNOWN', IOSTAT=errstat)
     IF ( errstat /= pge_errstat_ok ) THEN
        WRITE(*, '(2A)') modulename, ': Cannot open solar composite file!!!'
        pge_error_status = pge_errstat_error; RETURN
     END IF

     READ(lun, *) nx, nscpts(ch), swav, ewav, snorms(ch)
     IF (  (ch == 1 .AND. nx /= nxtrack / ncoadd) .OR. (ch == 2 .AND. nx /= nxtrack) ) THEN
        WRITE(*, '(2A)') modulename, ': Solar composite does not cover all xtrack positions!!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
     IF (nscpts(ch) > mscpt) THEN
        WRITE(*, '(2A)') modulename, ': Increase the dimension mscpt for solar composite!!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF

     DO j = 1, nscpts(ch)
        READ(lun, *) solcomp_wvl(ch, j), solcomp(ch, 1:nx, j)
     ENDDO     
     CLOSE(lun)         
  ENDDO

  ! Determine number of binning for different fitting windows
  DO iw = 1, numwin
     ch = band_selectors(iw)
     IF (ch == 1 .OR. .NOT. coadd_uv2) THEN
        nwbin(iw) = nxbin
     ELSE
        nwbin(iw) = nxbin * ncoadd
     ENDIF
  ENDDO
  
  DO ix = first_pix, last_pix
     fidx = 1
     DO iw = 1, numwin
        ch = band_selectors(iw)
        nbin = nwbin(iw);  iix = (ix - 1) * nbin
        npts1 = omi_nsolpix(iw, ix)
        lidx = fidx + npts1 - 1
        
        ! Check boundary levels
        IF ( omi_irradiance_wavl(fidx, ix) < solcomp_wvl(ch, 1) .OR. &
             omi_irradiance_wavl(lidx, ix) > solcomp_wvl(ch, nscpts(ch))) THEN
           WRITE(*, '(2A)') modulename, ': Reduce fitting window or switch off solar composite!!!'
           WRITE(*, '(3I5, 4F10.4)') ch, iw, ix, solcomp_wvl(ch, 1), solcomp_wvl(ch, nscpts(ch)), &
                omi_irradiance_wavl(fidx, ix), omi_irradiance_wavl(lidx, ix)
           pge_error_status = pge_errstat_error; RETURN
        ENDIF

        ! Replace Solar irradiance
        sidx = MINVAL(MAXLOC(solcomp_wvl(ch, 1:nscpts(ch)), MASK=(solcomp_wvl(ch, 1:nscpts(ch)) &
             <= omi_irradiance_wavl(fidx, ix))))
        eidx = MINVAL(MINLOC(solcomp_wvl(ch, 1:nscpts(ch)), MASK=(solcomp_wvl(ch, 1:nscpts(ch)) &
             >= omi_irradiance_wavl(lidx, ix)))) 
        npts = eidx - sidx + 1

        !print *, fidx, lidx, sidx, eidx
        !print *, solcomp_wvl(ch, sidx), solcomp_wvl(ch, eidx)
        !print *, omi_irradiance_wavl(fidx, ix), omi_irradiance_wavl(lidx, ix)

        tmpwvl(1:npts1) = omi_irradiance_wavl(fidx:lidx, ix)      
        omi_irradiance_spec(fidx:lidx, ix) = 0.0
           
        DO i = 1, nbin
           CALL interpolation (npts, solcomp_wvl(ch, sidx:eidx), solcomp(ch, iix + i, sidx:eidx), &
                npts1, tmpwvl(1:npts1), tmpspec(1:npts1), errstat )
           IF ( errstat > pge_errstat_warning ) THEN
              errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
           END IF
           omi_irradiance_spec(fidx:lidx, ix) = omi_irradiance_spec(fidx:lidx, ix) + tmpspec(1:npts1)
        ENDDO
        omi_irradiance_spec(fidx:lidx, ix) = omi_irradiance_spec(fidx:lidx, ix) / nbin
        omi_irradiance_spec(fidx:lidx, ix) = omi_irradiance_spec(fidx:lidx, ix) * snorms(ch) / omi_solnorm(ix)
     
        fidx = lidx + 1
     ENDDO

     ! Replace ring source spectrum
     IF (omi_solring_ndiv(ix) > 0) THEN  ! two channels
        DO iw = 1, 2
           ch = iw       
           IF (ch == 1) THEN
              fidx = solring_lin(ix); lidx = omi_solring_ndiv(ix)
              IF (omi_solspec_ring(wvl_idx, fidx, ix) <  solcomp_wvl(ch, nscpts(ch))) THEN
                 sidx = MINVAL(MINLOC(omi_solspec_ring(wvl_idx, fidx:lidx, ix), &
                      MASK=( omi_solspec_ring(wvl_idx, fidx:lidx, ix) >= solcomp_wvl(ch, 1:nscpts(ch)) ))) + fidx - 1
                 fidx = sidx; solring_lin(ix) = sidx
              ENDIF
           ELSE
              fidx = omi_solring_ndiv(ix)+1; lidx = solring_uin(ix) 
           ENDIF

           npts1 = lidx - fidx + 1
           tmpwvl(1:npts1) = omi_solspec_ring(wvl_idx, fidx:lidx, ix)
           
           IF ( tmpwvl(1) < solcomp_wvl(ch, 1) .OR. tmpwvl(npts1) > solcomp_wvl(ch, nscpts(ch))) THEN
              WRITE(*, '(2A)') modulename, ': Solar Composite does not cover Ring source spectrum !!!'
              pge_error_status = pge_errstat_error; RETURN
           ENDIF
           
           sidx = MINVAL(MAXLOC(solcomp_wvl(ch, 1:nscpts(ch)), MASK=(solcomp_wvl(ch, 1:nscpts(ch)) &
                <= tmpwvl(1))))
           eidx = MINVAL(MINLOC(solcomp_wvl(ch, 1:nscpts(ch)), MASK=(solcomp_wvl(ch, 1:nscpts(ch)) &
                >= tmpwvl(npts1)))) 
           npts = eidx - sidx + 1

           IF (ch == 1 ) THEN   
              nbin = nxbin
           ELSE
              nbin = nxbin * ncoadd
           ENDIF
           iix = (ix - 1) * nbin

           omi_solspec_ring(spc_idx, fidx:lidx, ix) = 0.0
           DO i = 1, nbin
              CALL interpolation (npts, solcomp_wvl(ch, sidx:eidx), solcomp(ch, iix + i, sidx:eidx), &
                   npts1, tmpwvl(1:npts1), tmpspec(1:npts1), errstat )
              IF ( errstat > pge_errstat_warning ) THEN
                 errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
              END IF
              omi_solspec_ring(spc_idx, fidx:lidx, ix) = omi_solspec_ring(spc_idx, fidx:lidx, ix) + tmpspec(1:npts1) 
           ENDDO
           omi_solspec_ring(spc_idx, fidx:lidx, ix) = omi_solspec_ring(spc_idx, fidx:lidx, ix) &
                / nbin * snorms(ch) / omi_solnorm(ix) 
        ENDDO
     ELSE                                ! one channel, and no coadding
        fidx = solring_lin(ix); lidx = solring_uin(ix)
        npts1 = lidx - fidx + 1
        tmpwvl(1:npts1) = omi_solspec_ring(wvl_idx, fidx:lidx, ix)
        ch = band_selectors(1)

        IF ( tmpwvl(1) < solcomp_wvl(ch, 1) .OR. tmpwvl(npts1) > solcomp_wvl(ch, nscpts(ch))) THEN
           WRITE(*, '(2A)') modulename, ': Solar Composite does not cover Ring source spectrum !!!'
           pge_error_status = pge_errstat_error; RETURN
        ENDIF

        sidx = MINVAL(MAXLOC(solcomp_wvl(ch, 1:nscpts(ch)), MASK=(solcomp_wvl(ch, 1:nscpts(ch)) &
             <= tmpwvl(1))))
        eidx = MINVAL(MINLOC(solcomp_wvl(ch, 1:nscpts(ch)), MASK=(solcomp_wvl(ch, 1:nscpts(ch)) &
             >= tmpwvl(npts1)))) 
        npts = eidx - sidx + 1

        nbin = nxbin; iix = (ix - 1) * nbin
        omi_solspec_ring(spc_idx, fidx:lidx, ix) = 0.0
        DO i = 1, nbin
           CALL interpolation (npts, solcomp_wvl(ch, sidx:eidx), solcomp(ch, iix + i, sidx:eidx), &
                npts1, tmpwvl(1:npts1), tmpspec(1:npts1), errstat )
           IF ( errstat > pge_errstat_warning ) THEN
              errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
           END IF
           omi_solspec_ring(spc_idx, fidx:lidx, ix) = omi_solspec_ring(spc_idx, fidx:lidx, ix) + tmpspec(1:npts1) 
        ENDDO
        omi_solspec_ring(spc_idx, fidx:lidx, ix) = &
             omi_solspec_ring(spc_idx, fidx:lidx, ix) / nbin * snorms(ch) / omi_solnorm(ix) 
     ENDIF
     
     ! Replace solar irradiance around 370 nm region
     ch = 2
     IF ( omi_solspecr(wvl_idx, 1, ix) < solcomp_wvl(ch, 1) .OR. &
          omi_solspecr(wvl_idx, nrefl, ix) > solcomp_wvl(ch, nscpts(ch))) THEN
        WRITE(*, '(2A)') modulename, ': Solar composite does not cover 370 nm!!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
     sidx = MINVAL(MAXLOC(solcomp_wvl(ch, 1:nscpts(ch)), MASK=(solcomp_wvl(ch, 1:nscpts(ch)) &
          <= omi_solspecr(wvl_idx, 1, ix))))
     eidx = MINVAL(MINLOC(solcomp_wvl(ch, 1:nscpts(ch)), MASK=(solcomp_wvl(ch, 1:nscpts(ch)) &
          >= omi_solspecr(wvl_idx, nrefl, ix)))) 
     npts = eidx - sidx + 1

     tmpwvl(1:nrefl) = omi_solspecr(wvl_idx, 1:nrefl, ix)

     IF (.NOT. coadd_uv2) THEN
        nbin = nxbin
     ELSE
        nbin = nxbin * ncoadd
     ENDIF
     iix = (ix - 1) * nbin
     omi_solspecr(spc_idx, 1:nrefl, ix) = 0.0
     
     DO i = 1, nbin
        CALL interpolation (npts, solcomp_wvl(ch, sidx:eidx), solcomp(ch, iix + i, sidx:eidx), &
             nrefl, tmpwvl(1:nrefl), tmpspec(1:nrefl), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
        omi_solspecr(spc_idx, 1:nrefl, ix) = omi_solspecr(spc_idx, 1:nrefl, ix) + tmpspec(1:nrefl)
     ENDDO
     omi_solspecr(spc_idx, 1:nrefl, ix)  = omi_solspecr(spc_idx, 1:nrefl, ix) / nbin * snorms(ch)
  ENDDO
  
  RETURN
END SUBROUTINE replace_solar_irradiance


!! Extrac OMI radiances above high reflectivity clouds 
!SUBROUTINE extract_highcld_radiances (pge_error_status )
!
!  USE OMSAO_precision_module
!  USE OMSAO_parameters_module, ONLY: maxchlen
!  USE OMSAO_variables_module,  ONLY: l1b_rad_filename
!  USE OMSAO_omidata_module,    ONLY: orbc, ntimes_max, nswath, mswath, nxtrack_max, &
!       omi_radiance_swathname, nfxtrack, nwavel_max, nxtrack
!  USE OMSAO_errstat_module
!  USE hdfeos4_parameters
!  USE L1B_Reader_class
!  IMPLICIT NONE
!
!  ! ----------------
!  ! Output variables
!  ! ----------------
!  INTEGER, INTENT (OUT) :: pge_error_status
!
!
!  ! ---------------
!  ! Local variables
!  ! ---------------
!  TYPE (l1b_block_type)                                  :: omi_data_block
!  INTEGER (KIND=i4)                                      :: blockline, errstat, nline, nx, nw
!  INTEGER, PARAMETER                                     :: nbits = 16
!  INTEGER (KIND=i2), DIMENSION(nwavel_max, 0:nbits-1)    :: flgbits
!
!  REAL (KIND=r8)                                         :: tmp_time
!  REAL (KIND=r4)                                         :: tmp_alt, tmp_lat, tmp_lon
!  INTEGER (KIND=i2), DIMENSION(nxtrack_max)              :: tmp_height, tmp_geoflg
!  REAL (KIND=r4), DIMENSION (nwavel_max, nxtrack_max)    :: tmp_rspec, tmp_rprec, tmp_rwavl
!  INTEGER (KIND=I2), DIMENSION (nwavel_max, nxtrack_max) :: tmp_rqflg
! 
!  INTEGER, PARAMETER :: CLDUNIT=133, UV1UNIT=134, UV2UNIT=135
!  INTEGER            :: fline, lline, UVUNIT, is, iloop, iw, ix, i, j, jj
!  INTEGER, DIMENSION(mswath)                             :: nhcld
!  REAL (KIND = dp)                                       :: refl_thresh
!  REAL (KIND = dp), DIMENSION(nxtrack_max, 0:ntimes_max) :: uv2_r331, uv1_r331, r331
!  CHARACTER (LEN=maxchlen)                               :: r331fname, uv1hcld_fname, uv2hcld_fname
!  
!  ! Exteranl functions
!  INTEGER :: estat
!  
!  ! ------------------------------
!  ! Name of this module/subroutine
!  ! ------------------------------
!  CHARACTER (LEN=25), PARAMETER :: modulename = 'extract_highcld_radiances'
! 
!  pge_error_status = pge_errstat_ok;  errstat = omi_s_success
!
!  ! Read reflectivity (r331 from OMTO3)
!  r331fname='/data/dumbo/xliu/OMIHCLD/r331-o' // orbc // '_60S-60N.dat'
!  OPEN(UNIT=CLDUNIT, FILE=TRIM(ADJUSTL(r331fname)), STATUS='OLD', IOSTAT=errstat)
!  READ(CLDUNIT, *) nline, fline, lline
!  READ(CLDUNIT, *) ((uv2_r331(i, j), i = 1, nxtrack_max), j = fline, lline) 
!  CLOSE(CLDUNIT)
!
!  DO i = fline, lline
!     DO j = 1, nfxtrack
!        jj = j * 2 - 1 
!        uv1_r331(j, i) = (uv2_r331(jj, i) + uv2_r331(jj + 1, i) ) / 2.0
!     ENDDO
!  ENDDO
!
!  refl_thresh = 80.0
!  nhcld(1)=COUNT(mask=(uv1_r331(1:nfxtrack, fline:lline) >= refl_thresh))
!  nhcld(2)=COUNT(mask=(uv2_r331(1:nxtrack,  fline:lline) >= refl_thresh))
!  tmp_rspec (1:nwavel_max, 1:nxtrack) = 0.0
!  tmp_rprec (1:nwavel_max, 1:nxtrack) = 0.0
!  tmp_rqflg (1:nwavel_max, 1:nxtrack) = 0.0
!  tmp_rwavl (1:nwavel_max, 1:nxtrack) = 0.0
!
!  uv1hcld_fname='/data/dumbo/xliu/OMIHCLD/OMIL1BUV1-o' // orbc // '_hcld_60S-60N.dat' 
!  uv2hcld_fname='/data/dumbo/xliu/OMIHCLD/OMIL1BUV2-o' // orbc // '_hcld_60S-60N.dat'   
!  OPEN(UNIT=UV1UNIT, FILE=TRIM(ADJUSTL(uv1hcld_fname)), STATUS='UNKNOWN', IOSTAT=errstat) 
!  OPEN(UNIT=UV2UNIT, FILE=TRIM(ADJUSTL(uv2hcld_fname)), STATUS='UNKNOWN', IOSTAT=errstat) 
!  
!   DO is = 1, nswath
!
!     IF (is == 1) THEN
!        UVUNIT = UV1UNIT; r331 = uv1_r331
!     ENDIF
!     IF (is == 2) THEN
!        UVUNIT = UV2UNIT; r331 = uv2_r331
!     ENDIF
!
!     WRITE(UVUNIT, *) nhcld(is)
!
!     ! Open data block called 'omi_data_block' with default size of 100 lines
!     errstat = L1Br_OPEN ( omi_data_block, l1b_rad_filename, omi_radiance_swathname(is), nline )
!     IF ( errstat /= omi_s_success ) THEN
!        estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, 'L1Br_OPEN failed.', modulename, 0 ) ; STOP 1
!     END IF
!
!     errstat = L1Br_getSWdims ( omi_data_block, nXtrack_k=nx, nWavel_k=nw)
!     IF( errstat /=  omi_s_success ) THEN
!        estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, "L1Br_getSWdims failed.", modulename, 0); STOP 1
!     END IF
!
!     print *, nx, nw
!     
!     DO iloop = 0, nline - 1
!        
!        blockline = fline + iloop
!  
!        !errstat = L1Br_getGEOline ( omi_data_block, blockline,          &
!        !     Time_k                    = tmp_time,                      &
!        !     SpacecraftAltitude_k      = tmp_alt,                       &
!        !     SpacecraftLongitude_k     = tmp_lon,                       &
!        !     SpacecraftLatitude_k      = tmp_lat,                       &
!        !     TerrainHeight_k           = tmp_height(1:nx),              &
!        !     GroundPixelQualityFlags_k = tmp_geoflg(1:nx))
!        !
!        !IF ( errstat /= omi_s_success ) THEN
!        !   estat = OMI_SMF_setmsg (omsao_e_read_l1b_file, 'L1Br_getGEOline failed.', modulename, 0); STOP 1
!        !END IF
!                
!        ! Get radiances associated with wavelength range
!        errstat = L1Br_getSIGline ( omi_data_block, blockline,   &
!             Signal_k            = tmp_rspec(1:nw, :),           &
!             SignalPrecision_k   = tmp_rprec(1:nw, :),           &
!             PixelQualityFlags_k = tmp_rqflg(1:nw, :),           &
!             Wavelength_k        = tmp_rwavl(1:nw, :) )     
!       
!        DO ix = 1, nx
!           CALL convert_2bytes_to_16bits ( nbits, nw, tmp_rqflg(1:nw, ix), flgbits(1:nw, 0:nbits-1))
!           
!           IF (r331(ix, blockline) >= refl_thresh) THEN
!              WRITE(UVUNIT, '(3I5,F8.2)') ix, blockline, nw, r331(ix, blockline)
!              
!              DO iw = 1, nw
!                 WRITE(UVUNIT, '(F10.4,D14.6,1X,16I1)') tmp_rwavl(iw, ix), tmp_rspec(iw, ix), flgbits(iw, 0:nbits-1)
!              ENDDO
!           ENDIF
!        ENDDO
!        
!        IF( errstat /= omi_s_success ) THEN
!           estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, "L1Br_getSIGline failed.", modulename, 0); STOP 1
!        ENDIF
!     END DO     ! end iloop
!
!     CLOSE (UVUNIT)
!     
!  ENDDO ! end swath loop
!
!END SUBROUTINE extract_highcld_radiances
!
!SUBROUTINE extract_highcld_visradiances (pge_error_status )
!
!  USE OMSAO_precision_module
!  USE OMSAO_parameters_module, ONLY: maxchlen
!  USE OMSAO_variables_module,  ONLY: l1b_rad_filename
!  USE OMSAO_omidata_module,    ONLY: orbc, ntimes_max, nswath, mswath, nxtrack_max, &
!       omi_radiance_swathname, nfxtrack, nwavel_max, nxtrack
!  USE OMSAO_errstat_module
!  USE hdfeos4_parameters
!  USE L1B_Reader_class
!  IMPLICIT NONE
!
!  ! ----------------
!  ! Output variables
!  ! ----------------
!  INTEGER, INTENT (OUT) :: pge_error_status
!
!
!  ! ---------------
!  ! Local variables
!  ! ---------------
!  TYPE (l1b_block_type)                                  :: omi_data_block
!  INTEGER (KIND=i4)                                      :: blockline, errstat, nline, nx, nw
!  INTEGER, PARAMETER                                     :: nbits = 16
!  INTEGER (KIND=i2), DIMENSION(nwavel_max, 0:nbits-1)    :: flgbits
!
!  REAL (KIND=r8)                                         :: tmp_time
!  REAL (KIND=r4)                                         :: tmp_alt, tmp_lat, tmp_lon
!  INTEGER (KIND=i2), DIMENSION(nxtrack_max)              :: tmp_height, tmp_geoflg
!  REAL (KIND=r4), DIMENSION (nwavel_max, nxtrack_max)    :: tmp_rspec, tmp_rprec, tmp_rwavl
!  INTEGER (KIND=I2), DIMENSION (nwavel_max, nxtrack_max) :: tmp_rqflg
! 
!  INTEGER, PARAMETER :: CLDUNIT=133, VISUNIT=134, ORBUNIT=135
!  INTEGER            :: fline, lline, UVUNIT, is, iloop, iw, ix, i, j, jj, nhcld
!  REAL (KIND = dp)                                       :: refl_thresh
!  REAL (KIND = dp), DIMENSION(nxtrack_max, 0:ntimes_max) :: r331
!  CHARACTER (LEN=maxchlen)                               :: r331fname, hcld_fname
!  
!  ! Exteranl functions
!  INTEGER :: estat
!  
!  ! ------------------------------
!  ! Name of this module/subroutine
!  ! ------------------------------
!  CHARACTER (LEN=28), PARAMETER :: modulename = 'extract_highcld_visradiances'
! 
!  pge_error_status = pge_errstat_ok;  errstat = omi_s_success
!
!  OPEN(UNIT=ORBUNIT, FILE='temp_rad.dat', STATUS='OLD', IOSTAT=errstat)
!  DO is = 1, 55
!     READ(ORBUNIT, '(A)') l1b_rad_filename
!     i = INDEX(l1b_rad_filename, '-o') + 2
!     orbc = l1b_rad_filename(i : i + 5)
!
!     ! Read reflectivity (r331 from OMTO3)
!     r331fname='/data/dumbo/xliu/OMIHCLD/r331-o' // orbc // '_60S-60N.dat'
!     OPEN(UNIT=CLDUNIT, FILE=TRIM(ADJUSTL(r331fname)), STATUS='OLD', IOSTAT=errstat)
!     READ(CLDUNIT, *) nline, fline, lline
!     READ(CLDUNIT, *) ((r331(i, j), i = 1, nxtrack_max), j = fline, lline) 
!     CLOSE(CLDUNIT)
!     
!     refl_thresh = 80.0
!     nhcld=COUNT(mask=(r331(1:nxtrack,  fline:lline) >= refl_thresh))
!     print *, nhcld
!     
!     tmp_rspec (1:nwavel_max, 1:nxtrack) = 0.0
!     tmp_rprec (1:nwavel_max, 1:nxtrack) = 0.0
!     tmp_rqflg (1:nwavel_max, 1:nxtrack) = 0.0
!     tmp_rwavl (1:nwavel_max, 1:nxtrack) = 0.0
!     
!     hcld_fname='/data/dumbo/xliu/OMIHCLD/OMIL1BVIS-o' // orbc // '_hcld_60S-60N.dat' 
!     OPEN(UNIT=VISUNIT, FILE=TRIM(ADJUSTL(hcld_fname)), STATUS='UNKNOWN', IOSTAT=errstat) 
!     print *, TRIM(ADJUSTL(r331fname)), TRIM(ADJUSTL(hcld_fname))
!     
!     WRITE(VISUNIT, *) nhcld
!     
!     ! Open data block called 'omi_data_block' with default size of 100 lines
!     
!     errstat = L1Br_OPEN ( omi_data_block, l1b_rad_filename, 'Earth VIS Swath ', nline )
!     IF ( errstat /= omi_s_success ) THEN
!        estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, 'L1Br_OPEN failed.', modulename, 0 ) ; STOP 1
!     END IF
!     
!     errstat = L1Br_getSWdims ( omi_data_block, nXtrack_k=nx, nWavel_k=nw)
!     IF( errstat /=  omi_s_success ) THEN
!        estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, "L1Br_getSWdims failed.", modulename, 0); STOP 1
!     END IF
!     
!     DO iloop = 0, nline - 1
!        
!        blockline = fline + iloop
!        
!        !errstat = L1Br_getGEOline ( omi_data_block, blockline,          &
!        !     Time_k                    = tmp_time,                      &
!        !     SpacecraftAltitude_k      = tmp_alt,                       &
!        !     SpacecraftLongitude_k     = tmp_lon,                       &
!        !     SpacecraftLatitude_k      = tmp_lat,                       &
!        !     TerrainHeight_k           = tmp_height(1:nx),              &
!        !     GroundPixelQualityFlags_k = tmp_geoflg(1:nx))
!        !
!        !IF ( errstat /= omi_s_success ) THEN
!        !   estat = OMI_SMF_setmsg (omsao_e_read_l1b_file, 'L1Br_getGEOline failed.', modulename, 0); STOP 1
!        !END IF
!        
!        ! Get radiances associated with wavelength range
!        errstat = L1Br_getSIGline ( omi_data_block, blockline,   &
!             Signal_k            = tmp_rspec(1:nw, :),           &
!             SignalPrecision_k   = tmp_rprec(1:nw, :),           &
!             PixelQualityFlags_k = tmp_rqflg(1:nw, :),           &
!             Wavelength_k        = tmp_rwavl(1:nw, :) )     
!        
!        DO ix = 1, nx
!           CALL convert_2bytes_to_16bits ( nbits, nw, tmp_rqflg(1:nw, ix), flgbits(1:nw, 0:nbits-1))
!           
!           IF (r331(ix, blockline) >= refl_thresh) THEN
!              WRITE(VISUNIT, '(3I5,F8.2)') ix, blockline, nw, r331(ix, blockline)
!              
!              DO iw = 1, nw
!                 WRITE(VISUNIT, '(F10.4,D14.6,1X,16I1)') tmp_rwavl(iw, ix), tmp_rspec(iw, ix), flgbits(iw, 0:nbits-1)
!              ENDDO
!           ENDIF
!        ENDDO
!        
!        IF( errstat /= omi_s_success ) THEN
!           estat = OMI_SMF_setmsg ( omsao_e_read_l1b_file, "L1Br_getSIGline failed.", modulename, 0); STOP 1
!        ENDIF
!     END DO     ! end iloop
!     
!     CLOSE (VISUNIT)
!  ENDDO
!  CLOSE(ORBUNIT)
!     
!END SUBROUTINE extract_highcld_visradiances

SUBROUTINE reduce_irrad_resolution (spec, qflag, np, nx, which_slit, slwth, &
     samprate, dwav, retswav, retewav, swav, ewav, np_out, pge_error_status)

  USE OMSAO_precision_module 
  USE OMSAO_parameters_module, ONLY: max_spec_pts
  USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx
  USE OMSAO_variables_module,  ONLY: use_redfixwav, nredfixwav, redfixwav
  USE ozprof_data_module,      ONLY: pos_alb
  USE OMSAO_errstat_module
  
  IMPLICIT NONE
  
  ! ====================
  ! In/Output variables
  ! ====================
  INTEGER, INTENT(IN)                                 :: np, nx, which_slit
  INTEGER, INTENT(OUT)                                :: np_out, pge_error_status
  INTEGER (KIND=i2), DIMENSION(np, nx), INTENT(IN)    :: qflag                   
  REAL (KIND=dp), INTENT(IN)                          :: samprate, slwth, retswav, retewav
  REAL (KIND=dp), INTENT(INOUT)                       :: dwav
  REAL (KIND=dp), INTENT(INOUT)                       :: swav, ewav
  REAL (KIND=dp), DIMENSION(sig_idx, np, nx), INTENT(INOUT) :: spec
  
  ! ====================
  ! Local variables
  ! ====================
  INTEGER, PARAMETER     :: nmax = max_spec_pts
  INTEGER                :: i, j, ix, mslit, nf, nsamp, nsamp1, nslit, errstat, iwin, idx
  INTEGER, DIMENSION(nx) :: nmod
  INTEGER, DIMENSION( 3) :: sidx, eidx, nstep
  REAL (KIND=dp)         :: dlam0, slitsum, redsnr, dx, fwav, lwav
  REAL (KIND=dp), DIMENSION(nmax)      :: slit
  REAL (KIND=dp), DIMENSION(sig_idx, np, nx) :: specmod
  REAL (KIND=dp), DIMENSION(sig_idx, nmax)   :: fnspec

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=28), PARAMETER :: modulename = 'reduce_irradiance_resolution'

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg
 
  pge_error_status = pge_errstat_ok

  ! If slwth == 0, not allowed, unless that user provides a fixed wavelength grid
  IF (slwth == 0 .AND. .NOT. use_redfixwav) THEN
     WRITE(*, '(2A)') modulename, ': Zero Slit Width Not allowed!!!'
     pge_error_status = pge_errstat_error
  ENDIF

  ! Filter bad measurements, determine wavelength range for all x-track positions
  dlam0 = spec(1, 2, 1) - spec(1, 1, 1)
  ewav = MAXVAL(spec(1, np, :)) 
  DO ix = 1, nx
     j = 0
     
     DO i = 1, np 
        IF (spec(2, i, ix) > 0. .AND. spec(2, i, ix) <= 4.0E14 .AND. qflag(i, ix) == 0) THEN
           j = j + 1; specmod(:, j, ix) = spec(:, i, ix)
        ENDIF
     ENDDO
     nmod(ix) = j
     IF (specmod(1, j, ix) < ewav) ewav = specmod(1, j, ix)
  ENDDO
  swav = MAXVAL(specmod(1, 1, :))

  IF (slwth == 0 .AND. use_redfixwav) THEN

     j = 0
     DO i = 1, nredfixwav
        IF (redfixwav(i) >= swav .AND. redfixwav(i) <=ewav) THEN
           j = j + 1
           spec(1, j, 1:nx) = redfixwav(i)
        ENDIF
     ENDDO

  ELSE

     ! Establish fine wavelength scale (common to all across-track positions)
     nf = INT((ewav - swav) / dwav + 1)
     fnspec(1, 1) = swav
     DO i = 2, nf
        fnspec(1, i) = fnspec(1, i - 1) + dwav
     ENDDO
     
     IF (dwav > dlam0 .OR. dwav <= 0) THEN
        dwav = dlam0
     ENDIF
     
     IF (which_slit == 1) THEN               ! Symmetric Gaussian (slwth is hw1e)
        ! hw1e * 1.6551 = FWHM
        mslit = NINT(2.62826 * slwth / dwav) ! slit truncation (<0.1%)
        nsamp = INT(slwth * 1.66511 / samprate / dwav)
     ELSE IF (which_slit == 2) THEN          ! Triangular (slwth is FWHM)
        mslit = NINT(slwth / dwav)
        nsamp = INT(slwth / samprate / dwav)
     ENDIF
     nslit = mslit * 2 + 1
     nsamp1 = INT(dlam0 / dwav);  IF (nsamp1 < 1.0) nsamp1 = 1
     
     ! Set up slit function
     IF (which_slit == 1) THEN
        slit(1:nslit) = EXP( -((fnspec(1, 1:nslit)-fnspec(1, mslit+1)) / slwth)**2 )
     ELSE
        slit(1:nslit) = 1.0 - ABS(fnspec(1, 1:nslit)-fnspec(1, mslit+1)) / slwth
        WHERE (slit(1:nslit) < 0.0)
           slit(1:nslit) = 0.0
        ENDWHERE
     ENDIF
     slitsum = SUM(slit(1:nslit)); slit(1:nslit) = slit(1:nslit) / slitsum       ! Normalization
     redsnr = 1.0 / SQRT(slitsum / dlam0 * dwav)   ! Measurement error/noise reduction
     
     ! Convolve and Sample
     ! More sampling at the two ends of the selected spectral range
     ! Especially for the first and last 4 positions 
     ! This is to avoid extrapolation while keeping as many measurements as possible
     IF (retswav <= fnspec(1, mslit+1)) THEN
        i = mslit + 1 + 3 * nsamp1
     ELSE
        i = MAXVAL(MINLOC(fnspec(1, 1:nf), MASK=(fnspec(1, 1:nf) >= retswav + dlam0 * 3)))
     ENDIF
     sidx(1) = mslit+1; eidx(1) = i; nstep(1) = nsamp1
     
     IF (retewav >= fnspec(1, nf - mslit)) THEN
        i = nf - mslit - nsamp1 * 3
     ELSE
        i = MAXVAL(MAXLOC(fnspec(1, 1:nf), MASK=(fnspec(1, 1:nf) <= retewav - dlam0 * 3)))
     ENDIF
     sidx(2) = eidx(1) + nsamp1; eidx(2) = i-1; nstep(2) = nsamp
     sidx(3) = i + nsamp1; eidx(3) = nf - mslit; nstep(3) = nsamp1

     j = 0
     DO i = 1, nredfixwav
        IF (redfixwav(i) > fnspec(1, mslit + 1) .AND. redfixwav(i) < fnspec(1, nf - mslit)) THEN
           j = j + 1
           spec(1, j, 1:nx) = redfixwav(i)
        ENDIF
     ENDDO
  ENDIF

  IF (use_redfixwav) THEN
     ! Add 3 * 2 extra wavelengths for irradiance and 2 * 2 extra wavelengths for radiance
     ! Add 3 wavelengths at the beginning of a spectra region
     DO i = 1, j
        fwav = MAX(swav, retswav)
        IF (spec(1, i, 1) > fwav) THEN
           IF (i == 1) THEN
              dx = (spec(1, i, 1) - fwav) / 3.
           ELSE
              dx = (spec(1, i, 1) - MAX(spec(1, i-1, 1), fwav)) / 3.
           ENDIF
           IF (dx > dlam0) dx = dlam0
           spec(1, i+3:j+3, 1:nx) = spec(1, i:j, 1:nx)
           spec(1, i, 1:nx)   = spec(1, i+3, 1:nx) - dx * 3.0
           spec(1, i+1, 1:nx)   = spec(1, i+3, 1:nx) - dx * 2.0
           spec(1, i+2, 1:nx)   = spec(1, i+3, 1:nx) - dx * 1.0
           j = j + 3
           EXIT
        ENDIF
     ENDDO
     
     ! Add 3 wavelengths at the end of a spectra region
     DO i = j, 1, -1
        lwav = MIN(ewav, retewav)
        IF (spec(1, i, 1) < lwav) THEN
           IF (i == j) THEN
              dx = (lwav - spec(1, i, 1)) / 3.
           ELSE
              dx = ( MAX(spec(1, i+1, 1), lwav) - spec(1, i, 1)) / 3.
           ENDIF
           IF (dx > dlam0) dx = dlam0
           spec(1, i+4:j+3, 1:nx) = spec(1, i+1:j, 1:nx)
           spec(1, i+1, 1:nx)   = spec(1, i, 1:nx) + dx * 1.0
           spec(1, i+2, 1:nx)   = spec(1, i, 1:nx) + dx * 2.0
           spec(1, i+3, 1:nx)   = spec(1, i, 1:nx) + dx * 3.0
           j = j + 3
           EXIT
        ENDIF
     ENDDO

     ! Add a wavelength at the wavelength to derive initial cloud fraction
     IF (swav < pos_alb .AND. ewav > pos_alb) THEN
        DO i = j, 1, -1 
           IF (spec(1, i, 1) < pos_alb) THEN
              spec(1, i+1:j, 1:nx) = spec(1, i+2:j+1, 1:nx)
              spec(1, i+1, 1:nx) = pos_alb
              j = j + 1
              EXIT
           ENDIF
        ENDDO
     ENDIF

  ENDIF
  np_out = j

  ! Perform direct interpolation
  IF (slwth == 0 .AND. use_redfixwav) THEN

     DO ix = 1, nx
        DO i = 2, 3
           CALL interpolation (nmod(ix),  specmod(1, 1:nmod(ix), ix), specmod(i, 1:nmod(ix), ix), &
                np_out, spec(1, 1:np_out, ix),  spec(i, 1:np_out, ix), pge_error_status )
           IF ( pge_error_status > pge_errstat_warning ) THEN
              errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; RETURN
           END IF
        ENDDO
     ENDDO
  ELSE
     
     DO ix = 1, nx 
        ! Pre-interpolation
        DO i = 2, 3
           CALL interpolation (nmod(ix),  specmod(1, 1:nmod(ix), ix), specmod(i, 1:nmod(ix), ix), &
                nf, fnspec(1, 1:nf), fnspec(i, 1:nf), pge_error_status )
           IF ( pge_error_status > pge_errstat_warning ) THEN
              print *, ix, 'interpolation problem'
              errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; RETURN
           END IF
        ENDDO
        
        IF (.NOT. use_redfixwav) THEN
           j = 0
           DO iwin = 1, 3
              DO i = sidx(iwin), eidx(iwin), nstep(iwin)
                 j = j + 1
                 spec(1, j, ix) = fnspec(1, i)
                 spec(2, j, ix) = SUM(slit(1:nslit) * fnspec(2, i-mslit:i+mslit))
                 spec(3, j, ix) = SUM(slit(1:nslit) * fnspec(3, i-mslit:i+mslit)) * redsnr ! Reduce noise
              ENDDO
           ENDDO
        ELSE
           DO j = 1,  np_out
              idx = MAXVAL(MINLOC(fnspec(1, 1:nf), MASK=(fnspec(1, 1:nf) >= spec(1, j, ix))))
              IF ( ABS(fnspec(1, idx-1) - spec(1, j, ix)) < ABS(fnspec(1, idx) - spec(1, j, ix))) idx = idx - 1
              spec(2, j, ix) = SUM(slit(1:nslit) * fnspec(2, idx-mslit:idx+mslit))
              spec(3, j, ix) = SUM(slit(1:nslit) * fnspec(3, idx-mslit:idx+mslit)) * redsnr   ! Reduce noise              
           ENDDO
        ENDIF
     ENDDO

     np_out = j
  ENDIF
  !print *, nx
  !IF (nx == 30) THEN
  !   WRITE(*, '(2F10.4)') spec(1, np_out, 15)
  !   WRITE(*, '(2D14.6)') spec(2, np_out, 15)
  !ELSE
  !   WRITE(*, '(2F10.4)') spec(1, np_out, 29:30)
  !   WRITE(*, '(2D14.6)') SUM(spec(2, np_out, 29:30)) / 2.
  !ENDIF
  !IF (np_out > np) THEN
  !   WRITE(*, '(2A)') modulename, ': Improper sampling rate or slit width!!!'
  !   pge_error_status = pge_errstat_error
  !ENDIF

  RETURN
END SUBROUTINE reduce_irrad_resolution


SUBROUTINE reduce_rad_resolution (spec, qflag, np, nx, which_slit, slwth, &
     samprate, dwav, retswav, retewav, swav, ewav, np_out, pge_error_status)

  USE OMSAO_precision_module 
  USE OMSAO_parameters_module, ONLY: max_spec_pts
  USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx
  USE OMSAO_variables_module,  ONLY: use_redfixwav, nredfixwav, redfixwav
  USE ozprof_data_module,      ONLY: pos_alb
  USE OMSAO_errstat_module
  
  IMPLICIT NONE
  
  ! ====================
  ! In/Output variables
  ! ====================
  INTEGER, INTENT(IN)                              :: np, nx, which_slit
  INTEGER, INTENT(OUT)                             :: np_out, pge_error_status
  INTEGER (KIND=i2), DIMENSION(np, nx), INTENT(IN) :: qflag                   
  REAL (KIND=dp), INTENT(IN)                       :: samprate, slwth, dwav, swav, ewav, retswav, retewav
  REAL (KIND=dp), DIMENSION(sig_idx, np, nx), INTENT(INOUT) :: spec
  
  ! ====================
  ! Local variables
  ! ====================
  INTEGER, PARAMETER     :: nmax = max_spec_pts
  INTEGER                :: i, j, ix, mslit, nf, nsamp, nsamp1, nslit, errstat, fidx, lidx, iwin, idx
  INTEGER, DIMENSION(nx) :: nmod
  INTEGER, DIMENSION( 3) :: sidx, eidx, nstep
  REAL (KIND=dp)         :: dlam0, slitsum, redsnr, dx, fwav, lwav
  REAL (KIND=dp), DIMENSION(nmax)      :: slit
  REAL (KIND=dp), DIMENSION(sig_idx, np, nx) :: specmod
  REAL (KIND=dp), DIMENSION(sig_idx, nmax)   :: fnspec

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=28), PARAMETER :: modulename = 'reduce_irradiance_resolution'

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  pge_error_status = pge_errstat_ok

  ! If slwth == 0, not allowed, unless that user provides a fixed wavelength grid
  IF (slwth == 0 .AND. .NOT. use_redfixwav) THEN
     WRITE(*, '(2A)') modulename, ': Zero Slit Width Not allowed!!!'
     pge_error_status = pge_errstat_error
  ENDIF

  dlam0 = spec(1, 2, 1) - spec(1, 1, 1)  

  ! Filter bad measurements
  DO ix = 1, nx
     j = 0
     
     DO i = 1, np 
        IF (spec(2, i, ix) > 0. .AND. spec(2, i, ix) <= 4.0E14 .AND. qflag(i, ix) == 0) THEN
           j = j + 1; specmod(:, j, ix) = spec(:, i, ix)
        ENDIF
     ENDDO
     nmod(ix) = j
  ENDDO

  IF (slwth == 0 .AND. use_redfixwav) THEN

     j = 0
     DO i = 1, nredfixwav
        IF (redfixwav(i) >= swav .AND. redfixwav(i) <=ewav) THEN
           j = j + 1
           spec(1, j, 1:nx) = redfixwav(i)
        ENDIF
     ENDDO

  ELSE

     ! Establish fine wavelength scale (common to all across-track positions)
     nf = INT((ewav - swav) / dwav + 1)
     fnspec(1, 1) = swav
     DO i = 2, nf
        fnspec(1, i) = fnspec(1, i - 1) + dwav
     ENDDO
     
     IF (which_slit == 1) THEN               ! Symmetric Gaussian (slwth is hw1e)
        ! hw1e * 1.6551 = FWHM
        mslit = NINT(2.62826 * slwth / dwav) ! slit truncation (<0.1%)
        nsamp = INT(slwth * 1.66511 / samprate / dwav)
     ELSE IF (which_slit == 2) THEN          ! Triangular (slwth is FWHM)
        mslit = NINT(slwth / dwav)
        nsamp = INT(slwth / samprate / dwav)
     ENDIF
     nslit  = mslit * 2 + 1
     nsamp1 = INT(dlam0 / dwav);  IF (nsamp1 < 1.0) nsamp1 = 1
     
     ! Set up slit function
     IF (which_slit == 1) THEN
        slit(1:nslit) = EXP( -((fnspec(1, 1:nslit)-fnspec(1, mslit+1)) / slwth)**2 )
     ELSE
        slit(1:nslit) = 1.0 - ABS(fnspec(1, 1:nslit)-fnspec(1, mslit+1)) / slwth
        WHERE (slit(1:nslit) < 0.0)
           slit(1:nslit) = 0.0
        ENDWHERE
     ENDIF
     slitsum = SUM(slit(1:nslit)); slit(1:nslit) = slit(1:nslit) / slitsum ! Normalization
     redsnr  = 1.0 / SQRT(slitsum / dlam0 * dwav)                          ! Measurement error/noise reduction

     ! Convolve and Sample
     ! More sampling at the two ends of the selected spectral range
     ! Especially for the first and last 4 positions 
     ! This is to avoid extrapolation while keeping as many measurements as possible
     IF (retswav <= fnspec(1, mslit+1)) THEN
        i = mslit + 1 + 3 * nsamp1
     ELSE
        i = MAXVAL(MINLOC(fnspec(1, 1:nf), MASK=(fnspec(1, 1:nf) >= retswav + dlam0 * 3)))
     ENDIF
     sidx(1) = mslit+1; eidx(1) = i; nstep(1) = nsamp1
     
     IF (retewav >= fnspec(1, nf - mslit)) THEN
        i = nf - mslit - nsamp1 * 3
     ELSE
        i = MAXVAL(MAXLOC(fnspec(1, 1:nf), MASK=(fnspec(1, 1:nf) <= retewav - dlam0 * 3)))
     ENDIF
     sidx(2) = eidx(1) + nsamp1; eidx(2) = i-1; nstep(2) = nsamp
     sidx(3) = i + nsamp1; eidx(3) = nf - mslit; nstep(3) = nsamp1
     
     j = 0
     DO i = 1, nredfixwav
        IF (redfixwav(i) > fnspec(1, mslit + 1) .AND. redfixwav(i) < fnspec(1, nf - mslit)) THEN
           j = j + 1
           spec(1, j, 1:nx) = redfixwav(i)
        ENDIF
     ENDDO

  ENDIF

  IF (use_redfixwav) THEN
     ! Add 3 * 2 extra wavelengths for irradiance and 2 * 2 extra wavelengths for radiance
     ! Add 3 wavelengths at the beginning of a spectra region
     DO i = 1, j
        fwav = MAX(swav, retswav)
        IF (spec(1, i, 1) > fwav) THEN
           IF (i == 1) THEN
              dx = (spec(1, i, 1) - fwav) / 3.
           ELSE
              dx = (spec(1, i, 1) - MAX(spec(1, i-1, 1), fwav)) / 3.
           ENDIF
           IF (dx > dlam0) dx = dlam0
           spec(1, i+3:j+3, 1:nx) = spec(1, i:j, 1:nx)
           spec(1, i, 1:nx)   = spec(1, i+3, 1:nx) - dx * 3.0
           spec(1, i+1, 1:nx)   = spec(1, i+3, 1:nx) - dx * 2.0
           spec(1, i+2, 1:nx)   = spec(1, i+3, 1:nx) - dx * 1.0
           j = j + 3
           EXIT
        ENDIF
     ENDDO
     
     ! Add 3 wavelengths at the end of a spectra region
     DO i = j, 1, -1
        lwav = MIN(ewav, retewav)
        IF (spec(1, i, 1) < lwav) THEN
           IF (i == j) THEN
              dx = (lwav - spec(1, i, 1)) / 3.
           ELSE
              dx = ( MAX(spec(1, i+1, 1), lwav) - spec(1, i, 1)) / 3.
           ENDIF
           IF (dx > dlam0) dx = dlam0
           spec(1, i+4:j+3, 1:nx) = spec(1, i+1:j, 1:nx)
           spec(1, i+1, 1:nx)   = spec(1, i, 1:nx) + dx * 1.0
           spec(1, i+2, 1:nx)   = spec(1, i, 1:nx) + dx * 2.0
           spec(1, i+3, 1:nx)   = spec(1, i, 1:nx) + dx * 3.0
           j = j + 3
           EXIT
        ENDIF
     ENDDO

     ! Add a wavelength at the wavelength to derive initial cloud fraction
     IF (swav < pos_alb .AND. ewav > pos_alb) THEN
        DO i = j, 1, -1 
           IF (spec(1, i, 1) < pos_alb) THEN
              spec(1, i+1:j, 1:nx) = spec(1, i+2:j+1, 1:nx)
              spec(1, i+1, 1:nx) = pos_alb
              j = j + 1
              EXIT
           ENDIF
        ENDDO
     ENDIF

  ENDIF
  np_out = j
  
  IF (slwth == 0 .AND. use_redfixwav) THEN

     DO ix = 1, nx
        DO i = 2, 3
           CALL interpolation (nmod(ix),  specmod(1, 1:nmod(ix), ix), specmod(i, 1:nmod(ix), ix), &
                np_out, spec(1, 1:np_out, ix),  spec(i, 1:np_out, ix), pge_error_status )
           IF ( pge_error_status > pge_errstat_warning ) THEN
              errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; RETURN
           END IF
        ENDDO
     ENDDO

  ELSE
     DO ix = 1, nx 
        ! Pre-interpolation
        fidx = MAXVAL(MINLOC(fnspec(1, 1:nf), MASK=(fnspec(1, 1:nf) >= specmod(1, 1, ix))))
        lidx = MAXVAL(MAXLOC(fnspec(1, 1:nf), MASK=(fnspec(1, 1:nf) <= specmod(1, nmod(ix), ix))))
        fnspec(2:3, 1:fidx-1) = 0.0; fnspec(2:3, lidx+1:nf) = 0.0
        
        IF (nmod(ix) > np * 0.75) THEN
           DO i = 2, 3           
              CALL interpolation (nmod(ix),  specmod(1, 1:nmod(ix), ix), specmod(i, 1:nmod(ix), ix), &
                   lidx-fidx+1, fnspec(1, fidx:lidx), fnspec(i, fidx:lidx), pge_error_status )
              IF ( pge_error_status > pge_errstat_warning ) THEN
                 errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; RETURN
              END IF
           ENDDO
        ELSE
           fnspec(2:3, 1:nf) = 0.0
        ENDIF
        
        IF (.NOT. use_redfixwav) THEN
           ! Convolve and Sample
           j = 0
           DO iwin = 1, 3
              DO i = sidx(iwin), eidx(iwin), nstep(iwin)
                 j = j + 1
                 
                 spec(1, j, ix) = fnspec(1, i)
                 spec(2, j, ix) = SUM(slit(1:nslit) * fnspec(2, i-mslit:i+mslit))
                 spec(3, j, ix) = SUM(slit(1:nslit) * fnspec(3, i-mslit:i+mslit)) * redsnr ! Reduce noise
                 
                 IF ((i - mslit < fidx) .OR. (i + mslit > lidx)) THEN
                    spec(2:3, j, ix) = 0.0
                 ENDIF
              ENDDO
           ENDDO
        ELSE
           DO j = 1,  np_out
              idx = MAXVAL(MINLOC(fnspec(1, 1:nf), MASK=(fnspec(1, 1:nf) >= spec(1, j, ix))))
              IF ( ABS(fnspec(1, idx-1) - spec(1, j, ix)) < ABS(fnspec(1, idx) - spec(1, j, ix))) idx = idx - 1
              spec(2, j, ix) = SUM(slit(1:nslit) * fnspec(2, idx-mslit:idx+mslit))
              spec(3, j, ix) = SUM(slit(1:nslit) * fnspec(3, idx-mslit:idx+mslit)) * redsnr   ! Reduce noise              
           ENDDO
        ENDIF
     ENDDO

     np_out = j
  ENDIF
  
  !print *, np_out, retswav, retewav
  !print *, nx
  !IF (nx == 1) THEN
  !   WRITE(*, '(2F10.4)') spec(1, np_out, 1)
  !   WRITE(*, '(2D14.6)') spec(2, np_out, 1)
  !ELSE
  !   WRITE(*, '(2F10.4)') spec(1, np_out, 1:2)
  !   WRITE(*, '(2D14.6)') SUM(spec(2, np_out, 1:2)) / 2.
  !ENDIF
  !IF (np_out > np) THEN
  !   WRITE(*, '(2A)') modulename, ': Improper sampling rate or slit width!!!'
  !   pge_error_status = pge_errstat_error
  !ENDIF
  
  RETURN
END SUBROUTINE reduce_rad_resolution
