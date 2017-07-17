!xliu, 03/25/2011, add the XtrackQualityFlags variable for row anomaly
MODULE OMSAO_pixelcorner_module

  ! ============================================================================== !
  ! 1. Compute geolocation of OMI pixel corners from latitude and longitude fields !
  ! 2. Compute one single effective viewing geometry for each pixel since current  !
  !    current viewing geometry is only provided at pixel center                   !
  ! 3. Deal with binning along and across the track                                !
  ! ============================================================================== !

  USE OMSAO_precision_module
  USE OMSAO_errstat_module
  USE L1B_Reader_class
  USE OMSAO_variables_module, ONLY: l1b_rad_filename, verb_thresh_lev, the_utc, &
       TAI93At0ZOfGranule, TAI93StartOfGranule, GranuleYear, GranuleMonth,      &
       GranuleDay, GranuleHour, GranuleMinute, GranuleSecond, GranuleJDay,here_stop
  USE OMSAO_omidata_module,   ONLY: omi_radiance_swathname, nxtrack_max, ntimes_max, &
       nxbin, nybin, offset_line, zoom_mode, nswath, zoom_p1
  IMPLICIT NONE

  ! Public parameters
  REAL (KIND=r8), DIMENSION (1:nxtrack_max, 0:ntimes_max-1) :: omi_alllat,  omi_alllon
  REAL (KIND=r8), DIMENSION (0:nxtrack_max, 0:ntimes_max-1) :: omi_allelat, omi_allelon
  REAL (KIND=r8), DIMENSION (1:nxtrack_max, 0:ntimes_max-1) :: omi_allsza,  omi_allvza, &
       omi_allaza, omi_allsca
  REAL (KIND=r8), DIMENSION (0:nxtrack_max, 0:ntimes_max)   :: omi_allclat, omi_allclon
  REAL (KIND=r8), DIMENSION (0:ntimes_max-1)                :: omi_alltime
  REAL (KIND=r4), DIMENSION (0:ntimes_max-1)                :: omi_allSpcftLat,        &
       omi_allSpcftLon, omi_allSpcftAlt, omi_allSecondsInDay
  INTEGER (KIND=i2), DIMENSION(1:nxtrack_max,0:ntimes_max-1):: omi_allGeoFlg, omi_allHeight
  INTEGER (KIND=i1), DIMENSION(1:nxtrack_max,0:ntimes_max-1):: omi_allXTrackQFlg
  INTEGER (KIND=i2), DIMENSION(0:ntimes_max-1)              :: omi_allMflg
  INTEGER (KIND=i1), DIMENSION(0:ntimes_max-1)              :: omi_allNSPC

  ! ----------------
  ! Local Parameters
  ! ---------------------------------------------------------------------
  ! * Values for Pi (rad, deg) and Conversions between Degree and Radians
  ! ---------------------------------------------------------------------
  REAL (KIND=r8), PARAMETER, PRIVATE :: pi         = 3.14159265358979_r8  ! 2*ASIN(1.0_r8)
  REAL (KIND=r8), PARAMETER, PRIVATE :: pihalf     = 0.5_r8  * pi
  REAL (KIND=r8), PARAMETER, PRIVATE :: twopi      = 2.0_r8  * pi
  REAL (KIND=r8), PARAMETER, PRIVATE :: pi_deg     = 180.0_r8
  REAL (KIND=r8), PARAMETER, PRIVATE :: pihalf_deg =  90.0_r8
  REAL (KIND=r8), PARAMETER, PRIVATE :: twopi_deg  = 360.0_r8
  REAL (KIND=r8), PARAMETER, PRIVATE :: deg2rad    = pi / 180.0_r8
  REAL (KIND=r8), PARAMETER, PRIVATE :: rad2deg    = 180.0_r8 / pi
  ! ---------------------------------------------------------------------
  ! * Precison for DEG <-> RAD conversion - anything less than EPS
  !   is effectively ZERO.
  ! ---------------------------------------------------------------------
  REAL (KIND=r8), PARAMETER, PRIVATE :: eps = 1.0E-10_r8
  ! ---------------------------------------------------------------------
  REAL (KIND=r4), PARAMETER, PRIVATE :: rearth0 = 6378  ! equatorial radius
  REAL (KIND=r4), PARAMETER, PRIVATE :: minza = 0.0, maxza=90.0, minaza = -360., maxaza = 360.0


CONTAINS

  SUBROUTINE compute_pixel_corners ( ntimes, nxtrack, nl, pge_error_status )

    ! =======================================================
    ! Computes OMI pixel corner coordinates, start to finish:
    !
    !  * Reads L1b geolocation data
    !  * Computes the corners
    !  * Writes output to file
    ! =======================================================

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: ntimes, nxtrack
    INTEGER, INTENT (IN)           :: nl
    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (OUT) :: pge_error_status

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4), PARAMETER                      :: ntloop = 100
    INTEGER (KIND=i4)  :: errstat, iline, i, j, nline, nx, ny, sline, eline, ix, &
         nxtrack1, xoff, ysidx, yeidx, ymidx, xsidx, xeidx, xmidx, iy, estat, nbits, ndim, ntmp
    REAL (KIND=r8), DIMENSION (1:nxtrack, 0:ntimes-1) :: sza, vza, saza, vaza, aza
    REAL (KIND=r4), DIMENSION (1:nxtrack)             :: tmp_lat, tmp_lon, tmp_sza, &
         tmp_vza, tmp_saza, tmp_vaza
    INTEGER (KIND=i2), DIMENSION (1:nxtrack)          :: tmp_geoflg, tmp_height
    INTEGER (KIND=i1), DIMENSION (1:nxtrack)          :: tmp_xtrackqflg
    INTEGER (KIND=i1), DIMENSION (1:nxtrack_max)      :: tmp_xtrackqflg1
    INTEGER (KIND=i2)                                 :: tmp_mflg
    INTEGER (KIND=i1)                                 :: tmp_NSPC
    REAL (KIND=r4)                                    :: tmp_auralon, tmp_auralat, &
         tmp_auraalt, tmp_SecondsInDay  
    REAL (KIND=r8)                                    :: tmp_time                 
    TYPE (L1B_block_type)                             :: omi_data_block

    ! Exteranl functions
    INTEGER                       :: OMI_SMF_setmsg
    INTEGER (KIND=i4), EXTERNAL   :: PGS_TD_TAItoUTC, PGS_TD_UTCtoTAI, day_of_year

    CHARACTER( LEN = 28 ) :: GranuleDAY0Z
    CHARACTER (LEN=21)    :: modulename = 'compute_pixel_corners'

    pge_error_status = pge_errstat_ok
    errstat          = omi_s_success

    ! -----------------------------------------------------------
    ! Open data block (UV-1, if both are selected) 
    ! called 'omi_data_block' with nTimes lines
    ! -----------------------------------------------------------
    errstat = L1Br_open ( omi_data_block, l1b_rad_filename, omi_radiance_swathname(1))
    IF( errstat /= omi_s_success ) THEN
       estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, "L1Br_open failed.", modulename, 0 ); STOP 1
    END IF

    sline = offset_line;  eline  = offset_line + nl * nybin - 1
    
    IF (eline == sline) eline = sline + 1
    nline = eline - sline + 1

    IF (zoom_mode .AND. nswath == 1) THEN
       nxtrack1 = nxtrack / 2; xoff = (zoom_p1 - 1) / nxbin
    ELSE
       nxtrack1 = nxtrack; xoff = 0
    ENDIF

    errstat = L1Br_getGEOline ( omi_data_block, 0, Time_k = TAI93StartOfGranule)
    IF( errstat /= omi_s_success ) THEN
       estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
            "L1Br_getGEOline failed.", modulename, 0 ); STOP 1
    ENDIF
    estat = PGS_TD_TAItoUTC(TAI93StartOfGranule, the_utc)
    READ (the_utc, '(I4, 1x, I2, 1x, I2, 1x, I2, 1x, I2, 1x, F9.6)') GranuleYear, &
         GranuleMonth, GranuleDay, GranuleHour, GranuleMinute, GranuleSecond
    WRITE( GranuleDAY0Z, FMT = '(I4.4,A1,I2.2,A1,I2.2,A)' ) &
         GranuleYear, '-', GranuleMonth, '-', GranuleDay, 'T00:00:00.000Z'
    estat = PGS_TD_UTCtoTAI( GranuleDAY0Z, TAI93At0zOfGranule )
    GranuleJDay = day_of_year(GranuleYear, GranuleMonth, GranuleDay)

    ! ---------------------------------
    ! Read all Latitudes and Longitudes
    ! ---------------------------------
    DO iline = sline, eline

       errstat = L1Br_getGEOline ( omi_data_block, iline,           &
            Time_k                      = tmp_time,                 &
            SecondsInDay_k              = tmp_SecondsInDay,         &
            SpacecraftAltitude_k        = tmp_auraalt,              &
            SpacecraftLongitude_k       = tmp_auralon,              &
            SpacecraftLatitude_k        = tmp_auralat,              &
            TerrainHeight_k             = tmp_height(1:nxtrack),    &
            GroundPixelQualityFlags_k   = tmp_geoflg(1:nxtrack),    &
            XTrackQualityFlags_k        = tmp_xtrackqflg(1:nxtrack),&
            Latitude_k                  = tmp_lat (1:nxtrack),      & 
            Longitude_k                 = tmp_lon (1:nxtrack),      &             
            SolarZenithANgle_k          = tmp_sza (1:nxtrack),      & 
            SolarAzimuthAngle_k         = tmp_saza(1:nxtrack),      & 
            ViewingZenithAngle_k        = tmp_vza (1:nxtrack),      & 
            ViewingAzimuthAngle_k       = tmp_vaza(1:nxtrack))  

       omi_alltime(iline)               = tmp_time
       omi_allSecondsInDay(iline)       = tmp_SecondsInDay    
       omi_allSpcftAlt(iline)           = tmp_auraalt
       omi_allSpcftLon(iline)           = tmp_auralon
       omi_allSpcftLat(iline)           = tmp_auralat
       omi_allHeight(1:nxtrack1, iline) = tmp_height(1:nxtrack1)  
       omi_allGeoFlg(1:nxtrack1, iline) = tmp_geoflg(1:nxtrack1)  
       omi_allXTrackQFlg(1:nxtrack1, iline) = tmp_xtrackqflg(1:nxtrack1)  
       omi_alllat(1:nxtrack1, iline)    = tmp_lat(1:nxtrack1)   
       omi_alllon(1:nxtrack1, iline)    = tmp_lon(1:nxtrack1)   
       sza (1:nxtrack1, iline)          = tmp_sza(1:nxtrack1)   
       saza(1:nxtrack1, iline)          = tmp_saza(1:nxtrack1)  
       vza (1:nxtrack1, iline)          = tmp_vza(1:nxtrack1)   
       vaza(1:nxtrack1, iline)          = tmp_vaza(1:nxtrack1) 
      
       IF( errstat /= omi_s_success ) THEN
          estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
               "L1Br_getGEOline failed.", modulename, 0 ); STOP 1
       ENDIF
    ENDDO

       
    ! --------------------------
    ! Close data block structure
    ! --------------------------
    errstat = L1Br_close ( omi_data_block )
    IF( errstat /= omi_s_success ) THEN
       estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
            "L1Br_close failed.", modulename, 0 ); STOP 1
    ENDIF

    ! Read measurement quality flags and NumberSmallPixelColumns (only in UV2 swath)

    ! For XtrackQuality flags, always use that from UV2 because UV1 and UV2 are different
    ! and more pixels are filtered in UV2

    errstat = L1Br_open ( omi_data_block, l1b_rad_filename, 'Earth UV-2 Swath')
    IF( errstat /= omi_s_success ) THEN
       estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, "L1Br_open failed.", modulename, 0 ); STOP 1
    ENDIF

    nbits = 8; ndim = 1
    DO iline = sline, eline
       errstat = L1Br_getDATA ( omi_data_block, iline, MeasurementQualityFlags_k = tmp_mflg)
       IF( errstat /= omi_s_success ) THEN
          estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
               "L1Br_getData failed.", modulename, 0 ); STOP 1
       END IF
       omi_allMflg(iline) = tmp_mflg
       
       errstat = L1Br_getSIGline ( omi_data_block, iline, NumberSmallPixelColumns_k = tmp_NSPC)
       IF( errstat /= omi_s_success ) THEN
          estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
               "L1Br_getSIGline failed.", modulename, 0 ); STOP 1
       END IF
       omi_allNSPC(iline) = tmp_NSPC

       IF (nswath == 2) THEN
          errstat = L1Br_getGEOline ( omi_data_block, iline,           &
               XTrackQualityFlags_k = tmp_xtrackqflg1(1:nxtrack_max))

          ! Coadd the flags
          DO ix = 1, nxtrack1
             i = ix * 2 - 1; j = i + 1
             CALL coadd_byte_qflgs(nbits, ndim, tmp_xtrackqflg1(i), tmp_xtrackqflg1(j)) 
             omi_allXTrackQFlg(ix, iline) = tmp_xtrackqflg1(i)
          ENDDO
       ENDIF      
    ENDDO

    errstat = L1Br_close ( omi_data_block )
    IF( errstat /= omi_s_success ) THEN
       estat = OMI_SMF_setmsg ( omsao_e_open_l1b_file, &
            "L1Br_close failed.", modulename, 0 ); STOP 1
    ENDIF

    WHERE (omi_allXTrackQFlg(1:nxtrack1, sline:eline) == -127)
       omi_allXTrackQFlg(1:nxtrack1, sline:eline) = 0
    ENDWHERE
    
    ! -----------------------------------------------------
    ! Compute the corner coordinates/viewing geometry for
    !    the spatially coadded pixels
    ! -----------------------------------------------------


      CALL get_sphgeoview_corners (nxtrack1, nline,  &
         omi_alllon(1:nxtrack1, sline:eline), omi_alllat(1:nxtrack1, sline:eline),       &
         sza(1:nxtrack1, sline:eline), saza(1:nxtrack1, sline:eline),                    &
         vza(1:nxtrack1, sline:eline), vaza(1:nxtrack1, sline:eline), &

         omi_allclon(0:nxtrack1, sline:eline+1), omi_allclat(0:nxtrack1, sline:eline+1), &
         omi_allelon(0:nxtrack1, sline:eline),   omi_allelat(0:nxtrack1, sline:eline),   &
         omi_allsza(1:nxtrack1, sline:eline),    omi_allvza(1:nxtrack1, sline:eline),    &
         omi_allaza(1:nxtrack1, sline:eline),    omi_allsca(1:nxtrack1, sline:eline))

   !DO iline = sline, eline
   !     WRITE(*,'(1i5, 15f8.3)') iline,sza(1:3,iline),vza(1:3,iline),saza(1:3,iline), vaza(1:3, iline), omi_allaza (1:3, iline)
   !ENDDO

    
     
    ! Re-store the data
    nx = nxtrack1 / nxbin 
    i = sline; j = i + nl - 1
!   print * , omi_allsca(1:nx, I:J), omi_allelat(0:nx, I:J)

    omi_alllon (1+xoff:nx+xoff, 0:nl-1)  = omi_alllon (1:nx, i:j)
    omi_alllat (1+xoff:nx+xoff, 0:nl-1)  = omi_alllat (1:nx, i:j)
    !sza        (1+xoff:nx+xoff, 0:nl-1)  = sza        (1:nx, i:j)
    !saza       (1+xoff:nx+xoff, 0:nl-1)  = saza       (1:nx, i:j)
    !vza        (1+xoff:nx+xoff, 0:nl-1)  = vza        (1:nx, i:j)
    !vaza       (1+xoff:nx+xoff, 0:nl-1)  = vaza       (1:nx, i:j)
    omi_allclon(0+xoff:nx+xoff, 0:nl)    = omi_allclon(0:nx, i:j+1)
    omi_allclat(0+xoff:nx+xoff, 0:nl)    = omi_allclat(0:nx, i:j+1)
    omi_allelon(0+xoff:nx+xoff, 0:nl-1)  = omi_allelon(0:nx, i:j)
    omi_allelat(0+xoff:nx+xoff, 0:nl-1)  = omi_allelat(0:nx, i:j)
    omi_allsza (1+xoff:nx+xoff, 0:nl-1)  = omi_allsza (1:nx, i:j)
    omi_allvza (1+xoff:nx+xoff, 0:nl-1)  = omi_allvza (1:nx, i:j)
    omi_allaza (1+xoff:nx+xoff, 0:nl-1)  = omi_allaza (1:nx, i:j)
    omi_allsca (1+xoff:nx+xoff, 0:nl-1)  = omi_allsca (1:nx, i:j)


     !print * , omi_alllat(2, i:j) , i, j
 
    ! Need to get Time, SecondsInDay, Spacecfraft altitude/latitude/longitude, 
    ! GroundPixelQualityFlags, XTrackQualityFlags, and Terrain Height for the spatially coadded pixels
    DO iy = 0, nl - 1 
       ysidx = sline + iy * nybin 
       yeidx = ysidx + nybin - 1
       ymidx = ysidx + nybin / 2
       omi_alltime(iy)         = SUM(omi_alltime(ysidx:yeidx))      / nybin
       omi_allSecondsInDay(iy) = SUM(omi_allSecondsInDay(ysidx:yeidx)) / nybin
       omi_allSpcftAlt(iy)     = SUM(omi_allSpcftAlt(ysidx:yeidx))  / nybin

       ! Use those from the middle point (avoid dealing with polar, dateline regions)
       omi_allSpcftLat(iy)     = omi_allSpcftLat(ymidx)
       omi_allSpcftLon(iy)     = omi_allSpcftLon(ymidx)
       omi_allMflg(iy)         = omi_allMflg(ymidx)
       omi_allNSPC(iy)         = omi_allNSPC(ymidx)       
       
       DO ix = 1, nx 
          xsidx = (ix - 1) * nxbin + 1
          xeidx = xsidx + nxbin - 1
          xmidx = xsidx + nxbin / 2
          omi_allHeight(ix, iy)  = INT(SUM(1.0 * omi_allHeight(xsidx:xeidx, ysidx:yeidx)) &
               / (1.0 * nxbin * nybin))
          omi_allGeoFlg(ix, iy) = omi_allGeoFlg(xmidx, ymidx)
          omi_allXTrackQFlg(ix, iy) = omi_allXTrackQFlg(xmidx, ymidx)
       ENDDO
    ENDDO
    IF (xoff > 0) THEN
       omi_allGeoFlg(1+xoff:nx+xoff, 0:nl-1) = omi_allGeoFlg(1:nx, 0:nl-1)
       omi_allXTrackQFlg(1+xoff:nx+xoff, 0:nl-1) = omi_allXTrackQFlg(1:nx, 0:nl-1)
       omi_allHeight(1+xoff:nx+xoff, 0:nl-1) = omi_allHeight(1:nx, 0:nl-1)
    ENDIF
        
    !DO ix = 1, nxtrack
    !   WRITE(90, '(I5, 2F10.4)'), ix, omi_alllat(ix, 1), omi_alllon(ix, 1)       
    !ENDDO
    !DO ix = 0, nx
    !   WRITe(90, '(I5, 4F10.4)') ix, omi_allclat(ix, 1:2), omi_allclon(ix, 1:2)
    !ENDDO
    !stop
    
    RETURN
  END SUBROUTINE compute_pixel_corners


  SUBROUTINE sphergeom_baseline_comp ( a0, b0, gam0, c0 )
    ! -------------------------------------------------------
    ! Finds the lengh of the baseline of a spherical triangle
    ! -------------------------------------------------------
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8), INTENT (IN) :: a0, b0, gam0

    ! ---------------
    ! Output variable
    ! ---------------
    REAL (KIND=r8), INTENT (OUT) :: c0

    ! --------------
    ! Local variable
    ! --------------
    REAL (KIND=r8) :: tmp

    ! --------------------------------------------------
    ! Initialize output variables to keep compiler happy
    ! --------------------------------------------------
    c0 = 0.0_r8

    ! -------------------------------------------------
    ! Compute length of baseline between the two points
    ! -------------------------------------------------
    tmp = COS(a0) * COS(b0) + SIN(a0) * SIN(b0) * COS(gam0)
    IF ( ABS(tmp) < eps ) THEN
       c0 = pihalf
    ELSE
       c0 = ACOS(tmp)
    END IF

    RETURN
  END SUBROUTINE sphergeom_baseline_comp

  SUBROUTINE lonlat_to_pi ( lon, lat )

    IMPLICIT NONE

    ! ------------------
    ! Modified variables
    ! ------------------
    REAL (KIND=r8), INTENT (INOUT) :: lon, lat

    ! ------------------------------------
    ! Adjust longitude values to [-pi,+pi]
    ! ------------------------------------
    IF ( ABS(lon) > twopi ) lon = MOD(lon, twopi)
    IF ( lon >  pi ) lon = lon - twopi
    IF ( lon < -pi ) lon = lon + twopi
    ! ---------------------------------------
    ! Adjust latitude values to [-pi/2,+pi/2]
    ! ---------------------------------------
    IF ( ABS(lat) > pihalf ) lat = MOD(lat, pihalf)
    IF ( lat >  pihalf ) lat =   pi - lat
    IF ( lat < -pihalf ) lat = -(pi + lat)
    
    RETURN
  END SUBROUTINE lonlat_to_pi

  REAL (KIND=r8) FUNCTION angle_minus_twopi ( gamma0, pival ) RESULT ( gamma )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8), INTENT (IN) :: gamma0, pival

    IF ( gamma0 > pival ) THEN
       gamma = gamma0 - 2.0_r8 * pival !SIGN(2.0_r8*pival - gamma0, gamma0)
    ELSE IF ( gamma0 < -pival ) THEN
       gamma = gamma0 + 2.0_r8 * pival 
    ELSE
       gamma = gamma0
    END IF

    RETURN
  END FUNCTION angle_minus_twopi

  SUBROUTINE get_sphgeoview_corners (nxtrack, ntimes, lon, lat, sza, saza, vza, vaza, &
         clon, clat, elon, elat, esza, evza, eaza, esca)

    IMPLICIT NONE
    
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                                    INTENT(IN)    :: nxtrack, ntimes
    REAL    (KIND=r8), DIMENSION (1:nxtrack, 0:ntimes-1), INTENT(INOUT) :: lon, lat, sza, saza, vza, vaza
    REAL    (KIND=r8), DIMENSION (0:nxtrack, 0:ntimes),   INTENT(OUT)   :: clon, clat
    REAL    (KIND=r8), DIMENSION (0:nxtrack, 0:ntimes-1), INTENT(OUT)   :: elon, elat
    REAL    (KIND=r8), DIMENSION (1:nxtrack, 0:ntimes-1), INTENT(OUT)   :: esza, evza, eaza, esca     

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                                   :: i, j, jj, ix, mpix, nx, ny
    INTEGER                                             :: errstat
    !REAL    (KIND=r8)                                  :: lat1, lat2, phi1, phi2, dis
    !REAL    (KIND=r8)                                  :: a0, b0, c0, gam0, alp0, a, gam
    REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1) :: omixsize
    REAL    (KIND=r8), DIMENSION (0:nxtrack, 0:ntimes-1):: edsza, edsazm, edvza, edvazm
    REAL    (KIND=r8), DIMENSION (1:nxtrack)            :: tmpdisx, xsize, tmpxmid
    REAL    (KIND=r8), DIMENSION (0:nxtrack)            :: tmpx, tmpsza, tmpvza, tmpsaza, tmpvaza
    REAL    (KIND=r8), DIMENSION (0:ntimes-1)           :: tmpdisy
    REAL    (KIND=dp), DIMENSION (3)                    :: zen0, zen, sazm, vazm, relaza

    ! ------------------------------------------------------------------
    ! Convert geolocation to radians; do everything in R8 rather than R4
    ! ------------------------------------------------------------------
    lon = lon * deg2rad; lat = lat * deg2rad
    
    ! -------------------------
    ! Initialize some variables
    ! -------------------------
    clon   = -999.9 ; clat = -999.9
    elon   = -999.9 ; elat = -999.9
    esza   = -999.9 ; evza = -999.9; eaza = -999.9; esca = -999.9

    ! Perform interpolation across the track
    DO i = 0, ntimes - 1
       ! Compute the distances between two pixels: (x1 + x2) / 2.
       DO ix = 1, nxtrack - 1 
          tmpdisx(ix) = circle_rdis(lat(ix, i), lon(ix, i), lat(ix+1, i), lon(ix+1, i))
       ENDDO
       ! Compute the pixel size across the track
       ! Assume the center two pixels have equal pixel size (which causes about < 0.1 km error for UV-2)
       mpix = nxtrack / 2
       xsize(mpix) = tmpdisx(mpix) / 2.0;  xsize(mpix + 1) = xsize(mpix)
       
       DO ix = mpix -1, 1, -1
          xsize(ix) = tmpdisx(ix) - xsize(ix + 1)
       ENDDO
       DO ix = mpix + 2, nxtrack
          xsize(ix) = tmpdisx(ix - 1) - xsize(ix - 1)
       ENDDO
       omixsize(:, i) = xsize

       !!  This is to test SUBROUTINE sphergeom_intermediate
       !!  Works for both interpolation and extrapolation (with certain limitation) 
       !lat1 = -88.0 * deg2rad; lat2 = -88.0 * deg2rad
       !phi1 = 0.0  * deg2rad; phi2 = 180.0 * deg2rad 
       !dis  = circle_rdis(lat1, phi1, lat2, phi2)     
       !print *, dis
       !CALL sphergeom_intermediate ( lat1, phi1, lat2, phi2, dis, dis*0.8, a, gam )
       !WRITE(*, '(6D14.6)')  lat1, phi1, a, gam, lat2, phi2
       !print *, circle_rdis(lat1, phi1, a, gam)/dis
       !print *, circle_rdis(a, gam, lat2, phi2)/dis
  
       ! Perform interpolation 
       DO ix = 1, nxtrack - 1       
            here_stop = 0
          if (ix == 2) here_stop = 1
          CALL sphergeom_intermediate(lat(ix, i), lon(ix, i), lat(ix+1, i), lon(ix+1, i), &
               tmpdisx(ix), xsize(ix), elat(ix, i), elon(ix, i))
!          if ( ix == 2) write(*,'(10f8.4)') lat(ix, i), lon(ix, i),lat(ix+1,i), lon(ix+1, i), tmpdisx(ix), xsize(ix),elat(ix, i), elon(ix, i)
       ENDDO
        
       ix = 1
       CALL sphergeom_intermediate(lat(ix, i), lon(ix, i), lat(ix+1, i), lon(ix+1, i),   &
            tmpdisx(ix), -xsize(ix), elat(ix-1, i), elon(ix-1, i))
       ix = nxtrack
       CALL sphergeom_intermediate(lat(ix, i), lon(ix, i), lat(ix-1, i), lon(ix-1, i),   &
            tmpdisx(ix-1), -xsize(ix), elat(ix, i), elon(ix, i))   

       ! Compute viewing geometry for west and east edge
       ! Performal interpolation/extrapolation (2 points) along the spherical lines (good enough)
       tmpx = 0.0
       DO ix = 1, nxtrack
          tmpx(ix) = tmpx(ix-1) + xsize(ix)
       ENDDO
       tmpxmid(1:nxtrack) = (tmpx(0:nxtrack-1) + tmpx(1:nxtrack)) / 2.0
       CALL interpol(tmpxmid(1:nxtrack), sza(1:nxtrack, i),  nxtrack, tmpx(0:nxtrack), &
            tmpsza(0:nxtrack),  nxtrack+1, errstat)      
       CALL interpol(tmpxmid(1:nxtrack), saza(1:nxtrack, i), nxtrack, tmpx(0:nxtrack), &
            tmpsaza(0:nxtrack), nxtrack+1, errstat)
       CALL interpol(tmpxmid(1:nxtrack), vza(1:nxtrack, i),  nxtrack, tmpx(0:nxtrack), &
            tmpvza(0:nxtrack),  nxtrack+1, errstat)
       CALL interpol(tmpxmid(1:nxtrack), vaza(1:nxtrack, i), nxtrack, tmpx(0:nxtrack), &
            tmpvaza(0:nxtrack), nxtrack+1, errstat)

       ! Check center pixel
       DO ix = mpix-1, mpix + 1
          IF (tmpvaza(ix) < 0) THEN
             tmpvaza(ix) = -tmpvaza(ix-1); EXIT
          ENDIF
       ENDDO

       edsza(0:nxtrack, i) = tmpsza(0:nxtrack); edsazm(0:nxtrack, i) = tmpsaza(0:nxtrack)
       edvza(0:nxtrack, i) = tmpvza(0:nxtrack); edvazm(0:nxtrack, i) = tmpvaza(0:nxtrack)
    ENDDO

    !PRINT *
    !WRITE(*, *) 'Solar Angle'
    !WRITE(*, '(10F9.2)') esza(1:nxtrack, 0)
    !WRITE(*, *) 'View Angle'
    !WRITE(*, '(10F9.2)') evza(1:nxtrack, 0)
    !WRITE(*, *) 'Relative Azmimuthal Angle'
    !WRITE(*, '(10F9.2)') eaza(1:nxtrack, 0)
    !WRITE(*, *) 'Scattering Angle'
    !WRITE(*, '(10F9.2)') esca(1:nxtrack, 0)
 
    ! Perform interpolation along the track with a simipler but simpler (center) approach
    ! since the pixel size along the track does not vary much
    IF (ntimes == 1) tmpdisy(0) = 0.00212031  ! ~ 13.5 km
    DO ix = 0, nxtrack
       DO i = 0, ntimes - 2 
          tmpdisy(i) = circle_rdis(elat(ix, i), elon(ix, i), elat(ix, i+1), elon(ix, i+1))
          CALL sphergeom_intermediate(elat(ix, i), elon(ix, i), elat(ix, i+1), &
               elon(ix, i+1), tmpdisy(i), tmpdisy(i)*0.5, clat(ix, i+1), clon(ix, i+1))  
          !print *, i+1, clat(ix, i+1), clon(ix, i+1)
       ENDDO
 
       i = 0
       CALL sphergeom_intermediate(elat(ix, i), elon(ix, i), elat(ix, i+1),    &
            elon(ix, i+1), tmpdisy(i), -tmpdisy(i)*0.5, clat(ix, i), clon(ix, i))
       !print *, i, clat(ix, i), clon(ix, i)

       i = ntimes - 1       
       CALL sphergeom_intermediate(elat(ix, i), elon(ix, i), elat(ix, i-1),    &
            elon(ix, i-1), tmpdisy(i-1), -tmpdisy(i-1)*0.5, clat(ix, i+1), clon(ix, i+1))
       !print *, i+1, clat(ix, i+1), clon(ix, i+1)
    ENDDO

    ! Perform coadding
    IF (nxbin > 1 .OR. nybin > 1) THEN
       nx = nxtrack / nxbin    
       ny = ntimes  / nybin

       ! cornor coordinates (only need sampling)
       j = 0
       DO ix = 0, nxtrack, nxbin
          clon(j, :) = clon(ix, :)
          clat(j, :) = clat(ix, :)
          j = j + 1
       ENDDO

       j = 0
       DO i = 0, ntimes, nybin
          clon(0:nx, j) = clon(0:nx, i)
          clat(0:nx, j) = clat(0:nx, i)
          j = j + 1
       ENDDO

       ! edge coordinates (easy to be re-computed from cornor coordinates)
       DO ix = 0, nx
          DO i = 0, ny - 1 
             tmpdisy(i) = circle_rdis(clat(ix, i), clon(ix, i), clat(ix, i+1), clon(ix, i+1))
             CALL sphergeom_intermediate(clat(ix, i), clon(ix, i), clat(ix, i+1), clon(ix, i+1), &
                  tmpdisy(i), tmpdisy(i)*0.5, elat(ix, i), elon(ix, i))  
             !print *, tmpdisy(i), elat(ix, i)*rad2deg, elon(ix, i) * rad2deg
             !print *, clat(ix, i)*rad2deg, clon(ix, i)*rad2deg, clat(ix, i+1)*rad2deg, clon(ix, i+1)*rad2deg
          ENDDO
       ENDDO

       ! Center coordinates (computed from edge coordinates)
       DO ix = 1, nx
          DO i = 0, ny - 1 
            
             tmpdisx(ix) = circle_rdis(elat(ix-1, i), elon(ix-1, i), elat(ix, i), elon(ix, i))
             CALL sphergeom_intermediate(elat(ix-1, i), elon(ix-1, i), elat(ix, i), elon(ix, i), &
                  tmpdisx(ix), tmpdisx(ix)*0.5, lat(ix, i), lon(ix, i))
             !print *, elat(ix-1, i)*rad2deg, elon(ix-1, i)*rad2deg, elat(ix, i)*rad2deg, elon(ix, i)*rad2deg
          ENDDO
       ENDDO 

       ! Average edge viewing geometries along the track, sample along the track
       i = 0
       DO ix = 0, nxtrack, nxbin
          jj = 0
          DO j = 0, ntimes-1, nybin
             edsza (i, jj) = SUM(edsza (ix, j:j+nybin-1)) / nybin
             edsazm(i, jj) = SUM(edsazm(ix, j:j+nybin-1)) / nybin
             edvza (i, jj) = SUM(edvza (ix, j:j+nybin-1)) / nybin
             edvazm(i, jj) = SUM(edvazm(ix, j:j+nybin-1)) / nybin
             jj = jj + 1
          ENDDO
          i = i + 1
       ENDDO   

       ! Compute center viewing geometries (interpolate across the track)
       DO i = 0, ny-1
          tmpx = 0.0
          DO ix = 1, nx
             tmpx(ix) = tmpx(ix-1) + circle_rdis(elat(ix-1, i), elon(ix-1, i), elat(ix, i), elon(ix, i))
          ENDDO
          tmpxmid(1:nx) = (tmpx(0:nx-1) + tmpx(1:nx)) / 2.0

          CALL interpol(tmpx(0:nx), edsza (0:nx, i),  nx+1, tmpxmid(1:nx), sza (1:nx, i),  nx, errstat) 
          CALL interpol(tmpx(0:nx), edsazm(0:nx, i),  nx+1, tmpxmid(1:nx), saza(1:nx, i),  nx, errstat)  
     
          CALL interpol(tmpx(0:nx), edvza (0:nx, i),  nx+1, tmpxmid(1:nx), vza (1:nx, i),  nx, errstat) 
          CALL interpol(tmpx(0:nx), edvazm(0:nx, i),  nx+1, tmpxmid(1:nx), vaza(1:nx, i),  nx, errstat)      

          ! Check center pixel
          mpix = nx / 2
          DO ix = mpix - 1, mpix + 1
             IF (vaza(ix, i) < 0) THEN
                vaza(ix, i) = -vaza(ix-1, i); EXIT
             ENDIF
          ENDDO  
       ENDDO
    ELSE
       nx = nxtrack;       ny = ntimes
    ENDIF
    
    ! Now compute effective viewing geometry
    ! Compute effective viewing geometry for each pixel (at a certain atmosphere) 
    DO i = 0, ny-1
       DO ix = 1, nx
          IF ( sza(ix, i)  >= minza  .AND. sza(ix, i)  < maxza  .AND. &
               vza(ix, i)  >= minza  .AND. vza(ix, i)  < maxza  .AND. &
               saza(ix, i) >= minaza .AND. saza(ix, i) < maxaza .AND. &
               vaza(ix, i) >= minaza .AND. vaza(ix, i) < maxaza) THEN
             zen0(1) = edsza(ix-1, i) ; zen0(2) = sza(ix, i) ; zen0(3) = edsza(ix, i)
             zen(1)  = edvza(ix-1, i) ; zen(2)  = vza(ix, i) ; zen(3)  = edvza(ix, i)
             sazm(1) = edsazm(ix-1, i); sazm(2) = saza(ix, i); sazm(3) = edsazm(ix, i)
             vazm(1) = edvazm(ix-1, i); vazm(2) = vaza(ix, i); vazm(3) = edvazm(ix, i)
             relaza  = ABS(vazm - sazm)
             
             WHERE (vazm > 0) 
                zen = - zen
             ENDWHERE
             
             ! -180 < relaza < 180          
             WHERE (relaza > 180.0)
                relaza = 360.0 - relaza
             ENDWHERE
             
             CALL omi_angle_sat2toa (3, zen0, zen, relaza, esza(ix, i), evza(ix, i), eaza(ix, i), esca(ix, i) )
             !print *, zen0 * rad2deg
             !print *, zen  * rad2deg
             !print *, 180.0 - relaza * rad2deg
             !print *, esza(ix, i), evza(ix, i), eaza(ix, i), esca(ix, i)
             !STOP             
          ENDIF
       ENDDO
    ENDDO
    
    clon(0:nx, 0:ny)   = clon(0:nx, 0:ny)   * rad2deg    
    clat(0:nx, 0:ny)   = clat(0:nx, 0:ny)   * rad2deg
    lon(1:nx, 0:ny-1)  = lon(1:nx, 0:ny-1)  * rad2deg    
    lat(1:nx, 0:ny-1)  = lat(1:nx, 0:ny-1)  * rad2deg
    elon(0:nx, 0:ny-1) = elon(0:nx, 0:ny-1) * rad2deg
    elat(0:nx, 0:ny-1) = elat(0:nx, 0:ny-1) * rad2deg
    WHERE (clon(0:nx, 0:ny) > 180.0)
       clon(0:nx, 0:ny) = clon(0:nx, 0:ny) - 360.0
    ENDWHERE
    WHERE (clon(0:nx, 0:ny) < -180.0)
       clon(0:nx, 0:ny) = clon(0:nx, 0:ny) + 360.0
    ENDWHERE

    WHERE (lon(1:nx, 0:ny-1) > 180.0)
       lon(1:nx, 0:ny-1) = lon(1:nx, 0:ny-1) - 360.0
    ENDWHERE
    WHERE (lon(1:nx, 0:ny-1) < -180.0)
       lon(1:nx, 0:ny-1) = lon(1:nx, 0:ny-1) + 360.0
    ENDWHERE

    WHERE (elon(0:nx, 0:ny-1) > 180.0)
       elon(0:nx, 0:ny-1) = elon(0:nx, 0:ny-1) - 360.0
    ENDWHERE
    WHERE (elon(0:nx, 0:ny-1) < -180.0)
       elon(0:nx, 0:ny-1) = elon(0:nx, 0:ny-1) + 360.0
    ENDWHERE
   

    RETURN

  END SUBROUTINE get_sphgeoview_corners

  ! Compute spherical distance between two points
  ! Haversine Formula (from R.W. Sinnott, "virtue of the Haversine", 
  ! Sky and Telescope V68 (2), 1984, p159  
  ! lat1, lon1, lat2, lon2 in radians
  FUNCTION circle_dis(lat1, lon1, lat2, lon2) RESULT(dis)

    IMPLICIT NONE

    ! ----------------------
    ! Input/output variables
    ! -----------------------
    REAL (KIND=r8), INTENT (IN) :: lat1, lon1, lat2, lon2
    REAL (KIND=r8)              :: dis

    ! Local variable
    REAL (KIND=r8)              :: dlon, dlat, a, rdis, mlat
     
    dlat = lat2 - lat1;    dlon = lon2 - lon1; mlat = (lat1 + lat2) / 2.0

    a = MIN(1.0, SQRT( SIN(dlat/2.0)**2.0 + COS(lat1) * COS(lat2) * SIN(dlon/2.0)**2.0 )  )
    rdis = 2.0 * ASIN(a)                                    ! relative distaince in radiances
    dis = rdis * (rearth0 - 21.0 * SIN(mlat))               ! in the same unit of rearth0

    RETURN

  END FUNCTION circle_dis


  FUNCTION circle_rdis(lat1, lon1, lat2, lon2) RESULT(rdis)

    IMPLICIT NONE

    ! ----------------------
    ! Input/output variables
    ! -----------------------
    REAL (KIND=r8), INTENT (IN) :: lat1, lon1, lat2, lon2
    REAL (KIND=r8)              :: rdis

    ! Local variable
    REAL (KIND=r8)              :: dlon, dlat, a

    dlat = lat2 - lat1;    dlon = lon2 - lon1
    a = MIN(1.0, SQRT( SIN(dlat/2.0)**2.0 + COS(lat1) * COS(lat2) * SIN(dlon/2.0)**2.0 )  )
    rdis = 2.0 * ASIN(a)                                    ! relative distaince in radiances
 
    RETURN

  END FUNCTION circle_rdis

  ! This one works, gam is the longitude difference (-pi < gam0 < pi) between 2 and 1 (phi2 - phi1)
  SUBROUTINE sphergeom_intermediate ( lat1, lon1, lat2, lon2, c0, c, lat, lon )

    ! -----------------------------------------------------------------
    ! Finds the co-ordinates of C the baseline extended from two
    ! lon/lat points (A, B) on a sphere given the hypotenuse C_IN.
    ! ----------------------------------------------------------------
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8),    INTENT (IN) :: lat1, lat2, lon1, lon2, c0, c

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8),  INTENT (OUT)  :: lat, lon

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8)  :: x, y, z, tmp1, tmp2, frc, gamsign, theta, gam0, gam

    lat = 0.0_r8 ; lon = 0.0_r8    
    gam0 = angle_minus_twopi ( lon2 - lon1, pi )

    gamsign = ABS(gam0) / gam0;  gam0 = ABS(gam0)
    if ( gam0 == 0 ) gamsign = 1 ! jbak correction
    ! Get straight line (AB) segment fraction frc intercepted by the line from center to C
    ! If frc < 0, extrapolation, but it is limited to |c| < (180-c0)/2.0
    tmp1 = SIN(c)
    frc = tmp1 / (SIN(c0 - c) + tmp1)
    
    ! Work in Cartesian Coordinate 
    tmp1 = frc * COS(lat2); tmp2 = 1.0 - frc
    x = tmp2 *   COS(lat1) + tmp1 * COS(gam0)
    y = tmp1 *   SIN(gam0)
    z = tmp2 *   SIN(lat1) + frc * SIN(lat2)

    gam = ATAN(y/x)                          ! -90 < gam < 90
    IF (frc >= 0) THEN
       IF (gam < 0) gam = gam + pi           ! 0 <= gam <= 180
    ELSE
       IF (gam > 0) gam = gam - pi           ! -180 <= gam <= 0
    ENDIF
    gam = gamsign * gam                      ! Get correct sign
    lon = gam + lon1
   
 
    theta = ATAN (SQRT(x**2 + y**2) / z)     ! -90 < theta < 90
    IF (theta < 0) theta = theta + pi        ! 0 <= theta <= 180
    lat = pihalf - theta                     ! Convert to latitud
    RETURN
  END SUBROUTINE sphergeom_intermediate


END MODULE OMSAO_pixelcorner_module
