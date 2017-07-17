  ! *********************** Modification History *******************
  ! xiong liu, july 2003 (!xliu)
  ! 1. read solar/earthshine spectrum around 370 nm for initializing
  !    effective albedo into LIDORT
  ! 2. Read longitude and latitude at pixel center for initializing
  !    atmospheric profiles
  ! 4. Don't read data where radiance/irradiance is 0
  ! ****************************************************************

SUBROUTINE gome_read_el1data_sol (funit, file_read_stat, error)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,ONLY : maxchlen
  USE OMSAO_gome_data_module, ONLY : lm_gome_solspec, lm_gome_eshine,&
       n_gome_max_pts, n_gome_data_dim, gome_spec_missing, &
       n_gome_solpts, gome_solspec
  USE OMSAO_errstat_module,   ONLY : file_read_ok, file_read_failed, &
       file_read_missing, file_read_eof
  USE OMSAO_variables_module, ONLY : numwin, winpix, winlim, &
       band_selectors, nsolpix, sol_spec_ring, nsol_ring
  USE ozprof_data_module,     ONLY : div_sun

  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER,  INTENT (IN)       :: funit
  LOGICAL, INTENT (OUT)       :: error	
  INTEGER, INTENT (OUT)       :: file_read_stat

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                     :: i, j
  CHARACTER (LEN=maxchlen)    :: lastline
  REAL (KIND=dp)              :: solar_norm
  
  file_read_stat = file_read_ok ; error = .FALSE.

  ! ------------------------------------------
  ! Position cursor at start of solar spectrum
  ! ------------------------------------------
  CALL skip_to_filemark ( funit, lm_gome_solspec, lastline, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) RETURN

  ! Read solar spectrum
  CALL gome_read_el1solspec (funit, numwin, winpix(1:numwin, :), &
       winlim(1:numwin, :), band_selectors(1:numwin), n_gome_solpts, &
       gome_solspec, nsolpix(1:numwin), file_read_stat)  
  IF ( file_read_stat /= file_read_ok ) RETURN

  ! -----------------------------------
  ! Normalize solar irradiance spectrum
  ! -----------------------------------
  solar_norm = SUM (gome_solspec(2,1:n_gome_solpts)) / n_gome_solpts
  IF ( solar_norm <= 0.0 ) THEN 
     error = .TRUE.; RETURN
  ENDIF

  gome_solspec(2,1:n_gome_solpts) = gome_solspec(2,1:n_gome_solpts) / solar_norm
  gome_solspec(3,1:n_gome_solpts) = gome_solspec(3,1:n_gome_solpts) / solar_norm
  sol_spec_ring(2, 1:nsol_ring) = sol_spec_ring(2, 1:nsol_ring)     / solar_norm 
  div_sun = solar_norm

  ! ---------------------------------------------------------------------
  ! Position cursor after solar spectrum, before first ground pixel entry
  ! ---------------------------------------------------------------------
  CALL skip_to_filemark ( funit, lm_gome_eshine, lastline, file_read_stat )
  
  RETURN

END SUBROUTINE gome_read_el1data_sol


SUBROUTINE gome_read_el1data_rad ( funit, file_read_stat, error )

  USE OMSAO_precision_module
  USE OMSAO_errstat_module,    ONLY : &
       file_read_ok, file_read_failed, file_read_missing, file_read_eof
  USE OMSAO_parameters_module, ONLY : maxchlen, deg2rad, rad2deg
  USE OMSAO_gome_data_module,  ONLY : &
       azm_idx, zen_idx, zen0_idx, lat_idx, lon_idx, n_gome_ang, n_gome_geo,   &
       lm_gome_groundpix, gome_spec_missing, n_gome_max_pts, n_gome_data_dim,  &
       n_gome_radpts, gome_radspec, gome_curpix, gome_curqual, gome_curscan,   &
       gome_pixdate, gome_angles_wrtn, gome_angles_wrts, gome_geoloc, ers2_alt,&
       earth_curv, orbnum, gome_stpix, gome_endpix, gome_npix
  USE ozprof_data_module,      ONLY : div_rad, coadd_after_b1ab, b1ab_change
  USE OMSAO_variables_module,  ONLY : the_month, the_year, the_day, numwin,    &
       winlim, winpix, band_selectors, nradpix, npix_fitting, sza_atm, vza_atm,&
       aza_atm, the_sza_atm, the_aza_atm, the_vza_atm, the_sca_atm, the_lat,   &
       the_lon, the_lons, the_lats, nview, nloc, edgelons, edgelats, szamax

  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER,  INTENT (IN)    :: funit
  LOGICAL, INTENT (OUT)    :: error	
  INTEGER, INTENT (OUT)    :: file_read_stat

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                  :: i, j, n_coadd
  LOGICAL                  :: view_coadd
  CHARACTER (LEN=maxchlen) :: lastline
  CHARACTER (LEN=3)        :: monthc
  CHARACTER (LEN=2)        :: day
  CHARACTER (LEN=4)        :: year
  REAL (KIND=dp)           :: spec_norm, lon1, lon2, templon, templat, temp
  REAL (KIND=dp), DIMENSION(n_gome_ang):: sumsza, sumaza, sumvza
  REAL (KIND=dp)  :: sumtsza, sumtaza, sumtvza, sumalt, sumcurv
  REAL (KIND=dp), DIMENSION(lon_idx, n_gome_geo) :: sumloc
  CHARACTER (LEN=3), DIMENSION(12), PARAMETER :: months = (/'JAN', 'FEB', 'MAR', &
       'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/)

  file_read_stat = file_read_ok ; error = .FALSE.
  
  IF ((band_selectors(1) > 1) .OR. (.NOT. coadd_after_b1ab .AND. b1ab_change)) THEN
     view_coadd = .FALSE.
  ELSE
     view_coadd = .TRUE.
  ENDIF

  IF (npix_fitting == 0 .AND. orbnum == 0) THEN
     CALL skip_to_filemark ( funit, 'ERS Information', lastline, file_read_stat )
     IF ( file_read_stat /= file_read_ok ) RETURN
     READ (funit, *), i, i, orbnum
  ENDIF

  IF (.NOT. view_coadd) THEN
     ! ----------------------------------------------------------------
     ! Position cursor at start of GOME ground pixel entry
     ! And Extract GOME ground pixel number and scan type from LASTLINE
     ! -----------------------------------------------------------------
     CALL skip_to_filemark ( funit, lm_gome_groundpix, lastline, file_read_stat)
     IF ( file_read_stat /= file_read_ok ) RETURN
     READ (lastline(14:17), '(I4)') gome_curpix
     READ (lastline(19:19), '(I1)') gome_curqual
     READ (lastline(21:21), '(I1)') gome_curscan
     IF ( gome_curqual == gome_spec_missing) THEN
        file_read_stat = file_read_missing; RETURN
     END IF

     CALL gome_read_el1geo (funit, gome_pixdate, gome_angles_wrtn, &
          gome_angles_wrts, gome_geoloc, ers2_alt, earth_curv, file_read_stat)  
     IF ( ANY(gome_angles_wrtn(zen0_idx, 1:n_gome_ang) >= szamax) ) THEN
        file_read_stat = file_read_missing; RETURN
     END IF   

     CALL calc_view_geometry
     IF ( the_sza_atm < 0.0 .OR. the_sza_atm >= szamax) THEN
        file_read_stat = file_read_missing; RETURN
     END IF

     IF (gome_curscan == 3) THEN
        the_vza_atm = the_vza_atm + 1.955
     ELSE IF (gome_curscan == 0 .OR. gome_curscan == 2) THEN
        the_vza_atm = the_vza_atm - 1.23
     ENDIF
     
     gome_stpix = gome_curpix; gome_endpix = gome_curpix; gome_npix = 1
  ELSE
     n_coadd =   0;    sumsza=   0.0;   sumvza  = 0.0;   sumaza = 0.0
     sumtaza = 0.0;    sumtvza = 0.0;   sumtsza = 0.0
     sumalt =  0.0;    sumcurv = 0.0
     
     DO
        CALL skip_to_filemark (funit, lm_gome_groundpix, lastline, file_read_stat)
        IF ( file_read_stat /= file_read_ok ) RETURN
        READ (lastline(14:17), '(I4)') gome_curpix
        IF (n_coadd == 0) gome_stpix = gome_curpix

        READ (lastline(19:19), '(I1)') gome_curqual
        READ (lastline(21:21), '(I1)') gome_curscan
        CALL gome_read_el1geo (funit, gome_pixdate, gome_angles_wrtn, &
             gome_angles_wrts, gome_geoloc, ers2_alt, earth_curv, file_read_stat)
        IF ( ANY(gome_angles_wrtn(zen0_idx, 1:n_gome_ang) >= szamax) ) THEN
           file_read_stat = file_read_missing; RETURN
        END IF
        CALL calc_view_geometry
        IF ( n_coadd == 3 .AND. gome_curscan == 3) THEN 
           sumloc = gome_geoloc
        ENDIF
        IF ( the_sza_atm < 0.0 .OR. the_sza_atm >= szamax) THEN
           file_read_stat = file_read_missing; RETURN
        END IF

        !WRITE(*, '(4x, A, 4f8.2)') 'Effective SZA, VZA, AZA = ', &
        !     the_sza_atm, the_vza_atm, the_aza_atm, the_sca_atm  
       
        IF (gome_curscan == 3) THEN
           the_vza_atm = the_vza_atm + 1.955
        ELSE IF (gome_curscan == 0 .OR. gome_curscan == 2) THEN
           the_vza_atm = the_vza_atm - 1.23
        ENDIF
             
        n_coadd = n_coadd + 1
        
        sumtsza = sumtsza + 1.0 / COS(the_sza_atm * deg2rad)
        sumsza = sumsza   + 1.0 / COS(sza_atm * deg2rad)
        sumtvza = sumtvza + 1.0 / COS(the_vza_atm * deg2rad)
        sumvza = sumvza   + 1.0 / COS(vza_atm * deg2rad)     
        sumtaza = sumtaza + the_aza_atm
        sumaza = sumaza   + aza_atm  
        sumcurv = sumcurv +  earth_curv
        sumalt = sumalt   +  ers2_alt
        IF (gome_curqual /= gome_spec_missing) EXIT
     ENDDO

     gome_endpix = gome_curpix
     gome_npix = n_coadd

     IF (n_coadd > 1 .AND. gome_curscan == 3) THEN
        sumloc(:, 1) = gome_geoloc(:, 1)
        sumloc(:, 3) = gome_geoloc(:, 3) 
        sumloc(:, 5) = (gome_geoloc(:, 5) + sumloc(:, 5)) / 2.0D0
     ENDIF
     
     IF (n_coadd > 1)  THEN
        the_sza_atm = ACOS(1.0D0 / (sumtsza / n_coadd)) * rad2deg
        the_vza_atm = ACOS(1.0D0 / (sumtvza / n_coadd)) * rad2deg
        the_aza_atm = sumtaza / n_coadd
        sza_atm = ACOS(1.0D0 / (sumsza / n_coadd)) * rad2deg
        vza_atm = ACOS(1.0D0 / (sumvza / n_coadd)) * rad2deg
        aza_atm = sumaza / n_coadd
        aza_atm(2) = (aza_atm(1) + aza_atm(3)) / 2.0
        earth_curv = sumcurv / n_coadd
        ers2_alt =   sumalt  / n_coadd
        gome_geoloc = sumloc

        the_sca_atm = 180.0 - ACOS(COS(the_sza_atm * deg2rad) * COS(the_vza_atm &
             * deg2rad) + SIN(the_sza_atm * deg2rad)* SIN(the_vza_atm * deg2rad)&
             * COS(the_aza_atm * deg2rad)) * rad2deg
        
        WRITE(*, '(4x, A, I5)') 'nocadd = ', n_coadd
        !WRITE(*, '(4x, A5, 3f8.2)')  'SZA =', sza_atm
        !WRITE(*, '(4x, A5, 3f8.2)')  'VZA =', vza_atm
        !WRITE(*, '(4x, A5, 3f8.2)')  'AZA =', aza_atm        
        WRITE(*, '(4x, A, 4f8.2)') 'Effective SZA, VZA, AZA = ', &
             the_sza_atm, the_vza_atm, the_aza_atm, the_sca_atm       
     ENDIF
  ENDIF
  
  ! unable to be handled
  IF (gome_npix /= 1 .AND. gome_npix /= 8 .AND. gome_npix /= 40) THEN
     file_read_stat = file_read_failed; RETURN
  ENDIF

  nview = n_gome_ang
  nloc = n_gome_geo

  the_lats(1:nloc)  = gome_geoloc(lat_idx, 1:n_gome_geo)
  the_lons(1:nloc)  = gome_geoloc(lon_idx, 1:n_gome_geo)

  WHERE (the_lons(1:nloc) > 180.0) 
     the_lons(1:nloc) = the_lons(1:nloc) - 360.0
  ENDWHERE
  the_lon =  the_lons(nloc); the_lat = the_lats(nloc)   ! center cordinate
  
  ! Get left/right center of a pixel  
  IF (ABS(the_lons(1) - the_lons(2)) < 100.0) THEN   ! normal
     edgelons(1) = (the_lons(1) + the_lons(2)) / 2.0
  ELSE                                               ! across date line
     lon1 = the_lons(1); lon2 = the_lons(2)
     IF (lon1 > 0.0) THEN 
        lon2 = lon2 + 360.0
     ELSE
        lon1 = lon1 + 360.0
     ENDIF
     edgelons(1) = (lon1 + lon2) / 2.0
     IF (edgelons(1) > 180.0) edgelons(1) = edgelons(1) - 360.0 
  ENDIF
     
  IF (ABS(the_lons(3) - the_lons(4)) < 100.0) THEN   ! normal
     edgelons(2) = (the_lons(3) + the_lons(4)) / 2.0
  ELSE                                               ! across date line
     lon1 = the_lons(3); lon2 = the_lons(4)
     IF (lon1 > 0.0) THEN 
        lon2 = lon2 + 360.0
     ELSE
        lon1 = lon1 + 360.0
     ENDIF
     edgelons(2) = (lon1 + lon2) / 2.0
     IF (edgelons(2) > 180.0) edgelons(2) = edgelons(2) - 360.0 
  ENDIF

  edgelats(1) = (the_lats(1) + the_lats(2)) / 2.0
  edgelats(2) =  (the_lats(3) + the_lats(4)) / 2.0

  ! special handling for data inteprolation
  ! If pixel across date line, then edgelons [0-360], otherwise [-180W, 180E]
  IF (ABS(edgelons(1) - edgelons(2)) > 100.0) THEN                                            
     IF (edgelons(1) > 0.0) THEN 
        edgelons(2) = edgelons(2) + 360.0
     ELSE
        edgelons(1)= edgelons(1) + 360.0
     ENDIF
  ENDIF

  ! make sure that edgelon(1) < edgelon(2)
  IF (edgelons(1) > edgelons(2)) THEN
     templon = edgelons(2)
     edgelons(2) = edgelons(1)
     edgelons(1) = templon
     
     templat = edgelats(2)
     edgelats(2) = edgelats(1)
     edgelats(1) = templat
  ENDIF
  
  year = gome_pixdate(8:11);  READ(year,'(I4)') the_year
  day =  gome_pixdate (1:2);  READ(day, '(I2)') the_day
  monthc = gome_pixdate(4:6)
  
  DO i = 1, 12
     IF (monthc == months(i)) EXIT
  ENDDO
  the_month = i
    
  CALL gome_read_el1radspec (funit, numwin, winpix(1:numwin, :), &
       winlim(1:numwin, :), band_selectors(1:numwin), &
       n_gome_radpts, gome_radspec, nradpix(1:numwin), file_read_stat)
  IF (file_read_stat /= file_read_ok ) RETURN
  
  spec_norm = SUM(gome_radspec(2,1:n_gome_radpts)) / n_gome_radpts
  IF (spec_norm <= 0.0 ) THEN
     error = .TRUE.;  RETURN
  ENDIF

  gome_radspec(2,1:n_gome_radpts) = gome_radspec(2,1:n_gome_radpts) / spec_norm
  gome_radspec(3,1:n_gome_radpts) = gome_radspec(3,1:n_gome_radpts) / spec_norm
  div_rad = spec_norm

  !OPEN(77, file='sim_noise.dat', status='old')
  !DO i = 1, n_gome_radpts
  !   read(77, *) j, temp
  !   gome_radspec(1,i) = gome_radspec(1,i) + temp * 0.004
  !ENDDO
  !CLOSE(77)

  RETURN

END SUBROUTINE gome_read_el1data_rad


SUBROUTINE gome_read_el1geo ( &
     funit, gome_pixdate, gome_angles_wrtn, gome_angles_wrts, gome_geoloc, &
     ers2_alt, earth_curv, file_read_stat)

  USE OMSAO_precision_module
  USE OMSAO_errstat_module,     ONLY : file_read_ok
  USE OMSAO_gome_data_module,   ONLY : zen0_idx, azm0_idx, &
       zen_idx, azm_idx, lat_idx, lon_idx, n_gome_ang, n_gome_geo

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER, INTENT (IN) :: funit

  ! ================
  ! Output variables
  ! ================
  INTEGER,           INTENT (OUT) :: file_read_stat
  CHARACTER (LEN=*), INTENT (OUT) :: gome_pixdate
  REAL (KIND=dp),    INTENT (OUT) :: ers2_alt, earth_curv
  REAL (KIND=dp), DIMENSION (azm_idx,n_gome_ang), INTENT (OUT) :: gome_angles_wrtn, &
       gome_angles_wrts
  REAL (KIND=dp), DIMENSION (lon_idx,n_gome_geo), INTENT (OUT) :: gome_geoloc

  ! ===============
  ! Local variables
  ! ===============
  INTEGER :: ios, k
  REAL (KIND=dp), DIMENSION (2*n_gome_ang) :: tmp_arr
  REAL (KIND=dp), DIMENSION (2*n_gome_geo) :: tmp_arr2

  ! ---------------
  ! Read pixel date
  ! ---------------
  READ (UNIT=funit, FMT='(A)', IOSTAT=ios) gome_pixdate

  CALL gome_check_read_status ( ios, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) RETURN

  ! ---------------------------------------------------------------------
  ! Read GOME pixel viewing angles:
  ! Solar and line-of-sight angles with respect to North or to Spacecraft
  ! ---------------------------------------------------------------------
  READ (UNIT=funit, FMT=*, IOSTAT=ios) tmp_arr(1:2*n_gome_ang)
  DO k = 1, n_gome_ang
     gome_angles_wrtn(zen0_idx,k) = tmp_arr(2*k-1)
     gome_angles_wrtn(azm0_idx,k) = tmp_arr(2*k)
  END DO
  CALL gome_check_read_status ( ios, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) RETURN

  READ (UNIT=funit, FMT=*, IOSTAT=ios) tmp_arr(1:2*n_gome_ang)
  DO k = 1, n_gome_ang
     gome_angles_wrtn(zen_idx, k) = tmp_arr(2*k-1)
     gome_angles_wrtn(azm_idx, k) = tmp_arr(2*k)
  END DO
  READ (UNIT=funit, FMT=*, IOSTAT=ios) tmp_arr(1:2*n_gome_ang)
  DO k = 1, n_gome_ang
     gome_angles_wrts(zen0_idx,k) = tmp_arr(2*k-1)
     gome_angles_wrts(azm0_idx,k) = tmp_arr(2*k)
  END DO
  CALL gome_check_read_status ( ios, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) RETURN

  READ (UNIT=funit, FMT=*, IOSTAT=ios) tmp_arr(1:2*n_gome_ang)
  DO k = 1, n_gome_ang
     gome_angles_wrts(zen_idx, k) = tmp_arr(2*k-1)
     gome_angles_wrts(azm_idx, k) = tmp_arr(2*k)
  END DO
  CALL gome_check_read_status ( ios, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) RETURN

  ! ---------------------------------------
  ! Read ERS-2 Altitude and Earth Curvature
  ! ---------------------------------------
  READ (UNIT=funit, FMT=*, IOSTAT=ios) ers2_alt, earth_curv
  CALL gome_check_read_status ( ios, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) RETURN

  ! -------------------------------------
  ! Read LAT/LON specifications for pixel
  ! -------------------------------------
  READ (UNIT=funit, FMT=*, IOSTAT=ios) tmp_arr2(1:2*n_gome_geo)
  DO k = 1, n_gome_geo
      gome_geoloc(lat_idx,k) = tmp_arr2(2*k-1)
      gome_geoloc(lon_idx,k) = tmp_arr2(2*k)
  END DO

  CALL gome_check_read_status ( ios, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) RETURN

  RETURN
END SUBROUTINE gome_read_el1geo


SUBROUTINE gome_read_el1solspec (funit, nwin, winpix, winlim, bds, &
     ngome, gspec, nsolpix, file_read_stat)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : maxchlen, maxwin, max_fit_pts, mrefl
  USE OMSAO_errstat_module,    ONLY : file_read_ok, file_read_failed, file_read_eof
  USE OMSAO_gome_data_module,  ONLY : n_gome_data_dim, n_gome_max_pts
  USE OMSAO_variables_module,  ONLY : sol_spec_ring, nsol_ring, sol_identifier
  USE ozprof_data_module,      ONLY : sun_posr, sun_specr, degcorr, &
       degfname, biascorr, biasfname, corr_unit, nrefl

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,        INTENT (IN)                          :: funit, nwin
  INTEGER, DIMENSION(nwin, 2), INTENT(IN)              :: winpix
  INTEGER, DIMENSION(nwin), INTENT(IN)                 :: bds


  ! ================
  ! Output variables
  ! ================
  INTEGER, INTENT (OUT)                 :: file_read_stat, ngome
  INTEGER, DIMENSION(nwin), INTENT(OUT) :: nsolpix
  REAL (KIND=dp), DIMENSION(nwin, 2), INTENT(OUT) :: winlim
  REAL (KIND=dp), DIMENSION(n_gome_data_dim, max_fit_pts), INTENT(OUT):: gspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER, PARAMETER          :: maxw = maxwin - 1, nuar = 2911
  REAL (KIND=dp), DIMENSION (maxw, n_gome_max_pts, n_gome_data_dim)  :: solspec
  INTEGER, DIMENSION (maxw)   :: numtpix = (/695, 841, 1024, 1024/), numpix
  INTEGER                     :: ios, i, j, k, totpix, ch, fch, lch, iwin, &
       fidx, lidx, sfidx, slidx, errstat, sidx, eidx
  REAL (KIND=dp)              :: swav, ewav, albcoe, factor, prewav
  CHARACTER (LEN=9)           :: tmpchar
  CHARACTER (LEN=maxchlen)    :: lastline

  INTEGER, PARAMETER          :: maxdeg = 500
  INTEGER                     :: ndeg, nwave, yr, mon, day, utc, ntot, nord
  REAL (KIND=dp)              :: curtm, frac
  INTEGER, DIMENSION(12)      :: ndays = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
  REAL (KIND=dp), DIMENSION(maxdeg)                 :: degtms
  REAL (KIND=dp), DIMENSION(n_gome_max_pts)         :: coeff
  REAL (KIND=dp), DIMENSION(10, n_gome_max_pts)     :: del
  REAL (KIND=dp), DIMENSION(maxdeg, 20)             :: corr
  REAL (KIND=dp), DIMENSION(20)                     :: bcorr
  REAL (KIND=dp), DIMENSION(nwin)                   :: woff
  
  
  ! ----------------------------------
  ! Initialization of output variables
  ! ----------------------------------
  file_read_stat = file_read_ok
  
  ! Determine the first and last channel to read
  fch = bds(1)
  IF (winpix(1, 1) - 40 < 0) fch = fch-1
  fch = MAX(1, fch)
  
  lch = bds(nwin)
  IF (winpix(nwin, 2) + 40 > numtpix(lch) ) lch = lch + 1
  lch = MAX(lch, 2)  ! must have channel for surface albedo
  lch = MIN(lch, 4)
  
  ! Read all the data, channel 1, 2, 3, 4
  numpix = 0; solspec = 0.0
  DO i = 1, maxw
     CALL skip_to_filemark ( funit, 'Channel', lastline, file_read_stat )
     IF ( file_read_stat /= file_read_ok ) THEN
        STOP 'Read Solar: Cannot find header!!!'
     END IF
     
     READ(lastline, *) tmpchar, ch, swav, ewav, totpix
     !IF (ch < fch) CYCLE
     numpix(ch) = totpix

     DO j = 1, totpix
        READ (UNIT=funit, FMT=*, IOSTAT = ios) solspec(ch, j, 1:n_gome_data_dim)
        IF (ios /= 0 ) THEN
           WRITE(*, *) 'Read Solar: Read Error: Channel', ch
           STOP
        ENDIF
     ENDDO

     IF (ch >= lch) EXIT
  ENDDO

  ! Get data for fitting
  prewav = 0.0
  ngome = 0; gspec = 0.0
  DO i = 1, nwin
     ch = bds(i)
     IF (numpix(ch) == 0) THEN
        WRITE(*, *) 'No pixels are read from channel ', ch; STOP
     ENDIF
     nsolpix(i) = 0
     IF (ch /= 1) THEN
        DO j = winpix(i, 1), winpix(i, 2)
           IF (solspec(ch, j, 2) /= 0.0 .AND. solspec(ch, j, 1) &
                > prewav + 0.1 .AND. solspec(ch, j, 5) == 0) THEN
              ngome = ngome + 1; nsolpix(i) = nsolpix(i) + 1
              gspec(:, ngome) = solspec(ch, j, :)
              prewav = gspec(1, ngome)
           ENDIF
        ENDDO
     ELSE
        DO j = winpix(i, 1), winpix(i, 2)
!           IF (solspec(ch, j, 2) /= 0.0 .AND.  j >= 231 &  ! > 264.0 nm
!                .AND. (j >= 393 .OR. j <= 276)          &  ! no 269.0-282.0 nm 
!                .AND. (j >= 456 .OR. j <= 411)          &  ! no 284-289 nm  
!                .AND. solspec(ch, j, 1) > prewav) THEN
           IF (solspec(ch, j, 2) /= 0.0 .AND.  j >= 231 &  ! > 264.0 nm
                .AND. (j >= 410 .OR. j <= 300)          &  ! no 276.0-282.0 nm 
                .AND. (j >= 448 .OR. j <= 410)          &  ! no 284-289 nm  
                .AND. solspec(ch, j, 1) > prewav + 0.1 &
                .AND. solspec(ch, j, 5) == 0) THEN      
              ngome = ngome + 1; nsolpix(i) = nsolpix(i) + 1
              gspec(:, ngome) = solspec(ch, j, :)
              prewav = gspec(1, ngome)
           ENDIF
        ENDDO
     ENDIF
     IF (i == 1) THEN
        winlim(i, 1) = gspec(1,1) 
     ELSE 
        winlim(i, 1) = gspec(1, SUM(nsolpix(1:i-1))+1)
     ENDIF

     winlim(i, 2) = gspec(1, ngome)
  ENDDO
  IF (ngome > max_fit_pts) THEN
     WRITE(*, *) 'Read Solar: Increase max_fit_pts!!!'; STOP
  ELSE IF (ngome == 0) THEN
     WRITE(*, *) 'Read Radiance: NO radiances are read!!!'; STOP
  ENDIF

  ! get data for ring calculation
  nsol_ring = nsolpix(1) + 40
  sol_spec_ring(1:2, 41:41+nsolpix(1)) = gspec(1:2, 1:nsolpix(1))

  j = 41
  DO 
     ch = bds(1)
     DO i = winpix(ch, 1)-1, 1, -1
        IF (solspec(ch, i, 2) /= 0.0 .AND. solspec(ch, i, 1) &
             < sol_spec_ring(1, j) - 0.1 .AND. solspec(ch, i, 5) == 0 ) THEN
           j = j - 1
           sol_spec_ring(1:2, j) = solspec(ch, i, 1:2)
           IF (j == 1) EXIT
        ENDIF
     ENDDO
     IF (j == 1) EXIT

     ch = ch - 1
     IF (ch < 1) EXIT

     DO i = numpix(ch), 1, -1
        IF (solspec(ch, i, 2) /= 0.0 .AND. solspec(ch, i, 1) &
             < sol_spec_ring(1, j) - 0.1 .AND. solspec(ch, i, 5) == 0) THEN
           j = j - 1
           sol_spec_ring(1:2, j) = solspec(ch, i, 1:2)
           IF (j == 1) EXIT
        ENDIF
     ENDDO
     IF (j == 1) EXIT
  ENDDO

  iwin = 1; ch = bds(iwin); i = winpix(1, 2) + 1
  DO 
     IF (solspec(ch, i, 2) /= 0.0 .AND. solspec(ch, i, 1) > &
          sol_spec_ring(1, nsol_ring) .AND. solspec(ch, i, 1) <= gspec(1, ngome) &
          .AND. i >= winpix(iwin, 1) - 40 .AND. i <= winpix(iwin, 2) + 40 &
          .AND. solspec(ch, i, 5) == 0) THEN
        IF (iwin /= nwin) THEN
           IF (solspec(ch, i, 1) < solspec(bds(iwin+1), winpix(iwin+1, 1), 1)) THEN
              nsol_ring = nsol_ring + 1
              sol_spec_ring(1:2, nsol_ring) = solspec(ch, i, 1:2)
           ENDIF
        ELSE IF (iwin == nwin) THEN
           nsol_ring = nsol_ring + 1
           sol_spec_ring(1:2, nsol_ring) = solspec(ch, i, 1:2)
        ENDIF
     ENDIF
     IF (i == numtpix(ch) .AND. iwin < nwin) THEN
        iwin = iwin + 1; ch = bds(iwin); i=1
     ELSE IF ( i == numtpix(ch)) THEN
        EXIT
     ELSE
        i = i + 1
     ENDIF
  ENDDO

  j = nsol_ring
  nsol_ring = nsol_ring + 40
  DO 
     ch = bds(nwin)
     DO i = winpix(nwin, 2)+1, numpix(ch)
        IF (solspec(ch, i, 2) /= 0.0 .AND. solspec(ch, i, 1)&
             > sol_spec_ring(1, j) + 0.1 .AND. solspec(ch, i, 5) == 0) THEN
           j = j + 1
           sol_spec_ring(1:2, j) = solspec(ch, i, 1:2)
           IF (j == nsol_ring) EXIT
        ENDIF
     ENDDO
     IF (j == nsol_ring) EXIT
     
     ch = ch + 1
     IF (ch >= maxw) EXIT

     DO i = 1, numpix(ch)
        IF (solspec(ch, i, 2) /= 0.0 .AND. solspec(ch, i, 1) &
             > sol_spec_ring(1, j) + 0.1 .AND. solspec(ch, i, 5) == 0) THEN
           j = j + 1
           sol_spec_ring(1:2, j) = solspec(ch, i, 1:2)
           IF (j == nsol_ring) EXIT
        ENDIF
     ENDDO
     IF (j == nsol_ring) EXIT
  ENDDO

  IF (nsol_ring >= max_fit_pts) THEN
     WRITE(*, *) 'Read Solar: Number of points in sol_spec_ring exceeds max_fit_pts!!!'
     STOP
  ENDIF
  
  ! Get data for surface albedo at 370.2 nm +/- 20 pixels
  j = 0
  DO i = 500, numpix(2)
     IF (solspec(2, i, 2) /= 0.0 .AND. solspec(2, i, 5) == 0) THEN
        j = j + 1
        sun_posr(j) = solspec(2, i, 1); sun_specr(j) = solspec(2, i, 2)
        IF (j == nrefl) EXIT
     ENDIF
  ENDDO

  WRITE(*, *) 'End Of Reading Irradiance Spectrum'
  DO i = 1, nwin
     WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', i, winlim(i,1), winlim(i,2), nsolpix(i)
     IF (nsolpix(i) <= 20) THEN
        WRITE(*, '(A,f8.3,A3,f8.3)') ' Not enough points in window: ', winlim(i,1),&
             ' - ', winlim(i,2)
        STOP
     ENDIF
  ENDDO

  ! check the units of GOME solar spectrum
  IF (gspec(2, 1) .LE. 1.0E5) THEN    ! in mw/m^2/nm
     factor = 5.0341176e+08
     gspec(2, 1:ngome) = gspec(2, 1:ngome) * factor * gspec(1, 1:ngome)
     gspec(3, 1:ngome) = gspec(3, 1:ngome) * factor * gspec(1, 1:ngome)
     sun_specr(1:nrefl) = sun_specr(1:nrefl) * factor * sun_posr(1:nrefl)
     sol_spec_ring(2, 1:nsol_ring) = sol_spec_ring(2, 1:nsol_ring) * &
          factor * sol_spec_ring(1, 1:nsol_ring)
  ENDIF

  IF (biascorr) THEN
     IF (nwin /=2) THEN
        WRITE(*, *) 'Bias LUT is not ready for nwin = ', nwin
        STOP
     ENDIF

     OPEN(UNIT=corr_unit, FILE=ADJUSTL(TRIM(biasfname)), STATUS='OLD')
     READ(corr_unit, *) 
     READ(corr_unit, *) 
     READ(corr_unit, *)
     READ(corr_unit, *) woff
     READ(corr_unit, *) nord     
     READ(corr_unit, *) (bcorr(i), i = 1, 2 * nord)   
     CLOSE(corr_unit) 

     sidx = 1
     DO i = 1, nwin 
        j = (i - 1) * nord + 1
        eidx = sidx + nsolpix(i) - 1
        del(1, 1:nsolpix(i)) = gspec(1, sidx:eidx) - woff(i)
        DO k = 2, nord-1
           del(k, 1:nsolpix(i)) = del(k-1, 1:nsolpix(i)) * del(1, 1:nsolpix(i))
        ENDDO
     
        coeff(1:nsolpix(i)) = bcorr(j)
        DO k = 1, nord - 1
           coeff(1:nsolpix(i)) =  coeff(1:nsolpix(i)) + bcorr(j + k) * del(k, 1:nsolpix(i)) 
        ENDDO
              
        gspec(2, sidx:eidx) = gspec(2, sidx:eidx) * coeff(1:nsolpix(i))
        gspec(3, sidx:eidx) = gspec(3, sidx:eidx) * coeff(1:nsolpix(i))
        sidx = eidx + 1
     ENDDO
  ENDIF

  IF (degcorr) THEN
     IF (nwin /= 2) THEN
        WRITE(*, *) 'Degradation LUT is not ready for nwin = ', nwin
        STOP
     ENDIF

     OPEN(UNIT=corr_unit, FILE=ADJUSTL(TRIM(degfname)), STATUS='OLD')
     READ(corr_unit, *) 
     READ(corr_unit, *) ndeg
     IF (ndeg > maxdeg) THEN 
        WRITE(*, *) 'Need to increase maxdeg to at least: ', ndeg;  STOP
     ENDIF
     READ(corr_unit, *) degtms(1:ndeg)
     READ(corr_unit, '(A)')
     READ(corr_unit, *) woff
     READ(corr_unit, *) nord
     READ(corr_unit, '(A)')
     READ(corr_unit, *) ((corr(i, j), j = 1, 2 * nord+1), i = 1, ndeg)     
     CLOSE(corr_unit) 
   
     READ(sol_identifier, '(I1, I2, I2, I3)') yr, mon, day, utc
     IF (yr >= 5) THEN 
        yr = yr + 1990
     ELSE 
        yr = yr + 2000
     ENDIF
     IF (yr >= 2004 .OR. (yr == 2003 .AND. mon >= 6)) THEN 
        WRITE(*, *) 'Degradation cannot be performed since June 2003!!!'; STOP
     ENDIF

     IF ( MOD(yr, 4) == 0) THEN
        ntot = 366;   ndays(2) = 29
     ELSE
        ntot = 365;  ndays(2) = 28
     ENDIF
     
     curtm = yr + 1.0 * (SUM(ndays(1:mon-1)) + day - 1 + utc / 240.0) / ntot

     ! find two days that bracket curtm
     lidx = MINVAL(MINLOC(degtms(1:ndeg), MASK = (degtms(1:ndeg) >= curtm)))
     IF (lidx > 1) THEN
        fidx = lidx - 1
        frac = 1.0 - (curtm - degtms(fidx)) / (degtms(lidx) - degtms(fidx))
     ELSE
        lidx = 2; fidx = 1
        frac = 1.0
     ENDIF
     
     sidx = 1
     DO i = 1, nwin 
        j = (i - 1) * nord + 1
        eidx = sidx + nsolpix(i) - 1
        
        del(1, 1:nsolpix(i)) = gspec(1, sidx:eidx) - woff(i)
        DO k = 2, nord - 1
           del(k, 1:nsolpix(i)) = del(k-1, 1:nsolpix(i)) * del(1, 1:nsolpix(i))
        ENDDO
        
        coeff(1:nsolpix(i)) = frac * corr(fidx, j) + (1.0 - frac) * corr(lidx, j)
        DO k = 1, nord - 1
           coeff(1:nsolpix(i)) = coeff(1:nsolpix(i)) + (frac * corr(fidx, j + k) + &
                (1.0 - frac) * corr(lidx, j + k)) * del(k, 1:nsolpix(i))
        ENDDO
        
        gspec(2, sidx:eidx) = gspec(2, sidx:eidx) * coeff(1:nsolpix(i))
        gspec(3, sidx:eidx) = gspec(3, sidx:eidx) * coeff(1:nsolpix(i))
        sidx = eidx + 1
     ENDDO

     j = 2 * nord + 1
     coeff(1) = frac * corr(fidx, j) + (1.0 - frac) * corr(lidx, j)
     sun_specr(1:nrefl) = sun_specr(1:nrefl) * coeff(1)  
  ENDIF
  
  CALL gome_check_read_status ( ios, file_read_stat)
  
  RETURN
END SUBROUTINE gome_read_el1solspec


SUBROUTINE gome_read_el1radspec (funit, nwin, winpix, winlim, bds, &
     ngome, gspec, nradpix, file_read_stat)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : maxchlen, maxwin, max_fit_pts, mrefl
  USE OMSAO_errstat_module,    ONLY : file_read_ok, file_read_failed, &
       file_read_eof, file_read_missing
  USE OMSAO_gome_data_module,  ONLY : n_gome_data_dim, n_gome_max_pts
  USE OMSAO_variables_module,  ONLY : b1ab_div_wav, nsolpix
  USE ozprof_data_module,      ONLY : rad_posr, rad_specr, nrefl

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,        INTENT (IN)             :: funit, nwin
  INTEGER, DIMENSION(nwin, 2), INTENT(IN) :: winpix
  INTEGER, DIMENSION(nwin), INTENT(IN)    :: bds


  ! ================
  ! Output variables
  ! ================
  INTEGER, INTENT (OUT)                 :: file_read_stat, ngome
  INTEGER, DIMENSION(nwin), INTENT(OUT) :: nradpix
  REAL (KIND=dp), DIMENSION(nwin, 2), INTENT(INOUT) :: winlim
  REAL (KIND=dp), DIMENSION(n_gome_data_dim, max_fit_pts), INTENT(OUT):: gspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER, PARAMETER          :: maxw = maxwin - 1
  REAL (KIND=dp), DIMENSION (maxw, n_gome_max_pts, n_gome_data_dim)  :: radspec
  INTEGER, DIMENSION (maxw) :: numtpix = (/695, 841, 1024, 1024/), numpix
  INTEGER                     :: ios, i, j, totpix, ch, fch, lch, offset
  REAL (KIND=dp)              :: integt, swav, ewav, prewav
  CHARACTER (LEN=9)           :: tmpchar
  CHARACTER (LEN=5)           :: tmpchar1
  CHARACTER (LEN=1)           :: subch
  CHARACTER (LEN=maxchlen)    :: lastline
 
  ! ----------------------------------
  ! Initialization of output variables
  ! ----------------------------------
  file_read_stat = file_read_ok
  
  ! Determine the first and last channel to read
  fch = bds(1)
  IF (winpix(1, 1) - 20 < 0) fch = fch-1
  fch = MAX(1, fch)
  
  lch = bds(nwin)
  IF (winpix(nwin, 2) + 20 > numtpix(lch) ) lch = lch + 1
  lch = MAX(lch, 2)  ! must have channel for surface albedo
  lch = MIN(lch, 4)

  b1ab_div_wav = 307.2
 
  ! Read all the data 
  radspec = 0.0; numpix = 0
  DO i = 1, maxw + 2   ! +2 because of 1a, 1b, 2a, 2b
     CALL skip_to_filemark ( funit, 'Band', lastline, file_read_stat )
     IF ( file_read_stat /= file_read_ok ) THEN
        STOP 'Read Radiance: Cannot find header!!!'
     END IF
     
     READ(lastline, '(A9,F9.5,2F9.4,I5)') tmpchar, integt, swav, ewav, totpix
     READ(tmpchar, '(A5,I1,A1)') tmpchar1, ch, subch

     IF (ch == 1 .AND. subch == 'b') b1ab_div_wav = swav - 0.01
     !IF (ch < fch) CYCLE

     IF ( (ch == 1 .OR. ch == 2) .AND. subch == 'b' ) THEN
        offset = numtpix(ch) - totpix
     ELSE
        offset = 0
     ENDIF
     numpix(ch) = numpix(ch) + totpix
        
     DO j = 1, totpix
        READ (UNIT=funit, FMT=*, IOSTAT = ios) &
             radspec(ch, offset+j, 1:n_gome_data_dim)
        IF (ios /= 0 ) THEN
           WRITE(*, *) 'Read Radiance: Read Error: Channel', ch
           STOP
        ENDIF
     ENDDO

     IF ( (ch >= lch) .AND. (subch == 'b' .OR. subch == ' ')) EXIT
  ENDDO

  ! Get data for fitting
  ngome = 0 ; gspec = 0.0
  prewav = 0.0
  DO i = 1, nwin
     ch = bds(i)
     IF (numpix(ch) == 0) THEN
        WRITE(*, *) 'No pixels are read from channel ', ch; STOP
     ENDIF
     nradpix(i) = 0
     IF (ch /= 1) THEN
        DO j = winpix(i, 1), winpix(i, 2)
           IF (radspec(ch, j, 2) /= 0.0 .AND. (radspec(ch, j, 1) &
                > prewav + 0.1 .OR. ngome == 0) .AND. radspec(ch, j, 5) == 0) THEN
              ngome = ngome + 1; nradpix(i) = nradpix(i) + 1
              gspec(:, ngome) = radspec(ch, j, :)
              prewav = gspec(1, ngome)
           ENDIF
        ENDDO
     ELSE
        DO j = winpix(i, 1), winpix(i, 2)
!           IF (radspec(ch, j, 2) /= 0.0 .AND.  j >= 231 &  ! > 264.0 nm
!                .AND. (j >= 393 .OR. j <= 276)          &  ! no 269.0-282.0 nm 
!                .AND. (j >= 456 .OR. j <= 411)          &  ! no 284-289 nm
!                .AND. radspec(ch, j, 1) > gspec(1, ngome) + 0.1) THEN
           IF (radspec(ch, j, 2) /= 0.0 .AND.  j >= 231 &  ! > 264.0 nm
                .AND. (j >= 410 .OR. j <= 300)          &  ! no 276.0-282.0 nm 
                .AND. (j >= 448 .OR. j <= 410)          &  ! no 284-289 nm
                .AND. (ngome == 0 .OR. radspec(ch, j, 1) > prewav + 0.1) &
                .AND. radspec(ch, j, 5) == 0) THEN
              ngome = ngome + 1; nradpix(i) = nradpix(i) + 1
              gspec(:, ngome) = radspec(ch, j, :)
              prewav = gspec(1, ngome)
           ENDIF
        ENDDO
     ENDIF
     IF (i == 1) THEN
        winlim(i, 1) = gspec(1,1) 
     ELSE
        winlim(i, 1) = gspec(1, SUM(nradpix(1:i-1))+1)
     ENDIF
     winlim(i, 2) = gspec(1, ngome)
  ENDDO
  IF (ngome > max_fit_pts) THEN
     WRITE(*, *) 'Read Radiance: Increase max_fit_pts!!!'; STOP
  ELSE IF (ngome == 0) THEN
     WRITE(*, *) 'Read Radiance: NO radiances are read!!!'; STOP
  ENDIF
 
  ! Get data for surface albedo at 370.2 nm +/- 20 pixels
  j = 0
  DO i = 500, numpix(2)
     IF (radspec(2, i, 2) /= 0.0 .AND. radspec(2, i, 5) == 0) THEN
        j = j + 1
        rad_posr(j) = radspec(2, i, 1); rad_specr(j) = radspec(2, i, 2)
        IF (j == nrefl) EXIT
     ENDIF
  ENDDO

  WRITE(*, *) 'End of Reading Radiance Spectrum'
  DO i = 1, nwin
     WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', i, winlim(i,1), winlim(i,2), nradpix(i)
     IF (nradpix(i) <= 20) THEN
        WRITE(*, '(A,f8.3,A3,f8.3)') 'Not enough points in window: ', winlim(i,1),&
             ' - ', winlim(i,2)
        STOP
     ENDIF
     
     IF (nradpix(i) > nsolpix(i)) THEN
        WRITE(*, *) 'Inconsistent pixels between I and F0 in channel ', bds(i)
        file_read_stat = file_read_missing; RETURN
     ENDIF
  ENDDO
  
  CALL gome_check_read_status ( ios, file_read_stat)
  
  RETURN
END SUBROUTINE gome_read_el1radspec

! Purpose: read and average data for ozone profile retrieval
! 1. Read channel 1 and all channel 2 pixels corresponding to channel 1
! 2. Average channel 2 pixels to channel 1 spatial resolution

SUBROUTINE gome_read_el1data_radavg ( funit, file_read_stat, error )

  USE OMSAO_precision_module
  USE OMSAO_errstat_module,    ONLY : &
       file_read_ok, file_read_failed, file_read_missing, file_read_eof
  USE OMSAO_parameters_module, ONLY : maxchlen, deg2rad, rad2deg
  USE OMSAO_gome_data_module,  ONLY : &
       azm_idx, zen_idx, zen0_idx, lat_idx, lon_idx, n_gome_ang, n_gome_geo,   &
       lm_gome_groundpix, gome_spec_missing, n_gome_max_pts, n_gome_data_dim,  &
       n_gome_radpts, gome_radspec, gome_curpix, gome_curqual, gome_curscan,   &
       gome_pixdate, gome_angles_wrtn, gome_angles_wrts, gome_geoloc, ers2_alt,&
       earth_curv, orbnum, gome_stpix, gome_endpix, gome_npix
  USE ozprof_data_module,      ONLY : div_rad, coadd_after_b1ab, b1ab_change
  USE OMSAO_variables_module,  ONLY : the_month, the_year, the_day, numwin,    &
       winlim, winpix, band_selectors, nradpix, npix_fitting, sza_atm, vza_atm,&
       aza_atm, the_sza_atm, the_aza_atm, the_vza_atm, the_sca_atm, the_lat,   &
       the_lon, the_lons, the_lats, nview, nloc, edgelons, edgelats, b1ab_div_wav, &
       rad_identifier, szamax

  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER,  INTENT (IN)    :: funit
  LOGICAL, INTENT (OUT)    :: error	
  INTEGER, INTENT (OUT)    :: file_read_stat

  ! ===============
  ! Local variables
  ! ===============
  INTEGER, PARAMETER       :: maxnpix = 40, nband = 2
  INTEGER                  :: i, j, n_coadd
  LOGICAL                  :: view_coadd, valid_pix
  CHARACTER (LEN=maxchlen) :: lastline
  CHARACTER (LEN=3)        :: monthc
  CHARACTER (LEN=2)        :: day
  CHARACTER (LEN=4)        :: year
  REAL (KIND=dp)           :: spec_norm, lon1, lon2, templon, templat
  REAL (KIND=dp), DIMENSION(n_gome_ang) :: sumsza, sumaza, sumvza
  REAL (KIND=dp)                        :: sumtsza, sumtaza, sumtvza, sumalt, sumcurv
  REAL (KIND=dp), DIMENSION (lon_idx, n_gome_geo) :: sumloc
  CHARACTER (LEN=3), DIMENSION(12), PARAMETER    :: months = (/'JAN', 'FEB', 'MAR', &
       'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/)
  ! 2: channel 1 and channel 2
  REAL (KIND=dp), DIMENSION (maxnpix, nband, n_gome_max_pts, n_gome_data_dim) :: radspec
  INTEGER,        DIMENSION (maxnpix, nband)                             :: nwpos, stpos

  file_read_stat = file_read_ok ; error = .FALSE.
 
  IF ((band_selectors(1) > 1) .OR. (.NOT. coadd_after_b1ab .AND. b1ab_change)) THEN
     view_coadd = .FALSE.
  ELSE
     view_coadd = .TRUE.
  ENDIF

  IF (npix_fitting == 0 .AND. orbnum == 0) THEN
     CALL skip_to_filemark ( funit, 'ERS Information', lastline, file_read_stat )
     IF ( file_read_stat /= file_read_ok ) RETURN
     READ (funit, *), i, i, orbnum
  ENDIF

  IF (.NOT. view_coadd) THEN
     ! ----------------------------------------------------------------
     ! Position cursor at start of GOME ground pixel entry
     ! And Extract GOME ground pixel number and scan type from LASTLINE
     ! -----------------------------------------------------------------
     CALL skip_to_filemark ( funit, lm_gome_groundpix, lastline, file_read_stat)
     IF ( file_read_stat /= file_read_ok ) RETURN
     READ (lastline(14:17), '(I4)') gome_curpix
     READ (lastline(19:19), '(I1)') gome_curqual
     READ (lastline(21:21), '(I1)') gome_curscan

     ! No retrieval for back scan
     IF ( gome_curqual == gome_spec_missing) THEN
        file_read_stat = file_read_missing; RETURN
     END IF

     CALL gome_read_el1geo (funit, gome_pixdate, gome_angles_wrtn, &
          gome_angles_wrts, gome_geoloc, ers2_alt, earth_curv, file_read_stat)  
     IF ( ANY(gome_angles_wrtn(zen0_idx, 1:n_gome_ang) >= szamax) ) THEN
        file_read_stat = file_read_missing; RETURN
     END IF  

     CALL calc_view_geometry
     IF ( the_sza_atm < 0.0 .OR. the_sza_atm >= szamax) THEN
        file_read_stat = file_read_missing; RETURN
     END IF

     ! correct for scanning asymmetry
     IF (gome_curscan == 3) THEN
        the_vza_atm = the_vza_atm + 1.955
     ELSE IF (gome_curscan == 0 .OR. gome_curscan == 2) THEN
        the_vza_atm = the_vza_atm - 1.23
     ENDIF
     
     gome_stpix = gome_curpix; gome_endpix = gome_curpix; gome_npix = 1
     CALL gome_read_el1radspec (funit, numwin, winpix(1:numwin, :), &
          winlim(1:numwin, :), band_selectors(1:numwin), &
          n_gome_radpts, gome_radspec, nradpix(1:numwin), file_read_stat)
     IF (file_read_stat /= file_read_ok ) RETURN

     ! No retrieval for backscan (geographically overlapped)
     IF ( gome_curscan == 3) THEN
        file_read_stat = file_read_missing; RETURN
     END IF
  ELSE
     n_coadd =   0;    sumsza=   0.0;   sumvza  = 0.0;   sumaza = 0.0
     sumtaza = 0.0;    sumtvza = 0.0;   sumtsza = 0.0
     sumalt =  0.0;    sumcurv = 0.0

     DO
        CALL skip_to_filemark (funit, lm_gome_groundpix, lastline, file_read_stat)
        IF ( file_read_stat /= file_read_ok ) RETURN
        READ (lastline(14:17), '(I4)') gome_curpix
        READ (lastline(19:19), '(I1)') gome_curqual
        READ (lastline(21:21), '(I1)') gome_curscan

        CALL gome_read_el1geo (funit, gome_pixdate, gome_angles_wrtn, &
             gome_angles_wrts, gome_geoloc, ers2_alt, earth_curv, file_read_stat)
        IF ( ANY(gome_angles_wrtn(zen0_idx, 1:n_gome_ang) >= szamax) ) THEN
           file_read_stat = file_read_missing; RETURN
        END IF

        ! First pixel:
        ! 1. scan = 0  2. sza le 90.0  3. # bands=3 (b1b,b2a,b2b) 4: n_coadd=0

        valid_pix = .FALSE.
        IF (gome_curqual >= 3 .AND. n_coadd <40 .AND. &
             ALL(gome_angles_wrtn(zen0_idx, 1:n_gome_ang) < szamax) ) THEN
        
           IF (n_coadd == 0 .AND. gome_curscan == 0) THEN
              n_coadd = n_coadd + 1; gome_stpix = gome_curpix; valid_pix = .TRUE.
           ELSE IF (n_coadd > 0) THEN
              n_coadd = n_coadd + 1; valid_pix = .TRUE.
           ENDIF
           ! ELSE: have not found the first pixel yet
           
           IF (valid_pix) THEN
              CALL calc_view_geometry
              IF ( the_sza_atm < 0.0 .OR. the_sza_atm >= szamax) THEN
                 file_read_stat = file_read_missing; RETURN
              END IF
              IF (gome_curscan == 3) THEN
                 the_vza_atm = the_vza_atm + 1.955
              ELSE IF (gome_curscan == 0 .OR. gome_curscan == 2) THEN
                 the_vza_atm = the_vza_atm - 1.23
              ENDIF
              
              IF ( n_coadd == 4 .AND. gome_curscan == 3) THEN 
                 sumloc = gome_geoloc
              ENDIF
              
              sumtsza = sumtsza + 1.0 / COS(the_sza_atm * deg2rad)
              sumsza = sumsza   + 1.0 / COS(sza_atm * deg2rad)
              sumtvza = sumtvza + 1.0 / COS(the_vza_atm * deg2rad)
              sumvza = sumvza   + 1.0 / COS(vza_atm * deg2rad)     
              sumtaza = sumtaza + the_aza_atm
              sumaza = sumaza   + aza_atm  
              sumcurv = sumcurv +  earth_curv
              sumalt = sumalt   +  ers2_alt

              ! read level 1 data
              CALL read_all_gomerad(funit, nband, gome_curqual, radspec(n_coadd, :, :, :), &
                   nwpos(n_coadd, :), stpos(n_coadd, :), file_read_stat)
              IF (file_read_stat /= file_read_ok ) RETURN              
           ENDIF
        ENDIF

        ! exit
        IF (gome_curqual == 4) EXIT
     ENDDO

     ! check if valid data is found
     IF (gome_curqual == 4 .AND. gome_curscan == 3 .AND. valid_pix .AND. &
          ((n_coadd == 8  .AND. gome_curpix - gome_stpix == 7) .OR. &
          (n_coadd == 40 .AND. gome_curpix - gome_stpix == 39))) THEN

        sumloc(:, 1) = gome_geoloc(:, 1)
        sumloc(:, 3) = gome_geoloc(:, 3) 
        sumloc(:, 5) = (gome_geoloc(:, 5) + sumloc(:, 5)) / 2.0D0

        ! Check if read data are ok, should not happen
        IF ( ABS(sumloc(1, 1) - sumloc(1, 2)) > 5.0      .OR. &   ! avoid Himalayas data gap
             ANY(nwpos(2:n_coadd, 2)   /= nwpos(1, 2))   .OR. &
             ANY(stpos(2:n_coadd, 2)   /= stpos(1, 2))   .OR. &
             ANY(nwpos(2:n_coadd-1, 1) /= nwpos(1, 1))   .OR. &
             ANY(stpos(2:n_coadd-1, 1) /= stpos(1, 1))   .OR. &
             nwpos(n_coadd, 1) <= nwpos(1, 1)            .OR. & 
             stpos(n_coadd, 1) >= stpos(1, 1)) THEN
           file_read_stat = file_read_missing; RETURN
        ENDIF
        b1ab_div_wav = radspec(n_coadd, 1, stpos(1, 1)-1, 1)

         !Average and subtract data
        CALL avg_and_subtract (radspec(1:n_coadd, :, :, :), n_coadd, nband, stpos(1, 1), &
             nwpos(n_coadd, :), numwin, winpix(1:numwin, :), winlim(1:numwin, :),        &
             band_selectors(1:numwin), n_gome_radpts, gome_radspec, nradpix(1:numwin), &
             file_read_stat)
        IF ( file_read_stat /= file_read_ok ) RETURN

        gome_endpix = gome_curpix; gome_npix = n_coadd
        nview = n_gome_ang; nloc = n_gome_geo
        
        the_sza_atm = ACOS(1.0D0 / (sumtsza / n_coadd)) * rad2deg
        the_vza_atm = ACOS(1.0D0 / (sumtvza / n_coadd)) * rad2deg
        the_aza_atm = sumtaza / n_coadd
        sza_atm = ACOS(1.0D0 / (sumsza / n_coadd)) * rad2deg
        vza_atm = ACOS(1.0D0 / (sumvza / n_coadd)) * rad2deg
        aza_atm = sumaza / n_coadd
        aza_atm(2) = (aza_atm(1) + aza_atm(3)) / 2.0
        earth_curv = sumcurv / n_coadd
        ers2_alt =   sumalt  / n_coadd
        gome_geoloc = sumloc
        
        the_sca_atm = 180.0 - ACOS(COS(the_sza_atm * deg2rad) * COS(the_vza_atm &
             * deg2rad) + SIN(the_sza_atm * deg2rad)* SIN(the_vza_atm * deg2rad)&
             * COS(the_aza_atm * deg2rad)) * rad2deg
        
        WRITE(*, *) 'nocadd = ', n_coadd
        !WRITE(*, '(4x, A5, 3f8.2)')  'SZA =', sza_atm
        !WRITE(*, '(4x, A5, 3f8.2)')  'VZA =', vza_atm
        !WRITE(*, '(4x, A5, 3f8.2)')  'AZA =', aza_atm
        WRITE(*, '(4x, A, 4f8.2)') 'Effective SZA, VZA, AZA = ', &
             the_sza_atm, the_vza_atm, the_aza_atm, the_sca_atm
        WRITE(*, *) ' ' 
     ELSE
        file_read_stat = file_read_missing; RETURN
     ENDIF
  ENDIF

  nview = n_gome_ang
  nloc = n_gome_geo
        
  the_lats(1:nloc)  = gome_geoloc(lat_idx, 1:n_gome_geo)
  the_lons(1:nloc)  = gome_geoloc(lon_idx, 1:n_gome_geo)
  
  WHERE (the_lons(1:nloc) > 180.0) 
     the_lons(1:nloc) = the_lons(1:nloc) - 360.0
  ENDWHERE
  the_lon =  the_lons(nloc); the_lat = the_lats(nloc)   ! center cordinate
  
  ! Get left/right center of a pixel  
  IF (ABS(the_lons(1) - the_lons(2)) < 100.0) THEN   ! normal
     edgelons(1) = (the_lons(1) + the_lons(2)) / 2.0
  ELSE                                               ! across date line
     lon1 = the_lons(1); lon2 = the_lons(2)
     IF (lon1 > 0.0) THEN 
        lon2 = lon2 + 360.0
     ELSE
        lon1 = lon1 + 360.0
     ENDIF
     edgelons(1) = (lon1 + lon2) / 2.0
     IF (edgelons(1) > 180.0) edgelons(1) = edgelons(1) - 360.0 
  ENDIF
  
  IF (ABS(the_lons(3) - the_lons(4)) < 100.0) THEN   ! normal
     edgelons(2) = (the_lons(3) + the_lons(4)) / 2.0
  ELSE                                               ! across date line
     lon1 = the_lons(3); lon2 = the_lons(4)
     IF (lon1 > 0.0) THEN 
        lon2 = lon2 + 360.0
     ELSE
        lon1 = lon1 + 360.0
     ENDIF
     edgelons(2) = (lon1 + lon2) / 2.0
     IF (edgelons(2) > 180.0) edgelons(2) = edgelons(2) - 360.0 
  ENDIF
  
  edgelats(1) = (the_lats(1) + the_lats(2)) / 2.0
  edgelats(2) =  (the_lats(3) + the_lats(4)) / 2.0
  
  ! special handling for data inteprolation
  ! If pixel across date line, then edgelons [0-360], otherwise [-180W, 180E]
  IF (ABS(edgelons(1) - edgelons(2)) > 100.0) THEN                                            
     IF (edgelons(1) > 0.0) THEN 
        edgelons(2) = edgelons(2) + 360.0
     ELSE
        edgelons(1)= edgelons(1) + 360.0
     ENDIF
  ENDIF
  
  ! make sure that edgelon(1) < edgelon(2)
  IF (edgelons(1) > edgelons(2)) THEN
     templon = edgelons(2)
     edgelons(2) = edgelons(1)
     edgelons(1) = templon
     
     templat = edgelats(2)
     edgelats(2) = edgelats(1)
     edgelats(1) = templat
  ENDIF
        
  year = gome_pixdate(8:11);  READ(year,'(I4)') the_year
  day =  gome_pixdate (1:2);  READ(day, '(I2)') the_day
  monthc = gome_pixdate(4:6)
  
  DO i = 1, 12
     IF (monthc == months(i)) EXIT
  ENDDO
  the_month = i
  
  spec_norm = SUM(gome_radspec(2,1:n_gome_radpts)) / n_gome_radpts
  IF (spec_norm <= 0.0 ) THEN
     error = .TRUE.;  RETURN
  ENDIF

  gome_radspec(2,1:n_gome_radpts) = gome_radspec(2,1:n_gome_radpts) / spec_norm
  gome_radspec(3,1:n_gome_radpts) = gome_radspec(3,1:n_gome_radpts) / spec_norm
  div_rad = spec_norm

  RETURN

END SUBROUTINE gome_read_el1data_radavg


SUBROUTINE read_all_gomerad(funit, nb, nsub, radspec, numpix, stpos, file_read_stat)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : maxchlen
  USE OMSAO_errstat_module,    ONLY : file_read_ok, file_read_failed, file_read_eof
  USE OMSAO_gome_data_module,  ONLY : n_gome_data_dim, n_gome_max_pts

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,        INTENT (IN)          :: funit, nb, nsub

  ! ================
  ! Output variables
  ! ================
  INTEGER, INTENT (OUT)                :: file_read_stat
  INTEGER, DIMENSION (nb), INTENT(OUT) :: numpix, stpos
  REAL (KIND=dp), DIMENSION (nb, n_gome_max_pts, n_gome_data_dim), INTENT(OUT) :: radspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER, DIMENSION (4)      :: numtpix = (/695, 841, 1024, 1024/)
  INTEGER                     :: ios, i, j, totpix, ch, offset
  REAL (KIND=dp)              :: integt, swav, ewav
  CHARACTER (LEN=9)           :: tmpchar
  CHARACTER (LEN=5)           :: tmpchar1
  CHARACTER (LEN=1)           :: subch
  CHARACTER (LEN=maxchlen)    :: lastline
 
  ! ----------------------------------
  ! Initialization of output variables
  ! ----------------------------------
  file_read_stat = file_read_ok
   
  ! Read all the data 
  radspec = 0.0; numpix = 0; stpos = 0
  DO i = 1, nsub 
     CALL skip_to_filemark ( funit, 'Band', lastline, file_read_stat )
     IF ( file_read_stat /= file_read_ok ) THEN
        STOP 'Read Radiance: Cannot find header!!!'
     END IF
     
     READ(lastline, '(A9,F9.5,2F9.4,I5)') tmpchar, integt, swav, ewav, totpix
     READ(tmpchar, '(A5,I1,A1)') tmpchar1, ch, subch

     IF ( (ch == 1 .OR. ch == 2) .AND. subch == 'b' ) THEN
        offset = numtpix(ch) - totpix
        IF (stpos(ch) == 0) stpos(ch) = offset + 1
     ELSE
        offset = 0
        IF (stpos(ch) == 0) stpos(ch) = 1
     ENDIF
     numpix(ch) = numpix(ch) + totpix
        
     DO j = 1, totpix
        READ (UNIT=funit, FMT=*, IOSTAT = ios) &
             radspec(ch, offset+j, 1:n_gome_data_dim)
        IF (ios /= 0 ) THEN
           WRITE(*, *) 'Read Radiance: Read Error: ', tmpchar
           STOP
        ENDIF
     ENDDO
     !PRINT *, radspec(ch, 1, 1), radspec(ch, 626, 1), radspec(ch, 695, 1)
  ENDDO

  CALL gome_check_read_status ( ios, file_read_stat)
  
  RETURN

END SUBROUTINE read_all_gomerad


SUBROUTINE avg_and_subtract (radspec, n_coadd, nb, b1b_stpix, &
     numpix, nwin, winpix, winlim, bds, ngome, gspec, nradpix, file_read_stat)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : maxchlen, maxwin, max_fit_pts, mrefl
  USE OMSAO_gome_data_module,  ONLY : n_gome_data_dim, n_gome_max_pts
  USE OMSAO_errstat_module,    ONLY : file_read_ok, file_read_failed, &
       file_read_eof, file_read_missing
  USE OMSAO_variables_module,  ONLY : nsolpix
  USE ozprof_data_module,      ONLY : rad_posr, rad_specr, nrefl

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER, INTENT (IN)                     :: nb, n_coadd, nwin, b1b_stpix
  INTEGER, DIMENSION (nwin, 2), INTENT(IN) :: winpix
  INTEGER, DIMENSION (nwin), INTENT(IN)    :: bds
  INTEGER, DIMENSION (nb), INTENT(IN)      :: numpix
  REAL (KIND=dp), DIMENSION (n_coadd, nb, n_gome_max_pts, n_gome_data_dim), &
       INTENT(IN)                          :: radspec
  
  ! ================
  ! Output variables
  ! ================
  INTEGER, INTENT (OUT)                    :: ngome, file_read_stat
  INTEGER, DIMENSION(nwin), INTENT(OUT)    :: nradpix
  REAL (KIND=dp), DIMENSION(nwin, 2), INTENT(INOUT) :: winlim
  REAL (KIND=dp), DIMENSION(n_gome_data_dim, max_fit_pts), INTENT(OUT):: gspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER         :: ch, i, j
  REAL (KIND=dp) :: prewav

  file_read_stat = file_read_ok

  ! Get data for fitting
  ngome = 0 ; gspec = 0.0
  prewav = 0.0
  DO i = 1, nwin
     ch = bds(i)
     IF (ch > nb) THEN
        file_read_stat = file_read_failed; RETURN
     ENDIF
     IF (numpix(ch) == 0) THEN
        WRITE(*, *) 'No pixels are read from channel ', ch
        file_read_stat = file_read_missing; RETURN
     ENDIF
     nradpix(i) = 0
     IF (ch /= 1) THEN
        DO j = winpix(i, 1), winpix(i, 2)
           IF (ALL(radspec(1:n_coadd, ch, j, 2) /= 0.0) .AND. &
                (ngome == 0 .OR. radspec(1, ch, j, 1) > prewav + 0.1) &
                .AND. ALL(radspec(1:n_coadd, ch, j, 5) == 0)) THEN
              ngome = ngome + 1; nradpix(i) = nradpix(i) + 1
              prewav = gspec(1, ngome)
              gspec(1, ngome) = SUM(radspec(1:n_coadd, ch, j, 1)) / n_coadd
              gspec(2, ngome) = SUM(radspec(1:n_coadd, ch, j, 2)) / n_coadd
              gspec(3, ngome) = SUM(radspec(1:n_coadd, ch, j, 3)) / n_coadd / SQRT(1.0 * n_coadd)
              gspec(4, ngome) = SUM(radspec(1:n_coadd, ch, j, 4)) / n_coadd
              gspec(5, ngome) = SUM(radspec(1:n_coadd, ch, j, 5)) / n_coadd
           ENDIF
        ENDDO
     ELSE   ! channel 1
        DO j = winpix(i, 1), winpix(i, 2) 
           IF (j < b1b_stpix .AND. radspec(n_coadd, ch, j, 2) /= 0.0 &
                .AND. j >= 231 &  ! > 264.0 nm
                .AND. (j >= 410 .OR. j <= 300)          &  ! no 276.0-282.0 nm 
                .AND. (j >= 448 .OR. j <= 410)          &  ! no 284-289 nm
                .AND. (ngome == 0 .OR. radspec(n_coadd, ch, j, 1) > prewav + 0.1) &
                .AND. radspec(n_coadd, ch, j, 5) == 0) THEN
              ngome = ngome + 1; nradpix(i) = nradpix(i) + 1
              gspec(:, ngome) = radspec(n_coadd, ch, j, :)
              prewav = gspec(1, ngome)
           ELSE IF (j>= b1b_stpix .AND. ALL(radspec(1:n_coadd, ch, j, 2) /= 0.0)  &
                .AND. (j >= 410 .OR. j <= 300)          &  ! no 276.0-282.0 nm 
                .AND. (j >= 448 .OR. j <= 410)          &  ! no 284-289 nm
                .AND. (ngome == 0 .OR. radspec(1, ch, j, 1) > prewav + 0.1) &
                .AND. ALL(radspec(1:n_coadd, ch, j, 5) == 0)) THEN
              ngome = ngome + 1; nradpix(i) = nradpix(i) + 1
              prewav = gspec(1, ngome)
              gspec(1, ngome) = SUM(radspec(1:n_coadd, ch, j, 1)) / n_coadd
              gspec(2, ngome) = SUM(radspec(1:n_coadd, ch, j, 2)) / n_coadd
              gspec(3, ngome) = SUM(radspec(1:n_coadd, ch, j, 3)) / n_coadd / SQRT(1.0 * n_coadd)
              gspec(4, ngome) = SUM(radspec(1:n_coadd, ch, j, 4)) / n_coadd
              gspec(5, ngome) = SUM(radspec(1:n_coadd, ch, j, 5)) / n_coadd
           ENDIF
        ENDDO
     ENDIF

     IF (i == 1) THEN
        winlim(i, 1) = gspec(1,1) 
     ELSE
        winlim(i, 1) = gspec(1, SUM(nradpix(1:i-1))+1)
     ENDIF
     winlim(i, 2) = gspec(1, ngome)
  ENDDO

  IF (ngome > max_fit_pts) THEN
     WRITE(*, *) 'Read Radiance: Increase max_fit_pts!!!'
     file_read_stat = file_read_failed; RETURN
  ELSE IF (ngome == 0) THEN
     WRITE(*, *) 'Read Radiance: NO radiances are read!!!'
     file_read_stat = file_read_failed; RETURN
  ENDIF
  
  ! Get data for surface albedo at 370.2 nm +/- 20 pixels
  j = 0
  DO i = 500, numpix(2)
     IF (ALL(radspec(1:n_coadd, 2, i, 2) /= 0.0) .AND. ALL(radspec(1:n_coadd, 2, i, 5) == 0)) THEN
        j = j + 1
        rad_posr(j)  = SUM(radspec(1:n_coadd, 2, i, 1)) / n_coadd
        rad_specr(j) = SUM(radspec(1:n_coadd, 2, i, 2)) / n_coadd
        IF (j == nrefl) EXIT
     ENDIF
  ENDDO

  WRITE(*, *) 'End of Reading Radiance Spectrum'
  DO i = 1, nwin
     WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', i, winlim(i,1), winlim(i,2), nradpix(i)
     IF (nradpix(i) <= 20) THEN
        WRITE(*, '(A,f8.3,A3,f8.3)') 'Not enough points in window: ', winlim(i,1),&
             ' - ', winlim(i,2)
        file_read_stat = file_read_missing; RETURN
     ENDIF
     
     IF (nradpix(i) > nsolpix(i)) THEN
        WRITE(*, *) 'Inconsistent pixels between I and F0 in channel ', bds(i)
        file_read_stat = file_read_missing; RETURN
     ENDIF
  ENDDO
  
  RETURN

END SUBROUTINE avg_and_subtract

! Purpose: read all the data corresponding to a gome channel 1
! Use channel 1 (large spectral resolution) togehter with each individual channel 2
! 1. Check whether there are data that are already read
! 2. If all data are already read, then constitude the spectrum one by one: 
!    channel 1 (same for all) + channel 2 individual
! 3. If there is no more data for the channel 1, read all the data 
!    corresponding to a new channel 1 spectra

SUBROUTINE gome_read_el1data_allradch1 ( funit, file_read_stat, error )

  USE OMSAO_precision_module
  USE OMSAO_errstat_module,    ONLY : &
       file_read_ok, file_read_failed, file_read_missing, file_read_eof
  USE OMSAO_parameters_module, ONLY : maxchlen, deg2rad, rad2deg, maxwin
  USE OMSAO_gome_data_module,  ONLY : &
       azm_idx, zen_idx, zen0_idx, lat_idx, lon_idx, n_gome_ang, n_gome_geo,   &
       lm_gome_groundpix, gome_spec_missing, n_gome_max_pts, n_gome_data_dim,  &
       n_gome_radpts, gome_radspec, gome_curpix, gome_curqual, gome_curscan,   &
       gome_pixdate, gome_angles_wrtn, gome_angles_wrts, gome_geoloc, ers2_alt,&
       earth_curv, orbnum, gome_stpix, gome_endpix, gome_npix, allsza, allvza, &
       allaza, allsca, thecurpix, allgeoloc, allpixdate, allspec, allcurscan,  &
       allcurpix, nwpos, stpos, maxc1c2r
  USE ozprof_data_module,      ONLY : div_rad 
  USE OMSAO_variables_module,  ONLY : the_month, the_year, the_day, numwin,    &
       winlim, winpix, band_selectors, nradpix, npix_fitting, sza_atm, vza_atm,&
       aza_atm, the_sza_atm, the_aza_atm, the_vza_atm, the_sca_atm, the_lat,   &
       the_lon, the_lons, the_lats, nview, nloc, edgelons, edgelats,           &
       b1ab_div_wav, avgsza, avgaza, avgvza, avgsca , szamax 

  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER,  INTENT (IN)    :: funit
  LOGICAL, INTENT (OUT)    :: error	
  INTEGER, INTENT (OUT)    :: file_read_stat

  ! ===============
  ! Local variables
  ! ===============
  INTEGER, PARAMETER       :: maxnpix = 40, nband = 2
  INTEGER                  :: i, j, nch2
  LOGICAL                  :: valid_pix
  CHARACTER (LEN=maxchlen) :: lastline
  CHARACTER (LEN=3)        :: monthc
  CHARACTER (LEN=2)        :: day
  CHARACTER (LEN=4)        :: year
  REAL (KIND=dp)           :: spec_norm, lon1, lon2, templon, templat
  REAL (KIND=dp)           :: sumsza, sumaza, sumvza
  CHARACTER (LEN=3), DIMENSION(12), PARAMETER :: months = (/'JAN', 'FEB', 'MAR', &
       'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/)
 
  file_read_stat = file_read_ok ; error = .FALSE.
  IF (npix_fitting == 0 .AND. orbnum == 0) THEN
     thecurpix = 0
     CALL skip_to_filemark ( funit, 'ERS Information', lastline, file_read_stat )
     IF ( file_read_stat /= file_read_ok ) RETURN
     READ (funit, *), i, i, orbnum
  ENDIF

  nch2    =   0;    sumsza=   0.0;   sumvza  = 0.0;   sumaza = 0.0
  DO WHILE (thecurpix == 0)
     CALL skip_to_filemark (funit, lm_gome_groundpix, lastline, file_read_stat)
     IF ( file_read_stat /= file_read_ok ) RETURN
     READ (lastline(14:17), '(I4)') gome_curpix
     READ (lastline(19:19), '(I1)') gome_curqual
     READ (lastline(21:21), '(I1)') gome_curscan

     CALL gome_read_el1geo (funit, gome_pixdate, gome_angles_wrtn, &
          gome_angles_wrts, gome_geoloc, ers2_alt, earth_curv, file_read_stat)
     IF ( ANY(gome_angles_wrtn(zen0_idx, 1:n_gome_ang) >= szamax) ) THEN
        file_read_stat = file_read_missing; RETURN
     END IF

     ! First pixel:
     ! 1. scan = 0  2. sza le 90.0  3. # bands=3 (b1b,b2a,b2b) 4: nch2=0

     valid_pix = .FALSE.
     IF (gome_curqual >= 3 .AND. nch2 < 8 .AND. &
          ALL(gome_angles_wrtn(zen0_idx, 1:n_gome_ang) < szamax) ) THEN
     
        IF (nch2 == 0 .AND. gome_curscan == 0) THEN
           nch2 = nch2 + 1; gome_stpix = gome_curpix; valid_pix = .TRUE.
        ELSE IF (nch2 > 0) THEN
           nch2 = nch2 + 1; valid_pix = .TRUE.
        ENDIF
        ! ELSE: have not found the first pixel yet
        
        IF (valid_pix) THEN
           CALL calc_view_geometry
           
           IF ( the_sza_atm < 0.0 .OR. the_sza_atm >= szamax) THEN
              file_read_stat = file_read_missing; RETURN
           END IF
           IF (gome_curscan == 3) THEN
              the_vza_atm = the_vza_atm + 1.955
           ELSE IF (gome_curscan == 0 .OR. gome_curscan == 2) THEN
              the_vza_atm = the_vza_atm - 1.23
           ENDIF

           allgeoloc(nch2, :, :) = gome_geoloc
           allpixdate(nch2)   = gome_pixdate           
           allsza(nch2)       = the_sza_atm
           allvza(nch2)       = the_vza_atm
           allaza(nch2)       = the_aza_atm
           allsca(nch2)       = the_sca_atm
           allcurpix(nch2)    = gome_curpix
           allcurscan(nch2)   = gome_curscan
           
           sumsza = sumsza + 1.0 / COS(the_sza_atm * deg2rad)
           sumvza = sumvza + 1.0 / COS(the_vza_atm * deg2rad)    
           sumaza = sumaza + the_aza_atm

           ! read level 1 data
           CALL read_all_gomerad(funit, nband, gome_curqual, allspec(nch2, 1:nband, :, :), &
                nwpos(nch2, 1:nband), stpos(nch2, 1:nband), file_read_stat)

           IF (file_read_stat /= file_read_ok ) RETURN              
        ENDIF
     ENDIF

     ! finish reading all channels corresponding to one GOME 1 channel
     IF (gome_curqual == 4) THEN
        ! check if valid data is found
        IF (gome_curqual == 4 .AND. gome_curscan == 3 .AND. valid_pix .AND. &
           ((nch2 == maxc1c2r  .AND. gome_curpix - gome_stpix == maxc1c2r-1))) THEN
           
           ! Check if read data are ok, should not happen
           IF ( ANY(nwpos(2:nch2, 2)   /= nwpos(1, 2))   .OR. &
                ANY(stpos(2:nch2, 2)   /= stpos(1, 2))   .OR. &
                ANY(nwpos(2:nch2-1, 1) /= nwpos(1, 1))   .OR. &
                ANY(stpos(2:nch2-1, 1) /= stpos(1, 1))   .OR. &
                nwpos(nch2, 1) <= nwpos(1, 1)            .OR. & 
                stpos(nch2, 1) >= stpos(1, 1)) THEN
              file_read_stat = file_read_missing; RETURN
           ENDIF
           b1ab_div_wav = allspec(nch2, 1, stpos(1, 1)-1, 1)
           gome_endpix = gome_curpix; gome_npix = nch2
           nview = n_gome_ang; nloc = n_gome_geo
           
           avgsza = ACOS(1.0D0 / (sumsza / nch2)) * rad2deg
           avgvza = ACOS(1.0D0 / (sumvza / nch2)) * rad2deg
           avgaza = sumaza / nch2
           
           avgsca = 180.0 - ACOS(COS(the_sza_atm * deg2rad) * COS(the_vza_atm &
                * deg2rad) + SIN(the_sza_atm * deg2rad)* SIN(the_vza_atm * deg2rad)&
                * COS(the_aza_atm * deg2rad)) * rad2deg

           thecurpix = 1

           DO i = 1, maxc1c2r - 1 
              allspec(i, 1, 1:stpos(1, 1)-1, :) = allspec(maxc1c2r, 1, 1:stpos(1, 1)-1, :)
           ENDDO

           !DO i = 1, maxc1c2r
           !   WRITE(*, '(2I5, 4f10.4)') allcurpix(i), allcurscan(i), allsza(i), allvza(i), allaza(i), allsca(i)
           !ENDDO
           !WRITE(*, '(4f10.4)') avgsza, avgvza, avgaza, avgsca
           !WRITE(*, '(6I8)') ((allcurpix(i), allcurscan(i), nwpos(i, 1:2), stpos(i, 1:2)), i = 1, nch2)

           EXIT  ! finish reading all the data corresponding to a gome pixel
        ELSE ! not valid data
           file_read_stat = file_read_missing; RETURN
        ENDIF       
     ENDIF
  ENDDO

  gome_geoloc  = allgeoloc(thecurpix,  :, :)
  gome_pixdate = allpixdate(thecurpix)
  gome_curpix  = allcurpix(thecurpix)
  gome_curscan = allcurscan(thecurpix)
  the_sza_atm  = allsza(thecurpix)
  the_vza_atm  = allvza(thecurpix)
  the_aza_atm  = allaza(thecurpix)
  the_sca_atm  = allsca(thecurpix)
  gome_stpix   = gome_curpix
  gome_endpix  = gome_curpix
  gome_npix    = 1
          
  the_lats(1:nloc)  = gome_geoloc(lat_idx, 1:n_gome_geo)
  the_lons(1:nloc)  = gome_geoloc(lon_idx, 1:n_gome_geo)
     
  WHERE (the_lons(1:nloc) > 180.0) 
     the_lons(1:nloc) = the_lons(1:nloc) - 360.0
  ENDWHERE
  the_lon =  the_lons(nloc); the_lat = the_lats(nloc)   ! center cordinate
     
  ! Get left/right center of a pixel  
  IF (ABS(the_lons(1) - the_lons(2)) < 100.0) THEN   ! normal
     edgelons(1) = (the_lons(1) + the_lons(2)) / 2.0
  ELSE                                               ! across date line
     lon1 = the_lons(1); lon2 = the_lons(2)
     IF (lon1 > 0.0) THEN 
        lon2 = lon2 + 360.0
     ELSE
        lon1 = lon1 + 360.0
     ENDIF
     edgelons(1) = (lon1 + lon2) / 2.0
     IF (edgelons(1) > 180.0) edgelons(1) = edgelons(1) - 360.0 
  ENDIF
  
  IF (ABS(the_lons(3) - the_lons(4)) < 100.0) THEN   ! normal
     edgelons(2) = (the_lons(3) + the_lons(4)) / 2.0
  ELSE                                               ! across date line
     lon1 = the_lons(3); lon2 = the_lons(4)
     IF (lon1 > 0.0) THEN 
        lon2 = lon2 + 360.0
     ELSE
        lon1 = lon1 + 360.0
     ENDIF
     edgelons(2) = (lon1 + lon2) / 2.0
     IF (edgelons(2) > 180.0) edgelons(2) = edgelons(2) - 360.0 
  ENDIF
  
  edgelats(1) = (the_lats(1) + the_lats(2)) / 2.0
  edgelats(2) =  (the_lats(3) + the_lats(4)) / 2.0
         
  ! special handling for data inteprolation
  ! If pixel across date line, then edgelons [0-360], otherwise [-180W, 180E]
  IF (ABS(edgelons(1) - edgelons(2)) > 100.0) THEN                                            
     IF (edgelons(1) > 0.0) THEN 
        edgelons(2) = edgelons(2) + 360.0
     ELSE
        edgelons(1)= edgelons(1) + 360.0
     ENDIF
  ENDIF
  
  ! make sure that edgelon(1) < edgelon(2)
  IF (edgelons(1) > edgelons(2)) THEN
     templon = edgelons(2)
     edgelons(2) = edgelons(1)
     edgelons(1) = templon
     
     templat = edgelats(2)
     edgelats(2) = edgelats(1)
     edgelats(1) = templat
  ENDIF
  
  year = gome_pixdate(8:11);  READ(year,'(I4)') the_year
  day =  gome_pixdate (1:2);  READ(day, '(I2)') the_day
  monthc = gome_pixdate(4:6)
     
  DO i = 1, 12
     IF (monthc == months(i)) EXIT
  ENDDO       
  the_month = i

  CALL avg_and_subtract (allspec(thecurpix, 1:nband, :, :), 1, nband, stpos(1, 1), &
       nwpos(maxc1c2r, 1:nband), numwin, winpix(1:numwin, :), winlim(1:numwin, :),        &
       band_selectors(1:numwin), n_gome_radpts, gome_radspec, nradpix(1:numwin), &
       file_read_stat)

  IF ( file_read_stat /= file_read_ok ) RETURN

  spec_norm = SUM(gome_radspec(2,1:n_gome_radpts)) / n_gome_radpts
  IF (spec_norm <= 0.0 ) THEN
     error = .TRUE.;  RETURN
  ENDIF

  gome_radspec(2, 1:n_gome_radpts) = gome_radspec(2, 1:n_gome_radpts) / spec_norm
  gome_radspec(3, 1:n_gome_radpts) = gome_radspec(3, 1:n_gome_radpts) / spec_norm
  div_rad = spec_norm
   
  thecurpix = thecurpix + 1

  ! neglect backscan pixels since they are overlapped with forward scans
  IF (allcurscan(thecurpix) == 3 .AND. thecurpix /= maxc1c2r) THEN
     thecurpix = thecurpix + 1
  ELSE IF (thecurpix == maxc1c2r) THEN
     thecurpix = 0
  ENDIF

  ! include backscan pixels as well
  !IF (thecurpix > maxc1c2r) thecurpix = 0
  
  RETURN

END SUBROUTINE gome_read_el1data_allradch1


