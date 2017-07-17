SUBROUTINE scia_read_el1data_sol (funit, file_read_stat, error )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,ONLY : maxchlen
  USE OMSAO_gome_data_module, ONLY : lm_gome_solspec, lm_gome_eshine,&
       n_gome_max_pts, n_gome_data_dim, gome_spec_missing, n_gome_solpts, gome_solspec
  USE OMSAO_errstat_module,   ONLY : file_read_ok, file_read_failed, file_read_missing, &
       file_read_eof
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
  CALL skip_to_filemark ( funit, 'Solar Irradiance', lastline, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) RETURN

  ! Read solar spectrum
  CALL scia_read_el1solspec (funit, numwin, winpix(1:numwin, :), &
       winlim(1:numwin, :), band_selectors(1:numwin), n_gome_solpts, &
       gome_solspec(1:3, :), nsolpix(1:numwin), file_read_stat)  
  IF ( file_read_stat /= file_read_ok ) RETURN

  ! -----------------------------------
  ! Normalize solar irradiance spectrum
  ! -----------------------------------
  ! Rough Correction for SCIAMACHY measurements (probably in irradiances)
  ! for data before Dec 24, 2003, a correction of 0.8566 is made
  ! for data after Dec 24., 2003, a correction of 1.1078 is made

  gome_solspec(2,1:n_gome_solpts) = gome_solspec(2,1:n_gome_solpts)  * 0.8566
  sol_spec_ring(2, 1:nsol_ring) = sol_spec_ring(2, 1:nsol_ring) * 0.8566

  solar_norm = SUM (gome_solspec(2,1:n_gome_solpts)) / n_gome_solpts
  IF ( solar_norm <= 0.0 ) THEN 
     error = .TRUE.; RETURN
  ENDIF
 
  gome_solspec(2,1:n_gome_solpts) = gome_solspec(2,1:n_gome_solpts) / solar_norm
  gome_solspec(3,1:n_gome_solpts) = gome_solspec(3,1:n_gome_solpts) * gome_solspec(2,1:n_gome_solpts)
  sol_spec_ring(2, 1:nsol_ring)   = sol_spec_ring(2, 1:nsol_ring)   / solar_norm 

  div_sun = solar_norm

  ! ---------------------------------------------------------------------
  ! Position cursor after solar spectrum, before first ground pixel entry
  ! ---------------------------------------------------------------------
  CALL skip_to_filemark ( funit, lm_gome_eshine, lastline, file_read_stat )
  
  RETURN
END SUBROUTINE scia_read_el1data_sol


SUBROUTINE scia_read_el1data_rad ( funit, file_read_stat, error )

  USE OMSAO_precision_module
  USE OMSAO_errstat_module,    ONLY : &
       file_read_ok, file_read_failed, file_read_missing, file_read_eof
  USE OMSAO_parameters_module, ONLY : maxchlen, deg2rad, rad2deg
  USE OMSAO_gome_data_module,  ONLY : &
       azm_idx, zen_idx, zen0_idx, lat_idx, lon_idx, n_gome_ang, n_gome_geo,   &
       lm_gome_groundpix, gome_spec_missing, n_gome_max_pts, n_gome_data_dim,  &
       n_gome_radpts, gome_radspec, gome_curpix, gome_curqual, gome_curscan,   &
       gome_pixdate, gome_angles_wrtn, gome_angles_wrts, gome_geoloc, ers2_alt,&
       earth_curv, orbnum, gome_stpix, gome_endpix, gome_npix, gome_integt
  USE ozprof_data_module,      ONLY : div_rad, scia_coadd
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
  INTEGER                  :: i, j, n_coadd, npts, stpix, endpix, fidx, lidx
  INTEGER, PARAMETER       :: nwp = 1048, raderr_unit =33
  CHARACTER (LEN=maxchlen) :: lastline
  CHARACTER (LEN=3)        :: monthc
  CHARACTER (LEN=2)        :: day
  CHARACTER (LEN=4)        :: year
  REAL (KIND=dp)           :: spec_norm, lon1, lon2, templon, templat, temp, integt
  REAL (KIND=dp), DIMENSION(n_gome_ang):: sumsza, sumaza, sumvza
  REAL (KIND=dp)  :: sumtsza, sumtaza, sumtvza, sumalt, sumcurv
  REAL (KIND=dp), DIMENSION(lon_idx, n_gome_geo) :: sumloc
  CHARACTER (LEN=3), DIMENSION(12), PARAMETER :: months = (/'JAN', 'FEB', 'MAR', &
       'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/)

  REAL (KIND=dp), DIMENSION(nwp), SAVE :: slope, intercept
  LOGICAL                        , SAVE :: first = .TRUE.

  file_read_stat = file_read_ok ; error = .FALSE.
  
  IF (npix_fitting == 0 .AND. orbnum == 0) THEN
     orbnum = 9999
  ENDIF

  ! ----------------------------------------------------------------
  ! Position cursor at start of GOME ground pixel entry
  ! And Extract GOME ground pixel number and scan type from LASTLINE
  ! -----------------------------------------------------------------
  CALL skip_to_filemark ( funit, 'Ground Pixel', lastline, file_read_stat)
  IF ( file_read_stat /= file_read_ok ) RETURN
  lastline = lastline(13:30)
  READ(lastline, *) gome_curpix, gome_curqual, gome_curscan  
  IF ( gome_curqual == gome_spec_missing) THEN
     file_read_stat = file_read_missing; RETURN
  END IF
  READ (funit, '(A39)') lastline
  gome_pixdate = lastline(1:22)
  READ(lastline, *) integt, integt, integt
  gome_integt = integt

  READ (funit, *) gome_angles_wrtn(1:2, :)
  READ (funit, *) gome_angles_wrtn(3:4, :)
  READ (funit, *) gome_angles_wrts(1:2, :)
  READ (funit, *) gome_angles_wrts(3:4, :)
  
  ! decide if there is coadding
  IF (gome_angles_wrts(1, 1) /= gome_angles_wrtn(1, 1)) THEN
     scia_coadd = .TRUE.
     !file_read_stat = file_read_missing; RETURN
  ELSE
     scia_coadd = .FALSE.
  ENDIF
  gome_angles_wrts = gome_angles_wrtn
  !PRINT *, gome_angles_wrtn(1, :)
  
  READ (funit, *) ers2_alt, earth_curv
  READ (funit, *) gome_geoloc

  WHERE(gome_geoloc(lon_idx, :) < 0.0)
     gome_geoloc(lon_idx, :) = gome_geoloc(lon_idx, :) + 360.0
  ENDWHERE
  
  IF ( ANY(gome_angles_wrtn(zen0_idx, 1:n_gome_ang) >= szamax) ) THEN
     file_read_stat = file_read_missing; RETURN
  END IF   

  CALL calc_view_geometry
  IF ( the_sza_atm < 0.0 .OR. the_sza_atm >= szamax) THEN
     file_read_stat = file_read_missing; RETURN
  END IF
   
  gome_stpix = gome_curpix; gome_endpix = gome_curpix; gome_npix = 1
  
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
  
  !   year = gome_pixdate(8:11);  READ(year,'(I4)') the_year
  !   day =  gome_pixdate (1:2);  READ(day, '(I2)') the_day
  !   monthc = gome_pixdate(4:6)
  !   
  !   DO i = 1, 12
  !      IF (monthc == months(i)) EXIT
  !   ENDDO
  !   the_month = i

  year   = gome_pixdate(1:4);   READ(year,'(I4)') the_year
  day    =  gome_pixdate (7:8); READ(day, '(I2)') the_day
  monthc = gome_pixdate(5:6);   READ(monthc, '(I2)') the_month
   
  CALL scia_read_el1radspec (funit, numwin, winpix(1:numwin, :), &
       winlim(1:numwin, :), band_selectors(1:numwin), &
       n_gome_radpts, gome_radspec(1:3, :), nradpix(1:numwin), file_read_stat)
  IF (file_read_stat /= file_read_ok ) RETURN

  ! gome_radspec(3,1:n_gome_radpts) = gome_radspec(3,1:n_gome_radpts) * gome_radspec(2,1:n_gome_radpts) 
  ! Need to get the relative radiance error
  IF (first) THEN
     OPEN(UNIT = raderr_unit, FILE = 'INP/scia_relerr_ozprof.dat', status='old')
     READ(raderr_unit, *)
     READ(raderr_unit, *)
     
     DO i = 1, nwp - 1 
        READ(raderr_unit, *) intercept(i), intercept(i), slope(i) 
     ENDDO
     CLOSE (raderr_unit)
     first = .FALSE.     
  ENDIF

  npts = 0
  DO i = 1, numwin
     npts = npts + winpix(i, 2) - winpix(i, 1) + 1
  ENDDO
  
  IF (npts == n_gome_radpts) THEN
     stpix = 1
     DO i = 1, numwin
        fidx = winpix(i, 1); lidx = winpix(i, 2)
        endpix = stpix + lidx - fidx 
        gome_radspec(3, stpix:endpix) = 10.0 ** (intercept(fidx:lidx) + slope(fidx:lidx) &
             * LOG10(gome_radspec(2, stpix:endpix) * integt)) * 2.0
        stpix = endpix + 1
     ENDDO
     
     !WRITE(90, '(f9.4, 2D12.4)') (gome_radspec(1:3, i), i = 1, n_gome_radpts)
  ELSE
     WRITE(*, *) 'Inconsistent Pixels. Could not get radiance error!!!'; STOP
  ENDIF
  
  spec_norm = SUM(gome_radspec(2,1:n_gome_radpts)) / n_gome_radpts
  IF (spec_norm <= 0.0 ) THEN
     error = .TRUE.;  RETURN
  ENDIF

  gome_radspec(2, 1:n_gome_radpts) = gome_radspec(2, 1:n_gome_radpts) / spec_norm
  gome_radspec(3, 1:n_gome_radpts) = gome_radspec(2, 1:n_gome_radpts) * gome_radspec(3, 1:n_gome_radpts)

  IF (scia_coadd) THEN
     WHERE (gome_radspec(1, 1:n_gome_radpts) > 303.53)
        gome_radspec(3, 1:n_gome_radpts) = gome_radspec(3, 1:n_gome_radpts) / 1.414
     ENDWHERE
  ENDIF
  div_rad = spec_norm
  
  RETURN
  
END SUBROUTINE scia_read_el1data_rad



SUBROUTINE scia_read_el1solspec (funit, nwin, winpix, winlim, bds, &
     ngome, gspec, nsolpix, file_read_stat)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : maxchlen, maxwin, max_fit_pts
  USE OMSAO_errstat_module,    ONLY : file_read_ok, file_read_failed, file_read_eof
  USE OMSAO_gome_data_module,  ONLY : n_gome_data_dim, n_gome_max_pts
  USE OMSAO_variables_module,  ONLY : sol_spec_ring, nsol_ring, scnwrt
  USE ozprof_data_module,      ONLY : sun_posr, sun_specr, nrefl, degcorr

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
  REAL (KIND=dp), DIMENSION(nwin, 2), INTENT(OUT)        :: winlim
  REAL (KIND=dp), DIMENSION(3, max_fit_pts), INTENT(OUT) :: gspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER, PARAMETER          :: maxw = 1
  REAL (KIND=dp), DIMENSION (2048, 3)            :: tmpspec
  REAL (KIND=dp), DIMENSION (maxw, n_gome_max_pts, 3)  :: solspec
  INTEGER, DIMENSION (maxw)   :: numtpix = (/1048/), numpix
  INTEGER                     :: ios, i, j, totpix, ch, fch, lch, iwin, &
       fidx, lidx, sfidx, slidx, errstat
  REAL (KIND=dp)              :: swav, ewav, albcoe, factor
  REAL (KIND=dp), DIMENSION(2):: coeff, wavs0
  CHARACTER (LEN=9)           :: tmpchar
  CHARACTER (LEN=maxchlen)    :: lastline
  
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
 
  ! Read all the data
  DO i = 1, 2048 
     READ (UNIT=funit, FMT=*, IOSTAT = ios) tmpspec(i, :)
     IF (ios /= 0 ) THEN
        WRITE(*, *) 'Read Solar Error!!!'
        STOP
     ENDIF
  ENDDO
  solspec(1, 1:290, :)    = tmpspec(553:842, :)
  solspec(1, 291:1048, :) = tmpspec(1101:1858, :)
  numpix(1) = 1048
  
  solspec(1, 264:1000, :) = solspec(1, 305:1041, :)
  numpix(1) = 1000

  ! Get data for fitting
  ngome = 0; gspec = 0.0
  DO i = 1, nwin
     ch = bds(i)
     IF (numpix(ch) == 0) THEN
        WRITE(*, *) 'No pixels are read from channel ', ch; STOP
     ENDIF
     nsolpix(i) = 0
     DO j = winpix(i, 1), winpix(i, 2)
        IF (solspec(ch, j, 2) /= 0.0 .AND. (solspec(ch, j, 1) > gspec(1, ngome) + 0.1 .OR. ngome == 0)) THEN
           ngome = ngome + 1; nsolpix(i) = nsolpix(i) + 1
           gspec(:, ngome) = solspec(ch, j, :)
        ENDIF
     ENDDO

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
             < sol_spec_ring(1, j) - 0.1) THEN
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
             < sol_spec_ring(1, j) - 0.1) THEN
           j = j - 1
           sol_spec_ring(1:2, j) = solspec(ch, i, 1:2)
           IF (j == 1) EXIT
        ENDIF
     ENDDO
     IF (j == 1) EXIT
  ENDDO

  iwin = 2; ch = bds(iwin); i = winpix(1, 2) + 1
  DO 
     IF (solspec(ch, i, 2) /= 0.0 .AND. solspec(ch, i, 1) > &
          sol_spec_ring(1, nsol_ring) .AND. solspec(ch, i, 1) <= gspec(1, ngome) &
          .AND. i >= winpix(iwin, 1) - 40 .AND. i <= winpix(iwin, 2) + 40 ) THEN
        IF (iwin /= nwin .AND. solspec(ch, i, 1) < &
             solspec(bds(iwin+1), winpix(iwin+1, 1), 1)) THEN
           nsol_ring = nsol_ring + 1
           sol_spec_ring(1:2, nsol_ring) = solspec(ch, i, 1:2)
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
             > sol_spec_ring(1, j) + 0.1) THEN
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
             > sol_spec_ring(1, j) + 0.1) THEN
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
  DO i = 783, numpix(1)
     IF (solspec(1, i, 2) /= 0.0 ) THEN
        j = j + 1
        sun_posr(j) = solspec(1, i, 1); sun_specr(j) = solspec(1, i, 2)
        IF (j == nrefl) EXIT
     ENDIF
  ENDDO

  IF (scnwrt) WRITE(*, *) 'End Of Reading Irradiance Spectrum'
  DO i = 1, nwin
     IF (scnwrt) WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', i, winlim(i,1), winlim(i,2), nsolpix(i)
     IF (nsolpix(i) <= 20) THEN
        IF (scnwrt) WRITE(*, '(A,f8.3,A3,f8.3)') ' Not enough points in window: ', winlim(i,1),&
             ' - ', winlim(i,2)
        STOP
     ENDIF
  ENDDO
 
  CALL gome_check_read_status ( ios, file_read_stat)
  
  RETURN
END SUBROUTINE scia_read_el1solspec


SUBROUTINE scia_read_el1radspec (funit, nwin, winpix, winlim, bds, &
     ngome, gspec, nradpix, file_read_stat)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : maxchlen, maxwin, max_fit_pts
  USE OMSAO_errstat_module,    ONLY : file_read_ok, file_read_failed, &
       file_read_eof, file_read_missing
  USE OMSAO_gome_data_module,  ONLY : n_gome_data_dim, n_gome_max_pts
  USE OMSAO_variables_module,  ONLY : b1ab_div_wav, nsolpix, scnwrt
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
  REAL (KIND=dp), DIMENSION(3, max_fit_pts), INTENT(OUT):: gspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER, PARAMETER          :: maxw = 1
  REAL (KIND=dp), DIMENSION (maxw, n_gome_max_pts, 3)  :: radspec
  INTEGER, DIMENSION (maxw) :: numtpix = (/1048/), numpix
  INTEGER                     :: ios, i, j, totpix, ch, fch, lch, offset
  REAL (KIND=dp)              :: integt, swav, ewav
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

  b1ab_div_wav = 311.0
 
  ! Read all the data 
  radspec = 0.0; numpix = 0
  DO i = 1, maxw  
     CALL skip_to_filemark (funit, 'Wavelength Radiance', lastline, file_read_stat )
     IF ( file_read_stat /= file_read_ok ) THEN
        STOP 'Read Radiance: Cannot find header!!!'
     END IF
     lastline = lastline(20:100)
     READ(lastline, *) swav, ewav, totpix
     ch = 1
     numpix(ch) = numpix(ch) + totpix
        
     DO j = 1, totpix
        READ (UNIT=funit, FMT=*, IOSTAT = ios) radspec(ch, j, 1:3)
        IF (ios /= 0 ) THEN
           WRITE(*, *) 'Read Radiance: Read Error: Channel', ch
           STOP
        ENDIF
     ENDDO
  ENDDO

  radspec(1, 264:1000, :) = radspec(1, 305:1041, :)
  numpix(1) = 1000

  ! Get data for fitting
  ngome = 0 ; gspec = 0.0

  DO i = 1, nwin
     ch = bds(i)
     IF (numpix(ch) == 0) THEN
        WRITE(*, *) 'No pixels are read from channel ', ch; STOP
     ENDIF
     nradpix(i) = 0

     DO j = winpix(i, 1), winpix(i, 2)
        IF (radspec(ch, j, 2) /= 0.0 .AND. (radspec(ch, j, 1) &
             > gspec(1, ngome) + 0.1 .OR. ngome == 0)) THEN
           ngome = ngome + 1; nradpix(i) = nradpix(i) + 1
           gspec(:, ngome) = radspec(ch, j, :)
        ENDIF
     ENDDO
     
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
  DO i = 783, numpix(1)
     IF (radspec(1, i, 2) /= 0.0) THEN
        j = j + 1
        rad_posr(j) = radspec(1, i, 1); rad_specr(j) = radspec(1, i, 2)
        IF (j == nrefl) EXIT
     ENDIF
  ENDDO

  IF (scnwrt) WRITE(*, *) 'End of Reading Radiance Spectrum'
  DO i = 1, nwin
     IF (scnwrt) WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', i, winlim(i,1), winlim(i,2), nradpix(i)
     IF (nradpix(i) <= 20) THEN
        IF (scnwrt) WRITE(*, '(A,f8.3,A3,f8.3)') 'Not enough points in window: ', winlim(i,1),&
             ' - ', winlim(i,2)
        STOP
     ENDIF
     
     IF (nradpix(i) > nsolpix(i)) THEN
        IF (scnwrt) WRITE(*, *) 'Inconsistent pixels between I and F0 in channel ', bds(i)
        file_read_stat = file_read_missing; RETURN
     ENDIF
  ENDDO
  
  CALL gome_check_read_status ( ios, file_read_stat)
  
  RETURN
END SUBROUTINE scia_read_el1radspec
