

! These subroutines read GOME-2 level 1-B data and are based on 
! routines by Trevor Beck (NOAA), which are based on routines
! from Eumetsat.
! They have been modified to work in a BOREAS program that is generalized for 
! both GOME-1 and GOME-2 data analysis.
!
! CRN (18/11/2008) This version modified to work with Xiong Liu's ozone profile
! retrieval code which goes by pixel number and not wavelength. It has also been
! updated to read in multiple bands.
! (08/01/2009) Updated to correct for xtrack indexing offset in Level 1B data files,
! reprocessed version 0.
! (04/05/09) Updated again to work with Xiong's newer OMI code with limit inputs
!  in wavelengths, not pixels, and with VLIDORT.
!
! Caroline Nowlan, cnowlan@cfa.harvard.edu, 18-Jun-2008


!*******************************************************************
SUBROUTINE gome2_read_l1b_sol ( fptr, file_read_stat, error )
!*******************************************************************

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxchlen, maxwin, max_fit_pts, mrefl
  USE BOREAS_gome2_data_module, ONLY: g2_smr, read_g2_viadr_smr
  USE OMSAO_gome_data_module, ONLY: gome_radspec, gome_solspec, n_gome_solpts
  USE OMSAO_variables_module, ONLY : numwin, winpix, winlim, &
       band_selectors, nsolpix, sol_spec_ring, nsol_ring, sring_fidx, sring_lidx
  USE OMSAO_errstat_module,    ONLY: &
       file_read_ok, file_read_failed, file_read_missing, file_read_eof
  USE file_handling  ! GOME-2 struct for fptr, the file pointer struct.
  USE ozprof_data_module,     ONLY : pos_alb, toms_fwhm, &
       div_sun, sun_posr, sun_specr, nrefl  


  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  !fptr has  intent(inout) because fptr has current byte offset which is modified for each read !CTB
  TYPE(fileptr), INTENT(INOUT) :: fptr 

  ! ================
  ! Output variables
  ! ================
  LOGICAL, INTENT (OUT) :: error
  INTEGER, INTENT (OUT) :: file_read_stat

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp)   :: swav, ewav, solar_norm
  INTEGER :: i, j      
  INTEGER :: status 
  INTEGER :: ichan 
  INTEGER :: fch, lch, ch, iwin, nwin
  INTEGER, PARAMETER          :: maxw = maxwin - 1
  INTEGER, DIMENSION (maxw)   :: numtpix = (/1024, 1024, 1024, 1024/), numpix


  file_read_stat = file_read_ok ; error = .FALSE.

  ! ----------------------------
  ! Read the GOME-2 solar spectrum
  ! ----------------------------


  !  Read a single viadr smr reference solar record
  CALL read_g2_viadr_smr(fptr, 1, g2_smr, status )
  CALL check_status( status, "gome2_read_l1b CALLed read_g2_viadr_smr")

  n_gome_solpts = 0
  DO i = 1, numwin
     swav = winlim(i,1)
     ewav = winlim(i,2)
     nsolpix(i) = 0
     DO j = 1, SIZE(g2_smr%LAMBDA_SMR(:,band_selectors(i)))
        IF ((g2_smr%LAMBDA_SMR(j,band_selectors(i))  .GE. swav ) .AND. &  
             (g2_smr%LAMBDA_SMR(j,band_selectors(i))  .LE. ewav )) THEN
           gome_solspec(1,n_gome_solpts+1)= g2_smr%LAMBDA_SMR(j,band_selectors(i))
           gome_solspec(2,n_gome_solpts+1)= g2_smr%SMR(j,band_selectors(i))
           gome_solspec(3,n_gome_solpts+1)= g2_smr%E_SMR(j,band_selectors(i))
           gome_solspec(4,n_gome_solpts+1)= g2_smr%E_REL_SUN(j,band_selectors(i)) * 100.0
           nsolpix(i) = nsolpix(i) + 1
           n_gome_solpts = n_gome_solpts+1
           winpix(i,2) = j
           winpix(i,1) = winpix(i,2) - nsolpix(i) + 1
        ENDIF
     ENDDO
  ENDDO


  ! Determine the first and last channel to read
  fch = band_selectors(1)
  IF (winpix(1, 1) - 40 < 0) fch = fch-1
  fch = MAX(1, fch)

  lch = band_selectors(numwin)
  IF (winpix(numwin, 2) + 40 > numtpix(lch) ) lch = lch + 1
  lch = MAX(lch, 2)  ! must have channel for surface albedo
  lch = MIN(lch, 4)


  ! -----------------------------
  ! Get data for Ring calculation
  ! -----------------------------

!!$  ch = 4
!!$  do j=1,1024
!!$     write(*,*) g2_smr%LAMBDA_SMR(j,ch), g2_smr%SMR(j,ch), g2_smr%E_SMR(j,ch)
!!$  enddo
!!$  STOP

  nsol_ring = nsolpix(1) + 40
  sol_spec_ring(1:2, 41:41+nsolpix(1)) = gome_solspec(1:2, 1:nsolpix(1))

  nwin = numwin
  j = 41
  DO
     ch = band_selectors(1)
     DO i = winpix(1,1)-1, 1, -1
        IF (g2_smr%SMR(i,ch) /= 0.0 &
             .AND. g2_smr%LAMBDA_SMR(i,ch) < sol_spec_ring(1,j) - 0.08) THEN
           j = j - 1
           sol_spec_ring(1,j) = g2_smr%LAMBDA_SMR(i,ch)
           sol_spec_ring(2,j) = g2_smr%SMR(i,ch)
           IF (j == 1) EXIT
        ENDIF
     ENDDO
     IF (j == 1) EXIT

     ch = ch - 1
     IF (ch < 1) EXIT

     DO i = 1024, 1, -1
        IF (g2_smr%SMR(i,ch) /= 0.0 &
             .AND. g2_smr%LAMBDA_SMR(i,ch) < sol_spec_ring(1,j) - 0.08) THEN
           j = j - 1
           sol_spec_ring(1,j) = g2_smr%LAMBDA_SMR(i,ch)
           sol_spec_ring(2,j) = g2_smr%SMR(i,ch)
           IF (j == 1) EXIT
        ENDIF
     ENDDO
     IF (j == 1) EXIT
  ENDDO

  iwin = 1; ch = band_selectors(iwin); i = winpix(1, 2) + 1
  DO
     IF (g2_smr%SMR(i,ch) /= 0.0 &
          .AND. g2_smr%LAMBDA_SMR(i,ch) > sol_spec_ring(1, nsol_ring) &
          .AND. i >= winpix(iwin,1) -40 &
          .AND. i <= winpix(iwin,2) + 40) THEN
        IF (iwin /= nwin) THEN
           IF (g2_smr%LAMBDA_SMR(i,ch) < g2_smr%LAMBDA_SMR(winpix(iwin+1,1),band_selectors(iwin+1))) THEN
              nsol_ring = nsol_ring + 1
              sol_spec_ring(1, nsol_ring) = g2_smr%LAMBDA_SMR(i,ch)
              sol_spec_ring(2, nsol_ring) = g2_smr%SMR(i,ch)
           ENDIF
        ELSE IF (iwin == nwin) THEN
           nsol_ring = nsol_ring + 1
           sol_spec_ring(1, nsol_ring) = g2_smr%LAMBDA_SMR(i,ch)
           sol_spec_ring(2, nsol_ring) = g2_smr%SMR(i,ch) 
        ENDIF
     ENDIF
     IF (i == numtpix(ch) .AND. iwin < nwin) THEN
        iwin = iwin + 1; ch = band_selectors(iwin); i = 1
     ELSE IF (i == numtpix(ch)) THEN
        EXIT
     ELSE
        i = i + 1
     ENDIF
  ENDDO



  j = nsol_ring
  nsol_ring = nsol_ring + 40
  DO
     ch = band_selectors(iwin) 
     DO i = winpix(nwin,2)+1, 1024
        IF (g2_smr%SMR(i,ch) /= 0.0 &
             .AND. g2_smr%LAMBDA_SMR(i,ch) > sol_spec_ring(1,j) + 0.08) THEN
           j = j + 1
           sol_spec_ring(1,j) = g2_smr%LAMBDA_SMR(i,ch)
           sol_spec_ring(2,j) = g2_smr%SMR(i,ch)
           IF (j == nsol_ring) EXIT
        ENDIF
     ENDDO

     IF (j == nsol_ring) EXIT

     ch = ch + 1
     IF (ch >= 4) EXIT

     DO i = 1, 1024
        IF (g2_smr%SMR(i,ch) /= 0.0 &
             .AND. g2_smr%LAMBDA_SMR(i,ch) > sol_spec_ring(1,j) + 0.08) THEN
           j = j + 1
           sol_spec_ring(1,j) = g2_smr%LAMBDA_SMR(i,ch)
           sol_spec_ring(2,j) = g2_smr%SMR(i,ch)
           IF (j == nsol_ring) EXIT
        ENDIF
     ENDDO
     IF (j == nsol_ring) EXIT
  ENDDO

  sring_fidx = 1
  sring_lidx = nsol_ring


  ! Get data for surface albedo at 370.2 nm +/- 20 pixels, or thereabouts
  ch = 2
  swav = MAXVAL ( MINLOC ( g2_smr%LAMBDA_SMR(:,ch), MASK = &
          ( g2_smr%LAMBDA_SMR(:,ch) > pos_alb - toms_fwhm * 1.4) ))
  ewav = MAXVAL ( MAXLOC ( g2_smr%LAMBDA_SMR(:,ch), MASK = &
          ( g2_smr%LAMBDA_SMR(:,ch) < pos_alb + toms_fwhm * 1.4) ))
  nrefl = ewav - swav + 1
  sun_posr(1:nrefl)  = g2_smr%LAMBDA_SMR(swav:ewav,ch)
  sun_specr(1:nrefl) = g2_smr%SMR(swav:ewav,ch)

  IF ( file_read_stat /= file_read_ok ) RETURN

  ! -----------------------------------
  ! Normalize solar irradiance spectrum
  ! -----------------------------------
  solar_norm = SUM ( gome_solspec(2,1:n_gome_solpts) ) / REAL ( n_gome_solpts, KIND=dp )
  IF ( solar_norm <= 0.0_dp ) solar_norm = 1.0_dp

  gome_solspec(2,1:n_gome_solpts) = gome_solspec(2,1:n_gome_solpts) / solar_norm
  gome_solspec(3,1:n_gome_solpts) = gome_solspec(3,1:n_gome_solpts) / solar_norm
  sol_spec_ring(2,1:nsol_ring)  = sol_spec_ring(2, 1:nsol_ring)   / solar_norm 
  div_sun = solar_norm


  RETURN

END SUBROUTINE gome2_read_l1b_sol




!*************************************************************
SUBROUTINE gome2_read_l1b_rad ( &
     fptr, iscan, ixtrack, file_read_stat, error )
  !*************************************************************


  USE OMSAO_precision_module
  USE OMSAO_errstat_module

  USE OMSAO_parameters_module, ONLY: maxchlen
  USE OMSAO_errstat_module
  USE OMSAO_variables_module, ONLY : numwin, winpix, winlim
  USE OMSAO_gome_data_module, ONLY: gome_radspec, gome_solspec, n_gome_radpts, &
       gome_curpix, gome_curqual, gome_curscan, gome_geoloc, lat_idx, lon_idx, &
       zen0_idx, azm0_idx, zen_idx, azm_idx, gome_angles_wrtn, gome_angles_wrts, &
       ers2_alt, gome_pixdate, obs_month, n_gome_ang, n_gome_geo, earth_curv, &
       gome_endpix, gome_npix
  USE BOREAS_gome2_data_module, ONLY: n_gome2_scans, gome2_lastscan, &
       n_gome2_ev, gome2_max_negrad, integration_time, gome2_utc_day, gome2_utc_millisec, &
       gome2_utc_year, gome2_curr_ixtrack, gome2_curr_nxtrack, xtrack_index, &
       xtrack_indices, gome2_relazm, cloud, gome2_surf_elevation, &
       gome2_cloud_press, gome2_cloud_frac, pol_ss, pcd_basic, g2_earthview_nextscan, &
       glint_flg
  USE file_handling, ONLY: fileptr
  USE gome_mdr_1b_earthshine_module
  USE control_input, ONLY: bands
  USE ozprof_data_module,      ONLY : pos_alb, toms_fwhm, div_rad, nrefl, &
       coadd_after_b1ab, b1ab_change, rad_posr, rad_specr, the_cfrac, the_ctp, &
       the_snowice, hres_i0
  USE OMSAO_variables_module,  ONLY : the_month, the_year, the_day, numwin,    &
       winlim, winpix, band_selectors, nradpix, npix_fitting, sza_atm, vza_atm,&
       aza_atm, the_sza_atm, the_aza_atm, the_vza_atm, the_sca_atm, the_lat,   &
       the_lon, the_lons, the_lats, the_surfalt, nview, nloc, edgelons, &
       edgelats, szamax, scnwrt

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                              INTENT (IN) :: iscan, ixtrack

  ! ================
  ! Output variables
  ! ================
  LOGICAL, INTENT (INOUT) :: error
  INTEGER, INTENT (INOUT) :: file_read_stat

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp)                      :: swav, ewav, spec_norm
  INTEGER                             :: i, j, g, x1, x2, k
  TYPE(gome_mdr_1b_earthshine_struct) :: g2_earthview
  INTEGER                             :: status
  INTEGER                             :: pge_error_status
  TYPE(fileptr), INTENT(INOUT)        :: fptr 
  INTEGER                             :: actual_idx
  INTEGER, DIMENSION(1)               :: int_index
  INTEGER(KIND=i4)                    :: mon, day, jday, hour, min, sec
  CHARACTER(LEN=25)                   :: datestring
  INTEGER                             :: spix, epix, ch
  LOGICAL                             :: view_coadd, fix_xtrack_idx_shift
  REAL (KIND=dp)                      :: lon1, lon2, templon, templat, temp
  INTEGER, DIMENSION(10)              :: n_coadd
  INTEGER                             :: min_num_recs, max_coadd



  ! FROM GOME2 READER, PASSED THROUGH USE FILE_HANDLING module.
  ! TYPE short_cds_time  
  !     INTEGER (KIND=4) :: day
  !     INTEGER (KIND=8):: msec
  ! END TYPE short_cds_time
  TYPE(short_cds_time) :: utc_time


  file_read_stat = file_read_ok ; error = .FALSE.



  ! check for EOF
  IF ( iscan > n_gome2_scans ) THEN
     file_read_stat=file_read_eof
     RETURN
  ENDIF



  ! There is an indexing problem in how the data is stored in Eumetsat's L1B files.  (This 
  ! alternatively could be somewhere in a Eumetsat reader -- but I get the same
  ! problem using the NOAA/Eumetsat g2_l1b_reader and S&T/ESA BEAT IDL software.)  
  ! It so far affects all processed data version N_0 and R_0.  In order to treat this
  ! issue, this subroutine reads a scan and saves it to use in the next scan, so that
  ! the scan read in uses pixels 2:32, but gets pixel 1 from previous scan
  ! (the problem: all pixels are shifted up one index, so that 
  ! the last scan pixel, often #32 for Band 2B, appears as pixel #1 in the next scan).
  ! (CRN 08-Jan-2009)

  ! Flag index shift problem from L1B, check version in filename.
  IF ((INDEX(fptr%filename, '_R_O_') /= 0) .AND. &
       (INDEX(fptr%filename, '_N_O_') /= 0)) THEN
     fix_xtrack_idx_shift = .FALSE.
  ELSE
     fix_xtrack_idx_shift = .TRUE.
  ENDIF

  ! Also read PMD bands
  bands(7:8) = 1

  IF (fix_xtrack_idx_shift) THEN

     IF ( gome2_lastscan /= iscan ) THEN 

        g2_earthview = g2_earthview_nextscan

        IF ( (gome_curpix == 1) .AND. (ixtrack == 1) ) THEN
           CALL read_g2_earthshine(fptr, iscan, g2_earthview, status )
        ENDIF
        
        IF (iscan .LT. n_gome2_scans) THEN  
           CALL read_g2_earthshine(fptr, iscan+1, g2_earthview_nextscan, status )
           IF ( status .NE. 0 ) THEN 
              file_read_stat=file_read_eof
              error=.TRUE.
              WRITE(*,*) "got error from read, returning EOF"
              RETURN
           END IF
           gome2_lastscan = iscan
           CALL shift_g2_earthview_xtrack(g2_earthview, g2_earthview_nextscan)
        ELSE   !Last scan will not be corrected, so stop retrieval in parent routine
           file_read_stat = file_read_eof  
           RETURN
        ENDIF
     END IF

  ELSE

     IF ( gome2_lastscan /= iscan ) THEN 
        CALL read_g2_earthshine(fptr, iscan+1, g2_earthview_nextscan, status )
        IF ( status .NE. 0 ) THEN 
           file_read_stat=file_read_eof
           error=.TRUE.
           WRITE(*,*) "got error from read, returning EOF"
           RETURN
        END IF
        gome2_lastscan = iscan
     END IF

  ENDIF



!!$  write(*,*) 'bands'
!!$  write(*,*) g2_earthview_nextscan%BAND_PP(1:10,1:10)%RAD
!!$
  open(67)
  do j = 1, 15
     write(67,'(F10.4, 256D16.7)') g2_earthview%WAVELENGTH_PS(j), g2_earthview%BAND_PS(j,1:256)%RAD
  enddo
  close(67)



!!$
!!$
!!$  do j = 1, g2_earthview%NUM_RECS(7)
!!$     write(*,*) g2_earthview%WAVELENGTH_PP(j), g2_earthview%BAND_PP(j)%RAD
!!$  enddo
   


!!$
!!$     write(*,*) '1A', g2_earthview%WAVELENGTH_1A(1)
!!$     write(*,*) '1B', g2_earthview%WAVELENGTH_1B(1)
!!$     write(*,*) '2A', g2_earthview%WAVELENGTH_2A(1)
!!$     write(*,*) '2B', g2_earthview%WAVELENGTH_2B(1)
!!$
!!$
!!$  
!!$ 
!!$     open(76)
!!$     do i=1,g2_earthview%REC_LENGTH(1)
!!$        write(76,*) g2_earthview%WAVELENGTH_1a(i), g2_earthview%band_1a(i,ixtrack)%RAD, &
!!$             g2_earthview%band_1a(i,ixtrack)%ERR_RAD
!!$     enddo
!!$     do i=1,g2_earthview%REC_LENGTH(2)
!!$        write(76,*) g2_earthview%WAVELENGTH_1b(i), g2_earthview%band_1b(i,ixtrack)%RAD, &
!!$             g2_earthview%band_1b(i,ixtrack)%ERR_RAD
!!$     enddo   
!!$     close(76)
!!$     
!!$     open(77)
!!$     do i=1,g2_earthview%REC_LENGTH(3)
!!$        write(77,*) g2_earthview%WAVELENGTH_2a(i), g2_earthview%band_2a(i,ixtrack)%RAD, &
!!$             g2_earthview%band_2a(i,ixtrack)%ERR_RAD
!!$     enddo
!!$     do i=1,g2_earthview%REC_LENGTH(4)
!!$        write(77,*) g2_earthview%WAVELENGTH_2b(i), g2_earthview%band_2b(i,ixtrack)%RAD, &
!!$             g2_earthview%band_2b(i,ixtrack)%ERR_RAD
!!$     enddo   
!!$     close(77)

!!$  
!!$  do i=1,g2_earthview%REC_LENGTH(4)
!!$     write(*,*) g2_earthview%WAVELENGTH_2b(i), g2_earthview%band_2b(i,12)%RAD, &
!!$          g2_earthview%band_2b(i,12)%ERR_RAD
!!$  enddo
     

  ! If integration times are not consistent between bands, the bands and the ground pixels
  ! must be averaged to determine a consistent ground pixel for the retrieval.
  ! Ie, Band 1A may have fewer pixels than Band 2B as it has a longer integration time.  
  n_coadd   = g2_earthview%NUM_RECS * bands / gome2_curr_nxtrack
  max_coadd = MAXVAL(g2_earthview%NUM_RECS, MASK = (bands == 1) ) / gome2_curr_nxtrack


  ! Flag co-adding (Xiong's method)
  IF ((band_selectors(1) > 1) .OR. (.NOT. coadd_after_b1ab .AND. b1ab_change) ) THEN
     view_coadd = .FALSE.
  ELSE
     view_coadd = .TRUE.
  ENDIF

  ! ---------------------
  ! Read in spectral data
  ! ---------------------

  j = 1
  n_gome_radpts = 0
  DO i = 1, numwin
     ! Simple case: only uses one band.
     ! Otherwise, need to look at if there are different integration times in band
     SELECT CASE (band_selectors(i))
     CASE (1)  
        ! Band 1A
        IF (bands(1) == 1 .AND. bands(2) == 0) THEN 
           spix = winpix(i,1)
           epix = winpix(i,2)
           nradpix(i) = epix - spix + 1
           IF (n_coadd(1) .GT. 1) THEN
              x1 = (ixtrack - 1) * n_coadd(1) + 1
              x2 = x1 + n_coadd(1) - 1
              k = 0
              DO g = spix, epix
                 gome_radspec(1,j+k) = g2_earthview%WAVELENGTH_1a(g)
                 gome_radspec(2,j+k) = SUM(g2_earthview%band_1a(g,x1:x2)%RAD)     /n_coadd(1)
                 gome_radspec(3,j+k) = SQRT( SUM(g2_earthview%band_1a(g,x1:x2)%ERR_RAD**2)/n_coadd(1) )
                 gome_radspec(4,j+k) = SUM(g2_earthview%band_1a(g,x1:x2)%STOKES_FRACTION)/n_coadd(1)
                 k = k + 1
              ENDDO
           ELSE
              gome_radspec(1,j:j+nradpix(i)-1)= g2_earthview%WAVELENGTH_1a(spix:epix)
              gome_radspec(2,j:j+nradpix(i)-1)= g2_earthview%band_1a(spix:epix,ixtrack)%RAD
              gome_radspec(3,j:j+nradpix(i)-1)= g2_earthview%band_1a(spix:epix,ixtrack)%ERR_RAD
              gome_radspec(4,j:j+nradpix(i)-1)= g2_earthview%band_1a(spix:epix,ixtrack)%STOKES_FRACTION  
           ENDIF
        ENDIF
        ! Band 1B
        IF (bands(1) == 0 .AND. bands(2) == 1) THEN
           spix = winpix(i,1) - (1024 - g2_earthview%REC_LENGTH(2))
           epix = winpix(i,2) - (1024 - g2_earthview%REC_LENGTH(2))
           nradpix(i) = epix - spix + 1
           IF (n_coadd(2) .GT. 1) THEN
              x1 = (ixtrack - 1) * n_coadd(2) + 1
              x2 = x1 + n_coadd(2) - 1
              k = 0
              DO g = spix, epix
                 gome_radspec(1,j+k) = g2_earthview%WAVELENGTH_1b(g)
                 gome_radspec(2,j+k) = SUM(g2_earthview%band_1b(g,x1:x2)%RAD)     /n_coadd(2)
                 gome_radspec(3,j+k) = SQRT( SUM(g2_earthview%band_1b(g,x1:x2)%ERR_RAD**2)/n_coadd(2) )
                 gome_radspec(4,j+k) = SUM(g2_earthview%band_1b(g,x1:x2)%STOKES_FRACTION)/n_coadd(2)
                 k = k + 1
              ENDDO
           ELSE
              gome_radspec(1,j:j+nradpix(i)-1)= g2_earthview%WAVELENGTH_1b(spix:epix)
              gome_radspec(2,j:j+nradpix(i)-1)= g2_earthview%band_1b(spix:epix,ixtrack)%RAD
              gome_radspec(3,j:j+nradpix(i)-1)= g2_earthview%band_1b(spix:epix,ixtrack)%ERR_RAD
              gome_radspec(4,j:j+nradpix(i)-1)= g2_earthview%band_1b(spix:epix,ixtrack)%STOKES_FRACTION
           ENDIF
        ENDIF
     CASE (2)
        ! Band 2A
        IF (bands(3) == 1 .AND. bands(4) == 0 ) THEN
           spix = winpix(i,1)
           epix = winpix(i,2)
           nradpix(i) = epix - spix + 1
           IF (n_coadd(3) .GT. 1) THEN
              x1 = (ixtrack - 1) * n_coadd(3) + 1
              x2 = x1 + n_coadd(3) - 1
              k = 0
              DO g = spix, epix
                 gome_radspec(1,j+k) = g2_earthview%WAVELENGTH_2a(g)
                 gome_radspec(2,j+k) = SUM(g2_earthview%band_2a(g,x1:x2)%RAD)     /n_coadd(3)
                 gome_radspec(3,j+k) = SQRT( SUM(g2_earthview%band_2a(g,x1:x2)%ERR_RAD**2)/n_coadd(3) )
                 gome_radspec(4,j+k) = SUM(g2_earthview%band_2a(g,x1:x2)%STOKES_FRACTION)/n_coadd(3)
                 k = k + 1
              ENDDO
           ELSE
              gome_radspec(1,j:j+nradpix(i)-1)= g2_earthview%WAVELENGTH_2a(spix:epix)
              gome_radspec(2,j:j+nradpix(i)-1)= g2_earthview%band_2a(spix:epix,ixtrack)%RAD
              gome_radspec(3,j:j+nradpix(i)-1)= g2_earthview%band_2a(spix:epix,ixtrack)%ERR_RAD
              gome_radspec(4,j:j+nradpix(i)-1)= g2_earthview%band_2a(spix:epix,ixtrack)%STOKES_FRACTION
           ENDIF
        ENDIF
        ! Band 2B
        IF (bands(3) == 0 .AND. bands(4) == 1) THEN
           spix = winpix(i,1) - (1024 - g2_earthview%REC_LENGTH(4))
           epix = winpix(i,2) - (1024 - g2_earthview%REC_LENGTH(4))
           nradpix(i) = epix - spix + 1
           IF (n_coadd(4) .GT. 1) THEN
              x1 = (ixtrack - 1) * n_coadd(4) + 1
              x2 = x1 + n_coadd(4) - 1
              k = 0
              DO g = spix, epix
                 gome_radspec(1,j+k) = g2_earthview%WAVELENGTH_2b(g)
                 gome_radspec(2,j+k) = SUM(g2_earthview%band_2b(g,x1:x2)%RAD)     /n_coadd(4)
                 gome_radspec(3,j+k) = SQRT( SUM(g2_earthview%band_2b(g,x1:x2)%ERR_RAD**2)/n_coadd(4) )
                 gome_radspec(4,j+k) = SUM(g2_earthview%band_2b(g,x1:x2)%STOKES_FRACTION)/n_coadd(4)
                 k = k + 1
              ENDDO
           ELSE
              gome_radspec(1,j:j+nradpix(i)-1)= g2_earthview%WAVELENGTH_2b(spix:epix)
              gome_radspec(2,j:j+nradpix(i)-1)= g2_earthview%band_2b(spix:epix,ixtrack)%RAD
              gome_radspec(3,j:j+nradpix(i)-1)= g2_earthview%band_2b(spix:epix,ixtrack)%ERR_RAD
              gome_radspec(4,j:j+nradpix(i)-1)= g2_earthview%band_2b(spix:epix,ixtrack)%STOKES_FRACTION
           ENDIF
        ENDIF
     CASE (3)
        spix = winpix(i,1) - (1024 - g2_earthview%REC_LENGTH(5))
        epix = winpix(i,2) - (1024 - g2_earthview%REC_LENGTH(5))
        nradpix(i) = epix - spix + 1
        IF (n_coadd(5) .GT. 1) THEN
           x1 = (ixtrack - 1) * n_coadd(5) + 1
           x2 = x1 + n_coadd(5) - 1
           k = 0
           DO g = spix, epix
              gome_radspec(1,j+k) = g2_earthview%WAVELENGTH_3(g)
              gome_radspec(2,j+k) = SUM(g2_earthview%band_3(g,x1:x2)%RAD)     /n_coadd(5)
              gome_radspec(3,j+k) = SQRT( SUM(g2_earthview%band_3(g,x1:x2)%ERR_RAD**2)/n_coadd(5) )
              gome_radspec(4,j+k) = SUM(g2_earthview%band_3(g,x1:x2)%STOKES_FRACTION)/n_coadd(5)
              k = k + 1
           ENDDO
        ELSE
           gome_radspec(1,j:j+nradpix(i)-1)= g2_earthview%WAVELENGTH_3(spix:epix)
           gome_radspec(2,j:j+nradpix(i)-1)= g2_earthview%band_3(spix:epix,ixtrack)%RAD
           gome_radspec(3,j:j+nradpix(i)-1)= g2_earthview%band_3(spix:epix,ixtrack)%ERR_RAD
           gome_radspec(4,j:j+nradpix(i)-1)= g2_earthview%band_3(spix:epix,ixtrack)%STOKES_FRACTION
        ENDIF
     CASE (4) 
        spix = winpix(i,1) - (1024 - g2_earthview%REC_LENGTH(6))
        epix = winpix(i,2) - (1024 - g2_earthview%REC_LENGTH(6))
        nradpix(i) = epix - spix + 1
        IF (n_coadd(6) .GT. 1) THEN
           x1 = (ixtrack - 1) * n_coadd(6) + 1
           x2 = x1 + n_coadd(6) - 1
           k = 0
           DO g = spix, epix
              gome_radspec(1,j+k) = g2_earthview%WAVELENGTH_4(g)
              gome_radspec(2,j+k) = SUM(g2_earthview%band_4(g,x1:x2)%RAD)     /n_coadd(6)
              gome_radspec(3,j+k) = SQRT( SUM(g2_earthview%band_4(g,x1:x2)%ERR_RAD**2)/n_coadd(6) )
              gome_radspec(4,j+k) = SUM(g2_earthview%band_4(g,x1:x2)%STOKES_FRACTION)/n_coadd(6)
              k = k + 1
           ENDDO
        ELSE
           gome_radspec(1,j:j+nradpix(i)-1)= g2_earthview%WAVELENGTH_4(spix:epix)
           gome_radspec(2,j:j+nradpix(i)-1)= g2_earthview%band_4(spix:epix,ixtrack)%RAD
           gome_radspec(3,j:j+nradpix(i)-1)= g2_earthview%band_4(spix:epix,ixtrack)%ERR_RAD
           gome_radspec(4,j:j+nradpix(i)-1)= g2_earthview%band_4(spix:epix,ixtrack)%STOKES_FRACTION
        ENDIF
     CASE DEFAULT
        WRITE(*, '(a)') &
             'ERROR...unknown METOP-a GOME-2 band selected:subr gome2_read_l1b_rad ', band_selectors(i)
     END SELECT

     j = j + nradpix(i)
     n_gome_radpts = n_gome_radpts + nradpix(i)

  ENDDO


  ! Normalize radiance values
  spec_norm = SUM(gome_radspec(2,1:n_gome_radpts)) / REAL ( n_gome_radpts, KIND=dp )
  IF ( spec_norm <= 0.0_dp ) spec_norm = 1.0_dp
  gome_radspec(2,1:n_gome_radpts) = gome_radspec(2,1:n_gome_radpts) / spec_norm
  gome_radspec(3,1:n_gome_radpts) = gome_radspec(3,1:n_gome_radpts) / spec_norm
  div_rad = spec_norm




  ! Get data for surface albedo at 370.2 nm +/- 20 pixels from Band 2B
  ch = 2
  swav = MAXVAL ( MINLOC ( g2_earthview%WAVELENGTH_2b, MASK = &
          ( g2_earthview%WAVELENGTH_2b > pos_alb - toms_fwhm * 1.4) ))
  ewav = MAXVAL ( MAXLOC ( g2_earthview%WAVELENGTH_2b, MASK = &
          ( g2_earthview%WAVELENGTH_2b < pos_alb + toms_fwhm * 1.4) ))
  x1 = (ixtrack - 1) * n_coadd(4) + 1
  x2 = x1 + n_coadd(4) - 1
  rad_posr(1:nrefl)  = g2_earthview%WAVELENGTH_2b(swav:ewav)

  DO j = 1, nrefl
     rad_specr(j) = SUM(g2_earthview%band_2b(swav+j-1,x1:x2)%RAD)/n_coadd(4)
  ENDDO

 
  ! ---------------------
  ! Read geolocation data
  ! ---------------------

  ! Find index for reading geolocation data in case of co-adding.
  min_num_recs = MINVAL(g2_earthview%NUM_RECS, MASK=(bands==1))
  actual_idx = 0
  DO i = 1, 6
     IF ( (bands(i) == 1) .AND. (g2_earthview%NUM_RECS(i) == min_num_recs) ) THEN
        actual_idx = g2_earthview%geo_earth_actual%int_index(i)
     ENDIF
  ENDDO
  actual_idx = actual_idx + 1 ! increment because fortran arrays defined to start at 1 and not 0.

  ! Read in geolocation and integration information
  integration_time         = g2_earthview%geo_earth_actual%unique_int(actual_idx)
  gome_geoloc(lat_idx,5)   = g2_earthview%geo_earth_actual%CENTRE_ACTUAL(ixtrack,actual_idx)%latitude
  gome_geoloc(lon_idx,5)   = g2_earthview%geo_earth_actual%CENTRE_ACTUAL(ixtrack,actual_idx)%longitude
  gome_geoloc(lat_idx,1:4) = g2_earthview%geo_earth_actual%CORNER_ACTUAL(ixtrack,:,actual_idx)%latitude
  gome_geoloc(lon_idx,1:4) = g2_earthview%geo_earth_actual%CORNER_ACTUAL(ixtrack,:,actual_idx)%longitude

  ! Angle convention is slightly different than for GOME-1, so viewing zenith angle is actually wrt S/C 
  ! in GOME-1 coordinates, while others are wrt North.
  gome_angles_wrtn(zen0_idx,:) = g2_earthview%geo_earth_actual%SOLAR_ZENITH_ACTUAL (ixtrack,:,actual_idx)
  gome_angles_wrtn(azm0_idx,:) = g2_earthview%geo_earth_actual%SOLAR_AZIMUTH_ACTUAL(ixtrack,:,actual_idx)
  gome_angles_wrts(zen_idx,:) = g2_earthview%geo_earth_actual%SAT_ZENITH_ACTUAL (ixtrack,:,actual_idx)
  gome_angles_wrtn(azm_idx,:) = g2_earthview%geo_earth_actual%SAT_AZIMUTH_ACTUAL(ixtrack,:,actual_idx)
  gome_angles_wrtn(zen_idx,:) = 180.0 - ABS(gome_angles_wrts(zen_idx,:))
  gome2_relazm = 180.0 + gome_angles_wrtn(azm0_idx,2) -  gome_angles_wrtn(azm_idx,2)
  IF ( gome2_relazm .GT.  180.0 )  gome2_relazm = gome2_relazm - 360.0
  IF ( gome2_relazm .LT. -180.0 )  gome2_relazm = gome2_relazm + 360.0



  pcd_basic = g2_earthview%pcd_basic
  cloud = g2_earthview%cloud

  ! Only the lat/lon and angle information above is included for ground pixels 
  ! of different integration times.  Most information that varies across-track
  ! appears in arrays of 32 for each scan.  Need to average if doing co-adding.

  IF (view_coadd) THEN
     max_coadd = 32/gome2_curr_nxtrack
     x1 = (ixtrack - 1) * max_coadd + 1
     x2 = x1 + max_coadd - 1
     ers2_alt = 0.001 * SUM(g2_earthview%GEO_BASIC%SATELLITE_ALTITUDE(x1:x2))/max_coadd  ! convert to kilometre
     gome2_utc_day = g2_earthview%GEO_BASIC%UTC_TIME(x1)%day
     gome2_utc_millisec = SUM(g2_earthview%GEO_BASIC%UTC_TIME(x1:x2)%msec)/max_coadd
     gome2_surf_elevation = SUM(g2_earthview%GEO_EARTH%SURFACE_ELEVATION(x1:x2))/max_coadd
     gome_curscan = ixtrack
     gome2_curr_ixtrack = ixtrack
     xtrack_indices = xtrack_index
     gome2_cloud_press = SUM(cloud%FIT_1(x1:x2))/max_coadd
     gome2_cloud_frac  = SUM(cloud%FIT_2(x1:x2))/max_coadd
     pol_ss%wl_pol_ss  = SUM(g2_earthview%pol_ss(x1:x2)%wl_pol_ss) /max_coadd
     pol_ss%p_pol_ss   = SUM(g2_earthview%pol_ss(x1:x2)%p_pol_ss)  /max_coadd
     pol_ss%chi_pol_ss = SUM(g2_earthview%pol_ss(x1:x2)%chi_pol_ss)/max_coadd
     pol_ss%q_pol_ss   = SUM(g2_earthview%pol_ss(x1:x2)%q_pol_ss)  /max_coadd
     pol_ss%u_pol_ss   = SUM(g2_earthview%pol_ss(x1:x2)%u_pol_ss)  /max_coadd
     the_snowice       = SUM(cloud%FIT_MODE(x1:x2)*100)/max_coadd
  ELSE
     ers2_alt = 0.001 *  g2_earthview%GEO_BASIC%SATELLITE_ALTITUDE(ixtrack)  ! convert to kilometre
     utc_time = g2_earthview%GEO_BASIC%UTC_TIME(ixtrack)	
     gome2_utc_day = utc_time%day
     gome2_utc_millisec = utc_time%msec
     gome2_surf_elevation = g2_earthview%GEO_EARTH%SURFACE_ELEVATION(ixtrack)
     !gome_curqual=-2.71828 !FIXME    !probably no exact analogue for gome2 ? 
     gome_curscan = ixtrack
     gome2_curr_ixtrack = ixtrack
     xtrack_indices = xtrack_index
     !IF ( gome2_curr_nxtrack < 32 )  xtrack_indices(gome2_curr_nxtrack+1:size(xtrack_indices))=0
     gome2_cloud_press = cloud%FIT_1(ixtrack)
     gome2_cloud_frac = cloud%FIT_2(ixtrack)
     pol_ss = g2_earthview%pol_ss(ixtrack)
     the_snowice = cloud%FIT_MODE(ixtrack)*100
  ENDIF

  ! Snow/ice flag given as 0 or 1, but changed to 0 (snow/ice free) or 103 (snow) for consistency with OMI.
  ! Co-added pixels are correctly represented by an ice fraction.
  IF (the_snowice == 100) the_snowice = 103

  the_ctp = gome2_cloud_press
  the_cfrac = gome2_cloud_frac
  glint_flg = pcd_basic%F_SUNGLINT

  earth_curv = 0.001 * g2_earthview%GEO_EARTH%EARTH_RADIUS

  CALL epoch_time_to_utc(gome2_utc_day, gome2_utc_millisec, gome2_utc_year, mon, & 
       day, jday, hour, min, sec, gome_pixdate)

  the_year = gome2_utc_year
  the_month = mon
  the_day = day

  obs_month = mon



  !Something weird is going on when I call epoch_time_to_utc --> messes up n_gome_radpts
  ! and resets it to 32!!
  !n_gome_radpts = winpix(1,2) - winpix(1,1) + 1

  n_gome_radpts = 0
  DO i = 1, numwin
     n_gome_radpts = n_gome_radpts + nradpix(i)
  ENDDO


  ! This is the last pixel used in the retrieval.  Note: if we are doing pixel-by-pixel,
  ! it just represents the current pixel number.  This is required for Xiong's ascii writing routine
  ! (not important for HDF writing).
  gome_endpix = gome_curpix
  gome_npix = 1



  ! -----------------------------------------
  ! Some extra processing of geolocation data
  ! -----------------------------------------

  ! Calculate angles with respect to S/C for consistent processing with GOME-1 data
   
  CALL calc_view_geometry
  IF ( the_sza_atm < 0.0 .OR. the_sza_atm >= szamax) THEN
     file_read_stat = file_read_missing; RETURN
  END IF

  the_surfalt = 0.001 * gome2_surf_elevation

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




  ! -----------------------------------------------------------
  ! We are at the end of the current ground pixel input. RETURN
  ! -----------------------------------------------------------
  RETURN

END SUBROUTINE gome2_read_L1B_rad


  


!******************************************************
FUNCTION n_gome2_xtrack(irecord, fptr) RESULT (nXtrack)
  !******************************************************

  USE control_input, ONLY: bands
  USE OMSAO_variables_module,  ONLY : numwin, winpix, band_selectors
  USE file_handling
  USE gome_mdr_1b_earthshine_module

  IMPLICIT NONE

  TYPE(fileptr), INTENT(INOUT) :: fptr
  INTEGER,       INTENT(INOUT)    :: irecord

  INTEGER :: nXtrack
  INTEGER :: i, iband
  INTEGER :: record_size
  INTEGER :: status
  TYPE(gome_mdr_1b_earthshine_struct) :: g2_ev_record

  ! Read record size from header
  record_size = earthview_grh(irecord)%record_size


  ! If we are using only Bands 2B, 3, or 4, can do a fast method by looking at header size.
  IF ( (bands(1) == 0 .AND. bands(2) == 0 .AND. bands(3) == 0) .AND. &
       (SUM(bands) == 1) ) THEN
     SELECT CASE (record_size)
     CASE( 263506) 
        nXtrack = 4 
     CASE( 283594) 
        nXtrack = 4 
     CASE ( 325618 )
        nXtrack = 4
     CASE (480130  )
        nXtrack = 8 
     CASE ( 1438774 )
        nXtrack = 32 
     CASE DEFAULT
        iband = MAXLOC(bands, 1)
        write (*, '(A, I10)') &
             'Cannot determine number of xtracks from record size in bytes.: Revert to slower method: ', record_size
        nXtrack=n_xtrack_per_band(irecord, iband, fptr) 
     END SELECT
  ELSE
     CALL read_g2_earthshine(fptr, irecord, g2_ev_record, status, .TRUE.)
     nXtrack = MINVAL(g2_ev_record%NUM_RECS, MASK = (bands == 1) )
  ENDIF


  RETURN

END FUNCTION n_gome2_xtrack


!*********************************
SUBROUTINE which_gome2_bands(fptr)
  !*********************************

  ! ----------------------------------------------------------------------
  ! This routine is used to determine which bands are being used
  ! in the SAOPROF code.  It should only be used once at the beginning of 
  ! an orbit, so that we don't read in each ground pixel an extra time.
  ! ----------------------------------------------------------------------

  USE control_input, ONLY: bands
  USE OMSAO_variables_module,  ONLY : numwin, winpix, band_selectors
  USE file_handling
  USE gome_mdr_1b_earthshine_module

  IMPLICIT NONE

  TYPE(fileptr), INTENT(INOUT) :: fptr

  INTEGER :: i, status, irecord
  TYPE(gome_mdr_1b_earthshine_struct) :: g2_earthview

 
  irecord = 1
  CALL read_g2_earthshine(fptr, irecord, g2_earthview, status )


  ! Determine which GOME-2 bands are being used
  bands(1:10)=(/0,0,0,0,0,0,0,0,0,0/)
  DO i = 1, numwin
     IF (band_selectors(i) == 1) THEN
        IF ( (winpix(i,1) < g2_earthview%REC_LENGTH(1)) .AND. &
             (winpix(i,2) < g2_earthview%REC_LENGTH(1)) ) THEN
           bands(1) = 1  !Band 1A
        ELSEIF ( (winpix(i,1) >= g2_earthview%REC_LENGTH(1)) .AND. &
             (winpix(i,2) < 1024 ) ) THEN
           bands(2) = 1  !Band 1B
        ELSEIF ( (winpix(i,1) < g2_earthview%REC_LENGTH(1)) .AND. &
             (winpix(i,2) > g2_earthview%REC_LENGTH(1) ) ) THEN
           bands(1) = 1  !Bands 1A and 1B
           bands(2) = 1
        ENDIF
     ELSEIF (band_selectors(i) == 2) THEN
        IF ( (winpix(i,1) < g2_earthview%REC_LENGTH(3)) .AND. &
             (winpix(i,2) < g2_earthview%REC_LENGTH(3) ) ) THEN
           bands(3) = 1  !Band 1A
        ELSEIF ( (winpix(i,1) >= g2_earthview%REC_LENGTH(3)) .AND. &
             (winpix(i,2) < 1024 ) ) THEN
           bands(4) = 1  !Band 1B
        ELSEIF ( (winpix(i,1) < g2_earthview%REC_LENGTH(3) ) .AND. &
             (winpix(i,2) > g2_earthview%REC_LENGTH(3) ) ) THEN
           bands(3) = 1  !Bands 1A and 1B
           bands(4) = 1
        ENDIF
     ELSEIF (band_selectors(i) == 3) THEN
        bands(5) = 1
     ELSEIF (band_selectors(i) == 4) THEN
        bands(6) = 1
     ENDIF
  ENDDO

  RETURN
  
END SUBROUTINE which_gome2_bands


!************************************************************************
SUBROUTINE shift_g2_earthview_xtrack(ev, ev_next)
!************************************************************************

  USE gome_mdr_1b_earthshine_module

  IMPLICIT NONE

  TYPE(gome_mdr_1b_earthshine_struct), INTENT(INOUT) :: ev
  TYPE(gome_mdr_1b_earthshine_struct), INTENT(IN)    :: ev_next

  INTEGER :: i, j, n

  ev%band_1a(:,1:31) = ev%band_1a(:,2:32)
  ev%band_1b(:,1:31) = ev%band_1b(:,2:32)
  ev%band_2a(:,1:31) = ev%band_2a(:,2:32)
  ev%band_2b(:,1:31) = ev%band_2b(:,2:32)
  ev%band_3(:,1:31)  = ev%band_3(:,2:32)
  ev%band_4(:,1:31)  = ev%band_4(:,2:32)
  !ev%band_pp(:,1:248)  = ev%band_pp(:,9:256)
  !ev%band_ps(:,1:248)  = ev%band_ps(:,9:256)
  ev%geo_earth_actual%CENTRE_ACTUAL       (1:31,:)   = ev%geo_earth_actual%CENTRE_ACTUAL(2:32,:)
  ev%geo_earth_actual%CORNER_ACTUAL       (1:31,:,:) = ev%geo_earth_actual%CORNER_ACTUAL(2:32,:,:)
  ev%geo_earth_actual%SOLAR_ZENITH_ACTUAL (1:31,:,:) = ev%geo_earth_actual%SOLAR_ZENITH_ACTUAL(2:32,:,:) 
  ev%geo_earth_actual%SOLAR_AZIMUTH_ACTUAL(1:31,:,:) = ev%geo_earth_actual%SOLAR_AZIMUTH_ACTUAL(2:32,:,:)
  ev%geo_earth_actual%SAT_ZENITH_ACTUAL   (1:31,:,:) = ev%geo_earth_actual%SAT_ZENITH_ACTUAL(2:32,:,:)
  ev%geo_earth_actual%SAT_AZIMUTH_ACTUAL  (1:31,:,:) = ev%geo_earth_actual%SAT_AZIMUTH_ACTUAL(2:32,:,:)
  ev%GEO_BASIC%SATELLITE_ALTITUDE(1:31) = ev%GEO_BASIC%SATELLITE_ALTITUDE(2:32)
  ev%GEO_BASIC%UTC_TIME(1:31) = ev%GEO_BASIC%UTC_TIME(2:32) 
  ev%GEO_EARTH%SURFACE_ELEVATION(1:31) = ev%GEO_EARTH%SURFACE_ELEVATION(2:32)
  ev%cloud%FIT_MODE(1:31) = ev%cloud%FIT_MODE(2:32)
  ev%cloud%FAIL_FLAG(1:31) = ev%cloud%FAIL_FLAG(2:32)
  ev%cloud%FIT_1(1:31) = ev%cloud%FIT_1(2:32)
  ev%cloud%FIT_2(1:31) = ev%cloud%FIT_2(2:32)
  ev%cloud%E_FIT_1(1:31) = ev%cloud%E_FIT_1(2:32)
  ev%cloud%E_FIT_2(1:31) = ev%cloud%E_FIT_2(2:32)
  ev%cloud%FINAL_CHI_SQUARE(1:31) = ev%cloud%FINAL_CHI_SQUARE(2:32)
  ev%cloud%CLOUD_ALBEDO(1:31) = ev%cloud%CLOUD_ALBEDO(2:32)
  ev%cloud%SURFACE_ALBEDO(1:31,:) = ev%cloud%SURFACE_ALBEDO(2:32,:)
  ev%cloud%SURFACE_PRESSURE(1:31) = ev%cloud%SURFACE_PRESSURE(2:32)
  ev%pol_ss(1:31) = ev%pol_ss(2:32)


  DO i = 1, 6
     n = ev%NUM_RECS(i)
     j = ev%geo_earth_actual%int_index(i) + 1 ! Increment 1 as Fortran array indexing starts at 1, not 0
     IF (i==1) ev%band_1a(:,n) = ev_next%band_1a(:,1)
     IF (i==2) ev%band_1b(:,n) = ev_next%band_1b(:,1)
     IF (i==3) ev%band_2a(:,n) = ev_next%band_2a(:,1)
     IF (i==4) ev%band_2b(:,n) = ev_next%band_2b(:,1)
     IF (i==5) ev%band_3(:,n)  = ev_next%band_3(:,1)
     IF (i==6) ev%band_4(:,n)  = ev_next%band_4(:,1)
     !IF (i==7) ev%band_pp(:,(n-7):n) = ev_next%band_pp(:,1:8)
     !IF (i==8) ev%band_ps(:,(n-7):n)  = ev_next%band_ps(:,1:8)
     ev%geo_earth_actual%CENTRE_ACTUAL       (n,j)   = ev_next%geo_earth_actual%CENTRE_ACTUAL(1,j)
     ev%geo_earth_actual%CORNER_ACTUAL       (n,:,j) = ev_next%geo_earth_actual%CORNER_ACTUAL(1,:,j)
     ev%geo_earth_actual%SOLAR_ZENITH_ACTUAL (n,:,j) = ev_next%geo_earth_actual%SOLAR_ZENITH_ACTUAL(1,:,j) 
     ev%geo_earth_actual%SOLAR_AZIMUTH_ACTUAL(n,:,j) = ev_next%geo_earth_actual%SOLAR_AZIMUTH_ACTUAL(1,:,j)
     ev%geo_earth_actual%SAT_ZENITH_ACTUAL   (n,:,j) = ev_next%geo_earth_actual%SAT_ZENITH_ACTUAL(1,:,j)
     ev%geo_earth_actual%SAT_AZIMUTH_ACTUAL  (n,:,j) = ev_next%geo_earth_actual%SAT_AZIMUTH_ACTUAL(1,:,j)
     ev%GEO_BASIC%SATELLITE_ALTITUDE(32) = ev_next%GEO_BASIC%SATELLITE_ALTITUDE(1)
     ev%GEO_BASIC%UTC_TIME(32) = ev_next%GEO_BASIC%UTC_TIME(1) 
     ev%GEO_EARTH%SURFACE_ELEVATION(32) = ev_next%GEO_EARTH%SURFACE_ELEVATION(1)
     ev%cloud%FIT_MODE(32) = ev_next%cloud%FIT_MODE(1)
     ev%cloud%FAIL_FLAG(32) = ev_next%cloud%FAIL_FLAG(1)
     ev%cloud%FIT_1(32) = ev_next%cloud%FIT_1(1)
     ev%cloud%FIT_2(32) = ev_next%cloud%FIT_2(1)
     ev%cloud%E_FIT_1(32) = ev_next%cloud%E_FIT_1(1)
     ev%cloud%E_FIT_2(32) = ev_next%cloud%E_FIT_2(1)
     ev%cloud%FINAL_CHI_SQUARE(32) = ev_next%cloud%FINAL_CHI_SQUARE(1)
     ev%cloud%CLOUD_ALBEDO(32) = ev_next%cloud%CLOUD_ALBEDO(1)
     ev%cloud%SURFACE_ALBEDO(32,:) = ev_next%cloud%SURFACE_ALBEDO(1,:)
     ev%cloud%SURFACE_PRESSURE(32) = ev_next%cloud%SURFACE_PRESSURE(1)
     ev%pol_ss(32) = ev_next%pol_ss(1)
  ENDDO


  RETURN

END SUBROUTINE shift_g2_earthview_xtrack



