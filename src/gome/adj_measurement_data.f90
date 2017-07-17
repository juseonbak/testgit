  ! *********************** Modification History *******************
  ! xiong liu, July 2003
  ! 1. Obtain the effective viewing geometry by averaging those at 
  !    points A, B, and C, which are used for LIDORT calculation 
  ! 2. Use actual measurement weight if use_meas_sig
  ! 3. *********** currently, no downweight is used **************
  ! ****************************************************************
SUBROUTINE adj_solar_data ( error )

  USE OMSAO_precision_module
  USE OMSAO_indices_module,         ONLY: wvl_idx, spc_idx, sig_idx
  USE OMSAO_parameters_module,      ONLY: downweight, normweight
  USE OMSAO_variables_module,       ONLY: curr_sol_spec, n_irrad_wvl, use_meas_sig
  USE OMSAO_gome_data_module,       ONLY: n_gome_solpts, gome_solspec

  IMPLICIT NONE

  ! ================
  ! Output variables
  ! ================
  LOGICAL, INTENT (OUT) :: error

  ! ===============
  ! Local variables
  ! ===============
  INTEGER :: i

  error = .FALSE.

  ! -----------------------------------------------------------
  ! Assign number of irradiance wavelengths to generic variable
  ! -----------------------------------------------------------
  n_irrad_wvl = n_gome_solpts

  ! -----------------------------------------------
  ! Assign irradiance spectrum to generic variables
  ! -----------------------------------------------
  curr_sol_spec(wvl_idx,1:n_irrad_wvl) = gome_solspec(1,1:n_irrad_wvl)
  curr_sol_spec(spc_idx,1:n_irrad_wvl) = gome_solspec(2,1:n_irrad_wvl)

  IF (use_meas_sig) THEN
     curr_sol_spec(sig_idx, 1:n_irrad_wvl) = gome_solspec(3, 1:n_irrad_wvl) 
  ELSE
     curr_sol_spec(sig_idx, 1:n_irrad_wvl) = normweight
  END IF


  RETURN
END SUBROUTINE adj_solar_data


SUBROUTINE adj_earthshine_data ( error )

  ! ************************************************
  !
  !   Read solar spectrum and all radiance spectra
  !
  ! ************************************************

  USE OMSAO_precision_module
  USE OMSAO_indices_module,      ONLY : wvl_idx, spc_idx, sig_idx, solar_idx
  USE OMSAO_parameters_module,   ONLY : deg2rad, rad2deg, downweight, normweight, pi
  USE OMSAO_variables_module,    ONLY : &
       zatmos, sza_atm, vza_atm, aza_atm, amfgeo, szamax, &
       curr_rad_spec, amf, have_amf, sol_zen_eff, n_rad_wvl, &
       the_sza_atm, the_vza_atm, the_aza_atm, the_sca_atm, use_meas_sig, &
       the_month, the_year, the_day, the_lon, the_lat, &
       numwin, nradpix, nsolpix, the_surfalt, actspec_rad, database, refidx, refwvl, &
       radwavcal_freq, yn_varyslit, slit_rad, wavcal_sol, radnhtrunc, refnhextra, &
       nradpix_sav, n_radwvl_sav, radwvl_sav, gome2_idx, instrument_idx
  USE OMSAO_gome_data_module,    ONLY : &
       los_idx, n_gome_ang, n_gome_radpts, azm0_idx, azm_idx, zen0_idx, &
       zen_idx, earth_curv, ers2_alt, n_gome_geo, gome_radspec, &
       gome_angles_wrtn, gome_angles_wrts, gome_geoloc, gome_solspec, gome_curpix
  USE ozprof_data_module,        ONLY : div_rad, div_sun, ozprof_flag, &
       the_cfrac, the_ctp, the_cod, the_orig_cfr, the_orig_ctp, the_orig_cod, &
       aerosol, which_aerosol, &
       radcalwrt, l1l2inp_unit, scale_aod, scaled_aod, do_simu, the_fixalb, &
       do_lambcld, lambcld_refl, has_glint, glintprob, calunit
  USE BOREAS_gome2_data_module,  ONLY : glint_flg
  USE OMSAO_errstat_module 


  IMPLICIT NONE


  ! ================
  ! Output variables
  ! ================
  LOGICAL, INTENT (OUT) :: error

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                                 :: i, j, k, ipix, isub, fidx, lidx, &
       nhtrunc, ntrunc, ntrunc1, idxoff
  REAL (KIND=dp), DIMENSION (n_gome_ang)  :: zen0_n, zen_s, rel_azm
  LOGICAL                                 :: redo_database
  INTEGER                                 :: pge_error_status

  error = .FALSE.

  ! ---------------------------------------------------------
  ! Assign number of radiance wavelengths to generic variable
  ! ---------------------------------------------------------
  n_rad_wvl = n_gome_radpts

  ! Make sure that radiance 

  ! =====================================================
  ! Complete array of wavelengths, radiances, and weights
  ! =====================================================
  curr_rad_spec(wvl_idx,1:n_rad_wvl) = gome_radspec(wvl_idx,1:n_rad_wvl)
  curr_rad_spec(spc_idx,1:n_rad_wvl) = gome_radspec(spc_idx,1:n_rad_wvl)

  IF (use_meas_sig) THEN
     curr_rad_spec(sig_idx, 1:n_rad_wvl) = gome_radspec(3, 1:n_rad_wvl) 
  ELSE
     curr_rad_spec(sig_idx, 1:n_rad_wvl) = normweight
  END IF

  ! Special treatments for sea glint
  IF (instrument_idx /= gome2_idx) glint_flg = 0
  has_glint = .FALSE.; glintprob = 0.0
  ! Land-water flag=1: >=8 not used, else contain water
  IF (glint_flg == 1) THEN
     has_glint = .TRUE.
     CALL SUNGLINT_PROBABILITY (the_sza_atm, the_vza_atm, the_aza_atm, glintprob)
  ENDIF

  IF (do_simu .AND. .NOT. radcalwrt) THEN
     OPEN(UNIT=l1l2inp_unit, FILE='INP/sim.inp', STATUS='unknown')
     READ(l1l2inp_unit, *) the_sza_atm, the_vza_atm, the_aza_atm, the_fixalb, the_surfalt, &
          the_cfrac, the_ctp, the_cod, the_lon, the_lat, the_month, the_day, which_aerosol, &
          scaled_aod, do_lambcld, lambcld_refl
     IF (which_aerosol < 0 ) THEN
        aerosol = .FALSE.
        scale_aod = .FALSE.
        scaled_aod = 0.0
     ELSE
        aerosol = .TRUE.
        scale_aod = .TRUE.
     ENDIF

     IF (the_cfrac == 0.0 .OR. the_cod == 0) THEN
        the_ctp = 0.0; the_cod = 0.0; the_cfrac = 0
     ENDIF
     
     IF (do_lambcld ) THEN
        the_cod = 0.0
     ENDIF

     CLOSE (l1l2inp_unit)
  ENDIF

  ! These properties may be slightly modified later (save them)
  the_orig_cfr = the_cfrac
  the_orig_ctp = the_ctp
  the_orig_cod = the_cod


  ! Make sure that reference spectra has  more wavelengths than
  ! irradiance interpolation and shifting      
  nhtrunc = radnhtrunc
  ntrunc = nhtrunc * 2; ntrunc1 = ntrunc + 1

  fidx = 1
  DO i = 1, numwin
     lidx = fidx + nradpix(i) - ntrunc1
     curr_rad_spec(1:sig_idx, fidx:lidx) = curr_rad_spec(1:sig_idx, fidx + nhtrunc : lidx + nhtrunc)
     !strayspec(1:2, fidx:lidx) = strayspec(1:2, fidx + nhtrunc : lidx + nhtrunc)  
     !strayspec(2, fidx:lidx) = strayspec(2, fidx + nhtrunc : lidx + nhtrunc)    
     IF (lidx  < n_rad_wvl - ntrunc1 ) THEN
        curr_rad_spec(1:sig_idx, lidx+1:n_rad_wvl - ntrunc) =  &
             curr_rad_spec(1:sig_idx, lidx+ntrunc1:n_rad_wvl)
        !strayspec(1:2, lidx+1:n_rad_wvl - ntrunc) =  strayspec(1:2, lidx+ntrunc1:n_rad_wvl)
        !strayspec(2, lidx+1:n_rad_wvl - ntrunc) =  strayspec(2, lidx+ntrunc1:n_rad_wvl)
     ENDIF
     
     nradpix(i) = nradpix(i) - ntrunc; fidx = lidx + 1; n_rad_wvl = n_rad_wvl - ntrunc
  ENDDO
 
  ! save the original grids for later obtain ozone cross section
  n_radwvl_sav = n_rad_wvl; nradpix_sav = nradpix
  radwvl_sav(1:n_rad_wvl) = curr_rad_spec(wvl_idx, 1:n_rad_wvl)  

!!$  redo_database = .FALSE.
!!$  IF (gome_curpix >= 1) THEN
!!$     IF (omi_nwav_rad(currpix, currloop) /= omi_nwav_rad(currpix, currloop-1)) THEN
!!$        redo_database = .TRUE.
!!$     ELSE 
!!$        IF (ANY(radwind(1:n_rad_wvl, currpix, currloop) - &
!!$             radwind(1:n_rad_wvl, currpix, currloop-1) /= 0)) redo_database = .TRUE.
!!$     ENDIF
!!$  ENDIF


  IF ( (gome_curpix == 1) .OR. MOD(gome_curpix, radwavcal_freq) == 0 ) THEN           
     ! --------------------------------------------------
     ! Perform earthshine radiance wavelength calibration
     ! --------------------------------------------------
     PRINT *, 'Performing radiance wavelength calibration'
           
     IF (yn_varyslit) THEN 
        IF (slit_rad) CALL rad_fit_vary(calunit, n_gome_radpts, &
             curr_rad_spec(wvl_idx:sig_idx, 1:n_gome_radpts), error)
              
        IF (.NOT. slit_rad .OR. wavcal_sol) CALL rad_wavcal_vary (calunit, &
             n_gome_radpts, curr_rad_spec(wvl_idx:sig_idx,1:n_gome_radpts), error )
     ELSE 
        CALL radiance_wavcal (n_gome_radpts, &
             curr_rad_spec(wvl_idx:sig_idx,1:n_gome_radpts), error )
     END IF
  ELSE   ! adjust wavelength positions to the calibrated ones           
     CALL rad_shisqu (n_gome_radpts, &
          curr_rad_spec(wvl_idx:sig_idx,1:n_gome_radpts) )          
  ENDIF



  IF ( (gome_curpix == 1) .OR. MOD (gome_curpix, radwavcal_freq) == 0 .OR. redo_database) THEN 
     ! --------------------------------------------------------------
     ! Spline data bases, compute undersampling spectrum, and prepare
     ! reference spectra for fitting.
     ! --------------------------------------------------------------
     CALL prepare_databases ( n_rad_wvl, curr_rad_spec(wvl_idx,1:n_rad_wvl), pge_error_status )
     IF ( pge_error_status >= pge_errstat_error ) RETURN
  ENDIF 


!!$  fidx = 1
!!$  DO i = 1, numwin
!!$     lidx = fidx + nradpix(i) - 1
!!$     refidx(fidx:lidx) = (/(j, j = fidx + 2 + (i - 1) * 4, lidx + 2 + (i - 1) * 4)/)
!!$     fidx = lidx  + 1
!!$  ENDDO

  fidx = 1
  DO i = 1, numwin
     lidx = fidx + nradpix(i) - 1
     idxoff = refnhextra + (i - 1) * 2 * refnhextra
     refidx(fidx:lidx) = (/(j, j = fidx + idxoff, lidx + idxoff)/)
     fidx = lidx  + 1
  ENDDO

  IF (radcalwrt) THEN
     !actspec_rad(1:n_rad_wvl) = curr_rad_spec(spc_idx, 1:n_rad_wvl)
     actspec_rad(1:n_rad_wvl) = curr_rad_spec(spc_idx, 1:n_rad_wvl) / &
          database(solar_idx, refidx(1:n_rad_wvl)) * div_rad / div_sun
  ENDIF

      

  ! ---------------------------------------------------------------------------
  ! Compute shape factor AMF from look-up tables. If the effective solar zenith
  ! angle is >= 90 deg set "HAVE_AMF = .FALSE." and write slant columns to file
  ! ---------------------------------------------------------------------------
  IF (.NOT. ozprof_flag) THEN
     CALL amf_conversion ( n_gome_ang, sza_atm, vza_atm, sol_zen_eff, &
          amf, amfgeo, have_amf )
     
     !   Decide whether to process
     IF ( ALL(sza_atm(1:n_gome_ang) <= szamax) ) THEN
        
        ! Calculate the geometric air mass factor. This can easily be
        ! changed to a look-up of the amf using LIDORT and the
        ! assumption of stratospheric BrO (already implemented
        ! elsewhere), but the difference is minor and until the issue of
        ! widespread tropospheric BrO is resolved, its use is probably
        ! counter-productive.
        amfgeo = ( &
             SUM ( 1.0 / ABS (COS (deg2rad * sza_atm  (1:n_gome_ang))) ) +  &
             SUM ( 1.0 / ABS (COS (deg2rad * vza_atm(1:n_gome_ang))) ) ) &
             / REAL(n_gome_ang, KIND=dp)
     END IF     
  ENDIF

  RETURN
END SUBROUTINE adj_earthshine_data


! Note: need to update the floor noise threshold values
SUBROUTINE adj_rad_sig (radspec, solspec, np)

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx
  USE OMSAO_parameters_module, ONLY: downweight, maxwin
  USE OMSAO_variables_module,  ONLY: numwin, nradpix, band_selectors, &
       do_bandavg, n_band_avg, n_band_samp, the_month, the_year, the_day
  USE ozprof_data_module,      ONLY: div_rad, div_sun, use_lograd, &
       use_flns, b1ab_change, corr_unit
  
  IMPLICIT NONE
  
  ! ====================
  ! In/Output variables
  ! ====================
  INTEGER, INTENT(IN)                             :: np
  REAL (KIND=dp), DIMENSION(3, np), INTENT(IN)    :: solspec
  REAL (KIND=dp), DIMENSION(3, np), INTENT(INOUT) :: radspec
  
  ! ====================
  ! Local variables
  ! ====================
  INTEGER                           :: i, iwin, fidx, lidx
  REAL (KIND=dp), DIMENSION(np)     :: relsig, normrad, sig
  REAL (KIND=dp), DIMENSION(maxwin) :: floor_noise = &
       (/0.0075, 0.005, 0.004, 0.004, 0.004/)
  REAL (KIND=dp), SAVE, DIMENSION(187) :: fudge_factor
  REAL (KIND=dp)                       :: avgmerr
  LOGICAL                              :: use_fudge_factor = .TRUE.
  LOGICAL,SAVE                         :: first = .TRUE.
  
  ! I/F
  normrad = radspec(spc_idx,:) / solspec(spc_idx,:) * (div_rad / div_sun)  
   
  ! Relative photon noise (I)
  relsig = radspec(sig_idx,:) / radspec(spc_idx,:)
 
  ! the estimated noise after changing of bounday is too high
  IF (b1ab_change) THEN  
     IF (use_fudge_factor) THEN
        IF (first) THEN
           OPEN(UNIT=corr_unit, FILE='INP/b1a1b_merr_fudge_factor.dat', STATUS='OLD')
           READ(corr_unit, *) fudge_factor
           CLOSE(corr_unit)
           first = .FALSE.
        ENDIF
        relsig(1:187) = relsig(1:187) / fudge_factor
     ELSE
        WHERE(radspec(1, :) < 307.2 .AND. radspec(1, :) > 282.6)
           relsig = relsig / SQRT(8.0)
        ENDWHERE
     ENDIF
  ENDIF

  relsig = SQRT( relsig ** 2.0 + (solspec(sig_idx,:) / solspec(spc_idx,:)) ** 2.0)
  
  IF (use_flns) THEN
     fidx = 1
     DO iwin = 1, numwin
        lidx = fidx + nradpix(iwin)
        IF (lidx - fidx + 1 > 1) THEN
           !WHERE (relsig(fidx:lidx) < floor_noise(band_selectors(iwin)))               
              relsig(fidx:lidx) = floor_noise(band_selectors(iwin))
           !END WHERE
        ENDIF

        fidx = lidx + 1
     ENDDO
  END IF
 
 IF (use_lograd) THEN        ! error in logarithmic radiance
    sig = LOG(1.0 + relsig)  ! ~relative error
 ELSE
    sig = relsig * normrad   ! absolute measurement error in I/F
 ENDIF
 
 ! Special for edge pixels (2 pixels at each end by assigning a much larger weight)
 !fidx = 1
 !DO iwin =1, numwin
 !   lidx = fidx + nradpix(iwin) - 1
 !WRITE(*, '(3i5,2f8.2)') iwin, fidx, lidx, winlim(iwin,1), winlim(iwin, 2)
 
 !   IF (lidx - fidx >= 10) THEN
 !      sig(fidx:fidx+1) = downweight; sig(lidx-1:lidx) = downweight
 !   ELSE 
 !      WRITE(*, *) 'adj_rad_sig: # wavelengths <= 10 in win ', &
 !           iwin, nradpix(iwin)
 !      STOP
 !   END IF

 !   fidx = lidx + 1
 !END DO
 radspec(sig_idx, :) = sig

 !DO i =1, np
 !   write(90, '(4d14.6)')  radspec(1, i), radspec(2, i), relsig(i), sig(i)
 !ENDDO
 
 RETURN
END SUBROUTINE adj_rad_sig


!SUBROUTINE reduce_resolution (radspec, solspec, np, nnp)
!
!  USE OMSAO_precision_module
!  USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx, hwe_idx, asy_idx
!  USE OMSAO_variables_module,  ONLY: slitfit, nslit
!  
!  IMPLICIT NONE
!  
!  ! ====================
!  ! In/Output variables
!  ! ====================
!  INTEGER, INTENT(IN)                             :: np
!  INTEGER, INTENT(OUT)                            :: nnp
!  REAL (KIND=dp), DIMENSION(3, np), INTENT(INOUT) :: solspec
!  REAL (KIND=dp), DIMENSION(3, np), INTENT(INOUT) :: radspec
!  
!  ! ====================
!  ! Local variables
!  ! ====================
!  INTEGER                           :: i
!  REAL (KIND=dp), DIMENSION(np)     :: specmod
!  REAL (KIND=dp), DIMENSION(3, np)  :: nsolspec
!  REAL (KIND=dp), DIMENSION(3, np)  :: nradspec
!
!  ! slit width 0.8 nm FWHM, 0.5-nm hw1e
!  slitfit(1:nslit, hwe_idx, 1) = 0.6
!  slitfit(1:nslit, asy_idx, 1) = 0.0
!  
!  CALL asym_gauss_vary(radspec(wvl_idx, 1:np), radspec(spc_idx, 1:np), specmod(1:np), np)
!  radspec(spc_idx, 1:np) = specmod(1:np)
!
!  CALL asym_gauss_vary(radspec(wvl_idx, 1:np), radspec(sig_idx, 1:np), specmod(1:np), np)
!  radspec(sig_idx, 1:np) = specmod(1:np)
!
!  CALL asym_gauss_vary(solspec(wvl_idx, 1:np), solspec(spc_idx, 1:np), specmod(1:np), np)
!  solspec(spc_idx, 1:np) = specmod(1:np)
!
!  CALL asym_gauss_vary(solspec(wvl_idx, 1:np), solspec(sig_idx, 1:np), specmod(1:np), np)
!  solspec(sig_idx, 1:np) = specmod(1:np)
!
!  nnp = (np + 1) / 2
!
!  DO i = 1, nnp
!     nsolspec(:, i) = solspec(:, i * 2 - 1)
!     nradspec(:, i) = radspec(:, i * 2 - 1)
!  ENDDO
!
!  radspec(:, 1:nnp) = nradspec(:, 1:nnp)
!  solspec(:, 1:nnp) = nsolspec(:, 1:nnp)
!
!  radspec(3, 1:nnp) = radspec(3, 1:nnp) / 2.0
!  solspec(3, 1:nnp) = solspec(3, 1:nnp) / 2.0
!
! 
! RETURN
!END SUBROUTINE reduce_resolution


SUBROUTINE calc_view_geometry

  ! ************************************************
  !
  !   Read solar spectrum and all radiance spectra
  !
  ! ************************************************

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,   ONLY: deg2rad, rad2deg, pi
  USE OMSAO_variables_module,    ONLY: zatmos, sza_atm, vza_atm,   &
       aza_atm, amfgeo, szamax, sol_zen_eff, the_sza_atm,          &
       the_vza_atm, the_aza_atm, the_sca_atm, &
       instrument_idx, gome_idx, scia_idx, gome2_idx, omi_idx
  USE OMSAO_gome_data_module, ONLY: los_idx, n_gome_ang,           &
       azm0_idx, azm_idx, zen0_idx, zen_idx, earth_curv, ers2_alt, &
       n_gome_geo, gome_angles_wrtn, gome_angles_wrts
  
  IMPLICIT NONE

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                                :: i, j
  REAL (KIND=dp), DIMENSION(n_gome_ang)  :: zen0_n, zen_s, rel_azm, rel_azm_tmp



  ! Get relative azimuthal angle
  ! If rel_azm > 180 then rel_azm=360-rel_azm because -50 is same as 50
  ! GOME-1 and GOME-2 calculations are different for rel_azm because
  ! GOME-1 azimuth is referenced to satellite and GOME-2 to ground.  If
  ! move GOME-1 reference to surface, get 180 deg turnaround. 
  IF (instrument_idx == gome2_idx) THEN
     rel_azm_tmp = ABS(gome_angles_wrtn(azm0_idx,1:n_gome_ang) - &
          gome_angles_wrtn(azm_idx,1:n_gome_ang))
     WHERE (rel_azm_tmp .LE. 180.0) rel_azm = 180.0 - rel_azm_tmp
     WHERE (rel_azm_tmp .GT. 180.0) rel_azm = rel_azm_tmp - 180.0
  ELSE
     rel_azm = ABS( gome_angles_wrtn(azm0_idx,1:n_gome_ang) - &
          gome_angles_wrtn(azm_idx,1:n_gome_ang))
     WHERE (rel_azm > 180.0)
        rel_azm = 360.0 - rel_azm
     END WHERE
  ENDIF


  ! If satellite azimuthal angle wrtn < 180, then negative view angle 
  ! Negative means on the other side of solar zenith angle
  ! The use of negative angle is to obtain the effective view angle
  WHERE(gome_angles_wrtn(azm_idx, 1:n_gome_ang) < 180.0)
     gome_angles_wrts(zen_idx, 1:n_gome_ang) =  &
          -gome_angles_wrts(zen_idx, 1:n_gome_ang) 
  END WHERE
  
  ! -----------------------------------------
  ! Convert angles from degrees to radians
  ! -----------------------------------------
  zen0_n(1:n_gome_ang) = deg2rad * gome_angles_wrtn(zen0_idx,1:n_gome_ang)
  zen_s (1:n_gome_ang) = deg2rad * gome_angles_wrts(zen_idx, 1:n_gome_ang)
  rel_azm              = deg2rad * rel_azm

  ! --------------------------------------------------------------------
  ! Perform correction of TOA; angles are overwritten with corrected values
  ! --------------------------------------------------------------------
  IF (instrument_idx == scia_idx) THEN
     CALL scia_angle_sat2toa (earth_curv, ers2_alt, zatmos, n_gome_ang, zen0_n, zen_s, &
          rel_azm, the_sza_atm, the_vza_atm, the_aza_atm, the_sca_atm )
  ELSE IF (instrument_idx == gome2_idx ) THEN
     CALL gome2_angle_sat2toa (earth_curv, ers2_alt, zatmos, n_gome_ang, zen0_n, zen_s, &
          rel_azm, the_sza_atm, the_vza_atm, the_aza_atm, the_sca_atm )
  ELSE 
     CALL angle_sat2toa (earth_curv, ers2_alt, zatmos, n_gome_ang, zen0_n, zen_s, &
          rel_azm, the_sza_atm, the_vza_atm, the_aza_atm, the_sca_atm )
  ENDIF


  ! ------------------------------
  ! Convert angles back to degrees
  ! ------------------------------
  zen0_n = rad2deg*zen0_n; zen_s = rad2deg*zen_s; rel_azm = rad2deg*rel_azm 
  vza_atm = zen_s  ;  sza_atm = zen0_n; aza_atm = rel_azm

  !WRITE(*, '(A5, 3f8.2)')  'SZA =', zen0_n 
  !WRITE(*, '(A5, 3f8.2)')  'VZA =', zen_s 
  !WRITE(*, '(A5, 3f8.2)')  'AZA =', rel_azm
  !WRITE(*, '(A, 4f8.2)') 'Effective SZA, VZA, AZA, SCA = ', &
  !     the_sza_atm, the_vza_atm, the_aza_atm, the_sca_atm
  !WRITE(*, *) ' '


  RETURN
END SUBROUTINE calc_view_geometry

SUBROUTINE SPIKE_DETECT_CORRECT1(ns, fitspec, simrad)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : fitwavs, fitweights, currspec, poly_order
  USE ozprof_data_module,     ONLY : use_lograd
  
  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN)                             :: ns
  REAL (KIND=dp), INTENT(INOUT), DIMENSION (ns)   :: simrad
  REAL (KIND=dp), INTENT(INOUT), DIMENSION (ns)   :: fitspec
  
  ! =======================
  ! local variables
  ! =======================
  INTEGER                        :: i, j, nspike, ncorr, approach, iter
  REAL (KIND=dp), DIMENSION(ns)  :: dfthresh, diff, mratio, sratio, &
       reldf, reldv, relavg
  REAL (KIND=dp)                 :: rms, thresh
  
  IF (use_lograd) THEN
     simrad = EXP(simrad); fitspec = EXP(fitspec)
  ENDIF
  !WRITE(90, '(f8.3, 2d14.6)') ((fitwavs(i), fitspec(i), simrad(i)), i=1, ns)
  
  sratio(1:ns-1) = simrad(1:ns-1)  /  simrad(2:ns) ; sratio(ns) = sratio(ns-1)
  dfthresh(1:ns-1) = 2.0 * SQRT(fitweights(1:ns-1)**2.0 + fitweights(2:ns)**2.0)
  rms = SUM(dfthresh(1:ns-1)) / (ns-1.0D0)
  
  nspike = 0; ncorr = 0; approach = 1
  
  DO iter = 1, 1     
     mratio(1:ns-1) = fitspec(1:ns-1) / fitspec(2:ns);  mratio(ns) = mratio(ns-1) 
     diff(1:ns-1)   = mratio(1:ns-1) - sratio(1:ns-1)
     
     DO i = ns - 2, 1, - 1 
        thresh = MAX(dfthresh(i), rms)
        IF (ABS(diff(i)) > thresh ) THEN        
           nspike = nspike + 1
        ENDIF
        
        IF (diff(i) < -thresh ) THEN
           ncorr = ncorr + 1
           fitspec (i+1)  = fitspec(i+1)  * mratio(i) / sratio(i)
           currspec(i+1)  = currspec(i+1) * mratio(i) / sratio(i)
           mratio(i) = fitspec(i)/fitspec(i+1); mratio(i+1) = fitspec(i+1)/fitspec(i+2)
           diff(i+1) = mratio(i+1) - sratio(i+1)
        ELSE IF  (diff(i) > thresh) THEN
           ncorr = ncorr + 1
           fitspec(i)  = fitspec(i)  * sratio(i) / mratio(i)
           currspec(i) = currspec(i) * sratio(i) / mratio(i)
           mratio(i-1) = fitspec(i-1) / fitspec(i)
           diff(i-1)   = mratio(i-1) - sratio(i-1)
        ENDIF
     ENDDO
     
     WRITE(*, *) 'Number of spikes = ', nspike, thresh
     WRITE(*, *) 'Number of corrections = ', ncorr
  ENDDO
  
  !WRITE(91, '(f8.3, 2d14.6)') ((fitwavs(i), fitspec(i), simrad(i)), i=1, ns)
  IF (use_lograd) THEN
     simrad = LOG(simrad); fitspec = LOG(fitspec)
  ENDIF
  
  RETURN
  
END SUBROUTINE SPIKE_DETECT_CORRECT1

SUBROUTINE ROUGH_SPIKE_DETECT1(ns, waves, rad, sol)

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN)                             :: ns
  REAL (KIND=dp), INTENT(INOUT), DIMENSION (ns)   :: rad
  REAL (KIND=dp), INTENT(IN), DIMENSION (ns)      :: sol, waves

  ! =======================
  ! local variables
  ! =======================
  INTEGER                       :: i, j, nspike
  REAL (KIND=dp), PARAMETER     :: thresh = 0.15
  REAL (KIND=dp)                :: diff, oldratio, newratio, newratio2
  REAL (KIND=dp), DIMENSION(ns) :: normrad

  !WRITE(90, '(f8.3, 2d14.6)') ((waves(i), rad(i), sol(i)), i=1, ns)

  normrad = rad / sol
  nspike  = 0           
  DO i = ns - 5, 1, - 1      ! only detect spike below 305 nm

     IF (waves(i) > 305.0D0) CYCLE

     oldratio = normrad(i) * normrad(i+2) / (normrad(i+1)**2.0)
     diff = oldratio-1.0

     IF (diff > thresh ) THEN   !  This pixel got large error
        nspike = nspike + 1
        newratio = normrad(i+1) * normrad(i+3) / (normrad(i+2)**2.0)
        
        IF (ABS(newratio - 1.0) < 0.5 * thresh) THEN
           !write(*, *) i, oldratio, normrad(i)           
           normrad(i) = normrad(i+1)**2.0 * newratio / normrad(i+2)
           rad(i) = rad(i) * newratio / oldratio
           !write(*, *) i, newratio, normrad(i)
        ELSE
           !write(*, *) i, oldratio, normrad(i), normrad(i+1)
           newratio2 = normrad(i+2) * normrad(i+4) / (normrad(i+3)**2.0)
           rad(i+1) =  rad(i+1) * newratio2 / newratio
           rad(i) = rad(i) * newratio2 / oldratio
           !write(*, *) i, newratio, normrad(i), normrad(i+1)
        ENDIF
     ENDIF
  ENDDO
  WRITE(*, *) 'Number of spikes = ', nspike
  
  RETURN
END SUBROUTINE ROUGH_SPIKE_DETECT1

