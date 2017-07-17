
  ! *********************** Modification History *******************
  ! xiong liu, July 2003
  ! 1. Add vgr, vgl, hwl, hwr for voigt profile shape
  ! 2. Add one argument to subroutine undersample
  ! ****************************************************************

SUBROUTINE prepare_databases (n_rad_wvl, curr_rad_wvl, pge_error_status )

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: phase, have_undersampling, n_refwvl, &
       refwvl, database, database_shiwf, do_bandavg, nradpix, numwin,    &
       lo_radbnd, up_radbnd, n_refspec_pts, refspec_orig_data, fitvar_rad_str, &
       curr_sol_spec, n_irrad_wvl, i0sav, refidx_sav, database_save, &
       n_refwvl_sav, refwvl_sav, curr_sol_spec, nsolpix, refsol_idx, &
       radnhtrunc, refnhextra, sol_spec_ring
  USE OMSAO_indices_module,   ONLY: max_calfit_idx, max_rs_idx, mxs_idx, &
       ring1_idx, ring_idx, comm_idx, wvl_idx, spc_idx, solar_idx, com1_idx, &
       us1_idx, us2_idx, no2_t1_idx, so2_idx, bro_idx, hcho_idx, shift_offset
  USE OMSAO_errstat_module,   ONLY: pge_errstat_error
  USE ozprof_data_module,     ONLY: ring_on_line, ozprof_flag, do_tracewf, &
       do_simu, radcalwrt

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,                               INTENT (IN) :: n_rad_wvl
  REAL (KIND=dp), DIMENSION (n_rad_wvl), INTENT (IN) :: curr_rad_wvl

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER, INTENT (OUT) :: pge_error_status

  ! Local variable
  INTEGER        :: ntemp, i, j, npts, fidx, lidx, k, choice, fidx1, lidx1, &
       nextra, nhsolextra, idxoff
  REAL (KIND=dp) :: deltlam

  ! ensure reference spectra having larger wave range than required
  ! and avoid interpolation failure to get garbage values
  ! in gome_pge_fitting_process, special processing is made to enusre
  ! n_rad_wvl is smaller than n_irrad_wvl by radnhextra * 2 in each fitting window
  nextra    = refnhextra * 2
  nhsolextra = radnhtrunc - refnhextra
  n_refwvl = n_rad_wvl + numwin * nextra

  ! Obtain reference position in solar spectra
  fidx = 1
  DO i = 1, numwin
     lidx = fidx + nradpix(i) + nextra - 1
     idxoff = nhsolextra * 2 * i - nhsolextra
     refsol_idx(fidx:lidx) = (/(j, j = fidx + idxoff, lidx + idxoff)/)
     fidx = lidx  + 1
  ENDDO

  fidx = 1; fidx1 = 1
  DO i = 1, numwin
     j = i - 1
     deltlam = curr_sol_spec(wvl_idx, fidx1 + 1) - curr_sol_spec(wvl_idx, fidx1)
     lidx = fidx + nradpix(i) + nextra - 1 ; lidx1 = fidx + nsolpix(i) - 1
     !IF (do_simu .AND. .NOT. radcalwrt) THEN
     !   deltlam = 0.01
     !   refwvl(fidx)   = curr_rad_wvl(fidx - j * nextra) - deltlam * 2.0
     !   refwvl(fidx+1) = curr_rad_wvl(fidx - j * nextra) - deltlam
     !ELSE
     !ENDIF
  
     !IF (do_simu .AND. .NOT. radcalwrt) THEN
     !   deltlam = 0.01
     !   refwvl(lidx-1) = curr_rad_wvl(lidx - i * nextra) + deltlam
     !   refwvl(lidx)   = curr_rad_wvl(lidx - i * nextra) + deltlam * 2.0
     !ELSE
     !ENDIF

     IF (refnhextra > 0) refwvl(fidx:fidx+refnhextra-1) = curr_sol_spec(wvl_idx, refsol_idx(fidx:fidx+refnhextra-1))
     IF (refnhextra > 0) refwvl(lidx-refnhextra+1:lidx) = curr_sol_spec(wvl_idx, refsol_idx(lidx-refnhextra+1:lidx))
     refwvl(fidx+refnhextra:lidx-refnhextra) = curr_rad_wvl(fidx - j * nextra : lidx - i * nextra)
 
     fidx = lidx + 1; fidx1 = lidx1 + 1
  ENDDO

  fidx = 1
  DO i = 1, numwin
     lidx = fidx + nradpix(i) - 1
     idxoff = refnhextra + (i - 1) * nextra
     refidx_sav(fidx:lidx) = (/(j, j = fidx + idxoff, lidx + idxoff)/)
     fidx = lidx  + 1
  ENDDO
  refwvl_sav(1:n_refwvl) = refwvl(1:n_refwvl); n_refwvl_sav = n_refwvl

 
  ! Initialize database
  database = 0.0; database_shiwf = 0.0

  ! --------------------------------------
  ! Calculate the splined fitting database
  ! --------------------------------------
  CALL prepare_refspecs (n_refwvl, refwvl(1:n_refwvl), pge_error_status)
  IF ( pge_error_status >= pge_errstat_error ) RETURN
  i0sav(1:n_refwvl) = database(solar_idx, 1:n_refwvl)

  ! ----------------------------------------------------------------------------
  ! Normalize high-resolution solar reference spectrum with actual solar spectra
  ! ----------------------------------------------------------------------------
  !CALL normalize_solar_refspec (n_refwvl, refwvl(1:n_refwvl), i0sav(1:n_refwvl), pge_error_status)
  !IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! ---------------------------------------------------------
  ! Spline external reference spectra to common radiance grid
  ! ---------------------------------------------------------
  CALL dataspline ( n_refwvl, refwvl(1:n_refwvl), pge_error_status)
  IF ( pge_error_status >= pge_errstat_error) RETURN

  ! ----------------------------------------------------------
  ! Calculate the undersampled spectrum
  ! -----------------------------------------------------------
  CALL undersample (n_refwvl, refwvl(1:n_refwvl), phase, pge_error_status)
  IF ( pge_error_status >= pge_errstat_error ) RETURN


  IF (do_bandavg) THEN
     DO i = 1, max_rs_idx 
        IF ( (i == ring_idx .OR. i == ring1_idx ) &
             .AND. ozprof_flag .AND. ring_on_line) CYCLE 

        IF ( i == comm_idx .OR. i == com1_idx) CYCLE
        
        IF (n_refspec_pts(i) > 0 ) THEN
           CALL avg_band_refspec(refwvl(1:n_refwvl), database(i,1:n_refwvl), &
                n_refwvl, ntemp, pge_error_status)
           IF ( pge_error_status >= pge_errstat_error ) RETURN
        ENDIF
        
        j = shift_offset + i
        IF (n_refspec_pts(i) > 0 .AND. lo_radbnd(j) < up_radbnd(j)) THEN
           CALL avg_band_refspec(refwvl(1:n_refwvl), database_shiwf(i,1:n_refwvl), &
                n_refwvl, ntemp, pge_error_status)    

           IF ( pge_error_status >= pge_errstat_error ) RETURN
        ENDIF
     ENDDO
     
     CALL avg_band_refspec(refwvl(1:n_refwvl), refwvl(1:n_refwvl), &
          n_refwvl, ntemp, pge_error_status)
     IF ( pge_error_status >= pge_errstat_error ) RETURN
     n_refwvl = ntemp  
  ENDIF
 
!  DO i = 1, 2
!  
!     IF (i == 1) THEN 
!        k = comm_idx
!     ELSE 
!        k = com1_idx
!     ENDIF
!     j = shift_offset + k
!
!     npts = n_refspec_pts(k)
!     
!     IF (lo_radbnd(j) < up_radbnd(j) .AND. npts > 4) THEN
!        CALL bspline1(refspec_orig_data(k,1:npts,wvl_idx), refspec_orig_data(k,1:npts, spc_idx), &
!             npts, refwvl(3:n_refwvl-2), database(k, 3:n_refwvl-2), database_shiwf(k, 3:n_refwvl-2), &
!             n_refwvl - 4, pge_error_status)
!        IF (pge_error_status < 0) THEN
!           WRITE(*, *) 'Prepare databases: BSPLINE1 error, errstat = ', pge_error_status; STOP
!        ENDIF
!     ELSE IF (npts > 4) THEN 
!        CALL bspline(refspec_orig_data(k,1:npts,wvl_idx), refspec_orig_data(k,1:npts, spc_idx), &
!             npts, refwvl(3:n_refwvl-2), database(k, 3:n_refwvl-2),  n_refwvl-4, pge_error_status)
!        IF (pge_error_status < 0) THEN
!           WRITE(*, *) 'Prepare databases: BSPLINE error, errstat = ', pge_error_status; STOP
!        ENDIF
!       
!     END IF        
!  ENDDO

  ! save it for fitting weighting function
  IF (do_tracewf) database_save = database

  !DO i = 1, n_refwvl
  !   WRITE(90, *) refwvl(i), database(us1_idx, i), database(us2_idx, i)
  !ENDDO
  !STOP
    
  RETURN
END SUBROUTINE prepare_databases

