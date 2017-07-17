SUBROUTINE prepare_refspecs (n_radpts, curr_rad_wvl, pge_error_status )

  ! ***********************************************************
  !
  !   Calculate the splined fitting database.
  !   Note that the undersampled spectrum has just been done.
  !   Finish filling in database array.
  !
  ! ***********************************************************

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: max_rs_idx, solar_idx, ring_idx, &
       us1_idx, us2_idx, wvl_idx, spc_idx, max_calfit_idx, mxs_idx, ring1_idx
  USE OMSAO_parameters_module, ONLY: maxchlen
  USE OMSAO_variables_module,  ONLY: curr_sol_spec, yn_doas, rad_winwav_idx, &
       n_refspec_pts, yn_smooth, database, database_shiwf, &
       n_irrad_wvl, lo_radbnd, up_radbnd, refsol_idx
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! *******************************************************************
  ! CAREFUL: Assumes that radiance and solar wavelength arrays have the
  ! same number of points. That must not be the case if we read in a
  ! general EL1 file. Examine and adjust! (tpk, note to himself)
  ! *******************************************************************

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,                              INTENT (IN)    :: n_radpts
  REAL (KIND=dp), DIMENSION (n_radpts), INTENT (IN)    :: curr_rad_wvl
  INTEGER,                              INTENT (INOUT) :: pge_error_status

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER :: i, j, npts, ll_rad, lu_rad, errstat
  REAL  (KIND=dp), DIMENSION (n_radpts) :: spline_sun

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=16), PARAMETER :: modulename = 'prepare_refspecs'


  ! Spline irradiance spectrum onto radiance grid

  ! PROBLEM: first wavlength in POS is smaller than first wavelength in
  !          curr_sol_spec(wvl_idx,*).
  ! SOLUTION: Don't include POS(1) in the interpolation, and assign 
  !           SPLINE_SUN(1) = SPEC_SUN(1)

  errstat = pge_errstat_ok
 
  CALL interpolation (n_irrad_wvl, curr_sol_spec(wvl_idx,1:n_irrad_wvl),        &
       curr_sol_spec(spc_idx,1:n_irrad_wvl), n_radpts, curr_rad_wvl(1:n_radpts), &
       spline_sun(1:n_radpts), errstat )  
  database(solar_idx, 1:n_radpts) = spline_sun(1:n_radpts)
     !    print  * , spline_sun(1:n_radpts)
  !database(solar_idx, 1:n_radpts) = curr_sol_spec(spc_idx,refsol_idx(1:n_radpts))

  ! =================================================================
  ! Note that the UNDERSAMPLING spectrum has already been assigned to
  ! DATABASE(us1/2_idx,*) in the UNDERSPEC routine.
  ! =================================================================

  ! ---------------------------------------------------------------------
  ! Set up Ring spectrum for high pass filtering of divided spectrum
  ! (afterward, add solar spectrum back in order to divide by solar
  ! spectrum with altered wavelength calibration in subroutine spectrum).
  ! ---------------------------------------------------------------------
  !IF ( yn_doas ) database(ring_idx, 1:n_radpts) = &
  !     database(ring_idx, 1:n_radpts) / spline_sun(1:n_radpts)
  !
  !! ---------------------------------------------------------
  !! For the DOAS case, high-pass filter the reference spectra
  !! ---------------------------------------------------------
  !IF ( yn_doas ) THEN
  !   ll_rad = rad_winwav_idx(2) ; lu_rad = rad_winwav_idx(3)
  !   CALL subtract_cubic (curr_rad_wvl, n_radpts, 1, n_radpts)
  !   database(ring_idx, 1:n_radpts) = &
  !        database(ring_idx, 1:n_radpts) * spline_sun(1:n_radpts)
  !END IF

  !IF ( yn_smooth ) THEN   ! 1/16, 1/4, 3/8, 1/4, 1/16
  !   DO i = 1, max_rs_idx
  !      IF (n_refspec_pts(i) > 0) &
  !         database(i, 3:n_radpts-2) = 0.375  * database(i, 3:n_radpts-2) + &
  !         0.25 * (database(i, 4:n_radpts-1) + database(i, 2:n_radpts-3)) + &
  !         0.0625 * (database(i, 5:n_radpts)   + database(i, 1:n_radpts-4))
  !  
  !      j = max_calfit_idx + max_rs_idx * mxs_idx + i
  !      IF (lo_radbnd(j) < up_radbnd(j) .AND. n_refspec_pts(i) > 0) &
  !        database_shiwf(j,3:n_radpts-2) = 0.375 * &
  !        database_shiwf(j,3:n_radpts-2) + 0.25*(database_shiwf(j,4:n_radpts-1)   &
  !        + database_shiwf(j,2:n_radpts-3))+0.0625*(database_shiwf(j,5:n_radpts)&
  !        + database_shiwf(j, 1:n_radpts-4))
  !   END DO
  !
  !END IF

  pge_error_status = MAX ( errstat, pge_error_status )

  RETURN
END SUBROUTINE prepare_refspecs


SUBROUTINE load_omi_comres(errstat)
  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: max_fit_pts, maxchlen
  USE OMSAO_variables_module,  ONLY: refdbdir, n_refwvl, refwvl, database, &
       nradpix, numwin, refidx, lo_radbnd, up_radbnd, currpix, fitvar_rad_str, &
       mask_fitvar_rad, n_fitvar_rad
  USE OMSAO_indices_module,    ONLY: com1_idx, comm_idx, comfidx, cm1fidx
  USE OMSAO_errstat_module,    ONLY: pge_errstat_error, pge_errstat_ok
  
  INTEGER, INTENT (OUT)                                :: errstat
  
  CHARACTER (LEN=15), PARAMETER                        :: modulename = 'load_omi_comres'
  INTEGER, PARAMETER                                   :: nx = 30, lun = 12
  INTEGER                                              :: ix, itemp, i, fidx, lidx
  CHARACTER (LEN=maxchlen)                             :: comres_fname
  
  LOGICAL,                                        SAVE :: first = .TRUE.
  REAL (KIND=dp), DIMENSION (nx, max_fit_pts, 2), SAVE :: comres
  INTEGER, DIMENSION(nx, 3),                      SAVE :: npts
    
  errstat = pge_errstat_ok
  
  IF (first) THEN
     comres = 0.D0
     
     !comres_fname = ADJUSTL(TRIM(refdbdir)) // 'OMI_jul11tropopt1com_comres.dat'
     !comres_fname = ADJUSTL(TRIM(refdbdir)) // 'OMI_hres8s48lncncrs1_L0400-0599_comres.dat'
     comres_fname = ADJUSTL(TRIM(refdbdir)) // 'OMI_hresjul11-30S-30N_comres.dat'

     OPEN (UNIT=lun, FILE=TRIM(ADJUSTL(comres_fname)), STATUS='UNKNOWN', IOSTAT=errstat)
     IF ( errstat /= pge_errstat_ok ) THEN
        WRITE(*, '(2A)') modulename, ': Cannot open common-mode residual file!!!'
        errstat = pge_errstat_error; RETURN
     END IF
     
     READ (lun, *) 
     READ (lun, *)
     DO ix = 1, nx
        READ (lun, *) itemp, npts(ix, 1:3)
        DO i = 1, npts(ix, 1)
           READ (lun, *) comres(ix, i, 1:2)
        ENDDO
     ENDDO
     CLOSE (lun)
     
     first = .FALSE.
  ENDIF
  
  IF (comfidx > 0) THEN
     database(comm_idx, 1:n_refwvl) = 0.D0
     fidx = 1; lidx = fidx + nradpix(1) - 1
     database(comm_idx, refidx(fidx:lidx)) = comres(currpix, fidx:lidx, 2)
  ENDIF

  IF ( cm1fidx > 0 ) THEN
     database(com1_idx, 1:n_refwvl) = 0.D0
     fidx = nradpix(1) + 1; lidx = fidx + nradpix(2) - 1
     database(com1_idx, refidx(fidx:lidx)) = comres(currpix, fidx:lidx, 2)
  ENDIF

  RETURN
  
END SUBROUTINE load_omi_comres
