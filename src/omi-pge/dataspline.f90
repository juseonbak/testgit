SUBROUTINE dataspline ( n_radwvl, curr_rad_wvl, errstat)

  USE OMSAO_precision_module
  USE OMSAO_indices_module,     ONLY: max_calfit_idx, mxs_idx, max_rs_idx, &
       solar_idx, wvl_idx, spc_idx, refspec_strings, comm_idx, com1_idx, &
       ring1_idx, ring_idx, so2_idx, hcho_idx, bro_idx, no2_t1_idx, shift_offset, &
       fsl_idx, rsl_idx, bro2_idx, so2v_idx, o2o2_idx
  USE OMSAO_parameters_module,  ONLY: maxchlen, max_spec_pts, zerospec_string
  USE OMSAO_variables_module,   ONLY: n_refspec_pts, refspec_orig_data,    &
       refspec_fname, database, yn_varyslit, database_shiwf,  which_slit, &
       lo_radbnd, up_radbnd, nslit, nslit_rad, nslit_sol, slitwav, slitwav_sol, &
       slitwav_rad, slitfit, solslitfit, radslitfit,  slit_rad, refspec_norm, & 
       instrument_idx, omps_idx, omi_idx
  USE OMSAO_slitfunction_module
  !USE OMPS_slit_module
  USE ozprof_data_module,       ONLY: ring_convol, ring_on_line       
  USE OMSAO_errstat_module
  IMPLICIT NONE

  INTEGER,                              INTENT (IN)  :: n_radwvl
  REAL (KIND=dp), DIMENSION (n_radwvl), INTENT (IN)  :: curr_rad_wvl
  INTEGER,                              INTENT (OUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER :: imin, imax, i, j, fidx, lidx, stidx, ni0, idx, npts
  REAL (KIND=dp), DIMENSION (max_spec_pts) :: specmod, lowresi0
  REAL (KIND=dp)                           :: frefw, lrefw, scalex
  LOGICAL                                  :: get_lresi0 = .FALSE.

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  CHARACTER (LEN=11), PARAMETER :: modulename = 'dataspline'

  errstat = pge_errstat_ok

  ! ---------------------------------------------------------------------
  ! Load results into the database array. The order of the spectra is
  ! determined by the molecule indices in OMSAO_indices_module. Only 
  ! those spectra that are requested for read-in are interpolated.
  ! Note that this might be different from the actual fitting parameters
  ! to be varied in the fit - a non-zero but constant fitting parameter
  ! will still require a spectrum to make a contribution.
  !
  ! If the original reference spectrum does not cover the wavelength
  ! range of the current radiance wavelength, we interpolate only the
  ! part that is covered and set the rest to Zero.
  ! ------------------ ---------------------------------------------------
 
  ! ------------------------
  ! Spline Reference Spectra
  ! ------------------------
  IF (slit_rad) THEN   ! use derived slit from solar for ring effect
     nslit  = nslit_sol; slitwav = slitwav_sol; slitfit = solslitfit
  END IF
  
  ! 1: solar_idx, obtained latter in prepare_refspecs
  ! 2: ring_idx, done on line
  IF (ring_on_line) THEN
     stidx = ring_idx + 1
  ELSE 
     stidx = ring_idx
  ENDIF
  
  ! Perform solar i0 effect or convole spectra with slit functions
  ! Should be used for spectra not used in VLIDORT calculation or 
  ! when effective cross sections are used
  ni0 = n_refspec_pts(solar_idx)
  DO idx = stidx, max_rs_idx
          
     IF ( n_refspec_pts(idx) >= 3 .AND. INDEX(TRIM(ADJUSTL(refspec_fname(idx))), &
          zerospec_string ) == 0) THEN
        npts = n_refspec_pts(idx)  ! Define short-hand
        specmod(1:npts) = refspec_orig_data(idx,1:npts,spc_idx)

      !  print *, idx, TRIM(ADJUSTL(refspec_fname(idx))), n_radwvl
        !IF (idx == o2o2_idx) refspec_orig_data(idx,1:npts,wvl_idx) = refspec_orig_data(idx,1:npts,wvl_idx)

        IF (idx == bro_idx .OR. idx == bro2_idx .OR. idx == hcho_idx .OR. idx == no2_t1_idx &
             .OR. idx == so2_idx .OR. idx == so2v_idx .OR. idx == o2o2_idx ) THEN
        !IF (idx == bro_idx .OR. idx == hcho_idx .OR. idx == no2_t1_idx) THEN
           IF (idx == bro_idx .OR. idx == bro2_idx) THEN
              scalex = 2.0E13
           ELSE IF (idx == o2o2_idx) THEN
              scalex = 2.6E33
           ELSE
              scalex = 5.0E15
           ENDIF

           scalex = scalex * refspec_norm(idx)
            
           CALL CORRECT_I0EFFECT(refspec_orig_data(idx,1:npts,wvl_idx), &
                specmod(1:npts), npts, refspec_orig_data(solar_idx,1:ni0,wvl_idx), &
                refspec_orig_data(solar_idx,1:ni0,spc_idx), ni0, scalex, get_lresi0, &
                errstat, lowresi0(1:npts))
            
           
           IF (errstat == pge_errstat_error) RETURN
           
        ! ----------------------------------------------------------------------------
        ! Call interpolation and check returned error status. WARNING status indicates
        ! missing parts of the interpolated spectrum, while ERROR status indicates a
        ! more serious condition that requires termination.
        ! ----------------------------------------------------------------------------  
    
        ELSE IF ((idx /= comm_idx .AND. idx /= com1_idx .AND. idx /= ring_idx .AND. idx /= ring1_idx ) .OR. &
             (idx == ring_idx .AND. ring_convol) .OR. (idx == ring1_idx .AND. ring_convol)) THEN
                
           IF (.NOT. yn_varyslit) THEN
              IF (which_slit == 0) THEN
                 CALL gauss_multi (refspec_orig_data(idx,1:npts,wvl_idx), &
                      refspec_orig_data(idx,1:npts,spc_idx), specmod(1:npts), npts)
              ELSE IF (which_slit == 1) THEN
                 CALL asym_gauss_multi (refspec_orig_data(idx,1:npts,wvl_idx), &
                      refspec_orig_data(idx,1:npts,spc_idx), specmod(1:npts), npts)
              ELSE IF (which_slit == 2) THEN
                 CALL asym_voigt_multi (refspec_orig_data(idx,1:npts,wvl_idx), &
                      refspec_orig_data(idx,1:npts,spc_idx), specmod(1:npts), npts)
              ELSE IF (which_slit == 3) THEN
                 CALL triangle_multi (refspec_orig_data(idx,1:npts,wvl_idx),   &
                      refspec_orig_data(idx,1:npts,spc_idx), specmod(1:npts), npts)
              ELSE IF (which_slit == 4) THEN 
                 CALL super_gauss_multi (refspec_orig_data(idx,1:npts,wvl_idx), &
                      refspec_orig_data(idx,1:npts,spc_idx), specmod(1:npts), npts)
              ELSE IF (which_slit == 5) THEN
                SELECT CASE (instrument_idx)
                CASE (omi_idx)
                 CALL omislit_multi (refspec_orig_data(idx,1:npts,wvl_idx),    &
                      refspec_orig_data(idx,1:npts,spc_idx), specmod(1:npts), npts)          
                CASE (omps_idx)
                 !CALL ompsslit_multi (refspec_orig_data(idx,1:npts,wvl_idx),    &
                 !     refspec_orig_data(idx,1:npts,spc_idx), specmod(1:npts), npts)                   
                END SELECT
              END IF
           ELSE 
              IF (which_slit == 0) THEN
                 CALL gauss_vary (refspec_orig_data(idx,1:npts,wvl_idx),   &
                      refspec_orig_data(idx,1:npts,spc_idx),specmod(1:npts), npts)
              ELSE IF (which_slit == 1) THEN
                 CALL asym_gauss_vary (refspec_orig_data(idx,1:npts,wvl_idx),   &
                      refspec_orig_data(idx,1:npts,spc_idx),specmod(1:npts), npts)
              ELSE IF (which_slit == 2) THEN
                 CALL asym_voigt_vary (refspec_orig_data(idx,1:npts,wvl_idx),   &
                      refspec_orig_data(idx,1:npts,spc_idx), specmod(1:npts), npts)
              ELSE IF (which_slit == 3) THEN
                 CALL triangle_vary (refspec_orig_data(idx,1:npts,wvl_idx),     &
                      refspec_orig_data(idx,1:npts,spc_idx), specmod(1:npts), npts)
              ELSE IF (which_slit == 4) THEN
                 CALL super_gauss_vary (refspec_orig_data(idx,1:npts,wvl_idx),   &
                      refspec_orig_data(idx,1:npts,spc_idx),specmod(1:npts), npts)
              ELSE IF (which_slit == 5) THEN
                SELECT CASE (instrument_idx)
                CASE (omi_idx)
                 CALL omislit_vary (refspec_orig_data(idx,1:npts,wvl_idx),      &
                      refspec_orig_data(idx,1:npts,spc_idx), specmod(1:npts), npts)
                CASE (omps_idx)
                 !CALL ompsslit_vary (refspec_orig_data(idx,1:npts,wvl_idx),      &
                 !     refspec_orig_data(idx,1:npts,spc_idx), specmod(1:npts), npts)
                END SELECT
              END IF
          ENDIF
        ENDIF
        
        j = shift_offset + idx

        frefw = refspec_orig_data(idx,1, wvl_idx)
        lrefw = refspec_orig_data(idx,npts, wvl_idx)
        fidx = MINVAL(MINLOC(curr_rad_wvl, MASK=(curr_rad_wvl >= &
             frefw + 0.1 .AND. curr_rad_wvl <= lrefw - 0.1)))
        lidx = MINVAL(MAXLOC(curr_rad_wvl, MASK=(curr_rad_wvl >= &
             frefw + 0.1 .AND. curr_rad_wvl <= lrefw - 0.1)))
        
        IF (lidx > fidx .AND. lidx > 0 .AND. fidx > 0) THEN 
           IF (lo_radbnd(j) < up_radbnd(j)) THEN

              ! xliu: 02/28/2009
              ! Save convolved but at original reference wavelength grid for further interpolation 
              ! Necessary when fitting a wavelength shift 
              refspec_orig_data(idx,1:npts,3) =  specmod(1:npts)

              CALL bspline1(refspec_orig_data(idx,1:npts,wvl_idx), specmod(1:npts),&
                   npts, curr_rad_wvl(fidx:lidx), database(idx, fidx:lidx), &
                   database_shiwf(idx, fidx:lidx), lidx - fidx + 1, errstat)
                
              IF (errstat < 0) THEN
                 WRITE(*, *) modulename, ': BSPLINE1 error, errstat = ', errstat
                 errstat = pge_errstat_error
              ENDIF
           ELSE
              CALL bspline(refspec_orig_data(idx,1:npts,wvl_idx), specmod(1:npts), &
                   npts, curr_rad_wvl(fidx:lidx), database(idx, fidx:lidx), &
                   lidx - fidx + 1, errstat)
              IF (errstat < 0) THEN
                 WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
                 errstat = pge_errstat_error
              ENDIF
           END IF
        ENDIF
     END IF

     IF (slit_rad .AND. i == 2 ) THEN   ! use derived slit for other reference spectra
        nslit  = nslit_rad; slitwav = slitwav_rad; slitfit = radslitfit
     END IF

  END DO

  RETURN
END SUBROUTINE dataspline

SUBROUTINE CORRECT_I0EFFECT(refwav, refspec, nref, i0wav, i0, ni0, scalex, get_lresi0, errstat, lowresi0)
  USE OMSAO_precision_module     
  USE OMSAO_errstat_module
  USE OMSAO_variables_module,   ONLY: yn_varyslit, which_slit, instrument_idx , omi_idx, omps_idx
  USE OMSAO_slitfunction_module
  !USE OMPS_slit_module
  IMPLICIT NONE

  INTEGER,                          INTENT (IN)    :: nref, ni0
  INTEGER,                          INTENT (OUT)   :: errstat
  REAL (KIND=dp),                   INTENT (IN)    :: scalex
  REAL (KIND=dp), DIMENSION (nref), INTENT (IN)    :: refwav
  REAL (KIND=dp), DIMENSION (nref), INTENT (INOUT) :: refspec
  REAL (KIND=dp), DIMENSION (ni0), INTENT  (IN)    :: i0wav, i0
  REAL (KIND=dp), DIMENSION (nref), INTENT(OUT)    :: lowresi0
  LOGICAL, INTENT(IN)                              :: get_lresi0
  
  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER                          :: i, j, fidx, lidx, npts, fidx1, lidx1
  REAL (KIND=dp), DIMENSION (nref) :: specmod, newi0, abspec, abspecmod, abspec1, abspecmod1, frac
  REAL (KIND=dp)                   :: frefw, lrefw
  CHARACTER (LEN=16), PARAMETER    :: modulename = 'CORRECT_I0EFFECT'

  ! T: weight by irradiance (works better for Hartley bands) F: I0 effect (works better for longer wavelengths)
  LOGICAL, PARAMETER               :: weight_irrad = .TRUE.


  errstat = pge_errstat_ok

  frefw = MAXVAL([refwav(1), i0wav(1)])
  lrefw = MINVAL([refwav(nref), i0wav(ni0)])

  fidx = MAXVAL(MINLOC(refwav, MASK = (refwav >= frefw)))
  lidx = MAXVAL(MAXLOC(refwav, MASK = (refwav <= lrefw)))
  npts = lidx - fidx + 1

  ! Interpolate i0 to refwav positions
  CALL bspline(i0wav, i0, ni0, refwav(fidx:lidx), newi0(fidx:lidx), npts, errstat)
  IF (errstat < 0) THEN
     print *, i0wav(1), i0wav(ni0)
     print *, refwav(fidx), refwav(lidx)
     WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
     errstat = pge_errstat_error;      RETURN
  ENDIF

 IF (weight_irrad) THEN
     abspec(fidx:lidx)  = newi0(fidx:lidx) * refspec(fidx:lidx)
     !abspec1(fidx:lidx) = newi0(fidx:lidx) * EXP(-refspec(fidx:lidx) * scalex)
  ELSE
     abspec(fidx:lidx) = newi0(fidx:lidx) * EXP(-refspec(fidx:lidx) * scalex)
  ENDIF

  IF (.NOT. yn_varyslit) THEN
     IF (which_slit == 0) THEN
        CALL gauss_multi (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        !CALL gauss_multi (refwav(fidx:lidx), abspec1(fidx:lidx), abspecmod1(fidx:lidx), npts)
        CALL gauss_multi (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
     ELSE IF (which_slit == 1) THEN
        CALL asym_gauss_multi (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        !CALL asym_gauss_multi (refwav(fidx:lidx), abspec1(fidx:lidx), abspecmod1(fidx:lidx), npts)
        CALL asym_gauss_multi (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
     ELSE IF (which_slit == 2) THEN
        CALL asym_voigt_multi (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        !CALL asym_voigt_multi (refwav(fidx:lidx), abspec1(fidx:lidx), abspecmod1(fidx:lidx), npts)
        CALL asym_voigt_multi (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
     ELSE IF (which_slit == 3) THEN
        CALL triangle_multi   (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        !CALL triangle_multi (refwav(fidx:lidx), abspec1(fidx:lidx), abspecmod1(fidx:lidx), npts)
        CALL triangle_multi   (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
     ELSE IF (which_slit == 4) THEN 
        CALL super_gauss_multi (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        CALL super_gauss_multi (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
     ELSE IF (which_slit == 5) THEN
        SELECT CASE (instrument_idx)
         CASE (omi_idx)
          CALL omislit_multi    (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
          !CALL omislit_multi (refwav(fidx:lidx), abspec1(fidx:lidx), abspecmod1(fidx:lidx), npts)
          CALL omislit_multi    (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)          
         CASE (omps_idx)
         ! CALL ompsslit_multi    (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
         ! CALL ompsslit_multi    (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
         END SELECT
     END IF
  ELSE 
     IF (which_slit == 0) THEN
        CALL gauss_vary (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        !CALL gauss_vary (refwav(fidx:lidx), abspec1(fidx:lidx), abspecmod1(fidx:lidx), npts)
        CALL gauss_vary (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
     ELSE IF (which_slit == 1) THEN
        CALL asym_gauss_vary (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        !CALL asym_gauss_vary (refwav(fidx:lidx), abspec1(fidx:lidx), abspecmod1(fidx:lidx), npts)
        CALL asym_gauss_vary (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
     ELSE IF (which_slit == 2) THEN
        CALL asym_voigt_vary (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        !CALL asym_voigt_vary (refwav(fidx:lidx), abspec1(fidx:lidx), abspecmod1(fidx:lidx), npts)
        CALL asym_voigt_vary (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
     ELSE IF (which_slit == 3) THEN
        CALL triangle_vary   (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        !CALL triangle_vary   (refwav(fidx:lidx), abspec1(fidx:lidx), abspecmod1(fidx:lidx), npts)
        CALL triangle_vary   (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
     ELSE IF (which_slit == 4) THEN 
        CALL super_gauss_vary (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        !CALL gauss_vary (refwav(fidx:lidx), abspec1(fidx:lidx), abspecmod1(fidx:lidx), npts)
        CALL super_gauss_vary (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
     ELSE IF (which_slit == 5) THEN 
       SELECT CASE (instrument_idx) 
       CASE (omi_idx) 
        CALL omislit_vary    (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        CALL omislit_vary    (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
       CASE (omps_idx)
        !CALL ompsslit_vary    (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        !CALL ompsslit_vary    (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
       END SELECT
     ENDIF
  ENDIF

!  ! < 300 nm use abspec, > 310 nm use abspec1
!  ! merge in between
!  fidx1 = MAXVAL(MINLOC(refwav, MASK = (refwav >= 300)))  
!  lidx1 = MAXVAL(MAXLOC(refwav, MASK = (refwav <= 310)))
!  !print *, fidx, fidx1, lidx1, lidx
!  !print *, refwav(fidx), refwav(lidx), refwav(fidx1), refwav(lidx1)
!  IF (fidx1 > fidx) THEN
!     refspec(fidx:fidx1-1) = abspecmod(fidx:fidx1-1) / specmod(fidx:fidx1-1)
!  ENDIF
!  IF (lidx1 < lidx) THEN
!     refspec(lidx1+1:lidx) = -LOG(abspecmod1(lidx1+1:lidx) / specmod(lidx1+1:lidx)) / scalex
!  ENDIF
!  IF (lidx1 > fidx1) THEN
!     frac(fidx1:lidx1) = 1.0 - (refwav(fidx1:lidx1) - 300.) / 10.
!     refspec(fidx1:lidx1) =  abspecmod(fidx1:lidx1) / specmod(fidx1:lidx1) * frac(fidx1:lidx1) - &
!          LOG(abspecmod1(fidx1:lidx1) / specmod(fidx1:lidx1)) / scalex * ( 1.0 - frac(fidx1:lidx1))
!  ENDIF
  
  IF (weight_irrad) THEN
     refspec(fidx:lidx) = abspecmod(fidx:lidx) / specmod(fidx:lidx)
     !refspec(fidx:lidx) = - LOG(abspecmod1(fidx:lidx) / specmod(fidx:lidx)) / scalex
  ELSE
     refspec(fidx:lidx) = - LOG(abspecmod(fidx:lidx) / specmod(fidx:lidx)) / scalex
  ENDIF

  IF ( get_lresi0 ) THEN
     lowresi0(fidx:lidx) = specmod(fidx:lidx)
     lowresi0(1:fidx-1) = 0.0; lowresi0(lidx+1:nref) = 0.0
  ENDIF

  RETURN

END SUBROUTINE CORRECT_I0EFFECT


SUBROUTINE CORRECT_COADDEFFECT(refwav, refspec, i0, nref, choice, nout, errstat)
  USE OMSAO_precision_module     
  USE OMSAO_errstat_module
  IMPLICIT NONE

  INTEGER,                          INTENT (IN)    :: nref, choice
  INTEGER,                          INTENT (OUT)   :: errstat, nout
  REAL (KIND=dp), DIMENSION (nref), INTENT (IN)    :: refwav, i0
  REAL (KIND=dp), DIMENSION (nref), INTENT (INOUT) :: refspec

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER                          :: i, j
  REAL (KIND=dp), DIMENSION (nref) :: newi0, abspec
  REAL (KIND=dp)                   :: frefw, lrefw, scalex
  CHARACTER (LEN=19), PARAMETER    :: modulename = 'CORRECT_COADDEFFECT'

  errstat = pge_errstat_ok
  scalex = 0.08              !  Arbitarily assummed
  newi0 = i0
  
  IF (choice == 1 .OR. choice == 2) THEN 
     abspec = i0 * EXP(-refspec * scalex)
  ENDIF
 
  IF (choice == 1) THEN 
     CALL avg_band_refspec(refwav, abspec, nref, nout, errstat)
     IF ( errstat >= pge_errstat_error ) RETURN
     CALL avg_band_refspec(refwav, newi0, nref, nout, errstat)
     IF ( errstat >= pge_errstat_error ) RETURN
  ELSE IF (choice == 2) THEN     
     CALL avg_band_ozcrs(refwav, abspec, nref, nout, errstat)
     IF ( errstat >= pge_errstat_error ) RETURN
     CALL avg_band_ozcrs(refwav, newi0, nref, nout, errstat)
     IF ( errstat >= pge_errstat_error ) RETURN
  ELSE IF (choice == 3) THEN
     CALL avg_band_refspec(refwav, refspec, nref, nout, errstat)
     IF ( errstat >= pge_errstat_error ) RETURN
  ENDIF
     

  IF (choice == 1 .OR. choice == 2) THEN
     refspec(1:nout) = -LOG(abspec(1:nout) / newi0(1:nout)) / scalex
  ENDIF

  RETURN

END SUBROUTINE CORRECT_COADDEFFECT

SUBROUTINE append_solring(nspec, n1, n2, wav, spec, nsol, solwav, solspec, errstat)
  
  USE OMSAO_precision_module
  USE OMSAO_errstat_module
       
  IMPLICIT NONE  

  ! ================================
  ! Input and Output variables
  ! =================================
  INTEGER, INTENT(IN)                             :: nspec, n1, n2, nsol
  REAL (KIND=dp), DIMENSION(nsol),     INTENT(IN) :: solwav, solspec
  REAL (KIND=dp), DIMENSION(nspec), INTENT(INOUT) :: wav, spec
  INTEGER,                            INTENT(OUT) :: errstat

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp)    :: delw, temp      
  INTEGER           :: i, fidx, lidx, nref

  INTEGER, EXTERNAL :: OMI_SMF_setmsg

  ! --------------------------------
  ! Name of this subroutine/module
  ! --------------------------------
  CHARACTER (LEN=14), PARAMETER :: modulename = 'append_solring'

  errstat = pge_errstat_ok
  
  IF (n1 > 1) THEN
     delw = wav(n1 + 1) - wav(n1)
     DO i = n1-1, 1, -1
        wav(i) = wav(i+1) - delw
     ENDDO
  ENDIF

  IF (n2 < nspec) THEN
     delw = wav(n2) - wav(n2-1)
     DO i = n2 + 1, nspec
        wav(i) = wav(i - 1) + delw
     ENDDO
  ENDIF

  IF (n1 > 1) THEN     
     ! Perform interpolatioon
     fidx = MINVAL(MAXLOC(solwav, MASK=(solwav < wav(1))))
     lidx = MINVAL(MINLOC(solwav, MASK=(solwav > wav(n1))))
     nref = lidx - fidx + 1
     IF (nref <= 4) THEN
        lidx = lidx + 2
        fidx = fidx - 2
     ENDIF 
     IF (fidx < 1 .OR. lidx > nsol) THEN 
        WRITE(*, *) modulename, ': Increase wavelength range of solar reference!!!'
        errstat = pge_errstat_error; RETURN
     ENDIF
     temp = spec(n1)

     CALL interpolation (nref, solwav(fidx:lidx), solspec(fidx:lidx), &
          n1,  wav(1:n1), spec(1:n1), errstat )

     IF ( errstat > pge_errstat_warning ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0); RETURN
     END IF

     spec(1:n1-1) = spec(1:n1-1) * temp / spec(n1)
  ENDIF

  IF (n2 < nspec) THEN  
     ! Perform interpolatioon
     fidx = MINVAL(MAXLOC(solwav, MASK=(solwav < wav(n2))))
     lidx = MINVAL(MINLOC(solwav, MASK=(solwav > wav(nspec))))
     nref = lidx - fidx + 1
     IF (nref <= 4) THEN
        lidx = lidx + 2; fidx = fidx -2
     ENDIF 

     IF (fidx < 1 .OR. lidx > nsol) THEN 
        WRITE(*, *) modulename, 'Increase wavelength range of solar reference!!!'
        errstat = pge_errstat_error; RETURN
     ENDIF

     temp = spec(n2)
     CALL interpolation (nref, solwav(fidx:lidx), solspec(fidx:lidx), &
          nspec-n2+1,  wav(n2:nspec), spec(n2:nspec), errstat )
     IF ( errstat > pge_errstat_warning ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0); RETURN
     END IF
     spec(n2+1:nspec) = spec(n2+1:nspec) * temp / spec(n2)
  ENDIF
     
  RETURN
  
END SUBROUTINE append_solring

!SUBROUTINE CORRECT_I0EFFECT(refwav, refspec, nref, i0wav, i0, ni0, scalex, hw1e, e_asym, errstat)
!  USE OMSAO_precision_module     
!  USE OMSAO_errstat_module
!  IMPLICIT NONE
!
!  INTEGER,                          INTENT (IN)    :: nref, ni0
!  INTEGER,                          INTENT (OUT)   :: errstat
!  REAL (KIND=dp),                   INTENT (IN)    :: scalex, hw1e, e_asym
!  REAL (KIND=dp), DIMENSION (nref), INTENT (IN)    :: refwav
!  REAL (KIND=dp), DIMENSION (nref), INTENT (INOUT) :: refspec
!  REAL (KIND=dp), DIMENSION (ni0), INTENT  (IN)    :: i0wav, i0
!
!  ! ---------------
!  ! Local variables
!  ! ---------------
!  INTEGER                          :: i, j, fidx, lidx, npts
!  REAL (KIND=dp), DIMENSION (nref) :: specmod, newi0, abspec, abspecmod
!  REAL (KIND=dp)                   :: frefw, lrefw
!  CHARACTER (LEN=16), PARAMETER    :: modulename = 'CORRECT_I0EFFECT'
!  LOGICAL                          :: do_exact_solari0corr = .FALSE.
!
!  errstat = pge_errstat_ok
!
!  frefw = MAXVAL([refwav(1), i0wav(1)])
!  lrefw = MINVAL([refwav(nref), i0wav(ni0)])
!
!  fidx = MAXVAL(MINLOC(refwav, MASK = (refwav >= frefw)))
!  lidx = MAXVAL(MAXLOC(refwav, MASK = (refwav <= lrefw)))
!  npts = lidx - fidx + 1
!
!  ! Interpolate i0 to refwav positions
!  CALL interpolation ( ni0, i0wav, i0, npts, refwav(fidx:lidx), newi0(fidx:lidx), errstat)
!
!  IF (errstat < 0) THEN
!     WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
!     errstat = pge_errstat_error;      RETURN
!  ENDIF
!
!  IF (do_exact_solari0corr) THEN
!     ! Solar-I0 effect, has to specify scalex
!     abspec(fidx:lidx) = newi0(fidx:lidx) * EXP(-refspec(fidx:lidx) * scalex)
!  ELSE 
!     ! Simply weight cross sections by solar irradiance
!     ! no need to specify scalex  seems to work slightly better
!     abspec(fidx:lidx)  = newi0(fidx:lidx) * refspec(fidx:lidx)
!  ENDIF
!  
!  CALL asym_gauss (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts, hw1e, e_asym)
!  CALL asym_gauss (refwav(fidx:lidx), newi0(fidx:lidx),    specmod(fidx:lidx), npts, hw1e, e_asym)
!  
!  IF (do_exact_solari0corr) THEN
!     refspec(fidx:lidx) = - LOG(abspecmod(fidx:lidx) / specmod(fidx:lidx)) / scalex
!  ELSE
!     refspec(fidx:lidx) = abspecmod(fidx:lidx) / specmod(fidx:lidx)
!  ENDIF
!
!  RETURN
!
!END SUBROUTINE CORRECT_I0EFFECT


SUBROUTINE normalize_solar_refspec ( n_radwvl, curr_rad_wvl, solar_spec, errstat)

  USE OMSAO_precision_module
  USE OMSAO_indices_module,     ONLY: solar_idx, wvl_idx, spc_idx
  USE OMSAO_variables_module,   ONLY: n_refspec_pts, refspec_orig_data, &
       yn_varyslit, which_slit, refspec_norm, solar_refspec, instrument_idx,omps_idx, omi_idx
  USE OMSAO_parameters_module,  ONLY: max_spec_pts
  USE ozprof_data_module,       ONLY: div_sun
  USE OMSAO_slitfunction_module
  !USE OMPS_slit_module        
  USE OMSAO_errstat_module
  IMPLICIT NONE

  INTEGER,                              INTENT (IN)  :: n_radwvl
  REAL (KIND=dp), DIMENSION (n_radwvl), INTENT (IN)  :: curr_rad_wvl, solar_spec
  INTEGER,                              INTENT (OUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER                                  :: fidx, lidx, ni0, i
  REAL (KIND=dp), DIMENSION (max_spec_pts) :: wave, specmod, ratio
  REAL (KIND=dp), DIMENSION (n_radwvl)     :: solar_spec0, ratio0
  REAL (KIND=dp)                           :: frefw, lrefw

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  CHARACTER (LEN=22), PARAMETER :: modulename = 'normalize_solar_refspec'

  errstat = pge_errstat_ok

  ! Convole high-resolution solar reference spectrum
  ni0  = n_refspec_pts(solar_idx)
  wave = refspec_orig_data(solar_idx,1:ni0,wvl_idx)
  
  IF (.NOT. yn_varyslit) THEN
     IF (which_slit == 0) THEN
        CALL gauss_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
     ELSE IF (which_slit == 1) THEN
        CALL asym_gauss_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
     ELSE IF (which_slit == 2) THEN
        CALL asym_voigt_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
     ELSE IF (which_slit == 3) THEN
        CALL triangle_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
     ELSE IF (which_slit == 4) THEN 
        CALL super_gauss_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
     ELSE IF (which_slit == 5) THEN
       SELECT CASE (instrument_idx)
       CASE(omi_idx)
        CALL omislit_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
       CASE(omps_idx)
       ! CALL ompsslit_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
       END SELECT
     END IF
  ELSE 
     IF (which_slit == 0) THEN
        CALL gauss_vary (wave(1:ni0), solar_refspec(1:ni0),specmod(1:ni0), ni0)
     ELSE IF (which_slit == 1) THEN
        CALL asym_gauss_vary (wave(1:ni0), solar_refspec(1:ni0),specmod(1:ni0), ni0)
     ELSE IF (which_slit == 2) THEN
        CALL asym_voigt_vary (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
     ELSE IF (which_slit == 3) THEN
        CALL triangle_vary (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
     ELSE IF (which_slit == 4) THEN
        CALL super_gauss_vary (wave(1:ni0), solar_refspec(1:ni0),specmod(1:ni0), ni0)
     ELSE IF (which_slit == 5) THEN
      SELECT CASE (instrument_idx)
      CASE (omi_idx)
        CALL omislit_vary (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
      CASE (omps_idx)
       ! CALL ompsslit_vary (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
      END SELECT
     ENDIF
  ENDIF

  CALL bspline(wave(1:ni0), specmod(1:ni0), ni0, curr_rad_wvl, solar_spec0, n_radwvl, errstat)
  ratio0 = solar_spec / solar_spec0 * div_sun / refspec_norm(solar_idx)
  
  frefw = curr_rad_wvl(1)
  lrefw = curr_rad_wvl(n_radwvl)
  fidx = MINVAL(MINLOC(wave(1:ni0), MASK=(wave(1:ni0) >= frefw)))
  lidx = MINVAL(MAXLOC(wave(1:ni0), MASK=(wave(1:ni0) <= lrefw)))
   
  CALL bspline(curr_rad_wvl, ratio0, n_radwvl, wave(fidx:lidx), ratio(fidx:lidx), lidx-fidx + 1, errstat)
  IF (errstat < 0) THEN
     WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
     errstat = pge_errstat_error
  ENDIF
  IF (fidx > 1) ratio(1:fidx-1) = ratio(fidx)
  IF (fidx < ni0) ratio(lidx+1:ni0) = ratio(lidx)

  refspec_orig_data(solar_idx, 1:ni0, spc_idx) = &
       refspec_orig_data(solar_idx, 1:ni0, spc_idx) * ratio(1:ni0)
  
  RETURN
END SUBROUTINE normalize_solar_refspec


SUBROUTINE SOLAR_INTERPOLATION ( n_in, x_in, y_in, n_out, x_out, y_out, errstat )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: max_spec_pts
  USE OMSAO_variables_module,   ONLY: n_refspec_pts, refspec_orig_data, refspec_fname, &
       yn_varyslit, which_slit, refspec_norm, refsol_idx, instrument_idx,omi_idx, omps_idx
  USE OMSAO_indices_module, ONLY: solar_idx, wvl_idx, spc_idx
  USE OMSAO_slitfunction_module
  !USE OMPS_slit_module    
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_warning, pge_errstat_error
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,                           INTENT (IN) :: n_in, n_out
  REAL (KIND=dp), DIMENSION (n_in),  INTENT (IN) :: x_in, y_in
  REAL (KIND=dp), DIMENSION (n_out), INTENT (IN) :: x_out

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER,                           INTENT (INOUT) :: errstat
  REAL (KIND=dp), DIMENSION (n_out), INTENT (OUT)   :: y_out

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER                                  :: ni0, i, errstat1, errstat2
  REAL (KIND=dp), DIMENSION (max_spec_pts) :: i0wav, i0spec, i0specmod
  REAL (KIND=dp), DIMENSION (n_in)         :: i0in
  REAL (KIND=dp), DIMENSION (n_out)        :: i0out
  CHARACTER (LEN=19), PARAMETER    :: modulename = 'SOLAR_INTERPOLATION'

  errstat = pge_errstat_ok

  ni0 = n_refspec_pts(solar_idx)       
  i0wav (1:ni0) = refspec_orig_data(solar_idx, 1:ni0, wvl_idx)
  i0spec(1:ni0) = refspec_orig_data(solar_idx, 1:ni0, spc_idx)

  ! Convolve high resolution solar reference with slit functions
  IF (.NOT. yn_varyslit) THEN
     IF (which_slit == 0) THEN
        CALL gauss_multi (i0wav(1:ni0), i0spec(1:ni0), i0specmod(1:ni0), ni0)
     ELSE IF (which_slit == 1) THEN
        CALL asym_gauss_multi (i0wav(1:ni0), i0spec(1:ni0), i0specmod(1:ni0), ni0)
     ELSE IF (which_slit == 2) THEN
        CALL asym_voigt_multi (i0wav(1:ni0), i0spec(1:ni0), i0specmod(1:ni0), ni0)
     ELSE IF (which_slit == 3) THEN
        CALL triangle_multi   (i0wav(1:ni0), i0spec(1:ni0), i0specmod(1:ni0), ni0)
     ELSE IF (which_slit == 4) THEN 
        CALL super_gauss_multi (i0wav(1:ni0), i0spec(1:ni0), i0specmod(1:ni0), ni0)
     ELSE IF (which_slit == 5) THEN
       SELECT CASE(instrument_idx)
       CASE (omi_idx)
        CALL omislit_multi  (i0wav(1:ni0), i0spec(1:ni0), i0specmod(1:ni0), ni0)
       CASE (omps_idx)
       ! CALL ompsslit_multi  (i0wav(1:ni0), i0spec(1:ni0), i0specmod(1:ni0), ni0)
       END SELECT
     ENDIF
  ELSE 
     IF (which_slit == 0) THEN
        CALL gauss_vary (i0wav(1:ni0), i0spec(1:ni0), i0specmod(1:ni0), ni0)
     ELSE IF (which_slit == 1) THEN
        CALL asym_gauss_vary (i0wav(1:ni0), i0spec(1:ni0), i0specmod(1:ni0), ni0)
     ELSE IF (which_slit == 2) THEN
        CALL asym_voigt_vary (i0wav(1:ni0), i0spec(1:ni0), i0specmod(1:ni0), ni0)
     ELSE IF (which_slit == 3) THEN
        CALL triangle_vary   (i0wav(1:ni0), i0spec(1:ni0), i0specmod(1:ni0), ni0)
     ELSE IF (which_slit == 4) THEN 
        CALL super_gauss_vary (i0wav(1:ni0), i0spec(1:ni0), i0specmod(1:ni0), ni0)
     ELSE IF (which_slit == 5) THEN
        SELECT CASE (instrument_idx)
        CASE (omi_idx)
        CALL omislit_vary (i0wav(1:ni0), i0spec(1:ni0), i0specmod(1:ni0), ni0)
        CASE (omps_idx)
        !CALL ompsslit_vary (i0wav(1:ni0), i0spec(1:ni0), i0specmod(1:ni0), ni0) 
        END SELECT
     ENDIF
  ENDIF

  !WRITE(90, '(F10.4, 2D14.5)') ((i0wav(i), i0spec(i), i0specmod(i)), i=1, ni0)

  CALL interpolation(ni0, i0wav(1:ni0), i0specmod(1:ni0), n_in, x_in, i0in, errstat1)
  CALL interpolation(ni0, i0wav(1:ni0), i0specmod(1:ni0), n_out, x_out, i0out, errstat2)

  y_out = y_in(refsol_idx(1:n_out)) * i0out(1:n_out) / i0in(refsol_idx(1:n_out))

 !WRITE(91, '(F10.4, 2D14.5)') ((x_in(i), y_in(i), i0in(i)), i=1, n_in)
 !WRITE(92, '(F10.4, D14.5, I5)') ((x_out(i), i0out(i), refsol_idx(i)), i=1, n_out)

 errstat = MAX(errstat1, errstat2)
  
  RETURN
  
END SUBROUTINE SOLAR_INTERPOLATION
