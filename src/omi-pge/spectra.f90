  ! *********************** Modification History **************
  ! xliu: 
  ! 1. Replace call to asym_gauss to voigt_gauss
  ! 2. Add variable vgr_idx, vgl_idx, hwr_idx, hwl_idx
  ! **********************************************************

SUBROUTINE spectrum_solar ( &
     npoints, nfitvar, smooth, sol_wav_avg, locwvl, fit, fitvar, doas )

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: vgl_idx, vgr_idx, hwr_idx, hwl_idx, &
       wvl_idx, spc_idx,      &
       max_rs_idx, solar_idx, &
       bl0_idx, bl1_idx, bl2_idx, bl3_idx, sc0_idx, sc1_idx, &
       sc2_idx, sc3_idx, sin_idx,  hwe_idx, asy_idx, shi_idx, squ_idx, &
       hwn_idx, ops_idx
  USE OMSAO_parameters_module, ONLY: max_spec_pts
  USE OMSAO_variables_module,  ONLY: n_refspec_pts, refspec_orig_data, &
       phase, fitwavs, lo_sunbnd, up_sunbnd, fitvar_sol, mask_fitvar_sol,&
       which_slit, fixslitcal, yn_varyslit, currspec,&
       instrument_idx, omps_idx, omi_idx, solwinfit
  USE OMSAO_slitfunction_module
  !USE OMPS_slit_module
  USE OMSAO_errstat_module

  IMPLICIT NONE


  INTEGER,                            INTENT (INOUT) :: npoints, nfitvar
  LOGICAL,                            INTENT (INOUT) :: smooth, doas
  REAL (KIND=dp),                     INTENT (IN)    :: sol_wav_avg
  REAL (KIND=dp), DIMENSION (nfitvar),INTENT (INOUT) :: fitvar
  REAL (KIND=dp), DIMENSION (npoints),INTENT (INOUT) :: locwvl, fit


  REAL (KIND=dp), DIMENSION (npoints) :: del, sunspec_ss, tempspec_ss

  ! =======================================
  ! Variable declarations for IMPLICIT NONE
  ! =======================================
  INTEGER :: i, j, npts, errstat, ref_fidx, ref_lidx

  ! =======================================
  ! Shorthands for solar reference spectrum
  ! =======================================
  REAL (KIND=dp), DIMENSION (n_refspec_pts(solar_idx)) :: kppos, kpspec,&
       kpspec_gauss

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=14), PARAMETER :: modulename = 'spectrum_solar'

  errstat = pge_errstat_ok

  fitvar_sol(mask_fitvar_sol(1:nfitvar)) = fitvar(1:nfitvar)

  npts = n_refspec_pts(solar_idx)
  ref_fidx = MINVAL(MINLOC(refspec_orig_data(solar_idx, 1:npts, wvl_idx), &
       MASK = (refspec_orig_data(solar_idx, 1:npts, wvl_idx) >= fitwavs(1) - 2.0)))
  
  ref_lidx = MINVAL(MAXLOC(refspec_orig_data(solar_idx, 1:npts, wvl_idx), &
       MASK=(refspec_orig_data(solar_idx, 1:npts, wvl_idx)<=fitwavs(npoints)+2.0)))

  npts = ref_lidx - ref_fidx + 1
  kppos (1:npts) = refspec_orig_data(solar_idx, ref_fidx:ref_lidx, wvl_idx)
  kpspec(1:npts) = refspec_orig_data(solar_idx, ref_fidx:ref_lidx, spc_idx)

  !     Spectrum calculation for both fitting and non-fitting cases.

  !     Calculate the spectrum:
  !     First do the shift and squeeze. Shift by FITVAR(SHI_IDX), squeeze by
  !     1 + FITVAR(SQU_IDX); do in absolute sense, to make it easy to back-convert
  !     OMI data.

  kppos(1:npts) = kppos(1:npts) * (1.0 + fitvar_sol(squ_idx)) + fitvar_sol(shi_idx)

  ! =============================================
  ! Broadening and re-sampling of solar spectrum:
  ! =============================================
  IF (which_slit == 5) THEN
     ! solar calibration, don't use variable slit (unknown)
    SELECT CASE (instrument_idx)
    CASE (omi_idx)   
     IF (.NOT. yn_varyslit ) THEN
        CALL omislit_multi (kppos, kpspec, kpspec_gauss, npts)
     ELSE
        CALL omislit_vary  (kppos, kpspec, kpspec_gauss, npts)
     END IF
    CASE (omps_idx)
      IF (.NOT. yn_varyslit .or. fixslitcal) THEN 
       ! CALL ompsslit(kppos, kpspec, kpspec_gauss, npts, fitvar_sol(ops_idx))
      ELSE
       ! CALL ompsslit_vary (kppos, kpspec, kpspec_gauss, npts)
      ENDIF
    END SELECT
  ELSE IF (which_slit == 4) THEN
 
     IF (.not. yn_varyslit .or. fixslitcal) THEN 
        CALL super_gauss (kppos, kpspec, kpspec_gauss,npts,fitvar_sol(hwe_idx), fitvar_sol(hwn_idx))
     ELSE 
        CALL super_gauss_vary (kppos, kpspec, kpspec_gauss, npts)
     ENDIF
  ELSE IF (which_slit == 3) THEN
     IF (.NOT. yn_varyslit .OR. fixslitcal ) THEN
        CALL triangle (kppos, kpspec, kpspec_gauss, npts, &
             fitvar_sol(hwe_idx), fitvar_sol(asy_idx))
     ELSE  
        CALL triangle_vary (kppos, kpspec, kpspec_gauss, npts)
     END IF
  ELSE IF (which_slit == 2) THEN
     IF (.NOT. yn_varyslit .OR. fixslitcal ) THEN
        CALL asym_voigt (kppos, kpspec, kpspec_gauss, npts, &
             fitvar_sol(vgl_idx), fitvar_sol(vgr_idx), &
             fitvar_sol(hwl_idx), fitvar_sol(hwr_idx) )
     ELSE  
        CALL asym_voigt_vary (kppos, kpspec, kpspec_gauss, npts)
     END IF
  ELSE IF (which_slit == 1) THEN
     IF (.NOT. yn_varyslit .OR. fixslitcal ) THEN
        CALL asym_gauss (kppos, kpspec, kpspec_gauss, npts, &
             fitvar_sol(hwe_idx), fitvar_sol(asy_idx))
     ELSE
        CALL asym_gauss_vary (kppos, kpspec, kpspec_gauss, npts)
     END IF
  ELSE
     IF (.NOT. yn_varyslit .OR. fixslitcal ) THEN
        CALL gauss (kppos, kpspec, kpspec_gauss, npts, fitvar_sol(hwe_idx))
     ELSE
        CALL gauss_vary (kppos, kpspec, kpspec_gauss, npts)
     END IF
  ENDIF

  ! ------------------------------------------------------
  ! Re-sample the solar reference spectrum to the OMI grid
  ! ------------------------------------------------------
  CALL interpolation ( &
       npts, kppos(1:npts), kpspec_gauss(1:npts), &
       npoints, locwvl(1:npoints), sunspec_ss(1:npoints), errstat )
  IF ( errstat > pge_errstat_warning ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
  END IF

     
  ! Add up the contributions, with solar intensity as FITVAR_SOL (SIN_IDX),
  ! to include possible linear and Beer's law forms.  Do these as 
  ! linear-Beer's-linear. In order to do DOAS we only need to be careful
  ! to include just linear contributions, since I already high-pass 
  ! filtered them.

  ! -----------
  !  Doing BOAS
  ! -----------
  fit(1:npoints) = fitvar_sol(sin_idx) * sunspec_ss(1:npoints)
  
  ! ----------------
  ! Add the scaling.
  ! ----------------
  del(1:npoints) = locwvl(1:npoints) - sol_wav_avg
 
  fit(1:npoints) = fit(1:npoints) * ( &
       fitvar_sol(sc0_idx)                                               + &
       fitvar_sol(sc1_idx) * del(1:npoints)                              + &
       fitvar_sol(sc2_idx) * del(1:npoints)*del(1:npoints)               + &
       fitvar_sol(sc3_idx) * del(1:npoints)*del(1:npoints)*del(1:npoints) )
 
  ! ------------------------
  ! Add baseline parameters.
  ! ------------------------
  fit(1:npoints) = fit(1:npoints) + &
       fitvar_sol(bl0_idx)                                               + &
       fitvar_sol(bl1_idx) * del(1:npoints)                              + &
       fitvar_sol(bl2_idx) * del(1:npoints)*del(1:npoints)               + &
       fitvar_sol(bl3_idx) * del(1:npoints)*del(1:npoints)*del(1:npoints)

  !WRITE(*, '(6d12.4)') locwvl(1), locwvl(npoints), currspec(1), &
  !       currspec(npoints), fit(1), fit(npoints)
  RETURN
END SUBROUTINE spectrum_solar

SUBROUTINE spectrum_earthshine ( &
     npoints, n_fitvar, smooth, rad_wav_avg, locwvl, fit, fitvar, database, doas )

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: &
       max_rs_idx, max_calfit_idx, solar_idx, ring_idx, ad1_idx, &
       lbe_idx, ad2_idx, mxs_idx, wvl_idx, &
       bl0_idx, bl1_idx, bl2_idx, bl3_idx, sc0_idx, sc1_idx, sc2_idx, &
       sc3_idx, sin_idx, hwe_idx, asy_idx, shi_idx, squ_idx
  USE OMSAO_parameters_module, ONLY: max_spec_pts
  USE OMSAO_variables_module,  ONLY: curr_rad_spec, curr_sol_spec, fitvar_rad, &
       mask_fitvar_rad, n_refwvl, refwvl
  USE OMSAO_errstat_module

  IMPLICIT NONE


  ! ===============
  ! Input variables
  ! ===============
  LOGICAL,                                             INTENT (IN) :: smooth, doas
  INTEGER,                                             INTENT (IN) :: npoints, n_fitvar
  REAL (KIND=dp),                                      INTENT (IN) :: rad_wav_avg
  REAL (KIND=dp), DIMENSION (n_fitvar),                INTENT (IN) :: fitvar
  REAL (KIND=dp), DIMENSION (npoints),                 INTENT (IN) :: locwvl
  REAL (KIND=dp), DIMENSION (max_rs_idx,max_spec_pts), INTENT (IN) :: database

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=dp), DIMENSION (npoints), INTENT (OUT) :: fit

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                              :: i, j, idx, errstat
  REAL (KIND=dp), DIMENSION (npoints)  :: del, sunspec_ss
  REAL (KIND=dp), DIMENSION (n_refwvl) :: sunpos_ss

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=8), PARAMETER :: modulename = 'spectrum'

  !     Spectrum calculation for both fitting and non-fitting cases.

  !     Calculate the spectrum:
  !     First do the shift and squeeze. Shift by FITVAR(SHI_IDX), squeeze by
  !     1 + FITVAR(SQU_IDX); do in absolute sense, to make it easy to back-convert
  !     OMI data.


  errstat = pge_errstat_ok

  ! -----------------------------------------------------------------------------------
  ! First, we have to undo the compression of the FITVAR_RAD array. This compression
  ! is performed in the RADIANCE_FIT subroutine and accelerates the fitting process,
  ! because ELSUNC has to handle less indices. But here we require the original layout,
  ! otherwise the index assingment is screwed.
  ! -----------------------------------------------------------------------------------
  DO i = 1, n_fitvar
     idx = mask_fitvar_rad(i)
     fitvar_rad(idx) = fitvar(i)
  END DO

  sunpos_ss(1:n_refwvl) = refwvl(1:n_refwvl) * (1.0_dp + fitvar_rad(squ_idx)) + &
       fitvar_rad(shi_idx)

  ! Broadening and re-sampling of solar spectrum:
  !       Re-sample the solar reference spectrum to the radiance grid
  
  CALL interpolation ( &
       n_refwvl, sunpos_ss(1:n_refwvl), database(solar_idx,1:n_refwvl), &
       npoints, locwvl(1:npoints), sunspec_ss(1:npoints), errstat )

  IF ( errstat > pge_errstat_warning ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
  END IF

  !     Add up the contributions, with solar intensity as FITVAR_RAD(sin_idx), trace
  !     species beginning at FITVAR_RAD(SQU_IDX+1), to include possible linear and
  !     Beer's law forms.  Do these as linear-Beer's-linear. In order to
  !     do DOAS I only need to be careful to include just linear
  !     contributions, since I already high-pass filtered them.

  ! ==================================================================
  ! For BOAS or any wavelength calibration, we have the following line
  ! ==================================================================
  fit(1:npoints) = fitvar_rad(sin_idx) * sunspec_ss(1:npoints)

  !     DOAS here - the spectrum to be fitted needs to be re-defined:
  IF ( doas ) THEN

     i = max_calfit_idx + (ring_idx-1)*mxs_idx + ad1_idx

     fit(1:npoints) = &
          ! For DOAS, FITVAR_RAD(SIN_IDX) should == 1., and not be varied
          fitvar_rad(sin_idx) * LOG ( sunspec_ss(1:npoints) ) + &
          ! Ring adjustment
          fitvar_rad(i) * (database(ring_idx,3:n_refwvl-2) / sunspec_ss (1:npoints))

     DO j = 1, max_rs_idx
        IF ( j /= solar_idx .AND. j /= ring_idx ) THEN !.AND. database(j,1:npoints) /= 0.0) THEN
           i = max_calfit_idx + (j-1)*mxs_idx + ad1_idx
           fit(1:npoints) = fit(1:npoints) + fitvar_rad(i) * database(j,3:n_refwvl-2)
        END IF
     END DO

  ELSE

     ! Initial add-on contributions.
     DO j = 1, max_rs_idx
        IF ( j /= solar_idx ) THEN !.AND. database(j,1:npoints) /= 0.0_dp ) THEN
           i = max_calfit_idx + (j-1)*mxs_idx + ad1_idx
           fit(1:npoints) = fit(1:npoints) + fitvar_rad(i) * database(j,3:n_refwvl-2)
        END IF
     END DO
     ! Beer's law contributions.
     DO j = 1, max_rs_idx
        IF ( j /= solar_idx ) THEN !.AND. database(j,1:npoints) /= 0.0_dp ) THEN
           i = max_calfit_idx + (j-1)*mxs_idx + lbe_idx
           fit(1:npoints) = fit(1:npoints) * EXP(-fitvar_rad(i)*database(j,3:n_refwvl-2))
           
        END IF
     END DO
     ! Final add-on contributions.
     DO j = 1, max_rs_idx
        IF ( j /= solar_idx ) THEN !.AND. database(j,1:npoints) /= 0.0_dp ) THEN
           i = max_calfit_idx + (j-1)*mxs_idx + ad2_idx
           fit(1:npoints) = fit(1:npoints) + fitvar_rad(i) * database(j,3:n_refwvl-2)
        END IF
     END DO
  END IF

  ! Add the scaling.
  del(1:npoints) = locwvl(1:npoints) - rad_wav_avg
  fit(1:npoints) = fit(1:npoints) * ( &
       fitvar_rad(sc0_idx)                                               + &
       fitvar_rad(sc1_idx) * del(1:npoints)                              + &
       fitvar_rad(sc2_idx) * del(1:npoints)*del(1:npoints)               + &
       fitvar_rad(sc3_idx) * del(1:npoints)*del(1:npoints)*del(1:npoints) )
  
  ! Add baseline parameters.
  fit(1:npoints) = fit(1:npoints) + &
       fitvar_rad(bl0_idx)                                               + &
       fitvar_rad(bl1_idx) * del(1:npoints)                              + &
       fitvar_rad(bl2_idx) * del(1:npoints)*del(1:npoints)               + &
       fitvar_rad(bl3_idx) * del(1:npoints)*del(1:npoints)*del(1:npoints)

  RETURN
END SUBROUTINE spectrum_earthshine


