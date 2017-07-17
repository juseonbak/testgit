  ! *********************** Modification History ************
  ! xliu: 
  ! 1. Add subroutine radiance_wavcal_vary
  ! 2. Print fitting variables at the end of radiance_wavcal
  ! 3. Add vgl, vgr, hwl, hwr, and corresponding indices
  ! ********************************************************

SUBROUTINE radiance_wavcal ( n_rad_wvl, curr_rad_spec, error)

  USE OMSAO_precision_module
  USE OMSAO_indices_module,     ONLY: max_calfit_idx, shi_idx, squ_idx, &
       wvl_idx, spc_idx, sig_idx,  hwe_idx, hwr_idx, asy_idx, vgr_idx,  &
       vgl_idx, hwl_idx, bl0_idx, bl1_idx, bl2_idx, bl3_idx, sc0_idx,   &
       sc1_idx, sc2_idx, sc3_idx
  USE OMSAO_variables_module,   ONLY: n_fitvar_sol, fitvar_sol, &
       fitvar_sol_saved,  mask_fitvar_sol, lo_sunbnd, up_sunbnd, wincal_wav, &
       winlim, radwinfit, solwinfit, fixslitcal, fitwavs, fitweights, &
       currspec, which_slit, sol_wav_avg, fitvar_sol_init, numwin, nradpix, scnwrt
       
  IMPLICIT NONE

  ! ===================
  ! IN/Output variables
  ! ===================
  INTEGER, INTENT (IN)                                        :: n_rad_wvl
  REAL(KIND=dp), DIMENSION(sig_idx, n_rad_wvl), INTENT(INOUT) :: curr_rad_spec
  LOGICAL, INTENT (OUT)                                       :: error

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION (n_rad_wvl)        :: allwaves, del
  REAL (KIND=dp), DIMENSION(max_calfit_idx, 2) :: tmp_fitvar
  REAL (KIND=dp)                               :: tmpwave             
  INTEGER :: i, iwin, fidx, lidx, n_fit_pts, slit_unit, solfit_exval
  LOGICAL, SAVE :: wrt_to_screen, wrt_to_file, slitcal, first = .TRUE.   

  IF (first) THEN
     fixslitcal = .FALSE.   
     slitcal = .FALSE.;  wrt_to_file = .FALSE.; slit_unit = 1000
     
     IF (scnwrt) THEN
        wrt_to_screen = .TRUE.
     ELSE
        wrt_to_screen = .FALSE.
     ENDIF
     
     ! fix variables for hwe, asy, vgr, vgl, hwr, hwl
     fitvar_sol(hwe_idx:asy_idx) = 0_dp
     lo_sunbnd(hwe_idx:asy_idx)  = 0_dp
     up_sunbnd(hwe_idx:asy_idx)  = 0_dp
     
     fitvar_sol(vgl_idx:hwr_idx) = 0_dp
     lo_sunbnd(vgl_idx:hwr_idx)  = 0_dp
     up_sunbnd(vgl_idx:hwr_idx)  = 0_dp
     
     ! find the locations of actually used fitting variables
     n_fitvar_sol = 0
     DO i = 1, max_calfit_idx
        IF (lo_sunbnd(i) < up_sunbnd(i) ) THEN
           n_fitvar_sol =  n_fitvar_sol + 1
           mask_fitvar_sol(n_fitvar_sol) = i
        END IF
     END DO
     first = .FALSE.
  ENDIF
  
  error = .FALSE.
  allwaves = curr_rad_spec(wvl_idx, 1:n_rad_wvl)
  fidx = 1
  DO iwin = 1, numwin
     
     n_fit_pts = nradpix(iwin)
     lidx = fidx + n_fit_pts - 1

     IF (scnwrt) WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', i, winlim(iwin,1), &
          winlim(iwin,2), nradpix(iwin)
     
     fitwavs   (1:n_fit_pts) = curr_rad_spec(wvl_idx, fidx:lidx)
     currspec  (1:n_fit_pts) = curr_rad_spec(spc_idx, fidx:lidx)
     fitweights(1:n_fit_pts) = curr_rad_spec(sig_idx, fidx:lidx)
     
     fitvar_sol = fitvar_sol_init
     fitvar_sol(hwe_idx:hwr_idx) = solwinfit(iwin, hwe_idx:hwr_idx, 1)
     CALL cal_fit_one (n_fit_pts, n_fitvar_sol, wrt_to_screen, wrt_to_file, &
          slitcal, slit_unit, wincal_wav(iwin), &
          radwinfit(iwin, 1:max_calfit_idx, 1:2), solfit_exval)

     IF (solfit_exval < 0) THEN
        WRITE(*, *) 'Radiance_wavcal: not converge for window: ', iwin
        error = .TRUE.
        fitvar_sol(shi_idx) = 0.0; fitvar_sol(squ_idx) = 0.0
     ENDIF
            
     ! =================================
     ! Shift and squeeze solar spectrum.
     ! =================================
     ! fitvar_sol is updated in solar_fit_one through common module variables
     curr_rad_spec(wvl_idx,fidx:lidx) = (curr_rad_spec(wvl_idx,fidx:lidx) - &
          fitvar_sol(shi_idx)) / (1.0 + fitvar_sol(squ_idx))

     fidx = lidx + 1
     
  END DO
  
  RETURN
END SUBROUTINE radiance_wavcal
