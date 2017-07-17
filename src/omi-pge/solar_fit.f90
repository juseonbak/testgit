  ! *********************** Modification History ********
  ! xliu: 
  ! 1. Add subroutine solar_fit_vary
  ! 2. Print fitting variables at the end of solar_fit
  ! 3. Add indices for voigt function
  ! *****************************************************

SUBROUTINE solar_fit (error)

  USE OMSAO_precision_module
  USE OMSAO_indices_module,     ONLY: max_calfit_idx, shi_idx, squ_idx,&
       wvl_idx, spc_idx, sig_idx, hwr_idx, hwe_idx, vgl_idx, asy_idx, ops_idx,hwn_idx
  USE OMSAO_variables_module,   ONLY: n_irrad_wvl, curr_sol_spec, nsolpix, sol_spec_ring, nsol_ring, &
       currpix,numwin, fitwavs,fitweights, currspec, &
       n_fitvar_sol, fitvar_sol, fitvar_sol_saved,  fitvar_sol_init, &
       mask_fitvar_sol, lo_sunbnd, up_sunbnd,lo_sunbnd_init, up_sunbnd_init, &
       which_slit, instrument_idx, omps_idx, fixslitcal,  &
       wavcal, wavcal_sol, slit_fname, write_is_slit, here_stop, calfitspec,&
       wincal_wav, solwinfit, sol_wav_avg, slit_ins_idx
       
  IMPLICIT NONE

  ! ================
  ! Output variables
  ! ================
  LOGICAL,            INTENT (OUT)             :: error

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION (n_irrad_wvl)      :: allwaves, del
  REAL (KIND=dp), DIMENSION(max_calfit_idx, 2) :: tmp_fitvar
  REAL (KIND=dp)                               :: tmpwave             
  INTEGER         :: i, iwin, fidx, lidx, n_fit_pts,  &
       solfit_exval, sfidx, slidx, finter, linter
  INTEGER, SAVE   :: slit_unit
  LOGICAL, SAVE   :: wrt_to_screen, wrt_to_file, slitcal
  LOGICAL, SAVE   :: first = .TRUE.

  IF (first) THEN
     wrt_to_screen = .false.
     fixslitcal = .TRUE.; slitcal = .TRUE.
     wrt_to_file = .false. ; slit_unit = 1000
     fitvar_sol = fitvar_sol_init
     lo_sunbnd  = lo_sunbnd_init
     up_sunbnd  = up_sunbnd_init
     ! find the locations of actually used fitting variables
     IF (which_slit == slit_ins_idx) THEN
        fitvar_sol(hwe_idx:asy_idx) = 0_dp
        lo_sunbnd(hwe_idx:asy_idx)  = 0_dp
        up_sunbnd(hwe_idx:asy_idx)  = 0_dp
        
        fitvar_sol(vgl_idx:hwn_idx) = 0_dp
        lo_sunbnd(vgl_idx:hwn_idx)  = 0_dp
        up_sunbnd(vgl_idx:hwn_idx)  = 0_dp
        fitvar_sol(ops_idx) = 1_dp
        slitcal = .FALSE.
        IF (instrument_idx == omps_idx) slitcal =.true.
     ELSE  
        fitvar_sol(ops_idx) = 0_dp
        lo_sunbnd(ops_idx) = 0_dp 
        up_sunbnd(ops_idx) = 0_dp
     ENDIF
      
     n_fitvar_sol = 0
     DO i = 1, max_calfit_idx
        IF (lo_sunbnd(i) < up_sunbnd(i) ) THEN
           n_fitvar_sol =  n_fitvar_sol + 1
           mask_fitvar_sol(n_fitvar_sol) = i            
        END IF
     END DO
  ENDIF

 
  error = .FALSE.
  allwaves = curr_sol_spec(wvl_idx, 1:n_irrad_wvl)
   
  fidx = 1
  
  DO iwin = 1, numwin    

     n_fit_pts = nsolpix(iwin)
     lidx = fidx + n_fit_pts - 1
     fitwavs   (1:n_fit_pts) = curr_sol_spec(wvl_idx, fidx:lidx)
     currspec  (1:n_fit_pts) = curr_sol_spec(spc_idx, fidx:lidx)
     fitweights(1:n_fit_pts) = curr_sol_spec(sig_idx, fidx:lidx)
       

     IF (wrt_to_screen) WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', iwin, fitwavs(1), &
          fitwavs(n_fit_pts), nsolpix(iwin)
     

     calfitspec(iwin, 1, 1:n_fit_pts) = fitwavs   (1:n_fit_pts)
     calfitspec(iwin, 2, 1:n_fit_pts) = currspec  (1:n_fit_pts)
      
     
     CALL cal_fit_one (n_fit_pts, n_fitvar_sol, wrt_to_screen, wrt_to_file,&
          slitcal, slit_unit, wincal_wav(iwin), &
          solwinfit(iwin,1:max_calfit_idx, 1:2), solfit_exval)
     calfitspec(iwin, 3, 1:n_fit_pts) = currspec(1:n_fit_pts)
     IF (solfit_exval < 0) THEN
        WRITE(*, *) 'Solar_fit: solar calibration not converge for window: ', iwin
        error = .TRUE.; RETURN
     END IF
              
     ! =================================
     ! Shift and squeeze solar spectrum.
     ! =================================
     ! fitvar_sol is updated in solar_fit_one through common module variables
    ! IF (wavcal .or. wavcal_sol) allwaves(fidx:lidx) = (allwaves(fidx:lidx) - fitvar_sol(shi_idx)) / (1.0 + fitvar_sol(squ_idx))  

     fidx = lidx + 1
  END DO
  curr_sol_spec(wvl_idx, 1:n_irrad_wvl) = allwaves

  RETURN

END SUBROUTINE solar_fit

SUBROUTINE solar_wavcal (error)

  USE OMSAO_precision_module
  USE OMSAO_indices_module,     ONLY: max_calfit_idx, shi_idx, squ_idx,&
       wvl_idx, spc_idx, sig_idx, hwr_idx, hwe_idx, vgl_idx, asy_idx,ops_idx,hwn_idx
  USE OMSAO_variables_module,   ONLY:    n_irrad_wvl, curr_sol_spec, nsol_ring,sol_spec_ring, &
       numwin,nsolpix, fitwavs,fitweights, currspec,  &
       n_fitvar_sol, fitvar_sol, fitvar_sol_saved,  fitvar_sol_init, &
       mask_fitvar_sol, lo_sunbnd, up_sunbnd,lo_sunbnd_init, up_sunbnd_init, &
       fixslitcal,   calfitspec, wincal_wav,  solwinfit

  IMPLICIT NONE

  ! ================
  ! Output variables
  ! ================
  LOGICAL,            INTENT (OUT)             :: error

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION (n_irrad_wvl)      :: allwaves, del
  REAL (KIND=dp), DIMENSION(max_calfit_idx, 2) :: tmp_fitvar
  REAL (KIND=dp)                               :: tmpwave
  INTEGER         :: i, iwin, fidx, lidx, n_fit_pts,  &
       solfit_exval, sfidx, slidx, finter, linter
  INTEGER, SAVE   :: slit_unit
  LOGICAL, SAVE   :: wrt_to_screen, wrt_to_file, slitcal
  LOGICAL, SAVE   :: first = .TRUE.
   
  IF (first) THEN
      wrt_to_screen = .TRUE.
      fixslitcal = .TRUE.  ; slitcal = .false.
      wrt_to_file = .false. ; slit_unit = 1000
      fitvar_sol = fitvar_sol_init
      lo_sunbnd  = lo_sunbnd_init; up_sunbnd  = up_sunbnd_init

      !fitvar_sol(hwe_idx:asy_idx) = 0_dp
      lo_sunbnd(hwe_idx:asy_idx)  = 0_dp
      up_sunbnd(hwe_idx:asy_idx)  = 0_dp

      !fitvar_sol(vgl_idx:ops_idx) = 0_dp
      lo_sunbnd(vgl_idx:ops_idx)  = 0_dp
      up_sunbnd(vgl_idx:ops_idx)  = 0_dp

      n_fitvar_sol = 0
     DO i = 1, max_calfit_idx
        IF (lo_sunbnd(i) < up_sunbnd(i) ) THEN
           n_fitvar_sol =  n_fitvar_sol + 1
           mask_fitvar_sol(n_fitvar_sol) = i
        END IF
     END DO
  ENDIF
  error = .FALSE.
  allwaves = curr_sol_spec(wvl_idx, 1:n_irrad_wvl)
  fidx = 1

  DO iwin = 1, numwin

     n_fit_pts = nsolpix(iwin)
     lidx = fidx + n_fit_pts - 1
     fitwavs   (1:n_fit_pts) = curr_sol_spec(wvl_idx, fidx:lidx)
     currspec  (1:n_fit_pts) = curr_sol_spec(spc_idx, fidx:lidx)
     fitweights(1:n_fit_pts) = curr_sol_spec(sig_idx, fidx:lidx)
     fitvar_sol(1:max_calfit_idx) = solwinfit(iwin,1:max_calfit_idx, 1)

     IF (wrt_to_screen) WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', iwin, fitwavs(1),&
          fitwavs(n_fit_pts), nsolpix(iwin)

     calfitspec(iwin, 1, 1:n_fit_pts) = fitwavs   (1:n_fit_pts)
     calfitspec(iwin, 2, 1:n_fit_pts) = currspec  (1:n_fit_pts)

     CALL cal_fit_one (n_fit_pts, n_fitvar_sol, wrt_to_screen, wrt_to_file,&
          slitcal, slit_unit, wincal_wav(iwin), &
          solwinfit(iwin,1:max_calfit_idx, 1:2), solfit_exval)

     calfitspec(iwin, 3, 1:n_fit_pts) = currspec(1:n_fit_pts)
     IF (solfit_exval < 0) THEN
        WRITE(*, *) 'Solar_fit: solar calibration not converge for window: ',iwin
        error = .TRUE.; RETURN
     END IF                                 
     allwaves(fidx:lidx) = (allwaves(fidx:lidx) - fitvar_sol(shi_idx))/(1.0 +fitvar_sol (squ_idx))
     fidx = lidx + 1
   END DO
   curr_sol_spec(wvl_idx, 1::n_irrad_wvl) = allwaves
   RETURN
END SUBROUTINE solar_wavcal
 
