SUBROUTINE solar_fit_vary (slit_unit, error )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : maxchlen
  USE OMSAO_indices_module,    ONLY : max_calfit_idx, shi_idx, squ_idx, wvl_idx,&
       spc_idx, sig_idx, hwe_idx, hwr_idx, vgr_idx, vgl_idx, asy_idx, hwl_idx
  USE OMSAO_variables_module,  ONLY : curr_sol_spec, n_fitvar_sol, fitvar_sol,  &
       mask_fitvar_sol, lo_sunbnd, up_sunbnd, n_irrad_wvl, slitwav_sol,         &
       nslit_sol, slit_fit_pts, n_slit_step, smooth_slit,slit_redo, slit_fname, &
       nsolpix, winlim, poly_order, solslitfit, which_slit, wavcal_sol,         &
       fixslitcal, fitweights, fitwavs, currspec, fitvar_sol_saved, numwin,     &
       nslit, slitwav, slitfit, sol_spec_ring, nsol_ring, currpixchar, &
       fitvar_sol_init,lo_sunbnd_init, up_sunbnd_init, scnwrt, wavcal
  USE OMSAO_errstat_module
       

  IMPLICIT NONE

  ! ================
  ! Input/OUTPUT variables
  ! ================
  INTEGER, INTENT(IN)    :: slit_unit
  LOGICAL, INTENT (OUT)  :: error

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp)                          :: tmpwave
  REAL (KIND=dp), DIMENSION (n_irrad_wvl) :: allwaves, locshi, locsqu, locspec
  REAL (KIND=dp), DIMENSION (max_calfit_idx, 2) :: tmp_fitvar
  INTEGER :: npoints, i, j, iwin, fidx, lidx, errstat = pge_errstat_ok, nstep,   &
       islit, fpos, lpos, n_fit_pts, fslit, lslit, ios, finter, linter,          &
       npoly, solfit_exval, slidx, sfidx
  CHARACTER(LEN=maxchlen)                       :: tmpchar, fname
  LOGICAL :: wrt_to_screen, wrt_to_file, slitcal,  calfname_exist = .TRUE.

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=14), PARAMETER :: modulename = 'solar_fit_vary'

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  fixslitcal = .TRUE.; error = .FALSE.
  wrt_to_screen = .FALSE.; wrt_to_file = .TRUE.; slitcal=.TRUE.
  allwaves = curr_sol_spec(wvl_idx, 1:n_irrad_wvl)

  ! Determine if file exists or not
  fname = TRIM(ADJUSTL(slit_fname)) // currpixchar // '.dat'
  INQUIRE (FILE=TRIM(ADJUSTL(fname)), EXIST=calfname_exist)

  ! Open a file for saving/reading the results
  OPEN (UNIT=slit_unit, FILE=TRIM(ADJUSTL(fname)), STATUS='UNKNOWN', IOSTAT=errstat)
  IF (errstat /= pge_errstat_ok) THEN
     WRITE(*, *) modulename, ' : cannot open slit_fname', TRIM(ADJUSTL(fname)) 
     error = .TRUE.; RETURN
  END IF
     
  IF (slit_redo .OR. .NOT. calfname_exist) THEN
     IF (which_slit == 2) then
        WRITE(slit_unit, '(A)') &
        ' wave    vgl     vgr     hwl    hwr   shift    squeeze     sin    rms   exval'
     ELSE
        WRITE(slit_unit, '(A)') &
        ' wave    hwe     asym  shift    squeeze      sin      rms   exval '
     END IF

     lo_sunbnd  = lo_sunbnd_init
     up_sunbnd  = up_sunbnd_init
     fitvar_sol_saved = fitvar_sol_init 
     n_fitvar_sol = 0
     DO i = 1, max_calfit_idx
        IF (lo_sunbnd(i) < up_sunbnd(i) ) THEN
           n_fitvar_sol =  n_fitvar_sol + 1
           mask_fitvar_sol(n_fitvar_sol) = i
        END IF
     END DO
     
     islit  = 0              ! number of sucessful calibrations      
     fidx = 1                ! first pixel
     DO iwin = 1, numwin    

        IF (scnwrt) WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', iwin, winlim(iwin,1), &
             winlim(iwin,2), nsolpix(iwin)

        lidx = fidx + nsolpix(iwin) - 1   ! ending index in this window

        ! do each fit using points (fpos:lpos)
        fpos = fidx
        
        DO WHILE (fpos < lidx)

           DO i = fpos+1, fpos + slit_fit_pts - 1
              ! until end of window or sudden gap of > 2.0 um
              IF (i > lidx) EXIT
              IF (allwaves(i) - allwaves(i-1) > 2.0) EXIT
           ENDDO
           lpos = i - 1; npoints = lpos - fpos + 1

           IF (npoints < slit_fit_pts / 2) THEN
              ! Either until the end with not enough points or gap behind
              fpos = lpos + 1; CYCLE
           ENDIF
           
           fitwavs   (1:npoints) = curr_sol_spec(wvl_idx, fpos:lpos)
           currspec  (1:npoints) = curr_sol_spec(spc_idx, fpos:lpos)
           fitweights(1:npoints) = curr_sol_spec(sig_idx, fpos:lpos)
           fitvar_sol =  fitvar_sol_saved
           CALL cal_fit_one (npoints, n_fitvar_sol, wrt_to_screen, wrt_to_file,&
                slitcal, slit_unit, tmpwave, tmp_fitvar, solfit_exval)
           
           IF (solfit_exval > 0) THEN
              islit = islit + 1
              slitwav_sol(islit) = tmpwave
              solslitfit(islit, 1:max_calfit_idx, 1:2) = &
                   tmp_fitvar(1:max_calfit_idx, 1:2)
           END IF

           fpos = fpos + n_slit_step          
        END DO
        
        fidx = lidx + 1       
     END DO
     nslit_sol = islit 
     
  ELSE  
     READ (slit_unit, '(A)') tmpchar
     islit = 0

     DO 
        islit = islit + 1
        IF (which_slit == 2) THEN
           READ (slit_unit, *, IOSTAT = ios) slitwav_sol(islit), &
                solslitfit(islit, vgl_idx, 1), solslitfit(islit, vgr_idx, 1), &
                solslitfit(islit, hwl_idx, 1), solslitfit(islit, hwr_idx, 1), &
                solslitfit(islit, shi_idx, 1), solslitfit(islit, squ_idx, 1)
        ELSE
           READ (slit_unit, *, IOSTAT = ios) slitwav_sol(islit), &
                solslitfit(islit, hwe_idx, 1), solslitfit(islit, asy_idx, 1), &
                solslitfit(islit, shi_idx, 1), solslitfit(islit, squ_idx, 1)
        END IF
        
        IF (ios < 0) THEN
           islit = islit - 1; nslit_sol = islit; EXIT  ! end of file
        ELSE IF (ios > 0) THEN
           WRITE(*, *) modulename, ': read slit file error'
           error = .TRUE. ; RETURN 
        END IF
     END DO

  END IF
  CLOSE(slit_unit)

  ! Disable smoothing slit funciton (only smooth shift)
  IF (smooth_slit ) THEN

     fidx = 1
     DO iwin = 1, numwin
     
        lidx = fidx + nsolpix(iwin) - 1     
        fslit = MINVAL(MINLOC(slitwav_sol(1:nslit_sol), &
             MASK=(slitwav_sol(1:nslit_sol) >= winlim(iwin, 1))))
        lslit = MINVAL(MAXLOC(slitwav_sol(1:nslit_sol), &
             MASK=(slitwav_sol(1:nslit_sol) <= winlim(iwin, 2))))
     
        IF (lslit > fslit + 3) THEN           ! use a 4-order polynomial
           poly_order = 3
  
           DO i = hwe_idx, hwr_idx
              IF (solslitfit(fslit,i,1)==0. .AND. solslitfit(lslit, i, 1)==0.0) CYCLE
              IF (wavcal_sol .AND. (i == shi_idx .OR. i == squ_idx)) CYCLE
              
              npoly = lslit - fslit + 1
              locspec(1:npoly) = solslitfit(fslit:lslit, i, 1)
              CALL subtract_poly_meas (slitwav_sol(fslit:lslit), npoly, &
                   locspec(1:npoly), 1, npoly)
              solslitfit(fslit:lslit, i,1) = solslitfit(fslit:lslit, i, 1)&
                   - locspec(1:npoly) 

              !DO j = 3, npoly - 2
              !  locspec(j) = (solslitfit(fslit+j-1, i, 1)*6.0 + solslitfit
              !  (fslit+j-2, i, 1)*4.0 + solslitfit(fslit+j-3, i, 1) + 
              !  solslitfit(fslit+j, i, 1)*4.0 + solslitfit(fslit+j+1, i, 1)) / 15.0
              !END DO
              !solslitfit(fslit:lslit, i, 1) =  locspec(1:npoly)
          END DO             
        ELSE                                 ! use average
           DO i = hwe_idx, hwr_idx
              solslitfit(fslit:lslit, i, 1) = SUM(solslitfit(fslit:lslit, i, 1)) / &
                   (lslit - fslit + 1.0)
           END DO
        ENDIF
     
        fidx = lidx + 1
     END DO
  END IF

  ! to be used in voigt.f90 or gauss.f90
  nslit = nslit_sol; slitwav = slitwav_sol; slitfit = solslitfit

  OPEN (unit=77, file='slit_after_smooth.dat')
  WRITE(77, *) nslit_sol
  DO i = 1, nslit_sol
     WRITE(77, '(3d14.6)') slitwav_sol(i), solslitfit(i, hwe_idx, 1), &
          solslitfit(i, shi_idx, 1)
  END DO
  CLOSE (77)
 
  ! squeeze and shift wavelength position
  IF (wavcal_sol .OR. .NOT. wavcal) RETURN     ! done in next iteration
  
  fslit = 1; fidx = 1; 
   DO iwin = 1, numwin
     
     lidx = fidx + nsolpix(iwin) - 1     
     lslit = MINVAL(MAXLOC(slitwav_sol, MASK=(slitwav_sol(1:nslit_sol)&
          <= winlim(iwin, 2))))

     IF (lslit < fslit + 3) THEN
        locshi(fidx:lidx) = solslitfit(fslit, shi_idx, 1)
        locsqu(fidx:lidx) = solslitfit(fslit, squ_idx, 1)
     ELSE
        finter = MINVAL(MINLOC(allwaves, MASK=(allwaves > slitwav_sol(fslit))))
        linter = MINVAL(MAXLOC(allwaves, MASK=(allwaves < slitwav_sol(lslit)))) 
        IF (finter == 0) CYCLE

        CALL interpolation (lslit-fslit+1,  slitwav_sol(fslit:lslit), &
             solslitfit(fslit:lslit, shi_idx, 1),   linter-finter+1, &
             allwaves(finter:linter), locshi(finter:linter),errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
     
        CALL interpolation (lslit-fslit+1, slitwav_sol(fslit:lslit),  &
             solslitfit(fslit:lslit, squ_idx, 1),  linter-finter+1,  &
             allwaves(finter:linter), locsqu(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
        
        IF (finter > fidx ) THEN
           locshi(fidx:finter-1)=solslitfit(fslit, shi_idx, 1)
           locsqu(fidx:finter-1)=solslitfit(fslit, squ_idx, 1)
        END IF
        
        IF (linter < lidx)  THEN
           locshi(linter+1:lidx)  = solslitfit(lslit, shi_idx, 1)
           locsqu(linter+1:lidx) =  solslitfit(lslit, squ_idx, 1)
        END IF
     END IF

     allwaves(fidx:lidx) = (allwaves(fidx:lidx) - locshi(fidx:lidx) ) &
          / ( 1.0 + locsqu(fidx:lidx))
     fidx = lidx + 1; fslit = lslit + 1     
  END DO

  fidx = 1; sfidx  =1
  DO i = 1, numwin
     lidx  = fidx + nsolpix(i) - 1
     IF (i == numwin) THEN 
        slidx = nsol_ring
     ELSE
        slidx = MINVAL(MAXLOC(sol_spec_ring(1, 1:nsol_ring), &
             MASK=(sol_spec_ring(1, 1:nsol_ring) < curr_sol_spec(wvl_idx, lidx+1))))
     ENDIF

     finter = MINVAL(MAXLOC(sol_spec_ring(1, 1:nsol_ring), &
             MASK=(sol_spec_ring(1, 1:nsol_ring) == curr_sol_spec(wvl_idx, fidx))))
     linter = finter + nsolpix(i) - 1

     sol_spec_ring(1, finter:linter) = allwaves(fidx:lidx)
     IF (finter > sfidx) sol_spec_ring(1, sfidx:finter-1) = &
          (sol_spec_ring(1,sfidx:finter-1) - locshi(fidx)) / (1.0 + locsqu(fidx))
     IF (linter < slidx) sol_spec_ring(1, linter+1:slidx) = &
          (sol_spec_ring(1, linter+1:slidx) - locshi(lidx)) / (1.0 + locsqu(lidx))
    
     fidx = lidx+ 1
     sfidx= slidx + 1
  ENDDO

  curr_sol_spec(wvl_idx, 1:n_irrad_wvl) = allwaves
   
  RETURN
END SUBROUTINE solar_fit_vary
