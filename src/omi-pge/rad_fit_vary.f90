SUBROUTINE rad_fit_vary (rslit_unit, n_rad_wvl, curr_rad_spec, error )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : maxchlen
  USE OMSAO_indices_module,    ONLY : max_calfit_idx, shi_idx, squ_idx,&
       wvl_idx, spc_idx, sig_idx, hwe_idx, hwr_idx, vgr_idx, vgl_idx,  &
       asy_idx, hwl_idx
  USE OMSAO_variables_module,  ONLY : n_fitvar_sol, fitvar_sol,        &
       mask_fitvar_sol, lo_sunbnd, up_sunbnd, slitwav_rad, nslit_rad,  &
       slit_fit_pts, n_slit_step, smooth_slit,  slit_redo, rslit_fname,&
       nradpix, winlim, poly_order, radslitfit,which_slit, wavcal_sol, &
       fixslitcal,fitweights, fitwavs, currspec, fitvar_sol_saved, &
       fitvar_sol_init,lo_sunbnd_init, up_sunbnd_init, numwin, &
       nslit, slitwav, slitfit, currpixchar, scnwrt, wavcal
  USE OMSAO_errstat_module
       

  IMPLICIT NONE

  ! ================
  ! Input/OUTPUT variables
  ! ================
  INTEGER, INTENT(IN)              :: rslit_unit
  LOGICAL, INTENT (OUT)            :: error
  INTEGER, INTENT (IN)             :: n_rad_wvl
  REAL (KIND=dp), DIMENSION (sig_idx,n_rad_wvl), INTENT (INOUT) :: curr_rad_spec

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp)                             :: tmpwave
  REAL (KIND=dp), DIMENSION (n_rad_wvl)      :: allwaves, locshi, locsqu, locspec
  REAL (KIND=dp), DIMENSION (max_calfit_idx, 2) :: tmp_fitvar
  INTEGER                                    :: npoints, i, iwin, fidx,  j, &
       lidx, errstat = pge_errstat_ok, nstep, islit, fpos, lpos, n_fit_pts, fslit, &
       lslit, ios, finter, linter, npoly, solfit_exval
  CHARACTER(LEN=maxchlen)                    :: tmpchar, fname
  LOGICAL        :: calfname_exist
  LOGICAL, SAVE  :: wrt_to_screen, wrt_to_file, slitcal, first = .TRUE.

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=12), PARAMETER :: modulename = 'rad_fit_vary'

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  IF (first) THEN
     fixslitcal = .TRUE.;  wrt_to_screen = .FALSE.; wrt_to_file = .FALSE.; slitcal=.TRUE.
     fitvar_sol = fitvar_sol_init
     
     ! need to restore back the initial fitting variables and bounds
     IF (wavcal_sol) THEN
        lo_sunbnd  = lo_sunbnd_init
        up_sunbnd  = up_sunbnd_init
        fitvar_sol_saved = fitvar_sol_init     
     END IF
     
     ! find the number of actual used fitting variables
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
  allwaves = curr_rad_spec(wvl_idx, 1 : n_rad_wvl) 
 
  ! Determine if file exists or not
  fname = TRIM(ADJUSTL(rslit_fname)) // currpixchar // '.dat'
  INQUIRE (FILE=TRIM(ADJUSTL(fname)), EXIST=calfname_exist)

!  ! Open a file for saving/reading the results
!  OPEN (UNIT=rslit_unit, FILE=TRIM(ADJUSTL(fname)), STATUS='UNKNOWN', IOSTAT=errstat)
!  IF (errstat /= pge_errstat_ok) THEN
!     WRITE(*, *) modulename, ' : cannot open rslit_fname', TRIM(ADJUSTL(fname))
!     error = .TRUE.; RETURN
!  END IF
     
  IF (slit_redo .OR. .NOT. calfname_exist) THEN
!     IF (which_slit == 2) then
!        WRITE(rslit_unit, '(A)') &
!        ' wave    vgl     vgr     hwl    hwr   shift    squeeze      sin     rms   exval'
!     ELSE
!        WRITE(rslit_unit, '(A)') &
!        ' wave    hwe     asym  shift    squeeze      sin      rms   exval '
!     END IF

     islit  = 0               ! number of sucessful calibrations      
     fidx = 1                 ! first pixel
     DO iwin = 1, numwin     

        IF (scnwrt) WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', iwin, winlim(iwin,1), &
             winlim(iwin,2), nradpix(iwin)
        lidx = fidx + nradpix(iwin) - 1

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

           npoints = lpos - fpos + 1
           fitwavs   (1:npoints) = curr_rad_spec(wvl_idx, fpos:lpos)
           currspec  (1:npoints) = curr_rad_spec(spc_idx, fpos:lpos)
           fitweights(1:npoints) = curr_rad_spec(sig_idx, fpos:lpos)
           fitvar_sol =  fitvar_sol_saved
           CALL cal_fit_one (npoints, n_fitvar_sol, wrt_to_screen, wrt_to_file,&
                slitcal, rslit_unit, tmpwave, tmp_fitvar, solfit_exval)
           
           IF (solfit_exval > 0) THEN
              islit = islit + 1;   slitwav_rad(islit) = tmpwave
              radslitfit(islit, 1:max_calfit_idx, 1:2) = &
                   tmp_fitvar(1:max_calfit_idx, 1:2)
           END IF
           
           fpos = fpos + n_slit_step           
        END DO
        
        fidx = lidx + 1
     END DO

     nslit_rad = islit 
  ELSE

     READ (rslit_unit, '(A)') tmpchar

     islit = 0
     DO 
        islit = islit + 1
        IF (which_slit == 2) THEN
           READ (rslit_unit, *, IOSTAT = ios) slitwav_rad(islit), &
                radslitfit(islit, vgl_idx, 1), radslitfit(islit, vgr_idx, 1), &
                radslitfit(islit, hwl_idx, 1), radslitfit(islit, hwr_idx, 1), &
                radslitfit(islit, shi_idx, 1), radslitfit(islit, squ_idx, 1)
        ELSE
           READ (rslit_unit, *, IOSTAT = ios) slitwav_rad(islit), &
                radslitfit(islit, hwe_idx, 1), radslitfit(islit, asy_idx, 1), &
                radslitfit(islit, shi_idx, 1), radslitfit(islit, squ_idx, 1)
        END IF

        IF (ios < 0) THEN
           islit = islit - 1; nslit_rad = islit; EXIT  ! end of file
        ELSE IF (ios > 0) THEN
           WRITE(*, *) modulename, ': read slit file error'
           error = .TRUE. ; RETURN
        END IF

     END DO

  END IF
  !CLOSE(rslit_unit)

  IF (smooth_slit) THEN

     fidx = 1
     DO iwin = 1, numwin
     
        lidx = fidx + nradpix(iwin) - 1
     
        fslit = MINVAL(MINLOC(slitwav_rad(1:nslit_rad), &
             MASK=(slitwav_rad(1:nslit_rad) >= winlim(iwin, 1))))
        lslit = MINVAL(MAXLOC(slitwav_rad(1:nslit_rad), &
             MASK=(slitwav_rad(1:nslit_rad) <= winlim(iwin, 2))))
     
        ! use a 4-order polynomial
        IF (lslit > fslit + 3) THEN
           poly_order = 4
           DO i = hwe_idx, hwr_idx
              IF (radslitfit(fslit,i,1)==0. .AND. radslitfit(lslit,i,1)==0.0) CYCLE
              IF (wavcal_sol .AND. (i == shi_idx .OR. i == squ_idx)) CYCLE
              
              npoly = lslit - fslit + 1
              locspec(1:npoly) = radslitfit(fslit:lslit, i, 1)
              CALL subtract_poly_meas ( slitwav_rad(fslit:lslit), npoly, &
                   locspec(1:npoly), 1, npoly )
              radslitfit(fslit:lslit, i, 1) = radslitfit(fslit:lslit, i, 1) &
                   - locspec(1:npoly) 

              !DO j = 3, npoly - 2
              !  locspec(j) = (radslitfit(fslit+j-1, i, 1)*6_dp + &
              !    radslitfit(fslit+j-2, i, 1)*4_dp &
              !  + radslitfit(fslit+j-3, i, 1) + radslitfit(fslit+j, i, 1)*4_dp &
              !  + radslitfit(fslit+j+1, i, 1)) / 15_dp
              !END DO
             radslitfit(fslit:lslit, i, 1) =  locspec(1:npoly)

           END DO
        ELSE
           DO i = hwe_idx, hwr_idx
              radslitfit(fslit:lslit, i, 1) = SUM(radslitfit(fslit:lslit, i, 1)) / &
                   (lslit - fslit + 1.0)
           END DO
        ENDIF
     
        fidx = lidx + 1
     END DO
  END IF

  ! to be used in voigt.f90 or gauss.f90
  nslit = nslit_rad; slitwav = slitwav_rad; slitfit = radslitfit

  !OPEN (unit=77, file='rslit_after_smooth.dat')
  !WRITE(77, *) nslit_rad
  !DO i = 1, nslit_rad
  !   WRITE(77, '(3d14.6)') slitwav_rad(i), radslitfit(i, hwe_idx, 1), &
  !        radslitfit(i, shi_idx, 1)
  !END DO
  !CLOSE (77)

  ! squeeze and shift wavelength position

  IF (wavcal_sol .OR. .NOT. wavcal) RETURN     ! do wavelength registration in next iteration

  fslit = 1; fidx = 1; 
  DO iwin = 1, numwin

     lidx = fidx + nradpix(iwin) - 1
     
     lslit = MINVAL(MAXLOC(slitwav_rad, MASK=(slitwav_rad(1:nslit_rad) &
          <= winlim(iwin, 2))))

     IF (lslit < fslit + 3) THEN
        locshi(fidx:lidx) = radslitfit(fslit, shi_idx, 1)
        locsqu(fidx:lidx) = radslitfit(fslit, squ_idx, 1)
     ELSE

        finter = MINVAL(MINLOC(allwaves, MASK = &
             (allwaves > slitwav_rad(fslit))))
        linter = MINVAL(MAXLOC(allwaves, MASK = &
             (allwaves < slitwav_rad(lslit))))
        CALL interpolation (lslit-fslit+1, slitwav_rad(fslit:lslit), &
             radslitfit(fslit:lslit, shi_idx, 1), linter-finter+1, &
             allwaves(finter:linter), locshi(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
     
        CALL interpolation ( lslit-fslit+1, slitwav_rad(fslit:lslit), &
             radslitfit(fslit:lslit, squ_idx, 1), linter-finter+1, &
             allwaves(finter:linter), locsqu(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
        
        IF (finter > fidx ) THEN
           locshi(fidx:finter-1)=radslitfit(fslit, shi_idx, 1)
           locsqu(fidx:finter-1)=radslitfit(fslit, squ_idx, 1)
        END IF
        
        IF (linter < lidx)  THEN
           locshi(linter+1:lidx)  = radslitfit(lslit, shi_idx, 1)
           locsqu(linter+1:lidx) =  radslitfit(lslit, squ_idx, 1)
        END IF
     END IF

     allwaves(fidx:lidx) = (allwaves(fidx:lidx) - locshi(fidx:lidx)) &
          / (1.0 + locsqu(fidx:lidx))
     fidx = lidx + 1; fslit = lslit + 1     
  END DO

  curr_rad_spec(wvl_idx, 1:n_rad_wvl) = allwaves
  
  RETURN
END SUBROUTINE rad_fit_vary
