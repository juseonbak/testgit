SUBROUTINE rad_wavcal_vary (wavcal_unit, n_rad_wvl, curr_rad_spec, error )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : maxchlen
  USE OMSAO_indices_module,    ONLY : max_calfit_idx, shi_idx, &
       squ_idx, wvl_idx, spc_idx, sig_idx, hwe_idx, hwr_idx, vgr_idx, &
       vgl_idx, asy_idx, hwl_idx, bl0_idx, bl1_idx, bl2_idx, bl3_idx, &
       sc0_idx, sc1_idx, sc2_idx, sc3_idx
  USE OMSAO_variables_module,  ONLY : n_fitvar_sol, fitvar_sol, &
       mask_fitvar_sol, lo_sunbnd, up_sunbnd, sswav_rad, nwavcal_rad, &
       wavcal_fit_pts, n_wavcal_step, smooth_slit,  wavcal_redo, &
       wavcal_fname, nradpix, winlim, poly_order, radwavfit, slitfit, &
       fixslitcal, which_slit, debug_boreas, fitweights, fitwavs, &
       currspec, fitvar_sol_saved, nslit, slitwav, fitvar_sol_init, &
       numwin, curr_sol_spec, currpixchar, scnwrt
  USE OMSAO_errstat_module
       

  IMPLICIT NONE

  ! ================
  ! Input/OUTPUT variables
  ! ================
  INTEGER, INTENT (IN)  :: n_rad_wvl
  INTEGER, INTENT (IN)  :: wavcal_unit
  LOGICAL, INTENT(OUT)  :: error
  REAL (KIND=dp), DIMENSION (sig_idx,n_rad_wvl), &
       INTENT (INOUT) :: curr_rad_spec

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp)                           :: tmpwave
  REAL (KIND=dp), DIMENSION (6, n_rad_wvl) :: tmpslit 
  REAL (KIND=dp), DIMENSION (n_rad_wvl) :: allwaves, locshi, locsqu, locspec
  REAL (KIND=dp), DIMENSION (max_calfit_idx, 2) :: tmp_fitvar
  INTEGER  :: npoints, i, iwin, fidx,  j, lidx, errstat = pge_errstat_ok, &
       nstep, iwavcal, fpos, lpos, n_fit_pts, fwavcal, lwavcal, ios, finter,&
       linter, npoly, solfit_exval
  INTEGER, DIMENSION(6)   :: slitind = (/hwe_idx, asy_idx, vgl_idx, &
       vgr_idx, hwl_idx, hwr_idx/)
  CHARACTER(LEN=maxchlen) :: tmpchar, fname
  LOGICAL :: wrt_to_screen, wrt_to_file, slitcal, calfname_exist = .TRUE.

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=14), PARAMETER :: modulename = 'rad_wavcal_vary'

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  error = .FALSE. ;  wrt_to_screen = .FALSE.
  wrt_to_file = .FALSE.; slitcal = .FALSE.

  ! add more variables (for scale and base line terms)
  ! fitvar_sol(bl0_idx) = 2.4626E-1; fitvar_sol(bl1_idx) = -1.716E-2
  ! fitvar_sol(bl2_idx) = -8.86E-05
  ! lo_sunbnd(bl0_idx:bl2_idx) = -1.0D+99; up_sunbnd(bl0_idx:bl2_idx) = 1.0D+99

  allwaves = curr_rad_spec(wvl_idx, 1 : n_rad_wvl)
  fitvar_sol_saved  =  fitvar_sol_init
  fitvar_sol = fitvar_sol_init

  IF (which_slit /= 4) THEN
     tmpslit = 0.0
     fpos = MINVAL(MINLOC(allwaves, MASK=(allwaves >= slitwav(1))))
     lpos = MINVAL(MAXLOC(allwaves, MASK=(allwaves <= slitwav(nslit))))
     IF (fpos > 1)  THEN
        DO i = 1, 6
           tmpslit(i, 1:fpos-1) = slitfit(1, slitind(i), 1)
        ENDDO
     ENDIF
     IF (lpos < n_rad_wvl)  THEN
        DO i = 1, 6
           tmpslit(i, lpos+1:n_rad_wvl) = slitfit(nslit, slitind(i), 1)
        ENDDO
     ENDIF
     DO i = 1, 6
        IF (slitfit(1, slitind(i), 1) /= 0.0) &
             CALL BSPLINE(slitwav(1:nslit), slitfit(1:nslit, slitind(i), 1), &
             nslit, allwaves(fpos:lpos), tmpslit(i, fpos:lpos), lpos-fpos+1,errstat)
        IF (errstat < 0) THEN
           WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat; STOP
        ENDIF
     ENDDO
  ENDIF

  ! Determine if file exists or not
  fname = TRIM(ADJUSTL(wavcal_fname)) // currpixchar // '.dat'
  INQUIRE (FILE=TRIM(ADJUSTL(fname)), EXIST=calfname_exist)

  ! Open a file for saving/reading the results
!  OPEN (UNIT=wavcal_unit, FILE=TRIM(ADJUSTL(fname)), STATUS='UNKNOWN', IOSTAT=errstat)
!  IF (errstat /= pge_errstat_ok) THEN
!     WRITE(*, *) modulename, ' : Cannot open wavcal_fname', TRIM(ADJUSTL(fname))
!     error = .TRUE.; RETURN
!  END IF

  IF (wavcal_redo .OR. .NOT. calfname_exist) THEN

!     WRITE(wavcal_unit, '(A)') ' wave  shift    squeeze     sin      rms    exval'

     ! fix variables for hwe, asy, vgr, vgl, hwr, hwl
     fitvar_sol(hwe_idx:asy_idx) = 0.0
     lo_sunbnd(hwe_idx:asy_idx) = 0.0
     up_sunbnd(hwe_idx:asy_idx) = 0.0
     fitvar_sol(vgl_idx:hwr_idx) = 0.0
     lo_sunbnd(vgl_idx:hwr_idx) = 0.0
     up_sunbnd(vgl_idx:hwr_idx) = 0.0

     ! find the number of actual used fitting variables     
     n_fitvar_sol = 0
     DO i = 1, max_calfit_idx
        IF (lo_sunbnd(i) < up_sunbnd(i) ) THEN
           n_fitvar_sol =  n_fitvar_sol + 1
           mask_fitvar_sol(n_fitvar_sol) = i
        END IF
     END DO

     iwavcal  = 0               ! number of sucessful calibrations      
     fidx = 1                   ! first pixel
     DO iwin = 1, numwin  
        lidx = fidx + nradpix(iwin) - 1     
        IF (scnwrt) WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', iwin, winlim(iwin,1), &
             winlim(iwin,2), nradpix(iwin)
        
        ! do each fit using points (fpos:lpos)
        fpos = fidx
        
        DO WHILE (fpos < lidx)

           DO i = fpos+1, fpos + wavcal_fit_pts - 1
              ! until end of window or sudden gap of > 2.0 um
              IF (i > lidx) EXIT
              IF (allwaves(i) - allwaves(i-1) > 2.0) EXIT
           ENDDO
           lpos = i - 1; npoints = lpos - fpos + 1
           
           IF (npoints < wavcal_fit_pts / 2) THEN
              ! Either until the end with not enough points or gap behind
              fpos = lpos + 1; CYCLE
           ENDIF

               ! use variable slit width for solar wavelength calibration
           fitvar_sol = fitvar_sol_saved
           !fixslitcal = .FALSE.; fitvar_sol(slitind) = 0.0
           fixslitcal = .TRUE. 
           IF (which_slit /= 4) THEN
              fitvar_sol(slitind) = tmpslit(:, (fpos+lpos)/2)
              lo_sunbnd(slitind) = fitvar_sol(slitind)
              up_sunbnd(slitind) = fitvar_sol(slitind)
           ENDIF

           fitwavs   (1:npoints) = curr_rad_spec(wvl_idx, fpos:lpos)
           currspec  (1:npoints) = curr_rad_spec(spc_idx, fpos:lpos)
           fitweights(1:npoints) = curr_rad_spec(sig_idx, fpos:lpos)
           
           CALL cal_fit_one (npoints, n_fitvar_sol, wrt_to_screen, wrt_to_file,&
                slitcal, wavcal_unit, tmpwave, tmp_fitvar, solfit_exval)

           IF (solfit_exval > 0) THEN
              iwavcal = iwavcal + 1
              sswav_rad(iwavcal) = tmpwave
              radwavfit(iwavcal, 1:max_calfit_idx, 1:2) = tmp_fitvar(1:max_calfit_idx, 1:2)
           END IF

           fpos = fpos + n_wavcal_step           
        END DO
          
        fidx = lidx + 1
     END DO

     nwavcal_rad = iwavcal 
  ELSE

     READ (wavcal_unit, '(A)') tmpchar
     iwavcal = 0
     DO 
        iwavcal = iwavcal + 1
        READ (wavcal_unit, *, IOSTAT = ios) sswav_rad(iwavcal), &
             radwavfit(iwavcal, shi_idx, 1), radwavfit(iwavcal, squ_idx, 1)

        IF (ios < 0) THEN
           iwavcal = iwavcal - 1; nwavcal_rad = iwavcal; EXIT  ! end of file
        ELSE IF (ios > 0) THEN
           WRITE(*, *) modulename, ': read slit file error'
           error = .TRUE. ; RETURN
        END IF
     END DO

  END IF

  !CLOSE(wavcal_unit)

  ! smooth shi, squ for window 3 and window 4
  IF (smooth_slit) THEN
     fidx = 1
     DO iwin = 1,numwin
        lidx = fidx + nradpix(iwin) - 1 
     
        fwavcal = MINVAL(MINLOC(sswav_rad(1:nwavcal_rad), &
             MASK=(sswav_rad(1:nwavcal_rad) >= winlim(iwin, 1))))
        lwavcal = MINVAL(MAXLOC(sswav_rad(1:nwavcal_rad), &
             MASK=(sswav_rad(1:nwavcal_rad) <= winlim(iwin, 2))))
     
        ! use a 4-order polynomial
        IF (lwavcal > fwavcal + 3) THEN
           IF (iwin == 1) THEN
              poly_order = 4
           ELSE
              poly_order = 3
           ENDIF
           DO i = shi_idx, squ_idx
              IF (radwavfit(fwavcal, i, 1) == 0 .AND. &
                   radwavfit(lwavcal, i, 1) == 0.0) CYCLE
              
              npoly = lwavcal - fwavcal + 1
              locspec(1:npoly) = radwavfit(fwavcal:lwavcal, i, 1)
              CALL subtract_poly_meas ( sswav_rad(fwavcal:lwavcal), npoly,&
                   locspec(1:npoly), 1, npoly )
              radwavfit(fwavcal:lwavcal, i, 1) = radwavfit(fwavcal:lwavcal, i, 1)&
                   - locspec(1:npoly)
           END DO
        ELSE
           DO i = shi_idx, squ_idx
              radwavfit(fwavcal:lwavcal,i,1)=SUM(radwavfit(fwavcal:lwavcal, i, 1))/ &
                   (lwavcal-lwavcal + 1)
           END DO
        END IF
     
        fidx = lidx + 1
     END DO
  END IF

  ! squeeze and shift wavelength position
  fwavcal = 1; fidx = 1; 
  DO iwin = 1, numwin

     lidx = fidx + nradpix(iwin) - 1
     lwavcal = MINVAL(MAXLOC(sswav_rad(1:nwavcal_rad), &
          MASK=(sswav_rad(1:nwavcal_rad) <= winlim(iwin, 2))))

     IF (lwavcal < fwavcal + 3) THEN
        locshi(fidx:lidx) = radwavfit(fwavcal, shi_idx, 1)
        locsqu(fidx:lidx) = radwavfit(fwavcal, squ_idx, 1)
     ELSE
 
        finter = MINVAL(MINLOC(allwaves, MASK = &
             (allwaves > sswav_rad(fwavcal))))  ! finter >=fidx
        linter = MINVAL(MAXLOC(allwaves, MASK = &
             (allwaves < sswav_rad(lwavcal))))  ! linter <=lidx

        IF (finter == 0)  CYCLE

        CALL interpolation (  lwavcal-fwavcal+1, sswav_rad(fwavcal:lwavcal), &
             radwavfit(fwavcal:lwavcal, shi_idx, 1), linter-finter+1,&
             allwaves(finter:linter), locshi(finter:linter),   errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
     
        CALL interpolation ( lwavcal-fwavcal+1,   sswav_rad(fwavcal:lwavcal), &
             radwavfit(fwavcal:lwavcal, squ_idx, 1), linter-finter+1, &
             allwaves(finter:linter), locsqu(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
        
        IF (finter - fidx > 0) THEN
           locshi(fidx:finter-1)=radwavfit(fwavcal, shi_idx, 1)
           locsqu(fidx:finter-1)=radwavfit(fwavcal, squ_idx, 1)
        END IF
        
        IF (linter < lidx)  THEN
           locshi(linter+1:lidx)  = radwavfit(lwavcal, shi_idx, 1)
           locsqu(linter+1:lidx) =  radwavfit(lwavcal, squ_idx, 1)
        END IF

     END IF
        
     allwaves(fidx:lidx) = (allwaves(fidx:lidx) - locshi(fidx:lidx) ) &
          / ( 1.0 + locsqu(fidx:lidx))
     
     fidx = lidx + 1; fwavcal = lwavcal + 1     
  END DO

  curr_rad_spec(wvl_idx, 1:n_rad_wvl)  = allwaves
  
  !OPEN (unit=77, file='rwavcal_after_smooth.dat')
  !write(77, *) nwavcal_rad
  !DO i = 1, nwavcal_rad
  !   WRITE(77, '(3d14.6)') sswav_rad(i), radwavfit(i, shi_idx, 1)
  !END DO
  !CLOSE (77)

  RETURN
END SUBROUTINE rad_wavcal_vary
