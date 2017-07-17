SUBROUTINE rad_shisqu (n_rad_wvl, curr_rad_spec)

  USE OMSAO_precision_module,     ONLY: dp
  USE OMSAO_indices_module,       ONLY: wvl_idx, sig_idx, shi_idx, squ_idx
  USE OMSAO_variables_module,     ONLY: nradpix, winlim, radwavfit, &
       nwavcal_rad, sswav_rad, numwin
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ================
  ! Input parameters
  ! ================
  INTEGER, INTENT (IN) :: n_rad_wvl

  ! ===================
  ! Modified parameters
  ! ===================
  REAL (KIND=dp), DIMENSION(sig_idx,n_rad_wvl), INTENT(INOUT) :: curr_rad_spec
  
  ! ===============
  ! Local variables
  ! ===============
  INTEGER ::  errstat = pge_errstat_ok, fidx, lidx, linter, finter,&
       fwavcal, lwavcal, i, iwin
  REAL (KIND=dp), DIMENSION (n_rad_wvl) ::  allwaves, locshi, locsqu

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=10), PARAMETER :: modulename = 'rad_shisqu'

  allwaves(1:n_rad_wvl) = curr_rad_spec(wvl_idx, 1:n_rad_wvl)

  ! squeeze and shift wavelength position
  fidx = 1; 
  fwavcal = fwavcal + 1
  DO iwin = 1, numwin

     lidx =  fidx + nradpix(iwin) - 1
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

        CALL interpolation (lwavcal-fwavcal+1, sswav_rad(fwavcal:lwavcal),&
             radwavfit(fwavcal:lwavcal, shi_idx, 1), linter-finter+1,&
             allwaves(finter:linter), locshi(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0); STOP 1
        END IF
     
        CALL interpolation (lwavcal-fwavcal+1, sswav_rad(fwavcal:lwavcal), &
             radwavfit(fwavcal:lwavcal, squ_idx, 1), linter-finter+1, &
             allwaves(finter:linter), locsqu(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0); STOP 1
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
        
     allwaves(fidx:lidx) = (allwaves(fidx:lidx) - locshi(fidx:lidx)) &
          / ( 1.0 + locsqu(fidx:lidx)) 
     fidx = lidx + 1; fwavcal = lwavcal + 1     
  END DO

  curr_rad_spec(wvl_idx, 1:n_rad_wvl) = allwaves

  RETURN
END SUBROUTINE rad_shisqu
