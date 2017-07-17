! ============================================================
! Convolve a spectrum with left/right asymmetric gaussian/voigt
! slit width with four parameter, al/alr, hwel/hwer
! al/ar must >= 0, hwer/hwel >=0
! hwer/hwel = 0 : no convolution
! al/ar = 0 : gaussian; > 0 : voigt
! ============================================================ 


SUBROUTINE asym_voigt(wvlarr, specarr, specmod, npoints, vgl, vgr, hwl, hwr)
 
  USE OMSAO_precision_module
  !USE OMSAO_variables_module, ONLY : lo_sunbnd, up_sunbnd
  !USE OMSAO_indices_module,   ONLY : vgr_idx, vgl_idx

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,        INTENT (IN) :: npoints
  REAL (KIND=dp), INTENT (IN) :: hwl, hwr, vgr, vgl
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN) :: wvlarr, specarr

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=dp), DIMENSION (npoints), INTENT (OUT) :: specmod

  ! ===============
  ! Local variables
  ! ===============
  INTEGER :: i, j, j1, j2, num_slit, mslit, sslit, eslit, odd
  REAL (KIND=dp) :: delwvl, slitsum,  maxslit, voigt, left_norm,&
       right_norm, norm, vg, hw1e
  REAL (KIND=dp), DIMENSION (npoints) :: slit, locwvl, locsli, xxx
  INTEGER,        DIMENSION (npoints) :: idx

  ! --------------------------------------------------------
  ! Initialize output variable (default for "no convolution"
  ! --------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  
  ! -----------------------------------------------
  ! No convolution if Halfwidth @ 1/e is 0
  ! -----------------------------------------------
  IF ( hwl == 0.0 .AND. hwr == 0.0 ) RETURN

  ! --------------------------------------------------------------
  ! Find the number of spectral points that fall within a Gaussian
  ! slit function with values >= 0.01. Remember that we have an
  ! asymmetric Gaussian, so we create a wavelength array symmetric
  ! around 0. The spacing is provided by the equidistant WVLARR.
  ! --------------------------------------------------------------
  left_norm =  voigt(0.0, vgl, hwl)
  right_norm = voigt(0.0, vgr, hwr)

  delwvl = wvlarr(2) - wvlarr(1) ; locsli = 0.0
  DO i = -npoints/2, npoints/2
     j = MIN ( i + npoints/2 + 1, npoints )
     locwvl(j) = delwvl * REAL(i,KIND=dp)
     IF (locwvl(j) < 0) THEN
        vg = vgl; hw1e = hwl; norm = left_norm
     ELSE
        vg = vgr; hw1e = hwr; norm = right_norm
     END IF
     locsli(j) = voigt(locwvl(j), vg, hw1e) / norm
     !WRITE(*, *) locwvl(j), locsli(j)
  END DO

  ! ------------------------------------------------------------------------
  ! Find the array entries that mark the region of SLIT >= 0.0001*MAX(SLIT).
  ! We start by initializing the auxilliary IDX as IDX(n) = n. We then set
  ! all those entries to 0 that are outside the accepted minimum value of 
  ! the slit function. The start and end indices simply follow as the 
  ! minimum and maximum values of IDX where IDX /= 0.
  ! ------------------------------------------------------------------------
  maxslit = MAXVAL(locsli(1:npoints))
  sslit = 0 ; eslit = 0 ; idx = (/ (i, i = 1, npoints) /)
  WHERE ( locsli < 0.01*maxslit )
     idx = 0
  END WHERE
  sslit = MINVAL( idx, MASK = (idx /= 0) )
  eslit = MAXVAL( idx, MASK = (idx /= 0) )

  ! ----------------------------------------------------------
  ! Compute number of slit function points, and the final slit
  ! function array; the latter will be normalized to 1.
  ! ----------------------------------------------------------
  num_slit = eslit - sslit + 1
  slit(1:num_slit) = locsli(sslit:eslit)
  slitsum = SUM(slit(1:num_slit))
  IF ( slitsum > 0.0 ) slit(1:num_slit) = slit(1:num_slit) / slitsum

  ! --------------------------------------------------------
  ! Find the maximum of the slit. Note that since we may be 
  ! with an asymmetric Gaussian, the maximum must not always
  ! be in the center of the Gaussian function.
  ! --------------------------------------------------------
  mslit = MAXVAL( MAXLOC ( slit(1:num_slit) ) )

  !WRITE(*, *) sslit, eslit, num_slit, mslit

  ! ---------------------------------------------------------------
  ! Convolve spectrum. First do the middle part, where we have full
  ! overlap coverage of the slit function. Again, remember the
  ! asymmetry of the Gaussian, which makes impossible a simple
  ! 50-50 division of the summation interval.
  ! ---------------------------------------------------------------

  ! Make a local copy of the NUM_SLIT spectrum points to be convolved
  ! with the slit function. The spectrum points to be convolved are
  ! arranged such that the updated index corresponds to the maximum
  ! of the slit function (MSLIT). For simplicity we reflect the
  ! spectrum at the array end points.

  ! ----------------------------------------------------
  ! Loop over all points of the spectrum to be convolved
  ! ----------------------------------------------------
  DO i = 1, npoints
     !locsli(1:num_slit) = 0.0

     ! First do the right half of the slit function
     DO j = mslit, num_slit
        j1 = i+(j-mslit) ; IF ( j1 > npoints ) j1 = npoints - MOD(j1, npoints)
        !locsli(j) = specarr(j1)
        idx(j) = j1
     END DO

     ! Now the left half of the slit function
     DO j = mslit-1, 1, -1
        j2 = i+(j-mslit) ; IF ( j2 < 1 ) j2 = ABS(j2) + 2
        !locsli(j) = specarr(j2)
        idx(j) = j2
     END DO

     specmod(i) = DOT_PRODUCT(slit(1:num_slit), specarr(idx(1:num_slit)))
     
  END DO
  !WRITE(*, *) 'num_slit = ', num_slit
  !DO i = 1, num_slit
  !   WRITE(*, '(2d14.6)') slit(i), locsli(i)
  !END DO
  !STOP

  RETURN
END SUBROUTINE asym_voigt



FUNCTION voigt (x, a, hw1e) RESULT (v)
 
 ! x:    Distance from the center     (+/-x gives the same result)
 ! a:    Lorentz HWHM / xg            (>= 0, = 0 for gaussian)
 ! hw1e: Gaussian half width at 1 / e (>=0, =0 for no concolution)

 ! ---------------------------------------------------------------------------
 ! Kelly Chance      May 31, 1994
 ! Evaluation of the Voigt function, v, and its first and second
 ! derivatives, dv and d2v, using algorithm 363 of the Collected
 ! Algorithms from CACM.  Also see "Efficient computation of the complex
 ! error function", W. Gautschi, SIAM J. Numer. Anal. 7, 187-198 (1970).
 ! The accuracy is stated as "10 decimal places after the decimal point,
 ! or better."  Certification shows it to be 10 significant figures or
 ! better in most circumstances.
 !
 ! The Voigt function is computed for an array of x values and any positive
 ! value of a, where x = (sigma - sigma0) / xg, xg = 4.30140E-7 * sigma0 *
 ! sqrt (T/m); a = (Lorentz HWHM)/xg.  This routine gives a Voigt profile
 ! normalized to area xg (in the Doppler limit the line center is always
 ! 1/sqrt (pi).
 ! The present version closely follows that of D.G. Hummer, including the
 ! approximation for a near 0.  The second derivative has been added here.
 ! ---------------------------------------------------------------------------

  USE OMSAO_precision_module
  IMPLICIT NONE
 
  ! ===============
  ! Input variables
  ! ===============
  REAL (KIND=dp), INTENT(IN) :: x, a, hw1e

  ! ===============
  ! Result variables
  ! ===============
  REAL (KIND=dp)             :: v

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                    :: capn, nu, in, n, np1
  REAL (KIND=dp)             :: lamda, r1, r2, t1, t2, c, &
       s, s1, s2, xnew, sfac, absx, h, h2, anew, hnew
  LOGICAL                    :: b
  REAL (KIND=dp), PARAMETER  :: sqrtpi =1.77245385090551, &
       twooverpi = 0.63661977236758,  fouroverpi = 1.27323954473516

  
  !WRITE(*, *) a, x, hw1e
  anew = a  ! ABS(a); 
  hnew =  ABS(hw1e); 

  IF (x /= 0.0) THEN
     xnew = ABS(x) / hnew
  ELSE 
     xnew = x
  END IF

  IF (anew == 0.0) THEN                ! a = 0., guassian

     v = EXP(- xnew * xnew ) / sqrtpi

  ELSE                                 ! voigt

     sfac = 1.0 - anew / 4.29     
     absx = xnew

     IF ((anew < 4.29) .AND. (absx < 5.33)) THEN
        s = sfac * SQRT (1.0 - (xnew / 5.33)**2)
        h = 1.6 * s
        h2 = 2.0 * h
        capn = INT (6.0 + 23.0 * s)
        lamda = h2 ** capn
        nu = INT(9.0 + 21.0 * s)
     ELSE
        h = 0.0; lamda = 0.0
        capn = 0
        nu = 8
     END IF

     b = ((h == 0.0) .OR. (lamda == 0.0))
     r1 = 0.0
     r2 = 0.0
     s1 = 0.0
     s2 = 0.0
     n = nu

     DO in = 1, nu + 1
        np1 = n + 1
        t1 = anew + h + REAL (np1, KIND=dp) * r1
        t2 = absx - REAL (np1, KIND=dp) * r2
        c = 0.5 / (t1 * t1 + t2 * t2)
        r1 = c * t1
        r2 = c * t2
        IF ((h > 0.0) .AND. (n <= capn)) THEN
           t1 = lamda + s1
           s1 = r1 * t1 - r2 * s2
           s2 = r2 * t1 + r1 * s2
           lamda = lamda / h2
        END IF
        n = n - 1
     END DO

     IF (b) THEN
        v = twooverpi * r1
     ELSE
        v = twooverpi * s1
     END IF
  END IF

  RETURN
END FUNCTION voigt



SUBROUTINE asym_voigt_multi (wvlarr, specarr, specmod, npoints)
 
  USE OMSAO_precision_module
  USE OMSAO_variables_module,  ONLY  : solwinfit, winlim, numwin
  USE OMSAO_indices_module,    ONLY  : hwl_idx, hwr_idx, vgr_idx, vgl_idx

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                             INTENT (IN) :: npoints
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN) :: wvlarr, specarr

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=dp), DIMENSION (npoints), INTENT (OUT) :: specmod

  ! ===============
  ! Local variables
  ! ===============
  INTEGER        :: i, j, j1, j2, num_slit, mslit, sslit, eslit, odd
  REAL (KIND=dp) :: delwvl, slitsum,  maxslit, voigt, upbnd, left_norm, &
       right_norm, norm, vg, hw1e, fidx, lidx, iwin, vgr, vgl, hwr, hwl
  REAL (KIND=dp), DIMENSION (npoints) :: slit, locwvl, locsli, xxx
  INTEGER,        DIMENSION (npoints) :: idx

  ! --------------------------------------------------------
  ! Initialize output variable (default for "no convolution"
  ! --------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  

  fidx = 1
  DO iwin = 1, numwin

     hwl = solwinfit(iwin, hwl_idx, 1); hwr = solwinfit(iwin, hwr_idx, 1)
     vgr = solwinfit(iwin, vgr_idx, 1); vgl = solwinfit(iwin, vgl_idx, 1)

     IF (iwin < numwin) THEN
        upbnd = (winlim(iwin, 2) + winlim(iwin+1, 1))/2.0 
        lidx = MINVAL(MAXLOC(wvlarr, MASK=(wvlarr <= upbnd )))
     ELSE 
        lidx = npoints
     END IF

     ! -----------------------------------------------
     ! No convolution if Halfwidth @ 1/e is 0
     ! -----------------------------------------------
     IF ( hwl == 0.0 .AND. hwr == 0.0 ) RETURN
     
     ! --------------------------------------------------------------
     ! Find the number of spectral points that fall within a Gaussian
     ! slit function with values >= 0.01. Remember that we have an
     ! asymmetric Gaussian, so we create a wavelength array symmetric
     ! around 0. The spacing is provided by the equidistant WVLARR.
     ! --------------------------------------------------------------
     left_norm =  voigt(0.0, vgl, hwl)
     right_norm = voigt(0.0, vgr, hwr)

     delwvl = wvlarr(2) - wvlarr(1) ; locsli = 0.0
     DO i = -npoints/2, npoints/2
        j = MIN ( i + npoints/2 + 1, npoints )
        locwvl(j) = delwvl * REAL(i,KIND=dp)
        IF (locwvl(j) < 0) THEN
           vg = vgl; hw1e = hwl; norm = left_norm
        ELSE
           vg = vgr; hw1e = hwr; norm = right_norm
        END IF
        locsli(j) = voigt(locwvl(j), vg, hw1e) / norm
     END DO

     ! ------------------------------------------------------------------------
     ! Find the array entries that mark the region of SLIT >= 0.0001*MAX(SLIT).
     ! We start by initializing the auxilliary IDX as IDX(n) = n. We then set
     ! all those entries to 0 that are outside the accepted minimum value of 
     ! the slit function. The start and end indices simply follow as the 
     ! minimum and maximum values of IDX where IDX /= 0.
     ! ------------------------------------------------------------------------
     maxslit = MAXVAL(locsli(1:npoints))
     sslit = 0 ; eslit = 0 ; idx = (/ (i, i = 1, npoints) /)
     WHERE ( locsli < 0.01*maxslit )
        idx = 0
     END WHERE
     sslit = MINVAL( idx, MASK = (idx /= 0) )
     eslit = MAXVAL( idx, MASK = (idx /= 0) )
     
     ! ----------------------------------------------------------
     ! Compute number of slit function points, and the final slit
     ! function array; the latter will be normalized to 1.
     ! ----------------------------------------------------------
     num_slit = eslit - sslit + 1
     slit(1:num_slit) = locsli(sslit:eslit)
     slitsum = SUM(slit(1:num_slit))
     IF ( slitsum > 0.0 ) slit(1:num_slit) = slit(1:num_slit) / slitsum
     
     ! --------------------------------------------------------
     ! Find the maximum of the slit. Note that since we may be 
     ! with an asymmetric Gaussian, the maximum must not always
     ! be in the center of the Gaussian function.
     ! --------------------------------------------------------
     mslit = MAXVAL( MAXLOC ( slit(1:num_slit) ) )
     
     ! ---------------------------------------------------------------
     ! Convolve spectrum. First do the middle part, where we have full
     ! overlap coverage of the slit function. Again, remember the
     ! asymmetry of the Gaussian, which makes impossible a simple
     ! 50-50 division of the summation interval.
     ! ---------------------------------------------------------------
     
     ! Make a local copy of the NUM_SLIT spectrum points to be convolved
     ! with the slit function. The spectrum points to be convolved are
     ! arranged such that the updated index corresponds to the maximum
     ! of the slit function (MSLIT). For simplicity we reflect the
     ! spectrum at the array end points.
     
     ! ----------------------------------------------------
     ! Loop over all points of the spectrum to be convolved
     ! ----------------------------------------------------
     DO i = fidx, lidx
        !locsli(1:num_slit) = 0.0
        
        ! First do the right half of the slit function
        DO j = mslit, num_slit
           j1 = i+(j-mslit)
           IF ( j1 > npoints ) j1 = npoints - MOD(j1, npoints)
           !locsli(j) = specarr(j1)
           idx(j) = j1
        END DO
     
        ! Now the left half of the slit function
        DO j = mslit-1, 1, -1
           j2 = i+(j-mslit) ; IF ( j2 < 1 ) j2 = ABS(j2) + 2
           !locsli(j) = specarr(j2)
           idx(j) = j2
        END DO
        
        specmod(i) = DOT_PRODUCT(slit(1:num_slit), specarr(idx(1:num_slit)))               
     END DO
     fidx = lidx + 1
  END DO

  RETURN
END SUBROUTINE asym_voigt_multi


SUBROUTINE asym_voigt_vary (wvlarr, specarr, specmod, npoints)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : slitwav, slitfit, nslit, winlim, &
       numwin, debug_boreas
  USE OMSAO_indices_module,   ONLY : vgr_idx, vgl_idx, hwr_idx, hwl_idx
  USE OMSAO_errstat_module
  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                             INTENT (IN) :: npoints
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN) :: wvlarr, specarr

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=dp), DIMENSION (npoints), INTENT (OUT) :: specmod

  ! ===============
  ! Local variables
  ! ===============
  INTEGER        :: i, j, j1, j2, num_slit, mslit, sslit, eslit, odd,&
       errstat, fidx, lidx, fslit, lslit, finter, linter, iwin
  REAL (KIND=dp) :: delwvl, slitsum,  maxslit, hwr, hwl, hw1e, &
       vg, voigt, vgr, vgl, left_norm, right_norm, norm
  REAL (KIND=dp), DIMENSION (npoints) :: slit, locwvl, locsli, xxx, &
       lochwl, lochwr, locvgl, locvgr
  INTEGER,        DIMENSION (npoints) :: idx

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=16), PARAMETER :: modulename = 'asym_voigt_vary'

  ! --------------------------------------------------------
  ! Initialize output variable (default for "no convolution"
  ! --------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  
  ! -----------------------------------------------
  ! No convolution if Halfwidth @ 1/e is 0
  ! -----------------------------------------------
  IF ( MAXVAL(slitfit(1:nslit, hwr_idx, 1)) == 0.0  ) RETURN

  ! most spread scenario for GOME
  hwr = 0.50; hwl= 0.50; hw1e = 0.50;
  vgr = 0.2;  vgl=0.20; vg    = 0.20;

  ! -------------------------------------------------------------------
  ! Find the maximum number of spectral points that fall within a Voigt
  ! slit function with values >= 0.01. 
  ! -------------------------------------------------------------------
  delwvl = wvlarr(2) - wvlarr(1) ; locsli = 0.0
  locsli(1) = voigt(0.0, vg, hw1e)
  DO i = 2, npoints/2
     locwvl(i) = delwvl * (i - 1.0)
     locsli(i) = voigt(locwvl(i), vg, hw1e)
     IF (locsli(i) < 0.01 * locsli(1)) EXIT
  END DO

  ! ------------------------------------------------------------------------
  ! Find the array entries that mark the region of SLIT >= 0.0001*MAX(SLIT).
  ! We start by initializing the auxilliary IDX as IDX(n) = n. We then set
  ! all those entries to 0 that are outside the accepted minimum value of 
  ! the slit function. The start and end indices simply follow as the 
  ! minimum and maximum values of IDX where IDX /= 0.
  ! ------------------------------------------------------------------------
  !DO i = 1, npoints / 2
  !   IF (locsli(i) < 0.01 * locsli(1)) EXIT
  !END DO

  mslit = i; num_slit = 2 * mslit + 1

  ! get an array of slit variables (lochwr, lochwl, locvgr, locvgl) &
  ! for each wavelength position
  lochwr = 0.; lochwl=0.; locvgr = 0.; locvgl = 0.

  fidx = 1; lidx= npoints
     
  fslit = MINVAL(MINLOC(slitwav(1:nslit), &
       MASK=(slitwav(1:nslit) >= wvlarr(fidx) )))
  lslit = MINVAL(MAXLOC(slitwav(1:nslit), &
       MASK=(slitwav(1:nslit) <= wvlarr(lidx) )))
  
  ! no slit between fidx:lidx, should never happen
  IF (fslit<=0 .OR. fslit>nslit .OR. lslit<=0 &
       .OR. lslit>nslit .OR. fslit>lslit) THEN  
     WRITE(*, *) modulename, ' Not slit available for this window!!!'
     STOP
  ENDIF
      
  IF (lslit < fslit + 3) THEN  ! extrapolate, use the nearest value
     locvgr(fidx:lidx) = SUM(slitfit(fslit:lslit,vgr_idx,1))/(lslit-fslit+1)
     locvgl(fidx:lidx) = SUM(slitfit(fslit:lslit,vgl_idx,1))/(lslit-fslit+1)
     lochwr(fidx:lidx) = SUM(slitfit(fslit:lslit,hwr_idx,1))/(lslit-fslit+1)
     lochwl(fidx:lidx) = SUM(slitfit(fslit:lslit,hwl_idx,1))/(lslit-fslit+1)
  ELSE                         ! at least 4 slits for interpolation
    
     finter = MINVAL(MINLOC(wvlarr, MASK=(wvlarr >= slitwav(fslit))))
     linter = MINVAL(MAXLOC(wvlarr, MASK=(wvlarr <= slitwav(lslit))))
     
     !WRITE(*, *) lslit, fslit, finter, linter
     !WRITE(*, *) wvlarr(finter), wvlarr(linter), slitwav(fslit), slitwav(lslit)

     CALL interpolation (lslit-fslit+1, slitwav(fslit:lslit), &
          slitfit(fslit:lslit, vgr_idx, 1), linter-finter+1, &
          wvlarr(finter:linter), locvgr(finter:linter), errstat )
     IF (errstat > pge_errstat_warning ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0); STOP 1
     END IF
     
     CALL interpolation (lslit-fslit+1, slitwav(fslit:lslit),&
          slitfit(fslit:lslit, vgl_idx,1), linter-finter+1,  &
          wvlarr(finter:linter), locvgl(finter:linter), errstat)
     IF ( errstat > pge_errstat_warning ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0); STOP 1
     END IF
     
     CALL interpolation (lslit-fslit+1, slitwav(fslit:lslit), &
          slitfit(fslit:lslit, hwr_idx, 1),  linter-finter+1, &
          wvlarr(finter:linter), lochwr(finter:linter), errstat )
     IF ( errstat > pge_errstat_warning ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0); STOP 1
     END IF

     
     CALL interpolation (lslit-fslit+1, slitwav(fslit:lslit), &
          slitfit(fslit:lslit, hwl_idx, 1), linter-finter+1, &
          wvlarr(finter:linter), lochwl(finter:linter), errstat )
     IF ( errstat > pge_errstat_warning ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0); STOP 1
     END IF

     IF (finter > fidx) THEN
        locvgr(fidx:finter-1)=slitfit(fslit, vgr_idx, 1)
        locvgl(fidx:finter-1)=slitfit(fslit, vgl_idx, 1)
        lochwr(fidx:finter-1)=slitfit(fslit, hwr_idx, 1)
        lochwl(fidx:finter-1)=slitfit(fslit, hwl_idx, 1)
     END IF
     
     IF (linter < lidx)  THEN
        locvgr(linter+1:lidx)=slitfit(lslit, vgr_idx, 1)
        locvgl(linter+1:lidx)=slitfit(lslit, vgl_idx, 1)
        lochwr(linter+1:lidx)=slitfit(lslit, hwr_idx, 1)
        lochwl(linter+1:lidx)=slitfit(lslit, hwl_idx, 1)
     END IF
       
  END IF
  
  DO i = 1, npoints
     ! get the hw1e and e_asym at this point
     vgl =  locvgl(i);   vgr =  locvgr(i);  
     hwl = lochwl(i);    hwr =  lochwr(i); 

     !locsli(1:num_slit) = 0.0
     !slit(1:num_slit) = 0.0

     ! First do the right half of the slit function
     right_norm = voigt(0.0, vgr, hwr)
     DO j = mslit, num_slit
        j1 = i+(j-mslit) 
        IF ( j1 > npoints ) j1 = npoints - MOD(j1, npoints)
        !locsli(j) = specarr(j1)
        delwvl = (wvlarr(j1)-wvlarr(i))
        slit(j) = voigt(delwvl, vgr, hwr) / right_norm
        idx(j) = j1
     END DO
 
     ! Now the left half of the slit function
     left_norm = voigt(0.0, vgl, hwl)
     DO j = mslit-1, 1, -1
        j2 =  i+(j-mslit)
        IF ( j2 < 1 ) j2 = ABS(j2) + 2
        !locsli(j) = specarr(j2)
        delwvl = (wvlarr(j2)-wvlarr(i))
        slit(j) = voigt(delwvl, vgl, hwl) / left_norm
        idx(j) = j2
     END DO

     ! Need to normalize the slit function
     specmod(i) = DOT_PRODUCT(slit(1:num_slit), specarr(idx(1:num_slit))) &
          / SUM(slit(1:num_slit))
  END DO

  RETURN
END SUBROUTINE asym_voigt_vary

!xliu, 10/29/2009, add the capability to convole multiple spectrum simultaneously
SUBROUTINE asym_voigt_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)
  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : solwinfit, winlim, numwin, slit_trunc_limit
  USE OMSAO_indices_module,   ONLY : hwl_idx, hwr_idx, vgr_idx, vgl_idx
  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                        INTENT (IN)         :: nc, nf, nspec
  REAL (KIND=dp), DIMENSION (nf), INTENT (IN)         :: fwave
  REAL (KIND=dp), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
  REAL (KIND=dp), DIMENSION (nc), INTENT (IN)         :: cwave
  REAL (KIND=dp), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER        :: i, j, iwin, iw, fidx, fidxc, lidx, lidxc, midx, sidx, eidx
  REAL (KIND=dp) :: temp, hwl, hwr, vgr, vgl, vg, hw, lnorm, rnorm, voigt, ssum
  REAL (KIND=dp), DIMENSION (nf) :: slit
  
  !slit(:) = 0.0; cspec(:, :)= 0.0
  fidx = 1; fidxc = 1
  DO iwin = 1, numwin

     IF (iwin == numwin) THEN
        lidx = nf; lidxc = nc
     ELSE
        temp = (winlim(iwin, 2) + winlim(iwin + 1, 1)) / 2.0
        lidx =  MINVAL(MAXLOC(fwave, MASK=(fwave <= temp)))
        lidxc = MINVAL(MAXLOC(cwave, MASK=(cwave <= temp)))
     ENDIF    
     hwl = solwinfit(iwin, hwl_idx, 1)
     hwr = solwinfit(iwin, hwr_idx, 1)
     vgl = solwinfit(iwin, vgl_idx, 1)
     vgr = solwinfit(iwin, vgr_idx, 1)
     lnorm = voigt(0.0, vgl, hwl)
     rnorm = voigt(0.0, vgr, hwr)
     
     DO i = fidxc, lidxc
        ! Find the closest pixel
        midx = MINVAL(MAXLOC(fwave(fidx:lidx), MASK=(fwave(fidx:lidx) <= cwave(i)))) + fidx

        DO j = midx, fidx, -1
           slit(j) = voigt(fwave(j)-cwave(i), vgl, hwl) / lnorm
           IF (slit(j) <= slit_trunc_limit .OR. j == 1) EXIT
        ENDDO
        sidx = j

        DO j = midx, lidx
           slit(j) = voigt(fwave(j)-cwave(i), vgr, hwr) / rnorm
           IF (slit(j) <= slit_trunc_limit .OR. j == nf) EXIT
        ENDDO
        eidx = j
        
        ssum = SUM(slit(sidx:eidx))
        DO j = 1, nspec
           cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx)) / ssum
        ENDDO
     ENDDO
     
     fidx = lidx + 1; fidxc = lidxc + 1
  ENDDO

  RETURN

END SUBROUTINE asym_voigt_f2c

!xliu, 10/29/2009, add the capability to convole multiple spectrum simultaneously
SUBROUTINE asym_voigt_vary_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)
  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : winlim, numwin, slit_trunc_limit, slitwav, slitfit, nslit
  USE OMSAO_indices_module,   ONLY : hwl_idx, hwr_idx, vgr_idx, vgl_idx
  USE OMSAO_errstat_module
  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                        INTENT (IN)         :: nc, nf, nspec
  REAL (KIND=dp), DIMENSION (nf), INTENT (IN)         :: fwave
  REAL (KIND=dp), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
  REAL (KIND=dp), DIMENSION (nc), INTENT (IN)         :: cwave
  REAL (KIND=dp), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER        :: i, j, iwin, iw, fidx, fidxc, lidx, lidxc, &
       midx, sidx, eidx, fslit,lslit, finter, linter, errstat
  REAL (KIND=dp) :: temp, hwl, hwr, vgr, vgl, vg, hw, lnorm, rnorm, voigt, ssum
  REAL (KIND=dp), DIMENSION (nf) :: slit
  REAL (KIND=dp), DIMENSION (nc) :: lochwl, lochwr, locvgl, locvgr

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=18), PARAMETER :: modulename = 'asym_voigt_vary_f2c'

  errstat = pge_errstat_ok
  
  !slit(:) = 0.0; cspec(:, :)= 0.0
  fidx = 1; fidxc = 1
  DO iwin = 1, numwin
     IF (iwin == numwin) THEN
        lidx = nf; lidxc = nc
     ELSE
        temp = (winlim(iwin, 2) + winlim(iwin + 1, 1)) / 2.0
        lidx =  MINVAL(MAXLOC(fwave, MASK=(fwave <= temp)))
        lidxc = MINVAL(MAXLOC(cwave, MASK=(cwave <= temp)))
     ENDIF 
               
     fslit = MINVAL(MINLOC(slitwav(1:nslit), &
          MASK=(slitwav(1:nslit) >= cwave(fidxc))))
     lslit = MINVAL(MAXLOC(slitwav(1:nslit), &
          MASK=(slitwav(1:nslit) <= cwave(lidxc))))

     ! no slit between fidx:lidx, should never happen
     IF (fslit <= 0) THEN
        fslit = lslit
     ELSE IF (lslit <= 0) THEN
        lslit = fslit
     ENDIF

     IF (fslit > nslit .OR. lslit > nslit) THEN  
        WRITE(*, *) fslit, lslit, cwave(fidxc), cwave(lidxc), slitwav(1), &
             slitwav(nslit), nslit, cwave(1), cwave(nc)
        WRITE(*, *) modulename, ': Not slit available for this window!!!'
        STOP
     ENDIF

     IF (lslit < fslit + 3) THEN  ! extrapolate, use the nearest value
        lochwl(fidxc:lidxc) = SUM(slitfit(fslit:lslit, hwl_idx, 1))/(lslit-fslit+1)
        lochwr(fidxc:lidxc) = SUM(slitfit(fslit:lslit, hwr_idx, 1))/(lslit-fslit+1)
        locvgl(fidxc:lidxc) = SUM(slitfit(fslit:lslit, vgl_idx, 1))/(lslit-fslit+1)
        locvgr(fidxc:lidxc) = SUM(slitfit(fslit:lslit, vgr_idx, 1))/(lslit-fslit+1)
     ELSE
        ! finter <= linter here
        finter = MINVAL(MINLOC(cwave, MASK = (cwave >= slitwav(fslit))))
        linter = MINVAL(MAXLOC(cwave, MASK = (cwave <= slitwav(lslit))))
        
        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, hwl_idx, 1), &
             linter-finter+1, cwave(finter:linter), lochwl(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF

        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, hwr_idx, 1), &
             linter-finter+1, cwave(finter:linter), lochwr(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF

        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, vgl_idx, 1), &
             linter-finter+1, cwave(finter:linter), locvgl(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF

        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, vgr_idx, 1), &
             linter-finter+1, cwave(finter:linter), locvgr(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
        
        IF (finter > fidxc) THEN
           lochwl(fidxc:finter-1)=slitfit(fslit, hwl_idx, 1)
           lochwr(fidxc:finter-1)=slitfit(fslit, hwr_idx, 1)
           locvgl(fidxc:finter-1)=slitfit(fslit, vgl_idx, 1)
           locvgr(fidxc:finter-1)=slitfit(fslit, vgr_idx, 1)
        END IF
        
        IF (linter < lidxc)  THEN
           lochwl(linter+1:lidxc)=slitfit(lslit, hwl_idx, 1)
           lochwr(linter+1:lidxc)=slitfit(lslit, hwr_idx, 1)
           locvgl(linter+1:lidxc)=slitfit(lslit, vgl_idx, 1)
           locvgr(linter+1:lidxc)=slitfit(lslit, vgr_idx, 1)
        END IF
     ENDIF
   
     DO i = fidxc, lidxc

        hwl = lochwl(i); hwr = lochwr(i)
        vgl = locvgl(i); vgr = locvgr(i)
        lnorm = voigt(0.0, vgl, hwl)
        rnorm = voigt(0.0, vgr, hwr)

        ! Find the closest pixel
        midx = MINVAL(MAXLOC(fwave(fidx:lidx), MASK=(fwave(fidx:lidx) <= cwave(i)))) + fidx

        DO j = midx, fidx, -1
           slit(j) = voigt(fwave(j)-cwave(i), vgl, hwl) / lnorm
           IF (slit(j) <= slit_trunc_limit .OR. j == 1) EXIT
        ENDDO
        sidx = j

        DO j = midx, lidx
           slit(j) = voigt(fwave(j)-cwave(i), vgr, hwr) / rnorm
           IF (slit(j) <= slit_trunc_limit .OR. j == nf) EXIT
        ENDDO
        eidx = j
        
        ssum = SUM(slit(sidx:eidx))
        DO j = 1, nspec
           cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx)) / ssum
        ENDDO
     ENDDO
     
     fidx = lidx + 1; fidxc = lidxc + 1
  ENDDO

  RETURN

END SUBROUTINE asym_voigt_vary_f2c

