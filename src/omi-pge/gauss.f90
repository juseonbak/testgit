!  Implement both symmetric gaussian and asymmetric gaussian slit function
!  If only symmetric gaussian is used, then use symmetric gaussian (faster)

! =========================================================================
!
! Convolves input spectrum with an asymmetric Gaussian slit function of
! specified HW1E (half-width at 1/e intensity) and asymmetry factor E_ASYM.
!
! The asymetric Gaussian g(x) is defined as
!                   _                                  _
!                  |               x^2                  |
!      g(x) =  EXP | - -------------------------------- |
!                  |_   (hw1e * (1 + SIGN(x)*e_aym))^2 _|
!
! g(x) bmwinecomes symmetric for E_ASYM = 0.
!
! =========================================================================

SUBROUTINE asym_gauss (wvlarr, specarr, specmod, npoints, hw1e, e_asym)


  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : slit_trunc_limit
  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                             INTENT (IN) :: npoints
  REAL (KIND=dp),                      INTENT (IN) :: hw1e, e_asym
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN) :: wvlarr, specarr

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=dp), DIMENSION (npoints), INTENT (OUT) :: specmod

  ! ===============
  ! Local variables
  ! ===============
  INTEGER         :: i, j, j1, j2, num_slit, mslit
  REAL (KIND=dp)  :: delwvl, slitsum,  maxslit, rsw, lsw
  REAL (KIND=dp), DIMENSION (npoints) :: slit, locwvl
  INTEGER,        DIMENSION (npoints) :: idx

  ! --------------------------------------------------------
  ! Initialize output variable (default for "no convolution"
  ! --------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  
  ! -----------------------------------------------
  ! No Gaussian convolution if Halfwidth @ 1/e is 0
  ! -----------------------------------------------
  IF ( hw1e == 0.0 ) RETURN

  ! --------------------------------------------------------------
  ! Find the number of spectral points that fall within a Gaussian
  ! slit function with values >= 0.001. Remember that we have an
  ! asymmetric Gaussian, so we create a wavelength array symmetric
  ! around 0. The spacing is provided by the equidistant WVLARR.
  ! --------------------------------------------------------------
  rsw = -(hw1e * (1.0 + ABS(e_asym))) ** 2.0
  delwvl = wvlarr(2) - wvlarr(1)
  mslit = NINT( SQRT( LOG(slit_trunc_limit) * rsw) / delwvl) 
  num_slit = mslit * 2 + 1; mslit = mslit + 1

  IF (num_slit > npoints) THEN
     mslit = (npoints-1) / 2 + 1;    num_slit = mslit * 2 - 1
     !WRITE(*, *) 'Too small a slit fitting window!!!'; RETURN
  ENDIF

  ! left side
  lsw = -(hw1e * (1.0 - e_asym)) ** 2.0
  j = -mslit + 1
  DO i = 1, mslit - 1
     locwvl(i) = delwvl * j
     j = j + 1
  ENDDO
  slit(1:mslit-1) = EXP(locwvl(1:mslit-1)**2 / lsw)
  locwvl(mslit) = 0.0; slit(mslit) = 1.0

  ! Right side
  rsw = -(hw1e * (1.0 + e_asym)) ** 2.0
  j = 1
  DO i = mslit+1, num_slit
     locwvl(i) = delwvl * j
     j = j + 1
  ENDDO
  slit(mslit+1:num_slit) = EXP(locwvl(mslit+1:num_slit)**2 / rsw)

  ! Normalization
  slitsum = SUM(slit(1:num_slit))
  slit(1:num_slit) = slit(1:num_slit) / slitsum

  ! ---------------------------------------------------------------
  ! Convolve spectrum. First do the middle part, where we have full
  ! overlap coverage of the slit function. Again, remember the
  ! asymmetry of the Gaussian, which makes impossible a simple
  ! 50-50 division of the summation interval.
  ! ---------------------------------------------------------------

  ! Make a local copy of the NSLIT spectrum points to be convolved
  ! with the slit function. The spectrum points to be convolved are
  ! arranged such that the updated index corresponds to the maximum
  ! of the slit function (MSLIT). For simplicity we reflect the
  ! spectrum at the array end points.

  ! ----------------------------------------------------
  ! Loop over all points of the spectrum to be convolved
  ! ----------------------------------------------------
  DO i = 1, npoints
     ! First do the right half of the slit function
     DO j = mslit, num_slit
        j1 = i+j-mslit ; IF ( j1 > npoints ) j1 = npoints - MOD(j1, npoints)
        idx(j) = j1
     END DO

     ! Now the left half of the slit function
     DO j = mslit-1, 1, -1
        j2 = i+j-mslit ; IF ( j2 < 1 ) j2 = ABS(j2) + 2
        idx(j) = j2
     END DO

     specmod(i) = DOT_PRODUCT(slit(1:num_slit), specarr(idx(1:num_slit)))
  END DO

  RETURN
END SUBROUTINE asym_gauss

SUBROUTINE asym_gauss_multi (wvlarr, specarr, specmod, npoints)

  USE OMSAO_precision_module
  USE OMSAO_variables_module,  ONLY : solwinfit, winlim, numwin, slit_trunc_limit
  USE OMSAO_indices_module,    ONLY : hwe_idx, asy_idx

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
  INTEGER                             :: i, j, j1, j2, num_slit, mslit, fidx, lidx, iwin
  REAL (KIND=dp)                      :: delwvl, slitsum,  maxslit, hw1e, e_asym, rsw, lsw
  REAL (KIND=dp), DIMENSION (npoints) :: slit, locwvl, upbnd
  INTEGER,        DIMENSION (npoints) :: idx
  REAL (KIND=dp), EXTERNAL            :: signdp

  ! --------------------------------------------------------
  ! Initialize output variable (default for "no convolution"
  ! --------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  
  fidx = 1
  DO iwin = 1, numwin

     hw1e = solwinfit(iwin, hwe_idx, 1); e_asym = solwinfit(iwin, asy_idx, 1)
     
     IF (iwin < numwin) THEN
        upbnd = (winlim(iwin, 2) + winlim(iwin+1, 1))/2.0 
        lidx = MINVAL(MAXLOC(wvlarr, MASK=(wvlarr <= upbnd )))
     ELSE 
        lidx = npoints
     END IF

     ! -----------------------------------------------
     ! No Gaussian convolution if Halfwidth @ 1/e is 0
     ! -----------------------------------------------
     IF ( hw1e == 0.0 ) RETURN

     ! --------------------------------------------------------------
     ! Find the number of spectral points that fall within a Gaussian
     ! slit function with values >= 0.001. Remember that we have an
     ! asymmetric Gaussian, so we create a wavelength array symmetric
     ! around 0. The spacing is provided by the equidistant WVLARR.
     ! --------------------------------------------------------------
     rsw = -(hw1e * (1.0 + ABS(e_asym))) ** 2.0
     delwvl = wvlarr(fidx+1) - wvlarr(fidx)
     mslit = NINT( SQRT( LOG(slit_trunc_limit) * rsw) / delwvl) 
     num_slit = mslit * 2 + 1; mslit = mslit + 1
     
     IF (num_slit > npoints) THEN
        mslit = (npoints-1) / 2 + 1;    num_slit = mslit * 2 - 1
     ENDIF
     
     ! left side
     lsw = -(hw1e * (1.0 - e_asym)) ** 2.0
     j = -mslit + 1
     DO i = 1, mslit - 1
        locwvl(i) = delwvl * j
        j = j + 1
     ENDDO
     slit(1:mslit-1) = EXP(locwvl(1:mslit-1)**2 / lsw)
     locwvl(mslit) = 0.0; slit(mslit) = 1.0

     ! Right side
     rsw = -(hw1e * (1.0 + e_asym)) ** 2.0
     j = 1
     DO i = mslit+1, num_slit
        locwvl(i) = delwvl * j
        j = j + 1
     ENDDO
     slit(mslit+1:num_slit) = EXP(locwvl(mslit+1:num_slit)**2 / rsw)
     
     ! Normalization
     slitsum = SUM(slit(1:num_slit))
     slit(1:num_slit) = slit(1:num_slit) / slitsum
     
     ! ---------------------------------------------------------------
     ! Convolve spectrum. First do the middle part, where we have full
     ! overlap coverage of the slit function. Again, remember the
     ! asymmetry of the Gaussian, which makes impossible a simple
     ! 50-50 division of the summation interval.
     ! ---------------------------------------------------------------
     
     ! Make a local copy of the NSLIT spectrum points to be convolved
     ! with the slit function. The spectrum points to be convolved are
     ! arranged such that the updated index corresponds to the maximum
     ! of the slit function (MSLIT). For simplicity we reflect the
     ! spectrum at the array end points.
     
     ! ----------------------------------------------------
     ! Loop over all points of the spectrum to be convolved
     ! ----------------------------------------------------
     DO i = fidx, lidx
        ! First do the right half of the slit function
        DO j = mslit, num_slit
           j1 = i + j - mslit
           IF ( j1 > npoints ) j1 = npoints - MOD(j1, npoints)
           idx(j) = j1
        END DO

        ! Now the left half of the slit function
        DO j = mslit-1, 1, -1
           j2 = i+ j - mslit ; IF ( j2 < 1 ) j2 = ABS(j2) + 2          
           idx(j) = j2
        END DO
       
        specmod(i) = DOT_PRODUCT(slit(1:num_slit), specarr(idx(1:num_slit)))
     END DO

     fidx = lidx + 1
  END DO

  RETURN
END SUBROUTINE asym_gauss_multi


SUBROUTINE asym_gauss_vary (wvlarr, specarr, specmod, npoints)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : slitwav, slitfit, nslit, winlim, numwin, slit_trunc_limit
  USE OMSAO_indices_module,  ONLY  : hwe_idx, asy_idx
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
  INTEGER                             :: i, j, j1, j2, num_slit, mslit, &
       errstat,  fidx, lidx, fslit, lslit, finter, linter, iwin
  REAL (KIND=dp)                      :: delwvl, slitsum,  maxslit, hw1e, e_asym, lsw, rsw, upbnd
  REAL (KIND=dp), DIMENSION (npoints) :: slit, locwvl, lochwe, locasy
  INTEGER,        DIMENSION (npoints) :: idx

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=15), PARAMETER :: modulename = 'asym_gauss_vary'


  !WRITE(*, *) 'nslit = ', nslit
  !WRITE(*, '(10f8.3)') slitwav(1:nslit), slitfit(1:nslit, hwe_idx, 1)

  ! --------------------------------------------------------
  ! Initialize output variable (default for "no convolution"
  ! --------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  

  ! -----------------------------------------------
  ! No Gaussian convolution if Halfwidth @ 1/e is 0
  ! -----------------------------------------------
  hw1e   = MAXVAL(slitfit(1:nslit, hwe_idx, 1))
  e_asym = MAXVAL(ABS(slitfit(1:nslit, asy_idx, 1)))

  IF ( hw1e == 0.0 ) RETURN

  ! --------------------------------------------------------------
  ! Find the number of spectral points that fall within a Gaussian
  ! slit function with values >= 0.01. Remember that we have an
  ! asymmetric Gaussian, so we create a wavelength array symmetric
  ! around 0. The spacing is provided by the equidistant WVLARR.
  ! --------------------------------------------------------------  
  rsw = -(hw1e * (1.0 + ABS(e_asym))) ** 2.0
  delwvl = wvlarr(2) - wvlarr(1)
  mslit = NINT( SQRT( LOG(slit_trunc_limit) * rsw) / delwvl) 
  num_slit = mslit * 2 + 1; mslit = mslit + 1

  IF (num_slit > npoints) THEN
     mslit = (npoints-1) / 2 + 1;    num_slit = mslit * 2 - 1
     !WRITE(*, *) 'Too small a slit fitting window!!!'; RETURN
  ENDIF
 
  ! ---------------------------------------------------------------
  ! Convolve spectrum. First do the middle part, where we have full
  ! overlap coverage of the slit function. Again, remember the
  ! asymmetry of the Gaussian, which makes impossible a simple
  ! 50-50 division of the summation interval.
  ! ---------------------------------------------------------------

  ! Make a local copy of the NSLIT spectrum points to be convolved
  ! with the slit function. The spectrum points to be convolved are
  ! arranged such that the updated index corresponds to the maximum
  ! of the slit function (MSLIT). For simplicity we reflect the
  ! spectrum at the array end points.

  errstat = pge_errstat_ok

  ! get an array of slit variables (lochwe, locasy) for each wavelength position
  lochwe = 0.0; locasy=0.0

  fidx = 1; lidx = npoints
  DO iwin = 1, numwin
     IF (iwin < numwin) THEN
        upbnd = (winlim(iwin, 2) + winlim(iwin+1, 1))/2.0 
        lidx = MINVAL(MAXLOC(wvlarr, MASK=(wvlarr <= upbnd )))
     ELSE 
        lidx = npoints
     END IF
               
     fslit = MINVAL(MINLOC(slitwav(1:nslit), &
          MASK=(slitwav(1:nslit) >= wvlarr(fidx))))
     lslit = MINVAL(MAXLOC(slitwav(1:nslit), &
          MASK=(slitwav(1:nslit) <= wvlarr(lidx))))
 
     ! no slit between fidx:lidx, should never happen
     IF (fslit <= 0) THEN
        fslit = lslit
     ELSE IF (lslit <= 0) THEN
        lslit = fslit
     ENDIF

     IF (fslit > nslit .OR. lslit > nslit) THEN  
        WRITE(*, *) fslit, lslit, wvlarr(fidx), wvlarr(lidx), slitwav(1), &
             slitwav(nslit), nslit, wvlarr(1), wvlarr(npoints)
        WRITE(*, *) modulename, ': Not slit available for this window!!!'
        STOP
     ENDIF
     
     IF (lslit < fslit + 3) THEN  ! extrapolate, use the nearest value
        lochwe(fidx:lidx) = SUM(slitfit(fslit:lslit, hwe_idx, 1))/(lslit-fslit+1)
        locasy(fidx:lidx) = SUM(slitfit(fslit:lslit, asy_idx, 1))/(lslit-fslit+1) 
     ELSE
        ! finter <= linter here
        finter = MINVAL(MINLOC(wvlarr, MASK = (wvlarr >= slitwav(fslit))))
        linter = MINVAL(MAXLOC(wvlarr, MASK = (wvlarr <= slitwav(lslit))))
        
        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, hwe_idx, 1), &
             linter-finter+1, wvlarr(finter:linter), lochwe(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
        
        CALL interpolation (lslit-fslit+1, slitwav(fslit:lslit), &
             slitfit(fslit:lslit, asy_idx, 1), linter-finter+1, wvlarr(finter:linter),&
             locasy(finter:linter), errstat)
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
        
        IF (finter > fidx) THEN
           lochwe(fidx:finter-1)=slitfit(fslit, hwe_idx, 1)
           locasy(fidx:finter-1)=slitfit(fslit, asy_idx, 1)
        END IF
        
        IF (linter < lidx)  THEN
           lochwe(linter+1:lidx)=slitfit(lslit, hwe_idx, 1)
           locasy(linter+1:lidx)=slitfit(lslit, asy_idx, 1)
        END IF
     ENDIF
     fidx = lidx + 1
  ENDDO
    
  DO i = 1, npoints
     ! get the hw1e and e_asym at this point
     hw1e = lochwe(i); e_asym = locasy(i)
     lsw = -(hw1e*(1.0 - e_asym))**2.0
     rsw = -(hw1e*(1.0 + e_asym))**2.0
     
     ! First do the right half of the slit function
     DO j = mslit, num_slit
        j1 = i + j - mslit
        IF (j1 > npoints) j1 = npoints - MOD(j1, npoints)
        locwvl(j) = (wvlarr(j1)-wvlarr(i)) ** 2 / rsw
        idx(j) = j1
     END DO
        
     ! Now the left half of the slit function
     DO j = mslit-1, 1, -1
        j2 = i + j - mslit
        IF (j2 < 1) j2 = ABS(j2) + 2
        locwvl(j) = (wvlarr(j2)-wvlarr(i)) ** 2 / lsw
        idx(j) = j2
     END DO
     slit(1:num_slit) = EXP(locwvl(1:num_slit))
         
     specmod(i) = DOT_PRODUCT(slit(1:num_slit), specarr(idx(1:num_slit))) &
          / SUM(slit(1:num_slit))
  END DO

  RETURN
END SUBROUTINE asym_gauss_vary


! =========================================================================
!
! Convolves input spectrum with an asymmetric Gaussian slit function of
! specified HW1E (half-width at 1/e intensity)
!
! The :q!symetric Gaussian g(x) is defined as
!                   _                _
!                  |         x^2      |
!      g(x) =  EXP | - -------------- |
!                  |_     hw1e ^2    _|
!  FWHM = 2.0 * sqrt(ln(2.0)) * hw1e = 1.66551 * hw1e
! =========================================================================

SUBROUTINE gauss (wvlarr, specarr, specmod, npoints, hw1e)
  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : slit_trunc_limit
  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                             INTENT (IN)    :: npoints
  REAL (KIND=dp),                      INTENT (IN)    :: hw1e
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN)    :: wvlarr, specarr

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=dp), DIMENSION (npoints), INTENT (OUT) :: specmod

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                             :: nhi, nlo, i, j, num_slit
  REAL (KIND=dp)                      :: emult, delwvl, slitsum, slit0
  REAL (KIND=dp), DIMENSION (npoints) :: slit

  ! ----------------------------------------------------------------
  ! Initialization of output variable (default for "no convolution")
  ! ----------------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  
  ! --------------------------------------
  ! No convolution if halfwidth @ 1/e is 0
  ! --------------------------------------
  IF ( hw1e == 0.0 ) RETURN

  emult  = -1.0 / ( hw1e * hw1e )
  delwvl = wvlarr(2) - wvlarr(1)
 
  !  Calculate slit function values out to 0.001 times x0 value,
  !     normalize so that sum = 1.

  slitsum = 1.0 ;  slit0 = 1.0
  i = 1  ;  num_slit = 0
  DO WHILE ( num_slit <= npoints )
     slit (i) = EXP (emult * (delwvl * i)**2)
     slitsum = slitsum + 2.0 * slit (i)
     IF (slit (i) <= slit_trunc_limit ) EXIT 
     i = i + 1
  ENDDO
  num_slit = i
  
  slit0 = slit0 / slitsum
  slit(1:num_slit) = slit(1:num_slit) / slitsum

  ! Convolve spectrum.  reflect at endpoints.
  ! Doesn't look right
  specmod(1:npoints) = slit0 * specarr(1:npoints)
  DO i = 1, npoints
     DO j = 1, num_slit
        nlo = i - j 
        IF (nlo < 1) nlo = -nlo + 2 
        nhi = i + j
        IF ( nhi > npoints ) nhi = npoints - MOD(nhi, npoints)
        specmod(i) = specmod(i) + slit(j) * ( specarr(nlo) + specarr(nhi) )
     END DO
  END DO

  RETURN
END SUBROUTINE gauss

SUBROUTINE gauss_multi (wvlarr, specarr, specmod, npoints)

  USE OMSAO_precision_module
  USE OMSAO_variables_module,  ONLY : solwinfit, winlim, numwin, slit_trunc_limit
  USE OMSAO_indices_module,    ONLY : hwe_idx

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
  INTEGER                             :: i, j, ii, j1, j2, num_slit, mslit, fidx, lidx, iwin
  REAL (KIND=dp)                      :: delwvl, slitsum, hw1e, sw
  REAL (KIND=dp), DIMENSION (npoints) :: slit, locwvl, upbnd

  ! --------------------------------------------------------
  ! Initialize output variable (default for "no convolution"
  ! --------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  
  fidx = 1
  DO iwin = 1, numwin

     hw1e = solwinfit(iwin, hwe_idx, 1)
     
     IF (iwin < numwin) THEN
        upbnd = (winlim(iwin, 2) + winlim(iwin+1, 1))/2.0 
        lidx = MINVAL(MAXLOC(wvlarr, MASK=(wvlarr <= upbnd )))
     ELSE 
        lidx = npoints
     END IF

     ! -----------------------------------------------
     ! No Gaussian convolution if Halfwidth @ 1/e is 0
     ! -----------------------------------------------
     IF ( hw1e == 0.0 ) RETURN

     ! --------------------------------------------------------------
     ! Find the number of spectral points that fall within a Gaussian
     ! slit function with values >= 0.001. Remember that we have an
     ! asymmetric Gaussian, so we create a wavelength array symmetric
     ! around 0. The spacing is provided by the equidistant WVLARR.
     ! --------------------------------------------------------------
     sw = - hw1e ** 2.0
     delwvl = wvlarr(fidx+1) - wvlarr(fidx)
     mslit = NINT( SQRT( LOG(slit_trunc_limit) * sw) / delwvl) 
     num_slit = mslit * 2 + 1
     
     IF (num_slit > npoints) THEN
        mslit = (npoints-1) / 2 ;    num_slit = mslit * 2 + 1
     ENDIF
     
     ! only compute half of the slit shape (say left side)
     j = -mslit
     DO i = 1, mslit
        locwvl(i) = delwvl * j
        j = j + 1
     ENDDO
     slit(1:mslit) = EXP(locwvl(1:mslit)**2 / sw)

     ! Normalization
     slitsum = SUM(slit(1:mslit)) * 2.0 + 1.0   ! (center value is 1.0)
     slit(1:mslit) = slit(1:mslit) / slitsum
          
     ! Make a local copy of the NSLIT spectrum points to be convolved
     ! with the slit function. The spectrum points to be convolved are
     ! arranged such that the updated index corresponds to the maximum
     ! of the slit function (MSLIT). For simplicity we reflect the
     ! spectrum at the array end points.
     
     ! ----------------------------------------------------
     ! Loop over all points of the spectrum to be convolved
     ! ----------------------------------------------------
     DO i = fidx, lidx
        ! Center value
        specmod(i) = specarr(i) / slitsum  ! (i.e., 1.0 / slitsum for center)

        ! Add both left and side (symmetric)
        DO j = 1, mslit
           ii = mslit + 1 - j
           j1 = i + ii; IF ( j1 > npoints ) j1 = npoints - MOD(j1, npoints)
           j2 = i - ii; IF ( j2 < 1 ) j2 = ABS(j2) + 2 
           specmod(i) = specmod(i) + slit(j) * ( specarr(j1) + specarr(j2))
        ENDDO
     END DO

     fidx = lidx + 1
  END DO

  RETURN
END SUBROUTINE gauss_multi

SUBROUTINE gauss_vary (wvlarr, specarr, specmod, npoints)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : slitwav, slitfit, nslit, winlim, numwin, slit_trunc_limit
  USE OMSAO_indices_module,  ONLY  : hwe_idx
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
  INTEGER                             :: i, ii, j, j1, j2, num_slit, mslit, &
       errstat,  fidx, lidx, fslit, lslit, finter, linter, iwin
  REAL (KIND=dp)                      :: delwvl, slitsum,  hw1e, sw, upbnd
  REAL (KIND=dp), DIMENSION (npoints) :: slit, lochwe

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=10), PARAMETER :: modulename = 'gauss_vary'


  !WRITE(*, *) 'nslit = ', nslit
  !WRITE(*, '(10f8.3)') slitwav(1:nslit), slitfit(1:nslit, hwe_idx, 1)

  ! --------------------------------------------------------
  ! Initialize output variable (default for "no convolution"
  ! --------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  

  ! -----------------------------------------------
  ! No Gaussian convolution if Halfwidth @ 1/e is 0
  ! -----------------------------------------------
  hw1e   = MAXVAL(slitfit(1:nslit, hwe_idx, 1))

  IF ( hw1e == 0.0 ) RETURN

  ! --------------------------------------------------------------
  ! Find the number of spectral points that fall within a Gaussian
  ! slit function with values >= 0.01. Remember that we have an
  ! asymmetric Gaussian, so we create a wavelength array symmetric
  ! around 0. The spacing is provided by the equidistant WVLARR.
  ! --------------------------------------------------------------  
  sw = -hw1e ** 2.0
  delwvl = wvlarr(2) - wvlarr(1)
  mslit = NINT( SQRT( LOG(slit_trunc_limit) * sw) / delwvl) 
  num_slit = mslit * 2 + 1

  IF (num_slit > npoints) THEN
     mslit = (npoints-1) / 2 + 1
  ENDIF
 
  ! ---------------------------------------------------------------
  ! Convolve spectrum. First do the middle part, where we have full
  ! overlap coverage of the slit function. Again, remember the
  ! asymmetry of the Gaussian, which makes impossible a simple
  ! 50-50 division of the summation interval.
  ! ---------------------------------------------------------------

  ! Make a local copy of the NSLIT spectrum points to be convolved
  ! with the slit function. The spectrum points to be convolved are
  ! arranged such that the updated index corresponds to the maximum
  ! of the slit function (MSLIT). For simplicity we reflect the
  ! spectrum at the array end points.

  errstat = pge_errstat_ok

  ! get an array of slit variables lochwe for each wavelength position
  lochwe = 0.0

  fidx = 1; lidx = npoints
  DO iwin = 1, numwin
     IF (iwin < numwin) THEN
        upbnd = (winlim(iwin, 2) + winlim(iwin+1, 1))/2.0 
        lidx = MINVAL(MAXLOC(wvlarr, MASK=(wvlarr <= upbnd )))
     ELSE 
        lidx = npoints
     END IF
               
     fslit = MINVAL(MINLOC(slitwav(1:nslit), &
          MASK=(slitwav(1:nslit) >= wvlarr(fidx))))
     lslit = MINVAL(MAXLOC(slitwav(1:nslit), &
          MASK=(slitwav(1:nslit) <= wvlarr(lidx))))
 
     ! no slit between fidx:lidx, should never happen
     IF (fslit <= 0) THEN
        fslit = lslit
     ELSE IF (lslit <= 0) THEN
        lslit = fslit
     ENDIF

     IF (fslit > nslit .OR. lslit > nslit) THEN  
        WRITE(*, *) fslit, lslit, wvlarr(fidx), wvlarr(lidx), slitwav(1), &
             slitwav(nslit), nslit, wvlarr(1), wvlarr(npoints)
        WRITE(*, *) modulename, ': Not slit available for this window!!!'
        STOP
     ENDIF
     
     IF (lslit < fslit + 3) THEN  ! extrapolate, use the nearest value
        lochwe(fidx:lidx) = SUM(slitfit(fslit:lslit, hwe_idx, 1))/(lslit-fslit+1)
     ELSE
        ! finter <= linter here
        finter = MINVAL(MINLOC(wvlarr, MASK = (wvlarr >= slitwav(fslit))))
        linter = MINVAL(MAXLOC(wvlarr, MASK = (wvlarr <= slitwav(lslit))))
        
        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, hwe_idx, 1), &
             linter-finter+1, wvlarr(finter:linter), lochwe(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
        
        IF (finter > fidx) THEN
           lochwe(fidx:finter-1)=slitfit(fslit, hwe_idx, 1)
        END IF
        
        IF (linter < lidx)  THEN
           lochwe(linter+1:lidx)=slitfit(lslit, hwe_idx, 1)
        END IF
     ENDIF
     fidx = lidx + 1
  ENDDO
    
  DO i = 1, npoints
     ! Get center value
     specmod(i) = specarr(i); slitsum = 1.0

     ! Compute slit contribution from both sides
     sw   = -lochwe(i) ** 2.0
     DO j = 1, mslit
        ii = mslit + 1 - j
        j1 = i + ii; IF ( j1 > npoints ) j1 = npoints - MOD(j1, npoints)
        j2 = i - ii; IF ( j2 < 1 ) j2 = ABS(j2) + 2 
        slit(j) = EXP( (wvlarr(j1)-wvlarr(i)) ** 2 / sw)
        slitsum = slitsum + 2.0 * slit(j)
        specmod(i) = specmod(i) + slit(j) * ( specarr(j1) + specarr(j2))
     ENDDO
     
     specmod(i) = specmod(i) / slitsum
  END DO

  RETURN
END SUBROUTINE gauss_vary

! =========================================================================
!
! Convolves input spectrum with an super Gaussian slit function 
! where the exponent of the Gaussian is a variable N and will give a Gaussian for N=2. 
! For big N the function will describe a more rectangular distribution, 
! while for small N it fits to a distribution with long tails on both sides. 
!
! The symetric super Gaussian g(x) is defined as
!                               _               _
!              K               |      abs(x)^k   |
!      g(x) =-----         EXP | - ------------- |
!            2hw1e gam(1/k)    |_     hw1e ^k   _|
!  FWHM = 2.0 * sqrt(ln(2.0)) * hw1e = 1.66551 * hw1e
! =========================================================================
SUBROUTINE super_gauss (wvlarr, specarr, specmod, npoints, hw1e, power)
  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : slit_trunc_limit
  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                             INTENT (IN)    :: npoints
  REAL (KIND=dp),                      INTENT (IN)    :: hw1e, power
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN)    :: wvlarr, specarr

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=dp), DIMENSION (npoints), INTENT (OUT)   :: specmod

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                             :: nhi, nlo, i, j, num_slit
  REAL (KIND=dp)                      :: emult, delwvl, slitsum, slit0, coef
  REAL (KIND=dp), DIMENSION (npoints) :: slit

  ! ----------------------------------------------------------------
  ! Initialization of output variable (default for "no convolution")
  ! ----------------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  
  ! --------------------------------------
  ! No convolution if halfwidth @ 1/e is 0
  ! --------------------------------------
  IF ( hw1e == 0.0 .or. power ==0.0) STOP

  
  delwvl = wvlarr(2) - wvlarr(1)

  !  Calculate slit function values out to 0.001 times x0 value,
  !     normalize so that sum = 1.
  coef   = power / ( 2*hw1e  *gamma(power)) 
  coef   = 1
  emult  = -1.0/(hw1e)**power
  slitsum = coef ! this initial value equal to center of slit
  slit0   = coef
  i = 1  ;  num_slit = 0
  DO WHILE ( num_slit <= npoints )
     slit (i) = coef*EXP (emult*(abs(delwvl*i))**power)
     slitsum = slitsum + 2.0 * slit (i)
     IF (slit (i) <= slit_trunc_limit ) EXIT 
     i = i + 1
  ENDDO
  num_slit = i

  slit0 = slit0 / slitsum
  slit(1:num_slit) = slit(1:num_slit) / slitsum
 
  ! Convolve spectrum. reflect at endpoints.
  ! Doesn't look right
  specmod(1:npoints) = slit0 * specarr(1:npoints)
  DO i = 1, npoints
     DO j = 1, num_slit
        nlo = i - j 
        IF (nlo < 1) nlo = -nlo + 2 
        nhi = i + j
        IF ( nhi > npoints ) nhi = npoints - MOD(nhi, npoints)
        specmod(i) = specmod(i) + slit(j) * ( specarr(nlo) + specarr(nhi) )
     END DO
  END DO

  RETURN
END SUBROUTINE super_gauss


SUBROUTINE super_gauss_multi (wvlarr, specarr, specmod, npoints)

  USE OMSAO_precision_module
  USE OMSAO_variables_module,  ONLY : solwinfit, winlim, numwin, slit_trunc_limit
  USE OMSAO_indices_module,    ONLY : hwe_idx, hwn_idx

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
  INTEGER                             :: i, j, ii, j1, j2, num_slit, mslit, fidx, lidx, iwin
  REAL (KIND=dp)                      :: delwvl, slitsum, hw1e, power, sw, coef
  REAL (KIND=dp), DIMENSION (npoints) :: slit, locwvl, upbnd

  ! --------------------------------------------------------
  ! Initialize output variable (default for "no convolution"
  ! --------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  
  fidx = 1
  DO iwin = 1, numwin

     hw1e  = solwinfit(iwin, hwe_idx, 1)
     power = solwinfit(iwin, hwn_idx, 1)

     IF (iwin < numwin) THEN
        upbnd = (winlim(iwin, 2) + winlim(iwin+1, 1))/2.0 
        lidx = MINVAL(MAXLOC(wvlarr, MASK=(wvlarr <= upbnd )))
     ELSE 
        lidx = npoints
     END IF

     ! -----------------------------------------------
     ! No Gaussian convolution if Halfwidth @ 1/e is 0
     ! -----------------------------------------------
     IF ( hw1e == 0.0 .or. power == 0.0 ) STOP

     ! --------------------------------------------------------------
     ! Find the number of spectral points that fall within a Gaussian
     ! slit function with values >= 0.001. Remember that we have an
     ! asymmetric Gaussian, so we create a wavelength array symmetric
     ! around 0. The spacing is provided by the equidistant WVLARR.
     ! --------------------------------------------------------------
     sw = - hw1e ** 2.0
     delwvl = wvlarr(fidx+1) - wvlarr(fidx)
     mslit = NINT( SQRT( LOG(slit_trunc_limit) * sw) / delwvl) 
     num_slit = mslit * 2 + 1
     
     IF (num_slit > npoints) THEN
        mslit = (npoints-1) / 2 ;    num_slit = mslit * 2 + 1
     ENDIF
         
     ! only compute half of the slit shape (say left side)
     j = -mslit
     DO i = 1, mslit
        locwvl(i) = delwvl * j
        j = j + 1
     ENDDO
     coef   = power / ( 2*hw1e  *gamma(power))
     coef   = 1.0
     slit (1:mslit) = coef*EXP (-1.*(abs(locwvl(1:mslit))/hw1e)**power)


     ! Normalization
     slitsum = SUM(slit(1:mslit)) * 2.0 + coef   ! (center value is 1.0)
     slit(1:mslit) = slit(1:mslit) / slitsum
          
     ! Make a local copy of the NSLIT spectrum points to be convolved
     ! with the slit function. The spectrum points to be convolved are
     ! arranged such that the updated index corresponds to the maximum
     ! of the slit function (MSLIT). For simplicity we reflect the
     ! spectrum at the array end points.
      
     ! ----------------------------------------------------
     ! Loop over all points of the spectrum to be convolved
     ! ----------------------------------------------------
     DO i = fidx, lidx
        ! Center value
        specmod(i) = specarr(i)* coef/ slitsum  ! (i.e., 1.0 / slitsum for center)

        ! Add both left and side (symmetric)
        DO j = 1, mslit
           ii = mslit + 1 - j
           j1 = i + ii; IF ( j1 > npoints ) j1 = npoints - MOD(j1, npoints)
           j2 = i - ii; IF ( j2 < 1 ) j2 = ABS(j2) + 2 
           specmod(i) = specmod(i) + slit(j) * ( specarr(j1) + specarr(j2))
        ENDDO
     END DO

     fidx = lidx + 1
  END DO

  RETURN
END SUBROUTINE super_gauss_multi

SUBROUTINE super_gauss_vary (wvlarr, specarr, specmod, npoints) ! need to check in detail

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : slitwav, slitfit, nslit, winlim, numwin, slit_trunc_limit
  USE OMSAO_indices_module,  ONLY  : hwe_idx,hwn_idx
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
  INTEGER                             :: i, ii, j, j1, j2, num_slit, mslit, &
       errstat,  fidx, lidx, fslit, lslit, finter, linter, iwin
  REAL (KIND=dp)                      :: delwvl, slitsum,  hw1e, sw, upbnd, POWER, coeff
  REAL (KIND=dp), DIMENSION (npoints) :: slit, lochwe, lochwn

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=10), PARAMETER :: modulename = 'gauss_vary'


  !WRITE(*, *) 'nslit = ', nslit
  !WRITE(*, '(10f8.3)') slitwav(1:nslit), slitfit(1:nslit, hwe_idx, 1)

  ! --------------------------------------------------------
  ! Initialize output variable (default for "no convolution"
  ! --------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  

  ! -----------------------------------------------
  ! No Gaussian convolution if Halfwidth @ 1/e is 0
  ! -----------------------------------------------
  hw1e   = MAXVAL(slitfit(1:nslit, hwe_idx, 1))
  power   = MAXVAL(slitfit(1:nslit, hwn_idx, 1))

  IF ( hw1e == 0.0 ) RETURN

  ! --------------------------------------------------------------
  ! Find the number of spectral points that fall within a Gaussian
  ! slit function with values >= 0.01. Remember that we have an
  ! asymmetric Gaussian, so we create a wavelength array symmetric
  ! around 0. The spacing is provided by the equidistant WVLARR.
  ! --------------------------------------------------------------  
  sw = -hw1e ** 2.0
  delwvl = wvlarr(2) - wvlarr(1)
  mslit = NINT( SQRT( LOG(slit_trunc_limit) * sw) / delwvl) 
  num_slit = mslit * 2 + 1

  IF (num_slit > npoints) THEN
     mslit = (npoints-1) / 2 + 1
  ENDIF
 
  ! ---------------------------------------------------------------
  ! Convolve spectrum. First do the middle part, where we have full
  ! overlap coverage of the slit function. Again, remember the
  ! asymmetry of the Gaussian, which makes impossible a simple
  ! 50-50 division of the summation interval.
  ! ---------------------------------------------------------------

  ! Make a local copy of the NSLIT spectrum points to be convolved
  ! with the slit function. The spectrum points to be convolved are
  ! arranged such that the updated index corresponds to the maximum
  ! of the slit function (MSLIT). For simplicity we reflect the
  ! spectrum at the array end points.

  errstat = pge_errstat_ok

  ! get an array of slit variables lochwe for each wavelength position
  lochwe = 0.0

  fidx = 1; lidx = npoints
  DO iwin = 1, numwin
     IF (iwin < numwin) THEN
        upbnd = (winlim(iwin, 2) + winlim(iwin+1, 1))/2.0 
        lidx = MINVAL(MAXLOC(wvlarr, MASK=(wvlarr <= upbnd )))
     ELSE 
        lidx = npoints
     END IF
               
     fslit = MINVAL(MINLOC(slitwav(1:nslit), &
          MASK=(slitwav(1:nslit) >= wvlarr(fidx))))
     lslit = MINVAL(MAXLOC(slitwav(1:nslit), &
          MASK=(slitwav(1:nslit) <= wvlarr(lidx))))
 
     ! no slit between fidx:lidx, should never happen
     IF (fslit <= 0) THEN
        fslit = lslit
     ELSE IF (lslit <= 0) THEN
        lslit = fslit
     ENDIF

     IF (fslit > nslit .OR. lslit > nslit) THEN  
        WRITE(*, *) fslit, lslit, wvlarr(fidx), wvlarr(lidx), slitwav(1), &
             slitwav(nslit), nslit, wvlarr(1), wvlarr(npoints)
        WRITE(*, *) modulename, ': Not slit available for this window!!!'
        STOP
     ENDIF
     
     IF (lslit < fslit + 3) THEN  ! extrapolate, use the nearest value
        lochwe(fidx:lidx) = SUM(slitfit(fslit:lslit, hwe_idx, 1))/(lslit-fslit+1)
        lochwn(fidx:lidx) = SUM(slitfit(fslit:lslit, hwn_idx, 1))/(lslit-fslit+1)
     ELSE
        ! finter <= linter here
        finter = MINVAL(MINLOC(wvlarr, MASK = (wvlarr >= slitwav(fslit))))
        linter = MINVAL(MAXLOC(wvlarr, MASK = (wvlarr <= slitwav(lslit))))
        
        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, hwe_idx, 1), &
             linter-finter+1, wvlarr(finter:linter), lochwe(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
        
        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, hwn_idx, 1), &
             linter-finter+1, wvlarr(finter:linter), lochwn(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF

        IF (finter > fidx) THEN
           lochwe(fidx:finter-1)=slitfit(fslit, hwe_idx, 1)
           lochwn(fidx:finter-1)=slitfit(fslit, hwn_idx, 1)
        END IF
        
        IF (linter < lidx)  THEN
           lochwe(linter+1:lidx)=slitfit(lslit, hwe_idx, 1)
           lochwn(linter+1:lidx)=slitfit(lslit, hwn_idx, 1)
        END IF


     ENDIF
     fidx = lidx + 1
  ENDDO
    
  DO i = 1, npoints
     ! Get center value
     coeff = lochwn(i) / ( 2 *lochwe(i) *gamma(lochwn(i)))
     specmod(i) = specarr(i)
     slitsum = coeff

     ! Compute slit contribution from both sides
     sw    = -lochwe(i) ** 2.0
     DO j = 1, mslit
        ii = mslit + 1 - j
        j1 = i + ii; IF ( j1 > npoints ) j1 = npoints - MOD(j1, npoints)
        j2 = i - ii; IF ( j2 < 1 ) j2 = ABS(j2) + 2 
        !slit(j) = EXP( (wvlarr(j1)-wvlarr(i)) ** 2 / sw)
        slit(j) = coeff*EXP (-1.*(abs(wvlarr(j1)-wvlarr(i))/lochwe(i))**lochwn(i))

        slitsum = slitsum + 2.0 * slit(j)
        specmod(i) = specmod(i) + slit(j) * ( specarr(j1) + specarr(j2))
     ENDDO
     
     specmod(i) = specmod(i)*coeff / slitsum !normalized to center value
  END DO
  !stop ! not tested
  RETURN
END SUBROUTINE super_gauss_vary

SUBROUTINE gauss_uneven (wvlarr, specarr, npoints, nwin, slw, lbnd, ubnd)
  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : slit_trunc_limit
  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                             INTENT (IN)    :: npoints, nwin
  REAL (KIND=dp), DIMENSION(nwin),     INTENT (IN)    :: slw, lbnd, ubnd
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN)    :: wvlarr
  REAL (KIND=dp), DIMENSION (npoints), INTENT (INOUT) :: specarr

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                             :: i, iw, j1, j2, fidx, lidx
  REAL (KIND=dp)                      :: hw1e, wmid
  REAL (KIND=dp), DIMENSION (npoints) :: specmod, slit

  ! ----------------------------------------------------------------
  ! Initialization of output variable (default for "no convolution")
  ! ----------------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  slit(1:npoints) = 0.0
  
  fidx = 1
  DO iw = 1, nwin
     IF (iw == nwin) THEN
        lidx = npoints
     ELSE
        wmid = (ubnd(iw) + lbnd(iw+1)) / 2.0
        lidx = MINVAL(MAXLOC(wvlarr, MASK=(wvlarr <= wmid )))
     ENDIF
     hw1e = slw(iw)
     !print *, iw, fidx, lidx, hw1e
     !print *, wvlarr(fidx), wvlarr(lidx), slit_trunc_limit
     
     IF (hw1e /= 0.0) THEN
        DO i = fidx, lidx
           DO j1 = i, 1, -1
              slit(j1) = EXP(-((wvlarr(i)-wvlarr(j1))/hw1e)**2.)
              IF (slit(j1) < slit_trunc_limit .OR. j1 == 1) EXIT
           ENDDO
           
           DO j2 = i, npoints
              slit(j2) = EXP(-((wvlarr(i)-wvlarr(j2))/hw1e)**2.)
              IF (slit(j2) < slit_trunc_limit .OR. j2 == npoints) EXIT
           ENDDO

           !WRITE(*, '(3I5, 2D14.5)') j1, i, j2, MINVAL(slit(j1:j2)), MAXVAL(slit(j1:j2))

           specarr(i) = SUM(specmod(j1:j2)*slit(j1:j2)) / SUM(slit(j1:j2))
        ENDDO
     ENDIF

     fidx = lidx + 1
  ENDDO

  !DO i = 1, npoints
  !   WRITE(90, *) wvlarr(i), specmod(i), specarr(i)
  !ENDDO
  !STOP
  
  RETURN

END SUBROUTINE gauss_uneven

!xliu, 10/29/2009, add the capability to convole multiple spectrum simultaneously
SUBROUTINE convol_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)
  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : yn_varyslit, which_slit,instrument_idx,&
                                     omi_idx, omps_idx
  USE OMSAO_slitfunction_module
!  USE OMPS_slit_module
  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                        INTENT (IN)         :: nc, nf, nspec
  REAL (KIND=dp), DIMENSION (nf), INTENT (IN)         :: fwave
  REAL (KIND=dp), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
  REAL (KIND=dp), DIMENSION (nc), INTENT (IN)         :: cwave
  REAL (KIND=dp), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

  IF (.NOT. yn_varyslit) THEN
     IF (which_slit == 0) THEN
        CALL gauss_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
     ELSE IF (which_slit == 1) THEN
        CALL asym_gauss_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
     ELSE IF (which_slit == 2) THEN
        CALL asym_voigt_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
     ELSE IF (which_slit == 3) THEN
        CALL triangle_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
     ELSE IF (which_slit == 4) THEN
        CALL super_gauss_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
     ELSE IF (which_slit == 5) THEN 
       SELECT CASE (instrument_idx)
       CASE (omi_idx); CALL omislit_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
       !CASE (omps_idx); CALL ompsslit_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
       END SELECT
     END IF
  ELSE 
     IF (which_slit == 0) THEN
        CALL gauss_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
     ELSE IF (which_slit == 1) THEN
        CALL asym_gauss_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
     ELSE IF (which_slit == 2) THEN
        CALL asym_voigt_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
     ELSE IF (which_slit == 3) THEN
        CALL triangle_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
     ELSE IF (which_slit == 4) THEN 
        CALL super_gauss_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
     ELSE IF (which_slit == 5) THEN 
       SELECT CASE (instrument_idx)  
       CASE (omi_idx) ;  CALL omislit_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
       !CASE (omps_idx) ; call ompsslit_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
       END SELECT
     ENDIF
  ENDIF

  RETURN

END SUBROUTINE convol_f2c

!xliu, 10/22/2009, change the way of finding start & end positions on high resolution
!      grid for slit convolution
!xliu, 10/29/2009, add the capability to convole multiple spectrum simultaneously
SUBROUTINE gauss_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)
  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : solwinfit, winlim, numwin, slit_trunc_limit
  USE OMSAO_indices_module,   ONLY : hwe_idx
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
  INTEGER        :: i, j, iwin, iw, fidx, fidxc, lidx, lidxc, midx, sidx, eidx, nhalf
  REAL (KIND=dp) :: temp, hw1esq, dfw, ssum
  REAL (KIND=dp), DIMENSION (nf) :: slit
  
  !slit(:) = 0.0; cspec(:, :)= 0.0
  fidx = 1; fidxc = 1

  dfw  = fwave(2) - fwave(1) !xliu, 10/22/2009

  DO iwin = 1, numwin

     IF (iwin == numwin) THEN
        lidx = nf; lidxc = nc
     ELSE
        temp = (winlim(iwin, 2) + winlim(iwin + 1, 1)) / 2.0
        lidx =  MINVAL(MAXLOC(fwave, MASK=(fwave <= temp)))
        lidxc = MINVAL(MAXLOC(cwave, MASK=(cwave <= temp)))
     ENDIF    
     hw1esq = solwinfit(iwin, hwe_idx, 1)**2
     nhalf  = CEILING(solwinfit(iwin, hwe_idx, 1) / dfw * SQRT(-LOG(slit_trunc_limit))) !xliu, 10/22/2009

     !print *, iwin, fidx, lidx, fidxc, lidxc, nf, nc
     !print *, fwave(fidx), fwave(lidx), cwave(fidxc), cwave(lidxc)

     DO i = fidxc, lidxc
        ! Find the closest pixel
        midx = MINVAL(MAXLOC(fwave(fidx:lidx), MASK=(fwave(fidx:lidx) <= cwave(i)))) + fidx
        !DO j = midx, fidx, -1
        !   slit(j) = EXP(-(cwave(i) - fwave(j))**2 / hw1esq )
        !   IF (slit(j) <= slit_trunc_limit .OR. j == 1) EXIT
        !ENDDO
        !sidx = j
        !
        !DO j = midx+1, lidx
        !   slit(j) = EXP(-(cwave(i) - fwave(j))**2 / hw1esq )
        !   IF (slit(j) <= slit_trunc_limit .OR. j == nf) EXIT
        !ENDDO
        !eidx = j

        !xliu, 10/22/2009, replace above with following
        sidx = MAX(midx - nhalf, 1)
        eidx = MIN(nf, midx + nhalf)
        slit(sidx:eidx) = EXP(-(cwave(i) - fwave(sidx:eidx))**2 / hw1esq )
                 
        ssum = SUM(slit(sidx:eidx))
        DO j = 1, nspec
           cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx)) / ssum
        ENDDO
     ENDDO
     
     fidx = lidx + 1; fidxc = lidxc + 1
  ENDDO

  RETURN

END SUBROUTINE gauss_f2c

!xliu, 10/22/2009, change the way of finding start & end positions on high resolution
!      grid for slit convolution
!xliu, 10/29/2009, add the capability to convole multiple spectrum simultaneously
SUBROUTINE asym_gauss_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)
  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : solwinfit, winlim, numwin, slit_trunc_limit
  USE OMSAO_indices_module,   ONLY : hwe_idx, asy_idx
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
  INTEGER        :: i, j, iwin, iw, fidx, fidxc, lidx, lidxc, midx, sidx, eidx, nhalf1, nhalf2
  REAL (KIND=dp) :: temp, hw1esq, lasy, rasy, dfw, ssum
  REAL (KIND=dp), DIMENSION (nf) :: slit
  
  !slit(:) = 0.0; cspec(:, :)= 0.0
  dfw  = fwave(2) - fwave(1) !xliu, 10/22/2009

  fidx = 1; fidxc = 1
  DO iwin = 1, numwin

     IF (iwin == numwin) THEN
        lidx = nf; lidxc = nc
     ELSE
        temp = (winlim(iwin, 2) + winlim(iwin + 1, 1)) / 2.0
        lidx =  MINVAL(MAXLOC(fwave, MASK=(fwave <= temp)))
        lidxc = MINVAL(MAXLOC(cwave, MASK=(cwave <= temp)))
     ENDIF    
     hw1esq = solwinfit(iwin, hwe_idx, 1)**2
     lasy = hw1esq * (1.0 - solwinfit(iwin, asy_idx, 1))**2   
     rasy = hw1esq * (1.0 + solwinfit(iwin, asy_idx, 1))**2 

     !xliu, 10/22/2009
     nhalf1  = CEILING( SQRT(lasy) / dfw * SQRT(-LOG(slit_trunc_limit))) 
     nhalf2  = CEILING( SQRT(rasy) / dfw * SQRT(-LOG(slit_trunc_limit))) 
     
     DO i = fidxc, lidxc
        ! Find the closest pixel
        midx = MINVAL(MAXLOC(fwave(fidx:lidx), MASK=(fwave(fidx:lidx) <= cwave(i)))) + fidx
        
        !DO j = midx, fidx, -1
        !   slit(j) = EXP(-(cwave(i) - fwave(j))**2 / lasy )
        !   IF (slit(j) <= slit_trunc_limit .OR. j == 1) EXIT
        !ENDDO
        !sidx = j
        !
        !DO j = midx+1, lidx
        !   slit(j) = EXP(-(cwave(i) - fwave(j))**2 / rasy )
        !   IF (slit(j) <= slit_trunc_limit .OR. j == nf) EXIT
        !ENDDO
        !eidx = j

        !xliu, 10/22/2009, replace above with following
        sidx = MAX(midx - nhalf1, 1)
        eidx = MIN(nf, midx + nhalf2)
        slit(sidx:midx) = EXP(-(cwave(i) - fwave(sidx:midx))**2 / lasy )
        slit(midx+1:eidx) = EXP(-(cwave(i) - fwave(midx+1:eidx))**2 / rasy )
        
        ssum = SUM(slit(sidx:eidx))
        DO j = 1, nspec
           cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx)) / ssum
        ENDDO
     ENDDO
     
     fidx = lidx + 1; fidxc = lidxc + 1
  ENDDO

  RETURN

END SUBROUTINE asym_gauss_f2c

SUBROUTINE super_gauss_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)
  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : solwinfit, winlim, numwin, slit_trunc_limit
  USE OMSAO_indices_module,   ONLY : hwe_idx, hwn_idx
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
  INTEGER        :: i, j, iwin, iw, fidx, fidxc, lidx, lidxc, midx, sidx, eidx, nhalf
  REAL (KIND=dp) :: power, temp, hw1e, dfw, ssum, coeff
  REAL (KIND=dp), DIMENSION (nf) :: slit
  
  !slit(:) = 0.0; cspec(:, :)= 0.0
  fidx = 1; fidxc = 1
  dfw  = fwave(2) - fwave(1) 

  DO iwin = 1, numwin

     IF (iwin == numwin) THEN
        lidx = nf; lidxc = nc
     ELSE
        temp = (winlim(iwin, 2) + winlim(iwin + 1, 1)) / 2.0
        lidx =  MINVAL(MAXLOC(fwave, MASK=(fwave <= temp)))
        lidxc = MINVAL(MAXLOC(cwave, MASK=(cwave <= temp)))
     ENDIF    
     power  = solwinfit(iwin, hwn_idx, 1)
     hw1e   = solwinfit(iwin, hwe_idx, 1)
     nhalf  = CEILING(solwinfit(iwin, hwe_idx, 1) / dfw * SQRT(-LOG(slit_trunc_limit))) !xliu, 10/22/2009
     coeff  = power / ( 2 * hw1e *gamma(power))
     coeff = 1.0
     DO i = fidxc, lidxc
        ! Find the closest pixel
        midx = MINVAL(MAXLOC(fwave(fidx:lidx), MASK=(fwave(fidx:lidx) <= cwave(i)))) + fidx
        sidx = MAX(midx - nhalf, 1)
        eidx = MIN(nf, midx + nhalf)
        slit (sidx:eidx) =coeff* EXP (-1.*(abs(cwave(i) - fwave(sidx:eidx))/hw1e)**power)         
        ssum = SUM(slit(sidx:eidx))
        DO j = 1, nspec
           cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx))/ ssum
        ENDDO
     ENDDO
     
     fidx = lidx + 1; fidxc = lidxc + 1
  ENDDO

  RETURN

END SUBROUTINE super_gauss_f2c

!xliu, 10/22/2009, change the way of finding start & end positions on high resolution
!      grid for slit convolution
!xliu, 10/29/2009, add the capability to convole multiple spectrum simultaneously
SUBROUTINE gauss_vary_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)
  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : winlim, numwin, slit_trunc_limit, slitwav, slitfit, nslit
  USE OMSAO_indices_module,   ONLY : hwe_idx
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
       midx, sidx, eidx, fslit,lslit, finter, linter, errstat, nhalf
  REAL (KIND=dp) :: temp, hw1esq, dfw, ssum
  REAL (KIND=dp), DIMENSION (nf) :: slit
  REAL (KIND=dp), DIMENSION (nc) :: lochwe

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=13), PARAMETER :: modulename = 'gauss_vary_f2c'

  errstat = pge_errstat_ok
  
  !slit(:) = 0.0; cspec(:, :)= 0.0
  dfw  = fwave(2) - fwave(1) !xliu, 10/22/2009

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
        lochwe(fidxc:lidxc) = SUM(slitfit(fslit:lslit, hwe_idx, 1))/(lslit-fslit+1)
     ELSE
        ! finter <= linter here
        finter = MINVAL(MINLOC(cwave, MASK = (cwave >= slitwav(fslit))))
        linter = MINVAL(MAXLOC(cwave, MASK = (cwave <= slitwav(lslit))))
        
        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, hwe_idx, 1), &
             linter-finter+1, cwave(finter:linter), lochwe(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
        
        IF (finter > fidxc) THEN
           lochwe(fidxc:finter-1)=slitfit(fslit, hwe_idx, 1)
        END IF
        
        IF (linter < lidxc)  THEN
           lochwe(linter+1:lidxc)=slitfit(lslit, hwe_idx, 1)
        END IF
     ENDIF
   
     DO i = fidxc, lidxc

        hw1esq = lochwe(i)**2
        nhalf  = CEILING(lochwe(i) / dfw * SQRT(-LOG(slit_trunc_limit))) !xliu, 10/22/2009

        ! Find the closest pixel
        midx = MINVAL(MAXLOC(fwave(fidx:lidx), MASK=(fwave(fidx:lidx) <= cwave(i)))) + fidx

        !DO j = midx, fidx, -1
        !   slit(j) = EXP(-(cwave(i) - fwave(j))**2 / hw1esq )
        !   IF (slit(j) <= slit_trunc_limit .OR. j == 1) EXIT
        !ENDDO
        !sidx = j
        !
        !DO j = midx+1, lidx
        !   slit(j) = EXP(-(cwave(i) - fwave(j))**2 / hw1esq )
        !   IF (slit(j) <= slit_trunc_limit .OR. j == nf) EXIT
        !ENDDO
        !eidx = j

        !xliu, 10/22/2009, replace above with following
        sidx = MAX(midx - nhalf, 1)
        eidx = MIN(nf, midx + nhalf)
        slit(sidx:eidx) = EXP(-(cwave(i) - fwave(sidx:eidx))**2 / hw1esq )
        
        ssum = SUM(slit(sidx:eidx))
        DO j = 1, nspec
           cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx)) / ssum
        ENDDO
     ENDDO
     
     fidx = lidx + 1; fidxc = lidxc + 1
  ENDDO

  RETURN

END SUBROUTINE gauss_vary_f2c

!xliu, 10/22/2009, change the way of finding start & end positions on high resolution
!      grid for slit convolution
!xliu, 10/29/2009, add the capability to convole multiple spectrum simultaneously
SUBROUTINE asym_gauss_vary_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)
  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : winlim, numwin, slit_trunc_limit, slitwav, slitfit, nslit
  USE OMSAO_indices_module,   ONLY : hwe_idx, asy_idx
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
       midx, sidx, eidx, fslit,lslit, finter, linter, errstat, nhalf1, nhalf2
  REAL (KIND=dp) :: temp, hw1esq, lasy, rasy, dfw, ssum
  REAL (KIND=dp), DIMENSION (nf) :: slit
  REAL (KIND=dp), DIMENSION (nc) :: lochwe, locasy

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=18), PARAMETER :: modulename = 'asym_gauss_vary_f2c'

  errstat = pge_errstat_ok
  
  !slit(:) = 0.0; cspec(:, :)= 0.0
  dfw  = fwave(2) - fwave(1) !xliu, 10/22/2009
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
        lochwe(fidxc:lidxc) = SUM(slitfit(fslit:lslit, hwe_idx, 1))/(lslit-fslit+1)
        locasy(fidxc:lidxc) = SUM(slitfit(fslit:lslit, asy_idx, 1))/(lslit-fslit+1)
     ELSE
        ! finter <= linter here
        finter = MINVAL(MINLOC(cwave, MASK = (cwave >= slitwav(fslit))))
        linter = MINVAL(MAXLOC(cwave, MASK = (cwave <= slitwav(lslit))))
        
        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, hwe_idx, 1), &
             linter-finter+1, cwave(finter:linter), lochwe(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF

        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, asy_idx, 1), &
             linter-finter+1, cwave(finter:linter), locasy(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
        
        IF (finter > fidxc) THEN
           lochwe(fidxc:finter-1)=slitfit(fslit, hwe_idx, 1)
           locasy(fidxc:finter-1)=slitfit(fslit, asy_idx, 1)
        END IF
        
        IF (linter < lidxc)  THEN
           lochwe(linter+1:lidxc)=slitfit(lslit, hwe_idx, 1)
           locasy(linter+1:lidxc)=slitfit(lslit, asy_idx, 1)
        END IF
     ENDIF
   
     DO i = fidxc, lidxc

        hw1esq = lochwe(i)**2
        lasy = hw1esq * (1.0 - locasy(i))**2   
        rasy = hw1esq * (1.0 + locasy(i))**2 

        !xliu, 10/22/2009
        nhalf1  = CEILING( SQRT(lasy) / dfw * SQRT(-LOG(slit_trunc_limit))) 
        nhalf2  = CEILING( SQRT(rasy) / dfw * SQRT(-LOG(slit_trunc_limit))) 

        ! Find the closest pixel
        midx = MINVAL(MAXLOC(fwave(fidx:lidx), MASK=(fwave(fidx:lidx) <= cwave(i)))) + fidx

        !DO j = midx, fidx, -1
        !   slit(j) = EXP(-(cwave(i) - fwave(j))**2 / lasy )
        !   IF (slit(j) <= slit_trunc_limit .OR. j == 1) EXIT
        !ENDDO
        !sidx = j
        !
        !DO j = midx+1, lidx
        !   slit(j) = EXP(-(cwave(i) - fwave(j))**2 / rasy )
        !   IF (slit(j) <= slit_trunc_limit .OR. j == nf) EXIT
        !ENDDO
        !eidx = j

        !xliu, 10/22/2009, replace above with following
        sidx = MAX(midx - nhalf1, 1)
        eidx = MIN(nf, midx + nhalf2)
        slit(sidx:midx) = EXP(-(cwave(i) - fwave(sidx:midx))**2 / lasy )
        slit(midx+1:eidx) = EXP(-(cwave(i) - fwave(midx+1:eidx))**2 / rasy )
        
        ssum = SUM(slit(sidx:eidx))
        DO j = 1, nspec
           cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx)) / ssum
        ENDDO
     ENDDO
     
     fidx = lidx + 1; fidxc = lidxc + 1
  ENDDO

  RETURN

END SUBROUTINE asym_gauss_vary_f2c

SUBROUTINE super_gauss_vary_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)
  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : winlim, numwin, slit_trunc_limit, slitwav, slitfit, nslit
  USE OMSAO_indices_module,   ONLY : hwe_idx, hwn_idx
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
       midx, sidx, eidx, fslit,lslit, finter, linter, errstat, nhalf
  REAL (KIND=dp) :: temp, hw1e, power, dfw, ssum, coeff
  REAL (KIND=dp), DIMENSION (nf) :: slit
  REAL (KIND=dp), DIMENSION (nc) :: lochwe, lochwn

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=13), PARAMETER :: modulename = 'super_gauss_vary_f2c'

  errstat = pge_errstat_ok 
  dfw  = fwave(2) - fwave(1) 
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
        lochwe(fidxc:lidxc) = SUM(slitfit(fslit:lslit, hwe_idx, 1))/(lslit-fslit+1)
        lochwn(fidxc:lidxc) = SUM(slitfit(fslit:lslit, hwn_idx, 1))/(lslit-fslit+1)
     ELSE
        ! finter <= linter here
        finter = MINVAL(MINLOC(cwave, MASK = (cwave >= slitwav(fslit))))
        linter = MINVAL(MAXLOC(cwave, MASK = (cwave <= slitwav(lslit))))
        
        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, hwe_idx, 1), &
             linter-finter+1, cwave(finter:linter), lochwe(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF

        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, hwn_idx, 1), &
             linter-finter+1, cwave(finter:linter), lochwn(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
        
        IF (finter > fidxc) THEN
           lochwe(fidxc:finter-1)=slitfit(fslit, hwe_idx, 1)
           lochwn(fidxc:finter-1)=slitfit(fslit, hwn_idx, 1)
        END IF
        
        IF (linter < lidxc)  THEN
           lochwe(linter+1:lidxc)=slitfit(lslit, hwe_idx, 1)
           lochwn(linter+1:lidxc)=slitfit(lslit, hwn_idx, 1)
        END IF
     ENDIF
   
     DO i = fidxc, lidxc

        hw1e   = lochwe(i)
        power  = lochwn(i)
        coeff  = power / ( 2 * hw1e *gamma(power))
        nhalf  = CEILING(lochwe(i) / dfw * SQRT(-LOG(slit_trunc_limit))) !xliu, 10/22/2009

        ! Find the closest pixel
        midx = MINVAL(MAXLOC(fwave(fidx:lidx), MASK=(fwave(fidx:lidx) <= cwave(i)))) + fidx

        !xliu, 10/22/2009, replace above with following
        sidx = MAX(midx - nhalf, 1)
        eidx = MIN(nf, midx + nhalf)
        slit(sidx:eidx) = coeff*EXP (-1.*(abs(cwave(i) - fwave(sidx:eidx))/hw1e)**power)         
     
        ssum = SUM(slit(sidx:eidx))
        DO j = 1, nspec
           cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx)) / ssum
        ENDDO
     ENDDO
     
     fidx = lidx + 1; fidxc = lidxc + 1
  ENDDO

  RETURN

END SUBROUTINE super_gauss_vary_f2c

