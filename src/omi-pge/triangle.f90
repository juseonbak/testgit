!  convolves input spectrum with trainglular slit function of specified fwhm
SUBROUTINE triangle (wvlarr, specarr, specmod, npoints, fwhm)
 
  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                             INTENT (IN)    :: npoints
  REAL (KIND=dp),                      INTENT (IN)    :: fwhm
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN)    :: wvlarr, specarr

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=dp), DIMENSION (npoints), INTENT (OUT) :: specmod

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                             :: nhi, nlo, i, j, num_slit
  REAL (KIND=dp)                      :: slitsum, slit0, delwvl
  REAL (KIND=dp), DIMENSION (npoints) :: slit

  ! ----------------------------------------------------------------
  ! Initialization of output variable (default for "no convolution")
  ! ----------------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  
  ! --------------------------------------
  ! No convolution if halfwidth @ 1/e is 0
  ! --------------------------------------
  IF ( fwhm == 0.0 ) RETURN

  delwvl = wvlarr(2) - wvlarr(1)     ! assume evenly-spaced

  ! Apply slit function convolution
  ! Calculate slit function values out to 0.00 times x0 value,
  ! Normalize so that sum = 1.

  slitsum = 1.0 ;  slit0 = 1.0
  i = 1  ;  num_slit = 0
  numslit: DO WHILE ( num_slit <= npoints )
     slit (i) = 1.0 - delwvl * i / fwhm
     slitsum = slitsum + 2.0 * slit (i)
     IF (slit(i) / slit0 < 0.0 ) EXIT
     i = i + 1
  END DO numslit
  num_slit = i
  
  IF (i /= npoints) num_slit = num_slit - 1
  
  slit0 = slit0 / slitsum
  slit(1:num_slit) = slit(1:num_slit) / slitsum

  ! Convolve spectrum.  reflect at endpoints.
  specmod(1:npoints) = slit0 * specarr(1:npoints)
  DO i = 1, npoints
     DO j = 1, num_slit
        nlo = MAXVAL( (/ 1, ABS(i - j) /) )
        nhi = i + j
        IF ( nhi > npoints ) nhi = MAXVAL ( (/ 1, 2*npoints-nhi /) )
        specmod(i) = specmod(i) + slit(j) * ( specarr(nlo) + specarr(nhi) )
     END DO
  END DO

  RETURN
END SUBROUTINE triangle


SUBROUTINE triangle_multi (wvlarr, specarr, specmod, npoints)

  USE OMSAO_precision_module
  USE OMSAO_variables_module,  ONLY : solwinfit, winlim, numwin
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
  INTEGER                             :: i, j, j1, j2, num_slit, mslit, sslit, &
       eslit, odd, fidx, lidx, iwin
  REAL (KIND=dp)                      :: delwvl, slitsum,  maxslit, fwhm, e_asym
  REAL (KIND=dp), DIMENSION (npoints) :: slit, locwvl, locsli, xxx, upbnd
  INTEGER,        DIMENSION (npoints) :: idx
  REAL (KIND=dp), EXTERNAL            :: signdp

  ! --------------------------------------------------------
  ! Initialize output variable (default for "no convolution"
  ! --------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  
  fidx = 1
  DO iwin = 1, numwin

     fwhm = solwinfit(iwin, hwe_idx, 1)
     
     IF (iwin < numwin) THEN
        upbnd = (winlim(iwin, 2) + winlim(iwin+1, 1))/2.0 
        lidx = MINVAL(MAXLOC(wvlarr, MASK=(wvlarr <= upbnd )))
     ELSE 
        lidx = npoints
     END IF

     ! -----------------------------------------------
     ! No Gaussian convolution if Halfwidth @ 1/e is 0
     ! -----------------------------------------------
     IF ( fwhm == 0.0 ) RETURN

     ! --------------------------------------------------------------
     ! Find the number of spectral points that fall within a Gaussian
     ! slit function with values >= 0.01. Remember that we have an
     ! asymmetric Gaussian, so we create a wavelength array symmetric
     ! around 0. The spacing is provided by the equidistant WVLARR.
     ! --------------------------------------------------------------
     delwvl = wvlarr(2) - wvlarr(1) ; locsli = 0.0
     DO i = -npoints/2, npoints/2
        j = MIN ( i + npoints/2 + 1, npoints )
        locwvl(j) = ABS(delwvl * REAL(i,KIND=dp))
        locsli(j) = 1.0 - ABS(locwvl(j)) / fwhm 
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
     WHERE ( locsli < 0.0)
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
     ! Find the maximum of the slit
     ! --------------------------------------------------------
     mslit = MAXVAL( MAXLOC ( slit(1:num_slit) ) )

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
        locsli(1:num_slit) = 0.0
        ! First do the right half of the slit function
        DO j = mslit, num_slit
           j1 = i+(j-mslit)
           IF ( j1 > npoints ) j1 = npoints - MOD(j1, npoints)
           locsli(j) = specarr(j1)
           idx(j) = j1
        END DO
        ! Now the left half of the slit function
        DO j = mslit-1, 1, -1
           j2 = i-((mslit-1)-j) - 1 ; IF ( j2 < 1 ) j2 = ABS(j2) + 2
           locsli(j) = specarr(j2)
           idx(j) = j2
        END DO
        specmod(i) = DOT_PRODUCT(slit(1:num_slit), locsli(1:num_slit))
     END DO

     fidx = lidx + 1
  END DO

  RETURN
END SUBROUTINE triangle_multi



SUBROUTINE triangle_vary (wvlarr, specarr, specmod, npoints)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : slitwav, slitfit, nslit, winlim, numwin
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
  INTEGER                             :: i, j, j1, j2, num_slit, mslit, sslit, &
       eslit, odd, errstat,  fidx, lidx, fslit, lslit, finter, linter, iwin
  REAL (KIND=dp)                      :: delwvl, slitsum,  maxslit, fwhm
  REAL (KIND=dp), DIMENSION (npoints) :: slit, locwvl, locsli, xxx, lochwe
  INTEGER,        DIMENSION (npoints) :: idx

  REAL (KIND=dp) :: signdp
  EXTERNAL signdp

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=13), PARAMETER :: modulename = 'triangle_vary'


  !WRITE(*, *) 'nslit = ', nslit
  !WRITE(*, '(10f8.3)') slitwav(1:nslit), slitfit(1:nslit, hwe_idx, 1)

  ! --------------------------------------------------------
  ! Initialize output variable (default for "no convolution"
  ! --------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  
  ! -----------------------------------------------
  ! No Gaussian convolution if Halfwidth @ 1/e is 0
  ! -----------------------------------------------
  i = MAXVAL(MAXLOC(slitfit(1:nslit, hwe_idx, 1)))
  fwhm = slitfit(i, hwe_idx, 1)

  IF ( fwhm == 0.0 ) RETURN

  ! --------------------------------------------------------------
  ! Find the number of spectral points that fall within a Gaussian
  ! slit function with values >= 0.01. Remember that we have an
  ! asymmetric Gaussian, so we create a wavelength array symmetric
  ! around 0. The spacing is provided by the equidistant WVLARR.
  ! --------------------------------------------------------------
     
  delwvl = wvlarr(2) - wvlarr(1) ; locsli = 0.0
  DO i = -npoints/2, npoints/2
     j = MIN ( i + npoints/2 + 1, npoints )
     locwvl(j) = ABS(delwvl * REAL(i,KIND=dp))
     locsli(j) = 1.0 - locwvl(j) / fwhm 
  END DO

  ! ------------------------------------------------------------------------
  ! Find the array entries that mark the region of SLIT >= 0.0001*MAX(SLIT).
  ! We start by initializing the auxilliary IDX as IDX(n) = n. We then set
  ! all those entries to 0 that are outside the accepted minimum value of 
  ! the slit function. The start and end indices simply follow as the 
  ! minimum and maximum values of IDX where IDX /= 0.
  ! ------------------------------------------------------------------------
  maxslit = MAXVAL(locsli(1:npoints))
  mslit = MAXVAL( MAXLOC ( locsli(1:npoints) ) )

  sslit = 0 ; eslit = 0 ; idx = (/ (i, i = 1, npoints) /)
  WHERE ( locsli < 0.0 )
     idx = 0
  END WHERE
  sslit = MINVAL( idx, MASK = (idx /= 0) )
  eslit = MAXVAL( idx, MASK = (idx /= 0) )
 
  ! ----------------------------------------------------------
  ! Compute number of slit function points, just based on the 
  ! maximum slit width
  ! ----------------------------------------------------------
  IF (mslit - sslit > eslit - mslit) THEN
     eslit = mslit + (mslit - sslit)
  ELSE 
     sslit = mslit - (eslit - mslit)
  END IF
  num_slit = eslit - sslit + 1
  mslit = num_slit / 2 + 1 
  
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
  lochwe = 0.0

  fidx = 1; lidx = npoints
          
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
    
  DO i = 1, npoints
     ! get the fwhm and e_asym at this point
     fwhm = lochwe(i)
     
     locsli(1:num_slit) = 0.0
     slit(1:num_slit)   = 0.0
     ! First do the right half of the slit function
     DO j = mslit, num_slit
        j1 = i+(j-mslit)
        IF (j1 > npoints) j1 = npoints - MOD(j1, npoints)
        locsli(j) = specarr(j1)
        delwvl = (wvlarr(j1)-wvlarr(i))
        slit(j) = 1.0 - ABS(delwvl) / fwhm
        IF (slit(j) < 0.0) slit(j) = 0.0
        idx(j) = j1
     END DO
        
     ! Now the left half of the slit function
     DO j = mslit-1, 1, -1
        j2 = i-((mslit-1)-j) - 1
        IF(j2 < 1) j2 = ABS(j2) + 2
        locsli(j) = specarr(j2)
        delwvl = (wvlarr(j2)-wvlarr(i))
        slit(j) = 1.0 - ABS(delwvl) / fwhm
        IF (slit(j) < 0.0) slit(j) = 0.0
        idx(j) = j2
     END DO
     
     specmod(i) = DOT_PRODUCT(slit(1:num_slit), locsli(1:num_slit)) &
          / SUM(slit(1:num_slit))
  END DO

  RETURN
END SUBROUTINE triangle_vary


SUBROUTINE triangle_uneven (wvlarr, specarr, npoints, nwin, slw, lbnd, ubnd)
  
  USE OMSAO_precision_module
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
  REAL (KIND=dp)                      :: fwhm, wmid
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
     fwhm = slw(iw)
     !print *, iw, fidx, lidx, hw1e
     !print *, wvlarr(fidx), wvlarr(lidx), slit_trunc_limit
     
     IF (fwhm /= 0.0) THEN
        DO i = fidx, lidx
           DO j1 = i, 1, -1
              slit(j1) = 1.0 - ABS(wvlarr(i)-wvlarr(j1))/fwhm
              IF (slit(j1) < 0 .OR. j1 == 1) EXIT
           ENDDO
           IF (slit(j1) < 0) slit(j1) = 0.0
           
           DO j2 = i, npoints
              slit(j2) = 1.0 - ABS(wvlarr(i)-wvlarr(j1))/fwhm
              IF (slit(j2) < 0 .OR. j2 == npoints) EXIT
           ENDDO
           IF (slit(j2) < 0) slit(j2) = 0.0

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

END SUBROUTINE triangle_uneven

!xliu, 10/29/2009, add the capability to convole multiple spectrum simultaneously
SUBROUTINE triangle_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)
  
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
  INTEGER        :: i, j, iwin, iw, fidx, fidxc, lidx, lidxc, midx, sidx, eidx
  REAL (KIND=dp) :: temp, fwhm, ssum
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
     fwhm = solwinfit(iwin, hwe_idx, 1)

     DO i = fidxc, lidxc
        ! Find the closest pixel
        midx = MINVAL(MAXLOC(fwave(fidx:lidx), MASK=(fwave(fidx:lidx) <= cwave(i)))) + fidx

        DO j = midx, fidx, -1
           slit(j) = 1.0 - ABS(fwave(j)-cwave(i)) / fwhm
           IF (slit(j) <= 0.0 .OR. j == 1) EXIT
        ENDDO
        IF (slit(j) < 0.0 ) slit(j) = 0.0
        sidx = j

        DO j = midx, lidx
           slit(j) = 1.0 - ABS(fwave(j)-cwave(i)) / fwhm
           IF (slit(j) <= 0.0 .OR. j == nf) EXIT
        ENDDO
        IF (slit(j) < 0.0 ) slit(j) = 0.0
        eidx = j
        
        ssum = SUM(slit(sidx:eidx))
        DO j = 1, nspec
           cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx)) / ssum
        ENDDO
     ENDDO
     
     fidx = lidx + 1; fidxc = lidxc + 1
  ENDDO

  RETURN

END SUBROUTINE triangle_f2c

!xliu, 10/29/2009, add the capability to convole multiple spectrum simultaneously
SUBROUTINE triangle_vary_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)
  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : winlim, numwin, slitwav, slitfit, nslit
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
       midx, sidx, eidx, fslit,lslit, finter, linter, errstat
  REAL (KIND=dp) :: temp, fwhm, ssum
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

        fwhm = lochwe(i)

        ! Find the closest pixel
        midx = MINVAL(MAXLOC(fwave(fidx:lidx), MASK=(fwave(fidx:lidx) <= cwave(i)))) + fidx

        DO j = midx, fidx, -1
           slit(j) = 1.0 - ABS(fwave(j)-cwave(i)) / fwhm
           IF (slit(j) <= 0.0 .OR. j == 1) EXIT
        ENDDO
        IF (slit(j) < 0.0) slit(j) = 0.0
        sidx = j

        DO j = midx, lidx
           slit(j) = 1.0 - ABS(fwave(j)-cwave(i)) / fwhm
           IF (slit(j) <= 0.0 .OR. j == nf) EXIT
        ENDDO
        IF (slit(j) < 0.0) slit(j) = 0.0
        eidx = j
        
        ssum = SUM(slit(sidx:eidx))
        DO j = 1, nspec
           cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx)) / ssum
        ENDDO
     ENDDO
     
     fidx = lidx + 1; fidxc = lidxc + 1
  ENDDO

  RETURN

END SUBROUTINE triangle_vary_f2c

