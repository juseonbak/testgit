MODULE OMSAO_omicloud_module
 
  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxchlen
  USE OMSAO_indices_module
  USE OMSAO_omidata_module,    ONLY: nxtrack_max, nlines_max, ntimes_max, ncoadd, nfxtrack, &
       nxbin, nybin, offset_line, zoom_mode, zoom_p1
  USE OMSAO_variables_module,  ONLY: l2_cld_filename, coadd_uv2
  USE OMSAO_he5_module
  USE OMSAO_errstat_module
  IMPLICIT NONE

  ! -------------------------------
  ! TYPE declaration for OMI clouds
  ! -------------------------------
  INTEGER (KIND=i4), PARAMETER                                  :: nbit = 16
  TYPE, PUBLIC :: OMI_CloudBlock
     REAL (KIND=r4) :: CFRmissing, CFRscale, CFRoffset, CTPmissing, CTPscale, CTPoffset
     REAL (KIND=r4),   DIMENSION (nxtrack_max,0:ntimes_max-1)   :: cfr
     REAL (KIND=r4),   DIMENSION (nxtrack_max,0:ntimes_max-1)   :: ctp
     REAL (KIND=r4),   DIMENSION (nxtrack_max,0:ntimes_max-1)   :: ai
     INTEGER(KIND=i1), DIMENSION (nxtrack_max, 0:ntimes_max-1)  :: qflags
     CHARACTER (LEN=maxchlen)                                   :: PGEversion
  END TYPE OMI_CloudBlock  
  TYPE (OMI_CloudBlock),   PUBLIC :: OMIL2_clouds

  ! ----------------------------------------------
  ! Names of various HE5 fields to read from files
  ! ----------------------------------------------
  CHARACTER (LEN=19), PARAMETER, PRIVATE :: omicld_cfrac1_field  = 'CloudFractionforO3'
  CHARACTER (LEN=22), PARAMETER, PRIVATE :: omicld_cfracd_field  = 'CloudFractionPrecision'
  CHARACTER (LEN=19), PARAMETER, PRIVATE :: omicld_cpres1_field  = 'CloudPressureforO3'
  CHARACTER (LEN=22), PARAMETER, PRIVATE :: omicld_cpresd_field  = 'CloudPressurePrecision'
  CHARACTER (LEN=13), PARAMETER, PRIVATE :: omicld_cfrac_field   = 'CloudFraction'
  CHARACTER (LEN=13), PARAMETER, PRIVATE :: omicld_cpres_field   = 'CloudPressure'
  CHARACTER (LEN=27), PARAMETER, PRIVATE :: omicld_qflag1_field   = 'ProcessingQualityFlagsforO3'
  CHARACTER (LEN=22), PARAMETER, PRIVATE :: omicld_qflag_field   = 'ProcessingQualityFlags'
  CHARACTER (LEN=14), PARAMETER, PRIVATE :: omicld_ai_field      = 'UVAerosolIndex'
  
CONTAINS
  
  SUBROUTINE read_omicldrr_clouds (nt, nx, nl, pge_error_status)
    
    ! ----------------------
    ! Input/Output variables
    ! ----------------------  
    INTEGER (KIND=i4), INTENT (IN) :: nt, nx, nl
    INTEGER, INTENT (OUT)          :: pge_error_status
    
    ! -----------------------
    ! Local variables
    ! -----------------------
    INTEGER   (KIND=i4)      :: locerrstat, swath_id, swath_file_id, nt_loc, nx_loc, &
         ix, iix, i, ii, j, k, nline, nbx, nbin, sline, eline, nx1, xoff, xoff0
    REAL      (KIND=r4)      :: scale_cfr, offset_cfr, missval_cfr, scale_ctp, &
         offset_ctp, missval_ctp, scfr, scfr1, scfr0, tmpsum
    CHARACTER (LEN=maxchlen) :: swath_name
    REAL    (KIND=r4), DIMENSION (nx, 0:nt-1)           :: cfr, ctp, ai
    INTEGER (KIND=i2), DIMENSION (nx, 0:nt-1)           :: qflag
    INTEGER (KIND=i2), DIMENSION (nx, 0:nt-1, 0:nbit-1) :: flgbits

    ! External functions
    INTEGER                  :: OMI_SMF_setmsg
    INTEGER                  :: estat
    
    ! ----------------------
    ! Name of the subroutine
    ! ----------------------
    CHARACTER (LEN=20), PARAMETER :: modulename = 'read_omicldrr_clouds'
    
    locerrstat = pge_errstat_ok
    CALL he5_init_input_file ( l2_cld_filename, swath_name,  swath_id, swath_file_id, &
         nt_loc, nx_loc, locerrstat)
    
    IF ( locerrstat /= pge_errstat_ok ) THEN
       estat = OMI_SMF_setmsg (omsao_e_open_l2cld_file, modulename, 'OML2 Cloud access failed', 0) ; STOP 1
    END IF

    ! check if the dimensions are consistent with OMI L1B data
    IF (nt_loc /= nt .OR. nx_loc /= nx) THEN
       print *, 'OMICLDRR: ', nt_loc, nx_loc
       print *, 'OMIL2 : ', nt, nx
       WRITE(*, '(A)') 'Inconsistent dimensions between OMIL1B and OMIL2!!!'
       pge_error_status = pge_errstat_error; RETURN
    ENDIF

    OMIL2_clouds%PGEversion     = '0.9.5.0'
    scale_cfr = 1.0; offset_cfr = 0.0; missval_cfr = 0.0
    scale_ctp = 1.0; offset_ctp = 0.0; missval_ctp = 0
    OMIL2_clouds%CFRmissing = missval_cfr
    OMIL2_clouds%CFRscale   = scale_cfr
    OMIL2_clouds%CFRoffset  = offset_cfr
    OMIL2_clouds%CTPmissing = missval_ctp
    OMIL2_clouds%CTPscale   = scale_ctp
    OMIL2_clouds%CTPoffset  = offset_ctp

    ! Read all the data for the swath
    he5_start_2d = (/ 0, 0 /) ; he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nx, nt /)
    
    ! ----------------------------------------------------------------
    ! Read cloud fraction and cloud top pressure, and check for error.
    ! Eventually we may read the cloud uncertainties also, but for the
    ! first version we stick with just the basic cloud products.
    ! ----------------------------------------------------------------
    locerrstat = HE5_SWrdfld (swath_id, omicld_cfrac1_field,   &
         he5_start_2d, he5_stride_2d, he5_edge_2d, cfr(1:nx, 0:nt-1) )
    IF ( locerrstat /= pge_errstat_ok ) THEN
       estat = OMI_SMF_setmsg (omsao_e_open_l2cld_file, modulename, 'OML2 Cloud CFR access failed', 0); STOP 1
    END IF
    
    locerrstat = HE5_SWrdfld (swath_id, omicld_cpres1_field,   &
         he5_start_2d, he5_stride_2d, he5_edge_2d, ctp(1:nx, 0:nt-1) )
    IF ( locerrstat /= pge_errstat_ok ) THEN
       estat = OMI_SMF_setmsg (omsao_e_open_l2cld_file, modulename, 'OML2 Cloud CTP access failed', 0); STOP 1
    END IF

    locerrstat = HE5_SWrdfld (swath_id, omicld_qflag1_field,   &
         he5_start_2d, he5_stride_2d, he5_edge_2d, qflag(1:nx, 0:nt-1) )
    IF ( locerrstat /= pge_errstat_ok ) THEN
       estat = OMI_SMF_setmsg (omsao_e_open_l2cld_file, modulename, 'OML2 Cloud Quality Flags access failed', 0); STOP 1
    END IF

    !locerrstat = HE5_SWrdfld (swath_id, omicld_ai_field,   &
    !     he5_start_2d, he5_stride_2d, he5_edge_2d, ai(1:nx, 0:nt-1) )
    !IF ( locerrstat /= pge_errstat_ok ) THEN
    !   estat = OMI_SMF_setmsg (omsao_e_open_l2cld_file, modulename, 'OML2 Cloud UV aerosol index access failed', 0); STOP 1
    !END IF

    ! ------------------------------------
    ! Save cloud arrays in OMI Cloud Block
    ! ------------------------------------
    IF (zoom_mode) THEN
       nx1 = nx / 2
    ELSE 
       nx1 = nx
    ENDIF

    cfr(1:nx1,0:nt-1) = cfr(1:nx1,0:nt-1) * OMIL2_clouds%CFRscale + OMIL2_clouds%CFRoffset
    ctp(1:nx1,0:nt-1) = ctp(1:nx1,0:nt-1) * OMIL2_clouds%CTPscale + OMIL2_clouds%CTPoffset
    ai(1:nx1, 0:nt-1) = 0.0
    flgbits = 0
    sline = offset_line;  nline = nl * nybin; eline  = offset_line + nline - 1
    DO ix = 1, nx1
       !CALL convert_2bytes_to_16bits (nbit, nline, qflag(ix, sline:eline), flgbits(ix, sline:eline, 0:nbit-1))
       CALL convert_2bytes_to_16bits (nbit, nt, qflag(ix, 0:nt-1), flgbits(ix, 0:nt-1, 0:nbit-1))
    ENDDO

    nbin = nxbin; IF (coadd_uv2) nbin = nbin * ncoadd
    nbx = nx1 / nbin
    IF (zoom_mode .AND. MOD(zoom_p1, 2) == 0) THEN
       nbx = nbx - 1; xoff0 = 1 
    ELSE
       xoff0 = 0
    ENDIF

    OMIL2_clouds%cfr   (1:nbx, 0:nl-1) = 0.0
    OMIL2_clouds%ctp   (1:nbx, 0:nl-1) = 0.0
    OMIL2_clouds%ai    (1:nbx, 0:nl-1) = 0.0
    OMIL2_clouds%qflags(1:nbx, 0:nl-1) = 0

    ! Fill in cloud top pressure values for bad pixels (interpolation/extrapolation)
    ! 0   - failed convergence check
    ! 1   - solar zenith angle, lat., lon., out of range (SZA > 88 deg) 
    ! 2   - cloud pressure less than low range of table
    ! 3   - cloud pressure greater than surface pressure
    ! 4   - matrix inversion failed
    ! 5   - snow/ice (if second byte of GroundPixelQualityFlags >= 50 and <= 130)
    ! 6   - reflectivity < 0 or > 1.0
    ! 7   - bad radiances detected
    ! 8   - aerosol index flag
    ! 9   - radiance PixelQuality error
    ! 10  - radiance PixelQuality warning
    ! 11  - irradiance PixelQuality error
    ! 12  - irradiance PixelQuality warning
    ! 13  - effective surface pressure retrieved because cloud fraction < 0.05
    ! 14  - missing data
    ! 15  - geolocation error
    qflag(1:nx1, :)  = flgbits(1:nx1, :, 0) + flgbits(1:nx1, :, 1) + flgbits(1:nx1, :, 2) &
         + flgbits(1:nx1, :, 3) + flgbits(1:nx1, :, 4) + flgbits(1:nx1, :, 6) + flgbits(1:nx1, :, 7) + &
         flgbits(1:nx1, :, 9) + flgbits(1:nx1, :, 11) + flgbits(1:nx1, :, 13) + &
         flgbits(1:nx1, :, 14) + flgbits(1:nx1, :, 15) 

    ! Fill in cloud top pressure values for bad pixels
    CALL fill_in_ctp(nx1, nt, ctp(1:nx1, 0:nt-1), qflag(1:nx1, 0:nt-1))

    DO ix = 1, nbx
       DO i = 0, nl - 1        
          iix = (ix - 1) * nbin + 1 + xoff0
          ii  = i * nybin + sline
          
          scfr = 0.0; scfr1 = 0.0; scfr0 = 0.0; tmpsum = 0.0
          DO j = iix, iix + nbin - 1
             DO k = ii, ii + nybin - 1
                IF (cfr(j, k) >= 0.0) THEN
                   OMIL2_clouds%cfr(ix, i) = OMIL2_clouds%cfr(ix, i) + cfr(j, k)
                   OMIL2_clouds%ai (ix, i) = OMIL2_clouds% ai(ix, i) + ai(j, k)
                   scfr1 = scfr1 + 1.0
                ENDIF

                IF ( ctp(j, k) > 0.0 .AND. cfr(j, k) >= 0.0 ) THEN                 
                   OMIL2_clouds%ctp(ix, i) = OMIL2_clouds%ctp(ix, i) + LOG(ctp(j, k)) * cfr(j, k)
                   tmpsum = tmpsum + LOG(ctp(j, k)) 
                   scfr = scfr + cfr(j, k); scfr0 = scfr0 + 1.0
                ENDIF
             ENDDO
          ENDDO
             
          IF (scfr /= 0.0) THEN       ! Weighted by Cloud Fraction
             OMIL2_clouds%ctp(ix, i) = EXP(OMIL2_clouds%ctp(ix, i) / scfr)
          ELSE IF (scfr0 > 0.0 ) THEN ! Simple average if cloud fraction is all zero
             OMIL2_clouds%ctp(ix, i) = EXP(tmpsum / scfr0)
          ELSE
             OMIL2_clouds%ctp(ix, i) = 0.0
          ENDIF
          
          IF (scfr1 /= 0.0) THEN
             OMIL2_clouds%cfr(ix, i) = OMIL2_clouds%cfr(ix, i) / scfr1
             OMIL2_clouds%ai (ix, i) = OMIL2_clouds%ai (ix, i) / scfr1
          ELSE
             OMIL2_clouds%cfr(ix, i) = 0.0
             OMIL2_clouds%ai (ix, i) = 0.0
          ENDIF
       ENDDO
    ENDDO
 
    DO ix = 1, nbx
       DO i = 0, nl - 1
          IF (OMIL2_clouds%ctp(ix, i) == 0.0 ) THEN
             OMIL2_clouds%qflags(ix, i) = 10          ! Bad results (should not be used)
          ELSE
             OMIL2_clouds%qflags(ix, i)  = 0          ! Good results
          ENDIF
       ENDDO
    ENDDO

    IF (zoom_mode) THEN
       xoff = NINT(( 1.0 * zoom_p1 - 1 ) / nbin)
       OMIL2_clouds%cfr   (1+xoff:nbx+xoff, 0:nl-1) = OMIL2_clouds%cfr   (1:nbx, 0:nl-1)
       OMIL2_clouds%ctp   (1+xoff:nbx+xoff, 0:nl-1) = OMIL2_clouds%ctp   (1:nbx, 0:nl-1)
       OMIL2_clouds%ai    (1+xoff:nbx+xoff, 0:nl-1) = OMIL2_clouds%ai    (1:nbx, 0:nl-1)
       OMIL2_clouds%qflags(1+xoff:nbx+xoff, 0:nl-1) = OMIL2_clouds%qflags(1:nbx, 0:nl-1)
    ENDIF

   RETURN
 END SUBROUTINE read_omicldrr_clouds

 SUBROUTINE read_omicldo2_clouds (nt, nx, nl, pge_error_status)
    
    ! ----------------------
    ! Input/Output variables
    ! ----------------------  
    INTEGER (KIND=i4), INTENT (IN) :: nt, nx, nl
    INTEGER, INTENT (OUT)          :: pge_error_status
    
    ! -----------------------
    ! Local variables
    ! -----------------------
    INTEGER   (KIND=i4)      :: locerrstat, swath_id, swath_file_id, nt_loc, nx_loc, &
         ix, iix, i, ii, j, k, nline, nbx, nbin, sline, eline, nx1, xoff, xoff0
    REAL      (KIND=r4)      :: scale_cfr, offset_cfr, missval_cfr, scale_ctp, &
         offset_ctp, missval_ctp, scfr, scfr1, scfr0, tmpsum
    CHARACTER (LEN=maxchlen) :: swath_name
    REAL (KIND=r4),    DIMENSION (nx, 0:nt-1)           :: cfr, ctp
    INTEGER (KIND=i2), DIMENSION (nx, 0:nt-1)           :: ctp1
    INTEGER (KIND=i2), DIMENSION (nx, 0:nt-1)           :: qflag
    INTEGER (KIND=i2), DIMENSION (nx, 0:nt-1, 0:nbit-1) :: flgbits

    ! External functions
    INTEGER                  :: OMI_SMF_setmsg
    INTEGER                  :: estat
    
    ! ----------------------
    ! Name of the subroutine
    ! ----------------------
    CHARACTER (LEN=20), PARAMETER :: modulename = 'read_omicldo2_clouds'
    
    locerrstat = pge_errstat_ok
    CALL he5_init_input_file ( l2_cld_filename, swath_name,  swath_id, swath_file_id, &
         nt_loc, nx_loc, locerrstat)
    
    IF ( locerrstat /= pge_errstat_ok ) THEN
       estat = OMI_SMF_setmsg (omsao_e_open_l2cld_file, modulename, 'OML2 Cloud access failed', 0) ; STOP 1
    END IF

    ! check if the dimensions are consistent with OMI L1B data
    !IF (nt_loc /= nt .OR. nx_loc /= nx) THEN
    !  print *, 'OMICLDO2: ', nt_loc, nx_loc
    !  print *, 'OMIL2 : ', nt, nx
    !  print *, l2_cld_filename
    !  WRITE(*, '(A)') 'Inconsistent dimensions between OMIL1B and OMIL2!!!'
    !  !pge_error_status = pge_errstat_error; RETURN
    !ENDIF

    OMIL2_clouds%PGEversion     = '1.0.1.0'
    scale_cfr = 1.0; offset_cfr = 0.0; missval_cfr = 0.0
    scale_ctp = 1.0; offset_ctp = 0.0; missval_ctp = 0
    OMIL2_clouds%CFRmissing = missval_cfr
    OMIL2_clouds%CFRscale   = scale_cfr
    OMIL2_clouds%CFRoffset  = offset_cfr
    OMIL2_clouds%CTPmissing = missval_ctp
    OMIL2_clouds%CTPscale   = scale_ctp
    OMIL2_clouds%CTPoffset  = offset_ctp

    ! Read all the data for the swath
    he5_start_2d = (/ 0, 0 /) ; he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nx, nt /)
    
    ! ----------------------------------------------------------------
    ! Read cloud fraction and cloud top pressure, and check for error.
    ! Eventually we may read the cloud uncertainties also, but for the
    ! first version we stick with just the basic cloud products.
    ! ----------------------------------------------------------------
    locerrstat = HE5_SWrdfld (swath_id, omicld_cfrac_field,   &
         he5_start_2d, he5_stride_2d, he5_edge_2d, cfr(1:nx, 0:nt-1) )
    IF ( locerrstat /= pge_errstat_ok ) THEN
       estat = OMI_SMF_setmsg (omsao_e_open_l2cld_file, modulename, 'OML2 Cloud CFR access failed', 0); STOP 1
    END IF

    locerrstat = HE5_SWrdfld (swath_id, omicld_cpres_field,   &
         he5_start_2d, he5_stride_2d, he5_edge_2d, ctp1(1:nx, 0:nt-1) )
    IF ( locerrstat /= pge_errstat_ok ) THEN
       estat = OMI_SMF_setmsg (omsao_e_open_l2cld_file, modulename, 'OML2 Cloud CTP access failed', 0); STOP 1
    END IF
    ctp = ctp1 * 1.0

    locerrstat = HE5_SWrdfld (swath_id, omicld_qflag_field,   &
         he5_start_2d, he5_stride_2d, he5_edge_2d, qflag(1:nx, 0:nt-1) )
    IF ( locerrstat /= pge_errstat_ok ) THEN
       estat = OMI_SMF_setmsg (omsao_e_open_l2cld_file, modulename, 'OML2 Cloud Quality Flags access failed', 0); STOP 1
    END IF

    !write(*, *) cfr(1:nx, 500)
    !write(*, *) ctp(1:nx, 500)

    IF (zoom_mode) THEN
       nx1 = nx / 2
    ELSE 
       nx1 = nx
    ENDIF
    
    ! ------------------------------------
    ! Save cloud arrays in OMI Cloud Block
    ! ------------------------------------
    cfr(1:nx1,0:nt-1) = cfr(1:nx1,0:nt-1) * OMIL2_clouds%CFRscale + OMIL2_clouds%CFRoffset
    ctp(1:nx1,0:nt-1) = ctp(1:nx1,0:nt-1) * OMIL2_clouds%CTPscale + OMIL2_clouds%CTPoffset

    flgbits = 0

    sline = offset_line;  nline = nl * nybin; eline  = offset_line + nline - 1
    DO ix = 1, nx1
       !CALL convert_2bytes_to_16bits (nbit, nline, qflag(ix, sline:eline), flgbits(ix, sline:eline, 0:nbit-1))
       CALL convert_2bytes_to_16bits (nbit, nt, qflag(ix, 0:nt-1), flgbits(ix, 0:nt-1, 0:nbit-1))
    ENDDO

    nbin = nxbin; IF (coadd_uv2) nbin = nbin * ncoadd
    nbx = nx1 / nbin
    IF (zoom_mode .AND. MOD(zoom_p1, 2) == 0) THEN
       nbx = nbx - 1; xoff0 = 1 
    ELSE
       xoff0 = 0
    ENDIF
    
    OMIL2_clouds%cfr   (1:nbx, 0:nl-1) = 0.0
    OMIL2_clouds%ctp   (1:nbx, 0:nl-1) = 0.0
    OMIL2_clouds%ai    (1:nbx, 0:nl-1) = 0.0
    OMIL2_clouds%qflags(1:nbx, 0:nl-1) = 0

    ! 0   - Solar irradiance warning
    ! 1   - Earth radiance missing flag 
    ! 2   - Earth radiance error flag
    ! 3   - Earth radiance warning flag
    ! 4   - No Snow/Ice Data Flag
    ! 5   - DOAS Fit Error Flag
    ! 6   - DOAS Fit Warning Flag
    ! 7   - Cloud Fraction Missing
    ! 8   - Cloud Fraction Warning
    ! 9   - Cloud Pressure Missing Flag
    ! 10  - Cloud Pressure Warning 
    ! 11  - Extrapolation Error Flag
    ! 12  - 15 Reserved

    !Where does this come from??? 
    !qflag(1:nx1, :)  = flgbits(1:nx1, :, 0) + flgbits(1:nx1, :, 1) + flgbits(1:nx1, :, 2) &
    !     + flgbits(1:nx1, :, 4) + flgbits(1:nx1, :, 6) + flgbits(1:nx1, :, 7) + &
    !     flgbits(1:nx1, :, 9) + flgbits(1:nx1, :, 11) + flgbits(1:nx1, :, 14) + flgbits(1:nx1, :, 15) 

    qflag(1:nx1, :)  = flgbits(1:nx1, :, 1) + flgbits(1:nx1, :, 2) &
         + flgbits(1:nx1, :, 5) + flgbits(1:nx1, :, 7)  &
         + flgbits(1:nx1, :, 9) + flgbits(1:nx1, :, 11) 
        
    ! Fill in cloud top pressure values for bad pixels
    CALL fill_in_ctp(nx1, nt, ctp(1:nx1, 0:nt-1), qflag(1:nx1, 0:nt-1))
 
    DO ix = 1, nbx
       DO i = 0, nl - 1        
          iix = (ix - 1) * nbin + 1 + xoff0
          ii  = i * nybin + sline
          
          scfr = 0.0; scfr1 = 0.0; scfr0 = 0.0; tmpsum = 0.0
          DO j = iix, iix + nbin - 1
             DO k = ii, ii + nybin - 1
                IF (cfr(j, k) >= 0.0) THEN
                   OMIL2_clouds%cfr(ix, i) = OMIL2_clouds%cfr(ix, i) + cfr(j, k)
                   scfr1 = scfr1 + 1.0
                ENDIF

                IF ( ctp(j, k) > 0.0 .AND. cfr(j, k) >= 0.0 ) THEN                 
                   OMIL2_clouds%ctp(ix, i) = OMIL2_clouds%ctp(ix, i) + LOG(ctp(j, k)) * cfr(j, k)
                   tmpsum = tmpsum + LOG(ctp(j, k)) 
                   scfr = scfr + cfr(j, k); scfr0 = scfr0 + 1.0
                ENDIF
             ENDDO
          ENDDO
             
          IF (scfr /= 0.0) THEN       ! Weighted by Cloud Fraction
             OMIL2_clouds%ctp(ix, i) = EXP(OMIL2_clouds%ctp(ix, i) / scfr)
          ELSE IF (scfr0 > 0.0 ) THEN ! Simple average if cloud fraction is all zero
             OMIL2_clouds%ctp(ix, i) = EXP(tmpsum / scfr0)
          ELSE
             OMIL2_clouds%ctp(ix, i) = 0.0
          ENDIF
          
          IF (scfr1 /= 0.0) THEN
             OMIL2_clouds%cfr(ix, i) = OMIL2_clouds%cfr(ix, i) / scfr1
          ELSE
             OMIL2_clouds%cfr(ix, i) = 0.0
          ENDIF
       ENDDO
    ENDDO
   
    DO ix = 1, nbx
       DO i = 0, nl - 1
          IF (OMIL2_clouds%ctp(ix, i) == 0.0 ) THEN
             OMIL2_clouds%qflags(ix, i) = 10          ! Bad results (should not be used)
          ELSE
             OMIL2_clouds%qflags(ix, i)  = 1          ! Good results
          ENDIF
       ENDDO
    ENDDO

    IF (zoom_mode) THEN
       xoff = NINT(( 1.0 * zoom_p1 - 1 ) / nbin)
       OMIL2_clouds%cfr   (1+xoff:nbx+xoff, 0:nl-1) = OMIL2_clouds%cfr   (1:nbx, 0:nl-1)
       OMIL2_clouds%ctp   (1+xoff:nbx+xoff, 0:nl-1) = OMIL2_clouds%ctp   (1:nbx, 0:nl-1)
       OMIL2_clouds%ai    (1+xoff:nbx+xoff, 0:nl-1) = OMIL2_clouds%ai    (1:nbx, 0:nl-1)
       OMIL2_clouds%qflags(1+xoff:nbx+xoff, 0:nl-1) = OMIL2_clouds%qflags(1:nbx, 0:nl-1)
    ENDIF
 
   RETURN
 END SUBROUTINE read_omicldo2_clouds


 SUBROUTINE fill_in_ctp(nx, nt, ctp, qflag)
   ! ----------------------
    ! Input/Output variables
    ! ----------------------  
    INTEGER (KIND=i4), INTENT (IN)                           :: nt, nx
    INTEGER (KIND=i2), DIMENSION (nx, 0:nt-1), INTENT(INOUT) :: qflag
    REAL    (KIND=r4), DIMENSION (nx, 0:nt-1), INTENT(INOUT) :: ctp
    
    ! -----------------------
    ! Local variables
    ! -----------------------
    INTEGER   (KIND=i4)      :: ix, i, j, k, nline, fidx, lidx, sidx, eidx
    REAL      (KIND=r4)      :: frac, dis

    ! If cloud top pressure is too small, then flag that retrieval
    WHERE (ctp < 90.0 ) 
       qflag = qflag + 1
    ENDWHERE
    
    ! If cloud top pressure is too large, then reset it to 1013 mb (i.e., surface pressure)
    WHERE (ctp > 1013.0 ) 
       ctp = 1013.0
    ENDWHERE
    
    ! Reset all flagged pixels to zero 
    WHERE (qflag >= 1) 
       ctp = 0.0
    ENDWHERE
    !PRINT *, 1.D0 * COUNT(ctp == 0.0) / (1.D0 * nx * nt * 1.0)
    
    ! Fill in flagged pixels
    DO i = 1, nt - 2   ! Fill in along the track
       DO ix = 1, nx 
          IF ( ctp(ix, i) == 0.0 .AND. qflag(ix, i-1) == 0 .AND. qflag(ix, i+1) == 0 ) &
               ctp(ix, i) = ( ctp(ix, i-1) + ctp(ix, i+1) ) / 2.0
       ENDDO
    ENDDO
    
    DO i = 0, nt - 1   ! Fill in across the track 
       DO ix = 2, nx - 1 
          IF ( ctp(ix, i) == 0.0 .AND. qflag(ix - 1, i) == 0 .AND. qflag(ix + 1, i) == 0 ) &
               ctp(ix, i) = ( ctp(ix - 1, i) + ctp(ix + 1, i) ) / 2.0
       ENDDO
    ENDDO
    !PRINT *, 1.D0 * COUNT(ctp == 0.0) / (1.D0 * nx * nt * 1.0)
    
    ! Linear interpolation along the track
    DO ix = 1, nx      
       DO i = 0, nt - 1
          IF (ctp(ix, i) > 0.0) EXIT
       ENDDO
       fidx = i
       
       DO i = nt-1, 0, -1
          IF (ctp(ix, i) > 0.0) EXIT
       ENDDO
       lidx = i
       
       IF (fidx >= lidx ) CYCLE
       
       i = fidx + 1
       DO WHILE ( i <= lidx )
          IF (ctp(ix, i) == 0.0 ) THEN
             sidx = i - 1; i = i + 1
             
             eidx = sidx - 1
             DO WHILE (i <= lidx) 
                IF ( ctp(ix, i) > 0.0 ) THEN
                   eidx = i; i = i + 1; EXIT
                ELSE
                   i = i + 1
                ENDIF
             ENDDO
             
             dis = REAL((eidx - sidx), KIND=dp)
             IF (dis <= 12) THEN
                DO j = sidx + 1, eidx - 1
                   frac = 1.0 - REAL((j - sidx), KIND=dp) / dis
                   ctp(ix, j) = frac * ctp(ix, sidx) + (1.0 - frac) * ctp(ix, eidx)
                ENDDO
             ENDIF
          ELSE
             i = i + 1
          ENDIF
       ENDDO
    ENDDO
    !print *, 1.D0 * COUNT(ctp == 0.0) / (1.D0 * nx * nt * 1.0)

    ! Linear interpolation across the track
    DO i = 0, nt - 1     
       DO ix = 1, nx
          IF (ctp(ix, i) > 0.0) EXIT
       ENDDO
       fidx = ix

       DO ix = nx, 1, -1
          IF (ctp(ix, i) > 0.0) EXIT
       ENDDO
       lidx = ix

       IF (fidx >= lidx ) CYCLE
       
       ix = fidx + 1
       DO WHILE ( ix <= lidx )
          IF (ctp(ix, i) == 0.0 ) THEN
             sidx = ix - 1; ix = ix + 1
             
             eidx = sidx  - 1
             DO WHILE (ix <= lidx) 
                IF ( ctp(ix, i) > 0.0 ) THEN
                   eidx = ix; ix = ix + 1; EXIT
                ELSE
                   ix = ix + 1
                ENDIF
             ENDDO
             
             dis = REAL((eidx - sidx), KIND=dp)
             IF (dis <= 6) THEN
                DO j = sidx + 1, eidx - 1
                   frac = 1.0 - REAL((j - sidx), KIND=dp) / dis
                   ctp(j, i) = frac * ctp(sidx, i) + (1.0 - frac) * ctp(eidx, i)
                ENDDO
             ENDIF
          ELSE
             ix = ix + 1
          ENDIF
       ENDDO
    ENDDO
    !print *, 1.D0 * COUNT(ctp == 0.0) / (1.D0 * nx * nt)  

    RETURN

 END SUBROUTINE fill_in_ctp

END MODULE OMSAO_omicloud_module
