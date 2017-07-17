SUBROUTINE omi_pge_fitting_process ( l2_hdf_flag,  pge_error_status )

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: ring_idx, shi_idx, squ_idx,&
       wvl_idx, spc_idx, sig_idx, hwe_idx, hwr_idx, vgr_idx, vgl_idx,  &
       asy_idx, hwl_idx, ops_idx, hwn_idx
  USE OMSAO_variables_module, ONLY: coadd_uv2, pixnum_lim, linenum_lim,npix_fitting, npix_fitted, l2_filename, &
                                    currpix, currline, currloop, n_fitvar_rad, fitvar_rad_saved,mask_fitvar_rad,&
                                    ozabs_convl, so2crs_convl,wavcal,  which_slit, scnwrt, curr_x, curr_y, &
                                    use_solcomp, use_backup,  reduce_resolution, l2_cld_filename, currtime, &
                                    the_lons, the_lats, slit_ins_idx,write_is_slit, the_sza_atm,&
                                    slit_fname, swavcal_fname, rslit_fname, wavcal_fname, slit_rad
  USE ozprof_data_module,      ONLY:calunit,ozwrtint,  l2funit, ozwrtint_fname,ozwrtint_unit, the_cfrac, &
       lcurve_fname, lcurve_write, lcurve_unit,  algorithm_name,algorithm_version, which_cld, num_iter, salbedo, ntp, avg_kernel
  USE OMSAO_pixelcorner_module
  USE OMSAO_omicloud_module
  USE OMSAO_slitfunction_module
  USE OMSAO_omidata_module,    ONLY: nlines_max, ntimes, ntimes_loop, nxtrack, nfxtrack, &
       omi_radpix_errstat,omi_solpix_errstat, omi_exitval, omi_fitvar, omi_initval,&
       ncoadd, nxbin, nybin, offset_line, zoom_mode, zoom_p1, zoom_p2
  USE he5_output_module
  USE OMSAO_errstat_module
  IMPLICIT NONE

  ! -----------------------
  ! Input/Output variables
  ! ----------------------
  INTEGER, INTENT (IN)  :: l2_hdf_flag
  INTEGER, INTENT (OUT) :: pge_error_status

  ! -------------------------
  ! Local variables (for now)
  ! -------------------------
  INTEGER :: initval, errstat,exval
  INTEGER :: first_line, last_line,first_pix, last_pix,curr_fitted_line,  nxcoadd,   &
              sline, eline, i,iline
  REAL (KIND=dp), DIMENSION(3)    :: fitcol
  REAL (KIND=dp), DIMENSION(3, 2) :: dfitcol
  REAL (KIND=dp)     :: fitcol_avg, rms_avg, dfitcol_avg, drel_fitcol_avg, rms,avgres
  REAL (KIND=dp), DIMENSION(24) :: dfs 
  CHARACTER (LEN=7)  :: cldtype, tmp_zoom
  LOGICAL            :: reduce_resolution_save
  INTEGER, SAVE :: slit_unit
  INTEGER, PARAMETER          :: n_slitvar = 6
  CHARACTER (len=100), DIMENSION(n_slitvar),  PARAMETER ::slitname=(/'shi','squ','hwle','asym','ops','hwn'/)
  INTEGER, DIMENSION (n_slitvar), PARAMETER :: slitidx=(/shi_idx, squ_idx,hwe_idx,asy_idx, ops_idx, hwn_idx/)
  CHARACTER (LEN=100) :: fname
 
  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=23), PARAMETER :: modulename = 'omi_pge_fitting_process'
  
  pge_error_status = pge_errstat_ok


  ! Initialize swath names 
  CALL omi_set_fitting_parameters ( pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN
 
  ! Initialize ntimes, nxtrack, nfxtrack (UV-1)
  CALL omi_read_radiance_paras (pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN
        
  ! xtrack positions to be processed (start from the actual binned position)
  IF (linenum_lim(2) >= ntimes)  linenum_lim(2) = ntimes

  IF (zoom_mode) THEN
     pixnum_lim(1) = MAX(NINT(1.0 * (zoom_p1 + 1) / ncoadd), pixnum_lim(1))
     pixnum_lim(2) = MIN(zoom_p2 / ncoadd, pixnum_lim(2))
     IF ( MOD(pixnum_lim(2) - pixnum_lim(1) + 1, nxbin) /= 0) THEN
        WRITE(*, '(A,2I4)') 'Incorrect across track binning option: ', pixnum_lim(1:2)
        pge_error_status = pge_errstat_error; RETURN
     ENDIF

     IF (pixnum_lim(1) > pixnum_lim(2)) THEN
        WRITE(*, *) 'This is a zoom in mode orbit!!!'
        WRITE(*, '(A, I2, A, I2)') 'Invalid pixel selection, must be between ', &
             (zoom_p1 / ncoadd) * ncoadd + 1, ' and ', (zoom_p2 / ncoadd ) * ncoadd 
        pge_error_status = pge_errstat_error; RETURN
     ENDIF  
  ENDIF

  first_pix  = CEILING(1.0 * pixnum_lim(1) / nxbin)
  last_pix  = NINT(1.0 * pixnum_lim(2) / nxbin )
  nxcoadd = nxbin * ncoadd

  ! Line number, starting from zero and keep track of offset
  offset_line = linenum_lim(1) - 1; first_line = 1
  last_line   = INT ((linenum_lim(2) - linenum_lim(1) + 1.0) / nybin) 

  ! load omi slit parameters  (Need to use it while getting coadded irradiance data)
  IF (which_slit == slit_ins_idx) THEN
     CALL load_slitpars (pge_error_status)
     IF ( pge_error_status >= pge_errstat_error ) RETURN
     IF (scnwrt) WRITE(*, '(A)') 'Finish loading OMI slit parameters !!!'
  ENDIF

  ! Read OMI irradiances
  reduce_resolution_save = reduce_resolution; reduce_resolution = .FALSE.
  CALL omi_read_irradiance_data (calunit, nxcoadd, first_pix, last_pix, pge_error_status ) 
  IF ( pge_error_status >= pge_errstat_error ) THEN
     use_backup = .TRUE.
     CALL omi_read_irradiance_data (calunit, nxcoadd, first_pix, last_pix, pge_error_status )  
     IF ( pge_error_status >= pge_errstat_error ) RETURN   
  ENDIF
  IF (scnwrt) WRITE(*, '(A)') 'Finish reading irradiances!!!'
  
  IF (use_solcomp) THEN
     CALL replace_solar_irradiance(calunit, nxcoadd, first_pix, last_pix, pge_error_status)
     IF ( pge_error_status >= pge_errstat_error ) RETURN
     IF (scnwrt) WRITE(*, '(A)') 'Finish replacing irradiances with solar composite!!!'
  ENDIF

  ! Calibrate OMI irradiances for all cross track pixels and save results
  IF (reduce_resolution_save) reduce_resolution = .TRUE.

  CALL omi_irrad_cross_calibrate (first_pix, last_pix, pge_error_status)
   
  IF ( pge_error_status >= pge_errstat_error ) RETURN
  IF (scnwrt) WRITE(*, '(A)') 'Finish calibrating irradiances!!!'
  IF (write_is_slit) THEN 
       slit_unit = 1004
       fname = slit_fname
       if (wavcal_sol) fname = swavcal_fname
         OPEN (slit_unit,file=TRIM(ADJUSTL(fname)), status='unknown')
        WRITE(slit_unit,'(5i5)') which_slit, last_pix - first_pix + 1 , 0, n_slitvar
      CALL omps_write_cali (1, slit_unit, n_slitvar, slitidx, slitname, first_pix, last_pix, first_line, last_line)
      IF (wavcal .or. slit_rad ) THEN
      ELSE
          WRITE(*,*)'Save slit/wav calibration in',trim(fname)
          STOP
      ENDIF
      CLOSE(slit_unit)
  ENDIF
   
  
  IF (reduce_resolution_save) THEN
     CALL omi_read_irradiance_data (calunit, nxcoadd, first_pix, last_pix, pge_error_status ) 
     IF ( pge_error_status >= pge_errstat_error ) RETURN  
     IF (scnwrt) WRITE(*, '(A)') 'Finish reading irradiances (with reduced resolution)!!!'      
  ENDIF
   
       
  ! Compute spatial pixel corners and effective viewing geometry
  CALL compute_pixel_corners ( ntimes, nfxtrack, last_line, pge_error_status)

  IF ( pge_error_status >= pge_errstat_error ) RETURN
  IF (scnwrt) WRITE(*, '(A)') 'Finish computing pixel corners!!!'

  ! Determine OMCLDRR or OMCLDO2
  i = INDEX(l2_cld_filename, '-o') - 22
  cldtype = l2_cld_filename(i : i + 6)
  IF (cldtype == 'OMCLDO2') THEN
     IF (which_cld == 0) which_cld = 1
     IF (which_cld == 3) which_cld = 4
  ENDIF

  IF (which_cld == 4 .OR. which_cld == 1) THEN
     CALL read_omicldo2_clouds (ntimes, nxtrack, last_line, pge_error_status)
     IF ( pge_error_status >= pge_errstat_error ) RETURN
     IF (scnwrt) WRITE(*, '(A)') 'Finish reading omi L2 o2 clouds!!!'
  ELSE IF (which_cld == 3 .OR. which_cld == 0) THEN
     CALL read_omicldrr_clouds (ntimes, nxtrack, last_line, pge_error_status)
     IF ( pge_error_status >= pge_errstat_error ) RETURN
     IF (scnwrt) WRITE(*, '(A)') 'Finish reading omi L2 RR clouds!!!'
  ENDIF
  
  ! Initialize fitting statistics
  npix_fitting = 0       ! number of pixels (failure + success)
  npix_fitted  = 0       ! number of successfully fitted pixels   
  fitcol_avg   = 0.0;  rms_avg = 0.0; dfitcol_avg = 0.0; drel_fitcol_avg = 0.0

  ! Open output files
  IF (lcurve_write) THEN
     OPEN(UNIT=lcurve_unit, file=TRIM(ADJUSTL(lcurve_fname)), status='unknown', IOSTAT=errstat)
     IF ( errstat /= pge_errstat_ok ) THEN
        WRITE(*, *) modulename, ': Cannot open lcurve file!!!'
        pge_error_status = pge_errstat_error; RETURN
     END IF
  ENDIF     
  IF (ozwrtint) THEN
     OPEN(UNIT=ozwrtint_unit, file=TRIM(ADJUSTL(ozwrtint_fname)), status='unknown', IOSTAT=errstat)
     IF ( errstat /= pge_errstat_ok ) THEN
        WRITE(*, *) modulename, ': Cannot open intermediate output file!!!'
        pge_error_status = pge_errstat_error; RETURN
     END IF
  ENDIF
  
  IF (l2_hdf_flag == 0) THEN
     OPEN (UNIT=l2funit, FILE=TRIM(ADJUSTL(l2_filename)), STATUS='UNKNOWN', IOSTAT=errstat)

     IF ( errstat /= pge_errstat_ok ) THEN
        WRITE(*, *) modulename, ': Cannot open output file!!!'
    print * , trim(adjustl(l2_filename))
        pge_error_status = pge_errstat_error; RETURN
     END IF
     CALL timestamp(currtime)
     WRITE(l2funit, '(3A,1x,A27,A10,I5,A10,I5)') TRIM(ADJUSTL(algorithm_name)), ', ', &
          TRIM(ADJUSTL(algorithm_version)), currtime, ' xbin = ', nxbin, ' ybin = ', nybin
  ELSE 
     CALL He5_L2WrtInit (first_pix, last_pix, first_line, last_line, errstat )
     IF ( errstat /= pge_errstat_ok ) THEN
        WRITE(*, *) modulename, ' : Cannot create HE5 output file!!!'
        pge_error_status = pge_errstat_error; RETURN
     END IF

     CALL He5_L2SetGeoFields (first_pix, last_pix, errstat )
     IF ( errstat /= pge_errstat_ok ) THEN
        WRITE(*, *) modulename, ' : Cannot write geolocation fields!!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
  ENDIF
  
  ! loop through each OMI data block
  ! 1. read each block
  ! 2. perform calibration for each block (middle line)
  ! 3. perform retrievals from first_pix (all lines within the block) to last_pix

  ! these pixels with exitval >= 0 will be used by subsequent retrievals as initial values
  omi_exitval = -10 ! no retrievals yet 
  omi_fitvar  = 0.0
  WRITE(*, '(2a5, 6a8,2a5)') 'LINE','XPix', 'Lat','SZA','ALB','CFR','RMS','avgres', '#','exval'


  OMIBlock: DO iline = 0, last_line-1, nlines_max
     ntimes_loop = nlines_max

     IF ( iline + ntimes_loop > last_line )   ntimes_loop = last_line - iline

     ! Actually lines in OMI data
     sline = offset_line + iline * nybin
     eline = sline + ntimes_loop * nybin - 1     
  
     ! Get NTIMES_LOOP radiance lines (with effective viewing geometry)
     
     CALL omi_read_radiance_lines (iline, ntimes_loop, sline, eline, first_pix, last_pix, nxcoadd, pge_error_status)
     IF ( pge_error_status >= pge_errstat_error ) RETURN
     IF (scnwrt) WRITE(*, '(A,I4,A,I4)') 'Finishing reading radiances for lines: ', sline + 1, ' - ', eline+1

     ! Perform calibration for the middle scan line and apply to the other scan lines
    
     IF (wavcal) THEN
        CALL omi_rad_cross_calibrate (first_pix, last_pix, pge_error_status)
        IF ( pge_error_status >= pge_errstat_error ) RETURN
        IF (scnwrt) WRITE(*, '(A,I4,A,I4)') 'Finishing calibrating radiances for lines: ', sline+1, ' - ', eline+1
     ENDIf
        
      
     ! loop through each xtrack position
     ! 1. Prepare databases, adjust radiances
     ! 2. Process all pixels at this poistion
        
     XtrackPix: DO currpix = first_pix, last_pix
        ! Need to convolve high-resolution ozone absorption cross section (for this position)
        ! Once the xsection is convolved, it will be set to false in ROUTINE getabs_crs
        ozabs_convl = .TRUE.; so2crs_convl = .TRUE.


        ! Load/adjust irradiances and slit calibration parameters
        CALL omi_adj_solar_data (pge_error_status)
        IF ( pge_error_status >= pge_errstat_error ) CYCLE

        curr_fitted_line = 0
             
        YfitLine: DO currloop = 0, ntimes_loop - 1

           omi_exitval(currpix, currloop) = -10
           currline = iline + currloop  
   
           IF (omi_radpix_errstat(currpix, currloop) == pge_errstat_error .OR. &
                omi_solpix_errstat(currpix) == pge_errstat_error ) &
                omi_exitval(currpix, currloop) = -9
           ! print *,'pge:', omi_radpix_errstat(currpix, currloop), omi_solpix_errstat(currpix)
           ! Load/adjust radiances/geolocations fields for a particular pixel
           ! Prepare databases for the first pixel (ifitline == 1)   
           IF (omi_exitval(currpix, currloop) == -10) THEN
              CALL omi_adj_earthshine_data (curr_fitted_line, pge_error_status)
              IF ( pge_error_status >= pge_errstat_error ) omi_exitval(currpix, currloop) = -9
           ENDIF
               curr_x = (currpix-1) * nxbin  + 1 
               curr_y = currline * nybin + offset_line + 1
           IF (omi_exitval(currpix, currloop) == -10) THEN
              IF (scnwrt) WRITE(*, '(A,I5,A10,I5, A10, I5)')       'OMI Pixel: Line = ', &
                   curr_y, ' XPix = ', curr_x, &
                   ' Loop = ', currloop
              initval = omi_initval(currpix, currloop)
              CALL specfit_ozprof (initval, fitcol, dfitcol, rms, exval)
              omi_exitval(currpix, currloop) = exval  ! Sotre exit status for current pixel
              omi_fitvar(currpix, currloop, 1:n_fitvar_rad) &
                                            = fitvar_rad_saved(mask_fitvar_rad(1:n_fitvar_rad))
           ELSE
              exval = -9
           ENDIF
           !CALL timestamp(currtime)
             
           ! Write retrievals
           IF (l2_hdf_flag == 0) THEN
              IF (exval > -9) CALL omi_write_intermed (l2funit, fitcol, dfitcol,  rms, avgres, exval)
           ELSE 
              CALL He5_L2SetDataFields (currpix, first_pix, last_pix, currloop, currline, &
                   ntimes_loop, exval, fitcol, dfitcol, pge_error_status )
              IF ( pge_error_status >= pge_errstat_error ) RETURN
           ENDIF
              
           WRITE(*, '(2I5,6f8.2,2i5)') &
                    curr_y, curr_x, the_lats(5),the_sza_atm, salbedo, the_cfrac, rms, avgres,num_iter, exval

           IF ( exval >= 0 .AND. fitcol(1) > 0.0 .AND. dfitcol(1, 1) >= 0.0 ) THEN         
              ! ----------------------------------------------------------------------
              ! Some general statistics on the average fitted column and uncertainty.
              ! Again, we make sure that only "good" fits are included in the average.
              ! ----------------------------------------------------------------------              
              fitcol_avg       = fitcol_avg + fitcol(1)
              rms_avg          = rms_avg + rms
              dfitcol_avg      = dfitcol_avg + dfitcol(1, 1)
              drel_fitcol_avg  = drel_fitcol_avg + dfitcol(1, 1) / fitcol(1)
              npix_fitted      = npix_fitted + 1
              curr_fitted_line = curr_fitted_line + 1
           ENDIF

           npix_fitting        = npix_fitting + 1          
        ENDDO YfitLine
     ENDDO XtrackPix
  ENDDO OMIBlock
  
  ! ------------
  ! Final output
  ! ------------
  IF ( npix_fitted == 0) npix_fitted = 1
  IF (scnwrt) CALL write_final(fitcol_avg, rms_avg, dfitcol_avg, drel_fitcol_avg, npix_fitted)
  IF (scnwrt) WRITE(*, '(2(A,I5))') 'Number of pixels = ', &
       npix_fitting, '   Number of fitted pixels = ', npix_fitted

  ! -----------------------------------------
  ! Close L1 radiance file and L2 output file
  ! -----------------------------------------
  IF (l2_hdf_flag == 0) THEN
     CALL timestamp(currtime)
     WRITE(l2funit, '(A27)') currtime 
     CLOSE ( l2funit )
  ENDIF

  IF (lcurve_write) CLOSE (lcurve_unit)
  IF (ozwrtint)     CLOSE (ozwrtint_unit)

  RETURN
END SUBROUTINE omi_pge_fitting_process
