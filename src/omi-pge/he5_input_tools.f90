SUBROUTINE he5_init_input_file ( &
     file_name, swath_name, swath_id, swath_file_id, ntimes_aux, nxtrack_aux, he5stat )

  !------------------------------------------------------------------------------
  ! This subroutine initializes the HE5 input files for for
  ! prefitted BrO and O3 that are used in OMHCHO
  !
  ! Input:
  !   file_name         - Name of HE5 input file
  !   swath_name        - Name of existing swath in file
  !   swath_file_id_inp - id number for HE5 input file (required for closing it)
  !   swath_id_inp      - id number for swath (required for reading from swath)
  !
  ! Output:
  !   ntimes_aux.............nTimes  as given in product file
  !   nxtrack_aux............nXtrack as given in product file
  !   he5stat................OMI_E_SUCCESS if everything went well
  !
  ! No Swath ID Variables passed through MODULE.
  !
  !------------------------------------------------------------------------------

  USE OMSAO_indices_module,    ONLY: pge_hcho_idx, pge_bro_idx, max_calfit_idx
  USE OMSAO_parameters_module, ONLY: maxchlen, i4_missval, r4_missval, r8_missval, vb_lev_default
  USE OMSAO_he5_module
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=19), PARAMETER :: modulename = 'he5_init_input_file'
 
  ! ---------------
  ! Input variables
  ! ---------------
  CHARACTER (LEN=*), INTENT (IN) :: file_name

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER   (KIND=i4),      INTENT (INOUT) :: he5stat
  INTEGER   (KIND=i4),      INTENT (OUT)   :: swath_id, swath_file_id, nxtrack_aux, ntimes_aux
  CHARACTER (LEN=maxchlen), INTENT (OUT)   :: swath_name

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4)                      :: ndim, nsep
  INTEGER   (KIND=i4)                      :: i, errstat, swlen, iend, istart, astat
  INTEGER   (KIND=i4),      DIMENSION(0:9) :: dim_array, dim_seps
  CHARACTER (LEN=maxchlen)                 :: dim_chars, tmp_char

  REAL (KIND=r8), DIMENSION (:,:,:), ALLOCATABLE :: o3fit_cols, o3fit_dcols
  REAL (KIND=r8), DIMENSION (:,:),   ALLOCATABLE :: tmp_col

  errstat = pge_errstat_ok

  ! -----------------------------------------------------------
  ! Open HE5 output file and check SWATH_FILE_ID ( -1 if error)
  ! -----------------------------------------------------------
  swath_file_id = HE5_SWopen ( TRIM(ADJUSTL(file_name)), he5f_acc_rdonly )
  IF ( swath_file_id == he5_stat_fail ) THEN
     CALL error_check ( &
          0, 1, pge_errstat_error, OMSAO_E_HE5SWOPEN, modulename, vb_lev_default, he5stat )
     IF ( he5stat >= pge_errstat_error ) RETURN
  END IF

  ! ---------------------------------------------
  ! Check for existing HE5 swath and attach to it
  ! ---------------------------------------------
  errstat  = HE5_SWinqswath  ( TRIM(ADJUSTL(file_name)), swath_name, swlen )
  swath_id = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(swath_name)) )
  IF ( swath_id == he5_stat_fail ) THEN
     CALL error_check ( &
          0, 1, pge_errstat_error, OMSAO_E_HE5SWATTACH, modulename, vb_lev_default, he5stat )
     IF ( he5stat >= pge_errstat_error ) RETURN
  END IF
  
  ! ------------------------------
  ! Inquire about swath dimensions
  ! ------------------------------
  errstat  = HE5_SWinqdims  ( swath_id, dim_chars, dim_array(0:9) )

  ! -------------------------------------------------------------
  ! The number of recorded dimension is the number of non-zero
  ! entries - i.a.w., it't the index of the first ZERO minus ONE.
  ! -------------------------------------------------------------
  ndim = MINVAL( MINLOC( dim_array, MASK=(dim_array == 0) ) ) - 1

  ! -----------------------------------------
  ! Extract nTimes and nXtrack from the swath
  ! -----------------------------------------
  nxtrack_aux = 0  ;  ntimes_aux = 0
  dim_chars = TRIM(ADJUSTL(dim_chars))
  swlen = LEN_TRIM(ADJUSTL(dim_chars))
  istart = 1  ;  iend = 1

  ! ----------------------------------------------------------------------
  ! Find the positions of separators (commas, ",") between the dimensions.
  ! Add a "pseudo separator" at the end to fully automate the consecutive
  ! check for nTimes and nXtrack.
  ! ----------------------------------------------------------------------
  nsep = 0 ; dim_seps = 0 ; dim_seps(0) = 0
  getseps: DO i = 1, swlen 
     IF ( dim_chars(i:i) == ',' ) THEN
        nsep = nsep + 1
        dim_seps(nsep) = i
     END IF
  END DO getseps
  nsep = nsep + 1 ; dim_seps(nsep) = swlen+1

  ! --------------------------------------------------------------------
  ! Hangle along the NSEP indices until we have found the two dimensions
  ! we are interested in.
  ! --------------------------------------------------------------------
  getdims:DO i = 0, nsep-1
     istart = dim_seps(i)+1 ; iend = dim_seps(i+1)-1
     IF  ( dim_chars(istart:iend) == "nTimes"  ) ntimes_aux  = dim_array(i)
     IF  ( dim_chars(istart:iend) == "nXtrack" ) nxtrack_aux = dim_array(i)
     IF ( ntimes_aux > 0 .AND. nxtrack_aux > 0 ) EXIT getdims
  END DO getdims

  RETURN
END SUBROUTINE he5_init_input_file

!SUBROUTINE he5_read_prefit_columns ( &
!     swath_id, ntimes_mol, nxtrack_mol, iline, &
!     mcol_len,  molcol_field, col_mol, mdcol_len, &
!     moldcol_field, dcol_mol, yn_read_amf, errstat )
!
!  USE OMSAO_precision_module
!  USE OMSAO_indices_module,    ONLY: pge_hcho_idx, pge_bro_idx
!  USE OMSAO_parameters_module, ONLY: maxchlen
!  USE OMSAO_omidata_module
!  USE OMSAO_he5_module
!  USE OMSAO_errstat_module
!
!  IMPLICIT NONE
!
!  ! ------------------------------
!  ! Name of this module/subroutine
!  ! ------------------------------
!  CHARACTER (LEN=23), PARAMETER :: modulename = 'he5_read_prefit_columns'
!
!  ! ---------------
!  ! Input variables
!  ! ---------------
!  INTEGER   (KIND=i4),       INTENT (IN) :: iline, swath_id, ntimes_mol, nxtrack_mol
!  INTEGER   (KIND=i4),       INTENT (IN) :: mcol_len, mdcol_len
!  CHARACTER (LEN=mcol_len) , INTENT (IN) :: molcol_field
!  CHARACTER (LEN=mdcol_len), INTENT (IN) :: moldcol_field
!  LOGICAL,                   INTENT (IN) :: yn_read_amf
!
!  ! ---------------
!  ! Output variable
!  ! ---------------
!  INTEGER (KIND=i4),                                          INTENT (INOUT) :: errstat
!  REAL    (KIND=r8), DIMENSION (nxtrack_mol, 0:ntimes_mol-1), INTENT (OUT)   :: col_mol, dcol_mol
!
!  ! ---------------
!  ! Local variables
!  ! ---------------
!  INTEGER (KIND=i4) :: locerrstat, i
!  REAL    (KIND=r8), DIMENSION (nxtrack_mol, 0:ntimes_mol-1) :: amf
!
!
!  locerrstat = pge_errstat_ok
!
!  ! ----------------------------------------------------
!  ! Read current data block fitting output from HE5 file
!  ! ----------------------------------------------------
!  he5_start_2d = (/ 0, iline /) ; he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nxtrack_mol, ntimes_mol /)
!
!  ! -----------------------------
!  ! Column amount and uncertainty
!  ! -----------------------------
!  locerrstat = HE5_SWrdfld ( swath_id, TRIM(ADJUSTL(molcol_field)),     &
!       he5_start_2d, he5_stride_2d, he5_edge_2d, col_mol(1:nxtrack_mol,0:ntimes_mol-1) )
!  locerrstat = HE5_SWrdfld ( swath_id, TRIM(ADJUSTL(moldcol_field)),    &
!       he5_start_2d, he5_stride_2d, he5_edge_2d, dcol_mol(1:nxtrack_mol,0:ntimes_mol-1) )
!
!
!  ! --------------------------------------
!  ! Air Mass Factor, but only if requested
!  ! --------------------------------------
!  IF ( yn_read_amf ) THEN
!     locerrstat = HE5_SWrdfld ( swath_id, TRIM(ADJUSTL(amfmol_field)),     &
!          he5_start_2d, he5_stride_2d, he5_edge_2d, amf(1:nxtrack_mol,0:ntimes_mol-1) )
!     WHERE ( amf > 0.0_r8 )
!        col_mol  =  col_mol * amf
!        dcol_mol = dcol_mol * amf
!     END WHERE
!  END IF
!
!  ! ------------------
!  ! Check error status
!  ! ------------------
!  CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, &
!       OMSAO_E_HE5SWRDFLD, modulename, vb_lev_default, errstat )
!
!  RETURN
!END SUBROUTINE he5_read_prefit_columns

!FUNCTION he5_close_prefit_file ( swath_id, swath_file_id ) RESULT ( he5stat )
!
!  !------------------------------------------------------------------------------
!  ! This function detatches from the HE5 swath and closes the HE5 input file.
!  !
!  ! Input: 
!  !
!  !   swath_id       - ID of pre-fitted input swath      (for BrO or O3)
!  !   swath_file_id  - ID of pre-fitted input swath file (for BrO or O3)
!  !
!  !------------------------------------------------------------------------------
!
!  USE OMSAO_variables_module,  ONLY: verb_thresh_lev
!  USE OMSAO_he5_module
!  USE OMSAO_errstat_module
!
!  IMPLICIT NONE
!
!  ! ---------------------------------------
!  ! Name of this module/subroutine/function
!  ! ---------------------------------------
!  CHARACTER (LEN=21), PARAMETER :: modulename = 'he5_close_prefit_file'
!
!  ! ---------------
!  ! Input variables
!  ! ---------------
!  INTEGER (KIND=i4), INTENT (IN) :: swath_id, swath_file_id
!  
!  ! ---------------
!  ! Result variable
!  ! ---------------
!  INTEGER (KIND=i4) :: he5stat
!  
!  ! --------------
!  ! Local variable
!  ! --------------
!  INTEGER (KIND=i4) :: locerr
!
!
!  he5stat = pge_errstat_ok
!  locerr  = pge_errstat_ok
!
!  ! -----------------------------------------------
!  ! Detach from HE5 swath and close HE5 output file
!  ! -----------------------------------------------
!  locerr = HE5_SWdetach ( swath_id )
!  locerr = HE5_SWclose  ( swath_file_id )
!  CALL error_check ( locerr, HE5_STAT_OK, pge_errstat_warning, &
!       OMSAO_W_HE5SWCLOSE, modulename, vb_lev_default, he5stat )
!
!  RETURN
!END FUNCTION he5_close_prefit_file

!SUBROUTINE init_prefit_files (                                                   &
!     ntimes, nxtrack, yn_o3_prefit, yn_bro_prefit,                               &
!     o3fit_swath_id, o3fit_swath_file_id, brofit_swath_id, brofit_swath_file_id, &
!     errstat )
!
!  USE OMSAO_precision_module
!  USE OMSAO_variables_module, ONLY: o3_prefit_fname, bro_prefit_fname
!  USE OMSAO_errstat_module
!  USE OMSAO_he5_module,       ONLY: o3fit_swath_name, brofit_swath_name
!
!  IMPLICIT NONE
!
!  ! ---------------
!  ! Input variables
!  ! ---------------
!  INTEGER (KIND=i4), INTENT (IN) :: ntimes, nxtrack
!  ! --------------------------------
!  ! (Possible modification on error)
!  ! --------------------------------
!  LOGICAL, INTENT (INOUT) :: yn_o3_prefit, yn_bro_prefit
!
!  ! ----------------
!  ! Output variables
!  ! ----------------
!  INTEGER (KIND=i4), INTENT (OUT) :: &
!       o3fit_swath_id, o3fit_swath_file_id, brofit_swath_id, brofit_swath_file_id
!
!  ! ------------------
!  ! Modified variables
!  ! ------------------
!  INTEGER (KIND=i4), INTENT (INOUT) :: errstat
!
!  ! ---------------
!  ! Local variables
!  ! ---------------
!  INTEGER (KIND=i4) :: locerrstat, ntimes_o3, nxtrack_o3, ntimes_bro, nxtrack_bro
!  CHARACTER (LEN=17), PARAMETER :: modulename = 'init_prefit_files'
!
!  ! ---------------------------
!  ! Initialize output variables
!  ! ---------------------------
!  o3fit_swath_id  = -1 ; o3fit_swath_file_id  = -1
!  brofit_swath_id = -1 ; brofit_swath_file_id = -1
!
!  IF ( .NOT. ANY((/yn_o3_prefit, yn_bro_prefit/)) ) RETURN
!
!  ! ----------
!  ! O3 prefits
!  ! ----------
!  locerrstat = pge_errstat_ok
!  IF ( yn_o3_prefit ) THEN
!     CALL he5_init_input_file ( &
!          o3_prefit_fname, o3fit_swath_name, o3fit_swath_id, o3fit_swath_file_id, &
!          ntimes_o3, nxtrack_o3, errstat )
!     IF ( ntimes_o3 /= ntimes .OR. nxtrack_o3 /= nxtrack ) THEN
!        locerrstat = pge_errstat_error
!        CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_fatal, OMSAO_E_PREFITDIM, &
!             modulename//f_sep//"O3 access failed.", vb_lev_default, errstat )
!        yn_o3_prefit = .FALSE.
!     END IF
!  END IF
!
!  ! -----------
!  ! BrO prefits
!  ! -----------
!  locerrstat = pge_errstat_ok
!  IF ( yn_bro_prefit ) THEN
!     CALL he5_init_input_file ( &
!          bro_prefit_fname, brofit_swath_name, brofit_swath_id, &
!          brofit_swath_file_id, ntimes_bro, nxtrack_bro, locerrstat )
!     IF ( ntimes_bro /= ntimes .OR. nxtrack_bro /= nxtrack ) THEN
!        locerrstat = pge_errstat_error
!        CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_fatal, OMSAO_E_PREFITDIM, &
!             modulename//f_sep//"BrO access failed.", vb_lev_default, errstat )
!        yn_bro_prefit = .FALSE.
!     END IF
!  END IF
!
!  RETURN
!END SUBROUTINE init_prefit_files

!SUBROUTINE read_prefit_columns (                                      &
!     ntimes, nxtrack, nloop, nstart,                                  &
!     yn_o3_prefit,  o3fit_swath_id,  o3_prefit_col,  o3_prefit_dcol,  &
!     yn_bro_prefit, brofit_swath_id, bro_prefit_col, bro_prefit_dcol, &
!     errstat )
!
!  USE OMSAO_precision_module
!  USE OMSAO_errstat_module
!  USE OMSAO_indices_module,    ONLY: o3_t1_idx, o3_t2_idx, o3_t3_idx, bro_idx
!  USE OMSAO_parameters_module, ONLY: r8_missval
!  USE OMSAO_variables_module,  ONLY: o3_prefit_fname, bro_prefit_fname, refspecs_original
!  USE OMSAO_he5_module,        ONLY: o3_prefit_fields
!
!  IMPLICIT NONE
!
!  ! ---------------
!  ! Input variables
!  ! ---------------
!  LOGICAL,           INTENT (IN) :: yn_o3_prefit, yn_bro_prefit
!  INTEGER (KIND=i4), INTENT (IN) :: ntimes, nxtrack, nloop, nstart, o3fit_swath_id, brofit_swath_id
!
!  ! ----------------
!  ! Output variables
!  ! ----------------
!  REAL (KIND=r8), DIMENSION (o3_t1_idx:o3_t3_idx,nxtrack,0:nloop-1), INTENT (OUT) :: &
!       o3_prefit_col, o3_prefit_dcol
!  REAL (KIND=r8), DIMENSION (nxtrack,0:nloop-1),                     INTENT (OUT) :: &
!       bro_prefit_col, bro_prefit_dcol
!
!  ! ------------------
!  ! Modified variables
!  ! ------------------
!  INTEGER (KIND=i4), INTENT (INOUT) :: errstat
!
!  ! ---------------
!  ! Local variables
!  ! ---------------
!  INTEGER (KIND=i4) :: locerrstat, i
!  LOGICAL           :: yn_read_amf
!  CHARACTER (LEN=12), PARAMETER ::  col_str = "ColumnAmount"
!  CHARACTER (LEN=17), PARAMETER :: dcol_str = "ColumnUncertainty"
!  INTEGER  (KIND=i4), PARAMETER :: lcolstr = LEN(col_str), ldcolstr = LEN(dcol_str)
!  
!  CHARACTER (LEN=19), PARAMETER :: modulename = 'read_prefit_columns'
!
!
!  ! ---------------------------------------------
!  ! O3 prefitted columns and column uncertainties
!  ! ---------------------------------------------
!  yn_read_amf = .FALSE. ; locerrstat = pge_errstat_ok
!  IF ( yn_o3_prefit ) THEN
!     o3_prefit_col = 0.0_r8  ;  o3_prefit_dcol = 0.0_r8
!     DO i = o3_t1_idx, o3_t3_idx
!        CALL he5_read_prefit_columns ( &
!             o3fit_swath_id, nloop, nxtrack, nstart, &
!             LEN_TRIM(ADJUSTL(o3_prefit_fields(i,1))), TRIM(ADJUSTL(o3_prefit_fields(i,1))), &
!             o3_prefit_col (i,1:nxtrack,0:nloop-1),  &
!             LEN_TRIM(ADJUSTL(o3_prefit_fields(i,2))), TRIM(ADJUSTL(o3_prefit_fields(i,2))), &
!             o3_prefit_dcol(i,1:nxtrack,0:nloop-1), &
!             yn_read_amf, locerrstat )
!        errstat = MAX( errstat, locerrstat )
!
!        ! ----------------------------------------------------------------------
!        ! Multiply O3 columns with normalization factor to return to true values
!        ! ----------------------------------------------------------------------
!        WHERE ( o3_prefit_col (i,1:nxtrack,0:nloop-1) > r8_missval )
!           o3_prefit_col (i,1:nxtrack,0:nloop-1) = &
!                o3_prefit_col (i,1:nxtrack,0:nloop-1) * refspecs_original(i)%NormFactor
!           o3_prefit_dcol(i,1:nxtrack,0:nloop-1) = &
!                o3_prefit_dcol(i,1:nxtrack,0:nloop-1) * refspecs_original(i)%NormFactor
!        END WHERE
!
!     END DO
!  END IF
!
!  ! -----------------------------------------------
!  ! BrO prefitted columns and column uncertainties
!  ! -----------------------------------------------
!  yn_read_amf = .TRUE. ; locerrstat = pge_errstat_ok
!  IF ( yn_bro_prefit ) THEN
!     CALL he5_read_prefit_columns (                                 &
!          brofit_swath_id, nloop, nxtrack, nstart,                  &
!          lcolstr,   col_str, bro_prefit_col (1:nxtrack,0:nloop-1), &
!          ldcolstr, dcol_str, bro_prefit_dcol(1:nxtrack,0:nloop-1), &
!          yn_read_amf, locerrstat )
!     errstat = MAX( errstat, locerrstat )
!
!     ! -----------------------------------------------------------------------
!     ! Multiply BrO columns with normalization factor to return to true values
!     ! -----------------------------------------------------------------------
!     WHERE ( bro_prefit_col (1:nxtrack,0:nloop-1) > r8_missval )
!        bro_prefit_col (1:nxtrack,0:nloop-1) = &
!             bro_prefit_col (1:nxtrack,0:nloop-1) * refspecs_original(bro_idx)%NormFactor
!        bro_prefit_dcol(1:nxtrack,0:nloop-1) = &
!             bro_prefit_dcol(1:nxtrack,0:nloop-1) * refspecs_original(bro_idx)%NormFactor
!     END WHERE
!  END IF
!
!  RETURN
!END SUBROUTINE read_prefit_columns
