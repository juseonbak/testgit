MODULE he5_l2writer_class
  
    USE OMSAO_HE5_module
    USE OMSAO_errstat_module
    USE OMSAO_omidata_module, only : nlines_max
    IMPLICIT NONE

    INTEGER (KIND=2), PARAMETER, PUBLIC :: fill_int8    = -127 
    INTEGER (KIND=2), PARAMETER, PUBLIC :: fill_uint8   = 255
    INTEGER (KIND=2), PARAMETER, PUBLIC :: fill_int16   = -32767 
    INTEGER (KIND=2), PARAMETER, PUBLIC :: fill_uint16  = 65535
    REAL (KIND=4), PUBLIC               :: fill_float32
    REAL (KIND=8), PUBLIC               :: fill_float64 
    !REAL (KIND=4), PARAMETER, PUBLIC    :: fill_float32 = -1.0E+30_4
    !REAL (KIND=8), PARAMETER, PUBLIC    :: fill_float64 = -1.0E+30_8

    INTEGER (KIND=4), PARAMETER, PRIVATE :: NM = 2

    INTEGER (KIND=4), PARAMETER, PRIVATE :: zero = 0, one = 1, two = 2
    INTEGER (KIND=4), PARAMETER, PRIVATE :: three = 3, four = 4, five = 5
    
    !INTEGER (KIND=1), PUBLIC :: I1   !! dummy varaible
    !INTEGER (KIND=2), PUBLIC :: I2   !! dummy varaible
    INTEGER (KIND=4), EXTERNAL :: r4Fill, r8Fill
    INTEGER (KIND=4), EXTERNAL :: OMI_SMF_setmsg
  
    PUBLIC  :: L2_createFile
    PUBLIC  :: L2_setupSwath
    PUBLIC  :: L2_writeBlock
    PUBLIC  :: L2_copyWriteBlock
    PUBLIC  :: L2_setDFattrs
    PUBLIC  :: L2_defSWgeofields
    PUBLIC  :: L2_defSWdatafields
    PUBLIC  :: L2_newBlockW
    PUBLIC  :: L2_disposeBlockW

    CONTAINS
       FUNCTION L2_createFile( L2_file_LUN, L2_filename ) RESULT( status )
         INTEGER (KIND=4), INTENT(IN) :: L2_file_LUN
         !CHARACTER( LEN = * ), INTENT(OUT) :: L2_filename
         CHARACTER( LEN = * ), INTENT(IN) :: L2_filename
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
         INTEGER (KIND=4) :: numfiles
         INTEGER (KIND=4) :: version
         INTEGER (KIND=4) :: status, ierr
         INTEGER (KIND=4) :: SW_fileid
 
         !! Get number of L2 output files from PCF
         !! and make sure there is only one 

         !status = PGS_PC_getnumberoffiles( L2_file_LUN,  numfiles )
         !IF( status /= PGS_S_SUCCESS ) THEN
         !   WRITE( msg,'(A,I)' ) "can't get numfiles from PCF file at LUN =", &
         !                         L2_file_LUN
         !   ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "L2_createFile", zero )
         !   status = OZT_E_FAILURE
         !   RETURN
         !ELSE IF( numfiles /= one ) THEN
         !   WRITE( msg,'(A,I,A,I)' ) "numfiles =", numfiles, "not equal to 1"//&
         !                            " at LUN =", L2_file_LUN
         !   ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "L2_createFile", zero )
         !   status = OZT_E_FAILURE
         !   RETURN
         !ENDIF
         !
         !!! Get the L2 file name from the PCF.
         !version = 1
         !status = PGS_PC_getreference( L2_file_LUN, version, L2_filename )
         !IF( status /= PGS_S_SUCCESS ) THEN
         !   WRITE( msg,'(A,I)' ) "get filename from PCF file at LUN =", &
         !                         L2_file_LUN
         !   ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "L2_createFile", zero )
         !   status = OZT_E_FAILURE
         !   RETURN
         !ENDIF
         version = 1

         SW_fileid = he5_swopen( L2_filename, HE5F_ACC_TRUNC )
         
         IF( SW_fileid == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_swopen:"// TRIM(L2_filename) // " failed."
            ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, msg, "L2_createFile", zero )
            status = OMI_E_FAILURE
            RETURN 
         ELSE
            WRITE( msg,'(A)' ) "he5_swopen:"// TRIM(L2_filename) //&
                               " created succesfully."
            ierr = OMI_SMF_setmsg( OMI_S_SUCCESS, msg, "L2_createFile", &
                                   four )
         ENDIF
         status = he5_swclose( SW_fileid )
         IF (status == 0) THEN
            status = OMI_S_SUCCESS
         ELSE
            WRITE( msg,'(A)' ) "he5_swclose:"// TRIM(L2_filename) // " failed."
            ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, msg, "L2_createFile", zero )
            status = OMI_E_FAILURE
         ENDIF

         ierr = r4Fill( fill_float32 )
         ierr = r8Fill( fill_float64 )

         RETURN 
       END FUNCTION L2_createFile

       FUNCTION L2_setupSwath( L2_filename, dimList, dims, swath_name_LUN, &
                               L2_swathname, SW_fileid, SW_id) RESULT( status )
         CHARACTER( LEN = * ), INTENT(IN) :: L2_filename
         CHARACTER( LEN = * ), INTENT(IN) :: dimList
         !CHARACTER( LEN = * ), INTENT(OUT) :: L2_swathname
         CHARACTER( LEN = * ), INTENT(IN) :: L2_swathname
         INTEGER (KIND=4), INTENT(IN) :: swath_name_LUN
         INTEGER (KIND=4), DIMENSION(:), INTENT(IN) :: dims
         INTEGER (KIND=4) :: status, ierr, id, nc
         INTEGER (KIND=4) :: nDims
         CHARACTER( LEN = MAX_STR_LEN ), DIMENSION(NDIM_MAX) :: dimName
         INTEGER (KIND=4), DIMENSION(NDIM_MAX) :: dimSize
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
         INTEGER (KIND=4), INTENT(OUT) :: SW_fileid, SW_id
         
         SW_fileid = he5_swopen( L2_filename, HE5F_ACC_RDWR )
         IF( SW_fileid == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_swopen:"// TRIM(L2_filename) // " failed."
            ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, msg, "L2_setupSwath", zero )
            status = OMI_E_FAILURE
            RETURN 
         ENDIF

         !!! Get the L2 Swath name from the PCF.
         !status = PGS_PC_GetConfigData( swath_name_LUN, msg )
         !IF( status /= PGS_S_SUCCESS ) THEN
         !   WRITE( msg,'(A,I)' ) "get swathname from PCF file at LUN =", &
         !                        swath_name_LUN
         !   ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "L2_setupSwath", zero )
         !   status = OZT_E_FAILURE
         !   RETURN
         !ELSE
         !  nc = deQuote( msg )
         !  IF( nc <= 0 ) THEN
         !     WRITE( msg,'(A)' ) "Error in Input Swath Name:"// TRIM(msg) 
         !     ierr = OMI_SMF_setmsg( OZT_E_HDFEOS,msg,"L2_setupSwath",zero )
         !     status = OZT_E_FAILURE
         !     RETURN 
         !  ELSE
         !     L2_swathname = TRIM( msg )
         !  ENDIF 
         !ENDIF

         !! create the swath
         SW_id = he5_swcreate( SW_fileid, L2_swathname )
            
         IF( SW_id == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_swcreate:"// TRIM(L2_swathname) //&
                               " failed in file " // TRIM(L2_filename )
            ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, msg, "L2_setupSwath", zero )
            status = OMI_E_FAILURE
            RETURN 
         ELSE
            WRITE( msg,'(A)' ) "he5_swcreate:"// TRIM(L2_swathname) //&
                               " created succesfully in file " // &
                               TRIM(L2_filename)
            ierr = OMI_SMF_setmsg( OMI_S_SUCCESS, msg, "L2_setupSwath", &
                                   four )
         ENDIF

         nDims = EH_parsestrF( dimList, ',', dimName )
         IF( nDims > SIZE( dims ) ) THEN
            ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                   "input dimList and dims do not match", &
                                   "L2_setupSwath", zero )
            status = OMI_E_FAILURE
            RETURN 
         ENDIF  
         dimSize(1:nDims) = dims(1:nDims)

         DO id = 1, nDims
            status = he5_swdefdim( SW_id, TRIM( dimName(id) ), dimSize(id) )
            IF( status == -1 ) THEN
               WRITE( msg,'(A)' ) "he5_swdefdim "//TRIM(dimName(id))//"failed."
               ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, msg, "L2_setupSwath", &
                                      zero )
               status = OMI_E_FAILURE
               RETURN
            ENDIF
         ENDDO
         status = OMI_S_SUCCESS
         RETURN
       END FUNCTION L2_setupSwath

       FUNCTION L2_defSWgeofields( SW_id, gfList ) RESULT( status )
         INTEGER (KIND=4), INTENT(IN) :: SW_id
         TYPE (DFHE5_T), DIMENSION(:), INTENT(IN) :: gfList
         INTEGER (KIND=4) :: ig, NG, status, ierr
         INTEGER (KIND=4) :: ifld, nflds
         INTEGER (KIND=4), DIMENSION(HE5_DTSETRANKMAX) :: dims
         INTEGER (KIND=4), DIMENSION(HE5_FLDNUMBERMAX) :: rank, ntype
         CHARACTER( LEN = MAX_STR_LEN ) :: dimlist, maxdimlist
         CHARACTER( LEN = 4096 ) :: fieldlist
         CHARACTER( LEN = 1 ) :: delim = ','
         CHARACTER ( LEN = 256 ), DIMENSION(HE5_FLDNUMBERMAX) :: outstrs
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg

         nflds = HE5_SWinqgflds( SW_id, fieldlist, rank, ntype )
         nflds = EH_parsestrF( fieldlist, delim, outstrs ) 
     
         NG = SIZE( gfList )
         DO ig = 1, NG
           !! check to make sure the field does not exist before
           !! define it 
           DO ifld = 1, nflds
          
             IF( TRIM( outstrs(ifld) ) == TRIM(gfList(ig)%name) ) EXIT 
           ENDDO
             
           IF( ifld > nflds ) THEN
              status = he5_swdefgfld( SW_id,                    &
                                      gfList(ig)%name,          &
                                      gfList(ig)%dimnames, " ", &
                                      gfList(ig)%datatype,      &
                                      HE5_HDFE_NOMERGE )
              IF( status == -1 ) THEN
                 WRITE( msg,'(A)' ) "he5_swdefgfld failed for " // &
                                    TRIM(gfList(ig)%name)
                 ierr=OMI_SMF_setmsg( OMI_E_HDFEOS, msg, &
                                     "L2_defSWgeofields", zero )
                 status = OMI_E_FAILURE
                 RETURN
              ENDIF
              status = L2_setDFattrs( SW_id, gfList(ig) )
           ENDIF
         ENDDO
 
         status = OMI_S_SUCCESS 
         RETURN
       END FUNCTION L2_defSWgeofields

       FUNCTION L2_defSWdatafields( SW_id, dfList ) RESULT( status )
         INTEGER (KIND=4), INTENT(IN) :: SW_id
         TYPE (DFHE5_T), DIMENSION(:), INTENT(IN) :: dfList
         INTEGER (KIND=4) :: id, ND, status, ierr
         INTEGER (KIND=4) :: ifld, nflds
         INTEGER (KIND=4), DIMENSION(HE5_DTSETRANKMAX) :: dims
         CHARACTER( LEN = MAX_STR_LEN ) :: dimlist, maxdimlist
                  INTEGER (KIND=4), DIMENSION(HE5_FLDNUMBERMAX) :: rank, ntype
         CHARACTER( LEN = 4096 ) :: fieldlist
         CHARACTER( LEN = 1 ) :: delim = ','
         CHARACTER ( LEN = 256 ), DIMENSION(HE5_FLDNUMBERMAX) :: outstrs
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg

         nflds = HE5_SWinqdflds( SW_id, fieldlist, rank, ntype )


         nflds = EH_parsestrF( fieldlist, delim, outstrs ) 
         ND = SIZE( dfList )
         DO id = 1, ND
           !! check to make sure the field does not exist before
           !! define it 
           DO ifld = 1, nflds
             IF( TRIM( outstrs(ifld) ) == TRIM(dfList(id)%name) ) EXIT 
           ENDDO

           !PRINT *, ID, ND, dfList(id)%name
           IF( ifld > nflds ) THEN
              status = he5_swdefdfld( SW_id,                    &
                                   dfList(id)%name,          &
                                   dfList(id)%dimnames, " ", &
                                   dfList(id)%datatype,      &
                                   HE5_HDFE_NOMERGE )
              IF( status == -1 ) THEN
                 WRITE( msg,'(A)' ) "he5_swdefdfld failed for " // &
                                    TRIM(dfList(id)%name)
                 ierr=OMI_SMF_setmsg( OMI_E_HDFEOS, msg, &
                                     "L2_defSWdatafields", zero )
                 status = OMI_E_FAILURE
                 RETURN
              ENDIF
              status = L2_setDFattrs( SW_id, dfList(id) )
           ENDIF
         ENDDO
 
         status = OMI_S_SUCCESS 
         RETURN
       END FUNCTION L2_defSWdatafields
       
       FUNCTION L2_newBlockW( this, filename, swathname, &
                              SW_fileid, SW_id, fieldList, nL ) &
                              RESULT( status )
         CHARACTER( LEN = * ), INTENT(IN) :: filename, swathname
         TYPE (DFHE5_T), DIMENSION(:), INTENT(INOUT) :: fieldList 
         INTEGER (KIND=4), INTENT(IN), OPTIONAL :: nL
         INTEGER (KIND=4) :: status, ierr, id, rankID
         INTEGER (KIND=4) :: ntype
         CHARACTER( LEN = MAX_STR_LEN ) :: dimlist, maxdimlist
         TYPE (L2_generic_type), INTENT( OUT ) :: this
         INTEGER (KIND=4), INTENT(IN) :: SW_fileid, SW_id
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg

         INTEGER (KIND=4), EXTERNAL :: HE5Tget_size
         
         ! Nulliify the pointer
         Nullify(this%data)

         this%filename  = filename
         this%swathname = swathname
         this%sw_fid    = SW_fileid
         this%swathID   = SW_id

         !! get the dimensions associated with this swath
         this%nDims = HE5_SWinqdims( this%swathID, this%dimnames, &
                                     this%dimSizes )
         IF( this%nDims == -1 ) THEN
            WRITE( msg, '(A)' ) "HE5_SWinqdims:" // TRIM( this%swathname ) // &
                                ", ", TRIM(this%filename)// " failed."
            ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, msg, "L2_newBlock", zero )
            status = OMI_E_FAILURE
            RETURN
         ENDIF

         this%nFields = SIZE( fieldList )
         IF( this%nFields == 0 .OR. this%nFields > NFLDS_MAX ) THEN
            WRITE( msg, '(A,I2,A)' ) "number of fields =", this%nFields, &
                                     " out of range."
            ierr = OMI_SMF_setmsg( OMI_E_INPUT, msg, "L2_newBlockW", zero )
            status = OMI_E_FAILURE
            RETURN
         ENDIF
         
         !print *, 'Define new block: '
         !print *, this%nDims, this%nFields
         
         !! Find out the dimension size, the rank, and the datatype and
         !! and byte size for each of the field (geo or data).
         this%SumElmSize = 0
         DO id = 1, this%nFields
           !print *, id, fieldList(id)%name  
           this%fieldname(id) = fieldList(id)%name  
           ierr = HE5_SWfldinfo( this%swathID, this%fieldname(id), rankID, &
                                 this%dims(id,:), ntype, dimlist, maxdimlist )

           IF( ierr == - 1 ) THEN
              WRITE( msg, '(A)') "HE5_SWfldinfo failed for "// &
                                 TRIM( this%fieldname(id)) // ", "// &
                                 TRIM( this%swathname )
              ierr = OMI_SMF_setmsg( OMI_E_INPUT, msg, "L2_newBlock", zero )
              status = OMI_E_FAILURE
              RETURN
           ENDIF
            
           IF( rankID /= fieldList(id)%rank ) THEN
              WRITE( msg,'(A,2I3,A)' ) "Inconsistent rank ", rankID, &
                                        fieldList(id)%rank, &
                                        TRIM( this%fieldname(id) )
              ierr = OMI_SMF_setmsg( OMI_E_INPUT, msg, "L2_newBlockW", zero )
              status = OMI_E_FAILURE
              RETURN
           ELSE IF( rankID > MAXRANK ) THEN
              WRITE( msg,'(A,I,A)' ) "rank= ", rankID, " exceed 3 for " // &
                                     TRIM( this%fieldname(id) )
              ierr = OMI_SMF_setmsg( OMI_E_INPUT, msg, "L2_newBlockW", zero )
              status = OMI_E_FAILURE
              RETURN
           ELSE
              this%rank(id) = rankID
           ENDIF
           fieldList(id)%dims(1:rankID) = this%dims(id,1:rankID) 

           this%elmSize(id) = HE5Tget_size( fieldList(id)%datatype )
           IF( this%elmSize(id) <= 0 ) THEN
              WRITE( msg,'(A,I,A)' ) "datatype: ", fieldList(id)%datatype, &
                                     " error for " // TRIM( this%fieldname(id) )
              ierr = OMI_SMF_setmsg( OMI_E_INPUT, msg, "L2_newBlockW", zero )
              status = OMI_E_FAILURE
              RETURN
           ELSE
              this%SumElmSize = this%SumElmSize + this%elmSize(id)
           ENDIF
           !print *, id, this%fieldname(id), this%rank(id), this%dims(id, 1:rankID), this%elmSize(id)

           !! make sure every field has the same first (HDF5) dimension size
           !! Note that first HDF5 or C dimension is the last Fortran
           !! dimension.
           IF( id == 1 ) THEN
              this%nTotLine = this%dims(1,this%rank(1))
           ELSE
              IF( this%nTotLine /= this%dims(id,this%rank(id)) ) THEN
                 WRITE( msg, '(A,2I5)') "First dimension size not equal:", &
                        this%dims(id,this%rank(id)), this%nTotLine
                 ierr = OMI_SMF_setmsg( OMI_E_INPUT, msg, "L2_newBlockW", zero )
                 status = OMI_E_FAILURE
                 RETURN
              ENDIF
           ENDIF
         ENDDO

         !! the first C (HDF5) dimesnion is usually referred to as the line
         !! The block size is determined by the number of Lines in the block.
         IF( PRESENT(nL) ) THEN
            IF( nL < this%nTotLine ) THEN
               this%nLine = nL      !! use input nL (if present) and if it is
            ELSE                    !! smaller than the dim size.
               this%nLine = this%nTotLine
            ENDIF
         ELSE
            IF( this%nTotLine < nlines_max ) THEN
               this%nLine = this%nTotLine
            ELSE
               this%nLine = nlines_max   !! default nLine is 100
            ENDIF
         ENDIF
         
         !!calculate the line size and block size and the accumulated block size
         this%SumLineSize = 0
         this%accuLineSize(0) = 0
         this%accuBlkSize(0) = 0
         DO id = 1, this%nFields
           rankID = this%rank(id)
           IF( rankID == 1 ) THEN
               this%pixSize(id) = this%elmSize(id)
              this%lineSize(id) = this%elmSize(id)
           ELSE IF( rankID == 2 ) THEN
               this%pixSize(id) = this%elmSize(id)
              this%lineSize(id) = this%elmSize(id) * this%dims(id,1)
           ELSE IF( rankID == 3 ) THEN
              this%pixSize(id) = this%elmSize(id) * this%dims(id,1)
              this%lineSize(id) = this%pixSize(id) * this%dims(id,2)
           ELSE 
              this%pixSize(id) = this%elmSize(id) * this%dims(id,1) * this%dims(id,2)
              this%lineSize(id) = this%pixSize(id) * this%dims(id,3)
           ENDIF
           this%SumLineSize = this%SumLineSize + this%lineSize(id)
           this%blkSize(id) = this%lineSize(id) * this%nLine
           this%accuLineSize(id) = this%accuLineSize(id-1) + this%lineSize(id)
            this%accuBlkSize(id) =  this%accuBlkSize(id-1) +  this%blkSize(id)
         ENDDO

         IF( this%accuBlkSize(this%nFields) > 0 ) THEN
            IF( ASSOCIATED( this%data ) ) DEALLOCATE( this%data )
            ALLOCATE( this%data( this%accuBlkSize(this%nFields) ), STAT = ierr )
            IF( ierr /= zero ) THEN
               status = OMI_E_FAILURE
               ierr = OMI_SMF_setmsg( OMI_E_MEM_ALLOC, "this%data", &
                                      "L2_newBlockW", zero )
               RETURN
            ENDIF
         ENDIF
         this%iLine = -1
         this%eLine = -1
         this%initialized = .TRUE.
         status = OMI_S_SUCCESS

         RETURN
       END FUNCTION L2_newBlockW
         
       FUNCTION L2_writeBlock( this, iLine, nLw )  RESULT( status )
         TYPE (L2_generic_type), INTENT( INOUT ) :: this
         INTEGER (KIND = 4), INTENT( IN ) :: iLine
         INTEGER (KIND = 4), INTENT( IN ), OPTIONAL :: nLw
         INTEGER (KIND = 4) :: ierr, status
         INTEGER (KIND = 4) :: ii, jj, id, k, rank, nL
         INTEGER (KIND=4 ), DIMENSION(MAXRANK) :: start, stride, edge
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg

         IF( .NOT. this%initialized ) THEN
            ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                  "input block not initialized", &
                                  "L2_writeBlock", zero )
            status = OMI_E_FAILURE
            RETURN
         ENDIF

         IF( iLine < 0 .OR. iLine >= this%nTotLine ) THEN
            ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                 "iLine out of range", "L2_writeBlock", zero )
            status = OMI_E_FAILURE
            RETURN
         ELSE
            this%iLine = iLine
         ENDIF

         IF( (iLine + this%nLine) > this%nTotLine ) THEN
            nL = this%nTotLine - iLine
         ELSE
            nL = this%nLine
         ENDIF

         this%eLine = this%iLine + nL - 1
         IF( PRESENT( nLw ) ) THEN
            IF( nL > nLw ) nL = nLw
         ENDIF
         
         stride(:) = 1
         DO id = 1, this%nFields
           rank = this%rank(id)
           DO k = 1, rank
             IF( k == rank ) THEN
                start(k) = iLine
                edge(k)  = nL
             ELSE
                start(k) = 0
                edge(k)  = this%dims(id,k)
             ENDIF
           ENDDO

           ii = this%accuBlkSize(id-1)+1
           jj = this%accuBlkSize(id)
           status = HE5_swwrfld( this%swathID, this%fieldname(id), &
                                 start, stride, edge, this%data(ii:jj) )
           IF( status == -1 ) THEN
              WRITE( msg,'(A)' ) "Write "//TRIM(this%fieldname(id))//&
                                  " failed."
              ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, msg,"L2_writeBlock",zero )
              status = OMI_E_FAILURE
              RETURN
           ENDIF

         ENDDO
         status = OMI_S_SUCCESS
         RETURN
       END FUNCTION L2_writeBlock

       FUNCTION L2_copyWriteBlock( swathID, this, iLine, nLw ) RESULT( status )
         TYPE (L2_generic_type), INTENT( INOUT ) :: this
         INTEGER (KIND = 4), INTENT( IN ) :: iLine
         INTEGER (KIND = 4), INTENT( IN ), OPTIONAL :: nLw
         INTEGER (KIND = 4), INTENT( IN ) :: swathID
         INTEGER (KIND = 4) :: ierr, status
         INTEGER (KIND = 4) :: ii, jj, id, k, rank, nL
         INTEGER (KIND=4 ), DIMENSION(MAXRANK) :: start, stride, edge
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg

         IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                  "input block not initialized", &
                                  "L2_writeBlock", zero )
            status = OMI_E_FAILURE
            RETURN
         ENDIF

         IF( iLine < 0 .OR. iLine >= this%nTotLine ) THEN
            ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                 "iLine out of range", "L2_writeBlock", zero )
            status = OMI_E_FAILURE
            RETURN
         ELSE
            this%iLine = iLine
         ENDIF

         IF( (iLine + this%nLine) > this%nTotLine ) THEN
            nL = this%nTotLine - iLine
         ELSE
            nL = this%nLine
         ENDIF
         this%eLine = this%iLine + nL - 1
         IF( PRESENT( nLw ) ) THEN
            IF( nL > nLw ) nL = nLw
         ENDIF

         stride(:) = 1
         DO id = 1, this%nFields
           rank = this%rank(id)
           DO k = 1, rank
             IF( k == rank ) THEN
                start(k) = iLine
                edge(k)  = nL
             ELSE
                start(k) = 0
                edge(k)  = this%dims(id,k)
             ENDIF
           ENDDO

           ii = this%accuBlkSize(id-1)+1
           jj = this%accuBlkSize(id)
           status = HE5_swwrfld( swathID, this%fieldname(id), &
                                 start, stride, edge, this%data(ii:jj) )
           IF( status == -1 ) THEN
              WRITE( msg,'(A)' ) "Write "//TRIM(this%fieldname(id))//&
                                  " failed."
              ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, msg,"L2_writeBlock",zero )
              status = OMI_E_FAILURE
              RETURN
           ENDIF

         ENDDO
         status = OMI_S_SUCCESS
         RETURN
       END FUNCTION L2_copyWriteBlock

       SUBROUTINE L2_disposeBlockW( this )
         TYPE (L2_generic_type), INTENT( INOUT ) :: this
         INTEGER (KIND=4) :: ierr

         ! Nullify the pointer
         Nullify(this%data)

         this%nDims   = 0
         this%nFields = 0
         this%iLine   = -1
         this%eLine   = -1
         this%nLine   =  0
         IF( ASSOCIATED( this%data ) ) DEALLOCATE( this%data )
     !   IF( this%initialized ) THEN
     !      IF( this%swathID /= -1 ) THEN
     !         ierr = HE5_SWdetach( this%swathID )
     !         this%swathID = -1
     !      ENDIF
     !      IF( this%sw_fid /= -1 ) THEN
     !         ierr = HE5_SWclose( this%sw_fid )
     !         this%sw_fid = -1
     !      ENDIF
     !   ENDIF
         this%initialized = .FALSE.
       END SUBROUTINE L2_disposeBlockW

       FUNCTION L2_setDFattrs( SW_id, df ) RESULT( status )
         INTEGER (KIND=4), INTENT(IN) :: SW_id
         TYPE (DFHE5_T), INTENT(IN) :: df
         INTEGER (KIND=4) :: nc
         INTEGER (KIND=4) :: status, ierr
         INTEGER (KIND=1) :: range_int8(2)
         INTEGER (KIND=2) :: range_int16(2)
         REAL (KIND=4) :: range_float(2)
         REAL (KIND=8) :: range_double(2)
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg

         nc = LEN( TRIM(df%Units) ) 
         IF( nc > 0 ) THEN
           status = he5_swwrlattr( SW_id, df%name, "Units", &
                                   HE5T_NATIVE_CHAR, nc, df%Units )      
           IF( status == -1 ) THEN
              WRITE( msg,'(A)' ) " he5_swwrlattr: Units, " // TRIM(df%name )
              ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
           ENDIF
         ENDIF
                                
         nc = LEN( TRIM(df%LongName) ) 
         IF( nc > 0 ) THEN
           status = he5_swwrlattr( SW_id, df%name, "Title", &
                                   HE5T_NATIVE_CHAR, nc, df%LongName )      
           IF( status == -1 ) THEN
              WRITE( msg,'(A)' ) " he5_swwrlattr: LongName, " // TRIM(df%name )
              ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
           ENDIF
         ENDIF

         nc = LEN( TRIM(df%UniqueFieldDefinition) ) 
         IF( nc > 0 ) THEN
           status = he5_swwrlattr( SW_id, df%name, "UniqueFieldDefinition", &
                                   HE5T_NATIVE_CHAR, nc, &
                                   df%UniqueFieldDefinition )      
           IF( status == -1 ) THEN
              WRITE( msg,'(A)' ) " he5_swwrlattr: UniqueFieldDefinition, " // &
                     TRIM(df%name )
              ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
           ENDIF
         ENDIF

         nc = 1
         status = he5_swwrlattr( SW_id, df%name, "ScaleFactor", &
                                 HE5T_NATIVE_DOUBLE, nc, df%ScaleFactor )      
         IF( status == -1 ) THEN
            WRITE( msg,'(A)' ) " he5_swwrlattr: ScaleFactor, " // TRIM(df%name )
            ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
         ENDIF
         
         status = he5_swwrlattr( SW_id, df%name, "Offset", &
                                 HE5T_NATIVE_DOUBLE, nc, df%Offset )      
         IF( status == -1 ) THEN
            WRITE( msg,'(A)' ) " he5_swwrlattr: Offset, " // TRIM(df%name )
            ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
         ENDIF

         nc = 2
         SELECT CASE( df%datatype )
           CASE ( HE5T_NATIVE_INT8 )
             range_int8(1) = NINT(df%ValidRange_l)
             range_int8(2) = NINT(df%ValidRange_h)
             status = he5_swwrlattr( SW_id, df%name, "ValidRange", &
                                     HE5T_NATIVE_INT8, nc, range_int8 )      
             IF( status == -1 ) THEN
                WRITE( msg,'(A)' ) " he5_swwrlattr: ValidRange, " // &
                                   TRIM(df%name )
                ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
             ENDIF

             status = he5_swwrlattr( SW_id, df%name, "_FillValue", &
                                     HE5T_NATIVE_INT8, 1, fill_int8 )      
             status = he5_swwrlattr( SW_id, df%name, "MissingValue", &
                                     HE5T_NATIVE_INT8, 1, fill_int8 )      
             IF( status == -1 ) THEN
                WRITE( msg,'(A)' ) " he5_swwrlattr: MissingValue, " // &
                                   TRIM(df%name )
                ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
             ENDIF

           CASE ( HE5T_NATIVE_UINT8 )
             range_int8(1) = NINT(df%ValidRange_l)
             range_int8(2) = NINT(df%ValidRange_h)
             status = he5_swwrlattr( SW_id, df%name, "ValidRange", &
                                     df%datatype, nc, range_int8 )      
             IF( status == -1 ) THEN
                WRITE( msg,'(A)' ) " he5_swwrlattr: ValidRange, " // &
                                   TRIM(df%name )
                ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
             ENDIF
             status = he5_swwrlattr( SW_id, df%name, "_FillValue", &
                                     df%datatype, 1, fill_uint8 )      
             status = he5_swwrlattr( SW_id, df%name, "MissingValue", &
                                     df%datatype, 1, fill_uint8 )      
             IF( status == -1 ) THEN
                WRITE( msg,'(A)' ) " he5_swwrlattr: MissingValue, " // &
                                   TRIM(df%name )
                ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
             ENDIF
           CASE ( HE5T_NATIVE_INT16 )
             range_int16(1) = NINT(df%ValidRange_l)
             range_int16(2) = NINT(df%ValidRange_h)
             status = he5_swwrlattr( SW_id, df%name, "ValidRange", &
                                     df%datatype, nc, range_int16 )      
             IF( status == -1 ) THEN
                WRITE( msg,'(A)' ) " he5_swwrlattr: ValidRange, " // &
                                   TRIM(df%name )
                ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
             ENDIF
             status = he5_swwrlattr( SW_id, df%name, "_FillValue", &
                                     df%datatype, 1, fill_int16 )      
             status = he5_swwrlattr( SW_id, df%name, "MissingValue", &
                                     df%datatype, 1, fill_int16 )      
             IF( status == -1 ) THEN
                WRITE( msg,'(A)' ) " he5_swwrlattr: MissingValue, " // &
                                   TRIM(df%name )
                ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
             ENDIF
           CASE ( HE5T_NATIVE_UINT16 )
             range_int16(1) = NINT(df%ValidRange_l)
             range_int16(2) = NINT(df%ValidRange_h)
             status = he5_swwrlattr( SW_id, df%name, "ValidRange", &
                                     df%datatype, nc, range_int16 )      
             IF( status == -1 ) THEN
                WRITE( msg,'(A)' ) " he5_swwrlattr: ValidRange, " // &
                                   TRIM(df%name )
                ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
             ENDIF
             status = he5_swwrlattr( SW_id, df%name, "_FillValue", &
                                     df%datatype, 1, fill_uint16 )      
             status = he5_swwrlattr( SW_id, df%name, "MissingValue", &
                                     df%datatype, 1, fill_uint16 )      
             IF( status == -1 ) THEN
                WRITE( msg,'(A)' ) " he5_swwrlattr: MissingValue, " // &
                                   TRIM(df%name )
                ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
             ENDIF
           CASE ( HE5T_NATIVE_FLOAT )
             range_float(1) = df%ValidRange_l
             range_float(2) = df%ValidRange_h
             status = he5_swwrlattr( SW_id, df%name, "ValidRange", &
                                     df%datatype, nc, range_float )      
             IF( status == -1 ) THEN
                WRITE( msg,'(A)' ) " he5_swwrlattr: ValidRange, " // &
                                   TRIM(df%name )
                ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
             ENDIF
             ierr = r4Fill( fill_float32 )
             status = he5_swwrlattr( SW_id, df%name, "_FillValue", &
                                     df%datatype, 1, fill_float32 )      
             status = he5_swwrlattr( SW_id, df%name, "MissingValue", &
                                     df%datatype, 1, fill_float32 )      
             IF( status == -1 ) THEN
                WRITE( msg,'(A)' ) " he5_swwrlattr: ValidRange, " // &
                                   TRIM(df%name )
                ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
             ENDIF
           CASE ( HE5T_NATIVE_DOUBLE )
             range_double(1) = df%ValidRange_l
             range_double(2) = df%ValidRange_h
             status = he5_swwrlattr( SW_id, df%name, "ValidRange", &
                                     df%datatype, nc, range_double )      
             IF( status == -1 ) THEN
                WRITE( msg,'(A)' ) " he5_swwrlattr: ValidRange, " // &
                                   TRIM(df%name )
                ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
             ENDIF
             ierr = r8Fill( fill_float64 )
             status = he5_swwrlattr( SW_id, df%name, "_FillValue", &
                                     df%datatype, 1, fill_float64 )      
             status = he5_swwrlattr( SW_id, df%name, "MissingValue", &
                                     df%datatype, 1, fill_float64 )      
             IF( status == -1 ) THEN
                WRITE( msg,'(A)' ) " he5_swwrlattr: ValidRange, " // &
                                   TRIM(df%name )
                ierr=OMI_SMF_setmsg( OMI_E_HDFEOS,msg,"L2_setDFattrs",zero )
             ENDIF
           CASE DEFAULT 
            ierr = OMI_SMF_setmsg( OMI_W_GENERAL, "Valid range not set for "// &
                                   TRIM(df%name) // " because of datatype.",&
                                  "L2_setDFattrs", zero )
         END SELECT 

         status = OMI_S_SUCCESS
         RETURN
       END FUNCTION L2_setDFattrs
       
     END MODULE he5_l2writer_class
