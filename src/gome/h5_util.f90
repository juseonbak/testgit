! From Trevor Beck, NOAA
! Modifications: cnowlan@cfa.harvard.edu

MODULE h5_util

  USE HDF5, ONLY: HID_T, HSIZE_T, HSSIZE_T

  IMPLICIT NONE

  INTEGER, PARAMETER :: MAXRANK=4
  TYPE :: DSh5_T
     CHARACTER(LEN=80) :: name
     INTEGER(HID_T)    :: dataset_id 
     INTEGER(HID_T)    :: datatype  
     INTEGER(HID_T)    :: dataspace 
     INTEGER(HID_T)    :: memspace  
     INTEGER           :: rank
     INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
     !     INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: maxdims
  END TYPE DSh5_T
  CHARACTER(LEN=11), PARAMETER :: fname = "dummy" ! File name


CONTAINS

!**************************************************
  SUBROUTINE createDS( fid, datatype, npixels, ds )
!**************************************************

    USE HDF5
    IMPLICIT NONE
    INTEGER(HID_T),   INTENT(IN)    :: fid 
    INTEGER, INTENT(IN)             :: npixels
    INTEGER(HID_T),   INTENT(IN)    :: datatype  
    TYPE (DSh5_T),    INTENT(INOUT) :: ds 
    INTEGER                         :: error

    CALL h5open_f(error)

    IF( error < 0 ) THEN
       WRITE(*,* ) "h5open_f exit fail error  "
       call exit(1)
    ENDIF

    ds%datatype = datatype
    ds%dims(1)  = npixels

    CALL h5screate_simple_f(ds%rank, ds%dims, ds%dataspace, error )
    CALL check_msg(error, "createds called h5screate_simple_f on "//trim(ds%name))

    CALL h5dcreate_f(fid, ds%name, ds%datatype, ds%dataspace, ds%dataset_id, error)
    CALL check_msg(error , "createDS, filenameh5_util "//trim(ds%name))

  END SUBROUTINE createDS
  

!***********************************************************
  SUBROUTINE createDS_ndims( fid, datatype, rank, dims, ds )
!***********************************************************
   
    ! Creates dataspace of several dimensions
    ! CRN, 28-Nov-2008
    
    USE HDF5

    IMPLICIT NONE

    INTEGER (HID_T),     INTENT(IN)      :: fid 
    INTEGER, INTENT(IN)                  :: datatype, rank
    INTEGER, INTENT(IN), DIMENSION(rank) :: dims
    TYPE (DSh5_T),       INTENT(INOUT)   :: ds 

    INTEGER          :: error

    CALL h5open_f(error)

    IF( error < 0 ) THEN
       WRITE(*,* ) "h5open_f exit fail error  "
       call exit(1)
    ENDIF

    ds%datatype = datatype
    ds%rank     = rank
    ds%dims     = dims


    CALL h5screate_simple_f(ds%rank, ds%dims, ds%dataspace, error )
    CALL check_msg(error, "createds called h5screate_simple_f on "//trim(ds%name))

    CALL h5dcreate_f(fid, ds%name, ds%datatype, ds%dataspace, ds%dataset_id, error)
    CALL check_msg(error , "createDS, filenameh5_util "//trim(ds%name))

  END SUBROUTINE createDS_ndims

!*******************************
  SUBROUTINE close_hdf5(ds, fid)
!*******************************

    USE HDF5

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: fid
    TYPE (DSh5_T), INTENT(INOUT) :: ds
    INTEGER :: error

    CALL h5sclose_f(ds%dataspace, error)
    CALL h5sclose_f(ds%memspace, error)
    CALL h5dclose_f(ds%dataset_id, error)
    CALL h5fclose_f(fid, error)

  END SUBROUTINE close_hdf5




!**************************
  SUBROUTINE open_hdf5(fid)
!**************************

    USE HDF5
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: fid
    INTEGER :: error


    CALL h5open_f(error) 
    CALL check_msg(error, "interface open call h5open_f in subr open_hdf5")

    CALL h5fcreate_f(fname,H5F_ACC_RDWR_F, fid, error)
    CALL check_msg(error, "h5freate_f  in subr open_hdf5")

  END SUBROUTINE open_hdf5





!*************************************
  SUBROUTINE write_scalar_ik4(ds, val)
!*************************************

    USE HDF5
    IMPLICIT NONE
    TYPE (DSh5_T), INTENT(INOUT) :: ds
    INTEGER, INTENT(IN) :: val
    INTEGER :: error
    INTEGER(HID_T)  :: dataspace 
    INTEGER(HID_T) :: memspace 

    dataspace = H5S_ALL_F
    memspace = H5S_ALL_F
    CALL h5dwrite_f( ds%dataset_id, ds%datatype, val, ds%dims, error,  memspace, dataspace )
    CALL check_msg(error, "h5dwrite_f on  failed  subr write_scalar_ik4 "//TRIM(ds%name))
    CALL h5sclose_f(ds%dataspace, error)

  END SUBROUTINE write_scalar_ik4





!************************************
  SUBROUTINE write_scalar_r4(ds, val)
!************************************

    USE HDF5
    IMPLICIT NONE
    TYPE (DSh5_T), INTENT(INOUT) :: ds
    REAL(KIND=4),  INTENT(IN) :: val
    INTEGER :: error
    INTEGER(HID_T) :: dataspace 
    INTEGER(HID_T) :: memspace 

    dataspace = H5S_ALL_F
    memspace = H5S_ALL_F
    CALL h5dwrite_f( ds%dataset_id, ds%datatype, val, ds%dims, error,  memspace, dataspace )
    CALL check_msg(error, "h5dwrite_f on  failed  subr write_scalar_ik4 "//trim(ds%name))
    CALL h5sclose_f(ds%dataspace, error)



  END SUBROUTINE write_scalar_r4



!*******************************************
  SUBROUTINE writeone_ik4( ds, ipix, val )
!******************************************* 

    USE HDF5

    IMPLICIT NONE

    TYPE (DSh5_T),     INTENT(INOUT) :: ds
    !INTEGER,           INTENT(IN)    :: ipix, val
    INTEGER,           INTENT(IN)    :: ipix, val(1)

    INTEGER(HSIZE_T),  DIMENSION(1)  :: dimsm = (/1/)
    INTEGER(HSSIZE_T), DIMENSION(1)  :: ipixel
    INTEGER                          :: error
    INTEGER, PARAMETER               :: memrank = 1, npts = 1  

    ipixel = ipix

    CALL h5dget_space_f(ds%dataset_id, ds%dataspace, error)
    CALL check_msg(error, "h5dget_space_f failed, writeone_ik4 on "//ds%name)

    CALL h5screate_simple_f(memrank, dimsm, ds%memspace, error)
    CALL check_msg(error, "writeone_ik4 called h5screate_simple_f on "//trim(ds%name))

    CALL h5sselect_elements_f(ds%dataspace, H5S_SELECT_SET_F, ds%rank, npts, ipixel, error)
    CALL check_msg(error, "writeone_ik4 called h5sselect_elements_f on "//trim(ds%name))

    CALL H5dwrite_f(ds%dataset_id,ds%datatype, val, ds%dims, error, ds%memspace, ds%dataspace)
    CALL check_msg(error,"H5dwrite_f writeone_ik4 "//trim(ds%name ))

    CALL h5sclose_f(ds%memspace, error)
    CALL h5sclose_f(ds%dataspace, error)

  END SUBROUTINE writeone_ik4


!*****************************************************************
  subroutine write_string( fid, ds, group_id, string, file_bname ) 
!*****************************************************************

    USE HDF5

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)    :: string
    TYPE(DSh5_T),     INTENT(INOUT) :: ds
    INTEGER(KIND=4),  INTENT(IN)    :: group_id
    LOGICAL, INTENT(IN), OPTIONAL   :: file_bname
    INTEGER                         :: nchar
    INTEGER(HID_T),   INTENT(IN)    :: fid 
    INTEGER                         :: status
    INTEGER(HID_T)                  :: dspace_id     ! Dataspace identifier  
    INTEGER(HID_T)                  :: type_id       ! Datatype identifier
    INTEGER(HSIZE_T), DIMENSION(1)  :: dimstr
    CHARACTER(LEN=1), PARAMETER     :: dir_separator="/"
    INTEGER                         :: sep_position

    dimstr(1)=1

    sep_position=0
    nchar=len(trim(adjustl(string)))
    !
    ! If the input is a filename and the file_basename is desired:
    !
    if (( present( file_bname )) .AND. ( file_bname )) then 
       sep_position = scan(string, dir_separator, BACK=.TRUE.)
       nchar=nchar-sep_position
    endif
    CALL h5tcopy_f (H5T_NATIVE_CHARACTER, type_id, status)
    CALL h5tset_size_f (type_id, nchar, status)
    CALL h5screate_f (H5S_SCALAR_F, dspace_id, status) 
    call check(status)
    CALL h5dcreate_f(fid, ds%name, type_id, dspace_id,ds%dataset_id, status)
    CALL H5dwrite_f(ds%dataset_id,type_id, trim(adjustl(string(sep_position+1:sep_position+nchar))),dimstr, status)
    call check(status)

  end subroutine write_string




!****************************************
  SUBROUTINE writeone_r4( ds, ipix, val )
!**************************************** 

    USE HDF5

    IMPLICIT NONE

    TYPE (DSh5_T),     INTENT(INOUT) :: ds
    INTEGER,           INTENT(IN)    :: ipix
    !REAL(KIND=4),      INTENT(IN)    :: val
    REAL(KIND=4),      INTENT(IN)    :: val(1)

    INTEGER(HSIZE_T),   DIMENSION(1) :: dimsm = (/1/)
    INTEGER(HSSIZE_T),  DIMENSION(1) :: ipixel
    INTEGER                          :: nwrite
    INTEGER                          :: error
    INTEGER, PARAMETER               :: memrank = 1, npts = 1
    INTEGER                          :: i

    ipixel = ipix

    CALL h5dget_space_f(ds%dataset_id, ds%dataspace, error)
    CALL check_msg(error, "h5dget_space_f failed, writeone_r4 on "//ds%name)

    CALL h5screate_simple_f(memrank, dimsm, ds%memspace, error)
    CALL check_msg(error, "writeoneR4 called h5screate_simple_f on "//trim(ds%name))

    CALL h5sselect_elements_f(ds%dataspace, H5S_SELECT_SET_F, ds%rank, npts, ipixel, error)
    CALL check_msg(error, "writeone_r4 called h5sselect_elements_f on "//trim(ds%name))

    CALL H5dwrite_f(ds%dataset_id, ds%datatype, val, ds%dims, error, ds%memspace, ds%dataspace)
    CALL check_msg(error,"H5dwrite_f writeone_r4 "//trim(ds%name))

    CALL h5sclose_f(ds%memspace, error)
    CALL h5sclose_f(ds%dataspace, error)

  END SUBROUTINE writeone_r4


!******************************************
  subroutine writeone_r8( ds, ipix, val )
!****************************************** 

    USE HDF5

    IMPLICIT NONE

    TYPE (DSh5_T),     INTENT(INOUT) :: ds
    INTEGER,           INTENT(IN)    :: ipix
    REAL(KIND=8),      INTENT(IN)    :: val

    INTEGER(HSIZE_T),  DIMENSION(1)  :: dimsm = (/1/)
    INTEGER(HSSIZE_T), DIMENSION(1)  :: ipixel
    INTEGER                          :: nwrite
    INTEGER                          :: error
    INTEGER, PARAMETER               :: memrank = 1, npts = 1 
    INTEGER                          :: i

    ipixel = ipix

    CALL h5dget_space_f(ds%dataset_id, ds%dataspace, error)
    CALL check_msg(error, "h5dget_space_f failed, writeone_r8 on "//ds%name)

    CALL h5screate_simple_f(memrank, dimsm, ds%memspace, error)
    CALL check_msg(error, "writeone_r8 called h5screate_simple_f on "//trim(ds%name))

    CALL h5sselect_elements_f(ds%dataspace, H5S_SELECT_SET_F, ds%rank, npts, ipixel, error)
    CALL check_msg(error, "writeone_r8 called h5sselect_elements_f on "//trim(ds%name))

    CALL H5dwrite_f(ds%dataset_id,ds%datatype, val, ds%dims, error, ds%memspace, ds%dataspace)
    CALL check_msg(error,"H5dwrite_f writeone_r8 "//trim(ds%name ))

    CALL h5sclose_f(ds%memspace, error)
    CALL h5sclose_f(ds%dataspace, error)

  END SUBROUTINE writeone_r8



!******************************************
  SUBROUTINE writeN_ik4( ds, ipix, N, val )
!****************************************** 

    USE HDF5

    IMPLICIT NONE

    TYPE (DSh5_T), INTENT(INOUT)    :: ds
    INTEGER      , INTENT(IN)       :: ipix, N, val(N)

    INTEGER(HSIZE_T),  DIMENSION(1) :: dimsm
    INTEGER                         :: nwrite
    INTEGER                         :: error
    INTEGER, PARAMETER              :: memrank = 1  ! Rank of the dataset in memory
    INTEGER(HSSIZE_T), DIMENSION(N) :: coord ! Elements coordinates
    INTEGER                         :: i

    dimsm = N

    DO i=1, N 
       coord(i) = ipix + i -1 
    ENDDO
    nwrite = SIZE(coord)

    CALL h5dget_space_f(ds%dataset_id, ds%dataspace, error)
    CALL check_msg(error, "h5dget_space_f failed, writeN_ik4 on "//ds%name)

    CALL h5screate_simple_f(memrank, dimsm, ds%memspace, error)
    CALL check_msg(error, "writeN_ik4 called h5screate_simple_f on "//trim(ds%name))

    CALL h5sselect_elements_f(ds%dataspace, H5S_SELECT_SET_F, ds%rank, nwrite, coord, error)
    CALL check_msg(error, "writeN_ik4 called h5sselect_elements_f on "//trim(ds%name))

    CALL H5dwrite_f(ds%dataset_id,ds%datatype, val, ds%dims, error, ds%memspace, ds%dataspace)
    CALL check_msg(error,"H5dwrite_f writeN_ik4 "//trim(ds%name ))

  END SUBROUTINE writeN_ik4



!*******************************************
  SUBROUTINE writeN_r4( ds, ipix, N, val )
!******************************************* 

    USE HDF5

    IMPLICIT NONE

    TYPE (DSh5_T), INTENT(INOUT) :: ds
    INTEGER,       INTENT(IN)    :: ipix
    INTEGER,       INTENT(IN)    :: N
    REAL(KIND=4),  INTENT(IN)    :: val(N)

    INTEGER, PARAMETER :: RANK = 2
    INTEGER(HSIZE_T), DIMENSION(1)  :: dimsm
    INTEGER                         :: error
    INTEGER, PARAMETER              :: memrank = 1
    INTEGER(HSIZE_T), DIMENSION(RANK, N) :: coord ! Elements coordinates
    INTEGER                         :: i

    dimsm = N

    DO i = 1, N
       coord(1,i) = i
       coord(2,i) = ipix
    ENDDO

    CALL h5dget_space_f(ds%dataset_id, ds%dataspace, error)
    CALL check_msg(error, "h5dget_space_f failed, writeN_r4 on "//ds%name)

    CALL h5screate_simple_f(memrank, dimsm, ds%memspace, error)
    CALL check_msg(error, "writeN_r4 called h5screate_simple_f on "//trim(ds%name))

    CALL h5sselect_elements_f(ds%dataspace, H5S_SELECT_SET_F, ds%rank, N, coord, error)
    CALL check_msg(error, "writeN_rk4 called h5sselect_elements_f on "//trim(ds%name))

    CALL h5dwrite_f(ds%dataset_id, ds%datatype, val, ds%dims, error, ds%memspace, ds%dataspace)
    CALL check_msg(error,"H5dwrite_f writeN_r4 "//trim(ds%name ))

  END SUBROUTINE writeN_r4



!*******************************************
  SUBROUTINE writeN_r8( ds, ipix, N, val )
!******************************************* 

    USE HDF5

    IMPLICIT NONE

    TYPE (DSh5_T), INTENT(INOUT) :: ds
    INTEGER,       INTENT(IN)    :: ipix, N
    REAL(KIND=8),  INTENT(IN)    :: val(N)

    INTEGER, PARAMETER :: RANK = 2
    INTEGER(HSIZE_T), DIMENSION(1)  :: dimsm
    INTEGER                         :: error
    INTEGER, PARAMETER              :: memrank = 1
    INTEGER(HSIZE_T), DIMENSION(RANK, N) :: coord ! Elements coordinates
    INTEGER                         :: i

    dimsm = N

    DO i = 1, N
       coord(1,i) = i
       coord(2,i) = ipix
    ENDDO

    CALL h5dget_space_f(ds%dataset_id, ds%dataspace, error)
    CALL check_msg(error, "h5dget_space_f failed, writeN_r4 on "//ds%name)

    CALL h5screate_simple_f(memrank, dimsm, ds%memspace, error)
    CALL check_msg(error, "writeN_r4 called h5screate_simple_f on "//trim(ds%name))

    CALL h5sselect_elements_f(ds%dataspace, H5S_SELECT_SET_F, ds%rank, N, coord, error)
    CALL check_msg(error, "writeN_rk4 called h5sselect_elements_f on "//trim(ds%name))

    CALL h5dwrite_f(ds%dataset_id, ds%datatype, val, ds%dims, error, ds%memspace, ds%dataspace)
    CALL check_msg(error,"H5dwrite_f writeN_r4 "//trim(ds%name ))

  END SUBROUTINE writeN_r8





  subroutine write_field_r4( ds, array )
    !    use
    use hdf5
    implicit none
    TYPE(DSh5_T), intent(in) :: ds
    real(kind=4), intent(in), dimension(:) :: array
    INTEGER, parameter  :: memrank = 1  ! Rank of the dataset in memory
    INTEGER(HSIZE_T), DIMENSION(1) :: dimsm 
    integer :: status
    integer :: npixels
    INTEGER (HID_T)  :: dataspace 
    INTEGER(HID_T)   :: memspace  
    INTEGER(HID_T)   :: ds_id
    INTEGER(HID_T)   :: datatype
    INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
    dims      = ds%dims
    dataspace = ds%dataspace
    memspace  = ds%memspace
    datatype  = ds%datatype
    ds_id     = ds%dataset_id
    npixels   = size(array)
    dimsm     = (/npixels/)
    !    CALL h5dget_space_f(ds_id, dataspace, status)
    !    write(*,*) "hdf5 r4 one status2 after create_simple", status, trim(ds%name)
    !    call check(status)
    CALL h5screate_simple_f(memrank, dimsm , memspace, status)
    !    call check(status)
    CALL H5dwrite_f(ds_id,datatype, array,dims, status,memspace,dataspace)
    !    call check(status)
    CALL h5sclose_f(memspace,status )
    !    call check(status)
  end subroutine write_field_r4

  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !!
  !
  !
  !
  !
  !
  !
  subroutine read_field_r4(fieldName ,group_id, array )
    use hdf5
    implicit none
    CHARACTER(LEN=*)::  fieldName
    !    CHARACTER(LEN=255)::  fieldName
    real(kind=4), intent(inout), dimension(:) :: array
    integer,intent(in) :: group_id
    INTEGER(HID_T)   :: ds_id
    INTEGER(HSIZE_T), DIMENSION(1) :: dimsm 
    integer :: status
    integer :: npixels
    npixels   = size(array)
    dimsm     = (/npixels/)
    ds_id=group_id


    CALL h5dopen_f(group_id, trim(fieldName), ds_id, status ) 
    call check_msg(status, "h5dopen_f filename h5_util.f90 in read_field_r4 name: "//trim(fieldName))
    call h5dread_f(ds_id ,H5T_NATIVE_REAL ,array, dimsm, status)
    call check_msg(status, "h5dread_f filename h5_util.f90 in read_field_r4 name: "//trim(fieldName))
    CALL h5dclose_f(ds_id, status)
    call check_msg(status, "h5dclose_f filename h5_util.f90 in read_field_r4 name: "//trim(fieldName))



  end subroutine read_field_r4


  subroutine read_field_i4(fieldName ,group_id, array )
    use hdf5
    implicit none
    CHARACTER(LEN=*)::  fieldName
    !    CHARACTER(LEN=255)::  fieldName
    integer(kind=4), intent(inout), dimension(:) :: array
    integer(kind=4),intent(in) :: group_id
    INTEGER(HID_T)   :: ds_id
    INTEGER(HSIZE_T), DIMENSION(1) :: dimsm 
    integer :: status
    integer :: npixels
    npixels   = size(array)
    dimsm     = (/npixels/)
    ds_id=group_id


    CALL h5dopen_f(group_id, trim(fieldName), ds_id, status ) 
    call check_msg(status, "h5dopen_f filename  in read_field_r4 name: "//trim(fieldName))
    call h5dread_f(ds_id ,H5T_NATIVE_INTEGER ,array, dimsm, status)
    call check_msg(status, "h5dread_f filename  in read_field_r4 name: "//trim(fieldName))
    CALL h5dclose_f(ds_id, status)
    call check_msg(status, "h5dclose_f filename  in read_field_r4 name: "//trim(fieldName))



  end subroutine read_field_i4




  subroutine check( error ) 
    implicit none
    integer, intent(in) :: error 
    if ( error < 0 ) then 
       write(*,*) "Hdf5 Error, Exit Fail"
       call exit(1)
    endif
  end subroutine check






  subroutine check_msg( error, msg ) 
    implicit none
    integer, intent(in) :: error 
    character(LEN=*), intent(in) :: msg
    if ( error < 0 ) then 
       write(*,*) "Hdf5 Error, Exit Fail "//msg
       call exit(1)
    endif
  end subroutine check_msg





end module h5_util
