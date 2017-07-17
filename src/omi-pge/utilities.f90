FUNCTION int2string ( int_num, ni ) RESULT ( int_str )

  ! ===============================================================
  ! Converts INTEGER number INT_NUM to STRING INT_STR of length NI
  ! or the number of digits in INT_NUM, whatever is larger.
  ! ===============================================================

  USE OMSAO_precision_module, ONLY: i4
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),  INTENT (IN)  :: int_num, ni

  ! ---------------
  ! Result variable
  ! ---------------
  CHARACTER (LEN=*) :: int_str

  ! ------------------------------
  ! Local variables and parameters
  ! ------------------------------
  ! * Arrays containing indices for number strings in ASCII table
  INTEGER (KIND=i4),                   PARAMETER :: n = 10
  INTEGER (KIND=i4), DIMENSION(0:n-1)            :: aidx
  CHARACTER (LEN=1), DIMENSION(0:n-1), PARAMETER :: astr = (/ "0","1","2","3","4","5","6","7","8","9" /)
  ! * Temporary and loop variables
  INTEGER (KIND=i4)                              :: i, k, nd, tmpint, ld

  ! ----------------------------------------------------
  ! Compute the index entries of ASTR in the ASCII table
  ! ----------------------------------------------------
  aidx = IACHAR(astr)

  ! ----------------------------------------------------------
  ! Find the number of digits in INT_NUM. This is equal to the 
  ! truncated integer of LOG10(INT_NUM) plus 1.
  ! ----------------------------------------------------------
  SELECT CASE ( int_num )
  CASE ( 0:9 ) 
     nd = 1
  CASE ( 10: )
     nd = INT ( LOG10(REAL(int_num)) ) + 1
  CASE DEFAULT
     int_str = '?' ; RETURN
  END SELECT

  ! -------------------------------------------------------------
  ! We may want to create a string that is longer than the number
  ! of digits in INT_NUM. This will create leading "0"s.
  ! -------------------------------------------------------------
  nd = MAX ( nd, ni )

  ! ----------------------------------------------
  ! Convert the integer to string "digit by digit"
  ! ----------------------------------------------
  int_str = "" ; tmpint = int_num
  DO i = 1, nd
     ld = MOD ( tmpint, 10 )        ! Current last digit
     tmpint = ( tmpint - ld ) / 10  ! Remove current last digit from INT_STR
     k = nd - i + 1                 ! Position of current digit in INT_STR
     int_str(k:k) = ACHAR(aidx(ld)) ! Convert INTEGER digit to CHAR
  END DO

  RETURN
END FUNCTION int2string


SUBROUTINE utc_julian_date_and_time ( year, month, day, julday, hour, minute, second )

  ! -------------------------------------------------------------
  ! Converts current date and time to UTC time and "Julian" date,
  ! i.e., day of the year.
  ! -------------------------------------------------------------

  USE OMSAO_precision_module, ONLY: i4, r4
  IMPLICIT NONE

  ! ---------------
  ! RESULT
  ! ---------------
  INTEGER (KIND=i4), INTENT (OUT) :: year, month, day, julday, hour, minute, second

  !INTEGER (KIND=i4), DIMENSION(8), INTENT (OUT) :: utc_juldate
  
  ! ---------------
  ! Local variables
  ! ---------------
  ! * Arguments for DATE_AND_TIME, and position indices for Year, Month, Day, Zone, Hour, Minute, Seconds
  CHARACTER (LEN= 8)                 :: date
  CHARACTER (LEN=10)                 :: time
  CHARACTER (LEN= 5)                 :: zone
  INTEGER   (KIND=i4), DIMENSION (8) :: date_vector
  INTEGER   (KIND=i4), PARAMETER     :: y_idx=1, m_idx=2, d_idx=3, z_idx=4, hh_idx=5, mm_idx=6, ss_idx=7
  ! * Some MAX values for Minutes and Hours (not that we would expect those to change :-)
  INTEGER   (KIND=i4), PARAMETER     :: max_mi = 60, max_hr = 24, max_dy = 31, max_mo = 12
  ! * Other local variables
  INTEGER   (KIND=i4)                :: del_mm, del_hh, max_julday, yyy, mmm, ddd

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER (KIND=i4), EXTERNAL :: day_of_year

  ! ----------------------------------------------------------
  ! First copy the input variable to the output variable. That
  ! saves work on premature return.
  ! ----------------------------------------------------------
  CALL DATE_AND_TIME ( date, time, zone, date_vector )

  ! -----------------------------------------------
  ! Initialize all output variables (except JULDAY)
  ! -----------------------------------------------
  year   = date_vector(y_idx )
  month  = date_vector(m_idx )
  day    = date_vector(d_idx )
  hour   = date_vector(hh_idx)
  minute = date_vector(mm_idx)
  second = date_vector(ss_idx)

  ! ----------------------------------------------------------------------
  ! Next, find the day of the year. At this point, we also compute the
  ! maximum number of days per year. This will come in handy further down.
  ! ----------------------------------------------------------------------
  julday     = day_of_year ( year, month,  day    )
  max_julday = day_of_year ( year, max_mo, max_dy )

  ! -----------------------------------------------
  ! Check if there is any difference to UTC at all.
  ! -----------------------------------------------
  IF ( ABS(date_vector(z_idx)) == 0 ) RETURN

  ! ------------------------------------------------------------
  ! Apply difference in UTC time to local time vector. This is
  ! given in minutes, thus apply to minutes index. 
  !
  ! REMEMBER: This is the DIFFERENCE to UTC, so we have to
  !           SUBTRACT whatever difference from the current time.
  ! -------------------------------------------------------------
  del_mm = MOD ( date_vector(z_idx), max_mi )
  minute = minute - del_mm

  ! ----------------------------------------
  ! Compute hour difference, return if none.
  ! ----------------------------------------
  del_hh =   ( date_vector(z_idx) - del_mm ) / max_mi
  IF ( del_hh == 0 ) RETURN

  ! --------------------------------------------------------
  ! Apply hour difference and check for/apply day difference
  ! --------------------------------------------------------
  hour = hour - del_hh

  ! -------------------------------------------------------------
  ! If we reach here, we may have had a day difference. We could
  ! make our lives miserable by checking for all the months of 
  ! the year. But what we ultimately want is the day of the year
  ! (a proxy for Julian day), and that makes life a lot easier.
  ! -------------------------------------------------------------
  SELECT CASE ( hour )
  CASE ( 0:23 )
     RETURN
  CASE ( :-1 )
     hour = max_hr + hour
     day  = day - 1
  CASE ( 24: )
     hour = hour - max_hr
     day  = day  + 1
  END SELECT
  
  ! --------------------------------------------------------------
  ! Finally, if we ever reach here, we've had a year difference.
  ! Apply correction, compute the new day of the year, and return.
  ! --------------------------------------------------------------
  yyy = year ; mmm = month ; ddd = day                     ! save old values
  CALL year_month_day ( yyy, mmm, ddd, year, month, day )  ! overwrite with new ones
  julday = day_of_year ( year, month, day )                ! recompute day-of-year

  RETURN
END SUBROUTINE utc_julian_date_and_time

FUNCTION day_of_year ( year, month, day ) RESULT ( jday )

  USE OMSAO_precision_module, ONLY: i4
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: year, month, day

  ! ---------------
  ! Result variable
  ! ---------------
  INTEGER (KIND=i4) :: jday

  ! ------------------------------
  ! Local variables and parameters
  ! ------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_month = 12
  INTEGER (KIND=i4), DIMENSION (n_month), PARAMETER :: &
       days_per_month = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  ! -------------------------------------------------------
  ! First find day of the year in a regular (non leap) year
  ! -------------------------------------------------------
  SELECT CASE ( month )
  CASE ( 1 ) 
     jday = day
  CASE (2:12)
     jday = SUM ( days_per_month(1:month-1) ) + day
  CASE DEFAULT
     jday = -9999
  END SELECT

  ! -------------------------
  ! Now apply leap year rules
  ! ------------------------------------
  ! * Divisible by 4:   leap year
  ! * Divisible by 100: not a leap year
  ! * Divisible by 400: leap year
  ! ------------------------------------
  IF ( MOD(year,4) == 0 .AND. ( MOD(year,100) /= 0 .OR. MOD(year,400) == 0 ) ) jday = jday+1

  RETURN
END FUNCTION day_of_year

SUBROUTINE year_month_day ( year, month, day, newyear, newmonth, newday )

  USE OMSAO_precision_module, ONLY: i4
  IMPLICIT NONE

  ! ----------------------
  ! Input/Output variables
  ! ----------------------
  INTEGER (KIND=i4), INTENT (IN)    :: year, month, day
  INTEGER (KIND=i4), INTENT (INOUT) :: newyear, newmonth, newday

  ! ------------------------------
  ! Local variables and parameters
  ! ------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_month = 12
  INTEGER (KIND=i4), DIMENSION (n_month), PARAMETER :: &
       days_per_month    = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /), &
       ly_days_per_month = (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  INTEGER (KIND=i4), DIMENSION (n_month) :: dpm

  ! -------------------------
  ! Check for leap year
  ! ------------------------------------
  ! * Divisible by 4:   leap year
  ! * Divisible by 100: not a leap year
  ! * Divisible by 400: leap year
  ! ------------------------------------
  IF ( MOD(year,4) == 0 .AND. ( MOD(year,100) /= 0 .OR. MOD(year,400) == 0 ) ) THEN
     dpm = days_per_month
  ELSE
     dpm = ly_days_per_month
  END IF

  newyear = year  ;  newmonth = month  ;  newday = day

  ! -------------------------------------------------------
  ! If DAY is within bounds there's nothing to do
  ! (MONTH hasn't been adjusted yet in the calling program)
  ! -------------------------------------------------------
  IF ( day > 0 .AND. day <= dpm(month) ) RETURN

  ! ---------------------------------------------------------------------------
  ! If we reach here, DAY must be either too large or <= 0. The crucial
  ! variable now is MONTH, since it might lead to an additional change in YEAR.
  ! ---------------------------------------------------------------------------
  SELECT CASE ( month )
  CASE ( 2:11 )   ! The easiest case: change only MONTH and DAY
     IF ( day <= 0 ) THEN
        newmonth = month - 1
        newday   = dpm(newmonth) + day  ! DAY is negative
     ELSE
        newmonth = month + 1
        newday   = day - dpm(month)
     END IF
  CASE ( 1 )      ! January requires a YEAR change if DAY <= 0
     IF ( day <= 0 ) THEN
        newmonth = 12
        newday   = dpm(newmonth) + day  ! DAY is negative
        newyear  = year - 1
     ELSE
        newmonth = month + 1
        newday   = day - dpm(month)
     END IF
  CASE ( 12 )      ! December requires a YEAR change if DAY > 31
     IF ( day <= 0 ) THEN
        newmonth = month - 1
        newday   = dpm(newmonth) + day  ! DAY is negative
     ELSE
        newmonth = 1
        newday   = day - dpm(month)
        newyear  = year + 1
     END IF
  END SELECT

  RETURN
END SUBROUTINE year_month_day


CHARACTER(LEN=*) FUNCTION upper_case ( mixstring ) RESULT ( upcase )

  ! =============================================
  ! Function to convert strings to all upper case
  ! =============================================

  IMPLICIT NONE

  ! ------------------------------------------------------------
  CHARACTER (LEN=*), INTENT (IN) :: mixstring   ! Input string
  ! ------------------------------------------------------------

  ! -------------------------------------------------------------
  !CHARACTER (LEN=*)              :: upcase      ! Result string
  ! -------------------------------------------------------------

  ! -------------------------------------------------------------
  INTEGER :: slen, i                            ! Local variables
  ! -------------------------------------------------------------

  DO i = 1, LEN(mixstring)
     SELECT CASE ( ICHAR(mixstring(i:i)) )
     CASE ( 97:122 )
        upcase(i:i) = ACHAR(ICHAR(mixstring(i:i))-32)
     CASE DEFAULT
        upcase(i:i) = mixstring(i:i)
     END SELECT
  END DO

  RETURN
END FUNCTION upper_case


SUBROUTINE check_for_endofinput ( iostring, yn_eoi )

  USE OMSAO_indices_module, ONLY : eoi_str
  IMPLICIT NONE

  ! ==============
  ! Input variable
  ! ==============
  CHARACTER (LEN=*), INTENT (IN) :: iostring

  ! ===============
  ! Output variable
  ! ===============
  LOGICAL, INTENT (OUT) :: yn_eoi

  yn_eoi = .FALSE.
  IF ( TRIM(ADJUSTL(iostring)) == eoi_str ) yn_eoi = .TRUE.

  RETURN
END SUBROUTINE check_for_endofinput

SUBROUTINE gome_check_read_status ( ios, file_read_stat )

  USE OMSAO_errstat_module, ONLY : file_read_ok, file_read_failed, file_read_eof
  IMPLICIT NONE

  INTEGER, INTENT (IN)  :: ios
  INTEGER, INTENT (OUT) :: file_read_stat

  SELECT CASE ( ios )
  CASE ( :-1 )  ;  file_read_stat = file_read_eof
  CASE (   0 )  ;  file_read_stat = file_read_ok
  CASE DEFAULT  ;  file_read_stat = file_read_failed
  END SELECT

  RETURN
END SUBROUTINE gome_check_read_status

SUBROUTINE skip_to_filemark ( funit, lm_string, lastline, file_read_stat )

  USE OMSAO_variables_module, ONLY : maxchlen
  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,           INTENT (IN) :: funit
  CHARACTER (LEN=*), INTENT (IN) :: lm_string

  ! ================
  ! Output variables
  ! ================
  INTEGER,           INTENT (OUT) :: file_read_stat
  CHARACTER (LEN=*), INTENT (OUT) :: lastline

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                  :: iline, lmlen, ios
  CHARACTER (LEN=maxchlen) :: tmpline

  ! -------------------------------------------
  ! Determine the length of the string landmark
  ! -------------------------------------------
  lmlen = LEN(TRIM(ADJUSTL(lm_string)))

  ! -------------------------------------------------------
  ! Read lines in the file until we either find the string,
  ! reach the end of the file, or reading fails otherwise.
  ! ----------------------------------------------------
  ios = 0
  getlm: DO WHILE ( ios == 0 )
     READ (UNIT=funit, FMT='(A)', IOSTAT=ios) tmpline
     tmpline = TRIM(ADJUSTL(tmpline))
     IF ( ios /= 0 .OR. tmpline(1:lmlen) == lm_string ) EXIT getlm
  END DO getlm
  
  ! ---------------------------------------------------
  ! Return the last line read for the case that we need 
  ! to extract further information from it
  ! ---------------------------------------------------
  lastline = TRIM(ADJUSTL(tmpline))

  CALL gome_check_read_status ( ios, file_read_stat )

  RETURN
END SUBROUTINE skip_to_filemark


SUBROUTINE get_substring ( string, sstart, substring, nsubstring, eostring )

  IMPLICIT NONE

  ! ==========================================================================
  ! Given string STRING, extract first space- or comma-delimited substring
  ! beginning on or after position SSTART, and return its value, SUBSTRING and
  ! its length, NSUBSTRING; update STRING and SSTART to remove this substring.
  ! 
  ! F90 version of the original GET_TOKEN by J. Lavanigno
  ! ==========================================================================

  ! NSUBSTRING represents the token's length without trailing blanks. It's
  ! set to 0 if no substring was found in STRING: this can be used to determine
  ! when to stop looking.
  !
  ! This routine makes no attempt to detect or handle the case in which
  ! the token is bigger than the token buffer.  It's assumed that the
  ! caller will declare line and token to be the same size so that this
  ! won't ever happen.

  ! ==================
  ! Modified arguments.
  ! ==================
  CHARACTER (LEN = *), INTENT (INOUT) :: string
  INTEGER,             INTENT (INOUT) :: sstart

  ! =================
  ! Output arguments.
  ! =================
  CHARACTER (LEN =LEN(string)), INTENT (OUT) :: substring
  INTEGER,                      INTENT (OUT) :: nsubstring
  LOGICAL,                      INTENT (OUT) :: eostring

  ! ================
  ! Local arguments.
  ! ================
  CHARACTER :: char
  INTEGER   :: lstart, lend, nline
  

  ! ------------------
  ! Initialize outputs
  ! ------------------
  substring  = ' '  ;  nsubstring = 0  ;  eostring = .FALSE.

  nline = LEN(TRIM(ADJUSTL(string)) )


  ! ------------------------------------------------------
  ! We are working on character variables, i.e., positions 
  ! are always larger than 0
  ! ------------------------------------------------------
  IF ( sstart <= 0 ) sstart = 1

  ! ----------------------------------------------
  ! Check first whether we have to any work at all
  ! ----------------------------------------------
  IF ( sstart >= nline ) THEN
     eostring = .TRUE.  ;  RETURN
  END IF

  ! --------------------------------------------------------
  ! Find first character in line that's not a blank or comma
  ! --------------------------------------------------------
  findchar: DO lstart = sstart, nline
     char = string(lstart:lstart)
     SELECT CASE ( char )
     CASE ( ' ' )
        IF ( lstart == nline ) THEN
           sstart = nline + 1; RETURN  ! there are no further substrings in STRING
        END IF
     CASE ( ',' )
        IF ( lstart == nline ) THEN
           sstart = nline + 1; RETURN  ! there are no further substrings in STRING
        END IF
     CASE DEFAULT
        EXIT findchar
     END SELECT
  END DO findchar


  ! --------------------------------------------------
  ! Start of the next substring is at position LSTART.
  ! Now find separator that ends the substring.
  ! --------------------------------------------------
  findsep: DO lend = lstart + 1, nline
     char = string (lend:lend)
     If ( char == ' ' .OR. char == ',' ) EXIT findsep
  END DO findsep
  IF ( (lend == nline) .AND. (char /= ' ' .AND. char /= ',') ) lend = lend + 1
  lend = lend - 1

  ! --------------------
  ! Output the substring
  ! --------------------
  nsubstring = lend - lstart + 1

  substring(1:nsubstring) = string (lstart:lend)

  ! -----------------------------------------------------------------
  ! The next substring, if any, starts at least two characters beyond
  ! the end of thelast (we have to skip over the comma or space that 
  ! marks the substring's end).
  ! -----------------------------------------------------------------
  sstart = lend + 2

  ! ---------------------------------------------
  ! Final check, whether we have to any more work
  ! ---------------------------------------------
  IF ( sstart >= nline ) eostring = .TRUE.

  RETURN
END SUBROUTINE get_substring


SUBROUTINE string2index ( table, ntable, string, stridx )

  ! =====================================================
  ! Looks up STRING in character table TABLE of dimension
  ! NTABLE, and returns position STRIDX. Defaults to
  ! STRIDX = -1 if STRING is not found in TABLE.
  ! =====================================================

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                               INTENT (IN) :: ntable
  CHARACTER (LEN=*), DIMENSION (ntable), INTENT (IN) :: table
  CHARACTER (LEN=*),                     INTENT (IN) :: string

  ! ================
  ! Output variables
  ! ================
  INTEGER, INTENT (OUT) :: stridx

  ! ===============
  ! Local variables
  ! ===============
  INTEGER :: i

  stridx = -1

  getidx: DO i = 1,  ntable
     IF ( TRIM(ADJUSTL(string)) == TRIM(ADJUSTL(table(i))) ) THEN
        stridx = i
        EXIT getidx
     END IF
  END DO getidx

  RETURN
END SUBROUTINE string2index


REAL (KIND=KIND(1.0D0)) FUNCTION signdp ( x )

  USE OMSAO_precision_module, ONLY: dp
  IMPLICIT NONE

  REAL (KIND=dp), INTENT (IN) :: x

  signdp = 0.0_dp
  IF ( x < 0.0_dp ) THEN
     signdp = -1.0_dp
  ELSE
     signdp = +1.0_dp
  END IF

  RETURN
END FUNCTION signdp

SUBROUTINE interpolation ( n_in, x_in, y_in, n_out, x_out, y_out, pge_error_status )

  USE OMSAO_precision_module
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_warning, pge_errstat_error
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,                           INTENT (IN) :: n_in, n_out
  REAL (KIND=dp), DIMENSION (n_in),  INTENT (IN) :: x_in, y_in
  REAL (KIND=dp), DIMENSION (n_out), INTENT (IN) :: x_out

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER,                           INTENT (INOUT) :: pge_error_status
  REAL (KIND=dp), DIMENSION (n_out), INTENT (OUT) :: y_out

  ! --------------
  ! Local variable
  ! --------------
  INTEGER :: errstat, imin, imax, nloc

  errstat = pge_errstat_ok

  ! -------------------------------
  ! Initialize interpolation output
  ! -------------------------------
  y_out(1:n_out) = 0.0_dp

  ! --------------------------------------------------------------------------
  ! Find indices in radiance wavelength spectrum that cover reference spectrum
  ! --------------------------------------------------------------------------

  imin = 1; imax = n_out
  IF ( x_out(1) < x_in(1) ) &
       imin = MAXVAL ( MAXLOC ( x_out(1:n_out), MASK = x_out(1:n_out) < x_in(1) ) ) + 1
  IF ( x_out(n_out) > x_in(n_in) ) &
       imax = MINVAL ( MINLOC ( x_out(1:n_out), MASK = x_out(1:n_out) > x_in(n_in) ) )

  ! ------------------------------------------------------------------------------
  ! Check whether we have the whole wavelength range. If not, set PGE_ERROR_STATUS
  ! to WARNING level, which will trigger an error message in the calling module.
  ! We don't want to report errors in this subroutine, since it is so generic that
  ! the error message would not be very helpful.
  ! ------------------------------------------------------------------------------
  IF ( imin /= 1 .OR. imax /= n_out ) THEN
     pge_error_status = pge_errstat_warning !; RETURN
  ENDIF

  ! --------------------------------------------------------------------------
  ! Now that we know the first and last index to cover with the interpolation,
  ! we can treat all cases alike. Only we need make sure that the number of
  ! interpolation points is consistent. And if we don't have *any* points,
  ! then we must return without calling the interpolation routine.
  ! --------------------------------------------------------------------------
  nloc = imax - imin + 1

  SELECT CASE ( nloc )
  CASE ( :0 )
     RETURN
  CASE DEFAULT
     CALL ezspline_interpolation (                          &
          n_in, x_in (1:n_in),    y_in (1:n_in),            &
          nloc, x_out(imin:imax), y_out(imin:imax), errstat )
     ! -----------------------------------------------------------------
     ! If we have non-zero exit status, something must have gone wrong
     ! in the interpolation. Set PGE_ERROR_STATUS to ERROR in this case.
     ! -----------------------------------------------------------------
     IF ( errstat /= pge_errstat_ok ) THEN
        pge_error_status = pge_errstat_error
     ELSE
        pge_error_status = pge_errstat_ok
     ENDIF
  END SELECT

  RETURN
END SUBROUTINE interpolation
