! ================================================================================
!
! Collection of subroutines that relate to PGE identification and processing aids.
!
! ================================================================================

SUBROUTINE get_pge_ident ( tmpchar, pge_name, pge_idx, pge_errstat )

  ! ====================================================
  ! Find name and index of current PGE from input string
  ! ====================================================

  USE OMSAO_indices_module, ONLY: n_sao_pge, sao_pge_names, sao_pge_min_idx, sao_pge_max_idx
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! --------------
  ! Input variable
  ! --------------
  CHARACTER (LEN=*), INTENT (IN) :: tmpchar

  ! ----------------
  ! Output variables
  ! ----------------
  CHARACTER (LEN=*), INTENT (OUT) :: pge_name
  INTEGER,           INTENT (OUT) :: pge_idx, pge_errstat

  ! --------------
  ! Local variable
  ! --------------
  INTEGER :: i

  ! ---------------------------
  ! Initialize output variables
  ! ---------------------------
  pge_idx = -1  ;  pge_name(1:1) = '?' ; pge_errstat = pge_errstat_ok

  ! -------------------------------------------------
  ! Find name and index by looping over all SAO PGEs.
  ! Not very elegant but simple and effective.
  ! -------------------------------------------------
  getpge: DO i = sao_pge_min_idx, sao_pge_max_idx
     IF ( TRIM(ADJUSTL(tmpchar)) == TRIM(ADJUSTL(sao_pge_names(i))) ) THEN
        pge_name = sao_pge_names(i)  ;  pge_idx = i
        EXIT getpge
     END IF
  END DO getpge
  IF ( pge_idx == -1 .OR.  pge_name(1:1) == '?' ) pge_errstat = pge_errstat_error

  RETURN
END SUBROUTINE get_pge_ident


SUBROUTINE pge_error_message ( pge_errstat, ok_msg, warn_msg, err_msg )

  ! ====================================================================
  ! Check PGE error status and report appropriate message. STOP on error
  ! ====================================================================

  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,           INTENT (IN) :: pge_errstat
  CHARACTER (LEN=*), INTENT (IN) :: ok_msg, warn_msg, err_msg

  SELECT CASE ( pge_errstat )
  CASE ( pge_errstat_ok )
     WRITE (*, '(A)') TRIM(ADJUSTL(ok_msg))
  CASE ( pge_errstat_warning )
     WRITE (*, '(A,A)') 'WARNING: ', TRIM(ADJUSTL(warn_msg))
  CASE ( pge_errstat_error )
     WRITE (*, '(A,A)') 'ERROR: ', TRIM(ADJUSTL(err_msg))
     STOP 1
  END SELECT

  RETURN
END SUBROUTINE pge_error_message
