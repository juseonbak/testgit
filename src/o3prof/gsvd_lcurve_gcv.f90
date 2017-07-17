  ! ***************************************************************************
  ! Author:  xiong liu (xliu)
  ! Date  :  July 24, 2003
  ! Purpose: Implement GSVD with Lcurve/GCV analysis
  !          GSVD: from NETLIB/LAPACK
  !          LCURVE: formula by Schimpf and Schrier, 1997
  !          SUBROUTINE FMIN is used to find function minimum
  ! Modification History
  ! xliu, Jan 6, 2004: implement average kernel, contribution function, 
  !   convergence criterion, information content, dfs according to OE theory
  !   writing fitting process for each iteration
  ! ***************************************************************************
  !  **********************************************************
  !  THE ROUTINE INVERT SOMTIMES DOES NOT WORK (March 1, 2004)
  !  Need to use LAPACK Routine Instead
  !  ************************************************************

  !Note: FMIN often fails to find the global minimum, need to be improved
  !Using 2-D cubic-spline interpolation is better? Considered later.

SUBROUTINE gsvd_lcurve_gcv (ozwrtint, ozwrtint_unit, epsrel, num_iter, &
     n_rad_wvl, n_fitvar_rad, ptr_nump, ptr_b, g_spec, dyda, xap, xold, xname,&
     fopt, covar, conv, avg_kernel, contri, ozdfs, ozinfo, nchisq)
				      
  USE gsvd_data_module,       ONLY : utg, xfa, x, gcv, gamma, alpha, ptr_l
  USE ozprof_data_module,     ONLY : lcurve_gcv
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variable
  ! ---------------
  LOGICAL, INTENT(IN) :: ozwrtint
  INTEGER, INTENT(IN) :: ozwrtint_unit, num_iter, n_fitvar_rad, n_rad_wvl, ptr_nump
  REAL (KIND=dp), DIMENSION(ptr_nump, n_fitvar_rad), INTENT(IN)  :: ptr_b
  REAL (KIND=dp), DIMENSION(n_rad_wvl, n_fitvar_rad), INTENT(IN) :: dyda
  REAL (KIND=dp), DIMENSION(n_rad_wvl), INTENT(IN)               :: g_spec
  REAL (KIND=dp), DIMENSION(n_fitvar_rad), INTENT (IN)           :: xold ! previous iteration
  REAL (KIND=dp), DIMENSION(n_fitvar_rad), INTENT (IN)           :: xap  ! xa - xk 
  REAL (KIND=dp), INTENT(IN)                                     :: epsrel
  CHARACTER (LEN=6), DIMENSION (n_fitvar_rad), INTENT (IN)       :: xname

  ! ---------------
  ! Output variable
  ! ---------------
  REAL (KIND=dp), INTENT (OUT), DIMENSION(n_fitvar_rad)              :: fopt
  REAL (KIND=dp), DIMENSION(n_fitvar_rad, n_fitvar_rad), INTENT(OUT) :: covar, &
	avg_kernel
  REAL (KIND=dp), DIMENSION(n_fitvar_rad, n_rad_wvl), INTENT(OUT)    :: contri
  REAL (KIND=dp), INTENT (OUT)  :: ozdfs, ozinfo, nchisq
  LOGICAL, INTENT (OUT)         :: conv

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER :: i, j, k, ptr_k, gsvd_info, inv_info
  INTEGER, DIMENSION(n_fitvar_rad) :: iwork, ipiv, work_inv
  REAL (KIND=dp), DIMENSION(n_rad_wvl, n_fitvar_rad)    :: acpy
  REAL (KIND=dp), DIMENSION(ptr_nump, n_fitvar_rad)     :: bcpy
  REAL (KIND=dp), DIMENSION(n_fitvar_rad, n_fitvar_rad) :: xinv, q, xcpy
  REAL (KIND=dp), DIMENSION(n_fitvar_rad)               :: filter, beta
  REAL (KIND=dp), DIMENSION(n_rad_wvl)                  :: g
  REAL (KIND=dp), DIMENSION(n_rad_wvl, n_rad_wvl)       :: u
  REAL (KIND=dp), DIMENSION(ptr_nump, ptr_nump)         :: v
  REAL (KIND=dp), DIMENSION(MAX(3*n_fitvar_rad,n_rad_wvl,ptr_nump)+n_fitvar_rad) :: work
  REAL (KIND=dp) :: temp, gamsq, picard, lcurve_gcv_ratio, ptr_alpha1, &
       ptr_alpha2,  OptEng, OptRho, OptCurv, OptEng2, OptRho2, OptCurv2,&
       ochisq, delchi, gtemp, ptr_alpha
  INTEGER        ::  pge_error_status = pge_errstat_ok

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=15), PARAMETER :: modulename = 'gsvd_lcurve_gcv'

  pge_error_status = pge_errstat_ok

  ! copy dyda, and ptr_b
  acpy(1:n_rad_wvl, 1:n_fitvar_rad) = dyda(1:n_rad_wvl, 1:n_fitvar_rad)
  bcpy(1:ptr_nump, 1:n_fitvar_rad) = ptr_b(1:ptr_nump, 1:n_fitvar_rad)
  g(1:n_rad_wvl) = g_spec(1:n_rad_wvl)
             
  !	Perform the GSVD operation
  CALL dggsvd( 'U', 'N', 'Q', n_rad_wvl, n_fitvar_rad, ptr_nump, ptr_k, ptr_l, &
       acpy, n_rad_wvl, bcpy, ptr_nump, alpha(1:n_fitvar_rad), beta(1:n_fitvar_rad), &
       u(1:n_rad_wvl, 1:n_rad_wvl), n_rad_wvl, v, ptr_nump, q, n_fitvar_rad, &
       work, iwork, gsvd_info)

  !	print the results
  IF (gsvd_info /= 0) THEN
     WRITE(*, *) modulename, ' : Call DGGSVD Fails !!!'
     pge_error_status = pge_errstat_error; GO TO 999
  END IF

  !	compute x, and sort the generalized singular values in increasing order
  CALL sort_gsvd (acpy, bcpy, u(1:n_rad_wvl, 1:n_rad_wvl), q, alpha(1:n_fitvar_rad), &
       beta(1:n_fitvar_rad), n_rad_wvl, n_fitvar_rad, ptr_nump, ptr_k, ptr_l, iwork, &
       x(1:n_fitvar_rad, 1:n_fitvar_rad), gamma(1:n_fitvar_rad), pge_error_status)
  
  IF (pge_error_status == pge_errstat_error)  GO TO 999

  ! invert x
  !WRITE(90, '(1p20d14.5)') ((x(i,j),j=1, n_fitvar_rad), i=1, n_fitvar_rad)
  
  xinv = x(1:n_fitvar_rad, 1:n_fitvar_rad)
  CALL dgetrf(n_fitvar_rad, n_fitvar_rad, xinv, n_fitvar_rad, ipiv, inv_info )
  IF (inv_info /= 0)  THEN
     WRITE(*, *) modulename, ' :LU of x fails and inv_info = ', inv_info
     pge_error_status = pge_errstat_error; GO TO 999;
  ELSE
     CALL dgetri(n_fitvar_rad, xinv, n_fitvar_rad, ipiv, work_inv, n_fitvar_rad, inv_info)
     IF (inv_info /= 0)  THEN
        WRITE(*, *) modulename, ' :Inverse of x fails and inv_info = ', inv_info
        pge_error_status = pge_errstat_error; GO TO 999;
     END IF
  END IF
  !WRITE(90, '(1p20d14.5)') ((xinv(i,j),j=1, n_fitvar_rad), i=1, n_fitvar_rad)

  ! This routine is faster, but without error reports.
  !CALL inverse(x(1:n_fitvar_rad, 1:n_fitvar_rad), n_fitvar_rad, xinv)
  
  ! get x^(-1)*fa
  DO i =1, n_fitvar_rad
     xfa(i) = SUM (xinv(i, 1:n_fitvar_rad) * &
          (xap(1:n_fitvar_rad) - xold(1:n_fitvar_rad)))
  END DO

  ! get Ui^T * g
  DO i = 1, n_rad_wvl
      utg(i) = SUM (u(1:n_rad_wvl, i) * g(1:n_rad_wvl))
  END DO
  
  !	Find the best regularization parameter using L-curve ir Generalized Cross Validation
  !  Lcurve-approach is more stable, but could over regularize (i.e. not extracting
  !  maximum information). While the GCV approach is not stable, but usually get better
  !  regularization level if it works
  IF (lcurve_gcv == 1) THEN       ! use lcurve only
     gcv = .FALSE.
     CALL find_lcorner (ptr_alpha, OptRho, OptEng, OptCurv, pge_error_status)
  ELSE IF (lcurve_gcv == 2) THEN  ! use gcv only
     gcv = .TRUE.
     CALL find_lcorner (ptr_alpha,  OptRho, OptEng, OptCurv, pge_error_status)
  ELSE IF (lcurve_gcv == 3) THEN  ! use gcv only unless gcv/lcurve differs by 10, use lcurve
     gcv = .TRUE.
     CALL find_lcorner (ptr_alpha1,  OptRho, OptEng, OptCurv, pge_error_status)
     
     gcv = .FALSE.
     CALL find_lcorner (ptr_alpha2,  OptRho2, OptEng2, OptCurv2, pge_error_status)
     
     lcurve_gcv_ratio = ptr_alpha2 / ptr_alpha1
     IF ( lcurve_gcv_ratio >= 10.0 .OR. lcurve_gcv_ratio <= 0.1) THEN
        ptr_alpha = ptr_alpha2; OptRho = OptRho2; OptEng = OptEng2; OptCurv = OptCurv2
     END IF
  ELSE
     WRITE(*, *) modulename, ' : This option is not implemented for lcurve_gcv!!!'
     pge_error_status = pge_errstat_error; GO TO 999
  END IF
  
  IF (pge_error_status == pge_errstat_error) THEN
     WRITE(*, *) modulename, ' :compute best regularization parameter'; GO TO 999
  END IF
          

  !WRITE(*, *) 'ptr_alpha = ', ptr_alpha
  !	Compute the solution
  CALL get_gsvd_sol(ptr_alpha, fopt, ozdfs, ozinfo)

  ! check for convergence 
  ochisq = SUM(g_spec(1:n_rad_wvl)**2)     
  nchisq = 0.0
  DO i = 1, n_rad_wvl
     gtemp = SUM(dyda(i,1:n_fitvar_rad) * fopt(1:n_fitvar_rad))
     nchisq = nchisq + (g_spec(i) - gtemp) ** 2
  END DO
  delchi = ABS(nchisq - ochisq) / ochisq 
  conv = .FALSE. ; IF (delchi <= epsrel) conv = .TRUE. 
  
  ! write output if needes
  IF (ozwrtint) THEN

     WRITE(ozwrtint_unit, '(A, I5)')    'Iteration = ', num_iter
     WRITE(ozwrtint_unit, '(A, d14.6)') 'Old Chi   = ', ochisq
     WRITE(ozwrtint_unit, '(A, d14.6)') 'New Chi   = ', nchisq
     WRITE(ozwrtint_unit, '(A, d14.6, A1, d14.6)') 'Delchi/ limit value = ', &
          delchi, '/', epsrel
     
     ! Degrees of Freedom Noise and Signal, Information content (dfn,dfs,h):
     WRITE(ozwrtint_unit, '(A, d14.6)') 'DFN       = ', n_fitvar_rad - ozdfs
     WRITE(ozwrtint_unit, '(A, d14.6)')  'O3 DFS    = ', ozdfs
     WRITE(ozwrtint_unit, '(A, d14.6)')  'Information content = ', ozinfo
     
     ! Retrieved state, previous state, a priori
     WRITE(ozwrtint_unit,'(A)') '  Var       Current       Previous      Delta    Apriori'
     DO i = 1, n_fitvar_rad
        WRITE(ozwrtint_unit, '(A6, 4d14.6)') xname(i), &
			fopt(i) + xold(i), xold(i), fopt(i), xap(i)
     END DO
     
  END IF

  ! ======================= computer general diagonastics ====================
  DO i = 1, n_fitvar_rad
     IF (i <= ptr_l) THEN
        filter(i) = gamma(i) ** 2.0 /  (gamma(i)**2.0 + ptr_alpha**2.0) / alpha(i)
     ELSE
        filter(i) = 1.0
     ENDIF
  ENDDO
  
  ! contribution function
  DO i =1, n_fitvar_rad
     DO j = 1, n_rad_wvl
        contri(i, j) = SUM(x(i, 1:n_fitvar_rad) * filter(1:n_fitvar_rad) &
             * u(j, 1:n_fitvar_rad) )
     END DO
  END DO
  
  ! Covariance matrix
  ! If dyda, g_spec already is being multiplied by measurement error, the sqrt of
  ! diagonal of covariance actually includes measurement error, otherwise, one need to
  ! multiply the rms (don't need to divide by measurement error) of the residual 
  ! (notice difference between using radiance and logarithmetic of the radiance)
  DO i = 1, n_fitvar_rad
     DO j = 1, i
        covar(i, j) =  SUM(contri(i, 1:n_rad_wvl) * contri(j, 1:n_rad_wvl))
        covar(j, i) = covar(i, j)
     END DO
  END DO
  
  ! average kernel (be careful here, the average kernel also include other parameters)
  !dfs1 = 0.0       ! cross check for both average kernel and dfs
  DO i = 1, n_fitvar_rad
     DO j = 1, n_fitvar_rad
        avg_kernel(i, j) = SUM(contri(i, 1:n_rad_wvl) * dyda(1:n_rad_wvl, j))
     END DO
     !dfs1 = dfs1 + avg_kernel(i, i)
  END DO

!  WRITE(91, *) 'rkernel1 = ', n_fitvar_rad, n_fitvar_rad
!  DO i  =1 , n_fitvar_rad
!     WRITE(91, '(18d16.7)') (avg_kernel(i, j), j=1, n_fitvar_rad)
!  ENDDO
  
  RETURN

999 STOP 'Error occur in GSVD_LCURVE_GCV'
 
END SUBROUTINE gsvd_lcurve_gcv

! =========================================================================      
!	SORT U, Q, X, ALPHA and BETA so that the generalized singular 
!	values are ordered in increasing order
! ========================================================================= 
SUBROUTINE sort_gsvd (a, b, u, q, alpha, beta, m, n, p, k, l, iwork, x, &
     gamma, pge_error_status)

  USE OMSAO_precision_module
  USE OMSAO_errstat_module
  IMPLICIT NONE

  ! ---------------
  ! Input variable
  ! ---------------
  INTEGER, INTENT(IN)                         :: p, m, n, k, l
  REAL (KIND=dp), DIMENSION(m, n), INTENT(IN) :: a
  REAL (KIND=dp), DIMENSION(p, n), INTENT(IN) :: b 
  INTEGER, DIMENSION(n),           INTENT(IN) :: iwork

  ! ---------------
  ! Output variable
  ! ---------------
  REAL (KIND=dp), DIMENSION(n), INTENT(OUT)        :: gamma
  REAL (KIND=dp), DIMENSION(n, n), INTENT(OUT)     :: x
  INTEGER, INTENT(OUT)                             :: pge_error_status

  ! -----------------
  ! modified variable
  ! -----------------
  REAL (KIND=dp), DIMENSION(m, m), INTENT(INOUT) :: u
  REAL (KIND=dp), DIMENSION(n, n), INTENT(INOUT) :: q
  REAL (KIND=dp), DIMENSION(n),    INTENT(INOUT) :: alpha, beta

  ! -----------------
  ! Local variable
  ! -----------------
  REAL (KIND=dp), DIMENSION(k+l, k+l)             :: r
  REAL (KIND=dp)                                  :: sum
  INTEGER                                         :: i, j, info, jfirst, jlast

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=9), PARAMETER :: modulename = 'sort_gsvd'

  ! initialize r to 0.0_dp
  r = 0.0

  !	obtain X so that U' A X = (0 D1) and V' B X = (0 D2)
  !	get R      
  IF (m >= (k + l))  THEN 
     r(1:k+l, 1:k+l) = a(1:k+l, n-k-l+1:n)
  ELSE  ! m < k + l  
     r(1:m, 1:k+l) = a(1:m, n-k-l+1:n)
     r(m+1:k+l, m+1:k+l) = b(m-k+1:l, n+m-k-l+1:n)
  END IF

  CALL dtrtri( 'U', 'N', k + l, r,  k + l, info )

  IF (info /= 0) THEN
     WRITE(*, *) modulename, ' : Call dtrtri Fails !!!'
     pge_error_status = pge_errstat_error; RETURN
  END IF

  !	compute X
  jfirst = n - (k + l)
  jlast = jfirst + 1

  x(1:n, 1:jfirst) = q(1:n, 1:jfirst)
  x(1:n, jlast:n) = MATMUL(q(1:n, jlast:n), r)

  !	re-sort x, u, v, based on the index stored in iwork for alpha/beta
  !	save space, but less efficient, may need to be modified later
  ! Note: here could not be done in parallel
  DO i = 1, k + l
     IF ( i /= iwork (i) ) THEN
        CALL swap(alpha(i), alpha(iwork(i)))
        CALL swap(beta(i), beta(iwork(i)))
         
        DO j = 1, m
           CALL swap(u(j, i), u(j, iwork(i))) 
        END DO
        
        DO j = 1, n
           CALL swap(q(j, i), q(j, iwork(i))) 
           CALL swap(x(j, i), x(j, iwork(i))) 
        END DO
     END IF
  END DO

  DO i = 1, l
     gamma(i) = alpha(i) / beta (i)
  END DO

  RETURN

END SUBROUTINE sort_gsvd

! =========================================================================       
! Subroutine to get lcurve/gcv points including (lcurve curvature/GCV value) 
! for a given regularization parameter
! ========================================================================= 
SUBROUTINE scurv(rlamda, TheCurv, TheRho, TheEng)

  !	u: GSVD for the kernel matrix (n_rad_wvl x n_rad_wvl)
  !	g: remainding measurement vector ( n_rad_wvl )
  !	n_rad_wvl, n_fitvar_rad :  dimension of kernel (n_rad_wvl x n_fitvar_rad)
  !	ptr_l: effective rank of b (ptr_l <= ptr_nump)
  !	gamma: generalized singular values  (n_fitvar_rad)
  !	rlamda:    input regularization parameter 

  USE gsvd_data_module,       ONLY : utg, xfa, gcv, gamma, ptr_l, dp, &
       n_fitvar_rad, n_rad_wvl, alpha
  IMPLICIT NONE

  ! ---------------
  ! Input variable
  ! ---------------
  REAL (KIND=dp), INTENT (IN) :: rlamda

  ! ---------------
  ! Result variable
  ! ---------------	 
  REAL (KIND=dp), INTENT (OUT):: TheCurv, TheEng, TheRho

  ! ---------------
  ! Lcoal variables
  ! ---------------	
  INTEGER                     :: i  
  REAL (KIND=dp)              :: TheEngfd, picard, picard1, picard2, &
	sumf, denom, gamsq, lamsq

  !	Initialize TheRho, TheEng, TheEngfd
  TheRho = 0.0
  TheEng = 0.0
  TheEngfd = 0.0
  sumf = 0.0  
  lamsq = rlamda * rlamda   
   
  !	Compute TheRho, TheEng, TheEngfd
  IF (.NOT. gcv) THEN

     DO i = 1, ptr_l
        picard = (utg(i) - xfa(i) *alpha(i))**2.0
        gamsq = gamma(i) * gamma(i)
        denom = gamsq + lamsq         
        sumf = sumf  + gamsq / denom
        
        TheRho = TheRho + picard * (lamsq / denom) ** 2.0
        TheEng = TheEng + picard * gamsq / (denom ** 2.0 )
        TheEngfd = TheEngfd - 4.0 * rlamda * gamsq * picard / (denom ** 3.0)		    
     END DO

     TheRho = TheRho + SUM(utg(n_fitvar_rad+1:n_rad_wvl)**2.0)
 
     TheCurv =  TheEng * TheRho / TheEngfd * ( lamsq *   & 
          TheEngfd * TheRho + 2.0 * rlamda * TheEng * TheRho  + &
          lamsq ** 2.0 * TheEng * TheEngfd ) /                 &
          (lamsq ** 2.0 * TheEng * TheEng + TheRho * TheRho) ** 1.5

     TheRho = LOG10(TheRho)
     TheEng = LOG10(TheEng)

  ELSE      ! using GCV
     DO i = 1, ptr_l
        picard = (utg(i) - xfa(i) *alpha(i))**2.0
		gamsq = gamma(i) * gamma(i)
        denom = gamsq + lamsq         
        sumf = sumf  + gamsq / denom
        TheRho = TheRho + TheRho + picard * (lamsq / denom) ** 2.0
     END DO

     TheRho = TheRho + SUM(utg(n_fitvar_rad+1:n_rad_wvl)**2.0)
 
     TheCurv = TheRho * n_rad_wvl / ( 0.0 + n_rad_wvl - &
          n_fitvar_rad + ptr_l - sumf) ** 2.0
    
     TheRho = 0.0
     TheEng = 0.0
  END IF
  
  RETURN
END SUBROUTINE scurv

! ========================================================================= 
! Function to calculate curvature for a given regularization parameter
! ========================================================================= 
FUNCTION fcurv(rlamda) RESULT (TheCurv)

  !	u: GSVD for the kernel matrix (n_rad_wvl x n_rad_wvl)
  !	g: remainding measurement vector ( n_rad_wvl )
  !	n_rad_wvl, n_fitvar_rad :  dimension of kernel (n_rad_wvl x n_fitvar_rad)
  !	ptr_l: effective rank of b (ptr_l <= ptr_nump)
  !	gamma: generalized singular values  (n_fitvar_rad)
  !	rlamda:    input regularization parameter 

  USE gsvd_data_module,       ONLY : utg, xfa, gcv, gamma, alpha, ptr_l, dp, &
       n_fitvar_rad, n_rad_wvl
  IMPLICIT NONE

  ! ---------------
  ! Input variable
  ! ---------------
  REAL (KIND=dp), INTENT (IN) :: rlamda

  ! ---------------
  ! Result variable
  ! ---------------	 
  REAL (KIND=dp)              :: TheCurv

  ! ---------------
  ! Lcoal variables
  ! ---------------	
  INTEGER                     :: i, j
  REAL (KIND=dp)              :: TheRho, TheEng, TheEngfd, picard, &
       sumf, denom, gamsq, lamsq

  !	Initialize TheRho, TheEng, TheEngfd
  TheRho = 0.0
  TheEng = 0.0
  TheEngfd = 0.0
  sumf = 0.0
  lamsq = rlamda * rlamda   
   
  !	Compute TheRho, TheEng, TheEngfd
  IF (.NOT. gcv) THEN

     DO i = 1, ptr_l
	     picard = (utg(i) - xfa(i) *alpha(i))**2.0
         gamsq = gamma(i) * gamma(i)
         denom = gamsq + lamsq         
         sumf = sumf  + gamsq / denom

         TheRho = TheRho + picard * (lamsq / denom) ** 2.0
         TheEng = TheEng + picard * gamsq / (denom ** 2.0 )
         TheEngfd = TheEngfd - 4.0 * rlamda * gamsq * picard / (denom ** 3.0)		    
     END DO

     TheRho = TheRho + SUM(utg(n_fitvar_rad+1:n_rad_wvl)**2.0)
 
     TheCurv =  TheEng * TheRho / TheEngfd * ( lamsq *   & 
          TheEngfd * TheRho + 2.0 * rlamda * TheEng * TheRho  + &
          lamsq ** 2.0 * TheEng * TheEngfd ) /                 &
          (lamsq ** 2.0 * TheEng * TheEng + TheRho * TheRho) ** 1.5

  ELSE  ! using GCV

     DO i = 1, ptr_l
        picard = (utg(i) - xfa(i) *alpha(i))**2
		gamsq = gamma(i) * gamma(i)
        denom = gamsq + lamsq         
        sumf = sumf  + gamsq / denom
        TheRho = TheRho + TheRho + picard * (lamsq / denom) ** 2.0
     END DO

     TheRho = TheRho + SUM(utg(n_fitvar_rad+1:n_rad_wvl)**2)
 
     TheCurv = TheRho * n_rad_wvl / ( 0.0 + n_rad_wvl - &
          n_fitvar_rad + ptr_l - sumf) ** 2.0

  END IF
  
  RETURN

END FUNCTION fcurv

! ========================================================================= 
! Subroutine to find the points of maximum curvature for lcurve or
! minimize GCV value for GCV
! ========================================================================= 
SUBROUTINE find_lcorner(ptr_alpha, OptRho, OptEng, OptCurv, pge_error_status)

  !	u: GSVD for the kernel matrix (n_rad_wvl x n_rad_wvl)
  !	g: remainding measurement vector ( n_rad_wvl )
  !	n_rad_wvl, n_fitvar_rad :  dimension of kernel (n_rad_wvl x n_fitvar_rad)
  !	ptr_l: effective rank of b (ptr_l <= ptr_nump)
  !	gamma: generalized singular values  (n_fitvar_rad)

  USE gsvd_data_module,        ONLY : gamma, ptr_l
  USE OMSAO_parameters_module, ONLY : lcurve_tol
  USE ozprof_data_module,      ONLY : lcurve_write, lcurve_unit
  USE OMSAO_errstat_module
  IMPLICIT NONE

  ! ---------------
  ! Output variable
  ! ---------------	 
  INTEGER,        INTENT (OUT) :: pge_error_status
  REAL (KIND=dp), INTENT (OUT) :: ptr_alpha, OptEng, OptRho, OptCurv

  ! ---------------
  ! local variable
  ! ---------------	 
  INTEGER, PARAMETER           :: maxlam = 100
  INTEGER                      :: MinIn, i
  REAL (KIND=dp), DIMENSION(maxlam) :: curv, rho, eng, RegPar
  REAL (KIND=dp) :: CurvMin, gmax, gmin, DeltaLam, alam, blam, fmin, OptLam, fcurv

  EXTERNAL fcurv

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=12), PARAMETER :: modulename = 'find_lcorner'

  !	Check to make sure that ptr_l is smaller than maxlam
  IF (ptr_l > maxlam) THEN
     WRITE(*, *) modulename, ' : Find_lcorner: Increase maxlam to at least ',  ptr_l
     pge_error_status = pge_errstat_error; RETURN
  END IF

  ! find the maximum curvature at points of singular values
  ! reduce the range for searching the global maximum curvature
  ! this is done because that the FMIN may find the local maximum
  !curvature if the function (curvature) is not unimodal
  CurvMin= 10.E10
  gmax = LOG10 (gamma(ptr_l))
  gmin = LOG10 (gamma(1))
  DeltaLam = (gmax - gmin) / (maxlam - 1.0)
 
  DO i = 1, maxlam
     RegPar (i) = 10.0 ** ( gmin + DeltaLam * i)
     curv(i) = fcurv(RegPar(i))
     IF (curv (i) < CurvMin) THEN
        CurvMin = curv (i)
        MinIn = i
     END IF
  END DO
     
  IF ( MinIn == 1 ) THEN
     alam = RegPar (MinIn)
  ELSE
     alam = RegPar (MinIn - 1)
  END IF

  IF ( MinIn == maxlam ) THEN
     blam = RegPar (MinIn)
  ELSE
     blam = RegPar (MinIn + 1)
  END IF

  ! identify te corner using FMIN routine
  ptr_alpha = fmin (alam, blam, fcurv, lcurve_tol)
  
  ! compute L-curve by calculating for further examination
  IF (lcurve_write) THEN
     DO i = 1, maxlam
        CALL scurv(RegPar(i), curv(i), rho(i), eng(i))
     END DO
   
     CALL scurv(ptr_alpha, Optcurv, OptRho, OptEng)
   
     WRITE (lcurve_unit, *) 'L-curve/GCV points between GSVs: '
     DO i = 1, maxlam
        WRITE (lcurve_unit, '(4D16.7)') RegPar(i), rho(i), eng(i), curv(i)
     END DO

     WRITE (lcurve_unit, *) 'L-curve/GCV point of maximum curvature/minimum GCV value: '
     WRITE (lcurve_unit, '(4D16.7)') ptr_alpha, OptRho, OptEng, OptCurv
  END IF
   
  RETURN
END SUBROUTINE find_lcorner

! ====================
! swap two data items
! ====================
SUBROUTINE SWAP (c, d)
  
  USE OMSAO_precision_module
  IMPLICIT NONE

  REAL (KIND=dp), INTENT(INOUT) :: c, d
  REAL(KIND=dp) :: temp

  temp = c
  c = d
  d = temp

  RETURN
END SUBROUTINE SWAP


! ============================================================
!  Get solution and degree of freedom for signal, and shannon 
!  information content given a regularization parameter, 
! ============================================================

SUBROUTINE get_gsvd_sol(ptr_alpha, sol, dfs, infoh)

  USE OMSAO_precision_module
  USE gsvd_data_module,       ONLY : utg, xfa, x, gamma, alpha, ptr_l
  USE OMSAO_variables_module, ONLY : n_fitvar_rad, n_rad_wvl

  IMPLICIT NONE

  ! ---------------
  ! Input variable
  ! ---------------
  REAL (KIND=dp), INTENT(IN)                             :: ptr_alpha

  ! ---------------
  ! Output variable
  ! ---------------
  REAL (KIND=dp), INTENT (OUT), DIMENSION(n_fitvar_rad)  :: sol
  REAL (KIND=dp), INTENT(OUT)                            :: dfs, infoh

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER         :: i
  REAL (KIND=dp)  :: gamsq, lamsq, picard, filter, denom
 

  !	Inititalize the soultion
  sol(1:n_fitvar_rad) = 0.0;  dfs = 0.0; infoh = 0.0

  DO i = 1, n_fitvar_rad
     picard = utg(i)

     IF (i <= ptr_l) THEN
        gamsq  = gamma(i) ** 2.0
        lamsq = ptr_alpha ** 2.0
		denom = gamsq + lamsq
        filter = gamsq / denom
        dfs = dfs + filter
        infoh = infoh + 0.5 * LOG(1.0 + gamsq / lamsq)
        picard = picard * filter / alpha(i) + xfa(i) * lamsq / denom
     END IF
     
     sol(1:n_fitvar_rad) = sol(1:n_fitvar_rad)  + picard * &
		x(1:n_fitvar_rad, i)
  END DO
       
  RETURN

END SUBROUTINE get_gsvd_sol


! --------------------------------------------------------------------
SUBROUTINE inverse (a, n, b) ! Invert matrix by Gauss method
! --------------------------------------------------------------------
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL(8), DIMENSION(n, n), INTENT(IN)  :: a
REAL(8), DIMENSION(n, n), INTENT(OUT) :: b

! - - - Local Variables - - -
REAL(8) :: c, d, temp(n)
INTEGER :: i, j, k, m, imax(1), ipvt(n)
! - - - - - - - - - - - - - -

b = a
ipvt = (/ (i, i = 1, n) /)

DO k = 1,n
imax = MAXLOC(ABS(b(k:n,k)))
m = k-1+imax(1)

IF (m /= k) THEN
ipvt( (/m,k/) ) = ipvt( (/k,m/) )
b((/m,k/),:) = b((/k,m/),:)
END IF
d = 1/b(k,k)

temp = b(:,k)
DO j = 1, n
c = b(k,j)*d
b(:,j) = b(:,j)-temp*c
b(k,j) = c
END DO
b(:,k) = temp*(-d)
b(k,k) = d
END DO

b(:,ipvt) = b

END SUBROUTINE inverse
