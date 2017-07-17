
      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRI computes the inverse of a matrix using the LU factorization
*  computed by DGETRF.
*
*  This method inverts U and then computes inv(A) by solving the system
*  inv(A)*L = inv(U) for inv(A).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the factors L and U from the factorization
*          A = P*L*U as computed by DGETRF.
*          On exit, if INFO = 0, the inverse of the original matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimal performance LWORK >= N*NB, where NB is
*          the optimal blocksize returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
*                singular and its inverse could not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB,
     $                   NBMIN, NN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGEMV, DSWAP, DTRSM, DTRTRI, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NB = ILAENV( 1, 'DGETRI', ' ', N, -1, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRI', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Form inv(U).  If INFO > 0 from DTRTRI, then U is singular,
*     and the inverse is not computed.
*
      CALL DTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
      IF( INFO.GT.0 )
     $   RETURN
*
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = MAX( LDWORK*NB, 1 )
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DGETRI', ' ', N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = N
      END IF
*
*     Solve the equation inv(A)*L = inv(U) for inv(A).
*
      IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         DO 20 J = N, 1, -1
*
*           Copy current column of L to WORK and replace with zeros.
*
            DO 10 I = J + 1, N
               WORK( I ) = A( I, J )
               A( I, J ) = ZERO
   10       CONTINUE
*
*           Compute current column of inv(A).
*
            IF( J.LT.N )
     $         CALL DGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
   20    CONTINUE
      ELSE
*
*        Use blocked code.
*
         NN = ( ( N-1 ) / NB )*NB + 1
         DO 50 J = NN, 1, -NB
            JB = MIN( NB, N-J+1 )
*
*           Copy current block column of L to WORK and replace with
*           zeros.
*
            DO 40 JJ = J, J + JB - 1
               DO 30 I = JJ + 1, N
                  WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
                  A( I, JJ ) = ZERO
   30          CONTINUE
   40       CONTINUE
*
*           Compute current block column of inv(A).
*
            IF( J+JB.LE.N )
     $         CALL DGEMM( 'No transpose', 'No transpose', N, JB,
     $                     N-J-JB+1, -ONE, A( 1, J+JB ), LDA,
     $                     WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
            CALL DTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,
     $                  ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
   50    CONTINUE
      END IF
*
*     Apply column interchanges.
*
      DO 60 J = N - 1, 1, -1
         JP = IPIV( J )
         IF( JP.NE.J )
     $      CALL DSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
   60 CONTINUE
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DGETRI
*
      END


      SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRTRI computes the inverse of a real upper or lower triangular
*  matrix A.
*
*  This is the Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*
*  DIAG    (input) CHARACTER*1
*          = 'N':  A is non-unit triangular;
*          = 'U':  A is unit triangular.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
*               matrix is singular and its inverse can not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J, JB, NB, NN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTRMM, DTRSM, DTRTI2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check for singularity if non-unit.
*
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
         INFO = 0
      END IF
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DTRTRI', UPLO // DIAG, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
      ELSE
*
*        Use blocked code
*
         IF( UPPER ) THEN
*
*           Compute inverse of upper triangular matrix
*
            DO 20 J = 1, N, NB
               JB = MIN( NB, N-J+1 )
*
*              Compute rows 1:j-1 of current block column
*
               CALL DTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, ONE, A, LDA, A( 1, J ), LDA )
               CALL DTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
*
*              Compute inverse of current diagonal block
*
               CALL DTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
   20       CONTINUE
         ELSE
*
*           Compute inverse of lower triangular matrix
*
            NN = ( ( N-1 ) / NB )*NB + 1
            DO 30 J = NN, 1, -NB
               JB = MIN( NB, N-J+1 )
               IF( J+JB.LE.N ) THEN
*
*                 Compute rows j+jb:n of current block column
*
                  CALL DTRMM( 'Left', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA,
     $                        A( J+JB, J ), LDA )
                  CALL DTRSM( 'Right', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, -ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
*
*              Compute inverse of current diagonal block
*
               CALL DTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
   30       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DTRTRI
*
      END

      SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRTI2 computes the inverse of a real upper or lower triangular
*  matrix.
*
*  This is the Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A is unit triangular.
*          = 'N':  Non-unit triangular
*          = 'U':  Unit triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading n by n upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DTRMV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTI2', -INFO )
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Compute inverse of upper triangular matrix.
*
         DO 10 J = 1, N
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
*
*           Compute elements 1:j-1 of j-th column.
*
            CALL DTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA,
     $                  A( 1, J ), 1 )
            CALL DSCAL( J-1, AJJ, A( 1, J ), 1 )
   10    CONTINUE
      ELSE
*
*        Compute inverse of lower triangular matrix.
*
         DO 20 J = N, 1, -1
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
            IF( J.LT.N ) THEN
*
*              Compute elements j+1:n of j-th column.
*
               CALL DTRMV( 'Lower', 'No transpose', DIAG, N-J,
     $                     A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
               CALL DSCAL( N-J, AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of DTRTI2
*
      END

      SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORM2R overwrites the general real m by n matrix C with
*
*        Q * C  if SIDE = 'L' and TRANS = 'N', or
*
*        Q'* C  if SIDE = 'L' and TRANS = 'T', or
*
*        C * Q  if SIDE = 'R' and TRANS = 'N', or
*
*        C * Q' if SIDE = 'R' and TRANS = 'T',
*
*  where Q is a real orthogonal matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q' from the Left
*          = 'R': apply Q or Q' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply Q  (No transpose)
*          = 'T': apply Q' (Transpose)
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
*          The i-th column must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          DGEQRF in the first k columns of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If SIDE = 'L', LDA >= max(1,M);
*          if SIDE = 'R', LDA >= max(1,N).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQRF.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension
*                                   (N) if SIDE = 'L',
*                                   (M) if SIDE = 'R'
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      DOUBLE PRECISION   AII
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
*
*     NQ is the order of Q
*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORM2R', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN
*
      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) )
     $     THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
*
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
*
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
*
*           H(i) is applied to C(i:m,1:n)
*
            MI = M - I + 1
            IC = I
         ELSE
*
*           H(i) is applied to C(1:m,i:n)
*
            NI = N - I + 1
            JC = I
         END IF
*
*        Apply H(i)
*
         AII = A( I, I )
         A( I, I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ),
     $               LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
*
*     End of DORM2R
*
      END
      SUBROUTINE DORMR2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORMR2 overwrites the general real m by n matrix C with
*
*        Q * C  if SIDE = 'L' and TRANS = 'N', or
*
*        Q'* C  if SIDE = 'L' and TRANS = 'T', or
*
*        C * Q  if SIDE = 'R' and TRANS = 'N', or
*
*        C * Q' if SIDE = 'R' and TRANS = 'T',
*
*  where Q is a real orthogonal matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by DGERQF. Q is of order m if SIDE = 'L' and of order n
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q' from the Left
*          = 'R': apply Q or Q' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply Q  (No transpose)
*          = 'T': apply Q' (Transpose)
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension
*                               (LDA,M) if SIDE = 'L',
*                               (LDA,N) if SIDE = 'R'
*          The i-th row must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          DGERQF in the last k rows of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,K).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGERQF.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension
*                                   (N) if SIDE = 'L',
*                                   (M) if SIDE = 'R'
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, MI, NI, NQ
      DOUBLE PRECISION   AII
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
*
*     NQ is the order of Q
*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMR2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN
*
      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) )
     $     THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
*
      IF( LEFT ) THEN
         NI = N
      ELSE
         MI = M
      END IF
*
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
*
*           H(i) is applied to C(1:m-k+i,1:n)
*
            MI = M - K + I
         ELSE
*
*           H(i) is applied to C(1:m,1:n-k+i)
*
            NI = N - K + I
         END IF
*
*        Apply H(i)
*
         AII = A( I, NQ-K+I )
         A( I, NQ-K+I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( I, 1 ), LDA, TAU( I ), C, LDC,
     $               WORK )
         A( I, NQ-K+I ) = AII
   10 CONTINUE
      RETURN
*
*     End of DORMR2
*
      END
      SUBROUTINE DTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B,
     $                   LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV,
     $                   Q, LDQ, WORK, NCYCLE, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBQ, JOBU, JOBV
      INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N,
     $                   NCYCLE, P
      DOUBLE PRECISION   TOLA, TOLB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), Q( LDQ, * ), U( LDU, * ),
     $                   V( LDV, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DTGSJA computes the generalized singular value decomposition (GSVD)
*  of two real upper triangular (or trapezoidal) matrices A and B.
*
*  On entry, it is assumed that matrices A and B have the following
*  forms, which may be obtained by the preprocessing subroutine DGGSVP
*  from a general M-by-N matrix A and P-by-N matrix B:
*
*               N-K-L  K    L
*     A =    K ( 0    A12  A13 ) if M-K-L >= 0;
*            L ( 0     0   A23 )
*        M-K-L ( 0     0    0  )
*
*             N-K-L  K    L
*     A =  K ( 0    A12  A13 ) if M-K-L < 0;
*        M-K ( 0     0   A23 )
*
*             N-K-L  K    L
*     B =  L ( 0     0   B13 )
*        P-L ( 0     0    0  )
*
*  where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
*  upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
*  otherwise A23 is (M-K)-by-L upper trapezoidal.
*
*  On exit,
*
*              U'*A*Q = D1*( 0 R ),    V'*B*Q = D2*( 0 R ),
*
*  where U, V and Q are orthogonal matrices, Z' denotes the transpose
*  of Z, R is a nonsingular upper triangular matrix, and D1 and D2 are
*  ``diagonal'' matrices, which are of the following structures:
*
*  If M-K-L >= 0,
*
*                      K  L
*         D1 =     K ( I  0 )
*                  L ( 0  C )
*              M-K-L ( 0  0 )
*
*                    K  L
*         D2 = L   ( 0  S )
*              P-L ( 0  0 )
*
*                 N-K-L  K    L
*    ( 0 R ) = K (  0   R11  R12 ) K
*              L (  0    0   R22 ) L
*
*  where
*
*    C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
*    S = diag( BETA(K+1),  ... , BETA(K+L) ),
*    C**2 + S**2 = I.
*
*    R is stored in A(1:K+L,N-K-L+1:N) on exit.
*
*  If M-K-L < 0,
*
*                 K M-K K+L-M
*      D1 =   K ( I  0    0   )
*           M-K ( 0  C    0   )
*
*                   K M-K K+L-M
*      D2 =   M-K ( 0  S    0   )
*           K+L-M ( 0  0    I   )
*             P-L ( 0  0    0   )
*
*                 N-K-L  K   M-K  K+L-M
* ( 0 R ) =    K ( 0    R11  R12  R13  )
*            M-K ( 0     0   R22  R23  )
*          K+L-M ( 0     0    0   R33  )
*
*  where
*  C = diag( ALPHA(K+1), ... , ALPHA(M) ),
*  S = diag( BETA(K+1),  ... , BETA(M) ),
*  C**2 + S**2 = I.
*
*  R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored
*      (  0  R22 R23 )
*  in B(M-K+1:L,N+M-K-L+1:N) on exit.
*
*  The computation of the orthogonal transformation matrices U, V or Q
*  is optional.  These matrices may either be formed explicitly, or they
*  may be postmultiplied into input matrices U1, V1, or Q1.
*
*  Arguments
*  =========
*
*  JOBU    (input) CHARACTER*1
*          = 'U':  U must contain an orthogonal matrix U1 on entry, and
*                  the product U1*U is returned;
*          = 'I':  U is initialized to the unit matrix, and the
*                  orthogonal matrix U is returned;
*          = 'N':  U is not computed.
*
*  JOBV    (input) CHARACTER*1
*          = 'V':  V must contain an orthogonal matrix V1 on entry, and
*                  the product V1*V is returned;
*          = 'I':  V is initialized to the unit matrix, and the
*                  orthogonal matrix V is returned;
*          = 'N':  V is not computed.
*
*  JOBQ    (input) CHARACTER*1
*          = 'Q':  Q must contain an orthogonal matrix Q1 on entry, and
*                  the product Q1*Q is returned;
*          = 'I':  Q is initialized to the unit matrix, and the
*                  orthogonal matrix Q is returned;
*          = 'N':  Q is not computed.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  P       (input) INTEGER
*          The number of rows of the matrix B.  P >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and B.  N >= 0.
*
*  K       (input) INTEGER
*  L       (input) INTEGER
*          K and L specify the subblocks in the input matrices A and B:
*          A23 = A(K+1:MIN(K+L,M),N-L+1:N) and B13 = B(1:L,N-L+1:N)
*          of A and B, whose GSVD is going to be computed by DTGSJA.
*          See Further details.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, A(N-K+1:N,1:MIN(K+L,M) ) contains the triangular
*          matrix R or part of R.  See Purpose for details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
*          On entry, the P-by-N matrix B.
*          On exit, if necessary, B(M-K+1:L,N+M-K-L+1:N) contains
*          a part of R.  See Purpose for details.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,P).
*
*  TOLA    (input) DOUBLE PRECISION
*  TOLB    (input) DOUBLE PRECISION
*          TOLA and TOLB are the convergence criteria for the Jacobi-
*          Kogbetliantz iteration procedure. Generally, they are the
*          same as used in the preprocessing step, say
*              TOLA = max(M,N)*norm(A)*MAZHEPS,
*              TOLB = max(P,N)*norm(B)*MAZHEPS.
*
*  ALPHA   (output) DOUBLE PRECISION array, dimension (N)
*  BETA    (output) DOUBLE PRECISION array, dimension (N)
*          On exit, ALPHA and BETA contain the generalized singular
*          value pairs of A and B;
*            ALPHA(1:K) = 1,
*            BETA(1:K)  = 0,
*          and if M-K-L >= 0,
*            ALPHA(K+1:K+L) = diag(C),
*            BETA(K+1:K+L)  = diag(S),
*          or if M-K-L < 0,
*            ALPHA(K+1:M)= C, ALPHA(M+1:K+L)= 0
*            BETA(K+1:M) = S, BETA(M+1:K+L) = 1.
*          Furthermore, if K+L < N,
*            ALPHA(K+L+1:N) = 0 and
*            BETA(K+L+1:N)  = 0.
*
*  U       (input/output) DOUBLE PRECISION array, dimension (LDU,M)
*          On entry, if JOBU = 'U', U must contain a matrix U1 (usually
*          the orthogonal matrix returned by DGGSVP).
*          On exit,
*          if JOBU = 'I', U contains the orthogonal matrix U;
*          if JOBU = 'U', U contains the product U1*U.
*          If JOBU = 'N', U is not referenced.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U. LDU >= max(1,M) if
*          JOBU = 'U'; LDU >= 1 otherwise.
*
*  V       (input/output) DOUBLE PRECISION array, dimension (LDV,P)
*          On entry, if JOBV = 'V', V must contain a matrix V1 (usually
*          the orthogonal matrix returned by DGGSVP).
*          On exit,
*          if JOBV = 'I', V contains the orthogonal matrix V;
*          if JOBV = 'V', V contains the product V1*V.
*          If JOBV = 'N', V is not referenced.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V. LDV >= max(1,P) if
*          JOBV = 'V'; LDV >= 1 otherwise.
*
*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
*          On entry, if JOBQ = 'Q', Q must contain a matrix Q1 (usually
*          the orthogonal matrix returned by DGGSVP).
*          On exit,
*          if JOBQ = 'I', Q contains the orthogonal matrix Q;
*          if JOBQ = 'Q', Q contains the product Q1*Q.
*          If JOBQ = 'N', Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q. LDQ >= max(1,N) if
*          JOBQ = 'Q'; LDQ >= 1 otherwise.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)
*
*  NCYCLE  (output) INTEGER
*          The number of cycles required for convergence.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          = 1:  the procedure does not converge after MAXIT cycles.
*
*  Internal Parameters
*  ===================
*
*  MAXIT   INTEGER
*          MAXIT specifies the total loops that the iterative procedure
*          may take. If after MAXIT cycles, the routine fails to
*          converge, we return INFO = 1.
*
*  Further Details
*  ===============
*
*  DTGSJA essentially uses a variant of Kogbetliantz algorithm to reduce
*  min(L,M-K)-by-L triangular (or trapezoidal) matrix A23 and L-by-L
*  matrix B13 to the form:
*
*           U1'*A13*Q1 = C1*R1; V1'*B13*Q1 = S1*R1,
*
*  where U1, V1 and Q1 are orthogonal matrix, and Z' is the transpose
*  of Z.  C1 and S1 are diagonal matrices satisfying
*
*                C1**2 + S1**2 = I,
*
*  and R1 is an L-by-L nonsingular upper triangular matrix.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 40 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
*
      LOGICAL            INITQ, INITU, INITV, UPPER, WANTQ, WANTU, WANTV
      INTEGER            I, J, KCYCLE
      DOUBLE PRECISION   A1, A2, A3, B1, B2, B3, CSQ, CSU, CSV, ERROR,
     $                   GAMMA, RWK, SNQ, SNU, SNV, SSMIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAGS2, DLAPLL, DLARTG, DLASET, DROT,
     $                   DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      INITU = LSAME( JOBU, 'I' )
      WANTU = INITU .OR. LSAME( JOBU, 'U' )
*
      INITV = LSAME( JOBV, 'I' )
      WANTV = INITV .OR. LSAME( JOBV, 'V' )
*
      INITQ = LSAME( JOBQ, 'I' )
      WANTQ = INITQ .OR. LSAME( JOBQ, 'Q' )
*
      INFO = 0
      IF( .NOT.( INITU .OR. WANTU .OR. LSAME( JOBU, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( INITV .OR. WANTV .OR. LSAME( JOBV, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( INITQ .OR. WANTQ .OR. LSAME( JOBQ, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( N.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDU.LT.1 .OR. ( WANTU .AND. LDU.LT.M ) ) THEN
         INFO = -18
      ELSE IF( LDV.LT.1 .OR. ( WANTV .AND. LDV.LT.P ) ) THEN
         INFO = -20
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -22
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTGSJA', -INFO )
         RETURN
      END IF
*
*     Initialize U, V and Q, if necessary
*
      IF( INITU )
     $   CALL DLASET( 'Full', M, M, ZERO, ONE, U, LDU )
      IF( INITV )
     $   CALL DLASET( 'Full', P, P, ZERO, ONE, V, LDV )
      IF( INITQ )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
*
*     Loop until convergence
*
      UPPER = .FALSE.
      DO 40 KCYCLE = 1, MAXIT
*
         UPPER = .NOT.UPPER
*
         DO 20 I = 1, L - 1
            DO 10 J = I + 1, L
*
               A1 = ZERO
               A2 = ZERO
               A3 = ZERO
               IF( K+I.LE.M )
     $            A1 = A( K+I, N-L+I )
               IF( K+J.LE.M )
     $            A3 = A( K+J, N-L+J )
*
               B1 = B( I, N-L+I )
               B3 = B( J, N-L+J )
*
               IF( UPPER ) THEN
                  IF( K+I.LE.M )
     $               A2 = A( K+I, N-L+J )
                  B2 = B( I, N-L+J )
               ELSE
                  IF( K+J.LE.M )
     $               A2 = A( K+J, N-L+I )
                  B2 = B( J, N-L+I )
               END IF
*
               CALL DLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU,
     $                      CSV, SNV, CSQ, SNQ )
*
*              Update (K+I)-th and (K+J)-th rows of matrix A: U'*A
*
               IF( K+J.LE.M )
     $            CALL DROT( L, A( K+J, N-L+1 ), LDA, A( K+I, N-L+1 ),
     $                       LDA, CSU, SNU )
*
*              Update I-th and J-th rows of matrix B: V'*B
*
               CALL DROT( L, B( J, N-L+1 ), LDB, B( I, N-L+1 ), LDB,
     $                    CSV, SNV )
*
*              Update (N-L+I)-th and (N-L+J)-th columns of matrices
*              A and B: A*Q and B*Q
*
               CALL DROT( MIN( K+L, M ), A( 1, N-L+J ), 1,
     $                    A( 1, N-L+I ), 1, CSQ, SNQ )
*
               CALL DROT( L, B( 1, N-L+J ), 1, B( 1, N-L+I ), 1, CSQ,
     $                    SNQ )
*
               IF( UPPER ) THEN
                  IF( K+I.LE.M )
     $               A( K+I, N-L+J ) = ZERO
                  B( I, N-L+J ) = ZERO
               ELSE
                  IF( K+J.LE.M )
     $               A( K+J, N-L+I ) = ZERO
                  B( J, N-L+I ) = ZERO
               END IF
*
*              Update orthogonal matrices U, V, Q, if desired.
*
               IF( WANTU .AND. K+J.LE.M )
     $            CALL DROT( M, U( 1, K+J ), 1, U( 1, K+I ), 1, CSU,
     $                       SNU )
*
               IF( WANTV )
     $            CALL DROT( P, V( 1, J ), 1, V( 1, I ), 1, CSV, SNV )
*
               IF( WANTQ )
     $            CALL DROT( N, Q( 1, N-L+J ), 1, Q( 1, N-L+I ), 1, CSQ,
     $                       SNQ )
*
   10       CONTINUE
   20    CONTINUE
*
         IF( .NOT.UPPER ) THEN
*
*           The matrices A13 and B13 were lower triangular at the start
*           of the cycle, and are now upper triangular.
*
*           Convergence test: test the parallelism of the corresponding
*           rows of A and B.
*
            ERROR = ZERO
            DO 30 I = 1, MIN( L, M-K )
               CALL DCOPY( L-I+1, A( K+I, N-L+I ), LDA, WORK, 1 )
               CALL DCOPY( L-I+1, B( I, N-L+I ), LDB, WORK( L+1 ), 1 )
               CALL DLAPLL( L-I+1, WORK, 1, WORK( L+1 ), 1, SSMIN )
               ERROR = MAX( ERROR, SSMIN )
   30       CONTINUE
*
            IF( ABS( ERROR ).LE.MIN( TOLA, TOLB ) )
     $         GO TO 50
         END IF
*
*        End of cycle loop
*
   40 CONTINUE
*
*     The algorithm has not converged after MAXIT cycles.
*
      INFO = 1
      GO TO 100
*
   50 CONTINUE
*
*     If ERROR <= MIN(TOLA,TOLB), then the algorithm has converged.
*     Compute the generalized singular value pairs (ALPHA, BETA), and
*     set the triangular matrix R to array A.
*
      DO 60 I = 1, K
         ALPHA( I ) = ONE
         BETA( I ) = ZERO
   60 CONTINUE
*
      DO 70 I = 1, MIN( L, M-K )
*
         A1 = A( K+I, N-L+I )
         B1 = B( I, N-L+I )
*
         IF( A1.NE.ZERO ) THEN
            GAMMA = B1 / A1
*
*           change sign if necessary
*
            IF( GAMMA.LT.ZERO ) THEN
               CALL DSCAL( L-I+1, -ONE, B( I, N-L+I ), LDB )
               IF( WANTV )
     $            CALL DSCAL( P, -ONE, V( 1, I ), 1 )
            END IF
*
            CALL DLARTG( ABS( GAMMA ), ONE, BETA( K+I ), ALPHA( K+I ),
     $                   RWK )
*
            IF( ALPHA( K+I ).GE.BETA( K+I ) ) THEN
               CALL DSCAL( L-I+1, ONE / ALPHA( K+I ), A( K+I, N-L+I ),
     $                     LDA )
            ELSE
               CALL DSCAL( L-I+1, ONE / BETA( K+I ), B( I, N-L+I ),
     $                     LDB )
               CALL DCOPY( L-I+1, B( I, N-L+I ), LDB, A( K+I, N-L+I ),
     $                     LDA )
            END IF
*
         ELSE
*
            ALPHA( K+I ) = ZERO
            BETA( K+I ) = ONE
            CALL DCOPY( L-I+1, B( I, N-L+I ), LDB, A( K+I, N-L+I ),
     $                  LDA )
*
         END IF
*
   70 CONTINUE
*
*     Post-assignment
*
      DO 80 I = M + 1, K + L
         ALPHA( I ) = ZERO
         BETA( I ) = ONE
   80 CONTINUE
*
      IF( K+L.LT.N ) THEN
         DO 90 I = K + L + 1, N
            ALPHA( I ) = ZERO
            BETA( I ) = ZERO
   90    CONTINUE
      END IF
*
  100 CONTINUE
      NCYCLE = KCYCLE
      RETURN
*
*     End of DTGSJA
*
      END

      SUBROUTINE DGEQPF( M, N, A, LDA, JPVT, TAU, WORK, INFO )
*
*  -- LAPACK test routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  This routine is deprecated and has been replaced by routine DGEQP3.
*
*  DGEQPF computes a QR factorization with column pivoting of a
*  real M-by-N matrix A: A*P = Q*R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A. N >= 0
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the upper triangle of the array contains the
*          min(M,N)-by-N upper triangular matrix R; the elements
*          below the diagonal, together with the array TAU,
*          represent the orthogonal matrix Q as a product of
*          min(m,n) elementary reflectors.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  JPVT    (input/output) INTEGER array, dimension (N)
*          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
*          to the front of A*P (a leading column); if JPVT(i) = 0,
*          the i-th column of A is a free column.
*          On exit, if JPVT(i) = k, then the i-th column of A*P
*          was the k-th column of A.
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(n)
*
*  Each H(i) has the form
*
*     H = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i).
*
*  The matrix P is represented in jpvt as follows: If
*     jpvt(j) = i
*  then the jth column of P is the ith canonical unit vector.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ITEMP, J, MA, MN, PVT
      DOUBLE PRECISION   AII, TEMP, TEMP2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEQR2, DLARF, DLARFG, DORM2R, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DNRM2
      EXTERNAL           IDAMAX, DNRM2
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQPF', -INFO )
         RETURN
      END IF
*
      MN = MIN( M, N )
*
*     Move initial columns up front
*
      ITEMP = 1
      DO 10 I = 1, N
         IF( JPVT( I ).NE.0 ) THEN
            IF( I.NE.ITEMP ) THEN
               CALL DSWAP( M, A( 1, I ), 1, A( 1, ITEMP ), 1 )
               JPVT( I ) = JPVT( ITEMP )
               JPVT( ITEMP ) = I
            ELSE
               JPVT( I ) = I
            END IF
            ITEMP = ITEMP + 1
         ELSE
            JPVT( I ) = I
         END IF
   10 CONTINUE
      ITEMP = ITEMP - 1
*
*     Compute the QR factorization and update remaining columns
*
      IF( ITEMP.GT.0 ) THEN
         MA = MIN( ITEMP, M )
         CALL DGEQR2( M, MA, A, LDA, TAU, WORK, INFO )
         IF( MA.LT.N ) THEN
            CALL DORM2R( 'Left', 'Transpose', M, N-MA, MA, A, LDA, TAU,
     $                   A( 1, MA+1 ), LDA, WORK, INFO )
         END IF
      END IF
*
      IF( ITEMP.LT.MN ) THEN
*
*        Initialize partial column norms. The first n elements of
*        work store the exact column norms.
*
         DO 20 I = ITEMP + 1, N
            WORK( I ) = DNRM2( M-ITEMP, A( ITEMP+1, I ), 1 )
            WORK( N+I ) = WORK( I )
   20    CONTINUE
*
*        Compute factorization
*
         DO 40 I = ITEMP + 1, MN
*
*           Determine ith pivot column and swap if necessary
*
            PVT = ( I-1 ) + IDAMAX( N-I+1, WORK( I ), 1 )
*
            IF( PVT.NE.I ) THEN
               CALL DSWAP( M, A( 1, PVT ), 1, A( 1, I ), 1 )
               ITEMP = JPVT( PVT )
               JPVT( PVT ) = JPVT( I )
               JPVT( I ) = ITEMP
               WORK( PVT ) = WORK( I )
               WORK( N+PVT ) = WORK( N+I )
            END IF
*
*           Generate elementary reflector H(i)
*
            IF( I.LT.M ) THEN
               CALL DLARFG( M-I+1, A( I, I ), A( I+1, I ), 1, TAU( I ) )
            ELSE
               CALL DLARFG( 1, A( M, M ), A( M, M ), 1, TAU( M ) )
            END IF
*
            IF( I.LT.N ) THEN
*
*              Apply H(i) to A(i:m,i+1:n) from the left
*
               AII = A( I, I )
               A( I, I ) = ONE
               CALL DLARF( 'LEFT', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                     A( I, I+1 ), LDA, WORK( 2*N+1 ) )
               A( I, I ) = AII
            END IF
*
*           Update partial column norms
*
            DO 30 J = I + 1, N
               IF( WORK( J ).NE.ZERO ) THEN
                  TEMP = ONE - ( ABS( A( I, J ) ) / WORK( J ) )**2
                  TEMP = MAX( TEMP, ZERO )
                  TEMP2 = ONE + 0.05D0*TEMP*
     $                    ( WORK( J ) / WORK( N+J ) )**2
                  IF( TEMP2.EQ.ONE ) THEN
                     IF( M-I.GT.0 ) THEN
                        WORK( J ) = DNRM2( M-I, A( I+1, J ), 1 )
                        WORK( N+J ) = WORK( J )
                     ELSE
                        WORK( J ) = ZERO
                        WORK( N+J ) = ZERO
                     END IF
                  ELSE
                     WORK( J ) = WORK( J )*SQRT( TEMP )
                  END IF
               END IF
   30       CONTINUE
*
   40    CONTINUE
      END IF
      RETURN
*
*     End of DGEQPF
*
      END

      SUBROUTINE DLAPLL( N, X, INCX, Y, INCY, SSMIN )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      DOUBLE PRECISION   SSMIN
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  Given two column vectors X and Y, let
*
*                       A = ( X Y ).
*
*  The subroutine first computes the QR factorization of A = Q*R,
*  and then computes the SVD of the 2-by-2 upper triangular matrix R.
*  The smaller singular value of R is returned in SSMIN, which is used
*  as the measurement of the linear dependency of the vectors X and Y.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The length of the vectors X and Y.
*
*  X       (input/output) DOUBLE PRECISION array,
*                         dimension (1+(N-1)*INCX)
*          On entry, X contains the N-vector X.
*          On exit, X is overwritten.
*
*  INCX    (input) INTEGER
*          The increment between successive elements of X. INCX > 0.
*
*  Y       (input/output) DOUBLE PRECISION array,
*                         dimension (1+(N-1)*INCY)
*          On entry, Y contains the N-vector Y.
*          On exit, Y is overwritten.
*
*  INCY    (input) INTEGER
*          The increment between successive elements of Y. INCY > 0.
*
*  SSMIN   (output) DOUBLE PRECISION
*          The smallest singular value of the N-by-2 matrix A = ( X Y ).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   A11, A12, A22, C, SSMAX, TAU
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DLARFG, DLAS2
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.1 ) THEN
         SSMIN = ZERO
         RETURN
      END IF
*
*     Compute the QR factorization of the N-by-2 matrix ( X Y )
*
      CALL DLARFG( N, X( 1 ), X( 1+INCX ), INCX, TAU )
      A11 = X( 1 )
      X( 1 ) = ONE
*
      C = -TAU*DDOT( N, X, INCX, Y, INCY )
      CALL DAXPY( N, C, X, INCX, Y, INCY )
*
      CALL DLARFG( N-1, Y( 1+INCY ), Y( 1+2*INCY ), INCY, TAU )
*
      A12 = Y( 1 )
      A22 = Y( 1+INCY )
*
*     Compute the SVD of 2-by-2 Upper triangular matrix.
*
      CALL DLAS2( A11, A12, A22, SSMIN, SSMAX )
*
      RETURN
*
*     End of DLAPLL
*
      END

      SUBROUTINE DLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV,
     $                   SNV, CSQ, SNQ )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      LOGICAL            UPPER
      DOUBLE PRECISION   A1, A2, A3, B1, B2, B3, CSQ, CSU, CSV, SNQ,
     $                   SNU, SNV
*     ..
*
*  Purpose
*  =======
*
*  DLAGS2 computes 2-by-2 orthogonal matrices U, V and Q, such
*  that if ( UPPER ) then
*
*            U'*A*Q = U'*( A1 A2 )*Q = ( x  0  )
*                        ( 0  A3 )     ( x  x  )
*  and
*            V'*B*Q = V'*( B1 B2 )*Q = ( x  0  )
*                        ( 0  B3 )     ( x  x  )
*
*  or if ( .NOT.UPPER ) then
*
*            U'*A*Q = U'*( A1 0  )*Q = ( x  x  )
*                        ( A2 A3 )     ( 0  x  )
*  and
*            V'*B*Q = V'*( B1 0  )*Q = ( x  x  )
*                        ( B2 B3 )     ( 0  x  )
*
*  The rows of the transformed A and B are parallel, where
*
*    U = (  CSU  SNU ), V = (  CSV SNV ), Q = (  CSQ   SNQ )
*        ( -SNU  CSU )      ( -SNV CSV )      ( -SNQ   CSQ )
*
*  Z' denotes the transpose of Z.
*
*
*  Arguments
*  =========
*
*  UPPER   (input) LOGICAL
*          = .TRUE.: the input matrices A and B are upper triangular.
*          = .FALSE.: the input matrices A and B are lower triangular.
*
*  A1      (input) DOUBLE PRECISION
*  A2      (input) DOUBLE PRECISION
*  A3      (input) DOUBLE PRECISION
*          On entry, A1, A2 and A3 are elements of the input 2-by-2
*          upper (lower) triangular matrix A.
*
*  B1      (input) DOUBLE PRECISION
*  B2      (input) DOUBLE PRECISION
*  B3      (input) DOUBLE PRECISION
*          On entry, B1, B2 and B3 are elements of the input 2-by-2
*          upper (lower) triangular matrix B.
*
*  CSU     (output) DOUBLE PRECISION
*  SNU     (output) DOUBLE PRECISION
*          The desired orthogonal matrix U.
*
*  CSV     (output) DOUBLE PRECISION
*  SNV     (output) DOUBLE PRECISION
*          The desired orthogonal matrix V.
*
*  CSQ     (output) DOUBLE PRECISION
*  SNQ     (output) DOUBLE PRECISION
*          The desired orthogonal matrix Q.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   A, AUA11, AUA12, AUA21, AUA22, AVB11, AVB12,
     $                   AVB21, AVB22, B, C, CSL, CSR, D, R, S1, S2,
     $                   SNL, SNR, UA11, UA11R, UA12, UA21, UA22, UA22R,
     $                   VB11, VB11R, VB12, VB21, VB22, VB22R
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARTG, DLASV2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( UPPER ) THEN
*
*        Input matrices A and B are upper triangular matrices
*
*        Form matrix C = A*adj(B) = ( a b )
*                                   ( 0 d )
*
         A = A1*B3
         D = A3*B1
         B = A2*B1 - A1*B2
*
*        The SVD of real 2-by-2 triangular C
*
*         ( CSL -SNL )*( A B )*(  CSR  SNR ) = ( R 0 )
*         ( SNL  CSL ) ( 0 D ) ( -SNR  CSR )   ( 0 T )
*
         CALL DLASV2( A, B, D, S1, S2, SNR, CSR, SNL, CSL )
*
         IF( ABS( CSL ).GE.ABS( SNL ) .OR. ABS( CSR ).GE.ABS( SNR ) )
     $        THEN
*
*           Compute the (1,1) and (1,2) elements of U'*A and V'*B,
*           and (1,2) element of |U|'*|A| and |V|'*|B|.
*
            UA11R = CSL*A1
            UA12 = CSL*A2 + SNL*A3
*
            VB11R = CSR*B1
            VB12 = CSR*B2 + SNR*B3
*
            AUA12 = ABS( CSL )*ABS( A2 ) + ABS( SNL )*ABS( A3 )
            AVB12 = ABS( CSR )*ABS( B2 ) + ABS( SNR )*ABS( B3 )
*
*           zero (1,2) elements of U'*A and V'*B
*
            IF( ( ABS( UA11R )+ABS( UA12 ) ).NE.ZERO ) THEN
               IF( AUA12 / ( ABS( UA11R )+ABS( UA12 ) ).LE.AVB12 /
     $             ( ABS( VB11R )+ABS( VB12 ) ) ) THEN
                  CALL DLARTG( -UA11R, UA12, CSQ, SNQ, R )
               ELSE
                  CALL DLARTG( -VB11R, VB12, CSQ, SNQ, R )
               END IF
            ELSE
               CALL DLARTG( -VB11R, VB12, CSQ, SNQ, R )
            END IF
*
            CSU = CSL
            SNU = -SNL
            CSV = CSR
            SNV = -SNR
*
         ELSE
*
*           Compute the (2,1) and (2,2) elements of U'*A and V'*B,
*           and (2,2) element of |U|'*|A| and |V|'*|B|.
*
            UA21 = -SNL*A1
            UA22 = -SNL*A2 + CSL*A3
*
            VB21 = -SNR*B1
            VB22 = -SNR*B2 + CSR*B3
*
            AUA22 = ABS( SNL )*ABS( A2 ) + ABS( CSL )*ABS( A3 )
            AVB22 = ABS( SNR )*ABS( B2 ) + ABS( CSR )*ABS( B3 )
*
*           zero (2,2) elements of U'*A and V'*B, and then swap.
*
            IF( ( ABS( UA21 )+ABS( UA22 ) ).NE.ZERO ) THEN
               IF( AUA22 / ( ABS( UA21 )+ABS( UA22 ) ).LE.AVB22 /
     $             ( ABS( VB21 )+ABS( VB22 ) ) ) THEN
                  CALL DLARTG( -UA21, UA22, CSQ, SNQ, R )
               ELSE
                  CALL DLARTG( -VB21, VB22, CSQ, SNQ, R )
               END IF
            ELSE
               CALL DLARTG( -VB21, VB22, CSQ, SNQ, R )
            END IF
*
            CSU = SNL
            SNU = CSL
            CSV = SNR
            SNV = CSR
*
         END IF
*
      ELSE
*
*        Input matrices A and B are lower triangular matrices
*
*        Form matrix C = A*adj(B) = ( a 0 )
*                                   ( c d )
*
         A = A1*B3
         D = A3*B1
         C = A2*B3 - A3*B2
*
*        The SVD of real 2-by-2 triangular C
*
*         ( CSL -SNL )*( A 0 )*(  CSR  SNR ) = ( R 0 )
*         ( SNL  CSL ) ( C D ) ( -SNR  CSR )   ( 0 T )
*
         CALL DLASV2( A, C, D, S1, S2, SNR, CSR, SNL, CSL )
*
         IF( ABS( CSR ).GE.ABS( SNR ) .OR. ABS( CSL ).GE.ABS( SNL ) )
     $        THEN
*
*           Compute the (2,1) and (2,2) elements of U'*A and V'*B,
*           and (2,1) element of |U|'*|A| and |V|'*|B|.
*
            UA21 = -SNR*A1 + CSR*A2
            UA22R = CSR*A3
*
            VB21 = -SNL*B1 + CSL*B2
            VB22R = CSL*B3
*
            AUA21 = ABS( SNR )*ABS( A1 ) + ABS( CSR )*ABS( A2 )
            AVB21 = ABS( SNL )*ABS( B1 ) + ABS( CSL )*ABS( B2 )
*
*           zero (2,1) elements of U'*A and V'*B.
*
            IF( ( ABS( UA21 )+ABS( UA22R ) ).NE.ZERO ) THEN
               IF( AUA21 / ( ABS( UA21 )+ABS( UA22R ) ).LE.AVB21 /
     $             ( ABS( VB21 )+ABS( VB22R ) ) ) THEN
                  CALL DLARTG( UA22R, UA21, CSQ, SNQ, R )
               ELSE
                  CALL DLARTG( VB22R, VB21, CSQ, SNQ, R )
               END IF
            ELSE
               CALL DLARTG( VB22R, VB21, CSQ, SNQ, R )
            END IF
*
            CSU = CSR
            SNU = -SNR
            CSV = CSL
            SNV = -SNL
*
         ELSE
*
*           Compute the (1,1) and (1,2) elements of U'*A and V'*B,
*           and (1,1) element of |U|'*|A| and |V|'*|B|.
*
            UA11 = CSR*A1 + SNR*A2
            UA12 = SNR*A3
*
            VB11 = CSL*B1 + SNL*B2
            VB12 = SNL*B3
*
            AUA11 = ABS( CSR )*ABS( A1 ) + ABS( SNR )*ABS( A2 )
            AVB11 = ABS( CSL )*ABS( B1 ) + ABS( SNL )*ABS( B2 )
*
*           zero (1,1) elements of U'*A and V'*B, and then swap.
*
            IF( ( ABS( UA11 )+ABS( UA12 ) ).NE.ZERO ) THEN
               IF( AUA11 / ( ABS( UA11 )+ABS( UA12 ) ).LE.AVB11 /
     $             ( ABS( VB11 )+ABS( VB12 ) ) ) THEN
                  CALL DLARTG( UA12, UA11, CSQ, SNQ, R )
               ELSE
                  CALL DLARTG( VB12, VB11, CSQ, SNQ, R )
               END IF
            ELSE
               CALL DLARTG( VB12, VB11, CSQ, SNQ, R )
            END IF
*
            CSU = SNR
            SNU = CSR
            CSV = SNL
            SNV = CSL
*
         END IF
*
      END IF
*
      RETURN
*
*     End of DLAGS2
*
      END

      SUBROUTINE DLAS2( F, G, H, SSMIN, SSMAX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   F, G, H, SSMAX, SSMIN
*     ..
*
*  Purpose
*  =======
*
*  DLAS2  computes the singular values of the 2-by-2 matrix
*     [  F   G  ]
*     [  0   H  ].
*  On return, SSMIN is the smaller singular value and SSMAX is the
*  larger singular value.
*
*  Arguments
*  =========
*
*  F       (input) DOUBLE PRECISION
*          The (1,1) element of the 2-by-2 matrix.
*
*  G       (input) DOUBLE PRECISION
*          The (1,2) element of the 2-by-2 matrix.
*
*  H       (input) DOUBLE PRECISION
*          The (2,2) element of the 2-by-2 matrix.
*
*  SSMIN   (output) DOUBLE PRECISION
*          The smaller singular value.
*
*  SSMAX   (output) DOUBLE PRECISION
*          The larger singular value.
*
*  Further Details
*  ===============
*
*  Barring over/underflow, all output quantities are correct to within
*  a few units in the last place (ulps), even in the absence of a guard
*  digit in addition/subtraction.
*
*  In IEEE arithmetic, the code works correctly if one matrix element is
*  infinite.
*
*  Overflow will not occur unless the largest singular value itself
*  overflows, or is within a few ulps of overflow. (On machines with
*  partial overflow, like the Cray, overflow may occur if the largest
*  singular value is within a factor of 2 of overflow.)
*
*  Underflow is harmless if underflow is gradual. Otherwise, results
*  may correspond to a matrix modified by perturbations of size near
*  the underflow threshold.
*
*  ====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AS, AT, AU, C, FA, FHMN, FHMX, GA, HA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      FA = ABS( F )
      GA = ABS( G )
      HA = ABS( H )
      FHMN = MIN( FA, HA )
      FHMX = MAX( FA, HA )
      IF( FHMN.EQ.ZERO ) THEN
         SSMIN = ZERO
         IF( FHMX.EQ.ZERO ) THEN
            SSMAX = GA
         ELSE
            SSMAX = MAX( FHMX, GA )*SQRT( ONE+
     $              ( MIN( FHMX, GA ) / MAX( FHMX, GA ) )**2 )
         END IF
      ELSE
         IF( GA.LT.FHMX ) THEN
            AS = ONE + FHMN / FHMX
            AT = ( FHMX-FHMN ) / FHMX
            AU = ( GA / FHMX )**2
            C = TWO / ( SQRT( AS*AS+AU )+SQRT( AT*AT+AU ) )
            SSMIN = FHMN*C
            SSMAX = FHMX / C
         ELSE
            AU = FHMX / GA
            IF( AU.EQ.ZERO ) THEN
*
*              Avoid possible harmful underflow if exponent range
*              asymmetric (true SSMIN may not underflow even if
*              AU underflows)
*
               SSMIN = ( FHMN*FHMX ) / GA
               SSMAX = GA
            ELSE
               AS = ONE + FHMN / FHMX
               AT = ( FHMX-FHMN ) / FHMX
               C = ONE / ( SQRT( ONE+( AS*AU )**2 )+
     $             SQRT( ONE+( AT*AU )**2 ) )
               SSMIN = ( FHMN*C )*AU
               SSMIN = SSMIN + SSMIN
               SSMAX = GA / ( C+C )
            END IF
         END IF
      END IF
      RETURN
*
*     End of DLAS2
*
      END

      SUBROUTINE DLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   CSL, CSR, F, G, H, SNL, SNR, SSMAX, SSMIN
*     ..
*
*  Purpose
*  =======
*
*  DLASV2 computes the singular value decomposition of a 2-by-2
*  triangular matrix
*     [  F   G  ]
*     [  0   H  ].
*  On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the
*  smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and
*  right singular vectors for abs(SSMAX), giving the decomposition
*
*     [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]
*     [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].
*
*  Arguments
*  =========
*
*  F       (input) DOUBLE PRECISION
*          The (1,1) element of the 2-by-2 matrix.
*
*  G       (input) DOUBLE PRECISION
*          The (1,2) element of the 2-by-2 matrix.
*
*  H       (input) DOUBLE PRECISION
*          The (2,2) element of the 2-by-2 matrix.
*
*  SSMIN   (output) DOUBLE PRECISION
*          abs(SSMIN) is the smaller singular value.
*
*  SSMAX   (output) DOUBLE PRECISION
*          abs(SSMAX) is the larger singular value.
*
*  SNL     (output) DOUBLE PRECISION
*  CSL     (output) DOUBLE PRECISION
*          The vector (CSL, SNL) is a unit left singular vector for the
*          singular value abs(SSMAX).
*
*  SNR     (output) DOUBLE PRECISION
*  CSR     (output) DOUBLE PRECISION
*          The vector (CSR, SNR) is a unit right singular vector for the
*          singular value abs(SSMAX).
*
*  Further Details
*  ===============
*
*  Any input parameter may be aliased with any output parameter.
*
*  Barring over/underflow and assuming a guard digit in subtraction, all
*  output quantities are correct to within a few units in the last
*  place (ulps).
*
*  In IEEE arithmetic, the code works correctly if one matrix element is
*  infinite.
*
*  Overflow will not occur unless the largest singular value itself
*  overflows or is within a few ulps of overflow. (On machines with
*  partial overflow, like the Cray, overflow may occur if the largest
*  singular value is within a factor of 2 of overflow.)
*
*  Underflow is harmless if underflow is gradual. Otherwise, results
*  may correspond to a matrix modified by perturbations of size near
*  the underflow threshold.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   FOUR
      PARAMETER          ( FOUR = 4.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            GASMAL, SWAP
      INTEGER            PMAX
      DOUBLE PRECISION   A, CLT, CRT, D, FA, FT, GA, GT, HA, HT, L, M,
     $                   MM, R, S, SLT, SRT, T, TEMP, TSIGN, TT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Executable Statements ..
*
      FT = F
      FA = ABS( FT )
      HT = H
      HA = ABS( H )
*
*     PMAX points to the maximum absolute element of matrix
*       PMAX = 1 if F largest in absolute values
*       PMAX = 2 if G largest in absolute values
*       PMAX = 3 if H largest in absolute values
*
      PMAX = 1
      SWAP = ( HA.GT.FA )
      IF( SWAP ) THEN
         PMAX = 3
         TEMP = FT
         FT = HT
         HT = TEMP
         TEMP = FA
         FA = HA
         HA = TEMP
*
*        Now FA .ge. HA
*
      END IF
      GT = G
      GA = ABS( GT )
      IF( GA.EQ.ZERO ) THEN
*
*        Diagonal matrix
*
         SSMIN = HA
         SSMAX = FA
         CLT = ONE
         CRT = ONE
         SLT = ZERO
         SRT = ZERO
      ELSE
         GASMAL = .TRUE.
         IF( GA.GT.FA ) THEN
            PMAX = 2
            IF( ( FA / GA ).LT.DLAMCH( 'EPS' ) ) THEN
*
*              Case of very large GA
*
               GASMAL = .FALSE.
               SSMAX = GA
               IF( HA.GT.ONE ) THEN
                  SSMIN = FA / ( GA / HA )
               ELSE
                  SSMIN = ( FA / GA )*HA
               END IF
               CLT = ONE
               SLT = HT / GT
               SRT = ONE
               CRT = FT / GT
            END IF
         END IF
         IF( GASMAL ) THEN
*
*           Normal case
*
            D = FA - HA
            IF( D.EQ.FA ) THEN
*
*              Copes with infinite F or H
*
               L = ONE
            ELSE
               L = D / FA
            END IF
*
*           Note that 0 .le. L .le. 1
*
            M = GT / FT
*
*           Note that abs(M) .le. 1/macheps
*
            T = TWO - L
*
*           Note that T .ge. 1
*
            MM = M*M
            TT = T*T
            S = SQRT( TT+MM )
*
*           Note that 1 .le. S .le. 1 + 1/macheps
*
            IF( L.EQ.ZERO ) THEN
               R = ABS( M )
            ELSE
               R = SQRT( L*L+MM )
            END IF
*
*           Note that 0 .le. R .le. 1 + 1/macheps
*
            A = HALF*( S+R )
*
*           Note that 1 .le. A .le. 1 + abs(M)
*
            SSMIN = HA / A
            SSMAX = FA*A
            IF( MM.EQ.ZERO ) THEN
*
*              Note that M is very tiny
*
               IF( L.EQ.ZERO ) THEN
                  T = SIGN( TWO, FT )*SIGN( ONE, GT )
               ELSE
                  T = GT / SIGN( D, FT ) + M / T
               END IF
            ELSE
               T = ( M / ( S+T )+M / ( R+L ) )*( ONE+A )
            END IF
            L = SQRT( T*T+FOUR )
            CRT = TWO / L
            SRT = T / L
            CLT = ( CRT+SRT*M ) / A
            SLT = ( HT / FT )*SRT / A
         END IF
      END IF
      IF( SWAP ) THEN
         CSL = SRT
         SNL = CRT
         CSR = SLT
         SNR = CLT
      ELSE
         CSL = CLT
         SNL = SLT
         CSR = CRT
         SNR = SRT
      END IF
*
*     Correct signs of SSMAX and SSMIN
*
      IF( PMAX.EQ.1 )
     $   TSIGN = SIGN( ONE, CSR )*SIGN( ONE, CSL )*SIGN( ONE, F )
      IF( PMAX.EQ.2 )
     $   TSIGN = SIGN( ONE, SNR )*SIGN( ONE, CSL )*SIGN( ONE, G )
      IF( PMAX.EQ.3 )
     $   TSIGN = SIGN( ONE, SNR )*SIGN( ONE, SNL )*SIGN( ONE, H )
      SSMAX = SIGN( SSMAX, TSIGN )
      SSMIN = SIGN( SSMIN, TSIGN*SIGN( ONE, F )*SIGN( ONE, H ) )
      RETURN
*
*     End of DLASV2
*
      END


      SUBROUTINE DLAPMT( FORWRD, M, N, X, LDX, K )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      LOGICAL            FORWRD
      INTEGER            LDX, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            K( * )
      DOUBLE PRECISION   X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  DLAPMT rearranges the columns of the M by N matrix X as specified
*  by the permutation K(1),K(2),...,K(N) of the integers 1,...,N.
*  If FORWRD = .TRUE.,  forward permutation:
*
*       X(*,K(J)) is moved X(*,J) for J = 1,2,...,N.
*
*  If FORWRD = .FALSE., backward permutation:
*
*       X(*,J) is moved to X(*,K(J)) for J = 1,2,...,N.
*
*  Arguments
*  =========
*
*  FORWRD  (input) LOGICAL
*          = .TRUE., forward permutation
*          = .FALSE., backward permutation
*
*  M       (input) INTEGER
*          The number of rows of the matrix X. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix X. N >= 0.
*
*  X       (input/output) DOUBLE PRECISION array, dimension (LDX,N)
*          On entry, the M by N matrix X.
*          On exit, X contains the permuted matrix X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X, LDX >= MAX(1,M).
*
*  K       (input) INTEGER array, dimension (N)
*          On entry, K contains the permutation vector.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, II, IN, J
      DOUBLE PRECISION   TEMP
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.1 )
     $   RETURN
*
      DO 10 I = 1, N
         K( I ) = -K( I )
   10 CONTINUE
*
      IF( FORWRD ) THEN
*
*        Forward permutation
*
         DO 50 I = 1, N
*
            IF( K( I ).GT.0 )
     $         GO TO 40
*
            J = I
            K( J ) = -K( J )
            IN = K( J )
*
   20       CONTINUE
            IF( K( IN ).GT.0 )
     $         GO TO 40
*
            DO 30 II = 1, M
               TEMP = X( II, J )
               X( II, J ) = X( II, IN )
               X( II, IN ) = TEMP
   30       CONTINUE
*
            K( IN ) = -K( IN )
            J = IN
            IN = K( IN )
            GO TO 20
*
   40       CONTINUE
*
   50    CONTINUE
*
      ELSE
*
*        Backward permutation
*
         DO 90 I = 1, N
*
            IF( K( I ).GT.0 )
     $         GO TO 80
*
            K( I ) = -K( I )
            J = K( I )
   60       CONTINUE
            IF( J.EQ.I )
     $         GO TO 80
*
            DO 70 II = 1, M
               TEMP = X( II, I )
               X( II, I ) = X( II, J )
               X( II, J ) = TEMP
   70       CONTINUE
*
            K( J ) = -K( J )
            J = K( J )
            GO TO 60
*
   80       CONTINUE
*
   90    CONTINUE
*
      END IF
*
      RETURN
*
*     End of DLAPMT
*
      END

      SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEQR2 computes a QR factorization of a real m by n matrix A:
*  A = Q * R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, the elements on and above the diagonal of the array
*          contain the min(m,n) by n upper trapezoidal matrix R (R is
*          upper triangular if m >= n); the elements below the diagonal,
*          with the array TAU, represent the orthogonal matrix Q as a
*          product of elementary reflectors (see Further Details).
*
*  LDA     (input) INTEGER
*        The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*  and tau in TAU(i).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K
      DOUBLE PRECISION   AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQR2', -INFO )
         RETURN
      END IF
*
      K = MIN( M, N )
*
      DO 10 I = 1, K
*
*        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
*
         CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,
     $                TAU( I ) )
         IF( I.LT.N ) THEN
*
*           Apply H(i) to A(i:m,i+1:n) from the left
*
            AII = A( I, I )
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
            A( I, I ) = AII
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGEQR2
*
      END
      SUBROUTINE DGERQ2( M, N, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGERQ2 computes an RQ factorization of a real m by n matrix A:
*  A = R * Q.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, if m <= n, the upper triangle of the subarray
*          A(1:m,n-m+1:n) contains the m by m upper triangular matrix R;
*          if m >= n, the elements on and above the (m-n)-th subdiagonal
*          contain the m by n upper trapezoidal matrix R; the remaining
*          elements, with the array TAU, represent the orthogonal matrix
*          Q as a product of elementary reflectors (see Further
*          Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (M)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in
*  A(m-k+i,1:n-k+i-1), and tau in TAU(i).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K
      DOUBLE PRECISION   AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGERQ2', -INFO )
         RETURN
      END IF
*
      K = MIN( M, N )
*
      DO 10 I = K, 1, -1
*
*        Generate elementary reflector H(i) to annihilate
*        A(m-k+i,1:n-k+i-1)
*
         CALL DLARFG( N-K+I, A( M-K+I, N-K+I ), A( M-K+I, 1 ), LDA,
     $                TAU( I ) )
*
*        Apply H(i) to A(1:m-k+i-1,1:n-k+i) from the right
*
         AII = A( M-K+I, N-K+I )
         A( M-K+I, N-K+I ) = ONE
         CALL DLARF( 'Right', M-K+I-1, N-K+I, A( M-K+I, 1 ), LDA,
     $               TAU( I ), A, LDA, WORK )
         A( M-K+I, N-K+I ) = AII
   10 CONTINUE
      RETURN
*
*     End of DGERQ2
*
      END


      SUBROUTINE DGGSVD( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B,
     $                   LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK,
     $                   IWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBQ, JOBU, JOBV
      INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), Q( LDQ, * ), U( LDU, * ),
     $                   V( LDV, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGGSVD computes the generalized singular value decomposition (GSVD)
*  of an M-by-N real matrix A and P-by-N real matrix B:
*
*      U'*A*Q = D1*( 0 R ),    V'*B*Q = D2*( 0 R )
*
*  where U, V and Q are orthogonal matrices, and Z' is the transpose
*  of Z.  Let K+L = the effective numerical rank of the matrix (A',B')',
*  then R is a K+L-by-K+L nonsingular upper triangular matrix, D1 and
*  D2 are M-by-(K+L) and P-by-(K+L) "diagonal" matrices and of the
*  following structures, respectively:
*
*  If M-K-L >= 0,
*
*                      K  L
*         D1 =     K ( I  0 )
*                  L ( 0  C )
*              M-K-L ( 0  0 )
*
*                    K  L
*         D2 =   L ( 0  S )
*              P-L ( 0  0 )
*
*                  N-K-L  K    L
*    ( 0 R ) = K (  0   R11  R12 )
*              L (  0    0   R22 )
*
*  where
*
*    C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
*    S = diag( BETA(K+1),  ... , BETA(K+L) ),
*    C**2 + S**2 = I.
*
*    R is stored in A(1:K+L,N-K-L+1:N) on exit.
*
*  If M-K-L < 0,
*
*                    K M-K K+L-M
*         D1 =   K ( I  0    0   )
*              M-K ( 0  C    0   )
*
*                      K M-K K+L-M
*         D2 =   M-K ( 0  S    0  )
*              K+L-M ( 0  0    I  )
*                P-L ( 0  0    0  )
*
*                     N-K-L  K   M-K  K+L-M
*    ( 0 R ) =     K ( 0    R11  R12  R13  )
*                M-K ( 0     0   R22  R23  )
*              K+L-M ( 0     0    0   R33  )
*
*  where
*
*    C = diag( ALPHA(K+1), ... , ALPHA(M) ),
*    S = diag( BETA(K+1),  ... , BETA(M) ),
*    C**2 + S**2 = I.
*
*    (R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N), and R33 is stored
*    ( 0  R22 R23 )
*    in B(M-K+1:L,N+M-K-L+1:N) on exit.
*
*  The routine computes C, S, R, and optionally the orthogonal
*  transformation matrices U, V and Q.
*
*  In particular, if B is an N-by-N nonsingular matrix, then the GSVD of
*  A and B implicitly gives the SVD of A*inv(B):
*                       A*inv(B) = U*(D1*inv(D2))*V'.
*  If ( A',B')' has orthonormal columns, then the GSVD of A and B is
*  also equal to the CS decomposition of A and B. Furthermore, the GSVD
*  can be used to derive the solution of the eigenvalue problem:
*                       A'*A x = lambda* B'*B x.
*  In some literature, the GSVD of A and B is presented in the form
*                   U'*A*X = ( 0 D1 ),   V'*B*X = ( 0 D2 )
*  where U and V are orthogonal and X is nonsingular, D1 and D2 are
*  ``diagonal''.  The former GSVD form can be converted to the latter
*  form by taking the nonsingular matrix X as
*
*                       X = Q*( I   0    )
*                             ( 0 inv(R) ).
*
*  Arguments
*  =========
*
*  JOBU    (input) CHARACTER*1
*          = 'U':  Orthogonal matrix U is computed;
*          = 'N':  U is not computed.
*
*  JOBV    (input) CHARACTER*1
*          = 'V':  Orthogonal matrix V is computed;
*          = 'N':  V is not computed.
*
*  JOBQ    (input) CHARACTER*1
*          = 'Q':  Orthogonal matrix Q is computed;
*          = 'N':  Q is not computed.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and B.  N >= 0.
*
*  P       (input) INTEGER
*          The number of rows of the matrix B.  P >= 0.
*
*  K       (output) INTEGER
*  L       (output) INTEGER
*          On exit, K and L specify the dimension of the subblocks
*          described in the Purpose section.
*          K + L = effective numerical rank of (A',B')'.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, A contains the triangular matrix R, or part of R.
*          See Purpose for details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
*          On entry, the P-by-N matrix B.
*          On exit, B contains the triangular matrix R if M-K-L < 0.
*          See Purpose for details.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDA >= max(1,P).
*
*  ALPHA   (output) DOUBLE PRECISION array, dimension (N)
*  BETA    (output) DOUBLE PRECISION array, dimension (N)
*          On exit, ALPHA and BETA contain the generalized singular
*          value pairs of A and B;
*            ALPHA(1:K) = 1,
*            BETA(1:K)  = 0,
*          and if M-K-L >= 0,
*            ALPHA(K+1:K+L) = C,
*            BETA(K+1:K+L)  = S,
*          or if M-K-L < 0,
*            ALPHA(K+1:M)=C, ALPHA(M+1:K+L)=0
*            BETA(K+1:M) =S, BETA(M+1:K+L) =1
*          and
*            ALPHA(K+L+1:N) = 0
*            BETA(K+L+1:N)  = 0
*
*  U       (output) DOUBLE PRECISION array, dimension (LDU,M)
*          If JOBU = 'U', U contains the M-by-M orthogonal matrix U.
*          If JOBU = 'N', U is not referenced.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U. LDU >= max(1,M) if
*          JOBU = 'U'; LDU >= 1 otherwise.
*
*  V       (output) DOUBLE PRECISION array, dimension (LDV,P)
*          If JOBV = 'V', V contains the P-by-P orthogonal matrix V.
*          If JOBV = 'N', V is not referenced.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V. LDV >= max(1,P) if
*          JOBV = 'V'; LDV >= 1 otherwise.
*
*  Q       (output) DOUBLE PRECISION array, dimension (LDQ,N)
*          If JOBQ = 'Q', Q contains the N-by-N orthogonal matrix Q.
*          If JOBQ = 'N', Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q. LDQ >= max(1,N) if
*          JOBQ = 'Q'; LDQ >= 1 otherwise.
*
*  WORK    (workspace) DOUBLE PRECISION array,
*                      dimension (max(3*N,M,P)+N)
*
*  IWORK   (workspace/output) INTEGER array, dimension (N)
*          On exit, IWORK stores the sorting information. More
*          precisely, the following loop will sort ALPHA
*             for I = K+1, min(M,K+L)
*                 swap ALPHA(I) and ALPHA(IWORK(I))
*             endfor
*          such that ALPHA(1) <= ALPHA(2) <= ... <= ALPHA(N).
*
*  INFO    (output)INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = 1, the Jacobi-type procedure failed to
*                converge.  For further details, see subroutine DTGSJA.
*
*  Internal Parameters
*  ===================
*
*  TOLA    DOUBLE PRECISION
*  TOLB    DOUBLE PRECISION
*          TOLA and TOLB are the thresholds to determine the effective
*          rank of (A',B')'. Generally, they are set to
*                   TOLA = MAX(M,N)*norm(A)*MAZHEPS,
*                   TOLB = MAX(P,N)*norm(B)*MAZHEPS.
*          The size of TOLA and TOLB may affect the size of backward
*          errors of the decomposition.
*
*  Further Details
*  ===============
*
*  2-96 Based on modifications by
*     Ming Gu and Huan Ren, Computer Science Division, University of
*     California at Berkeley, USA
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            WANTQ, WANTU, WANTV
      INTEGER            I, IBND, ISUB, J, NCYCLE
      DOUBLE PRECISION   ANORM, BNORM, SMIN, TEMP, TOLA, TOLB, ULP, UNFL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGGSVP, DTGSJA, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      WANTU = LSAME( JOBU, 'U' )
      WANTV = LSAME( JOBV, 'V' )
      WANTQ = LSAME( JOBQ, 'Q' )
*
      INFO = 0
      IF( .NOT.( WANTU .OR. LSAME( JOBU, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTV .OR. LSAME( JOBV, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( WANTQ .OR. LSAME( JOBQ, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( P.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDU.LT.1 .OR. ( WANTU .AND. LDU.LT.M ) ) THEN
         INFO = -16
      ELSE IF( LDV.LT.1 .OR. ( WANTV .AND. LDV.LT.P ) ) THEN
         INFO = -18
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -20
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGSVD', -INFO )
         RETURN
      END IF
*
*     Compute the Frobenius norm of matrices A and B
*
      ANORM = DLANGE( '1', M, N, A, LDA, WORK )
      BNORM = DLANGE( '1', P, N, B, LDB, WORK )
*
*     Get machine precision and set up threshold for determining
*     the effective numerical rank of the matrices A and B.
*
      ULP = DLAMCH( 'Precision' )
      UNFL = DLAMCH( 'Safe Minimum' )
      TOLA = MAX( M, N )*MAX( ANORM, UNFL )*ULP
      TOLB = MAX( P, N )*MAX( BNORM, UNFL )*ULP
*
*     Preprocessing
*
      CALL DGGSVP( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA,
     $             TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, WORK,
     $             WORK( N+1 ), INFO )

*
*     Compute the GSVD of two upper "triangular" matrices
*
      CALL DTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, LDB,
     $             TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ,
     $             WORK, NCYCLE, INFO )
*
*     Sort the singular values and store the pivot indices in IWORK
*     Copy ALPHA to WORK, then sort ALPHA in WORK
*
      CALL DCOPY( N, ALPHA, 1, WORK, 1 )
      IBND = K + L
      DO 20 I = 1, IBND
*
*        Scan for smallest ALPHA(I)  !xliu
*
         ISUB = I
         SMIN = WORK( I )
         DO 10 J = I + 1, IBND
            TEMP = WORK( J )
            IF( TEMP.LT.SMIN ) THEN
               ISUB = J
               SMIN = TEMP
            END IF
   10    CONTINUE
         IF( ISUB.NE.I ) THEN
            WORK( ISUB ) = WORK( I )
            WORK( I ) = SMIN
            IWORK( I ) = ISUB
         ELSE
            IWORK( I ) =  I
         END IF
   20 CONTINUE
*
      RETURN
*
*     End of DGGSVD
*
      END

      SUBROUTINE DGGSVP( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB,
     $                   TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ,
     $                   IWORK, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          JOBQ, JOBU, JOBV
      INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P
      DOUBLE PRECISION   TOLA, TOLB
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGGSVP computes orthogonal matrices U, V and Q such that
*
*                   N-K-L  K    L
*   U'*A*Q =     K ( 0    A12  A13 )  if M-K-L >= 0;
*                L ( 0     0   A23 )
*            M-K-L ( 0     0    0  )
*
*                   N-K-L  K    L
*          =     K ( 0    A12  A13 )  if M-K-L < 0;
*              M-K ( 0     0   A23 )
*
*                 N-K-L  K    L
*   V'*B*Q =   L ( 0     0   B13 )
*            P-L ( 0     0    0  )
*
*  where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
*  upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
*  otherwise A23 is (M-K)-by-L upper trapezoidal.  K+L = the effective
*  numerical rank of the (M+P)-by-N matrix (A',B')'.  Z' denotes the
*  transpose of Z.
*
*  This decomposition is the preprocessing step for computing the
*  Generalized Singular Value Decomposition (GSVD), see subroutine
*  DGGSVD.
*
*  Arguments
*  =========
*
*  JOBU    (input) CHARACTER*1
*          = 'U':  Orthogonal matrix U is computed;
*          = 'N':  U is not computed.
*
*  JOBV    (input) CHARACTER*1
*          = 'V':  Orthogonal matrix V is computed;
*          = 'N':  V is not computed.
*
*  JOBQ    (input) CHARACTER*1
*          = 'Q':  Orthogonal matrix Q is computed;
*          = 'N':  Q is not computed.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  P       (input) INTEGER
*          The number of rows of the matrix B.  P >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and B.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, A contains the triangular (or trapezoidal) matrix
*          described in the Purpose section.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
*          On entry, the P-by-N matrix B.
*          On exit, B contains the triangular matrix described in
*          the Purpose section.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,P).
*
*  TOLA    (input) DOUBLE PRECISION
*  TOLB    (input) DOUBLE PRECISION
*          TOLA and TOLB are the thresholds to determine the effective
*          numerical rank of matrix B and a subblock of A. Generally,
*          they are set to
*             TOLA = MAX(M,N)*norm(A)*MAZHEPS,
*             TOLB = MAX(P,N)*norm(B)*MAZHEPS.
*          The size of TOLA and TOLB may affect the size of backward
*          errors of the decomposition.
*
*  K       (output) INTEGER
*  L       (output) INTEGER
*          On exit, K and L specify the dimension of the subblocks
*          described in Purpose.
*          K + L = effective numerical rank of (A',B')'.
*
*  U       (output) DOUBLE PRECISION array, dimension (LDU,M)
*          If JOBU = 'U', U contains the orthogonal matrix U.
*          If JOBU = 'N', U is not referenced.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U. LDU >= max(1,M) if
*          JOBU = 'U'; LDU >= 1 otherwise.
*
*  V       (output) DOUBLE PRECISION array, dimension (LDV,M)
*          If JOBV = 'V', V contains the orthogonal matrix V.
*          If JOBV = 'N', V is not referenced.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V. LDV >= max(1,P) if
*          JOBV = 'V'; LDV >= 1 otherwise.
*
*  Q       (output) DOUBLE PRECISION array, dimension (LDQ,N)
*          If JOBQ = 'Q', Q contains the orthogonal matrix Q.
*          If JOBQ = 'N', Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q. LDQ >= max(1,N) if
*          JOBQ = 'Q'; LDQ >= 1 otherwise.
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  TAU     (workspace) DOUBLE PRECISION array, dimension (N)
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (max(3*N,M,P))
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*
*  Further Details
*  ===============
*
*  The subroutine uses LAPACK subroutine DGEQPF for the QR factorization
*  with column pivoting to detect the effective numerical rank of the
*  a matrix. It may be replaced by a better rank determination strategy.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FORWRD, WANTQ, WANTU, WANTV
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEQPF, DGEQR2, DGERQ2, DLACPY, DLAPMT, DLASET,
     $                   DORG2R, DORM2R, DORMR2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      WANTU = LSAME( JOBU, 'U' )
      WANTV = LSAME( JOBV, 'V' )
      WANTQ = LSAME( JOBQ, 'Q' )
      FORWRD = .TRUE.
*
      INFO = 0
      IF( .NOT.( WANTU .OR. LSAME( JOBU, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTV .OR. LSAME( JOBV, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( WANTQ .OR. LSAME( JOBQ, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( N.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -10
      ELSE IF( LDU.LT.1 .OR. ( WANTU .AND. LDU.LT.M ) ) THEN
         INFO = -16
      ELSE IF( LDV.LT.1 .OR. ( WANTV .AND. LDV.LT.P ) ) THEN
         INFO = -18
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -20
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGSVP', -INFO )
         RETURN
      END IF
*
*     QR with column pivoting of B: B*P = V*( S11 S12 )
*                                           (  0   0  )
*
      DO 10 I = 1, N
         IWORK( I ) = 0
   10 CONTINUE

c      WRITE(*, *) 'Before call DGEQPF'
      CALL DGEQPF( P, N, B, LDB, IWORK, TAU, WORK, INFO )
c      WRITE(*, *) 'After call DGEQPF'
*
*     Update A := A*P
*
c      WRITE(*, *) 'Before call DLAPMT'
      CALL DLAPMT( FORWRD, M, N, A, LDA, IWORK )

*
*     Determine the effective rank of matrix B.
*
      L = 0
      DO 20 I = 1, MIN( P, N )
         IF( ABS( B( I, I ) ).GT.TOLB )
     $      L = L + 1
   20 CONTINUE
*
      IF( WANTV ) THEN
*
*        Copy the details of V, and form V.
*
         CALL DLASET( 'Full', P, P, ZERO, ZERO, V, LDV )
         IF( P.GT.1 )
     $      CALL DLACPY( 'Lower', P-1, N, B( 2, 1 ), LDB, V( 2, 1 ),
     $                   LDV )
         CALL DORG2R( P, P, MIN( P, N ), V, LDV, TAU, WORK, INFO )
      END IF
*
*     Clean up B
*
      DO 40 J = 1, L - 1
         DO 30 I = J + 1, L
            B( I, J ) = ZERO
   30    CONTINUE
   40 CONTINUE
      IF( P.GT.L )
     $   CALL DLASET( 'Full', P-L, N, ZERO, ZERO, B( L+1, 1 ), LDB )
*
      IF( WANTQ ) THEN
*
*        Set Q = I and Update Q := Q*P
*
         CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
         CALL DLAPMT( FORWRD, N, N, Q, LDQ, IWORK )
      END IF
*
      IF( P.GE.L .AND. N.NE.L ) THEN
*
*        RQ factorization of (S11 S12): ( S11 S12 ) = ( 0 S12 )*Z
*
         CALL DGERQ2( L, N, B, LDB, TAU, WORK, INFO )
*
*        Update A := A*Z'
*
         CALL DORMR2( 'Right', 'Transpose', M, N, L, B, LDB, TAU, A,
     $                LDA, WORK, INFO )
*
         IF( WANTQ ) THEN
*
*           Update Q := Q*Z'
*
            CALL DORMR2( 'Right', 'Transpose', N, N, L, B, LDB, TAU, Q,
     $                   LDQ, WORK, INFO )
         END IF
*
*        Clean up B
*
         CALL DLASET( 'Full', L, N-L, ZERO, ZERO, B, LDB )
         DO 60 J = N - L + 1, N
            DO 50 I = J - N + L + 1, L
               B( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE
*
      END IF
*
*     Let              N-L     L
*                A = ( A11    A12 ) M,
*
*     then the following does the complete QR decomposition of A11:
*
*              A11 = U*(  0  T12 )*P1'
*                      (  0   0  )
*
      DO 70 I = 1, N - L
         IWORK( I ) = 0
   70 CONTINUE
c      WRITE(*, *) 'Before call DGEQPF'
      CALL DGEQPF( M, N-L, A, LDA, IWORK, TAU, WORK, INFO )
c      WRITE(*, *) 'After call WANTU'
*
*     Determine the effective rank of A11
*
      K = 0
      DO 80 I = 1, MIN( M, N-L )
         IF( ABS( A( I, I ) ).GT.TOLA )
     $      K = K + 1
   80 CONTINUE
*
*     Update A12 := U'*A12, where A12 = A( 1:M, N-L+1:N )
*
      CALL DORM2R( 'Left', 'Transpose', M, L, MIN( M, N-L ), A, LDA,
     $             TAU, A( 1, N-L+1 ), LDA, WORK, INFO )
*


      IF( WANTU ) THEN
*
*        Copy the details of U, and form U
*
         CALL DLASET( 'Full', M, M, ZERO, ZERO, U, LDU )
         IF( M.GT.1 )
     $      CALL DLACPY( 'Lower', M-1, N-L, A( 2, 1 ), LDA, U( 2, 1 ),
     $                   LDU )
         CALL DORG2R( M, M, MIN( M, N-L ), U, LDU, TAU, WORK, INFO )
      END IF
*
      IF( WANTQ ) THEN
*
*        Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1
*
         CALL DLAPMT( FORWRD, N, N-L, Q, LDQ, IWORK )
      END IF
*
*     Clean up A: set the strictly lower triangular part of
*     A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0.
*
      DO 100 J = 1, K - 1
         DO 90 I = J + 1, K
            A( I, J ) = ZERO
   90    CONTINUE
  100 CONTINUE
      IF( M.GT.K )
     $   CALL DLASET( 'Full', M-K, N-L, ZERO, ZERO, A( K+1, 1 ), LDA )
*
      IF( N-L.GT.K ) THEN
*
*        RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1
*
         CALL DGERQ2( K, N-L, A, LDA, TAU, WORK, INFO )
*
         IF( WANTQ ) THEN
*
*           Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1'
*
            CALL DORMR2( 'Right', 'Transpose', N, N-L, K, A, LDA, TAU,
     $                   Q, LDQ, WORK, INFO )
         END IF
*
*        Clean up A
*
         CALL DLASET( 'Full', K, N-L-K, ZERO, ZERO, A, LDA )
         DO 120 J = N - L - K + 1, N - L
            DO 110 I = J - N + L + K + 1, K
               A( I, J ) = ZERO
  110       CONTINUE
  120    CONTINUE
*
      END IF
*
      IF( M.GT.K ) THEN
*
*        QR factorization of A( K+1:M,N-L+1:N )
*
         CALL DGEQR2( M-K, L, A( K+1, N-L+1 ), LDA, TAU, WORK, INFO )
*
         IF( WANTU ) THEN
*
*           Update U(:,K+1:M) := U(:,K+1:M)*U1
*
            CALL DORM2R( 'Right', 'No transpose', M, M-K, MIN( M-K, L ),
     $                   A( K+1, N-L+1 ), LDA, TAU, U( 1, K+1 ), LDU,
     $                   WORK, INFO )
         END IF
*
*        Clean up
*
         DO 140 J = N - L + 1, N
            DO 130 I = J - N + K + L + 1, M
               A( I, J ) = ZERO
  130       CONTINUE
  140    CONTINUE
*
      END IF
*
      RETURN
*
*     End of DGGSVP
*
      END


      subroutine dgauss (pos, spec, specmod, npoints, hw1e)

c     convolves input spectrum with gaussian slit function of specified
c     hw1e (half-width at 1/e intensity). Assumes input spectrum has
c     constant spacing.

      implicit double precision (a - h, o - z)
      parameter (maxpts = 2000)
      dimension pos (maxpts), spec (maxpts)
      dimension specmod (maxpts), slit (maxpts)

      if (hw1e .eq. 0) then
         write (*, *) ' No gaussian convolution applied!!!'
         do i = 1, npoints
            specmod (i) = spec (i)
         end do
         return
      end if

      emult = - 1.d0 / (hw1e**2)
      delpos = pos (2) - pos (1)

c     apply slit function convolution

c     calculate slit function values out to 0.001 times x0 value,
c     normalize so that sum = 1.

      slitsum = 1.d0
      slit0 = 1.d0
      do i = 1, 1000
         slit (i) = exp (emult * (delpos * i)**2)
         slitsum = slitsum + 2.d0 * slit (i)
         if (slit (i) / slit0 .le. 0.001) then
            nslit = i
c         write (*, *) ' nslit = ', nslit
            go to 40
         end if
      end do
 40   continue
      slit0 = slit0 / slitsum
      do i = 1, nslit
         slit (i) = slit (i) / slitsum
c       write (*, *) i, slit (i)
      end do

c     convolve spectrum.  reflect at endpoints.
      do i = 1, npoints
         specmod (i) = slit0 * spec (i)
         do j = 1, nslit
            nlo = i - j
            nhi = i + j
            if (nlo .lt. 0) nlo = - nlo
            if (nlo .eq. 0) nlo = 1
            if (nhi .gt. npoints) nhi = 2 * npoints - nhi
            specmod (i) = specmod (i) +
     *           slit (j) * (spec (nlo) + spec (nhi))
         end do
      end do

      return
      end

c  To get d1mach, mail netlib
c       send d1mach from core
      double precision function fmin(ax,bx,f,tol)
      double precision ax,bx,f,tol
c
c      an approximation  x  to the point where  f  attains a minimum  on
c  the interval  (ax,bx)  is determined.
c
c  input..
c
c  ax    left endpoint of initial interval
c  bx    right endpoint of initial interval
c  f     function subprogram which evaluates  f(x)  for any  x
c        in the interval  (ax,bx)
c  tol   desired length of the interval of uncertainty of the final
c        result (.ge.0.)
c
c  output..
c
c  fmin  abcissa approximating the point where  f  attains a
c        minimum
c
c      the method used is a combination of  golden  section  search  and
c  successive parabolic interpolation.  convergence is never much slower
c  than  that  for  a  fibonacci search.  if  f  has a continuous second
c  derivative which is positive at the minimum (which is not  at  ax  or
c  bx),  then  convergence  is  superlinear, and usually of the order of
c  about  1.324....
c      the function  f  is never evaluated at two points closer together
c  than  eps*abs(fmin)+(tol/3), where eps is  approximately  the  square
c  root  of  the  relative  machine  precision.   if   f   is a unimodal
c  function and the computed values of   f   are  always  unimodal  when
c  separated  by  at least  eps*abs(x)+(tol/3), then  fmin  approximates
c  the abcissa of the global minimum of  f  on the interval  ax,bx  with
c  an error less than  3*eps*abs(fmin)+tol.  if   f   is  not  unimodal,
c  then fmin may approximate a local, but perhaps non-global, minimum to
c  the same accuracy.
c      this function subprogram is a slightly modified  version  of  the
c  algol  60 procedure  localmin  given in richard brent, algorithms for
c  minimization without derivatives, prentice-hall, inc. (1973).
c
c
      double precision  a,b,c,d,e,eps,xm,p,q,r,tol1,t2,u,v,w,fu,fv,fw,
     2    fx,x,tol3
      double precision  dabs,dsqrt !,d1mach
c
c  c is the squared inverse of the golden ratio
      c=0.5d0*(3.0d0-dsqrt(5.0d0))
c
c  eps is approximately the square root of the relative machine
c  precision.
c
! the first statement can be changed to eps = 2.219209D-16
c   10 eps=d1mach(4)   
   10 eps=2.219209D-16 
      tol1=eps+1.0d0
      eps=dsqrt(eps)
c
      a=ax
      b=bx
      v=a+c*(b-a)
      w=v
      x=v
      e=0.0d0
      fx=f(x)
      fv=fx
      fw=fx
      tol3=tol/3.0d0
c
c  main loop starts here
c
   20 xm=0.5d0*(a+b)
      tol1=eps*dabs(x)+tol3
      t2=2.0d0*tol1
c
c  check stopping criterion
c
      if (dabs(x-xm).le.(t2-0.5d0*(b-a))) go to 190
      p=0.0d0
      q=0.0d0
      r=0.0d0
      if (dabs(e).le.tol1) go to 50
c
c  fit parabola
c
      r=(x-w)*(fx-fv)
      q=(x-v)*(fx-fw)
      p=(x-v)*q-(x-w)*r
      q=2.0d0*(q-r)
      if (q.le.0.0d0) go to 30
      p=-p
      go to 40
   30 q=-q
   40 r=e
      e=d
   50 if ((dabs(p).ge.dabs(0.5d0*q*r)).or.(p.le.q*(a-x))
     2          .or.(p.ge.q*(b-x))) go to 60
c
c  a parabolic-interpolation step
c
      d=p/q
      u=x+d
c
c  f must not be evaluated too close to ax or bx
c
      if (((u-a).ge.t2).and.((b-u).ge.t2)) go to 90
      d=tol1
      if (x.ge.xm) d=-d
      go to 90
c
c  a golden-section step
c
   60 if (x.ge.xm) go to 70
      e=b-x
      go to 80
   70 e=a-x
   80 d=c*e
c
c  f must not be evaluated too close to x
c
   90 if (dabs(d).lt.tol1) go to 100
      u=x+d
      go to 120
  100 if (d.le.0.0d0) go to 110
      u=x+tol1
      go to 120
  110 u=x-tol1
  120 fu=f(u)
c
c  update  a, b, v, w, and x
c
      if (fx.gt.fu) go to 140
      if (u.ge.x) go to 130
      a=u
      go to 140
  130 b=u
  140 if (fu.gt.fx) go to 170
      if (u.ge.x) go to 150
      b=x
      go to 160
  150 a=x
  160 v=w
      fv=fw
      w=x
      fw=fx
      x=u
      fx=fu
      go to 20
  170 if ((fu.gt.fw).and.(w.ne.x)) go to 180
      v=w
      fv=fw
      w=u
      fw=fu
      go to 20
  180 if ((fu.gt.fv).and.(v.ne.x).and.(v.ne.w)) go to 20
      v=u
      fv=fu
      go to 20
c
c  end of main loop
c
  190 fmin=x
      return
      end


c     Routine to perform Cholesky decomposition
c     From Numerical Recipes by W.H. Press [1992]
      SUBROUTINE choldc(a, n, np, p, errstat)
      INTEGER n,np
      DOUBLE PRECISION a(np, np),p(np)
      INTEGER i,j,k, errstat
      DOUBLE PRECISION sum

      errstat = 0
      do 13 i=1,n
        do 12 j=i,n
          sum=a(i,j)
          do 11 k=i-1,1,-1
            sum=sum-a(i,k)*a(j,k)
11        continue
          if(i.eq.j)then
            if(sum.le.0.d0) then
               errstat = 2
               return
             endif
            p(i)=sqrt(sum)
          else
            a(j,i)=sum/p(i)
          endif
12      continue
13    continue
      return
      END


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    routine to invert a matrix by Gaussian elimination
c     Ainv=inverse(A)    A(n,n)   Ainv(n,n)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine invert( A, Ainv, n, errstat)
       integer errstat, n
       double precision A(n,n), Ainv(n,n)
       double precision D(n,n*2)

       errstat = 0
c  initialize the reduction matrix
       n2 = 2*n
       do 1 i = 1,n
        do 2 j = 1,n
         D(i,j) = A(i,j)
         d(i,n+j) = 0.d0
 2      continue
        D(i,n+i) = 1.0d0
 1     continue

c  do the reduction 
       do 3 i = 1,n
        alpha = D(i,i)
        if(alpha .eq. 0.) then
           errstat = 2
           go to 300
         endif
        do 4 j = 1,n2
         D(i,j) = D(i,j)/alpha
 4      continue
        do 5 k = 1,n
         if((k-i).eq.0.d0) go to 5
         beta = D(k,i)
         do 6 j = 1,n2
          D(k,j) = D(k,j) - beta*D(i,j)
 6       continue
 5      continue
 3     continue

c  copy result into output matrix
       do 7 i = 1,n
        do 8 j = 1,n
         Ainv(i,j) = D(i,j+n)
 8      continue
 7     continue

       return
 300   print *,'*** ERROR: Singular matrix ***', i
       return
       end


      SUBROUTINE dsvdcmp1(a,m,n,mp,np,w,v)
      INTEGER m,mp,n,np,NMAX
      DOUBLE PRECISION a(mp,np),v(np,np),w(np)
      PARAMETER (NMAX=2500)
CU    USES dpythag
      INTEGER i,its,j,jj,k,l,nm
      DOUBLE PRECISION anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),dpythag
      g=0.0d0
      scale=0.0d0
      anorm=0.0d0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if(i.le.m)then
          do 11 k=i,m
            scale=scale+abs(a(k,i))
11        continue
          if(scale.ne.0.0d0)then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do 15 j=l,n
              s=0.0d0
              do 13 k=i,m
                s=s+a(k,i)*a(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
14            continue
15          continue
            do 16 k=i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scale=scale+abs(a(i,k))
17        continue
          if(scale.ne.0.0d0)then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0d0
              do 21 k=l,n
                s=s+a(j,k)*a(i,k)
21            continue
              do 22 k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.0d0)then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0d0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0d0
            v(j,i)=0.0d0
31        continue
        endif
        v(i,i)=1.0d0
        g=rv1(i)
        l=i
32    continue
      do 39 i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do 33 j=l,n
          a(i,j)=0.0d0
33      continue
        if(g.ne.0.0d0)then
          g=1.0d0/g
          do 36 j=l,n
            s=0.0d0
            do 34 k=l,m
              s=s+a(k,i)*a(k,j)
34          continue
            f=(s/a(i,i))*g
            do 35 k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
35          continue
36        continue
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0d0
38        continue
        endif
        a(i,i)=a(i,i)+1.0d0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0d0
          s=1.0d0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=dpythag(f,g)
            w(i)=h
            h=1.0d0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=w(k)
          if(l.eq.k)then
            if(z.lt.0.0d0)then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            goto 3
          endif
          if(its.eq.30) pause 'no convergence in svdcmp'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
          g=dpythag(f,1.0d0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0d0
          s=1.0d0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=dpythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
45          continue
            z=dpythag(f,h)
            w(j)=z
            if(z.ne.0.0d0)then
              z=1.0d0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0d0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      END

      FUNCTION dpythag1(a,b)
      DOUBLE PRECISION a,b,dpythag
      DOUBLE PRECISION absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        dpythag=absa*sqrt(1.0d0+(absb/absa)**2)
      else
        if(absb.eq.0.0d0)then
          dpythag=0.0d0
        else
          dpythag=absb*sqrt(1.0d0+(absa/absb)**2)
        endif
      endif
      return
      END


      SUBROUTINE dsvdcmp(a,m,n,mp,np,w,v)
      INTEGER m,mp,n,np,NMAX
      DOUBLE PRECISION a(mp,np),v(np,np),w(np)
      PARAMETER (NMAX=500)
CU    USES dpythag
      INTEGER i,its,j,jj,k,l,nm
      DOUBLE PRECISION anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),dpythag
      g=0.0d0
      scale=0.0d0
      anorm=0.0d0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if(i.le.m)then
          do 11 k=i,m
            scale=scale+abs(a(k,i))
11        continue
          if(scale.ne.0.0d0)then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do 15 j=l,n
              s=0.0d0
              do 13 k=i,m
                s=s+a(k,i)*a(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
14            continue
15          continue
            do 16 k=i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scale=scale+abs(a(i,k))
17        continue
          if(scale.ne.0.0d0)then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0d0
              do 21 k=l,n
                s=s+a(j,k)*a(i,k)
21            continue
              do 22 k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.0d0)then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0d0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0d0
            v(j,i)=0.0d0
31        continue
        endif
        v(i,i)=1.0d0
        g=rv1(i)
        l=i
32    continue
      do 39 i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do 33 j=l,n
          a(i,j)=0.0d0
33      continue
        if(g.ne.0.0d0)then
          g=1.0d0/g
          do 36 j=l,n
            s=0.0d0
            do 34 k=l,m
              s=s+a(k,i)*a(k,j)
34          continue
            f=(s/a(i,i))*g
            do 35 k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
35          continue
36        continue
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0d0
38        continue
        endif
        a(i,i)=a(i,i)+1.0d0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0d0
          s=1.0d0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=dpythag(f,g)
            w(i)=h
            h=1.0d0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=w(k)
          if(l.eq.k)then
            if(z.lt.0.0d0)then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            goto 3
          endif
          if(its.eq.30) then 
          print * , 'no convergence in svdcmp'; return
          !pause 'no convergence in svdcmp'
          ENDIF
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
          g=dpythag(f,1.0d0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0d0
          s=1.0d0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=dpythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
45          continue
            z=dpythag(f,h)
            w(j)=z
            if(z.ne.0.0d0)then
              z=1.0d0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0d0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0@.1Y..

      FUNCTION dpythag(a,b)
      DOUBLE PRECISION a,b,dpythag
      DOUBLE PRECISION absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        dpythag=absa*sqrt(1.0d0+(absb/absa)**2)
      else
        if(absb.eq.0.0d0)then
          dpythag=0.0d0
        else
          dpythag=absb*sqrt(1.0d0+(absa/absb)**2)
        endif
      endif
      return
      END
     
