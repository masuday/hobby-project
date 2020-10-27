module lapack_wrapper

!implicit none
integer,parameter,private :: rk=selected_real_kind(15,300)

contains

!
! Extracted from LAPACK95
! http://www.netlib.org/lapack95/
!
SUBROUTINE DSYEV_F95SUB( A, W, JOBZ, UPLO, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
!   USE LASUB_PRECISION, ONLY: WP => DP
!   USE LASUB_AUXMOD, ONLY: ERINFO, LSAME
!   USE F77_LAPACKSUB, ONLY: SYEV_F77 => LASUB_SYEV, ILAENV_F77 => ILAENV
!
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
   INTEGER,PARAMETER :: WP=selected_real_kind(15,300)
   INTEGER :: ILAENV
!  .. CHARACTER ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ, UPLO
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(INOUT) :: A(:,:)
   REAL(WP), INTENT(OUT) :: W(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!      LA_SYEV and LA_SYEVD compute all eigenvalues and, optionally, all
! eigenvectors of a real symmetric matrix A.
!      LA_HEEV and LA_HEEVD compute all eigenvalues and, optionally, all
! eigenvectors of a complex Hermitian matrix A.
!      LA_SYEVD and LA_HEEVD use a divide and conquer algorithm. If 
! eigenvectors are desired, they can be much faster than LA_SYEV and 
! LA_HEEV for large matrices but use more workspace.
! 
! =========
! 
!       SUBROUTINE LA_SYEV / LA_HEEV / LA_SYEVD / LA_HEEVD( A, W, &
!                       JOBZ=jobz, UPLO=uplo, INFO=info )
!           <type>(<wp>), INTENT(INOUT) :: A(:,:)
!           REAL(<wp>), INTENT(OUT) :: W(:)
!           CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ, UPLO
!           INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!           <type> ::= REAL | COMPLEX
!           <wp> ::= KIND(1.0) | KIND(1.0D0)
! 
! 
! Arguments
! =========
! 
! A      (input/output) REAL or COMPLEX square array, shape (:,:).
!        On entry, the matrix A.
!        If UPLO = 'U', the upper triangular part of A contains the upper
!        triangular part of the matrix A. If UPLO = 'L', the lower 
!        triangular part of A contains the lower triangular part of the
!        matrix A.
!        On exit:
!        If JOBZ = 'V', then the columns of A contain the orthonormal 
!        eigenvectors of the matrix A in the order of the eigenvalues.
!        If JOBZ = 'N', then the upper triangle (if UPLO = 'U') or the 
!        lower triangle (if UPLO = 'L') of A, including the diagonal, is
!        destroyed.
! W      (output) REAL array, shape (:) with size(W) = size(A,1).
!        The eigenvalues in ascending order.
! JOBZ   Optional (input) CHARACTER(LEN=1).
!        = 'N': Computes eigenvalues only;
!        = 'V': Computes eigenvalues and eigenvectors.
!        Default value: 'N'.
! UPLO   Optional (input) CHARACTER(LEN=1).
!        = 'U': Upper triangle of A is stored;
!        = 'L': Lower triangle of A is stored.
!        Default value: 'U'.
! INFO   Optional (output) INTEGER.
!        = 0: successful exit.
!        < 0: if INFO = -i, the i-th argument had an illegal value
!        > 0: if INFO = i, then i off-diagonal elements of an
!        intermediate tridiagonal form did not converge to zero.
!        If INFO is not present and an error occurs, then the program is 
!        terminated with an error message.
!-----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_SYEV'
   CHARACTER(LEN=6), PARAMETER :: BSNAME = 'DSYTRD'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LJOBZ, LUPLO
   INTEGER :: N, LINFO, LD, ISTAT, ISTAT1, LWORK, NB
!  .. LOCAL ARRAYS ..
   REAL(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC MAX, PRESENT
!  .. EXECUTABLE STATEMENTS ..
   N = SIZE( A, 1 ); LINFO = 0; ISTAT = 0; LD = MAX(1,N)
   IF( PRESENT(JOBZ) ) THEN
      LJOBZ = JOBZ
   ELSE
      LJOBZ = 'N'
   END IF
   IF( PRESENT(UPLO) ) THEN
      LUPLO = UPLO
   ELSE
      LUPLO = 'U'
   END IF
!  .. TEST THE ARGUMENTS
   IF( SIZE( A, 2 ) /= N .OR. N < 0 )THEN
      LINFO = -1
   ELSE IF( SIZE( W ) /= N )THEN
      LINFO = -2
   ELSE IF( .NOT.LSAME(LJOBZ,'N') .AND. .NOT.LSAME(LJOBZ,'V') )THEN
      LINFO = -3
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN
      LINFO = -4
   ELSE IF( N > 0 )THEN
!  .. DETERMINE THE WORKSPACE
      !NB = ILAENV_F77( 1, BSNAME, LUPLO, N, -1, -1, -1 )
      NB = ILAENV( 1, BSNAME, LUPLO, N, -1, -1, -1 )
      IF( NB <= 1 .OR. NB >= N )THEN
         NB = 1
      END IF
     LWORK = (2+NB)*N
      ALLOCATE(WORK(LWORK), STAT=ISTAT)
      IF( ISTAT /= 0 )THEN
         LWORK = 3*N-1
         ALLOCATE(WORK(LWORK), STAT=ISTAT)
         IF( ISTAT /= 0 ) THEN
            LINFO = - 100
         ELSE
            CALL ERINFO( -200, SRNAME, LINFO )
         ENDIF
      ENDIF
!
      IF( LINFO == 0 )THEN
         !CALL SYEV_F77( LJOBZ, LUPLO, N, A, LD, W, WORK, LWORK, LINFO )
         CALL DSYEV( LJOBZ, LUPLO, N, A, LD, W, WORK, LWORK, LINFO )
      ENDIF
      DEALLOCATE(WORK, STAT=ISTAT1)
   ENDIF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DSYEV_F95SUB

      LOGICAL FUNCTION LSAME( CA, CB )
!
!  PURPOSE
!  =======
!
!  LSAME  TESTS IF CA IS THE SAME LETTER AS CB REGARDLESS OF CASE.
!
!  PARAMETERS
!  ==========
!
!  CA      (INPUT) CHARACTER*1
!  CB      (INPUT) CHARACTER*1
!          CHARACTERS TO BE COMPARED.
!
!  .. SCALAR ARGUMENTS ..
      CHARACTER( LEN=1 ), INTENT(IN) :: CA, CB
!  .. PARAMETERS ..
      INTEGER, PARAMETER      :: IOFF=32
!  .. LOCAL SCALARS ..
      INTEGER                 :: INTA, INTB, ZCODE
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC                  ICHAR
!
!  .. EXECUTABLE STATEMENTS ..
!
!  TEST IF THE CHARACTERS ARE EQUAL
!
      LSAME = CA == CB
!
!  NOW TEST FOR EQUIVALENCE
!
      IF( .NOT.LSAME )THEN
!
!     USE 'Z' RATHER THAN 'A' SO THAT ASCII CAN BE DETECTED ON PRIME
!     MACHINES, ON WHICH ICHAR RETURNS A VALUE WITH BIT 8 SET.
!     ICHAR('A') ON PRIME MACHINES RETURNS 193 WHICH IS THE SAME AS
!     ICHAR('A') ON AN EBCDIC MACHINE.
!
         ZCODE = ICHAR( 'Z' )
!
         INTA = ICHAR( CA )
         INTB = ICHAR( CB )
!
         IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 )THEN
!
!        ASCII IS ASSUMED - ZCODE IS THE ASCII CODE OF EITHER LOWER OR
!        UPPER CASE 'Z'.
!
            IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
            IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
         ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 )THEN
!
!        EBCDIC IS ASSUMED - ZCODE IS THE EBCDIC CODE OF EITHER LOWER OR
!        UPPER CASE 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.                         &
!    &       INTA.GE.145 .AND. INTA.LE.153 .OR.                         &
     &       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.                         &
     &       INTB.GE.145 .AND. INTB.LE.153 .OR.                         &
     &       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
         ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 )THEN
!
!        ASCII IS ASSUMED, ON PRIME MACHINES - ZCODE IS THE ASCII CODE
!        PLUS 128 OF EITHER LOWER OR UPPER CASE 'Z'.
!
            IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
         ENDIF
         LSAME = INTA == INTB
      ENDIF
      END FUNCTION LSAME

      SUBROUTINE ERINFO(LINFO, SRNAME, INFO, ISTAT)
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. IMPLICIT STATEMENT ..
         IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
         CHARACTER( LEN = * ), INTENT(IN)              :: SRNAME
         INTEGER             , INTENT(IN)              :: LINFO
         INTEGER             , INTENT(OUT), OPTIONAL   :: INFO
         INTEGER             , INTENT(IN), OPTIONAL    :: ISTAT
!  .. EXECUTABLE STATEMENTS ..
!         IF( ( LINFO < 0 .AND. LINFO > -200 ) .OR.                     &
!    &       ( LINFO > 0 .AND. .NOT.PRESENT(INFO) ) )THEN
      IF( ( ( LINFO < 0 .AND. LINFO > -200 ) .OR. LINFO > 0 )           &
     &           .AND. .NOT.PRESENT(INFO) )THEN
        WRITE (*,*) 'Program terminated in LAPACK95 subroutine ',SRNAME
        WRITE (*,*) 'Error indicator, INFO = ',LINFO
        IF( PRESENT(ISTAT) )THEN
          IF( ISTAT /= 0 ) THEN
            IF( LINFO == -100 )THEN
              WRITE (*,*) 'The statement ALLOCATE causes STATUS = ',    &
     &                    ISTAT
            ELSE
              WRITE (*,*) 'LINFO = ', LINFO, ' not expected'
            END IF
          END IF   
        END IF
        STOP
         ELSE IF( LINFO <= -200 ) THEN
           WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
           WRITE(*,*) '*** WARNING, INFO = ', LINFO, ' WARNING ***'
           IF( LINFO == -200 )THEN
             WRITE(*,*)                                                 &
     &        'Could not allocate sufficient workspace for the optimum'
             WRITE(*,*)                                                 &
     &        'blocksize, hence the routine may not have performed as'
             WRITE(*,*) 'efficiently as possible'
         ELSE
           WRITE(*,*) 'Unexpected warning'
         END IF
           WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
        END IF
        IF( PRESENT(INFO) ) THEN
          INFO = LINFO
        END IF
      END SUBROUTINE ERINFO

end module lapack_wrapper
