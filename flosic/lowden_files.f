C UTEP Electronic Structure Lab (2020)
       SUBROUTINE LOWSIC(ISPN,NDH,M,OVER,HAM,FMAT,EVAL,SC1)
C   Mark Pederson 14 September 2001
C TO USE
C PLACE OVERLAP MATRIX IN OVER.
C NEW WAVEFUNCTIONS ARE RETURNED IN HAM...
       USE MPIDAT1,only : IRANK
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION HAM(NDH,NDH),OVER(NDH,NDH),FMAT(NDH,NDH)
     &,EVAL(NDH),SC1(NDH)
!      DIMENSION OVER2(NDH,NDH)
       REAL*8,allocatable :: OVER2(:,:)
       LOGICAL CHECKIT
       LOGICAL LOWEVAL
       LOGICAL AUTOCORR
       DATA CHECKIT/.FALSE./
       DATA LOWEVAL/.FALSE./
       DATA AUTOCORR/.FALSE./
       DATA ISPN_SAV/0/
       SAVE
C      IF(ISPN.EQ.ISPN_SAV)THEN
C        ISPN_SAV=ISPN
C      ELSE
C        RETURN
C      END IF
C      PRINT*,'HELLO'
       IF(CHECKIT)THEN
         PRINT*," NOTE IF CHECKIT=.TRUE. THERE IS N**4 SCALING."
         PRINT*," NORMAL CALLS SHOULD BE WITH CHECKIT=.FALSE."
       END IF
C THa: Check if we want to have very small LOWDEN 
C      eigenvalues autocorrected
C      this may be usefull during startup and bad
C      FO positions 
C      CAUTION: *may* lead to undefined results !
       LOWEVAL  = .FALSE.
       AUTOCORR = .FALSE.
       INQUIRE(FILE='LOWAUTO',EXIST=AUTOCORR)
           
C      OPEN(32,FILE='SCRLOW',FORM='UNFORMATTED')
C      WRITE(32)((OVER(J,I),J=1,M),I=1,M) 
       ALLOCATE(OVER2(NDH,NDH))
       FORALL(I=1:M, J=1:M)
         HAM(J,I) = 0.0D0
         OVER2(J,I)=OVER(J,I)
       END FORALL
       DO I=1,M
         HAM(I,I)=1.0D0
       END DO
C      DO I=1,M
C        DO J=1,M
C          HAM(J,I)=0.0D0
C          OVER2(J,I)=OVER(J,I)
C        END DO  
C        HAM(I,I)=1.0D0
C      END DO
       WRITE(6+IRANK,*)'CALLING DIAGGE:',NDH,M
       CALL DIAGGE_FO(NDH,M,OVER,HAM,EVAL,SC1,1)
       FORALL(I=1:M, J=1:M)
         FMAT(J,I)=OVER(J,I)
       END FORALL
           
C      DO I=1,M
C        DO J=1,M
C          FMAT(J,I)=OVER(J,I)
C        END DO
C      END DO
       WRITE(6+IRANK,*)' LOWDEN OVERLAP EIGENVALUES:'
       WRITE(6+IRANK,20)(EVAL(I),I=1,M)

       DO I=1,M
C        check for too small eigenvalues
C        if there are some, stop after printing them all out
         IF(EVAL(I).LT.1.0D-8) THEN
           LOWEVAL=.TRUE.
C          PRINT *, "> SET <"
         ENDIF
       END DO
       IF(LOWEVAL .eqv. .TRUE.) THEN
         WRITE(6+IRANK,*)
     &    "-----------------------------------------------"
         WRITE(6+IRANK,*)
     &    "> LOWSIC: LOWDEN EIGENVALUES TOO SMALL (< 1e-8)"
         DO I=1,M
C          PRINT *, '>> ', EVAL(I)
           IF(EVAL(I).LT.1.0D-8) THEN
             WRITE(6+IRANK,*) ">  EV: ", I , "=", EVAL(I)
             IF (AUTOCORR .eqv. .TRUE.) THEN
               EVAL(I) = 1.0D-8
               WRITE(6+IRANK,*) ">    corrected to", EVAL(I)
             END IF
           END IF
         END DO
         WRITE(6+IRANK,*)
     &    ">         BAD POSITIONS IN FRMORB ??           "
         WRITE(6+IRANK,*)
     &    "-----------------------------------------------"
         IF(.NOT.AUTOCORR) CALL STOPIT
       END IF
C      normalize
       DO I=1,M
         EVAL(I)=1.0D0/SQRT(EVAL(I))
       END DO
           
C CHECK FOR ORTHOGONALITY:
       IF(CHECKIT)THEN
         DO I=1,M
           DO J=1,M
             SC1(J)=0.0D0
             DO K=1,M
               SC1(J)=SC1(J)+OVER(K,I)*OVER(K,J)
             END DO 
           END DO
C          PRINT 20,(SC1(J),J=1,M)
         END DO
       END IF
C PERFORM BACK TRANSFORMATION:
       DO J=1,M
         DO I=1,M
           HAM(I,J)=0.0D0
           DO N=1,M
             HAM(I,J)=HAM(I,J)+OVER(J,N)*OVER(I,N)*EVAL(N)
           END DO
         END DO
       END DO
C RELOAD OVERLAP MATRIX:
C      REWIND(32)
C      READ(32)((OVER(J,I),J=1,M),I=1,M) 
C      CLOSE(32,STATUS='DELETE')
       DO J=1,M
         DO I=1,M
           OVER(J,I)=OVER2(J,I)
         END DO
       END DO
!YY
       DEALLOCATE(OVER2)
C CHECK AGAINST 3-ORDER LOWDEN:
       IF(CHECKIT)THEN
         WRITE(6+IRANK,*)'COMPARISON TO 3RD ORDER LOWDEN...'
         DO I=1,M
           DO J=1,M
             SC1(J)=-OVER(J,I)/2.0D0
           END DO
           SC1(I)=1.0D0
           DO J=1,M
             DO K=1,M
               IF(K.NE.J.AND.K.NE.I)THEN
                 SC1(J)=SC1(J)+OVER(K,J)*OVER(K,I)*(3./8.)
               END IF
             END DO 
           END DO
C          PRINT 20,(HAM(J,I),J=1,M)
C          PRINT 20,(SC1(J  ),J=1,M)
           PRINT*,' '
         END DO
C CHECK FOR ORTHOGONALITY:
C
C        PRINT*,'TRANSFORMED OVERLAP MATRIX:'
         DO I=1,M
           DO J=1,M
             SC1(J)=0.0D0
             DO K=1,M
               DO L=1,M
                 SC1(J)=SC1(J)+HAM(L,I)*HAM(K,J)*OVER(L,K)
               END DO
             END DO
           END DO
C          PRINT 20,(SC1(J),J=1,M)
         END DO
       END IF
 20    FORMAT(' ',5g15.6)
       END
