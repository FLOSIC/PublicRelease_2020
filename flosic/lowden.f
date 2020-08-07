C UTEP Electronic Structure Lab (2020)
           SUBROUTINE LOWDEN(NBAS,M)

C   Mark Pederson 14 September 2001
C TO USE
C PLACE OVERLAP MATRIX IN OVER.
C NEW WAVEFUNCTIONS ARE RETURNED IN HAM...
         use for_diag1
         IMPLICIT REAL*8 (A-H,O-Z)
         LOGICAL CHECKIT
         REAL*8 EVAL(NBAS)
         ALLOCATE(ASC1(M),STAT=IERR)
         IF(IERR.NE.0)THEN
           WRITE(6,*)'LOWDEN:ERROR ALLOCATING ASC1'
         ENDIF

           DATA CHECKIT/.FALSE./
           IF(CHECKIT)THEN
           PRINT*," NOTE THAT IF CHECKIT=.TRUE. THERE IS N**4 SCALING."
           PRINT*," NORMAL CALLS SHOULD BE WITH CHECKIT=.FALSE."
           END IF
              OPEN(32,FILE='SCRLOW',FORM='UNFORMATTED')
              WRITE(32)((AOVER(J,I),J=1,M),I=1,M) 
C              eps=1.0d-3
C               DO I=1,M
C                  OVER(I,I)=(1.0d0+eps)*OVER(I,I)
C               END DO
                  DO I=1,M
                  DO J=1,M
                  AHAM(J,I)=0.0D0
                  END DO  
                  AHAM(I,I)=1.0D0
                  END DO
           CALL DIAGGE(NBAS,M,AOVER,AHAM,EVAL,ASC1,1)
           PRINT*,' LOWDEN OVERLAP EIGENVALUES:'
           PRINT 20,(EVAL(I),I=1,M)
C                   DO I=1,M
C                   EVAL(I)=1.0D0/SQRT(EVAL(I))
C                   END DO
C CHECK FOR ORTHOGONALITY:
          IF(CHECKIT)THEN
               DO I=1,M
               DO J=1,M
                 ASC1(J)=0.0D0
                 DO K=1,M
                 ASC1(J)=ASC1(J)+AOVER(K,I)*AOVER(K,J)
                 END DO 
               END DO
               PRINT 20,(ASC1(J),J=1,M)
               END DO
           END IF
C PERFORM BACK TRANSFORMATION:
               DO J=1,M
                  DO I=1,M
                  AHAM(I,J)=0.0D0
                   DO N=1,M
C                   IF (EVAL(N).GT.1.0D-6) THEN
                   EV=1.0d0/SQRT(EVAL(N))
C                   HAM(I,J)=HAM(I,J)+OVER(J,N)*OVER(I,N)*EVAL(N)
                   AHAM(I,J)=AHAM(I,J)+AOVER(J,N)*AOVER(I,N)*EV
C                   END IF
                   END DO
                  END DO
               END DO
C RELOAD OVERLAP MATRIX:
              REWIND(32)
              READ(32)((AOVER(J,I),J=1,M),I=1,M) 
              DO I=1,M
                  ASC1(I)=0.0d0
                  DO J=1, M
                   DO K=1,M
                   ASC1(I)=ASC1(I)+AHAM(J,I)*AOVER(J,K)*AHAM(K,I)
                   END DO
                  END DO
               END DO
           PRINT*,' LOWDEN NORMALIZATION: ' 
           PRINT 20,(ASC1(I),I=1,M)
               DO I=1,M
                 ANORM=1.0d0/SQRT(ASC1(I))
                 DO J=1,M
                   AHAM(J,I)=AHAM(J,I)*ANORM
                 END DO
               END DO
              DO I=1,M
                  ASC1(I)=0.0d0
                  DO J=1, M
                   DO K=1,M
                   ASC1(I)=ASC1(I)+AHAM(J,I)*AOVER(J,K)*AHAM(K,I)
                   END DO
                  END DO
               END DO
           PRINT*,' AFTER NORMALIZATION: ' 
           PRINT 20,(ASC1(I),I=1,M)
              CLOSE(32,STATUS='DELETE')
C CHECK AGAINST 3-ORDER LOWDEN:
              IF(CHECKIT)THEN
               PRINT*,'COMPARISON TO 3RD ORDER LOWDEN...'
               DO I=1,M
                  DO J=1,M
                  ASC1(J)=-AOVER(J,I)/2.0D0
                  END DO
                  ASC1(I)=1.0D0
                  DO J=1,M
                  DO K=1,M
                  IF(K.NE.J.AND.K.NE.I)THEN
                  ASC1(J)=ASC1(J)+AOVER(K,J)*AOVER(K,I)*(3./8.)
                  END IF
                  END DO 
                  END DO
               PRINT 20,(AHAM(J,I),J=1,M)
               PRINT 20,(ASC1(J  ),J=1,M)
               PRINT*,' '
               END DO
C CHECK FOR ORTHOGONALITY:
C
               PRINT*,'TRANSFORMED OVERLAP MATRIX:'
               DO I=1,M
               DO J=1,M
                   ASC1(J)=0.0D0
                   DO K=1,M
                   DO L=1,M
                   ASC1(J)=ASC1(J)+AHAM(L,I)*AHAM(K,J)*AOVER(L,K)
                   END DO
                   END DO
               END DO
               PRINT 20,(ASC1(J),J=1,M)
               END DO
              END IF
 20            FORMAT(' ',5g15.6)
          DEALLOCATE(ASC1,STAT=IERR)
          IF(IERR.NE.0)THEN
            WRITE(6,*)'LOWDEN:ERROR DEALLOCATING ASC1'
          ENDIF


         END

