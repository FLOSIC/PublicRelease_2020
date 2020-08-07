C UTEP Electronic Structure Lab (2020)
C
C this defines an external potential of the form
C
C Vext= Sum_N Vext_i
C Vext_i= C(i)/r**M(i) * (x-Ax(i))**Nx(i)*exp(-Gx(i)*(x-Ax(i))**2) *
C                      * (y-Ay(i))**Ny(i)*exp(-Gy(i)*(y-Ay(i))**2) *
C                      * (z-Az(i))**Nz(i)*exp(-Gz(i)*(z-Az(i))**2)  
C
C===========EXTPOT==================================================
C NTERMS                # NUMBER OF TERMS'
C C M                   # PREFACTOR C(i), POWER M(i) of 1/r**M(i)')
C Nx Ny Nz              # POWER OF POLYNOMIALS NX(i), NY(i),NZ(i)')
C Ax Ay Az              # CENTER OF EXT POT Ax(i),Ay(i),Az(i)')
C gx gy gz              # EXP OF ENVELOPE Gx(i),Gy(i),Gz(i)')
C the last 4 lines are repeated for every NTERM
C==================================================================
C
      SUBROUTINE READEXT
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MAXT=10)
      COMMON/EXPOT/C(MAXT),M(MAXT),N(3,MAXT),A(3,MAXT),G(3,MAXT),NTERMS
      LOGICAL EXIST
      CHARACTER*1 AXIS(3)
      DATA AXIS/'X','Y','Z'/
      INQUIRE(FILE='EXTPOT',EXIST=EXIST)
      IF(.NOT.EXIST)RETURN
      OPEN(90,FILE='EXTPOT',FORM='FORMATTED')
      READ(90,*)NTERMS
      IF(NTERMS.GT.MAXT) THEN
       write(6,*)'MAXT IN EXTPOT TOO SMALL : ',MAXT,NTERMS
       CALL STOPIT
      ENDIF
      write(6,*)'EXTERNAL POTENTIAL FROM EXTPOT USED'
      DO I=1,NTERMS
       READ(90,*)C(I),M(I)
       READ(90,*)(N(J,I),J=1,3)
       READ(90,*)(A(J,I),J=1,3)
       READ(90,*)(G(J,I),J=1,3)
      END DO
      REWIND(90)
      WRITE(90,'(I4,23X,A)')NTERMS,'# NUMBER OF TERMS'
      DO I=1,NTERMS
       WRITE(90,10)C(I),M(I)
       WRITE(90,11)(N(J,I),J=1,3)
       WRITE(90,12)(A(J,I),J=1,3)
       WRITE(90,13)(G(J,I),J=1,3)
      END DO
      WRITE(90,*)
      WRITE(90,*)
      DO I=1,NTERMS
       WRITE(90,20)C(I),M(I),(AXIS(J),A(J,I),N(J,I),G(J,I),
     &  AXIS(J),A(J,I),J=1,3)
      END DO
 10   FORMAT(F8.4,I4,15X,'# PREFACTOR C(i), POWER M(i) of 1/r**M(i)')
 11   FORMAT(3I4,15X,'# POWER OF POLYNOMIALS NX(i), NY(i),NZ(i)')
 12   FORMAT(3F8.4,3X,'# CENTER OF EXT POT Ax(i),Ay(i),Az(i)')
 13   FORMAT(3F8.4,3X,'# EXP OF ENVELOPE Gx(i),Gy(i),Gz(i)')
 20   FORMAT(' +',F6.2,'/R**',I1,3(/,'*(',A1,'-',F6.2,')**',I1,
     &   '*EXP(-',F6.2,'*(',A1,'-',F6.2,')**2)'))
      CLOSE(90)
      RETURN
      END                                      
