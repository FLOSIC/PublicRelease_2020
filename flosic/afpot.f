C UTEP Electronic Structure Lab (2020)
C
C 
      SUBROUTINE AFPOT(NSPN,RHOG,POT,POTIN)
C AFPOT adds a spin dependent potential of the form
C V(s_up)=V(s_up)+c(s_up)*exp(-alpha*(r-a)**2)
C V(s_dn)=V(s_dn)+c(s_dn)*exp(-alpha(r-a)**2)
C at sites given in AFPOTSYM/AFPOT
C
      use mesh1,only : wmsh,R=>rmsh,nmsh

! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:33 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NSPN, MAXT, I, IPTS, ISPN, J, K, MSITES, NATOMS,
     & NTERMS
       REAL*8 :: SYMBOL , RHOG, POT, POTIN, A, AFTMP, AHELP, B, C, EPS,
     & G, RR, XX, YY, ZZ
      PARAMETER(MAXT=10)
C      COMMON/MESH/WMSH(MAX_PTS),R(3,MAX_PTS),NMSH
      DIMENSION C(2,MAXT),A(3,MAXT),G(MAXT)
      DIMENSION POT(*),B(3,MX_GRP),MSITES(1)
      DIMENSION POTIN(*),RHOG(MAX_PTS,KRHOG,MXSPN)
      DIMENSION AHELP(2)
      LOGICAL FIRST,EXIST
      DATA FIRST/.TRUE./
      DATA EPS/ 1.0D-7/
      SAVE
      IF(FIRST)THEN
        NTERMS=0
        INQUIRE(FILE='AFPOT',EXIST=EXIST)
       IF(EXIST)THEN
        OPEN(43,FILE='AFPOT')
        READ(43,*)NTERMS
         DO I=1,NTERMS
         READ(43,*)(A(J,I),J=1,3)
         READ(43,*)G(I),(C(ISPN,I),ISPN=1,NSPN)
         END DO
        CLOSE(43)
       END IF                     
C        FIRST=.FALSE.
      END IF
C
      DO I=1,NTERMS
      CALL GASITES(1,A(1,I),NATOMS,B,MSITES)
        DO J=1,NATOMS
        write(6,*)' AFPOT: ',I,(B(K,J),K=1,3),C(1,I),C(2,I),NSPN
          DO IPTS=1,NMSH
          XX=R(1,IPTS)-B(1,J)
          YY=R(2,IPTS)-B(2,J)
          ZZ=R(3,IPTS)-B(3,J)
          RR=XX*XX+YY*YY+ZZ*ZZ
          DO ISPN=1,NSPN
           AHELP(ISPN)=EXP(-G(I)*RR)*C(ISPN,I)
           POT(IPTS+(ISPN-1)*NMSH)=
     &     POT(IPTS+(ISPN-1)*NMSH)+AHELP(ISPN)
          END DO
          IF(NSPN.EQ.1) THEN
           AHELP(2)=AHELP(1)
          ENDIF
          AFTMP=AHELP(1)*(RHOG(IPTS,1,1)-RHOG(IPTS,1,NSPN))
          AFTMP=AFTMP+AHELP(2)*RHOG(IPTS,1,NSPN)
          IF(RHOG(IPTS,1,1).GT.1D-6) THEN
          POTIN(IPTS)=POTIN(IPTS)+AFTMP/RHOG(IPTS,1,1)
          ENDIF
          END DO
        END DO                    
      END DO                    
      RETURN
      END      
