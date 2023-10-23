C UTEP Electronic Structure Lab (2020)
C
C
CJK2
       SUBROUTINE RHOFFTPAR(ALPIV,ALPJV,AIV,AJV,DMTV,MUIV,MUJV,NCNTP)
       use formfakm
       use debug1
       use common2,only : NCNT, ISPN, NSPN
       use mpidat1,only : IRANK
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:58 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NCNTP, MUIV, MUJV, ICNT, IKPT, J, MUI, MUJ
       REAL*8 :: ALPIV , ALPJV, AIV, AJV, DMTV, ADDON1, ADDON2, SS
       SAVE       
       DIMENSION ALPIV(NCNTP),ALPJV(NCNTP)
       DIMENSION MUIV (NCNTP),MUJV (NCNTP)
       DIMENSION AIV(3,NCNTP),AJV(3,NCNTP)
       DIMENSION DMTV(10,10,2,NCNTP),SS(10,10,2)
C       COMMON/TMP2/ RKPT(3,MXKPT),RHOKPT(MXKPT,2,2),NKPT          
       LOGICAL FIRST
       DATA FIRST/.TRUE./
       write(6,*)'RFFT : ', NCNT,NKPT,irank
       DO ICNT=1,NCNT
        IF(DEBUG) THEN
        write(6,*)ALPIV(ICNT),ALPJV(ICNT)
        write(6,*)(AIV(j,ICNT),j=1,3),'  AI'
        write(6,*)(AjV(j,ICNT),j=1,3),'  Aj'
        write(6,*)MUIV (ICNT),MUJV (ICNT),'MU'
        ENDIF
        DO IKPT=1,NKPT
        CALL RHOFFT(ALPIV(ICNT),ALPJV(ICNT),AIV(1,ICNT),AJV(1,ICNT),
     &            SS,RKPT(1,IKPT)) 
        DO ISPN=1,NSPN
         ADDON1=0.0D0
         ADDON2=0.0D0
         DO MUJ=1,MUJV(ICNT)
          DO MUI=1,MUIV(ICNT)
           ADDON1=ADDON1+DMTV(MUI,MUJ,ISPN,ICNT)*SS(MUI,MUJ,1)
           ADDON2=ADDON2+DMTV(MUI,MUJ,ISPN,ICNT)*SS(MUI,MUJ,2)
          END DO
         END DO
         RHOKPT(IKPT,1,ISPN)=RHOKPT(IKPT,1,ISPN)+ADDON1
         RHOKPT(IKPT,2,ISPN)=RHOKPT(IKPT,2,ISPN)+ADDON2
         END DO                 !ISPN
        END DO                  !IKPT   
       END DO                   !ICNT
       RETURN
       END
