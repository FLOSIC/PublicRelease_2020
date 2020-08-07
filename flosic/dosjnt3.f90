! UTEP Electronic Structure Lab (2020)
      SUBROUTINE DOSJNT3(IWF)
        use debug1
        USE DOSJNT_MOD,only : H,SPTOT,SPDIP,SOS_FREQ,ESTEP,EALP, &
         SOS_POL,VFAC,ENJD,TEMP,NSPEC
        use common5,only : NWF,NWFS,EFERMI,EVLOCC
        use mpidat1,only : IRANK
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: IWF
        INTEGER :: IPTS,ISPN,JSPN,JWF
        REAL*8 :: BACKW,DIP,EDIFF,EDIFF27,ERG,EV1
        REAL*8 :: EV2,FACT,FCT,FUNC,OCCI,OCCJ
        REAL*8 :: TIME1,TIME2
        REAL*8,EXTERNAL :: FFERMI

        ISPN=1
        IF (IWF.GT.NWFS(1)) ISPN=2
        OCCI=FFERMI(EVLOCC(IWF),EFERMI(ISPN),TEMP)
        DO JWF=IWF,NWF
         JSPN=1
         IF (JWF.GT.NWFS(1)) JSPN=2
         OCCJ=FFERMI(EVLOCC(JWF),EFERMI(JSPN),TEMP)
         DIP=H(JWF,IWF,1)**2 &
            +H(JWF,IWF,2)**2 &
            +H(JWF,IWF,3)**2
         IF(ISPN.NE.JSPN)DIP=0.0D0
         EDIFF=EVLOCC(JWF)-EVLOCC(IWF)
         FCT=MAX(OCCI*(1.-OCCJ),OCCJ*(1.-OCCI))
!         write(IRANK+10,121)evlocc(iwf),evlocc(jwf),EDIFF, &
!                h(jwf,iwf,1),h(jwf,iwf,2),h(jwf,iwf,3),dip,OCCI,OCCJ
 121     format(1x,3(f9.4,1x), 6(f13.6))

         IF(FCT*DIP.GT.0.00001)THEN
             EV1=MIN(EVLOCC(IWF),EVLOCC(JWF))
             EV2=MAX(EVLOCC(IWF),EVLOCC(JWF))
             EDIFF=EV1-EV2
             SOS_POL = SOS_POL + DIP/EDIFF
         END IF
!
! LOOP OVER ENERGY GRID TO GET INTENSITIES
!
         BACKW= 1.0D0
         IF (IWF.EQ.JWF) BACKW= 0.0D0
         FACT = OCCI*(1.0d0-OCCJ)
         DO IPTS=1,NSPEC
          ERG=(IPTS-1)*ESTEP+ENJD
!          FUNC=VFAC*EXP(-EALP*(ERG+EDIFF)**2)
! LB: Use the following expression to get same numbers as excite
! VFAC is set to one
          FUNC=VFAC*EXP(-((ERG+EDIFF)/0.01)**2)
          SPTOT(IPTS)=SPTOT(IPTS)+FACT*FUNC
          SPDIP(IPTS)=SPDIP(IPTS)+FACT*FUNC*DIP
          IF(ERG<1E-5) THEN
            SOS_FREQ(IPTS) = 0.0
          ELSE
            SOS_FREQ(IPTS) = SOS_FREQ(IPTS)+FACT*DIP/ERG
          ENDIF
!          FUNC=BACKW*VFAC*EXP(-EALP*(ERG+EDIFF)**2)
! LB: Use the following expression to get same numbers as excite
! VFAC is set to one
          FUNC=BACKW*VFAC*EXP(-((ERG+EDIFF)/0.01)**2)
          SPTOT(IPTS)=SPTOT(IPTS)+OCCJ*(1.0D0-OCCI)*FUNC
          SPDIP(IPTS)=SPDIP(IPTS)+OCCJ*(1.0D0-OCCI)*FUNC*DIP
         END DO
        END DO
        RETURN
      END SUBROUTINE DOSJNT3
