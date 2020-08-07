C UTEP Electronic Structure Lab (2020)
! Conversion to implicit none.  Raja Zope Wed Nov  9 16:19:44 MST 2016
      SUBROUTINE WPARAMS
!       INCLUDE  'PARAMAS'  
        INCLUDE  'PARAMA2'  
! Conversion to implicit none.  Raja Zope Sun Oct 23 21:40:28 MDT 2016

       write(6,*)'LIST OF PARAMETERS USED FOR THIS ADVENTURE'
       write(6,*)'MXSPN=(',MXSPN,')'
       write(6,*)'MAX_PTS=(',MAX_PTS ,')'
       write(6,*)'MAX_FUSET=(',MAX_FUSET,')'
       write(6,*)'MAX_IDENT=(',MAX_IDENT,')'
       write(6,*)'MX_CNT=(','MX_CNT',')'
       write(6,*)'MXATMS=(',MXATMS,')'
       write(6,*)'MX_CNT=(',MX_CNT,')'
       write(6,*)'NMUPMAX=(',NMUPMAX,')'
       write(6,*)'MAX_BARE=(',MAX_BARE,')'
       write(6,*)'MAX_CON=(',MAX_CON,')'
       write(6,*)'MAX_OCC=(',MAX_OCC   ,')'
       write(6,*)'MAXUNSYM=(',MAXUNSYM,')'
       write(6,*)'NDH=(',NDH    ,')'
       write(6,*)'NDH_TOT=(',NDH_TOT,')'
       write(6,*)'MX_GRP=(',MX_GRP,')'
       write(6,*)'MAX_REP=(',MAX_REP,')'
       write(6,*)'MAXSLC=(',MAXSLC ,')'
       write(6,*)'LOCMAX=(',LOCMAX ,')'
       write(6,*)'ISMAX=(',ISMAX,')'
       write(6,*)'MAXSYMSALC=(',MAXSYMSALC,')'
       write(6,*)'MXRPSP=(',MXRPSP,')'
       write(6,*)'MXPTAB=(',MXPTAB,')'
       write(6,*)'MXLPSP=(',MXLPSP,')'
       write(6,*)'MXKPT=(',MXKPT,')'
       write(6,*)'MPBLOCK=(',MPBLOCK ,')'
       write(6,*)'MPBLOCKXC=(',MPBLOCKXC ,')'
       write(6,*)'MXPOISS=(',MXPOISS,')'
       write(6,*)'NSPEED=(',NSPEED ,')'
       write(6,*)'LDIM=(',LDIM,')'
       write(6,*)'CUTEXP=(',CUTEXP ,')'
       write(6,*)'MX_SPH=(',MX_SPH,')'
       write(6,*)'MAX_VIRT_PER_SYM=(',MAX_VIRT_PER_SYM,')'
       write(6,*)'MAXGFIT=(',MAXGFIT,')'
       write(6,*)'MAXLSQF=(',MAXLSQF    ,')'
       write(6,*)'MTEMP=(',MTEMP,')'
       write(6,*)'KRHOG=(',KRHOG,')'
       write(6,*)'C   DONE....'
       OPEN(70,FILE='PARASAV')
       REWIND(70)
       WRITE(70,*)'       IMPLICIT NONE'
       WRITE(70,*)'       INTEGER, PARAMETER :: dp = KIND(1.d0)'
       WRITE(70,*)'       INTEGER, PARAMETER :: MXSPN=    ',MXSPN
       WRITE(70,*)'       INTEGER, PARAMETER :: MAX_PTS=  ',MAX_PTS
       WRITE(70,*)'       INTEGER, PARAMETER :: MAX_FUSET=',MAX_FUSET
       WRITE(70,*)'       INTEGER, PARAMETER :: MAX_IDENT=',MAX_IDENT
       WRITE(70,*)'       INTEGER, PARAMETER :: MX_CNT= ',MX_CNT
       WRITE(70,*)'       INTEGER, PARAMETER :: MXATMS= MX_CNT'
       WRITE(70,*)'       INTEGER, PARAMETER :: NMUPMAX=  ',NMUPMAX
       WRITE(70,*)'       INTEGER, PARAMETER :: MAX_BARE=',MAX_BARE
       WRITE(70,*)'       INTEGER, PARAMETER :: MAX_CON= ',MAX_CON
       WRITE(70,*)'       INTEGER, PARAMETER :: MAX_OCC=',MAX_OCC
       WRITE(70,*)'       INTEGER, PARAMETER :: MAXUNSYM=',MAXUNSYM
       WRITE(70,*)'       INTEGER, PARAMETER :: NDH=  ',NDH
       WRITE(70,*)'       INTEGER, PARAMETER :: NDH_TOT= NDH*(NDH+1)/2'
       WRITE(70,*)'       INTEGER, PARAMETER :: MX_GRP=  ',MX_GRP
       WRITE(70,*)'       INTEGER, PARAMETER :: MAX_REP= ',MAX_REP
       WRITE(70,*)'       INTEGER, PARAMETER :: MAXSLC=  ',MAXSLC
       WRITE(70,*)'       INTEGER, PARAMETER :: LOCMAX=  ',LOCMAX
       WRITE(70,*)'       INTEGER, PARAMETER :: ISMAX=   ',ISMAX
       WRITE(70,*)'       INTEGER, PARAMETER :: MAXSYMSALC=',MAXSYMSALC
       WRITE(70,*)'       INTEGER, PARAMETER :: MXRPSP= ',MXRPSP
       WRITE(70,*)'       INTEGER, PARAMETER :: MXPTAB= ',MXPTAB
       WRITE(70,*)'       INTEGER, PARAMETER :: MXLPSP=   ',MXLPSP
       WRITE(70,*)'       INTEGER, PARAMETER :: MXKPT=    ',MXKPT
       WRITE(70,*)'       INTEGER, PARAMETER :: MPBLOCK=',MPBLOCK
       WRITE(70,*)'       INTEGER, PARAMETER :: MPBLOCKXC=',MPBLOCKXC
       WRITE(70,*)'       INTEGER, PARAMETER :: MXPOISS=',MXPOISS
       WRITE(70,*)'       INTEGER, PARAMETER :: NSPEED=  ',NSPEED
       WRITE(70,*)'       INTEGER, PARAMETER :: LDIM=     ',LDIM
       WRITE(70,*)'       REAL(dp),PARAMETER :: CUTEXP=',CUTEXP
       WRITE(70,*)'       INTEGER, PARAMETER :: MX_SPH=',MX_SPH
       WRITE(70,*)'       INTEGER, PARAMETER :: MAX_VIRT_PER_SYM=',
     &                           MAX_VIRT_PER_SYM
       WRITE(70,*)'       INTEGER, PARAMETER :: MAXGFIT= ',MAXGFIT
       WRITE(70,*)'       INTEGER, PARAMETER :: MAXLSQF= ',MAXLSQF
       WRITE(70,*)'       INTEGER, PARAMETER :: MTEMP=',MTEMP
       WRITE(70,*)'       INTEGER, PARAMETER :: KRHOG=   ',KRHOG
       CLOSE(70)
       RETURN
       END SUBROUTINE WPARAMS
