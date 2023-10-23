C UTEP Electronic Structure Lab (2020)
C
C ****************************************************************
C
        SUBROUTINE SETDFTYP(STRING)
C
C DIRK POREZAG, AUGUST 1997
C DEFINE DENSITY FUNCTIONAL
C
       use common2,only : IGGA, IDFTYP
       use global_inputs, only : libxc1
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:58 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, IEND, ILEN, ISTART, ISTGGA, IXC, NDFTYP, FUNC_ID
       REAL*8 :: SYMBOL 
         SAVE
         PARAMETER(NDFTYP=15)
         PARAMETER(ISTGGA=6)
         CHARACTER*35 STRING
         CHARACTER*6  NAMDF(0:NDFTYP)
         CHARACTER*25 NAMPR(0:NDFTYP)
         DATA NAMDF/'NONE  ','PERZUN','CEPERL','RPA   ','WIGNER',
     &              'GUNLUN','PW91  ','PBE   ','REVPBE','RPBE  ',
     &              'B88   ','SCAN  ','MSCAN', 'RSCAN', 'R2SCAN',
     &              'PBEIBP'/
         DATA NAMPR/'NONE                     ',
     &              'PERDEW-ZUNGER 81         ',
     &              'CEPERLEY                 ',
     &              'RPA                      ',
     &              'WIGNER                   ',
     &              'GUNNARSSON-LUNDQVIST     ',
     &              'PERDEW-WANG 91           ',
     &              'PERDEW-BURKE-ERNZERHOF 96',
     &              'REV PBE                  ',
     &              'R PBE                    ',
     &              'BECKE 88                 ',
     &              'SCAN                     ',
     &              'modSCAN                  ',
     &              'regularizedSCAN          ',
     &              'r2SCAN                   ',
     &              'PBE integration by parts '/
         CHARACTER*40 DEFSTR
         DATA DEFSTR/'[A.D. Becke, J.Chem.Phys.96, 2155, 1992]'/
C
C LOOP OVER EXCHANGE AND CORRELATION 
C
         IGGA(1)= -1
         IGGA(2)= -1

! YY. Libxc functional string should be family_kind_type e.g GGA_X_PBE*GGA_C_PBE
         call CHECK_INPUTS
         libxc1=.false.
         !use built in functionals
         if(libxc1) then
            print *, "LIBXC DISABLED IN CODE. TERMINATING"
            call stopit
         end if

         ISTART=1
         DO IXC=1,2
          IF      (STRING(ISTART:ISTART+2) .EQ. 'LDA') THEN
           IGGA(IXC)=0
          ELSE IF (STRING(ISTART:ISTART+2) .EQ. 'GGA') THEN
           IGGA(IXC)=1
          ELSE IF (STRING(ISTART:ISTART+3) .EQ. 'MGGA') THEN
           IGGA(IXC)=3
           ISTART=ISTART+1
          ELSE
           GOTO 900
          END IF
          IF (STRING(ISTART+3:ISTART+3) .NE. '-') GOTO 900
          ISTART=ISTART+4
          DO I=ISTART,35
           IF ((STRING(I:I) .EQ. '*') .OR. 
     &         (STRING(I:I) .EQ. ' ')) GOTO 10
          END DO
          GOTO 900
   10     IEND=I-1 
          ILEN=IEND-ISTART+1
          IF (ILEN .GT. 6) GOTO 900
          DO I=0,NDFTYP
           IF (STRING(ISTART:IEND) .EQ. NAMDF(I)(1:ILEN)) GOTO 20
          END DO
          GOTO 900
   20     IDFTYP(IXC)=I
C
C CHECK GGA AVAILABILITY
C
          IF ((IGGA(IXC) .EQ. 1) .AND. (IDFTYP(IXC) .NE. 0)
     &                           .AND. (IDFTYP(IXC) .LT. ISTGGA)) THEN
           write(6,*)'SETDFTYP: FUNCTIONAL ',STRING(ISTART:IEND)
           write(6,*)'IS NOT A GRADIENT-CORRECTED FUNCTIONAL'
           CALL STOPIT
          END IF
          ISTART=IEND+2
         END DO
         IF (IGGA(1) .LT. 0) GOTO 900
         IF (IGGA(2) .LT. 0) GOTO 900
C
C PRINT RESULTS
C
         IF (IDFTYP(1) .EQ. 0) THEN
          PRINT '(A)', 'EXCHANGE    FUNCTIONAL: NONE'
         ELSE IF (IGGA(1) .EQ. 0) THEN
          PRINT '(A)', 'EXCHANGE    FUNCTIONAL: LDA'
         ELSE IF (IGGA(1) .EQ. 3) THEN
          PRINT '(2A)', 'EXCHANGE    FUNCTIONAL: MGGA ',NAMPR(IDFTYP(1))
         ELSE
           IF(IDFTYP(1).EQ.10) THEN
            PRINT '(2A)','EXCHANGE FUNCTIONAL: GGA-B88 ',DEFSTR
           ELSE
            PRINT '(2A)','EXCHANGE    FUNCTIONAL: GGA ',NAMPR(IDFTYP(1))
           ENDIF
         END IF
         IF (IDFTYP(2) .EQ. 0) THEN
          PRINT '(A)', 'CORRELATION FUNCTIONAL: NONE'
         ELSE IF (IGGA(2) .EQ. 0) THEN
          PRINT '(2A)','CORRELATION FUNCTIONAL: LDA ',NAMPR(IDFTYP(2))
         ELSE IF (IGGA(2) .EQ. 3) THEN
          PRINT '(2A)','CORRELATION FUNCTIONAL: MGGA ',NAMPR(IDFTYP(2))
         ELSE
          PRINT '(2A)','CORRELATION FUNCTIONAL: GGA ',NAMPR(IDFTYP(2))
         END IF

C REMIND USER THAT MESH AND MIXING SETTING ARE TURNED ON (SEE
C check_input.F91 FILE)
         IF(IDFTYP(1).EQ.11 .OR. IDFTYP(2).EQ.11) THEN
          PRINT '(A)','> SCAN functional detected. <'
          PRINT '(A)','> SCAN mesh will be used. < '
         END IF
         IF(IGGA(1).EQ.3 .OR. IGGA(2).EQ.3) THEN
          PRINT '(A)','> Meta-GGA detected. <' 
          PRINT '(A)','> Hamiltonian mixing will be used. <'
         END IF

         RETURN
  900    write(6,*)'SETDFTYP: INVALID FUNCTIONAL TYPE:',STRING
         CALL STOPIT
        END
