C UTEP Electronic Structure Lab (2020)
      SUBROUTINE EZSTART
C
C     MODIFIED BY ULISES REVELES, JULY 2013.
C
C     ------------------------------------------------------------------
C
      use xmol,only : AU2ANG
      USE GLOBAL_INPUTS,only : SYMMETRY1
      IMPLICIT REAL*8 (A-H,O-Z)
C
      CHARACTER*30, NAME
C
      LOGICAL EXIST
      INTEGER STREXT
C
      DIMENSION R(5,675) 
      DIMENSION G(3,3,60),IND(60)
C
C
C     ------------------------------------------------------------------
C
C     --- OPEN AND REWIND CLUSTER FILE ---
C
      OPEN(80,FILE='CLUSTER')   
      REWIND(80)
C
      NAME = ' '
C
      READ(80,'(A30)') NAME
      IF (NAME(1:1).NE.'@') THEN
        WRITE(6,*)'IMPROPER CALL TO EZSTART'
        CALL STOPIT
      END IF
C
      IF (STREXT(NAME).GT.29) THEN 
        WRITE(6,*)'XMOL FILE MUST BE'
        WRITE(6,*)'LESS THAN 30 CHARACTERS'
      END IF
C
C     --- REMOVE FIRST CHARACTER OF LINE ---
C
      DO I=2,30
        NAME(I-1:I-1) = NAME(I:I)
      END DO
C
C     --- CLOSE CLUSTER FILE ---
C
      CLOSE (80)
C
C     --- OPEN FILE WITH COORDINATES ---
C
        OPEN(90,FILE=NAME)
C
        READ(90,*)NAT
        READ(90,*)
C
C     --- CHECK MAXIMUM NUMBER OF ATOMS ---
C
        IF (NAT.GT.675) THEN
          WRITE(6,*)'NAT TO BIG IN EZSTART'
        END IF
C
        do iat=1,nat
          read(90,*)nuc,(r(j,iat),j=1,3)
          print 80,(r(j,iat),j=1,3),nuc
          r(5,iat)=dfloat(nuc)
        end do
        CLOSE(90)
C DETERMINE IF WE WANT TO DETECT SYMMETRY OR NOT
      CALL CHECK_INPUTS
      IF(SYMMETRY1)THEN
        MAT=NAT
C
C     --- OPEN CLUSTER FILE ---
C
        OPEN(80,FILE='CLUSTER')
C
        WRITE(80,78)
        WRITE(80,79)
        WRITE(80,*)MAT,' Number of inequivalent atoms'
C
        DO IAT=1,NAT
          WRITE(80,80)(R(J,IAT)/AU2ANG,J=1,3),NINT(R(5,IAT))
        END DO
C
      ELSE
C
C     --- INTIALIZATION ---
C
      G(1:3,1:3,1:48) = 0.0D0
C
      G(1,1,1) = 1.0d0
      G(2,2,1) = 1.0d0
      G(3,3,1) = 1.0d0
      G(1,1,4) =-1.0d0
      G(2,2,4) = 1.0d0
      G(3,3,4) = 1.0d0
      G(1,2,2) = 1.0d0
      G(2,3,2) = 1.0d0
      G(3,1,2) = 1.0d0
      G(1,2,3) = 1.0d0
      G(2,1,3) = 1.0d0
      G(3,3,3) = 1.0d0
C
      MGP = 4
      DO ITIME=1,5
        NGP = MGP
        DO IGP=1,NGP
C
C         print 20,((g(i,j,igp),j=1,3),i=1,3)
C20       format(' ',3f12.4)
C         write(6,*)' '
C
          DO JGP=1,NGP
            MGP=MGP+1
            DO I=1,3
              DO J=1,3
                G(I,J,MGP) = 0.0D0
                DO K=1,3
                  G(I,J,MGP)=G(I,J,MGP)+G(I,K,IGP)*G(K,J,JGP)
                END DO
              END DO
            END DO
C
            no=0
            do kgp=1,mgp-1
              err=0.0d0
              do i=1,3
                do j=1,3
                  err=err+abs(g(i,j,mgp)-g(i,j,kgp))
                end do
              end do
C
              if(err.lt.0.001)no=1
            end do
C
            if(no.eq.1)mgp=mgp-1
C
          END DO
        END DO
        WRITE(6,*)itime,ngp,mgp
      END DO
      ngp=0
      do kgp=1,mgp
        kat=nat
        do iat=1,nat
          kat=kat+1
          do j=1,3
            r(j,kat)=0.0d0
            do k=1,3
              r(j,kat)=r(j,kat)+g(j,k,kgp)*r(k,iat)
            end do
          end do
          errmin=1.0d0
          do jat=1,nat
          dd=distance(r(1,kat),r(1,jat))
     &                    +abs(r(5,iat)-r(5,jat))
          errmin=min(dd,errmin)
c          errmin=min(distance(r(1,kat),r(1,jat)),errmin)
        end do
        if(errmin.lt.0.001)kat=kat-1
      end do
      if(kat.eq.nat)then
        ngp=ngp+1
        ind(ngp)=kgp
      end if
      end do
C
C     --- OPEN GRPMAT FILE ---
C
      OPEN(80,FILE='GRPMAT')
C
      write(80,*)ngp,' Number of group elements'
      do igp=1,ngp
        write(80,75)((g(i,j,ind(igp)),j=1,3),i=1,3)
        do i=1,3
          do j=1,3
            g(i,j,igp)=g(i,j,ind(igp))
          end do
        end do
        write(80,*)' '
      end do
C
 75   format(' ',3g15.6)
C
C     --- CLOSE GRPMAT FILE ---
C
      CLOSE (80)
C
      mat=nat+1
      do iat=1,nat
        r(4,iat)=1.0d0
        errmin=1.0d0
        do igp=1,ngp
          do j=1,3
            r(j,mat)=0.0d0
            do k=1,3
              r(j,mat)=r(j,mat)+g(j,k,igp)*r(k,iat)
            end do
          end do
          do jat=1,iat-1   
            dd=distance(r(1,mat),r(1,jat))
     &         +abs(r(5,jat)-r(5,iat))
            errmin=min(dd,errmin)
          end do
        end do
        if(errmin.lt.0.001)r(4,iat)=-1.0d0
      end do
      mat=0
      do iat=1,nat
        if(r(4,iat).gt.0.0d0)mat=mat+1
      end do

C
C     --- OPEN CLUSTER FILE ---
C
      OPEN(80,FILE='CLUSTER')
C
      WRITE(80,78)
      WRITE(80,79)
      WRITE(80,*)MAT,' Number of inequivalent atoms'
C
      DO IAT=1,NAT
        IF (R(4,IAT).GT.0.0d0) THEN
          WRITE(80,80)(R(J,IAT)/AU2ANG,J=1,3),NINT(R(5,IAT))
        END IF
      END DO
      ENDIF
C
      WRITE(80,81)CHRG,SPIN
      CLOSE(80)
C
C     --- FORMATS ---
C
 30   FORMAT(a20)
 78   FORMAT('GGA-PBE*GGA-PBE')
 79   FORMAT('GRP')
 80   FORMAT(' ',3f12.4,I5,' ALL UPO ')
 81   FORMAT(' ',2f12.4,' Net Charge and Moment')
C
C     ------------------------------------------------------------------
C
      END
