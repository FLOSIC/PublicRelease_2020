! UTEP Electronic Structure Lab (2020)

subroutine FODONMSHGEN(ITBEG,NCALC,IMESH,FOMSH_INFO)
use common2,only : RIDT,RCNT,IFUIDT,IFUCNT,NIDENT,NCNT
use FRM,only     : BFRM,RESULTS,LFRM,DEBDAX
INCLUDE  'PARAMA2'

integer, intent(in) :: ITBEG,NCALC,IMESH
real*8,  intent(out):: FOMSH_INFO
logical :: EXIST
integer :: IFUIDT_SAV,IFUCNT_SAV,NSPN, &
           LSPN,IFRM,IATOM
real*8 :: RIDT_SAV,RCNT_SAV,FOMSH_CUT,FOMSH_ISCORE,FOMSH_MDIST, &
          DFONUC

! These vars are needed for the SIC O(N) mesh code
DIMENSION    RIDT_SAV(3,MAX_IDENT)
DIMENSION    RCNT_SAV(3,MX_CNT)
DIMENSION    IFUIDT_SAV(MAX_IDENT)
DIMENSION    IFUCNT_SAV(MX_CNT)
DIMENSION    FOMSH_INFO(MXSPN,MAX_OCC)
LOGICAL  ::  FOMSH_ENABLE
INTEGER  ::  NIDENT_SAV, NCNT_SAV
CHARACTER*15 :: FOMSH_NAME
CHARACTER*5  :: VM_NAME


!C     THa: generate VMSH for every Fermi-orbital
!     we have to tweak and tmp store the following vars
!      COMMON/NUCLEI/RIDT(3,MAX_IDENT),RCNT(3,MX_CNT)
!    &  ,IFUIDT(MAX_IDENT),IFUCNT(MX_CNT),NIDENT,NCNT
INQUIRE(FILE='FODONMSH',EXIST=EXIST)
IF (EXIST) THEN
  FOMSH_ENABLE  = .TRUE.
  RIDT_SAV(:,:) = RIDT(:,:)
  RCNT_SAV(:,:) = RCNT(:,:)
  IFUIDT_SAV(:) = IFUIDT(:)
  IFUCNT_SAV(:) = IFUCNT(:)
  NIDENT_SAV    = NIDENT
  NCNT_SAV      = NCNT


  FOMSH_CUT    = 6.0D0     ! cut of for inclusion into vmsh
  FOMSH_ISCORE = 0.5D0     ! closer that that -> is core level
  FOMSH_MDIST  = 0.0D0
  VM_NAME      = 'VMOLD'

!C     THa: generate VMSH for every Fermi-orbital

! for every fermi orbital
!   calculate distance to all atoms
!   if distance is smaller than a given th:
!     include this atom in structure S
!   generate mesh for structure S

! this generates too much meshes (one for every fermi orbital)

! different approach:
!  check distance of FO to every atom
!  build list that: flags FO as core, holds id to closest atom
!  then build partial meshes for all atoms

  !PRINT *, 'NFNCT', NFNCT
  !DO IFNCT=1,NFNCT
  !  PRINT *, 'INFNCT', IFNCT, N_POS(IFNCT)
  !END DO
  !
  !
  !DO IFNCT=1,NCNT
  !  PRINT *, 'INFNCT', IFNCT, IFUIDT(IFNCT)
  !END DO
  print *,'test',NCNT_SAV,LFRM(1),LFRM(NSPN)

  DO LSPN=1,2
    DO IFRM=1,LFRM(LSPN)

      RIDT(:,:)  = 0.0D0
      RCNT(:,:)  = 0.0D0
      IFUIDT(:)  = 0.0D0
      IFUCNT(:)  = 0.0D0
      NIDENT     = 0
      NCNT       = 0

      FOMSH_MDIST = 1.0D99 ! set to some very large value
      DO IATOM=1,NCNT_SAV
        DFONUC = (BFRM(1,IFRM,LSPN) - RCNT_SAV(1,IATOM))**2 +  &
                 (BFRM(2,IFRM,LSPN) - RCNT_SAV(2,IATOM))**2 +  &
                 (BFRM(3,IFRM,LSPN) - RCNT_SAV(3,IATOM))**2  
        DFONUC = SQRT(DFONUC)
        ! assign distance to nearest core
        IF (DFONUC .LT. FOMSH_MDIST) FOMSH_MDIST = DFONUC

        PRINT *, LSPN, IFRM, IATOM, ':', DFONUC,   &
          IFUIDT_SAV(IATOM)

        IF (DFONUC .LT. FOMSH_CUT) THEN
          NIDENT          = NIDENT+1
          NCNT            = NCNT+1
          RCNT(1:3,NCNT) = RCNT_SAV(1:3,IATOM)
          RIDT(1:3,NCNT) = RIDT_SAV(1:3,IATOM)
          IFUIDT(NCNT)   = IFUIDT_SAV(IATOM)
          IFUCNT(NCNT)   = IFUCNT_SAV(IATOM)
        END IF

      END DO

      PRINT *, 'NIDENT, NCNT', NIDENT, NCNT
      print *,'FOMSH_MDIST',FOMSH_MDIST

! The commented section below does the equivalent of the section
! above. The commented section was moved from main.ftn within the
! orbital loop.
!            DO J=1,3
!              POSITION(J)=BFRM(J,IORBX,JSPNX)
!            END DO
!C FIND CLOSEST ATOM:
!            OPEN(83,FILE='XMOL.DAT')
!            REWIND(83)
!            DSTMIN=1.0D30
!            JXAT=0
!            READ(83,*)NXAT
!            READ(83,*)
!            DO IXAT=1,NXAT
!              READ(83,*)MMM,X1,Y1,Z1
!              DST=(POSITION(1)-X1/0.529177)**2
!     &           +(POSITION(2)-Y1/0.529177)**2
!     &           +(POSITION(3)-Z1/0.529177)**2
!              DST=SQRT(DST)
!              IF(DST.LT.DSTMIN)THEN
!                DSTMIN=DST
!                JXAT=IXAT
!              END IF
!            END DO
!            CLOSE(83)
!            FOMSH_ENABLE
!            CALL READMESH(JXAT,POSITION,DSTMIN)
!            CALL APOTNL_SIC(TOTQNUM,JSPNX,IORBX)
!            CALL READMESH(0,POSITION,DSTMIN)




      ! call the original mesh code, but only if this is
      ! not a core orbital
      ! (core orbitals will only get a simple radial mesh)
      FOMSH_INFO(LSPN,IFRM) = FOMSH_MDIST
      !YY Create mesh if FOD is far away from the nearest nucleus
      !print *,'ifstatement',FOMSH_MDIST,FOMSH_ISCORE
      IF (FOMSH_MDIST .GT. FOMSH_ISCORE) THEN
        ! assign name to mesh that is beeing to generated
        FOMSH_NAME(1:7)='VMSHFOD'
        WRITE(FOMSH_NAME(8:11),601)LSPN
        WRITE(FOMSH_NAME(12:15),601)IFRM
 601    format(i4.4)


        IF (IMESH .EQ. 1) THEN
         CALL VMESH(ITBEG,NCALC)
        ELSE
         CALL DVPMESH(ITBEG,NCALC)
        ENDIF
! testing. don't rename if not created. cmd
        print *,'calling rename',VM_NAME, ' ',FOMSH_NAME
        CALL RENAME(VM_NAME, FOMSH_NAME)
! 
      ELSE
        PRINT *, '-->', LSPN, IFUIDT,  &
                 'core state, use rad mesh'
      END IF

!     print *,'calling rename',VM_NAME, ' ',FOMSH_NAME
!     CALL RENAME(VM_NAME, FOMSH_NAME)
      ! save the
      !CALL STOPIT

    END DO
  END DO


  RIDT(:,:)  = RIDT_SAV(:,:)
  RCNT(:,:)  = RCNT_SAV(:,:)
  IFUIDT(:)  = IFUIDT_SAV(:)
  IFUCNT(:)  = IFUCNT_SAV(:)
  NIDENT     = NIDENT_SAV
  NCNT       = NCNT_SAV
END IF


END SUBROUTINE
