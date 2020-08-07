C UTEP Electronic Structure Lab (2020)
      SUBROUTINE DEFAULTKEY
C
C     SET FLAGS TO ITS DEFAULT VALUES.
C
C     BY ULISES REVELES, SEP. 2013.
C
C     ------------------------------------------------------------------
C
C     GLOBAL VARIABLES:
C
      use common2,only : CHNET,EGMIN,EGMAX,INPTYP,MAXDR,MOLDEN,OPTTYP,
     %                   PRTOPT,PRTPRI,PRTMOS,PRTPMAT,PRTSMAT,UNITS,
     %                   RTRUST,STEPTYP,SPNNET,NBO
C
      IMPLICIT NONE
C
C     ------------------------------------------------------------------
C     --- SET DEFAULT VALUES FOR CHARGE AND SPIN STATE ---
C
      CHNET = 0.0
      SPNNET = 0.0
C
C     --- DEFAULT GEOMETRY OPTIONS ---
C
      UNITS = 'ANGSTROMS'
C
C     --- DEFAULT OPTIMIZATION OPTIONS ---
C
      MAXDR   = 0.3
      RTRUST  = 0.1
C
      INPTYP  = 'CARTESIAN'
      OPTTYP  = 'CARTESIAN'
      STEPTYP = 'LEVENBERG'
C
C     --- DEFAULT PRINT OPTIONS ---
C
      MOLDEN = 'JMOL'
      PRTOPT = .FALSE.
      PRTPRI = .FALSE.
      PRTMOS = .FALSE.
      PRTPMAT = .FALSE.
      PRTSMAT = .FALSE.
C
C     --- DEFAULT CHARGE POPULATION ANAYSIS ---
C
      NBO = .TRUE.
C
C     --- DEFAULT ENERGY GAP FOR PRINTING ORBITALS IN MOLDEN ---
C
      EGMIN = 0.5
      EGMAX = 0.5
C
C     ------------------------------------------------------------------
C
      END
