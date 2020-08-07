C UTEP Electronic Structure Lab (2020)
c-----------------------------------------------------------------------------------------
c
c     SUBROUTINE ELE_INFO(ngrp,nspn,ele_up,ele_dn)
c       This subroutine will return number of electrons and number 
c    of symmetry operations. This information is needed in DIAGGE routine if
c    partial diagonalization option is to be used.
c   
      SUBROUTINE ELE_INFO(lgrp,lspn,ele_up,ele_dn)
       use common2,only : E_UP, E_DN, NSPN
       use common3,only : RMAT, NGRP
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:43 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       SAVE
       real*8  ele_up, ele_dn
       integer lgrp,lspn
       ele_up = e_up
       ele_dn = e_dn
       lgrp = ngrp
       lspn = nspn
       return
      end
