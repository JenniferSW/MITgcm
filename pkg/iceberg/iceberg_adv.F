#include "ICEBERG_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

CBOP
C !ROUTINE: ICEBERG_ADVDIFF

C !INTERFACE: ==========================================================
      SUBROUTINE ICEBERG_ADV(
     I                  myTime, myIter, myThid )

C !DESCRIPTION: \bv
C     *===========================================================*
C     | SUBROUTINE ICEBERG_ADVDIFF
C     | o driver for different advection routines
C     |   calls an adaption of gad_advection to call different
C     |   advection routines of pkg/generic_advdiff
C     | o the advection runs for every iceberg class separately in a loop
C     *===========================================================*
C \ev

C !USES: ===============================================================
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"
#include "DYNVARS.h"

C !INPUT/OUTPUT PARAMETERS: ===================================================
C     === Routine arguments ===
C     uc/vc     :: current ice velocity on C-grid;
C               :: C-Grid : Input only ; B-grid : Output only
C     myTime    :: current time in simulation
C     myIter    :: current iteration number in simulation
C     myThid    :: my Thread Id number
      _RL uc   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vc   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL myTime
      INTEGER myIter
      INTEGER myThid

C !LOCAL VARIABLES: ====================================================
C     === Local variables ===
C     i,j,bi,bj :: Loop counters
C     iCl        :: Loop counter for iceberg size classes
C     uTrans    :: volume transport, x direction
C     vTrans    :: volume transport, y direction
C     afx       :: horizontal advective flux, x direction
C     afy       :: horizontal advective flux, y direction
C     gFld      :: tendency of seaice field
C     xA,yA     :: "areas" of X and Y face of tracer cells
      INTEGER i, j, bi, bj, iCl, it
      INTEGER GAD_HEFF, GAD_AREA


      _RL uTrans    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL vTrans    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL afx       (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL afy       (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gFld      (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS xA        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS yA        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL recip_heff(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      CHARACTER*(18) diagLabel
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C TODO introduce GAD-Identifier for diagnostics
      GAD_HEFF = 1
      GAD_AREA = 2
C--   get surface velocities from MITgcm
C TODO calculate the velocity from momentum equation
      DO iCl = 1, ICEBERG_numClUsed
        CALL ICEBERG_MOM_EQ(uc,vc,iCl,myTime, myIter, myThid)

       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
c         DO j=1-OLy,sNy+OLy
c          DO i=1-OLx,sNx+OLx
c           uc(i,j,bi,bj)=uc(i,j,bi,bj)*10000000
c           vc(i,j,bi,bj)=vc(i,j,bi,bj)*10000000
c          ENDDO
c         ENDDO


C-    compute cell areas used by all tracers
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
c           xA(i,j,bi,bj) = 1. _d 0
c           yA(i,j,bi,bj) = 1. _d 0
           xA(i,j,bi,bj) = _dyG(i,j,bi,bj)*maskW(i,j,1,bi,bj)
           yA(i,j,bi,bj) = _dxG(i,j,bi,bj)*maskS(i,j,1,bi,bj)
          ENDDO
         ENDDO
C         End bi, bj loops
        ENDDO
       ENDDO




       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
C---   loops on tile indices bi,bj

#ifdef ALLOW_AUTODIFF_TAMC
C     Initialise for TAF
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           gFld(i,j)       = 0. _d 0
          ENDDO
         ENDDO
C
#endif /* ALLOW_AUTODIFF_TAMC */
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           recip_heff(i,j)    = 1. _d 0
          ENDDO
         ENDDO

C-    Calculate "volume transports" through tracer cell faces.

         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           uTrans(i,j) = uc(i,j,bi,bj)*xA(i,j,bi,bj)
           vTrans(i,j) = vc(i,j,bi,bj)*yA(i,j,bi,bj)
          ENDDO
         ENDDO
C C TODO customize ICEBERG_ADVECTION
C HINT GAD_HEFF is a tracer identifier to customize the diagnostics; at the moment it is not used for icebergs but
C       it coul be used to save the tendency for icebergs and fractional area

C ICEBERG_ADVECTION(tracerIdentity, advectionScheme, uFld, vFl, uTrans, vTrans
C                   iceFld, r_hFld, gFld, afx, afy, bi, bj, myTime, myIter, myThid)
C uFld, vFld - input velocity to use for iceberg
C uTrans, vTrans - input volume transport through grid boundaries
C iceFld - input iceberg distribution field
C r_hFld - ????? input reciprocal of ice thickness
C gFld - output tendency of icebergs used to do the advection
C afx, afy - horizontal advective flux

         CALL ICEBERG_ADVECTION(
     I        GAD_HEFF, ICEBERGadvScheme,
     I        uc(1-OLx,1-OLy,bi,bj), vc(1-OLx,1-OLy,bi,bj),
     I        uTrans, vTrans, iceberg_distr(1-OLx,1-OLy,bi,bj,iCl),
     I        recip_heff,
     O        gFld, afx, afy,
     I        bi, bj, myTime, myIter, myThid )
C     now do the "explicit" time step
         DO j=1,sNy
          DO i=1,sNx
              iceberg_distr(i,j,bi,bj,iCl) =
     &            iceberg_distr(i,j,bi,bj,iCl) + ICEBERG_deltaTtherm
     &            * gFld(i,j)
          ENDDO
         ENDDO


C fractional area



         CALL ICEBERG_ADVECTION(
     I        GAD_AREA, ICEBERGadvScheme,
     I        uc(1-OLx,1-OLy,bi,bj), vc(1-OLx,1-OLy,bi,bj),
     I        uTrans, vTrans, iceberg_area(1-OLx,1-OLy,bi,bj,iCl),
     I        recip_heff,
     O        gFld, afx, afy,
     I        bi, bj, myTime, myIter, myThid )

         DO j=1,sNy
          DO i=1,sNx
              iceberg_area(i,j,bi,bj,iCl) =
     &            iceberg_area(i,j,bi,bj,iCl) + ICEBERG_deltaTtherm
     &            * gFld(i,j)
           ENDDO
          ENDDO


C---   end bi,bj loops
        ENDDO
       ENDDO


C HINT Maybe diagnostics has to be moved into bi, bj loop



           WRITE(diagLabel,'(A,A)') 'IB_distr_',ICEBERG_ClLbl(iCl)
           CALL DIAGNOSTICS_FILL(iceberg_distr(1-Olx,1-Oly,1,1,iCl),
     &                            ICEBERG_ClLbl(iCl),0,Nr, 0, 1, 1,
     &                            myThid )



           WRITE(diagLabel,'(A,A)') 'IB_area_',ICEBERG_ClLbl(iCl)
           CALL DIAGNOSTICS_FILL(iceberg_area(1-Olx,1-Oly,1,1,iCl),
     &                            ICEBERG_ClLbl(iCl),0,Nr, 0, 1, 1,
     &                            myThid )

      ENDDO
      RETURN
      END

