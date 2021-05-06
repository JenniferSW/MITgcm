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
C     | this routine is based on the advection of seaice in the pkg seaice
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
C               :: C-Grid : Input only ; B-grid : Output only
C     myTime    :: current time in simulation
C     myIter    :: current iteration number in simulation
C     myThid    :: my Thread Id number
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
      INTEGER i, j, bi, bj, iCl
      INTEGER GAD_HEFF, GAD_AREA
      _RL tmpscal, recip_nCl
      _RL uTrans    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL vTrans    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL afx       (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL afy       (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gFld      (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS xA        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS yA        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL recip_heff(1-OLx:sNx+OLx,1-OLy:sNy+OLy)


      CHARACTER*(8) diagName
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C TODO introduce GAD-Identifier for diagnostics
      GAD_HEFF = 1
      GAD_AREA = 2

      DO iCl = 1, ICEBERG_numClUsed
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)

C-    compute cell areas used by all tracers
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
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
           uTrans(i,j) = iceberg_uvel(i,j,bi,bj,iCl)*xA(i,j,bi,bj)
           vTrans(i,j) = iceberg_vvel(i,j,bi,bj,iCl)*yA(i,j,bi,bj)
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
C r_hFld - input reciprocal of ice thickness
C gFld - output tendency of icebergs used to do the advection
C afx, afy - horizontal advective flux

C- advection of concentration
         CALL ICEBERG_ADVECTION(
     I        GAD_HEFF, ICEBERGadvScheme,
     I        iceberg_uvel(1-OLx,1-OLy,bi,bj,iCl),
     I        iceberg_vvel(1-OLx,1-OLy,bi,bj,iCl),
     I        uTrans, vTrans,
     I        iceberg_concentration(1-OLx,1-OLy,bi,bj,iCl),
     I        recip_heff,
     O        gFld, afx, afy,
     I        bi, bj, myTime, myIter, myThid )
C     now do the "explicit" time step
         DO j=1,sNy
          DO i=1,sNx
              iceberg_concentration(i,j,bi,bj,iCl) =
     &            maskC(i,j,1,bi,bj) * (
     &            iceberg_concentration(i,j,bi,bj,iCl)
     &            + ICEBERG_deltaT * gFld(i,j))
          ENDDO
         ENDDO


C_________________________________________________________________________________________________
C Fix Fields after Advection
C ------------------------------------------------------------------------------------------------
C- after advection fix possible negative values of effective iceberg thickness and fractional area
C followed the example of seaice_reg_ridge.F
#ifdef ALLOW_FIXING
          DO i = 1,sNy
            DO j = 1, sNx
             tmpscal = 0. _d 0
             tmpscal =
     &        MAX(-iceberg_concentration(i,j,bi,bj,iCl),0. _d 0)
             iceberg_concentration(i,j,bi,bj,iCl) =
     &          iceberg_concentration(i,j,bi,bj,iCl) + tmpscal
            ENDDO
          ENDDO

#endif

         CALL EXCH_3D_RL( iceberg_concentration,  iCl, myThid )

C---   end iCl bi,bj loops
        ENDDO !bi
       ENDDO !bj
      ENDDO !iCl
      RETURN
      END
