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
      INTEGER GAD_HEFF

c#ifdef SEAICE_ITD
c#endif /* SEAICE_ITD */
c#ifdef ALLOW_SITRACER
c#endif /* ALLOW_SITRACER */
      _RL uTrans    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL vTrans    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL afx       (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL afy       (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gFld      (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS xA        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS yA        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL recip_heff(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      GAD_HEFF = 1
C--   get surface velocities from MITgcm
C TODO calculate the velocity from momentum equation
C HINT at the moment the velocities are taken from the model as the surface velocities

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          uc(i,j,bi,bj)=uVel(i,j,1,bi,bj)*10000000
          vc(i,j,bi,bj)=vVel(i,j,1,bi,bj)*10000000
         ENDDO
        ENDDO

        PRINT *,'debug velocities u', uc
        PRINT*, 'debug velocities v', vc

C-    compute cell areas used by all tracers
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          xA(i,j,bi,bj) = 1. _d 0
          yA(i,j,bi,bj) = 1. _d 0


c          xA(i,j,bi,bj) = _dyG(i,j,bi,bj)*SIMaskU(i,j,bi,bj)
c          yA(i,j,bi,bj) = _dxG(i,j,bi,bj)*SIMaskV(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      IF ( SEAICEmultiDimAdvection ) THEN
C  TODO replace dummy variable also in ICEBERG.h

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
C TODO loop over size classes
        iCl = 1
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          uTrans(i,j) = uc(i,j,bi,bj)*xA(i,j,bi,bj)
          vTrans(i,j) = vc(i,j,bi,bj)*yA(i,j,bi,bj)
         ENDDO
        ENDDO
C C TODO customize ICEBERG_ADVECTION
C ICEBERG_ADVECTION(tracerIdentity, advectionScheme, uFld, vFl, uTrans, vTrans
C                   iceFld, r_hFld, gFld, afx, afy, bi, bj, myTime, myIter, myThid)
C uFld, vFld - input velocity to use for iceberg
C uTrans, vTrans - input volume transport through grid boundaries
C iceFld - input iceberg distribution field
C r_hFld - ????? input reciprocal of ice thickness
C gFld - output tendency of icebergs used to do the advection
C afx, afy - horizontal advective flux

C uc, vc filled with ocean surface velocity
C
C dummy variables - uTrans, vTrans, recip_heff, afx, afy



          CALL ICEBERG_ADVECTION(
     I         GAD_HEFF, ICEBERGadvScheme,
     I         uc(1-OLx,1-OLy,bi,bj), vc(1-OLx,1-OLy,bi,bj),
     I         uTrans, vTrans, iceberg_distr(1-OLx,1-OLy,1,bi,bj,iCl),
     I         recip_heff,
     O         gFld, afx, afy,
     I         bi, bj, myTime, myIter, myThid )
C     now do the "explicit" time step
          DO j=1,sNy
           DO i=1,sNx
            iceberg_distr(i,j,1,bi,bj,iCl) =
     &           iceberg_distr(i,j,1,bi,bj,iCl) + ICEBERG_deltaTtherm
     &           * gFld(i,j)

           ENDDO
          ENDDO

C        ENDIF


C---   end bi,bj loops
        ENDDO
       ENDDO
C HINT Maybe diagnostics has to be moved into bi, bj loop

       PRINT *, 'foobar'
       PRINT *, gFld(i,j)
       DO iCl = 1, ICEBERG_numClUsed
           CALL DIAGNOSTICS_FILL(iceberg_distr(1-Olx,1-Oly,1,1,1,iCl),
     &                            ICEBERG_ClLbl(iCl),0,Nr, 0, 1, 1,
     &                            myThid )
       ENDDO


C      ENDDO


C      ELSE
C       STOP 'undefined advection scheme'

C--   end if multiDimAdvection
      ENDIF

      RETURN
      END

