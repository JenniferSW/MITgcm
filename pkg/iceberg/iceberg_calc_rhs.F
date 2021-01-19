#include "ICEBERG_OPTIONS.h"

CBOP
C !ROUTINE: ICEBERG_CALC_RHS

C !INTERFACE: ==========================================================
      SUBROUTINE ICEBERG_CALC_RHS( myTime, myIter, myThid )

C !DESCRIPTION:
C Calculate custom tendency terms outside k-loop in DO_OCEANIC_PHYS

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"


C !INPUT PARAMETERS: ===================================================
C myTime    :: Current time in simulation
C myIter    :: Current time-step number
C myThid    :: my Thread Id number
      _RL     myTime
      INTEGER myIter, myThid

C !OUTPUT PARAMETERS: ==================================================


C !LOCAL VARIABLES: ====================================================
C bi,bj     :: Tile indices
      INTEGER bi, bj, i, j, k, iCl
      _RL maskU     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL maskV     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL mask_C     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CEOP

C--   Masking of iceberg - Mask cells with icebergs (1) and without (0)
       DO iCl = 1, ICEBERG_numClUsed
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO j=1-OLy+1,sNy+OLy
           DO i=1-OLx+1,sNx+OLx
C HINT iceberg_MaskU can be evaluated only by the added area, because when it is .gt. 0, then there is an iceberg anyway
            iceberg_MaskU(i,j,bi,bj) =
     &         iceberg_concentration(i,j,bi,bj,iCl)
     &         + iceberg_concentration(i-1,j,bi,bj,iCl)

            IF ( ( iceberg_MaskU(i,j,bi,bj) .GT. 0. _d 0 ) THEN
             iceberg_MaskU(i,j,bi,bj) = 1. _d 0
            ELSE
             iceberg_MaskU(i,j,bi,bj) = 0. _d 0
            ENDIF

            iceberg_MaskV(i,j,bi,bj) =
     &         iceberg_concentration(i,j,bi,bj,iCl)
     &         + iceberg_concentration(i,j-1,bi,bj,iCl)

            IF ( (iceberg_MaskV(i,j,bi,bj) .GT. 0. _d 0) THEN
             iceberg_MaskV(i,j,bi,bj) = 1. _d 0
            ELSE
             iceberg_MaskV(i,j,bi,bj) = 0. _d 0
            ENDIF

            iceberg_MaskC(i,j,bi,bj) =
     &         iceberg_concentration(i,j,bi,bj,iCl)
            IF ( (iceberg_MaskC(i,j,bi,bj) .GT. 0. _d 0) .AND.
     &           (mask_C(i,j,bi,bj) .GT. 0. _d 0) ) THEN
             iceberg_MaskC(i,j,bi,bj) = 1. _d 0
            ELSE
             iceberg_MaskC(i,j,bi,bj) = 0. _d 0
            ENDIF

           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO


C--   Calculate iceberg velocities
      DO iCl = 1, ICEBERG_numClUsed
        CALL ICEBERG_MOM_EQ(
     I              iCl,
     I              myTime, myIter, myThid)
      ENDDO


C--   Advect iceberg (effective height and fractional area)

      CALL ICEBERG_ADV(myTime, myIter, myThid)


C--   Melt iceberg concentration (melting is done considering the mechanism
C     from Thomas Rackows pyberg but using the mass balance from Martin and Adcroft (2011)
C     MA10 used different definitions of impact areas for the melting mechanisms than pyberg

c      CALL ICEBERG_MELTING(myTime,myIter,myThid)


CML#ifdef ALLOW_DIAGNOSTICS
CML      IF ( useDiagnostics ) THEN
CML         DO iCl = 1, ICEBERG_numClUsed
CML           CALL DIAGNOSTICS_FILL( iceberg_distr(1-Olx,1-Oly,1,1,iCl),
CML     &                            ICEBERG_ClLbl(iCl),0,Nr, 0, 1, 1,
CML     &                            myThid )
CML         ENDDO
CML
CML
CML      ENDIF
CML#endif /*ALLOW_DIAGNOSTICS*/




      RETURN
      END
