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
CEOP


      bi = 1
      bj = 1
      k = 1
c      i = 1
c      j = 1

C--   Calculate custom tendency terms ...

      CALL ICEBERG_ADV(myTime, myIter, myThid)

       
#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
         DO iCl = 1, ICEBERG_numClUsed
           CALL DIAGNOSTICS_FILL( iceberg_distr(1-Olx,1-Oly,1,1,iCl),
     &                            ICEBERG_ClLbl(iCl),0,Nr, 0, 1, 1,
     &                            myThid )
         ENDDO


      ENDIF
#endif /*ALLOW_DIAGNOSTICS*/

      


      RETURN
      END
