#include "ICEBERG_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: ICEBERG_DIAGNOSTICS_STATE

C !INTERFACE:
      SUBROUTINE ICEBERG_DIAGNOSTICS_STATE(
     I                     myTime, myIter, myThid )

C     !DESCRIPTION:
C     Fill-in diagnostics array for ICEBERG state variables

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
C  put following include within #ifdef in case this file is copied to local
C  code dir but pkg/iceberg is not compiled (missing file ICEBERG.h).
#include "ICEBERG_SIZE.h"
#ifdef ALLOW_ICEBERG
#include "ICEBERG.h"
#endif


C     !INPUT PARAMETERS:
C     myTime   :: Current time of simulation ( s )
C     myIter   :: Current iteration number in simulation
C     myThid   :: my Thread Id number
      _RL     myTime
      INTEGER myIter, myThid

      CHARACTER*(8) diagName

CEOP

#ifdef ALLOW_ICEBERG
C     !LOCAL VARIABLES:
c     iCl  ::  loop parameter for size classes

      INTEGER iCl


#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN


c     Loop over iceberg size classes and write them to diagnostics

       DO iCl = 1, ICEBERG_numClUsed
        WRITE(diagName,'(A8)') 'IBD'//ICEBERG_ClLbl(iCl)
        CALL DIAGNOSTICS_FILL( iceberg_distr(1-Olx,1-Oly,1,1,iCl),
     &                         diagName,0,Nr, 0, 1, 1,
     &                         myThid )
       ENDDO

       DO iCl = 1, ICEBERG_numClUsed
        WRITE(diagName,'(A8)') 'IBA'//ICEBERG_ClLbl(iCl)
        CALL DIAGNOSTICS_FILL( iceberg_area(1-Olx,1-Oly,1,1,iCl),
     &                         diagName,0,Nr, 0, 1, 1,
     &                         myThid )
       ENDDO

       DO iCl = 1, ICEBERG_numClUsed
        WRITE(diagName,'(A8)') 'IBH'//ICEBERG_ClLbl(iCl)
        CALL DIAGNOSTICS_FILL( iceberg_height(1-Olx,1-Oly,1,1,iCl),
     &                         diagName,0,Nr, 0, 1, 1,
     &                         myThid )
       ENDDO

      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_ICEBERG */

      RETURN
      END
