#include "ICEBERG_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: ICEBERG_INIT_FIXED

C !INTERFACE:
      SUBROUTINE ICEBERG_INIT_FIXED( myThid )

C     !DESCRIPTION:
C     Initialize ICEBERG variables that are kept fixed during the run.

C     !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"


C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     === Local variables ===




#ifdef ALLOW_MNC
      IF (useMNC) THEN
        CALL ICEBERG_MNC_INIT( myThid )
      ENDIF
#endif /* ALLOW_MNC */

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL ICEBERG_DIAGNOSTICS_INIT( myThid )
      ENDIF
#endif

      RETURN
      END
