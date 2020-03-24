#include "ICEBERG_OPTIONS.h"

CBOP
      SUBROUTINE ICEBERG_OUTPUT( myTime, myIter, myThid)

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "ICEBERG.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myTime :: Current time of simulation ( s )
C     myIter :: Iteration number
C     myThid :: my Thread Id number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP


#ifdef ALLOW_TIMEAVE
      CHARACTER*(MAX_LEN_MBUF) fn
      INTEGER bi, bj, kl
#endif
#ifdef ALLOW_MNC
      INTEGER i
      CHARACTER*(MAX_LEN_MBUF) mncf
      CHARACTER*(1) pf
#endif

#ifdef ALLOW_TIMEAVE

c make a control output for seeing if the modelis running
c Write a testfile with some text for every time step

      WRITE(fn,'(A5)') 'Done '

#endif /*ALLOW_TIMEAVE*/


      RETURN
      END
