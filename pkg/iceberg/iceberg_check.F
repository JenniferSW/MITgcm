#include "ICEBERG_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: ICEBERG_CHECK

C !INTERFACE:
      SUBROUTINE ICEBERG_CHECK( myThid )

C     !DESCRIPTION:
C     Check dependances with other packages

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"


C     !INPUT PARAMETERS:
C     myThid   :: my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_ICEBERG
C     !LOCAL VARIABLES:
C     msgBuf   :: Informational/error message buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      _BEGIN_MASTER(myThid)

       WRITE(msgBuf,'(A)') 'ICEBERG_CHECK: #define ICEBERG'
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT , 1)

C--  Print out some kee parameters :
       CALL WRITE_0D_L( ICEBERG_numClUsed, INDEX_NONE,
     &  'ICEBERG_numClUsed =', ' /* number of used size classes */')
       CALL WRITE_0D_L( ICEBERG_Iter0, INDEX_NONE,
     &  'ICEBERG_Iter0 =', ' /* Start time of iceberg simulation */')

C-   add some other as example:
c       CALL WRITE_0D_RL( myPa_param1, INDEX_NONE, ' myPa_param1 =',
c     &  ' /* user defined parameter */')


C--  Check parameters:

C-     Put stops here if some flags are wrongly chosen.
C-     For example, require this package to use tempStepping
c       IF ( .NOT.tempStepping .AND. myPa_applyTendT ) THEN
c        WRITE(msgBuf,'(A)') 'An example requirement is:'
c        CALL PRINT_ERROR( msgBuf , 1)
c        WRITE(msgBuf,'(A)') 'myPa_applyTendT needs tempStepping=.true.'
c        CALL PRINT_ERROR( msgBuf , 1)
c        CALL ALL_PROC_DIE( 0 )
c        STOP 'ABNORMAL END: S/R ICEBERG_CHECK'
c       ENDIF

      _END_MASTER(myThid)

#endif /* ALLOW_ICEBERG */

      RETURN
      END
