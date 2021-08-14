#include "ICEBERG_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

CBOP
C !ROUTINE: ICEBERG_GROUNDING

C !INTERFACE: ==========================================================
      SUBROUTINE ICEBERG_GROUNDING(nCl,
     I                  myTime, myIter, myThid)

C !DESCRIPTION: \bv
C     *===========================================================*
C     | SUBROUTINE ICEBERG_GROUNDING
C     | o checks if water is too shallow for iceberg drift
C     |   if iceberg hieght is larger than ocean depth the iceberg
C     |   velocity in this cell is set to zero
C     |   icebergs with larger height will stay in this cell until
C     |   they melted enough mass to get into a smaller size class.
C     |   If height of smaller size class is lower than ater depth
C     |   the iceberg will start drifting again.
C     *===========================================================*




C !USES: ===============================================================
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SURFACE.h"
#include "FFIELDS.h"
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"
#include "DYNVARS.h"
#ifdef ALLOW_EXF
# include "EXF_FIELDS.h"
# include "EXF_PARAM.h"
#endif
#ifdef ALLOW_SEAICE
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
#include "SEAICE.h"
#endif

C !INPUT/OUTPUT PARAMETERS: ===================================================
C     === Routine arguments ===
C     myTime    :: current time in simulation
C     myIter    :: current iteration number in simulation
C     myThid    :: my Thread Id number
      _RL myTime
      INTEGER myIter
      INTEGER myThid


C !LOCAL VARIABLES: ====================================================
C     === Local variables ===
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1 ,sNy
         DO i= 1,sNx
          IF (ICEBERG_height(nCl) .GT. R_low)
            iceberg_velocity(i,j,bi,bj,nCl) = 0. _d 0

         ENDDO ! end i loop
        ENDDO ! end j loop
       ENDDO ! end bi loop
      ENDDO ! end bj loop
      RETURN
      END
