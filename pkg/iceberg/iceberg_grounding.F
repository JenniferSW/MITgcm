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
      INTEGER i, j, bi, bj
      INTEGER nCl

      PRINT *, 'dbg grounding'
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=0 ,sNy
         DO i= 0,sNx
           IF (ICEBERG_height(nCl) .GT. ABS(R_low(i,j,bi,bj)).OR.
     &     ICEBERG_height(nCl) .GT. ABS(R_low(i-1,j,bi,bj))) THEN
            iceberg_uvel(i,j,bi,bj,nCl) = 0. _d 0
          ENDIF
          IF (ICEBERG_height(nCl) .GT. ABS(R_low(i,j,bi,bj)).OR.
     &     ICEBERG_height(nCl) .GT. ABS(R_low(i,j-1,bi,bj))) THEN
            iceberg_vvel(i,j,bi,bj,nCl) = 0. _d 0
          ENDIF
          CALL EXCH_UV_XY_RL(iceberg_uvel(1-OLx,1-OLy,1,1,nCl),
     &    iceberg_vvel(1-OLx,1-OLy,1,1,nCl),
     &    .TRUE.,myThid)
         ENDDO ! end i loop
        ENDDO ! end j loop
       ENDDO ! end bi loop
      ENDDO ! end bj loop
      RETURN
      END
