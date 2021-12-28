#include "ICEBERG_OPTIONS.h"
C     Without ALLOW_SEAICE, this would never work anyway.
C     Needs to be checked in iceberg_check.F
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

CBOP
C !ROUTINE: ICEBERG_ADVDIFF

C !INTERFACE: ==========================================================
      SUBROUTINE ICEBERG_CHANGE_SALINITY( nCl,
     I                  myTime, myIter, myThid)

C !DESCRIPTION: \bv
C     *===========================================================*
C     | SUBROUTINE ICEBERG_CHANGE_SALINITY
C     | o in this routine the salinity is changed due to meltwater input
C     |
C     *===========================================================*
C     | written by Jennifer S. Wesselbaum, https://github.com/JenniferSW
C     *===========================================================*
C
C !USES: ===============================================================
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SURFACE.h"
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"
#include "DYNVARS.h"

C !INPUT/OUTPUT PARAMETERS: ===================================================
C     === Routine arguments ===
C     myTime    :: current time in simulation
C     myIter    :: current iteration number in simulation
C     myThid    :: my Thread Id number
C     nCl       :: number of class the velocity is calculated for
      _RL myTime
      INTEGER myIter
      INTEGER myThid
      INTEGER nCl

C !LOCAL VARIABLES: ====================================================
C     === Local variables ===
C     i,j,bi,bj :: Loop counters
      INTEGER i, j, bi, bj, k
      _RL Vmelt
      _RL FW
      _RL VFW
      _RL FW_melt
      _RL cell_volume


CEOP
       k = 1
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1,sNy
          DO i=1,sNx
           IF(ICEBERG_MaskC(i,j,bi,bj) .NE. 0. _d 0)THEN
C Evaluate Vmelt for this grid cell
            Vmelt = 0. _d 0
c            DO nCl = 1, ICEBERG_numClUsed
              Vmelt = Vmelt + iceberg_Vloss(i,j,bi,bj,nCl)
c            ENDDO
            cell_volume = rA(i,j,bi,bj)*drF(k)*hFacC(i,j,k,bi,bj)
C actual freshwater fraction
            FW = 1 - 1e-3 * salt(i,j,1,bi,bj)
C actual freshwater volume
            VFW = cell_volume * FW
C new freshwater fraction
            FW_melt = (VFW + Vmelt)/(cell_volume + Vmelt)
C new salinity
            salt(i,j,1,bi,bj) = (1 - FW_melt) *1e3
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       _EXCH_XYZ_RL( salt, myThid )
      RETURN
      END
