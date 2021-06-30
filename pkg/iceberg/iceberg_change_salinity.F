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
      INTEGER i, j, bi, bj, k, l
      _RL Vmelt
      _RL FW
      _RL VFW
      _RL FW_melt
      _RL cell_volume
      _RL iceberg_maxlayer


CEOP

C Initialise local variables
      Vmelt = 0. _d 0
      FW = 0. _d 0
      VFW = 0. _d 0
      FW_melt = 0. _d 0
      cell_volume = 0. _d 0
      iceberg_maxlayer = 0. _d 0


       k = 1
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1,sNy
          DO i=1,sNx
           IF(usingPCoords)THEN
             Vmelt = 0. _d 0
             IF(maskC(i,j,k,bi,bj) .EQ. 1. _d 0)THEN
C Evaluate Vmelt for this grid cell
              Vmelt = iceberg_Vloss(i,j,bi,bj,nCl)

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


           ELSEIF (usingZCoords)THEN
C get Vmelt for actual size class in this grid cell
             Vmelt = iceberg_Vloss(i,j,bi,bj,nCl)
C get iceberg_maxlayer
             DO l =1,Nr
               IF(ABS(rF(l)) .GT. iceberg_height(nCl))THEN
C check if the depth at cell bottom is higher  than hieght of iceberg. Then the actual layer is the maximum layer
C the iceberg reaches into. The index of the bottom boundary of the layer is one higher than the layer index
                 iceberg_maxlayer = l-1
                 EXIT
               ENDIF
             ENDDO

C calculate freshwater volume per layer as follows
C  - get cell volume
C  - get actual freshwater volume in this cell
C  - calculate new freshwater fraction in this layer
C  - calculate new salinity in certain layer
C loop over all layers that are "filled" completely by iceberg height
             DO l = 1, iceberg_maxlayer -1
               cell_volume = rA(i,j,bi,bj)*drF(l)*hFacC(i,j,l,bi,bj)
               FW = 1 - 1e-3 * salt(i,j,l,bi,bj)
               VFW = cell_volume * FW
               Vmelt_cell = Vmelt * drF(l)/iceberg_height(nCl)
               FW_melt = (VFW + Vmelt_cell)/(cell_volume + Vmelt_cell)
               salt(i,j,l,bi,bj) = (1 - FW_melt) *1e3
             ENDDO
C last layer only partly covered by iceberg
             cell_volume = rA(i,j,bi,bj)*drF(iceberg_maxlayer)*hFacC(i,j,iceberg_maxlayer,bi,bj)
             FW = 1 - 1e-3 * salt(i,j,iceberg_maxlayer,bi,bj)
             VFW = cell_volume * FW
             Vmelt_cell = Vmelt * (iceberg_height(nCl) - rF(iceberg_maxlayer))/iceberg_height(nCl)
             FW_melt = (VFW + Vmelt_cell)/(cell_volume + Vmelt_cell)





           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       _EXCH_XYZ_RL( salt, myThid )
      RETURN
      END
