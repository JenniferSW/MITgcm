#include "ICEBERG_OPTIONS.h"
C--  File iceberg_tendency_apply.F: Routines to apply ICEBERG tendencies
C--   Contents
C--   o ICEBERG_TENDENCY_APPLY_U
C--   o ICEBERG_TENDENCY_APPLY_V
C--   o ICEBERG_TENDENCY_APPLY_T
C--   o ICEBERG_TENDENCY_APPLY_S

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: ICEBERG_TENDENCY_APPLY_U

C !INTERFACE:
      SUBROUTINE ICEBERG_TENDENCY_APPLY_U(
     U                     gU_arr,
     I                     iMin,iMax,jMin,jMax, k, bi, bj,
     I                     myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add iceberg tendency terms to U tendency.  Routine works for one
C     level at a time. Assumes that U and V tendencies are on C-Grid or
C     on A-grid.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
c#include "DYNVARS.h"
#include "ICEBERG.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gU_arr    :: the tendency array
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gU_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ICEBERG_TENDENCY
C     !LOCAL VARIABLES:
      INTEGER i, j

      IF ( myPa_applyTendU .AND. myPa_Tend_Cgrid ) THEN
        DO j=jMin,jMax
         DO i=iMin,iMax
          gU_arr(i,j) = gU_arr(i,j)
     &     + maskW(i,j,k,bi,bj) * myPa_TendVelU(i,j,k,bi,bj)
         ENDDO
        ENDDO
      ELSEIF ( myPa_applyTendU ) THEN
C--   Assume Agrid position:
        DO j=jMin,jMax
         DO i=iMin,iMax
          gU_arr(i,j) = gU_arr(i,j)
     &     + maskW(i,j,k,bi,bj)*( myPa_TendVelU(i-1,j,k,bi,bj)
     &                           +myPa_TendVelU( i ,j,k,bi,bj)
     &                          )*0.5 _d 0
         ENDDO
        ENDDO
      ENDIF

#endif /* ICEBERG_TENDENCY */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: ICEBERG_TENDENCY_APPLY_V

C !INTERFACE:
      SUBROUTINE ICEBERG_TENDENCY_APPLY_V(
     U                     gV_arr,
     I                     iMin,iMax,jMin,jMax, k, bi, bj,
     I                     myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add iceberg tendency terms to V tendency.  Routine works for one
C     level at a time. Assumes that U and V tendencies are on C-Grid or
C     on A-grid.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
c#include "DYNVARS.h"
#include "ICEBERG.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gV_arr    :: the tendency array
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gV_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ICEBERG_TENDENCY
C     !LOCAL VARIABLES:
      INTEGER i, j

      IF ( myPa_applyTendV .AND. myPa_Tend_Cgrid ) THEN
        DO j=jMin,jMax
         DO i=iMin,iMax
          gV_arr(i,j) = gV_arr(i,j)
     &     + maskS(i,j,k,bi,bj) * myPa_TendVelV(i,j,k,bi,bj)
         ENDDO
        ENDDO
      ELSEIF ( myPa_applyTendV ) THEN
C--   Assume Agrid position:
        DO j=jMin,jMax
         DO i=iMin,iMax
          gV_arr(i,j) = gV_arr(i,j)
     &     + maskS(i,j,k,bi,bj)*( myPa_TendVelV(i,j-1,k,bi,bj)
     &                           +myPa_TendVelV(i, j ,k,bi,bj)
     &                          )*0.5 _d 0
         ENDDO
        ENDDO
      ENDIF

#endif /* ICEBERG_TENDENCY */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: ICEBERG_TENDENCY_APPLY_T

C !INTERFACE:
      SUBROUTINE ICEBERG_TENDENCY_APPLY_T(
     U                     gT_arr,
     I                     iMin,iMax,jMin,jMax, k, bi, bj,
     I                     myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add iceberg tendency terms to T (theta) tendency.  Routine works
C     for one level at a time

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
c#include "DYNVARS.h"
#include "ICEBERG.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gT_arr    :: the tendency array
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ICEBERG_TENDENCY
C     !LOCAL VARIABLES:
      INTEGER i, j

      IF ( myPa_applyTendT ) THEN
        DO j=jMin,jMax
         DO i=iMin,iMax
          gT_arr(i,j) = gT_arr(i,j)
     &     + maskC(i,j,k,bi,bj) * myPa_TendScal1(i,j,k,bi,bj)
         ENDDO
        ENDDO
      ENDIF

#endif /* ICEBERG_TENDENCY */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: ICEBERG_TENDENCY_APPLY_S

C !INTERFACE:
      SUBROUTINE ICEBERG_TENDENCY_APPLY_S(
     U                     gS_arr,
     I                     iMin,iMax,jMin,jMax, k, bi, bj,
     I                     myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add iceberg tendency terms to S tendency.  Routine works for one
C     level at a time

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
c#include "DYNVARS.h"
#include "ICEBERG.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gS_arr    :: the tendency array
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gS_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ICEBERG_TENDENCY
C     !LOCAL VARIABLES:
      INTEGER i, j

      IF ( myPa_applyTendS ) THEN
        DO j=jMin,jMax
         DO i=iMin,iMax
          gS_arr(i,j) = gS_arr(i,j)
     &     + maskC(i,j,k,bi,bj) * myPa_TendScal2(i,j,k,bi,bj)
         ENDDO
        ENDDO
      ENDIF

#endif /* ICEBERG_TENDENCY */

      RETURN
      END
