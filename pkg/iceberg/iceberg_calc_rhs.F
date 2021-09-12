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
#include "GRID.h"

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
      _RL maskU     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL maskV     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL mask_C     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CEOP

C--   Masking of iceberg - Mask cells with icebergs (1) and without (0)
       DO iCl = 1, ICEBERG_numClUsed

        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO j=1-OLy+1,sNy+OLy
           DO i=1-OLx+1,sNx+OLx

            iceberg_MaskU(i,j,bi,bj) =
     &         iceberg_concentration(i,j,bi,bj,iCl)
     &         + iceberg_concentration(i-1,j,bi,bj,iCl)

            IF ( iceberg_MaskU(i,j,bi,bj) .GT. 0. _d 0 ) THEN
             iceberg_MaskU(i,j,bi,bj) = 1. _d 0
            ELSE
             iceberg_MaskU(i,j,bi,bj) = 0. _d 0
            ENDIF

            iceberg_MaskV(i,j,bi,bj) =
     &         iceberg_concentration(i,j,bi,bj,iCl)
     &         + iceberg_concentration(i,j-1,bi,bj,iCl)

            IF ( iceberg_MaskV(i,j,bi,bj) .GT. 0. _d 0) THEN
             iceberg_MaskV(i,j,bi,bj) = 1. _d 0
            ELSE
             iceberg_MaskV(i,j,bi,bj) = 0. _d 0
            ENDIF

            iceberg_MaskC(i,j,bi,bj) =
     &         iceberg_concentration(i,j,bi,bj,iCl)
            IF ( iceberg_MaskC(i,j,bi,bj) .GT. 0. _d 0 ) THEN
             iceberg_MaskC(i,j,bi,bj) = 1. _d 0
            ELSE
             iceberg_MaskC(i,j,bi,bj) = 0. _d 0
            ENDIF

           ENDDO
          ENDDO
         ENDDO
        ENDDO

#ifdef ALLOW_ICEBERG_DYNAMICS
C--   Calculate iceberg velocities
        CALL ICEBERG_MOM_EQ(
     I              iCl,
     I              myTime, myIter, myThid)

#endif
C--   Melt iceberg concentration
#ifdef ALLOW_ICEBERG_MELTING
        CALL ICEBERG_MELTING(myTime,myIter,myThid)
#endif

#ifdef ALLOW_ICEBERG_CALVING
        CALL ICEBERG_CALVING( iCl,
     I                  myTime, myIter, myThid)
#endif

#ifdef ALLOW_ICEBERG_DYNAMICS
C--   Advect iceberg (effective height and fractional area)
        CALL ICEBERG_ADV(myTime, myIter, myThid)
#endif
      ENDDO

C calculate new iceberg volume in grid cell after melting


      DO iCl = 1, ICEBERG_numClUsed
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1-OLy+1,sNy+OLy
          DO i=1-OLx+1,sNx+OLx

#ifdef ALLOW_ICEBERG_MELTING
            IF(iCl .EQ. 1 .AND. ICEBERG_numClUsed .GT. 1) THEN
            ! there is no transition from the smallest size class
            ! from the smallest size class the iceberg go to a non iceberg status
            ! This is only happening when there are minimum two size classes
            iceberg_concentration(i,j,bi,bj,iCl) =
     &         iceberg_concentration(i,j,bi,bj,iCl)
     &         - iceberg_Vloss(i,j,bi,bj,iCl)
     &         /(iceberg_Vsgl(iCl)-ICEBERG_minV) !dann ist kein Eisberg mehr da
     &         + iceberg_Vloss(i,j,bi,bj,iCl+1)
     &         /(iceberg_Vsgl(iCl+1)-iceberg_Vsgl(iCl))

            ELSE IF(iCl .EQ. ICEBERG_numClUsed
     &        .AND. ICEBERG_numClUsed .GT. 1) THEN
            ! There is no transition to the highest size class, so icebergs
            ! only go to the next smaller size class and new icebergs go into by calving
            ! This is only happening when there are minimum two size classes
               iceberg_concentration(i,j,bi,bj,iCl) =
     &         iceberg_concentration(i,j,bi,bj,iCl)
     &         - iceberg_Vloss(i,j,bi,bj,iCl)
     &         /(iceberg_Vsgl(iCl)-iceberg_Vsgl(iCl-1))
            ELSE IF(ICEBERG_numClUsed .GT. 2) THEN
            ! Icebergs go to the next lower size class and icebergs are added
            ! from next higher size class
            ! only works if there are more than two size classes
             iceberg_concentration(i,j,bi,bj,iCl) =
     &         iceberg_concentration(i,j,bi,bj,iCl)
     &         - iceberg_Vloss(i,j,bi,bj,iCl)
     &         /(iceberg_Vsgl(iCl)-iceberg_Vsgl(iCl-1))
     &         + iceberg_Vloss(i,j,bi,bj,iCl+1)
     &         /(iceberg_Vsgl(iCl+1)-iceberg_Vsgl(iCl))
            ELSE !ICEBERG_numClUsed is 1
            ! if there is only one size class there are no transitions to
            ! other size classes, only to non iceberg status
             iceberg_concentration(i,j,bi,bj,iCl) =
     &         iceberg_concentration(i,j,bi,bj,iCl)
     &         - iceberg_Vloss(i,j,bi,bj,iCl)
     &         /(iceberg_Vsgl(iCl)-0)
            ENDIF
C catch values that are too small or values that are negative because of too strong melting
#endif

#ifdef ALLOW_FIXING
            IF (iceberg_concentration(i,j,bi,bj,iCl).LE.1. _d -5) THEN
             iceberg_concentration(i,j,bi,bj,iCl) = 0. _d 0
            ENDIF
#endif

          ENDDO
         ENDDO
        ENDDO
       ENDDO

C Update calculated fields (iceberg concentration and iceberg volume loss
        CALL EXCH_XY_RL(iceberg_concentration(1-OLx,1-OLy,1,1,iCl),
     &   myThid)
        CALL EXCH_XY_RL(iceberg_Vloss(1-OLx,1-OLy,1,1,iCl),
     &   myThid)

#ifdef ALLOW_FRESHWATER_INPUT
       CALL ICEBERG_CHANGE_SALINITY( iCl,
     I     myTime, myIter, myThid)
#endif

      ENDDO !ENDDO over size classes

      RETURN
      END
