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
C HINT iceberg_MaskU can be evaluated only by the added area, because when it is .gt. 0, then there is an iceberg anyway
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
       ENDDO


C--   Calculate iceberg velocities
      DO iCl = 1, ICEBERG_numClUsed
        CALL ICEBERG_MOM_EQ(
     I              iCl,
     I              myTime, myIter, myThid)
      ENDDO


C--   Advect iceberg (effective height and fractional area)

      CALL ICEBERG_ADV(myTime, myIter, myThid)


C--   Melt iceberg concentration

      CALL ICEBERG_MELTING(myTime,myIter,myThid)

c            PRINT*, 'debug melt vloss',i,j, iceberg_Vloss(i,j,bi,bj,iCl)
C calculate new iceberg volume in grid cell
      DO iCl = 1, ICEBERG_numClUsed
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1-OLy+1,sNy+OLy
          DO i=1-OLx+1,sNx+OLx

            IF(iCl .EQ. 1 .AND. ICEBERG_numClUsed .GT. 1) THEN
            !Wenn ich in der kleinsten size class bin, dann kann nichts in die untere Klasse wandern
            ! gilt aber nur, wenn ich mind 2 size classes habe
            iceberg_concentration(i,j,bi,bj,iCl) =
     &         iceberg_concentration(i,j,bi,bj,iCl)
     &         - iceberg_Vloss(i,j,bi,bj,iCl)
     &         /(iceberg_Vsgl(iCl)-0) !dann ist kein Eisberg mehr da
     &         + iceberg_Vloss(i,j,bi,bj,iCl+1)
     &         /(iceberg_Vsgl(iCl+1)-iceberg_Vsgl(iCl))

            ELSE IF(iCl .EQ. ICEBERG_numClUsed
     &        .AND. ICEBERG_numClUsed .GT. 1) THEN
            !Wenn ich in der größten size class bin, dann kann nichts aus der höheren kommen
            ! geht nur, wenn ich mindestens 2 size classes habe
               iceberg_concentration(i,j,bi,bj,iCl) =
     &         iceberg_concentration(i,j,bi,bj,iCl)
     &         - iceberg_Vloss(i,j,bi,bj,iCl)
     &         /(iceberg_Vsgl(iCl)-iceberg_Vsgl(iCl-1))
            ELSE IF(ICEBERG_numClUsed .GT. 2) THEN
            ! ansonsten bin ich in den anderen Klassen und habe eine höhere und eine kleinere
            ! geht nur, wenn ich mind. 3 size classes habe, sonst greifen nur die oberen beiden Fälle
             iceberg_concentration(i,j,bi,bj,iCl) =
     &         iceberg_concentration(i,j,bi,bj,iCl)
     &         - iceberg_Vloss(i,j,bi,bj,iCl)
     &         /(iceberg_Vsgl(iCl)-iceberg_Vsgl(iCl-1))
     &         + iceberg_Vloss(i,j,bi,bj,iCl+1)
     &         /(iceberg_Vsgl(iCl+1)-iceberg_Vsgl(iCl))
            ELSE !ICEBERG_numClUsed is 1
            ! nichts kommt aus einer höheren size class
            ! es gibt keine kleinere size class
             iceberg_concentration(i,j,bi,bj,iCl) =
     &         iceberg_concentration(i,j,bi,bj,iCl)
     &         - iceberg_Vloss(i,j,bi,bj,iCl)
     &         /(iceberg_Vsgl(iCl)-0)
            ENDIF
C catch values that are too small or values that are negative because of too strong melting

            IF (iceberg_concentration(i,j,bi,bj,iCl).LE.1. _d -5) THEN
             iceberg_concentration(i,j,bi,bj,iCl) = 0. _d 0
            ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END
