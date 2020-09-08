
#include "ICEBERG_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: ICEBERG_INIT_VARIA

C !INTERFACE:
      SUBROUTINE ICEBERG_INIT_VARIA( myThid )

C     !DESCRIPTION:
C     Initialize ICEBERG variables

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid


      CHARACTER*(MAX_LEN_FNAM) msgBuf
CEOP

#ifdef ALLOW_ICEBERG

C     !LOCAL VARIABLES:
C     === Local variables ===
      CHARACTER*(MAX_LEN_FNAM) tmpInitialArea
      CHARACTER*(MAX_LEN_FNAM) tmpInitialHeight


      INTEGER i,j,k,bi,bj,iCl



      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

C     Initialize arrays in common blocks :
C     2-dim.

         DO j=1-Oly,sNy+OLy
          DO i=1-Olx,sNx+Olx
           DO iCl = 1,ICEBERG_numClUsed
            iceberg_distr(i,j,bi,bj,iCl) = 0. _d 0
            iceberg_area(i,j,bi,bj,iCl) = 0. _d 0
            iceberg_height(i,j,bi,bj,iCl) = 0. _d 0
           ENDDO
          ENDDO
         ENDDO


C- end bi,bj loops
       ENDDO
      ENDDO

      IF ( startTime.EQ.baseTime .AND. nIter0.EQ.0
     &                           .AND. pickupSuff.EQ.' ' ) THEN

C     Prevent loading of new data before everyone has finished initialisation
C     this is done within IO routines => no longer needed
c      _BARRIER


      IF (nIter0 .EQ. ICEBERG_Iter0) THEN
       DO iCl = 1, ICEBERG_numClUsed

        tmpInitialArea = ICEBERG_initialArea(iCl)
        tmpInitialHeight = ICEBERG_initialHeight(iCl)

        IF ( tmpInitialArea.EQ.' ' ) THEN
            WRITE(msgBuf,'(A)') 'Not enough input files for ICEBERG
     &      AREA. Make sure every iceberg class has one input
     &      file to be read in data.iceberg'
            CALL PRINT_ERROR(msgBuf,myThid)
            STOP 'ABNORMAL END: ICEBERG_INIT_VARIA'
        ENDIF

        IF ( tmpInitialHeight.EQ.' ' ) THEN
            WRITE(msgBuf,'(A)') 'Not enough input files for ICEBERG
     &      HEIGHT. Make sure every iceberg class has one input
     &      file to be read in data.iceberg'
            CALL PRINT_ERROR(msgBuf,myThid)
            STOP 'ABNORMAL END: ICEBERG_INIT_VARIA'
        ENDIF


        CALL READ_FLD_XY_RL( tmpInitialArea, ' ',
     &                       iceberg_area(1-OLx,1-OLy,1,1,iCl), 0,
     &                       myThid )
        CALL READ_FLD_XY_RL( tmpInitialHeight, ' ',
     &                       iceberg_height(1-OLx,1-OLy,1,1,iCl), 0,
     &                       myThid )


        _EXCH_XY_RL(iceberg_area , myThid)
        _EXCH_XY_RL(iceberg_area , myThid)

       ENDDO
      ENDIF

C Calculate effective iceberg thickness from fraction area and height

      DO iCl = 1,ICEBERG_numClUsed
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
          DO j=1-Oly,sNy+OLy
           DO i=1-Olx,sNx+Olx
            iceberg_distr(i,j,bi,bj,iCl) =
     &      iceberg_height(i,j,bi,bj,iCl) * iceberg_area(i,j,bi,bj,iCl)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO

C Initialize iceberg mass fields and masks


       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1-Oly,sNy+OLy
          DO i=1-Olx,sNx+Olx
            iceberg_MassC(i,j,bi,bj) = 0. _d 0
            iceberg_MassU(i,j,bi,bj) = 0. _d 0
            iceberg_MassV(i,j,bi,bj) = 0. _d 0

            iceberg_MaskC(i,j,bi,bj) = 0. _d 0
            iceberg_MaskU(i,j,bi,bj) = 0. _d 0
            iceberg_MaskV(i,j,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO
       ENDDO




      ELSE
C--   restart from a pickup:
         CALL ICEBERG_READ_PICKUP( nIter0, myThid )

C--   end start-from-iter-zero if/else block
      ENDIF

#endif /* ALLOW_ICEBERG */

      RETURN
      END
