
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

      CHARACTER*(MAX_LEN_FNAM) tmpInitialFile
CEOP

#ifdef ALLOW_ICEBERG

C     !LOCAL VARIABLES:
C     === Local variables ===
      INTEGER i,j,k,bi,bj,iCl

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

C     Initialize arrays in common blocks :
C     3-dim.
        DO k=1,Nr
         DO j=1-Oly,sNy+OLy
          DO i=1-Olx,sNx+Olx
#ifdef ICEBERG_3D_STATE
           DO iCl = 1,ICEBERG_numClUsed
            iceberg_distr(i,j,k,bi,bj,iCl) = 0. _d 0
           ENDDO
           myPa_StatScal2(i,j,k,bi,bj) = 0. _d 0
           myPa_StatVelU (i,j,k,bi,bj) = 0. _d 0
           myPa_StatVelV (i,j,k,bi,bj) = 0. _d 0
#endif /* ICEBERG_3D_STATE */
#ifdef ICEBERG_TENDENCY
           myPa_TendScal1(i,j,k,bi,bj) = 0. _d 0
           myPa_TendScal2(i,j,k,bi,bj) = 0. _d 0
           myPa_TendVelU (i,j,k,bi,bj) = 0. _d 0
           myPa_TendVelV (i,j,k,bi,bj) = 0. _d 0
#endif /* ICEBERG_TENDENCY */
          ENDDO
         ENDDO
        ENDDO

C     2-dim.
#ifdef ICEBERG_2D_STATE
        DO j=1-Oly,sNy+OLy
         DO i=1-Olx,sNx+Olx
          myPa_Surf1(i,j,bi,bj) = 0. _d 0
          myPa_Surf2(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
#endif /* ICEBERG_2D_STATE */

C- end bi,bj loops
       ENDDO
      ENDDO

      IF ( startTime.EQ.baseTime .AND. nIter0.EQ.0
     &                           .AND. pickupSuff.EQ.' ' ) THEN

C     Prevent loading of new data before everyone has finished initialisation
C     this is done within IO routines => no longer needed
c      _BARRIER

#ifdef ICEBERG_3D_STATE
c      IF (nIter0 .EQ.ICEBERG_Iter0) THEN
       DO iCl = 1, ICEBERG_numClUsed
        tmpInitialFile = ICEBERG_initialFile(iCl)
        IF ( tmpInitialFile .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( tmpInitialFile, ' ',
     &                         iceberg_distr(1-OLx,1-OLy,1,1,1,iCl), 0,
     &                         myThid )
         CALL EXCH_3D_RL( iceberg_distr(1-OLx,1-OLy,1,1,1,iCl),
     &                         Nr, myThid )
        ENDIF
       ENDDO
c      ENDIF
       IF ( myPa_Scal2File .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( myPa_Scal2File, ' ',
     &                         myPa_StatScal2, 0, myThid )
         CALL EXCH_3D_RL( myPa_StatScal2, Nr, myThid )
       ENDIF

       IF ( myPa_VelUFile .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( myPa_VelUFile, ' ',
     &                         myPa_StatVelU, 0, myThid )
       ENDIF
       IF ( myPa_VelVFile .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( myPa_VelVFile, ' ',
     &                         myPa_StatVelV, 0, myThid )
       ENDIF
       IF ( myPa_VelUFile .NE. ' ' .OR.
     &      myPa_VelVFile .NE. ' ' ) THEN
        IF ( myPa_StaV_Cgrid ) THEN
         CALL EXCH_UV_3D_RL( myPa_StatVelU, myPa_StatVelV,
     &                       .TRUE., Nr, myThid )
        ELSE
C-      Assume Agrid position:
         CALL EXCH_UV_AGRID_3D_RL( myPa_StatVelU, myPa_StatVelV,
     &                       .TRUE., Nr, myThid )
        ENDIF
       ENDIF
#endif /* ICEBERG_3D_STATE */

#ifdef ICEBERG_2D_STATE
       IF ( myPa_Surf1File .NE. ' ' ) THEN
         CALL READ_FLD_XY_RL( myPa_Surf1File,' ',
     &                        myPa_Surf1, 0, myThid )
         CALL EXCH_XY_RL( myPa_Surf1, myThid )
       ENDIF
       IF ( myPa_Surf2File .NE. ' ' ) THEN
         CALL READ_FLD_XY_RL( myPa_Surf2File,' ',
     &                        myPa_Surf2, 0, myThid )
         CALL EXCH_XY_RL( myPa_Surf2, myThid )
       ENDIF
#endif /* ICEBERG_2D_STATE */

      ELSE
C--   restart from a pickup:
         CALL ICEBERG_READ_PICKUP( nIter0, myThid )

C--   end start-from-iter-zero if/else block
      ENDIF

#endif /* ALLOW_ICEBERG */

      RETURN
      END
