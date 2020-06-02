
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
      CHARACTER*(MAX_LEN_FNAM) tmpInitialFile

      INTEGER i,j,k,bi,bj,iCl

      REAL bla1

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

#endif /* ICEBERG_3D_STATE */
#ifdef ICEBERG_TENDENCY

#endif /* ICEBERG_TENDENCY */
          ENDDO
         ENDDO
        ENDDO

C     2-dim.
#ifdef ICEBERG_2D_STATE

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
      IF (nIter0 .EQ.ICEBERG_Iter0) THEN
       DO iCl = 1, ICEBERG_numClUsed

        tmpInitialFile = ICEBERG_initialFile(iCl)
        IF ( tmpInitialFile .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( tmpInitialFile, ' ',
     &                         iceberg_distr, 0,
     &                         myThid )
c         _EXCH_XYZ_RL(iceberg_distr(1-OLx,1-OLy,1,1,1,iCl),myThid)
         CALL EXCH_XYZ_RL( iceberg_distr,
     &                         Nr, myThid )


         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j = 1-OLy,sNy+OLy
            DO i = 1-OLx,sNx+Olx
              bla1 = bla1 + iceberg_distr(i,j,1,bi,bj,iCl)

              IF (iceberg_distr(i,j,1,bi,bj,iCl) .NE. 0) THEN
               PRINT *, 'debug coor', i, j
              ENDIF
            ENDDO
           ENDDO
          ENDDO
         ENDDO
         PRINT *, 'debug iceberg_distr' , bla1

        ELSE
                WRITE(msgBuf,'(A)') 'Not enough input files for iceberg
     &           classes. Make sure every iceberg class has one input
     &           file read in data.iceberg'
                CALL PRINT_ERROR(msgBuf,myThid)
                STOP 'ABNORMAL END: ICEBERG_INIT_VARIA'
        ENDIF
       ENDDO
      ENDIF

#endif /* ICEBERG_3D_STATE */

#ifdef ICEBERG_2D_STATE

#endif /* ICEBERG_2D_STATE */

      ELSE
C--   restart from a pickup:
         CALL ICEBERG_READ_PICKUP( nIter0, myThid )

C--   end start-from-iter-zero if/else block
      ENDIF

#endif /* ALLOW_ICEBERG */

      RETURN
      END
