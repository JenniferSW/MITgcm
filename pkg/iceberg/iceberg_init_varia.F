
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
      CHARACTER*(MAX_LEN_FNAM) tmpInitialConcentration
      CHARACTER*(MAX_LEN_FNAM) tmpCalvFile


      INTEGER i,j,k,bi,bj,iCl



      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

C     Initialize arrays in common blocks :
C     2-dim.
         DO iCl = 1,ICEBERG_numClUsed
          DO j=1-Oly,sNy+OLy
           DO i=1-Olx,sNx+Olx

            iceberg_concentration(i,j,iCl,bi,bj) = 0. _d 0
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

        tmpInitialConcentration = ICEBERG_initConcentration(iCl)

         IF ( tmpInitialConcentration.EQ.' ' ) THEN
            WRITE(msgBuf,'(A)') 'Not enough input files for ICEBERG
     &      Concentration. Make sure every iceberg class has one input
     &      file to be read in data.iceberg'
            CALL PRINT_ERROR(msgBuf,myThid)
            STOP 'ABNORMAL END: ICEBERG_INIT_VARIA'
         ENDIF

         CALL READ_FLD_XY_RL( tmpInitialConcentration, ' ',
     &                iceberg_concentration(1-OLx,1-OLy,iCl,1,1),
     &                0, myThid )

        ENDDO
C        _EXCH_XY_RL(iceberg_concentration , myThid)
        CALL EXCH_3D_RL( iceberg_concentration, ICEBERG_numClUsed
     &    , myThid )


       ENDIF

      ELSE
C--   restart from a pickup:
         CALL ICEBERG_READ_PICKUP( nIter0, myThid )
C--   end start-from-iter-zero if/else block
      ENDIF

      DO iCl = 1, ICEBERG_numClUsed
        tmpCalvFile = ICEBERG_CalvFile(iCl)
        IF ( tmpCalvFile .EQ.' ' ) THEN
            WRITE(msgBuf,'(A)') 'Not enough input files for ICEBERG
     &      Calving Rates. Make sure every iceberg class has one input
     &      file to be read in data.iceberg'
            CALL PRINT_ERROR(msgBuf,myThid)
            STOP 'ABNORMAL END: ICEBERG_INIT_VARIA'
        ENDIF

        CALL READ_FLD_XY_RL( tmpCalvFile, ' ',
     &                iceberg_CalvRate(1-OLx,1-OLy,iCl,1,1),
     &                0, myThid )

         ENDDO
C        _EXCH_XY_RL(iceberg_CalvRate , myThid)
         CALL EXCH_3D_RL( iceberg_CalvRate, ICEBERG_numClUsed, myThid )


C calculate iceberg length from iceberg width
      DO iCl = 1,ICEBERG_numClUsed

      iceberg_length(iCl) = 1.5 * ICEBERG_width(iCl)
      iceberg_Vsgl(iCl) = ICEBERG_height(iCl)
     &   * ICEBERG_width(iCl) * iceberg_length(iCl)

C initialize iceberg velocity and mass loss
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
          DO j=1-Oly,sNy+OLy
           DO i=1-Olx,sNx+Olx
            iceberg_uvel(i,j,iCl,bi,bj) = 0. _d 0
            iceberg_vVel(i,j,iCl,bi,bj) = 0. _d 0

            iceberg_Vloss(i,j,iCl,bi,bj) = 0. _d 0
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

#endif /* ALLOW_ICEBERG */

      RETURN
      END
