#include "ICEBERG_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: ICEBERG_READPARMS

C !INTERFACE:
      SUBROUTINE ICEBERG_READPARMS( myThid )

C     !DESCRIPTION:
C     Initialize ICEBERG variables and constants.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"

C     !INPUT PARAMETERS:
      INTEGER myThid
CEOP

#ifdef ALLOW_ICEBERG

C     !LOCAL VARIABLES:
C     msgBuf     :: Informational/error message buffer
C     iUnit      :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit, iCl

      NAMELIST /ICEBERG_PARM01/
     &       ICEBERG_Iter0,
     &       ICEBERG_numClUsed,
C HINT if not used later ICEBERG_deltaTtherm has to be deleted

     &       ICEBERG_deltaTtherm,
     &       ICEBERGadvScheme,
     &       ICEBERG_initialFile

      IF (.NOT.useICEBERG) THEN
C-    pkg ICEBERG is not used
        _BEGIN_MASTER(myThid)
C-    Track pkg activation status:
C     print a (weak) warning if data.iceberg is found
         CALL PACKAGES_UNUSED_MSG( 'useICEBERG', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

      _BEGIN_MASTER(myThid)

C--   Default values for ICEBERG

C-    additional parameters:
      ICEBERG_Iter0      = 0
      ICEBERG_numClUsed = -1
      ICEBERG_deltaTtherm = dTtracerLev(1)
      ICEBERGadvScheme = 77

C-    file names for initial conditions:
      DO iCl = 1, ICEBERG_numCl
        ICEBERG_initialFile(iCl)   = ' '
      ENDDO

      WRITE(msgBuf,'(A)') 'ICEBERG_READPARMS: opening data.iceberg'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)
      CALL OPEN_COPY_DATA_FILE(
     I                     'data.iceberg', 'ICEBERG_READPARMS',
     O                     iUnit,
     I                     myThid )

C     Read parameters from open data file
      READ(UNIT=iUnit,NML=ICEBERG_PARM01)
      WRITE(msgBuf,'(A)')
     &    'ICEBERG_READPARMS: finished reading data.iceberg'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &    SQUEEZE_RIGHT , 1)
C     Close the open data file
#ifdef SINGLE_DISK_IO
      CLOSE(iUnit)
#else
      CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C Create parameters resutling from input parameters


C     If ICEBERG_numClUsed is not set the value from ICEBERG_SIZE is
C     used and error handling if no icebergs are used
      IF (ICEBERG_numClUsed .LT. 0) THEN
        ICEBERG_numClUsed = ICEBERG_numCl
      ELSE IF (ICEBERG_numClUsed .EQ. 0) THEN
              WRITE(msgBuf,'(A)') 'Seems like you do not use
     & icebergs in the model. Set ICEBERG_numClUsed > 0 or remove the iceberg
     & package from data.pkg.'
              STOP 'ABNORMAL END: ICEBERG_READPARMS'
      ENDIF

C     Error handling: Do not use more size classes than introduced
      IF (ICEBERG_numClUsed .GT. ICEBERG_numCl) THEN
              WRITE(msgBuf,'(A)') 'You try to use more size classes than
     & implemented in ICEBERG_SIZE.h. ICEBERG_numClUsed has to be
     & smaller or equal ICEBERG_numCl.'
              CALL PRINT_ERROR(msgBuf, myThid)
              STOP 'ABNORMAL END: ICEBERG_READPARMS'
      ENDIF


C     Create diagnostics labels
      IF (ICEBERG_numClUsed .GT. 0) THEN
         CALL ICEBERG_SET_CLASS_LABELS(ICEBERG_ClLbl, ICEBERG_numClUsed,
     C                    myThid)
      ENDIF


C     Make sure that we locally honor the global MNC on/off flag

#ifndef ALLOW_MNC
C     Fix to avoid running without getting any output:

#endif


      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_ICEBERG */

      RETURN
      END
