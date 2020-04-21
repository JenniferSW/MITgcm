#include "ICEBERG_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: ICEBERG_DIAGNOSTICS_INIT

C !INTERFACE:
      SUBROUTINE ICEBERG_DIAGNOSTICS_INIT( myThid )

C     !DESCRIPTION:
C     Initialize list of all available diagnostics

C     !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS
C     !LOCAL VARIABLES:
C     === Local variables ===
C     diagNum   :: diagnostics number in the (long) list of available diag.
C     diagMate  :: diag. mate number in the (long) list of available diag.
C     diagName  :: local short name (8c) of a diagnostics
C     diagCode  :: local parser field with characteristics of the diagnostics
C              cf head of S/R DIAGNOSTICS_INIT_EARLY or DIAGNOSTICS_MAIN_INIT
C     diagUnits :: local string (16c): physical units of a diagnostic field
C     diagTitle :: local string (80c): description of field in diagnostic
      INTEGER       diagNum
      INTEGER       diagMate
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

      INTEGER iCl
      iCl = 1
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c     IF ( useDiagnostics ) THEN
        DO iCl = 1, ICEBERG_numClUsed
          WRITE(diagName,'(A5,I1,A2)') 'IBCL_',iCl,'  '
          WRITE(diagTitle,'(A,I2)') 'Iceberg distribution for class'
     &          ,iCl
          diagUnits = 'none            '
          diagCode  = 'SM      MR      '
          CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I            diagName, diagCode, diagUnits, diagTitle, 0, myThid )
        ENDDO
c        diagName  = 'MYPaSta2'
c        diagTitle = 'Scalar state variable 2 for package ICEBERG'
c        diagUnits = 'unit            '
c        diagCode  = 'SM      MR      '
c        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
c     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

c        diagName  = 'MYPaStaU'
c        diagTitle = 'uVel-type state variable for package ICEBERG'
c        diagUnits = 'unit            '
c        diagCode  = 'UM      MR      '
c        diagMate  = diagNum + 2
c        IF ( myPa_StaV_Cgrid ) diagCode(2:2) = 'U'
c        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
c     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

c        diagName  = 'MYPaStaV'
c        diagTitle = 'vVel-type state variable for package ICEBERG'
c        diagUnits = 'unit            '
c        diagCode  = 'VM      MR      '
c        diagMate  = diagNum
c        IF ( myPa_StaV_Cgrid ) diagCode(2:2) = 'V'
c        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
c     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

c        diagName  = 'MYPaSur1'
c        diagTitle = 'Surface field 1 for package ICEBERG'
c        diagUnits = 'unit            '
c        diagCode  = 'SM      L1      '
c        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
c     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

c        diagName  = 'MYPaSur2'
c        diagTitle = 'Surface field 2 for package ICEBERG'
c        diagUnits = 'unit            '
c        diagCode  = 'SM      L1      '
c        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
c     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

c        diagName  = 'MYPadTdt'
c        diagTitle = ' Pot.Temp. tendency from package ICEBERG'
c        diagUnits = 'K/s             '
c        diagCode  = 'SMR     MR      '
c        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
c     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

c        diagName  = 'MYPadSdt'
c        diagTitle = ' Salinity  tendency from package ICEBERG'
c        diagUnits = 'psu/s           '
c        diagCode  = 'SMR     MR      '
c        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
c     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

c        diagName  = 'MYPadUdt'
c        diagTitle = 'u-Velocity tendency from package ICEBERG'
c        diagUnits = 'm/s^2           '
c        diagCode  = 'UMR     MR      '
c        diagMate  = diagNum + 2
c        IF ( myPa_Tend_Cgrid ) diagCode(2:2) = 'U'
c        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
c     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

c        diagName  = 'MYPadVdt'
c        diagTitle = 'v-Velocity tendency from package ICEBERG'
c        diagUnits = 'm/s^2           '
c        diagCode  = 'VMR     MR      '
c        diagMate  = diagNum
c        IF ( myPa_Tend_Cgrid ) diagCode(2:2) = 'V'
c        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
c     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

C--     an example of how to set "free levels option" diagnostics
c        diagName  = 'MYPaLev3'
c        diagTitle = '3-levels field for package ICEBERG'
c        diagUnits = 'unit            '
c        diagCode  = 'SM      MX      '
c        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
c     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
c        CALL DIAGNOSTICS_SETKLEV( diagName, 3, myThid )

c     ENDIF

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END
