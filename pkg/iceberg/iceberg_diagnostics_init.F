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
#include "PARAMS.h"
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
      CHARACTER*(18) diagLabel

      INTEGER iCl

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C TODO Fix Diagnostics Names to maximum of 8 Characters

      IF ( useDiagnostics ) THEN
C Diagnostics for iceberg distribution
       DO iCl = 1, ICEBERG_numCl
        WRITE(diagTitle,'(A,I2)') 'Iceberg distribution for class'
     &       ,iCl
        diagUnits = 'm               '
        diagCode  = 'SM      M1      '
        WRITE(diagLabel,'(A,A)') 'IB_distr_',ICEBERG_ClLbl(iCl)
        WRITE(diagName,'(A8)') 'IBD'//ICEBERG_ClLbl(iCl)
        PRINT *, 'debug diagnostics ' ,diagLabel, diagName
C     TODO fill unit and diagCode
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle,
     &       0, myThid )
       ENDDO

C Diagnostics for iceberg area
       DO iCl = 1, ICEBERG_numCl
        WRITE(diagTitle,'(A,I2)') 'Iceberg fraction area for class'
     &       ,iCl
        diagUnits = 'none            '
        diagCode  = 'SM      M1      '
        WRITE(diagLabel,'(A,A)') 'IB_area_',ICEBERG_ClLbl(iCl)
        WRITE(diagName,'(A8)') 'IBA'//ICEBERG_ClLbl(iCl)
C TODO fill unit and diagCode
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle,
     &       0, myThid )
       ENDDO


       DO iCl = 1, ICEBERG_numCl
        WRITE(diagTitle,'(A,I2)') 'Iceberg mean height for class'
     &       ,iCl
        diagUnits = 'm               '
        diagCode  = 'SM      M1      '
        WRITE(diagLabel,'(A,A)') 'IB_height_',ICEBERG_ClLbl(iCl)
        WRITE(diagName,'(A8)') 'IBH'//ICEBERG_ClLbl(iCl)
C TODO fill unit and diagCode
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle,
     &       0, myThid )
       ENDDO

      ENDIF

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END
