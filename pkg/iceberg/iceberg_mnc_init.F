#include "ICEBERG_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C     !ROUTINE: ICEBERG_MNC_INIT

C     !INTERFACE:
      SUBROUTINE ICEBERG_MNC_INIT( myThid )

C     !DESCRIPTION:
C     Initialize MNC lookup tables for Iceberg package.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"

C     !INPUT PARAMETERS:
      INTEGER myThid
      INTEGER iCl
      CHARACTER*8 ClassName
CEOP

#ifdef ALLOW_ICEBERG
#ifdef ALLOW_MNC

      IF (useMNC) THEN
c     fill in iceberg distribution in for using mnc
C TODO fill in border of classes as Attribute Text
        DO iCl = 1, ICEBERG_numClUsed
          WRITE(ClassName,'(A5,I2)') 'IBCL_',iCl
          CALL MNC_CW_ADD_VNAME(ClassName,
     &        'Cen_xy_Hn__L__t', 4,5,myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(ClassName,
     &        'units','[unit]', myThid)
          CALL MNC_CW_ADD_VATTR_TEXT(ClassName,'class border',
     &        'fill in border of classes',myThid)
        ENDDO

      ENDIF

#endif /* ALLOW_MNC */
#endif /* ALLOW_ICEBERG */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
