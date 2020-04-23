C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CBOP
C !ROUTINE: ICEBERG_SET_CLASS_LABELS

C !INTERFACE: =========================================================
      SUBROUTINE ICEBERG_SET_CLASS_LABELS(outLbl, numLbl, myThid)


C !Description:
C Creates Label for iceberg classes in order to name the diagnostics
C output.
C It is possible to write labels up to 99 size classes


c !USES: ==============================================================
      IMPLICIT NONE


C !INPUT PARAMETERS: ==================================================
C numLbl    :: Number of Labels that are needed (=ICEBERG_numCLUsed)
C myThid    :: my Thread Id number
      INTEGER numLbl
      INTEGER myThid




C !OUTPUT PARAMETERS: =================================================
C outLbl     :: Output Label for each size class (=ICEBERG_CLbl)
      CHARACTER*8 outLbl(numLbl)



C !LOCAL VARIABLES: ===================================================
      INTEGER iCl


CEOP
      DO iCl = 1, numLbL
        IF (iCl .LT. 10) THEN
c       if number of size class is lower than 10 a '0' has to be added
c       to the label
          WRITE(outLbl(iCl),'(A6,I1,A1)') 'IBCL_0', iCl, ' '

        ELSE
C       otherwise the label is just the number of the size class
          WRITE(outLbl(iCl),'(A5,I2,A1)') 'IBCL_', iCl, ' '
        
        ENDIF
      ENDDO

      RETURN
      END
