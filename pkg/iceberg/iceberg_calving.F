#include "ICEBERG_OPTIONS.h"
C     Without ALLOW_SEAICE, this would never work anyway.
C     Needs to be checked in iceberg_check.F
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

CBOP
C !ROUTINE: ICEBERG_ADVDIFF

C !INTERFACE: ==========================================================
      SUBROUTINE ICEBERG_CALVING( nCl,
     I                  myTime, myIter, myThid)

C !DESCRIPTION: \bv
C     *===========================================================*
C     | SUBROUTINE ICEBERG_CALVING
C     | o in this routine the calving of icebergs is done
C     |   The calving sites are given by an input file as well as the rates
C     *===========================================================*
C



C !USES: ===============================================================
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SURFACE.h"
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"
#include "DYNVARS.h"

C !INPUT/OUTPUT PARAMETERS: ===================================================
C     === Routine arguments ===
C     myTime    :: current time in simulation
C     myIter    :: current iteration number in simulation
C     myThid    :: my Thread Id number
C     nCl       :: number of class the velocity is calculated for
      _RL myTime
      INTEGER myIter
      INTEGER myThid
      INTEGER nCl

C !LOCAL VARIABLES: ====================================================
C     === Local variables ===
C     i,j,bi,bj :: Loop counters
C     addV       :: Volumen to be added to iceberg_concentration
C     u_temp
      INTEGER i, j, bi, bj
      _RL calved_iceberg
      _RL addV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)


CEOP
      calvRateInKgpersecond = .TRUE.
C TODO fix if statement: calRateInkgpersecond wird bei einigen Kernen vor der declaration aufgerufen
C ICEBERG_calvrate :: field with calving rates at calving spots matching the calving mask

       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1 ,sNy
          DO i= 1,sNx
          addV(i,j,bi,bj) = 0. _d 0
          calved_iceberg = 0. _d 0
           IF( iceberg_CalvRate(i,j,bi,bj,nCl) .GT. 0. _d 0) THEN
             IF(calvRateInKgpersecond .EQV. .TRUE.)THEN ! calving rate in kg/s
               addV(i,j,bi,bj) =
     &          iceberg_CalvRate(i,j,bi,bj,nCl)/ICEBERG_rhoIceberg
     &          * ICEBERG_deltaT
             ELSE ! calving rate in m3/s
               addV(i,j,bi,bj) = iceberg_CalvRate(i,j,bi,bj,nCl)
     &         * ICEBERG_deltaT
             ENDIF
             iceberg_Reservoir(i,j,bi,bj,nCl) =
     &         iceberg_Reservoir(i,j,bi,bj,nCl)
     &         + addV(i,j,bi,bj)/iceberg_Vsgl(nCl)
c if iceberg reservoir contains more than a full iceberg, then release it
            IF (iceberg_Reservoir(i,j,bi,bj,nCl) .GT. 1. _d 0) THEN
                calved_iceberg = AINT(iceberg_Reservoir(i,j,bi,bj,nCl))
            PRINT*, 'dbg calving before',
     &       iceberg_concentration(i,j,bi,bj,nCl),
     &       calved_iceberg,iceberg_Reservoir(i,j,bi,bj,nCl),
     &       iceberg_Reservoir(i,j,bi,bj,nCl) - calved_iceberg
                iceberg_concentration(i,j,bi,bj,nCl) =
     &            iceberg_concentration(i,j,bi,bj,nCl) + calved_iceberg
                iceberg_Reservoir(i,j,bi,bj,nCl) =
     &            iceberg_Reservoir(i,j,bi,bj,nCl) - calved_iceberg
            ENDIF

           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO



      RETURN
      END
