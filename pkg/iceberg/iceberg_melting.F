#include "ICEBERG_OPTIONS.h"
#include "SEAICE_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

CBOP
C !ROUTINE: ICEBERG_ADVDIFF

C !INTERFACE: ==========================================================
      SUBROUTINE ICEBERG_MELTING(
     I                  myTime, myIter, myThid)

C !DESCRIPTION: \bv
C     *===========================================================*
C     | SUBROUTINE ICEBERG_MELTING
C     | o calculate mass loss from melting terms on C-point
C     *===========================================================*
C
c   *------------
c   Melting consists of
c       Basal Melting
c       Wave erosion
c       Buoyant convection
c   *------------



C !USES: ===============================================================
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SURFACE.h"
#include "FFIELDS.h"
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"
#include "DYNVARS.h"
#ifdef ALLOW_EXF
# include "EXF_FIELDS.h"
# include "EXF_PARAM.h"
#endif
#ifdef ALLOW_SEAICE
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
#include "SEAICE.h"
#endif

C !INPUT/OUTPUT PARAMETERS: ===================================================
C     === Routine arguments ===
C     myTime    :: current time in simulation
C     myIter    :: current iteration number in simulation
C     myThid    :: my Thread Id number
      _RL myTime
      INTEGER myIter
      INTEGER myThid


C !LOCAL VARIABLES: ====================================================
C     === Local variables ===
C     i,j,bi,bj :: Loop counters
C     iCl       :: Loop counter for iceberg size classes
C     u_temp
      INTEGER i, j, bi, bj, iCl
      _RL factor
      _RL S
      _RL recip_nCl
      _RL norm_b, norm_aw
      _RL Mb(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL Mv(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL Me(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      factor = 86400. _d 0
      DO iCl = 1, ICEBERG_numClUsed
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)

         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx

              Mb(i,j) = 0. _d 0
              Mv(i,j) = 0. _d 0
              Me(i,j) = 0. _d 0

          ENDDO
         ENDDO

        DO j = 1-OLy, sNy+OLy
         DO i = 1-OLx, sNx+OLx
          IF(ICEBERG_MaskC(i,j,bi,bj) .EQ. 0. _d 0) THEN
           iceberg_Vloss(i,j,bi,bj,iCl) = 0. _d 0

          ELSE


C *-------- Basal convection --------*
C Mb = 0.58 |v_w - v_i|^0.8 (T_w - T_i)/(L_a^0.2) [m/d]
C T_w - water temperature
C T_i - surface temperature of ice (-4°C in Bigg et al 1997)
C L_a - berg's horizontal long axis [m]

            norm_b = sqrt(
     &       (0.5 * ((uVel(i,j,1,bi,bj)
     &       -iceberg_uvel(i,j,bi,bj,iCl))
     &       + (uVel(i+1,j,1,bi,bj)
     &       -iceberg_uvel(i+1,j,bi,bj,iCl))))**2
     &       +
     &       (0.5 * ((vVel(i,j,1,bi,bj)
     &       -iceberg_vvel(i,j,bi,bj,iCl))
     &       +(vVel(i,j+1,1,bi,bj)
     &       -iceberg_vvel(i,j+1,bi,bj,iCl))))**2
     &       )

            Mb(i,j) = 1/factor * 0.58 *
     &       norm_b**0.8
     &       * (SST(i,j,bi,bj) - ICEBERG_temp)
     &       /sqrt(iceberg_area(i,j,bi,bj,iCl)*rA(i,j,bi,bj)) **0.2
            PRINT*, 'debug Mb',i,j, Mb(i,j)/factor
     &       ,Mb(i,j)*1/factor



C *-------- wave ersosion --------*
C Me = 1/12 (1+cos(A^3*pi)(To+2)*S
C A :: seaice concentration
C S_s = 3/2* |u_a-u_o|^1/2 + 1/10 |u_a-u_o|
C S_s - seat state calculated as a function of wind speed according to marine beaufeort scale

            norm_aw = sqrt((UWIND(i,j,bi,bj)
     &       - 0.5*(uVel(i,j,1,bi,bj) + uVel(i-1,j,1,bi,bj)))**2
     &       + (VWIND(i,j,bi,bj)
     &       - 0.5*(vVel(i,j,1,bi,bj) + vVel(i,j-1,1,bi,bj)))**2)

            S = 3/2 * norm_aw**0.5
     &       + 1/10 * norm_aw


            Me(i,j) = 1/factor * 1/12
     &       * (1+ COS(AREA(i,j,bi,bj) **3*pi))
     &       * (SST(i,j,bi,bj) + 2)
     &       * S
            PRINT*, 'debug Me',i,j, Me(i,j)

C *-------- buoyant convection --------*
C Mv = 7.62*10^-3 T_w +1.29*10^-3 T_w^2
C T_w water temperature
            Mv(i,j) = 1/factor * (
     &       0.00762* SST(i,j,bi,bj)
     &       + 0.00129 * SST(i,j,bi,bj)**2)
            PRINT*, 'debug Mv',i,j, Mv(i,j), SST(i,j,bi,bj)



C *-------- mass balance of iceberg --------*
C rho dV/dt rho(-L^2Mb - 2TLMe -2TL Mv))
C iceberg_Vloss = basal melting (change of height)
C               + wave erosion (change of area)
C               + buoyant convection (change of area)

            iceberg_Vloss(i,j,bi,bj,iCl) =
     &       iceberg_area(i,j,bi,bj,iCl)*Mb(i,j) * iceberg_deltaT
     &       + iceberg_distr(i,j,bi,bj,iCl) * Me(i,j) * iceberg_deltaT
     &       + iceberg_distr(i,j,bi,bj,iCl) * Mv(i,j) * iceberg_deltaT

            PRINT*, 'debug melt vloss',i,j, iceberg_Vloss(i,j,bi,bj,iCl)
C calculate new iceberg volume in grid cell
            iceberg_distr(i,j,bi,bj,iCl) =
     &       iceberg_distr(i,j,bi,bj,iCl) - iceberg_Vloss(i,j,bi,bj,iCl)

C catch values that are too small or values that are negative because of too strong melting
            IF (iceberg_distr(i,j,bi,bj,iCl).LE.1. _d -5) THEN
             iceberg_distr(i,j,bi,bj,iCl) = 0. _d 0
             iceberg_area(i,j,bi,bj,iCl) = 0. _d 0
            ENDIF

C calculate new 'dimension' by subtraction of change of height form iceberg_height
            iceberg_height(i,j,bi,bj,iCl) =
     &       iceberg_height(i,j,bi,bj,iCl)
     &       - iceberg_area(i,j,bi,bj,iCl)*Mb(i,j) * iceberg_deltaT

C catch values that are too small or values that are negative because of too strong melting
            IF (iceberg_height(i,j,bi,bj,iCl).LE.1. _d -5) THEN
             iceberg_height(i,j,bi,bj,iCl) = 0. _d 0
             iceberg_area(i,j,bi,bj,iCl) = 0. _d 0
            ENDIF

C calculate new concentration of iceberg from new volume and new height
            iceberg_area(i,j,bi,bj,iCl) =
     &       iceberg_distr(i,j,bi,bj,iCl)/iceberg_height(i,j,bi,bj,iCl)

cC catch if iceberg concentration is now too low
c            recip_nCl = 1. _d 0 / float(ICEBERG_numClUsed)
c            IF (iceberg_distr(i,j,bi,bj,iCl).GT.0) THEN
c             iceberg_area(i,j,bi,bj,iCl)=
c     &         MAX(iceberg_area(i,j,bi,bj,iCl),1. _d -5*recip_nCl)
c            ENDIF

           ENDIF ! end if iceberg_maskC = 0
          ENDDO ! end i loop
         ENDDO ! end j loop
        ENDDO ! end bi loop
       ENDDO ! end bj loop
      ENDDO
      RETURN
      END
