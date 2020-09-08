#include "ICEBERG_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

CBOP
C !ROUTINE: ICEBERG_ADVDIFF

C !INTERFACE: ==========================================================
      SUBROUTINE ICEBERG_MOM_EQ(uc, vc , nCl,
     I                  MaskU,MaskV,
     I                  myTime, myIter, myThid)

C !DESCRIPTION: \bv
C     *===========================================================*
C     | SUBROUTINE ICEBERG_MOM_EQ
C     | o solve momentum equation with explicit Euler forward method
C     | 1) calculate parts of momentum equation
C     | 2) Sum parts of Momentum equation
C     | 3) use Euler forward method to calculate iceberg velocity in next time step
C     *===========================================================*
C
c   *------------
c   Momentum Equation consists of
c       Coriolis Force
c       Ocean drag
c       Atmosphere drag
c       Seaice drag
c       Wave radiation force
c       Pressure Gradient force
c   *------------



C !USES: ===============================================================
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"
#include "DYNVARS.h"
#ifdef ALLOW_EXF
# include "EXF_FIELDS.h"
# include "EXF_PARAM.h"
#endif

C !INPUT/OUTPUT PARAMETERS: ===================================================
C     === Routine arguments ===
C     uc/vc     :: current ice velocity on C-grid;
C     myTime    :: current time in simulation
C     myIter    :: current iteration number in simulation
C     myThid    :: my Thread Id number
C     nCl       :: number of class the velocity is calculated for
      _RL uc   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vc   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL MaskU (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL MaskV (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL myTime
      INTEGER myIter
      INTEGER myThid
      INTEGER nCl

C !LOCAL VARIABLES: ====================================================
C     === Local variables ===
C     i,j,bi,bj :: Loop counters
C     nCl       :: Loop counter for iceberg size classes
C     u_temp
      INTEGER i, j, bi, bj
      _RL utemp(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vtemp(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Fcu(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Fcv(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Fwu(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Fwv(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Fau(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Fav(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Fsu(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Fsv(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Fru(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Frv(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Fpu(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Fpv(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dvdt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dudt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

CEOP



        DO bj=myByLo(myThid),myByHi(myThid)
          DO bi=myBxLo(myThid),myBxHi(myThid)

           DO j = 1-OLy+1, sNy+OLy
            DO i = 1-OLx+1, sNx+OLx
              utemp(i,j,bi,bj) = uc(i,j,bi,bj)
              vtemp(i,j,bi,bj) = vc(i,j,bi,bj)

            ENDDO
           ENDDO

C Set up mass per unit area
C
            DO j=1-OLy+1,sNy+OLy
              DO i=1-OLx+1,sNx+OLx

                iceberg_MassC(i,j,bi,bj) =
     &                 ICEBERG_rhoSeaice * iceberg_distr(i,j,bi,bj,nCl)
                iceberg_MassU(i,j,bi,bj) = ICEBERG_rhoSeaice * 0.5 *
     &                 (iceberg_distr(i,j,bi,bj,nCl)
     &                 + iceberg_distr(i-1,j,bi,bj,nCl))
                iceberg_MassV(i,j,bi,bj) = ICEBERG_rhoSeaice * 0.5 *
     &                 (iceberg_distr(i,j,bi,bj,nCl)
     &                 + iceberg_distr(i,j-1,bi,bj,nCl))

C     HINT this Enddo maybe has to be shifted behind the calculation
              ENDDO
            ENDDO


c   *------- Coriolis Force -------*
C  Fc = -M * f * k x v_i
C _fCori  ::  Makro from model for Coriolis parameter at cell center
C           DO j=1-OLy+1,sNy+OLy
C              DO i=1-OLx+1,sNx+OLx
C                Fcu(i,j,bi,bj) = iceberg_MassU(i,j,bi,bj)
C     &            * _fCori
C     &            * vtemp(i,j,bi,bj)


C                Fcv(i,j,bi,bj) = - iceberg_MassV(i,j,bi,bj)
C     &            * _fCori
C     &            * utemp(i,j,bi,bj)
C              ENDDO
C            ENDDO


c   *------- Ocean Drag -------*
C   Fw = 1/2 * rho_ocean * C_ocean * A_ocean * |v_ocean - v_i| (v_ocean -v_i)
c    Fw = 1/2 * ICEBERG_rhoOcean * ICEBERG_OceanDrag * 0.5 (iceberg_area(i) + iceberg_area(i-1)
            DO j=1-OLy+1,sNy+OLy
              DO i=1-OLx+1,sNx+OLx
                Fwu(i,j,bi,bj) = 0.5 * ICEBERG_rhoOcean
     &            * ICEBERG_OceanDrag
     &            * 0.5 * (iceberg_area(i,j,bi,bj,nCl)
     &            + iceberg_area(i-1,j,bi,bj,nCl))
     &            * ABS( uVel(i,j,1,bi,bj)-utemp(i,j,bi,bj) )
     &            * (uVel(i,j,1,bi,bj)-utemp(i,j,bi,bj))



                Fwv(i,j,bi,bj) = 0.5 * ICEBERG_rhoOcean
     &            * ICEBERG_OceanDrag
     &            * 0.5 * (iceberg_area(i,j,bi,bj,nCl)
     &            + iceberg_area(i,j-1,bi,bj,nCl))
     &            * ABS( vVel(i,j,1,bi,bj)-vtemp(i,j,bi,bj) )
     &            * (vVel(i,j,1,bi,bj)-vtemp(i,j,bi,bj))
              ENDDO
            ENDDO

c   *------ Atmosphere Drag ------*
C   Fa = 1/2 * rho_atm * C_atm * A_atm * |v_atm - v_i| (v_atm -v_i)
c   Fa = 1/2 * ICEBERG_rhoAir * ICEBERG_AirDrag *

c HINT Same as Fwu and Fwv but instead of uVel and vVel we need WINDX and WINDY???
#ifdef ALLOW_EXF
            DO j=1-OLy+1,sNy+OLy
              DO i=1-OLx+1,sNx+OLx
                Fau(i,j,bi,bj) = 0.5 * ICEBERG_rhoAir
     &            * ICEBERG_AirDrag
     &            * 0.5 * ( iceberg_area(i,j,bi,bj,nCl)
     &            + iceberg_area(i-1,j,bi,bj,nCl) )
     &            * ABS( UWIND(i,j,bi,bj)-utemp(i,j,bi,bj) )
     &            * (UWIND(i,j,bi,bj)-utemp(i,j,bi,bj))
                IF (UWIND(i,j,bi,bj) .GT. 0) THEN

                 PRINT *, 'debug uwind', UWIND(i,j,bi,bj)
c              PRINT *, 'debug Vwind', Fau(i,j,bi,bj)
                ENDIF
                Fav(i,j,bi,bj) = 0.5 * ICEBERG_rhoAir
     &            * ICEBERG_AirDrag
     &            * 0.5 * (iceberg_area(i,j,bi,bj,nCl)
     &            + iceberg_area(i,j-1,bi,bj,nCl))
     &            * ABS( VWIND(i,j,bi,bj)-vtemp(i,j,bi,bj) )
     &            * (VWIND(i,j,bi,bj)-vtemp(i,j,bi,bj))
              ENDDO
            ENDDO

#endif /* ALLOW_EXF */


c   *------ Seaice Drag ------*
c   Fs = 1/2 * rho_sice * C_sice * A_sice * |v_sice - v_i| (v_sice - v_i
c   Fs = 1/2 * ICEBERG_rhoSeaice * ICEBERG_SeaiceDrag
C TODO define area for interaction with seaice
C TODO include seaice package and seaice velocity
C HINT have to identify if there is seaice (Mask) otherwise drag = 0


C           DO j=1-OLy+1,sNy+OLy
C              DO i=1-OLx+1,sNx+OLx
C                IF (Seaice_maskU .EQ. 0.0) THEN
C                 Fsu(i,j,bi,bj) = 0.5 * ICEBERG_rhoSeaice
C     &             * ICEBERG_SeaiceDrag
C     &             * 0.5 * ( iceberg_area(i,j,bi,bj,nCl)
C     &             + iceberg_area(i-1,j,bi,bj,nCl) )
C     &             * ABS( UICE(i,j,bi,bj)-utemp(i,j,bi,bj) )
C     &             * (UICE(i,j,bi,bj)-utemp(i,j,bi,bj))
C                 ELSE
C                  Fsu(i,j,bi,bj) = 0. _d 0
C                 ENDIF

C                 IF (Seaice_maskV .EQ. 0.0) THEN
C                  Fsv(i,j,bi,bj) = 0.5 * ICEBERG_rhoSeaice
C     &             * ICEBERG_SeaiceDrag
C     &             * 0.5 * (iceberg_area(i,j,bi,bj,nCl)
C     &             + iceberg_area(i,j-1,bi,bj,nCl))
C     &             * ABS( VICE(i,j,bi,bj)-vtemp(i,j,bi,bj) )
C     &             * (VICE(i,j,bi,bj)-vtemp(i,j,bi,bj))
C                 ELSE
C                  Fsv(i,j,bi,bj) = 0. _d 0
C                 ENDIF
C              ENDDO
C            ENDDO


c   *------ Wave Radiation Force ------*
c   Fr = 1/4 * rho_ocean * gravity * a * L * v_air/|v_air|


C           DO j=1-OLy+1,sNy+OLy
C              DO i=1-OLx+1,sNx+OLx
C                Fru(i,j,bi,bj) = 0.25 * ICEBERG_rhoOcean
C     &            * ICEBERG_gravity
C     &            * 0.5 * (0.02025 * ABS( UWIND(i,j,bi,bj) ) ** 2 )
C     &            + HINT length of iceberg ???
C     &            * UWIND(i,j,bi,bj)/ ABS( UWIND(i,j,bi,bj)


C                Frv(i,j,bi,bj) = 0.25 * ICEBERG_rhoOcean
C     &            * ICEBERG_gravity
C     &            * 0.5 * (0.02025 * ABS( VWIND(i,j,bi,bj) ) ** 2 )
C     &            + HINT length of iceberg ???
C     &            * UWIND(i,j,bi,bj)/ ABS( UWIND(i,j,bi,bj)
C              ENDDO
C            ENDDO



c   *------ Pressure Gradient Force ------*
c   Fp = -M * g * Nabla * sea surface slope
C           DO j=1-OLy+1,sNy+OLy
C              DO i=1-OLx+1,sNx+OLx
C                Fpu(i,j,bi,bj) = iceberg_MassU(i,j,bi,bj)
C     &            * ICEBERG_gravity
C     &            + HINT sea surface slope ???


C                Fpv(i,j,bi,bj) = iceberg_MassV(i,j,bi,bj)
C     &            * ICEBERG_gravity
C     &            + HINT sea surface slope ???
C              ENDDO
C            ENDDO




C   *------ Momentum Equation ------*
C   dv_i/dt = Coriolis + Fw + Fa + Fs + Fr + Fp
            DO j=1-OLy+1,sNy+OLy
              DO i=1-OLx+1,sNx+OLx
                dudt(i,j,bi,bj) = Fwu(i,j,bi,bj) + Fau(i,j,bi,bj)
                dvdt(i,j,bi,bj) = Fwv(i,j,bi,bj) + Fav(i,j,bi,bj)
              ENDDO
            ENDDO



C   *------ Euler Forward Step ------*
c   v_i[n+1] = v_i[n] + dv_i * dt
            DO j=1-OLy+1,sNy+OLy
              DO i=1-OLx+1,sNx+OLx
                IF (MaskU(i,j,bi,bj) .GT. 0. _d 0) THEN
                 uc(i,j,bi,bj) = utemp(i,j,bi,bj)
     &              + dudt(i,j,bi,bj)/ iceberg_MassU(i,j,bi,bj)
     &              * ICEBERG_deltaTtherm
c
                ELSE
                 uc(i,j,bi,bj) = 0

                ENDIF

                IF (MaskV(i,j,bi,bj) .GT. 0) THEN
                 vc(i,j,bi,bj) = vtemp(i,j,bi,bj)
     &              + dvdt(i,j,bi,bj)/ iceberg_MassV(i,j,bi,bj)
     &              * ICEBERG_deltaTtherm
c
                ELSE
                 vc(i,j,bi,bj) = 0
                ENDIF

              ENDDO
            ENDDO

          ENDDO
        ENDDO



      RETURN
      END