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
C     | o solve momentum equation with Euler forward method
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
      _RL duw
      _RL dvw
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
      _RL urhs(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vrhs(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL CDw(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL CDa(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL CDs(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

CEOP



        DO bj=myByLo(myThid),myByHi(myThid)
          DO bi=myBxLo(myThid),myBxHi(myThid)

           DO j = 1-OLy, sNy+OLy
            DO i = 1-OLx, sNx+OLx
              utemp(i,j,bi,bj) = uc(i,j,bi,bj)
              vtemp(i,j,bi,bj) = vc(i,j,bi,bj)

            ENDDO
           ENDDO

C Set up mass per unit area
C
            DO j=1-OLy,sNy+OLy-1
              DO i=1-OLx,sNx+OLx-1

                iceberg_MassC(i,j,bi,bj) =
     &                 ICEBERG_rhoSeaice * iceberg_distr(i,j,bi,bj,nCl)
                iceberg_MassU(i,j,bi,bj) = ICEBERG_rhoSeaice * 0.5 *
     &                 (iceberg_distr(i,j,bi,bj,nCl)
     &                 + iceberg_distr(i+1,j,bi,bj,nCl))
                iceberg_MassV(i,j,bi,bj) = ICEBERG_rhoSeaice * 0.5 *
     &                 (iceberg_distr(i,j,bi,bj,nCl)
     &                 + iceberg_distr(i,j+1,bi,bj,nCl))

C     HINT this Enddo maybe has to be shifted behind the calculation
              ENDDO
            ENDDO


c   *------- Coriolis Force -------*
C TODO adapt Coriolis with discretisation from seaice_evp
C  Fc = -M * f * k x v_i[n]
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
C   explicit: Fw = 1/2 * rho_ocean * C_ocean * A_ocean * |v_ocean - v_i| (v_ocean -v_i)
c   implicit: rhsw = 1/2* rho:ocean * C_ocean * A_ocean + |v_ocean - v_i[n]| * v_ocean

            DO j=1-OLy,sNy+OLy-1
              DO i=1-OLx,sNx+OLx-1
C u/vVel and u/vtemp are both located on u/v Points and difference is interpolated to C-Points
               duw = 0.5 * ((uVel(i,j,1,bi,bj)-utemp(i,j,bi,bj))
     &          + (uVel(i+1,j,1,bi,bj)-utemp(i+1,j,bi,bj)))
               dvw = 0.5 * ((vVel(i,j,1,bi,bj)-vtemp(i,j,bi,bj))
     &          +(vVel(i,j+1,1,bi,bj)-vtemp(i,j+1,bi,bj)))
c Calculate Norm on C-Point with constants for ocean stress
               CDw(i,j) = 0.5 * ICEBERG_rhoOcean
     &            * ICEBERG_OceanDrag
     &            * iceberg_area(i,j,bi,bj,nCl)
     &            * sqrt(duw**2+dvw**2)


C TODO Error appears when only running this section,
C when everything from the sqrt downwards is commented, then the error disappears, so it is the sqrt
C The values should not be able to get negative because of the square. Possibly there are NANs?
C Maybe Fortran cannot also not deal with sqrt(0)?

C Fwu/v are calculated on C-Point for implicit scheme
C iceberg_area is located on C-Points
c                Fwu(i,j,bi,bj) = CDw(i,j)
c     &            * 0.5 * ( uVel(i,j,1,bi,bj) + uVel(i-1,j,1,bi,bj) )

c                Fwv(i,j,bi,bj) = CDw(i,j)
c     &            * 0.5 * ( vVel(i,j,1,bi,bj) + vVel(i,j-1,1,bi,bj) )
              ENDDO
            ENDDO

c   *------ Atmosphere Drag ------*
C   Fa = 1/2 * rho_atm * C_atm * A_atm * |v_atm - v_i| (v_atm -v_i)
C   Assumed that v_i is much smaller than v_atm (v_i << v_atm)
c   rhsa = 1/2 * rho_atm * C_atm * A_atm * | v_atm| * v_atm

c HINT Same as Fwu and Fwv but instead of uVel and vVel we need WINDX and WINDY???
C TODO Fix discretization of norm and take out iceberg velocity
#ifdef ALLOW_EXF
c            DO j=1-OLy+1,sNy+OLy
c              DO i=1-OLx+1,sNx+OLx
C Uwind is located on southwest C-Grid U-point and Vwind is located on southwest C-Grid V-Point
c                dua = 0.5 * ((UWIND(i,j,bi,bj)-utemp(i,j,bi,bj))
c     &           + (UWIND(i-1,j,bi,bj)-utemp(i-1,j,bi,bj)))
c                dva = 0.5 * (VWIND(i,j,bi,bj)-vtemp(i,j,bi,bj))
c     &           +(VWIND(i,j-1,bi,bj)-vtemp(i,j-1,bi,bj)))
C Calculate Norm on C-Point with constant values for wind stress
c                CDa(i,j) = 0.5 * ICEBERG_rhoAir
c     &            * ICEBERG_AirDrag
c     &            * iceberg_area(i,j,bi,bj,nCl)
c     &            * sqrt(dua**2+dva**2)

c                Fau(i,j,bi,bj) = CDa(i,j)
c     &            * 0.5 * (UWIND(i,j,bi,bj) + UWIND(i-1,j,bi,bj))

c                Fav(i,j,bi,bj) = CDa(i,j)
c     &            * 0.5 * (VWIND(i,j,bi,bj) + VWIND(i,j-1,bi,bj))
c              ENDDO
c            ENDDO

#endif /* ALLOW_EXF */


c   *------ Seaice Drag ------*
c   Fs = 1/2 * rho_sice * C_sice * A_sice * |v_sice - v_i| (v_sice - v_i
c   Fs = 1/2 * ICEBERG_rhoSeaice * ICEBERG_SeaiceDrag
C TODO define area for interaction with seaice
C TODO include seaice package and seaice velocity
C HINT have to identify if there is seaice (Mask) otherwise drag = 0
C TODO fix discretisation for norm

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
C TODO fix discretization of the norm

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
C TODO change name from dvdt to rhs
            DO j=1-OLy+1,sNy+OLy
              DO i=1-OLx+1,sNx+OLx
                urhs(i,j,bi,bj) = Fwu(i,j,bi,bj)
c                + Fau(i,j,bi,bj)
                vrhs(i,j,bi,bj) = Fwv(i,j,bi,bj)
c                + Fav(i,j,bi,bj)
              ENDDO
            ENDDO



C   *------ Euler Forward Step ------*
c   v_i[n+1] = v_i[n] + dv_i * dt
C Calculated at C-Point
            DO j=1-OLy,sNy+OLy-1
             DO i=1-OLx,sNx+OLx-1
               IF(iceberg_MaskC(i,j,bi,bj) .EQ. 0.0) THEN
C     &             .AND. iceberg_MaskU(i+1,j,bi,bj) .EQ. 0.0) THEN
                  uc(i,j,bi,bj) = 0. _d 0
                ELSE
c                  PRINT * ,'debug mass',iceberg_MassC(i,j,bi,bj), i, j
c                  PRINT*,'dbg mass',iceberg_distr(i,j,bi,bj,nCl)
                  PRINT*,'mask', iceberg_MaskC(i,j,bi,bj), i, j
                  PRINT*,'mask+1',iceberg_MaskC(i+1,j,bi,bj), i,j
c                  uc(i,j,bi,bj) =
c     &              (iceberg_MassC(i,j,bi,bj)
c     &              * 0.5 * (utemp(i,j,bi,bj)+utemp(i+1,j,bi,bj))
c     &              - Fwu(i,j,bi,bj) * ICEBERG_deltaTtherm)
c     &              / (iceberg_MassC(i,j,bi,bj)
c     &              * ICEBERG_deltaTtherm * CDw(i,j))
                ENDIF

c                IF(iceberg_MaskV(i,j,bi,bj) .EQ. 0.0
c     &             .AND. iceberg_MaskV(i,j+1,bi,bj) .EQ. 0.0) THEN
c                  vc(i,j,bi,bj) = 0. _d 0
c                ELSE
c                 vc(i,j,bi,bj) =
c     &              (iceberg_MassC(i,j,bi,bj)
c     &              * 0.5 * (vtemp(i,j,bi,bj)+vtemp(i,j+1,bi,bj))
c     &              - Fwv(i,j,bi,bj) * ICEBERG_deltaTtherm)
c     &              / (iceberg_MassC(i,j,bi,bj)
c     &              * ICEBERG_deltaTtherm * CDw(i,j))
c                ENDIF

               ENDDO
             ENDDO

          ENDDO
        ENDDO



      RETURN
      END