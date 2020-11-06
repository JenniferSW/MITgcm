#include "ICEBERG_OPTIONS.h"
C     Without ALLOW_SEAICE, this would never work anyway.
C     Needs to be checked in iceberg_check.F
#include "SEAICE_OPTIONS.h"
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
#ifdef ALLOW_SEAICE
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
#include "SEAICE.h"
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
      _RL dua
      _RL dva
      _RL dus
      _RL dvs
      _RL norm_a
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
      _RL CDr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

CEOP



        DO bj=myByLo(myThid),myByHi(myThid)
          DO bi=myBxLo(myThid),myBxHi(myThid)

           DO j = 1-OLy, sNy+OLy
            DO i = 1-OLx, sNx+OLx
              utemp(i,j,bi,bj) = uc(i,j,bi,bj)
              vtemp(i,j,bi,bj) = vc(i,j,bi,bj)

              CDw(i,j) = 0. _d 0
              CDa(i,j) = 0. _d 0
              CDs(i,j) = 0. _d 0
              CDr(i,j) = 0. _d 0
            ENDDO
           ENDDO

C Set up mass per unit area
C interpolation from c to u/v point
C u-point -> 0.5*( c(i,j) + c(i-1,j))
C v-point -> 0.5*( c(i,j) + c(i,j-1))

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
C Question Do I have to do the EXCH_XY....?
              ENDDO
            ENDDO


c   *------- Coriolis Force -------*
C  Fc = -M * f * k x v_i[n]
C _fCori  ::  Makro from model for Coriolis parameter at cell center
C iceberg_mass and _fcori are already on C-points, only iceberg velocity is interpolated

C HINT section is working (compilable and run without error)
           DO j=1-OLy+1,sNy+OLy-1
             DO i=1-OLx+1,sNx+OLx-1
                Fcu(i,j,bi,bj)=
     &            0.5 _d 0 * (
     &            iceberg_MassC(i  ,j,bi,bj) * _fCori(i  ,j,bi,bj)
     &            * 0.5 _d 0*( vtemp(i  ,j,bi,bj)+vtemp(i  ,j+1,bi,bj) )
     &            + iceberg_MassC(i-1,j,bi,bj) * _fCori(i-1,j,bi,bj)
     &            * 0.5 _d 0*( vtemp(i-1,j,bi,bj)+vtemp(i-1,j+1,bi,bj) )
     &            )
                Fcv(i,j,bi,bj)=
     &            - 0.5 _d 0*(
     &            iceberg_MassC(i,j  ,bi,bj) * _fCori(i,j  ,bi,bj)
     &            * 0.5 _d 0*( utemp(i,j  ,bi,bj)+utemp(i+1,  j,bi,bj) )
     &            + iceberg_MassC(i,j-1,bi,bj) * _fCori(i,j-1,bi,bj)
     &            * 0.5 _d 0*( utemp(i,j-1,bi,bj)+utemp(i+1,j-1,bi,bj) )
     &            )
             ENDDO
           ENDDO

c   *------- Ocean Drag -------*
C   explicit: Fw = 1/2 * rho_ocean * C_ocean * A_ocean * |v_ocean - v_i| (v_ocean -v_i)
c   implicit: rhsw = 1/2* rho_ocean * C_ocean * A_ocean + |v_ocean - v_i[n]| * v_ocean

C calculate drag coefficient (later shifted to separate function)
            DO j=1-OLy,sNy+OLy-1
              DO i=1-OLx,sNx+OLx-1
               duw = 0.5 * ((uVel(i,j,1,bi,bj)-utemp(i,j,bi,bj))
     &          + (uVel(i+1,j,1,bi,bj)-utemp(i+1,j,bi,bj)))
               dvw = 0.5 * ((vVel(i,j,1,bi,bj)-vtemp(i,j,bi,bj))
     &          +(vVel(i,j+1,1,bi,bj)-vtemp(i,j+1,bi,bj)))

c Calculate Norm on C-Point with constants for ocean stress
               CDw(i,j) = 0.5 * ICEBERG_rhoOcean
     &            * ICEBERG_OceanDrag
     &            * iceberg_area(i,j,bi,bj,nCl)
     &            * sqrt(duw**2+dvw**2)
              ENDDO
            ENDDO

C calculate rhs of ocean drag force
            DO j=1-OLy+1,sNy+OLy
              DO i=1-OLx+1,sNx+OLx
                Fwu(i,j,bi,bj) = 0.5 * ( CDw(i,j) + CDw(i-1,j) )
     &            * uVel(i,j,1,bi,bj)

                Fwv(i,j,bi,bj) = 0.5 * ( CDw(i,j) + CDw(i,j-1) )
     &            * vVel(i,j,1,bi,bj)
              ENDDO
            ENDDO

c   *------ Atmosphere Drag ------*
C Uwind is located on southwest C-Grid U-point and Vwind is located on southwest C-Grid V-Point
C https://mitgcm.readthedocs.io/en/latest/phys_pkgs/exf.html?highlight=uwind#field-attributes

C   Fa = 1/2 * rho_atm * C_atm * A_atm * |v_atm - v_i| (v_atm -v_i)
C   Assumed that v_i is much smaller than v_atm (v_i << v_atm)
c   rhsa = 1/2 * rho_atm * C_atm * A_atm * | v_atm| * v_atm

#ifdef ALLOW_EXF
            DO j=1-OLy+1,sNy+OLy
              DO i=1-OLx+1,sNx+OLx

                dua = 0.5 * (UWIND(i,j,bi,bj)
     &           + UWIND(i-1,j,bi,bj))
                dva = 0.5 * (VWIND(i,j,bi,bj)
     &            +VWIND(i,j-1,bi,bj))
C Calculate Norm on C-Point with constant values for wind stress
                CDa(i,j) = 0.5 * ICEBERG_rhoAir
     &            * ICEBERG_AirDrag
     &            * iceberg_area(i,j,bi,bj,nCl)
     &            * sqrt(dua**2+dva**2)
              ENDDO
            ENDDO

            DO j=1-OLy+1,sNy+OLy
              DO i=1-OLx+1,sNx+OLx
                Fau(i,j,bi,bj) = 0.5 * ( CDa(i,j) + CDa(i-1,j) )
     &            * UWIND(i,j,bi,bj)
                Fav(i,j,bi,bj) = 0.5 * ( CDa(i,j) + CDa(i,j-1) )
     &            * VWIND(i,j,bi,bj)
              ENDDO
            ENDDO

#endif /* ALLOW_EXF */


c   *------ Seaice Drag ------*
c   Fs = 1/2 * rho_sice * C_sice * A_sice * |v_sice - v_i| (v_sice - v_i)
c   implicit: rhss = 1/2* rho_sice * C_sice * A_sice + |v_ice - v_i[n]| * v_ice
C A_sice = L*h = length of iceberg * seaice thickness
C length of iceberg = sqrt(iceberg_area * dxG**2)
C seaice thickness = HEFF * AREA

C TODO fix problems with MDSIO_BUFF_3D
#ifdef ALLOW_SEAICE
      IF (useSEAICE) THEN

           DO j=1-OLy+1,sNy+OLy
             DO i=1-OLx+1,sNx+OLx

               dus = 0.5 * ((UICE(i,j,bi,bj)-utemp(i,j,bi,bj))
     &          + (UICE(i-1,j,bi,bj)-utemp(i-1,j,bi,bj)))
               dvs = 0.5 * ((VICE(i,j,bi,bj)-vtemp(i,j,bi,bj))
     &          + (VICE(i,j-1,bi,bj)-vtemp(i,j-1,bi,bj)))
C Calculate Norm on C-Point with constant values for wind stress
               CDs(i,j) = 0.5 * ICEBERG_rhoSeaice
     &            * ICEBERG_SeaiceDrag
     &            * sqrt(iceberg_area(i,j,bi,bj,nCl)*dxG(i,j,bi,bj)**2) !ändern zur Fläche
     &            * HEFF(i,j,bi,bj)/(AREA(i,j,bi,bj)+ seaice_area_reg)
     &            * sqrt(dus**2+dvs**2)
             ENDDO
           ENDDO

           DO j=1-OLy+1,sNy+OLy
              DO i=1-OLx+1,sNx+OLx

                IF (seaiceMaskU(i,j,bi,bj) .EQ. 0.0) THEN
                  Fsu(i,j,bi,bj) = 0. _d 0
                ELSE
                  Fsu(i,j,bi,bj) = 0.5 * ( CDs(i,j) + CDs(i-1,j) )
     &             * UICE(i,j,bi,bj)
                ENDIF

                IF (seaiceMaskV(i,j,bi,bj) .EQ. 0.0) THEN
                  Fsv(i,j,bi,bj) = 0. _d 0
                ELSE
                  Fsv(i,j,bi,bj) = 0.5 * ( CDs(i,j) + CDs(i,j-1) )
     &             * VICE(i,j,bi,bj)
                ENDIF
             ENDDO
           ENDDO
      ENDIF

#endif /*ALLOW_SEAICE*/


c   *------ Wave Radiation Force ------*
c   Fr = 1/4 * rho_ocean * gravity * a * L * v_air/|v_air|
C TODO What happens when UWIND and VWIND = 0? Assumption = No waves????
C constanten als eigene Variable for schleife

           DO j=0,sNy+1
              DO i=0,sNx+1
                norm_a = sqrt((0.5*(UWIND(i,j,bi,bj)
     &                      + UWIND(i+1,j,bi,bj)))**2
     &                      + (0.5*(VWIND(i,j,bi,bj)
     &                      + VWIND(i,j+1,bi,bj)))**2)
                Cdr(i,j) = 0.25 * ICEBERG_rhoOcean
     &            * ICEBERG_gravity
     &            * (0.5 * 0.02025)**2 *norm_a !quadrat prüfen
     &            * sqrt(iceberg_area(i,j,bi,bj,nCl)*dxG(i,j,bi,bj)**2)
              ENDDO
           ENDDO

           DO j=1,sNy+1
              DO i=1,sNx+1

                Fru(i,j,bi,bj) = 0.5* ( CDr(i,j)+ CDr(i-1,j))
     &            * UWIND(i,j,bi,bj)

                Frv(i,j,bi,bj) = 0.5 * ( CDr(i,j) + CDr(i,j-1))
     &            * UWIND(i,j,bi,bj)
              ENDDO
            ENDDO



c   *------ Pressure Gradient Force ------*
c   Fp = -M * g * Nabla * sea surface height
C QUESTION tilt is phiSurf in seaice Package, do i have to use phiSurf oder EtaN?

           DO j=1,sNy+1
             DO i=1,sNx+1
                Fpu(i,j,bi,bj) =
     &           iceberg_Massc(i,j,bi,bj) * ICEBERG_gravity
     &         * etaN(i,j,bi,bj)
     &         - iceberg_MassC(i-1,j,bi,bj)* ICEBERG_gravity
     &         * etaN(i-1,j,bi,bj)

                Fpv(i,j,bi,bj) =
     &           iceberg_MassC(i,j,bi,bj) * ICEBERG_gravity
     &         * etaN(i,j,bi,bj)
     &         - iceberg_MassC(i-1,j,bi,bj) * ICEBERG_gravity
     &         * etaN(i-1,j,bi,bj)
             ENDDO
           ENDDO



C TODO uncomment this section when all parts are compilable
C   *------ Momentum Equation ------*
C   dv_i/dt = Coriolis + Fw + Fa + Fs + Fr + Fp
            DO j=1,sNy+1
              DO i=1,sNx+1
                urhs(i,j,bi,bj) = Fcu(i,j,bi,bj) + Fwu(i,j,bi,bj)
     &              + Fau(i,j,bi,bj) + Fsu(i,j,bi,bj)
     &              + Fru(i,j,bi,bj) + Fpu(i,j,bi,bj)
                vrhs(i,j,bi,bj) = Fcv(i,j,bi,bj) + Fwv(i,j,bi,bj)
     &              + Fav(i,j,bi,bj) + Fsv(i,j,bi,bj)
     &              + Frv(i,j,bi,bj) + Fpv(i,j,bi,bj)
              ENDDO
            ENDDO



C   *------ Euler Forward Step ------*
c   v_i[n+1] = v_i[n] + dv_i * dt
C Calculated at C-Point
            DO j=1,sNy+1
             DO i=1,sNx+1
c               IF(iceberg_MaskU(i,j,bi,bj) .EQ. 0.0) THEN
c                  uc(i,j,bi,bj) = 0. _d 0
c                ELSE
                  uc(i,j,bi,bj) = iceberg_MaskU(i,j,bi,bj) *
     &              (iceberg_MassU(i,j,bi,bj) * utemp(i,j,bi,bj)
     &              + urhs(i,j,bi,bj) * ICEBERG_deltaTtherm)
     &              / (1. - iceberg_MaskU(i,j,bi,bj)
     &              + iceberg_MassU(i,j,bi,bj)
     &              + ICEBERG_deltaTtherm * (
     &              0.5 * ( CDw(i,j) + CDw(i-1,j) )
     &              + 0.5 * ( CDs(i,j) + CDs(i-1,j) )
     &               ))
c                ENDIF
c                IF(iceberg_MaskV(i,j,bi,bj) .EQ. 0.0) THEN
c                  vc(i,j,bi,bj) = 0. _d 0
c                ELSE
                 vc(i,j,bi,bj) = iceberg_MaskV(i,j,bi,bj) *
     &              (iceberg_MassV(i,j,bi,bj) * vtemp(i,j,bi,bj)
     &              + vrhs(i,j,bi,bj)* ICEBERG_deltaTtherm)
     &              / (1. - iceberg_MaskV(i,j,bi,bj)
     &              + iceberg_MassV(i,j,bi,bj)
     &              + ICEBERG_deltaTtherm * (
     &              0.5 * ( CDw(i,j) + CDw(i,j-1) )
     &              + 0.5 * ( CDs(i,j) + CDs(i,j-1) )
     &               ))
c                ENDIF

               ENDDO
             ENDDO

          ENDDO
        ENDDO



      RETURN
      END
