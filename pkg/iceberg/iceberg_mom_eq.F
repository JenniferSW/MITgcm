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
      SUBROUTINE ICEBERG_MOM_EQ( nCl,
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
#include "SURFACE.h"
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
C     nCl       :: number of class the velocity is calculated for
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
      _RL CFL_test
      _RL duw , dvw
      _RL dua , dva
      _RL dus , dvs
      _RL norm_a, Lw
      _RL cds_const, cdr_const
      _RL form_a, form_w
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
      _RL cwave(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL phiSurf(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL massC(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL massU(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL massV(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

CEOP



        DO bj=myByLo(myThid),myByHi(myThid)
          DO bi=myBxLo(myThid),myBxHi(myThid)

           DO j = 1-OLy, sNy+OLy
            DO i = 1-OLx, sNx+OLx
              utemp(i,j,bi,bj) = iceberg_uvel(i,j,bi,bj,nCl)
              vtemp(i,j,bi,bj) = iceberg_vvel(i,j,bi,bj,nCl)

              phiSurf(i,j) = Bo_surf(i,j,bi,bj)*etaN(i,j,bi,bj) ! calculate sea surface height

              Fcu(i,j,bi,bj) = 0. _d 0
              Fcv(i,j,bi,bj) = 0. _d 0
              Fwu(i,j,bi,bj) = 0. _d 0
              Fwv(i,j,bi,bj) = 0. _d 0
              Fau(i,j,bi,bj) = 0. _d 0
              Fav(i,j,bi,bj) = 0. _d 0
              Fsu(i,j,bi,bj) = 0. _d 0
              Fsv(i,j,bi,bj) = 0. _d 0
              Fru(i,j,bi,bj) = 0. _d 0
              Frv(i,j,bi,bj) = 0. _d 0
              Fpu(i,j,bi,bj) = 0. _d 0
              Fpv(i,j,bi,bj) = 0. _d 0
              urhs(i,j,bi,bj) = 0. _d 0
              vrhs(i,j,bi,bj) = 0. _d 0
              CDw(i,j) = 0. _d 0 ! initialise coefficients
              CDa(i,j) = 0. _d 0
              CDs(i,j) = 0. _d 0
              CDr(i,j) = 0. _d 0
              cwave(i,j) = 0. _d 0
              massC(i,j) = 0. _d 0
              massU(i,j) = 0. _d 0
              massV(i,j) = 0. _d 0
            ENDDO
           ENDDO

C Set up mass per unit area
C interpolation from c to u/v point
C u-point -> 0.5*( c(i,j) + c(i-1,j))
C v-point -> 0.5*( c(i,j) + c(i,j-1))

C This is the mass of all icebergs at grid point
!           DO j=1-OLy+1,sNy+OLy
!             DO i=1-OLx+1,sNx+OLx
!               iceberg_MassC(i,j,bi,bj) =
!    &                 ICEBERG_rhoSeaice
!    &                 * iceberg_concentration(i,j,bi,bj,nCl)
!    &                 * ICEBERG_width(nCl) * iceberg_length(nCl)
!    &                 * ICEBERG_height(nCl)
!               iceberg_MassU(i,j,bi,bj) = ICEBERG_rhoSeaice * 0.5 *
!    &                 (iceberg_concentration(i,j,bi,bj,nCl)
!    &                 + iceberg_concentration(i-1,j,bi,bj,nCl))
!    &                 * ICEBERG_width(nCl) * iceberg_length(nCl)
!    &                 * ICEBERG_height(nCl)
!               iceberg_MassV(i,j,bi,bj) = ICEBERG_rhoSeaice * 0.5 *
!    &                 (iceberg_concentration(i,j,bi,bj,nCl)
!    &                 + iceberg_concentration(i,j-1,bi,bj,nCl))
!    &                 * ICEBERG_width(nCl) * iceberg_length(nCl)
!    &                 * ICEBERG_height(nCl)
!             ENDDO
!           ENDDO
!           _EXCH_XY_RL( iceberg_MassC, mythid )
!           _EXCH_XY_RL( iceberg_MassU, mythid )
!           _EXCH_XY_RL( iceberg_MassV, mythid )
            DO j=-1,sNy+1
              DO i=-1,sNx+1
                massC(i,j) =
     &                 ICEBERG_rhoSeaice
     &                 * iceberg_concentration(i,j,bi,bj,nCl)
     &                 * ICEBERG_width(nCl) * iceberg_length(nCl)
     &                 * ICEBERG_height(nCl)
                massU(i,j) = ICEBERG_rhoSeaice * 0.5 *
     &                 (iceberg_concentration(i,j,bi,bj,nCl)
     &                 + iceberg_concentration(i-1,j,bi,bj,nCl))
     &                 * ICEBERG_width(nCl) * iceberg_length(nCl)
     &                 * ICEBERG_height(nCl)
                massV(i,j) = ICEBERG_rhoSeaice * 0.5 *
     &                 (iceberg_concentration(i,j,bi,bj,nCl)
     &                 + iceberg_concentration(i,j-1,bi,bj,nCl))
     &                 * ICEBERG_width(nCl) * iceberg_length(nCl)
     &                 * ICEBERG_height(nCl)
              ENDDO
            ENDDO

c   *------- Coriolis Force -------*
C  Fc = -M * f * k x v_i[n]
C _fCori  ::  Makro from model for Coriolis parameter at cell center
C iceberg_mass and _fcori are already on C-points, only iceberg velocity is interpolated

#ifdef ALLOW_ICEBERG_CORIOLIS
           DO j=0,sNy+1
             DO i=0,sNx+1


                Fcu(i,j,bi,bj)=
     &            0.5 _d 0 * (
     &            massC(i  ,j) * _fCori(i  ,j,bi,bj)
     &            * 0.5 _d 0*( vtemp(i  ,j,bi,bj)+vtemp(i  ,j+1,bi,bj) )
     &            + massC(i-1,j) * _fCori(i-1,j,bi,bj)
     &            * 0.5 _d 0*( vtemp(i-1,j,bi,bj)+vtemp(i-1,j+1,bi,bj) )
     &            )
                Fcv(i,j,bi,bj)=
     &            - 0.5 _d 0*(
     &            massC(i,j) * _fCori(i,j  ,bi,bj)
     &            * 0.5 _d 0*( utemp(i,j  ,bi,bj)+utemp(i+1,  j,bi,bj) )
     &            + massC(i,j-1) * _fCori(i,j-1,bi,bj)
     &            * 0.5 _d 0*( utemp(i,j-1,bi,bj)+utemp(i+1,j-1,bi,bj) )
     &            )
             ENDDO
           ENDDO
#endif
c   *------- Ocean Drag -------*
C   explicit: Fw = 1/2 * rho_ocean * C_ocean * A_ocean * |v_ocean - v_i| (v_ocean -v_i)
c   implicit: rhsw = 1/2* rho_ocean * C_ocean * A_ocean + |v_ocean - v_i[n]| * v_ocean


#ifdef ALLOW_ICEBERG_OCEAN_DRAG
            DO j=-1,sNy+1
              DO i=-1,sNx+1
               duw = 0.5 _d 0 * ((uVel(i,j,1,bi,bj)-utemp(i,j,bi,bj))
     &          + (uVel(i+1,j,1,bi,bj)-utemp(i+1,j,bi,bj)))
               dvw = 0.5 _d 0 * ((vVel(i,j,1,bi,bj)-vtemp(i,j,bi,bj))
     &          +(vVel(i,j+1,1,bi,bj)-vtemp(i,j+1,bi,bj)))

c Calculate Norm on C-Point with constants for ocean stress
                IF (iceberg_MaskC(i,j,bi,bj) .EQ. 0. _d 0 ) THEN
                  CDw(i,j) = 0. _d 0
                ELSE
                 form_w = 0.5 _d 0 * ICEBERG_OceanDragform
     &              * ICEBERG_rhoIceberg/ICEBERG_rhoOcean
     &              * ICEBERG_height(nCl) * ICEBERG_width(nCl)

                 CDw(i,j) = ICEBERG_rhoOcean *
     &              (form_w
     &              + ICEBERG_OceanDrag
     &              * iceberg_length(nCl)*iceberg_width(nCl))
     &              * sqrt(duw**2+dvw**2)
                ENDIF
              ENDDO
            ENDDO
C calculate rhs of ocean drag force
            DO j=0,sNy+1
              DO i=0,sNx+1
                Fwu(i,j,bi,bj) = 0.5 * ( CDw(i,j) + CDw(i-1,j) )
     &            * uVel(i,j,1,bi,bj)

                Fwv(i,j,bi,bj) = 0.5 * ( CDw(i,j) + CDw(i,j-1) )
     &            * vVel(i,j,1,bi,bj)
              ENDDO
            ENDDO
#endif
c   *------ Atmosphere Drag ------*
C Uwind/Vwind is located on southwest C-Grid points

C   Fa = 1/2 * rho_atm * C_atm * A_atm * |v_atm - v_i| (v_atm -v_i)
C   Assumed that v_i is much smaller than v_atm (v_i << v_atm)
c   rhsa = 1/2 * rho_atm * C_atm * A_atm * | v_atm| * v_atm

#ifdef ALLOW_ICEBERG_ATMOSPHERE_DRAG
#ifdef ALLOW_EXF


            DO j=-1,sNy+1
              DO i=-1,sNx+1

C Calculate Norm on C-Point with constant values for wind stress
                IF (iceberg_MaskC(i,j,bi,bj) .EQ. 0. _d 0 ) THEN
                  CDa(i,j) = 0. _d 0
                ELSE
                  form_a = 0.5 * ICEBERG_AirDragform
     &              * (1-ICEBERG_rhoIceberg/ICEBERG_rhoOcean)
     &              * ICEBERG_height(nCl) * ICEBERG_width(nCl)

                  CDa(i,j) = ICEBERG_rhoAir *
     &              (form_a
     &              + ICEBERG_AirDrag
     &              * iceberg_length(nCl)*iceberg_width(nCl))
     &              * sqrt(UWIND(i,j,bi,bj)**2+VWIND(i,j,bi,bj)**2)
                ENDIF
                IF(massC(i,j) .NE. 0. _d 0)THEN
c                  PRINT *, 'debug forma', form_a
                ENDIF
              ENDDO
            ENDDO

            DO j=0,sNy+1
              DO i=0,sNx+1
                Fau(i,j,bi,bj) = 0.5 *
     &          ( CDa(i,j) * UWIND(i,j,bi,bj)
     &           + CDa(i-1,j)* UWIND(i-1,j,bi,bj))

                Fav(i,j,bi,bj) = 0.5 *
     &           ( CDa(i,j) * VWIND(i,j,bi,bj)
     &            + CDa(i,j-1) * VWIND(i,j-1,bi,bj))

              ENDDO
            ENDDO

#endif /* ALLOW_EXF */
#endif

c   *------ Seaice Drag ------*
c   Fs = 1/2 * rho_sice * C_sice * A_sice * |v_sice - v_i| (v_sice - v_i)
c   implicit: rhss = 1/2* rho_sice * C_sice * A_sice + |v_ice - v_i[n]| * v_ice
C A_sice = L*h = length of iceberg * seaice thickness
C seaice thickness = HEFF * AREA

C TODO fix problems with MDSIO_BUFF_3D, related to SEAICE_ITD
#ifdef ALLOW_ICEBERG_SEAICE_DRAG
#ifdef ALLOW_SEAICE
      IF (useSEAICE) THEN
           cds_const = 0.5 * ICEBERG_rhoSeaice
     &            * ICEBERG_SeaiceDrag
           DO j=-1,sNy+1
             DO i=-1,sNx+1

               dus = 0.5 * ((UICE(i,j,bi,bj)-utemp(i,j,bi,bj))
     &          + (UICE(i+1,j,bi,bj)-utemp(i+1,j,bi,bj)))
               dvs = 0.5 * ((VICE(i,j,bi,bj)-vtemp(i,j,bi,bj))
     &          + (VICE(i,j+1,bi,bj)-vtemp(i,j+1,bi,bj)))
C     Calculate Norm on C-Point with constant values for wind stress
C     AREA can be zero, so avoid division by regularisation stolen from seaice_growth
               CDs(i,j) = cds_const
     &            * iceberg_width(nCl)
     &            * HEFF(i,j,bi,bj)/(AREA(i,j,bi,bj)+ seaice_area_reg)
     &            * sqrt(dus**2+dvs**2)
             ENDDO
           ENDDO

           DO j=0,sNy+1
              DO i=0,sNx+1

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
#endif

c   *------ Wave Radiation Force ------*
c   Fr = 1/4 * rho_ocean * gravity * a^2 * L * v_air/|v_air|
C term conversion as shown in Rackow (2011)
#ifdef ALLOW_ICEBERG_WAVE_RADIATION
           cdr_const = 0.25 * ICEBERG_rhoOcean
     &            * ICEBERG_gravity * (0.5 * 0.02025)**2
            DO j=0,sNy+1
              DO i=0,sNx+1
                Lw = 0.32 *
     &             sqrt((UWIND(i,j,bi,bj)
     &             - 0.5*(uVel(i,j,1,bi,bj) + uVel(i-1,j,1,bi,bj)))**2
     &             + (VWIND(i,j,bi,bj)
     &             - 0.5*(vVel(i,j,1,bi,bj) + vVel(i,j-1,1,bi,bj)))**2)
     &              **2
                IF (uVel(i-1,j,1,bi,bj) .NE. 0. _d 0) THEN
                 PRINT*,'dbg uvel', uVel(i-1,j,1,bi,bj)
                ENDIF

                cwave(i,j) = 0.06 * MIN(MAX(0. _d 0,
     &           (iceberg_length(nCl)
     &            -0.125*Lw)/(0.25*Lw-0.125*Lw)),1. _d 0)
              ENDDO
            ENDDO


           DO j=0,sNy+1
              DO i=0,sNx+1
                norm_a = sqrt(UWIND(i,j,bi,bj)**2
     &                      +VWIND(i,j,bi,bj)**2)
               IF (iceberg_MaskC(i,j,bi,bj) .EQ. 0. _d 0 ) THEN
                Cdr(i,j) = 0. _d 0
               ELSE
                Cdr(i,j) = cdr_const * cwave(i,j)
     &            * (2*iceberg_length(nCl)*iceberg_width(nCl))
     &            /(iceberg_length(nCl)+iceberg_width(nCl))
     &            * norm_a**3
               ENDIF
              ENDDO
           ENDDO

           DO j=0,sNy+1
              DO i=0,sNx+1

                Fru(i,j,bi,bj) = 0.5 * ( CDr(i,j) + CDr(i-1,j))
     &            * 0.5 * ( UWIND(i,j,bi,bj) + UWIND(i-1,j,bi,bj) )

                Frv(i,j,bi,bj) = 0.5 * ( CDr(i,j) + CDr(i,j-1))
     &            * 0.5 * ( VWIND(i,j,bi,bj) + VWIND(i,j-1,bi,bj) )
              ENDDO
            ENDDO

#endif


c   *------ Pressure Gradient Force ------*
c   Fp = -M * g * Nabla * sea surface height
c   phiSurf = g * etaN
#ifdef ALLOW_ICEBERG_PRESSURE_GRADIENT
           DO j=0,sNy+1
             DO i=0,sNx+1
                Fpu(i,j,bi,bj) = - massU(i,j) *
     &           _recip_dxC(i,j,bi,bj) *
     &           (phiSurf(i,j)
     &         - phiSurf(i-1,j))

                Fpv(i,j,bi,bj) = - massV(i,j) *
     &          _recip_dyC(i,j,bi,bj) *
     &          (phiSurf(i,j)
     &         - phiSurf(i,j-1))
             ENDDO
           ENDDO

#endif

C   *------ Momentum Equation ------*
C   dv_i/dt = Coriolis + Fw + Fa + Fs + Fr + Fp
            DO j=0,sNy+1
              DO i=0,sNx+1
                urhs(i,j,bi,bj) = Fcu(i,j,bi,bj) + Fwu(i,j,bi,bj)
     &              + Fau(i,j,bi,bj) + Fsu(i,j,bi,bj)+ Fpu(i,j,bi,bj)
     &              + Fru(i,j,bi,bj)
                vrhs(i,j,bi,bj) = Fcv(i,j,bi,bj) + Fwv(i,j,bi,bj)
     &              + Fav(i,j,bi,bj) + Fsv(i,j,bi,bj)+ Fpv(i,j,bi,bj)
     &              + Frv(i,j,bi,bj)
              ENDDO
            ENDDO



C   *------ Euler Forward Step ------*
c   v_i[n+1] = v_i[n] + dv_i * dt
C Calculated at C-Point
            DO j=0,sNy+1
             DO i=0,sNx+1

                  iceberg_uvel(i,j,bi,bj,nCl) =
     &              iceberg_MaskU(i,j,bi,bj) * (
     &              (massU(i,j) * utemp(i,j,bi,bj)
     &              + urhs(i,j,bi,bj) * ICEBERG_deltaT)
     &              / (1. - iceberg_MaskU(i,j,bi,bj)
     &              + massU(i,j)
     &              + ICEBERG_deltaT * (
     &              0.5 * ( CDw(i,j) + CDw(i-1,j) )
     &              + 0.5 * ( CDs(i,j) + CDs(i-1,j) )
     &               )))

                 iceberg_vvel(i,j,bi,bj,nCl) =
     &              iceberg_MaskV(i,j,bi,bj) * (
     &              (massV(i,j) * vtemp(i,j,bi,bj)
     &              + vrhs(i,j,bi,bj) * ICEBERG_deltaT)
     &              / (1. - iceberg_MaskV(i,j,bi,bj)
     &              + massV(i,j)
     &              + ICEBERG_deltaT * (
     &              0.5 * ( CDw(i,j) + CDw(i,j-1) )
     &              + 0.5 * ( CDs(i,j) + CDs(i,j-1) )
     &               )))

               ENDDO
             ENDDO

          ENDDO
        ENDDO
        CALL EXCH_UV_XY_RL(iceberg_uvel(1-OLx,1-OLy,1,1,nCl),
     &    iceberg_vvel(1-OLx,1-OLy,1,1,nCl),
     &    .TRUE.,myThid)

      RETURN
      END
