#include "ICEBERG_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

CBOP
C !ROUTINE: ICEBERG_ADVDIFF

C !INTERFACE: ==========================================================
      SUBROUTINE ICEBERG_MOM_EQ(
     I                  myTime, myIter, myThid )

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

C !INPUT/OUTPUT PARAMETERS: ===================================================
C     === Routine arguments ===
C     uc/vc     :: current ice velocity on C-grid;
C     myTime    :: current time in simulation
C     myIter    :: current iteration number in simulation
C     myThid    :: my Thread Id number
      _RL uc   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vc   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL myTime
      INTEGER myIter
      INTEGER myThid

C !LOCAL VARIABLES: ====================================================
C     === Local variables ===
C     i,j,bi,bj :: Loop counters
C     iCl        :: Loop counter for iceberg size classes

      INTEGER i, j, bi, bj, iCl
      _RL u_temp(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL v_temp(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

CEOP

c   *------- Coriolis Force -------*
C  -M * f * k x v_i




c   *------- Ocean Drag -------*
C   Fw = 1/2 * rho_ocean * C_ocean * A_ocean * |v_ocean - v_i| (v_ocean -v_i)

c   *------ Atmosphere Drag ------*
C   Fa = 1/2 * rho_atm * C_atm * A_atm * |v_atm - v_i| (v_atm -v_i)

c   *------ Seaice Drag ------*
c   Fs = 1/2 * rho_sice * C_sice * A_sice * |v_sice - v_i| (v_sice - v_i


c   *------ Wave Radiation Force ------*
c   Fr = 1/4 * rho_ocean * gravity * a * L v_air/|v_air|

c   *------ Pressure Gradient Force ------*
c   Fp = -M * Nabla P/rho_ocean

C   *------ Momentum Equation ------*
C   dv_i/dt = Coriolis + Fw + Fa + Fs + Fr + Fp


C   *------ Euler Forward Step ------*
c   v_i[n+1] = v_i[n] + dv_i * dt




      END