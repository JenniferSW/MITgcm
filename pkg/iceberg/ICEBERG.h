#ifdef ALLOW_ICEBERG

C     Package flag
c      LOGICAL myPa_MNC
c      LOGICAL myPa_MDSIO
c      COMMON /ICEBERG_PACKAGE/
c     &                     myPa_MNC, myPa_MDSIO



C     ICEBERG parameters
c     ICEBERG_ITER0   :: iteration number to start iceberg drift
c     ICEBERG_numClUsed :: Number of icebergs used in the simulation

c     ICEBERGadvScheme  ::  Scheme for calculating the advection
C     ICEBERG_rhoAir  :: density of the air
c     ICEBERG_rhoOcean  :: density of sea water
c     ICEBERG_rhoSeaice  ::  density of sea ice
c     ICEBERG_AirDrag  ::  air drag coefficient
c     ICEBERG_OceanDrag ::  ocean drag coefficient
c     ICEBERG_SeaiceDrag  ::  sea ice drag coefficient
c     ICEEBERG_gravity  ::  gravitational acceleration
c     ICEBERG_temp      :: iceberg temperature (default - 4 deg C)
C

      INTEGER ICEBERG_Iter0
      INTEGER ICEBERG_numClUsed
      INTEGER ICEBERGadvScheme

      _RL ICEBERG_rhoAir
      _RL ICEBERG_rhoOcean
      _RL ICEBERG_rhoSeaice
      _RL ICEBERG_rhoIceberg
      _RL ICEBERG_AirDrag
      _RL ICEBERG_AirDragform
      _RL ICEBERG_OceanDrag
      _RL ICEBERG_OceanDragform
      _RL ICEBERG_SeaiceDrag
      _RL ICEBERG_gravity
      _RL ICEBERG_temp






C-    additional parameters:
c     ICEBERG_deltaT :: Seaice timestep for thermodynamic equations (s)
c     ICEBERG_ClLbl   :: Labels for diagnostics output
C
      _RL ICEBERG_deltaT

      CHARACTER*(8) ICEBERG_ClLbl(ICEBERG_numCl)




C-    file names for initial conditions:
C     ICEBERG_initConcentration   :: File with initial concentration for each size class

      CHARACTER*(MAX_LEN_FNAM) ICEBERG_initConcentration(ICEBERG_numCl)
      CHARACTER*(MAX_LEN_FNAM) ICEBERG_CalvFile(ICEBERG_numCl)

c      COMMON /ICEBERG_PARAMS_L/

C Common Block for integer parameters

      COMMON /ICEBERG_PARAMS_I/ ICEBERGadvScheme,
     &     ICEBERG_Iter0,
     &     ICEBERG_numClUsed

C Common block for real parameters
      COMMON /ICEBERG_PARAMS_R/ ICEBERG_deltaT,
     &     ICEBERG_rhoAir, ICEBERG_rhoOcean,
     &     ICEBERG_rhoSeaice, ICEBERG_rhoIceberg,
     &     ICEBERG_AirDrag, ICEBERG_AirDragform,
     &     ICEBERG_OceanDrag, ICEBERG_OceanDragform,
     &     ICEBERG_SeaiceDrag,
     &     ICEBERG_gravity, ICEBERG_temp



C Common Block for Character Variables and file names
      COMMON /ICEBERG_PARAMS_C/ ICEBERG_ClLbl,
     &       ICEBERG_initConcentration,
     &       ICEBERG_CalvFile

#ifdef ICEBERG_3D_STATE
C     ICEBERG (3D) fields

C- Common Block for 3D fields
c      COMMON /ICEBERG_STATE_3D/


#endif /* ICEBERG_3D_STATE */

/* QUESTION kann ifdef 2D STate weg? */
#ifdef ICEBERG_2D_STATE
C     Iceberg 2-dim. fields
C     iceberg_concentration :: field for iceberg concentration
C                         - number of icebergs per grid point
C     ICEBERG_width   ::  field for iceberg width for each size class
C     ICEBERG_height  ::  list of iceberg height for each size class
C     iceberg_Vsgl :: iceberg volume for one iceberg in each size class
C     iceberg_MassC   ::  iceberg Mass at cell center
C     iceberg_MassU   ::  iceberg Mass at U grid point
C     iceberg_MassV   ::  iceberg Mass at V grid point
C     iceberg_MaskC   ::  iceberg Mask to identify cell with icebergs and without (center)
C     iceberg_MaskU   ::  iceberg Mask to identify cell with and without icebergs (U point)
C     iceberg_MaskV   ::  iceberg Mask to identify cell with and without icebergs (V point)
C     iceberg_Vloss   ::  total loss of volume in a time step for every size class

/* TODO die Listen in eine eigene Kategorie fassen */

      _RL iceberg_concentration(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,
     &     ICEBERG_numCl)
      _RL iceberg_CalvRate(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,
     &     ICEBERG_numCl)
      _RL ICEBERG_width (ICEBERG_numCl)
      _RL iceberg_length (ICEBERG_numCl)
      _RL ICEBERG_height (ICEBERG_numCl)
      _RL iceberg_Vsgl(ICEBERG_numCl)

      _RL iceberg_MassC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL iceberg_MassU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL iceberg_MassV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL iceberg_MaskC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL iceberg_MaskU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL iceberg_MaskV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL iceberg_uvel(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,
     &     ICEBERG_numCl)
      _RL iceberg_vvel(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,
     &     ICEBERG_numCl)
      _RL iceberg_Vloss(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,
     &     ICEBERG_numCl)

      COMMON /ICEBERG_STATE_2D/
     &    iceberg_concentration, iceberg_CalvRate,
     &    ICEBERG_height,ICEBERG_width, iceberg_length,
     &    iceberg_Vsgl,
     &    iceberg_MassC, iceberg_MassU, iceberg_MassV,
     &    iceberg_MaskC, iceberg_MaskU, iceberg_MaskV,
     &    iceberg_uvel, iceberg_vvel,
     &    iceberg_Vloss
#endif /* ICEBERG_2D_STATE */

#ifdef ICEBERG_TENDENCY
c Tendency fields

C      COMMON /ICEBERG_TENDENCY/

#endif /* ICEBERG_TENDENCY */

#endif /* ALLOW_ICEBERG */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
