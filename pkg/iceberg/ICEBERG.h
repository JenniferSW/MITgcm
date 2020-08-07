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
C
c     ICEBERG_mass  ::  mass of iceberg in cell
c     ICEBERG_waveampl ::  wave amplitude (calculated f)

      INTEGER ICEBERG_Iter0
      INTEGER ICEBERG_numClUsed
      INTEGER ICEBERGadvScheme

      _RL ICEBERG_rhoAir
      _RL ICEBERG_rhoOcean
      _RL ICEBERG_rhoSeaice
      _RL ICEBERG_AirDrag
      _RL ICEBERG_OceanDrag
      _RL ICEBERG_SeaiceDrag
      _RL ICEBERG_gravity

/* TODO prove whether these have to be fields or not */
      _RL ICEBERG_mass
      _RL ICEBERG_waveampl



/* TODO see seaice package for wind or air velocity */
      _RL WINDX      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL WINDY      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)



C-    additional parameters:
c     ICEBERG_deltaTtherm :: Seaice timestep for thermodynamic equations (s)
c     ICEBERG_ClLbl   :: Labels for diagnostics output
C
      _RL ICEBERG_deltaTtherm

      CHARACTER*(8) ICEBERG_ClLbl(ICEBERG_numCl)




C-    file names for initial conditions:
C     ICEBERG_initialArea   :: File with initial fractional iceberg cover at cell centers
C     ICEBERG_initialHeight :: File with initial mean iceberg height at cell centers

      CHARACTER*(MAX_LEN_FNAM) ICEBERG_initialArea(ICEBERG_numCl)
      CHARACTER*(MAX_LEN_FNAM) ICEBERG_initialHeight(ICEBERG_numCl)

c      COMMON /ICEBERG_PARAMS_L/

C Common Block for integer parameters

      COMMON /ICEBERG_PARAMS_I/ ICEBERGadvScheme,
     &       ICEBERG_Iter0,
     &       ICEBERG_numClUsed

C Common block for real parameters
      COMMON /ICEBERG_PARAMS_R/ ICEBERG_deltaTtherm,
     &       ICEBERG_rhoAir, ICEBERG_rhoOcean, ICEBERG_rhoSeaice,
     &       ICEBERG_AirDrag, ICEBERG_OceanDrag, ICEBERG_SeaiceDrag,
     &       ICEBERG_gravity, ICEBERG_mass, ICEBERG_waveampl,
     &       WINDX, WINDY


C Common Block for Character Variables and file names
      COMMON /ICEBERG_PARAMS_C/ ICEBERG_ClLbl,
     &       ICEBERG_initialArea, ICEBERG_initialHeight

#ifdef ICEBERG_3D_STATE
C     ICEBERG (3D) fields

C- Common Block for 3D fields
c      COMMON /ICEBERG_STATE_3D/


#endif /* ICEBERG_3D_STATE */

/* QUESTION kann ifdef 2D STate weg? */
#ifdef ICEBERG_2D_STATE
C     Iceberg 2-dim. fields
C     iceberg_distr    :: field for iceberg distribution/effective height (m)
C                         - actual height is iceberg_distr/iceberg_area
C     iceberg_area     ::  field for fractional area covered by icebergs
C                         - here 0 is no cover and 1 is 100% cover
C     iceberg_height  ::  field for mean iceberg height in grid cell (m)
C     iceberg_MassC   ::  iceberg Mass at cell center
C     iceberg_MassU   ::  iceberg Mass at U grid point
C     iceberg_MassV   ::  iceberg Mass at V grid point

      _RL iceberg_distr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,
     &  	      nSy,ICEBERG_numCl)
      _RL iceberg_area(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,
     &  	      nSy,ICEBERG_numCl)
      _RL iceberg_height(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,
     &  	      nSy,ICEBERG_numCl)
      _RL iceberg_MassC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,
     &  	      nSy)
      _RL iceberg_MassU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,
     &  	      nSy)
      _RL iceberg_MassV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,
     &  	      nSy)


      COMMON /ICEBERG_STATE_2D/
     &    iceberg_distr,
     &    iceberg_area, iceberg_height

#endif /* ICEBERG_2D_STATE */

#ifdef ICEBERG_TENDENCY
c Tendency fields

C      COMMON /ICEBERG_TENDENCY/

#endif /* ICEBERG_TENDENCY */

#endif /* ALLOW_ICEBERG */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
