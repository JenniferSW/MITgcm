#ifdef ALLOW_ICEBERG

C     Package flag
c      LOGICAL myPa_MNC
c      LOGICAL myPa_MDSIO
c      COMMON /ICEBERG_PACKAGE/
c     &                     myPa_MNC, myPa_MDSIO



C     ICEBERG parameters
/* TODO sort parameters into ICEBERG and additional parameters */
c     ICEBERG_ITER0   :: iteration number to start iceberg drift
c     ICEBERG_numClUsed :: Number of icebergs used in the simulation
c     ICEBERG_deltaTtherm :: Seaice timestep for thermodynamic equations (s)
c     ICEBERGadvScheme  ::  Scheme for calculating the advection

      INTEGER ICEBERG_Iter0
      INTEGER ICEBERG_numClUsed
      INTEGER ICEBERGadvScheme


C-    additional parameters:
c     ICEBERG_ClLbl   :: Labels for diagnostics output
C
      _RL ICEBERG_deltaTtherm
      CHARACTER*(8) ICEBERG_ClLbl(ICEBERG_numCl)


/* TODO read in size class borders for each class that can be applied by class number
     Best to have a left and a right border and if one is not read it is the border
     of the previous of following size class or a default minimum or maximum
*/

C-    file names for initial conditions:
C     ICEBERG_initialFile   :: File with initial conditions for icebergs

      CHARACTER*(MAX_LEN_FNAM) ICEBERG_initialFile(ICEBERG_numCl)

c      COMMON /ICEBERG_PARAMS_L/

      COMMON /ICEBERG_PARAMS_I/ ICEBERGadvScheme,
     &       ICEBERG_Iter0,
     &       ICEBERG_numClUsed
      COMMON /ICEBERG_PARAMS_R/ ICEBERG_deltaTtherm

      COMMON /ICEBERG_PARAMS_C/ ICEBERG_ClLbl,
     &       ICEBERG_initialFile

#ifdef ICEBERG_3D_STATE
C     ICEBERG (3D) fields

C- Common Block for 3D fields
c      COMMON /ICEBERG_STATE_3D/


#endif /* ICEBERG_3D_STATE */

#ifdef ICEBERG_2D_STATE
C     Iceberg 2-dim. fields
C     iceberg_distr    :: field for iceberg distribution

      _RL iceberg_distr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,
     &  	      nSy,ICEBERG_numCl)
      COMMON /ICEBERG_STATE_2D/
     &    iceberg_distr

#endif /* ICEBERG_2D_STATE */

#ifdef ICEBERG_TENDENCY
c Tendency fields

C      COMMON /ICEBERG_TENDENCY/

#endif /* ICEBERG_TENDENCY */

#endif /* ALLOW_ICEBERG */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
