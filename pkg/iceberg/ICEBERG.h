#ifdef ALLOW_ICEBERG

C     Package flag
c      LOGICAL myPa_MNC
c      LOGICAL myPa_MDSIO
c      COMMON /ICEBERG_PACKAGE/
c     &                     myPa_MNC, myPa_MDSIO



C     ICEBERG parameters



C-    additional parameters:
c     ICEBERG_ITER0   :: iteration number to start iceberg drift
c     ICEBERG_numClUsed :: Number of icebergs used in the simulation
c     ICEBERG_ClLbl   :: Labels for diagnostics output

      INTEGER ICEBERG_Iter0
      INTEGER ICEBERG_numClUsed
      CHARACTER*(8) ICEBERG_ClLbl(ICEBERG_numCl)

C-    file names for initial conditions:
C     ICEBERG_initialFile   :: File with initial conditions for icebergs

      CHARACTER*(MAX_LEN_FNAM) ICEBERG_initialFile(ICEBERG_numCl)

C- Common Block for Logical Parameter
c      COMMON /ICEBERG_PARAMS_L/
c     &       myPa_StaV_Cgrid, myPa_Tend_Cgrid,

C- Common Block for Integer, Real and Character
      COMMON /ICEBERG_PARAMS_I/ ICEBERG_Iter0,
     &       ICEBERG_numClUsed
c      COMMON /ICEBERG_PARAMS_R/
      COMMON /ICEBERG_PARAMS_C/ ICEBERG_ClLbl,
     &       ICEBERG_initialFile

#ifdef ICEBERG_3D_STATE
C     ICEBERG (3D) fields
C     iceberg_distr    :: field for iceberg distribution

      _RL iceberg_distr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,
     &  	      nSy,ICEBERG_numCl)

C- Common Block for 3D fields
      COMMON /ICEBERG_STATE_3D/
     &    iceberg_distr

#endif /* ICEBERG_3D_STATE */
#ifdef ICEBERG_2D_STATE
C     Iceberg 2-dim. fields

c      COMMON /ICEBERG_STATE_2D/

#endif /* ICEBERG_2D_STATE */

#ifdef ICEBERG_TENDENCY
c Tendency fields

C      COMMON /ICEBERG_TENDENCY/

#endif /* ICEBERG_TENDENCY */

#endif /* ALLOW_ICEBERG */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
