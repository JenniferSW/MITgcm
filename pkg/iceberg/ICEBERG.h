#ifdef ALLOW_ICEBERG

C     Package flag
      LOGICAL myPa_MNC
      LOGICAL myPa_MDSIO
      COMMON /MYPA_ICEBERG/
     &                     myPa_MNC, myPa_MDSIO



C     ICEBERG parameters
      LOGICAL myPa_StaV_Cgrid
      LOGICAL myPa_Tend_Cgrid
      LOGICAL myPa_applyTendT
      LOGICAL myPa_applyTendS
      LOGICAL myPa_applyTendU
      LOGICAL myPa_applyTendV

C-    additional parameters:
c     ICEBERG_ITER0   :: iteration number to start iceberg drift
c     ICEBERG_numClUsed :: Number of icebergs used in the simulation
c     ICEBERG_ClLbl   :: Labels for diagnostics output

      LOGICAL myPa_doSwitch1
      LOGICAL myPa_doSwitch2
      INTEGER ICEBERG_Iter0
      INTEGER myPa_index2
      INTEGER ICEBERG_numClUsed
      _RL myPa_param1
      _RL myPa_param2
      CHARACTER*(8) ICEBERG_ClLbl(ICEBERG_numCl)
      CHARACTER*(MAX_LEN_FNAM) myPa_string2

C-    file names for initial conditions:
C     ICEBERG_initialFile   :: File with initial conditions for icebergs

      CHARACTER*(MAX_LEN_FNAM) ICEBERG_initialFile(ICEBERG_numCl)
      CHARACTER*(MAX_LEN_FNAM) myPa_Scal2File
      CHARACTER*(MAX_LEN_FNAM) myPa_VelUFile
      CHARACTER*(MAX_LEN_FNAM) myPa_VelVFile
      CHARACTER*(MAX_LEN_FNAM) myPa_Surf1File
      CHARACTER*(MAX_LEN_FNAM) myPa_Surf2File

      COMMON /MYPA_PARAMS_L/
     &       myPa_StaV_Cgrid, myPa_Tend_Cgrid,
     &       myPa_applyTendT, myPa_applyTendS,
     &       myPa_applyTendU, myPa_applyTendV,
     &       myPa_doSwitch1, myPa_doSwitch2
      COMMON /MYPA_PARAMS_I/ ICEBERG_Iter0, myPa_index2,
     &       ICEBERG_numClUsed
      COMMON /MYPA_PARAMS_R/ myPa_param1, myPa_param2
      COMMON /MYPA_PARAMS_C/ ICEBERG_ClLbl, myPa_string2,
     &       ICEBERG_initialFile, myPa_Scal2File,
     &       myPa_VelUFile,  myPa_VelVFile,
     &       myPa_Surf1File, myPa_Surf2File

#ifdef ICEBERG_3D_STATE
C     ICEBERG (3D) fields
C     iceberg_distr    :: field for iceberg distribution

      _RL iceberg_distr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,
     &  	      nSy,ICEBERG_numCl)
      _RL myPa_StatScal2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL myPa_StatVelU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL myPa_StatVelV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      COMMON /MYPA_STATE_3D/
     &    iceberg_distr, myPa_StatScal2,
     &    myPa_StatVelU,  myPa_StatVelV
#endif /* ICEBERG_3D_STATE */
#ifdef ICEBERG_2D_STATE
C     MYPA 2-dim. fields
      _RL myPa_Surf1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL myPa_Surf2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /MYPA_STATE_2D/
     &    myPa_Surf1, myPa_Surf2
#endif /* ICEBERG_2D_STATE */

#ifdef ICEBERG_TENDENCY
      _RL myPa_TendScal1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL myPa_TendScal2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL myPa_TendVelU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL myPa_TendVelV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      COMMON /MYPA_TENDENCY/
     &    myPa_TendScal1, mypa_TendScal2,
     &    myPa_TendVelU,  mypa_TendVelV
#endif /* ICEBERG_TENDENCY */

#endif /* ALLOW_ICEBERG */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
