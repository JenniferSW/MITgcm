#include "ICEBERG_OPTIONS.h"
#ifdef ALLOW_GENERIC_ADVDIFF
# include "GAD_OPTIONS.h"
#endif
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: ICEBERG_ADVECTION
C
C !INTERFACE: ==========================================================
      SUBROUTINE ICEBERG_ADVECTION(
     I     tracerIdentity,
     I     advectionScheme,
     I     uFld, vFld, uTrans, vTrans, iceFld, r_hFld,
     O     gFld, afx, afy,
     I     bi, bj, myTime, myIter, myThid)

C !DESCRIPTION:
C Calculates the tendency of a iceberg field due to advection.
C It uses the multi-dimensional method given in \ref{sect:multiDimAdvection}
C and can only be used for the non-linear advection schemes such as the
C direct-space-time method and flux-limiters.
C
C This routine is an adaption of the GAD_ADVECTION for 2D-fields.
C for Area, effective thickness or other "extensive" iceberg field,
C  the contribution iceFld*div(u) (that is present in gad_advection)
C  is not included here.
C
C The algorithm is as follows:
C \begin{itemize}
C \item{$\theta^{(n+1/2)} = \theta^{(n)}
C      - \Delta t \partial_x (u\theta^{(n)}) + \theta^{(n)} \partial_x u$}
C \item{$\theta^{(n+2/2)} = \theta^{(n+1/2)}
C      - \Delta t \partial_y (v\theta^{(n+1/2)}) + \theta^{(n)} \partial_y v$}
C \item{$G_\theta = ( \theta^{(n+2/2)} - \theta^{(n)} )/\Delta t$}
C \end{itemize}
C
C The tendency (output) is over-written by this routine.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"
#ifdef ALLOW_GENERIC_ADVDIFF
# include "GAD.h"
#endif
#ifdef ALLOW_AUTODIFF_TAMC
# include "tamc.h"
# include "tamc_keys.h"
# ifdef ALLOW_PTRACERS
#  include "PTRACERS_SIZE.h"
# endif
#endif
#ifdef ALLOW_EXCH2
#include "W2_EXCH2_SIZE.h"
#include "W2_EXCH2_TOPOLOGY.h"
#endif /* ALLOW_EXCH2 */
      LOGICAL extensiveFld
      PARAMETER ( extensiveFld = .TRUE. )

C !INPUT PARAMETERS: ===================================================
C  tracerIdentity  :: (input) tracer identifier
C  advectionScheme :: (input) advection scheme to use (Horizontal plane)
C  extensiveFld    :: indicates to advect an "extensive" type of ice field
C  uFld            :: (input) velocity, zonal component
C  vFld            :: (input) velocity, meridional component
C  uTrans,vTrans   :: (input) volume transports at U,V points
C  iceFld          :: (input) field of iceberg distribution
C  r_hFld          :: (input) reciprocal of ice-thickness (only used for "intensive"
C                     type of iceberg field)
C  bi,bj           ::  tile indices
C  myTime          :: current time
C  myIter          :: iteration number
C  myThid          :: my Thread Id number
      INTEGER tracerIdentity
      INTEGER advectionScheme
      _RL uFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL vFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL uTrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL vTrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL iceFld(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL r_hFld(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER bi,bj
      _RL myTime
      INTEGER myIter
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  gFld          :: tendency array
C  afx           :: horizontal advective flux, x direction
C  afy           :: horizontal advective flux, y direction
      _RL gFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL afx   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL afy   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)

#ifdef ALLOW_GENERIC_ADVDIFF
C !LOCAL VARIABLES: ====================================================
C  maskLocW      :: 2-D array for mask at West points
C  maskLocS      :: 2-D array for mask at South points
C  iMin,iMax,    :: loop range for called routines
C  jMin,jMax     :: loop range for called routines
C [iMin,iMax]Upd :: loop range to update iceberg field
C [jMin,jMax]Upd :: loop range to update iceberg field
C  i,j,k         :: loop indices
C  af            :: 2-D array for horizontal advective flux
C  localTij      :: 2-D array, temporary local copy of iceberg fld
C  calc_fluxes_X :: logical to indicate to calculate fluxes in X dir
C  calc_fluxes_Y :: logical to indicate to calculate fluxes in Y dir
C  interiorOnly  :: only update the interior of myTile, but not the edges
C  overlapOnly   :: only update the edges of myTile, but not the interior
C  nipass        :: number of passes in multi-dimensional method
C  ipass         :: number of the current pass being made
C  myTile        :: variables used to determine which cube face
C  nCFace        :: owns a tile for cube grid runs using
C                :: multi-dim advection.
C [N,S,E,W]_edge :: true if N,S,E,W edge of myTile is an Edge of the cube
C     msgBuf     :: Informational/error message buffer
      _RS maskLocW(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS maskLocS(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      INTEGER iMin,iMax,jMin,jMax
      INTEGER iMinUpd,iMaxUpd,jMinUpd,jMaxUpd
      INTEGER i,j,k
      _RL af      (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL localTij(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      LOGICAL calc_fluxes_X, calc_fluxes_Y, withSigns
      LOGICAL interiorOnly, overlapOnly
      INTEGER nipass,ipass
      INTEGER nCFace
      LOGICAL N_edge, S_edge, E_edge, W_edge
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#ifdef ALLOW_EXCH2
      INTEGER myTile
#endif
#ifdef ALLOW_DIAGNOSTICS
      CHARACTER*8 diagName
#endif

      INTEGER ioUnit
      _RL tmpFac

CEOP



C myBxLo, myBxHi are bounds of loop bi
C myByLo, myByHi are bounds of loop bj
C nTx number of threads per process in x
C nTy number of threads per process in y
C maxpass default =2
C ikey_dynamics - key for reading and writing data in routine dynamics
C igadkey ????
C these variables are from tamc.h which has to be customized to use it in the simulation and be put in the code dir


#ifdef ALLOW_AUTODIFF_TAMC
          act0 = tracerIdentity - 1
          max0 = maxpass
          act1 = bi - myBxLo(myThid)
          max1 = myBxHi(myThid) - myBxLo(myThid) + 1
          act2 = bj - myByLo(myThid)
          max2 = myByHi(myThid) - myByLo(myThid) + 1
          act3 = myThid - 1
          max3 = nTx*nTy
          act4 = ikey_dynamics - 1
          igadkey = (act0 + 1)
     &                      + act1*max0
     &                      + act2*max0*max1
     &                      + act3*max0*max1*max2
     &                      + act4*max0*max1*max2*max3

C Wenn die tracerNummer groeßer ist als die maximale Anzahl an Dingen, die uebergeben werden koennen --> FEHLER
          if (tracerIdentity.GT.maxpass) then
           WRITE(msgBuf,'(A,2I3)')
     &      'ICEBERG_ADVECTION: maxpass < tracerIdentity ',
     &      maxpass, tracerIdentity
           CALL PRINT_ERROR( msgBuf, myThid )
           STOP 'ABNORMAL END: S/R ICEBERG_ADVECTION'
          endif
#endif /* ALLOW_AUTODIFF_TAMC */


C--   Set up work arrays with valid (i.e. not NaN) values
C     These inital values do not alter the numerical results. They
C     just ensure that all memory references are to valid floating
C     point numbers. This prevents spurious hardware signals due to
C     uninitialised but inert locations.
#ifdef ALLOW_AUTODIFF_TAMC
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
C local copy of iceberg field
        localTij(i,j) = 0. _d 0
       ENDDO
      ENDDO
#endif

C--   Set tile-specific parameters for horizontal fluxes
C uceCubedSphereExchange has to be set to TRUE if any
C type of grid composed of interconnected individual faces is used
C including cube sphere and lat-lon cap grid
      IF (useCubedSphereExchange) THEN
C number of passes in multi-dimensional method
       nipass=3
#ifdef ALLOW_AUTODIFF_TAMC
       IF ( nipass.GT.maxcube ) THEN
           WRITE(msgBuf,'(A)')
     &      'ICEBERG_ADVECTION: maxcube needs to be =3; check tamc.h '
           CALL PRINT_ERROR( msgBuf, myThid )
           STOP 'ABNORMAL END: S/R ICEBERG_ADVECTION'
       ENDIF
#endif
#ifdef ALLOW_EXCH2
C Check if actual tile is at the edge of a cube sphere subdomain
C 1 if yes and 0 if not
C if yes the case requires special exchange and numerical
C handling for singularities at the eight corners of the cube
C
C W2_myTileList :: List of tiles owned by this process
C exch2_myFace :: contains number of subdomains of each tile
C exhc2_is_edge :: set to 1 if tile lies on the edge
       myTile = W2_myTileList(bi,bj)
       nCFace = exch2_myFace(myTile)
       N_edge = exch2_isNedge(myTile).EQ.1
       S_edge = exch2_isSedge(myTile).EQ.1
       E_edge = exch2_isEedge(myTile).EQ.1
       W_edge = exch2_isWedge(myTile).EQ.1
#else
C if exch2 is not used but cubed sphere scheme
       nCFace = bi
       N_edge = .TRUE.
       S_edge = .TRUE.
       E_edge = .TRUE.
       W_edge = .TRUE.
#endif
      ELSE
C if cube sphere or lat-lon cap not used
       nipass=2
       nCFace = bi
       N_edge = .FALSE.
       S_edge = .FALSE.
       E_edge = .FALSE.
       W_edge = .FALSE.
      ENDIF

C set loop range variables for called routines
      iMin = 1-OLx
      iMax = sNx+OLx
      jMin = 1-OLy
      jMax = sNy+OLy

      k = 1
C--   Start of k loop for horizontal fluxes
#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE iceFld =
CADJ &     comlev1_bibj_k_gadice, key=igadkey, byte=isbyte
#endif /* ALLOW_AUTODIFF_TAMC */

C     Content of CALC_COMMON_FACTORS, adapted for 2D fields
C--   Get temporary terms used by tendency routines



C--   Make local copy of sea-ice field and mask West & South
C from GAD_ADVECTION adapted for SIMask
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
         localTij(i,j)=iceFld(i,j)


#ifdef ALLOW_OBCS
C SImaskU is maskW at the Surface
C SIMaskU(i,j,bi,bj)  :: maskW(i,j,kSurface,bi,bj) with kSurface as Nr or 1
C SIMaskV(i,j,bi,bj)  :: maskS(i,j,kSurface,bi,bj)
C TODO define SImask
C maskInW :: West face 2-D Interior mask (i.e., zero on and beyond OB)
C         maskLocW(i,j) = SIMaskU(i,j,bi,bj)*maskInW(i,j,bi,bj)
         maskLocW(i,j) = maskW(i,j,1,bi,bj)*maskInW(i,j,bi,bj)
C         maskLocS(i,j) = SIMaskV(i,j,bi,bj)*maskInS(i,j,bi,bj)
         maskLocS(i,j) = maskS(i,j,1,bi,bj)*maskInS(i,j,bi,bj)

#else /* ALLOW_OBCS */
c         maskLocW(i,j) = SIMaskU(i,j,bi,bj)
         maskLocW(i,j) = maskW(i,j,1,bi,bj)
c         maskLocS(i,j) = SIMaskV(i,j,bi,bj)
         maskLocS(i,j) = maskS(i,j,1,bi,bj)
#endif /* ALLOW_OBCS */

       ENDDO
      ENDDO



#ifdef ALLOW_AUTODIFF_TAMC
C-     Initialise Advective flux in X & Y
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         afx(i,j) = 0.
         afy(i,j) = 0.
        ENDDO
       ENDDO
#endif



C--   Multiple passes for different directions on different tiles
C--   For cube need one pass for each of red, green and blue axes.
C from GAD_ADVECTION
      DO ipass=1,nipass
#ifdef ALLOW_AUTODIFF_TAMC
C passkey is used for storing data
       passkey = ipass + (igadkey-1)*maxpass
       IF (nipass .GT. maxpass) THEN
           WRITE(msgBuf,'(A,2I3)')
     &      'ICEBERG_ADVECTION: nipass > max[ass. check tamc.h ',
     &      nipass, maxpass
           CALL PRINT_ERROR( msgBuf, myThid )
           STOP 'ABNORMAL END: S/R ICEBERG_ADVECTION'
       ENDIF
#endif /* ALLOW_AUTODIFF_TAMC */

       interiorOnly = .FALSE.
       overlapOnly  = .FALSE.
       IF (useCubedSphereExchange) THEN
C--   CubedSphere : pass 3 times, with partial update of local si field
C calc directions and interior or overlap are defined (see manual fig 2.11)
        IF (ipass.EQ.1) THEN
         overlapOnly  = MOD(nCFace,3).EQ.0
         interiorOnly = MOD(nCFace,3).NE.0
         calc_fluxes_X = nCFace.EQ.6 .OR. nCFace.EQ.1 .OR. nCFace.EQ.2
         calc_fluxes_Y = nCFace.EQ.3 .OR. nCFace.EQ.4 .OR. nCFace.EQ.5
        ELSEIF (ipass.EQ.2) THEN
         overlapOnly  = MOD(nCFace,3).EQ.2
         calc_fluxes_X = nCFace.EQ.2 .OR. nCFace.EQ.3 .OR. nCFace.EQ.4
         calc_fluxes_Y = nCFace.EQ.5 .OR. nCFace.EQ.6 .OR. nCFace.EQ.1
        ELSE
         calc_fluxes_X = nCFace.EQ.5 .OR. nCFace.EQ.6
         calc_fluxes_Y = nCFace.EQ.2 .OR. nCFace.EQ.3
        ENDIF

       ELSE
C--   not CubedSphere
        calc_fluxes_X = MOD(ipass,2).EQ.1
        calc_fluxes_Y = .NOT.calc_fluxes_X
       ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   X direction


#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE localTij(:,:)  =
CADJ &     comlev1_bibj_k_gadice_pass, key=passkey, byte=isbyte
# ifndef DISABLE_MULTIDIM_ADVECTION
CADJ STORE af(:,:)  =
CADJ &     comlev1_bibj_k_gadice_pass, key=passkey, byte=isbyte
# endif
#endif /* ALLOW_AUTODIFF_TAMC */
C


C also from GAD_ADVECTION
       IF (calc_fluxes_X) THEN

C-     Do not compute fluxes if
C       a) needed in overlap only
C   and b) the overlap of myTile are not cube-face Edges
        IF ( .NOT.overlapOnly .OR. N_edge .OR. S_edge ) THEN

C-     Advective flux in X
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           af(i,j) = 0.
          ENDDO
         ENDDO

cph-exch2#ifndef ALLOW_AUTODIFF_TAMC
C-     Internal exchange for calculations in X
C if cubesphere and overlap only or first pass then fill the corners
         IF ( useCubedSphereExchange .AND.
     &      ( overlapOnly .OR. ipass.EQ.1 ) ) THEN
          CALL FILL_CS_CORNER_TR_RL( 1, .FALSE.,
     &                               localTij, bi,bj, myThid )
         ENDIF
cph-exch2#endif

#ifdef ALLOW_AUTODIFF_TAMC
# ifndef DISABLE_MULTIDIM_ADVECTION
CADJ STORE localTij(:,:)  =
CADJ &     comlev1_bibj_k_gadice_pass, key=passkey, byte=isbyte
# endif
#endif /* ALLOW_AUTODIFF_TAMC */

C Do calculation for advection flux in X with build in schemes
C defined as code in advectionScheme

C Upwind 1st order (1) or 2nd order DST (Lax-Wendroff) (2)
         IF ( advectionScheme.EQ.ENUM_UPWIND_1RST
     &        .OR. advectionScheme.EQ.ENUM_DST2 ) THEN
          CALL GAD_DST2U1_ADV_X( bi,bj,k, advectionScheme, .TRUE.,
     I         ICEBERG_deltaTtherm, uTrans, uFld, localTij,
     O         af, myThid )

C flux limiter
         ELSEIF (advectionScheme.EQ.ENUM_FLUX_LIMIT) THEN
          CALL GAD_FLUXLIMIT_ADV_X( bi,bj,k, .TRUE.,
     I         ICEBERG_deltaTtherm, uTrans, uFld, maskLocW, localTij,
     O         af, myThid )

C 3rd order DST (30)
         ELSEIF (advectionScheme.EQ.ENUM_DST3 ) THEN
          CALL GAD_DST3_ADV_X(      bi,bj,k, .TRUE.,
     I         ICEBERG_deltaTtherm, uTrans, uFld, maskLocW, localTij,
     O         af, myThid )

C 3rd order DST flux limiter (33)
         ELSEIF (advectionScheme.EQ.ENUM_DST3_FLUX_LIMIT ) THEN
          CALL GAD_DST3FL_ADV_X(    bi,bj,k, .TRUE.,
     I         ICEBERG_deltaTtherm, uTrans, uFld, maskLocW, localTij,
     O         af, myThid )
#ifndef ALLOW_AUTODIFF_TAMC

c 7th order one-step method (7)
         ELSEIF (advectionScheme.EQ.ENUM_OS7MP ) THEN
          CALL GAD_OS7MP_ADV_X(     bi,bj,k, .TRUE.,
     I         ICEBERG_deltaTtherm, uTrans, uFld, maskLocW, localTij,
     O         af, myThid )

C piecewiese parabolic "null" limiter (40)
         ELSEIF (advectionScheme.EQ.ENUM_PPM_NULL_LIMIT  .OR.
     &           advectionScheme.EQ.ENUM_PPM_MONO_LIMIT  .OR.
     &           advectionScheme.EQ.ENUM_PPM_WENO_LIMIT) THEN
          CALL GAD_PPM_ADV_X( advectionScheme, bi, bj, k , .TRUE.,
     I         ICEBERG_deltaTtherm, uFld, uTrans, localTij,
     O         af, myThid )

C piecewise quartic "null" limiter (40)
         ELSEIF (advectionScheme.EQ.ENUM_PQM_NULL_LIMIT  .OR.
     &           advectionScheme.EQ.ENUM_PQM_MONO_LIMIT  .OR.
     &           advectionScheme.EQ.ENUM_PQM_WENO_LIMIT) THEN
          CALL GAD_PQM_ADV_X( advectionScheme, bi, bj, k , .TRUE.,
     I         ICEBERG_deltaTtherm, uFld, uTrans, localTij,
     O         af, myThid )
#endif
         ELSE
           WRITE(msgBuf,'(A,I3,A)')
     &      'ICEBERG_ADVECTION: adv. scheme ', advectionScheme,
     &      ' incompatibale with multi-dim. adv.'
           CALL PRINT_ERROR( msgBuf, myThid )
           STOP 'ABNORMAL END: S/R ICEBERG_ADVECTION'
         ENDIF

C--   Advective flux in X : done
        ENDIF

cph-exch2#ifndef ALLOW_AUTODIFF_TAMC
C--   Internal exchange for next calculations in Y
        IF ( overlapOnly .AND. ipass.EQ.1 ) THEN
         CALL FILL_CS_CORNER_TR_RL( 2, .FALSE.,
     &                              localTij, bi,bj, myThid )
        ENDIF
cph-exch2#endif

C-     Update the local si field where needed:

C     update in overlap-Only
        IF ( overlapOnly ) THEN
         iMinUpd = 1-OLx+1
         iMaxUpd = sNx+OLx-1
C--   notes: these 2 lines below have no real effect (because recip_hFac=0
C            in corner region) but safer to keep them.
         IF ( W_edge ) iMinUpd = 1
         IF ( E_edge ) iMaxUpd = sNx

C extensiveFld :: indicates to advect an "extensive" type of ice field
C recip_rA is the reciprocal of the grid
C masinC :: cell center 2D interior mask
         IF ( S_edge .AND. extensiveFld ) THEN
          DO j=1-OLy,0
           DO i=iMinUpd,iMaxUpd
            localTij(i,j)=localTij(i,j)
     &         -ICEBERG_deltaTtherm*maskInC(i,j,bi,bj)
     &           *recip_rA(i,j,bi,bj)
     &           *(  af(i+1,j)-af(i,j)
     &            )
           ENDDO
          ENDDO
         ELSEIF ( S_edge ) THEN


          DO j=1-OLy,0
           DO i=iMinUpd,iMaxUpd
            localTij(i,j)=localTij(i,j)
     &         -ICEBERG_deltaTtherm*maskInC(i,j,bi,bj)
     &           *recip_rA(i,j,bi,bj)*r_hFld(i,j)
     &           *( (af(i+1,j)-af(i,j))
     &             -(uTrans(i+1,j)-uTrans(i,j))*iceFld(i,j)
     &            )
           ENDDO
          ENDDO
         ENDIF


         IF ( N_edge .AND. extensiveFld ) THEN


          DO j=sNy+1,sNy+OLy
           DO i=iMinUpd,iMaxUpd
            localTij(i,j)=localTij(i,j)
     &         -ICEBERG_deltaTtherm*maskInC(i,j,bi,bj)
     &           *recip_rA(i,j,bi,bj)
     &           *(  af(i+1,j)-af(i,j)
     &            )
           ENDDO
          ENDDO
         ELSEIF ( N_edge ) THEN
          DO j=sNy+1,sNy+OLy
           DO i=iMinUpd,iMaxUpd
            localTij(i,j)=localTij(i,j)
     &         -ICEBERG_deltaTtherm*maskInC(i,j,bi,bj)
     &           *recip_rA(i,j,bi,bj)*r_hFld(i,j)
     &           *( (af(i+1,j)-af(i,j))
     &             -(uTrans(i+1,j)-uTrans(i,j))*iceFld(i,j)
     &            )
           ENDDO
          ENDDO
         ENDIF
C--   keep advective flux (for diagnostics) - not in GAD_advection
C save advective flux from local parameter af to output parameter afx
         IF ( S_edge ) THEN
          DO j=1-OLy,0
           DO i=1-OLx+1,sNx+OLx
            afx(i,j) = af(i,j)
           ENDDO
          ENDDO
         ENDIF
         IF ( N_edge ) THEN
          DO j=sNy+1,sNy+OLy
           DO i=1-OLx+1,sNx+OLx
            afx(i,j) = af(i,j)
           ENDDO
          ENDDO
         ENDIF

        ELSE
C     do not only update the overlap (from GAD_advection)
         jMinUpd = 1-OLy
         jMaxUpd = sNy+OLy
         IF ( interiorOnly .AND. S_edge ) jMinUpd = 1
         IF ( interiorOnly .AND. N_edge ) jMaxUpd = sNy
         IF ( extensiveFld ) THEN
          DO j=jMinUpd,jMaxUpd
           DO i=1-OLx+1,sNx+OLx-1
            localTij(i,j)=localTij(i,j)
     &         -ICEBERG_deltaTtherm*maskInC(i,j,bi,bj)
     &           *recip_rA(i,j,bi,bj)
     &           *(  af(i+1,j)-af(i,j)
     &            )
           ENDDO
          ENDDO

         ELSE
          DO j=jMinUpd,jMaxUpd
           DO i=1-OLx+1,sNx+OLx-1
            localTij(i,j)=localTij(i,j)
     &         -ICEBERG_deltaTtherm*maskInC(i,j,bi,bj)
     &           *recip_rA(i,j,bi,bj)*r_hFld(i,j)
     &           *( (af(i+1,j)-af(i,j))
     &             -(uTrans(i+1,j)-uTrans(i,j))*iceFld(i,j)
     &            )
           ENDDO
          ENDDO
         ENDIF

C--   keep advective flux (for diagnostics)
         DO j=jMinUpd,jMaxUpd
           DO i=1-OLx+1,sNx+OLx
            afx(i,j) = af(i,j)
           ENDDO
         ENDDO

C-     end if/else update overlap-Only
        ENDIF

C--   End of X direction
       ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Y direction

#ifdef ALLOW_AUTODIFF_TAMC
# ifndef DISABLE_MULTIDIM_ADVECTION
CADJ STORE localTij(:,:)  =
CADJ &     comlev1_bibj_k_gadice_pass, key=passkey, byte=isbyte
CADJ STORE af(:,:)  =
CADJ &     comlev1_bibj_k_gadice_pass, key=passkey, byte=isbyte
# endif
#endif /* ALLOW_AUTODIFF_TAMC */

C from GAD_Advection
       IF (calc_fluxes_Y) THEN

C-     Do not compute fluxes if
C       a) needed in overlap only
C   and b) the overlap of myTile are not cube-face edges
        IF ( .NOT.overlapOnly .OR. E_edge .OR. W_edge ) THEN

C-     Advective flux in Y
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           af(i,j) = 0.
          ENDDO
         ENDDO

cph-exch2#ifndef ALLOW_AUTODIFF_TAMC
C-     Internal exchange for calculations in Y
C if cube sphere and overlap only or it is the first tile
         IF ( useCubedSphereExchange .AND.
     &      ( overlapOnly .OR. ipass.EQ.1 ) ) THEN
          CALL FILL_CS_CORNER_TR_RL( 2, .FALSE.,
     &                               localTij, bi,bj, myThid )
         ENDIF
cph-exch2#endif

#ifdef ALLOW_AUTODIFF_TAMC
#ifndef DISABLE_MULTIDIM_ADVECTION
CADJ STORE localTij(:,:)  =
CADJ &     comlev1_bibj_k_gadice_pass, key=passkey, byte=isbyte
#endif
#endif /* ALLOW_AUTODIFF_TAMC */

         IF ( advectionScheme.EQ.ENUM_UPWIND_1RST
     &        .OR. advectionScheme.EQ.ENUM_DST2 ) THEN
          CALL GAD_DST2U1_ADV_Y( bi,bj,k, advectionScheme, .TRUE.,
     I         ICEBERG_deltaTtherm, vTrans, vFld, localTij,
     O         af, myThid )

         ELSEIF (advectionScheme.EQ.ENUM_FLUX_LIMIT) THEN
          CALL GAD_FLUXLIMIT_ADV_Y( bi,bj,k, .TRUE.,
     I         ICEBERG_deltaTtherm, vTrans, vFld, maskLocS, localTij,
     O         af, myThid )
         ELSEIF (advectionScheme.EQ.ENUM_DST3 ) THEN
          CALL GAD_DST3_ADV_Y(      bi,bj,k, .TRUE.,
     I         ICEBERG_deltaTtherm, vTrans, vFld, maskLocS, localTij,
     O         af, myThid )
         ELSEIF (advectionScheme.EQ.ENUM_DST3_FLUX_LIMIT ) THEN
          CALL GAD_DST3FL_ADV_Y(    bi,bj,k, .TRUE.,
     I         ICEBERG_deltaTtherm, vTrans, vFld, maskLocS, localTij,
     O         af, myThid )
#ifndef ALLOW_AUTODIFF_TAMC
         ELSEIF (advectionScheme.EQ.ENUM_OS7MP ) THEN
          CALL GAD_OS7MP_ADV_Y(     bi,bj,k, .TRUE.,
     I         ICEBERG_deltaTtherm, vTrans, vFld, maskLocS, localTij,
     O         af, myThid )
         ELSEIF (advectionScheme.EQ.ENUM_PPM_NULL_LIMIT  .OR.
     &           advectionScheme.EQ.ENUM_PPM_MONO_LIMIT  .OR.
     &           advectionScheme.EQ.ENUM_PPM_WENO_LIMIT) THEN
          CALL GAD_PPM_ADV_Y( advectionScheme, bi, bj, k , .TRUE.,
     I         ICEBERG_deltaTtherm, vFld, vTrans, localTij,
     O         af, myThid )
         ELSEIF (advectionScheme.EQ.ENUM_PQM_NULL_LIMIT  .OR.
     &           advectionScheme.EQ.ENUM_PQM_MONO_LIMIT  .OR.
     &           advectionScheme.EQ.ENUM_PQM_WENO_LIMIT) THEN
          CALL GAD_PQM_ADV_Y( advectionScheme, bi, bj, k , .TRUE.,
     I         ICEBERG_deltaTtherm, vFld, vTrans, localTij,
     O         af, myThid )
#endif
         ELSE
           WRITE(msgBuf,'(A,I3,A)')
     &      'ICEBERG_ADVECTION: adv. scheme ', advectionScheme,
     &      ' incompatibale with multi-dim. adv.'
           CALL PRINT_ERROR( msgBuf, myThid )
           STOP 'ABNORMAL END: S/R ICEBERG_ADVECTION'
         ENDIF

C-     Advective flux in Y : done
        ENDIF

cph-exch2#ifndef ALLOW_AUTODIFF_TAMC
C-     Internal exchange for next calculations in X
        IF ( overlapOnly .AND. ipass.EQ.1 ) THEN
         CALL FILL_CS_CORNER_TR_RL( 1, .FALSE.,
     &                              localTij, bi,bj, myThid )
        ENDIF
cph-exch2#endif

C-     Update the local si field where needed:

C      update in overlap-Only
        IF ( overlapOnly ) THEN
         jMinUpd = 1-OLy+1
         jMaxUpd = sNy+OLy-1
C- notes: these 2 lines below have no real effect (because recip_hFac=0
C         in corner region) but safer to keep them.
         IF ( S_edge ) jMinUpd = 1
         IF ( N_edge ) jMaxUpd = sNy

         IF ( W_edge .AND. extensiveFld ) THEN
          DO j=jMinUpd,jMaxUpd
           DO i=1-OLx,0
            localTij(i,j)=localTij(i,j)
     &         -ICEBERG_deltaTtherm*maskInC(i,j,bi,bj)
     &           *recip_rA(i,j,bi,bj)
     &           *(  af(i,j+1)-af(i,j)
     &            )
           ENDDO
          ENDDO
         ELSEIF ( W_edge ) THEN
          DO j=jMinUpd,jMaxUpd
           DO i=1-OLx,0

            localTij(i,j)=localTij(i,j)
     &         -ICEBERG_deltaTtherm*maskInC(i,j,bi,bj)
     &           *recip_rA(i,j,bi,bj)*r_hFld(i,j)
     &           *( (af(i,j+1)-af(i,j))
     &             -(vTrans(i,j+1)-vTrans(i,j))*iceFld(i,j)
     &            )
           ENDDO
          ENDDO
         ENDIF
         IF ( E_edge .AND. extensiveFld ) THEN
          DO j=jMinUpd,jMaxUpd
           DO i=sNx+1,sNx+OLx
            localTij(i,j)=localTij(i,j)
     &         -ICEBERG_deltaTtherm*maskInC(i,j,bi,bj)
     &           *recip_rA(i,j,bi,bj)
     &           *(  af(i,j+1)-af(i,j)
     &            )
           ENDDO
          ENDDO
         ELSEIF ( E_edge ) THEN
          DO j=jMinUpd,jMaxUpd
           DO i=sNx+1,sNx+OLx
            localTij(i,j)=localTij(i,j)
     &         -ICEBERG_deltaTtherm*maskInC(i,j,bi,bj)
     &           *recip_rA(i,j,bi,bj)*r_hFld(i,j)
     &           *( (af(i,j+1)-af(i,j))
     &             -(vTrans(i,j+1)-vTrans(i,j))*iceFld(i,j)
     &            )
           ENDDO
          ENDDO
         ENDIF


C--   keep advective flux (for diagnostics)
         IF ( W_edge ) THEN
          DO j=1-OLy+1,sNy+OLy
           DO i=1-OLx,0
            afy(i,j) = af(i,j)
           ENDDO
          ENDDO
         ENDIF
         IF ( E_edge ) THEN
          DO j=1-OLy+1,sNy+OLy
           DO i=sNx+1,sNx+OLx
            afy(i,j) = af(i,j)
           ENDDO
          ENDDO
         ENDIF

        ELSE
C     do not only update the overlap
         iMinUpd = 1-OLx
         iMaxUpd = sNx+OLx
         IF ( interiorOnly .AND. W_edge ) iMinUpd = 1
         IF ( interiorOnly .AND. E_edge ) iMaxUpd = sNx
         IF ( extensiveFld ) THEN
          DO j=1-OLy+1,sNy+OLy-1
           DO i=iMinUpd,iMaxUpd

            localTij(i,j)=localTij(i,j)
     &         -ICEBERG_deltaTtherm*maskInC(i,j,bi,bj)
     &           *recip_rA(i,j,bi,bj)
     &           *(  af(i,j+1)-af(i,j)
     &            )
           ENDDO
          ENDDO
         ELSE

          DO j=1-OLy+1,sNy+OLy-1
           DO i=iMinUpd,iMaxUpd
            localTij(i,j)=localTij(i,j)
     &         -ICEBERG_deltaTtherm*maskInC(i,j,bi,bj)
     &           *recip_rA(i,j,bi,bj)*r_hFld(i,j)
     &           *( (af(i,j+1)-af(i,j))
     &             -(vTrans(i,j+1)-vTrans(i,j))*iceFld(i,j)
     &            )
           ENDDO
          ENDDO
         ENDIF
         
C--   keep advective flux (for diagnostics)
         DO j=1-OLy+1,sNy+OLy
           DO i=iMinUpd,iMaxUpd
            afy(i,j) = af(i,j)
           ENDDO
         ENDDO

C      end if/else update overlap-Only
        ENDIF

C--   End of Y direction
       ENDIF

C--   End of ipass loop
      ENDDO



C-    explicit advection is done ; store tendency in gFld:
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
        gFld(i,j)=(localTij(i,j)-iceFld(i,j))/ICEBERG_deltaTtherm
       ENDDO
      ENDDO


c TODO adapt diagnostics_fill for iceberg_distr
c#ifdef ALLOW_DIAGNOSTICS
c        IF ( useDiagnostics ) THEN
c         diagName = ICEBERG_ClLbl(tracerIdentity)
c         CALL DIAGNOSTICS_FILL(afx,diagName, k,1, 2,bi,bj, myThid)
c         diagName = 'ADVy'//diagSufx
c         CALL DIAGNOSTICS_FILL(afy,diagName, k,1, 2,bi,bj, myThid)
c        ENDIF
c#endif

#endif /* ALLOW_GENERIC_ADVDIFF */
      RETURN
      END
