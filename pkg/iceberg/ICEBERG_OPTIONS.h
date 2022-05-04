C CPP options file for ICEBERG
C Use this file for selecting options within package "Iceberg"

#ifndef ICEBERG_OPTIONS_H
#define ICEBERG_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_ICEBERG
C Place CPP define/undef flag here

C CPP flags for use of dynamics, melting, calving and freshwater input
#define ALLOW_ICEBERG_DYNAMICS
#define ALLOW_ICEBERG_MELTING
#define ALLOW_ICEBERG_CALVING
#define ALLOW_FRESHWATER_INPUT
#define ALLOW_GROUNDING

C CPP flags for the use of forces in dynamics
#define ALLOW_ICEBERG_CORIOLIS
#define ALLOW_ICEBERG_OCEAN_DRAG
#define ALLOW_ICEBERG_ATMOSPHERE_DRAG
#define ALLOW_ICEBERG_SEAICE_DRAG
#define ALLOW_ICEBERG_WAVE_RADIATION
#define ALLOW_ICEBERG_PRESSURE_GRADIENT

C to reduce memory storage, disable unused array with those CPP flags :
#define ICEBERG_3D_STATE
#define ICEBERG_2D_STATE
#define ICEBERG_TENDENCY


#undef MYPA_SPECIAL_COMPILE_OPTION1

#define MYPA_SPECIAL_COMPILE_OPTION2

#endif /* ALLOW_ICEBERG */
#endif /* ICEBERG_OPTIONS_H */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
