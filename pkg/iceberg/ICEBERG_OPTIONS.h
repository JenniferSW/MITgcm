C CPP options file for ICEBERG
C Use this file for selecting options within package "Iceberg"

#ifndef ICEBERG_OPTIONS_H
#define ICEBERG_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_ICEBERG
C Place CPP define/undef flag here

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
