#ifdef ALLOW_ICEBERG

CBOP
C     !ROUTINE: ICEBERG_SIZE.h
C     !INTERFACE:
C #include ICEBERG_SIZE.h


C     !DESPCRIPTION
C     Specifies number of size classes as it has to be defined before compiling
c     in order to initializes the fields for iceberg distribution and some other
C     variables relying on the size class number. (default is 1)


C    Number of size classes
      INTEGER ICEBERG_numCl
      PARAMETER(ICEBERG_numCl = 1)
CEOP
#endif
