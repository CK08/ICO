        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 14 08:53:47 2014
        MODULE SYMMLQ__genmod
          INTERFACE 
            SUBROUTINE SYMMLQ(N,B,R1,R2,V,W,X,Y,APROD,MSOLVE,CHECKA,    &
     &GOODB,PRECON,SHIFT,NOUT,ITNLIM,RTOL,ISTOP,ITN,ANORM,ACOND,RNORM,  &
     &YNORM)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: B(N)
              REAL(KIND=8) :: R1(N)
              REAL(KIND=8) :: R2(N)
              REAL(KIND=8) :: V(N)
              REAL(KIND=8) :: W(N)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: Y(N)
              EXTERNAL APROD
              EXTERNAL MSOLVE
              LOGICAL(KIND=4) :: CHECKA
              LOGICAL(KIND=4) :: GOODB
              LOGICAL(KIND=4) :: PRECON
              REAL(KIND=8) :: SHIFT
              INTEGER(KIND=4) :: NOUT
              INTEGER(KIND=4) :: ITNLIM
              REAL(KIND=8) :: RTOL
              INTEGER(KIND=4) :: ISTOP
              INTEGER(KIND=4) :: ITN
              REAL(KIND=8) :: ANORM
              REAL(KIND=8) :: ACOND
              REAL(KIND=8) :: RNORM
              REAL(KIND=8) :: YNORM
            END SUBROUTINE SYMMLQ
          END INTERFACE 
        END MODULE SYMMLQ__genmod
