        !COMPILER-GENERATED INTERFACE MODULE: Thu May 29 14:18:24 2014
        MODULE FINDSPAN__genmod
          INTERFACE 
            SUBROUTINE FINDSPAN(N,P,U,U_KNOT,S)
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(IN) :: P
              REAL(KIND=8), INTENT(IN) :: U
              REAL(KIND=8), INTENT(IN) :: U_KNOT(N)
              INTEGER(KIND=4), INTENT(OUT) :: S
            END SUBROUTINE FINDSPAN
          END INTERFACE 
        END MODULE FINDSPAN__genmod
