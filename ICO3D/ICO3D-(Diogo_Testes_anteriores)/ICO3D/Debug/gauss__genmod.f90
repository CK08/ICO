        !COMPILER-GENERATED INTERFACE MODULE: Tue Feb 18 09:01:08 2014
        MODULE GAUSS__genmod
          INTERFACE 
            SUBROUTINE GAUSS(DP,A,VS,V)
              INTEGER(KIND=4), INTENT(IN) :: DP
              REAL(KIND=8), INTENT(INOUT) :: A(DP,DP)
              REAL(KIND=8), INTENT(INOUT) :: VS(DP)
              REAL(KIND=8), INTENT(INOUT) :: V(DP)
            END SUBROUTINE GAUSS
          END INTERFACE 
        END MODULE GAUSS__genmod
