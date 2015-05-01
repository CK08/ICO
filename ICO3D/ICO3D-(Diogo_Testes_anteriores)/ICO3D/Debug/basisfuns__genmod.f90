        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 14 08:53:48 2014
        MODULE BASISFUNS__genmod
          INTERFACE 
            SUBROUTINE BASISFUNS(NB,S,U,P,U_KNOTS,BF)
              INTEGER(KIND=4), INTENT(IN) :: P
              INTEGER(KIND=4), INTENT(IN) :: NB
              INTEGER(KIND=4), INTENT(IN) :: S
              REAL(KIND=8), INTENT(IN) :: U
              REAL(KIND=8), INTENT(IN) :: U_KNOTS(NB)
              REAL(KIND=8), INTENT(OUT) :: BF(P+1,1)
            END SUBROUTINE BASISFUNS
          END INTERFACE 
        END MODULE BASISFUNS__genmod
