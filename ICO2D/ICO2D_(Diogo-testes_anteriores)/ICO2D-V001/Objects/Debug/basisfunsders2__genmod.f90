        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 03 16:45:32 2014
        MODULE BASISFUNSDERS2__genmod
          INTERFACE 
            SUBROUTINE BASISFUNSDERS2(NB,S,U,P,N,U_KNOT,DERS)
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(IN) :: P
              INTEGER(KIND=4), INTENT(IN) :: NB
              INTEGER(KIND=4), INTENT(IN) :: S
              REAL(KIND=8), INTENT(IN) :: U
              REAL(KIND=8), INTENT(IN) :: U_KNOT(NB)
              REAL(KIND=8), INTENT(OUT) :: DERS(N+1,NB)
            END SUBROUTINE BASISFUNSDERS2
          END INTERFACE 
        END MODULE BASISFUNSDERS2__genmod
