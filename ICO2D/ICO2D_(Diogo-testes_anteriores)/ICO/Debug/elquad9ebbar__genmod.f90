        !COMPILER-GENERATED INTERFACE MODULE: Thu May 29 14:18:25 2014
        MODULE ELQUAD9EBBAR__genmod
          INTERFACE 
            SUBROUTINE ELQUAD9EBBAR(NEL,KE,KBAR,EL_DDISP,FINTE)
              USE MOD_VARIABLES
              INTEGER(KIND=4) :: NEL
              REAL(KIND=8), INTENT(OUT) :: KE((P+1)*(Q+1)*NDS,(P+1)*(Q+1&
     &)*NDS)
              REAL(KIND=8), INTENT(OUT) :: KBAR(P*Q*NDS,P*Q*NDS)
              REAL(KIND=8), INTENT(IN) :: EL_DDISP((P+1)*(Q+1)*NDS,1)
              REAL(KIND=8), INTENT(OUT) :: FINTE((P+1)*(Q+1)*NDS,1)
            END SUBROUTINE ELQUAD9EBBAR
          END INTERFACE 
        END MODULE ELQUAD9EBBAR__genmod
