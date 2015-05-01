        !COMPILER-GENERATED INTERFACE MODULE: Thu May 29 14:38:37 2014
        MODULE ELQUAD4E__genmod
          INTERFACE 
            SUBROUTINE ELQUAD4E(NEL,NELP,KE,EL_DDISP,FINTE,INC)
              USE MOD_VARIABLES
              INTEGER(KIND=4) :: NEL
              INTEGER(KIND=4) :: NELP
              REAL(KIND=8), INTENT(OUT) :: KE((P+1)*(Q+1)*NDS,(P+1)*(Q+1&
     &)*NDS)
              REAL(KIND=8), INTENT(IN) :: EL_DDISP((P+1)*(Q+1)*NDS,1)
              REAL(KIND=8), INTENT(OUT) :: FINTE((P+1)*(Q+1)*NDS,1)
              INTEGER(KIND=4) :: INC
            END SUBROUTINE ELQUAD4E
          END INTERFACE 
        END MODULE ELQUAD4E__genmod
