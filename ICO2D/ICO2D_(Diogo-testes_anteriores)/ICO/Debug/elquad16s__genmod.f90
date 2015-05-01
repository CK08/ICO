        !COMPILER-GENERATED INTERFACE MODULE: Wed May 28 14:32:19 2014
        MODULE ELQUAD16S__genmod
          INTERFACE 
            SUBROUTINE ELQUAD16S(NEL,NELP,KE,EL_DDISP,FINTE,INC)
              USE MOD_VARIABLES
              INTEGER(KIND=4) :: NEL
              INTEGER(KIND=4) :: NELP
              REAL(KIND=8), INTENT(OUT) :: KE((P+1)*(Q+1)*NDS,(P+1)*(Q+1&
     &)*NDS)
              REAL(KIND=8), INTENT(IN) :: EL_DDISP((P+1)*(Q+1)*NDS,1)
              REAL(KIND=8), INTENT(OUT) :: FINTE((P+1)*(Q+1)*NDS,1)
              INTEGER(KIND=4) :: INC
            END SUBROUTINE ELQUAD16S
          END INTERFACE 
        END MODULE ELQUAD16S__genmod
