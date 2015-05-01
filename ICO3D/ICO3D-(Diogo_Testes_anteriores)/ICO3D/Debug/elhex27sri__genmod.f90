        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 14 08:53:54 2014
        MODULE ELHEX27SRI__genmod
          INTERFACE 
            SUBROUTINE ELHEX27SRI(NEL,KE,EL_DDISP,FINTE)
              USE MOD_VARIABLES
              INTEGER(KIND=4) :: NEL
              REAL(KIND=8), INTENT(OUT) :: KE((P+1)*(Q+1)*(W+1)*NDS,(P+1&
     &)*(Q+1)*(W+1)*NDS)
              REAL(KIND=8), INTENT(IN) :: EL_DDISP((P+1)*(Q+1)*(W+1)*NDS&
     &,1)
              REAL(KIND=8), INTENT(OUT) :: FINTE((P+1)*(Q+1)*(W+1)*NDS,1&
     &)
            END SUBROUTINE ELHEX27SRI
          END INTERFACE 
        END MODULE ELHEX27SRI__genmod
