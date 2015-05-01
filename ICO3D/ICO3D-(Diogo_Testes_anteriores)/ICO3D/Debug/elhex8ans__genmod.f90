        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 14 08:53:52 2014
        MODULE ELHEX8ANS__genmod
          INTERFACE 
            SUBROUTINE ELHEX8ANS(NEL,KE,EL_DDISP,FINTE,INC)
              USE MOD_VARIABLES
              INTEGER(KIND=4) :: NEL
              REAL(KIND=8), INTENT(OUT) :: KE((P+1)*(Q+1)*(W+1)*NDS,(P+1&
     &)*(Q+1)*(W+1)*NDS)
              REAL(KIND=8), INTENT(IN) :: EL_DDISP((P+1)*(Q+1)*(W+1)*NDS&
     &,1)
              REAL(KIND=8), INTENT(OUT) :: FINTE((P+1)*(Q+1)*(W+1)*NDS,1&
     &)
              INTEGER(KIND=4) :: INC
            END SUBROUTINE ELHEX8ANS
          END INTERFACE 
        END MODULE ELHEX8ANS__genmod
