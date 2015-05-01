        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 14 08:53:51 2014
        MODULE ELHEX27EAS_PW__genmod
          INTERFACE 
            SUBROUTINE ELHEX27EAS_PW(NEL,KE,EL_DDISP,FINTE,MATA,MATC,   &
     &MATV)
              USE MOD_VARIABLES
              INTEGER(KIND=4) :: NEL
              REAL(KIND=8), INTENT(OUT) :: KE((P+1)*(Q+1)*(W+1)*NDS,(P+1&
     &)*(Q+1)*(W+1)*NDS)
              REAL(KIND=8), INTENT(IN) :: EL_DDISP((P+1)*(Q+1)*(W+1)*NDS&
     &,1)
              REAL(KIND=8), INTENT(OUT) :: FINTE((P+1)*(Q+1)*(W+1)*NDS,1&
     &)
              REAL(KIND=8), INTENT(OUT) :: MATA((P+1)*(Q+1)*(W+1)*NDS,(P&
     &+1)*(Q+1)*(W+1)*NDS)
              REAL(KIND=8), INTENT(OUT) :: MATC((P+1)*(Q+1)*(W+1)*NDS,  &
     &NALPHA)
              REAL(KIND=8), INTENT(OUT) :: MATV(NALPHA,NALPHA)
            END SUBROUTINE ELHEX27EAS_PW
          END INTERFACE 
        END MODULE ELHEX27EAS_PW__genmod
