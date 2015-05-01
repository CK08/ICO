        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 14 08:53:50 2014
        MODULE ELHEX8BBAR__genmod
          INTERFACE 
            SUBROUTINE ELHEX8BBAR(NEL,KE,EL_DDISP,FINTE,MATA,MATC,MATM)
              USE MOD_VARIABLES
              INTEGER(KIND=4) :: NEL
              REAL(KIND=8), INTENT(OUT) :: KE((P+1)*(Q+1)*(W+1)*NDS,(P+1&
     &)*(Q+1)*(W+1)*NDS)
              REAL(KIND=8), INTENT(IN) :: EL_DDISP((P+1)*(Q+1)*(W+1)*NDS&
     &,1)
              REAL(KIND=8), INTENT(OUT) :: FINTE((P+1)*(Q+1)*(W+1)*NDS,1&
     &)
              REAL(KIND=8) :: MATA((P+1)*(Q+1)*(W+1)*NDS,(P+1)*(Q+1)*(W+&
     &1)*NDS)
              REAL(KIND=8) :: MATC((P+1)*(Q+1)*(W+1)*NDS,P*Q*W)
              REAL(KIND=8) :: MATM(P*Q*W,P*Q*W)
            END SUBROUTINE ELHEX8BBAR
          END INTERFACE 
        END MODULE ELHEX8BBAR__genmod
