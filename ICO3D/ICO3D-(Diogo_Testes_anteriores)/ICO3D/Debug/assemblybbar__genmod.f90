        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 14 09:09:31 2014
        MODULE ASSEMBLYBBAR__genmod
          INTERFACE 
            SUBROUTINE ASSEMBLYBBAR(NZE,MA,MC,MV)
              USE MOD_VARIABLES
              INTEGER(KIND=4), INTENT(IN) :: NZE
              REAL(KIND=8), INTENT(IN) :: MA((P+1)*(Q+1)*(W+1)*NDS,(P+1)&
     &*(Q+1)*(W+1)*NDS)
              REAL(KIND=8), INTENT(IN) :: MC((P+1)*(Q+1)*(W+1)*NDS,P*Q*W&
     &)
              REAL(KIND=8), INTENT(IN) :: MV(P*Q*W,P*Q*W)
            END SUBROUTINE ASSEMBLYBBAR
          END INTERFACE 
        END MODULE ASSEMBLYBBAR__genmod
