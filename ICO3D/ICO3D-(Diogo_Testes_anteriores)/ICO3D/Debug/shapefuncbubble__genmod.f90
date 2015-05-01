        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 14 08:53:54 2014
        MODULE SHAPEFUNCBUBBLE__genmod
          INTERFACE 
            SUBROUTINE SHAPEFUNCBUBBLE(NEL,XII,ETAI,ZETAI,R,DRDX,DETJ,  &
     &JAC,UDISP)
              USE MOD_VARIABLES
              INTEGER(KIND=4) :: NEL
              REAL(KIND=8), INTENT(IN) :: XII
              REAL(KIND=8), INTENT(IN) :: ETAI
              REAL(KIND=8), INTENT(IN) :: ZETAI
              REAL(KIND=8), INTENT(OUT) :: R((P+1)*(Q+1)*(W+1))
              REAL(KIND=8), INTENT(OUT) :: DRDX((P+1)*(Q+1)*(W+1),NDS)
              REAL(KIND=8), INTENT(OUT) :: DETJ
              REAL(KIND=8), INTENT(OUT) :: JAC(NDS,NDS)
              REAL(KIND=8), INTENT(IN) :: UDISP((P+1)*(Q+1)*(W+1)*NDS,1)
            END SUBROUTINE SHAPEFUNCBUBBLE
          END INTERFACE 
        END MODULE SHAPEFUNCBUBBLE__genmod
