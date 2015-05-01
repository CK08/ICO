        !COMPILER-GENERATED INTERFACE MODULE: Wed May 28 10:58:15 2014
        MODULE SHAPEFUNC2__genmod
          INTERFACE 
            SUBROUTINE SHAPEFUNC2(NEL,XII,ETAI,R,DRDX,JAC,DETJ,UDISP)
              USE MOD_VARIABLES
              INTEGER(KIND=4) :: NEL
              REAL(KIND=8), INTENT(IN) :: XII
              REAL(KIND=8), INTENT(IN) :: ETAI
              REAL(KIND=8), INTENT(OUT) :: R((P+1)*(Q+1))
              REAL(KIND=8), INTENT(OUT) :: DRDX((P+1)*(Q+1),NDS)
              REAL(KIND=8), INTENT(OUT) :: JAC(NDS,NDS)
              REAL(KIND=8), INTENT(OUT) :: DETJ
              REAL(KIND=8), INTENT(IN) :: UDISP((P+1)*(Q+1)*NDS,1)
            END SUBROUTINE SHAPEFUNC2
          END INTERFACE 
        END MODULE SHAPEFUNC2__genmod
