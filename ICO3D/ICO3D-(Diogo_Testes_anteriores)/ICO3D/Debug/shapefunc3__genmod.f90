        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb 20 09:46:58 2014
        MODULE SHAPEFUNC3__genmod
          INTERFACE 
            SUBROUTINE SHAPEFUNC3(NEL,XII,ETAI,ZETAI,R,DRDX,DRDXI,DRDXII&
     &,DETJ,JAC,DXDXI,UDISP)
              USE MOD_VARIABLES
              INTEGER(KIND=4) :: NEL
              REAL(KIND=8), INTENT(IN) :: XII
              REAL(KIND=8), INTENT(IN) :: ETAI
              REAL(KIND=8), INTENT(IN) :: ZETAI
              REAL(KIND=8), INTENT(OUT) :: R((P+1)*(Q+1)*(W+1))
              REAL(KIND=8), INTENT(OUT) :: DRDX((P+1)*(Q+1)*(W+1),NDS)
              REAL(KIND=8), INTENT(OUT) :: DRDXI((P+1)*(Q+1)*(W+1),NDS)
              REAL(KIND=8), INTENT(OUT) :: DRDXII((P+1)*(Q+1)*(W+1),NDS)
              REAL(KIND=8), INTENT(OUT) :: DETJ
              REAL(KIND=8), INTENT(OUT) :: JAC(NDS,NDS)
              REAL(KIND=8), INTENT(OUT) :: DXDXI(NDS,NDS)
              REAL(KIND=8), INTENT(IN) :: UDISP((P+1)*(Q+1)*(W+1)*NDS,1)
            END SUBROUTINE SHAPEFUNC3
          END INTERFACE 
        END MODULE SHAPEFUNC3__genmod
