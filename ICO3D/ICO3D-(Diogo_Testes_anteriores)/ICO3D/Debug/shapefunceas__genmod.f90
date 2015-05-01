        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb 20 09:46:58 2014
        MODULE SHAPEFUNCEAS__genmod
          INTERFACE 
            SUBROUTINE SHAPEFUNCEAS(NEL,NNODES,NDS,NELEMS,NSHPL,P,Q,W,  &
     &XII,ETAI,ZETAI,R,DRDX,DETJ,JAC,UDISP,POINTS,IEN)
              INTEGER(KIND=4) :: W
              INTEGER(KIND=4) :: Q
              INTEGER(KIND=4) :: P
              INTEGER(KIND=4) :: NSHPL
              INTEGER(KIND=4) :: NELEMS
              INTEGER(KIND=4) :: NDS
              INTEGER(KIND=4) :: NNODES
              INTEGER(KIND=4) :: NEL
              REAL(KIND=8), INTENT(IN) :: XII
              REAL(KIND=8), INTENT(IN) :: ETAI
              REAL(KIND=8), INTENT(IN) :: ZETAI
              REAL(KIND=8), INTENT(OUT) :: R((P+1)*(Q+1)*(W+1))
              REAL(KIND=8), INTENT(OUT) :: DRDX((P+1)*(Q+1)*(W+1),NDS)
              REAL(KIND=8), INTENT(OUT) :: DETJ
              REAL(KIND=8), INTENT(OUT) :: JAC(NDS,NDS)
              REAL(KIND=8), INTENT(IN) :: UDISP((P+1)*(Q+1)*(W+1)*NDS,1)
              REAL(KIND=8) :: POINTS(NNODES,NDS)
              INTEGER(KIND=4) :: IEN(NELEMS,NSHPL)
            END SUBROUTINE SHAPEFUNCEAS
          END INTERFACE 
        END MODULE SHAPEFUNCEAS__genmod
