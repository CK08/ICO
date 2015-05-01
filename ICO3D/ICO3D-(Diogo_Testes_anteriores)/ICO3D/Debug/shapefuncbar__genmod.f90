        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 14 08:53:48 2014
        MODULE SHAPEFUNCBAR__genmod
          INTERFACE 
            SUBROUTINE SHAPEFUNCBAR(NEL,XII,ETAI,ZETAI,INN,IEN,B_NET,   &
     &U_KNOT,V_KNOT,W_KNOT,P,Q,W,NCPX,NCPY,NCPZ,NDS,NNODES,NELEMS,R,DRDX&
     &,DETJ)
              INTEGER(KIND=4), INTENT(IN) :: NELEMS
              INTEGER(KIND=4), INTENT(IN) :: NNODES
              INTEGER(KIND=4), INTENT(IN) :: NDS
              INTEGER(KIND=4), INTENT(IN) :: NCPZ
              INTEGER(KIND=4), INTENT(IN) :: NCPY
              INTEGER(KIND=4), INTENT(IN) :: NCPX
              INTEGER(KIND=4), INTENT(IN) :: W
              INTEGER(KIND=4), INTENT(IN) :: Q
              INTEGER(KIND=4), INTENT(IN) :: P
              INTEGER(KIND=4) :: NEL
              REAL(KIND=8), INTENT(IN) :: XII
              REAL(KIND=8), INTENT(IN) :: ETAI
              REAL(KIND=8), INTENT(IN) :: ZETAI
              INTEGER(KIND=4), INTENT(IN) :: INN(NNODES,NDS)
              INTEGER(KIND=4), INTENT(IN) :: IEN(NELEMS,(P+1)*(Q+1)*(W+1&
     &))
              REAL(KIND=8), INTENT(IN) :: B_NET(NCPX,NCPY,NCPZ,NDS+1)
              REAL(KIND=8), INTENT(IN) :: U_KNOT(NCPX+P)
              REAL(KIND=8), INTENT(IN) :: V_KNOT(NCPY+Q)
              REAL(KIND=8), INTENT(IN) :: W_KNOT(NCPZ+W)
              REAL(KIND=8), INTENT(OUT) :: R((P+1)*(Q+1)*(W+1))
              REAL(KIND=8), INTENT(OUT) :: DRDX((P+1)*(Q+1)*(W+1),NDS)
              REAL(KIND=8), INTENT(OUT) :: DETJ
            END SUBROUTINE SHAPEFUNCBAR
          END INTERFACE 
        END MODULE SHAPEFUNCBAR__genmod
