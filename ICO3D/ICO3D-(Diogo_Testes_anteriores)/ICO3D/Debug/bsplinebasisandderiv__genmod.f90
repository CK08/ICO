        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 14 08:53:46 2014
        MODULE BSPLINEBASISANDDERIV__genmod
          INTERFACE 
            SUBROUTINE BSPLINEBASISANDDERIV(NCP,P,U,U_KNOT,N,DNDXI)
              INTEGER(KIND=4), INTENT(IN) :: P
              INTEGER(KIND=4), INTENT(IN) :: NCP
              REAL(KIND=8), INTENT(IN) :: U
              REAL(KIND=8), INTENT(IN) :: U_KNOT(NCP+P+1)
              REAL(KIND=8) :: N(P+1)
              REAL(KIND=8) :: DNDXI(P+1)
            END SUBROUTINE BSPLINEBASISANDDERIV
          END INTERFACE 
        END MODULE BSPLINEBASISANDDERIV__genmod