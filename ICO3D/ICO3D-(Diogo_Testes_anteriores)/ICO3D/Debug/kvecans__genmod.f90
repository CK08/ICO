        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 14 08:53:46 2014
        MODULE KVECANS__genmod
          INTERFACE 
            SUBROUTINE KVECANS(NCPX,P,U_KNOT,UK,UB)
              INTEGER(KIND=4), INTENT(IN) :: P
              INTEGER(KIND=4), INTENT(IN) :: NCPX
              REAL(KIND=8), INTENT(IN) :: U_KNOT(NCPX+P+1)
              REAL(KIND=8), INTENT(IN) :: UK
              REAL(KIND=8), INTENT(OUT) :: UB(NCPX+P+1+P)
            END SUBROUTINE KVECANS
          END INTERFACE 
        END MODULE KVECANS__genmod
