        !COMPILER-GENERATED INTERFACE MODULE: Wed May 28 10:58:14 2014
        MODULE ASSEMBLY2D__genmod
          INTERFACE 
            SUBROUTINE ASSEMBLY2D(NZE,KE,FINTE)
              USE MOD_VARIABLES
              INTEGER(KIND=4), INTENT(IN) :: NZE
              REAL(KIND=8), INTENT(IN) :: KE((P+1)*(Q+1)*NDS,(P+1)*(Q+1)&
     &*NDS)
              REAL(KIND=8), INTENT(IN) :: FINTE((P+1)*(Q+1)*NDS,1)
            END SUBROUTINE ASSEMBLY2D
          END INTERFACE 
        END MODULE ASSEMBLY2D__genmod
