        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 14 08:53:52 2014
        MODULE ASSEMBLY3D__genmod
          INTERFACE 
            SUBROUTINE ASSEMBLY3D(NZE,KE,FINTE)
              USE MOD_VARIABLES
              INTEGER(KIND=4), INTENT(IN) :: NZE
              REAL(KIND=8), INTENT(IN) :: KE((P+1)*(Q+1)*(W+1)*NDS,(P+1)&
     &*(Q+1)*(W+1)*NDS)
              REAL(KIND=8), INTENT(IN) :: FINTE((P+1)*(Q+1)*(W+1)*NDS,1)
            END SUBROUTINE ASSEMBLY3D
          END INTERFACE 
        END MODULE ASSEMBLY3D__genmod
