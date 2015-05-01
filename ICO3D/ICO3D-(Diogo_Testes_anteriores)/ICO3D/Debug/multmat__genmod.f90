        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 14 08:53:49 2014
        MODULE MULTMAT__genmod
          INTERFACE 
            SUBROUTINE MULTMAT(AMATR,BMATR,ROWA,COLA,COLB,CMATR)
              INTEGER(KIND=4) :: COLB
              INTEGER(KIND=4) :: COLA
              INTEGER(KIND=4) :: ROWA
              REAL(KIND=8) :: AMATR(ROWA,COLA)
              REAL(KIND=8) :: BMATR(COLA,COLB)
              REAL(KIND=8) :: CMATR(ROWA,COLB)
            END SUBROUTINE MULTMAT
          END INTERFACE 
        END MODULE MULTMAT__genmod
