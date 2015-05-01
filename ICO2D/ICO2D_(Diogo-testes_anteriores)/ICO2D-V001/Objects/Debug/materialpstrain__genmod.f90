        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 03 16:45:32 2014
        MODULE MATERIALPSTRAIN__genmod
          INTERFACE 
            SUBROUTINE MATERIALPSTRAIN(IPROPS,NODEDOF,NNODES,PROPS,B,   &
     &BLINES,DDISP,STRESS,STRAIN,DSTRAIN,DSTRESS,HARD,MATD)
              INTEGER(KIND=4), INTENT(IN) :: BLINES
              INTEGER(KIND=4), INTENT(IN) :: NNODES
              INTEGER(KIND=4), INTENT(IN) :: NODEDOF
              INTEGER(KIND=4), INTENT(IN) :: IPROPS
              REAL(KIND=8), INTENT(IN) :: PROPS(IPROPS)
              REAL(KIND=8), INTENT(IN) :: B(BLINES,NODEDOF*NNODES)
              REAL(KIND=8), INTENT(IN) :: DDISP(NODEDOF*NNODES,1)
              REAL(KIND=8), INTENT(INOUT) :: STRESS(1,1,3)
              REAL(KIND=8), INTENT(INOUT) :: STRAIN(1,1,3)
              REAL(KIND=8), INTENT(INOUT) :: DSTRAIN(1,1,3)
              REAL(KIND=8), INTENT(INOUT) :: DSTRESS(1,1,3)
              REAL(KIND=8), INTENT(INOUT) :: HARD(1,1)
              REAL(KIND=8), INTENT(INOUT) :: MATD(3,3)
            END SUBROUTINE MATERIALPSTRAIN
          END INTERFACE 
        END MODULE MATERIALPSTRAIN__genmod
