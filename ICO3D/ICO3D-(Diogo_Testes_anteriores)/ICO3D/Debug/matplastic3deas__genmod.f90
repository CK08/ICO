        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 14 08:53:46 2014
        MODULE MATPLASTIC3DEAS__genmod
          INTERFACE 
            SUBROUTINE MATPLASTIC3DEAS(IPROPS,NODEDOF,NNODES,PROPS,B,   &
     &BALPHA,BLINES,NALPHA,DDISP,DALPHA,STRESS,STRAIN,DSTRAIN,DSTRESS,  &
     &HARD,MATD)
              INTEGER(KIND=4), INTENT(IN) :: NALPHA
              INTEGER(KIND=4), INTENT(IN) :: BLINES
              INTEGER(KIND=4), INTENT(IN) :: NNODES
              INTEGER(KIND=4), INTENT(IN) :: NODEDOF
              INTEGER(KIND=4), INTENT(IN) :: IPROPS
              REAL(KIND=8), INTENT(IN) :: PROPS(IPROPS)
              REAL(KIND=8), INTENT(IN) :: B(BLINES,NODEDOF*NNODES)
              REAL(KIND=8), INTENT(IN) :: BALPHA(BLINES,NALPHA)
              REAL(KIND=8), INTENT(IN) :: DDISP(NODEDOF*NNODES,1)
              REAL(KIND=8), INTENT(IN) :: DALPHA(NALPHA,1)
              REAL(KIND=8), INTENT(INOUT) :: STRESS(1,1,6)
              REAL(KIND=8), INTENT(INOUT) :: STRAIN(1,1,6)
              REAL(KIND=8), INTENT(INOUT) :: DSTRAIN(1,1,6)
              REAL(KIND=8), INTENT(INOUT) :: DSTRESS(1,1,6)
              REAL(KIND=8), INTENT(INOUT) :: HARD(1,1)
              REAL(KIND=8), INTENT(INOUT) :: MATD(6,6)
            END SUBROUTINE MATPLASTIC3DEAS
          END INTERFACE 
        END MODULE MATPLASTIC3DEAS__genmod
