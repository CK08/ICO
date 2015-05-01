        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 14 08:53:47 2014
        MODULE GRVT__genmod
          INTERFACE 
            SUBROUTINE GRVT(ISHP,NDS,R,GWT,GCONST,GRAVDIR,DENSITY,FINTE)
              INTEGER(KIND=4), INTENT(IN) :: NDS
              INTEGER(KIND=4), INTENT(IN) :: ISHP
              REAL(KIND=8), INTENT(IN) :: R(ISHP)
              REAL(KIND=8), INTENT(IN) :: GWT
              REAL(KIND=8), INTENT(IN) :: GCONST
              INTEGER(KIND=4), INTENT(IN) :: GRAVDIR
              REAL(KIND=8), INTENT(IN) :: DENSITY
              REAL(KIND=8), INTENT(INOUT) :: FINTE(ISHP*NDS,1)
            END SUBROUTINE GRVT
          END INTERFACE 
        END MODULE GRVT__genmod
