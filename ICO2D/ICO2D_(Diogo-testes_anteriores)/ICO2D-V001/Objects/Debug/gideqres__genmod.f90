        !COMPILER-GENERATED INTERFACE MODULE: Wed May 28 10:58:16 2014
        MODULE GIDEQRES__genmod
          INTERFACE 
            SUBROUTINE GIDEQRES(INC,COORDI)
              USE MOD_VARIABLES
              INTEGER(KIND=4), INTENT(IN) :: INC
              REAL(KIND=8), INTENT(IN) :: COORDI(NCPX*NCPY,NDS)
            END SUBROUTINE GIDEQRES
          END INTERFACE 
        END MODULE GIDEQRES__genmod
