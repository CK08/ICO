        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 14 08:53:53 2014
        MODULE PRSSR__genmod
          INTERFACE 
            SUBROUTINE PRSSR(NEL,IE,IN,IC,JAC,XI,ETA,ZETA,WE,WN,WC,FINTE&
     &)
              USE MOD_VARIABLES
              INTEGER(KIND=4), INTENT(IN) :: NEL
              INTEGER(KIND=4) :: IE
              INTEGER(KIND=4) :: IN
              INTEGER(KIND=4) :: IC
              REAL(KIND=8) :: JAC(NDS,NDS)
              REAL(KIND=8), INTENT(IN) :: XI
              REAL(KIND=8), INTENT(IN) :: ETA
              REAL(KIND=8), INTENT(IN) :: ZETA
              REAL(KIND=8), INTENT(IN) :: WE
              REAL(KIND=8), INTENT(IN) :: WN
              REAL(KIND=8), INTENT(IN) :: WC
              REAL(KIND=8), INTENT(INOUT) :: FINTE(NSHPL*NDS,1)
            END SUBROUTINE PRSSR
          END INTERFACE 
        END MODULE PRSSR__genmod
