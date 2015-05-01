        !COMPILER-GENERATED INTERFACE MODULE: Thu May 29 14:18:25 2014
        MODULE BSPLINESTONURBS__genmod
          INTERFACE 
            SUBROUTINE BSPLINESTONURBS(NCP,PS,XIB,KNOT_VEC,RMM,DRMM,    &
     &DDRMM,KSPAN,CONNECT,RM,DRMDXI,D2RMDXI2)
              INTEGER(KIND=4) :: PS
              INTEGER(KIND=4) :: NCP
              REAL(KIND=8) :: XIB
              REAL(KIND=8) :: KNOT_VEC(NCP+PS+1)
              REAL(KIND=8) :: RMM(PS+1)
              REAL(KIND=8) :: DRMM(PS+1)
              REAL(KIND=8) :: DDRMM(PS+1)
              INTEGER(KIND=4) :: KSPAN
              INTEGER(KIND=4) :: CONNECT(NCP)
              REAL(KIND=8) :: RM(NCP)
              REAL(KIND=8) :: DRMDXI(NCP)
              REAL(KIND=8) :: D2RMDXI2(NCP)
            END SUBROUTINE BSPLINESTONURBS
          END INTERFACE 
        END MODULE BSPLINESTONURBS__genmod
