!----------------------------------------------------------------------------------------------
!
! Subroutine for precomputing the L matrix employed in the ANS method
!
! The analytical derivation of the L matrix was proposed by Josef Kiendl
!
!----------------------------------------------------------------------------------------------
subroutine PreCompANS()

    use mod_variables
    implicit none
    real(8)::pa,pb,r1,r2,r3,r4,r5,CA,CB
    
    allocate(MatL1(9,6),MatL2(9,6))
    allocate(MatL3(9,4))
    
    allocate(TyPtANS(16,nds))
    
    !ANS Tying Points -----
    CA = dsqrt(1.0d0/3.0d0)
    CB = dsqrt(3.0d0/5.0d0)
    
    TyPtANS = 0.0d0
    TyPtANS(1,1)  =    CA; TyPtANS(1,2)  =    CB;
    TyPtANS(2,1)  =   -CA; TyPtANS(2,2)  =    CB;
    TyPtANS(3,1)  =    CA; TyPtANS(3,2)  = 0.0d0;
    TyPtANS(4,1)  =   -CA; TyPtANS(4,2)  = 0.0d0;
    TyPtANS(5,1)  =    CA; TyPtANS(5,2)  =   -CB;
    TyPtANS(6,1)  =   -CA; TyPtANS(6,2)  =   -CB;

    TyPtANS(7,1)  =    CB; TyPtANS(7,2)  =    CA;
    TyPtANS(8,1)  = 0.0d0; TyPtANS(8,2)  =    CA;
    TyPtANS(9,1)  =   -CB; TyPtANS(9,2)  =    CA;
    TyPtANS(10,1) =    CB; TyPtANS(10,2) =   -CA;
    TyPtANS(11,1) = 0.0d0; TyPtANS(11,2) =   -CA;
    TyPtANS(12,1) =   -CB; TyPtANS(12,2) =   -CA;

    TyPtANS(13,1) =  CA; TyPtANS(13,2) =  CA
    TyPtANS(14,1) =  CA; TyPtANS(14,2) = -CA
    TyPtANS(15,1) = -CA; TyPtANS(15,2) =  CA
    TyPtANS(16,1) = -CA; TyPtANS(16,2) = -CA
    
    pa = 1.0d0/2.0d0 - dsqrt(9.0d0/20.0d0)
    pb = 1.0d0/2.0d0 + dsqrt(9.0d0/20.0d0)

    MatL1 = 0.0d0
    MatL1(1,5) = pa
    MatL1(1,6) = pb
    MatL1(2,5) = 0.5d0
    MatL1(2,6) = 0.5d0
    MatL1(3,5) = pb
    MatL1(3,6) = pa

    MatL1(4,3) = pa
    MatL1(4,4) = pb
    MatL1(5,3) = 0.5d0
    MatL1(5,4) = 0.5d0
    MatL1(6,3) = pb
    MatL1(6,4) = pa

    MatL1(7,1) = pa
    MatL1(7,2) = pb
    MatL1(8,1) = 0.5d0
    MatL1(8,2) = 0.5d0
    MatL1(9,1) = pb
    MatL1(9,2) = pa

    !-------------------------------------
    MatL2 = 0.0d0
    MatL2(1,3) = pa
    MatL2(1,6) = pb
    MatL2(2,2) = pa
    MatL2(2,5) = pb
    MatL2(3,1) = pa
    MatL2(3,4) = pb

    MatL2(4,3) = 0.5d0
    MatL2(4,6) = 0.5d0
    MatL2(5,2) = 0.5d0
    MatL2(5,5) = 0.5d0
    MatL2(6,1) = 0.5d0
    MatL2(6,4) = 0.5d0

    MatL2(7,3) = pb
    MatL2(7,6) = pa
    MatL2(8,2) = pb
    MatL2(8,5) = pa
    MatL2(9,1) = pb
    MatL2(9,4) = pa

    !-------------------------------------
    r1 = 14.0d0/5.0d0 - dsqrt(36.0/5.0d0)
    r2 = 14.0d0/5.0d0 + dsqrt(36.0/5.0d0)

    r3 = -4.0d0/5.0d0

    r4 = 1.0d0 - dsqrt(9.0/5.0d0)
    r5 = 1.0d0 + dsqrt(9.0/5.0d0)

    MatL3 = 0.0d0
    MatL3(1,1) = r1
    MatL3(1,2) = r3
    MatL3(1,3) = r3
    MatL3(1,4) = r2

    MatL3(2,1) = r4
    MatL3(2,2) = r5
    MatL3(2,3) = r4
    MatL3(2,4) = r5

    MatL3(3,1) = r3 
    MatL3(3,2) = r2
    MatL3(3,3) = r1
    MatL3(3,4) = r3

    MatL3(4,1) = r4
    MatL3(4,2) = r4
    MatL3(4,3) = r5
    MatL3(4,4) = r5

    MatL3(5,1) = 1.0d0
    MatL3(5,2) = 1.0d0
    MatL3(5,3) = 1.0d0
    MatL3(5,4) = 1.0d0

    MatL3(6,1) = r5
    MatL3(6,2) = r5
    MatL3(6,3) = r4
    MatL3(6,4) = r4

    MatL3(7,1) = r3
    MatL3(7,2) = r1
    MatL3(7,3) = r2
    MatL3(7,4) = r3

    MatL3(8,1) = r5
    MatL3(8,2) = r4
    MatL3(8,3) = r5
    MatL3(8,4) = r4

    MatL3(9,1) = r2
    MatL3(9,2) = r3
    MatL3(9,3) = r3
    MatL3(9,4) = r1

    MatL3 = MatL3/4.0d0

end subroutine