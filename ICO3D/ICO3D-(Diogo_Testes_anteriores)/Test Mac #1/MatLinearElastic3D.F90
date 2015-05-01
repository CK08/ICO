!----------------------------------------------------------------------------
!
! Linear Elastic material model for 3D elements
!
! Input: E - Young Modulus
!        ve - Poissons coefficient
!
! Output: MatD - Linear isotropic constitutive matrix for 3D applications 
!----------------------------------------------------------------------------
    subroutine MatLinearElastic3D(E,ve,MatD)
    
    implicit none !real*8(a-h,o-z)

    real(8),intent(IN)::E,ve
    real(8),dimension(6,6),intent(OUT)::MatD
    real(8)::d1,d2,d3
    
    matD=0.0d0
    
    d1=E*(1.0d0-ve)/((1.0d0-2.0d0*ve)*(1.0d0+ve))
    d2=ve/(1.0d0-ve)
    d3=(1.0d0-2.0d0*ve)/(2.0d0*(1.0d0-ve))
    
    matD(1,1)=1.0d0
    matD(1,2)=d2
    matD(1,3)=d2
    matD(2,1)=d2
    matD(2,2)=1.0d0
    matD(2,3)=d2 
    matD(3,1)=d2
    matD(3,2)=d2
    matD(3,3)=1.0d0
    matD(4,4)=d3
    matD(5,5)=d3
    matD(6,6)=d3
    

    matD=matD*d1
    
    end subroutine MatLinearElastic3D