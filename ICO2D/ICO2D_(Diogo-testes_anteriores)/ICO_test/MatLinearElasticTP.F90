    !-----------------------------------------------
    ! Linear Elastic material model for plane stress
    !-----------------------------------------------
    subroutine MatLinearElasticTP(E,ve,MatD)
    
    implicit real*8(a-h,o-z)

    real(8),intent(IN)::E,ve
    real(8),dimension(3,3),intent(OUT)::MatD
    real(8)::d1,d2,d3
    
    matD=0.0d0
    
    d1=E/(1.0d0-ve**2.0d0)
    d2=ve
    d3=(1.0d0-ve)/(2.0d0)
    
    matD(1,1)=1.0d0
    matD(1,2)=d2
    matD(2,1)=d2
    matD(2,2)=1.0d0
    matD(3,3)=d3

    matD=matD*d1
    
    end subroutine MatLinearElasticTP