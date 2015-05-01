!------------------------------------------------------------------------------
!
! Subroutine to calculate local axis based on the jacobian matrix
!
! Input: nds - problem dimension
!        jac - jacobian matrix
!
! Output: rconv - local axis
!
!------------------------------------------------------------------------------

subroutine local_axis(nds,jac,rconv)

    implicit none
    integer(4),intent(IN)::nds
    real(8),dimension(nds,nds),intent(IN)::jac
    real(8),dimension(nds,nds),intent(OUT)::rconv
    
    integer(4)::j
    real(8)::norm1,norm2,norm3
    real(8),dimension(nds,1)::r1,r2,r3
    
    r1=0.0d0
    r1(1,1) = jac(1,1)
    r1(2,1) = jac(1,2)
    r1(3,1) = jac(1,3)

    r2=0.0d0
    r2(1,1) = jac(2,1)
    r2(2,1) = jac(2,2)
    r2(3,1) = jac(2,3)

    !Normal vector
    r3=0.0d0
    call cross(r1,r2,r3)

    norm3=0.0d0
    do j=1,nds
        norm3=norm3 + r3(j,1)**2.0d0
    end do
    norm3=sqrt(norm3)

    r3=r3/norm3

    !Normalized tangent
    norm1=0.0d0
    norm2=0.0d0
    do j=1,nds
        norm1=norm1 + r1(j,1)**2.0d0
    end do
    norm1=sqrt(norm1)

    r1=r1/norm1

    !new local vector 2
    call cross(r3,r1,r2)

    rconv=0.0d0
    do j=1,nds
        rconv(j,1)=r1(j,1)
        rconv(j,2)=r2(j,1)
        rconv(j,3)=r3(j,1)
    end do

end subroutine