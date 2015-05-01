subroutine local_axis_2D(nds,jac,rconv)

    implicit none
    integer(4),intent(IN)::nds
    real(8),dimension(nds,nds),intent(IN)::jac
    real(8),dimension(nds,nds),intent(OUT)::rconv
    
    integer(4)::j
    real(8)::norm1,norm2,norm3
    real(8),dimension(3,1)::r1,r2,r3
    
    r1=0.0d0
    r1(1,1) = jac(1,1)
    r1(2,1) = jac(1,2)
    r1(3,1) = 0.0d0

    r2=0.0d0
    r2(1,1) = jac(2,1)
    r2(2,1) = jac(2,2)
    r2(3,1) = 0.0d0

    !Normal vector
    r3=0.0d0
    call cross(r1,r2,r3)

    norm3=0.0d0
    do j=1,3
        norm3=norm3 + r3(j,1)**2.0d0
    end do
    norm3=sqrt(norm3)

    r3=r3/norm3

    !Normalized tangent
    norm1=0.0d0
    norm2=0.0d0
    do j=1,3
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
    end do

end subroutine