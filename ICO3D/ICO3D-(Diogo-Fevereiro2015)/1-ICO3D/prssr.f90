!----------------------------------------------------------------------------------------------
!
! Subroutine compute pressure loads
!
! Warning: This subroutine may require further testing!
!
!----------------------------------------------------------------------------------------------

subroutine prssr(nel,ie,in,ic,jac,xi,eta,zeta,we,wn,wc,Finte)

    use Mod_Variables
    implicit none
    
    integer(4),intent(IN)::nel
    integer(4)::ie,in,ic
    real(8)::detj
    real(8),intent(IN)::xi,we
    real(8),intent(IN)::eta,wn
    real(8),intent(IN)::zeta,wc
    real(8),dimension(nds,nds)::jac
    
    real(8),dimension((p+1)*(q+1)*(w+1))::R
    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdx
    
    real(8),dimension(nshpl*nds,1),intent(INOUT)::Finte
    
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1)::zeros
    
    integer(4)::i,j,k,count
    character(3)::tag
    real(8)::mag
    
    real(8),dimension(3,1)::dir,vec1,vec2,vec3
    real(8),dimension(3,nshpl*nds)::Ni
    real(8),dimension(nshpl)::ShapeFunc
    real(8)::norm,det
    real(8),dimension(2,2)::mat
    
    zeros = 0.0d0
    
    vec1 = 0.0d0
    vec2 = 0.0d0
    vec3 = 0.0d0
    
    
    do i=1,nds
        vec1(i,1) = jac(i,1) 
        vec2(i,1) = jac(i,2)
        vec3(i,1) = jac(i,3)
    end do

    norm = 0.0d0
    
    Ni = 0.0d0
    
    tag = ftag(nel)
    mag = pressvec(nel) !/((p+1)*(q+1)*(w+1))
      
    jac = 0.0d0
    if(tag=='S1')then
        
        call ShapeFunc2(nel,1.0d0,eta,zeta,R,dRdx,detj,jac,zeros)
        
        do i=1,nds
            vec1(i,1) = jac(i,1) 
            vec2(i,1) = jac(i,2)
            vec3(i,1) = jac(i,3)
        end do
        
        !Surface normal vector
        call cross(vec2,vec3,dir)
        
        !Normal vector norm
        norm = 0.0d0
        do i=1,nds
            norm = norm + dir(i,1)*dir(i,1)
        end do
        
        norm = dsqrt(norm)
        
        mag = mag*norm
        
        !Normalized normal direction
        dir = dir/norm
        
        do i=1,(p+1)*(q+1)*(w+1)
            Ni(1,i*nds-2) = R(i)
            Ni(2,i*nds-1) = R(i)
            Ni(3,i*nds  ) = R(i)
        end do
        
        !Contribution to the internal forces
        Finte = Finte - matmul(transpose(Ni),dir)*mag/(ie)*wn*wc
        
        continue
        
    elseif(tag=='S2')then
        
        call ShapeFunc2(nel,xi,1.0d0,zeta,R,dRdx,detj,jac,zeros)
        
        do i=1,nds
            vec1(i,1) = jac(i,1) 
            vec2(i,1) = jac(i,2)
            vec3(i,1) = jac(i,3)
        end do
        
        !Surface normal vector
        call cross(vec3,vec1,dir)
        
        !Normal vector norm
        norm = 0.0d0
        do i=1,nds
            norm = norm + dir(i,1)*dir(i,1)
        end do
        
        norm = dsqrt(norm)
        
        mag = mag*norm
        
        !Normalized normal direction
        dir = dir/norm
        
        do i=1,(p+1)*(q+1)*(w+1)
            Ni(1,i*nds-2) = R(i)
            Ni(2,i*nds-1) = R(i)
            Ni(3,i*nds  ) = R(i)
        end do
        
        !Contribution to the internal forces
        Finte = Finte - matmul(transpose(Ni),dir)*mag/(in)*we*wc
        
    elseif(tag=='S3')then
        
        call ShapeFunc2(nel,-1.0d0,eta,zeta,R,dRdx,detj,jac,zeros)
        
        do i=1,nds
            vec1(i,1) = jac(i,1) 
            vec2(i,1) = jac(i,2)
            vec3(i,1) = jac(i,3)
        end do
        
        !Surface normal vector
        call cross(vec3,vec2,dir)
        
        !Normal vector norm
        norm = 0.0d0
        do i=1,nds
            norm = norm + dir(i,1)*dir(i,1)
        end do
        
        norm = dsqrt(norm)
        
        !Normalized normal direction
        dir = dir/norm
        
        mag = mag*norm
        
        do i=1,(p+1)*(q+1)*(w+1)
            Ni(1,i*nds-2) = R(i)
            Ni(2,i*nds-1) = R(i)
            Ni(3,i*nds  ) = R(i)
        end do
        
        !Contribution to the internal forces
        Finte = Finte - matmul(transpose(Ni),dir)*mag/(ie)*wn*wc 

    elseif(tag=='S4')then
        
        call ShapeFunc2(nel,xi,-1.0d0,zeta,R,dRdx,detj,jac,zeros)
        
        do i=1,nds
            vec1(i,1) = jac(i,1) 
            vec2(i,1) = jac(i,2)
            vec3(i,1) = jac(i,3)
        end do
        
        !Surface normal vector
        call cross(vec1,vec3,dir)
        
        !Normal vector norm
        norm = 0.0d0
        do i=1,nds
            norm = norm + dir(i,1)*dir(i,1)
        end do
        
        norm = dsqrt(norm)
        
        !Normalized normal direction
        dir = dir/norm
        
        mag = mag*norm
        
        do i=1,(p+1)*(q+1)*(w+1)
            Ni(1,i*nds-2) = R(i)
            Ni(2,i*nds-1) = R(i)
            Ni(3,i*nds  ) = R(i)
        end do
        
        !Contribution to the internal forces
        Finte = Finte - matmul(transpose(Ni),dir)*mag/(in)*we*wc
        
        continue

    elseif(tag=='S5')then
        
        call ShapeFunc2(nel,xi,eta,-1.0d0,R,dRdx,detj,jac,zeros)
        
        do i=1,nds
            vec1(i,1) = jac(i,1) 
            vec2(i,1) = jac(i,2)
            vec3(i,1) = jac(i,3)
        end do
        
        !Surfae normal vector
        call cross(vec2,vec1,dir)
        
        !Normal vector norm
        norm = 0.0d0
        do i=1,nds
            norm = norm + dir(i,1)*dir(i,1)
        end do
        
        norm = dsqrt(norm)
        
        !Normalized normal direction
        dir = dir/norm
        
        mag = mag*norm
        
        do i=1,(p+1)*(q+1)*(w+1)
            Ni(1,i*nds-2) = R(i)
            Ni(2,i*nds-1) = R(i)
            Ni(3,i*nds  ) = R(i)
        end do
        
        !Contribution to the internal forces
        Finte = Finte - matmul(transpose(Ni),dir)*mag/(ic)*we*wn
        
        continue
    

    elseif(tag=='S6')then
        
        call ShapeFunc2(nel,xi,eta,1.0d0,R,dRdx,detj,jac,zeros)
        
        do i=1,nds
            vec1(i,1) = jac(i,1) 
            vec2(i,1) = jac(i,2)
            vec3(i,1) = jac(i,3)
        end do
        
        !Surfae normal vector
        call cross(vec1,vec2,dir)
        
        !Normal vector norm
        norm = 0.0d0
        do i=1,nds
            norm = norm + dir(i,1)*dir(i,1)
        end do
        
        norm = dsqrt(norm)
        
        !Normalized normal direction
        dir = dir/norm
        
        mag = mag*norm
        
        do i=1,(p+1)*(q+1)*(w+1)
            Ni(1,i*nds-2) = R(i)
            Ni(2,i*nds-1) = R(i)
            Ni(3,i*nds  ) = R(i)
        end do
        
        !Contribution to the internal forces
        Finte = Finte - matmul(transpose(Ni),dir)*mag/(ic)*we*wn

    end if
    
    
end subroutine