!--------------------------------------------------------------------------
!   Subroutine to calculate derivatives of the basis functions 
!   of a B-Spline based on the algorithm A2.3 from "The NURBS book"
!    
!   Input:
!       nb - number of basis functions
!       s - knot span
!       u - data point
!       p - B-Spline degree
!       n - nth derivative
!       U_knot - knot vector
!   Output:
!       ders(n+1,p+1) - Derivative of the basis functions
!                       Line 1 contains the basis functions
!                       Line i + 1 contains the ith derivative
!--------------------------------------------------------------------------

subroutine BasisFunsDers(nb,s,u,p,n,U_knot,ders)
    
    implicit none
    
    integer(4),intent(IN)::nb,s,p,n
    real(8),intent(IN)::u
    real(8),dimension(nb),intent(IN)::u_knot
    real(8),dimension(n+1,p+1),intent(OUT)::ders
    
    integer(4)::j,k,r,s1,s2,rk,pk,j1,j2
    real(8)::saved,temp,d
    real(8),dimension(p+1,p+1)::ndu
    real(8),dimension(p+1)::left,right
    !real(8),dimension(2,p+1)::a
    real(8),dimension(p+1,p+1)::a
    
    ders = 0.0d0
    ndu = 0.0d0
    left = 0.0d0
    right = 0.0d0
    a = 0.0d0

    ndu(1,1) = 1.0;

    do j=1,p
        left(j+1) = u - u_knot(s+1-j)
        right(j+1) = u_knot(s+j) - u
        saved = 0.0
        
        do r=0,j-1
            ndu(j+1,r+1) = right(r+2)+left(j-r+1)
            temp = ndu(r+1,j)/ndu(j+1,r+1)
            
            ndu(r+1,j+1) = saved + right(r+2)*temp
            saved = left(j-r+1)*temp
        end do
        ndu(j+1,j+1) = saved
        continue
    end do

    !Load basis functions
    do j=0,p
        ders(1,j+1) = ndu(j+1,p+1)
    end do

    do r=0,p
       s1 = 0
       s2 = 1
       a(1,1) = 1.0
       
       do k=1,n
           d =0.0
           rk = r-k
           pk = p-k
           
           if(r>=k)then
               a(s2+1,1) = a(s1+1,1)/ndu(pk+2,rk+1)
               d = a(s2+1,1)*ndu(rk+1,pk+1)
           end if
           
           if(rk>= -1)then
               j1 = 1
           else
               j1 = -rk
           end if
           
           if((r-1) <= pk) then
               j2 = k-1
           else
               j2 = p-r
           end if
          
           do j=j1,j2
               a(s2+1,j+1) = (a(s1+1,j+1)-a(s1+1,j))/ndu(pk+2,rk+j+1)
               d = d + a(s2+1,j+1)*ndu(rk+j+1,pk+1)
           end do
            
           if(r<=pk) then
               a(s2+1,k+1) = -a(s1+1,k)/ndu(pk+2,r+1)
               d = d + a(s2+1,k+1)*ndu(r+1,pk+1)
           end if
           
           ders(k+1,r+1) = d
           j = s1
           s1 = s2
           s2 = j + 1
            
            
       end do
       
    end do

    r = p;
    do k = 1,n
        do j = 0,p
            ders(k+1,j+1) = ders(k+1,j+1)*r
        end do
        r = r*(p-k);
    end do

    return
    
    
end subroutine BasisFunsDers



subroutine BasisFunsDers2(nb,s,u,p,n,U_knot,ders)
    
    implicit none
    
    integer(4),intent(IN)::nb,s,p,n
    real(8),intent(IN)::u
    real(8),dimension(nb),intent(IN)::u_knot
    real(8),dimension(n+1,nb),intent(OUT)::ders
    
    integer(4)::j,k,r,s1,s2,rk,pk,j1,j2
    real(8)::saved,temp,d
    real(8),dimension(p+1,p+1)::ndu
    real(8),dimension(p+1)::left,right
    !real(8),dimension(2,p+1)::a
    real(8),dimension(p+1,p+1)::a
    
    ders = 0.0d0
    ndu = 0.0d0
    left = 0.0d0
    right = 0.0d0
    a = 0.0d0

    ndu(1,1) = 1.0;

    do j=1,p
        left(j+1) = u - u_knot(s+1-j)
        right(j+1) = u_knot(s+j) - u
        saved = 0.0
        
        do r=0,j-1
            ndu(j+1,r+1) = right(r+2)+left(j-r+1)
            temp = ndu(r+1,j)/ndu(j+1,r+1)
            
            ndu(r+1,j+1) = saved + right(r+2)*temp
            saved = left(j-r+1)*temp
        end do
        ndu(j+1,j+1) = saved
        continue
    end do

    !Load basis functions
    do j=0,p
        ders(1,j+1) = ndu(j+1,p+1)
    end do

    do r=0,p
       s1 = 0
       s2 = 1
       a(1,1) = 1.0
       
       do k=1,n
           d =0.0
           rk = r-k
           pk = p-k
           
           if(r>=k)then
               a(s2+1,1) = a(s1+1,1)/ndu(pk+2,rk+1)
               d = a(s2+1,1)*ndu(rk+1,pk+1)
           end if
           
           if(rk>= -1)then
               j1 = 1
           else
               j1 = -rk
           end if
           
           if((r-1) <= pk) then
               j2 = k-1
           else
               j2 = p-r
           end if
          
           do j=j1,j2
               a(s2+1,j+1) = (a(s1+1,j+1)-a(s1+1,j))/ndu(pk+2,rk+j+1)
               d = d + a(s2+1,j+1)*ndu(rk+j+1,pk+1)
           end do
            
           if(r<=pk) then
               a(s2+1,k+1) = -a(s1+1,k)/ndu(pk+2,r+1)
               d = d + a(s2+1,k+1)*ndu(r+1,pk+1)
           end if
           
           ders(k+1,r+1) = d
           j = s1
           s1 = s2
           s2 = j + 1
            
            
       end do
       
    end do

    r = p;
    do k = 1,n
        do j = 0,p
            ders(k+1,j+1) = ders(k+1,j+1)*r
        end do
        r = r*(p-k);
    end do

    return
    
    
end subroutine BasisFunsDers2