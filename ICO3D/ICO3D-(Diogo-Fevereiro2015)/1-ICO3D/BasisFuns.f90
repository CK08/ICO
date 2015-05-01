!--------------------------------------------------------------------------
!   Calculate basis functions of a B-Spline
!   Based on the algorithm A2.2 from "The NURBS book"
!    
!   Input:
!       s - knot span
!       u - data point
!       p - B-Spline degree
!       u_knots - knot vector
!   Output:
!       BF(p+1,1) - Non-zero basis functions
!--------------------------------------------------------------------------

subroutine BasisFuns(nb,s,u,p,u_knots,BF)
    
    implicit none
    
    integer(4),intent(IN)::s,p,nb
    real(8),intent(IN)::u
    real(8),dimension(nb),intent(IN)::u_knots
    real(8),dimension(p+1,1),intent(OUT)::BF
    
    integer(4)::i,j,r,k1
    real(8)::saved,temp
    real(8),dimension(p+1)::N,left,right
    

    N = 0.0d0
    BF = 0.0d0
    left = 0.0d0
    right = 0.0d0

    N(1) = 1.0;

    i = s;

    do j=1,p
        left(j+1) = u - u_knots(i+1-j)
        right(j+1) = u_knots(i+j) - u
        saved = 0.0
        
        do r=0,j-1
            temp = N(r+1)/(right(r+2)+left(j-r+1))
            N(r+1) = saved+right(r+2)*temp
            saved = left(j-r+1)*temp
        end do
        
        N(j+1) = saved;
        continue
        
    end do

    do k1=1,p+1
        BF(k1,1) = N(k1);
    end do

    return

    
    
    
end subroutine BasisFuns