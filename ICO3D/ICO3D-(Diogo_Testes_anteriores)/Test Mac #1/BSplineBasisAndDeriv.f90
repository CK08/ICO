!--------------------------------------------------------------------------------------
! Subroutne to calculate the basis functions and derivatives of a B-spline 
!
! From The NURBS book by Piegl and Tiller
!
! Input: ncp - number of control points
!        p - B-Spline degree
!        u - current data point
!        u_knot - knot vector
!
! Output: N - B-Splines basis funtions
!         dNdxi - B-Splines basis funtions first derivatives
!--------------------------------------------------------------------------------------

Subroutine BSplineBasisAndDeriv(ncp,p,u,u_knot,N,dNdxi)
    
    
    implicit none
    
    integer(4),intent(IN)::ncp,p
    real(8),intent(IN)::u
    real(8),dimension(ncp+p+1),intent(IN)::u_knot
    real(8),dimension(p+1,p+1)::dBF
    
    real(8),dimension(p+1)::N
    real(8),dimension(p+1)::dNdxi
    
    real(8),dimension(p+1,1)::BF
    
    integer(4)::i,nb,uspan,nder
    
    nb = ncp+p+1
    
    !Determine span index
    call FindSpan(nb,p,u,U_Knot,uspan)
    
    !Determine basis functions values
    !-call BasisFuns(nb,uspan,u,p,u_knot,BF);
    
    !Determine basis functions and derivatives
    !Line 1 of dBF contains the basis functions
    !Line i+1 contains the ith derivative of the basis functions
    nder = p
    call BasisFunsDers(nb,uspan,u,p,nder,u_knot,dBF)
    
    do i=1,p+1
        N(i) = dBF(1,i)
        if(p==0)then
            dNdxi(i) = 1.0d0
        else
            dNdxi(i) = dBF(2,i)
        endif
    end do
    
    
    continue
    
    
    
end subroutine BSplineBasisAndDeriv