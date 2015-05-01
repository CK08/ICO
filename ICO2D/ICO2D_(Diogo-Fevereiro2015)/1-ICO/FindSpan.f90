!--------------------------------------------------------------------------
!   Calculate knot span to which the data point u belongs
!   Based on the algorithm A2.1 from "The NURBS book"
!    
!   Input:
!       n = nknots - p - 1
!       p - B-Spline degree
!       u - data point
!       U_knot - knot vector
!   Output:
!       s - Knot span of data point u
!--------------------------------------------------------------------------

subroutine FindSpan(n,p,u,u_knot,s)
    
    implicit none
    
    integer(4),intent(IN)::n,p
    real(8),intent(IN)::u
    real(8),dimension(n),intent(IN)::u_knot
    
    integer(4),intent(OUT)::s
    
    integer(4)::low,high,mid

!Special Case - The point is the last knot in the knot span
    if(u == U_knot(n)) then
        s=n-p-1;
        goto 1
    end if

    low = p;
    high = n+1;
    mid = low + 1;

    do while(u .lt. U_knot(mid) .or. u >= U_knot(mid+1))
       if(u .lt. U_knot(mid)) then
           mid = mid - 1;
       else
           low = mid;
           mid = mid + 1;
       end if
    end do
        
    s = mid;
    
    1 continue
    
    return


end subroutine