!--------------------------------------------------------------------
!The subroutine solves the system of equations Av=vs
!--------------------------------------------------------------------
subroutine Gauss (dp,A,vs,v)
	implicit none

	integer, intent (in):: dp
	real*8, dimension (dp,dp), intent (inout):: A
	real*8, dimension (dp), intent (inout):: vs, v
	integer:: i, j, k
	real*8:: aux

	v=0.0
	do k=1, dp-1
		do i=k+1, dp
			aux=A(i,k)/A(k,k)
			do j=k+1, dp
				A(i,j)=A(i,j)-aux*A(k,j)
			end do
			vs(i)=vs(i)-aux*vs(k)
		end do
	end do

	do i=dp, 1, -1
		aux=0
		do j=i+1, dp
			aux=aux+A(i,j)*v(j)
		end do
		v(i)=(vs(i)-aux)/A(i,i)
	end do

end subroutine Gauss

!--------------------------------------------------------------------