!------------------------------------------------------------------
!Aprod  computes  y = A*x  for some matrix  A.
!------------------------------------------------------------------

subroutine Aprod ( n, x, y )

    use Mod_Variables

    implicit           double precision (a-h,o-z)
    double precision   x(n), y(n)

    y = 0.0d0
!    do i=1,n
!        do j=1,n
!            if(Kfeq(i,j) .ne. 0.0d0) y(i) = y(i) + Kfeq(i,j)*x(j)
!        end do
!    end do

    do i=1,runner
        y(indx(i,1)) = y(indx(i,1)) + Kfeq(indx(i,1),indx(i,2))*x(indx(i,2))
    end do

end subroutine
      
      
      
      
      subroutine Msolve( n, x, y )

      implicit           double precision (a-h,o-z)
      double precision   x(n), y(n)

!     ------------------------------------------------------------------
!     Msolve  solves  M*y = x  for some symmetric positive-definite
!     matrix  M.
!     This is a simple example for testing  SYMMLQ.
!     shiftm will be the same as shift in SYMMLQ.
!
!     If pertbn = 0, the preconditioner will be exact, so
!     SYMMLQ should require either one or two iterations,
!     depending on whether (A - shift*I) is positive definite or not.
!
!     If pertbn is nonzero, somewhat more iterations will be required.
!     ------------------------------------------------------------------

      intrinsic          abs, mod
      common    /mshift/ shiftm, pertm
        
      do 10 i = 1, n
         d    = i ! 1.1
         d    = d / n
         d    = abs( d - shiftm )
         if (mod(i,10) .eq. 0) d = d + pertm
         y(i) = x(i) / d
   10 continue

!     end of Msolve
      end