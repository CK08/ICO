!--------------------------------------------------------------
! Subroutine to calculate the inverse of a matrix
!
! Input: A(n,n) - Matrix to be inverted
!        n      - matrix dimension
!        B(n)   - Array of zeros 
!        ierr   - error indicator
! Output: A(n,n) - Inverse of the matrix used in the input
!--------------------------------------------------------------

SUBROUTINE GAUSSJ(A,N,B,IERR)
	IMPLICIT REAL(8) (A-H,O-Z)
	INTEGER::N
	REAL(8),DIMENSION(N,N)::A
	REAL(8),DIMENSION(N)::B
	INTEGER,PARAMETER::NMAX=50
	INTEGER,DIMENSION(N)::INDXC,INDXR,IPIV
	IERR=0
!	IF(N.GT.NMAX)THEN
!	IERR=2
!	RETURN
!	ENDIF
	DO 11 J=1,N
	IPIV(J)=0
	11 CONTINUE
	DO 22 I=1,N
	BIG=0.
	DO 13 J=1,N
	IF(IPIV(J).NE.1)THEN
	DO 12 K=1,N
	IF(IPIV(K).EQ.0) THEN
	IF(ABS(A(J,K)).GE.BIG)THEN
	BIG=ABS(A(J,K))
	IROW=J
	ICOL=K
	ENDIF
	ELSE IF (IPIV(K).GT.1) THEN
	IERR=1
	RETURN
	ENDIF
	12 CONTINUE
	ENDIF
	13 CONTINUE
	IPIV(ICOL)=IPIV(ICOL)+1
	IF(IROW.NE.ICOL) THEN
	DO 14 L=1,N
	DUM=A(IROW,L)
	A(IROW,L)=A(ICOL,L)
	A(ICOL,L)=DUM
	14 CONTINUE
	DUM=B(IROW)
	B(IROW)=B(ICOL)
	B(ICOL)=DUM
	ENDIF
	INDXR(I)=IROW
	INDXC(I)=ICOL
	IF(A(ICOL,ICOL).EQ.0.)THEN
	IERR=1
	ENDIF
	PIVINV=1./A(ICOL,ICOL)
	A(ICOL,ICOL)=1.
	DO 16 L=1,N
	A(ICOL,L)=A(ICOL,L)*PIVINV
	16 CONTINUE
	B(ICOL)=B(ICOL)*PIVINV
	DO 21 LL=1,N
	IF(LL.NE.ICOL)THEN
	DUM=A(LL,ICOL)
	A(LL,ICOL)=0.
	DO 18 L=1,N
	A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
	18 CONTINUE
	B(LL)=B(LL)-B(ICOL)*DUM
	ENDIF
	21 CONTINUE
	22 CONTINUE
	DO 24 L=N,1,-1
	IF(INDXR(L).NE.INDXC(L))THEN
	DO 23 K=1,N
	DUM=A(K,INDXR(L))
	A(K,INDXR(L))=A(K,INDXC(L))
	A(K,INDXC(L))=DUM
	23 CONTINUE
	ENDIF
	24 CONTINUE
END SUBROUTINE GAUSSJ

!----------------------------------------------------------------------------------------------
!Subroutine to calculate the determinant of a 3x3 matrix
!----------------------------------------------------------------------------------------------
    subroutine fdeterm(amatr,deter)
        
      real(8),dimension(3,3),intent(IN)::amatr
      real(8),intent(OUT)::deter

      deter = (amatr(1,1)*amatr(2,2)*amatr(3,3)) +&
     &                    (amatr(1,2)*amatr(2,3)*amatr(3,1)) +&
     &                    (amatr(1,3)*amatr(2,1)*amatr(3,2)) -&
     &                    (amatr(1,3)*amatr(2,2)*amatr(3,1)) -&
     &                    (amatr(1,1)*amatr(2,3)*amatr(3,2)) -&
     &                    (amatr(1,2)*amatr(2,1)*amatr(3,3))

      end subroutine
      
      
!------------------------------------------------------------------------------------------------------
!Subroutine to calculate the cross-product of two vectors
!------------------------------------------------------------------------------------------------------
subroutine cross(a,b,axb)
 
    implicit none
    integer,parameter :: wp=8
    real(wp),dimension(3) :: axb
    real(wp),dimension(3),intent(in) :: a
    real(wp),dimension(3),intent(in) :: b 
 
    axb(1) = a(2)*b(3) - a(3)*b(2)
    axb(2) = a(3)*b(1) - a(1)*b(3)
    axb(3) = a(1)*b(2) - a(2)*b(1)
    
    continue
    
end subroutine cross