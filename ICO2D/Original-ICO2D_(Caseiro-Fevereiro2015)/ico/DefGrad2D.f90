!----------------------------------------------------------------------------------------------
!
! Subroutine to calculate the deformation gradient matF of 2D problems
!
! Input: nnodes - number of nodes of the element
!        ndof - problem dimension
!        dNdc - shape function derivatives wrt xx
!        dNdc - shape function derivatives wrt yy
!        dNdc - shape function derivatives wrt zz
!        jacinv - inverse of the jacobian matrix (not used)
!        coords - nodal coordinates
!        disp - nodal displacement
!
! Output: MatF - deformation gradient
!
!----------------------------------------------------------------------------------------------
  subroutine DefGrad2D(nnodes,ndof,dNdc,dNde,jacinv,coords,disp,MatF)
    implicit none
    integer*4,intent(IN)::nnodes,ndof
    real*8,dimension(nnodes),intent(IN)::dNdc,dNde
    real*8,dimension(ndof,ndof),intent(IN)::jacinv
    real*8,dimension(nnodes,ndof),intent(IN)::coords
    real*8,dimension(ndof*nnodes,1),intent(IN)::disp
    real*8,dimension(3,3),intent(OUT)::MatF
    
    integer*4::j,k,l
    
    real*8,dimension(ndof,1)::dNdce,dNdxy
    
    
    matF=0.0d0
    do j=1,nnodes
        dNdce=0.0d0
        dNdxy=0.0d0
        dNdce(1,1)=dNdc(j)
        dNdce(2,1)=dNde(j)
        dNdxy=dNdce
        
        do k=1,ndof
            do l=1,ndof
                matF(k,l) = matF(k,l) + dNdxy(l,1)*(coords(j,k) + disp(j*2-(2-k),1))
            end do
        end do
        
    end do
    
    MatF(3,3) = 1.0d0

  end subroutine