!------------------------------------------------------------------------------------
!
! Polar Decomposition algorithm based on a 2D simplification
!
! Input: ndof - problem  dimension
!        matF - deformation gradient matrix
!
! Output: MatR - Rotation Matrix
! 
!------------------------------------------------------------------------------------
subroutine PolarDecomp2D(ndof,matF,matR)
    
    integer(4),intent(IN)::ndof
    real(8),dimension(ndof,ndof),intent(IN)::matF
    real(8),dimension(ndof,ndof),intent(OUT)::matR
    
    real(8),dimension(ndof,ndof)::umatr,uinve,rmatr,cmatr
    real(8),dimension(ndof,ndof)::verify1,verify2
    real(8),dimension(ndof)::temp1
    
    real(8)::in_c1,in_c2,in_u1,in_u2
    real(8)::Ident
    integer(4)::ilixo
    
    cmatr=matmul(transpose(matF),matF)
    
    !2d approximation for rotation operator
    in_c1 = cmatr(1,1) + cmatr(2,2)
    in_c2 = cmatr(1,1)*cmatr(2,2) - cmatr(1,2)*cmatr(2,1)
    in_u1 = dsqrt(in_c1 + 2.d0*dsqrt(in_c2))
    in_u2 = dsqrt(in_c2)
          
    do i = 1,3
        do j = 1,3
            if (i.eq.j) ident = 1.d0
            if (i.ne.j) ident = 0.d0
            umatr(i,j) = (1.d0/in_u1)*(in_u2*ident + cmatr(i,j))
        enddo
    enddo
    
    do i = 1,3
        do j = 1,3
            uinve(i,j) = umatr(i,j)
        enddo
    enddo
    
    do i = 1,3
        temp1(i) = 0.d0
    enddo
    
    ilixo = 0
    call gaussj (uinve,3,temp1,ilixo)
    
    rmatr=matmul(matf,uinve)
    
    matR=rmatr
    
    verify1=matmul(matR,transpose(matR))
    verify2=matmul(transpose(matR),matR)
     
    continue
          
    end subroutine PolarDecomp2D