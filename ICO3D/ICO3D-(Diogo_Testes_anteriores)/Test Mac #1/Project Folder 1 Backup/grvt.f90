!-----------------------------------------------------------------------------------------
!
! Subroutine to calculate gravitic load in the integration points
!
! NOTE: The direction of the load is considered in the opposite direction of the gravity
!
! Input: 
!       ishp: number of shape functions
!       nds: problem dimension
!       R: shape function array for th integration point
!       gwt: weight factor for the integration point
!       gconst: gravitational constant
!       gravdir: direction of the gravity
!       density: material density
!
! Output:
!       Finte: Updated internal force vector
!
!-----------------------------------------------------------------------------------------

subroutine grvt(ishp,nds,R,gwt,gconst,gravdir,density,Finte)
    
    implicit none
    
    integer(4),intent(IN)::ishp,nds
    real(8),dimension(ishp),intent(IN)::R
    real(8),intent(IN)::gwt,gconst,density
    integer(4),intent(IN)::gravdir
    
    real(8),dimension(ishp*nds,1),intent(INOUT)::Finte
    
    integer(4)::i,j,k
    real(8),dimension(nds,ishp*nds)::MatShp
    real(8),dimension(nds,1)::vec
    real(8),dimension(ishp*nds,1)::temp
    
    !Shape Function Matrix ---------
    MatShp = 0.0d0
    do i=1,ishp
        MatShp(1,i*3-2) = R(i)
        MatShp(2,i*3-1) = R(i)
        MatShp(3,i*3  ) = R(i)
    end do
    
    !Direction and magnitude of the gravitic load
    vec = 0.0d0
    vec(gravdir,1) = density*gwt*gconst*(-1.0d0)
    
    temp = 0.0d0
    temp = matmul(transpose(MatShp),vec)
    
    !Update Internal Forces ----
    do i=1,ishp*nds
        Finte(i,1) = Finte(i,1) - temp(i,1)
    end do
    
    continue
    

end subroutine
