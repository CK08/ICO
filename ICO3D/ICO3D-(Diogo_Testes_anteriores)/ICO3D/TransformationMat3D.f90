!------------------------------------------------------------------------------------------------------
!
! Subroutine to create the transformation matrix between configurations
!
! Input: MatIn : 3x3 transformation matrix
! 
! Output: MatTrans: 6x6 transformation matrix
!
!------------------------------------------------------------------------------------------------------
    subroutine TransformationMat3D(MatIn,MatTrans)
        implicit none
        
        integer(4)::k1,k2
        real(8),dimension(3,3),intent(IN)::MatIn
        real(8),dimension(6,6),intent(OUT)::MatTrans
                
        do k1=1,6
            do k2=1,6
                MatTrans(k1,k2)=0.0d0
            end do
        end do
    
        MatTrans(1,1)=MatIn(1,1)*MatIn(1,1)
        MatTrans(1,2)=MatIn(1,2)*MatIn(1,2)
        MatTrans(1,3)=MatIn(1,3)*MatIn(1,3)
        MatTrans(1,4)=MatIn(1,1)*MatIn(1,2)
        MatTrans(1,5)=MatIn(1,1)*MatIn(1,3)
        MatTrans(1,6)=MatIn(1,2)*MatIn(1,3)
                
        MatTrans(2,1)=MatIn(2,1)*MatIn(2,1)
        MatTrans(2,2)=MatIn(2,2)*MatIn(2,2)
        MatTrans(2,3)=MatIn(2,3)*MatIn(2,3)
        MatTrans(2,4)=MatIn(2,1)*MatIn(2,2)
        MatTrans(2,5)=MatIn(2,1)*MatIn(2,3)
        MatTrans(2,6)=MatIn(2,2)*MatIn(2,3)
                
        MatTrans(3,1)=MatIn(3,1)*MatIn(3,1)
        MatTrans(3,2)=MatIn(3,2)*MatIn(3,2)
        MatTrans(3,3)=MatIn(3,3)*MatIn(3,3)
        MatTrans(3,4)=MatIn(3,1)*MatIn(3,2)
        MatTrans(3,5)=MatIn(3,1)*MatIn(3,3)
        MatTrans(3,6)=MatIn(3,2)*MatIn(3,3)
               
        MatTrans(4,1)=2.0d0*MatIn(1,1)*MatIn(2,1)
        MatTrans(4,2)=2.0d0*MatIn(1,2)*MatIn(2,2)
        MatTrans(4,3)=2.0d0*MatIn(1,3)*MatIn(2,3)
        MatTrans(4,4)=MatIn(1,1)*MatIn(2,2) + MatIn(1,2)*MatIn(2,1)
        MatTrans(4,5)=MatIn(1,1)*MatIn(2,3) + MatIn(2,1)*MatIn(1,3)
        MatTrans(4,6)=MatIn(1,2)*MatIn(2,3) + MatIn(2,2)*MatIn(1,3)
                
        MatTrans(5,1)=2.0d0*MatIn(1,1)*MatIn(3,1)
        MatTrans(5,2)=2.0d0*MatIn(1,2)*MatIn(3,2)
        MatTrans(5,3)=2.0d0*MatIn(1,3)*MatIn(3,3)
        MatTrans(5,4)=MatIn(1,1)*MatIn(3,2) + MatIn(1,2)*MatIn(3,1)
        MatTrans(5,5)=MatIn(1,1)*MatIn(3,3) + MatIn(3,1)*MatIn(1,3)
        MatTrans(5,6)=MatIn(1,2)*MatIn(3,3) + MatIn(3,2)*MatIn(1,3)
                
        MatTrans(6,1)=2.0d0*MatIn(2,1)*MatIn(3,1)
        MatTrans(6,2)=2.0d0*MatIn(2,2)*MatIn(3,2)
        MatTrans(6,3)=2.0d0*MatIn(2,3)*MatIn(3,3)
        MatTrans(6,4)=MatIn(2,1)*MatIn(3,2) + MatIn(2,2)*MatIn(3,1)
        MatTrans(6,5)=MatIn(2,1)*MatIn(3,3) + MatIn(3,1)*MatIn(2,3)
        MatTrans(6,6)=MatIn(2,2)*MatIn(3,3) + MatIn(3,2)*MatIn(2,3)

    end subroutine TransformationMat3D
    
    
 