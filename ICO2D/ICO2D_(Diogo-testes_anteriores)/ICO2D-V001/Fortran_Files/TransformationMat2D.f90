!------------------------------------------------------------------------------------------------------
!
! Subroutine to create the transformation matrix between configurations
!
! Input: MatIn : 2x2 transformation matrix
! 
! Output: MatTrans: 3x3 transformation matrix
!
!------------------------------------------------------------------------------------------------------
    subroutine TransformationMat2D(MatIn,MatTrans)
        implicit none
        
        real(8),dimension(2,2),intent(IN)::MatIn
        real(8),dimension(3,3),intent(OUT)::MatTrans
                
        MatTrans = 0.0d0

        MatTrans=0.0d0
        MatTrans(1,1)=MatIn(1,1)*MatIn(1,1)
        MatTrans(1,2)=MatIn(1,2)*MatIn(1,2)
        MatTrans(1,3)=MatIn(1,1)*MatIn(1,2)
        MatTrans(2,1)=MatIn(2,1)*MatIn(2,1)
        MatTrans(2,2)=MatIn(2,2)*MatIn(2,2)
        MatTrans(2,3)=MatIn(2,1)*MatIn(2,2)
        MatTrans(3,1)=2.0d0*MatIn(1,1)*MatIn(2,1)
        MatTrans(3,2)=2.0d0*MatIn(1,2)*MatIn(2,2)
        MatTrans(3,3)=MatIn(1,1)*MatIn(2,2) + MatIn(2,1)*MatIn(1,2)

    end subroutine