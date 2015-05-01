!------------------------------------------------------------------------------------------------------------
!
! Subroutine to calculate the equivalent von Mises stress for plane stress and 3D cases
!
! Input: stress - stress array
!        sdim - number of stess components
! Output: SVM - Equivalent von Mises stress
!
!------------------------------------------------------------------------------------------------------------
subroutine vonMisesStress(stress,sdim,SVM)
    implicit none
    
    integer(4),intent(IN)::sdim
    real(8),dimension(sdim,1)::stress
    real(8),intent(OUT)::SVM
    real(8)::aux
    
    SVM=0.0d0
    
    if(sdim==3)then
        SVM=sqrt(stress(1,1)**2.0d0+stress(2,1)**2.0d0-stress(1,1)*stress(2,1)+3.0d0*stress(3,1)**2.0d0)
    elseif(sdim==6)then
        SVM=(stress(1,1)-stress(2,1))**2+(stress(2,1)-stress(3,1))**2+(stress(1,1)-stress(3,1))**2
        SVM=SVM+6.0d0*(stress(4,1)**2+stress(5,1)**2+stress(6,1)**2)
        SVM=SVM/2.0d0
        SVM=sqrt(SVM)
    end if

    continue

end subroutine