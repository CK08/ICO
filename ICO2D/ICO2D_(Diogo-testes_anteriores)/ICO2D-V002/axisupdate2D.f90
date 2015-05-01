!----------------------------------------------------------------------------------------------
!
! Axis Update subroutine
!
! Input: raux - Axis Matrix
!        matR - Rotation Matrix
!
! Output: rconv - Updated axis
!
!----------------------------------------------------------------------------------------------
subroutine AxisUpdate2D(raux,matr,rconv)
        implicit none
        
        real(8),dimension(2,2),intent(IN)::raux
        real(8),dimension(3,3),intent(IN)::matr
        real(8),dimension(2,2),intent(OUT)::rconv
        integer(4)::j
        real(8),dimension(3,1)::aux31,auxx,auxy,auxz
        
        do j=1,2
            aux31(j,1)=raux(j,1)
        end do
        aux31(3,1) = 0.0d0
        
        auxx=matmul(matR,aux31)
        
        do j=1,2
            aux31(j,1)=raux(j,2)
        end do
        aux31(3,1 )= 0.0d0
                
        auxy=matmul(matR,aux31)
        
        call cross(auxx,auxy,auxz)
                
        do j=1,2
            rconv(j,1)=auxx(j,1)
            rconv(j,2)=auxy(j,1)
        end do


    end subroutine