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
subroutine axisupdate(raux,matr,rconv)
        implicit none
        
        real(8),dimension(3,3),intent(IN)::raux,matr
        real(8),dimension(3,3),intent(OUT)::rconv
        integer(4)::j
        real(8),dimension(3,1)::aux31,auxx,auxy,auxz
        
        do j=1,3
            aux31(j,1)=raux(j,1)
        end do
        
        auxx=matmul(matR,aux31)
        
        do j=1,3
            aux31(j,1)=raux(j,2)
        end do
        
        auxy=matmul(matR,aux31)
        
        call cross(auxx,auxy,auxz)
                
        do j=1,3
            rconv(j,1)=auxx(j,1)
            rconv(j,2)=auxy(j,1)
            rconv(j,3)=auxz(j,1)
        end do


    end subroutine axisupdate