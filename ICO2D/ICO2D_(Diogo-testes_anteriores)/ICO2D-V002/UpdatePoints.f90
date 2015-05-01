!----------------------------------------------------------------------------------------------
!
! Subroutine to update the control points coordinates for geometric nonlinearity
!
!----------------------------------------------------------------------------------------------

subroutine UpdatePoints(nze,nzep)

    use Mod_Variables
    implicit none
    
    integer(4)::k1,k2,k3,loc_num
    integer(4),intent(IN)::nze,nzep
    
    
    if(nptch .gt. 1)then
        loc_num = 0
        do k1=1,(p+1)
            do k2=1,(q+1)
                loc_num = loc_num + 1
                do k3=1,nds
                    Points(IEN(nze,loc_num),k3) = Points0(IEN(nze,loc_num),k3) + u(MP_conn(nzep,loc_num)*nds-(nds-k3),1)
                end do
            end do   
        end do
    else
        loc_num = 0
        do k1=1,(p+1)
            do k2=1,(q+1)
                loc_num = loc_num + 1
                do k3=1,nds
                    Points(IEN(nze,loc_num),k3) = Points0(IEN(nze,loc_num),k3) + u(ien(nze,loc_num)*nds-(nds-k3),1)
                end do
            end do   
        end do
    end if
    
end subroutine