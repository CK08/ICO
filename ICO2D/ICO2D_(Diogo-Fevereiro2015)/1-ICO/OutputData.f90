!------------------------------------------------------------------------------------------------------------
!
! Subroutine to compute updated control point's coordinate for output pouposes
!
!------------------------------------------------------------------------------------------------------------

subroutine OutputMesh(nze,nzep,iptch)

    use Mod_Variables
    implicit none
    
    integer(4)::k1,k2,count,count2
    integer(4),intent(IN)::nze,nzep,iptch
    
    if(nptch .gt. 1)then
        count = 0
        count2 = (p+1)*(q+1) + 1
        do k1=1,q+1
            do k2=1,p+1
                count = count + 1
                count2 = count2 - 1
                PMesh(ien(nze,count),1) = PMesh0(ien(nze,count),1) + u(MP_conn(nzep,count)*2-1,1) + ddisp(MP_conn(nzep,count)*2-1,1)
                PMesh(ien(nze,count),2) = PMesh0(ien(nze,count),2) + u(MP_conn(nzep,count)*2  ,1) + ddisp(MP_conn(nzep,count)*2  ,1)
                continue
            end do
        end do
        
        count = 0
        do k2=1,ncpy
            do k1=1,ncpx
                count = count+1
                B_net_final(k1,k2,1:nds) = PMesh(count,:)
                B_net_final(k1,k2,nds+1) = B_net(k1,k2,nds+1)
            end do
        end do
        
        do k2=1,ncpy
            do k1=1,ncpx
                MP_b_net_final(iptch,k1,k2,:) = B_net_final(k1,k2,:)
            end do
        end do
        
        continue
        
    else    
        count = 0
        do k1=1,q+1
            do k2=1,p+1
                count = count + 1
                PMesh(ien(nze,count),1) = PMesh0(ien(nze,count),1) + u(ien(nze,count)*2-1,1) + ddisp(ien(nze,count)*2-1,1)
                PMesh(ien(nze,count),2) = PMesh0(ien(nze,count),2) + u(ien(nze,count)*2  ,1) + ddisp(ien(nze,count)*2  ,1)
                continue
            end do
        end do
        
        count = 0
        do k2=1,ncpy
            do k1=1,ncpx
                count = count+1
                B_net_final(k1,k2,1:nds) = PMesh(count,:)
                B_net_final(k1,k2,nds+1) = B_net(k1,k2,nds+1)
            end do
        end do
    end if
    
end subroutine