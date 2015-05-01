!-------------------------------------------------------------------------------------------------------
! Subroutine to assemble the global stiffness matrix and internal forces vector for multipatch
! 2D applications
!
! Input: nze - Number of the current element
!         Ke - Elemental stiffness matrix
! 
! Updates Kf - Global stiffness matrix
!-------------------------------------------------------------------------------------------------------

subroutine AssemblyPatch(nel,nze,Ke,finte)
    
    use Mod_Variables
    implicit none
    
    integer(4),intent(IN)::nel,nze
    real(8),dimension((p+1)*(q+1)*nds,(p+1)*(q+1)*nds),intent(IN)::Ke
    real(8),dimension((p+1)*(q+1)*nds,1),intent(IN)::finte
    real(8)::kmax,wp
    integer(4)::k1,k2,k3,i,j,k,count,ni,nj,loc_num
    
    do k1=1,(p+1)*(q+1)
        do k2=1,(p+1)*(q+1)
            Kf(MP_Conn(nze,k1)*2-1,MP_Conn(nze,k2)*2-1) = Kf(MP_Conn(nze,k1)*2-1,MP_Conn(nze,k2)*2-1) + ke(k1*2-1,k2*2-1);
            Kf(MP_Conn(nze,k1)*2-1,MP_Conn(nze,k2)*2  ) = Kf(MP_Conn(nze,k1)*2-1,MP_Conn(nze,k2)*2  ) + ke(k1*2-1,k2*2  );
            Kf(MP_Conn(nze,k1)*2  ,MP_Conn(nze,k2)*2-1) = Kf(MP_Conn(nze,k1)*2  ,MP_Conn(nze,k2)*2-1) + ke(k1*2  ,k2*2-1);
            Kf(MP_Conn(nze,k1)*2  ,MP_Conn(nze,k2)*2  ) = Kf(MP_Conn(nze,k1)*2  ,MP_Conn(nze,k2)*2  ) + ke(k1*2  ,k2*2  );
        end do   
    end do
    
    do k1=1,(p+1)*(q+1)
        fint(MP_Conn(nze,k1)*2-1,1) = fint(MP_Conn(nze,k1)*2-1,1) + finte(k1*2-1,1);
        fint(MP_Conn(nze,k1)*2  ,1) = fint(MP_Conn(nze,k1)*2  ,1) + finte(k1*2  ,1); 
    end do
    
!    ni = INN(IEN(nel,1),1)
!    nj = INN(IEN(nel,1),2)
!    
    loc_num = 0
    do k1=1,(p+1)
        do k2=1,(q+1)
            loc_num = loc_num + 1
            do k3=1,nds
                GCoords(MP_Conn(nze,loc_num),k3) = Points(IEN(nel,loc_num),k3) + u(MP_conn(nze,loc_num)*nds-(nds-k3),1)
            end do
            
            GCoords(MP_Conn(nze,loc_num),3) = Weights(IEN(nel,loc_num))

        end do   
    end do
    
    
    continue
    
end subroutine