!----------------------------------------------------------------------------------------------
!
! To perform the assembly of the global stiffness and internal force vector
! for a multipatch analusis
!
! Input: nze - number of the element in the current patch
!        nzep - global number of the element
!        Ke - elemental stiffness
!        finte - elemental internal force vector
!
!----------------------------------------------------------------------------------------------
subroutine AssemblyMP(nze,nzep,Ke,Finte)
    
    use Mod_Variables
    implicit none
    integer(4),intent(IN)::nze,nzep
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(IN)::Ke
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::Finte
    
    integer(4)::k1,k2
    
    do k1=1,(p+1)*(q+1)*(w+1)
        do k2=1,(p+1)*(q+1)*(w+1)
            Kf(MP_IEN(nzep,k1)*3-2,MP_IEN(nzep,k2)*3-2) = Kf(MP_IEN(nzep,k1)*3-2,MP_IEN(nzep,k2)*3-2) + ke(k1*3-2,k2*3-2);
            Kf(MP_IEN(nzep,k1)*3-2,MP_IEN(nzep,k2)*3-1) = Kf(MP_IEN(nzep,k1)*3-2,MP_IEN(nzep,k2)*3-1) + ke(k1*3-2,k2*3-1);
            Kf(MP_IEN(nzep,k1)*3-2,MP_IEN(nzep,k2)*3  ) = Kf(MP_IEN(nzep,k1)*3-2,MP_IEN(nzep,k2)*3  ) + ke(k1*3-2,k2*3  );
            
            Kf(MP_IEN(nzep,k1)*3-1,MP_IEN(nzep,k2)*3-2) = Kf(MP_IEN(nzep,k1)*3-1,MP_IEN(nzep,k2)*3-2) + ke(k1*3-1,k2*3-2);
            Kf(MP_IEN(nzep,k1)*3-1,MP_IEN(nzep,k2)*3-1) = Kf(MP_IEN(nzep,k1)*3-1,MP_IEN(nzep,k2)*3-1) + ke(k1*3-1,k2*3-1);
            Kf(MP_IEN(nzep,k1)*3-1,MP_IEN(nzep,k2)*3  ) = Kf(MP_IEN(nzep,k1)*3-1,MP_IEN(nzep,k2)*3  ) + ke(k1*3-1,k2*3  );
            
            Kf(MP_IEN(nzep,k1)*3  ,MP_IEN(nzep,k2)*3-2) = Kf(MP_IEN(nzep,k1)*3  ,MP_IEN(nzep,k2)*3-2) + ke(k1*3  ,k2*3-2);
            Kf(MP_IEN(nzep,k1)*3  ,MP_IEN(nzep,k2)*3-1) = Kf(MP_IEN(nzep,k1)*3  ,MP_IEN(nzep,k2)*3-1) + ke(k1*3  ,k2*3-1);
            Kf(MP_IEN(nzep,k1)*3  ,MP_IEN(nzep,k2)*3  ) = Kf(MP_IEN(nzep,k1)*3  ,MP_IEN(nzep,k2)*3  ) + ke(k1*3  ,k2*3  );
        end do   
    end do
    
    do k1=1,(p+1)*(q+1)*(w+1)
        fint(MP_IEN(nzep,k1)*3-2,1) = fint(MP_IEN(nzep,k1)*3-2,1) + finte(k1*3-2,1);
        fint(MP_IEN(nzep,k1)*3-1,1) = fint(MP_IEN(nzep,k1)*3-1,1) + finte(k1*3-1,1);
        fint(MP_IEN(nzep,k1)*3  ,1) = fint(MP_IEN(nzep,k1)*3  ,1) + finte(k1*3  ,1);
    end do
    
    continue

end subroutine