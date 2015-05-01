!-------------------------------------------------------------------------------------------------------
! Subroutine to assemble the various matrices in the B-Bar method for 3D applications
!
! Input: nze - Number of the current element
!        Ke - Elemental stiffness matrix
!        finte - Elemental internal force vector
! 
! Updates: Kf - Global stiffness matrix
!          Fint - Global internal force vector
!-------------------------------------------------------------------------------------------------------

subroutine AssemblyBBar(nze,MA,MC,MV)
    
    use Mod_Variables
    implicit none
    
    integer(4),intent(IN)::nze
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(IN)::MA
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,p*q*w),intent(IN)::MC
    real(8),dimension(p*q*w,p*q*w),intent(IN)::MV

    integer(4)::k1,k2,k3,i,j,k,count
    
    !Assemble Matrix A ---------------------------------------------------------------------------------------------
    do k1=1,(p+1)*(q+1)*(w+1)
        do k2=1,(p+1)*(q+1)*(w+1)
            MatAf(conn(nze,k1)*3-2,conn(nze,k2)*3-2) = MatAf(conn(nze,k1)*3-2,conn(nze,k2)*3-2) + MA(k1*3-2,k2*3-2);
            MatAf(conn(nze,k1)*3-2,conn(nze,k2)*3-1) = MatAf(conn(nze,k1)*3-2,conn(nze,k2)*3-1) + MA(k1*3-2,k2*3-1);
            MatAf(conn(nze,k1)*3-2,conn(nze,k2)*3  ) = MatAf(conn(nze,k1)*3-2,conn(nze,k2)*3  ) + MA(k1*3-2,k2*3  );
            
            MatAf(conn(nze,k1)*3-1,conn(nze,k2)*3-2) = MatAf(conn(nze,k1)*3-1,conn(nze,k2)*3-2) + MA(k1*3-1,k2*3-2);
            MatAf(conn(nze,k1)*3-1,conn(nze,k2)*3-1) = MatAf(conn(nze,k1)*3-1,conn(nze,k2)*3-1) + MA(k1*3-1,k2*3-1);
            MatAf(conn(nze,k1)*3-1,conn(nze,k2)*3  ) = MatAf(conn(nze,k1)*3-1,conn(nze,k2)*3  ) + MA(k1*3-1,k2*3  );
            
            MatAf(conn(nze,k1)*3  ,conn(nze,k2)*3-2) = MatAf(conn(nze,k1)*3  ,conn(nze,k2)*3-2) + MA(k1*3  ,k2*3-2);
            MatAf(conn(nze,k1)*3  ,conn(nze,k2)*3-1) = MatAf(conn(nze,k1)*3  ,conn(nze,k2)*3-1) + MA(k1*3  ,k2*3-1);
            MatAf(conn(nze,k1)*3  ,conn(nze,k2)*3  ) = MatAf(conn(nze,k1)*3  ,conn(nze,k2)*3  ) + MA(k1*3  ,k2*3  );
        end do   
    end do
    
    !Assemble Matrix C ---------------------------------------------------------------------------------------------
    do k1=1,(p+1)*(q+1)*(w+1)
        do k2=1,p*q*w
            MatCf(conn(nze,k1)*3-2,ienb(nze,k2)) = MatCf(conn(nze,k1)*3-2,ienb(nze,k2)) + MC(k1*3-2,k2);
            MatCf(conn(nze,k1)*3-1,ienb(nze,k2)) = MatCf(conn(nze,k1)*3-1,ienb(nze,k2)) + MC(k1*3-1,k2);
            MatCf(conn(nze,k1)*3  ,ienb(nze,k2)) = MatCf(conn(nze,k1)*3  ,ienb(nze,k2)) + MC(k1*3  ,k2);
        end do   
    end do
    
    !Assemble Matrix V ---------------------------------------------------------------------------------------------
    do k1=1,p*q*w
        do k2=1,p*q*w
            MatVf(ienb(nze,k1),ienb(nze,k2)) = MatVf(ienb(nze,k1),ienb(nze,k2)) + MV(k1,k2);
        end do   
    end do
    
    continue
    
end subroutine

