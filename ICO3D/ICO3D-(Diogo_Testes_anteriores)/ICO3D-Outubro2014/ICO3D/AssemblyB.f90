!-------------------------------------------------------------------------------------------------------
! Subroutine to assemble the global stiffness matrix and internal forces vector for 2D applications
!
! Input: nze - Number of the current element
!         Ke - Elemental stiffness matrix
! 
! Updates Kf - Global stiffness matrix
!-------------------------------------------------------------------------------------------------------

subroutine AssemblyB(nze,Ke,finte)
    
    use Mod_Variables
    !implicit none
    
    integer(4),intent(IN)::nze
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(IN)::Ke
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::finte
    integer(4),dimension(:,:),allocatable::penaldof
    real(8)::kmax,wp
    integer(4)::k1,k2,k3,i,j,k,count
    
    do k1=1,(p+1)*(q+1)*(w+1)

        do k2=1,6
            BBf(k2,conn(nze,k1)*3-2) = BBf(k2,conn(nze,k1)*3-2)+ BBare(k2,k1*3-2)
            BBf(k2,conn(nze,k1)*3-1) = BBf(k2,conn(nze,k1)*3-1)+ BBare(k2,k1*3-1)
            BBf(k2,conn(nze,k1)*3  ) = BBf(k2,conn(nze,k1)*3  )+ BBare(k2,k1*3  )
            
            BDf(k2,conn(nze,k1)*3-2) = BDf(k2,conn(nze,k1)*3-2)+ BDeve(k2,k1*3-2)
            BDf(k2,conn(nze,k1)*3-1) = BDf(k2,conn(nze,k1)*3-1)+ BDeve(k2,k1*3-1)
            BDf(k2,conn(nze,k1)*3  ) = BDf(k2,conn(nze,k1)*3  )+ BDeve(k2,k1*3  )
        
!            Kf(conn(nze,k1)*3-2,conn(nze,k2)*3-2) = Kf(conn(nze,k1)*3-2,conn(nze,k2)*3-2) + ke(k1*3-2,k2*3-2);
!            Kf(conn(nze,k1)*3-2,conn(nze,k2)*3-1) = Kf(conn(nze,k1)*3-2,conn(nze,k2)*3-1) + ke(k1*3-2,k2*3-1);
!            Kf(conn(nze,k1)*3-2,conn(nze,k2)*3  ) = Kf(conn(nze,k1)*3-2,conn(nze,k2)*3  ) + ke(k1*3-2,k2*3  );
!            
!            Kf(conn(nze,k1)*3-1,conn(nze,k2)*3-2) = Kf(conn(nze,k1)*3-1,conn(nze,k2)*3-2) + ke(k1*3-1,k2*3-2);
!            Kf(conn(nze,k1)*3-1,conn(nze,k2)*3-1) = Kf(conn(nze,k1)*3-1,conn(nze,k2)*3-1) + ke(k1*3-1,k2*3-1);
!            Kf(conn(nze,k1)*3-1,conn(nze,k2)*3  ) = Kf(conn(nze,k1)*3-1,conn(nze,k2)*3  ) + ke(k1*3-1,k2*3  );
!            
!            Kf(conn(nze,k1)*3  ,conn(nze,k2)*3-2) = Kf(conn(nze,k1)*3  ,conn(nze,k2)*3-2) + ke(k1*3  ,k2*3-2);
!            Kf(conn(nze,k1)*3  ,conn(nze,k2)*3-1) = Kf(conn(nze,k1)*3  ,conn(nze,k2)*3-1) + ke(k1*3  ,k2*3-1);
!            Kf(conn(nze,k1)*3  ,conn(nze,k2)*3  ) = Kf(conn(nze,k1)*3  ,conn(nze,k2)*3  ) + ke(k1*3  ,k2*3  );
!            
!!            Kf(conn(nze,k1)*2-1,conn(nze,k2)*2-1) = Kf(conn(nze,k1)*2-1,conn(nze,k2)*2-1) + ke(k1*2-1,k2*2-1);
!!            Kf(conn(nze,k1)*2-1,conn(nze,k2)*2  ) = Kf(conn(nze,k1)*2-1,conn(nze,k2)*2  ) + ke(k1*2-1,k2*2  );
!!            Kf(conn(nze,k1)*2  ,conn(nze,k2)*2-1) = Kf(conn(nze,k1)*2  ,conn(nze,k2)*2-1) + ke(k1*2  ,k2*2-1);
!!            Kf(conn(nze,k1)*2  ,conn(nze,k2)*2  ) = Kf(conn(nze,k1)*2  ,conn(nze,k2)*2  ) + ke(k1*2  ,k2*2  );
        end do
    end do
    
    continue
 
    
end subroutine
