!-------------------------------------------------------------------------------------------------------
! Subroutine to assemble the global stiffness matrix and internal forces vector for 2D applications
!
! Input: nze - Number of the current element
!         Ke - Elemental stiffness matrix
! 
! Updates Kf - Global stiffness matrix
!-------------------------------------------------------------------------------------------------------

subroutine Assembly2D(nze,Ke,finte)
    
    use Mod_Variables
    implicit none
    
    integer(4),intent(IN)::nze
    real(8),dimension((p+1)*(q+1)*nds,(p+1)*(q+1)*nds),intent(IN)::Ke
    real(8),dimension((p+1)*(q+1)*nds,1),intent(IN)::finte
    integer(4),dimension(:,:),allocatable::penaldof
    real(8)::kmax,wp
    integer(4)::k1,k2,i,j,k,count
    
    do k1=1,(p+1)*(q+1)
        do k2=1,(p+1)*(q+1)
            Kf(conn(nze,k1)*2-1,conn(nze,k2)*2-1) = Kf(conn(nze,k1)*2-1,conn(nze,k2)*2-1) + ke(k1*2-1,k2*2-1);
            Kf(conn(nze,k1)*2-1,conn(nze,k2)*2  ) = Kf(conn(nze,k1)*2-1,conn(nze,k2)*2  ) + ke(k1*2-1,k2*2  );
            Kf(conn(nze,k1)*2  ,conn(nze,k2)*2-1) = Kf(conn(nze,k1)*2  ,conn(nze,k2)*2-1) + ke(k1*2  ,k2*2-1);
            Kf(conn(nze,k1)*2  ,conn(nze,k2)*2  ) = Kf(conn(nze,k1)*2  ,conn(nze,k2)*2  ) + ke(k1*2  ,k2*2  );
        end do   
    end do
    
    do k1=1,(p+1)*(q+1)
        fint(conn(nze,k1)*2-1,1) = fint(conn(nze,k1)*2-1,1) + finte(k1*2-1,1);
        fint(conn(nze,k1)*2  ,1) = fint(conn(nze,k1)*2  ,1) + finte(k1*2  ,1); 
    end do
    
    !Assemble penalty
    count=0
    if (npenal .gt. 0)then
        allocate(penaldof(npenal,ipenal))
        count = 0
        do k1=1,npenal
            do j=1,ncpy
                do i =1,ncpx
                    count = count + 1
                    if(penalmeth(k1,1) == i .and. penalmeth(k1,2) == j)then
                        penaldof(k1,1) = count
                    elseif(penalmeth(k1,3) == i .and. penalmeth(k1,4) == j) then
                        penaldof(k1,2) = count
                    endif
                    continue
                end do
            end do
        end do
        
        kmax = 1.0d-15
        do k1=1,ncpx*ncpy*nds
            do k2=1,ncpx*ncpy*nds
                if (abs(Kf(k1,k2)) .gt. Kmax) Kmax = Kf(k1,k2)
            end do
        end do
        
        wp = Kmax*1.0d4
        
        do i=1,npenal
            
            !Stiffness matrix penalty -----------------------------------------------------------------------------
            do j=1,ipenal
                Kf(penaldof(i,j)*nds-1,penaldof(i,j)*nds-1) = Kf(penaldof(i,j)*nds-1,penaldof(i,j)*nds-1) + wp
                Kf(penaldof(i,j)*nds  ,penaldof(i,j)*nds  ) = Kf(penaldof(i,j)*nds  ,penaldof(i,j)*nds  ) + wp
            end do 
            
            Kf(penaldof(i,1)*nds-1,penaldof(i,2)*nds-1) = Kf(penaldof(i,1)*nds-1,penaldof(i,2)*nds-1) - wp
            Kf(penaldof(i,1)*nds  ,penaldof(i,2)*nds  ) = Kf(penaldof(i,1)*nds  ,penaldof(i,2)*nds  ) - wp
            Kf(penaldof(i,2)*nds-1,penaldof(i,1)*nds-1) = Kf(penaldof(i,2)*nds-1,penaldof(i,1)*nds-1) - wp
            Kf(penaldof(i,2)*nds  ,penaldof(i,1)*nds  ) = Kf(penaldof(i,2)*nds  ,penaldof(i,1)*nds  ) - wp
            
            
            !Internal forces penalty ---------------------------------------------------------------------------------------------------
            Fint(penaldof(i,1)*nds-1,1) = Fint(penaldof(i,1)*nds-1,1) + &
             & wp*(dDisp(penaldof(i,1)*nds-1,1) - dDisp(penaldof(i,2)*nds-1,1))
            Fint(penaldof(i,1)*nds  ,1) = Fint(penaldof(i,1)*nds  ,1) + &
             & wp*(dDisp(penaldof(i,1)*nds  ,1) - dDisp(penaldof(i,2)*nds  ,1))
            Fint(penaldof(i,2)*nds-1,1) = Fint(penaldof(i,2)*nds-1,1) + &
             & wp*(dDisp(penaldof(i,2)*nds-1,1) - dDisp(penaldof(i,1)*nds-1,1))
            Fint(penaldof(i,2)*nds  ,1) = Fint(penaldof(i,2)*nds  ,1) + &
             & wp*(dDisp(penaldof(i,2)*nds  ,1) - dDisp(penaldof(i,1)*nds  ,1))
                
        end do
        
    endif
    
    continue
    !(ien(nze,(p+1)*(q+1)*nds-k1+1)-1)
    
end subroutine
