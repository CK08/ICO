!----------------------------------------------------------------------------------------------
!
! Subroutine to eliminate lines and columns of the global stiffness matrix due to  fixed BC's
!
!----------------------------------------------------------------------------------------------

subroutine ApplyBCs()

    use Mod_Variables
    implicit none
    
    integer(4)::i,j,k,count
    
    if(nptch == 1)then
        count = 0
        do j=1,ncpy
            do i=1,ncpx
                do k=1,nds
                    count = count+1
                    if(BC(i,j,k)==0)then
                        Kf(count,:) = 0.0d0
                        Kf(:,count) = 0.0d0
                        !Fint(count,1) = 0.0d0
                        redcount(count) = 0
                        continue
                    endif    
                end do
            end do
        end do
        
        count = 0
        do i=1,tnodes*nds
            count = count + 1
            if(BCdof(i) == 0) then
                Kf(i,:) = 0.0d0
                Kf(:,i) = 0.0d0
                redcount(count) = 0
            end if
        end do
        
    else
        count = 0
        do i=1,tnodes*nds
            count = count + 1
            if(BCdof(i) == 0) then
                Kf(i,:) = 0.0d0
                Kf(:,i) = 0.0d0
                redcount(count) = 0
            end if
        end do
    end if

end subroutine