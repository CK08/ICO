!----------------------------------------------------------------------------------------------
!
! Subroutine to eliminate internal loads at fixed BC's and nodes with prescribed displacements
!
!----------------------------------------------------------------------------------------------

subroutine EliminateReac(iLM)

    use Mod_variables
    implicit none
    
    integer(4)::i,j,k,count
    integer(4),intent(IN)::iLM
    
    if (nptch==1)then
        count = 0
        do j=1,ncpy
            do i=1,ncpx
                do k=1,nds
                    count = count+1
                    if(BC(i,j,k)==0)then
                        Fint(count,1) = 0.0d0
                        continue
                    endif    
                end do
            end do
        end do
        
        count = 0
        do i=1,tnodes*nds
            count = count + 1
            if(BCdof(i) == 0) then
                Fint(i,1) = 0.0d0
                if(iLM .gt. 0) FintLM(i,1) = 0.0d0
            end if
        end do
        
    else
        count = 0
        do i=1,tnodes*nds
            count = count + 1
            if(BCdof(i) == 0) then
                Fint(i,1) = 0.0d0
                if(iLM .gt. 0) FintLM(i,1) = 0.0d0
            end if
        end do
    end if
    
    if(nptch==1)then
        do i=1,nnodes*nds
            if(ImpDisp(i,1) /= 0.0d0)then
                Fint(i,1)=0.0d0
            end if 
        end do
    else
        do i=1,tnodes*nds
            if(ImpDisp(i,1) /= 0.0d0)then
                Fint(i,1)=0.0d0
                if(iLM .gt. 0) FintLM(i,1)=0.0d0
            end if 
        end do
    end if

end subroutine