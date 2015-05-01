!---------------------------------------------------------------------------------------
!
! Subroutine to calculate the Hardening slope based on the inserted data points
!
! Input: nprops - number of properties
!        props - material properties
!        iprops - position where data starts
!        Syield - Current yield stress
!
! Output: H - Hardening parameter to be determined

subroutine MatIsoHardening(nprops,props,iprops,SYield,H)
    
    
    implicit none
    
    integer(4),intent(IN)::nprops,iprops
    real(8),dimension(nprops),intent(IN)::props
    real(8),intent(IN)::Syield
    real(8),intent(OUT)::H
    
    integer(4)::i
    real(8)::YieldUB,YieldLB,StrainLB,StrainUB
    
    H=0.0d0
    
    do i=iprops, nprops, 2
        if(syield >= props(i) .and. i .lt. nprops-1)then
            
            if(syield >= props(i+2)) goto 2
            
            !Lower Data Point ------
            YieldLB = props(i)
            StrainLB = props(i+1)
            
            !Higher Data Point -----
            YieldUB = props(i+2)
            StrainUB = props(i+3)
            
            H = (YieldUB - YieldLB)/(StrainUB - StrainLB)
            
            goto 1
            
        elseif(syield .gt. props(i) .and. i == nprops-1) then
            write(*,*)'**warning** - Yield Stress higher than inserted data'
            write(*,*)'              Please, insert additional data points.'
        
        elseif(syield .lt. props(i)) then
            write(*,*)'**warning** - Yielding not initiated: should not be using this routine'
            write(*,*)'              Please, check code!'
        
        end if
        
        2 continue
        
    end do
    write(*,*)'**warning** - End of hardening cycle without calculating H'
    write(*,*)'              Please check code or input data!'
    
    1 continue 

end subroutine MatIsoHardening
