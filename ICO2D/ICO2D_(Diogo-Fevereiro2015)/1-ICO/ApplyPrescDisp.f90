!----------------------------------------------------------------------------------------------
!
! Subroutine to apply prescribed displacement boundary conditions
!
!----------------------------------------------------------------------------------------------

subroutine ApplyPrescDisp

    use Mod_Variables
    implicit none
    
    integer(4)::i,k1
    
    if(nptch == 1)then
        do i=1,nnodes*nds
            if(ImpDisp(i,1)/=0.0d0)then
                do k1=1,nnodes*nds
                    if(i /= k1) then
                        Fext(k1,1)=Fext(k1,1) - Kf(k1,i)*ImpDisp(i,1)/(1.0d0*incmax)*(1.0d0-stpfrac(i))
                        Kf(i,k1)=0.0d0
                        Kf(k1,i)=0.0d0
                     else 
                        Fext(i,1)=Kf(i,i)*ImpDisp(i,1)/(1.0d0*incmax)*(1.0d0-stpfrac(i))  
                        continue
                    end if
                end do
            end if
        end do
    else
        do i=1,tnodes*nds
            if(ImpDisp(i,1)/=0.0d0)then
                do k1=1,tnodes*nds
                    if(i /= k1) then
                        Fext(k1,1)=Fext(k1,1) - Kf(k1,i)*(ImpDisp(i,1)/(1.0d0*incmax)*(1.0d0-stpfrac(i)))
                        Kf(i,k1)=0.0d0
                        Kf(k1,i)=0.0d0
                     else 
                        Fext(i,1)=Kf(i,i)*(ImpDisp(i,1)/(1.0d0*incmax)*(1.0d0-stpfrac(i)))
                        continue
                    end if
                end do
            end if
        end do
    end if
    
    
end subroutine