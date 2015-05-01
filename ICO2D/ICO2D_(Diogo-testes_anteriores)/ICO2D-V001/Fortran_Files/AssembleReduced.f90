!----------------------------------------------------------------------------------------------
!
! Subroutine to assemble reduced stiffness matrix and load vector
!
!----------------------------------------------------------------------------------------------

subroutine AssembleReduced()

    use Mod_Variables
    implicit none
    
    integer(4)::i,j,k1,k2
    
    Kfeq = 0.0d0
    FextEq = 0.0d0
    if(nptch==1)then
        k1 = 1
        k2 = 1
        do i=1,nnodes*nds 
            do j=1,nnodes*nds
                if (redcount(i)/= 0 .and. redcount(j)/= 0)then
                    Kfeq(k1,k2)=Kf(i,j)
                    FextEq(k1,1)=Fext(i,1)                
                    k2=k2+1
                    if(k2==nnodes*nds-nbc+1)then
                        k1=k1+1
                        k2=1
                    end if
                end if
            end do
        end do
    
    else
        k1 = 1
        k2 = 1
        do i=1,tnodes*nds 
            do j=1,tnodes*nds
                if (redcount(i)/= 0 .and. redcount(j)/= 0)then
                    Kfeq(k1,k2)=Kf(i,j)
                    FextEq(k1,1)=Fext(i,1)                
                    k2=k2+1
                    if(k2==tnodes*nds-nbc+1)then
                        k1=k1+1
                        k2=1
                    end if
                end if
            end do
        end do
        
    end if

end subroutine