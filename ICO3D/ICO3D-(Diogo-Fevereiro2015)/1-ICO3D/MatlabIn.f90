!---------------------------------------------------------------------------------------
!
! Subroutine to writte data of the undeformed geometry for the matlab plot subroutine 
!
!---------------------------------------------------------------------------------------

subroutine MatlabIn
    
    use Mod_Variables
    implicit none
    integer(4)::i,j,k,l
    
    if(nptch .gt. 1) then
        open(unit=7,file='NURBSInput.txt')
        
        write(7,FMT=100)nptch
        
        do i=1,nptch
            write(7,FMT=101)MP_ncpx(i), MP_ncpy(i), MP_ncpz(i)
        end do
        
        do i=1,nptch
            write(7,FMT=101)MP_p(i), MP_q(i), MP_w(i)
        end do
        
        do i=1,nptch
            
            write(7,FMT=102)MP_uknot(i,1:MP_ncpx(i)+MP_p(i)+1)
            write(7,FMT=102)MP_vknot(i,1:MP_ncpy(i)+MP_q(i)+1)
            write(7,FMT=102)MP_wknot(i,1:MP_ncpz(i)+MP_w(i)+1)
            
            do j=1,MP_ncpx(i)
                do k=1,MP_ncpy(i)
                    do l=1,MP_ncpz(i)
                        write(7,FMT=102) MP_b_net(i,j,k,l,:)
                    end do
                end do
            end do
        end do
        
        100 format(I4,',')
        101 format(3(I4,','))
        102 format(100(E,','))
        
        close(unit=7)
    else
        
        open(unit=7,file='NURBSInput.txt')
        
        write(7,FMT=200)nptch
        
        write(7,FMT=201)ncpx, ncpy, ncpz
        
        write(7,FMT=201)p, q, w
        
        write(7,FMT=202)u_knot(:)
        write(7,FMT=202)v_knot(:)
        write(7,FMT=202)w_knot(:)
            
        do j=1,ncpx
            do k=1,ncpy
                do l=1,ncpz
                    write(7,FMT=202) b_net(j,k,l,:)
                end do
            end do
        end do
        
        200 format(I4,',')
        201 format(3(I4,','))
        202 format(100(E,','))
        
        close(unit=7)
        
        continue
    
    end if
    
    continue

end subroutine