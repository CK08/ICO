!----------------------------------------------------------------------------------------------
!
! Subroutine to update the control lattice
!
!----------------------------------------------------------------------------------------------

subroutine Updt_bnet(iptch,nzep)
    
    use Mod_Variables
    implicit none
    integer(4),intent(IN)::iptch,nzep
    integer(4)::nnod,i,j,k,l,count,nel,k1
    real(8),dimension(:,:),allocatable::PT,PT0
    
    nel = nzep
    
    !Allocate order of patch
    p = MP_p(iptch)
    q = MP_q(iptch)
    w = MP_q(iptch)
    
    !Allocate number of control points
    ncpx = MP_ncpx(iptch)
    ncpy = MP_ncpy(iptch)
    ncpz = MP_ncpz(iptch)
    
    nnodes = ncpx*ncpy*ncpz
    nelems = (ncpx-p)*(ncpy-q)*(ncpz-w)
    nshpl = (p+1)*(q+1)*(q+1)
    
    !call gen_IEN_INN()
    
    allocate(PT(ncpx*ncpy*ncpz,nds))
    PT = 0.0d0
    
    allocate(PT0(ncpx*ncpy*ncpz,nds))
    PT0 = 0.0d0
    
    count = 0
    do k = 1,ncpz
        do j = 1,ncpy
            do i = 1,ncpx
                count = count + 1
                PT(count,1) = b_net(i,j,k,1)
                PT(count,2) = b_net(i,j,k,2)
                PT(count,3) = b_net(i,j,k,3)
            end do
        end do
    end do
    
    PT0 = PT
    
    do i=1, ncpx
        do j=1,ncpy
            do k=1,ncpz
                if((u_knot(i) .ne. u_knot(i+1)) .and. (v_knot(j) .ne. v_knot(j+1)) .and. (w_knot(k) .ne. w_knot(k+1))) then
                    
                    nel = nel + 1
                    count = 0
                    
                    do k1=1,(p+1)*(q+1)*(w+1)
                        count = count + 1
                        if(nlgeom==.false.)then
                            PT(IEN(nel-nzep,k1),1) = PT0(IEN(nel-nzep,k1),1) + u(MP_IEN(nel,count)*3-2,1) + &
                                    & ddisp(MP_IEN(nel,count)*3-2,1)
                            PT(IEN(nel-nzep,k1),2) = PT0(IEN(nel-nzep,k1),2) + u(MP_IEN(nel,count)*3-1,1) + &
                                    & ddisp(MP_IEN(nel,count)*3-1,1)
                            PT(IEN(nel-nzep,k1),3) = PT0(IEN(nel-nzep,k1),3) + u(MP_IEN(nel,count)*3  ,1) + &
                                    & ddisp(MP_IEN(nel,count)*3  ,1)
                        else
                            PT(IEN(nel-nzep,k1),1) = PT0(IEN(nel-nzep,k1),1) + ddisp(MP_IEN(nel,count)*3-2,1)
                            PT(IEN(nel-nzep,k1),2) = PT0(IEN(nel-nzep,k1),2) + ddisp(MP_IEN(nel,count)*3-1,1)
                            PT(IEN(nel-nzep,k1),3) = PT0(IEN(nel-nzep,k1),3) + ddisp(MP_IEN(nel,count)*3  ,1)
                        end if
                    end do
                    
                    continue
                    
                end if        
            end do
        end do
    end do    
    
    count = 0
    b_net_final = 0.0d0
    do k=1,ncpz
        do j = 1,ncpy
            do i = 1,ncpx
                count = count + 1
                b_net_final(i,j,k,1) = PT(count,1) 
                b_net_final(i,j,k,2) = PT(count,2)
                b_net_final(i,j,k,3) = PT(count,3)
                b_net_final(i,j,k,4) = Weight(count)
            end do
        end do
    end do
    

    do i=1,ncpx
        do j=1,ncpy
            do k=1,ncpz
                do l=1,nds+1
                    MP_b_net_final(iptch,i,j,k,l) = b_net_final(i,j,k,l) 
                end do
            end do
        end do
    end do
   
    deallocate(PT)
    
    continue
    
    
end subroutine
