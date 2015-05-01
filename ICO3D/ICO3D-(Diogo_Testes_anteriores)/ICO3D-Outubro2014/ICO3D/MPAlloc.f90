!----------------------------------------------------------------------------------------------
!
! Subroutine to allocate the variables in a multipatch analysis
!
!----------------------------------------------------------------------------------------------

subroutine MPAlloc(iptch,nzep)
    
    use Mod_Variables
    implicit none
    integer(4),intent(IN)::iptch,nzep
    integer(4)::nnod,i,j,k,l,count,nel,k1
    
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
    
    !Allocate knot vectors
    allocate(u_knot(ncpx+p+1),v_knot(ncpy+q+1),w_knot(ncpz+w+1))
    
    do i=1,ncpx+p+1
        u_knot(i) = MP_uknot(iptch,i)
    end do
    
    do i=1,ncpy+q+1
        v_knot(i) = MP_vknot(iptch,i)
    end do
    
    do i=1,ncpz+w+1
        w_knot(i) = MP_wknot(iptch,i)
    end do
    
    !Allocate control lattice
    allocate(b_net(ncpx,ncpy,ncpz,nds+1))
    allocate(b_net_final(ncpx,ncpy,ncpz,nds+1))
    
    do i=1,ncpx
        do j=1,ncpy
            do k=1,ncpz
                do l=1,nds+1
                    b_net(i,j,k,l) =  MP_b_net(iptch,i,j,k,l)
                    !b_net_final(i,j,k,l) =  MP_b_net_final(iptch,i,j,k,l)
                end do
            end do
        end do
    end do
    
    allocate(Points(ncpx*ncpy*ncpz,nds))
    allocate(Weight(ncpx*ncpy*ncpz))
    
    Weight = 0.0d0
    Points = 0.0d0
    
    nnod = 0
    Points = 0.0d0
    Weight = 0.0d0
    do k=1,ncpz
        do j = 1,ncpy
            do i = 1,ncpx
                nnod = nnod + 1
                Points(nnod,1) = b_net(i,j,k,1)
                Points(nnod,2) = b_net(i,j,k,2)
                Points(nnod,3) = b_net(i,j,k,3)
                Weight(nnod)   = b_net(i,j,k,4)
            end do
        end do
    end do
    
    !Allocate element data
    elcode = MP_elcode(iptch)
    npi = MP_npi(iptch)
    
    !allocate material props
    !allocate(props(iprops))
    props(:) = MP_props(iptch,:)
    
    call gen_IEN_INN()
    
!    count = (p+1)*(q+1)*(w+1) + 1
!    if(nlgeom==.true. .and. nptch .gt. 1)then
!        do i=1, ncpx+p
!            do j=1,ncpy+q
!                do k=1,ncpz+w
!                    if((u_knot(i) .ne. u_knot(i+1)) .and. (v_knot(j) .ne. v_knot(j+1)) .and. (w_knot(k) .ne. w_knot(k+1))) then
!                        
!                        nel = nel + 1
!                        
!                        count = (p+1)*(q+1)*(w+1) + 1
!                        
!                        do k1=1,(p+1)*(q+1)*(w+1)
!                            count = count - 1
!                            Points(k1,1) = Points(k1,1) + u(MP_IEN(nel,count)*3-2,1)
!                            Points(k1,2) = Points(k1,2) + u(MP_IEN(nel,count)*3-1,1)
!                            Points(k1,3) = Points(k1,3) + u(MP_IEN(nel,count)*3  ,1)
!                        end do
!                        
!                    end if        
!                end do
!            end do
!        end do    
!    end if
    
    
    continue
    
    
end subroutine