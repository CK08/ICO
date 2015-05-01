
subroutine MPAlloc(iptch)
    
    use Mod_Variables
    implicit none
    integer(4),intent(IN)::iptch
    integer(4)::nnod,i,j,k,count
    
    p = order(iptch,1)
    q = order(iptch,2)
    ncpx = ncpptch(iptch,1)
    ncpy = ncpptch(iptch,2)
    nnodes = ncpx*ncpy
    nelems = (ncpx-p)*(ncpy-q)
    nshpl = (p+1)*(q+1)
    
    allocate(u_knot(ncpx+p+1),v_knot(ncpy+q+1))
    u_knot(:) = MP_u_knot(iptch,:)
    v_knot(:) = MP_v_knot(iptch,:)
    
    allocate(b_net(ncpx,ncpy,nds+1))
    allocate(b_net_final(ncpx,ncpy,nds+1))
    b_net(:,:,:) =  MP_b_net(iptch,:,:,:)
    b_net_final(:,:,:) =  MP_b_net_final(iptch,:,:,:)
    
    allocate(Points(ncpx*ncpy,nds))
    allocate(Points0(ncpx*ncpy,nds))
    allocate(Weights(ncpx*ncpy))
    
    Weights = 0.0d0
    Points = 0.0d0
    Points0 = 0.0d0
    
    elcode = MP_elcode(iptch)
    npi = MP_npi(iptch)
    
    nnod = 0
    Points = 0.0d0
    Weights = 0.0d0
    do j = 1,ncpy
        do i = 1,ncpx
            nnod = nnod + 1
            Points(nnod,1) = b_net(i,j,1)
            Points(nnod,2) = b_net(i,j,2)
            Points0(nnod,1) = b_net(i,j,1)
            Points0(nnod,2) = b_net(i,j,2)
            Weights(nnod)  = b_net(i,j,3)
        end do
    end do  
    
    allocate(props(iprops))
    props(:) = MP_props(iptch,:)
    
    allocate(PMesh(ncpx*ncpy,nds))
    allocate(PMesh0(ncpx*ncpy,nds))
    PMesh = 0.0d0
    PMesh0 = 0.0d0
    
    count = 0
    do j=1,ncpy
        do i=1,ncpx
            count = count+1
            PMesh(count,:) = B_net(i,j,1:nds)
            PMesh0(count,:) = B_net(i,j,1:nds)
        end do
    end do

end subroutine  


subroutine MPDalloc(iptch)
    
    use Mod_Variables
    implicit none
    integer(4),intent(IN)::iptch
    integer(4)::i,j,k
    
    MP_b_net_final(iptch,:,:,:) = 0.0d0
    
    do i=1,ncpx
        do j=1,ncpy
            do k=1,nds+1
                MP_b_net_final(iptch,i,j,k) = b_net_final(i,j,k)
            end do
        end do
    end do
            
            
    deallocate(u_knot,v_knot,b_net,b_net_final)
    deallocate(Points,Weights,Points0)
    
    deallocate(props)
    
    deallocate(PMesh,PMesh0)
    
end subroutine   