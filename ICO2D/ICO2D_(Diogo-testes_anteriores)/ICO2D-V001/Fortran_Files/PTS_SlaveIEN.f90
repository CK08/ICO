subroutine PTS_SlaveIEN(  p,nnodes, nelems,nds,      u_knot,INN,IEN)
                        !p_s, ncp_s,ncp_s-p,nds,u_knot_slave,INN_s,IEN_s
    implicit none
    integer :: i, j, k, loop1,loop2,loop3,g, e, gtemp, ln
    
    integer(4), intent(IN):: p,nnodes,nelems,nds
    real(8), dimension(nnodes+p+1), intent(IN):: u_knot
    integer(4), dimension(nnodes,nds), intent(OUT) :: INN
    integer(4), dimension(nelems,p+1), intent(OUT) :: IEN
    
    IEN = 0
    INN = 0
    g   = 0
    e   = 0
    
    g = 0
    do i = 1,nnodes
        g = g + 1
        INN(g,1) = i
        if((u_knot(i) .ne. u_knot(i+1))) then
            if(i.gt.p) then
                e = e + 1
                do loop3 = 0,p
                    gtemp     = g - loop3
                    ln        = loop3 + 1
                    IEN(e,ln) = gtemp
                enddo 
            endif
        end if
    enddo 
    
    continue
    
    
    
end subroutine