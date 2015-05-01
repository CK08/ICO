subroutine GPTS_Initialisation()

    use Mod_variables
    implicit none
    
    integer(4)::icel,i,j
    
    icel = 0
    do i=1, ncp_s+p_s
        if(u_knot_slave(i) .ne. u_knot_slave(i+1)) icel = icel + 1
    end do   
    
    allocate(initial_gap(icel*SurfGP),final_gap(icel*SurfGP))
    initial_gap = 0.0d0
    final_gap = 0.0d0
    
    allocate(PTS_active(icel*SurfGP))
    PTS_active = 0
    
    allocate(Press(icel*SurfGP))
    Press = 0.0d0
    
    allocate(Fct(tnodes*nds,1))
    Fct = 0.0d0
   
   

end subroutine