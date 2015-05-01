subroutine PTS_CollocationPoints()
    
    use Mod_variables
    implicit none
    integer(4)::i,j
    
    !Calculate Greville Points on the parametric surface
    ngrv = ncp_s
    allocate(grev(ngrv),initial_gap(ngrv),final_gap(ngrv))
    
    allocate(CP_area(ngrv),PTS_Stress(ngrv))
    CP_area = 0.0d0
    PTS_Stress = 0.0d0
    
    allocate(PTS_active(ngrv))
    
    allocate(plot_xy(ngrv,2))
    plot_xy = 0.0d0
    
    grev = 0.0d0
    initial_gap = 0.0d0
    final_gap = 0.0d0
    
    PTS_active = 0
    
    grev = 0.0d0
    do i=1,ngrv
        do j=1,p_s
            grev(i) = grev(i) + u_knot_slave(i+j)/(1.0d0*p_s)
        end do
    end do
    
    !Initialize Lagrange Multipliers
    allocate(Lagrange(ncp_s))
    Lagrange = 0.0d0
    
    allocate(LagrangeConv(ncp_s))
    LagrangeConv = 0.0d0
    
    allocate(dLagrange(ncp_s))
    dLagrange = 0.0d0
    
    allocate(ddLagrange(ncp_s))
    ddLagrange = 0.0d0
    
    allocate(Lagmx(ncp_s))
    Lagmx = 0.0d0
    
    allocate(Lagmy(ncp_s))
    Lagmy = 0.0d0
    
    allocate(FintLMi(tnodes*nds,ncp_s))
    FintLMi = 0.0d0
    
    
    continue
    
end subroutine