!----------------------------------------------------------------------------------------------
!
! Subroutine allocate global variables
!
!----------------------------------------------------------------------------------------------

subroutine AllocGlobalVars()

    use mod_variables
    
    implicit none
    
    integer(4)::k1

    allocate(u_knot(ncpx+p+1),v_knot(ncpy+q+1),w_knot(ncpz+w+1))
    u_knot = 0
    v_knot = 0
    w_knot = 0
    
    allocate(b_net(ncpx,ncpy,ncpz,nds+1))
    b_net = 0.0d0
    
    allocate(b_net_final(ncpx,ncpy,ncpz,nds+1))
    b_net_final = 0.0d0
    
    !allocate(bw(ncpx,ncpy,ncpz,nds+1))
    !bw = 0.0d0

    allocate(iper(nnodes))
    do k1 = 1, nnodes
        iper(k1) = k1
    enddo
    
    allocate(bc(ncpx,ncpy,ncpz,nds),load(ncpx,ncpy,ncpz,nds))
    bc = 1
    load = 0.0d0
    
    allocate(dispBC(ncpx,ncpy,ncpz,nds))
    dispBC = 0.0d0

    allocate(ld_el(nelems,nds))
    ld_el = 0.0d0
    allocate(ey(nelems))
    ey = 0.0d0
    allocate(nu(nelems))
    nu = 0.0d0
    
    !allocate(Ke((p+1)*(q+1)*nds,(p+1)*(q+1)*nds)
    allocate(Kf(ncpx*ncpy*ncpz*nds,ncpx*ncpy*ncpz*nds))
    Kf = 0.0d0
    
    allocate(Fext(nnodes*nds,1))
    Fext = 0.0d0
    
    allocate(ImpDisp(nnodes*nds,1),StpFrac(nnodes*nds))
    ImpDisp = 0.0d0
    StpFrac = 0.0d0
    
    allocate(Fini(nnodes*nds,1))
    Fini = 0.0d0
    
    allocate(FInc(nnodes*nds,1))
    FInc = 0.0d0
    
    allocate(Res(nnodes*nds,1))
    Res = 0.0d0
    
    allocate(Fint(nnodes*nds,1))
    Fint = 0.0d0
    
    allocate(RedCount(nnodes*nds))
    redcount = 1
    
    allocate(PMesh(ncpx*ncpy*ncpz,nds))
    PMesh = 0.0d0
    
    allocate(Weight(ncpx*ncpy*ncpz))
    Weight = 0.0d0
    
    allocate(Points(ncpx*ncpy*ncpz,nds))
    Points = 0.0d0
    
    allocate(loaddof(nnodes*nds))
    loaddof = 0.0d0
    
    allocate(dispdof(nnodes*nds))
    dispdof = 0.0d0
    
    allocate(bcdof(tnodes*nds))
    bcdof = 1

end subroutine