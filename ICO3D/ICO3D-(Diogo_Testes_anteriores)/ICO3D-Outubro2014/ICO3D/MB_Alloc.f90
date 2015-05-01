!----------------------------------------------------------------------------------------------
!
! Subroutine to allocate varaiables for a contact analysis involving more than one
! contact pair
!
!----------------------------------------------------------------------------------------------

subroutine MB_alloc(ipair)
    
    use mod_variables
    use mod_MBvariables
    
    implicit none
    
    integer(4),intent(IN)::ipair
    integer(4)::i,j,ix,iy,ix2,iy2,is,im,isurf,k1,count
    integer(4)::imult,imult2
    
    if(npair .gt. 1)then
        
        p_slv = MB_p_slv(ipair)
        q_slv = MB_q_slv(ipair)
        
        p_mst = MB_p_mst(ipair)
        q_mst = MB_q_mst(ipair)
        
        MRigid = MB_MRigid(ipair)
        
        islave = MB_islave(ipair)
        
        ncpx_mst = MB_ncpx_mst(ipair)
        ncpy_mst = MB_ncpy_mst(ipair)
        
        ncpx_slv = MB_ncpx_slv(ipair)
        ncpy_slv = MB_ncpy_slv(ipair)
        
        !--------------------------------------------------------
        ! Knot vectors of the master surface
        !--------------------------------------------------------
        allocate(u_knot_mst(ncpx_mst+p_mst+1))
        u_knot_mst = 0.0d0
        do i=1,ncpx_mst+p_mst+1
            u_knot_mst(i) = MB_u_knot_mst(ipair,i)
        end do
        
        allocate(v_knot_mst(ncpy_mst+q_mst+1))
        v_knot_mst = 0.0d0
        do i=1,ncpy_mst+q_mst+1
            v_knot_mst(i) = MB_v_knot_mst(ipair,i)
        end do  
        
        !--------------------------------------------------------
        ! Knot vectors of the slave surface 
        !--------------------------------------------------------
        allocate(u_knot_slv(ncpx_slv+p_slv+1))
        u_knot_slv = 0.0d0
        do i=1,ncpx_slv+p_slv+1
            u_knot_slv(i) = MB_u_knot_slv(ipair,i)
        end do
        
        allocate(v_knot_slv(ncpy_slv+q_slv+1))
        v_knot_slv = 0.0d0
        do i=1,ncpy_slv+q_slv+1
            v_knot_slv(i) = MB_v_knot_slv(ipair,i)
        end do    
        
        !--------------------------------------------------------
        ! Surface connectivities
        !--------------------------------------------------------
        allocate(conn_mst(ncpx_mst*ncpy_mst))
        conn_mst = 0
        do i=1,ncpx_mst*ncpy_mst
            conn_mst(i) = MB_conn_mst(ipair,i)
        end do
        
        allocate(conn_slv(ncpx_slv*ncpy_slv))
        conn_slv = 0
        do i=1,ncpx_slv*ncpy_slv
            conn_slv(i) = MB_conn_slv(ipair,i)
        end do
        
        !--------------------------------------------------------
        ! Number of Greville points in each direction
        !--------------------------------------------------------
        ngrv_xi = ncpx_slv
        ngrv_eta = ncpy_slv
        
        !--------------------------------------------------------
        ! Contact collocation points in xi-direction
        !--------------------------------------------------------
        allocate(grv_xi(ngrv_xi))
        grv_xi = 0.0d0
        do i=1,ngrv_xi
            do j=1,p_slv
                grv_xi(i) = grv_xi(i) + u_knot_slv(i+j)/(1.0d0*p_slv)
            end do 
        end do
        
        !--------------------------------------------------------
        ! Contact collocation points in eta-direction
        !--------------------------------------------------------
        allocate(grv_eta(ngrv_eta))
        grv_eta = 0.0d0
        do i=1,ngrv_eta
            do j=1,q_slv
                grv_eta(i) = grv_eta(i) + v_knot_slv(i+j)/(1.0d0*q_slv)
            end do 
        end do
        
        if(ipair==1)then
            
            ix = 0
            iy = 0
            imult = 0
            do i=1,npair
                ix = ix + MB_ncpx_slv(i)
                iy = iy + MB_ncpy_slv(i)
                imult = imult + MB_ncpx_slv(i)*MB_ncpy_slv(i)
            end do
            
            ix2 = 0
            iy2 = 0
            imult2 = 0
            do i=1,npair
                ix2 = ix2 + MB_ncpx_mst(i)
                iy2 = iy2 + MB_ncpy_mst(i)
                imult2 = imult + MB_ncpx_mst(i)*MB_ncpy_mst(i)
            end do
        
            allocate(gap(imult))
            gap = 0.0d0
            
        end if 

        is = imult
        im = imult2
        incp = imult + imult2
        
        isurf = ncpx_slv*ncpy_slv + ncpx_mst*ncpy_mst
          
    else
        
        allocate(gap(ncpx_slv*ncpy_slv))
        gap = 0.0d0
        
        is = ncpx_slv*ncpy_slv
        im = ncpx_mst*ncpy_mst
        incp = ncpx_slv*ncpy_slv + ncpx_mst*ncpy_mst
        isurf = ncpx_slv*ncpy_slv + ncpx_mst*ncpy_mst
        
    end if
    
    allocate(Nslv(p_slv+1),dNslvdxi(p_slv+1),dNslv2dxi2(p_slv+1))
    allocate(Mslv(q_slv+1),dMslvdeta(q_slv+1),dMslv2deta2(q_slv+1))
    
    
    allocate(Rslv(ncpx_slv*ncpy_slv))
    allocate(dRslv(ncpx_slv*ncpy_slv,2))
    allocate(ddRslv(ncpx_slv*ncpy_slv,4))
    
    allocate(Rmst(ncpx_mst*ncpy_mst))
    allocate(dRmst(ncpx_mst*ncpy_mst,2))
    allocate(ddRmst(ncpx_mst*ncpy_mst,4))
    
    Nslv = 0.0d0
    dNslvdxi = 0.0d0
    dNslv2dxi2 = 0.0d0
    
    Mslv = 0.0d0
    dMslvdeta = 0.0d0
    dMslv2deta2 = 0.0d0
    
    Rslv = 0.0d0
    dRslv = 0.0d0
    ddRslv = 0.0d0    
    
    Rmst = 0.0d0
    dRmst = 0.0d0
    ddRmst = 0.0d0
    
    allocate(N(isurf*nds,1),Na(isurf*nds,1),Nb(isurf*nds,1))
    N = 0.0d0
    Na = 0.0d0
    Nb = 0.0d0
    
    allocate(Ta(isurf*nds,1),Tb(isurf*nds,1))
    Ta = 0.0d0
    Tb = 0.0d0
    
    allocate(Nhat(isurf*nds,2),That(isurf*nds,2),D(isurf*nds,2),Nbar(isurf*nds,2))
    Nhat = 0.0d0
    That = 0.0d0
    Nbar = 0.0d0
    
    allocate(Kgeo(isurf*nds,isurf*nds))
    Kgeo = 0.0d0
    
    allocate(Kg1(isurf*nds,isurf*nds))
    allocate(Kg2(isurf*nds,isurf*nds))
    allocate(Kg3(isurf*nds,isurf*nds))
    allocate(Kg4(isurf*nds,isurf*nds))
    Kg1 = 0.0d0
    Kg2 = 0.0d0
    Kg3 = 0.0d0
    Kg4 = 0.0d0
    
    allocate(PTS_conn(isurf))
    PTS_conn = 0
    do k1=1,ncpx_slv*ncpy_slv
        PTS_conn(k1) = conn_slv(k1)
    end do
    
    count = 0
    do k1=ncpx_slv*ncpy_slv + 1,isurf
        count = count + 1
        PTS_conn(k1) = conn_mst(count)
    end do
    
    if(ipair==1)then
        allocate(KLM(tnodes*nds+is,tnodes*nds+is))
        KLM = 0.0d0
    
        allocate(FLM(tnodes*nds,is))
        FLM = 0.0d0
        
        allocate(GLM(tnodes*nds+is,1))
        GLM = 0.0d0
        
        allocate(disptempLM(tnodes*nds+is))
        disptempLM = 1
        
        allocate(KeqLM(tnodes*nds+is-nbc,tnodes*nds+is-nbc))
        KeqLM = 0.0d0
        
        allocate(FeqLM(tnodes*nds+is-nbc,1))
        FeqLM = 0.0d0
        
        allocate(DispeqLM(tnodes*nds+is-nbc,1))
        DispeqLM = 0.0d0
        
        allocate(PTS_Active(incp))
        PTS_Active = 0
    end if
    
                
end subroutine