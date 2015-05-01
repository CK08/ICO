!----------------------------------------------------------------------------------------------
!
! Subroutine compute the collocation points for the Point-to-Segment contact algorithm
!
!----------------------------------------------------------------------------------------------

subroutine PTS_CollocationPoints()
    
    use mod_variables
    implicit none
    integer(4)::i,j,ix,iy,mmult
    
    if(npair == 1)then
        !Number of Greville points in each direction
        ngrv_xi = ncpx_slv
        ngrv_eta = ncpy_slv
        
        !Contact collocation points in xi-direction
        allocate(grv_xi(ngrv_xi))
        grv_xi = 0.0d0
        do i=1,ngrv_xi
            do j=1,p_slv
                grv_xi(i) = grv_xi(i) + u_knot_slv(i+j)/(1.0d0*p_slv)
            end do 
        end do
        
        !Contact collocation points in eta-direction
        allocate(grv_eta(ngrv_eta))
        grv_eta = 0.0d0
        do i=1,ngrv_eta
            do j=1,q_slv
                grv_eta(i) = grv_eta(i) + v_knot_slv(i+j)/(1.0d0*q_slv)
            end do 
        end do
        
        allocate(Lagrange(ncpx_slv*ncpy_slv))
        Lagrange = 0.0d0
        
        allocate(dLagrange(ncpx_slv*ncpy_slv))
        dLagrange = 0.0d0
        
        allocate(FintLM(tnodes*nds,1))
        FintLM = 0.0d0
    else
        
        ix = 0
        iy = 0
        mmult = 0   
        do i=1,npair
            ix = ix + MB_ncpx_slv(i)
            iy = iy + MB_ncpy_slv(i)
            mmult = mmult + MB_ncpx_slv(i)*MB_ncpy_slv(i)
        end do
        
        allocate(Lagrange(mmult))
        Lagrange = 0.0d0
        
        allocate(dLagrange(mmult))
        dLagrange = 0.0d0
        
        allocate(FintLM(tnodes*nds,1))
        FintLM = 0.0d0
        
    end if    
    continue

end subroutine