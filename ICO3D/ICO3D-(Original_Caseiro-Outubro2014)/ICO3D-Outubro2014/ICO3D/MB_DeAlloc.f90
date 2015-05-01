!----------------------------------------------------------------------------------------------
!
! Subroutine to deallocate varaiables for a contact analysis involving more than one
! contact pair
!
!----------------------------------------------------------------------------------------------

subroutine MB_DeAlloc(ipair)
    
    use mod_variables
    use mod_MBvariables
    
    implicit none
    
    integer(4),intent(IN)::ipair
    
    if(npair .gt. 1)then
        
        deallocate(u_knot_mst)
        deallocate(v_knot_mst)
        deallocate(u_knot_slv)
        deallocate(v_knot_slv)
        deallocate(conn_mst)
        deallocate(conn_slv)
        deallocate(grv_xi)
        deallocate(grv_eta) 
        deallocate(PTS_conn)   
    
    else
        
        deallocate(gap)
    
    end if
    
    deallocate(Nslv,dNslvdxi,dNslv2dxi2)
    deallocate(Mslv,dMslvdeta,dMslv2deta2)
        
    deallocate(Rslv,dRslv,ddRslv)
    deallocate(Rmst,dRmst,ddRmst)
    
    deallocate(N,Na,Nb,Ta,Tb)
    deallocate(Nhat,That,D,Nbar)
    deallocate(Kgeo,KG1,KG2,KG3,KG4)

end subroutine