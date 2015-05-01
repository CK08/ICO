!----------------------------------------------------------------------------------------------
!
! Subroutine to deallocate the variables in a multipatch analysis
!
!----------------------------------------------------------------------------------------------

subroutine MPDealloc(iptch)

    use mod_variables
    implicit none
    integer(4),intent(IN)::iptch
    
    deallocate(u_knot,v_knot,w_knot)
    deallocate(b_net)
    deallocate(b_net_final)
    
    deallocate(Points)
    deallocate(Weight)
    
end subroutine