!--------------------------------------------------------------------
!
! Subroutine to determine the connectivity arrays
! as presented in the book
! Isogeometric Analysis -Toward Integration of CAD and FEA
!
!--------------------------------------------------------------------

subroutine gen_ien_inn
    
    use Mod_Variables
    
    implicit none
    
    logical :: prnt
    integer :: i, j, k, loop1,loop2,loop3,g, e, gtemp, ln

    allocate(INN(nnodes,nds))
    allocate(IEN(nelems,(p+1)*(q+1)*(w+1)))

    !Initialize matrices and variables

    IEN = 0
    INN = 0
    g   = 0
    e   = 0

    do k = 1,ncpz
        do j = 1,ncpy           ! loop through control points in V direction
            do i = 1,ncpx       ! loop through control points in U direction
                g = g + 1       ! increment global function number
                INN(g,1) = i
                INN(g,2) = j
                INN(g,3) = k
                if((u_knot(i) .ne. u_knot(i+1)) .and. (v_knot(j) .ne. v_knot(j+1)) .and. (w_knot(k) .ne. w_knot(k+1))) then
                    if(i .gt. p .and. j .gt. q .and. k .gt. w) then
                        e = e + 1                       ! increment element number

                        do loop1 = 0,w
                            do loop2 = 0,q
                                do loop3 = 0,p
                                    gtemp     = g - ncpx*(loop1*ncpy + loop2) - loop3
                                    ln        = (p+1)*((q+1)*loop1 + loop2) + loop3 + 1
                                    IEN(e,ln) = gtemp
                                    continue
                                enddo ! loop3
                            enddo ! loop2
                        enddo !loop1

                    endif
                end if
            enddo ! i
        enddo ! j
    enddo !k
    
!    do i=1,nelems
!        g = (p+1)*(q+1)
!        do j=1,(p+1)*(q+1)
!            conn(i,j) = ien(i,g)
!            g = g - 1
!        end do
!    end do
    
    conn = ien
    
    continue
    
    
end subroutine gen_ien_inn    