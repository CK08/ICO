!--------------------------------------------------------------------
!
! Subroutine to determine the connectivity arrays for 
! reduced basis of the B-Bar method
!
!--------------------------------------------------------------------

subroutine gen_ien_inn_bar
    
    use Mod_Variables
    
    implicit none
    
    logical :: prnt
    integer :: i, j, k, loop1,loop2,loop3,g, e, gtemp, ln

    !allocate(INN(nnodes,nds))
    allocate(IENb(nelems,(p)*(q)*(w)))

    !Initialize matrices and variables

    IENb = 0
    !INN = 0
    g   = 0
    e   = 0

    do k = 1,ncpz-1
        do j = 1,ncpy-1           ! loop through control points in V direction
            do i = 1,ncpx-1      ! loop through control points in U direction
                g = g + 1       ! increment global function number
                !INN(g,1) = i
                !INN(g,2) = j
                !INN(g,3) = k
                if(i .gt. p-1 .and. j .gt. q-1 .and. k .gt. w-1) then
                    e = e + 1                       ! increment element number

                    do loop1 = 0,w-1
                        do loop2 = 0,q-1
                            do loop3 = 0,p-1
                                gtemp     = g - (ncpx-1)*(loop1*(ncpy-1) + loop2) - loop3
                                ln        = (p)*((q)*loop1 + loop2) + loop3 + 1
                                IENb(e,ln) = gtemp
                            enddo ! loop3
                        enddo ! loop2
                    enddo !loop1

                endif
            enddo ! i
        enddo ! j
    enddo !k
    
    
    continue
    
    
end subroutine gen_ien_inn_bar 