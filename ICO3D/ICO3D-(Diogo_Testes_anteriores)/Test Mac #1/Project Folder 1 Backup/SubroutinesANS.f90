subroutine kvecANS(ncpx,p,u_knot,uk,ub)

    integer(4),intent(IN)::ncpx,p
    real(8),intent(IN)::uk
    real(8),dimension(ncpx+p+1),intent(IN)::u_knot
    real(8),dimension(ncpx+p+1+p),intent(OUT)::ub
    logical::lsp,endcycly

    ub = 0.0d0
    count = p+2
    lsp = .false.
    endcycle = .false.
    do k2=p+2,ncpx+p+1
        if(u_knot(k2) .gt. uk .and. lsp==.false.) then
            ub(count) = uk
            count = count + 1
            ub(count) = uk
            count = count + 1
            lsp = .true.
            do k3=k2,ncpx+p+1
                ub(count) = u_knot(k3)
                count = count + 1
            end do
            endcycle = .true.
        elseif(endcycle==.false.)then
            ub(count) = u_knot(k2)
            count = count + 1
        end if
    end do

    continue

end subroutine