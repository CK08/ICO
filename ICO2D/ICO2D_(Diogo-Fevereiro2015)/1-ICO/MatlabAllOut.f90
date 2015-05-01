subroutine MatlabAllOut(inc,iter)

    use Mod_Variables
    implicit none
    integer(4),intent(IN)::inc,iter
    integer(4)::i,j,k,l

    if(inc == 1 .and. iter == 1)then
        open(unit=8,file='NURBSAllOut.txt', status = 'REPLACE')
        open(unit=7,file='GPCoordsAll.txt', status = 'REPLACE')
        open(unit=10,file='GPStressAll.txt', status = 'REPLACE')
        open(unit=9,file='GPStrainAll.txt', status = 'REPLACE')
    else
        open(unit=8,file='NURBSAllOut.txt', access = 'APPEND')
        open(unit=7,file='GPCoordsAll.txt', access = 'APPEND')
        open(unit=10,file='GPStressAll.txt', access = 'APPEND')
        open(unit=9,file='GPStrainAll.txt', access = 'APPEND')

    end if

    if(nptch .gt. 1)then
        !open(unit=8,file='NURBSAllOut.txt', access = 'APPEND')

        write(8,FMT=100)nptch

        do i=1,nptch
            write(8,FMT=101)ncpptch(i,1), ncpptch(i,2)
        end do

        do i=1,nptch
            write(8,FMT=101)order(i,1), order(i,2)
        end do

        do i=1,nptch

            write(8,FMT=102)MP_u_knot(i,1:ncpptch(i,1)+order(i,1)+1)
            write(8,FMT=102)MP_v_knot(i,1:ncpptch(i,2)+order(i,2)+1)

            do j=1,ncpptch(i,1)
                do k=1,ncpptch(i,2)
                    write(8,FMT=102) MP_b_net_final(i,j,k,:)
                end do
            end do
        end do

        100 format(I4,',')
        101 format(2(I4,','))
        102 format(100(E,','))

        close(unit=8)

    else
        !open(unit=8,file='NURBSAllOut.txt', status = 'REPLACE')

        write(8,FMT=200)nptch

        write(8,FMT=201)ncpx, ncpy

        write(8,FMT=201)p, q

        write(8,FMT=202)u_knot(:)
        write(8,FMT=202)v_knot(:)

        do j=1,ncpx
            do k=1,ncpy
                write(8,FMT=202) b_net_final(j,k,:)
            end do
        end do

        200 format(I4,',')
        201 format(2(I4,','))
        202 format(100(E,','))

        close(unit=8)

    end if

    !open(unit=7,file='GPCoordsAll.txt')

    write(7,FMT=203) telems*npi

    do i=1,telems*npi
        write(7,FMT=204) GP_Coords(i,:)
    end do

    203 format(I4,',')
    204 format(100(E,','))

    close(unit=7)

    !open(unit=10,file='GPStressAll.txt')
    !open(unit=9,file='GPStrainAll.txt')

    do i=1,telems
        do j=1,npi
            write(10,FMT=205) Stress_Conv(i,j,:)
            write(9,FMT=205) Strain_Conv(i,j,:)
        end do
    end do

    205 format(100(E,','))

    close(unit=10)
    close(unit=9)


end subroutine