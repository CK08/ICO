!------------------------------------------------------------------------------------
!
! Subroutine to export control points displacements to GiD
!
!------------------------------------------------------------------------------------

subroutine GiDRes(inc)
    
    use Mod_Variables
    implicit none
    
    integer(4),intent(IN)::inc
    integer(4)::k1,count,i,j,k
    
    open(unit=3,file='IGA.flavia.res')
    write(3,*)'GiD Post Results File 1.0'
    
    !Control knots displacements
    write(3,FMT=115)inc
    115 format('Result "Displacements" "Load Analysis"', I3, ' Vector OnNodes')
    if (nds==2) write(3,*)'ComponentNames "X-Disp", "Y-Disp"'
    !if (nds==3) write(13,*)'ComponentNames "X-Disp", "Y-Disp", "Z-Disp"'
    write(3,*)'Values'
    
    b_net_final = b_net
    
    count = 0
    do k=1,ncpz
        do j=1,ncpy
            do i=1,ncpx
                count = count + 1
                !Update Polygon Mesh
                do k1=1,nds
                    b_net_final(i,j,k,k1)=b_net_final(i,j,k,k1) + u(count*3-(nds-k1),1)
                end do
                
                PMesh(count,1) = PMesh(count,1) + u(count*3-2,1)
                PMesh(count,2) = PMesh(count,2) + u(count*3-1,1)
                PMesh(count,3) = PMesh(count,3) + u(count*3  ,1)
                !Write control knots displacements
                write(3,FMT=114)count, u(count*3-2,1), u(count*3-1,1), u(count*3,1)
                114 format(I5,3(E))
            end do
        end do
    end do
    
!    do k1=1,nnodes/4
!        if (nds==2)then
!            write(3,FMT=114)(k1*4-3), u((k1*4-3)*2-1,1), u((k1*4-3)*2,1)
!            write(3,FMT=114)(k1*4-1), u((k1*4-2)*2-1,1), u((k1*4-2)*2,1)
!            write(3,FMT=114)(k1*4  ), u((k1*4  )*2-1,1), u((k1*4  )*2,1)
!            write(3,FMT=114)(k1*4-2), u((k1*4-1)*2-1,1), u((k1*4-1)*2,1)
!        endif
!        !if (nds==3)write(13,FMT=114)k1, u(k1*3-2,1), u(k1*3-1,1),u(k1*3,1)
!        114 format(I3,3(E))
!    end do
    write(3,*)'End Values'
    
    close(3,status='keep')
    
    continue
    
end subroutine