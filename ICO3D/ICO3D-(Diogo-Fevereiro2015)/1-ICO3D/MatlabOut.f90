!---------------------------------------------------------------------------------------
!
! Subroutine to writte data of the deformed geometry for the matlab plot subroutine 
!
!---------------------------------------------------------------------------------------

subroutine MatlabOut

    use Mod_Variables
    implicit none
    integer(4)::i,j,k,l, count, nze,nnod
    integer(4)::k1,k2,k3,k4,k5,k6
    
    if(nptch .gt. 1)then
        open(unit=7,file='NURBSOutput.txt')
        
        write(7,FMT=100)nptch
        
        do i=1,nptch
            write(7,FMT=101)MP_ncpx(i), MP_ncpy(i), MP_ncpz(i)
        end do
        
        do i=1,nptch
            write(7,FMT=101)MP_p(i), MP_q(i), MP_w(i)
        end do
        
        do i=1,nptch
            
            write(7,FMT=102)MP_uknot(i,1:MP_ncpx(i)+MP_p(i)+1)
            write(7,FMT=102)MP_vknot(i,1:MP_ncpy(i)+MP_q(i)+1)
            write(7,FMT=102)MP_wknot(i,1:MP_ncpz(i)+MP_w(i)+1)
            
            do j=1,MP_ncpx(i)
                do k=1,MP_ncpy(i)
                    do l=1,MP_ncpz(i)
                        write(7,FMT=102) MP_b_net_final(i,j,k,l,:)
                    end do
                end do
            end do 
        end do
        
        close(unit=7)
        
        !
        open(unit=8,File='GCoords_Final.txt')
        write(8,FMT=100)tnodes
        
        do i=1,tnodes
            write(8,FMT=102)GCoords(i,1:nds+1)
        end do
        
        
        close(8)
        
        
        100 format(I4,',')
        101 format(3(I4,','))
        102 format(100(E,','))
    else
        open(unit=7,file='NURBSOutput.txt')
        
        write(7,FMT=200)nptch
        
        write(7,FMT=201)ncpx, ncpy, ncpz
        
        write(7,FMT=201)p, q, w
        
        write(7,FMT=202)u_knot(:)
        write(7,FMT=202)v_knot(:)
        write(7,FMT=202)w_knot(:)
        
        nnod = 0
        do k=1,ncpz
            do j=1,ncpy
                do i=1,ncpx
                    nnod = nnod + 1
                    b_net_final(i,j,k,nds+1) = weight(nnod)
                end do
            enddo
        enddo
      
        nnod = 0
        do k=1,ncpz
            do j=1,ncpy
                do i=1,ncpx 
                    nnod = nnod + 1
                    b_net_final(i,j,k,1) = Points(nnod,1)
                    b_net_final(i,j,k,2) = Points(nnod,2)
                    b_net_final(i,j,k,3) = Points(nnod,3)
                end do
            enddo
        enddo
        
        open(unit=8,File='GCoords_Final.txt')
        write(8,FMT=100)tnodes
        
        do i=1,tnodes
            write(8,FMT=108)Points(i,1),Points(i,2),Points(i,3),weight(i)
        end do
        
        close(8)
            
        do j=1,ncpx
            do k=1,ncpy
                do l=1,ncpz
                    write(7,FMT=202) b_net_final(j,k,l,:)
                end do
            end do
        end do
        
        
        108 format(100(E,','))
        200 format(I4,',')
        201 format(3(I4,','))
        202 format(100(E,','))
        
        close(unit=7)
        
        continue
        
        
    end if
    
!    open(unit=7,file='GPCoords.txt')
!        
!        write(7,FMT=203) telems*npi
!        
!        do i=1,telems*npi
!            write(7,FMT=204) GP_Coords(i,:)
!        end do
!        
!        203 format(I4,',')
!        204 format(100(E,','))
!        
!        close(unit=7)
!        
!        open(unit=8,file='GPStress.txt')
!        open(unit=9,file='GPStrain.txt')
!        
!        do i=1,telems
!            do j=1,npi
!                write(8,FMT=205) Stress_Conv(i,j,:)
!                write(9,FMT=205) Strain_Conv(i,j,:)
!            end do
!        end do
!        
!        205 format(100(E,','))
!        
!        close(unit=8)
!        close(unit=9)
!    
!    continue
end subroutine