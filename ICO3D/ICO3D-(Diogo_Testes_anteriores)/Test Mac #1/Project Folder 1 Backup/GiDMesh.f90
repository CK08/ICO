!------------------------------------------------------------------------------------
!
! Subroutine to create the mesh file for GiD, based on the control polygon
!
!------------------------------------------------------------------------------------

subroutine GiDMesh()
    
    use Mod_Variables
    implicit none
    
    integer(4)::i,j,k,nodpelem,count
    integer(4)::el1,el2,el3,el4,el5,el6,el7,el8
    character(16)::my_type
    
    open(unit=2, file='IGA.flavia.msh')

    !---------------------------------------------------------------    
    ! Ciclo comentado devido a problemas de compilacao no MAC: Diogo
    !        ideia: trocar as condicoes todas por algo do genero
    !               elcode(primeiras 3 letras) == 'Hex'
    !---------------------------------------------------------------
!    if(elcode == 'Hex8' .or. elcode == 'Hex8BBar' .or. elcode == 'Hex27' .or. elcode == 'Hex27BBar' .or. elcode == 'Hex64' .or. elcode == 'Hex8ANS' .or. elcode == 'Hex27ANS' .or. elcode == 'Hex27EAS_PW' .or. elcode == 'Hex27EAS' .or. elcode == 'Hex27SRI' .or. elcode == 'Hex27PVS' .or. elcode == 'Hex27PV' .or. elcode == 'Hex27_ProjVol') then
!        my_type = "Hexahedra"
!        nodpelem = 8
!    end if
    
    write(2,FMT=111)nds, my_type, nodpelem
    
    !Heading ----------------------------------------------------------------
    111 format('MESH "problem" dimension ', I1,' ElemType ',A13, ' Nnode ', I2)
    
    !Coordinates ---------------------------
    count=0
    write(2,*)'coordinates'
    
    do k=1,ncpz
        do j=1,ncpy
            do i=1,ncpx
                count = count+1
                !Create polygon mesh...
                PMesh(count,:) = B_net(i,j,k,1:nds)
                !...and write to output file
                write(2,FMT=112)count,B_net(i,j,k,1:nds)
                112 format(I4,5(E,1x))
            end do
        end do
    end do
    write(2,*)'end coordinates'
    
    
    !Elements ------------------------------
    count = 0
    write(2,*)'Elements'
    do k=1,ncpz-1    
        do i=1,ncpy-1
            do j=1,ncpx-1
                count = count + 1
                el1 = (i-1)*ncpx + (k-1)*ncpx*ncpy + j
                el2 = el1 + 1
                el3 = el2 + ncpx
                el4 = el3 - 1
                el5 = (i-1)*ncpx + (k-1)*ncpx*ncpy + j + ncpx*ncpy
                el6 = el5 + 1
                el7 = el6 + ncpx
                el8 = el7 - 1
                write(2,FMT=113)count,el1,el2,el3,el4,el5,el6,el7,el8
                113 format(21(I4,1x))
            end do
        end do
    end do
    write(2,*)'end elements'
    
    close(2,status='keep')
    
    continue
    
    
end subroutine
