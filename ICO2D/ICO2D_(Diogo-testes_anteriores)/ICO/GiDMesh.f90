!------------------------------------------------------------------------------------
!
! Subroutine to create the mesh file for GiD, based on the control polygon
!
!------------------------------------------------------------------------------------

subroutine GiDMesh()
    
    use Mod_Variables
    implicit none
    
    integer(4)::i,j,nodpelem,count
    integer(4)::el1,el2,el3,el4
    character(16)::my_type
    
    open(unit=2, file='IGA.flavia.msh')
    
    if(elcode == 'Quad4E' .or. elcode == 'Quad9E' .or. elcode == 'Quad16E' .or. &
    &  elcode == 'Quad4S' .or. elcode == 'Quad9S' .or. elcode == 'Quad16S' .or. &
    &  elcode == 'Quad9EBBar') then
        my_type = "Quadrilateral"
        nodpelem = 4
    endif
    
    write(2,FMT=111)nds, my_type, nodpelem
    
    !Heading ----------------------------------------------------------------
    111 format('MESH "problem" dimension ', I1,' ElemType ',A13, ' Nnode ', I2)
    
    !Coordinates ---------------------------
    count=0
    write(2,*)'coordinates'
    
    do j=1,ncpy
        do i=1,ncpx
            count = count+1
            !Create polygon mesh...
            PMesh(count,:) = B_net(i,j,1:nds)
            !...and write to output file
            write(2,FMT=112)count,B_net(i,j,1:nds)
            112 format(I4,5(E,1x))
        end do
    end do
    write(2,*)'end coordinates'
    
    
    !Elements ------------------------------
    count = 0
    write(2,*)'Elements'
        do i=1,ncpy-1
            do j=1,ncpx-1
                count = count + 1
                el1 = (i-1)*ncpx + j
                el2 = el1 + 1
                el3 = el2 + ncpx
                el4 = el3 - 1
                write(2,FMT=113)count,el1,el2,el3,el4
                113 format(21(I4,1x))
            end do
        end do
    write(2,*)'end elements'
    
    close(2,status='keep')
    
    continue
    
    
end subroutine