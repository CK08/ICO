!------------------------------------------------------------------------------------
!
! Subroutine to create the mesh file for GiD, based on the control polygon
!
!------------------------------------------------------------------------------------

subroutine GiDEqMesh(coord)
    
    use Mod_Variables
    implicit none
    
    integer(4)::i,j,l,k,k1,uind,nodpelem,count,nbx,nby,counter
    integer(4)::nelx,nely,el1,el2,el3,el4
    real(8)::ucheck,vcheck,den
    integer(4)::uspan,vspan,idx
    real(8),dimension(p+1,1)::BFx
    real(8),dimension(q+1,1)::BFy
    real(8),dimension(1,1,nds+1)::temp
    character(16)::my_type
    
    real(8),dimension(q+1)::Ny
    real(8),dimension(q+1)::dNydxi
    
    !Nodes coordinates (size overestimated)
    real(8),dimension(ncpx*ncpy,nds)::coord
    
    open(unit=2, file='IGAEq.flavia.msh')
    
    if(elcode == 'Quad4E' .or. elcode == 'Quad9E' .or. elcode == 'Quad4S' .or. elcode == 'Quad9S') then
        my_type = "Quadrilateral"
        nodpelem = 4
    endif
    
    write(2,FMT=111)nds, my_type, nodpelem
    
    !Heading ----------------------------------------------------------------
    111 format('MESH "problem" dimension ', I1,' ElemType ',A13, ' Nnode ', I2)
    
    !Coordinates ---------------------------
    counter=0
    write(2,*)'coordinates'
    
    nbx = ncpx+p+1
    nby = ncpy+q+1
    
    nelx = ncpx - p
    nely = ncpy - q
    
    coord = 0.0d0
!    do i=1,ncpy
!        do j=1,nds
!            coord(i,j) = b_net(1,i,j)
!        enddo
!    end do
    
    !--------------------------------------------------------------------------------
    !u = 0.0
    ucheck = 0.0d0
    call FindSpan(nbx,p,ucheck,U_Knot,uspan)  
    BFx = 0.0d0
    call BasisFuns(nbx,uspan,ucheck,p,u_knot,BFx) 
    
    vcheck = 0.0d0   
    call FindSpan(nby,q,vcheck,V_Knot,vspan)
    BFy = 0.0d0
    call BasisFuns(nby,vspan,vcheck,q,v_knot,BFy)
    
    counter = counter + 1
    den = 0.0;
    j = 1
    do l=0,q
        temp = 0.0d0
        do k=0,p
            do k1=1,nds+1
                temp(1,1,k1)  = temp(1,1,k1)  + BFx(k+1,1)*Bw(uspan-p+k, vspan-q+l, k1)
            enddo
            continue
        enddo
        
        den = den + temp(1,1,nds+1)*BFy(l+1,1);
        do k1=1,nds
            coord(counter,k1) = coord(counter,k1) + temp(1,1,k1)*BFy(l+1,1);
        enddo
    enddo
    
    do j=1,ncpy+q
        if(v_knot(j) .ne. v_knot(j+1))then
            vcheck = v_knot(j+1)
            
            call FindSpan(nby,q,vcheck,V_Knot,vspan)
    
            BFy = 0.0d0
            call BasisFuns(nby,vspan,vcheck,q,v_knot,BFy)

            counter = counter + 1
            den = 0.0;
            do l=0,q
                temp = 0.0d0 !zeros(1,1,3);
                do k=0,p
                    do k1=1,nds+1
                        temp(1,1,k1)  = temp(1,1,k1)  + BFx(k+1,1)*Bw(uspan-p+k, vspan-q+l, k1);
                    enddo
                    continue
                enddo
                
                den = den + temp(1,1,nds+1)*BFy(l+1,1);
                do k1=1,nds
                    coord(counter,k1) = coord(counter,k1) + temp(1,1,k1)*BFy(l+1,1);
                end do
            enddo
            
            do k1=1,nds
                coord(counter,k1) = coord(counter,k1)/den;
            end do
            
            continue
            
        endif
    end do
    
    !--------------------------------------------------------------------------------
    !Search u_knot
    !counter = ncpy
    do i=1,ncpx+p
        if(u_knot(i) .ne. u_knot(i+1))then
            ucheck = u_knot(i+1)
            
            call FindSpan(nbx,p,ucheck,U_Knot,uspan)
            
            BFx = 0.0d0
            call BasisFuns(nbx,uspan,ucheck,p,u_knot,BFx)
            
            
            vcheck = 0.0d0
                    
            call FindSpan(nby,q,vcheck,V_Knot,vspan)
    
            BFy = 0.0d0
            call BasisFuns(nby,vspan,vcheck,q,v_knot,BFy)
            
            counter = counter + 1
            den = 0.0;
            j = 1
            do l=0,q
                temp = 0.0d0
                do k=0,p
                    do k1=1,nds+1
                        temp(1,1,k1)  = temp(1,1,k1)  + BFx(k+1,1)*Bw(uspan-p+k, vspan-q+l, k1)
                    enddo
                    continue
                enddo
                
                den = den + temp(1,1,nds+1)*BFy(l+1,1);
                do k1=1,nds
                    coord(counter,k1) = coord(counter,k1) + temp(1,1,k1)*BFy(l+1,1);
                enddo
            enddo
            
            do k1=1,nds
                coord(counter,k1) = coord(counter,k1)/den;
            enddo
            
            do j=1,ncpy+q
                if(v_knot(j) .ne. v_knot(j+1))then
                    vcheck = v_knot(j+1)
                    
                    call FindSpan(nby,q,vcheck,V_Knot,vspan)
            
                    BFy = 0.0d0
                    call BasisFuns(nby,vspan,vcheck,q,v_knot,BFy)

                    counter = counter + 1
                    den = 0.0;
                    do l=0,q
                        temp = 0.0d0 !zeros(1,1,3);
                        do k=0,p
                            do k1=1,nds+1
                                temp(1,1,k1)  = temp(1,1,k1)  + BFx(k+1,1)*Bw(uspan-p+k, vspan-q+l, k1);
                            enddo
                            continue
                        enddo
                        
                        den = den + temp(1,1,nds+1)*BFy(l+1,1);
                        do k1=1,nds
                            coord(counter,k1) = coord(counter,k1) + temp(1,1,k1)*BFy(l+1,1);
                        end do
                    enddo
                    
                    do k1=1,nds
                        coord(counter,k1) = coord(counter,k1)/den;
                    end do
                    
                    continue
                    
                endif
            end do
            
            continue
        end if

    end do
    !--------------------------------------------------------------------------------
    
    do j=1,ncpx*ncpy
        write(2,FMT=112)j,coord(j,:)
        112 format(I4,5(E,1x))
    end do

    write(2,*)'end coordinates'
    
    !Elements ------------------------------
    count = 0
    write(2,*)'Elements'
    
    do i=1,nelx
        do j = 1, nely
        count = count + 1
        idx = i + nely*(i-1) + (j-1)
        el1 = idx
        el2 = idx + 1
        el3 = idx + nely + 2
        el4 = idx + nely + 1
        
        write(2,FMT=113)count,el1,el2,el3,el4
        113 format(21(I4,1x))
        end do
    end do
    
    write(2,*)'end elements'
    
    close(2,status='keep')
    
    continue
    
    
end subroutine