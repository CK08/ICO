subroutine GiDEqRes(inc,coordi)
    
    use Mod_Variables
    implicit none
    
    integer(4),intent(IN)::inc
    integer(4)::k1,count,i,j,nelx,nely,nbx,nby,l,k,counter
    character(16)::my_type
    real(8),dimension(ncpx*ncpy,nds),intent(IN)::coordi
    real(8),dimension(ncpx*ncpy,nds)::coord,du
    real(8),dimension(q+1)::Ny
    real(8),dimension(q+1)::dNydxi
    real(8)::ucheck,vcheck,den
    integer(4)::uspan,vspan,idx
    real(8),dimension(p+1,1)::BFx
    real(8),dimension(q+1,1)::BFy
    real(8),dimension(1,1,nds+1)::temp
    
    
    !To Do: verificar a ordem dos pontos dos GP para o caso do quadrilatero 3x3
    
    if(elcode == 'Quad4E' .or. elcode == 'Quad4S') then
        my_type = "Quadrilateral"
        npi = 4 
    elseif(elcode == 'Quad9E' .or. elcode == 'Quad9S') then
        my_type = "Quadrilateral"
        npi = 9 
    endif
    
    open(unit=3,file='IGAEq.flavia.res')
    
    !Heading ---
    write(3,FMT=99)
    99 format('GiD Post Results File 1.0')
    write(3,*)''
    
    !Final polygon shape 
    count = 0
    do j=1,ncpy
        do i=1,ncpx
            count = count+1
            B_net_final(i,j,1:nds) = PMesh(count,:)
            B_net_final(i,j,nds+1) = B_net(i,j,nds+1)
        end do
    end do
    
    
    nbx = ncpx+p+1
    nby = ncpy+q+1
    
    nelx = ncpx - p
    nely = ncpy - q
    
    coord = 0.0d0

    !--------------------------------------------------------------------------------
    !Determine position of the equivalent nodes in the deformed configuration
    !--------------------------------------------------------------------------------
    !u = 0.0
    counter = 0
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
                temp(1,1,k1)  = temp(1,1,k1)  + BFx(k+1,1)*b_net_final(uspan-p+k, vspan-q+l, k1)
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
                        temp(1,1,k1)  = temp(1,1,k1)  + BFx(k+1,1)*b_net_final(uspan-p+k, vspan-q+l, k1);
                    enddo
                    continue
                enddo
                
                den = den + temp(1,1,nds+1)*BFy(l+1,1);
                do k1=1,nds
                    coord(counter,k1) = coord(counter,k1) + temp(1,1,k1)*BFy(l+1,1);
                end do
            enddo
            
            do k1=1,nds
                coord(counter,k1) = coord(counter,k1) !/den;
            end do
            
            continue
            
        endif
    end do
    
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
                        temp(1,1,k1)  = temp(1,1,k1)  + BFx(k+1,1)*b_net_final(uspan-p+k, vspan-q+l, k1)
                    enddo
                    continue
                enddo
                
                den = den + temp(1,1,nds+1)*BFy(l+1,1);
                do k1=1,nds
                    coord(counter,k1) = coord(counter,k1) + temp(1,1,k1)*BFy(l+1,1);
                enddo
            enddo
            
            do k1=1,nds
                coord(counter,k1) = coord(counter,k1) !/den;
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
                                temp(1,1,k1)  = temp(1,1,k1)  + BFx(k+1,1)*b_net_final(uspan-p+k, vspan-q+l, k1);
                            enddo
                            continue
                        enddo
                        
                        den = den + temp(1,1,nds+1)*BFy(l+1,1);
                        do k1=1,nds
                            coord(counter,k1) = coord(counter,k1) + temp(1,1,k1)*BFy(l+1,1);
                        end do
                    enddo
                    
                    do k1=1,nds
                        coord(counter,k1) = coord(counter,k1) !/den;
                    end do
                    
                    continue
                    
                endif
            end do
            
            continue
        end if

    end do
    
    !--------------------------------------------------------------------------------
    !Displacement increment of the equivalent nodes
    !--------------------------------------------------------------------------------
    du = 0.0d0
    du = coord-coordi
    
    
    
    
    !Gauss Point Data -----------
    write(3,FMT=100)my_type
    100 format('GaussPoints "GP" ElemType ', A13,' "problem"')
    write(3,FMT=101)npi
	101 format('    Number Of Gauss Points: ', I2)
	write(3,FMT=102)
	102 format('    Natural Coordinates: internal')
    write(3,FMT=103) 
    103 format('End GaussPoints')
    write(3,*)''
    
    !Displacement data ----
    write(3,FMT=115)inc
    115 format('Result "Displacements" "Load Analysis"', I3, ' Vector OnNodes')
    if (nds==2) write(3,FMT=116)
    116 format('ComponentNames "X-Disp", "Y-Disp"')
    !if (nds==3) write(13,*)'ComponentNames "X-Disp", "Y-Disp", "Z-Disp"'
    write(3,*)'Values'
    
    do j=1,ncpx*ncpy
        write(3,FMT=117)j, du(j,:)
        117 format(I4,5(E,1x))
    end do
    
    write(3,FMT=118)
    118 format('End Values')
    write(3,*)''
    
    
    
    
    !Stress data
    write(3,FMT=104)inc
    104 format('Result "Stress" "Load Analysis"', I3, ' PlainDeformationMatrix OnGaussPoints "GP"')
    write(3,FMT=105)
    105 format('ComponentNames "Sx", "Sy", "Sxy", "Szz"')
    write(3,FMT=106)
    106 format('Values')
    
    do i=1,nelems
        write(3,FMT=107)i, stress(i,1,1), stress(i,1,2), stress(i,1,3), 0.0d0
        107 format(I6,6(E,1x))
        do j=2,npi
            write(3,FMT=108)stress(i,j,1), stress(i,j,2), stress(i,j,3), 0.0d0
            108 format('      '6(E,1x))
        end do
    end do
    
    write(3,FMT=109)
    109 format('End Values')
    
    
    
    
    
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
!    write(3,*)'End Values'
    
    close(3,status='keep')
    
    continue
    
end subroutine