!-------------------------------------------------------------------
!
! Subroutine to read the input file
!
!-------------------------------------------------------------------

Subroutine ReadInFile(FileName)

    use Mod_Variables
    
    implicit none
    
    character*256,intent(IN)::FileName
    character*256::Line
    
    integer(4)::i,j,k,l,k1,k2,nnod,iptch
    integer(4)::int1,int2,int3,nload
    
    integer(4)::count1,count2,nnodesmax
    
    logical::pstrain, pstress
    !character::elcode
    
    real(8)::temp1
    real(8),allocatable,dimension(:,:)::matD
    
    integer(4),allocatable,dimension(:,:)::Grid,Mat
    
    integer(4)::pmax,qmax,ncpxmax,ncpymax
    
    integer(4),dimension(:),allocatable::sizex,sizey
    integer(4)::sumx,sumy
    
    open(unit=1,file=FileName)
    
    !Logical allocation
    pstrain = .false.
    pstress = .false.
    contactPTS = .false.
    contactGPTS = .false.
    Multistep = .false.
    
    nstp = 1
    
    do while(Line .ne. '*end')
        read(1,*)Line
        
        continue
        
        if(line=='*begin')then
            read(1,*)nds                    !Number of spatial directions
            if(nds==2) then
                
                nptch = 1
                
                read(1,*)p,q         !Degrees of curves in u,v,w
                read(1,*)ncpx,ncpy     !Number of control points in each direction
                read(1,*)closed_u, closed_v !Cloded knot vectors - 1: Yes, 0: No:
            
                nnodes = ncpx*ncpy          !Total number of nodes
                nelems = (ncpx-p)*(ncpy-q)  !Total number of elements
                telems = nelems
                nshpl = (p+1)*(q+1)         !Number of local shape functions
                
                tnodes = nnodes
                
                allocate(u_knot(ncpx+p+1))
                u_knot = 0
                allocate(v_knot(ncpy+q+1))
                v_knot = 0
                
                allocate(b_net(ncpx,ncpy,nds+1))
                b_net = 0.0d0
                
                allocate(b_net_final(ncpx,ncpy,nds+1))
                b_net_final = 0.0d0
                
                allocate(bw(ncpx,ncpy,nds+1))
                bw = 0.0d0
      
                allocate(iper(nnodes))
                do k1 = 1, nnodes
                    iper(k1) = k1
                enddo
                
                allocate(bc(ncpx,ncpy,nds),load(ncpx,ncpy,nds))
                bc = 1
                load = 0.0d0
                
                allocate(dispBC(ncpx,ncpy,nds))
                dispBC = 0.0d0

                allocate(ld_el(nelems,nds))
                ld_el = 0.0d0
                allocate(ey(nelems))
                ey = 0.0d0
                allocate(nu(nelems))
                nu = 0.0d0
                
                !allocate(Ke((p+1)*(q+1)*nds,(p+1)*(q+1)*nds)
                allocate(Kf(ncpx*ncpy*nds,ncpx*ncpy*nds))
                Kf = 0.0d0
                
                allocate(KCont(ncpx*ncpy*nds,ncpx*ncpy*nds))
                KCont = 0.0d0
                
                allocate(Fext(nnodes*nds,1))
                Fext = 0.0d0
                
                allocate(ImpDisp(nnodes*nds,1),StpFrac(nnodes*nds))
                ImpDisp = 0.0d0
                StpFrac = 0.0d0
                
                allocate(Fini(nnodes*nds,1))
                Fini = 0.0d0
                
                allocate(FInc(nnodes*nds,1))
                FInc = 0.0d0
                
                allocate(Res(nnodes*nds,1))
                Res = 0.0d0
                
                allocate(Fint(nnodes*nds,1))
                Fint = 0.0d0
                
                allocate(RedCount(nnodes*nds))
                redcount = 1
                
                allocate(PMesh(ncpx*ncpy,nds))
                PMesh = 0.0d0
                
                allocate(Points(ncpx*ncpy,nds))
                allocate(Points0(ncpx*ncpy,nds))
                allocate(Weights(ncpx*ncpy))
                
                allocate(dispdof(tnodes*nds))
                dispdof = 0.0d0
                
                allocate(dispdofold(tnodes*nds))
                dispdofold = 0.0d0
                
                allocate(loaddof(tnodes*nds))
                loaddof = 0.0d0
                
                allocate(distload(telems,4))
                distload = 0.0d0
                
                allocate(bcdof(tnodes*nds))
                bcdof = 1
                
                Weights = 0.0d0
                Points = 0.0d0
                Points0 = 0.0d0
                
                ndi = nds
                nshr = 1
                ntens = ndi + nshr
                                
            endif !nds=2
        
        elseif(line=='*begin_MP')then
            read(1,*)nptch
            read(1,*)nds                    !Number of spatial directions
            
            allocate(order(nptch,nds),ncpptch(nptch,nds),statptch(nptch,nds))
            
            !Read order of each patch
            pmax = 0
            qmax = 0
            do i=1,nptch
                read(1,*)order(i,1),order(i,2)
                if(order(i,1) .gt. pmax) pmax = order(i,1)
                if(order(i,2) .gt. qmax) qmax = order(i,2)
            end do
            
            !Read number of control points in each direction
            ncpxmax = 0
            ncpymax = 0
            do i=1,nptch
                read(1,*)ncpptch(i,1),ncpptch(i,2)  
                if(ncpptch(i,1) .gt. ncpxmax) ncpxmax = ncpptch(i,1)
                if(ncpptch(i,2) .gt. ncpymax) ncpymax = ncpptch(i,2)
            end do
            
            !Open-close status
            do i=1,nptch
                read(1,*)statptch(i,1),statptch(i,2)      
            end do
            
            read(1,*)tnodes
            
            !nnodes = ncpx*ncpy          !Total number of nodes
            !nelems = (ncpx-p)*(ncpy-q)  !Total number of elements
            !nshpl = (p+1)*(q+1)         !Number of local shape functions
            
            allocate(MP_u_knot(nptch,ncpxmax+pmax+1))
            MP_u_knot = 0
            allocate(MP_v_knot(nptch,ncpymax+qmax+1))
            MP_v_knot = 0
            
            allocate(MP_b_net(nptch,ncpxmax,ncpymax,nds+1))
            MP_b_net = 0.0d0
            
            allocate(MP_b_net_final(nptch,ncpxmax,ncpymax,nds+1))
            MP_b_net_final = 0.0d0
            
            allocate(dispdof(tnodes*nds))
            dispdof = 0.0d0
           
            !allocate(Ke((p+1)*(q+1)*nds,(p+1)*(q+1)*nds)
            allocate(Kf(tnodes*nds,tnodes*nds))
            Kf = 0.0d0
            
            allocate(Kcont(tnodes*nds,tnodes*nds))
            Kcont = 0.0d0
            
            allocate(Fext(tnodes*nds,1))
            Fext = 0.0d0
            
            allocate(ImpDisp(tnodes*nds,1),StpFrac(tnodes*nds))
            ImpDisp = 0.0d0
            StpFrac = 0.0d0
!            
            allocate(Fini(tnodes*nds,1))
            Fini = 0.0d0
            
            allocate(FInc(tnodes*nds,1))
            FInc = 0.0d0
            
            allocate(Res(tnodes*nds,1))
            Res = 0.0d0
            
            allocate(Fint(tnodes*nds,1))
            Fint = 0.0d0
            
            allocate(RedCount(tnodes*nds))
            redcount = 1
            
            allocate(loaddof(tnodes*nds))
            loaddof = 0.0d0
            
            allocate(GCoords(tnodes,nds+1))
            GCoords = 0.0d0
            
            allocate(bcdof(tnodes*nds))
            bcdof = 1
            
            allocate(dispdofold(tnodes*nds))
            dispdofold = 0.0d0
!            
!            allocate(PMesh(ncpx*ncpy,nds))
!            PMesh = 0.0d0
!            
!!!            allocate(Points(ncpx*ncpy,nds))
!!!            allocate(Weights(ncpx*ncpy))
!!!            
!!!            Weights = 0.0d0
!!!            Points = 0.0d0
            
            ndi = nds
            nshr = 1
            ntens = ndi + nshr
        
        elseif(line=='*knots')then
            if(nptch==1)then
                read (1,*) (u_knot(i),i=1,ncpx+p+1)
                read (1,*) (v_knot(i),i=1,ncpy+q+1)
            else
                do i=1,nptch
                     read (1,*) (MP_u_knot(i,j),j=1,ncpptch(i,1)+order(i,1)+1)
                     read (1,*) (MP_v_knot(i,j),j=1,ncpptch(i,2)+order(i,2)+1)
                end do
            end if
            
            continue
            
        elseif(line=='*element')then
            if (nptch==1)then
                read(1,*)line
                if(line == 'Quad4E')then
                    npi = 4
                    pstrain = .true.
                    elcode = 'Quad4E'
                elseif(line == 'Quad4S')then
                    npi = 4
                    pstress = .true.
                    elcode = 'Quad4S'
                elseif(line == 'Quad9E')then
                    npi = 9
                    pstrain = .true.
                    elcode = 'Quad9E'
                elseif(line == 'Quad9EBBar')then
                    npi = 9
                    pstrain = .true.
                    elcode = 'Quad9EBBar'
                elseif(line == 'Quad9S')then
                    npi = 9
                    pstress = .true.
                    elcode = 'Quad9S'
                elseif(line == 'Quad16E')then
                    npi = 16
                    pstrain = .true.
                    elcode = 'Quad16E'
                elseif(line == 'Quad16S')then
                    npi = 16
                    pstress = .true.
                    elcode = 'Quad16S'    
                else
                    write(*,*)'ERROR! - Element key does not exit'
                endif
            else
                allocate(MP_elcode(nptch))
                allocate(MP_npi(nptch))
                do i=1,nptch
                    read(1,*)line
                    if(line == 'Quad4E')then
                        MP_npi(i) = 4
                        pstrain = .true.
                        MP_elcode(i) = 'Quad4E'
                    elseif(line == 'Quad4S')then
                        MP_npi(i) = 4
                        pstress = .true.
                        MP_elcode(i) = 'Quad4S'
                    elseif(line == 'Quad9E')then
                        MP_npi(i) = 9
                        pstrain = .true.
                        MP_elcode(i) = 'Quad9E'
                    elseif(line == 'Quad9EBBar')then
                        MP_npi(i) = 9
                        pstrain = .true.
                        MP_elcode(i) = 'Quad9EBBar'
                    elseif(line == 'Quad9S')then
                        MP_npi(i) = 9
                        pstress = .true.
                        MP_elcode(i) = 'Quad9S'
                    elseif(line == 'Quad16E')then
                        MP_npi(i) = 16
                        pstrain = .true.
                        MP_elcode(i) = 'Quad16E'
                    elseif(line == 'Quad16S')then
                        MP_npi(i) = 16
                        pstress = .true.
                        MP_elcode(i) = 'Quad16S'    
                    else
                        write(*,*)'ERROR! - Element key does not exit'
                    endif
                end do
            end if
            
            continue
            
        elseif(line=='*bnet')then
            
            if (nptch==1)then
                nnod = 0
                do i = 1,ncpx
                    do j = 1,ncpy
                        read (1,*) (b_net(i,j,l),l=1,nds+1)
                        nnod = nnod + 1
                  enddo
                enddo
                
                nnod = 0
                Points = 0.0d0
                Weights = 0.0d0
                do j = 1,ncpy
                    do i = 1,ncpx
                        nnod = nnod + 1
                        Points(nnod,1) = b_net(i,j,1)
                        Points(nnod,2) = b_net(i,j,2)
                        Weights(nnod)  = b_net(i,j,3)
                  enddo
                enddo
                
                Points0 = Points
                  
            else
                nnod = 0
                do iptch = 1, nptch
                    do i = 1,ncpptch(iptch,1)
                        do j = 1,ncpptch(iptch,2)
                            read (1,*) (MP_b_net(iptch,i,j,l),l=1,nds+1)
                            nnod = nnod + 1
                        end do
                    end do
                end do
            end if
            
            continue
            
        elseif(line == '*penalty')then
            read(1,*)npenal
            read(1,*)ipenal
            
            allocate(penalmeth(npenal,ipenal*nds))
            
            do i=1,npenal
                read (1,*) (penalmeth(i,j),j=1,nds*ipenal)
            end do
            
            continue
        
        elseif(line == '*boundary')then
            
!            allocate(dispeqv(nnodes*nds-nbc,1))
!            dispeqv = 0.0d0
            
            read(1,*)nbc
            do k1=1,nbc
                read(1,*)int1,int2,int3
                bc(int1,int2,int3) = 0
            end do
            
!            allocate(KfEq(nnodes*nds-nbc,nnodes*nds-nbc))
!            KfEq = 0.0d0
!            
!            allocate(KfEqInv(nnodes*nds-nbc,nnodes*nds-nbc))
!            KfEqInv = 0.0d0
!            
!            allocate(FextEq(nnodes*nds-nbc,1))
!            FextEq = 0.0d0
            
!            allocate(dDisp(nnodes*nds,1))
!            dDisp=0.0d0
!    
!            allocate(dddisp(nnodes*nds,1))
!            ddDisp=0.0d0
!            
!            allocate(u(nnodes*nds,1))
!            u=0.0d0
            
        elseif(line == '*load')then
            read(1,*)nload
            do k1=1,nload
                read(1,*)int1,int2,int3,temp1
                load(int1,int2,int3) = temp1
            end do
        
        elseif(line=='*material')then
            
            if(nptch==1)then
                read(1,*)iprops
                allocate(props(iprops))
                read(1,*)(props(i),i=1,iprops) 
            else
                read(1,*)iprops
                allocate(MP_props(nptch,iprops))
                
                do i=1,nptch
                     read(1,*)(MP_props(i,j),j=1,iprops) 
                end do
                
            end if
            
            continue  
        
        elseif(line=='*displacement')then
            read(1,*)ndisp
            do k1=1,ndisp
                read(1,*)int1,int2,int3,temp1
                dispBC(int1,int2,int3) = temp1
            end do
        
        elseif(line=='*bcdof') then
            read(1,*)nbc
            
!            allocate(bcdof(tnodes*nds))
!            bcdof = 1
            
            do k1=1,nbc
                read(1,*)int1,int2
                bcdof((int1-1)*nds+int2) = 0
            end do
            
!            allocate(dispeqv(nnodes*nds-nbc,1))
!            dispeqv = 0.0d0
            
            !if(nptch .gt. 1)then
!                allocate(KfEq(tnodes*nds-nbc,tnodes*nds-nbc))
!                KfEq = 0.0d0
!                
!                allocate(KfEqInv(tnodes*nds-nbc,tnodes*nds-nbc))
!                KfEqInv = 0.0d0
!                
!                allocate(FextEq(tnodes*nds-nbc,1))
!                FextEq = 0.0d0
                
!                allocate(dDisp(tnodes*nds,1))
!                dDisp=0.0d0
!        
!                allocate(dddisp(tnodes*nds,1))
!                ddDisp=0.0d0
!                
!                allocate(u(tnodes*nds,1))
!                u=0.0d0
                
!                allocate(dispeqv(tnodes*nds-nbc,1))
!                dispeqv = 0.0d0
            !end if
            
            continue
            
        elseif(line == '*loaddof')then
            read(1,*)nload
            
            do k1=1,nload
                read(1,*)int1,int2,temp1
                loaddof((int1-1)*nds+int2) = temp1
            end do
            
            continue
        
        elseif(line=='*dispdof')then
            read(1,*)ndisp
            
            do k1=1,ndisp
                read(1,*)int1,int2,temp1
                dispdof((int1-1)*nds+int2) = temp1
            end do
            
        elseif(nptch .gt.1 .and. line=='*MP_conn')then
            
            allocate(MP_nnodes(nptch))
            allocate(MP_nelems(nptch))
            
            telems = 0
            npi = 0
            nnodesmax = 0
            do iptch=1,nptch
                p = order(iptch,1)
                q = order(iptch,2)
                ncpx = ncpptch(iptch,1)
                ncpy = ncpptch(iptch,2)
                
                !MP_nnodes(iptch) = ncpx*ncpy
                !if(MP_nnodes(iptch) .gt. nnodesmax) nnodesmax = MP_nnodes(iptch)
                
                MP_nnodes(iptch) = (p+1)*(q+1)
                !if(MP_nnodes(iptch) .gt. (p+1)*(q+1)) nnodesmax = (p+1)*(q+1)
                nnodesmax = (p+1)*(q+1)
                
                count1 = 0
                do i=1,ncpx+p
                    do j=1,ncpy+q
                        if((MP_u_knot(iptch,i) .ne. MP_u_knot(iptch,i+1)) .and. (MP_v_knot(iptch,j) .ne. MP_v_knot(iptch,j+1)))then
                            count1 = count1 + 1
                        end if
                    end do
                end do
                
                MP_nelems(iptch) = count1
                
                telems = telems + MP_nelems(iptch)
                
                if(MP_npi(iptch) .gt. npi) npi = MP_npi(iptch)
            end do
            
            allocate(MP_conn(telems,nnodesmax))
            
            count1 = 1
            count2 = 1
            do i=1,telems
                
                read(1,*)int1,(MP_conn(i,j),j=1,MP_nnodes(count2)) 
                
                count1 = count1 + 1
                
                if(MP_nelems(count2) .lt. count1) then
                    count1 = 1
                    count2 = count2 + 1
                end if
                
                continue
                
            end do
            
            continue
        
!        elseif(nptch .gt.1 .and. line=='*MP_Grid')then
!            
!            read(1,*)int1,int2
!            allocate(Grid(int1,int2))
!            do i=1,int1
!                read(1,*)(Grid(i,j),j=1,int2)
!            end do
!            
!            allocate(sizex(int1+1),sizey(int2+1))
!            
!            sizey = 0
!            sumy = 0
!            do i=1,int2
!                sizey(i+1) = sizey(i) + ncpptch(Grid(1,i),1)
!                sumy = sumy + ncpptch(Grid(1,i),1)
!            end do
!            
!            sizex = 0
!            sumx = 0
!            do i=1,int1
!                sizex(i+1) = sizex(i) + ncpptch(Grid(i,1),2)
!                sumx = sumx + ncpptch(Grid(i,1),2)
!            end do
!            
!            allocate(Mat(sumx,sumy))
!            Mat = 0
!            
!            do i=1,int1
!                do j=1,int2
!                    Mat(sizex(i)+1:sizex(i+1),sizey(j)+1:sizey(j+1)) = Grid(i,j)
!                end do
!            end do
!            
!            do i=1,int1
!                do j=1,int2
!                    Mat(sizex(i)+1:sizex(i+1),sizey(j)+1:sizey(j+1)) = Grid(i,j)
!                end do
!            end do
!            
!            
!            continue
        
        
        elseif(line=='*ContactPTS')then
            ContactPTS = .true.
            
            do i=1,2
                read(1,*)line
                if(line=='*Master')then
                    
                    read(1,*)p_m
                    read(1,*)ncp_m
                    
                    allocate(u_knot_master(ncp_m+p_m+1))
                    u_knot_master = 0
                    
                    read (1,*) (u_knot_master(j),j=1,ncp_m+p_m+1)
                    
                    allocate(b_master(ncp_m,nds+1))
                    b_master = 0.0d0
                    
                    allocate(conn_master(ncp_m))
                    conn_master = 0
                    
                    read (1,*) (conn_master(j),j=1,ncp_m)
                    
                    continue
                    
                elseif(line=='*Slave')then
                    
                    !read(1,*)islv
                    read(1,*)p_s
                    read(1,*)ncp_s
                    
                    allocate(u_knot_slave(ncp_s+p_s+1))
                    u_knot_slave = 0
                    
                    read (1,*) (u_knot_slave(j),j=1,ncp_s+p_s+1)
                    
                    allocate(b_slave(ncp_s,nds+1))
                    b_slave = 0.0d0
                    
                    allocate(conn_slave(ncp_s))
                    conn_slave = 0
                    
                    read (1,*) (conn_slave(j),j=1,ncp_s)
                    
                    open(unit=113,file='ContactData.txt')
                    write(113,*)p_s
                    write(113,*)ncp_s
                    write(113,FMT=58)u_knot_slave
                    write(113,FMT=59)conn_slave
                    close(113)
                    
                    58 format(150(E,',',1x))
                    59 format(150(I4,',',1x))
                    
                    continue
                    
                end if
            end do
        
        elseif(line=='*ContactGPTS')then
            ContactGPTS = .true.
            
            read(1,*)SurfGP
            
            do i=1,2
                read(1,*)line
                if(line=='*Master')then
                    
                    read(1,*)p_m
                    read(1,*)ncp_m
                    
                    allocate(u_knot_master(ncp_m+p_m+1))
                    u_knot_master = 0
                    
                    read (1,*) (u_knot_master(j),j=1,ncp_m+p_m+1)
                    
                    allocate(b_master(ncp_m,nds+1))
                    b_master = 0.0d0
                    
                    allocate(conn_master(ncp_m))
                    conn_master = 0
                    
                    read (1,*) (conn_master(j),j=1,ncp_m)
                    
                    continue
                    
                elseif(line=='*Slave')then
                    
                    !read(1,*)islv
                    read(1,*)p_s
                    read(1,*)ncp_s
                    
                    allocate(u_knot_slave(ncp_s+p_s+1))
                    u_knot_slave = 0
                    
                    read (1,*) (u_knot_slave(j),j=1,ncp_s+p_s+1)
                    
                    allocate(b_slave(ncp_s,nds+1))
                    b_slave = 0.0d0
                    
                    allocate(conn_slave(ncp_s))
                    conn_slave = 0
                    
                    read (1,*) (conn_slave(j),j=1,ncp_s)
                    
                    open(unit=113,file='ContactData.txt')
                    write(113,*)p_s
                    write(113,*)ncp_s
                    write(113,FMT=60)u_knot_slave
                    write(113,FMT=61)conn_slave
                    close(113)
                    
                    60 format(150(E,',',1x))
                    61 format(150(I4,',',1x))
                    
                    continue
                    
                end if
            end do
            
            continue
        
        elseif(line=='*nlgeom')then
            nlgeom = .true.
            
        elseif(line=='*Multistep')then
            
            Multistep = .true.
            
            read(1,*)nstp  
            
            allocate(STP_nbc(nstp),STP_ndisp(nstp),STP_nload(nstp),STP_dload(nstp))
            allocate(STP_itermax(nstp),STP_incmax(nstp))
            STP_nbc = 0
            STP_ndisp = 0
            STP_nload = 0
            STP_dload = 0
            STP_itermax = 10
            STP_incmax = 1
            
            allocate(STP_bcdof(nstp,tnodes*nds))
            STP_bcdof = 1
            
            allocate(STP_dispdof(nstp,tnodes*nds))
            STP_dispdof = 0.0d0
            
            allocate(STP_distload(nstp,telems,4))
            STP_distload = 0.0d0
            
            allocate(STP_loaddof(nstp,tnodes*nds))
            STP_loaddof = 0.0d0
            
            i=1
            do while(i <= nstp)
                read(1,*)line
                
                if(line=='*endstep') then
                    i = i + 1
                
                elseif(line=='*bcdof')then
                    
                    read(1,*)STP_nbc(i)
                    
                    do k1=1,STP_nbc(i)
                        read(1,*)int1,int2
                        STP_bcdof(i,(int1-1)*nds+int2) = 0
                    end do
                
                elseif(line=='*dispdof')then  
                    
                    read(1,*)STP_ndisp(i)
            
                    do k1=1,STP_ndisp(i)
                        read(1,*)int1,int2,temp1
                        STP_dispdof(i,(int1-1)*nds+int2) = temp1
                    end do    
                
                elseif(line == '*loaddof')then
                    
                    read(1,*)STP_nload(i)
                    
                    do k1=1,STP_nload(i)
                        read(1,*)int1,int2,temp1
                        STP_loaddof(i,(int1-1)*nds+int2) = temp1
                    end do
                
                elseif(line == '*distload')then
                    read(1,*)STP_dload(i)
                    
                    do k1=1,STP_dload(i)
                        read(1,*)int1,int2,temp1
                        STP_distload(i,int1,int2) = temp1
                    end do    
                    
                    continue
                    
                elseif(line=='*increment')then
                    read(1,*)STP_incmax(i)
                    
                elseif(line=='*iteration')then
                    read(1,*)STP_itermax(i)
                    
                    continue
                    
!                    if(nptch .gt. 1)then
!                        allocate(KfEq(tnodes*nds-nbc,tnodes*nds-nbc))
!                        KfEq = 0.0d0
!                        
!                        allocate(KfEqInv(tnodes*nds-nbc,tnodes*nds-nbc))
!                        KfEqInv = 0.0d0
!                        
!                        allocate(FextEq(tnodes*nds-nbc,1))
!                        FextEq = 0.0d0
!                        
!                        allocate(dDisp(tnodes*nds,1))
!                        dDisp=0.0d0
!                
!                        allocate(dddisp(tnodes*nds,1))
!                        ddDisp=0.0d0
!                        
!                        allocate(u(tnodes*nds,1))
!                        u=0.0d0
!                        
!                        allocate(dispeqv(tnodes*nds-nbc,1))
!                        dispeqv = 0.0d0
!                    end if

                
                end if
                
               continue 
                
            end do
            
            continue
            
        endif !line
    
    end do
    
    allocate(KfEq(tnodes*nds-nbc,tnodes*nds-nbc))
    KfEq = 0.0d0
    
    allocate(KfEqInv(tnodes*nds-nbc,tnodes*nds-nbc))
    KfEqInv = 0.0d0
    
    allocate(FextEq(tnodes*nds-nbc,1))
    FextEq = 0.0d0
    
    allocate(dispeqv(tnodes*nds-nbc,1))
    dispeqv = 0.0d0
    
    if(nptch .gt. 1)then
!        telems = 0
!        npi = 0
!        do iptch=1,nptch
!            p = order(iptch,1)
!            q = order(iptch,2)
!            ncpx = ncpptch(iptch,1)
!            ncpy = ncpptch(iptch,2)
!
!            nelems = (ncpx-p)*(ncpy-q)
!            
!            telems = telems + nelems
!            
!            if(MP_npi(iptch) .gt. npi) npi = MP_npi(iptch)
!        end do

        
        allocate(stress(telems,npi,ntens),dstress(telems,npi,ntens))
        stress = 0.0d0
        dstress = 0.0d0
        
        allocate(strain(telems,npi,ntens),dstrain(telems,npi,ntens))
        strain = 0.0d0
        dstrain = 0.0d0
        
        allocate(stress_conv(telems,npi,ntens),strain_conv(telems,npi,ntens))
        stress_conv = 0.0d0
        strain_conv = 0.0d0
        
        allocate(hard(telems,npi),hard_conv(telems,npi))
        hard = 0.0d0
        hard_conv = 0.0d0
        
        allocate(laxis(telems,npi,nds,nds),laxis_conv(telems,npi,nds,nds))
        laxis = 0.0d0
        laxis_conv = 0.0d0
        
        allocate(FintLM(tnodes*nds,1),FintLMConv(tnodes*nds,1))
        FintLM = 0.0d0
        
        allocate(distload(telems,4))
        distload = 0.0d0
        
        continue
        
    else
        allocate(stress(nelems,npi,ntens),dstress(nelems,npi,ntens))
        stress = 0.0d0
        dstress = 0.0d0
        
        allocate(strain(nelems,npi,ntens),dstrain(nelems,npi,ntens))
        strain = 0.0d0
        dstrain = 0.0d0
        
        allocate(stress_conv(nelems,npi,ntens),strain_conv(nelems,npi,ntens))
        stress_conv = 0.0d0
        strain_conv = 0.0d0
        
        allocate(hard(nelems,npi),hard_conv(nelems,npi))
        hard = 0.0d0
        hard_conv = 0.0d0
        
        allocate(TmatD(nelems,npi,ntens,ntens))
        allocate(matD(ntens,ntens))
        TmatD = 0.0d0
        
        allocate(laxis(nelems,npi,nds,nds),laxis_conv(nelems,npi,nds,nds))
        laxis = 0.0d0
        laxis_conv = 0.0d0
        
!        allocate(distload(telems,4))
!        distload = 0.0d0
        
        if (pstrain==.true.) call MatLinearElasticDP(props(1),props(2),MatD)
        if (pstress==.true.) call MatLinearElasticTP(props(1),props(2),MatD)
        
        do k1=1,nelems
            do k2=1,npi
                TmatD(k1,k2,:,:) = matD(:,:)
            end do
        end do
    end if
    
    allocate(GP_coords(telems*npi,nds))
    GP_coords = 0.0d0
    
    close(1,status='keep')
    
    continue
    
end