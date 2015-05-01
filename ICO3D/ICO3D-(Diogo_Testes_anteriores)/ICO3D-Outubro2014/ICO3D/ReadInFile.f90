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
    
    integer(4)::i,j,k,l,nnod,count,counter,ipair
    integer(4)::int1,int2,int3,int4,nload
    
    logical::pstrain, pstress
    !character::elcode
    
    real(8)::temp1
    integer(4)::nmax,el1
    real(8),allocatable,dimension(:,:)::matD
    
    
    
    
    integer(4)::iptch,k1,k2,k3,k4,k5,k6,k7,k8,nnod1,nnod2
    integer :: loop1,loop2,loop3,g, e, gtemp, ln
    
    open(unit=1,file=FileName)
    
    !Logical allocation
    pstrain = .false.
    pstress = .false.
    gravity = .false.
    pressure = .false.
    
    ntens = 6
    ndi = 3
    nshr = 3
    
    do while(Line .ne. '*end')
        !write(*,*)Line
        read(1,*)Line
        if(line=='*Contact_PTS' .or. line=='*MultiBody_Contact_PTS')then
            PTS = .true.
        end if
    end do
    
    rewind(unit=1)
    read(1,*)Line
    close(1)
    
    
    open(unit=1,file=FileName)
    
    
    do while(Line .ne. '*end')
        !write(*,*)Line
        read(1,*)Line
        
        if(line=='*begin')then
            
            nptch = 1
            
            read(1,*)nds                    !Number of spatial directions
            read(1,*)p,q,w         !Degrees of curves in u,v,w
            read(1,*)ncpx,ncpy,ncpz     !Number of control points in each direction
            read(1,*)closed_u, closed_v, closed_w !Cloded knot vectors - 1: Yes, 0: No:
        
            nnodes = ncpx*ncpy*ncpz              !Total number of nodes
            nelems = (ncpx-p)*(ncpy-q)*(ncpz-w)  !Total number of elements
            
            nshpl = (p+1)*(q+1)*(w+1)            !Number of local shape functions
            
            tnodes = nnodes
            
            call AllocGlobalVars()
            
            ndi = nds
            nshr = 3
            ntens = ndi + nshr
        
        
        elseif(line=='*beginMP')then
            read(1,*)nptch
            read(1,*)nds
            
            allocate(MP_nds(nptch))
            allocate(MP_p(nptch),MP_q(nptch),MP_w(nptch))
            allocate(MP_ncpx(nptch),MP_ncpy(nptch),MP_ncpz(nptch))
            allocate(MP_closed_u(nptch),MP_closed_v(nptch),MP_closed_w(nptch))
            
            pmax = 0
            qmax = 0
            wmax = 0
            do i=1,nptch
                read(1,*)MP_p(i),MP_q(i),MP_w (i)
                if(MP_p(i) .gt. pmax) pmax = MP_p(i)
                if(MP_q(i) .gt. qmax) qmax = MP_q(i)
                if(MP_w(i) .gt. wmax) wmax = MP_w(i)
            end do
            
            ncpxmax = 0
            ncpymax = 0
            ncpzmax = 0
            do i=1,nptch
                read(1,*)MP_ncpx(i),MP_ncpy(i),MP_ncpz(i)
                if(MP_ncpx(i) .gt. ncpxmax) ncpxmax = MP_ncpx(i)
                if(MP_ncpy(i) .gt. ncpymax) ncpymax = MP_ncpy(i)
                if(MP_ncpz(i) .gt. ncpzmax) ncpzmax = MP_ncpz(i)
            end do
            
            do i=1,nptch
                read(1,*)MP_closed_u(i), MP_closed_v(i), MP_closed_w (i)
            end do
            
            nnodes = 0
            do i=1,nptch
                nnodes = nnodes + MP_ncpx(i)*MP_ncpy(i)*MP_ncpz(i)
            end do
            
            nelems = 0
            do i=1,nptch
                nelems = nelems + (MP_ncpx(i)-MP_p(i))*(MP_ncpy(i)-MP_q(i))*(MP_ncpz(i)-MP_w(i))
            end do
            
            nshpl = 0
            do i=1,nptch
                nshpl = nshpl + (MP_p(i)+1)*(MP_q(i)+1)*(MP_w(i)+1)
            end do                    
            
            continue
            
        elseif(line=='*knots')then
            
            if(nptch==1)then
                read (1,*) (u_knot(i),i=1,ncpx+p+1)
                read (1,*) (v_knot(i),i=1,ncpy+q+1)
                read (1,*) (w_knot(i),i=1,ncpz+w+1)
                
                nelems = 0
                do k = 1,ncpz
                    do j = 1,ncpy           
                        do i = 1,ncpx
                            if((u_knot(i) .ne. u_knot(i+1)) .and. &
                            &  (v_knot(j) .ne. v_knot(j+1)) .and. &
                            &  (w_knot(k) .ne. w_knot(k+1))) then  
                                nelems = nelems + 1     
                            end if
                        enddo
                    enddo
                enddo
                
            else
                
                allocate(MP_uknot(nptch,ncpxmax+pmax+1))
                allocate(MP_vknot(nptch,ncpymax+qmax+1))
                allocate(MP_wknot(nptch,ncpzmax+wmax+1))
                
                do i=1,nptch
                    read (1,*) (MP_uknot(i,j),j=1,MP_ncpx(i)+MP_p(i)+1)
                end do
                
                do i=1,nptch
                    read (1,*) (MP_vknot(i,j),j=1,MP_ncpy(i)+MP_q(i)+1)
                end do
                
                do i=1,nptch
                    read (1,*) (MP_wknot(i,j),j=1,MP_ncpz(i)+MP_w(i)+1)
                end do
                
                nelems = 0
                do iptch=1,nptch
                    do k = 1,ncpz
                        do j = 1,ncpy           
                            do i = 1,ncpx
                                if((MP_uknot(iptch,i) .ne. MP_uknot(iptch,i+1)) .and. &
                                &  (MP_vknot(iptch,j) .ne. MP_vknot(iptch,j+1)) .and. &
                                &  (MP_wknot(iptch,k) .ne. MP_wknot(iptch,k+1))) then
                                    nelems = nelems + 1     
                                end if
                            end do
                        end do
                    end do
                end do
            end if
        
        elseif(line=='*element')then

            npi = 0
            do i=1,nptch
                read(1,*)line
                if(line == 'Hex8')then
                    npi = 8
                    !pstrain = .true.
                    elcode = 'Hex8'
                elseif(line == 'Hex8BBar')then
                    npi = 8
                    !pstrain = .true.
                    elcode = 'Hex8BBar'
                elseif(line == 'Hex8ANS')then
                    npi = 8
                    !pstrain = .true.
                    elcode = 'Hex8ANS'    
                elseif(line == 'Hex27')then
                    npi = 27
                    !pstrain = .true.
                    elcode = 'Hex27'
                elseif(line == 'Hex27EAS')then
                    npi = 27
                    !pstrain = .true.
                    elcode = 'Hex27EAS'
                 elseif(line == 'Hex27SRI')then
                    npi = 27
                    !pstrain = .true.
                    elcode = 'Hex27SRI'
                elseif(line == 'Hex27BBar')then
                    npi = 27
                    !pstrain = .true.
                    elcode = 'Hex27BBar'
                elseif(line == 'Hex27EAS_PW')then
                    npi = 27
                    !pstrain = .true.
                    nalpha = 12
                    elcode = 'Hex27EAS_PW'
                elseif(line == 'Hex27ANS_PW')then
                    npi = 27
                    !pstrain = .true.
                    nalpha = 0
                    elcode = 'Hex27ANS_PW'
                elseif(line == 'Hex27ANS')then
                    npi = 27
                    !pstrain = .true.
                    elcode = 'Hex27ANS'
                    
                    !Precompute required matrices
                    call PreCompANS()
                    
                elseif(line == 'Hex27R')then
                    npi = 8
                    !pstrain = .true.
                    elcode = 'Hex27R'
                elseif(line == 'Hex27PVS')then
                    npi = 27
                    !pstrain = .true.
                    elcode = 'Hex27PVS'
                elseif(line == 'Hex27PV')then
                    npi = 27
                    !pstrain = .true.
                    elcode = 'Hex27PV'    
                elseif(line == 'Hex64')then
                    npi = 64
                    !pstrain = .true.
                    elcode = 'Hex64'
                elseif(line == 'Hex27_Teste')then
                    npi = 27
                    !pstrain = .true.
                    elcode = 'Hex27_Teste'
                elseif(line == 'HexRed')then
                    npi = 27
                    !pstrain = .true.
                    elcode = 'HexRed'
                elseif(line == 'Hex27_ProjVol')then
                    npi = 27
                    !pstrain = .true.
                    elcode = 'Hex27_ProjVol' 
                else
                    write(*,*)'ERROR! - Element key does not exit'
                endif
            
                if(npi .gt. npimax) npimax=npi
            
                if(nptch .gt. 1)then
                    if(i==1)then
                        allocate(MP_npi(nptch),MP_elcode(nptch))
                    end if
                    
                    MP_npi(i) = npi
                    MP_elcode(i) = elcode
                end if
                
                continue
                
            end do
            
            
            !Allocations for projection method ---------------------------------------
!            if(elcode=='Hex8BBar' .or. elcode=='Hex27BBar')then
!                allocate(MatAf(ncpx*ncpy*ncpz*nds,ncpx*ncpy*ncpz*nds))
!                MatAf = 0.0d0
!                allocate(MatCf(ncpx*ncpy*ncpz*nds,(ncpx-1)*(ncpy-1)*(ncpz-1)))
!                MatCf = 0.0d0
!                allocate(MatCfT((ncpx-1)*(ncpy-1)*(ncpz-1),ncpx*ncpy*ncpz*nds))
!                MatCfT = 0.0d0
!                allocate(MatVf((ncpx-1)*(ncpy-1)*(ncpz-1),(ncpx-1)*(ncpy-1)*(ncpz-1)))
!                MatVf = 0.0d0
!                
!                allocate(Kbar(ncpx*ncpy*ncpz*nds,ncpx*ncpy*ncpz*nds))
!                Kbar = 0.0d0
!                
!                allocate(fbar(ncpx*ncpy*ncpz*nds,1))
!                fbar = 0.0d0
!                
!            elseif(elcode=='Hex27EAS_PW')then
!                allocate(MatAf(ncpx*ncpy*ncpz*nds,ncpx*ncpy*ncpz*nds))
!                MatAf = 0.0d0
!                allocate(MatCf(ncpx*ncpy*ncpz*nds,nalpha*nelems))
!                MatCf = 0.0d0
!                allocate(MatCfT(nalpha*nelems,ncpx*ncpy*ncpz*nds))
!                MatCfT = 0.0d0
!                allocate(MatVf(nalpha*nelems,nalpha*nelems))
!                MatVf = 0.0d0
!                
!                allocate(diag(nalpha*nelems))
!                diag = 0.0d0
!                
!                continue
!                
!            elseif(elcode=='Hex27ANS_PW')then
!                allocate(MatAf(ncpx*ncpy*ncpz*nds,ncpx*ncpy*ncpz*nds))
!                MatAf = 0.0d0
!                allocate(MatCf(ncpx*ncpy*ncpz*nds,ncpx*ncpy*ncpz*nds))
!                MatCf = 0.0d0
!                allocate(MatCfT(ncpx*ncpy*ncpz*nds,ncpx*ncpy*ncpz*nds))
!                MatCfT = 0.0d0
!                allocate(MatVf(ncpx*ncpy*ncpz*nds,ncpx*ncpy*ncpz*nds))
!                MatVf = 0.0d0
!            endif
            
            continue
        
        elseif(line=='*bnet')then
            
            if(nptch==1)then
                nnod = 0
                do i = 1,ncpx
                    do j = 1,ncpy
                        do k=1, ncpz
                            read (1,*) (b_net(i,j,k,l),l=1,nds+1)
                            nnod = nnod + 1
                        end do
                    enddo
                enddo
              
                nnod = 0
                do k=1,ncpz
                    do j=1,ncpy
                        do i=1,ncpx
                            nnod = nnod + 1
                            weight(nnod) = b_net(i,j,k,nds+1)
                        end do
                    enddo
                enddo
              
                nnod = 0
                do k=1,ncpz
                    do j=1,ncpy
                        do i=1,ncpx 
                            nnod = nnod + 1
                            Points(nnod,1) = b_net(i,j,k,1)
                            Points(nnod,2) = b_net(i,j,k,2)
                            Points(nnod,3) = b_net(i,j,k,3)
                        end do
                    enddo
                enddo
            
            else
            
                allocate(MP_b_net(nptch,ncpxmax,ncpymax,ncpzmax,nds+1))
                MP_b_net = 0.0d0
                allocate(MP_b_net_final(nptch,ncpxmax,ncpymax,ncpzmax,nds+1))
                MP_b_net_final = 0.0d0
                
                nnod = 0
                do iptch=1,nptch
                    do i = 1,MP_ncpx(iptch)
                        do j = 1,MP_ncpy(iptch)
                            do k=1, MP_ncpz(iptch)
                                read (1,*) (MP_b_net(iptch,i,j,k,l),l=1,nds+1)
                                nnod = nnod + 1
                            end do
                        enddo
                    enddo
                end do
                
                if(PTS == .false.)then
                    !Check for coincident nodes is multiple patches
                    coin = 0
                    nnod1 = 0
                    do iptch=1,nptch-1
                        do k3 = 1, MP_ncpz(iptch)
                            do k2 = 1,MP_ncpy(iptch)
                                do k1 = 1,MP_ncpx(iptch)
                                    nnod1 = nnod1 + 1
                                    nnod2 = 0
                                 !--------------------------
                                    do k6 = 1, MP_ncpz(iptch+1)
                                        do k5 = 1,MP_ncpy(iptch+1)
                                            do k4 = 1,MP_ncpx(iptch+1)
                                              nnod2 = nnod2 + 1
                                              !--------------------------  
                                                if(MP_b_net(iptch,k1,k2,k3,1) == MP_b_net(iptch+1,k4,k5,k6,1)) then
                                                    if(MP_b_net(iptch,k1,k2,k3,2) == MP_b_net(iptch+1,k4,k5,k6,2)) then
                                                        if(MP_b_net(iptch,k1,k2,k3,3) == MP_b_net(iptch+1,k4,k5,k6,3)) then
                                                            if(MP_b_net(iptch,k1,k2,k3,4) == MP_b_net(iptch+1,k4,k5,k6,4)) then
                                                                coin = coin + 1
                                                            end if
                                                        end if
                                                    end if
                                                end if
                                              !--------------------------  
                                            end do
                                        end do
                                    end do
                                  !--------------------------
                                end do
                            enddo
                        enddo  
                    end do
                    
                    allocate(CoinCP(nptch,coin))
                    CoinCP = 0
                    
    !                allocate(GPoints(tnodes-coin,nds+1))
    !                GPoints = 0.0d0
                    
                    counter = 0
                    coin = 0
                    do iptch=1,nptch-1
                        ln = MP_ncpx(iptch)*MP_ncpy(iptch)*MP_ncpz(iptch)
                        
                        if(iptch .gt. 1)then
                            nnod1 = MP_ncpx(iptch-1)*MP_ncpy(iptch-1)*MP_ncpz(iptch-1)
                        else
                            nnod1 = 0
                        end if
                        
                        do k3 = 1, MP_ncpz(iptch)
                            do k2 = 1,MP_ncpy(iptch)
                                do k1 = 1,MP_ncpx(iptch)
                                    nnod1 = nnod1 + 1
                                    
                                    nnod2 = 0
                                    do i=1,iptch
                                        ln = MP_ncpx(i)*MP_ncpy(i)*MP_ncpz(i)
                                        nnod2 = nnod2 + ln
                                    end do
                                 !--------------------------
                                    do k6 = 1, MP_ncpz(iptch+1)
                                        do k5 = 1,MP_ncpy(iptch+1)
                                            do k4 = 1,MP_ncpx(iptch+1)
                                              nnod2 = nnod2 + 1
                                              !--------------------------  
                                                if(MP_b_net(iptch,k1,k2,k3,1) == MP_b_net(iptch+1,k4,k5,k6,1)) then
                                                    if(MP_b_net(iptch,k1,k2,k3,2) == MP_b_net(iptch+1,k4,k5,k6,2)) then
                                                        if(MP_b_net(iptch,k1,k2,k3,3) == MP_b_net(iptch+1,k4,k5,k6,3)) then
                                                            if(MP_b_net(iptch,k1,k2,k3,4) == MP_b_net(iptch+1,k4,k5,k6,4)) then
                                                                coin = coin + 1
                                                                CoinCP(iptch,coin) = nnod1
                                                                CoinCP(iptch+1,coin) = nnod2
                                                            end if
                                                        end if
                                                    end if
                                                else
    !                                                counter = counter + 1
    !                                                GPoints(counter,1) =     
                                                    
                                                        
                                                end if
                                              !--------------------------  
                                            end do
                                        end do
                                    end do
                                  !--------------------------
                                end do
                            enddo
                        enddo
                    end do
                    
                    end if
                    
                    tnodes = nnodes - coin
                    
                    !Generate global connectivities
                    allocate(MP_INN(nnodes,nds))
                    allocate(MP_IEN(nelems,(pmax+1)*(qmax+1)*(wmax+1)))
                    allocate(MP_IEN_Temp(nelems,(pmax+1)*(qmax+1)*(wmax+1)))
                    
                    MP_IEN = 0
                    MP_INN = 0
                    g   = 0
                    e   = 0
                    
                    do iptch=1,nptch
                        do k = 1,MP_ncpz(iptch)
                            do j = 1,MP_ncpy(iptch)     ! loop through control points in V direction
                                do i = 1,MP_ncpx(iptch)       ! loop through control points in U direction
                                    g = g + 1       ! increment global function number
                                    MP_INN(g,1) = i
                                    MP_INN(g,2) = j
                                    MP_INN(g,3) = k
                                    if(i .gt. MP_p(iptch) .and. j .gt. MP_q(iptch) .and. k .gt. MP_w(iptch)) then
                                        e = e + 1                       ! increment element number
                                        do loop1 = 0,MP_w(iptch)
                                            do loop2 = 0,MP_q(iptch)
                                                do loop3 = 0,MP_p(iptch)
                                                    gtemp     = g - MP_ncpx(iptch)*(loop1*MP_ncpy(iptch) + loop2) - loop3
                                                    ln        = (MP_p(iptch)+1)*((MP_q(iptch)+1)*loop1 + loop2) + loop3 + 1
                                                    MP_IEN(e,ln) = gtemp
                                                    continue
                                                enddo ! loop3
                                            enddo ! loop2
                                        enddo !loop1

                                    endif
                                enddo ! i
                            enddo ! j
                        enddo !k
                    end do
                    
                    if(coin .gt. 0)then
                        !Replace coincident control points
                        MP_IEN_Temp = MP_IEN
                        el1 = 1
                        
                        do iptch=1,nptch-1
                            do i=1,coin
                            
                                if(CoinCP(iptch,i) .gt. 0 .and. CoinCP(iptch+1,i) .gt. 0)then
                                    
                                    do j=1,nelems
                                        do k=1,(pmax+1)*(qmax+1)*(wmax+1)
                                            if(MP_IEN(j,k) == CoinCP(iptch+1,i)) then
                                                MP_IEN(j,k)=CoinCP(iptch,i)
                                                continue
                                            end if
                                        end do
                                    end do

                                end if
                            end do
                            
                            el1 = el1 + (MP_ncpx(iptch)-MP_p(iptch))*(MP_ncpy(iptch)-MP_q(iptch))*(MP_ncpz(iptch)-MP_w(iptch))
                            
                        end do
                        
                        !MP_IEN = MP_IEN_Temp
                        
                        nmax = MP_ncpx(1)*MP_ncpy(1)*MP_ncpz(1)
                        el1 = (MP_ncpx(1)-MP_p(1))*(MP_ncpy(1)-MP_q(1))*(MP_ncpz(1)-MP_w(1))
                        
                        do j=el1+1,nelems
                            do k=(pmax+1)*(qmax+1)*(wmax+1),1,-1
                                if(MP_IEN(j,k) .gt. nmax) then
                                    nmax = nmax + 1
                                    
                                    do k1=el1+1,nelems
                                        do k2=1,(pmax+1)*(qmax+1)*(wmax+1)
                                            if (MP_IEN(j,k) == MP_IEN(k1,k2) .and. j .ne. k1 .and. k .ne. k2)then
                                                MP_IEN(k1,k2) = nmax
                                            end if
                                        end do
                                    end do
                                    
                                    MP_IEN(j,k) = nmax
                                    
                                    continue
                                end if
                            end do
                            
                            continue
                            
                        end do
                    end if
                !end if !PTS == .false.
                continue
                    
                if(nptch .gt. 1)then
                    allocate(Kf(tnodes*nds,tnodes*nds))
                    Kf = 0.0d0
                    
                    allocate(Fext(tnodes*nds,1))
                    Fext = 0.0d0
                    
                    allocate(ImpDisp(tnodes*nds,1),StpFrac(tnodes*nds))
                    ImpDisp = 0.0d0
                    StpFrac = 0.0d0
                    
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
                    
                    allocate(PMesh(tnodes,nds))
                    PMesh = 0.0d0
                    
                    allocate(bcdof(tnodes*nds))
                    bcdof = 1
!                    allocate(Weight(tnodes))
!                    Weight = 0.0d0
                    
!                    allocate(Points(tnodes,nds))
!                    Points = 0.0d0
                    allocate(GCoords(tnodes,nds+1))
                    GCoords = 0.0d0
            
                    allocate(dispdof(tnodes*nds))
                    dispdof = 0.0d0
                    
                    allocate(loaddof(tnodes*nds))
                    loaddof = 0.0d0
                end if 
                
                continue
            
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
            read(1,*)nbc
            do k1=1,nbc
                read(1,*)int1,int2,int3,int4
                bc(int1,int2,int3,int4) = 0
            end do
            
!            allocate(KfEq(tnodes*nds-nbc,tnodes*nds-nbc))
!            KfEq = 0.0d0
!            
!            allocate(dispeqv(tnodes*nds-nbc,1))
!            dispeqv = 0.0d0
!            
!            !allocate(KfEqInv(nnodes*nds-nbc,nnodes*nds-nbc))
!            !KfEqInv = 0.0d0
!            
!            allocate(FextEq(tnodes*nds-nbc,1))
!            FextEq = 0.0d0
!            
!            allocate(dDisp(tnodes*nds,1))
!            dDisp=0.0d0
!    
!            allocate(dddisp(tnodes*nds,1))
!            ddDisp=0.0d0
!            
!            allocate(u(tnodes*nds,1))
!            u=0.00d
            
            continue
        
        elseif(line == '*bcdof')then
            read(1,*)nbc
            
            !allocate(bcdof(nbc))
            
            do k1=1,nbc
                read(1,*)int1,int2
                bcdof((int1-1)*nds+int2) = 0
            end do
            
            count = 0
            do k=1, ncpz
                do j = 1,ncpy
                    do i = 1,ncpx
                        do k1=1,nds
                            count = count + 1
                            if(bcdof(count)==0) bc(i,j,k,k1) = 0
                        end do
                    end do
                enddo
            enddo
            
            continue
            
        elseif(line == '*load')then
            read(1,*)nload
            do k1=1,nload
                read(1,*)int1,int2,int3,int4,temp1
                load(int1,int2,int3,int4) = temp1
            end do
            
            continue
        
        elseif(line == '*loaddof')then
            read(1,*)nload
            
            do k1=1,nload
                read(1,*)int1,int2,temp1
                loaddof((int1-1)*nds+int2) = temp1
            end do
            
!            count = 0
!            do k=1, ncpz
!                do j = 1,ncpy
!                    do i = 1,ncpx
!                        do k1=1,nds
!                            count = count + 1
!                            load(i,j,k,k1) = loaddof(count)
!                        end do
!                    end do
!                enddo
!            enddo
            
            continue
        elseif(line == '*dispdof')then
            read(1,*)ndisp
 
            do k1=1,ndisp
                read(1,*)int1,int2,temp1
                dispdof((int1-1)*nds+int2) = temp1
            end do
            
            continue
            
        elseif(line=='*material')then
            
            do i=1,nptch
                read(1,*)iprops
                
                if(i==1) allocate(props(iprops))
                
                read(1,*)(props(j),j=1,iprops)
                
                if(nptch .gt. 1)then
                    if(i==1)then
                        allocate(MP_iprops(nptch),MP_props(nptch,iprops))
                    end if
                    
                    MP_iprops(i)=iprops
                    MP_props(i,:) = props(:)
                    
                end if
                 
            end do
            
            continue  
        
        elseif(line=='*displacement')then
            read(1,*)ndisp
            do k1=1,ndisp
                read(1,*)int1,int2,int3,int4,temp1
                dispBC(int1,int2,int3,int4) = temp1
            end do
            
            continue
        
        elseif(line=='*gravity')then
            gconst = 0.0d0
            gravdir = 0
            read(1,*)gconst, gravdir
            gravity = .true.
        
        elseif(line=='*pressure') then
            
            allocate(pressvec(nelems),ftag(nelems))
            pressvec = 0.0d0
            ftag=''
            
            read(1,*) ipress
            do k1=1,ipress
                read(1,*)int1,tempc,temp1
                pressvec(int1) = temp1
                ftag(int1) = tempc
            end do
            
            pressure = .true.
            
            continue
       
        elseif(line=='*NLGeom') then
            NLGeom = .true.
        
        elseif(line=='*Increments')then
            read(1,*)incmax
        
        elseif(line=='*Iterations')then
            read(1,*)itermax
            
        elseif(line=='*Contact_PTS')then
            PTS = .true.
            npair = 1
            
            read(1,*)islave
            
            do k1=1,2
                read(1,*)line
                
                if(line=='*Master')then
                    
                    MRigid = .false.
                    
                    read(1,*)p_mst, q_mst
                    read(1,*)ncpx_mst, ncpy_mst
                    
                    allocate(u_knot_mst(ncpx_mst+p_mst+1))
                    u_knot_mst = 0.0d0
                    read (1,*) (u_knot_mst(i),i=1,ncpx_mst+p_mst+1)
                    
                    allocate(v_knot_mst(ncpy_mst+q_mst+1))
                    v_knot_mst = 0.0d0
                    read (1,*) (v_knot_mst(i),i=1,ncpy_mst+q_mst+1)
                    
                    allocate(Conn_mst(ncpx_mst*ncpy_mst))
                    Conn_mst = 0
                    
                    read(1,*)(Conn_mst(i),i=1,ncpx_mst*ncpy_mst)
                    
                    continue
                
                elseif(line=='*RigidMaster')then
                    
                    MRigid = .true.
                    
                    read(1,*)p_mst, q_mst
                    read(1,*)ncpx_mst, ncpy_mst
                    
                    allocate(u_knot_mst(ncpx_mst+p_mst+1))
                    u_knot_mst = 0.0d0
                    read (1,*) (u_knot_mst(i),i=1,ncpx_mst+p_mst+1)
                    
                    allocate(v_knot_mst(ncpy_mst+q_mst+1))
                    v_knot_mst = 0.0d0
                    read (1,*) (v_knot_mst(i),i=1,ncpy_mst+q_mst+1)
                    
                    allocate(Conn_mst(ncpx_mst*ncpy_mst))
                    Conn_mst = 0
                    
                    read(1,*)(Conn_mst(i),i=1,ncpx_mst*ncpy_mst)
                    
                    continue
                    
                elseif(line=='*Slave')then
                    
                    read(1,*)p_slv, q_slv
                    read(1,*)ncpx_slv, ncpy_slv
                    
                    allocate(u_knot_slv(ncpx_slv+p_slv+1))
                    u_knot_slv = 0.0d0
                    read (1,*) (u_knot_slv(i),i=1,ncpx_slv+p_slv+1)
                    
                    allocate(v_knot_slv(ncpy_slv+q_slv+1))
                    v_knot_slv = 0.0d0
                    read (1,*) (v_knot_slv(i),i=1,ncpy_slv+q_slv+1)
                    
                    allocate(Conn_slv(ncpx_slv*ncpy_slv))
                    Conn_slv = 0
                    
                    read(1,*)(Conn_slv(i),i=1,ncpx_slv*ncpy_slv)
                
                    continue
                
                end if
            end do
            continue
        
        elseif(line=='*MultiBody_Contact_PTS')then
            
            PTS = .true.
            
            read(1,*)npair
            
            allocate(MB_islave(npair))
            MB_islave = 0
            
            do k1=1,npair
                read(1,*)MB_islave(k1)
            end do
            
            allocate(MB_MRigid(npair))
            MB_MRigid = .false.
            
            allocate(MB_p_mst(npair),MB_q_mst(npair))
            MB_p_mst = 0
            MB_q_mst = 0
            
            allocate(MB_ncpx_mst(npair),MB_ncpy_mst(npair))
            MB_ncpx_mst = 0
            MB_ncpy_mst = 0
            
            allocate(MB_u_knot_mst(npair,250),MB_v_knot_mst(npair,250))
            MB_u_knot_mst = 0.0d0
            MB_v_knot_mst = 0.0d0
            
            allocate(MB_conn_mst(npair,1000))
            MB_conn_mst = 0
            
            allocate(MB_p_slv(npair),MB_q_slv(npair))
            MB_p_slv = 0
            MB_q_slv = 0
            
            allocate(MB_ncpx_slv(npair),MB_ncpy_slv(npair))
            MB_ncpx_slv = 0
            MB_ncpy_slv = 0
            
            allocate(MB_u_knot_slv(npair,250),MB_v_knot_slv(npair,250))
            MB_u_knot_slv = 0.0d0
            MB_v_knot_slv = 0.0d0
            
            allocate(MB_conn_slv(npair,1000))
            MB_conn_slv = 0
            
            do ipair=1,npair
            
                do k1=1,2
                    read(1,*)line
                    
                    if(line=='*Master')then
                        
                        MB_MRigid(ipair) = .false.
                        
                        read(1,*)MB_p_mst(ipair), MB_q_mst(ipair)
                        read(1,*)MB_ncpx_mst(ipair), MB_ncpy_mst(ipair)
                        
                        read (1,*) (MB_u_knot_mst(ipair,i),i=1,MB_ncpx_mst(ipair)+MB_p_mst(ipair)+1)
                        
                        read (1,*) (MB_v_knot_mst(ipair,i),i=1,MB_ncpy_mst(ipair)+MB_q_mst(ipair)+1)
                        
                        read(1,*)(MB_Conn_mst(ipair,i),i=1,MB_ncpx_mst(ipair)*MB_ncpy_mst(ipair))
                        
                        continue
                    
                    elseif(line=='*RigidMaster')then
                        
                        MB_MRigid(ipair) = .true.
                        
                        read(1,*)MB_p_mst(ipair), MB_q_mst(ipair)
                        read(1,*)MB_ncpx_mst(ipair), MB_ncpy_mst(ipair)
                        
                        read (1,*) (MB_u_knot_mst(ipair,i),i=1,MB_ncpx_mst(ipair)+MB_p_mst(ipair)+1)
                        
                        read (1,*) (MB_v_knot_mst(ipair,i),i=1,MB_ncpy_mst(ipair)+MB_q_mst(ipair)+1)
                        
                        read(1,*)(MB_Conn_mst(ipair,i),i=1,MB_ncpx_mst(ipair)*MB_ncpy_mst(ipair))
                        
                        continue
                        
                    elseif(line=='*Slave')then
                        
                        read(1,*)MB_p_slv(ipair), MB_q_slv(ipair)
                        read(1,*)MB_ncpx_slv(ipair), MB_ncpy_slv(ipair)
                        
                        read (1,*) (MB_u_knot_slv(ipair,i),i=1,MB_ncpx_slv(ipair)+MB_p_slv(ipair)+1)
                        
                        read (1,*) (MB_v_knot_slv(ipair,i),i=1,MB_ncpy_slv(ipair)+MB_q_slv(ipair)+1)
                        
                        read(1,*)(MB_Conn_slv(ipair,i),i=1,MB_ncpx_slv(ipair)*MB_ncpy_slv(ipair))
                    
                        continue
                    
                    end if
                end do
                
            end do !ipair
        endif !line
    
    end do
    
    allocate(stress(nelems,npimax,ntens),dstress(nelems,npimax,ntens))
    stress = 0.0d0
    dstress = 0.0d0
    
    allocate(strain(nelems,npimax,ntens),dstrain(nelems,npimax,ntens))
    strain = 0.0d0
    dstrain = 0.0d0
    
    allocate(stress_conv(nelems,npimax,ntens),strain_conv(nelems,npimax,ntens))
    stress_conv = 0.0d0
    strain_conv = 0.0d0
    
    allocate(hard(nelems,npimax),hard_conv(nelems,npimax))
    hard = 0.0d0
    hard_conv = 0.0d0
    
    allocate(TmatD(nelems,npimax,ntens,ntens))
    allocate(matD(ntens,ntens))
    TmatD = 0.0d0
    
    allocate(laxis(nelems,npimax,nds,nds))
    allocate(laxisconv(nelems,npimax,nds,nds))
    laxis = 0.0d0
    laxisconv = 0.0d0
    
!    call MatLinearElastic3D(props(1),props(2),MatD)
!    
!    do k1=1,nelems
!        do k2=1,npi
!            TmatD(k1,k2,:,:) = matD(:,:)
!        end do
!    end do
    
    allocate(KfEq(tnodes*nds-nbc,tnodes*nds-nbc))
    KfEq = 0.0d0
    
    allocate(KfEqInv(tnodes*nds-nbc,tnodes*nds-nbc))
    KfEqInv = 0.0d0
    
    allocate(FextEq(tnodes*nds-nbc,1))
    FextEq = 0.0d0
    
    allocate(dDisp(tnodes*nds,1))
    dDisp=0.0d0

    allocate(dddisp(tnodes*nds,1))
    ddDisp=0.0d0
    
    allocate(u(tnodes*nds,1))
    u=0.0d0
    
    allocate(dispeqv(tnodes*nds-nbc,1))
    dispeqv = 0.0d0
    
!    allocate(coordi(tnodes,nds))
    
    close(1,status='keep')
    
    continue
    
end