!--------------------------------------------------------------------------------------------------
!
!        %%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%
!              %%          %%               %%             %%
!             %%          %%               %%             %%
!            %%          %%               %%             %%
!           %%          %%               %%             %%
!          %%          %%               %%             %%
!         %%          %%               %%             %%
! %%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%
!
!--------------------------------------------------------------------------------------------------
!
! Isogeometric COde (ICO) for 2D applications
!    
!--------------------------------------------------------------------------------------------------    
    program ICO
    
    use ifport
    use Mod_Variables
    
    implicit none
    
    character*256::FileName
    integer(4)::i,j,k,nze,k1,k2,k3
    integer(4)::iptch,nzep,iLM,loc_num,istp
    integer(4)::count,iter,inc,count2
    real(8)::da
    
    real(8), dimension(:,:), allocatable::Ke,KBar,dispin,finte
    real(8)::absRes,absDisp
    
    real(8)::RF1,RF2
    real(8)::real1,real2
    
    real(8),dimension(:,:), allocatable::coordi
    
    call ICOHeader()
    
    !write(*,*)'Input file: '
    !read(*,*)FileName
    
    FileName = 'Contact1g'
    !FileName = 'ContactConv1b'
    !FileName = 'Multistep2'
    !FileName = 'Ironing_p2'
    !FileName = 'PtS_HalfRing_8el_p2q2'
    !FileName = 'SingleQuad'
    !FileName = 'PtS_Punch_100el_p2q2'
    !FileName = 'PtS_PunchFull_392el_p2q2'
    !FileName = 'Hertz1_p2q2'
    !FileName = 'Hertz_16x8_Mesh3_p2q2'
    !FileName = 'Hertz_32x16_p2q2'
    !FileName = 'ContactPT_Crisfield'
    !FileName = 'PtS_Tube_40el_p2q2'
    !FileName = 'Hertz_16x16_Mesh1_Bathe_p2q2'
    !FileName = 'FlatPunch_16x16_Mesh3_p2q2'
    FileName = 'SingleElementQ2'
    
    
    open(file='reac.txt', unit=123)
    
    incmax = 1
    itermax = 100
    
    nlgeom = .false.
    
    !------------------------------------------------------------------------------
    !Read input file
    !------------------------------------------------------------------------------
    call ReadInFile(FileName)
    
    !------------------------------------------------------------------------------
    !Data for underformed mesh (Matlab)
    !------------------------------------------------------------------------------
    call MatlabIn()
    
    !------------------------------------------------------------------------------
    !Initialisation of some required vectors/matrices
    !------------------------------------------------------------------------------
    call AllocVectors()
    
    !Initialise number of Lagrange Multipliers
    iLM = 0
    
    !------------------------------------------------------------------------------
    ! Determine Collocation Points for PTS algorithm
    !------------------------------------------------------------------------------
    if(ContactPTS == .true.)then
        call PTS_CollocationPoints()
    elseif(ContactGPTS == .true.)then
        call GPTS_Initialisation()
    end if
    
    !------------------------------------------------------------------------------
    ! Step cycle
    !------------------------------------------------------------------------------
    do istp=1,nstp
    
        !------------------------------------------------------------------------------
        ! Multistep allocations
        !------------------------------------------------------------------------------
        if(multistep==.true.)then
            call STPAlloc(istp)
        end if
        
        !------------------------------------------------------------------------------
        ! Multipatch allocations
        !------------------------------------------------------------------------------
        if(nptch == 1)then 
            
            if(istp==1)then
                allocate(coordi(ncpx*ncpy,nds))
                coordi = 0.0d0
                
                allocate(PMesh0(ncpx*ncpy,nds))
                PMesh0 = 0.0d0
            end if

            !------------------------------------------------------------------------------
            !Allocate coordinates of the mesh at the begining of the step
            !------------------------------------------------------------------------------
            count = 0
            do j=1,ncpy
                do i=1,ncpx
                    count = count+1
                    PMesh0(count,:) = B_net(i,j,1:nds)
                end do
            end do
            
            !------------------------------------------------------------------------------
            !Assemble the external forces vector
            !------------------------------------------------------------------------------
            count = 0
            do j=1,ncpy
                do i =1,ncpx
                    do k =1,nds
                        count = count + 1
                        Fext(count,1) = load(i,j,k)
                        ImpDisp(count,1) = dispBC(i,j,k)
                    end do
                end do
            end do
            
            Fini = Fext
            Fext = 0.0d0
            Finc=Fini/(1.0d0*incmax)
            
            !ImpDisp = 0.0d0
            ImpDisp(:,1) = ImpDisp(:,1) + dispdof(:)
        
            Fini(:,1) = Fini(:,1) + loaddof(:)
            Fext = 0.0d0
            Finc=Fini/(1.0d0*incmax)
            
        else
            ImpDisp = 0.0d0
            ImpDisp(:,1) = dispdof(:)
        
            Fini(:,1) = loaddof(:)
            Fext = 0.0d0
            Finc=Fini/(1.0d0*incmax)
        end if
        
        !Some Arrays Initializations
        !SEnergyConv = 0.0d0
        !FintLM = 0.0d0
        !FintLMConv = 0.0d0
        
        !------------------------------------------------------------------------------
        ! Increment cycle
        !------------------------------------------------------------------------------
        inc = 0
        do inc=1,incmax
            
            stpfrac=0.0d0
            dDisp=0.0d0
            ddDisp=0.0d0
            
            !Get converged values
            SEnergy = SEnergyConv
            
            !Increment the external load
            !Fext=Fext+Finc
            
            Fext = Finc
            
            !KCont = 0.0d0
            
            !------------------------------------------------------------------------------
            ! Iteration cycle
            !------------------------------------------------------------------------------
            iter = 0
            do while(iter <= itermax)
                iter = iter + 1
                
                !Use only converged variables in the iteration cycle
                Stress=Stress_Conv
                Strain=Strain_Conv
                hard=hard_conv
                laxis=laxis_conv
                !Lagrange = LagrangeConv
                
                !Some initialisations
                NZEp = 0
                GP_coords = 0.0d0
                !if(iter== 1) then
                !if(inc == 1 .and. iter == 1 .and. istp==1)then
                if(inc == 1 .and. iter == 1)then
                    
                    !------------------------------------------------------------------------------
                    ! Patch cycle
                    !------------------------------------------------------------------------------
                    Fint = 0.0d0
                    Kf = 0.0d0
                    do iptch=1,nptch
                        
                        if(nptch .gt. 1)then
                            
                            call MPAlloc(iptch)
                            
                            allocate(Ke((p+1)*(q+1)*nds,(p+1)*(q+1)*nds))
                            allocate(KBar(p*q*nds,p*q*nds))
                            Ke = 0.0d0
                            KBar = 0.0d0
                            
                            allocate(dispin((p+1)*(q+1)*nds,1))
                            dispin = 0.0d0
                            allocate(Finte((p+1)*(q+1)*nds,1))
                            
                            allocate(coordi(ncpx*ncpy,nds))
                            coordi = 0.0d0
                            
                            !Generate connectivity arrays
                            call gen_ien_inn()
                            
                            continue
                            
                        elseif(nptch == 1 .and. inc == 1 .and. iter == 1)then
                            allocate(Ke((p+1)*(q+1)*nds,(p+1)*(q+1)*nds))
                            allocate(KBar(p*q*nds,p*q*nds))
                            Ke = 0.0d0
                            KBar = 0.0d0
                            
                            allocate(dispin((p+1)*(q+1)*nds,1))
                            dispin = 0.0d0
                            allocate(Finte((p+1)*(q+1)*nds,1))
                            
                            !Generate connectivity arrays
                            call gen_ien_inn()
                        end if
                    
                        !------------------------------------------------------------------------------
                        ! Element cycle
                        !------------------------------------------------------------------------------
                        nze = 0
                        do i=1, ncpx+p
                            do j=1,ncpy+q
                                
                                !Check if it is a non-zero element
                                if((u_knot(i) .ne. u_knot(i+1)) .and. (v_knot(j) .ne. v_knot(j+1))) then
                                    nzep = nzep + 1
                                    nze = nze + 1
                                    
                                    SelectCase(Elcode)
                                        case('Quad4E')
                                            call ElQuad4E(nze,nzep,Ke,dispin,Finte,inc)
                                        case('Quad4S')
                                            call ElQuad4S(nze,nzep,Ke,dispin,Finte,inc)
                                        case('Quad9E')
                                            call ElQuad9E(nze,nzep,Ke,dispin,Finte,inc)
                                        case('Quad9S')
                                            call ElQuad9S(nze,nzep,Ke,dispin,Finte,inc)
                                        case('Quad16E')
                                            call ElQuad16E(nze,nzep,Ke,dispin,Finte,inc)
                                        case('Quad16S')
                                            call ElQuad16S(nze,nzep,Ke,dispin,Finte,inc)
!                                        case('Quad9EBBar')
!                                            call ElQuad9EBBar(nze,Ke,KBar,dispin,Finte)
                                    EndSelect
                                    
                                    if(nptch == 1)then
                                        call Assembly2D(nze,Ke,finte)
                                    else
                                        call AssemblyPatch(nze,nzep,Ke,finte)
                                    end if
                                   
                                    continue
                                    
                                endif
                                
                            end do
                        end do
                        
                        SEnergy = 0.0d0
                    
                        if(nptch .gt. 1)then
                            !call AssemblyPatch(nzep,Ke,finte)
                            call MPDalloc(iptch)
                            
                            deallocate(Ke,Kbar,dispin)
                            deallocate(Finte)
                            deallocate(coordi)
                            deallocate(INN,IEN,conn)

                        else
                            !deallocate(Ke,KBar,dispin,Finte,coordi)
                        end if
                        
                    end do !ipatch
                end if !inc == 1 .and. iter == 1
                
                
                !------------------------------------------------------------------------------
                !Obtained updated global coordinates (for contact purpouses)
                !------------------------------------------------------------------------------
                nzep = 0
                do iptch=1,nptch
                    if(nptch .gt. 1)then
                        call MPAlloc(iptch)
                        
                        allocate(Ke((p+1)*(q+1)*nds,(p+1)*(q+1)*nds))
                        allocate(KBar(p*q*nds,p*q*nds))
                        Ke = 0.0d0
                        KBar = 0.0d0
                        
                        allocate(dispin((p+1)*(q+1)*nds,1))
                        dispin = 0.0d0
                        allocate(Finte((p+1)*(q+1)*nds,1))
                        
                        allocate(coordi(ncpx*ncpy,nds))
                        coordi = 0.0d0
                        
                        !Generate connectivity arrays
                        call gen_ien_inn()
                        
                        continue
                    end if

                    if(nptch .gt. 1)then
                        nze = 0
                        do i=1, ncpx+p
                            do j=1,ncpy+q
                                
                                !Check if it is a non-zero element
                                if((u_knot(i) .ne. u_knot(i+1)) .and. (v_knot(j) .ne. v_knot(j+1))) then
                                    nzep = nzep + 1
                                    nze = nze + 1
                                    
                                    loc_num = 0
                                    do k1=1,(p+1)
                                        do k2=1,(q+1)
                                            loc_num = loc_num + 1
                                            do k3=1,nds
                                                GCoords(MP_Conn(nzep,loc_num),k3) = Points0(IEN(nze,loc_num),k3) + u(MP_conn(nzep,loc_num)*nds-(nds-k3),1)
                                            end do
                                            
                                            GCoords(MP_Conn(nzep,loc_num),3) = Weights(IEN(nze,loc_num))
                                            
                                        end do   
                                    end do
                                    
                                endif
                                
                            end do
                        end do
                    
                        !call AssemblyPatch(nzep,Ke,finte)
                        call MPDalloc(iptch)
                        
                        deallocate(Ke,Kbar,dispin)
                        deallocate(Finte)
                        deallocate(coordi)
                        deallocate(INN,IEN,conn)

                    else
                        !deallocate(Ke,KBar,dispin,Finte,coordi)
                    end if
                    
                    if(nptch .gt. 1)then
                        
                        !Slave Segment IEN array
                        if(iter==1 .and. inc==1 .and. istp==1 .and. iptch==1)then
                            allocate(Points_slv(ncp_s,nds))
                            allocate(Weights_slv(ncp_s))
                            
                            allocate(INN_s(ncp_s,nds))
                            allocate(IEN_s(ncp_s-p,p+1))
                            call PTS_SlaveIEN(p_s,ncp_s,ncp_s-p,nds,u_knot_slave,INN_s,IEN_s)
                        end if
                        
                        !Slave Segment control points coordinates and weights
                        do k1=1, ncp_s
                            do k2 = 1, nds+1
                                b_slave(k1,k2) = GCoords(conn_slave(k1),k2)
                            end do
                        end do
                        
                        count = 0
                        Points_slv = 0.0d0
                        Weights_slv = 0.0d0
                        do i = 1,ncp_s
                            count = count + 1
                            Points_slv(count,1) = b_slave(i,1)
                            Points_slv(count,2) = b_slave(i,2)
                            Weights_slv(count)  = b_slave(i,3)
                        enddo
                        
                        !call PTS_SlaveSegment(iter,inc,istp,iptch)
                        
                    end if
                    
                end do !ipatch
                
                
                KT = Kf
                
                !-----------------------------------------------------------------
                !Apply zero-displacement (Dirichelet) boundary conditions
                !-----------------------------------------------------------------
                call ApplyBCs()
                
                !-----------------------------------------------------------------
                !Apply prescibed displacement conditions
                !-----------------------------------------------------------------
                call ApplyPrescDisp()
                
                !-----------------------------------------------------------------
                !Determine Gap and Contact Contributions
                !-----------------------------------------------------------------
                !if(ContactPTS == .true.)then
                if(ContactPTS == .true. .and. iter .gt. 1)then
                    call PTS_ComputeGap(iter,inc,istp,iLM)
                    
                elseif(ContactGPTS == .true. .and. iter .gt. 1)then
                !elseif(ContactGPTS == .true. )then
                    call GPTS_Contact(iter,inc,istp)
                      
                else
                    !-----------------------------------------------------------------
                    !Assemble equivalent reduced matrices
                    !-----------------------------------------------------------------
                    call AssembleReduced()
                    
                    !-----------------------------------------------------------------
                    !Solve system of equations
                    !-----------------------------------------------------------------
                    dispeqv = 0.0d0
                    if(nptch == 1)then
                        call Gauss (nnodes*nds-nbc,Kfeq,FextEq,dispeqv)
                    else
                        call Gauss (tnodes*nds-nbc,Kfeq,FextEq,dispeqv)
                    end if
                    
                    if(nptch == 1)then
                        j=1
                        ddDisp = 0.0d0
                        do i=1,nnodes*nds
                            if(redcount(i)==1)then
                                ddDisp(i,1)=dispeqv(j,1)
                                j=j+1
                            end if
                        end do
                    else
                        j=1
                        ddDisp = 0.0d0
                        do i=1,tnodes*nds
                            if(redcount(i)==1)then
                                ddDisp(i,1)=dispeqv(j,1)
                                j=j+1
                            end if
                        end do
                    end if
                
                end if !ContactPTS == .true.
                
                !-----------------------------------------------------------------
                ! Increment of the displacement increment
                !-----------------------------------------------------------------
                dDisp = dDisp + ddDisp
                continue
                
                
                !-----------------------------------------------------------------
                ! Check for contact
                !-----------------------------------------------------------------
!                if(ContactPTS == .true.)then
!                    if(istp .gt. 1 .and. inc == 1 .and. iter==1)then
!                        continue
!                    else
!                       call CtPTS(iter,inc,istp,iLM)
!                    end if
!                end if
                
                !-----------------------------------------------------------------
                ! Compute reaction forces
                !-----------------------------------------------------------------
                count = 0
                do i=1,tnodes*nds
                    count = count + 1
                    if(BCdof(i) == 0) then
                        do j=1,tnodes*nds
                            Reac(i,1) = Reac(i,1) + KT(i,j)*dddisp(j,1)
                        end do
                    end if
                end do
                
                continue
                
                !-----------------------------------------------------------------
                ! Patch Cycle (Second)
                !-----------------------------------------------------------------
                NZE = 0
                NZEp = 0
                Kf = 0.0d0 !KCont
                Fint = 0.0d0
                SEnergy = 0.0d0
                GP_coords = 0.0d0
                do iptch=1,nptch
                    
                    if(nptch .gt. 1)then
                        call MPAlloc(iptch)
                        
                        allocate(Ke((p+1)*(q+1)*nds,(p+1)*(q+1)*nds))
                        allocate(KBar(p*q*nds,p*q*nds))
                        Ke = 0.0d0
                        KBar = 0.0d0
                        
                        allocate(dispin((p+1)*(q+1)*nds,1))
                        dispin = 0.0d0
                        allocate(Finte((p+1)*(q+1)*nds,1))
                        
                        allocate(coordi(ncpx*ncpy,nds))
                        coordi = 0.0d0
                        
                        !Generate connectivity arrays
                        call gen_ien_inn()
                        
                        continue
                        
                    else
                        count = 0
                        PMesh = 0.0d0
                        do j=1,ncpy
                            do i=1,ncpx
                                count = count+1
                                PMesh(count,:) = B_net(i,j,1:nds)
                                PMesh0(count,:) = B_net(i,j,1:nds)
                            end do
                        end do
                    end if
                    
                    !-----------------------------------------------------------------
                    ! Elements Cycle (Second)
                    !-----------------------------------------------------------------
                    nze = 0
                    do i=1, ncpx+p
                        do j=1,ncpy+q
                            
                            !Check if it is a non-zero element
                            if((u_knot(i) .ne. u_knot(i+1)) .and. (v_knot(j) .ne. v_knot(j+1))) then
                                nzep = nzep + 1
                                nze = nze + 1
                                
                                !Elemental displacement -------------
                                if(nptch==1)then
                                    do k1=1,(p+1)*(q+1)
                                        dispin(k1*2-1,1) = dDisp(conn(nze,k1)*2-1,1)
                                        dispin(k1*2  ,1) = dDisp(conn(nze,k1)*2  ,1) 
                                    end do
                                else
                                    do k1=1,(p+1)*(q+1)
                                        dispin(k1*2-1,1) = dDisp(MP_conn(nzep,k1)*2-1,1)
                                        dispin(k1*2  ,1) = dDisp(MP_conn(nzep,k1)*2  ,1) 
                                    end do
                                end if
                                
                                !Update Points coordinates for geometric nonlinear analysis
                                if(nlgeom==.true.) call UpdatePoints(nze,nzep)
                                
                                !Create deformed mesh for output purpouses
                                call OutputMesh(nze,nzep,iptch)

                                !Select element subroutine
                                SelectCase(Elcode)
                                    case('Quad4E')
                                        call ElQuad4E(nze,nzep,Ke,dispin,Finte,inc)
                                    case('Quad4S')
                                        call ElQuad4S(nze,nzep,Ke,dispin,Finte,inc)
                                    case('Quad9E')
                                        call ElQuad9E(nze,nzep,Ke,dispin,Finte,inc)
                                    case('Quad9S')
                                        call ElQuad9S(nze,nzep,Ke,dispin,Finte,inc)
                                    case('Quad16E')
                                        call ElQuad16E(nze,nzep,Ke,dispin,Finte,inc)
                                    case('Quad16S')
                                        call ElQuad16S(nze,nzep,Ke,dispin,Finte,inc)
!                                    case('Quad9EBBar')
!                                        call ElQuad9EBBar(nze,Ke,KBar,dispin,Finte)
                                EndSelect
                                
                                !Assemble global matrices
                                if(nptch == 1)then
                                    call Assembly2D(nze,Ke,finte)
                                else
                                    call AssemblyPatch(nze,nzep,Ke,finte)
                                end if
                                
                            endif
                            
                        end do
                    end do
                    
                    if(nptch .gt. 1)then
                        !call AssemblyPatch(nzep,Ke,finte)
                        call MPDalloc(iptch)
                        
                        deallocate(Ke,Kbar,dispin)
                        deallocate(Finte)
                        deallocate(coordi)
                        deallocate(INN,IEN,conn)

                    else
                        !deallocate(Ke,KBar,dispin,Finte,coordi)
                    end if
                    
                end do !iptch
                
                !-----------------------------------------------------------------
                ! Eliminate internal forces at fixed boundaries and at points 
                ! with prescribed displacement
                !-----------------------------------------------------------------
                call EliminateReac(iLM)
                
                continue
                
                
                !-----------------------------------------------------------------
                ! Contact Forces
                !-----------------------------------------------------------------
!                if(ContactPTS == .true. .and. iter .gt. 1)then
!                    call PTS_ContactForces(iter,inc)
!                else
!                    !FintLM = 0.0d0
!                    !Lagrange = 0.0d0
!                    continue
!                end if
                
                
                !-----------------------------------------------------------------
                ! Evaluate Residual Forces
                !-----------------------------------------------------------------
                !Residual = Fext(n)- Fint
                Res=0.0d0
                if(iLM == 0)then
                    !if(ContactPTS == .true. .and. iter .gt. 1)then
                    if(ContactPTS == .true.)then
                        Res=Fini/(1.0d0*incmax)*(1.0d0*inc) - Fint - FintLM !- FintLMConv
                    elseif(ContactGPTS == .true.)then
                        Res=Fini/(1.0d0*incmax)*(1.0d0*inc) - Fint - Fct
                    else
                        Res=Fini/(1.0d0*incmax)*(1.0d0*inc)-Fint
                    end if
                else
                    !if(ContactPTS == .true. .and. iter .gt. 1)then
                    if(ContactPTS == .true.)then
                        Res=Fini/(1.0d0*incmax)*(1.0d0*inc) - Fint - FintLM !- FintLMConv
                    elseif(ContactGPTS == .true.)then
                        Res=Fini/(1.0d0*incmax)*(1.0d0*inc) - Fint - Fct
                    else
                        Res=Fini/(1.0d0*incmax)*(1.0d0*inc)-Fint
                    end if
                end if
                
                
                continue
                
                !-----------------------------------------------------------------
                ! Percentage of applied prescribed displacements
                !-----------------------------------------------------------------
                if(nptch==1)then
                    do i=1,ncpx*ncpy*nds
                        if(ImpDisp(i,1) /= 0.0d0)then
                            stpfrac(i)=(u(i,1)+dDisp(i,1))/(ImpDisp(i,1)/(1.0d0*incmax)*(1.0d0*inc))
                            continue
                        end if
                    end do
                else
                    do i=1,tnodes*nds
                        if(ImpDisp(i,1) /= 0.0d0)then
                            stpfrac(i)=(u(i,1)+dDisp(i,1))/(ImpDisp(i,1)/(1.0d0*incmax)*(1.0d0*inc))
                            continue
                        end if
                    end do
                end if
                
                !-----------------------------------------------------------------
                ! Fictional step to activate geometric non-linearity
                !-----------------------------------------------------------------
                if(iter==1 .and. inc==1 .and. nlgeom==.true.) then
                    Fext=Res
                    continue
                    goto 80
                end if
                
                !-----------------------------------------------------------------
                ! Convergence Check
                !-----------------------------------------------------------------
                call ConvCheck (inc,AbsRes,AbsDisp)

                !-----------------------------------------------------------------
                ! Contact stress at slave segment (PTS)
                !-----------------------------------------------------------------
                if(ContactPTS == .true.) call PTS_SlaveSegment(iter,inc,istp,iptch)
                           
                !-----------------------------------------------------------------
                !Output relevant data
                !-----------------------------------------------------------------
                call MatlabOut()
                
                continue
                
                
                if(absRes .lt. RToler .and. absdisp .lt. RToler) then
                    !Fext = Res
                    
                    stpfrac = 0.0d0
                    
                    !Residual after convergence -----
                    !Resconv=Res
                    hard_conv=hard
                    
                    !Update converged stress
                    Stress_Conv=Stress_Conv + dStress
                    Strain_Conv=Strain_Conv + dStrain
                    
                    !Update converged displacement -----
                    u=u+dDisp
                    
                    !Update Strain energy -----
                    SEnergyConv = SEnergyConv + SEnergy
                    
                    !Update converged local axis -------
                    laxis_conv=laxis
                    
                    FintLMConv = FintLM
                    
                    !Update right Cauchy-Green tensor
                    !StoreCConv = StoreC
                    
                    !Write information to matlab output
                    call MatlabOut()
                    
!                    real1 =  1.0d+15
!                    real2 = -1.0d+15
!                    
!                    do i=1,sum(MP_nelems(:))
!                        do j=1,npi
!                            if(Stress_Conv(i,j,2) .lt. real1) real1 = Stress_Conv(i,j,2)
!                            if(Stress_Conv(i,j,2) .gt. real2) real2 = Stress_Conv(i,j,2)
!                        end do
!                    end do

!                    RF1 = 0.0d0
!                    RF2 = 0.0d0
!                    !do i=1,14
!                    do i=1,38
!                        RF1 = RF1 + Reac(i*2-1,1)
!                        RF2 = RF2 + Reac(i*2  ,1)
!                    end do
                    
!                    write(123,FMT=11)RF1, RF2, sumLag, u(2,1)
!                    11 format(6(E,1x))
                    !LagrangeConv = Lagrange
                    
                    continue
                    
                    !FintLM = 0.0d0
                    !Lagrange = 0.0d0
                    
                        
                    goto 90
                else
                    
                    Fext = Res
                    !deallocate(FintLM,FintLMConv)
                    
                    continue
                    !goto 10
                end if  
            
                80 continue
            
            end do !iteration cycle
            
            write(*,FMT=10)itermax
            10 format('Analysis failed to converge after ', I4, ' iterations')
            goto 999
            
        90 continue
        
        end do !Increment cycle
    
        !------------------------------------------------------------------------------
        ! Multistep deallocations
        !------------------------------------------------------------------------------
        call STPDalloc(istp)
      
        !deallocate(Lagrange)
      
        continue
    end do !Step Cyle
    
    999 continue
    
    continue
    
    close(123)
    
    end program ICO

