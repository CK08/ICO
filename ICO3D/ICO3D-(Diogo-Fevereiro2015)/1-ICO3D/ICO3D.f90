    program ICO3D
    
    !use ifport
    use Mod_Variables
    
    implicit none
    
    character*256::FileName
    integer(4)::i,j,k,l,nze,nzep,k1,k2,k3,k4,kt,errorflag,loc_num
    integer(4)::count,iter,inc
    integer(4)::isolver
    integer(4)::nred,iptch
    
    !Iterative Solver Variables ----
    !-----------------------------------------------------------------------------------------------------------------------------------------------
    ! Esta declaracao foi comentada devido a problemas de compilacao (nao sao necessarias para a analise - teste de solver2 e 3 pelo Caseiro): Diogo  
    !----------------------------------------------------------------------------------------------------------------------------------------------- 
    !external::aprod,msolve
    real(8),dimension(:),allocatable::r1,r2,vs,ws,ys
    integer(4)::istop,itn
    real(8)::Anorm,Rnorm,ynorm,Acond,STol
    integer(4)::nout,itnlim
    real(8)::shift
    logical::checkA,goodb,precon
    
    !Variables for solver3
    !------------------------------------------------------------------------------------------------------------------------
    ! Variables for solver3 (Nota: como apenas se está a usar o solver de Gauss, estas variaveis nao sao necessarias - Diogo)
    !------------------------------------------------------------------------------------------------------------------------
!    real(8),dimension(:,:),allocatable::coef
!    real(8),dimension(:),allocatable::AR
!    integer(4),dimension(:),allocatable::IA,JA
!    real(8),dimension(:),allocatable::wksp,rparm
!    integer(4),dimension(:),allocatable::colnz,iwksp,iparm
!    integer(4)::mnz,countx,county
    
    real(8)::da,tst1,tst2,tst3
    
    real(8), dimension(:,:), allocatable::Ke,dispin,finte
    real(8)::absRes,absDisp,sumRes,sumFext,sumd,sumdt
    
    real(8), dimension(:,:), allocatable::MatA,MatC,MatV
    real(8), dimension(:), allocatable::temp1
    
    !real(8),dimension(:,:), allocatable::coordi
    real(8),dimension(:,:), allocatable::MMult
    
    !------------------------------------------------------------------------------
    ! Variavies usadas para o modulo ifport (calculo do tempo de computacao): Diogo 
    !------------------------------------------------------------------------------
    real(4)::CPUt
!    real(4),dimension(2)::atime
    cput = 0 !Em MAC ainda nao ha modulo para calculo do CPUt !!!!!!!!!!!!!!!!!!!!!!!!!!
    !------------------------------------------------------------------
    ! Write data to screen and file
    !------------------------------------------------------------------
    call ScreenData(1,inc,iter,AbsRes,absdisp,cput)
    
    read(*,*)FileName
    
    
    !FileName = 'BernoulliStrBeam_4el_p2q2w2'
    !FileName = 'CCP_16el_p2q2w2'
    !FileName = 'CookPlast_8el_p1q1w1'
    !FileName = 'NLG_MBending_10el_p2q2w2'
    !FileName = 'FullSphere_32el_p2q2w2'
    !FileName = 'CylShellStrip_32el_p2q2w2'
    !FileName = 'CylShellStrip_ML5_p2q2w2'
    !FileName = 'TwoHex'
    !FileName = 'BCs'
    !FileName = 'CSP_VL2_p2q2w2'
    !FileName = 'MeshDist_e03_p2q2w2'
    !FileName = 'CSP_4el_p2q2w2'
    !FileName = 'PTS_5'
    !FileName = 'PTS_RigidPunch_10x10_MeshII_p2q2w2'
    !FileName = 'PTS_ContactPT'
    !FileName = 'PTS_MultiSurf1'
    !FileName = 'PTS_Stamping'
    !FileName = 'PTS_Cylinder'
    !FileName = 'PTS_Cylinder3'
    !FileName = 'PTS_Extrusion'
    
    !FileName = 'BernoulliStrBeam_16el_p2q2w2'
    !FileName = 'NLG_Ring_16x2x1_p2q2w2'
    !FileName = 'FullHemisphere'
    
    
    itermax = 100
    incmax = 1
    nlgeom=.false.
    
    eigen = .false.
    
    isolver = 1
    
!    open(unit=999,file='Basis.txt')
!    close(999)
    
    !------------------------------------------------------------------------------
    !Read input file
    !------------------------------------------------------------------------------
    call ScreenData(2,inc,iter,AbsRes,absdisp,cput)
    
    call ReadInFile(FileName)
    
    call ScreenData(3,inc,iter,AbsRes,absdisp,cput)
    
    !------------------------------------------------------------------------------
    !Data for underformed mesh (Matlab)
    !------------------------------------------------------------------------------
    call MatlabIn()
    
    !------------------------------------------------------------------------------
    ! Write information to screen
    !------------------------------------------------------------------------------
    call ScreenData(4,inc,iter,AbsRes,absdisp,cput)
    
    open(unit=9,file='Results.txt', status = 'REPLACE')
    write(9,*)''
    write(9,*)'ICO3D OUTPUT FILE'
    write(9,*)''
    close(9)
    
    open(unit=17,file='ReactionForces.txt', status = 'REPLACE')
    write(17,*)''
    write(17,*)'Reaction Forces'
    write(17,*)''
    close(17)
    
    !------------------------------------------------------------------------------
    ! Allocate some local variables
    !------------------------------------------------------------------------------
    if(nptch == 1)then
        allocate(Ke((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds))
        Ke = 0.0d0
        allocate(dispin((p+1)*(q+1)*(w+1)*nds,1))
        dispin = 0.0d0
        allocate(Finte((p+1)*(q+1)*(w+1)*nds,1))
        
        !------------------------------------------------------------------------------
        !Generate connectivity arrays
        !------------------------------------------------------------------------------
        call gen_ien_inn()
        call gen_ien_inn_bar()
    end if
    
    allocate(Reac(tnodes*nds))
    Reac = 0.0d0
    allocate(KTot(tnodes*nds,tnodes*nds))
    KTot = 0.0d0
    
    
    !coordi = 0.0d0
    
!    if(elcode == 'Hex8BBar' .or. elcode == 'Hex27BBar')then
!        allocate(MatA((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds))
!        MatA = 0.0d0
!        allocate(MatC((p+1)*(q+1)*(w+1)*nds,p*q*w))
!        MatC = 0.0d0
!        allocate(MatV(p*q*w,p*q*w))
!        MatV = 0.0d0
!        allocate(temp1((ncpx-1)*(ncpy-1)*(ncpz-1)))
!        Temp1 = 0.0d0
!    elseif(elcode == 'Hex27EAS_PW')then
!        allocate(MatA((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds))
!        MatA = 0.0d0
!        allocate(MatC((p+1)*(q+1)*(w+1)*nds,nalpha))
!        MatC = 0.0d0
!        allocate(MatV(nalpha,nalpha))
!        MatV = 0.0d0
!        allocate(temp1(nalpha))
!        Temp1 = 0.0d0
!    elseif(elcode == 'Hex27ANS_PW')then
!        allocate(MatA((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds))
!        MatA = 0.0d0
!        allocate(MatC((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds))
!        MatC = 0.0d0
!        allocate(MatV((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds))
!        MatV = 0.0d0
!        allocate(temp1((p+1)*(q+1)*(w+1)*nds))
!        Temp1 = 0.0d0
!    end if
    
    !Generate control polygon in GiD
    !call GiDMesh()
    
!    !Generate equivalent FEM mesh
!    call GiDEqMesh(coordi)
    
    !------------------------------------------------------------------------------
    ! Assemble the external forces vector
    !------------------------------------------------------------------------------
    if(nptch==1)then
        count = 0
        do k=1,ncpz
            do j=1,ncpy
                do i =1,ncpx
                    do l =1,nds
                        count = count + 1
                        Fext(count,1) = load(i,j,k,l) + loaddof(count)
                        ImpDisp(count,1) = dispBC(i,j,k,l) + dispdof(count)
                    end do
                end do
            end do
        end do
    else
        do i=1,tnodes*nds
            Fext(i,1) =  loaddof(i)
            ImpDisp(i,1) = dispdof(i)
        end do
    end if
    
    Fini = Fext
    Fext = 0.0d0
    Finc=Fini/(1.0d0*incmax)
    
    SEnergyConv = 0.0d0
    
    call ScreenData(5,inc,iter,AbsRes,absdisp,cput)
    
    !------------------------------------------------------------------------------
    ! Determine Collocation Points for PTS algorithm
    !------------------------------------------------------------------------------
    if(PTS == .true.)then
        call PTS_CollocationPoints()
    end if
    
    !------------------------------------------------------------------------------
    ! Increment cycle
    !------------------------------------------------------------------------------
    inc = 0
    do inc=1,incmax
        
        call ScreenData(6,inc,iter,AbsRes,absdisp,cput)
        
        stpfrac=0.0d0
        dDisp=0.0d0
        ddDisp=0.0d0
        
        SEnergy = SEnergyConv
        
        Fext=Fext+Finc
        
        !------------------------------------------------------------------------------
        ! Iteration cycle
        !------------------------------------------------------------------------------
        iter = 0
        do while(iter <= itermax)
            iter = iter + 1
            
            !------------------------------------------------------------------------------
            ! Recover converged variables
            !------------------------------------------------------------------------------
            Stress = Stress_Conv
            Strain = Strain_Conv
            hard = hard_conv
            laxis = laxisconv
            SEnergy = SEnergyConv

            !------------------------------------------------------------------------------
            ! Elastic Step 
            !------------------------------------------------------------------------------
            nzep = 0
            if(iter == 1 .and. inc==1)then
            
                !------------------------------------------------------------------------------
                ! Patch cycle
                !------------------------------------------------------------------------------
                do iptch=1,nptch
                    
                    !------------------------------------------------------------------------------
                    ! Multipatch allocations
                    !------------------------------------------------------------------------------
                    if(nptch .gt. 1)then
                        call MPAlloc(iptch,nzep)
                        
                        allocate(Ke((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds))
                        Ke = 0.0d0
                        
                        allocate(dispin((p+1)*(q+1)*(w+1)*nds,1))
                        dispin = 0.0d0
                        
                        allocate(Finte((p+1)*(q+1)*(w+1)*nds,1))
                        Finte = 0.0d0
    
                    end if
                    
                    !------------------------------------------------------------------------------
                    ! Compute global coordinates
                    !------------------------------------------------------------------------------
                    if(nptch .gt. 1)then
                        nze = 0
                        do i=1, ncpx+p
                            do j=1,ncpy+q
                                do k=1,ncpz+w
                                
                                    if((u_knot(i) .ne. u_knot(i+1)) .and. (v_knot(j) .ne. v_knot(j+1))&
                                        & .and. (w_knot(k) .ne. w_knot(k+1))) then
                                        nzep = nzep + 1
                                        nze = nze + 1
                                        
                                        loc_num = 0
                                        
                                        do k1=1,(p+1)
                                            do k2=1,(q+1)
                                                do k3=1,(w+1)
                                                    loc_num = loc_num + 1
                                                    
                                                    do k4=1,nds
                                                        GCoords(MP_IEN(nzep,loc_num),k4) = Points(IEN(nze,loc_num),k4)
                                                    end do
                                                
                                                    GCoords(MP_IEN(nzep,loc_num),nds+1) = Weight(IEN(nze,loc_num))
                                                
                                                end do
                                            end do   
                                        end do
                                        
                                    endif
                                end do   
                            end do
                        end do
                    end if
                    
                    nzep = nzep - nze
                    NZE = 0
                    
                    !------------------------------------------------------------------------------
                    ! Element cycle
                    !------------------------------------------------------------------------------
                    GP_coords = 0.0d0                    
                    do i=1, ncpx+p
                        do j=1,ncpy+q
                            do k=1,ncpz+w
                            
                            !------------------------------------------------------------------------------
                            ! Check if it is a non-zero element
                            !------------------------------------------------------------------------------
                            if((u_knot(i) .ne. u_knot(i+1)) .and. (v_knot(j) .ne. v_knot(j+1)) .and.&
                                &  (w_knot(k) .ne. w_knot(k+1))) then
                                if(nptch .gt. 1)then
                                    nzep = nzep + 1
                                else
                                    nzep = nze + 1
                                end if
                                
                                nze = nze + 1
                                
                                SelectCase(Elcode)
                                    case('Hex8')
                                        call ElHex8(nze,nzep,Ke,dispin,Finte,inc)
!                                    case('Hex8BBar')
!                                        call ElHex8BBar(nze,Ke,dispin,Finte,MatA,MatC,MatV)
!                                    case('Hex8ANS')
!                                        call ElHex8ANS(nze,Ke,dispin,Finte,inc)
                                    case('Hex27')
                                        call ElHex27(nze,nzep,Ke,dispin,Finte,inc)
!                                    case('Hex27SRI')
!                                        call ElHex27SRI(nze,Ke,dispin,Finte)
!                                    case('Hex27EAS')
!                                        call ElHex27EAS(nze,Ke,dispin,Finte,inc)
                                    case('Hex27ANS')
                                        call ElHex27ANS(nze,nzep,Ke,dispin,Finte,inc)
                                    case('Hex27R')
                                        call ElHex27R(nze,nzep,Ke,dispin,Finte,inc)
!                                    case('Hex27PVS')
!                                        call ElHex27PVS(nze,Ke,dispin,Finte,inc)
!                                    case('Hex27PV')
!                                        call ElHex27PV(nze,Ke,dispin,Finte,inc)
!                                    case('Hex27BBar')
!                                        call ElHex27BBar(nze,Ke,dispin,Finte,MatA,MatC,MatV)
    !                                case('Hex27EAS_PW')
    !                                    call ElHex27EAS_PW(nze,Ke,dispin,Finte,MatA,MatC,MatV)
    !                                case('Hex27ANS_PW')
    !                                    call ElHex27ANS_PW(nze,Ke,dispin,Finte,MatA,MatC,MatV)
    !                                case('Hex27_Teste')
    !                                    call ElHex27_Teste(nze,Ke,dispin,Finte,inc)
!                                    case('HexRed')
!                                        call ElHexRed(nze,Ke,dispin,Finte,inc)
                                    case('Hex27_ProjVol')
                                        call ElHex27_ProjVol(nze,nzep,Ke,dispin,Finte,inc)
                                    case('Hex64')
                                        call ElHex64(nze,nzep,Ke,dispin,Finte,inc)
                                EndSelect
                                
                                !------------------------------------------------------------------------------
                                ! Assemble contributions to global matrices/arrays
                                !------------------------------------------------------------------------------
                                if(nptch == 1)then
                                    if(elcode=='Hex27BBar' .or. elcode=='Hex8BBar')then
                                        call AssemblyBBar(nze,MatA,MatC,MatV)
        !                            elseif(elcode=='Hex27EAS_PW')then
        !                                call AssemblyEAS(nze,MatA,MatC,MatV)
        !                            elseif(elcode=='Hex27ANS_PW')then
        !                                call AssemblyANS(nze,MatA,MatC,MatV)
                                    else
                                        call Assembly3D(nze,Ke,finte)
                                    end if
                                else
                                    
                                    call AssemblyMP(nze,nzep,Ke,Finte)
                                    
                                end if
                                
                                continue
                                
                            endif
                            
                            end do
                        end do
                    end do
                    
                    !------------------------------------------------------------------------------
                    ! Write Global coordinates at the first increment
                    !------------------------------------------------------------------------------
                    if(iter == 1 .and. inc == 1 .and. nptch .gt. 1)then
                        open(unit=8,File='GCoords.txt')
                        write(8,FMT=100)tnodes
                        
                        do i=1,tnodes
                            write(8,FMT=102)GCoords(i,1:nds+1)
                        end do
                        
                        100 format(I4,',')
                        101 format(3(I4,','))
                        102 format(100(E,','))
                        
                        close(8)
                    end if
                    
                    SEnergy = 0.0d0
                    
                    !------------------------------------------------------------------------------
                    ! Multipatch variables deallocation
                    !------------------------------------------------------------------------------
                    if(nptch .gt. 1)then

                        call MPDealloc(iptch)
                        
                        deallocate(Ke,dispin)
                        deallocate(Finte)
                        deallocate(INN,IEN)

                    else
                        !deallocate(Ke,KBar,dispin,Finte,coordi)
                    end if
                end do    
            
            end if !Elastic Step
                
            
            !Lumped approach for projected space methodologies -------------------------------------
!            temp1 = 0.0d0
!            do k1=1,(ncpx-1)*(ncpy-1)*(ncpz-1)
!                temp1(k1) = sum(MatVf(:,k1))
!            end do
!    
!            MatVf = 0.0d0
!            do k1=1,(ncpx-1)*(ncpy-1)*(ncpz-1)
!                MatVf(k1,k1) = temp1(k1)
!            end do
            
            !Inversion of the projected space stiffness matrix -------------------------------------
!            if(elcode=='Hex27BBar' .or. elcode=='Hex8BBar')then
!                temp1 = 0.0d0
!                call gaussj(MatVf,(ncpx-1)*(ncpy-1)*(ncpz-1),temp1,errorflag)
!            
!                do i=1,ncpx*ncpy*ncpz*nds
!                    do j=1,(ncpx-1)*(ncpy-1)*(ncpz-1)
!                        MatCfT(j,i) = MatCf(i,j)
!                    end do
!                end do
!            
!                Kf = MatAf + matmul(matmul(matCf,MatVf),MatCfT)
!            
!            elseif(elcode=='Hex27EAS_PW')then
!                temp1 = 0.0d0
!                
!                do j=1,nalpha*nelems
!                    diag(j) = sum(MatVf(:,j))
!                end do
!                
!                MatVf = 0.0d0
!                do j=1,nalpha*nelems
!                     MatVf(j,j) = 1.0d0/diag(j)
!                end do
!                
!                
!                !call gaussj(MatVf,nalpha,temp1,errorflag)
!            
!                do i=1,ncpx*ncpy*ncpz*nds
!                    do j=1,nalpha*nelems
!                        MatCfT(j,i) = MatCf(i,j)
!                    end do
!                end do
!            
!                Kf = MatAf - matmul(matmul(matCf,MatVf),MatCfT)
!                
!                continue
!                
!            elseif(elcode=='Hex27ANS_PW')then
!                temp1 = 0.0d0
!                call gaussj(MatVf,ncpx*ncpy*ncpz*nds,temp1,errorflag)
!            
!                do i=1,ncpx*ncpy*ncpz*nds
!                    do j=1,nalpha
!                        MatCfT(j,i) = MatCf(i,j)
!                    end do
!                end do
!            
!                Kf = MatAf + matmul(matmul(matCf,MatVf),MatCfT)
!                
!                continue
!            end if
            
            continue
            
!            if (Eigen==.true.)then
!                allocate(vcp(nnodes*nds,nnodes*nds))
!                vcp = 0.0d0
!                allocate(vlp(nnodes*nds))
!                vlp = 0.0d0
!                call VECP23(nnodes*nds,Kf,vlp,vcp,ierror)
!            end if
            
            !------------------------------------------------------------------------------
            ! Store full stiffness matrix
            !------------------------------------------------------------------------------
            KTot = Kf
            
            !------------------------------------------------------------------------------
            ! Apply zero-displacement boundary conditions
            !------------------------------------------------------------------------------
            if(nptch == 1)then
                count = 0
                do k=1,ncpz
                    do j=1,ncpy
                        do i=1,ncpx
                            do l=1,nds
                                count = count+1
                                if(BC(i,j,k,l)==0)then
                                    Kf(count,:) = 0.0d0
                                    Kf(:,count) = 0.0d0
                                    !Fint(count,1) = 0.0d0
                                    redcount(count) = 0
                                    continue
                                endif    
                            end do
                        end do
                    end do
                end do
            
            else
            
                count = 0
                do i=1,tnodes*nds
                    count = count + 1
                    if(BCdof(i) == 0) then
                        Kf(i,:) = 0.0d0
                        Kf(:,i) = 0.0d0
                        redcount(count) = 0
                    end if
                end do
            end if
    
            !------------------------------------------------------------------------------
            ! Apply prescribed displacements
            !------------------------------------------------------------------------------
            do i=1,tnodes*nds
                if(ImpDisp(i,1)/=0.0d0)then
                    do k1=1,tnodes*nds
                        if(i /= k1) then
                            !Fext(k1,1)=Fext(k1,1) - Kf(k1,i)*ImpDisp(i,1)/(1.0d0*incmax)*(1.0d0*inc)*(1.0d0-stpfrac(i))
                            Fext(k1,1)=Fext(k1,1) - Kf(k1,i)*ImpDisp(i,1)/(1.0d0*incmax)*(1.0d0-stpfrac(i))
                            Kf(i,k1)=0.0d0
                            Kf(k1,i)=0.0d0
                         else 
                            !Fext(i,1)=Kf(i,i)*ImpDisp(i,1)/(1.0d0*incmax)*(1.0d0*inc)*(1.0d0-stpfrac(i))
                            Fext(i,1)=Kf(i,i)*ImpDisp(i,1)/(1.0d0*incmax)*(1.0d0-stpfrac(i))  
                            continue
                        end if
                    end do
                end if
            end do
            
            
            if(PTS ==.true. .and. iter .gt. 1)then
                
                !------------------------------------------------------------------------------
                ! Point-to-segment contact module
                !------------------------------------------------------------------------------
                call PTS_ComputeGap()
            else
                
                !------------------------------------------------------------------------------
                ! Assemble equivalent reduced matrices
                !------------------------------------------------------------------------------
                k1 = 1
                k2 = 1
                Kfeq = 0.0d0
                FextEq = 0.0d0
                do i=1,tnodes*nds 
                    do j=1,tnodes*nds
                        if (redcount(i)/= 0 .and. redcount(j)/= 0)then
                            Kfeq(k1,k2)=Kf(i,j)
                            FextEq(k1,1)=Fext(i,1)                
                            k2=k2+1
                            if(k2==tnodes*nds-nbc+1)then
                                k1=k1+1
                                k2=1
                            end if
                        end if
                    end do

                    continue
                    
                
                end do
                
                !------------------------------------------------------------------------------
                ! Solve system of equations
                !------------------------------------------------------------------------------
                dispeqv = 0.0d0
                call Gauss (tnodes*nds-nbc,Kfeq,FextEq,dispeqv)

                continue

                j=1
                ddDisp = 0.0d0
                do i=1,tnodes*nds
                    if(redcount(i)==1)then
                        ddDisp(i,1)=dispeqv(j,1)
                        j=j+1
                    end if
                end do
            
                !------------------------------------------------------------------------------
                !Increment of the displacement increment
                !------------------------------------------------------------------------------
                dDisp = dDisp + ddDisp
                continue
            end if
            
            !Bypass for the B-Bar elements
!            if(elcode=='Hex27BBar' .or. elcode=='Hex8BBar')then
!                goto 999
!            end if
            
            !-----------------------------------------------------------------
            ! Compute reaction forces
            !-----------------------------------------------------------------
            count = 0
            do i=1,tnodes*nds
                count = count + 1
                if(BCdof(i) == 0) then
                    do j=1,tnodes*nds
                        Reac(i) = Reac(i) + KTot(i,j)*dddisp(j,1)
                    end do
                end if
            end do
            
            
            !goto 999
            continue

            nzep = 0
            Kf = 0.0d0
            Fint = 0.0d0
            SEnergy = 0.0d0
            MatAf = 0.0d0
            MatCf = 0.0d0
            MatVf = 0.0d0
            
            !------------------------------------------------------------------------------
            ! Patch cycle
            !------------------------------------------------------------------------------
            do iptch=1,nptch
                
                if(nptch .gt. 1)then
                    call MPAlloc(iptch,nzep)
                    
                    !Update b_net for output purpouses
                    call Updt_bnet(iptch,nzep)
                    
                    !Allocate elemental arrays
                    allocate(Ke((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds))
                    Ke = 0.0d0
                    
                    allocate(dispin((p+1)*(q+1)*(w+1)*nds,1))
                    dispin = 0.0d0
                    
                    allocate(Finte((p+1)*(q+1)*(w+1)*nds,1))
                    Finte = 0.0d0

                end if
                              
                if(nptch .gt. 1)then

                    nze = 0
                    count = (p+1)*(q+1)*(w+1) + 1
                    do i=1, ncpx+p
                        do j=1,ncpy+q
                            do k=1,ncpz+w
                            
                                !Check if it is a non-zero element
                                if((u_knot(i) .ne. u_knot(i+1)) .and. (v_knot(j) .ne. v_knot(j+1))&
                                        &  .and. (w_knot(k) .ne. w_knot(k+1))) then
                                    nzep = nzep + 1
                                    nze = nze + 1
                                    count = (p+1)*(q+1)*(w+1) + 1
                                    
                                    loc_num = 0
                                    do k1=1,p+1
                                        do k2=1,q+1
                                            do k3=1,w+1
                                                count = count - 1
                                                loc_num = loc_num + 1
                                                do k4=1,nds
                                                    GCoords(MP_IEN(nzep,loc_num),k4) = Points(IEN(nze,loc_num),k4)&
                                                        & + ddisp(MP_IEN(nzep,loc_num)*nds - nds + k4,1)
                                                end do
                                            end do
                                            GCoords(MP_IEN(nzep,loc_num),nds+1) = Weight(IEN(nze,loc_num))
                                        end do   
                                    end do
                                    
                                endif
                            end do   
                        end do
                    end do
                end if
                
                nzep = nzep - nze
                nze = 0            
                GP_coords = 0.0d0
                !------------------------------------------------------------------------------
                ! Element cycle
                !------------------------------------------------------------------------------
                do i=1, ncpx+p
                    do j=1,ncpy+q
                        do k=1,ncpz+w
                        
                            !Check if it is a non-zero element
                            if((u_knot(i) .ne. u_knot(i+1)) .and. (v_knot(j) .ne. v_knot(j+1))&
                                &  .and. (w_knot(k) .ne. w_knot(k+1))) then
                                
                                Finte = 0.0d0
                                
                                if(nptch .gt. 1)then
                                    nzep = nzep + 1
                                else
                                    nzep = nze + 1
                                end if
                                
                                nze = nze + 1
                                
                                if(nptch==1)then
                                    do k1=1,(p+1)*(q+1)*(w+1)
                                        dispin(k1*3-2,1) = dDisp(IEN(nze,k1)*3-2,1)
                                        dispin(k1*3-1,1) = dDisp(IEN(nze,k1)*3-1,1)
                                        dispin(k1*3  ,1) = dDisp(IEN(nze,k1)*3  ,1)
                                    end do
                                else
                                    do k1=1,(p+1)*(q+1)*(w+1)
                                        dispin(k1*3-2,1) = dDisp(MP_IEN(nzep,k1)*3-2,1)
                                        dispin(k1*3-1,1) = dDisp(MP_IEN(nzep,k1)*3-1,1)
                                        dispin(k1*3  ,1) = dDisp(MP_IEN(nzep,k1)*3  ,1)
                                    end do
                                end if
                                
                                SelectCase(Elcode)
                                    case('Hex8')
                                        call ElHex8(nze,nzep,Ke,dispin,Finte,inc)
!                                    case('Hex8BBar')
!                                        call ElHex8BBar(nze,Ke,dispin,Finte,MatA,MatC,MatV)
!                                    case('Hex8ANS')
!                                        call ElHex8ANS(nze,Ke,dispin,Finte,inc)
                                    case('Hex27')
                                        call ElHex27(nze,nzep,Ke,dispin,Finte,inc)
!                                    case('Hex27SRI')
!                                        call ElHex27SRI(nze,Ke,dispin,Finte)
!                                    case('Hex27EAS')
!                                        call ElHex27EAS(nze,Ke,dispin,Finte,inc)
!                                    case('Hex27PVS')
!                                        call ElHex27PVS(nze,Ke,dispin,Finte,inc)
!                                    case('Hex27PV')
!                                        call ElHex27PV(nze,Ke,dispin,Finte,inc)
                                    case('Hex27ANS')
                                        call ElHex27ANS(nze,nzep,Ke,dispin,Finte,inc)
                                    case('Hex27R')
                                        call ElHex27R(nze,nzep,Ke,dispin,Finte,inc)
!                                    case('Hex27BBar')
!                                        call ElHex27BBar(nze,Ke,dispin,Finte,MatA,MatC,MatV)
    !                                case('Hex27EAS_PW')
    !                                    call ElHex27EAS_PW(nze,Ke,dispin,Finte,MatA,MatC,MatV)
    !                                case('Hex27_Teste')
    !                                    call ElHex27_Teste(nze,Ke,dispin,Finte,inc)
!                                    case('HexRed')
!                                        call ElHexRed(nze,Ke,dispin,Finte,inc)
                                    case('Hex27_ProjVol')
                                        call ElHex27_ProjVol(nze,nzep,Ke,dispin,Finte,inc)
                                    case('Hex64')
                                        call ElHex64(nze,nzep,Ke,dispin,Finte,inc)
                                EndSelect
                                
                                !------------------------------------------------------------------------------
                                ! Assemble contributions to global matrices/arrays
                                !------------------------------------------------------------------------------
                                if(nptch == 1)then
                                    if(elcode=='Hex27BBar' .or. elcode=='Hex8BBar')then
                                        call AssemblyBBar(nze,MatA,MatC,MatV)
        !                            elseif(elcode=='Hex27EAS_PW')then
        !                                call AssemblyEAS(nze,MatA,MatC,MatV)
        !                            elseif(elcode=='Hex27ANS_PW')then
        !                                call AssemblyANS(nze,MatA,MatC,MatV)
                                    else
                                        call Assembly3D(nze,Ke,finte)
                                    end if
                                else
                                    
                                    call AssemblyMP(nze,nzep,Ke,Finte)
                                    
                                end if
                               
                                continue
                                
                            endif
                            
                        end do
                    end do
                end do
                
                if(nptch .gt. 1)then

                    call MPDealloc(iptch)
                    
                    deallocate(Ke,dispin)
                    deallocate(Finte)
                    deallocate(INN,IEN)

                else
                    !deallocate(Ke,KBar,dispin,Finte,coordi)
                end if
                
            end do
            
            !------------------------------------------------------------------------------
            !Eliminate reactions at fixed boundaries
            !------------------------------------------------------------------------------
            if(nptch==1)then
                count = 0
                do k=1,ncpz
                    do j=1,ncpy
                        do i=1,ncpx
                            do l=1,nds
                                count = count+1
                                if(BC(i,j,k,l)==0)then
                                    Fint(count,1) = 0.0d0
                                    continue
                                endif    
                            end do
                        end do
                    end do
                end do
            else
                do i=1,tnodes*nds
                    if(BCdof(i) == 0) then
                        Fint(i,1) = 0.0d0
                        if(PTS == .true.) FintLM(i,1) = 0.0d0
                    end if
                end do
            end if
            
            continue
            
            !------------------------------------------------------------------------------
            ! Eliminate reactions at points with prescribed displacement
            !------------------------------------------------------------------------------
            do i=1,tnodes*nds
                if(ImpDisp(i,1) /= 0.0d0)then
                    fint(i,1)=0.0d0
                    if(PTS == .true.)FintLM(i,1)=0.0d0
                end if
            end do
            
            !------------------------------------------------------------------------------
            ! Evaluate Residual Forces
            !------------------------------------------------------------------------------
            ! Residual = Fext(n)- Fint
            Res=0.0d0
            if(PTS == .false.)then
                Res = Fini/(1.0d0*incmax)*(1.0d0*inc) - Fint
            else
                Res = Fini/(1.0d0*incmax)*(1.0d0*inc) - Fint - FintLM
            end if
            
            continue
            
            !------------------------------------------------------------------------------
            ! Update prescribed displacements fraction
            !------------------------------------------------------------------------------
            do i=1,tnodes*nds
                if(ImpDisp(i,1) /= 0.0d0)then
                    stpfrac(i)=(u(i,1)+dDisp(i,1))/(ImpDisp(i,1)/(1.0d0*incmax)*(1.0d0*inc))
                    continue
                end if
            end do
            
            !------------------------------------------------------------------------------
            ! Fictional step to activate geometric non-linearity
            !------------------------------------------------------------------------------
            if(iter==1 .and. inc==1 .and. nlgeom==.true.) then
                Fext=Res
                continue
                goto 80
            end if
            
            !------------------------------------------------------------------------------
            ! Convergence Check
            !------------------------------------------------------------------------------
            absRes=0.0d0
            sumRes=0.0d0
            sumFext=0.0d0
            sumd=0.0d0
            sumdt=0.0d0
            do i=1,tnodes*nds
                sumRes = sumRes + Res(i,1)**2.0d0
                sumFext = sumFext + (Fini(i,1)/(1.0d0*incmax)*(1.0d0*inc))**2.0d0

                sumd= sumd + dddisp(i,1)**2.0d0
                sumdt= sumdt + ddisp(i,1)**2.0d0
            end do
            
            if(sumFext==0.0d0 .and. sumdt == 0.0d0 .and. sumres==0.0d0) then
                call ScreenData(7,inc,iter,AbsRes,absdisp,cput)
                absres = 0.0d0
                absdisp = 0.0d0
            
            elseif(sumFext==0.0d0)then
                !absres=0.0d0
                sumFext =1.0d0
                if(sumdt==0.0d0) then
                    absdisp=999.0d0
                else
                    absRes=sqrt(sumRes)/sqrt(sumFext)*100.0d0
                    absdisp=sqrt(sumd)/sqrt(sumdt)*100.0d0
                end if
            
            else
                absRes=sqrt(sumRes)/sqrt(sumFext)*100.0d0
                absdisp=sqrt(sumd)/sqrt(sumdt)*100.0d0
            end if

            call ScreenData(12,inc,iter,AbsRes,absdisp,cput)
            
            !-----------------------------------------------------------------
            ! Output relevant data
            !-----------------------------------------------------------------
            call MatlabOut()
            
            continue
            
            if(absRes .lt. RToler .and. absdisp .lt. RToler) then
                
                999 continue
                
                Fext=Res
                
                stpfrac = 0.0d0
                
                !-----------------------------------------------------------------
                ! Store converged variables
                !-----------------------------------------------------------------
                hard_conv = hard
                Stress_Conv = Stress
                Strain_Conv = Strain
                SEnergyConv = SEnergy
                laxisconv = laxis
        
                !-----------------------------------------------------------------
                ! Update displacement field
                !-----------------------------------------------------------------
                u = u + dDisp
                
                !-----------------------------------------------------------------
                ! Update CP's coordinates for geometry nonlinear analyses
                !-----------------------------------------------------------------
                if(nlgeom == .true.) then
                    if(nptch == 1)then
                        count = 0
                        do k=1,ncpz
                            do j=1,ncpy
                                do i=1,ncpx 
                                    count = count + 1
                                    Points(count,1) = Points(count,1) + dDisp(count*3-2,1)
                                    Points(count,2) = Points(count,2) + dDisp(count*3-1,1)
                                    Points(count,3) = Points(count,3) + dDisp(count*3  ,1)
                                end do
                            enddo
                        enddo
                    else
                        MP_B_Net = MP_B_Net_Final
                    end if
                end if
                
                !-----------------------------------------------------------------
                ! Write displacement field to output file
                !-----------------------------------------------------------------
                open(unit=9,file='Results.txt', access = 'APPEND')
                write(9,*)''
                write(9,*)''
                write(9,*)'----------------------------'
                write(9,*)'INCREMENT ', inc
                write(9,*)'----------------------------'
                write(9,*)''
                write(9,*)'Displacement'
                write(9,*)''
                do i=1, tnodes
                    write(9,11)i,u(i*nds-2,1),u(i*nds-1,1),u(i*nds,1)
                end do
                
                !-----------------------------------------------------------------
                ! Write stress field to output file
                !-----------------------------------------------------------------
                write(9,*)''
                write(9,*)'Stress'
                write(9,*)''
                do i=1, nze
                    do j=1, npi
                        write(9,12)i,j,Stress(i,j,:)
                    end do
                end do
                
                !-----------------------------------------------------------------
                ! Write strain field to output file
                !-----------------------------------------------------------------
                write(9,*)''
                write(9,*)'Strain'
                write(9,*)''
                do i=1, nze
                    do j=1, npi
                        write(9,12)i,j,Strain(i,j,:)
                    end do
                end do
                
                !-----------------------------------------------------------------
                ! Write reaction forces to output file
                !-----------------------------------------------------------------
                open(unit=7,file='ReactionForces.txt', access = 'APPEND')
                write(7,*)''
                write(7,*)''
                write(7,*)'----------------------------'
                write(7,*)'INCREMENT ', inc
                write(7,*)'----------------------------'
                write(7,*)''
                write(7,*)'Reaction Forces'
                write(7,*)''
                
                tst1 = 0.0d0
                tst2 = 0.0d0
                tst3 = 0.0d0
                do i=1,tnodes
                    tst1 = tst1 + Reac(i*nds-2)
                    tst2 = tst2 + Reac(i*nds-1)
                    tst3 = tst3 + Reac(i*nds  )
                    write(7,11)i,Reac(i*nds-2),Reac(i*nds-1),Reac(i*nds)
                end do
                
                write(7,FMT=13)tst1, tst2, tst3
                
                
                11 format(1(I5,1x),6(E,1x))
                12 format(2(I5,1x),6(E,1x))
                13 format('Sum ',6(E,1x))
                close(9)
                close(7)
                
                !-----------------------------------------------------------------
                ! Write Gaus Points coordinates
                !-----------------------------------------------------------------
                open(unit=7,file='GPCoords.txt', access = 'APPEND')
                write(7,*)''
                write(7,*)''
                write(7,*)'----------------------------'
                write(7,*)'INCREMENT ', inc
                write(7,*)'----------------------------'
                write(7,*)''
                write(7,*)'Gauss Points Coordinates'
                write(7,*)''
                
                count=0                
                do i=1, nelems                
                        do j=1,npi
                                count = count + 1                               
                                write(7,FMT=14)i, j, GP_coords(count,1), GP_coords(count,2), GP_coords(count,3)
                        end do
                end do
                14 format(2(I4,','),3(E,','))
                close(7)
                
                call MatlabOut()
                
                call ScreenData(8,inc,iter,AbsRes,absdisp,cput)
                
                !pause
                
                continue
                    
                goto 90
            else
                
                !-----------------------------------------------------------------
                ! If maximum iteration number is reached terminate the analysis
                !-----------------------------------------------------------------
                Fext=Res
                if(iter == itermax) then
                    call ScreenData(9,inc,iter,AbsRes,absdisp,cput)
                    goto 70
                end if    
            end if  
            
         80 continue
            
        end do !iteration cycle
    
    90 continue
    
    end do !increment cycle 
    
    call ScreenData(10,inc,iter,AbsRes,absdisp,cput)
    
    70 continue
    
    !-----------------------------------------------------------------
    ! CPU time of the analysis
    !-----------------------------------------------------------------
    !CPUt=etime(atime)
    call ScreenData(11,inc,iter,AbsRes,absdisp,cput)

    continue
    
    
    end program ICO3D

