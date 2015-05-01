    program ICO3D
    
    !use ifport
    use Mod_Variables
    
    implicit none
    
    character*256::FileName
    integer(4)::i,j,k,l,nze,k1,k2,kt,errorflag
    integer(4)::count,iter,inc
    integer(4)::isolver
    integer(4)::nred
    
    !Iterative Solver Variables ----
    real(8),dimension(:),allocatable::r1,r2,vs,ws,ys
    external::aprod,msolve
    integer(4)::istop,itn
    real(8)::Anorm,Rnorm,ynorm,Acond,STol
    integer(4)::nout,itnlim
    real(8)::shift
    logical::checkA,goodb,precon
    
    !Variables for solver3
    real(8),dimension(:,:),allocatable::coef
    real(8),dimension(:),allocatable::AR
    integer(4),dimension(:),allocatable::IA,JA
    real(8),dimension(:),allocatable::wksp,rparm
    integer(4),dimension(:),allocatable::colnz,iwksp,iparm
    integer(4)::mnz,countx,county
    
    
    real(8)::da,tst
    
    real(8), dimension(:,:), allocatable::Ke,dispin,finte
    real(8)::absRes,absDisp,sumRes,sumFext,sumd,sumdt
    
    real(8), dimension(:,:), allocatable::MatA,MatC,MatV
    real(8), dimension(:), allocatable::temp1
    
    real(8),dimension(:,:), allocatable::coordi
    real(8),dimension(:,:), allocatable::MMult
    
    real(4)::CPUt
    real(4),dimension(2)::atime
    
    
    write(*,*)'!------------------------------------------------------------------------------'
    write(*,*)'!                                                                             !'
    write(*,*)'!      %%%%%%%%%%%%%%%%   %%%%%%%%%%%%   %%%%%%%%%%%%%%%  %%%%%%%%%  %%%%     !'
    write(*,*)'!            %%          %%             %%           %%         %%  %%  %%    !'
    write(*,*)'!           %%          %%             %%           %%         %%  %%    %%   !'
    write(*,*)'!          %%          %%             %%           %%   %%%%%%%%  %%      %%  !'
    write(*,*)'!         %%          %%             %%           %%         %%  %%      %%   !'
    write(*,*)'!        %%          %%             %%           %%         %%  %%    %%      !'
    write(*,*)'!       %%          %%             %%           %%         %%  %%  %%         !'
    write(*,*)'! %%%%%%%%%%%%%%   %%%%%%%%%%%%   %%%%%%%%%%%%%%%  %%%%%%%%%  %%%%            !'
    write(*,*)'!                                                                             !'
    write(*,*)'!------------------------------------------------------------------------------'
    write(*,*)'!                                                                             !'
    write(*,*)'!                   Isogeometric COde (ICO) for 3D applications               !'
    write(*,*)'!                                                                             !'
    write(*,*)'!------------------------------------------------------------------------------' 
    
    write(*,*)''
    write(*,*)'Input file name (including file extension): '

    !read(*,*)FileName
    
    !FileName = 'BernoulliStrBeam_16el_p2q2w2'
    !FileName = 'CCP_16el_p2q2w2'
    !FileName = 'CookPlast_8el_p1q1w1'
    !FileName = 'NLG_MBending_10el_p2q2w2'
    !FileName = 'SLo_8el_p2q2w2'
    !FileName = 'CylShellStrip_32el_p2q2w2'
    !FileName = 'CylShellStrip_ML5_p2q2w2'
    !FileName = 'TwoHex'
    !FileName = 'BCs'
    !FileName = 'CSP_VL2_p2q2w2'
    !FileName = 'MeshDist_e03_p2q2w2'
    FileName = 'FullSphere_4el_p2q2w2'
    !FileName = 'CurvedBeam_16el_p2q2w2'
    !FileName = 'SLo_2el_p2q2w2'
    
    !FileName = 'PressurePlate_8el_p2q2w2'
    !FileName = 'SingleHex1'
    
    itermax = 25
    incmax = 5
    nlgeom=.false.
    
    eigen = .false.
    
    isolver = 1
    
!    open(unit=999,file='Basis.txt')
!    close(999)
    
    write(*,*)''
    write(*,*)'Reading Input File'
    
    call ReadInFile(FileName)
    
    write(*,*)'...'
    write(*,*)'Done'
    
    !Obtain the weighted coordinates for the control points
!    Bw = B_NET
!    do i = 1,nds
!        Bw(:,:,:,i) = Bw(:,:,:,i)*B_net(:,:,:,nds+1)
!    enddo
    
    write(*,*)''
    write(*,*)'Generating Data Structure'
    
    open(unit=9,file='Results.txt', status = 'REPLACE')
    write(9,*)''
    write(9,*)'ICO3D OUTPUT FILE'
    write(9,*)''
    close(9)
    
    allocate(Ke((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds))
    Ke = 0.0d0
    allocate(dispin((p+1)*(q+1)*(w+1)*nds,1))
    dispin = 0.0d0
    allocate(Finte((p+1)*(q+1)*(w+1)*nds,1))
    
    allocate(coordi(ncpx*ncpy*ncpz,nds))
    !coordi = 0.0d0
    
    if(elcode == 'Hex8BBar' .or. elcode == 'Hex27BBar')then
        allocate(MatA((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds))
        MatA = 0.0d0
        allocate(MatC((p+1)*(q+1)*(w+1)*nds,p*q*w))
        MatC = 0.0d0
        allocate(MatV(p*q*w,p*q*w))
        MatV = 0.0d0
        allocate(temp1((ncpx-1)*(ncpy-1)*(ncpz-1)))
        Temp1 = 0.0d0
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
    end if
    
    !Generate connectivity arrays
    call gen_ien_inn()
    call gen_ien_inn_bar()
    
    !Generate control polygon in GiD
    call GiDMesh()
    
!    !Generate equivalent FEM mesh
!    call GiDEqMesh(coordi)
    
    !Assemble the external forces vector
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
    
    Fini = Fext
    Fext = 0.0d0
    Finc=Fini/(1.0d0*incmax)
    
    SEnergyConv = 0.0d0
    
    write(*,*)'...'
    write(*,*)'Done'
    
    write(*,*)''
    write(*,*)'!------------------------------------------------------------------------------'
    write(*,*)'!                               Begin Analysis                                !'
    write(*,*)'!------------------------------------------------------------------------------' 
    
    !Begin incremental cycle
    inc = 0
    do inc=1,incmax
        
        write(*,*)''
        write(*,*)'!------------------------------------------------------------------------------'
        write(*,FMT=60)inc
        write(*,*)'!------------------------------------------------------------------------------' 
        60 format(' ! Increment ', I2)
        
        stpfrac=0.0d0
        dDisp=0.0d0
        ddDisp=0.0d0
        
        SEnergy = SEnergyConv
        
        Fext=Fext+Finc
        
        iter = 0
        do while(iter <= itermax)
            iter = iter + 1
            
            !Use only converged variables in the iteration cycle
            Stress = Stress_Conv
            Strain = Strain_Conv
            hard = hard_conv
            laxis = laxisconv
            SEnergy = SEnergyConv
            !alpha=alphaconv
            !StoreC = StoreCConv
            
            !Elastic Step ---------------------------------------------
            if(inc == 1 .and. iter == 1)then
                NZE = 0
                do i=1, ncpx+p
                    do j=1,ncpy+q
                        do k=1,ncpz+w
                        
                        !Check if it is a non-zero element
                        if((u_knot(i) .ne. u_knot(i+1)) .and. (v_knot(j) .ne. v_knot(j+1)) .and. (w_knot(k) .ne. w_knot(k+1))) then
                            nze = nze + 1
                            
                            SelectCase(Elcode)
                                case('Hex8')
                                    call ElHex8(nze,Ke,dispin,Finte,inc)
                                case('Hex8BBar')
                                    call ElHex8BBar(nze,Ke,dispin,Finte,MatA,MatC,MatV)
                                case('Hex8ANS')
                                    call ElHex8ANS(nze,Ke,dispin,Finte,inc)
                                case('Hex27')
                                    call ElHex27(nze,Ke,dispin,Finte,inc)
                                case('Hex27SRI')
                                    call ElHex27SRI(nze,Ke,dispin,Finte)
                                case('Hex27EAS')
                                    call ElHex27EAS(nze,Ke,dispin,Finte,inc)
                                case('Hex27ANS')
                                    call ElHex27ANS(nze,Ke,dispin,Finte,inc)
                                case('Hex27R')
                                    call ElHex27R(nze,Ke,dispin,Finte,inc)
                                case('Hex27PVS')
                                    call ElHex27PVS(nze,Ke,dispin,Finte,inc)
                                case('Hex27PV')
                                    call ElHex27PV(nze,Ke,dispin,Finte,inc)
                                case('Hex27BBar')
                                    call ElHex27BBar(nze,Ke,dispin,Finte,MatA,MatC,MatV)
!                                case('Hex27EAS_PW')
!                                    call ElHex27EAS_PW(nze,Ke,dispin,Finte,MatA,MatC,MatV)
!                                case('Hex27ANS_PW')
!                                    call ElHex27ANS_PW(nze,Ke,dispin,Finte,MatA,MatC,MatV)
!                                case('Hex27_Teste')
!                                    call ElHex27_Teste(nze,Ke,dispin,Finte,inc)
                                case('HexRed')
                                    call ElHexRed(nze,Ke,dispin,Finte,inc)
                                case('Hex27_ProjVol')
                                    call ElHex27_ProjVol(nze,Ke,dispin,Finte,inc)
                                case('Hex64')
                                    call ElHex64(nze,Ke,dispin,Finte,inc)
                            EndSelect
                            
                            if(elcode=='Hex27BBar' .or. elcode=='Hex8BBar')then
                                call AssemblyBBar(nze,MatA,MatC,MatV)
!                            elseif(elcode=='Hex27EAS_PW')then
!                                call AssemblyEAS(nze,MatA,MatC,MatV)
!                            elseif(elcode=='Hex27ANS_PW')then
!                                call AssemblyANS(nze,MatA,MatC,MatV)
                            else
                                call Assembly3D(nze,Ke,finte)
                            end if
                            
                            continue
                            
                        endif
                        
                        end do
                    end do
                end do
                SEnergy = 0.0d0
                
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
            if(elcode=='Hex27BBar' .or. elcode=='Hex8BBar')then
                temp1 = 0.0d0
                call gaussj(MatVf,(ncpx-1)*(ncpy-1)*(ncpz-1),temp1,errorflag)
            
                do i=1,ncpx*ncpy*ncpz*nds
                    do j=1,(ncpx-1)*(ncpy-1)*(ncpz-1)
                        MatCfT(j,i) = MatCf(i,j)
                    end do
                end do
            
                Kf = MatAf + matmul(matmul(matCf,MatVf),MatCfT)
            
            elseif(elcode=='Hex27EAS_PW')then
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
                
            elseif(elcode=='Hex27ANS_PW')then
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
            end if
            
            continue
            
!            if (Eigen==.true.)then
!                allocate(vcp(nnodes*nds,nnodes*nds))
!                vcp = 0.0d0
!                allocate(vlp(nnodes*nds))
!                vlp = 0.0d0
!                call VECP23(nnodes*nds,Kf,vlp,vcp,ierror)
!            end if
            
            !Apply zero-displacement boundary conditions
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
            
            do i=1,nnodes*nds
                if(ImpDisp(i,1)/=0.0d0)then
                    do k1=1,nnodes*nds
                        if(i /= k1) then
                            Fext(k1,1)=Fext(k1,1) - Kf(k1,i)*ImpDisp(i,1)/(1.0d0*incmax)*(1.0d0*inc)*(1.0d0-stpfrac(i))
                            Kf(i,k1)=0.0d0
                            Kf(k1,i)=0.0d0
                         else 
                            Fext(i,1)=Kf(i,i)*ImpDisp(i,1)/(1.0d0*incmax)*(1.0d0*inc)*(1.0d0-stpfrac(i))  
                            continue
                        end if
                    end do
                end if
            end do
            
            
            !Assemble equivalent reduced matrices
            k1 = 1
            k2 = 1
            do i=1,nnodes*nds 
                do j=1,nnodes*nds
                    if (redcount(i)/= 0 .and. redcount(j)/= 0)then
                        Kfeq(k1,k2)=Kf(i,j)
                        FextEq(k1,1)=Fext(i,1)                
                        k2=k2+1
                        if(k2==nnodes*nds-nbc+1)then
                            k1=k1+1
                            k2=1
                        end if
                    end if
                end do

                continue
            
            end do
            
            !Solve system of equations
            if(isolver==1)then
                dispeqv = 0.0d0
                call Gauss (nnodes*nds-nbc,Kfeq,FextEq,dispeqv)
                
                continue
                
            elseif(isolver==2)then
                
                nred = nnodes*nds-nbc
                if(inc==1 .and. iter==1)then
                    allocate(r1(nred),r2(nred),vs(nred),ws(nred),ys(nred))
                    allocate(diag(nred))
                end if
                
                r1 = 0.0d0
                r2 = 0.0d0
                vs = 0.0d0
                ws = 0.0d0
                ys = 0.0d0
                
                runner = 0
                nlin = nred
                ncol = nred
                do i=1,nred
                    do j=1,nred
                        if(Kfeq(i,j).ne.0.0d0)runner = runner + 1
                    end do
                end do
                123 continue
                
                if(inc==1 .and. iter==1)allocate(indx(runner,2))
                
                indx = 0
                
                count = 0
                do i=1,nred
                    do j=1,nred
                        if(Kfeq(i,j).ne.0.0d0)then
                            count = count + 1
                            indx(count,1) = i
                            indx(count,2) = j
                        end if
                    end do
                end do
                
!                j=1
!                ddDisp = 0.0d0
!                do i=1,nnodes*nds
!                    if(redcount(i)==1)then
!                        dispeqv(j,1) = ddDisp(i,1)
!                        j=j+1
!                    end if
!                end do
                
                checkA = .false.
                goodb = .false.
                precon = .false.
                shift = 0.0d0
                nout = 0
                itnlim = 5000
                
                STol = 1.0d-10
                
                dispeqv = 0.0d0
                
                call SYMMLQ(nred,FextEq,r1,r2,vs,ws,dispeqv,ys,&  
                &           Aprod,Msolve,checkA,goodb, &
                &           precon,shift,nout,itnlim,STol,&
                &           istop, itn, Anorm, Acond, rnorm, ynorm)
                
                continue
            elseif(isolver==3)then
                
                nred = nnodes*nds-nbc
                
                if(inc==1 .and. iter==1)then
                    allocate(colnz(nred))
                    allocate(iwksp(nred))
                    iwksp = 0
                    allocate(wksp(4*nred+2*500*nred))
                    wksp = 0.0d0
                    allocate(iparm(12))
                    allocate(rparm(12))
                end if
                
                colnz = 0
                do k1=1,nred
                    do k2=1,nred
                        if(kfeq(k1,k2) .ne. 0.0d0) colnz(k1) = colnz(k1) + 1
                    end do
                end do
                
                mnz = 0
                do k1=1,nred
                    mnz = mnz + colnz(k1)
                end do
                
                if(inc==1 .and. iter==1)then
                    allocate(Ar(mnz))
                    allocate(IA(nred))
                    allocate(JA(mnz))
                end if
                
                countx = 0
                county = 1
                IA(1) = 1
                do k1=1,nred
                    
                    if (k1 .gt. 1)then
                        county = county + 1
                        IA(county) = countx+1
                    end if
                    
                    do k2=1,nred
                        if(kfeq(k1,k2) .ne. 0.0d0) then
                            countx = countx + 1
                            Ar(countx) = kfeq(k1,k2)
                            JA(countx) = k2
                        end if
                    end do
                    
                end do
                
                iparm(1) = 100
                iparm(2) = 0
                iparm(3) = 0
                iparm(4) = 6
                iparm(5) = 0
                iparm(6) = 1
                iparm(7) = 1
                iparm(8) = 900
                iparm(9) = -1
                iparm(10) = 0
                iparm(11) = 0
                iparm(12) = 0
                
                rparm(1) = 1.0d-6
                rparm(2) = 0.0d0
                rparm(3) = 0.0d0
                rparm(4) = 0.75d0
                rparm(5) = 1.0d0
                rparm(6) = 0.0d0
                rparm(7) = 0.25d0
                rparm(8) = 100.0d0/16.0d0
                rparm(9) = 0.0d0
                rparm(10) = 0.0d0
                rparm(11) = 0.0d0
                rparm(12) =  0.0d0
                
            end if

            j=1
            ddDisp = 0.0d0
            do i=1,nnodes*nds
                if(redcount(i)==1)then
                    ddDisp(i,1)=dispeqv(j,1)
                    j=j+1
                end if
            end do
            
            !Increment in the displacement increment ------ 
            dDisp = dDisp + ddDisp
            continue
            
            !Bypass for the B-Bar elements
            if(elcode=='Hex27BBar' .or. elcode=='Hex8BBar')then
                goto 999
            end if
            
            !goto 999
            continue

            NZE = 0
            Kf = 0.0d0
            Fint = 0.0d0
            SEnergy = 0.0d0
            MatAf = 0.0d0
            MatCf = 0.0d0
            MatVf = 0.0d0
            
            do i=1, ncpx+p
                do j=1,ncpy+q
                    do k=1,ncpz+w
                    
                        !Check if it is a non-zero element
                        if((u_knot(i) .ne. u_knot(i+1)) .and. (v_knot(j) .ne. v_knot(j+1)) .and. (w_knot(k) .ne. w_knot(k+1))) then
                            nze = nze + 1
                            
                            do k1=1,(p+1)*(q+1)*(w+1)
                                dispin(k1*3-2,1) = dDisp(IEN(nze,k1)*3-2,1)
                                dispin(k1*3-1,1) = dDisp(IEN(nze,k1)*3-1,1)
                                dispin(k1*3  ,1) = dDisp(IEN(nze,k1)*3  ,1)
                            end do
                            
                            SelectCase(Elcode)
                                case('Hex8')
                                    call ElHex8(nze,Ke,dispin,Finte,inc)
                                case('Hex8BBar')
                                    call ElHex8BBar(nze,Ke,dispin,Finte,MatA,MatC,MatV)
                                case('Hex8ANS')
                                    call ElHex8ANS(nze,Ke,dispin,Finte,inc)
                                case('Hex27')
                                    call ElHex27(nze,Ke,dispin,Finte,inc)
                                case('Hex27EAS')
                                    call ElHex27EAS(nze,Ke,dispin,Finte,inc)
                                case('Hex27PVS')
                                    call ElHex27PVS(nze,Ke,dispin,Finte,inc)
                                case('Hex27PV')
                                    call ElHex27PV(nze,Ke,dispin,Finte,inc)
                                case('Hex27ANS')
                                    call ElHex27ANS(nze,Ke,dispin,Finte,inc)
                                case('Hex27R')
                                    call ElHex27R(nze,Ke,dispin,Finte,inc)
                                case('Hex27BBar')
                                    call ElHex27BBar(nze,Ke,dispin,Finte,MatA,MatC,MatV)
!                                case('Hex27EAS_PW')
!                                    call ElHex27EAS_PW(nze,Ke,dispin,Finte,MatA,MatC,MatV)
!                                case('Hex27_Teste')
!                                    call ElHex27_Teste(nze,Ke,dispin,Finte,inc)
                                case('HexRed')
                                    call ElHexRed(nze,Ke,dispin,Finte,inc)
                                case('Hex27_ProjVol')
                                    call ElHex27_ProjVol(nze,Ke,dispin,Finte,inc)
                                case('Hex64')
                                    call ElHex64(nze,Ke,dispin,Finte,inc)
                            EndSelect
                            
                            if(elcode=='Hex27BBar' .or. elcode=='Hex8BBar')then
                                call AssemblyBBar(nze,MatA,MatC,MatV)
!                            elseif(elcode=='Hex27EAS_PW')then
!                                call AssemblyEAS(nze,MatA,MatC,MatV)
!                            elseif(elcode=='Hex27ANS_PW')then
!                                call AssemblyANS(nze,MatA,MatC,MatV)
                            else
                                call Assembly3D(nze,Ke,finte)
                            end if
                           
                            continue
                            
                        endif
                        
                    end do
                end do
            end do

            !Eliminate reactions at fixed boundaries
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
            
            continue
            
            !Eliminate reactions at points with prescribed displacement
            do i=1,nnodes*nds
                if(ImpDisp(i,1) /= 0.0d0)then
                    fint(i,1)=0.0d0
                end if 
            end do
            
            
            !Evaluate Residual Forces
            !Residual = Fext(n)- Fint
            Res=0.0d0
            Res=Fini/(1.0d0*incmax)*(1.0d0*inc)-Fint 
            
            continue
            
            !Fictional step to activate geometric non-linearity ------------------------------------
            if(iter==1 .and. inc==1 .and. nlgeom==.true.) then
                Fext=Res
                continue
                goto 80
            end if
            
            !Convergence Check ---------------------------------------------------------------------
            absRes=0.0d0
            sumRes=0.0d0
            sumFext=0.0d0
            sumd=0.0d0
            sumdt=0.0d0
            do i=1,ncpx*ncpy*ncpz*nds
                sumRes = sumRes + Res(i,1)**2.0d0
                sumFext = sumFext + (Fini(i,1)/(1.0d0*incmax)*(1.0d0*inc))**2.0d0

                sumd= sumd + dddisp(i,1)**2.0d0
                sumdt= sumdt + ddisp(i,1)**2.0d0
            end do
            
            if(sumFext==0.0d0 .and. sumdt == 0.0d0 .and. sumres==0.0d0) then
                write(*,*)'No load or displacement applied to the model'
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
            
            do i=1,ncpx*ncpy*ncpz*nds
                if(ImpDisp(i,1) /= 0.0d0)then
                    stpfrac(i)=(u(i,1)+dDisp(i,1))/(ImpDisp(i,1)/(1.0d0*incmax)*(1.0d0*inc))
                    continue
                end if
            end do
            
            !tst = sum(res(:,1))
            
            continue
            
            if(iter==20)then
                continue
            end if
            
!            if(iter==2) then
!                absRes = 0.0d0
!                absdisp = 0.0d0
!            end if
            
            if(absRes .lt. RToler .and. absdisp .lt. RToler) then
                
                999 continue
                
                Fext=Res
                
                stpfrac = 0.0d0
                
                !Residual after convergence -----
                !Resconv=Res
                hard_conv = hard
                
                !Update converged stress
                !Stress_Conv = Stress_Conv + dStress
                !Strain_Conv = Strain_Conv + dStrain
                
                Stress_Conv = Stress
                Strain_Conv = Strain
        
                !Update converged displacement -----
                u = u + dDisp
                
                !Update Enhanced Variables
                !alphaconv=alpha
                
                !Update Strain energy -----
                SEnergyConv = SEnergy
                
                !Update Control Points Coordinates
                if(nlgeom == .true.) then
                    
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
                end if
                
                open(unit=9,file='Results.txt', access = 'APPEND')
                write(9,*)''
                write(9,*)''
                write(9,*)'----------------------------'
                write(9,*)'INCREMENT ', inc
                write(9,*)'----------------------------'
                write(9,*)''
                write(9,*)'Stress'
                write(9,*)''
                do i=1, nze
                    do j=1, npi
                        write(9,12)i,j,Stress(i,j,:)
                    end do
                end do
                
                write(9,*)''
                write(9,*)'Strain'
                write(9,*)''
                do i=1, nze
                    do j=1, npi
                        write(9,12)i,j,Strain(i,j,:)
                    end do
                end do
                12 format(2(I5,1x),6(E,1x))
                close(9)
                
                !Update converged local axis -------
                laxisconv = laxis
                
                !Update right Cauchy-Green tensor
                !StoreCConv = StoreC
                
                !Output displacement in the control polygon
                call GiDRes(inc)
                
                !Output stress to equivalent FEM mesh
                !call GiDEqRes(inc,coordi)
                
                write(*,*)''
                write(*,*)' ! Solution has converged at iteration', iter
                write(*,FMT=65)AbsRes
                write(*,FMT=66)AbsDisp
             65 format('  -> The absolute residual of the forces is R=        ', E)
             66 format('  -> The absolute residual of the displacement is Rd= ', E)
                write(*,*)''
                
                continue
                    
                goto 90
            else
                Fext=Res
                continue
                if(iter == itermax) then
                
                    write(*,*)''
                    write(*,*)' ! Solution has not converged after iteration', iter
                    write(*,*)' ! Teminating analysis'
                    write(*,*)'!------------------------------------------------------------------------------'
                    
                    goto 70
                    
                end if    
                !goto 10
            end if  
            
         80 continue
            
        end do !iteration cycle
    
    90 continue
    
    end do !increment cycle 
    
    write(*,*)''
    write(*,*)'!------------------------------------------------------------------------------'
    write(*,*)'!                              Analysis Complete                              !'
    write(*,*)'!------------------------------------------------------------------------------'
    
    70 continue
    
    !CPUt=etime(atime)
    
    continue
    
    
    end program ICO3D

