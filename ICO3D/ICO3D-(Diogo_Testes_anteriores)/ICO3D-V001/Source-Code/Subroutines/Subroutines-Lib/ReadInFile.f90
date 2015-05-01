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
    
    integer(4)::i,j,k,l,k1,k2,nnod,count
    integer(4)::int1,int2,int3,int4,nload
    
    logical::pstrain, pstress
    !character::elcode
    
    real(8)::temp1
    real(8),allocatable,dimension(:,:)::matD
    
    integer(4),allocatable,dimension(:)::bcdof
    
    
    open(unit=1,file=FileName)
    
    !Logical allocation
    pstrain = .false.
    pstress = .false.
    gravity = .false.
    pressure = .false.
    
    
    do while(Line .ne. '*end')
        !write(*,*)Line
        read(1,*)Line
        
        if(line=='*begin')then
            read(1,*)nds                    !Number of spatial directions
            read(1,*)p,q,w         !Degrees of curves in u,v,w
            read(1,*)ncpx,ncpy,ncpz     !Number of control points in each direction
            read(1,*)closed_u, closed_v, closed_w !Cloded knot vectors - 1: Yes, 0: No:
        
            !-----------------------------------------------
            ! Allocate Variables with the ReadInFile(*begin) 
            !-----------------------------------------------
            nnodes = ncpx*ncpy*ncpz              !Total number of nodes
            nelems = (ncpx-p)*(ncpy-q)*(ncpz-w)  !Total number of elements
            nshpl = (p+1)*(q+1)*(w+1)            !Number of local shape functions
            
!                integer(4),allocatable::u_knot(:),v_knot(:),b_net(:,:,:)
!                integer(4),allocatable::ibc(:,:),iper(:)
!                real(8),allocatable::dir_bc(:,:),ld_el(:,:)
            
            allocate(u_knot(ncpx+p+1),v_knot(ncpy+q+1),w_knot(ncpz+w+1))
            u_knot = 0.0d0
            v_knot = 0.0d0
            w_knot = 0.0d0
            
            allocate(b_net(ncpx,ncpy,ncpz,nds+1))
            b_net = 0.0d0
            
            allocate(b_net_final(ncpx,ncpy,ncpz,nds+1))
            b_net_final = 0.0d0
            
            !allocate(bw(ncpx,ncpy,ncpz,nds+1))
            !bw = 0.0d0
  
            allocate(iper(nnodes))
            do k1 = 1, nnodes
                iper(k1) = k1
            enddo
            
            allocate(bc(ncpx,ncpy,ncpz,nds),load(ncpx,ncpy,ncpz,nds))
            bc = 1
            load = 0.0d0
            
            allocate(dispBC(ncpx,ncpy,ncpz,nds))
            dispBC = 0.0d0

            allocate(ld_el(nelems,nds))
            ld_el = 0.0d0
            allocate(ey(nelems))
            ey = 0.0d0
            allocate(nu(nelems))
            nu = 0.0d0
            
            !allocate(Ke((p+1)*(q+1)*nds,(p+1)*(q+1)*nds)
            allocate(Kf(ncpx*ncpy*ncpz*nds,ncpx*ncpy*ncpz*nds))
            Kf = 0.0d0
            
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
            
            allocate(PMesh(ncpx*ncpy*ncpz,nds))
            PMesh = 0.0d0
            
            allocate(Weight(ncpx*ncpy*ncpz))
            Weight = 0.0d0
            
            allocate(Points(ncpx*ncpy*ncpz,nds))
            Points = 0.0d0
            
            allocate(loaddof(nnodes*nds))
            loaddof = 0.0d0
            
            allocate(dispdof(nnodes*nds))
            dispdof = 0.0d0
            
            ndi = nds
            nshr = 3
            ntens = ndi + nshr
                                
        
        elseif(line=='*knots')then
            read (1,*) (u_knot(i),i=1,ncpx+p+1)
            read (1,*) (v_knot(i),i=1,ncpy+q+1)
            read (1,*) (w_knot(i),i=1,ncpz+w+1)
        
        elseif(line=='*element')then
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
            
            !Allocations for projection method ---------------------------------------
            if(elcode=='Hex8BBar' .or. elcode=='Hex27BBar')then
                allocate(MatAf(ncpx*ncpy*ncpz*nds,ncpx*ncpy*ncpz*nds))
                MatAf = 0.0d0
                allocate(MatCf(ncpx*ncpy*ncpz*nds,(ncpx-1)*(ncpy-1)*(ncpz-1)))
                MatCf = 0.0d0
                allocate(MatCfT((ncpx-1)*(ncpy-1)*(ncpz-1),ncpx*ncpy*ncpz*nds))
                MatCfT = 0.0d0
                allocate(MatVf((ncpx-1)*(ncpy-1)*(ncpz-1),(ncpx-1)*(ncpy-1)*(ncpz-1)))
                MatVf = 0.0d0
                
                allocate(Kbar(ncpx*ncpy*ncpz*nds,ncpx*ncpy*ncpz*nds))
                Kbar = 0.0d0
                
                allocate(fbar(ncpx*ncpy*ncpz*nds,1))
                fbar = 0.0d0
                
            elseif(elcode=='Hex27EAS_PW')then
                allocate(MatAf(ncpx*ncpy*ncpz*nds,ncpx*ncpy*ncpz*nds))
                MatAf = 0.0d0
                allocate(MatCf(ncpx*ncpy*ncpz*nds,nalpha*nelems))
                MatCf = 0.0d0
                allocate(MatCfT(nalpha*nelems,ncpx*ncpy*ncpz*nds))
                MatCfT = 0.0d0
                allocate(MatVf(nalpha*nelems,nalpha*nelems))
                MatVf = 0.0d0
                
                allocate(diag(nalpha*nelems))
                diag = 0.0d0
                
                continue
                
            elseif(elcode=='Hex27ANS_PW')then
                allocate(MatAf(ncpx*ncpy*ncpz*nds,ncpx*ncpy*ncpz*nds))
                MatAf = 0.0d0
                allocate(MatCf(ncpx*ncpy*ncpz*nds,ncpx*ncpy*ncpz*nds))
                MatCf = 0.0d0
                allocate(MatCfT(ncpx*ncpy*ncpz*nds,ncpx*ncpy*ncpz*nds))
                MatCfT = 0.0d0
                allocate(MatVf(ncpx*ncpy*ncpz*nds,ncpx*ncpy*ncpz*nds))
                MatVf = 0.0d0
            endif
            
            continue
        
        elseif(line=='*bnet')then
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
            
            allocate(KfEq(nnodes*nds-nbc,nnodes*nds-nbc))
            KfEq = 0.0d0
            
            allocate(dispeqv(nnodes*nds-nbc,1))
            dispeqv = 0.0d0
            
            !allocate(KfEqInv(nnodes*nds-nbc,nnodes*nds-nbc))
            !KfEqInv = 0.0d0
            
            allocate(FextEq(nnodes*nds-nbc,1))
            FextEq = 0.0d0
            
            allocate(dDisp(nnodes*nds,1))
            dDisp=0.0d0
    
            allocate(dddisp(nnodes*nds,1))
            ddDisp=0.0d0
            
            allocate(u(nnodes*nds,1))
            u=0.0d0
            
            continue
        
        elseif(line == '*bcdof')then
            read(1,*)nbc
            
            allocate(bcdof(nbc))
            
            do k1=1,nbc
                read(1,*)int1,int2
                bcdof(k1) = (int1-1)*nds+int2
            end do
            
            count = 0
            do k=1, ncpz
                do j = 1,ncpy
                    do i = 1,ncpx
                        do k1=1,nds
                            count = count + 1
                            do k2=1,nbc
                                if(bcdof(k2)==count) bc(i,j,k,k1) = 0
                            end do
                        end do
                    end do
                enddo
            enddo
            
            allocate(KfEq(nnodes*nds-nbc,nnodes*nds-nbc))
            KfEq = 0.0d0
            
            allocate(KfEqInv(nnodes*nds-nbc,nnodes*nds-nbc))
            KfEqInv = 0.0d0
            
            allocate(FextEq(nnodes*nds-nbc,1))
            FextEq = 0.0d0
            
            allocate(dDisp(nnodes*nds,1))
            dDisp=0.0d0
    
            allocate(dddisp(nnodes*nds,1))
            ddDisp=0.0d0
            
            allocate(u(nnodes*nds,1))
            u=0.0d0
            
            allocate(dispeqv(nnodes*nds-nbc,1))
            dispeqv = 0.0d0
            
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
            read(1,*)iprops
            allocate(props(iprops))
            read(1,*)(props(i),i=1,iprops) 
            
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
            
        endif !line
    
    end do
    
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
    
    allocate(laxis(nelems,npi,nds,nds))
    allocate(laxisconv(nelems,npi,nds,nds))
    laxis = 0.0d0
    laxisconv = 0.0d0
    
    call MatLinearElastic3D(props(1),props(2),MatD)
    
    do k1=1,nelems
        do k2=1,npi
            TmatD(k1,k2,:,:) = matD(:,:)
        end do
    end do
    
    close(1,status='keep')
    
    continue
    
end
