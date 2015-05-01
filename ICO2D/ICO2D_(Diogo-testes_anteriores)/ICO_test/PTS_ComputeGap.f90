subroutine PTS_ComputeGap(iter,inc,istep,iLM)

    use Mod_Variables
    implicit none
    
    integer(4),intent(IN)::iter,inc,istep
    integer(4),intent(OUT)::iLM
    
    integer(4)::i,j,k,l,m,k1,k2,k3,k4
    integer(4)::temp1,count,count1,count2,kspan,ilixo
    integer(4),parameter::knewton = 50
    integer(4),dimension(ncp_m+1)::ctconn
    integer(4),dimension(ncp_s+ncp_m)::smconn
    integer(4),dimension(ncp_s)::isct
    real(8)::cm,c1,c2,c3,c4
    
    real(8)::xib,lm,dlm,sumres,rsd,gn,a11,b11,gt,num,tst
    real(8)::xibi,xibf,diff
    real(8),dimension(p_s+1)::Rss,dRss,ddRss
    real(8),dimension(p_m+1)::Rmm,dRmm,ddRmm
    real(8),dimension(ncp_s)::Rs,dRsdxi,d2Rsdxi2
    real(8),dimension(ncp_m)::Rm,dRmdxi,d2Rmdxi2,sumRm
    real(8),dimension(nds,1)::tang,dtang,nor,dpos,aux21
    real(8),dimension(nds,1)::xs,xm,dxm,vec,dvec,vec0
    real(8),dimension(nds,1)::gap,x11,x21,gapi,gapf
    real(8),dimension((ncp_m+ncp_s)*nds,1)::Ns,N0s,Ts,T0s
    real(8),dimension((ncp_m+ncp_s)*nds,(ncp_m+ncp_s)*nds)::Kdelta
    real(8),dimension(:,:),allocatable::FLM,FLMStr,KLM,LMStr,GLM,GLMStr,Ki
    real(8),dimension(:,:),allocatable::KfeqLM,FeqLM,dispLM
    real(8),dimension(:),allocatable::disptempLM
    real(8),dimension(:,:),allocatable::ImpDispLM
    
    real(8),dimension(nds+1,nds+1)::Kc
    real(8),dimension(nds+1,1)::Fc
    
    real(8),dimension(ncp_s*nds+1,ncp_s*nds+1)::Kcc
    real(8),dimension(ncp_s*nds+1,1)::Fcc
    real(8),dimension(ncp_s)::str_xib,Fs
    real(8),dimension(ncp_m)::Fm
    
    real(8)::gti,gtf,xib0,rsd0,sumt,sumd,sumdd,cfv
    
    real(8),dimension(:),allocatable::temp
    real(8),dimension(:,:),allocatable::n_str
    
    !Contact Forces
    real(8)::suml,sumr,dx,dy,sumdl
    real(8),dimension(ncp_s)::CF,bl,br,dl,Rl,Rr,norm,normx,normy
    real(8),dimension(ncp_s,2)::Fslv
    real(8),dimension(ncp_m,2)::Fmst
    real(8),dimension(ncp_m,ncp_m)::MM
    logical::tag,gaplog
    
    real(8),dimension(:),allocatable::tempinv
    
    allocate(KLM(tnodes*nds+ngrv,tnodes*nds+ngrv))
    KLM = 0.0d0
    allocate(GLM(tnodes*nds+ngrv,1))
    GLM = 0.0d0
    allocate(FLM(tnodes*nds,ngrv))
    FLM = 0.0d0
    
    allocate(n_str((ncp_m+ncp_s)*nds,ngrv))
    n_str = 0.0d0
    
    iLM = 0
    
    sumRm = 0.0d0
    
    if(istep==1 .and. iter==1 .and. inc==1)then
        continue
        dLagrange = 0.0d0
    elseif(iter == 1)then
        !Lagrange = 0.0d0 !LagrangeConv
        dLagrange = 0.0d0
    end if
    
    do i=1,ngrv
        
        !------------------------------------------------------------------
        ! LAGRANGE MULTIPLIER METHOD - Compute gap
        !------------------------------------------------------------------
        
        !Compute B-Splines basis functions...
        call BSplineBasisAndDeriv2(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan)
        
        !... and convert to NURBS basis, including vanishing terms 
        call BSplinesToNURBS(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
        
        !Slave (Greville) point in physical slave surface
        xs = 0.0d0
        do j=1,ncp_s
            xs(1,1) = xs(1,1) + Rs(j)*(GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*nds-1),1))
            xs(2,1) = xs(2,1) + Rs(j)*(GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*nds  ),1))
        end do
        
        plot_xy(i,1) = xs(1,1)
        plot_xy(i,2) = xs(2,1)
        
        continue
        
        !------------------------------------------------------------------
        ! Contact check (maybe)
        !------------------------------------------------------------------
!        xibi = 0.0d0
!        !Compute B-Splines basis functions...
!        call BSplineBasisAndDeriv2(ncp_m,p_m,xibi,u_knot_master,Rmm,dRmm,ddRmm,kspan)
!        
!        !... and convert to NURBS basis, including vanishing terms 
!        call BSplinesToNURBS(ncp_m,p_m,xibi,u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
!        
!        !Tangent vector in the master segment and its derivative
!        tang = 0.0d0
!        dtang = 0.0d0
!        do k=1,ncp_m
!            tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*nds-1,1))
!            tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*nds  ,1))
!        end do
!        
!        !Length of the master segment and its derivative
!        lm  = dsqrt( tang(1,1)*tang(1,1) + tang(2,1)*tang(2,1))
!        
!        !Normalise tangent vectors
!        tang = tang/lm
!        
!        !Possible contact check
!        x11(1,1) = GCoords(conn_master(1),1) + ddisp(conn_master(1)*nds-1,1)
!        x11(2,1) = GCoords(conn_master(1),2) + ddisp(conn_master(1)*nds  ,1)
!        gapi = xs - x11
!        gti = gapi(1,1)*tang(1,1) + gapi(2,1)*tang(2,1)
!        
!        xibf = 1.0d0
!        !Compute B-Splines basis functions...
!        call BSplineBasisAndDeriv2(ncp_m,p_m,xibf,u_knot_master,Rmm,dRmm,ddRmm,kspan)
!        
!        !... and convert to NURBS basis, including vanishing terms 
!        call BSplinesToNURBS(ncp_m,p_m,xibf,u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
!        
!        !Tangent vector in the master segment
!        tang = 0.0d0
!        dtang = 0.0d0
!        do k=1,ncp_m
!            tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
!            tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
!        end do
!        
!        !Length of the master segment
!        lm  = dsqrt( tang(1,1)*tang(1,1) + tang(2,1)*tang(2,1))
!        
!        !Normalise tangent vector
!        tang = tang/lm
        
        !------------------------------------------------------------------
        ! Secant Method to determine the Closest Point Projection
        !------------------------------------------------------------------
        xib0 = 0.0d0
        xib = 0.123581d0 ! De onde surge este valor??????? 
        do j=1,knewton
            if(j==1)then
                !Compute B-Splines basis functions...
                call BSplineBasisAndDeriv2(ncp_m,p_m,xib0,u_knot_master,Rmm,dRmm,ddRmm,kspan)
                
               !... and convert to NURBS basis, including vanishing terms 
                call BSplinesToNURBS(ncp_m,p_m,xib0,u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
                
                !Closest Point Projection and its derivative
                xm = 0.0d0
                do k=1,ncp_m
                    xm(1,1)  = xm(1,1)  + Rm(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*nds-1,1))
                    xm(2,1)  = xm(2,1)  + Rm(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*nds  ,1))
                    continue
                end do
                
                !Tangent vector in the master segment and its derivative
                tang = 0.0d0
                do k=1,ncp_m
                    tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*nds-1,1))
                    tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*nds  ,1))
                end do
                
                !Length of the master segment and its derivative
                lm  = dsqrt(tang(1,1)*tang(1,1) + tang(2,1)*tang(2,1))
                
                !Normalise tangent vectors
                tang = tang/lm
                
                continue
                
                !Auxiliary vectors
                vec = xs - xm
                rsd0 = vec(1,1)*tang(1,1)+vec(2,1)*tang(2,1)
                
            end if
            
            !Compute B-Splines basis functions...
            call BSplineBasisAndDeriv2(ncp_m,p_m,xib,u_knot_master,Rmm,dRmm,ddRmm,kspan)
                
            !... and convert to NURBS basis, including vanishing terms 
            call BSplinesToNURBS(ncp_m,p_m,xib,u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
            
            !Closest Point Projection and its derivative
            xm = 0.0d0
            dxm = 0.0d0
            do k=1,ncp_m
                xm(1,1)  = xm(1,1)  + Rm(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*nds-1,1))
                xm(2,1)  = xm(2,1)  + Rm(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*nds  ,1))
                continue
            end do
            
            !Tangent vector in the master segment and its derivative
            tang = 0.0d0
            do k=1,ncp_m
                tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*nds-1,1))
                tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*nds  ,1))
            end do
            
            !Length of the master segment and its derivative
            lm  = dsqrt(tang(1,1)*tang(1,1) + tang(2,1)*tang(2,1))
            
            !Normalise tangent vectors
            tang = tang/lm
            
            continue
            
            !Auxiliary vectors
            vec = xs - xm
            
            rsd = vec(1,1)*tang(1,1) + vec(2,1)*tang(2,1)
            
            tst = rsd*(xib - xib0)/(rsd - rsd0)
            
            xibi = xib
            
            xib = xib - rsd*(xib - xib0)/(rsd - rsd0)
            
            ! In the limits
            if(xib .gt. 1.0d0) xib = 1.0d0
            if(xib .lt. 0.0d0) xib = 0.0d0
            
            xib0 = xibi
            rsd0 = rsd
            
            if(abs(rsd) .lt. 1.0d-12) goto 1
            
            continue
            
        end do
        
        write(*,FMT=16)i
        16 format('Warning: Failed to determine CPP for Greville point ',I3)
        
        1 continue
        
        !Compute normal vector from cross product {tang}x{0 0 -1}
        nor(1,1) = -tang(2,1)
        nor(2,1) =  tang(1,1)
        
        !Compute initial normal gap
        initial_gap(i) = final_gap(i) 
        final_gap(i) = vec(1,1)*nor(1,1) + vec(2,1)*nor(2,1)
        
        !---------------------------------------------
        ! Commented by: Diogo
        ! Problem compiling ISNAN with IMPLICIT NONE 
        !---------------------------------------------
        !gaplog = ISNAN (final_gap(i))
        !if(gaplog == .true.) final_gap(i) = 100.0d0
        
        if(final_gap(i) <= -1.0d-8) then
            PTS_active(i) = 1
        else
            PTS_active(i) = 0
        end if
        
        !----------------------------------------------------------------------------------------------------------------------------------
        ! Tese da Matzen: serve para quando temos condicao de "descolar" um contacto, no caso de compressao e libertacao
        ! Ter atencao que o MP_props(1,1) nao esta bem para todos os casos. Deve ser a propriedade 1 do master, que nem sempre e o patch 1
        !----------------------------------------------------------------------------------------------------------------------------------
        cfv = -1.0d0*Lagrange(i) - MP_props(1,1)*final_gap(i)
        
        if(cfv .gt. 0.0d0)then
            PTS_active(i) = 1
        else
            PTS_active(i) = 0
        end if
        
!        if(Lagrange(i)==0.0d0)then
!            if(final_gap(i) <= -1.0d-10) then
!                PTS_active(i) = 1
!            else
!                PTS_active(i) = 0
!            end if
!        end if
        
        !If contact has occured
        Kdelta = 0.0d0
        if(PTS_active(i) == 1)then
            
            gn = final_gap(i)
            
            !Metric of the boundary
            a11 = lm*lm
            
            !Second derivative of the Closest Point Projection
            dpos = 0.0d0
            do k=1,ncp_m
                dpos(1,1) = dpos(1,1) + d2Rmdxi2(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*nds-1,1))
                dpos(2,1) = dpos(2,1) + d2Rmdxi2(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*nds  ,1))
                continue
            end do
            
            !Curvature of the boundary
            b11 = dpos(1,1)*nor(1,1) + dpos(2,1)*nor(2,1)
        
            !Auxiliary Matrices
            Ns = 0.0d0
            Ts = 0.0d0
            N0s = 0.0d0
            !T0s = 0.0d0
            
            ! Slave components of the vectors 
            do j=1,ncp_s
                Ns(j*nds-1,1) = Rs(j)*nor(1,1)
                Ns(j*nds  ,1) = Rs(j)*nor(2,1)
                
                Ts(j*nds-1,1) = Rs(j)*tang(1,1)
                Ts(j*nds  ,1) = Rs(j)*tang(2,1)
            end do
            
            ! Master components of the vectors 
            count = 0
            do j=ncp_s+1,ncp_s+ncp_m
                
                count = count + 1
                
                Ns(j*nds-1,1) = -1.0d0*Rm(count)*nor(1,1)
                Ns(j*nds  ,1) = -1.0d0*Rm(count)*nor(2,1)
                
                N0s(j*nds-1,1) = dRmdxi(count)*nor(1,1)
                N0s(j*nds  ,1) = dRmdxi(count)*nor(2,1)
                
                Ts(j*nds-1,1) = -1.0d0*Rm(count)*tang(1,1)
                Ts(j*nds  ,1) = -1.0d0*Rm(count)*tang(2,1)
                
                ! Esta variavel nao esta a ser utilizada
                n_str(j*nds-1,i) = -1.0d0*Rm(count)
                n_str(j*nds  ,i) = -1.0d0*Rm(count)
                
                continue
            end do
            
            !Coefficients for the contact stiffness matrices
            cm = a11 - gn*b11
            c1 = -1.0d0*lm/cm - b11*lm*gn/(cm*cm) + b11*gn/(cm*lm) + b11*b11*gn*gn/(lm*cm*cm)
            c2 = -1.0d0*lm/cm - b11*lm*gn/(cm*cm) + b11*gn/(cm*lm) + b11*b11*gn*gn/(lm*cm*cm)
            c3 = -2.0d0*gn/cm - b11*gn*gn/(cm*cm) + gn/(lm*lm) + 2.0d0*b11*gn*gn/(lm*lm*cm) + b11*b11*gn*gn*gn/(lm*lm*cm*cm)
            c4 = -1.0d0*b11*lm*lm/(cm*cm) + b11*b11*gn/(cm*cm)
            
            !Linearized variation of the gap
            ! O Kdelta deve estar sempre activo - esta condicao de GNL nao devia existir: Diogo
            if(nlgeom==.false.)then
                Kdelta = 0.0d0
            else
                Kdelta = 0.0d0
                Kdelta = c1*matmul(N0s,transpose( Ts)) + &
                &        c2*matmul( Ts,transpose(N0s)) + &
                &        c3*matmul(N0s,transpose(N0s)) + &
                &        c4*matmul( Ts,transpose( Ts))
                
                Kdelta = Kdelta*Lagrange(i) !*-1.0d0
            end if
            
            iLM = iLM + 1
            
            continue 
            
            !------------------------------------------------------------------
            ! LAGRANGE MULTIPLIER METHOD - Normal Contact
            !------------------------------------------------------------------
            
            ! Esta variavel nao esta a ser utilizada: Diogo 
!            ctconn = 0
!            ctconn(1) = conn_slave(i)
!            do j=1,ncp_m
!                ctconn(j+1)=conn_master(j)
!            end do
            
            do j=1,ncp_s
                smconn(j)=conn_slave(j)
            end do
            
            count = 0
            do j=ncp_s+1,ncp_m+ncp_s
                count = count + 1
                smconn(j)=conn_master(count)
                continue
            end do
            
            !Deformable-rigid contact
!            do k1=1,ncp_s
!                do k2=1,ncp_s
!                    KLM(conn_slave(k1)*2-1,conn_slave(k2)*2-1) = KLM(conn_slave(k1)*2-1,conn_slave(k2)*2-1) + Kdelta(k1*2-1,k2*2-1)
!                    KLM(conn_slave(k1)*2-1,conn_slave(k2)*2  ) = KLM(conn_slave(k1)*2-1,conn_slave(k2)*2  ) + Kdelta(k1*2-1,k2*2  )
!                    KLM(conn_slave(k1)*2  ,conn_slave(k2)*2-1) = KLM(conn_slave(k1)*2  ,conn_slave(k2)*2-1) + Kdelta(k1*2  ,k2*2-1)
!                    KLM(conn_slave(k1)*2  ,conn_slave(k2)*2  ) = KLM(conn_slave(k1)*2  ,conn_slave(k2)*2  ) + Kdelta(k1*2  ,k2*2  )
!                end do
!                
!                KLM(conn_slave(k1)*2-1,tnodes*nds+i) =  KLM(conn_slave(k1)*2-1,tnodes*nds+i) + Ns(k1*2-1,1)
!                KLM(conn_slave(k1)*2  ,tnodes*nds+i) =  KLM(conn_slave(k1)*2  ,tnodes*nds+i) + Ns(k1*2  ,1)
!                
!                KLM(tnodes*nds+i,conn_slave(k1)*2-1) =  KLM(tnodes*nds+i,conn_slave(k1)*2-1) + Ns(k1*2-1,1)
!                KLM(tnodes*nds+i,conn_slave(k1)*2  ) =  KLM(tnodes*nds+i,conn_slave(k1)*2  ) + Ns(k1*2  ,1)
!                
!                FLM(conn_slave(k1)*nds-1,i) = FLM(conn_slave(k1)*nds-1,i) + Ns(k1*nds-1,1)
!                FLM(conn_slave(k1)*nds  ,i) = FLM(conn_slave(k1)*nds  ,i) + Ns(k1*nds  ,1)
!                
!            end do
            
            !Deformable-Deformable contact
            do k1=1,ncp_s+ncp_m
                do k2=1,ncp_s+ncp_m
                    KLM(smconn(k1)*nds-1,smconn(k2)*nds-1) = KLM(smconn(k1)*nds-1,smconn(k2)*nds-1) + Kdelta(k1*nds-1,k2*nds-1)
                    KLM(smconn(k1)*nds-1,smconn(k2)*nds  ) = KLM(smconn(k1)*nds-1,smconn(k2)*nds  ) + Kdelta(k1*nds-1,k2*nds  )
                    KLM(smconn(k1)*nds  ,smconn(k2)*nds-1) = KLM(smconn(k1)*nds  ,smconn(k2)*nds-1) + Kdelta(k1*nds  ,k2*nds-1)
                    KLM(smconn(k1)*nds  ,smconn(k2)*nds  ) = KLM(smconn(k1)*nds  ,smconn(k2)*nds  ) + Kdelta(k1*nds  ,k2*nds  )
                end do
                
                KLM(smconn(k1)*nds-1,tnodes*nds+i) =  KLM(smconn(k1)*nds-1,tnodes*nds+i) + Ns(k1*nds-1,1)
                KLM(smconn(k1)*nds  ,tnodes*nds+i) =  KLM(smconn(k1)*nds  ,tnodes*nds+i) + Ns(k1*nds  ,1)
                
                KLM(tnodes*nds+i,smconn(k1)*nds-1) =  KLM(tnodes*nds+i,smconn(k1)*nds-1) + Ns(k1*nds-1,1)
                KLM(tnodes*nds+i,smconn(k1)*nds  ) =  KLM(tnodes*nds+i,smconn(k1)*nds  ) + Ns(k1*nds  ,1)
                
                FLM(smconn(k1)*nds-1,i) = FLM(smconn(k1)*nds-1,i) + Ns(k1*nds-1,1)
                FLM(smconn(k1)*nds  ,i) = FLM(smconn(k1)*nds  ,i) + Ns(k1*nds  ,1)
            end do

            
            !Contact Residual
            GLM(tnodes*nds+i,1) = -1.0d0*final_gap(i)
            
            continue

        else
            
            99 continue
            
            iLM = iLM + 1
            KLM(tnodes*nds+i,:) = 0.0d0
            KLM(:,tnodes*nds+i) = 0.0d0
            KLM(tnodes*nds+i,tnodes*nds+i) = 1.0d0
            !FLM(:,i) = 0.0d0
            !Lagrange(i) = 0.0d0
            
            do j=1,ncp_s
                smconn(j)=conn_slave(j)
            end do
            count = 0
            do j=ncp_s+1,ncp_m+ncp_s
                count = count + 1
                smconn(j)=conn_master(count)
                continue
            end do
            
            do k1=1,ncp_s+ncp_m
                FintLM(smconn(k1)*nds-1,1) = 0.0d0
                FintLM(smconn(k1)*nds  ,1) = 0.0d0
            end do
              
        end if !gn .lt. 0.0d0
        
    end do !Greville points cycle
    
    ! Solve system of equations in the contact subroutine
    if(iLM .gt. 0)then
             
        !if(iter==1 .and. inc == 1)then
        allocate(KfeqLM(tnodes*nds-nbc+ngrv,tnodes*nds-nbc+ngrv))
        allocate(FeqLM(tnodes*nds-nbc+ngrv,1),dispLM(tnodes*nds-nbc+ngrv,1))
        allocate(disptempLM(tnodes*nds+ngrv),impdispLM(tnodes*nds+ngrv,1))
        !end do
        
        !Store Contact Matrix ----
        KCont = 0.0d0
        KCont(1:tnodes*nds,1:tnodes*nds) =  KCont(1:tnodes*nds,1:tnodes*nds) + KLM(1:tnodes*nds,1:tnodes*nds)
        
        !Total Stiffness matrix -----
        KLM(1:tnodes*nds,1:tnodes*nds) = KLM(1:tnodes*nds,1:tnodes*nds) + Kf(1:tnodes*nds,1:tnodes*nds)
        KT(1:tnodes*nds,1:tnodes*nds) = KT(1:tnodes*nds,1:tnodes*nds) + KLM(1:tnodes*nds,1:tnodes*nds)
        
        do j=1,ncp_s
            smconn(j)=conn_slave(j)
        end do
        
        count = 0
        do j=ncp_s+1,ncp_m+ncp_s
            count = count + 1
            smconn(j)=conn_master(count)
            continue
        end do
        
        GLM(1:tnodes*nds,1) = GLM(1:tnodes*nds,1) + Fext(1:tnodes*nds,1)
        
        disptempLM = 1.0d0
        disptempLM(1:tnodes*nds) = redcount(1:tnodes*nds)
        
        ImpDispLM = 0.0d0
        ImpDispLM(1:tnodes*nds,1) = ImpDisp(1:tnodes*nds,1)
        do i=1,ngrv
            if(PTS_active(i) == 0)then
                ImpDispLM(tnodes*nds+i,1) = -1.0d0*Lagrange(i)
            end if
        end do
        
        !Prescribed Lagrange Multipliers -----
        do i=tnodes*nds+1,tnodes*nds+iLM
            if(ImpDispLM(i,1)/=0.0d0)then
                do k1=tnodes*nds,tnodes*nds+iLM
                    if(i /= k1) then
                        GLM(k1,1) = GLM(k1,1) - KLM(k1,i)*ImpDispLM(i,1)
                        KLM(i,k1) = 0.0d0
                        KLM(k1,i) = 0.0d0
                     else 
                        GLM(i,1) = KLM(i,i)*ImpDispLM(i,1)
                        continue
                    end if
                end do
            end if
        end do
        
        k=1
        m=1
        KfeqLM = 0.0d0
        FEqLM = 0.0d0
        do i=1,tnodes*nds + ngrv
            do j=1,tnodes*nds + ngrv
                if (disptempLM(i)/=0.0d0 .and. disptempLM(j)/=0.0d0)then
                    KfeqLM(k,m)=KLM(i,j)
                    FEqLM(k,1)=GLM(i,1)                
                    m=m+1
                    if(m==tnodes*nds-nBC+ngrv+1)then
                        k=k+1
                        m=1
                    end if
                end if
            end do
        end do
        
        continue
        
        temp1=tnodes*nds-nBC+ngrv
        dispLM = 0.0d0
        call Gauss (temp1,KfeqLM,FEqLM,dispLM)
        
        ddDisp = 0.0d0
        j=1
        do i=1,tnodes*nds
            if(redcount(i)==1.0d0)then
                ddDisp(i,1)=dispLM(j,1)
                j=j+1
            end if
        end do
        
        do i=1,ngrv
            
            if(dispLM(tnodes*nds-nbc+i,1) .gt. 0.0d0)then
                if(dispLM(tnodes*nds-nbc+i,1) .gt. 0.0d0 .and. (abs(dispLM(tnodes*nds-nbc+i,1)) .gt. abs(Lagrange(i)))) then
                    dispLM(tnodes*nds-nbc+i,1) = -1.0d0*Lagrange(i)
                    PTS_active(i) = 0
                end if
            end if
            
            Lagrange(i) = Lagrange(i) + dispLM(tnodes*nds-nbc+i,1)
            
            if(Lagrange(i) .gt. 0.0d0)then
                Lagrange(i) = 0.0d0
            end if
            
            if(Lagrange(i) == 0.0d0)then
                FLM(:,i) = 0.0d0
            end if
            
            dLagrange(i) = dispLM(tnodes*nds-nbc+i,1)
            
            LagrangeConv(i) = LagrangeConv(i) + dispLM(tnodes*nds-nbc+i,1)
            
        end do
         
        if(sum(Lagrange) == 0.0d0)then
            FintLM = 0.0d0
        end if
        
        sumLag = sum(Lagrange)
        
        sumRm = sumRm/sum(sumRm)
        
        !Eliminate Contact Forces in master control points outside contact area
!        do k1=1,ncp_m
!            if(sumRm(k1) == 0.0d0)then
!                FintLM(conn_master(k1)*nds-1,1) = 0.0d0
!                FintLM(conn_master(k1)*nds  ,1) = 0.0d0
!            end if
!        end do
        
        !Right-hand side -----------------------------------------------------------
        FintLM = 0.0d0
        count = 0
        do i=1,ngrv
!            do k1=1,ncp_s+ncp_m
!                FintLM(smconn(k1)*nds-1,1) = FintLM(smconn(k1)*nds-1,1) + FLM(smconn(k1)*nds-1,i)*Lagrange(i)
!                FintLM(smconn(k1)*nds  ,1) = FintLM(smconn(k1)*nds  ,1) + FLM(smconn(k1)*nds  ,i)*Lagrange(i)
!            end do

            FintLM(:,1) = FintLM(:,1) + FLM(:,i)*Lagrange(i)

        end do

    end if
    
    
    continue
    
end subroutine


















!subroutine PTS_ComputeGap(iter,inc,istep,iLM)
!
!    use Mod_Variables
!    implicit none
!    
!    integer(4),intent(IN)::iter,inc,istep
!    integer(4),intent(OUT)::iLM
!    
!    integer(4)::i,j,k,l,m,k1,k2,k3,k4
!    integer(4)::temp1,count,count1,count2,kspan,ilixo
!    integer(4),parameter::knewton = 50
!    integer(4),dimension(ncp_m+1)::ctconn
!    integer(4),dimension(ncp_s)::isct
!    real(8)::cm,c1,c2,c3,c4
!    
!    real(8)::xib,lm,dlm,sumres,rsd,gn,a11,b11,gt,num,tst
!    real(8)::xibi,xibf,diff
!    real(8),dimension(p_s+1)::Rss,dRss,ddRss
!    real(8),dimension(p_m+1)::Rmm,dRmm,ddRmm
!    real(8),dimension(ncp_s)::Rs,dRsdxi,d2Rsdxi2
!    real(8),dimension(ncp_m)::Rm,dRmdxi,d2Rmdxi2,sumRm
!    real(8),dimension(nds,1)::tang,dtang,nor,dpos,aux21
!    real(8),dimension(nds,1)::xs,xm,dxm,vec,dvec,vec0
!    real(8),dimension(nds,1)::gap,x11,x21,gapi,gapf
!    real(8),dimension((ncp_m+1)*nds,1)::Ns,N0s,Ts,T0s
!    real(8),dimension((ncp_m+1)*nds,(ncp_m+1)*nds)::Kdelta
!    real(8),dimension(:,:),allocatable::FLM,FLMStr,KLM,LMStr,GLM,GLMStr,Ki
!    real(8),dimension(:,:),allocatable::KfeqLM,FeqLM,dispLM
!    real(8),dimension(:),allocatable::disptempLM
!    real(8),dimension(:,:),allocatable::ImpDispLM
!    
!    real(8),dimension(nds+1,nds+1)::Kc
!    real(8),dimension(nds+1,1)::Fc
!    
!    real(8),dimension(ncp_s*nds+1,ncp_s*nds+1)::Kcc
!    real(8),dimension(ncp_s*nds+1,1)::Fcc
!    real(8),dimension(ncp_s)::str_xib,Fs
!    real(8),dimension(ncp_m)::Fm
!    
!    real(8)::gti,gtf,xib0,rsd0,sumt,sumd,sumdd,cfv
!    
!    real(8),dimension(:),allocatable::temp
!    real(8),dimension(:,:),allocatable::n_str
!    
!    !Contact Forces
!    real(8)::suml,sumr,dx,dy,sumdl
!    real(8),dimension(ncp_s)::CF,bl,br,dl,Rl,Rr,norm,normx,normy
!    real(8),dimension(ncp_s,2)::Fslv
!    real(8),dimension(ncp_m,2)::Fmst
!    real(8),dimension(ncp_m,ncp_m)::MM
!    logical::tag,gaplog
!    
!    real(8),dimension(:),allocatable::tempinv
!    
!    allocate(KLM(tnodes*nds+ngrv,tnodes*nds+ngrv))
!    KLM = 0.0d0
!    allocate(GLM(tnodes*nds+ngrv,1))
!    GLM = 0.0d0
!    allocate(FLM(tnodes*nds,ngrv))
!    FLM = 0.0d0
!    
!    allocate(n_str((ncp_m)*nds,ngrv))
!    n_str = 0.0d0
!    
!    iLM = 0
!    
!    sumRm = 0.0d0
!    
!    if(istep==1 .and. iter==1 .and. inc==1)then
!        continue
!        dLagrange = 0.0d0
!    elseif(iter == 1)then
!        !Lagrange = 0.0d0 !LagrangeConv
!        dLagrange = 0.0d0
!    end if
!    
!    do i=1,ngrv
!        
!        !------------------------------------------------------------------
!        ! LAGRANGE MULTIPLIER METHOD - Compute gap
!        !------------------------------------------------------------------
!        
!        !Compute B-Splines basis functions...
!        call BSplineBasisAndDeriv2(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan)
!        
!        !... and convert to NURBS basis, including vanishing terms 
!        call BSplinesToNURBS(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
!        
!        !Slave (Greville) point in physical slave surface
!        xs = 0.0d0
!        do j=1,ncp_s
!            xs(1,1) = xs(1,1) + Rs(j)*(GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*2-1),1))
!            xs(2,1) = xs(2,1) + Rs(j)*(GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*2  ),1))
!        end do
!        
!        plot_xy(i,1) = xs(1,1)
!        plot_xy(i,2) = xs(2,1)
!        
!        continue
!        
!        !------------------------------------------------------------------
!        ! Contact check (maybe)
!        !------------------------------------------------------------------
!        xibi = 0.0d0
!        !Compute B-Splines basis functions...
!        call BSplineBasisAndDeriv2(ncp_m,p_m,xibi,u_knot_master,Rmm,dRmm,ddRmm,kspan)
!        
!        !... and convert to NURBS basis, including vanishing terms 
!        call BSplinesToNURBS(ncp_m,p_m,xibi,u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
!        
!        !Tangent vector in the master segment and its derivative
!        tang = 0.0d0
!        dtang = 0.0d0
!        do k=1,ncp_m
!            tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
!            tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
!        end do
!        
!        !Length of the master segment and its derivative
!        lm  = dsqrt( tang(1,1)*tang(1,1) + tang(2,1)*tang(2,1))
!        
!        !Normalise tangent vectors
!        tang = tang/lm
!        
!        !Possible contact check
!        x11(1,1) = GCoords(conn_master(1),1) + ddisp(conn_master(1)*2-1,1)
!        x11(2,1) = GCoords(conn_master(1),2) + ddisp(conn_master(1)*2  ,1)
!        gapi = xs - x11
!        gti = gapi(1,1)*tang(1,1) + gapi(2,1)*tang(2,1)
!        
!        xibf = 1.0d0
!        !Compute B-Splines basis functions...
!        call BSplineBasisAndDeriv2(ncp_m,p_m,xibf,u_knot_master,Rmm,dRmm,ddRmm,kspan)
!        
!        !... and convert to NURBS basis, including vanishing terms 
!        call BSplinesToNURBS(ncp_m,p_m,xibf,u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
!        
!        !Tangent vector in the master segment
!        tang = 0.0d0
!        dtang = 0.0d0
!        do k=1,ncp_m
!            tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
!            tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
!        end do
!        
!        !Length of the master segment
!        lm  = dsqrt( tang(1,1)*tang(1,1) + tang(2,1)*tang(2,1))
!        
!        !Normalise tangent vector
!        tang = tang/lm
!        
!        !------------------------------------------------------------------
!        ! Secant Method to determine the Closest Point Projection
!        !------------------------------------------------------------------
!        xib0 = 0.0d0
!        xib = 0.5d0
!        do j=1,knewton
!            if(j==1)then
!                !Compute B-Splines basis functions...
!                call BSplineBasisAndDeriv2(ncp_m,p_m,xib0,u_knot_master,Rmm,dRmm,ddRmm,kspan)
!                
!               !... and convert to NURBS basis, including vanishing terms 
!                call BSplinesToNURBS(ncp_m,p_m,xib0,u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
!                
!                !Closest Point Projection and its derivative
!                xm = 0.0d0
!                dxm = 0.0d0
!                do k=1,ncp_m
!                    xm(1,1)  = xm(1,1)  + Rm(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
!                    xm(2,1)  = xm(2,1)  + Rm(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
!                    continue
!                end do
!                
!                !Tangent vector in the master segment and its derivative
!                tang = 0.0d0
!                do k=1,ncp_m
!                    tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
!                    tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
!                end do
!                
!                !Length of the master segment and its derivative
!                lm  = dsqrt( tang(1,1)*tang(1,1) + tang(2,1)*tang(2,1))
!                
!                !Normalise tangent vectors
!                tang = tang/lm
!                
!                continue
!                
!                !Auxiliary vectors
!                vec = xs - xm
!                rsd0 = vec(1,1)*tang(1,1)+vec(2,1)*tang(2,1)
!                
!            end if
!            
!            !Compute B-Splines basis functions...
!            call BSplineBasisAndDeriv2(ncp_m,p_m,xib,u_knot_master,Rmm,dRmm,ddRmm,kspan)
!                
!            !... and convert to NURBS basis, including vanishing terms 
!            call BSplinesToNURBS(ncp_m,p_m,xib,u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
!            
!            !Closest Point Projection and its derivative
!            xm = 0.0d0
!            dxm = 0.0d0
!            do k=1,ncp_m
!                xm(1,1)  = xm(1,1)  + Rm(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
!                xm(2,1)  = xm(2,1)  + Rm(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
!                continue
!            end do
!            
!            !Tangent vector in the master segment and its derivative
!            tang = 0.0d0
!            do k=1,ncp_m
!                tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
!                tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
!            end do
!            
!            !Length of the master segment and its derivative
!            lm  = dsqrt( tang(1,1)*tang(1,1) + tang(2,1)*tang(2,1))
!            
!            !Normalise tangent vectors
!            tang = tang/lm
!            
!            continue
!            
!            !Auxiliary vectors
!            vec = xs - xm
!            
!            rsd = vec(1,1)*tang(1,1)+vec(2,1)*tang(2,1)
!            
!            tst = rsd*(xib - xib0)/(rsd - rsd0)
!            
!            xibi = xib
!            
!            xib = xib - rsd*(xib - xib0)/(rsd - rsd0)
!            
!            if(xib .gt. 1.0d0) xib = 1.0d0
!            if(xib .lt. 0.0d0) xib = 0.0d0
!            
!            xib0 = xibi
!            rsd0 = rsd
!            
!            if(abs(rsd) .lt. 1.0d-10) goto 1
!            
!            continue
!            
!        end do
!        
!        write(*,FMT=16)i
!        16 format('Warning: Failed to determine initial CPP for Greville point ',I3)
!        
!        1 continue
!        
!        str_xib(i) = xib
!        sumRm = sumRm + Rm
!        
!        !Compute normal vector from cross product {tang}x{0 0 -1}
!        nor(1,1) = -tang(2,1)
!        nor(2,1) =  tang(1,1)
!        
!        !Compute initial normal gap
!        initial_gap(i) = final_gap(i) 
!        final_gap(i) = vec(1,1)*nor(1,1)+vec(2,1)*nor(2,1)
!        
!        gaplog = ISNAN (final_gap(i))
!        if(gaplog == .true.) final_gap(i) = 100.0d0
!        
!        if(final_gap(i) <= -1.0d-8) then
!            PTS_active(i) = 1
!        else
!            PTS_active(i) = 0
!        end if
!        
!        cfv = -1.0d0*Lagrange(i) - MP_props(1,1)*final_gap(i)
!        
!        !if(cfv .gt. 0.0d0 .or. final_gap(i) <= 0.0d0)then
!        if(cfv .gt. 0.0d0)then
!            PTS_active(i) = 1
!        else
!            PTS_active(i) = 0
!        end if
!        
!        if(Lagrange(i)==0.0d0)then
!            if(final_gap(i) <= -1.0d-10) then
!                PTS_active(i) = 1
!            else
!                PTS_active(i) = 0
!            end if
!        end if
!        
!        !If contact has occured
!        Kdelta = 0.0d0
!        if(PTS_active(i) == 1)then
!            
!            !Metric of the boundary
!            a11 = lm*lm
!            
!            !Second derivative of the Closest Point Projection
!            dpos = 0.0d0
!            do k=1,ncp_m
!                dpos(1,1) = dpos(1,1) + d2Rmdxi2(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
!                dpos(2,1) = dpos(2,1) + d2Rmdxi2(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
!                continue
!            end do
!            
!            !Curvature of the boundary
!            b11 = dpos(1,1)*nor(1,1) + dpos(2,1)*nor(2,1)
!        
!            !Auxiliary Matrices
!            Ns = 0.0d0
!            Ts = 0.0d0
!            N0s = 0.0d0
!            !T0s = 0.0d0
!            
!            do j=1,nds
!                Ns(j,1) = nor(j,1)
!                Ts(j,1) = tang(j,1)
!            end do
!            
!            do j=1,ncp_m
!                Ns(nds+j*2-1,1) = -1.0d0*Rm(j)*nor(1,1)
!                Ns(nds+j*2  ,1) = -1.0d0*Rm(j)*nor(2,1)
!                
!                N0s(nds+j*2-1,1) = dRmdxi(j)*nor(1,1)
!                N0s(nds+j*2  ,1) = dRmdxi(j)*nor(2,1)
!                
!                Ts(nds+j*2-1,1) = -1.0d0*Rm(j)*tang(1,1)
!                Ts(nds+j*2  ,1) = -1.0d0*Rm(j)*tang(2,1)
!                
!                n_str(j*2-1,i) = -1.0d0*Rm(j)
!                n_str(j*2  ,i) = -1.0d0*Rm(j)
!                
!                continue
!            end do
!            
!            !Coefficients for the contact stiffness matrices
!            cm = a11 - gn*b11
!            c1 = -1.0d0*lm/cm - b11*lm*gn/(cm*cm) + b11*gn/(cm*lm) + b11*b11*gn*gn/(lm*cm*cm)
!            c2 = -1.0d0*lm/cm - b11*lm*gn/(cm*cm) + b11*gn/(cm*lm) + b11*b11*gn*gn/(lm*cm*cm)
!            c3 = -2.0d0*gn/cm - b11*gn*gn/(cm*cm) + gn/(lm*lm) + 2.0d0*b11*gn*gn/(lm*lm*cm) + b11*b11*gn*gn*gn/(lm*lm*cm*cm)
!            c4 = -1.0d0*b11*lm*lm/(cm*cm) + b11*b11*gn/(cm*cm)
!            
!            !Linearized variation of the gap
!            if(nlgeom==.false.)then
!                Kdelta = 0.0d0
!            else
!                Kdelta = 0.0d0
!                Kdelta = c1*matmul(N0s,transpose( Ts)) + &
!                &        c2*matmul( Ts,transpose(N0s)) + &
!                &        c3*matmul(N0s,transpose(N0s)) + &
!                &        c4*matmul( Ts,transpose( Ts))
!                
!                Kdelta = Kdelta*Lagrange(i)*-1.0d0
!                
!            end if
!            
!            Kdelta = 0.0d0
!            
!            iLM = iLM + 1
!            
!            continue 
!            
!            !------------------------------------------------------------------
!            ! LAGRANGE MULTIPLIER METHOD - Normal Contact
!            !------------------------------------------------------------------
!            ctconn = 0
!            ctconn(1) = conn_slave(i)
!            do j=1,ncp_m
!                ctconn(j+1)=conn_master(j)
!            end do
!            
!            !------------------------------------------------------------------
!            !Contact contribution from the Greville points to the slave points
!            !------------------------------------------------------------------
!            Kc = 0.0d0
!            Fc = 0.0d0
!            
!            Kc(1,1) = Kdelta(1,1)
!            Kc(1,2) = Kdelta(1,2)
!            Kc(2,1) = Kdelta(2,1)
!            Kc(2,2) = Kdelta(2,2)
!            
!            Kc(1,3) = Ns(1,1)
!            Kc(2,3) = Ns(2,1)
!            
!            Kc(3,1) = Ns(1,1)
!            Kc(3,2) = Ns(2,1)
!            
!            Fc(1,1) = Ns(1,1)
!            Fc(2,1) = Ns(2,1)
!            
!            
!!            count1 = 0
!!            count2 = 0
!!            Kcc = 0.0d0
!!            Fcc = 0.0d0
!!            do k1=1,ncp_s
!!                do k2=1,ncp_s
!!                    count1 = (k1-1)*nds
!!                    do k3=1,nds
!!                        count1 = count1+1
!!                        count2 = (k2-1)*nds
!!                        do k4=1,nds
!!                            count2=count2 + 1
!!                            Kcc(count1,count2)= Kc(k3,k4)*Rs(k1)*Rs(k2)
!!                            continue
!!                        end do
!!                        
!!                         Kcc(ncp_s*nds+1,count1) = Kc(nds+1,k3)*Rs(k1)
!!                         Kcc(count1,ncp_s*nds+1) = Kc(k3,nds+1)*Rs(k1)
!!                            
!!                         Fcc(count1,1) = Fc(k3,1)*Rs(k1)
!!                         continue
!!                    end do
!!                end do
!!            end do
!            
!            
!            Kcc = 0.0d0
!            Fcc = 0.0d0
!            count1 = 1
!            do k1=1,2*ncp_s,2
!                count2 = 0
!                do k2=1,2*ncp_s,2
!                    count2 = count2 + 1
!                    Kcc(k1  ,k2  ) = Kc(1,1)*Rs(count1)*Rs(count2)
!                    Kcc(k1+1,k2  ) = Kc(2,1)*Rs(count1)*Rs(count2)
!                    Kcc(k1  ,k2+1) = Kc(1,2)*Rs(count1)*Rs(count2)
!                    Kcc(k1+1,k2+1) = Kc(2,2)*Rs(count1)*Rs(count2)
!                    
!                    !Kcc(2*ncp_s+1,k2  ) = Kc(3,1)*Rs(count2) !*Rs(ncp_s)
!                    !Kcc(2*ncp_s+1,k2+1) = Kc(3,2)*Rs(count2) !*Rs(ncp_s)
!                    continue
!                end do
!                
!                Kcc(2*ncp_s+1,k1  ) = Kc(3,1)*Rs(count1) !*Rs(ncp_s)
!                Kcc(2*ncp_s+1,k1+1) = Kc(3,2)*Rs(count1) !*Rs(ncp_s)
!                
!                Kcc(k1  ,2*ncp_s+1) = Kc(1,3)*Rs(count1) !*Rs(ncp_s)
!                Kcc(k1+1,2*ncp_s+1) = Kc(2,3)*Rs(count1) !*Rs(ncp_s)
!                
!                Fcc(k1  ,1) = Fc(1,1)*Rs(count1)
!                Fcc(k1+1,1) = Fc(2,1)*Rs(count1)
!                
!                count1 = count1 + 1
!                
!                continue
!                
!            end do
!            
!            do k1=1,ncp_s
!                do k2=1,ncp_s
!                    KLM(conn_slave(k1)*2-1,conn_slave(k2)*2-1) = KLM(conn_slave(k1)*2-1,conn_slave(k2)*2-1) + Kcc(k1*2-1,k2*2-1)
!                    KLM(conn_slave(k1)*2-1,conn_slave(k2)*2  ) = KLM(conn_slave(k1)*2-1,conn_slave(k2)*2  ) + Kcc(k1*2-1,k2*2  )
!                    KLM(conn_slave(k1)*2  ,conn_slave(k2)*2-1) = KLM(conn_slave(k1)*2  ,conn_slave(k2)*2-1) + Kcc(k1*2  ,k2*2-1)
!                    KLM(conn_slave(k1)*2  ,conn_slave(k2)*2  ) = KLM(conn_slave(k1)*2  ,conn_slave(k2)*2  ) + Kcc(k1*2  ,k2*2  )
!                end do
!                
!                KLM(conn_slave(k1)*2-1,tnodes*nds+i) =  KLM(conn_slave(k1)*2-1,tnodes*nds+i) + Kcc(k1*2-1,ncp_s*nds+1)
!                KLM(conn_slave(k1)*2  ,tnodes*nds+i) =  KLM(conn_slave(k1)*2  ,tnodes*nds+i) + Kcc(k1*2  ,ncp_s*nds+1)
!                
!                KLM(tnodes*nds+i,conn_slave(k1)*2-1) =  KLM(tnodes*nds+i,conn_slave(k1)*2-1) + Kcc(ncp_s*nds+1,k1*2-1)
!                KLM(tnodes*nds+i,conn_slave(k1)*2  ) =  KLM(tnodes*nds+i,conn_slave(k1)*2  ) + Kcc(ncp_s*nds+1,k1*2  )
!                
!                FLM(conn_slave(k1)*nds-1,i) = FLM(conn_slave(k1)*nds-1,i) + Fcc(k1*nds-1,1)
!                FLM(conn_slave(k1)*nds  ,i) = FLM(conn_slave(k1)*nds  ,i) + Fcc(k1*nds  ,1)
!                
!            end do
!            
!            continue
!
!            !------------------------------------------------------------------
!            !Contact contribution to the master control points
!            !------------------------------------------------------------------
!!            do k1=2,ncp_m+1
!!                do k2=2,ncp_m+1
!!                    KLM(ctconn(k1)*2-1,ctconn(k2)*2-1) = KLM(ctconn(k1)*2-1,ctconn(k2)*2-1) + Kdelta(k1*2-1,k2*2-1)
!!                    KLM(ctconn(k1)*2-1,ctconn(k2)*2  ) = KLM(ctconn(k1)*2-1,ctconn(k2)*2  ) + Kdelta(k1*2-1,k2*2  )
!!                    KLM(ctconn(k1)*2  ,ctconn(k2)*2-1) = KLM(ctconn(k1)*2  ,ctconn(k2)*2-1) + Kdelta(k1*2  ,k2*2-1)
!!                    KLM(ctconn(k1)*2  ,ctconn(k2)*2  ) = KLM(ctconn(k1)*2  ,ctconn(k2)*2  ) + Kdelta(k1*2  ,k2*2  )
!!                end do
!!            end do
!!            
!!            do k1=2,ncp_m+1
!!                KLM(ctconn(k1)*2-1,tnodes*nds+i) = KLM(ctconn(k1)*2-1,tnodes*nds+i) + Ns(k1*2-1,1)
!!                KLM(ctconn(k1)*2  ,tnodes*nds+i) = KLM(ctconn(k1)*2  ,tnodes*nds+i) + Ns(k1*2  ,1)
!!                
!!                KLM(tnodes*nds+i,ctconn(k1)*2-1) = KLM(tnodes*nds+i,ctconn(k1)*2-1) + Ns(k1*2-1,1)
!!                KLM(tnodes*nds+i,ctconn(k1)*2  ) = KLM(tnodes*nds+i,ctconn(k1)*2  ) + Ns(k1*2  ,1)
!!            end do
!!            
!!            do k1=2,ncp_m+1
!!                FLM(ctconn(k1)*nds-1,i) = FLM(ctconn(k1)*nds-1,i) + Ns(k1*2-1,1)
!!                FLM(ctconn(k1)*nds  ,i) = FLM(ctconn(k1)*nds  ,i) + Ns(k1*2  ,1)
!!            end do
!            
!            GLM(tnodes*nds+i,1) = -1.0d0*final_gap(i)
!            
!            continue
!        
!        else
!            
!            99 continue
!            
!            iLM = iLM + 1
!            KLM(tnodes*nds+i,:) = 0.0d0
!            KLM(:,tnodes*nds+i) = 0.0d0
!            KLM(tnodes*nds+i,tnodes*nds+i) = 1.0d0
!            !FLM(:,i) = 0.0d0
!            !Lagrange(i) = 0.0d0
!            
!            ctconn = 0
!            ctconn(1) = conn_slave(i)
!            do j=1,ncp_m
!                ctconn(j+1)=conn_master(j)
!            end do
!            
!            do k1=1,ncp_s
!                FintLM(conn_slave(k1)*nds-1,1) = 0.0d0
!                FintLM(conn_slave(k1)*nds  ,1) = 0.0d0
!            end do
!              
!!            do k1=2,ncp_m+1
!!                FintLM(ctconn(k1)*nds-1,1) = 0.0d0
!!                FintLM(ctconn(k1)*nds  ,1) = 0.0d0
!!            end do
!            
!        end if !gn .lt. 0.0d0
!        
!    end do !Greville points cycle
!    
!    
!    if(iLM .gt. 0)then
!             
!        !if(iter==1 .and. inc == 1)then
!        allocate(KfeqLM(tnodes*nds-nbc+ngrv,tnodes*nds-nbc+ngrv))
!        allocate(FeqLM(tnodes*nds-nbc+ngrv,1),dispLM(tnodes*nds-nbc+ngrv,1))
!        allocate(disptempLM(tnodes*nds+ngrv),impdispLM(tnodes*nds+ngrv,1))
!        !end do
!        
!        !Store Contact Matrix ----
!        KCont = 0.0d0
!        KCont(1:tnodes*nds,1:tnodes*nds) =  KCont(1:tnodes*nds,1:tnodes*nds) + KLM(1:tnodes*nds,1:tnodes*nds)
!        
!        !Total Stiffness matrix -----
!        KLM(1:tnodes*nds,1:tnodes*nds) = KLM(1:tnodes*nds,1:tnodes*nds) + Kf(1:tnodes*nds,1:tnodes*nds)
!        KT(1:tnodes*nds,1:tnodes*nds) = KT(1:tnodes*nds,1:tnodes*nds) + KLM(1:tnodes*nds,1:tnodes*nds)
!        
!        !Right-hand side ------
!!        FintLM = 0.0d0
!!        count = 0
!!        do i=1,ngrv
!!                
!!            count = count + 1
!!            
!!            ctconn = 0
!!            ctconn(1) = conn_slave(i)
!!            do j=1,ncp_m
!!                ctconn(j+1)=conn_master(j)
!!            end do
!!            
!!            do k1=1,ncp_s
!!                FintLM(conn_slave(k1)*nds-1,1) = FintLM(conn_slave(k1)*nds-1,1) + FLM(conn_slave(k1)*nds-1,i)*Lagrange(i)
!!                FintLM(conn_slave(k1)*nds  ,1) = FintLM(conn_slave(k1)*nds  ,1) + FLM(conn_slave(k1)*nds  ,i)*Lagrange(i)
!!            end do
!!              
!!            do k1=2,ncp_m+1
!!                FintLM(ctconn(k1)*nds-1,1) = FintLM(ctconn(k1)*nds-1,1) + FLM(ctconn(k1)*nds-1,i)*Lagrange(i)
!!                FintLM(ctconn(k1)*nds  ,1) = FintLM(ctconn(k1)*nds  ,1) + FLM(ctconn(k1)*nds  ,i)*Lagrange(i)
!!            end do
!!            
!!            continue
!!        end do
!        
!        continue
!    
!        GLM(1:tnodes*nds,1) = GLM(1:tnodes*nds,1) + Fext(1:tnodes*nds,1) !- FintLM(1:tnodes*nds,1)
!        !GLM(1:tnodes*nds,1) = GLM(1:tnodes*nds,1) + Fini(1:tnodes*nds,1)/(1.0d0*incmax)*(1.0d0*inc) - Fint(1:tnodes*nds,1) - FintLM(1:tnodes*nds,1)
!        
!        disptempLM = 1.0d0
!        disptempLM(1:tnodes*nds) = redcount(1:tnodes*nds)
!        
!        ImpDispLM = 0.0d0
!        ImpDispLM(1:tnodes*nds,1) = ImpDisp(1:tnodes*nds,1)
!        do i=1,ngrv
!            if(PTS_active(i) == 0)then
!                ImpDispLM(tnodes*nds+i,1) = -1.0d0*Lagrange(i)
!            end if
!        end do
!        
!!        count = 0
!!        do i=1,tnodes*nds
!!            count = count + 1
!!            if(BCdof(i) == 0) then
!!                KLM(i,:) = 0.0d0
!!                KLM(:,i) = 0.0d0
!!                redcount(count) = 0
!!            end if
!!        end do
!        
!        !Prescribed Lagrange Multipliers -----
!        do i=tnodes*nds+1,tnodes*nds+iLM
!            if(ImpDispLM(i,1)/=0.0d0)then
!                do k1=tnodes*nds,tnodes*nds+iLM
!                    if(i /= k1) then
!                        GLM(k1,1) = GLM(k1,1) - KLM(k1,i)*ImpDispLM(i,1)
!                        KLM(i,k1) = 0.0d0
!                        KLM(k1,i) = 0.0d0
!                     else 
!                        GLM(i,1) = KLM(i,i)*ImpDispLM(i,1)
!                        continue
!                    end if
!                end do
!            end if
!        end do
!        
!        
!        
!        k=1
!        m=1
!        KfeqLM = 0.0d0
!        FEqLM = 0.0d0
!        do i=1,tnodes*nds + ngrv
!            do j=1,tnodes*nds + ngrv
!                if (disptempLM(i)/=0.0d0 .and. disptempLM(j)/=0.0d0)then
!                    KfeqLM(k,m)=KLM(i,j)
!                    FEqLM(k,1)=GLM(i,1)                
!                    m=m+1
!                    if(m==tnodes*nds-nBC+ngrv+1)then
!                        k=k+1
!                        m=1
!                    end if
!                end if
!            end do
!        end do
!        
!        continue
!        
!        temp1=tnodes*nds-nBC+ngrv
!        dispLM = 0.0d0
!        call Gauss (temp1,KfeqLM,FEqLM,dispLM)
!        
!!        allocate(tempinv(temp1))
!!        tempinv = 0.0d0
!!        ilixo = 0
!!        call gaussj(KfeqLM,temp1,tempinv,ilixo)
!!        deallocate(tempinv)
!!        dispLM = matmul(KfeqLM,FEqLM)
!        
!        !dDisp = dDisp - ddDisp
!        
!        ddDisp = 0.0d0
!        j=1
!        do i=1,tnodes*nds
!            if(redcount(i)==1.0d0)then
!                ddDisp(i,1)=dispLM(j,1)
!                j=j+1
!            end if
!        end do
!        
!        do i=1,ngrv
!            
!            if(dispLM(tnodes*nds-nbc+i,1) .gt. 0.0d0)then
!                if(dispLM(tnodes*nds-nbc+i,1) .gt. 0.0d0 .and. (abs(dispLM(tnodes*nds-nbc+i,1)) .gt. abs(Lagrange(i)))) then
!                    dispLM(tnodes*nds-nbc+i,1) = -1.0d0*Lagrange(i)
!                    PTS_active(i) = 0
!                end if
!            end if
!            
!            Lagrange(i) = Lagrange(i) + dispLM(tnodes*nds-nbc+i,1)
!            
!            if(Lagrange(i) .gt. 0.0d0)then
!                Lagrange(i) = 0.0d0
!            end if
!            
!            if(Lagrange(i) == 0.0d0)then
!                FLM(:,i) = 0.0d0
!            end if
!            
!            dLagrange(i) = dispLM(tnodes*nds-nbc+i,1)
!            
!            LagrangeConv(i) = LagrangeConv(i) + dispLM(tnodes*nds-nbc+i,1)
!            
!        end do
!         
!        
!!        !Compute Contact Stress based on the area associated with the collocation point
!!        PTS_Stress = 0.0d0
!!        do i=1,ngrv
!!            !PTS_Stress(i) = PTS_Stress(i) + dLagrange(i)/CP_area(i)
!!            PTS_Stress(i) = Lagrange(i)/CP_area(i)
!!        end do 
!         
!               
!!        FS = 0.0d0
!!        do i=1,ncp_s
!!            do j=1, ngrv
!!                !Compute B-Splines basis functions...
!!                call BSplineBasisAndDeriv2(ncp_s,p_s,grev(j),u_knot_slave,Rss,dRss,ddRss,kspan)
!!                !... and convert to NURBS basis, including vanishing terms 
!!                call BSplinesToNURBS(ncp_s,p_s,grev(j),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
!!                
!!                Fs(i) = Fs(i) + Rs(i)*Lagrange(j)
!!            end do
!!        end do
!!        
!!        Fm = 0.0d0
!!        do i=1,ncp_m
!!            do j=1, ngrv
!!                call BSplineBasisAndDeriv2(ncp_m,p_m,str_xib(j),u_knot_master,Rmm,dRmm,ddRmm,kspan)
!!                call BSplinesToNURBS(ncp_m,p_m,str_xib(j),u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
!!                
!!                Fm(i) = Fm(i) + Rm(i)*Lagrange(j)
!!                
!!                tang = 0.0d0
!!                do k=1,ncp_m
!!                    tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
!!                    tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
!!                end do
!!            
!!                !Length of the master segment and its derivative
!!                lm  = dsqrt( tang(1,1)*tang(1,1) + tang(2,1)*tang(2,1))
!!            
!!                !Normalise tangent vectors
!!                tang = tang/lm
!!            end do
!!        end do
!        
!        if(sum(Lagrange) == 0.0d0)then
!            FintLM = 0.0d0
!        end if
!        
!        sumLag = sum(Lagrange)
!        
!        sumRm = sumRm/sum(sumRm)
!        
!        !Eliminate Contact Forces in master control points outside contact area
!        do k1=1,ncp_m
!            if(sumRm(k1) == 0.0d0)then
!                FintLM(conn_master(k1)*nds-1,1) = 0.0d0
!                FintLM(conn_master(k1)*nds  ,1) = 0.0d0
!            end if
!        end do
!        
!        !Right-hand side -----------------------------------------------------------
!        FintLM = 0.0d0
!        count = 0
!        do i=1,ngrv
!                
!            count = count + 1
!            
!            ctconn = 0
!            ctconn(1) = conn_slave(i)
!            do j=1,ncp_m
!                ctconn(j+1)=conn_master(j)
!            end do
!            
!            do k1=1,ncp_s
!                FintLM(conn_slave(k1)*nds-1,1) = FintLM(conn_slave(k1)*nds-1,1) + FLM(conn_slave(k1)*nds-1,i)*Lagrange(i)
!                FintLM(conn_slave(k1)*nds  ,1) = FintLM(conn_slave(k1)*nds  ,1) + FLM(conn_slave(k1)*nds  ,i)*Lagrange(i)
!            end do
!            
!!            do k1=2,ncp_m+1
!!                FintLM(ctconn(k1)*nds-1,1) = FintLM(ctconn(k1)*nds-1,1) + FLM(ctconn(k1)*nds-1,i)*Lagrange(i)
!!                FintLM(ctconn(k1)*nds  ,1) = FintLM(ctconn(k1)*nds  ,1) + FLM(ctconn(k1)*nds  ,i)*Lagrange(i)
!!            end do
!            
!            continue
!        end do
!        
!        continue
!
!!        !Right-hand side ------
!!        !FintLM = 0.0d0
!!        count = 0
!!        do i=1,ngrv
!!                
!!            count = count + 1
!!            
!!            ctconn = 0
!!            ctconn(1) = conn_slave(i)
!!            do j=1,ncp_m
!!                ctconn(j+1)=conn_master(j)
!!            end do
!!            
!!            do k1=1,ncp_s
!!                FintLM(conn_slave(k1)*nds-1,1) = FintLM(conn_slave(k1)*nds-1,1) + FLM(conn_slave(k1)*nds-1,i)*dLagrange(i)
!!                FintLM(conn_slave(k1)*nds  ,1) = FintLM(conn_slave(k1)*nds  ,1) + FLM(conn_slave(k1)*nds  ,i)*dLagrange(i)
!!            end do
!!              
!!            continue
!!        end do
!!        
!!        continue
!!        
!!        normx = 0.0d0
!!        normy = 0.0d0
!!        do i=1,ngrv
!!            normx(i) = FintLM(conn_slave(i)*nds-1,1)
!!            normy(i) = FintLM(conn_slave(i)*nds  ,1)
!!        end do
!!        
!!        Lagmx = 0.0d0
!!        Lagmy = 0.0d0
!!        do i=1,ngrv
!!            
!!            !if (PTS_active(i)==1 .and. iter .gt. 1)then
!!            if (PTS_active(i)==1)then
!!                !Compute B-Splines basis functions...
!!                call BSplineBasisAndDeriv2(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan)
!!                !... and convert to NURBS basis, including vanishing terms 
!!                call BSplinesToNURBS(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
!!                
!!                do j=1,ncp_s
!!                    Lagmx(i) = Lagmx(i) - normx(i)*Rs(j)
!!                    Lagmy(i) = Lagmy(i) - normy(i)*Rs(j)
!!                end do
!!        
!!            end if
!!        end do
!!        
!!        continue
!!        
!!        !Right-hand side ------
!!        !FintLM = 0.0d0
!!        count = 0
!!        do i=1,ngrv
!!                
!!            count = count + 1
!!            ctconn = 0
!!            ctconn(1) = conn_slave(i)
!!            do j=1,ncp_m
!!                ctconn(j+1)=conn_master(j)
!!            end do
!!              
!!            do k1=2,ncp_m+1
!!                FintLM(ctconn(k1)*nds-1,1) = 0.0d0
!!                FintLM(ctconn(k1)*nds  ,1) = 0.0d0
!!            end do
!!        end do
!!        
!!        continue
!!        
!!        do i=1,ngrv
!!                
!!            count = count + 1
!!            
!!            ctconn = 0
!!            ctconn(1) = conn_slave(i)
!!            do j=1,ncp_m
!!                ctconn(j+1)=conn_master(j)
!!            end do
!!              
!!            do k1=2,ncp_m+1
!!                FintLM(ctconn(k1)*nds-1,1) = FintLM(ctconn(k1)*nds-1,1) - n_str((k1-1)*nds-1,i)*Lagmx(i)
!!                FintLM(ctconn(k1)*nds  ,1) = FintLM(ctconn(k1)*nds  ,1) - n_str((k1-1)*nds  ,i)*Lagmy(i)
!!            end do
!!            
!!            continue
!!        end do
!!        
!!        continue
!        
!        
!!        Fslv = 0.0d0
!!        do i=1,ngrv
!!            call BSplineBasisAndDeriv2(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan)
!!            call BSplinesToNURBS(      ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
!!            
!!            do j=1,ncp_s
!!                Fslv(i,1) = Fslv(i,1) + FintLM(conn_slave(j)*nds-1,1)*Rs(j)
!!                Fslv(i,2) = Fslv(i,2) + FintLM(conn_slave(j)*nds  ,1)*Rs(j)
!!            end do
!!            
!!        end do
!!        
!!        do k1=1,ncp_m
!!            FintLM(conn_master(k1)*nds-1,1) = 0.0d0
!!            FintLM(conn_master(k1)*nds  ,1) = 0.0d0
!!        end do
!!        
!!        do i=1,ngrv
!!                
!!            count = count + 1
!!            
!!            ctconn = 0
!!            ctconn(1) = conn_slave(i)
!!            do j=1,ncp_m
!!                ctconn(j+1)=conn_master(j)
!!            end do
!!            
!!            call BSplineBasisAndDeriv2(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan)
!!            call BSplinesToNURBS(      ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
!!            
!!            call BSplineBasisAndDeriv2(ncp_m,p_m,str_xib(i),u_knot_master,Rmm,dRmm,ddRmm,kspan)
!!            call BSplinesToNURBS(ncp_m,p_m,str_xib(i),u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
!!            
!!            do k1=2,ncp_m+1
!!                FintLM(ctconn(k1)*nds-1,1) = FintLM(ctconn(k1)*nds-1,1) - Fslv(i,1)*Rm(k1-1)
!!                FintLM(ctconn(k1)*nds  ,1) = FintLM(ctconn(k1)*nds  ,1) - Fslv(i,2)*Rm(k1-1)
!!            end do
!!            
!!            continue
!!        end do
!        
!!        Fmst = 0.0d0
!!        do i=1,ngrv
!!            
!!            call BSplineBasisAndDeriv2(ncp_m,p_m,str_xib(i),u_knot_master,Rmm,dRmm,ddRmm,kspan)
!!            call BSplinesToNURBS(ncp_m,p_m,str_xib(i),u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
!!            
!!            do j=1,ncp_m
!!                Fmst(j,1) = Fmst(j,1) - Fslv(i,1)*Rm(j)
!!                Fmst(j,1) = Fmst(j,1) - Fslv(i,1)*Rm(j)
!!            end do
!!        end do
!!        
!!        continue
!        
!        !FintLM = 0.0d0
!        
!!        count = 0
!!        do i=1,ngrv
!!                
!!            count = count + 1
!!            
!!            ctconn = 0
!!            !ctconn(1) = conn_slave(i)
!!            do j=1,ncp_m
!!                ctconn(j+1)=conn_master(j)
!!            end do
!!            
!!            do k1=1,ncp_s
!!                FintLM(conn_slave(k1)*nds-1,1) = FintLM(conn_slave(k1)*nds-1,1) + FLM(conn_slave(k1)*nds-1,i)*dispLM(tnodes*nds-nbc+i,1)
!!                FintLM(conn_slave(k1)*nds  ,1) = FintLM(conn_slave(k1)*nds  ,1) + FLM(conn_slave(k1)*nds  ,i)*dispLM(tnodes*nds-nbc+i,1)
!!            end do
!!              
!!            do k1=2,ncp_m+1
!!                FintLM(ctconn(k1)*nds-1,1) = FintLM(ctconn(k1)*nds-1,1) + FLM(ctconn(k1)*nds-1,i)*dispLM(tnodes*nds-nbc+i,1)
!!                FintLM(ctconn(k1)*nds  ,1) = FintLM(ctconn(k1)*nds  ,1) + FLM(ctconn(k1)*nds  ,i)*dispLM(tnodes*nds-nbc+i,1)
!!            end do
!!            
!!            continue
!!        end do
!!        
!!        continue
!!        
!!        !dDisp = dDisp + ddDisp
!!        continue
!
!    end if
!    
!    
!    continue
!    
!end subroutine
