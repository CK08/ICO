subroutine PTS_ComputeGapPenal(iter,inc,istep,iLM)

    use Mod_Variables
    implicit none
    
    integer(4),intent(IN)::iter,inc,istep
    integer(4),intent(OUT)::iLM
    
    integer(4)::i,j,k,l,m,k1,k2,k3,k4
    integer(4)::temp1,count,count1,count2,kspan,ilixo
    integer(4),parameter::knewton = 50
    integer(4),dimension(ncp_m+1)::ctconn
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
    real(8),dimension((ncp_m+1)*nds,1)::Ns,N0s,Ts,T0s
    real(8),dimension((ncp_m+1)*nds,(ncp_m+1)*nds)::Kdelta
    real(8),dimension(:,:),allocatable::FLM,FLMStr,KLM,LMStr,GLM,GLMStr,Ki
    real(8),dimension(:,:),allocatable::KfeqLM,FeqLM,dispLM
    real(8),dimension(:),allocatable::disptempLM
    real(8),dimension(:,:),allocatable::ImpDispLM
    
    real(8),dimension(nds,nds)::Kc
    real(8),dimension(nds,1)::Fc
    
    real(8),dimension(ncp_s*nds,ncp_s*nds)::Kcc
    real(8),dimension(ncp_s*nds,1)::Fcc
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
    logical::tag
    
    real(8)::ep
    real(8),dimension((ncp_m+1)*nds,1)::FPenal
    real(8),dimension((ncp_m+1)*nds,(ncp_m+1)*nds)::KPenal
    
    
    allocate(KLM(tnodes*nds,tnodes*nds))
    KLM = 0.0d0
    allocate(GLM(tnodes*nds,1))
    GLM = 0.0d0
    
    allocate(FLM(tnodes*nds,1))
    FLM = 0.0d0
    
    allocate(n_str((ncp_m)*nds,1))
    n_str = 0.0d0
    
    iLM = 0
    
    ep = 1.0d5
    
    sumRm = 0.0d0
    
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
            xs(1,1) = xs(1,1) + Rs(j)*(GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*2-1),1))
            xs(2,1) = xs(2,1) + Rs(j)*(GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*2  ),1))
        end do
        
        continue
        
        !------------------------------------------------------------------
        ! Contact check (maybe)
        !------------------------------------------------------------------
        xibi = 0.0d0
        !Compute B-Splines basis functions...
        call BSplineBasisAndDeriv2(ncp_m,p_m,xibi,u_knot_master,Rmm,dRmm,ddRmm,kspan)
        
        !... and convert to NURBS basis, including vanishing terms 
        call BSplinesToNURBS(ncp_m,p_m,xibi,u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
        
        !Tangent vector in the master segment and its derivative
        tang = 0.0d0
        dtang = 0.0d0
        do k=1,ncp_m
            tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
            tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
        end do
        
        !Length of the master segment and its derivative
        lm  = dsqrt( tang(1,1)*tang(1,1) + tang(2,1)*tang(2,1))
        
        !Normalise tangent vectors
        tang = tang/lm
        
        !Possible contact check
        x11(1,1) = GCoords(conn_master(1),1) + ddisp(conn_master(1)*2-1,1)
        x11(2,1) = GCoords(conn_master(1),2) + ddisp(conn_master(1)*2  ,1)
        gapi = xs - x11
        gti = gapi(1,1)*tang(1,1) + gapi(2,1)*tang(2,1)
        
        xibf = 1.0d0
        !Compute B-Splines basis functions...
        call BSplineBasisAndDeriv2(ncp_m,p_m,xibf,u_knot_master,Rmm,dRmm,ddRmm,kspan)
        
        !... and convert to NURBS basis, including vanishing terms 
        call BSplinesToNURBS(ncp_m,p_m,xibf,u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
        
        !Tangent vector in the master segment
        tang = 0.0d0
        dtang = 0.0d0
        do k=1,ncp_m
            tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
            tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
        end do
        
        !Length of the master segment
        lm  = dsqrt( tang(1,1)*tang(1,1) + tang(2,1)*tang(2,1))
        
        !Normalise tangent vector
        tang = tang/lm
        
        !------------------------------------------------------------------
        ! Secant Method to determine the Closest Point Projection
        !------------------------------------------------------------------
        xib0 = 0.0d0
        xib = 0.5d0
        do j=1,knewton
            if(j==1)then
                !Compute B-Splines basis functions...
                call BSplineBasisAndDeriv2(ncp_m,p_m,xib0,u_knot_master,Rmm,dRmm,ddRmm,kspan)
                
               !... and convert to NURBS basis, including vanishing terms 
                call BSplinesToNURBS(ncp_m,p_m,xib0,u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
                
                !Closest Point Projection and its derivative
                xm = 0.0d0
                dxm = 0.0d0
                do k=1,ncp_m
                    xm(1,1)  = xm(1,1)  + Rm(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
                    xm(2,1)  = xm(2,1)  + Rm(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
                    continue
                end do
                
                !Tangent vector in the master segment and its derivative
                tang = 0.0d0
                do k=1,ncp_m
                    tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
                    tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
                end do
                
                !Length of the master segment and its derivative
                lm  = dsqrt( tang(1,1)*tang(1,1) + tang(2,1)*tang(2,1))
                
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
                xm(1,1)  = xm(1,1)  + Rm(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
                xm(2,1)  = xm(2,1)  + Rm(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
                continue
            end do
            
            !Tangent vector in the master segment and its derivative
            tang = 0.0d0
            do k=1,ncp_m
                tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
                tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
            end do
            
            !Length of the master segment and its derivative
            lm  = dsqrt( tang(1,1)*tang(1,1) + tang(2,1)*tang(2,1))
            
            !Normalise tangent vectors
            tang = tang/lm
            
            continue
            
            !Auxiliary vectors
            vec = xs - xm
            
            rsd = vec(1,1)*tang(1,1)+vec(2,1)*tang(2,1)
            
            tst = rsd*(xib - xib0)/(rsd - rsd0)
            
            xibi = xib
            
            xib = xib - rsd*(xib - xib0)/(rsd - rsd0)
            
            if(xib .gt. 1.0d0) xib = 1.0d0
            if(xib .lt. 0.0d0) xib = 0.0d0
            
            xib0 = xibi
            rsd0 = rsd
            
            if(abs(rsd) .lt. 1.0d-10) goto 1
            
            continue
            
        end do
        
        write(*,FMT=16)i
        16 format('Warning: Failed to determine initial CPP for Greville point ',I3)
        
        1 continue
        
        str_xib(i) = xib
        sumRm = sumRm + Rm
        
        !Compute normal vector from cross product {tang}x{0 0 -1}
        nor(1,1) = -tang(2,1)
        nor(2,1) =  tang(1,1)
        
        !Compute initial normal gap
        initial_gap(i) = final_gap(i) 
        final_gap(i) = vec(1,1)*nor(1,1)+vec(2,1)*nor(2,1)
        
        if(final_gap(i) <= -1.0d-8) then
            PTS_active(i) = 1
        else
            PTS_active(i) = 0
        end if
        
        cfv = -ep*final_gap(i) - MP_props(1,1)*final_gap(i)
        
!        if(cfv .gt. 0.0d0)then
!            PTS_active(i) = 1
!        else
!            PTS_active(i) = 0
!        end if
        
        !If contact has occured
        Kdelta = 0.0d0
        if(PTS_active(i) == 1)then
            
            !Metric of the boundary
            a11 = lm*lm
            
            !Second derivative of the Closest Point Projection
            dpos = 0.0d0
            do k=1,ncp_m
                dpos(1,1) = dpos(1,1) + d2Rmdxi2(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
                dpos(2,1) = dpos(2,1) + d2Rmdxi2(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
                continue
            end do
            
            !Curvature of the boundary
            b11 = dpos(1,1)*nor(1,1) + dpos(2,1)*nor(2,1)
        
            !Auxiliary Matrices
            Ns = 0.0d0
            Ts = 0.0d0
            N0s = 0.0d0
            !T0s = 0.0d0
            
            do j=1,nds
                Ns(j,1) = nor(j,1)
                Ts(j,1) = tang(j,1)
            end do
            
            do j=1,ncp_m
                Ns(nds+j*2-1,1) = -1.0d0*Rm(j)*nor(1,1)
                Ns(nds+j*2  ,1) = -1.0d0*Rm(j)*nor(2,1)
                
                N0s(nds+j*2-1,1) = dRmdxi(j)*nor(1,1)
                N0s(nds+j*2  ,1) = dRmdxi(j)*nor(2,1)
                
                Ts(nds+j*2-1,1) = -1.0d0*Rm(j)*tang(1,1)
                Ts(nds+j*2  ,1) = -1.0d0*Rm(j)*tang(2,1)
                
                continue
            end do
            
            !Coefficients for the contact stiffness matrices
            cm = a11 - gn*b11
            c1 = -1.0d0*lm/cm - b11*lm*gn/(cm*cm) + b11*gn/(cm*lm) + b11*b11*gn*gn/(lm*cm*cm)
            c2 = -1.0d0*lm/cm - b11*lm*gn/(cm*cm) + b11*gn/(cm*lm) + b11*b11*gn*gn/(lm*cm*cm)
            c3 = -2.0d0*gn/cm - b11*gn*gn/(cm*cm) + gn/(lm*lm) + 2.0d0*b11*gn*gn/(lm*lm*cm) + b11*b11*gn*gn*gn/(lm*lm*cm*cm)
            c4 = -1.0d0*b11*lm*lm/(cm*cm) + b11*b11*gn/(cm*cm)
            
            KPenal = 0.0d0
            KPenal = ep*matmul(Ns,transpose(Ns))
            
            FPenal = ep*Ns
            
            
            iLM = iLM + 1
            
            continue 
            
            !------------------------------------------------------------------
            ! LAGRANGE MULTIPLIER METHOD - Normal Contact
            !------------------------------------------------------------------
            ctconn = 0
            ctconn(1) = conn_slave(i)
            do j=1,ncp_m
                ctconn(j+1)=conn_master(j)
            end do
            
            !------------------------------------------------------------------
            !Contact contribution from the Greville points to the slave points
            !------------------------------------------------------------------
            Kc = 0.0d0
            Fc = 0.0d0
            
            Kc(1,1) = KPenal(1,1)
            Kc(1,2) = KPenal(1,2)
            Kc(2,1) = KPenal(2,1)
            Kc(2,2) = KPenal(2,2)
            
!            Kc(1,3) = Ns(1,1)
!            Kc(2,3) = Ns(2,1)
!            
!            Kc(3,1) = Ns(1,1)
!            Kc(3,2) = Ns(2,1)
!            
            Fc(1,1) = Ns(1,1)*ep*final_gap(i)
            Fc(2,1) = Ns(2,1)*ep*final_gap(i)
            
            
            Kcc = 0.0d0
            Fcc = 0.0d0
            count1 = 1
            do k1=1,2*ncp_s,2
                count2 = 0
                do k2=1,2*ncp_s,2
                    count2 = count2 + 1
                    Kcc(k1  ,k2  ) = Kc(1,1)*Rs(count1)*Rs(count2)
                    Kcc(k1+1,k2  ) = Kc(2,1)*Rs(count1)*Rs(count2)
                    Kcc(k1  ,k2+1) = Kc(1,2)*Rs(count1)*Rs(count2)
                    Kcc(k1+1,k2+1) = Kc(2,2)*Rs(count1)*Rs(count2)
                    
                    !Kcc(2*ncp_s+1,k2  ) = Kc(3,1)*Rs(count2) !*Rs(ncp_s)
                    !Kcc(2*ncp_s+1,k2+1) = Kc(3,2)*Rs(count2) !*Rs(ncp_s)
                    continue
                end do
                
                Fcc(k1  ,1) = Fc(1,1)*Rs(count1)
                Fcc(k1+1,1) = Fc(2,1)*Rs(count1)
                
                count1 = count1 + 1
                
                continue
                
            end do
            
            
            continue
            
            
            do k1=1,ncp_s
                do k2=1,ncp_s
                    KLM(conn_slave(k1)*2-1,conn_slave(k2)*2-1) = KLM(conn_slave(k1)*2-1,conn_slave(k2)*2-1) + Kcc(k1*2-1,k2*2-1)
                    KLM(conn_slave(k1)*2-1,conn_slave(k2)*2  ) = KLM(conn_slave(k1)*2-1,conn_slave(k2)*2  ) + Kcc(k1*2-1,k2*2  )
                    KLM(conn_slave(k1)*2  ,conn_slave(k2)*2-1) = KLM(conn_slave(k1)*2  ,conn_slave(k2)*2-1) + Kcc(k1*2  ,k2*2-1)
                    KLM(conn_slave(k1)*2  ,conn_slave(k2)*2  ) = KLM(conn_slave(k1)*2  ,conn_slave(k2)*2  ) + Kcc(k1*2  ,k2*2  )
                end do
                
                FLM(conn_slave(k1)*nds-1,1) = FLM(conn_slave(k1)*nds-1,1) + Fcc(k1*nds-1,1)
                FLM(conn_slave(k1)*nds  ,1) = FLM(conn_slave(k1)*nds  ,1) + Fcc(k1*nds  ,1)
                
            end do
            
            continue

            !------------------------------------------------------------------
            !Contact contribution to the master control points
            !------------------------------------------------------------------
            KPenal(1,1) = 0.0d0
            KPenal(2,1) = 0.0d0
            KPenal(1,2) = 0.0d0
            KPenal(2,2) = 0.0d0
            
            do k1=1,ncp_m+1
                do k2=1,ncp_m+1
                    KLM(ctconn(k1)*2-1,ctconn(k2)*2-1) = KLM(ctconn(k1)*2-1,ctconn(k2)*2-1) + KPenal(k1*2-1,k2*2-1)
                    KLM(ctconn(k1)*2-1,ctconn(k2)*2  ) = KLM(ctconn(k1)*2-1,ctconn(k2)*2  ) + KPenal(k1*2-1,k2*2  )
                    KLM(ctconn(k1)*2  ,ctconn(k2)*2-1) = KLM(ctconn(k1)*2  ,ctconn(k2)*2-1) + KPenal(k1*2  ,k2*2-1)
                    KLM(ctconn(k1)*2  ,ctconn(k2)*2  ) = KLM(ctconn(k1)*2  ,ctconn(k2)*2  ) + KPenal(k1*2  ,k2*2  )
                end do
            end do
            
            do k1=2,ncp_m+1
                FLM(ctconn(k1)*nds-1,1) = FLM(ctconn(k1)*nds-1,1) + Ns(k1*2-1,1)*ep*final_gap(i)
                FLM(ctconn(k1)*nds  ,1) = FLM(ctconn(k1)*nds  ,1) + Ns(k1*2  ,1)*ep*final_gap(i)
            end do
            
            continue
        
        else
            
            99 continue
            
            iLM = iLM + 1
            !KLM = 0.0d0
            
            !FintLM = 0.0d0
              
!            do k1=2,ncp_m+1
!                FintLM(ctconn(k1)*nds-1,1) = 0.0d0
!                FintLM(ctconn(k1)*nds  ,1) = 0.0d0
!            end do
            
        end if !gn .lt. 0.0d0
        
    end do !Greville points cycle
    
    
    if(iLM .gt. 0)then
             

        !Store Contact Matrix ----
        !KCont = 0.0d0
        !KCont(1:tnodes*nds,1:tnodes*nds) =  KCont(1:tnodes*nds,1:tnodes*nds) + KLM(1:tnodes*nds,1:tnodes*nds)
        
        !Total Stiffness matrix -----
        !KLM(1:tnodes*nds,1:tnodes*nds) = KLM(1:tnodes*nds,1:tnodes*nds) + Kf(1:tnodes*nds,1:tnodes*nds)
        !KT(1:tnodes*nds,1:tnodes*nds) = KT(1:tnodes*nds,1:tnodes*nds) + KLM(1:tnodes*nds,1:tnodes*nds)
        
        KF = KF + KLM
        
        FintLM = FLM !*-1.0d0
        
        !-----------------------------------------------------------------
        !Assemble equivalent reduced matrices
        !-----------------------------------------------------------------
        Kfeq = 0.0d0
        FextEq = 0.0d0
        k1 = 1
        k2 = 1
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
        end do

        
        !-----------------------------------------------------------------
        !Solve system of equations
        !-----------------------------------------------------------------
        dispeqv = 0.0d0
        call Gauss (tnodes*nds-nbc,Kfeq,FextEq,dispeqv)

        j=1
        ddDisp = 0.0d0
        do i=1,tnodes*nds
            if(redcount(i)==1)then
                ddDisp(i,1)=dispeqv(j,1)
                j=j+1
            end if
        end do

    end if
    
    
    continue
    
end subroutine