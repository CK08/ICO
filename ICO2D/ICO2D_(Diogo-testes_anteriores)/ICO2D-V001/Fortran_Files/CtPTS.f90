!-------------------------------------------------------------------
!
! Subroutine for the point-to-segment contact procedure
!
! Based on the paper by Matzen et al. "A point to segment contat 
! formulation for isogeometric, NURBS based finite elements", 
! CMAME, 255:27-39, 2013
! 
! IMPORTANT NOTE: The master surface control points must be given
! in a counter-clockwise fashion
!
!-------------------------------------------------------------------
subroutine CtPTS(iter,inc,istep,iLM)
    
    use Mod_Variables
    implicit none
    
    integer(4),intent(IN)::iter,inc,istep
    integer(4),intent(OUT)::iLM
    integer(4)::i,j,k,l,m,k1,k2,k3,k4
    integer(4)::temp1,count,count1,count2,kspan,ilixo
    !integer(4)::ngrv
    integer(4),parameter::knewton = 50
    integer(4),dimension(ncp_m+1)::ctconn
    integer(4),dimension(ncp_s)::isct
    real(8)::cm,c1,c2,c3,c4
    
    real(8)::xib,lm,dlm,sumres,rsd,gn,a11,b11,gt,num,tst
    real(8)::xibi,xibf,diff
    !real(8),dimension(:),allocatable::grev
    real(8),dimension(p_s+1)::Rss,dRss,ddRss
    real(8),dimension(p_m+1)::Rmm,dRmm,ddRmm
    real(8),dimension(ncp_s)::Rs,dRsdxi,d2Rsdxi2
    real(8),dimension(ncp_m)::Rm,dRmdxi,d2Rmdxi2
    real(8),dimension(nds,1)::tang,dtang,nor,dpos,aux21
    real(8),dimension(nds,1)::xs,xm,dxm,vec,dvec,vec0
    real(8),dimension(nds,1)::gap,x11,x21,gapi,gapf
    real(8),dimension((ncp_m+1)*nds,1)::Ns,N0s,Ts,T0s
    real(8),dimension((ncp_m+1)*nds,(ncp_m+1)*nds)::Kdelta
    real(8),dimension(:,:),allocatable::FLM,FLMStr,KLM,LMStr,GLM,GLMStr,Ki
    real(8),dimension(:,:),allocatable::KfeqLM,FeqLM,dispLM
    real(8),dimension(:),allocatable::disptempLM
    !real(8),dimension(:),allocatable::initial_gap,final_gap
    
    real(8),dimension(nds+1,nds+1)::Kc
    real(8),dimension(nds+1,1)::Fc
    
    real(8),dimension(ncp_s*nds+1,ncp_s*nds+1)::Kcc
    real(8),dimension(ncp_s*nds+1,1)::Fcc
    
    real(8)::gti,gtf,xib0,rsd0,sumt,sumd,sumdd
    
    real(8),dimension(:),allocatable::temp
    
    !Contact Forces
    real(8)::suml,sumr,dx,dy,sumdl
    real(8),dimension(ncp_s)::CF,bl,br,dl,Rl,Rr
    logical::tag
    
    
    
    
    !if(istep==1 .and. iter==1 .and. inc==1)then
    !    allocate(Lagrange(ncp_s))
    !    Lagrange = 0.0d0
    !    allocate(LagrangeConv(ncp_s))
    !    LagrangeConv = 0.0d0
    !elseif(iter == 1)then
    !    Lagrange = LagrangeConv
    !end if
    
    !Calculate Greville Points on the parametric surface
    !ngrv = ncp_s
    !allocate(grev(ngrv),initial_gap(ngrv),final_gap(ngrv))
    
    allocate(KLM(tnodes*nds+ngrv,tnodes*nds+ngrv))
    KLM = 0.0d0
    allocate(GLM(tnodes*nds+ngrv,1))
    GLM = 0.0d0
    allocate(FLM(tnodes*nds,ngrv))
    FLM = 0.0d0
    
    if(firstinc==.true.)then
        allocate(xibs(ngrv))
        xibs = 0.0d0
        firstinc = .false.
    end if
    
    grev = 0.0d0
    do i=1,ngrv
        do j=1,p_s
            grev(i) = grev(i) + u_knot_slave(i+j)/(1.0d0*p_s)
        end do
    end do
    
    !Set the number of Lagrange Multipliers to zero
    iLM = 0
    
    !Reset info about contact poits
    isct = 0
    
    !Reset initial gap
    initial_gap = 0.0d0
    final_gap = 0.0d0
    
    !Slave (Greville) points cycle
    do i=1,ngrv
        
        !------------------------------------------------------------------
        ! LAGRANGE MULTIPLIER METHOD - Compute initial gap
        !------------------------------------------------------------------
        
        !Compute B-Splines basis functions...
        call BSplineBasisAndDeriv2(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan)
        
        !... and convert to NURBS basis, including vanishing terms 
        call BSplinesToNURBS(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
        
        !Slave (Greville) point in physical slave surface
        xs = 0.0d0
        do j=1,ncp_s
            xs(1,1) = xs(1,1) + Rs(j)*(GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*2-1),1) - ddDisp((conn_slave(j)*2-1),1))
            xs(2,1) = xs(2,1) + Rs(j)*(GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*2  ),1) - ddDisp((conn_slave(j)*2  ),1))
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
            tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1) - ddDisp(conn_master(k)*2-1,1))
            tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1) - ddDisp(conn_master(k)*2  ,1))
        end do
        
        !Length of the master segment and its derivative
        lm  = dsqrt( tang(1,1)*tang(1,1) + tang(2,1)*tang(2,1))
        
        !Normalise tangent vectors
        tang = tang/lm
        
        !Possible contact check
        x11(1,1) = GCoords(conn_master(1),1) + ddisp(conn_master(1)*2-1,1) - ddDisp(conn_master(1)*2-1,1)
        x11(2,1) = GCoords(conn_master(1),2) + ddisp(conn_master(1)*2  ,1) - ddDisp(conn_master(1)*2  ,1)
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
            tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1) - ddDisp(conn_master(k)*2-1,1))
            tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1) - ddDisp(conn_master(k)*2  ,1))
        end do
        
        !Length of the master segment
        lm  = dsqrt( tang(1,1)*tang(1,1) + tang(2,1)*tang(2,1))
        
        !Normalise tangent vector
        tang = tang/lm
        
!        x21(1,1) = GCoords(conn_master(ncp_m),1) + ddisp(conn_master(ncp_m)*2-1,1) - ddDisp(conn_master(ncp_m)*2-1,1)
!        x21(2,1) = GCoords(conn_master(ncp_m),2) + ddisp(conn_master(ncp_m)*2  ,1) - ddDisp(conn_master(ncp_m)*2  ,1)
!        gapf = xs - x21
!        gtf = gapf(1,1)*tang(1,1) + gapf(2,1)*tang(2,1)
!        
!        if(gti .lt. 0.0d0 .and. gtf .lt. 0.0d0 .or. gti .gt. 0.0d0 .and. gtf .gt. 0.0d0)then
!            write(*,FMT=15)i
!            15 format('Greville point',I3,' not in contact')
!            goto 99
!        end if
        
        continue

        
        
!        !------------------------------------------------------------------
!        ! Newton iterative cycle to determine Closest Point Projection
!        !------------------------------------------------------------------
!
!!        xib = 0.0d0
!        !xib = xibs(i)
!        rsd = 0.0d0
!        do j=1,knewton
!            
!            !Determine basis functions and derivatives for the initial parametric coordinate 
!            !of the Closest Point Projection, including vanishing terms
!            call BSplineBasisAndDeriv2(ncp_m,p_m,xib,u_knot_master,Rmm,dRmm,ddRmm,kspan)
!            
!            Rm       = 0.0d0
!            dRmdxi   = 0.0d0
!            d2Rmdxi2 = 0.0d0
!            count = 0
!            do k=kspan-p_m,kspan-p_m+p_m
!                count = count + 1
!                Rm(k)       =   Rmm(count)
!                dRmdxi(k)   =  dRmm(count)
!                d2Rmdxi2(k) = ddRmm(count)
!            end do
!            
!            !Closest Point Projection and its derivative
!            xm = 0.0d0
!            dxm = 0.0d0
!            do k=1,ncp_m
!                xm(1,1)  = xm(1,1)  + Rm(k)*(    GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1) - ddDisp(conn_master(k)*2-1,1))
!                xm(2,1)  = xm(2,1)  + Rm(k)*(    GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1) - ddDisp(conn_master(k)*2  ,1))
!                dxm(1,1) = dxm(1,1) + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1) - ddDisp(conn_master(k)*2-1,1))
!                dxm(2,1) = dxm(2,1) + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1) - ddDisp(conn_master(k)*2  ,1))
!                continue
!            end do
!            
!            !Tangent vector in the master segment and its derivative
!            tang = 0.0d0
!            dtang = 0.0d0
!            do k=1,ncp_m
!                tang(1,1)  = tang(1,1)  + dRmdxi(k)*(  GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1) - ddDisp(conn_master(k)*2-1,1))
!                tang(2,1)  = tang(2,1)  + dRmdxi(k)*(  GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1) - ddDisp(conn_master(k)*2  ,1))
!                dtang(1,1) = dtang(1,1) + d2Rmdxi2(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1) - ddDisp(conn_master(k)*2-1,1))
!                dtang(2,1) = dtang(2,1) + d2Rmdxi2(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1) - ddDisp(conn_master(k)*2  ,1))
!                continue
!            end do
!            
!            !Length of the master segment and its derivative
!            lm  = dsqrt( tang(1,1)*tang(1,1) + tang(2,1)*tang(2,1))
!            dlm = dsqrt(dtang(1,1)*dtang(1,1)+ dtang(2,1)*dtang(2,1))
!            
!            !Normalise tangent vectors
!            tang = tang/lm
!            dtang = dtang/dlm
!            
!            continue
!            
!            !Special case in which the surface has no curvature
!            if(abs(dlm) .lt. 1.0d-7)then
!                
!                x11(1,1) = GCoords(conn_master(1),1) + ddisp(conn_master(1)*2-1,1) - ddDisp(conn_master(1)*2-1,1)
!                x11(2,1) = GCoords(conn_master(1),2) + ddisp(conn_master(1)*2  ,1) - ddDisp(conn_master(1)*2  ,1)
!                
!                gap = xs - x11
!                
!                gt = gap(1,1)*tang(1,1) + gap(2,1)*tang(2,1)
!                        
!                xib = gt/lm
!                
!                call BSplineBasisAndDeriv2(ncp_m,p_m,xib,u_knot_master,Rmm,dRmm,ddRmm,kspan)
!                
!                Rm       = 0.0d0
!                dRmdxi   = 0.0d0
!                d2Rmdxi2 = 0.0d0
!                count = 0
!                do k=kspan-p_m,kspan-p_s+p_m
!                    count = count + 1
!                    Rm(k)       =   Rmm(count)
!                    dRmdxi(k)   =  dRmm(count)
!                    d2Rmdxi2(k) = ddRmm(count)
!                end do
!                
!                continue
!            end if !dlm == 0.0d0
!            
!            !Auxiliary vectors
!            vec = xs - xm
!            !vec = vec/(dsqrt(vec(1,1)*vec(1,1)+vec(2,1)*vec(2,1)))
!            
!            dvec = xs - dxm
!            !dvec = dvec/(dsqrt(dvec(1,1)*dvec(1,1)+dvec(2,1)*dvec(2,1)))
!            
!            if(abs(dlm) .lt. 1.0d-7) goto 1
!            
!            !Residual
!            rsd = vec(1,1)*tang(1,1)+vec(2,1)*tang(2,1)
!            num = dvec(1,1)*dtang(1,1)+dvec(2,1)*dtang(2,1)
!            
!!            aux21 = vec/(dsqrt(vec(1,1)*vec(1,1)+vec(2,1)*vec(2,1))) !-tang
!!            diff = dsqrt(aux21(1,1)*aux21(1,1)+aux21(2,1)*aux21(2,1))
!!            
!            tst = rsd/num
!            
!            !Update Parametric Coordinate of the Closest Point Projection 
!            if(rsd .gt. 0.0d0 .and. num .lt. 0.0d0) then !For a convex master surface
!                xib = xib - rsd/num
!            elseif(rsd .gt. 0.0d0 .and. num .gt. 0.0d0)then !For a concave master surface
!                xib = xib + rsd/num
!            elseif(rsd .lt. 0.0d0 .and. num .lt. 0.0d0)then
!                xib = xib - rsd/num
!            else
!                xib = xib + rsd/num
!            end if
!            
!            if(xib .gt. 1.0d0) xib = 1.0d0
!            if(xib .lt. 0.0d0) xib = 0.0d0
!            
!            if(abs(rsd) .lt. 1.0d-8) goto 1
!            
!            continue
!        end do !Newton iterative cycle
        

        !------------------------------------------------------------------
        ! Secant Method to determine the initial Closest Point Projection
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
                    xm(1,1)  = xm(1,1)  + Rm(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1) - ddDisp(conn_master(k)*2-1,1))
                    xm(2,1)  = xm(2,1)  + Rm(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1) - ddDisp(conn_master(k)*2  ,1))
                    continue
                end do
                
                !Tangent vector in the master segment and its derivative
                tang = 0.0d0
                do k=1,ncp_m
                    tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1) - ddDisp(conn_master(k)*2-1,1))
                    tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1) - ddDisp(conn_master(k)*2  ,1))
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
                xm(1,1)  = xm(1,1)  + Rm(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1) - ddDisp(conn_master(k)*2-1,1))
                xm(2,1)  = xm(2,1)  + Rm(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1) - ddDisp(conn_master(k)*2  ,1))
                continue
            end do
            
            !Tangent vector in the master segment and its derivative
            tang = 0.0d0
            do k=1,ncp_m
                tang(1,1)  = tang(1,1)  + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1) - ddDisp(conn_master(k)*2-1,1))
                tang(2,1)  = tang(2,1)  + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1) - ddDisp(conn_master(k)*2  ,1))
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
        
        !Compute normal vector from cross product {tang}x{0 0 -1}
        nor(1,1) = -tang(2,1)
        nor(2,1) =  tang(1,1)
        
        !Compute initial normal gap
        initial_gap(i) = vec(1,1)*nor(1,1)+vec(2,1)*nor(2,1)
        continue
        
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
            xs(1,1) = xs(1,1) + Rs(j)*(GCoords(conn_slave(j),1) + ddisp(conn_slave(j)*2-1,1))
            xs(2,1) = xs(2,1) + Rs(j)*(GCoords(conn_slave(j),2) + ddisp(conn_slave(j)*2  ,1))
        end do
        
        continue
        
!        !Newton iterative cycle to determine Closest Point Projection
!        xib = 0.0d0
!        !xib = xibs(i)
!        rsd = 0.0d0
!        do j=1,knewton
!            
!            !Determine basis functions and derivatives for the current parametric coordinate 
!            !of the Closest Point Projection, including vanishing terms
!            call BSplineBasisAndDeriv2(ncp_m,p_m,xib,u_knot_master,Rmm,dRmm,ddRmm,kspan)
!            
!            Rm       = 0.0d0
!            dRmdxi   = 0.0d0
!            d2Rmdxi2 = 0.0d0
!            count = 0
!            do k=kspan-p_m,kspan-p_m+p_m
!                count = count + 1
!                Rm(k)       =   Rmm(count)
!                dRmdxi(k)   =  dRmm(count)
!                d2Rmdxi2(k) = ddRmm(count)
!            end do
!            
!            !Closest Point Projection and its derivative
!            xm = 0.0d0
!            dxm = 0.0d0
!            do k=1,ncp_m
!                 xm(1,1) =  xm(1,1) + Rm(k)*(    GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
!                 xm(2,1) =  xm(2,1) + Rm(k)*(    GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
!                dxm(1,1) = dxm(1,1) + dRmdxi(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
!                dxm(2,1) = dxm(2,1) + dRmdxi(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
!                continue
!            end do
!            
!            !Tangent vector in the master segment and its derivative
!            tang = 0.0d0
!            dtang = 0.0d0
!            do k=1,ncp_m
!                 tang(1,1) =  tang(1,1) + dRmdxi(k)*(  GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
!                 tang(2,1) =  tang(2,1) + dRmdxi(k)*(  GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
!                dtang(1,1) = dtang(1,1) + d2Rmdxi2(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
!                dtang(2,1) = dtang(2,1) + d2Rmdxi2(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
!                continue
!            end do
!            
!            !Length of the master segment and its derivative
!            lm  = dsqrt(tang(1,1)*tang(1,1)+tang(2,1)*tang(2,1))
!            dlm = dsqrt(dtang(1,1)*dtang(1,1)+dtang(2,1)*dtang(2,1))
!            
!            !Normalise tangent vectors
!            tang = tang/lm
!            dtang = dtang/dlm
!            
!            !Special case in which the surface has no curvature
!            if(abs(dlm) .lt. 1.0d-7)then
!                
!                x11(1,1) = GCoords(conn_master(1),1) + ddisp(conn_master(1)*nds-1,1)
!                x11(2,1) = GCoords(conn_master(1),2) + ddisp(conn_master(1)*nds  ,1)
!                
!                gap = xs - x11
!                
!                gt = gap(1,1)*tang(1,1) + gap(2,1)*tang(2,1)
!                        
!                xib = gt/lm
!                
!                call BSplineBasisAndDeriv2(ncp_m,p_m,xib,u_knot_master,Rmm,dRmm,ddRmm,kspan)
!            
!                Rm       = 0.0d0
!                dRmdxi   = 0.0d0
!                d2Rmdxi2 = 0.0d0
!                count = 0
!                do k=kspan-p_m,kspan-p_m+p_m
!                    count = count + 1
!                    Rm(k)       =   Rmm(count)
!                    dRmdxi(k)   =  dRmm(count)
!                    d2Rmdxi2(k) = ddRmm(count)
!                end do
!                
!                continue
!            end if !dlm == 0.0d0
!            
!            !Auxiliary vectors
!            vec = xs - xm
!            dvec = xs - dxm
!            
!            if(abs(dlm) .lt. 1.0d-7) goto 2
!            
!            !Residual
!            rsd = vec(1,1)*tang(1,1) + vec(2,1)*tang(2,1)
!            num = dvec(1,1)*dtang(1,1) + dvec(2,1)*dtang(2,1)
!            
!            !Update Parametric Coordinate of the Closest Point Projection 
!            if(rsd .gt. 0.0d0 .and. num .lt. 0.0d0) then !For a convex master surface
!                xib = xib - rsd/num
!            elseif(rsd .gt. 0.0d0 .and. num .gt. 0.0d0)then !For a concave master surface
!                xib = xib + rsd/num
!            elseif(rsd .lt. 0.0d0 .and. num .lt. 0.0d0)then
!                xib = xib - rsd/num
!            else
!                xib = xib + rsd/num
!            end if
!            
!            
!            if(xib .gt. 1.0d0) xib = 1.0d0
!            if(xib .lt. 0.0d0) xib = 0.0d0
!            
!            if(abs(rsd) .lt. 1.0d-8) goto 2
!            
!            continue
!        end do !Newton iterative cycle
        
        
        !------------------------------------------------------------------
        ! Secant Method to determine the final Closest Point Projection
        !------------------------------------------------------------------
        xib0 = 0.0d0
        xib = 0.9d0
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
                    xm(1,1)  = xm(1,1)  + Rm(k)*(    GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
                    xm(2,1)  = xm(2,1)  + Rm(k)*(    GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
                    continue
                end do
                
                !Tangent vector in the master segment and its derivative
                tang = 0.0d0
                do k=1,ncp_m
                    tang(1,1)  = tang(1,1)  + dRmdxi(k)*(  GCoords(conn_master(k),1) + ddisp(conn_master(k)*2-1,1))
                    tang(2,1)  = tang(2,1)  + dRmdxi(k)*(  GCoords(conn_master(k),2) + ddisp(conn_master(k)*2  ,1))
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
            
            if(abs(rsd - rsd0) .lt. 1.0d-12)then
                xib = xib - rsd*(xib - xib0)
                rsd = 1.0d-15
            else    
                xib = xib - rsd*(xib - xib0)/(rsd - rsd0)
            end if
           
            xib0 = xibi
            rsd0 = rsd
            
            if(xib .gt. 1.0d0) xib = 1.0d0
            if(xib .lt. 0.0d0) xib = 0.0d0
            
            if(abs(rsd) .lt. 1.0d-10) goto 2
            
            continue
            
        end do
        
        write(*,FMT=17)i
        17 format('Warning: Failed to determine final CPP for Greville point ',I3)
        
        2 continue
        
        !Compute normal vector from cross product {tang}x{0 0 -1}
        nor(1,1) = -tang(2,1)
        nor(2,1) =  tang(1,1)
        
        !Compute normal gap
        gn = vec(1,1)*nor(1,1)+vec(2,1)*nor(2,1)
        
        final_gap(i) = gn
        
        !dl(i)=lm
        continue
        
        !If contact has occured
        Kdelta = 0.0d0
        if(final_gap(i) .lt. -1.0d-10)then
            
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
            T0s = 0.0d0
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
            
            !Linearized variation of the gap
            if(nlgeom==.false.)then
                Kdelta = 0.0d0
            else
                Kdelta = 0.0d0
                Kdelta = c1*matmul(N0s,transpose( Ts)) + &
                &        c2*matmul( Ts,transpose(N0s)) + &
                &        c3*matmul(N0s,transpose(N0s)) + &
                &        c4*matmul( Ts,transpose( Ts))
                
                Kdelta = Kdelta*Lagrange(i)
                
            end if
            
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
            do k1=1,1
                do k2=1,1
                    Kc(k1*2-1,k2*2-1) = Kc(k1*2-1,k2*2-1) + Kdelta(k1*2-1,k2*2-1)
                    Kc(k1*2-1,k2*2  ) = Kc(k1*2-1,k2*2  ) + Kdelta(k1*2-1,k2*2  )
                    Kc(k1*2  ,k2*2-1) = Kc(k1*2  ,k2*2-1) + Kdelta(k1*2  ,k2*2-1)
                    Kc(k1*2  ,k2*2  ) = Kc(k1*2  ,k2*2  ) + Kdelta(k1*2  ,k2*2  )
                end do
                
                Kc(k1*2-1,3) = Ns(k1*2-1,1)
                Kc(k1*2  ,3) = Ns(k1*2  ,1)
                
                Kc(3,k1*2-1) = Ns(k1*2-1,1)
                Kc(3,k1*2  ) = Ns(k1*2  ,1)
                
                Fc(k1*nds-1,1) = Ns(k1*2-1,1)
                Fc(k1*nds  ,1) = Ns(k1*2  ,1)
            end do
            
            count1 = 0
            count2 = 0
            Kcc = 0.0d0
            Fcc = 0.0d0
            do k1=1,ncp_s
                do k2=1,ncp_s
                    count1 = (k1-1)*nds
                    do k3=1,nds
                        count1 = count1+1
                        count2 = (k2-1)*nds
                        do k4=1,nds
                            count2=count2 + 1
                            Kcc(count1,count2)= Kc(k3,k4)*Rs(k1)*Rs(k2)
                            continue
                        end do
                        
                         Kcc(ncp_s*nds+1,count1) = Kc(nds+1,k3)*Rs(k1)
                         Kcc(count1,ncp_s*nds+1) = Kc(k3,nds+1)*Rs(k1)
                            
                         Fcc(count1,1) = Fc(k3,1)*Rs(k1)
                         continue
                    end do
                end do
            end do
            
            
            do k1=1,ncp_s
                do k2=1,ncp_s
                    KLM(conn_slave(k1)*2-1,conn_slave(k2)*2-1) = KLM(conn_slave(k1)*2-1,conn_slave(k2)*2-1) + Kcc(k1*2-1,k2*2-1)
                    KLM(conn_slave(k1)*2-1,conn_slave(k2)*2  ) = KLM(conn_slave(k1)*2-1,conn_slave(k2)*2  ) + Kcc(k1*2-1,k2*2  )
                    KLM(conn_slave(k1)*2  ,conn_slave(k2)*2-1) = KLM(conn_slave(k1)*2  ,conn_slave(k2)*2-1) + Kcc(k1*2  ,k2*2-1)
                    KLM(conn_slave(k1)*2  ,conn_slave(k2)*2  ) = KLM(conn_slave(k1)*2  ,conn_slave(k2)*2  ) + Kcc(k1*2  ,k2*2  )
                end do
                
                KLM(conn_slave(k1)*2-1,tnodes*nds+i) =  Kcc(k1*2-1,ncp_s*nds+1)
                KLM(conn_slave(k1)*2  ,tnodes*nds+i) =  Kcc(k1*2  ,ncp_s*nds+1)
                
                KLM(tnodes*nds+i,conn_slave(k1)*2-1) =  Kcc(ncp_s*nds+1,k1*2-1)
                KLM(tnodes*nds+i,conn_slave(k1)*2  ) =  Kcc(ncp_s*nds+1,k1*2  )
                
                FLM(conn_slave(k1)*nds-1,i) = Fcc(k1*nds-1,1)
                FLM(conn_slave(k1)*nds  ,i) = Fcc(k1*nds  ,1)
                
            end do
            
            continue

            !------------------------------------------------------------------
            !Contact contribution to the master points
            !------------------------------------------------------------------
            do k1=2,ncp_m+1
                do k2=2,ncp_m+1
                    KLM(ctconn(k1)*2-1,ctconn(k2)*2-1) = KLM(ctconn(k1)*2-1,ctconn(k2)*2-1) + Kdelta(k1*2-1,k2*2-1)
                    KLM(ctconn(k1)*2-1,ctconn(k2)*2  ) = KLM(ctconn(k1)*2-1,ctconn(k2)*2  ) + Kdelta(k1*2-1,k2*2  )
                    KLM(ctconn(k1)*2  ,ctconn(k2)*2-1) = KLM(ctconn(k1)*2  ,ctconn(k2)*2-1) + Kdelta(k1*2  ,k2*2-1)
                    KLM(ctconn(k1)*2  ,ctconn(k2)*2  ) = KLM(ctconn(k1)*2  ,ctconn(k2)*2  ) + Kdelta(k1*2  ,k2*2  )
                end do
            end do
            
            do k1=2,ncp_m+1
                KLM(ctconn(k1)*2-1,tnodes*nds+i) = Ns(k1*2-1,1)
                KLM(ctconn(k1)*2  ,tnodes*nds+i) = Ns(k1*2  ,1)
                
                KLM(tnodes*nds+i,ctconn(k1)*2-1) = Ns(k1*2-1,1)
                KLM(tnodes*nds+i,ctconn(k1)*2  ) = Ns(k1*2  ,1)
            end do
            
            do k1=2,ncp_m+1
                FLM(ctconn(k1)*nds-1,i) = FLM(ctconn(k1)*nds-1,i) + Ns(k1*2-1,1)
                FLM(ctconn(k1)*nds  ,i) = FLM(ctconn(k1)*nds  ,i) + Ns(k1*2  ,1)
            end do
            
            GLM(tnodes*nds+i,1) = -1.0d0*initial_gap(i)
            
            continue
        
        else
            
            99 continue
            
            iLM = iLM + 1
            KLM(tnodes*nds+i,tnodes*nds+i) = 1.0d0
        
            
        end if !gn .lt. 0.0d0
        
        
        
        !xibs(i) = xib
        
    end do !Greville points cycle
    
    
    if(iLM .gt. 0)then
             
        !if(iter==1 .and. inc == 1)then
        allocate(KfeqLM(tnodes*nds-nbc+ngrv,tnodes*nds-nbc+ngrv))
        allocate(FeqLM(tnodes*nds-nbc+ngrv,1),dispLM(tnodes*nds-nbc+ngrv,1))
        allocate(disptempLM(tnodes*nds+ngrv))
        !end do
        
        !KCont = 0.0d0
        KCont(1:tnodes*nds,1:tnodes*nds) = KLM(1:tnodes*nds,1:tnodes*nds)
        KLM(1:tnodes*nds,1:tnodes*nds) = KLM(1:tnodes*nds,1:tnodes*nds) + Kf(1:tnodes*nds,1:tnodes*nds)
        
        KT(1:tnodes*nds,1:tnodes*nds) = KT(1:tnodes*nds,1:tnodes*nds) + KLM(1:tnodes*nds,1:tnodes*nds)
        
        GLM(1:tnodes*nds,1) = GLM(1:tnodes*nds,1) + Fext(1:tnodes*nds,1)
        
        disptempLM = 1.0d0
        disptempLM(1:tnodes*nds) = redcount(1:tnodes*nds)
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
        
!        count = 0
!        do k1=1,tnodes*nds+iLM-nbc
!            do k2=1,tnodes*nds+iLM-nbc
!                tst = KfeqLM(k1,k2)-KfeqLM(k2,k1)
!                if(abs(tst) .gt. 1.0d-14)then
!                    count = count + 1
!                    continue
!                end if
!            end do
!        end do
!        
!        continue
        

!        do k1=1,tnodes*nds-nBC+iLM
!            do k2=1,tnodes*nds-nBC+iLM
!                if(abs(KfeqLM(k1,k2)) .lt. 1.0d-10) KfeqLM(k1,k2) = 0.0d0
!            end do
!            if(abs(FEqLM(k1,1)) .lt. 1.0d-10) FEqLM(k1,1) = 0.0d0
!        end do
        
        temp1=tnodes*nds-nBC+ngrv
        dispLM = 0.0d0
        call Gauss (temp1,KfeqLM,FEqLM,dispLM)
        
        dDisp = dDisp - ddDisp
        
        ddDisp = 0.0d0
        j=1
        do i=1,tnodes*nds
            if(redcount(i)==1.0d0)then
                ddDisp(i,1)=dispLM(j,1)
                j=j+1
            end if
        end do
        
        do i=1,ngrv
            Lagrange(i) = Lagrange(i) + dispLM(tnodes*nds-nbc+i,1)
            if(Lagrange(i) .gt. 0.0d0) then
                !Lagrange(i) = 0.0d0
!                dispLM(tnodes*nds-nbc+i,1) = 0.0d0
             end if
!            if(dispLM(tnodes*nds-nbc+i,1) >= 0.0d0) then
!            !if(final_gap(i) .gt. 0.0d0) then
!                dispLM(tnodes*nds-nbc+i,1) = 0.0d0
!                Lagrange(i) = 0.0d0
!            end if
        end do
        
        count = 0
        do i=1,ngrv
                
            count = count + 1
            
            ctconn = 0
            !ctconn(1) = conn_slave(i)
            do j=1,ncp_m
                ctconn(j+1)=conn_master(j)
            end do
            
            do k1=1,ncp_s
                FintLM(conn_slave(k1)*nds-1,1) = FintLM(conn_slave(k1)*nds-1,1) + FLM(conn_slave(k1)*nds-1,i)*dispLM(tnodes*nds-nbc+i,1)
                FintLM(conn_slave(k1)*nds  ,1) = FintLM(conn_slave(k1)*nds  ,1) + FLM(conn_slave(k1)*nds  ,i)*dispLM(tnodes*nds-nbc+i,1)
            end do
              
            do k1=2,ncp_m+1
                FintLM(ctconn(k1)*nds-1,1) = FintLM(ctconn(k1)*nds-1,1) + FLM(ctconn(k1)*nds-1,i)*dispLM(tnodes*nds-nbc+i,1)
                FintLM(ctconn(k1)*nds  ,1) = FintLM(ctconn(k1)*nds  ,1) + FLM(ctconn(k1)*nds  ,i)*dispLM(tnodes*nds-nbc+i,1)
            end do
            
            continue
        end do
        
        continue
        
        dDisp = dDisp + ddDisp
        continue
        
        
        !------------------------------------------------------------------
        !Contact Stress
        !------------------------------------------------------------------
!        CF = 0.0d0
!        do i=1,ncp_s
!            do j=1, ngrv
!                !Compute B-Splines basis functions...
!                call BSplineBasisAndDeriv2(ncp_s,p_s,grev(j),u_knot_slave,Rss,dRss,ddRss,kspan)
!                !... and convert to NURBS basis, including vanishing terms 
!                call BSplinesToNURBS(ncp_s,p_s,grev(j),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
!                
!                CF(i) = CF(i) + Rs(i)*Lagrange(j)
!            end do
!        end do
!        
!        !CF = CF*(1.0d0*ngrv)
!        
!        sumdl = 0.0d0
!        bl = 0.0d0
!        br = 0.0d0
!        CF = 0.0d0
!        tag = .false.
!        do i=1,ngrv
!            
!            if(i==1)then
!                bl(i) = 0.0d0
!                br(i) = (grev(i) + grev(i+1))*(2.0d0/(ngrv*p_s))
!                if(p_s==1)then
!                    br(i) = (grev(i) + grev(i+1))/2.0d0
!                end if
!            elseif(i==ngrv)then
!                bl(i) = br(i-1)
!                br(i) = 1.0d0
!            else
!                bl(i) = br(i-1)
!                br(i) = grev(i) - br(i-1) + (grev(i) + grev(i+1))*(2.0d0/(ngrv*p_s))
!                if(p_s==1)then
!                    br(i) = grev(i) - br(i-1) + (grev(i) + grev(i+1))/2.0d0
!                end if
!            end if
!            
!            call BSplineBasisAndDeriv2(ncp_s,p_s,bl(i),u_knot_slave,Rss,dRss,ddRss,kspan)
!            call BSplinesToNURBS(      ncp_s,p_s,bl(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rl,dRsdxi,d2Rsdxi2)
!            
!            dx = 0.0d0
!            dy = 0.0d0
!            do j=1,ncp_s
!                dx = dx + dRsdxi(j)*(GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*2-1),1))
!                dy = dy + dRsdxi(j)*(GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*2  ),1))
!            end do
!            
!            suml = 0.0d0
!            suml = dsqrt(dx*dx+dy*dy)*bl(i)
!            
!            call BSplineBasisAndDeriv2(ncp_s,p_s,br(i),u_knot_slave,Rss,dRss,ddRss,kspan)
!            call BSplinesToNURBS(      ncp_s,p_s,br(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rr,dRsdxi,d2Rsdxi2)
!            
!            
!            dx = 0.0d0
!            dy = 0.0d0
!            do j=1,ncp_s
!                dx = dx + dRsdxi(j)*(GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*2-1),1))
!                dy = dy + dRsdxi(j)*(GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*2  ),1))
!            end do
!            
!            sumr = 0.0d0
!            sumr = dsqrt(dx*dx+dy*dy)*br(i)
!            
!            dl(i) = sumr - suml
!            
!            call BSplineBasisAndDeriv2(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan)
!            call BSplinesToNURBS(      ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
!
!            
!            !dl(i) =2.0d0*3.1416d0*10.0d0/4.0d0*(br(i)-bl(i)) 
!            
!            CF(i) = Lagrange(i)/dl(i)
!            
!            sumdl = sumdl + dl(i)
!            
!!            if(i == 1 .or. i == ngrv)then
!!                CF(i) = Lagrange(i)*(2.0d0*ngrv/p_s)*dl(i)
!!            else
!!                CF(i) = Lagrange(i)*(ngrv/2.0d0/p_s)*dl(i)
!!            end if
!!            
!!            if(i==1 .or. i==ngrv)then
!!                CF(i) = CF(i)*3.0d0
!!            end if
!            
!            
!            continue
!            
!            
!        end do
!        
!        continue






!        CF = 0.0d0
!        do i=1,ncp_s
!            do j=1, ngrv
!                !Compute B-Splines basis functions...
!                call BSplineBasisAndDeriv2(ncp_s,p_s,grev(j),u_knot_slave,Rss,dRss,ddRss,kspan)
!                !... and convert to NURBS basis, including vanishing terms 
!                call BSplinesToNURBS(ncp_s,p_s,grev(j),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
!                
!                CF(i) = CF(i) + Rs(i)*Lagrange(j)
!            end do
!        end do
!        
!        sumdl = 0.0d0
!        
!        bl = 0.0d0
!        br = 0.0d0
!        dl = 0.0d0
!        do i=1,ncp_s
!            bl(i) = u_knot_slave(i)
!            br(i) = u_knot_slave(i+p_s+1)
!            
!            br(i) = 0.0d0
!            
!!            do j=1,i
!!                do k=1,i+p_s+1
!!                    br(i) = br(i) + Rs(i)
!!                end do
!!            end do
!!            
!!            continue
!!            
!            sumr = 0.0d0
!            suml = 0.0d0
!            br(i) = 0.0d0
!            do k=1,i
!                !do j=i,i+p_s+1
!                    call BSplineBasisAndDeriv2(ncp_s,p_s,u_knot_slave(i),u_knot_slave,Rss,dRss,ddRss,kspan)
!                    call BSplinesToNURBS(      ncp_s,p_s,u_knot_slave(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
!                    
!                    suml = suml + Rs(i)
!                    
!                    call BSplineBasisAndDeriv2(ncp_s,p_s,u_knot_slave(i+p_s+1),u_knot_slave,Rss,dRss,ddRss,kspan)
!                    call BSplinesToNURBS(      ncp_s,p_s,u_knot_slave(i+p_s+1),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
!                    
!                    sumr = sumr + Rs(i)
!                    
!                    continue
!                    
!                !end do
!            end do
!            br(i) = abs(sumr-suml)
!            
!
!            call BSplineBasisAndDeriv2(ncp_s,p_s,bl(i),u_knot_slave,Rss,dRss,ddRss,kspan)
!            call BSplinesToNURBS(      ncp_s,p_s,bl(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rl,dRsdxi,d2Rsdxi2)
!            
!!            suml = 0.0d0
!!            do j=1,ncp_s
!!                dx = GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*2-1),1)
!!                dy = GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*2  ),1)
!!                suml = suml + Rl(j)*dsqrt(dx*dx+dy*dy)
!!            end do
!
!            suml = 0.0d0
!            dx = 0.0d0
!            dy = 0.0d0
!            do j=1,ncp_s
!                dx = dx + Rl(j)*(GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*2-1),1))
!                dy = dy + Rl(j)*(GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*2  ),1))
!            end do
!            
!            suml = dsqrt(dx*dx+dy*dy)
!            
!            call BSplineBasisAndDeriv2(ncp_s,p_s,br(i),u_knot_slave,Rss,dRss,ddRss,kspan)
!            call BSplinesToNURBS(      ncp_s,p_s,br(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rr,dRsdxi,d2Rsdxi2)
!            
!            sumr = 0.0d0
!            dx = 0.0d0
!            dy = 0.0d0
!            do j=1,ncp_s
!                dx = dx + Rr(j)*(GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*2-1),1))
!                dy = dy + Rr(j)*(GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*2  ),1))
!            end do
!            
!            sumr = dsqrt(dx*dx+dy*dy)
!            
!            
!            dl(i) = sumr-suml
!            
!            CF(i) = CF(i)/dl(i)
!            
!            sumdl = sumdl + dl(i)/2.0d0
!            
!            continue
!            
!        end do
!        
!        continue
        
!        bl = 0.0d0
!        br = 0.0d0
!        dl = 0.0d0
!        do i=1,ngrv-1
!        
!            call BSplineBasisAndDeriv2(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan)
!            call BSplinesToNURBS(      ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rl,dRsdxi,d2Rsdxi2)
!            
!            suml = 0.0d0
!            dx = 0.0d0
!            dy = 0.0d0
!            do j=1,ncp_s
!                dx = dx + Rl(j)*(GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*2-1),1))
!                dy = dy + Rl(j)*(GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*2  ),1))
!            end do
!            
!            suml = dsqrt(dx*dx+dy*dy)
!            
!            call BSplineBasisAndDeriv2(ncp_s,p_s,grev(i+1),u_knot_slave,Rss,dRss,ddRss,kspan)
!            call BSplinesToNURBS(      ncp_s,p_s,grev(i+1),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rr,dRsdxi,d2Rsdxi2)
!            
!            sumr = 0.0d0
!            dx = 0.0d0
!            dy = 0.0d0
!            do j=1,ncp_s
!                dx = dx + Rr(j)*(GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*2-1),1))
!                dy = dy + Rr(j)*(GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*2  ),1))
!            end do
!            
!            sumr = dsqrt(dx*dx+dy*dy)
!            
!            dl(i) = sumr - suml
!            
!            !CF(i) = CF(i)/dl(i)
!            
!            continue
!            
!        end do
!        
!        do i=1,ngrv-1
!            CF(i) = CF(i)/dl(i)
!        end do
!        
!        
!        continue
        
        
!        sumr = 0.0d0
!        do i=1,ncp_s
!            call BSplineBasisAndDeriv2(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan)
!            call BSplinesToNURBS(      ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
!            
!            dx = 0.0d0
!            dy = 0.0d0
!            do j=1,ncp_s
!                dx = dx + dRsdxi(j)*(GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*2-1),1))
!                dy = dy + dRsdxi(j)*(GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*2  ),1))
!            end do
!            
!            dl(i) = dsqrt(dx*dx+dy*dy)
!            
!            CF(i) = CF(i)/dl(i)
!            
!            sumr = sumr + dl(i)
!            
!            continue
!            
!        end do
!        
!        
!        continue
        
        
!        dx = 0.0d0
!        dy = 0.0d0
!        dl = 0.0d0
!        do i=1,ngrv
!            
!            !Compute B-Splines basis functions...
!            call BSplineBasisAndDeriv2(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan)
!            !... and convert to NURBS basis, including vanishing terms 
!            call BSplinesToNURBS(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
!            
!            
!            do j=1,ncp_s
!                dx(i) = dx(i) +  dRsdxi(j)*(GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*2-1),1))
!                dy(i) = dy(i) +  dRsdxi(j)*(GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*2  ),1))
!            end do
!            
!            dl(i) = dsqrt(dx(i)*dx(i)+dy(i)*dy(i))
!            
!            CF(i) = CF(i)/dl(i)
!            
!            continue
!            
!        end do
        
        continue
    
    end if
    
    continue

end subroutine