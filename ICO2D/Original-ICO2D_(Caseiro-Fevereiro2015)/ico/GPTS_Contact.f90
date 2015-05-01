subroutine GPTS_Contact(iter,inc,istep)
    
    use Mod_Variables
    implicit none
    
    integer(4),intent(IN)::iter,inc,istep
    
    integer(4),parameter::knewton = 50
    
    integer(4)::i,j,k,k1,k2,k3,loc_num
    integer(4)::icel,iGP,ni,kspan,count,ipos
    
    integer(4),dimension(ncp_s+ncp_m)::smconn
    
    real(8)::xi_s,det1,det2,det
    real(8)::sumtot,sumxi
    real(8)::xib,xib0,lm,rsd,rsd0,xibi
    real(8)::cfv,gn
    
    real(8)::en
    
    real(8),dimension(nds,1)::tang,dtang,nor,vec
    
    real(8),dimension(nds,1)::xs,xm,dpos
    real(8),dimension(SurfGP)::GP,wg
    real(8),dimension(p_s+1)::N
    real(8),dimension(p_s+1)::dNdxi
    real(8),dimension(p_s+1)::R,dRdxi
    real(8),dimension(2)::dxdxi
    
    real(8),dimension(p_s+1)::Rss,dRss,ddRss
    real(8),dimension(p_m+1)::Rmm,dRmm,ddRmm
    real(8),dimension(ncp_s)::Rs,dRsdxi,d2Rsdxi2
    real(8),dimension(ncp_m)::Rm,dRmdxi,d2Rmdxi2
    
    real(8),dimension((ncp_m+ncp_s)*nds,1)::Ns,N0s,Ts,T0s
    
    real(8)::a11,k11,m11
    real(8),dimension((ncp_m+ncp_s)*nds,1)::MatD, NBar
    
    logical::gaplog
    
    real(8),dimension(:,:),allocatable::KPenal, FPenal, KGeo
    
    en = 1.0d4
    
    allocate(KPenal((ncp_s+ncp_m)*nds,(ncp_s+ncp_m)*nds))
    KPenal = 0.0d0
    
    allocate(KGeo((ncp_s+ncp_m)*nds,(ncp_s+ncp_m)*nds))
    Kgeo = 0.0d0
    
    allocate(FPenal((ncp_s+ncp_m)*nds,1))
    FPenal = 0.0d0
    
    
    !Contact Element Cycle (Slave surface elements)
    icel = 0
    do i=1, ncp_s+p_s
        if(u_knot_slave(i) .ne. u_knot_slave(i+1)) then
        
            icel = icel + 1
        
            !Gauss points coordinates and weights
            call gauleg(SurfGP, GP, wg)
            
            !Gauss point cycle
            do iGP=1, SurfGP
            
                !Parametric coordinate on the slave surface
                ni = INN_s(IEN_s(icel,1),1)
                xi_s  = ((u_knot_slave(ni+1) - u_knot_slave(ni))*GP(iGP)  + (u_knot_slave(ni+1) + u_knot_slave(ni)))/2.0d0   
                
                !Determinant of the Jacobian for the parent-parametric transformation
                det1 = (U_knot_slave(ni+1)-U_knot_slave(ni))/2.0d0
                
                !B-Spline Basis for the slave element
                call BSplineBasisAndDeriv(ncp_s,p_s,xi_s,u_knot_slave,N,dNdxi)
                
                !NURBS Basis of the slave surface
                loc_num = 0
                R = 0.0d0
                dRdxi = 0.0d0
                sumtot = 0.0d0
                sumxi = 0.0d0
                do j=0,p_s
                    loc_num = loc_num + 1
                    R(loc_num) = N(p_s+1-j)*Weights_slv(IEN_s(icel,loc_num))
                    sumtot = sumtot + R(loc_num)
                    
                    dRdxi(loc_num) = dNdxi(p_s+1-j)*Weights_slv(IEN_s(icel,loc_num))
                    sumxi = sumxi + dRdxi(loc_num)
                end do
                
                do loc_num=1,(p+1)
                    dRdxi(loc_num) = (dRdxi(loc_num)*sumtot - R(loc_num)*sumxi )/(sumtot*sumtot)
                    R(loc_num) = R(loc_num)/sumtot
                end do
                
                loc_num = 0
                dxdxi = 0.0d0
                do j=0,p_s
                    loc_num = loc_num + 1
                    dxdxi(1) = dxdxi(1) + dRdxi(loc_num)*Points_slv(IEN_s(icel,loc_num),1)
                    dxdxi(2) = dxdxi(2) + dRdxi(loc_num)*Points_slv(IEN_s(icel,loc_num),2)
                end do
                
                !Determinant of the Jacobian for the parametric-global transformation
                det2 = dsqrt(dxdxi(1)*dxdxi(1)+dxdxi(2)*dxdxi(2))
                
                !Determinant of the Jacobian for the parent-global transformation
                det = det1*det2
                
                !Compute B-Splines basis functions...
                call BSplineBasisAndDeriv2(ncp_s,p_s,xi_s,u_knot_slave,Rss,dRss,ddRss,kspan)
        
                !... and convert to NURBS basis, including vanishing terms 
                call BSplinesToNURBS(ncp_s,p_s,xi_s,u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
                
                !Slave point in physical slave surface
                xs = 0.0d0
                do j=1,ncp_s
                    xs(1,1) = xs(1,1) + Rs(j)*(GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*nds-1),1))
                    xs(2,1) = xs(2,1) + Rs(j)*(GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*nds  ),1))
                end do
                
                !------------------------------------------------------------------
                ! Secant Method to determine the Closest Point Projection
                !------------------------------------------------------------------
                xib0 = 0.0d0
                xib = 0.123581d0
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
                        !tang = tang/lm
                        
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
                    !tang = tang/lm
                    
                    continue
                    
                    !Auxiliary vectors
                    vec = xs - xm
                    
                    rsd = vec(1,1)*tang(1,1) + vec(2,1)*tang(2,1)

                    xibi = xib
                    
                    xib = xib - rsd*(xib - xib0)/(rsd - rsd0)
                    
                    if(xib .gt. 1.0d0) xib = 1.0d0
                    if(xib .lt. 0.0d0) xib = 0.0d0
                    
                    xib0 = xibi
                    rsd0 = rsd
                    
                    if(abs(rsd) .lt. 1.0d-12) goto 1
                    
                    continue
                    
                end do
                
                write(*,FMT=16)iGP
                16 format('Warning: Failed to determine CPP for contact point ',I3)
                
                1 continue
                
                !Compute normal vector from cross product {tang}x{0 0 -1}
                nor(1,1) = -tang(2,1)/lm
                nor(2,1) =  tang(1,1)/lm
                
                ipos = SurfGP*(icel-1) + iGP
                
                !Compute initial normal gap
                initial_gap(ipos) = final_gap(ipos)
                final_gap(ipos) = vec(1,1)*nor(1,1) + vec(2,1)*nor(2,1)
                
                gaplog = ISNAN (final_gap(ipos))
                if(gaplog == .true.) final_gap(ipos) = 100.0d0
                
                if(final_gap(ipos) <= 0.0d0) then
                    PTS_active(ipos) = 1
                else
                    PTS_active(ipos) = 0
                end if
                
                if(final_gap(ipos) <= 0.0d0) then
                    Press(ipos) = en*final_gap(ipos)
                else
                    Press(ipos) = 0.0d0
                end if
                
                cfv = -1.0d0*Press(ipos) - MP_props(1,1)*final_gap(ipos)
                
                if(cfv .gt. 0.0d0)then
                    PTS_active(ipos) = 1
                else
                    PTS_active(ipos) = 0
                end if
                
                !If contact has occured
                if(PTS_active(ipos) == 1)then
                    
                    gn = final_gap(ipos)
                    
                    m11 = lm*lm
            
                    !Second derivative of the Closest Point Projection
                    dpos = 0.0d0
                    do k=1,ncp_m
                        dpos(1,1) = dpos(1,1) + d2Rmdxi2(k)*(GCoords(conn_master(k),1) + ddisp(conn_master(k)*nds-1,1))
                        dpos(2,1) = dpos(2,1) + d2Rmdxi2(k)*(GCoords(conn_master(k),2) + ddisp(conn_master(k)*nds  ,1))
                        continue
                    end do
                    
                    !Curvature of the boundary
                    k11 = dpos(1,1)*nor(1,1) + dpos(2,1)*nor(2,1)
                    
                    a11 = m11-gn*k11
                    
                    Ns = 0.0d0
                    Ts = 0.0d0
                    N0s = 0.0d0
                    do j=1,ncp_s
                        Ns(j*nds-1,1) = Rs(j)*nor(1,1)
                        Ns(j*nds  ,1) = Rs(j)*nor(2,1)
                        
                        Ts(j*nds-1,1) = Rs(j)*tang(1,1)
                        Ts(j*nds  ,1) = Rs(j)*tang(2,1)
                    end do
                    
                    count = 0
                    do j=ncp_s+1,ncp_s+ncp_m
                        
                        count = count + 1
                        
                        Ns(j*nds-1,1) = -1.0d0*Rm(count)*nor(1,1)
                        Ns(j*nds  ,1) = -1.0d0*Rm(count)*nor(2,1)
                        
                        N0s(j*nds-1,1) = -1.0d0*dRmdxi(count)*nor(1,1)
                        N0s(j*nds  ,1) = -1.0d0*dRmdxi(count)*nor(2,1)
                        
                        Ts(j*nds-1,1) = -1.0d0*Rm(count)*tang(1,1)
                        Ts(j*nds  ,1) = -1.0d0*Rm(count)*tang(2,1)
                        continue
                    end do
                    
                    !Geometric contribution ---------------------------------------------------------------------------------------------------------------------------------------------
                    matD = (Ts-gn*N0s)/a11
                    Nbar = N0s - matD*k11
                    
                    Kgeo = Kgeo + ( gn/m11*matmul(Nbar,transpose(Nbar)) + matmul(matD,transpose(N0s)) + matmul(N0s,transpose(matD)) - matmul(matD,transpose(matD))*k11 )*det*wg(iGP)*gn
                    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    
                    KPenal = KPenal + matmul(Ns,transpose(Ns))*det*wg(iGP)
                    
                    FPenal = FPenal + en*gn*Ns*det*wg(iGP)
                    
                    continue
                        
                end if
                
                continue
            
            end do !Gauss point cycle
            
            
            
           continue 
            
            
        end if
    end do !Surface element cycle
    
    KPenal = KPenal*en + Kgeo*en
    
    
    !------------------------------------------------------------------------------        
    !Assembly contributions to global stiffness matrix and Contact Forces
    !------------------------------------------------------------------------------
    
    !Deformable-Rigid contact
!    Fct = 0.0d0
!    do k1=1,ncp_s
!        do k2=1,ncp_s
!            Kf(conn_slave(k1)*2-1,conn_slave(k2)*2-1) = Kf(conn_slave(k1)*2-1,conn_slave(k2)*2-1) + KPenal(k1*2-1,k2*2-1)
!            Kf(conn_slave(k1)*2-1,conn_slave(k2)*2  ) = Kf(conn_slave(k1)*2-1,conn_slave(k2)*2  ) + KPenal(k1*2-1,k2*2  )
!            Kf(conn_slave(k1)*2  ,conn_slave(k2)*2-1) = Kf(conn_slave(k1)*2  ,conn_slave(k2)*2-1) + KPenal(k1*2  ,k2*2-1)
!            Kf(conn_slave(k1)*2  ,conn_slave(k2)*2  ) = Kf(conn_slave(k1)*2  ,conn_slave(k2)*2  ) + KPenal(k1*2  ,k2*2  )
!        end do
!        
!        Fct(conn_slave(k1)*nds-1,1) = Fct(conn_slave(k1)*nds-1,1) + FPenal(k1*nds-1,1)
!        Fct(conn_slave(k1)*nds  ,1) = Fct(conn_slave(k1)*nds  ,1) + FPenal(k1*nds  ,1)
!    end do
    
    !Deformable-Deformable contact
    do j=1,ncp_s
        smconn(j)=conn_slave(j)
    end do
    
    count = 0
    do j=ncp_s+1,ncp_m+ncp_s
        count = count + 1
        smconn(j)=conn_master(count)
        continue
    end do
    
    Fct = 0.0d0
    do k1=1,ncp_s+ncp_m
        do k2=1,ncp_s+ncp_m
            Kf(smconn(k1)*2-1,smconn(k2)*2-1) = Kf(smconn(k1)*2-1,smconn(k2)*2-1) + KPenal(k1*2-1,k2*2-1)
            Kf(smconn(k1)*2-1,smconn(k2)*2  ) = Kf(smconn(k1)*2-1,smconn(k2)*2  ) + KPenal(k1*2-1,k2*2  )
            Kf(smconn(k1)*2  ,smconn(k2)*2-1) = Kf(smconn(k1)*2  ,smconn(k2)*2-1) + KPenal(k1*2  ,k2*2-1)
            Kf(smconn(k1)*2  ,smconn(k2)*2  ) = Kf(smconn(k1)*2  ,smconn(k2)*2  ) + KPenal(k1*2  ,k2*2  )
        end do
        
        Fct(smconn(k1)*nds-1,1) = Fct(smconn(k1)*nds-1,1) + FPenal(k1*nds-1,1)
        Fct(smconn(k1)*nds  ,1) = Fct(smconn(k1)*nds  ,1) + FPenal(k1*nds  ,1)
    end do
    
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
                FextEq(k1,1)=Fext(i,1) - Fct(i,1)
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
    
    continue

end subroutine