!------------------------------------------------------------------------------------------------------------------------
!
! Point-to-Segment contact algorithm module
!
!------------------------------------------------------------------------------------------------------------------------

subroutine PTS_ComputeGap()

    use mod_variables
    use mod_MBvariables
    implicit none
    integer(4),parameter::knewton=50
    
    integer(4)::i,j,k,l,k1,k2,k3,k4,k5,kn,igrv,ior,ipair,itot
    integer(4)::ispan,jspan,count,count2,itemp
    real(8)::sumt,sumxi,sumeta,sumxixi,sumetaeta,sumxieta,sumetaxi
    
    real(8),dimension(nds,1)::xs,xm,txi,teta,xab
    real(8)::xib,etab,xib0,etab0,xibi,etabi
    real(8)::rsd_xi,rsd_eta,rsd0_xi,rsd0_eta,rsd
    real(8)::gn,cfv
    real(8)::nor_xi,nor_eta, norm
    real(8),dimension(nds,1)::vec,nor
    real(8),dimension(nds-1,nds-1)::Mab,MabInv,Kab,Aab,AabInv
    real(8),dimension(nds-1)::temp1
    
    real(8),dimension(1,1)::conv1,conv2
    real(8),dimension(2,2)::met,ddF,aux
    real(8),dimension(2,1)::dxm,new,dnew,dF
    real(8),dimension(3,4)::ddxm

    !------------------------------------------------------------------------------
    ! Contact pair cycle
    !------------------------------------------------------------------------------
    igrv = 0
    do ipair=1,npair
        
        call MB_Alloc(ipair) 
            
        !------------------------------------------------------------------------------
        ! Slave points cycle
        !------------------------------------------------------------------------------
        do j=1,ngrv_eta
            do i=1,ngrv_xi 
            
                igrv = igrv + 1
                
                !------------------------------------------------------------------------------
                ! Slave NURBS surface basis, including vanishing terms
                !------------------------------------------------------------------------------
                call SurfaceBasis(ncpx_slv,p_slv,grv_xi(i),u_knot_slv, &
                &                 ncpy_slv,q_slv,grv_eta(j),v_knot_slv, &
                &                 GCoords, conn_slv, tnodes,nds, &
                &                 Rslv,dRslv,ddRslv)
                
                !------------------------------------------------------------------------------
                ! Slave point in the physical space
                !------------------------------------------------------------------------------
                xs = 0.0d0
                count = 0
                do k1=1,ncpx_slv*ncpy_slv
                    count = count + 1
                    do k3=1,nds
                        xs(k3,1) = xs(k3,1) + Rslv(count)*GCoords(conn_slv(count),k3)
                    end do
                end do
                
                !------------------------------------------------------------------------------
                ! Closest Point Projection (CPP) for direction xi
                !------------------------------------------------------------------------------
                xib = 0.5d0
                etab = 0.5d0
                new = 0.5d0
                do kn=1,knewton
                    
                    !------------------------------------------------------------------------------
                    ! Master NURBS surface basis, including vanishing terms
                    !------------------------------------------------------------------------------
                    call SurfaceBasis(ncpx_mst,p_mst,xib,u_knot_mst, &
                    &                 ncpy_mst,q_mst,etab,v_knot_mst, &
                    &                 GCoords, conn_mst, tnodes,nds, &
                    &                 Rmst,dRmst,ddRmst)
                    
                    !------------------------------------------------------------------------------
                    ! CPP in the physical space
                    !------------------------------------------------------------------------------
                    xm = 0.0d0
                    count = 0
                    do k1=1,ncpx_mst*ncpy_mst
                        count = count + 1
                        do k3=1,nds
                            xm(k3,1) = xm(k3,1) + Rmst(count)*GCoords(conn_mst(count),k3)
                        end do
                    end do
                    
                    !------------------------------------------------------------------------------
                    ! Tangent vectors in the CPP
                    !------------------------------------------------------------------------------
                    txi = 0.0d0
                    teta = 0.0d0
                    count = 0
                    do k2=1,ncpx_mst*ncpy_mst
                        count = count + 1
                        do k3=1,nds
                             txi(k3,1) =  txi(k3,1) + dRmst(count,1)*GCoords(conn_mst(count),k3)
                            teta(k3,1) = teta(k3,1) + dRmst(count,2)*GCoords(conn_mst(count),k3)
                        end do
                    end do

                    !------------------------------------------------------------------------------
                    ! Hessian
                    !------------------------------------------------------------------------------
                    ddxm = 0.0d0
                    count = 0
                    do k1=1,ncpx_mst*ncpy_mst
                        count = count + 1
                        do k3=1,nds
                            ddxm(k3,1) = ddxm(k3,1) + ddRmst(count,1)*GCoords(conn_mst(count),k3)
                            ddxm(k3,2) = ddxm(k3,2) + ddRmst(count,2)*GCoords(conn_mst(count),k3)
                            ddxm(k3,3) = ddxm(k3,3) + ddRmst(count,3)*GCoords(conn_mst(count),k3)
                            ddxm(k3,4) = ddxm(k3,4) + ddRmst(count,4)*GCoords(conn_mst(count),k3)
                        end do
                    end do
                    
                    !------------------------------------------------------------------------------
                    ! Metric tensor
                    !------------------------------------------------------------------------------
                    Mab = 0.0d0
                    do k1=1,nds
                        Mab(1,1) = Mab(1,1) +  txi(k1,1)* txi(k1,1)
                        Mab(2,1) = Mab(2,1) + teta(k1,1)* txi(k1,1)
                        Mab(1,2) = Mab(1,2) +  txi(k1,1)*teta(k1,1)
                        Mab(2,2) = Mab(2,2) + teta(k1,1)*teta(k1,1)
                    end do
                    
                    !------------------------------------------------------------------------------
                    ! Auxilliary vector
                    !------------------------------------------------------------------------------
                    vec = 0.0d0
                    vec = xs - xm
                    
                    !------------------------------------------------------------------------------
                    ! First derivative
                    !------------------------------------------------------------------------------
                    dF = 0.0d0
                    do k1=1,nds
                        dF(1,1) = dF(1,1) +  txi(k1,1)* vec(k1,1)
                        dF(2,1) = dF(2,1) + teta(k1,1)* vec(k1,1)
                    end do
                    
                    dF = -2.0d0*dF
                    
                    !------------------------------------------------------------------------------
                    ! Second derivative
                    !------------------------------------------------------------------------------
                    aux = 0.0d0
                    do k1=1,nds
                        aux(1,1) = aux(1,1) + ddxm(k1,1)* vec(k1,1)
                        aux(1,2) = aux(1,2) + ddxm(k1,2)* vec(k1,1)
                        aux(2,1) = aux(2,1) + ddxm(k1,3)* vec(k1,1)
                        aux(2,2) = aux(2,2) + ddxm(k1,4)* vec(k1,1)
                    end do
                    
                    ddF = 0.0d0
                    ddF = 2.0d0*(Mab-aux) 
                    
                    !------------------------------------------------------------------------------
                    ! Inverse
                    !------------------------------------------------------------------------------
                    aux = ddF
                    
                    ddF = 0.0d0
                    ddF(1,1) = aux(2,2)
                    ddF(1,2) = -1.0d0*aux(1,2)
                    ddF(2,1) = -1.0d0*aux(2,1)
                    ddF(2,2) = aux(1,1)
                    
                    ddF = 1.0d0/(aux(1,1)*aux(2,2)-aux(1,2)*aux(2,1))*ddF
                    
                    !------------------------------------------------------------------------------
                    ! New point and parametric coordinates
                    !------------------------------------------------------------------------------
                    dnew = matmul(ddF,dF)
                    
                    new = new - dnew
                    
                    if(new(1,1) .gt. 1.0d0) new(1,1) = 1.0d0
                    if(new(2,1) .gt. 1.0d0) new(2,1) = 1.0d0
                    if(new(1,1) .lt. 0.0d0) new(1,1) = 0.0d0
                    if(new(2,1) .lt. 0.0d0) new(2,1) = 0.0d0
                    
                    xib = new(1,1)
                    etab = new(2,1)
                    
                    !------------------------------------------------------------------------------
                    ! Check if no change of the variables values
                    !------------------------------------------------------------------------------
                    conv1 = matmul(transpose(dnew),dnew)
                    
                    !------------------------------------------------------------------------------
                    ! Check if convergence was obtained
                    !------------------------------------------------------------------------------
                    conv2 = matmul(transpose(dF),dF)
     
                    !------------------------------------------------------------------------------
                    ! Exit conditions
                    !------------------------------------------------------------------------------
                    if(abs(conv1(1,1)) .lt. 1.0d-12 .or. abs(conv2(1,1)) .lt. 1.0d-12) goto 1
                    
                    continue
                
                end do
                
                PTS_active(igrv) = 0
                goto 10
                
                1 continue
                
                continue

                !------------------------------------------------------------------------------
                ! Cross product to determine normal to the master surface
                !------------------------------------------------------------------------------
                nor_xi = 0.0d0
                nor_eta = 0.0d0
                do k1=1,nds
                    nor_xi=nor_xi + txi(k1,1)*txi(k1,1)
                    nor_eta=nor_eta + teta(k1,1)*teta(k1,1)
                end do
                
                nor_xi = dsqrt(nor_xi)
                nor_eta = dsqrt(nor_eta)
                
                call CrossProduct(txi,teta,nor)
                
                norm = 0.0d0
                do k1=1,nds
                    norm=norm + nor(k1,1)*nor(k1,1)
                end do
                
                norm = dsqrt(norm)
                
                nor = nor/norm
                
                !------------------------------------------------------------------------------
                ! Cumpute normal gap
                !------------------------------------------------------------------------------
                gn = 0.0d0
                do k1=1,nds
                    gn = gn +  (xs(k1,1) - xm(k1,1))*nor(k1,1)
                end do
                
                gap(igrv) = gn
                
                continue
                
                !------------------------------------------------------------------------------
                ! Check if slave point is active
                !------------------------------------------------------------------------------
                if(gap(igrv) .lt. 0.0d0)then
                    PTS_Active(igrv) = 1
                end if
                
                cfv = -1.0d0*Lagrange(igrv) - MP_props(islave,1)*gap(igrv)
                
                if(cfv .gt. 0.0d0)then
                    PTS_Active(igrv) = 1
                else
                    PTS_Active(igrv) = 0
                end if
                
                10 continue
                
                !------------------------------------------------------------------------------
                ! If contact is active
                !------------------------------------------------------------------------------
                if(PTS_Active(igrv) == 1)then
                
                    !------------------------------------------------------------------------------
                    ! Metric tensor
                    !------------------------------------------------------------------------------
                    Mab = 0.0d0
                    do k1=1,nds
                        Mab(1,1) = Mab(1,1) +  txi(k1,1)* txi(k1,1)
                        Mab(2,1) = Mab(2,1) + teta(k1,1)* txi(k1,1)
                        Mab(1,2) = Mab(1,2) +  txi(k1,1)*teta(k1,1)
                        Mab(2,2) = Mab(2,2) + teta(k1,1)*teta(k1,1)
                    end do
                    
                    MabInv = 0.0d0
                    temp1 = 0.0d0
                    MabInv = Mab
                    call gaussj(MabInv,nds-1,temp1,ior)
                    
                    !------------------------------------------------------------------------------
                    ! Curvature tensor
                    !------------------------------------------------------------------------------
                    Kab = 0.0d0
                    count2 = 0
                    do k4=1,nds-1
                        do k5=1,nds-1
                            
                            count2 = count2 + 1
                            
                            xab = 0.0d0
                            count = 0
                            do k1=1,ncpx_mst*ncpy_mst
                                count = count + 1
                                do k3=1,nds
                                    xab(k3,1) = xab(k3,1) + ddRmst(count,count2)*GCoords(conn_mst(count),k3)
                                end do
                            end do
                            
                            do k1=1,nds
                                Kab(k4,k5) = Kab(k4,k5) + xab(k1,1)*nor(k1,1)
                            end do
                        end do
                    end do
                    
                    !------------------------------------------------------------------------------
                    ! Auxiliary tensor
                    !------------------------------------------------------------------------------
                    Aab = 0.0d0
                    do k1=1, nds-1
                        do k2=1,nds-1
                            Aab(k1,k2) = Mab(k1,k2) - gn*Kab(k1,k2)
                        end do
                    end do
                    
                    AabInv = 0.0d0
                    temp1 = 0.0d0
                    AabInv = Aab
                    call gaussj(AabInv,nds-1,temp1,ior)
                    
                    !------------------------------------------------------------------------------
                    ! Auxiliary arrays
                    !------------------------------------------------------------------------------
                    N = 0.0d0
                    Na = 0.0d0
                    Nb = 0.0d0
                    Ta = 0.0d0
                    Tb = 0.0d0
                    do k1=1,ncpx_slv*ncpy_slv
                        N(k1*nds-2,1) = Rslv(k1)*nor(1,1)
                        N(k1*nds-1,1) = Rslv(k1)*nor(2,1)
                        N(k1*nds  ,1) = Rslv(k1)*nor(3,1)
                        
                        Ta(k1*nds-2,1) = Rslv(k1)*txi(1,1)
                        Ta(k1*nds-1,1) = Rslv(k1)*txi(2,1)
                        Ta(k1*nds  ,1) = Rslv(k1)*txi(3,1)
                        
                        Tb(k1*nds-2,1) = Rslv(k1)*teta(1,1)
                        Tb(k1*nds-1,1) = Rslv(k1)*teta(2,1)
                        Tb(k1*nds  ,1) = Rslv(k1)*teta(3,1)
                    end do
                    
                    count = 0
                    do k1=ncpx_slv*ncpy_slv+1, ncpx_slv*ncpy_slv + ncpx_mst*ncpy_mst
                        
                        count = count + 1
                        
                        N(k1*nds-2,1) = -1.0d0*Rmst(count)*nor(1,1)
                        N(k1*nds-1,1) = -1.0d0*Rmst(count)*nor(2,1)
                        N(k1*nds  ,1) = -1.0d0*Rmst(count)*nor(3,1)
                        
                        Na(k1*nds-2,1) = -1.0d0*dRmst(count,1)*nor(1,1)
                        Na(k1*nds-1,1) = -1.0d0*dRmst(count,1)*nor(2,1)
                        Na(k1*nds  ,1) = -1.0d0*dRmst(count,1)*nor(3,1)
                        
                        Nb(k1*nds-2,1) = -1.0d0*dRmst(count,2)*nor(1,1)
                        Nb(k1*nds-1,1) = -1.0d0*dRmst(count,2)*nor(2,1)
                        Nb(k1*nds  ,1) = -1.0d0*dRmst(count,2)*nor(3,1)
                        
                        Ta(k1*nds-2,1) = -1.0d0*Rmst(count)*txi(1,1)
                        Ta(k1*nds-1,1) = -1.0d0*Rmst(count)*txi(2,1)
                        Ta(k1*nds  ,1) = -1.0d0*Rmst(count)*txi(3,1)
                        
                        Tb(k1*nds-2,1) = -1.0d0*Rmst(count)*teta(1,1)
                        Tb(k1*nds-1,1) = -1.0d0*Rmst(count)*teta(2,1)
                        Tb(k1*nds  ,1) = -1.0d0*Rmst(count)*teta(3,1)
                    end do
                    
                    Nhat = 0.0d0
                    That = 0.0d0
                    do k1=1,(ncpx_slv*ncpy_slv + ncpx_mst*ncpy_mst)*nds
                        Nhat(k1,1) = Na(k1,1)
                        Nhat(k1,2) = Nb(k1,1)
                        
                        That(k1,1) = Ta(k1,1)
                        That(k1,2) = Tb(k1,1)
                    end do
                    
                    D = 0.0d0
                    D = matmul((That-gn*Nhat),AabInv)
                    
                    Nbar = 0.0d0
                    Nbar = Nhat - matmul(D,Kab)
                    
                    !------------------------------------------------------------------------------
                    ! Geometric stiffness
                    !------------------------------------------------------------------------------
                    Kgeo = 0.0d0
                    
                    if(Lagrange(igrv) .gt. 0.0d0)then
                        KG1 = gn*matmul(matmul(Nbar,MabInv),transpose(Nbar))
                        KG2 = matmul(D,transpose(Nhat))
                        KG3 = matmul(Nhat,transpose(D))
                        KG4 = matmul(matmul(D,Kab),transpose(D))
                        Kgeo = (KG1 + KG2 + KG3 -KG4)*Lagrange(igrv)
                    end if
                    !Kgeo = 0.0d0
                    
    !                Kgeo = gn*matmul(matmul(Nbar,MabInv),transpose(Nbar)) + &
    !                &      matmul(D,transpose(Nhat)) + &
    !                &      matmul(Nhat,transpose(D)) - &
    !                &      matmul(matmul(D,Kab),transpose(D))
    !                Kgeo = Kgeo*Lagrange(igrv)
                    
                    !------------------------------------------------------------------------------
                    ! Assemble contact contributions to stiffness matrix
                    !------------------------------------------------------------------------------
                    If(MRigid == .true.)then
                        do k1=1,ncpx_slv*ncpy_slv
                        
                            KLM(conn_slv(k1)*nds-2,tnodes*nds+igrv) = KLM(conn_slv(k1)*nds-2,tnodes*nds+igrv) + N(k1*nds-2,1)
                            KLM(conn_slv(k1)*nds-1,tnodes*nds+igrv) = KLM(conn_slv(k1)*nds-1,tnodes*nds+igrv) + N(k1*nds-1,1)
                            KLM(conn_slv(k1)*nds  ,tnodes*nds+igrv) = KLM(conn_slv(k1)*nds  ,tnodes*nds+igrv) + N(k1*nds  ,1)
                            
                            KLM(tnodes*nds+igrv,conn_slv(k1)*nds-2) = KLM(tnodes*nds+igrv,conn_slv(k1)*nds-2) + N(k1*nds-2,1)
                            KLM(tnodes*nds+igrv,conn_slv(k1)*nds-1) = KLM(tnodes*nds+igrv,conn_slv(k1)*nds-1) + N(k1*nds-1,1)
                            KLM(tnodes*nds+igrv,conn_slv(k1)*nds  ) = KLM(tnodes*nds+igrv,conn_slv(k1)*nds  ) + N(k1*nds  ,1)
                            
                            continue
                        end do
                    else
                        do k1=1,ncpx_slv*ncpy_slv + ncpx_mst*ncpy_mst
                            do k2=1,ncpx_slv*ncpy_slv + ncpx_mst*ncpy_mst
                                KLM(PTS_conn(k1)*nds-2,PTS_conn(k2)*nds-2) = &
                                    & KLM(PTS_conn(k1)*nds-2,PTS_conn(k2)*nds-2) + Kgeo(k1*nds-2,k2*nds-2)
                                KLM(PTS_conn(k1)*nds-2,PTS_conn(k2)*nds-1) = &
                                    & KLM(PTS_conn(k1)*nds-2,PTS_conn(k2)*nds-1) + Kgeo(k1*nds-2,k2*nds-1)
                                KLM(PTS_conn(k1)*nds-2,PTS_conn(k2)*nds  ) = &
                                    & KLM(PTS_conn(k1)*nds-2,PTS_conn(k2)*nds  ) + Kgeo(k1*nds-2,k2*nds  )
                                
                                KLM(PTS_conn(k1)*nds-1,PTS_conn(k2)*nds-2) = &
                                    & KLM(PTS_conn(k1)*nds-1,PTS_conn(k2)*nds-2) + Kgeo(k1*nds-1,k2*nds-2)
                                KLM(PTS_conn(k1)*nds-1,PTS_conn(k2)*nds-1) = &
                                    & KLM(PTS_conn(k1)*nds-1,PTS_conn(k2)*nds-1) + Kgeo(k1*nds-1,k2*nds-1)
                                KLM(PTS_conn(k1)*nds-1,PTS_conn(k2)*nds  ) = &
                                    & KLM(PTS_conn(k1)*nds-1,PTS_conn(k2)*nds  ) + Kgeo(k1*nds-1,k2*nds  )
                                
                                KLM(PTS_conn(k1)*nds  ,PTS_conn(k2)*nds-2) = &
                                    & KLM(PTS_conn(k1)*nds  ,PTS_conn(k2)*nds-2) + Kgeo(k1*nds  ,k2*nds-2)
                                KLM(PTS_conn(k1)*nds  ,PTS_conn(k2)*nds-1) = &
                                    & KLM(PTS_conn(k1)*nds  ,PTS_conn(k2)*nds-1) + Kgeo(k1*nds  ,k2*nds-1)
                                KLM(PTS_conn(k1)*nds  ,PTS_conn(k2)*nds  ) = &
                                    & KLM(PTS_conn(k1)*nds  ,PTS_conn(k2)*nds  ) + Kgeo(k1*nds  ,k2*nds  )
                            end do
                            
                            KLM(PTS_conn(k1)*nds-2,tnodes*nds+igrv) = &
                                    & KLM(PTS_conn(k1)*nds-2,tnodes*nds+igrv) + N(k1*nds-2,1)
                            KLM(PTS_conn(k1)*nds-1,tnodes*nds+igrv) = &
                                    & KLM(PTS_conn(k1)*nds-1,tnodes*nds+igrv) + N(k1*nds-1,1)
                            KLM(PTS_conn(k1)*nds  ,tnodes*nds+igrv) = &
                                    & KLM(PTS_conn(k1)*nds  ,tnodes*nds+igrv) + N(k1*nds  ,1)
                            
                            KLM(tnodes*nds+igrv,PTS_conn(k1)*nds-2) = &
                                    & KLM(tnodes*nds+igrv,PTS_conn(k1)*nds-2) + N(k1*nds-2,1)
                            KLM(tnodes*nds+igrv,PTS_conn(k1)*nds-1) = &
                                    & KLM(tnodes*nds+igrv,PTS_conn(k1)*nds-1) + N(k1*nds-1,1)
                            KLM(tnodes*nds+igrv,PTS_conn(k1)*nds  ) = &
                                    & KLM(tnodes*nds+igrv,PTS_conn(k1)*nds  ) + N(k1*nds  ,1)
                        end do
                    end if
                   
                    !------------------------------------------------------------------------------
                    ! Build contact forces auxiliary array
                    !------------------------------------------------------------------------------
                    If(MRigid == .true.)then
                        do k1=1,ncpx_slv*ncpy_slv     
                            FLM(conn_slv(k1)*nds-2,igrv) = FLM(conn_slv(k1)*nds-2,igrv) + N(k1*nds-2,1)
                            FLM(conn_slv(k1)*nds-1,igrv) = FLM(conn_slv(k1)*nds-1,igrv) + N(k1*nds-1,1)
                            FLM(conn_slv(k1)*nds  ,igrv) = FLM(conn_slv(k1)*nds  ,igrv) + N(k1*nds  ,1)
                        end do
                    else
                        do k1=1,ncpx_slv*ncpy_slv + ncpx_mst*ncpy_mst
                            FLM(PTS_conn(k1)*nds-2,igrv) = FLM(PTS_conn(k1)*nds-2,igrv) + N(k1*nds-2,1)
                            FLM(PTS_conn(k1)*nds-1,igrv) = FLM(PTS_conn(k1)*nds-1,igrv) + N(k1*nds-1,1)
                            FLM(PTS_conn(k1)*nds  ,igrv) = FLM(PTS_conn(k1)*nds  ,igrv) + N(k1*nds  ,1)
                        end do
                    end if
                    !------------------------------------------------------------------------------
                    ! Right-hand side vector containing gap
                    !------------------------------------------------------------------------------
                    GLM(tnodes*nds+igrv,1) = -1.0d0*gap(igrv)
                    
                    continue
                else
                    
                    !------------------------------------------------------------------------------
                    ! If not in contact put 1 in the diagonal of the contact stiffness
                    !------------------------------------------------------------------------------
                    KLM(tnodes*nds+igrv,:) = 0.0d0
                    KLM(:,tnodes*nds+igrv) = 0.0d0
                    KLM(tnodes*nds+igrv,tnodes*nds+igrv) = 1.0d0
                    
                end if !contact active
                
                continue
                
            end do !ngrv_eta
        end do !ngrv_xi
        
        call MB_DeAlloc(ipair)
        
    end do !ipair
    
    
    !------------------------------------------------------------------------------
    ! Global stiffness matrix
    !------------------------------------------------------------------------------
    Ktot(1:tnodes*nds,1:tnodes*nds) = Ktot(1:tnodes*nds,1:tnodes*nds) + KLM(1:tnodes*nds,1:tnodes*nds)
    
    KLM(1:tnodes*nds,1:tnodes*nds) = KLM(1:tnodes*nds,1:tnodes*nds) + Kf(1:tnodes*nds,1:tnodes*nds)
    Kf(1:tnodes*nds,1:tnodes*nds)  = KLM(1:tnodes*nds,1:tnodes*nds)
    
    !------------------------------------------------------------------------------
    ! Global force vector
    !------------------------------------------------------------------------------
    GLM(1:tnodes*nds,1) = GLM(1:tnodes*nds,1) + Fext(1:tnodes*nds,1)
    
    if(npair .gt. 1)then
        itot = 0   
        do i=1,npair
            itot = itot + MB_ncpx_slv(i)*MB_ncpy_slv(i)
        end do
    else
        itot = ncpx_slv*ncpy_slv
    end if
    
    !------------------------------------------------------------------------------
    ! Assemble equivalent reduced matrices
    !------------------------------------------------------------------------------
    disptempLM = 1
    disptempLM(1:tnodes*nds) = redcount(1:tnodes*nds)
    KeqLM = 0.0d0
    FeqLM = 0.0d0
    k1 = 1
    k2 = 1
    do i=1,tnodes*nds + itot 
        do j=1,tnodes*nds + itot
            if (disptempLM(i)/= 0 .and. disptempLM(j)/= 0)then
                KeqLM(k1,k2)=KLM(i,j)
                FeqLM(k1,1)=GLM(i,1)                
                k2=k2+1
                if(k2==tnodes*nds - nbc + 1 + itot)then
                    k1=k1+1
                    k2=1
                end if
            end if
        end do
    end do
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    open(unit=1,file='KLM.txt')
!    write(1,*)tnodes*nds - nbc + itot
!    do k1=1,tnodes*nds - nbc + itot
!        do k2=1,tnodes*nds - nbc + itot
!            write(1,*)KeqLM(k1,k2)
!        end do
!    end do
!    close(1)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    !------------------------------------------------------------------------------
    ! Solve system of equations
    !------------------------------------------------------------------------------
    itemp = tnodes*nds - nbc + itot
    call Gauss (itemp,KeqLM,FeqLM,dispeqLM)
    
    !------------------------------------------------------------------------------
    ! Assemble global displacement increment increment
    !------------------------------------------------------------------------------
    ddDisp = 0.0d0
    j=1
    do i=1,tnodes*nds
        if(redcount(i)==1)then
            ddDisp(i,1)=dispeqLM(j,1)
            j=j+1
        end if
    end do
    
    !------------------------------------------------------------------------------
    ! Increment in the displacement increment
    !------------------------------------------------------------------------------
    dDisp = dDisp + ddDisp
    continue
    
    
    !------------------------------------------------------------------------------
    ! Strore Lagrange Multiplier
    !------------------------------------------------------------------------------
    do i=1,itot

        dLagrange(i) = dispeqLM(tnodes*nds-nbc+i,1)
        
        Lagrange(i) = Lagrange(i) + dLagrange(i)
        
        if(Lagrange(i) .gt. 0.0d0)then
            Lagrange(i) = 0.0d0
        end if
  
    end do
    
    !------------------------------------------------------------------------------
    ! Right-hand side contribution
    !------------------------------------------------------------------------------
    FintLM = 0.0d0
    do i=1,itot
        FintLM(:,1) = FintLM(:,1) + FLM(:,i)*Lagrange(i)
    end do
    
    !------------------------------------------------------------------------------
    ! Deallocate contact variables
    !------------------------------------------------------------------------------
    if(npair .gt. 1)then
        deallocate(gap)
    end if  
      
    deallocate(KLM)
    deallocate(FLM)
    deallocate(GLM)
    deallocate(disptempLM)
    deallocate(KeqLM)
    deallocate(FeqLM)
    deallocate(DispeqLM)
    deallocate(PTS_Active)
    
    if(npair==1)then
        deallocate(PTS_conn)
    end if
    
    continue
    
end subroutine






















            
            
!            !------------------------------------------------------------------------------
!            ! Closest Point Projection (CPP) for direction xi
!            !------------------------------------------------------------------------------
!            xib = 1.0d0
!            xib0 = 0.0d0
!            rsd_xi = 0.0d0
!            rsd0_xi = 0.0d0
!            do kn=1,knewton
!                if(kn==1)then
!                    !------------------------------------------------------------------------------
!                    ! Master NURBS surface basis, including vanishing terms
!                    !------------------------------------------------------------------------------
!                    call SurfaceBasis(ncpx_mst,p_mst,xib0,u_knot_mst, &
!                    &                 ncpy_mst,q_mst,etab0,v_knot_mst, &
!                    &                 GCoords, conn_mst, tnodes,nds, &
!                    &                 Rmst,dRmst,ddRmst)
!                    
!                    !------------------------------------------------------------------------------
!                    ! CPP in the physical space
!                    !------------------------------------------------------------------------------
!                    xm = 0.0d0
!                    count = 0
!                    do k2=1,ncpy_mst
!                        do k1=1,ncpx_mst
!                            count = count + 1
!                            do k3=1,nds
!                                xm(k3,1) = xm(k3,1) + Rmst(count)*GCoords(conn_mst(count),k3)
!                            end do
!                        end do
!                    end do
!                    
!                    !------------------------------------------------------------------------------
!                    ! Tangent vectors in the CPP
!                    !------------------------------------------------------------------------------
!                    txi = 0.0d0
!                    count = 0
!                    do k2=1,ncpy_mst
!                        do k1=1,ncpx_mst
!                            count = count + 1
!                            do k3=1,nds
!                                txi(k3,1)  =  txi(k3,1) + dRmst(count,1)*GCoords(conn_mst(count),k3)
!                            end do
!                        end do
!                    end do
!                    
!                    !------------------------------------------------------------------------------
!                    ! Auxilliary vector
!                    !------------------------------------------------------------------------------
!                    vec = xs - xm
!                    
!                    !------------------------------------------------------------------------------
!                    ! Initial residual
!                    !------------------------------------------------------------------------------
!                    rsd0_xi = 0.0d0
!                    do k1=1,nds
!                        rsd0_xi  = rsd0_xi +  txi(k1,1)*vec(k1,1)
!                    end do
!                end if
!                
!                !------------------------------------------------------------------------------
!                ! Master NURBS surface basis, including vanishing terms
!                !------------------------------------------------------------------------------
!                call SurfaceBasis(ncpx_mst,p_mst,xib,u_knot_mst, &
!                &                 ncpy_mst,q_mst,etab,v_knot_mst, &
!                &                 GCoords, conn_mst, tnodes,nds, &
!                &                 Rmst,dRmst,ddRmst)
!                
!                !------------------------------------------------------------------------------
!                ! CPP in the physical space
!                !------------------------------------------------------------------------------
!                xm = 0.0d0
!                count = 0
!                do k2=1,ncpy_mst
!                    do k1=1,ncpx_mst
!                        count = count + 1
!                        do k3=1,nds
!                            xm(k3,1) = xm(k3,1) + Rmst(count)*GCoords(conn_mst(count),k3)
!                        end do
!                    end do
!                end do
!                
!                !------------------------------------------------------------------------------
!                ! Tangent vectors in the CPP
!                !------------------------------------------------------------------------------
!                txi = 0.0d0
!                count = 0
!                do k2=1,ncpy_mst
!                    do k1=1,ncpx_mst
!                        count = count + 1
!                        do k3=1,nds
!                            txi(k3,1)  =  txi(k3,1) + dRmst(count,1)*GCoords(conn_mst(count),k3)
!                        end do
!                    end do
!                end do
!                
!                !------------------------------------------------------------------------------
!                ! Auxilliary vector
!                !------------------------------------------------------------------------------
!                vec = 0.0d0
!                vec = xs - xm
!                
!                !------------------------------------------------------------------------------
!                ! Residual
!                !------------------------------------------------------------------------------
!                rsd_xi = 0.0d0
!                do k1=1,nds
!                    rsd_xi  = rsd_xi  +  txi(k1,1)*vec(k1,1)
!                end do
!                
!                xibi = xib
!                
!                xib  =  xib - rsd_xi*( xib -  xib0)/(rsd_xi - rsd0_xi)
!                
!                
!                if(xib .gt. 1.0d0) xib = 1.0d0
!                if(xib .lt. 0.0d0) xib = 0.0d0
!                
!                xib0 = xibi
!                
!                rsd0_xi = rsd_xi
!            
!                if(abs(rsd_xi) .lt. 1.0d-10) goto 1
!                
!                continue
!            
!            end do
!            
!            1 continue
!            
!            continue
!            
!            !------------------------------------------------------------------------------
!            ! Closest Point Projection (CPP) for direction eta
!            !------------------------------------------------------------------------------
!            etab = 0.1d0
!            etab0 = 0.0d0
!            rsd_eta = 0.0d0
!            rsd0_eta = 0.0d0
!            do kn=1,knewton
!                if(kn==1)then
!                    !------------------------------------------------------------------------------
!                    ! Master NURBS surface basis, including vanishing terms
!                    !------------------------------------------------------------------------------
!                    call SurfaceBasis(ncpx_mst,p_mst,xib,u_knot_mst, &
!                    &                 ncpy_mst,q_mst,etab0,v_knot_mst, &
!                    &                 GCoords, conn_mst, tnodes,nds, &
!                    &                 Rmst,dRmst,ddRmst)
!                    
!                    !------------------------------------------------------------------------------
!                    ! CPP in the physical space
!                    !------------------------------------------------------------------------------
!                    xm = 0.0d0
!                    count = 0
!                    do k2=1,ncpy_mst
!                        do k1=1,ncpx_mst
!                            count = count + 1
!                            do k3=1,nds
!                                xm(k3,1) = xm(k3,1) + Rmst(count)*GCoords(conn_mst(count),k3)
!                            end do
!                        end do
!                    end do
!                    
!                    !------------------------------------------------------------------------------
!                    ! Tangent vectors in the CPP
!                    !------------------------------------------------------------------------------
!                    teta = 0.0d0
!                    count = 0
!                    do k2=1,ncpy_mst
!                        do k1=1,ncpx_mst
!                            count = count + 1
!                            do k3=1,nds
!                                teta(k3,1) = teta(k3,1) + dRmst(count,2)*GCoords(conn_mst(count),k3)
!                            end do
!                        end do
!                    end do
!                    
!                    !------------------------------------------------------------------------------
!                    ! Auxilliary vector
!                    !------------------------------------------------------------------------------
!                    vec = xs - xm
!                    
!                    !------------------------------------------------------------------------------
!                    ! Initial residual
!                    !------------------------------------------------------------------------------
!                    rsd0_eta = 0.0d0
!                    do k1=1,nds
!                        rsd0_eta = rsd0_xi + teta(k1,1)*vec(k1,1)
!                    end do
!                end if
!                
!                !------------------------------------------------------------------------------
!                ! Master NURBS surface basis, including vanishing terms
!                !------------------------------------------------------------------------------
!                call SurfaceBasis(ncpx_mst,p_mst,xib,u_knot_mst, &
!                &                 ncpy_mst,q_mst,etab,v_knot_mst, &
!                &                 GCoords, conn_mst, tnodes,nds, &
!                &                 Rmst,dRmst,ddRmst)
!                
!                !------------------------------------------------------------------------------
!                ! CPP in the physical space
!                !------------------------------------------------------------------------------
!                xm = 0.0d0
!                count = 0
!                do k2=1,ncpy_mst
!                    do k1=1,ncpx_mst
!                        count = count + 1
!                        do k3=1,nds
!                            xm(k3,1) = xm(k3,1) + Rmst(count)*GCoords(conn_mst(count),k3)
!                        end do
!                    end do
!                end do
!                
!                !------------------------------------------------------------------------------
!                ! Tangent vectors in the CPP
!                !------------------------------------------------------------------------------
!                teta = 0.0d0
!                count = 0
!                do k2=1,ncpy_mst
!                    do k1=1,ncpx_mst
!                        count = count + 1
!                        do k3=1,nds
!                            teta(k3,1) = teta(k3,1) + dRmst(count,2)*GCoords(conn_mst(count),k3)
!                        end do
!                    end do
!                end do
!                
!                !------------------------------------------------------------------------------
!                ! Auxilliary vector
!                !------------------------------------------------------------------------------
!                vec = 0.0d0
!                vec = xs - xm
!                
!                !------------------------------------------------------------------------------
!                ! Residual
!                !------------------------------------------------------------------------------
!                rsd_eta = 0.0d0
!                do k1=1,nds
!                    rsd_eta = rsd_eta + teta(k1,1)*vec(k1,1)
!                end do
!                
!                etabi = etab
!                
!                
!                etab = etab - rsd_eta*(etab - etab0)/(rsd_eta - rsd0_eta)
!                
!                if(etab .gt. 1.0d0) etab = 1.0d0
!                if(etab .lt. 0.0d0) etab = 0.0d0
!                
!                etab0 = etabi
!                
!                rsd0_eta = rsd_eta
!            
!                if(abs(rsd_eta) .lt. 1.0d-10) goto 2
!                
!                continue
!            
!            end do
!            
!            2 continue
