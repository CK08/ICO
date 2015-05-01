subroutine PTS_ComputeGap_testD(iter,inc,istep,iLM)

    use Mod_Variables ! Usar as variaveis de m*#rd@ que o Caseiro definiu 
    implicit none ! Porque é um menino, ai e tal para nao dar stress com as variaveis de sistema ou la o que e
    
    integer(4),intent(IN)::iter,inc,istep
    integer(4),intent(OUT)::iLM
    
    !------------------------------------------------------------------------
    ! Variaveis para calculos                                               ! 
    !------------------------------------------------------------------------
    integer(4)::i,j,k,l,m,k1,k2,k3,k4
    integer(4)::temp1,count,count1,count2,kspan,ilixo
    integer(4),parameter::knewton = 50
    !integer(4),dimension(ncp_m+1)::ctconn
    integer(4),dimension(ncp_s+ncp_m)::smconn
    integer(4),dimension(ncp_s)::isct
    real(8)::cm,c1,c2,c3,c4
    
    real(8)::xib,lm,dlm,sumres,rsd,gn,a11,b11,gt,num,tst,xib_new
    real(8)::xibi,xibf,diff
    real(8),dimension(p_s+1)::Rss,dRss,ddRss
    real(8),dimension(p_m+1)::Rmm,dRmm,ddRmm
    real(8),dimension(ncp_s)::Rs,dRsdxi,d2Rsdxi2
    real(8),dimension(ncp_m)::Rm,dRmdxi,d2Rmdxi2,sumRm
    real(8),dimension(nds,1)::tang,dtang,nor,dpos,aux21,d2xdxi2
    real(8),dimension(nds,1)::xs,dxs,xm,dxm,d2xm,vec,dvec,vec0
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
    
    if(iter==1) dLagrange = 0.0d0
    
    ! Ciclo a cada ponto slave Greville
    do i=1,ngrv
            !------------------------------------------------!
            ! Determine the slave point physical coordenates !
            !------------------------------------------------!
            
            ! Calcular as Basis Functions e as derivadas
            call BSplineBasisAndDeriv2(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan) 
            ! Converter para o formato com os termos zero
            call BSplinesToNURBS(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
            
            ! Alocar as coordenadas fisicas do ponto slave
            xs = 0.0d0 
            do j=1,ncp_s
                    xs(1,1) = xs(1,1) + Rs(j)*(GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*nds-1),1))
                    xs(2,1) = xs(2,1) + Rs(j)*(GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*nds  ),1))
                    dxs(1,1) = dxs(1,1) + dRsdxi(j)*(GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*nds-1),1))
                    dxs(2,1) = dxs(2,1) + dRsdxi(j)*(GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*nds  ),1))
            end do        
            !-----------------------------------------------------!
            ! Clossest Point Projection (CPP) of the master point !
            !                    Newton Method                    !     
            !-----------------------------------------------------!
            xib0 = 0.0d0 
            xib = 0.123581d0
            do j=1,knewton             
                    ! 1a passagem: Definir um CPP inicial
                    if(j==1)then
                            ! Calcular as Basis Functions e Derivadas 
                            call BSplineBasisAndDeriv2(ncp_m,p_m,xib0,u_knot_master,Rmm,dRmm,ddRmm,kspan)
                            ! Converter para o formato com vanishing terms
                            call BSplinesToNURBS(ncp_m,p_m,xib0,u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
                            ! Calcular as coordenadas fisicas do ponto master e as derivadas 
                            xm = 0.0d0 
                            dxm = 0.0d0 
                            do k=1,ncp_m
                                    xm(1,1) = xm(1,1) + Rm(k)*(GCoords(conn_master(k),1) + dDisp((conn_master(k)*nds-1),1))
                                    xm(2,1) = xm(2,1) + Rm(k)*(GCoords(conn_master(k),2) + dDisp((conn_master(k)*nds  ),1))
                                    dxm(1,1) = dxm(1,1) + dRmdxi(k)*(GCoords(conn_master(k),1) + dDisp((conn_master(k)*nds-1),1))
                                    dxm(2,1) = dxm(2,1) + dRmdxi(k)*(GCoords(conn_master(k),2) + dDisp((conn_master(k)*nds  ),1))
                            end do        
                            ! Calcular o comprimento do segmento
                            lm = dsqrt(dxm(1,1)*dxm(1,1) + dxm(2,1)*dxm(2,1))  
                            ! Calcular o vector tangente 
                            tang = dxm/lm
                            ! Calculos auxiliares (diferenca espacial entre slave e master, e residual 0 )
                            vec = xs - xm
                            rsd0 = vec(1,1)*tang(1,1) + vec(2,1)*tang(2,1)               
                    end if
                    ! Calcular as Basis Functions e Derivadas 
                    call BSplineBasisAndDeriv2(ncp_m,p_m,xib,u_knot_master,Rmm,dRmm,ddRmm,kspan)
                    ! Converter para o formato com vanishing terms 
                    call BSplinesToNURBS(ncp_m,p_m,xib,u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
                    ! Calcular as coordenadas fisicas do ponto master e as derivadas 
                    xm = 0.0d0 
                    dxm = 0.0d0 
                    d2xm = 0.0d0
                    do k=1,ncp_m
                            xm(1,1) = xm(1,1) + Rm(k)*(GCoords(conn_master(k),1) + dDisp((conn_master(k)*nds-1),1))
                            xm(2,1) = xm(2,1) + Rm(k)*(GCoords(conn_master(k),2) + dDisp((conn_master(k)*nds  ),1))
                            dxm(1,1) = dxm(1,1) + dRmdxi(k)*(GCoords(conn_master(k),1) + dDisp((conn_master(k)*nds-1),1))
                            dxm(2,1) = dxm(2,1) + dRmdxi(k)*(GCoords(conn_master(k),2) + dDisp((conn_master(k)*nds  ),1))
                            d2xm(1,1) = d2xm(1,1) + d2Rmdxi2(k)*(GCoords(conn_master(k),1) + dDisp((conn_master(k)*nds-1),1))
                            d2xm(2,1) = d2xm(2,1) + d2Rmdxi2(k)*(GCoords(conn_master(k),2) + dDisp((conn_master(k)*nds  ),1))
                    end do              
                    ! Calcular o comprimento do segmento
                    lm  = dsqrt(dxm(1,1)*dxm(1,1) + dxm(2,1)*dxm(2,1))
                    dlm = dsqrt(d2xm(1,1)*d2xm(1,1) + d2xm(2,1)*d2xm(2,1)) 
                    ! Calcular o vector tangente 
                    tang  = dxm/lm  
                    dtang = dxm/dlm       
                    ! Calcular as novas coordenadas parametricas 
                    vec = xs - xm
                    rsd = vec(1,1)*tang(1,1) + vec(2,1)*tang(2,1)
                    
                    dvec = (dxs - dxm)
                    tst = dsqrt(dvec(1,1)*dtang(1,1) + dvec(2,1)*dtang(2,1)) 
                    
                    xib_new = xib 
                    xib = xib - rsd*(xib-xib0)/(rsd-rsd0) 
                    ! Nos limites de XIb         
                    if(xib .gt. 1.0d0) xib = 1.0d0 
                    if(xib .lt. 0.0d0) xib = 0.0d0
                    ! Alocar variaveis para uma proxima iteracao
                    xib0 = xib_new 
                    rsd0 = rsd         
                    ! Verificar a viabilidade dos novos pontos
                    if(abs(rsd) .lt. 1.0d-8) goto 1
            end do
        
            ! Error case: CPP not found
            write(*,FMT=16)i
            16 format('Warning: Failed to determine CPP for Greville point ',I3)
            ! Continue: no error found
            1 continue 
            
            ! Calcular a normal do master: tangent x {0 0 -1}
            nor(1,1) = -tang(2,1)
            nor(2,1) =  tang(1,1)
            ! Alocar o gap anterior e calcular o novo gap final 
            initial_gap(i) = final_gap(i)
            final_gap(i) = vec(1,1)*nor(1,1) + vec(2,1)*nor(2,1)
            ! Verificacao de contacto atraves do final gap e da condicao de separacao (Matzen dissertation) 
            if (final_gap(i) <= -1.0d-8) then 
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
            !-----------------------!
            ! Caso exista contacto: !
            !-----------------------!
            Kdelta = 0.0d0 
            if(PTS_active(i) == 1) then
                    gn = final_gap(i)
                    do k=1,ncp_m
                            d2xdxi2(1,1) = d2xdxi2(1,1) + d2Rmdxi2(k)*((GCoords(conn_master(k),1)) + ddisp(conn_master(k)*nds-1,1))
                            d2xdxi2(2,1) = d2xdxi2(2,1) + d2Rmdxi2(k)*((GCoords(conn_master(k),2)) + ddisp(conn_master(k)*nds  ,1))
                    end do  
                    call KLM_Matzen(i,iLM,lm,gn,nor,tang,d2xdxi2,Rs,Rm,dRmdxi,a11,b11,Ns,N0s,Ts,T0s,Kdelta,FLM,KLM,GLM)
                    
                    ! Alocar o Residuo de contacto GLM (=-1.0d0*final_gap(i))    
                    GLM(tnodes*nds+i,1) = -1.0d0*final_gap(i)
            else
            !---------------------------!             
            ! Caso nao exista contacto: !
            !---------------------------!
                                
                   call KLM_NoContact(i,iLM,KLM,FLM)
            end if        
                                   

    end do !End of Greville cycle 
    
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
    
    end subroutine
