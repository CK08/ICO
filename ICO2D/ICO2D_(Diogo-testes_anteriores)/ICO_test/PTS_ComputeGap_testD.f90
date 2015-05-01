subroutine PTS_ComputeGap(iter,inc,istep,iLM)

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
    
    if(iter==1) dLagrange = 0.0d0
    
    ! Ciclo a cada ponto slave Greville
    do i=1, ngrv
            !------------------------------------------------!
            ! Determine the slave point physical coordenates !
            !------------------------------------------------!
            
            ! Calcular as Basis Functions e as derivadas
            call BSplineBasisAndDeriv2(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan) 
            ! Converter para o formato com os termos zero
            call BSplineToNURBS(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
            ! Alocar as coordenadas fisicas do ponto slave
            xs = 0.0d0 
            do j=1,ncp_s
                    xs(1,1) = xs(1,1) + Rs(j)*(GCoords(conn_slave(j),1) + dDisp((conn_slave(j)*nds-1),1))
                    xs(1,1) = xs(2,1) + Rs(j)*(GCoords(conn_slave(j),2) + dDisp((conn_slave(j)*nds  ),1))
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
                            call BSplineToNURBS(ncp_m,p_m,xib0,u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
                            ! Calcular as coordenadas fisicas do ponto master e as derivadas 
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
                    call BSplineToNURBS(ncp_m,p_m,xib,u_knot_master,Rmm,dRmm,ddRmm,kspan,conn_master,Rm,dRmdxi,d2Rmdxi2)
                    ! Calcular as coordenadas fisicas do ponto master e as derivadas 
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
                    ! Calcular as novas coordenadas parametricas 
                    vec = xs - xm
                    rsd = vec(1,1)*tang(1,1) + vec(2,1)*tang(2,1)
                    xib_new = xib - rsd*((xib-xib0)/(rsd-rsd0))
                    ! Alocar variaveis para uma proxima iteracao
                    xib0 = xib 
                    xib = xib_new
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
                    call KLM_Matzen(iLM,lm,gn,nor,tang,d2xdxi2,Rs,Rm,dRmdxi,a11,b11,Ns,N0s,Ts,T0s,Kdelta,FLM,KLM,GLM)
                    
            else
            !---------------------------!             
            ! Caso nao exista contacto: !
            !---------------------------!
                                
                    ! Incrementar o contador iLM  
                    
                    ! Contruir a matrix KLM com diagonal =1 e resto 0
                    
                    ! Construir connectividades (caso deformavel-deformavel)
                    
                    ! Alocar no vector de FintLM valor 0 para os dois graus do ponto de Greville, i
            end if        
            ! Alocar o Residuo de contacto GLM (=-1.0d0*final_gap(i))
                            
    end do !End of Greville cycle 
    
    
    
    end subroutine
