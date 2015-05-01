!subroutine ElHex27EAS(nel,Ke,el_ddisp,Finte,inc)
!    
!    use Mod_Variables
!    implicit none
!    
!    integer(4)::nel,cpi,count,inc
!    
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::el_ddisp
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(OUT)::Ke
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(OUT)::Finte
!    
!    integer(4),parameter::npi_xi= 3, npi_eta = 3, npi_zeta = 3
!    integer(4),parameter:: nlph = 15
!    
!    integer(4)::i,j,k,ni,nj,nk,nel_nza,k1
!    real(8)::xi,eta,zeta
!    real(8),dimension(npi_xi)::e,we
!    real(8),dimension(npi_eta)::n,wn
!    real(8),dimension(npi_zeta)::c,wc
!    
!    real(8),dimension((p+1)*(q+1)*(w+1))::R,R_mid
!    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdx,dRdx_mid
!
!    real(8)::detJ,detj_mid,gwt
!    
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob,Bloc
!    
!    real(8),dimension(6,6)::matD
!    real(8),dimension(ntens,1)::tempstr,deform
!    
!    real(8)::density
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1)::Ftest
!    
!    real(8),dimension(nds,nds)::jac
!    
!    !Geometric Nonlinearity
!    integer(4)::inodes,ilixo
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1)::tdisp,mdisp,tzero
!    real(8),dimension(nds,nds)::rconv,rconv_mid,jacinv,raux
!    real(8),dimension(nds,nds)::jac_mid,jacinv_mid
!    real(8),dimension(nds,nds)::matF,matF_mid,matR,matR_mid
!    real(8),dimension((p+1)*(q+1)*(w+1))::dRde,dRdn,dRdc
!    real(8),dimension(nds)::temp1
!    real(8),dimension((p+1)*(q+1)*(w+1),nds)::nodes
!    real(8),dimension(nds,nds)::Temp33_mid,Temp33
!    real(8),dimension(nds,nds)::Trans_mid,Trans
!    real(8),dimension(nds*2,nds*2)::TCL_mid,TGL_mid
!    real(8),dimension(nds*2,nds*2)::TCL,TGL
!    real(8),dimension(9,(p+1)*(q+1)*(w+1)*nds)::BNL
!    real(8),dimension(9,9)::MSTR
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds)::KNLG
!    real(8),dimension(nds,1)::dRdenc,dRdxyz
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob_mid,Bloc_mid
!    
!    real(8)::detj0
!    real(8),dimension(nlph)::tempa,diagon
!    real(8),dimension(nds*2,nds*2)::MatT0
!    real(8),dimension(6,nlph)::Malpha,Balpha,MAT
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,nlph)::Kua
!    real(8),dimension(nlph,nlph)::Kaa,Kaainv
!    
!    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdxi,dRdxii
!    real(8),dimension(nds,nds)::dxdxi
!    real(8)::xib,etab,zetab,summa
!    
!    inodes = (p+1)*(q+1)*(w+1)
!    density = props(3)
!    
!    tzero = 0.0d0
!    
!    Mat = 0.0d0
!    
!    cpi = 0
!
!    we = 1.0d0
!    wn = 1.0d0
!    wc = 1.0d0
!    call gauleg(npi_xi, e, we)
!    call gauleg(npi_eta, n, wn)
!    call gauleg(npi_zeta, c, wc)
!    
!    !#####
!    if(nlgeom==.true.)then
!        tdisp = el_ddisp
!        mdisp = el_ddisp/2.0d0
!    else
!        tdisp = 0.0d0
!        mdisp= 0.0d0
!        
!        TGL = 0.0d0
!        do i=1,6
!            TGL(i,i) = 1.0d0
!        end do
!    end if
!    !#####
!    
!    !Transformation Matrix at the centre of the element
!    !call ShapeFunc2(nel,0.0d0,0.0d0,0.0d0,R,dRdx,detj0,jac,tzero)
!    call ShapeFunc3(nel,0.0d0,0.0d0,0.0d0,R,dRdx,dRdxi,dRdxii,detj0,jac,dxdxi,tzero)
!    
!    ni = INN(IEN(nel,1),1)
!    nj = INN(IEN(nel,1),2)
!    nk = INN(IEN(nel,1),3)
!    
!    !detj0=detj0/((U_knot(ni+1)-U_knot(ni))*(V_knot(nj+1)-V_knot(nj))*(W_knot(nk+1)-W_knot(nk))/8.0d0)
!    
!    temp1 = 0.0d0
!    jacinv = 0.0d0
!    jacinv = jac
!    call gaussj(jacinv,nds,temp1,ilixo)
!    
!    call TransformationMat3D(JacInv,MatT0)
!    
!    
!    !Gauss points cycle
!    Ke = 0.0d0
!    KNLG = 0.0d0
!    BNL = 0.0d0
!    MSTR = 0.0d0
!    
!    Kua = 0.0d0
!    Kaa = 0.0d0
!    
!    Finte = 0.0d0
!    do i=1,npi_xi
!        do j=1,npi_eta
!            do k=1,npi_zeta
!            
!                cpi = cpi + 1
!                MatD(:,:) = TmatD(nel,cpi,:,:)
!                
!                xi = e(i)
!                eta = n(j)
!                zeta = c(k)
!                
!                call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tzero)
!                
!                !Weight factor
!                gwt = we(i)*wn(j)*wc(k)*detj
!                
!                
!                !#########
!                if(nlgeom==.true.)then
!                    
!                    !Local Axis ---------------------------
!                    if(inc==1)then
!                        call local_axis(nds,jac,rconv)
!                    else
!                        rconv = laxis(nel,cpi,:,:)
!                    end if
!                    
!                    !-------------------------------------------------------------------------
!                    !Mid-Configuration
!                    !-------------------------------------------------------------------------
!                    matF_mid = 0.0d0
!                    do k1=1,inodes
!                        dRde(k1) = dRdx(k1,1)
!                        dRdn(k1) = dRdx(k1,2)
!                        dRdc(k1) = dRdx(k1,3)
!                    end do
!                    
!                    jacinv = 0.0d0
!                    jacinv = jac
!                    call gaussj(jacinv,nds,temp1,ilixo)
!                    
!                    count = 0
!                    do k1=1,inodes
!                        count = count + 1
!                        nodes(k1,1) = Points(IEN(nel,count),1) !+ u(IEN(nel,count)*3-2,1)
!                        nodes(k1,2) = Points(IEN(nel,count),2) !+ u(IEN(nel,count)*3-1,1)
!                        nodes(k1,3) = Points(IEN(nel,count),3) !+ u(IEN(nel,count)*3  ,1)
!                    end do
!                    
!                    !Deformation Gradient for mid-configuration
!                    call DefGrad3D(inodes,nds,dRde,dRdn,dRdc,jacinv,nodes,mdisp,MatF_mid)
!                    
!                    !Polar Decomposition for mid configuration
!                    matR_mid=0.0d0
!                    call PolarDecomp2D(nds,matF_mid,matR_mid)
!                    
!                    !Update local reference system to the mid configuration
!                    call axisupdate(rconv,matR_mid,rconv_mid)
!                    
!                    !Global to natural jacobian in mid configuration
!                    call ShapeFunc2(nel,xi,eta,zeta,R_mid,dRdx_mid,detj_mid,jac_mid,mdisp)
!                    
!                    jacinv_mid = 0.0d0
!                    temp1 = 0.0d0
!                    jacinv_mid = jac_mid
!                    call gaussj(jacinv_mid,nds,temp1,ilixo)
!                    
!                    !Transformation matrices mid configuration 
!                    Temp33_mid=0.0d0
!                    Temp33_mid= matmul(transpose(rconv_mid),jacinv_mid)
!                    Trans_mid=0.0d0
!                    Trans_mid=transpose(rconv_mid)
!                    
!                    !Natural to local transformation matrix for mid configuration
!                    call TransformationMat3D(Temp33_mid,TCL_mid)
!                    
!                    !Global to local transformation matrix for mid configuration
!                    call TransformationMat3D(Trans_mid,TGL_mid)
!                    
!                    Bglob_mid = 0.0d0
!                    do k1=1,(p+1)*(q+1)*(w+1)
!                        Bglob_mid(1,k1*3-2) = dRdx_mid(k1,1)
!                        Bglob_mid(2,k1*3-1) = dRdx_mid(k1,2)
!                        Bglob_mid(3,k1*3  ) = dRdx_mid(k1,3)
!                        Bglob_mid(4,k1*3-2) = dRdx_mid(k1,2)
!                        Bglob_mid(4,k1*3-1) = dRdx_mid(k1,1)
!                        Bglob_mid(5,k1*3-2) = dRdx_mid(k1,3)
!                        Bglob_mid(5,k1*3  ) = dRdx_mid(k1,1)
!                        Bglob_mid(6,k1*3-1) = dRdx_mid(k1,3)
!                        Bglob_mid(6,k1*3  ) = dRdx_mid(k1,2)
!                    end do
!                    
!                    Bloc_mid=0.0d0
!                    Bloc_mid=matmul(TGL_mid,Bglob_mid)
!                    !-------------------------------------------------------------------------
!                    !End-Configuration
!                    !-------------------------------------------------------------------------
!!                    matF = 0.0d0
!!                    do k1=1,inodes
!!                        dRde(k1) = dRdx(k1,1)
!!                        dRdn(k1) = dRdx(k1,2)
!!                        dRdc(k1) = dRdx(k1,3)
!!                    end do
!!                    
!!                    jacinv = 0.0d0
!!                    jacinv = jac
!!                    call gaussj(jacinv,nds,temp1,ilixo)
!!                    
!!                    count = 0
!!                    do k1=1,inodes
!!                        count = count + 1
!!                        nodes(k1,1) = Points(IEN(nel,count),1)
!!                        nodes(k1,2) = Points(IEN(nel,count),2)
!!                        nodes(k1,3) = Points(IEN(nel,count),3)
!!                    end do
!                    
!                    !Deformation Gradient for mid-configuration
!                    call DefGrad3D(inodes,nds,dRde,dRdn,dRdc,jacinv,nodes,tdisp,MatF)
!                    
!                    !Polar Decomposition for mid configuration
!                    matR = 0.0d0
!                    call PolarDecomp2D(nds,matF,matR)
!                    
!                    !Update local reference system to the mid configuration
!                    raux = rconv
!                    call axisupdate(rconv,matR,rconv)
!                    
!                    !Global to natural jacobian in mid configuration
!                    call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tdisp)
!                    
!                    jacinv = 0.0d0
!                    temp1 = 0.0d0
!                    jacinv = jac
!                    call gaussj(jacinv,nds,temp1,ilixo)
!                    
!                    !Transformation matrices mid configuration 
!                    Temp33 = 0.0d0
!                    Temp33 = matmul(transpose(rconv),jacinv)
!                    Trans = 0.0d0
!                    Trans = transpose(rconv)
!                    
!                    !Natural to local transformation matrix for mid configuration
!                    call TransformationMat3D(Temp33,TCL)
!                    
!                    !Global to local transformation matrix for mid configuration
!                    call TransformationMat3D(Trans,TGL)
!                    
!                    !Update Weight factor
!                    gwt = we(i)*wn(j)*wc(k)*detj
!                    
!                    continue
!                    
!                end if
!                !############
!                
!                !Good Results with this one for volumetric locking!
!!                Malpha = 0.0d0
!!                Malpha(1,1) = 3.0d0*xi**2.0d0 - 1.0d0
!!                Malpha(2,2) = 3.0d0*eta**2.0d0 - 1.0d0
!!                Malpha(3,3) = 3.0d0*zeta**2.0d0 - 1.0d0
!!                
!!                Malpha(1,4) = xi
!!                Malpha(2,5) = eta
!!                Malpha(3,6) = zeta
!!                
!!                Malpha(1,7) = xi*eta
!!                Malpha(1,8) = xi*zeta
!!                
!!                Malpha(2,9) = xi*eta
!!                Malpha(2,10) = eta*zeta
!!                
!!                Malpha(3,11) = xi*zeta
!!                Malpha(3,12) = eta*zeta
!!                
!!                Malpha(1,13) = xi*eta*zeta
!!                Malpha(2,14) = xi*eta*zeta
!!                Malpha(3,15) = xi*eta*zeta
!                ! --------------------------------------------------------
!                xib = xi !(xi+1.0d0)/2.0d0
!                etab = eta !(eta+1.0d0)/2.0d0
!                zetab = zeta !(zeta+1.0d0)/2.0d0
!                
!!                ni = INN(IEN(nel,1),1)
!!                nj = INN(IEN(nel,1),2)
!!                nk = INN(IEN(nel,1),3)
!!                
!!                xib   = ((u_knot(ni+1) - u_knot(ni))*xi   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                etab  = ((v_knot(nj+1) - v_knot(nj))*eta  + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                zetab = ((w_knot(nk+1) - w_knot(nk))*zeta + (w_knot(nk+1) + w_knot(nk)))/2.0d0
!                
!                Malpha = 0.0d0
!                Malpha(1,1) = 3.0d0*xib**2.0d0 - 1.0d0
!                Malpha(2,2) = 3.0d0*etab**2.0d0 - 1.0d0
!                Malpha(3,3) = 3.0d0*zetab**2.0d0 - 1.0d0
!                
!                Malpha(1,4) = xib
!                Malpha(2,5) = etab
!                Malpha(3,6) = zetab
!                
!                Malpha(1,7) = xib*etab
!                Malpha(1,8) = xib*zetab
!                
!                Malpha(2,9) = xib*etab
!                Malpha(2,10) = etab*zetab
!                
!                Malpha(3,11) = xib*zetab
!                Malpha(3,12) = etab*zetab
!                
!                Malpha(1,13) = xib*etab*zetab
!                Malpha(2,14) = xib*etab*zetab
!                Malpha(3,15) = xib*etab*zetab
!                     
!!                Malpha(4,16) = 3.0d0*xi*xi - 1.0d0
!!                Malpha(4,17) = 3.0d0*eta*eta - 1.0d0
!!                Malpha(5,18) = 3.0d0*xi*xi - 1.0d0
!!                Malpha(5,19) = 3.0d0*zeta*zeta - 1.0d0
!!                Malpha(6,20) = 3.0d0*eta*eta - 1.0d0
!!                Malpha(6,21) = 3.0d0*zeta*zeta - 1.0d0
!!                
!!                Malpha(4,22) = xi
!!                Malpha(4,23) = eta
!!                Malpha(5,24) = xi
!!                Malpha(5,25) = zeta
!!                Malpha(6,26) = eta
!!                Malpha(6,27) = zeta
!!
!!                
!!                Malpha(4,28) = xi*zeta
!!                Malpha(4,29) = eta*zeta
!!                Malpha(5,30) = xi*eta
!!                Malpha(5,31) = eta*zeta
!!                Malpha(6,32) = xi*eta
!!                Malpha(6,33) = xi*zeta
!!                
!!                Malpha(5,15) = xi*eta
!!                Malpha(5,16) = eta*zeta
!!                
!!                Malpha(6,17) = xi*eta
!!                Malpha(6,18) = xi*zeta
!!                
!!                Malpha(4,13) = xi*eta
!!                Malpha(5,14) = xi*zeta
!!                Malpha(6,15) = eta*zeta
!!                
!!                Malpha(4,16) = xi*eta*zeta
!!                Malpha(5,17) = xi*eta*zeta
!!                Malpha(6,18) = xi*eta*zeta
!!                
!!                Malpha(1,1) = -xi*(1.0-eta*eta)*(1.0-zeta*zeta);
!!                Malpha(2,2) = -eta*(1.0-xi*xi)*(1.0-zeta*zeta);
!!                Malpha(3,3) = -zeta*(1.0-xi*xi)*(1.0-eta*eta);
!!                
!!                Malpha(1,4) = 2.0d0*xi*eta*(1.0-zeta*zeta);
!!                Malpha(2,4) = 2.0d0*xi*eta*(1.0-zeta*zeta);
!!                Malpha(3,4) = 2.0d0*xi*eta*(1.0-zeta*zeta);
!!                
!!                Malpha(1,5) = 2.0d0*xi*zeta*(1.0-eta*eta);
!!                Malpha(2,5) = 2.0d0*xi*zeta*(1.0-eta*eta);
!!                Malpha(3,5) = 2.0d0*xi*zeta*(1.0-eta*eta);
!!                
!!                Malpha(1,6) = 2.0d0*eta*zeta*(1.0-xi*xi);
!!                Malpha(2,6) = 2.0d0*eta*zeta*(1.0-xi*xi);
!!                Malpha(3,6) = 2.0d0*eta*zeta*(1.0-xi*xi);
!!            
!!                Malpha(1,7) = 4.0*xi*eta*zeta;
!!                Malpha(2,8) = 4.0*xi*eta*zeta;
!!                Malpha(3,9) = 4.0*xi*eta*zeta;
!!           
!!                Malpha(1,10) = -(1.0-eta*eta)*(1.0-zeta*zeta);
!!                Malpha(2,11) = -(1.0-xi*xi)*(1.0-zeta*zeta);
!!                Malpha(3,12) = -(1.0-xi*xi)*(1.0-eta*eta);
!!                
!!                Malpha(4,1) = -xi*(1.0-eta*eta)*(1.0-zeta*zeta);
!!                Malpha(4,2) = -eta*(1.0-xi*xi)*(1.0-zeta*zeta);
!!                Malpha(5,3) = -xi*(1.0-eta*eta)*(1.0-zeta*zeta);
!!                Malpha(5,4) = -zeta*(1.0-xi*xi)*(1.0-eta*eta);
!!                Malpha(6,5) = -eta*(1.0-xi*xi)*(1.0-zeta*zeta);
!!                Malpha(6,6) = -zeta*(1.0-xi*xi)*(1.0-eta*eta);
!!                
!!                Malpha(4,19) = 2.0d0*xi*zeta*(1.0-eta*eta);
!!                Malpha(4,20) = 2.0d0*eta*zeta*(1.0-xi*xi);
!!                Malpha(5,21) = 2.0d0*xi*eta*(1.0-zeta*zeta);
!!                Malpha(5,22) = 2.0d0*eta*zeta*(1.0-xi*xi);
!!                Malpha(6,23) = 2.0d0*xi*eta*(1.0-zeta*zeta);
!!                Malpha(6,24) = 2.0d0*xi*zeta*(1.0-eta*eta);
!                
!                Balpha = 0.0d0
!                Balpha = matmul(MatT0,Malpha)*detj0/detj   !*props(1)/(3.0d0*(1.0d0-2.0d0*0.499999d0));
!                
!                Mat = Mat + Balpha*gwt
!                
!                Bglob = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
!                    
!                    Bglob(1,k1*3-2) = dRdx(k1,1)
!                    Bglob(2,k1*3-1) = dRdx(k1,2)
!                    Bglob(3,k1*3  ) = dRdx(k1,3)
!                    
!                    Bglob(4,k1*3-2) = dRdx(k1,2)
!                    Bglob(4,k1*3-1) = dRdx(k1,1)
!                    Bglob(5,k1*3-2) = dRdx(k1,3)
!                    Bglob(5,k1*3  ) = dRdx(k1,1)
!                    Bglob(6,k1*3-1) = dRdx(k1,3)
!                    Bglob(6,k1*3  ) = dRdx(k1,2)
!                    
!                end do
!                
!                !Local strain-displacement matrix at end configuration
!                Bloc=0.0d0
!                Bloc=matmul(TGL,Bglob)
!                
!                if(nlgeom == .false.) Bloc_mid = Bloc
!                
!                call MatPlastic3D(iprops,nds,(p+1)*(q+1)*(w+1),props,Bloc_mid,6,el_ddisp,stress(nel,cpi,:),strain(nel,cpi,:),dstrain(nel,cpi,:),dstress(nel,cpi,:),hard(nel,cpi),matD)
!                
!                !Store Constitutive Matrix
!                TmatD(nel,cpi,:,:) = MatD(:,:)
!                
!                tempstr = 0.0d0
!                do k1=1,ntens
!                    tempstr(k1,1) = stress(nel,cpi,k1)
!                    deform(k1,1) = dstrain(nel,cpi,k1)
!                end do
!                
!                !#########
!                if (nlgeom == .true.) then
!                
!                    BNL = 0.0d0
!                    do k1=1,inodes
!                        dRdenc=0.0d0
!                        dRdenc(1,1)=dRdx(k1,1) !
!                        dRdenc(2,1)=dRdx(k1,2) !
!                        dRdenc(3,1)=dRdx(k1,3) !
!                        
!                        dRdxyz = 0.0d0
!                        dRdxyz = matmul(trans,dRdenc)
!
!                        BNL(1,k1*3-2)=dRdxyz(1,1)
!                        BNL(2,k1*3-2)=dRdxyz(2,1)
!                        BNL(3,k1*3-2)=dRdxyz(3,1)
!                        BNL(4,k1*3-1)=dRdxyz(1,1)
!                        BNL(5,k1*3-1)=dRdxyz(2,1)
!                        BNL(6,k1*3-1)=dRdxyz(3,1)
!                        BNL(7,k1*3  )=dRdxyz(1,1)
!                        BNL(8,k1*3  )=dRdxyz(2,1)
!                        BNL(9,k1*3  )=dRdxyz(3,1)
!                    end do 
!                 
!                    MSTR=0.0d0
!                    do k1=1,nds
!                        MSTR(k1*3-2,k1*3-2) = tempstr(1,1)
!                        MSTR(k1*3-2,k1*3-1) = tempstr(4,1)
!                        MSTR(k1*3-2,k1*3  ) = tempstr(5,1)
!                        MSTR(k1*3-1,k1*3-2) = tempstr(4,1)
!                        MSTR(k1*3-1,k1*3-1) = tempstr(2,1)
!                        MSTR(k1*3-1,k1*3  ) = tempstr(6,1)
!                        MSTR(k1*3  ,k1*3-2) = tempstr(5,1)
!                        MSTR(k1*3  ,k1*3-1) = tempstr(6,1)
!                        MSTR(k1*3  ,k1*3  ) = tempstr(3,1)
!                    end do
!                    
!                    KNLG = KNLG + matmul(matmul(transpose(BNL),MSTR),BNL)*gwt
!                    
!                    !Export local axis ---------------
!                    laxis(nel,cpi,:,:) = rconv
! 
!                end if
!                !###########
!                
!                if(nlgeom == .false.) then
!                    Ke = Ke + matmul(matmul(transpose(Bloc),matD),Bloc)*gwt
!                else
!                    Ke = Ke + matmul(matmul(transpose(Bloc),matD),Bloc)*gwt + matmul(matmul(transpose(BNL),MSTR),BNL)*gwt
!                end if
!                
!                
!                Kua = Kua + matmul(matmul(transpose(Bloc),MatD),Balpha)*gwt
!                Kaa = Kaa + matmul(matmul(transpose(Balpha),MatD),Balpha)*gwt
!                
!                
!                
!                Finte = Finte + matmul(transpose(Bloc),tempstr)*gwt
!                
!                
!                !Gravitic loads -------------------------------------------------------------
!                if (gravity==.true.)then
!                    call grvt(nshpl,nds,R,gwt,gconst,gravdir,density,Finte)
!                end if
!                
!                !Pressure loads -------------------------------------------------------------
!                if(pressure==.true.) then
!                    do k1=1,nelems
!                        if(k1==nel .and. abs(pressvec(k1)) .gt. 0.0d0) then
!                            call prssr(nel,npi_xi,npi_eta,npi_zeta,jac,we(i),wn(j),wc(k),R,Finte)
!                        end if
!                    end do
!                end if
!                
!                !Strain Energy --------------------------------------------------------------
!                do k1=1,ntens
!                    SEnergy = SEnergy + tempstr(k1,1)*deform(k1,1)*gwt/2.0d0
!                end do
!                
!                continue
!                
!            end do  
!        end do
!    end do
!    
!!    diagon = 0.0d0
!!    do i=1,nlph
!!        diagon(i) = sum(Kaa(:,i))
!!    end do
!!    
!!    Kaa = 0.0d0
!!    do i=1,nlph
!!        Kaa(i,i) = diagon(i)
!!        !Kaainv(i,i) = 1.0d0/diagon(i)
!!    end do
!
!    summa = 0.0d0
!    summa = sum(Mat)
!    
!    tempa = 0.0d0
!    Kaainv = 0.0d0
!    Kaainv = Kaa
!    call gaussj(Kaainv,nlph,tempa,ilixo)
!    
!    Ke = Ke - matmul(matmul(Kua,Kaainv),transpose(Kua))
!    
!    continue
!
!end subroutine
!
!
!
!
!!------------------------------------------------------------------------------------------------------------
!!
!! Second order hexahedral element with EAS
!!
!! Input: nel - Element number
!!        el_ddisp - Element displacement
!!        inc - Increment number
!!
!! Output: Ke - Elementar stiffness matrix
!!         Finte - Elementar internal forces
!!------------------------------------------------------------------------------------------------------------
!
!subroutine ElHex27EAS(nel,Ke,el_ddisp,Finte,inc)
!    
!    use Mod_Variables
!    implicit none
!    
!    integer(4)::nel,cpi,count,inc
!    
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::el_ddisp
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(OUT)::Ke
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(OUT)::Finte
!    
!    integer(4),parameter::npi_xi= 3, npi_eta = 3, npi_zeta = 3
!    
!    integer(4)::i,j,k,ni,nj,nk,nel_nza,k1
!    real(8)::xi,eta,zeta
!    real(8),dimension(npi_xi)::e,we
!    real(8),dimension(npi_eta)::n,wn
!    real(8),dimension(npi_zeta)::c,wc
!    
!    real(8),dimension((p+1)*(q+1)*(w+1))::R,R_mid
!    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdx,dRdx_mid
!
!    real(8)::detJ,detj_mid,gwt
!    
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob,Bloc
!    
!    real(8),dimension(6,6)::matD
!    real(8),dimension(ntens,1)::tempstr,deform
!    
!    real(8)::density
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1)::Ftest
!    
!    real(8),dimension(nds,nds)::jac
!    
!    !Geometric Nonlinearity
!    integer(4)::inodes,ilixo
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1)::tdisp,mdisp,tzero
!    real(8),dimension(nds,nds)::rconv,rconv_mid,jacinv,raux
!    real(8),dimension(nds,nds)::jac_mid,jacinv_mid
!    real(8),dimension(nds,nds)::matF,matF_mid,matR,matR_mid
!    real(8),dimension((p+1)*(q+1)*(w+1))::dRde,dRdn,dRdc
!    real(8),dimension(nds)::temp1
!    real(8),dimension((p+1)*(q+1)*(w+1),nds)::nodes
!    real(8),dimension(nds,nds)::Temp33_mid,Temp33
!    real(8),dimension(nds,nds)::Trans_mid,Trans
!    real(8),dimension(nds*2,nds*2)::TCL_mid,TGL_mid
!    real(8),dimension(nds*2,nds*2)::TCL,TGL
!    real(8),dimension(9,(p+1)*(q+1)*(w+1)*nds)::BNL
!    real(8),dimension(9,9)::MSTR
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds)::KNLG
!    real(8),dimension(nds,1)::dRdenc,dRdxyz
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob_mid,Bloc_mid
!    
!    integer(4)::ii,jj,kk,k2,k3,loc_num
!    real(8)::CA,CB,detjA
!    real(8),dimension(nds,nds)::jacA,dxdxiA,dxdxi
!    real(8),dimension(npi,nds)::NatPoints
!    real(8),dimension(6)::ShpNat
!    real(8),dimension(npi)::Ra
!    real(8),dimension(npi,3)::DerShpNat,dRadx,dRAdxi,dRAdxii
!    real(8),dimension(npi,3)::dRdxi,dRdxii
!    real(8),dimension(24,nds)::TyPt
!    real(8)::xib,etab,zetab
!    
!    real(8),dimension(p+1)::Nb
!    real(8),dimension(q+1)::Mb
!    real(8),dimension(w+1)::Lb
!    real(8),dimension(p+1)::dNbdxi
!    real(8),dimension(q+1)::dMbdeta
!    real(8),dimension(w+1)::dLbdzeta
!    
!    real(8),dimension(2*(p+1))::ub
!    real(8),dimension(2*(q+1))::vb
!    
!    real(8),dimension(2*p)::ubr
!    real(8),dimension(2*q)::vbr
!    
!    real(8),dimension(p)::Nbr
!    real(8),dimension(q)::Mbr
!    real(8),dimension(p)::dNbrdxi
!    real(8),dimension(q)::dMbrdeta
!    
!    real(8),dimension(nds*2,nds*2)::TGC,TGC2
!    
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BNat,BANS,BNG,BNGold
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BConv
!    
!    real(8)::shp,suma,uk,vk,wk,wei,sumw,vol,A,B
!    
!    real(8),dimension(16)::shpANS
!    
!    real(8),dimension(3,3)::jact
!    
!    real(8),dimension(6,6)::MTP
!    real(8),dimension(4,4)::MTPr
!    
!    real(8),dimension(1,6)::M,MMult
!    
!    real(8),dimension(1,4)::Mr,MMultr
!    real(8),dimension(4)::temp4
!    real(8),dimension(6)::temp6
!    
!    logical::lsp,endcycle
!    logical::com1,com2,com3,com4,com5,com6
!    
!    integer(4),parameter:: ialpha = 3
!    real(8)::yc,zc,dNBde,dNBdn,dNBdc,detj0
!    real(8),dimension(6,ialpha)::BEAS,MEAS
!    real(8),dimension(6,6)::T0
!    real(8),dimension(ialpha)::tempa
!    real(8),dimension(ialpha,ialpha)::Kaa
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,ialpha)::Kua
!    
!    ShpNat = 0.0d0
!    derShpNat = 0.0d0
!    
!    suma = 0.0d0
!    
!    inodes = (p+1)*(q+1)*(w+1)
!    density = props(3)
!    
!    tzero = 0.0d0
!    
!    cpi = 0
!
!    we = 1.0d0
!    wn = 1.0d0
!    wc = 1.0d0
!    call gauleg(npi_xi, e, we)
!    call gauleg(npi_eta, n, wn)
!    call gauleg(npi_zeta, c, wc)
!    
!    call ShapeFunc3(nel,0.0d0,0.0d0,0.0d0,R,dRdx,dRdxi,dRdxii,detj0,jac,dxdxi,tzero)
!    temp1 = 0.0d0
!    jacinv = 0.0d0
!    jacinv = jac
!    call gaussj(jacinv,nds,temp1,ilixo)
!    
!    call TransformationMat3D(JacInv,T0)
!    
!    !#####
!    if(nlgeom==.true.)then
!        tdisp = el_ddisp
!        mdisp = el_ddisp/2.0d0
!    else
!        tdisp = 0.0d0
!        mdisp= 0.0d0
!        
!        TGL = 0.0d0
!        do i=1,6
!            TGL(i,i) = 1.0d0
!        end do
!    end if
!    !#####
!    
!    !Gauss points cycle
!    Ke = 0.0d0
!    Kua = 0.0d0
!    Kaa = 0.0d0
!    KNLG = 0.0d0
!    BNL = 0.0d0
!    MSTR = 0.0d0
!    
!    vol = 0.0d0
!    
!    Finte = 0.0d0
!    do k=1,npi_zeta
!        do j=1,npi_eta
!            do i=1,npi_xi
!                
!                !suma = 0.0d0
!                
!                cpi = cpi + 1
!                MatD(:,:) = TmatD(nel,cpi,:,:)
!                
!                xi = e(i)
!                eta = n(j)
!                zeta = c(k)
!                
!                dNBde=-1.0d0*xi*(1.0d0-eta**2.0d0)*(1.0d0-zeta**2.0d0)
!                dNBdn=-1.0d0*eta*(1.0d0-xi**2.0d0)*(1.0d0-zeta**2.0d0)
!                dNBdc=-1.0d0*zeta*(1.0d0-xi**2.0d0)*(1.0d0-eta**2.0d0)
!                
!                BEAS = 0.0d0
!                MEAS = 0.0d0
!                
!                call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tzero)
!                
!                !Weight factor
!                gwt = we(i)*wn(j)*wc(k)*detj
!                
!                
!                !#########
!                if(nlgeom==.true.)then
!                    
!                    !Local Axis ---------------------------
!                    if(inc==1)then
!                        call local_axis(nds,jac,rconv)
!                    else
!                        rconv = laxis(nel,cpi,:,:)
!                    end if
!                    
!                    !-------------------------------------------------------------------------
!                    !Mid-Configuration
!                    !-------------------------------------------------------------------------
!                    matF_mid = 0.0d0
!                    do k1=1,inodes
!                        dRde(k1) = dRdx(k1,1)
!                        dRdn(k1) = dRdx(k1,2)
!                        dRdc(k1) = dRdx(k1,3)
!                    end do
!                    
!                    jacinv = 0.0d0
!                    jacinv = jac
!                    call gaussj(jacinv,nds,temp1,ilixo)
!                    
!                    count = 0
!                    do k1=1,inodes
!                        count = count + 1
!                        nodes(k1,1) = Points(IEN(nel,count),1)
!                        nodes(k1,2) = Points(IEN(nel,count),2)
!                        nodes(k1,3) = Points(IEN(nel,count),3)
!                    end do
!                    
!                    !Deformation Gradient for mid-configuration
!                    call DefGrad3D(inodes,nds,dRde,dRdn,dRdc,jacinv,nodes,mdisp,MatF_mid)
!                    
!                    !Polar Decomposition for mid configuration
!                    matR_mid=0.0d0
!                    call PolarDecomp2D(nds,matF_mid,matR_mid)
!                    
!                    !Update local reference system to the mid configuration
!                    call axisupdate(rconv,matR_mid,rconv_mid)
!                    
!                    !Global to natural jacobian in mid configuration
!                    call ShapeFunc2(nel,xi,eta,zeta,R_mid,dRdx_mid,detj_mid,jac_mid,mdisp)
!                    
!                    jacinv_mid = 0.0d0
!                    temp1 = 0.0d0
!                    jacinv_mid = jac_mid
!                    call gaussj(jacinv_mid,nds,temp1,ilixo)
!                    
!                    !Transformation matrices mid configuration 
!                    Temp33_mid=0.0d0
!                    Temp33_mid= matmul(transpose(rconv_mid),jacinv_mid)
!                    Trans_mid=0.0d0
!                    Trans_mid=transpose(rconv_mid)
!                    
!                    !Natural to local transformation matrix for mid configuration
!                    call TransformationMat3D(Temp33_mid,TCL_mid)
!                    
!                    !Global to local transformation matrix for mid configuration
!                    call TransformationMat3D(Trans_mid,TGL_mid)
!                    
!                    Bglob_mid = 0.0d0
!                    do k1=1,(p+1)*(q+1)*(w+1)
!                        Bglob_mid(1,k1*3-2) = dRdx_mid(k1,1)
!                        Bglob_mid(2,k1*3-1) = dRdx_mid(k1,2)
!                        Bglob_mid(3,k1*3  ) = dRdx_mid(k1,3)
!                        Bglob_mid(4,k1*3-2) = dRdx_mid(k1,2)
!                        Bglob_mid(4,k1*3-1) = dRdx_mid(k1,1)
!                        Bglob_mid(5,k1*3-2) = dRdx_mid(k1,3)
!                        Bglob_mid(5,k1*3  ) = dRdx_mid(k1,1)
!                        Bglob_mid(6,k1*3-1) = dRdx_mid(k1,3)
!                        Bglob_mid(6,k1*3  ) = dRdx_mid(k1,2)
!                    end do
!                    
!                    Bloc_mid=0.0d0
!                    Bloc_mid=matmul(TGL_mid,Bglob_mid)
!                    
!                    !-------------------------------------------------------------------------
!                    !End-Configuration
!                    !-------------------------------------------------------------------------
!                    !Deformation Gradient for mid-configuration
!                    call DefGrad3D(inodes,nds,dRde,dRdn,dRdc,jacinv,nodes,tdisp,MatF)
!                    
!                    !Polar Decomposition for mid configuration
!                    matR = 0.0d0
!                    call PolarDecomp2D(nds,matF,matR)
!                    
!                    !Update local reference system to the mid configuration
!                    raux = rconv
!                    call axisupdate(rconv,matR,rconv)
!                    
!                    !Global to natural jacobian in mid configuration
!                    call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tdisp)
!                    
!                    jacinv = 0.0d0
!                    temp1 = 0.0d0
!                    jacinv = jac
!                    call gaussj(jacinv,nds,temp1,ilixo)
!                    
!                    !Transformation matrices mid configuration 
!                    Temp33 = 0.0d0
!                    Temp33 = matmul(transpose(rconv),jacinv)
!                    Trans = 0.0d0
!                    Trans = transpose(rconv)
!                    
!                    !Natural to local transformation matrix for mid configuration
!                    call TransformationMat3D(Temp33,TCL)
!                    
!                    !Global to local transformation matrix for mid configuration
!                    call TransformationMat3D(Trans,TGL)
!                    
!                    !Update Weight factor
!                    gwt = we(i)*wn(j)*wc(k)*detj
!                    
!                    continue
!                    
!                end if
!                !############
!                
!                call ShapeFunc3(nel,xi,eta,zeta,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,tzero)
!                
!                jac = transpose(jac)
!                
!                jacinv = 0.0d0
!                temp1 = 0.0d0
!                jacinv = jac
!                call gaussj(jacinv,nds,temp1,ilixo)
!                
!                !###################################################################
!                !Corrigir no caso de nlgeom
!                call TransformationMat3D(jac,TGC)
!                call TransformationMat3D(jacinv,TCL)
!                !###################################################################
!                
!                continue
!                
!                Bglob = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
!                    
!                    Bglob(1,k1*3-2) = dRdx(k1,1)
!                    Bglob(2,k1*3-1) = dRdx(k1,2)
!                    Bglob(3,k1*3  ) = dRdx(k1,3)
!                    
!                    Bglob(4,k1*3-2) = dRdx(k1,2)
!                    Bglob(4,k1*3-1) = dRdx(k1,1)
!                    Bglob(5,k1*3-2) = dRdx(k1,3)
!                    Bglob(5,k1*3  ) = dRdx(k1,1)
!                    Bglob(6,k1*3-1) = dRdx(k1,3)
!                    Bglob(6,k1*3  ) = dRdx(k1,2)
!                    
!                end do
!                
!!                BConv = 0.0d0
!!                do k1=1,(p+1)*(q+1)*(w+1)
!!                    do k2 =1,nds
!!                        BConv(1,(k1-1)*3+k2) = dRdxii(k1,1)*jac(1,k2)
!!                        BConv(2,(k1-1)*3+k2) = dRdxii(k1,2)*jac(2,k2)
!!                        BConv(3,(k1-1)*3+k2) = dRdxii(k1,3)*jac(3,k2)
!!                        
!!                        BConv(4,(k1-1)*3+k2) = dRdxii(k1,1)*jac(2,k2) + dRdxii(k1,2)*jac(1,k2)
!!                        BConv(5,(k1-1)*3+k2) = dRdxii(k1,1)*jac(3,k2) + dRdxii(k1,3)*jac(1,k2)
!!                        BConv(6,(k1-1)*3+k2) = dRdxii(k1,2)*jac(3,k2) + dRdxii(k1,3)*jac(2,k2)
!!                    end do
!!                end do
!                
!                !Local strain-displacement matrix at end configuration
!                Bloc=0.0d0
!                
!                !BConv=matmul(TGC,Bglob)
!                !BNG = matmul(TGC,Bglob)
!                
!                BNG = BConv
!                
!                !BNGold = BNG
!                
!                !Let's Try the ANS ideia ---------------------------------------------------------------
!                CA = dsqrt(1.0d0/3.0d0)
!                CB = dsqrt(3.0d0/5.0d0)
!                
!                !Tying Points Coordinates
!                TyPt = zeta
!                TyPt(1,1)  =    CA; TyPt(1,2)  =    CB;
!                TyPt(2,1)  =   -CA; TyPt(2,2)  =    CB;
!                TyPt(3,1)  =    CA; TyPt(3,2)  = 0.0d0;
!                TyPt(4,1)  =   -CA; TyPt(4,2)  = 0.0d0;
!                TyPt(5,1)  =    CA; TyPt(5,2)  =   -CB;
!                TyPt(6,1)  =   -CA; TyPt(6,2)  =   -CB;
!                
!                TyPt(7,1)  =    CB; TyPt(7,2)  =    CA;
!                TyPt(8,1)  = 0.0d0; TyPt(8,2)  =    CA;
!                TyPt(9,1)  =   -CB; TyPt(9,2)  =    CA;
!                TyPt(10,1) =    CB; TyPt(10,2) =   -CA;
!                TyPt(11,1) = 0.0d0; TyPt(11,2) =   -CA;
!                TyPt(12,1) =   -CB; TyPt(12,2) =   -CA;
!                
!                TyPt(13,1) =  CA; TyPt(13,2) =  CA; TyPt(13,3) =  zeta
!                TyPt(14,1) =  CA; TyPt(14,2) = -CA; TyPt(14,3) =  zeta
!                TyPt(15,1) = -CA; TyPt(15,2) =  CA; TyPt(15,3) =  zeta
!                TyPt(16,1) = -CA; TyPt(16,2) = -CA; TyPt(16,3) =  zeta
!                
!!                TyPt(17,1) =     CA; TyPt(17,2) =  1.0d0; TyPt(17,3) =  zeta
!!                TyPt(18,1) =    -CA; TyPt(18,2) =  1.0d0; TyPt(18,3) =  zeta
!!                TyPt(19,1) =     CA; TyPt(19,2) = -1.0d0; TyPt(19,3) =  zeta
!!                TyPt(20,1) =    -CA; TyPt(20,2) = -1.0d0; TyPt(20,3) =  zeta
!!                TyPt(21,1) =  1.0d0; TyPt(21,2) =     CA; TyPt(21,3) =  zeta
!!                TyPt(22,1) =  1.0d0; TyPt(22,2) =    -CA; TyPt(22,3) =  zeta
!!                TyPt(23,1) = -1.0d0; TyPt(23,2) =     CA; TyPt(23,3) =  zeta
!!                TyPt(24,1) = -1.0d0; TyPt(24,2) =    -CA; TyPt(24,3) =  zeta
!                
!                com1 = .true.
!                com2 = .true.
!                !com3 = .false.
!                com4 = .true.
!                com5 = .true.
!                com6 = .true.
!                
!!                if(com1 == .true.) BNG(1,:) = 0.0d0
!!                if(com2 == .true.) BNG(2,:) = 0.0d0
!!                !if(com3 == .true.) BNG(3,:) = 0.0d0
!!                if(com4 == .true.) BNG(4,:) = 0.0d0
!!                if(com5 == .true.) BNG(5,:) = 0.0d0
!!                if(com6 == .true.) BNG(6,:) = 0.0d0
!                
!                ni = INN(IEN(nel,1),1)
!                nj = INN(IEN(nel,1),2)
!                nk = INN(IEN(nel,1),3)
!    
!                xib   = ((u_knot(ni+1) - u_knot(ni))*xi   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                etab  = ((v_knot(nj+1) - v_knot(nj))*eta  + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                zetab = ((w_knot(nk+1) - w_knot(nk))*zeta + (w_knot(nk+1) + w_knot(nk)))/2.0d0
!                
!                sumw = 0.0d0
!                
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + R(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !-------------------------------------------------------------------------------------------
!
!                !Bezier knot vectors -----------------------------------------------
!                ub = 1.0d0
!                vb = 1.0d0
!                ubr = 1.0d0
!                vbr = 1.0d0
!                
!                do k1=1,p+1
!                    ub(k1) = 0.0d0
!                end do
!                do k1=1,p
!                    ubr(k1) = 0.0d0
!                end do
!                do k1=1,q+1
!                    vb(k1) = 0.0d0
!                end do
!                do k1=1,q
!                    vbr(k1) = 0.0d0
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)
!                   
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                
!                
!                !First TP set cycle ------------------------------------------------
!                !Mass Matrix -------------------------------------------------------
!                MTP = 0.0d0
!                do k1=1,6
!                    uk = (TyPt(k1,1) + 1.0d0)/2.0d0
!                    vk = (TyPt(k1,2) + 1.0d0)/2.0d0
!                    
!                    Nbr = 0.0d0
!                    call BSplineBasisAndDeriv(ncpx+p,p-1,uk,ubr,Nbr,dNbrdxi)
!                    Mb = 0.0d0
!                    call BSplineBasisAndDeriv(ncpy+q,q  ,vk, vb, Mb,dMbdeta)
!                    
!                    count = 0
!                    do k2=0,q
!                        do k3=0,p-1
!                            count = count + 1
!                            MTP(k1,count) = Nbr(p-k3)*Mb(q+1-k2)
!                        end do
!                    end do      
!                end do
!                
!                !IP Mass Matrix ----------------------------------------------------
!                xib  = ( xi + 1.0d0)/2.0d0
!                etab = (eta + 1.0d0)/2.0d0
!                
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q  ,etab, vb, Mb,dMbdeta)
!                    
!                count = 0
!                do k2=0,q
!                    do k3=0,p-1
!                        count = count + 1
!                        M(1,count) = Nbr(p-k3)*Mb(q+1-k2)
!                    end do
!                end do 
!                
!                temp6 = 0.0d0
!                call gaussj(MTP,6,temp6,ilixo)
!                
!                MMult = matmul(M,MTP)
!                
!                do k1=1,6
!                    
!                    dNBde=-1.0d0*TyPt(k1,1)*(1.0d0-TyPt(k1,2)**2.0d0)*(1.0d0-TyPt(k1,3)**2.0d0)
!                    dNBdn=-1.0d0*TyPt(k1,2)*(1.0d0-TyPt(k1,1)**2.0d0)*(1.0d0-TyPt(k1,3)**2.0d0)
!                    dNBdc=-1.0d0*TyPt(k1,3)*(1.0d0-TyPt(k1,1)**2.0d0)*(1.0d0-TyPt(k1,2)**2.0d0)
!                    
!                    MEAS = 0.0d0
!                    MEAS(1,1) = dNBde
!                    MEAS(2,2) = dNBdn
!                    MEAS(3,3) = dNBdc
!                    MEAS(4,1) = dNBdn
!                    MEAS(4,2) = dNBde
!                    MEAS(5,1) = dNBdc
!                    MEAS(5,3) = dNBde
!                    MEAS(6,2) = dNBdc
!                    MEAS(6,3) = dNBdn
!                    
!                    BEAS = BEAS + matmul(T0,MEAS)*detj0/detj*MMult(1,k1)
!                    
!                end do
!                
!                !Second TP set cycle -----------------------------------------------
!                !Mass Matrix -------------------------------------------------------
!                MTP = 0.0d0
!                do k1=1,6
!                    uk = (TyPt(k1+6,1) + 1.0d0)/2.0d0
!                    vk = (TyPt(k1+6,2) + 1.0d0)/2.0d0
!                    
!                    Nb = 0.0d0
!                    call BSplineBasisAndDeriv(ncpx+p,p  ,uk, ub,  Nb,  dNbdxi)
!                    Mbr = 0.0d0
!                    call BSplineBasisAndDeriv(ncpy+q,q-1,vk,vbr, Mbr,dMbrdeta)
!                    
!                    count = 0
!                    do k2=0,q-1
!                        do k3=0,p
!                            count = count + 1
!                            MTP(k1,count) = Nb(p+1-k3)*Mbr(q-k2)
!                        end do
!                    end do      
!                end do
!                
!                temp6 = 0.0d0
!                call gaussj(MTP,6,temp6,ilixo)
!                continue
!                
!                !IP Mass Matrix ----------------------------------------------------
!                xib  = ( xi + 1.0d0)/2.0d0
!                etab = (eta + 1.0d0)/2.0d0
!                
!                call BSplineBasisAndDeriv(ncpx+p,p  , xib, ub, Nb,  dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                    
!                count = 0
!                M = 0.0d0
!                do k2=0,q-1
!                    do k3=0,p
!                        count = count + 1
!                        M(1,count) = Nb(p+1-k3)*Mbr(q-k2)
!                    end do
!                end do 
!                
!                MMult = 0.0d0
!                MMult = matmul(M,MTP)
!                continue
!                
!                do k1=1,6
!                    dNBde=-1.0d0*TyPt(k1+6,1)*(1.0d0-TyPt(k1+6,2)**2.0d0)*(1.0d0-TyPt(k1+6,3)**2.0d0)
!                    dNBdn=-1.0d0*TyPt(k1+6,2)*(1.0d0-TyPt(k1+6,1)**2.0d0)*(1.0d0-TyPt(k1+6,3)**2.0d0)
!                    dNBdc=-1.0d0*TyPt(k1+6,3)*(1.0d0-TyPt(k1+6,1)**2.0d0)*(1.0d0-TyPt(k1+6,2)**2.0d0)
!                    
!                    MEAS = 0.0d0
!                    MEAS(1,1) = dNBde
!                    MEAS(2,2) = dNBdn
!                    MEAS(3,3) = dNBdc
!                    MEAS(4,1) = dNBdn
!                    MEAS(4,2) = dNBde
!                    MEAS(5,1) = dNBdc
!                    MEAS(5,3) = dNBde
!                    MEAS(6,2) = dNBdc
!                    MEAS(6,3) = dNBdn
!                    
!                    BEAS = BEAS + matmul(T0,MEAS)*detj0/detj*MMult(1,k1)
!                end do
!                
!                
!                !Third TP set cycle ------------------------------------------------
!                !Mass Matrix -------------------------------------------------------
!                MTPr = 0.0d0
!                do k1=1,4
!                    uk = (TyPt(k1+12,1) + 1.0d0)/2.0d0
!                    vk = (TyPt(k1+12,2) + 1.0d0)/2.0d0
!                    
!                    Nbr = 0.0d0
!                    call BSplineBasisAndDeriv(ncpx+p,p-1,uk,ubr, Nbr, dNbrdxi)
!                    Mbr = 0.0d0
!                    call BSplineBasisAndDeriv(ncpy+q,q-1,vk,vbr, Mbr,dMbrdeta)
!                    
!                    count = 0
!                    do k2=0,q-1
!                        do k3=0,p-1
!                            count = count + 1
!                            MTPr(k1,count) = Nbr(p-k3)*Mbr(q-k2)
!                        end do
!                    end do      
!                end do
!                
!                temp4 = 0.0d0
!                call gaussj(MTPr,4,temp4,ilixo)
!                
!                
!                !IP Mass Matrix ----------------------------------------------------
!                xib  = ( xi + 1.0d0)/2.0d0
!                etab = (eta + 1.0d0)/2.0d0
!                
!                call BSplineBasisAndDeriv(ncpx+p,p-1, xib,ubr,Nbr, dNbrdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                    
!                count = 0
!                Mr = 0.0d0
!                do k2=0,q-1
!                    do k3=0,p-1
!                        count = count + 1
!                        Mr(1,count) = Nbr(p-k3)*Mbr(q-k2)
!                    end do
!                end do 
!                
!                
!                MMultr = matmul(Mr,MTPr)
!                continue
!                
!                do k1=1,4
!                    dNBde=-1.0d0*TyPt(k1+12,1)*(1.0d0-TyPt(k1+12,2)**2.0d0)*(1.0d0-TyPt(k1+12,3)**2.0d0)
!                    dNBdn=-1.0d0*TyPt(k1+12,2)*(1.0d0-TyPt(k1+12,1)**2.0d0)*(1.0d0-TyPt(k1+12,3)**2.0d0)
!                    dNBdc=-1.0d0*TyPt(k1+12,3)*(1.0d0-TyPt(k1+12,1)**2.0d0)*(1.0d0-TyPt(k1+12,2)**2.0d0)
!                    
!                    MEAS = 0.0d0
!                    MEAS(1,1) = dNBde
!                    MEAS(2,2) = dNBdn
!                    MEAS(3,3) = dNBdc
!                    MEAS(4,1) = dNBdn
!                    MEAS(4,2) = dNBde
!                    MEAS(5,1) = dNBdc
!                    MEAS(5,3) = dNBde
!                    MEAS(6,2) = dNBdc
!                    MEAS(6,3) = dNBdn
!                    
!                    BEAS = BEAS + matmul(T0,MEAS)*detj0/detj*MMult(1,k1)
!                end do
!                !-------------------------------------------------------------------------------------------
!
!                ! Strain-Displacement Matrix in local frame
!                Bloc = BGlob
!
!                if(nlgeom == .false.) Bloc_mid = Bloc
!                
!                call MatPlastic3D(iprops,nds,(p+1)*(q+1)*(w+1),props,Bloc_mid,6,el_ddisp,stress(nel,cpi,:),strain(nel,cpi,:),dstrain(nel,cpi,:),dstress(nel,cpi,:),hard(nel,cpi),matD)
!                
!                !Store Constitutive Matrix
!                TmatD(nel,cpi,:,:) = MatD(:,:)
!                
!                tempstr = 0.0d0
!                do k1=1,ntens
!                    tempstr(k1,1) = stress(nel,cpi,k1)
!                    deform(k1,1) = dstrain(nel,cpi,k1)
!                end do
!                
!                if (nlgeom == .true.) then
!                
!                    BNL = 0.0d0
!                    do k1=1,inodes
!                        dRdenc=0.0d0
!                        dRdenc(1,1)=dRdx(k1,1) !
!                        dRdenc(2,1)=dRdx(k1,2) !
!                        dRdenc(3,1)=dRdx(k1,3) !
!                        
!                        dRdxyz = 0.0d0
!                        dRdxyz = matmul(trans,dRdenc)
!
!                        BNL(1,k1*3-2)=dRdxyz(1,1)
!                        BNL(2,k1*3-2)=dRdxyz(2,1)
!                        BNL(3,k1*3-2)=dRdxyz(3,1)
!                        BNL(4,k1*3-1)=dRdxyz(1,1)
!                        BNL(5,k1*3-1)=dRdxyz(2,1)
!                        BNL(6,k1*3-1)=dRdxyz(3,1)
!                        BNL(7,k1*3  )=dRdxyz(1,1)
!                        BNL(8,k1*3  )=dRdxyz(2,1)
!                        BNL(9,k1*3  )=dRdxyz(3,1)
!                    end do 
!                 
!                    MSTR=0.0d0
!                    do k1=1,nds
!                        MSTR(k1*3-2,k1*3-2) = tempstr(1,1)
!                        MSTR(k1*3-2,k1*3-1) = tempstr(4,1)
!                        MSTR(k1*3-2,k1*3  ) = tempstr(5,1)
!                        MSTR(k1*3-1,k1*3-2) = tempstr(4,1)
!                        MSTR(k1*3-1,k1*3-1) = tempstr(2,1)
!                        MSTR(k1*3-1,k1*3  ) = tempstr(6,1)
!                        MSTR(k1*3  ,k1*3-2) = tempstr(5,1)
!                        MSTR(k1*3  ,k1*3-1) = tempstr(6,1)
!                        MSTR(k1*3  ,k1*3  ) = tempstr(3,1)
!                    end do
!                    
!                    KNLG = KNLG + matmul(matmul(transpose(BNL),MSTR),BNL)*gwt
!                    
!                    !Export local axis ---------------
!                    laxis(nel,cpi,:,:) = rconv
! 
!                end if
!                
!                if(nlgeom == .false.) then
!                    Ke = Ke + matmul(matmul(transpose(Bloc),matD),Bloc)*gwt
!                else
!                    Ke = Ke + matmul(matmul(transpose(Bloc),matD),Bloc)*gwt + matmul(matmul(transpose(BNL),MSTR),BNL)*gwt
!                end if
!                
!                
!                Kua = Kua + matmul(matmul(transpose(Bloc),matD),BEAS)*gwt
!                Kaa  =Kaa + matmul(matmul(transpose(BEAS),matD),BEAS)*gwt
!                
!                
!                
!                
!                Finte = Finte + matmul(transpose(Bloc),tempstr)*gwt
!                
!                
!                !Gravitic loads -------------------------------------------------------------
!                if (gravity==.true.)then
!                    call grvt(nshpl,nds,R,gwt,gconst,gravdir,density,Finte)
!                end if
!                
!                !Pressure loads -------------------------------------------------------------
!                if(pressure==.true.) then
!                    do k1=1,nelems
!                        if(k1==nel .and. abs(pressvec(k1)) .gt. 0.0d0) then
!                            call prssr(nel,npi_xi,npi_eta,npi_zeta,jac,we(i),wn(j),wc(k),R,Finte)
!                        end if
!                    end do
!                end if
!                
!                !Strain Energy --------------------------------------------------------------
!                do k1=1,ntens
!                    SEnergy = SEnergy + tempstr(k1,1)*deform(k1,1)*gwt/2.0d0
!                end do
!                
!                nj = INN(IEN(nel,1),2)
!    
!                yc  = ((v_knot(nj+1) - v_knot(nj))*eta + (v_knot(nj+1) + v_knot(nj)))/2.0d0*10.0
!                zc = ((w_knot(nk+1) - w_knot(nk))*zeta + (w_knot(nk+1) + w_knot(nk)))/2.0d0
!                
!                continue
!                
!            end do  
!        end do
!    end do
!    
!    tempa=0.0d0
!    call gaussj(kaa,ialpha,tempa,ilixo)
!    
!    Ke = Ke  - matmul(matmul(Kua,Kaa),transpose(Kua))
!    
!    
!    
!    
!    
!    continue
!
!end subroutine ElHex27EAS






!------------------------------------------------------------------------------------------------------------
!
! Second order hexahedral element with EAS
!
! Input: nel - Element number
!        el_ddisp - Element displacement
!        inc - Increment number
!
! Output: Ke - Elementar stiffness matrix
!         Finte - Elementar internal forces
!------------------------------------------------------------------------------------------------------------

subroutine ElHex27EAS(nel,Ke,el_ddisp,Finte,inc)
    
    use Mod_Variables
    implicit none
    
    integer(4)::nel,cpi,count,inc
    
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::el_ddisp
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(OUT)::Ke
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(OUT)::Finte
    
    integer(4),parameter::npi_xi= 3, npi_eta = 3, npi_zeta = 3
    
    integer(4)::i,j,k,ni,nj,nk,nel_nza,k1
    real(8)::xi,eta,zeta
    real(8),dimension(npi_xi)::e,we
    real(8),dimension(npi_eta)::n,wn
    real(8),dimension(npi_zeta)::c,wc
    
    real(8),dimension((p+1)*(q+1)*(w+1))::R,R_mid
    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdx,dRdx_mid

    real(8)::detJ,detj_mid,gwt
    
    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob,Bloc
    
    real(8),dimension(6,6)::matD
    real(8),dimension(ntens,1)::tempstr,deform
    
    real(8)::density
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1)::Ftest
    
    real(8),dimension(nds,nds)::jac
    
    !Geometric Nonlinearity
    integer(4)::inodes,ilixo
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1)::tdisp,mdisp,tzero
    real(8),dimension(nds,nds)::rconv,rconv_mid,jacinv,raux
    real(8),dimension(nds,nds)::jac_mid,jacinv_mid
    real(8),dimension(nds,nds)::matF,matF_mid,matR,matR_mid
    real(8),dimension((p+1)*(q+1)*(w+1))::dRde,dRdn,dRdc
    real(8),dimension(nds)::temp1
    real(8),dimension((p+1)*(q+1)*(w+1),nds)::nodes
    real(8),dimension(nds,nds)::Temp33_mid,Temp33
    real(8),dimension(nds,nds)::Trans_mid,Trans
    real(8),dimension(nds*2,nds*2)::TCL_mid,TGL_mid
    real(8),dimension(nds*2,nds*2)::TCL,TGL
    real(8),dimension(9,(p+1)*(q+1)*(w+1)*nds)::BNL
    real(8),dimension(9,9)::MSTR
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds)::KNLG
    real(8),dimension(nds,1)::dRdenc,dRdxyz
    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob_mid,Bloc_mid
    
    integer(4)::ii,jj,kk,k2,k3,loc_num
    real(8)::CA,CB,detjA
    real(8),dimension(nds,nds)::jacA,dxdxiA,dxdxi
    real(8),dimension(npi,nds)::NatPoints
    real(8),dimension(6)::ShpNat
    real(8),dimension(npi)::Ra
    real(8),dimension(npi,3)::DerShpNat,dRadx,dRAdxi,dRAdxii
    real(8),dimension(npi,3)::dRdxi,dRdxii
    real(8),dimension(24,nds)::TyPt
    real(8)::xib,etab,zetab
    
    real(8),dimension(p+1)::Nb
    real(8),dimension(q+1)::Mb
    real(8),dimension(w+1)::Lb
    real(8),dimension(p+1)::dNbdxi
    real(8),dimension(q+1)::dMbdeta
    real(8),dimension(w+1)::dLbdzeta
    
    real(8),dimension(2*(p+1))::ub
    real(8),dimension(2*(q+1))::vb
    
    real(8),dimension(2*p)::ubr
    real(8),dimension(2*q)::vbr
    
    real(8),dimension(p)::Nbr
    real(8),dimension(q)::Mbr
    real(8),dimension(p)::dNbrdxi
    real(8),dimension(q)::dMbrdeta
    
    real(8),dimension(nds*2,nds*2)::TGC,TGC2
    
    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BNat,BANS,BNG,BNGold
    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BConv,Bdil,Bdev
    
    real(8)::shp,suma,uk,vk,wk,wei,sumw,vol,A,B
    
    real(8),dimension(16)::shpANS
    
    real(8),dimension(3,3)::jact
    
    real(8),dimension(6,6)::MTP
    real(8),dimension(4,4)::MTPr
    
    real(8),dimension(1,6)::M,MMult
    
    real(8),dimension(1,4)::Mr,MMultr
    real(8),dimension(4)::temp4
    real(8),dimension(6)::temp6
    
    logical::lsp,endcycle
    logical::com1,com2,com3,com4,com5,com6
    
    integer(4),parameter:: ialpha = 3
    real(8)::yc,zc,dNBde,dNBdn,dNBdc,detj0
    real(8),dimension(6,ialpha)::BEAS,MEAS
    real(8),dimension(6,6)::T0
    real(8),dimension(ialpha)::tempa
    real(8),dimension(ialpha,ialpha)::Kaa
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,ialpha)::Kua
    
    ShpNat = 0.0d0
    derShpNat = 0.0d0
    
    suma = 0.0d0
    
    inodes = (p+1)*(q+1)*(w+1)
    density = props(3)
    
    tzero = 0.0d0
    
    cpi = 0

    we = 1.0d0
    wn = 1.0d0
    wc = 1.0d0
    call gauleg(npi_xi, e, we)
    call gauleg(npi_eta, n, wn)
    call gauleg(npi_zeta, c, wc)
    
    call ShapeFunc3(nel,0.0d0,0.0d0,0.0d0,R,dRdx,dRdxi,dRdxii,detj0,jac,dxdxi,tzero)
    temp1 = 0.0d0
    jacinv = 0.0d0
    jacinv = jac
    call gaussj(jacinv,nds,temp1,ilixo)
    
    call TransformationMat3D(JacInv,T0)
    
    !#####
    if(nlgeom==.true.)then
        tdisp = el_ddisp
        mdisp = el_ddisp/2.0d0
    else
        tdisp = 0.0d0
        mdisp= 0.0d0
        
        TGL = 0.0d0
        do i=1,6
            TGL(i,i) = 1.0d0
        end do
    end if
    !#####
    
    !Gauss points cycle
    Ke = 0.0d0
    Kua = 0.0d0
    Kaa = 0.0d0
    KNLG = 0.0d0
    BNL = 0.0d0
    MSTR = 0.0d0
    
    vol = 0.0d0
    
    Finte = 0.0d0
    do k=1,npi_zeta
        do j=1,npi_eta
            do i=1,npi_xi
                
                !suma = 0.0d0
                
                cpi = cpi + 1
                MatD(:,:) = TmatD(nel,cpi,:,:)
                
                xi = e(i)
                eta = n(j)
                zeta = c(k)
                
                dNBde=-1.0d0*xi*(1.0d0-eta**2.0d0)*(1.0d0-zeta**2.0d0)
                dNBdn=-1.0d0*eta*(1.0d0-xi**2.0d0)*(1.0d0-zeta**2.0d0)
                dNBdc=-1.0d0*zeta*(1.0d0-xi**2.0d0)*(1.0d0-eta**2.0d0)
                
                BEAS = 0.0d0
                MEAS = 0.0d0
                
                call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tzero)
                
                !Weight factor
                gwt = we(i)*wn(j)*wc(k)*detj
                
                
                !#########
                if(nlgeom==.true.)then
                    
                    !Local Axis ---------------------------
                    if(inc==1)then
                        call local_axis(nds,jac,rconv)
                    else
                        rconv = laxis(nel,cpi,:,:)
                    end if
                    
                    !-------------------------------------------------------------------------
                    !Mid-Configuration
                    !-------------------------------------------------------------------------
                    matF_mid = 0.0d0
                    do k1=1,inodes
                        dRde(k1) = dRdx(k1,1)
                        dRdn(k1) = dRdx(k1,2)
                        dRdc(k1) = dRdx(k1,3)
                    end do
                    
                    jacinv = 0.0d0
                    jacinv = jac
                    call gaussj(jacinv,nds,temp1,ilixo)
                    
                    count = 0
                    do k1=1,inodes
                        count = count + 1
                        nodes(k1,1) = Points(IEN(nel,count),1)
                        nodes(k1,2) = Points(IEN(nel,count),2)
                        nodes(k1,3) = Points(IEN(nel,count),3)
                    end do
                    
                    !Deformation Gradient for mid-configuration
                    call DefGrad3D(inodes,nds,dRde,dRdn,dRdc,jacinv,nodes,mdisp,MatF_mid)
                    
                    !Polar Decomposition for mid configuration
                    matR_mid=0.0d0
                    call PolarDecomp2D(nds,matF_mid,matR_mid)
                    
                    !Update local reference system to the mid configuration
                    call axisupdate(rconv,matR_mid,rconv_mid)
                    
                    !Global to natural jacobian in mid configuration
                    call ShapeFunc2(nel,xi,eta,zeta,R_mid,dRdx_mid,detj_mid,jac_mid,mdisp)
                    
                    jacinv_mid = 0.0d0
                    temp1 = 0.0d0
                    jacinv_mid = jac_mid
                    call gaussj(jacinv_mid,nds,temp1,ilixo)
                    
                    !Transformation matrices mid configuration 
                    Temp33_mid=0.0d0
                    Temp33_mid= matmul(transpose(rconv_mid),jacinv_mid)
                    Trans_mid=0.0d0
                    Trans_mid=transpose(rconv_mid)
                    
                    !Natural to local transformation matrix for mid configuration
                    call TransformationMat3D(Temp33_mid,TCL_mid)
                    
                    !Global to local transformation matrix for mid configuration
                    call TransformationMat3D(Trans_mid,TGL_mid)
                    
                    Bglob_mid = 0.0d0
                    do k1=1,(p+1)*(q+1)*(w+1)
                        Bglob_mid(1,k1*3-2) = dRdx_mid(k1,1)
                        Bglob_mid(2,k1*3-1) = dRdx_mid(k1,2)
                        Bglob_mid(3,k1*3  ) = dRdx_mid(k1,3)
                        Bglob_mid(4,k1*3-2) = dRdx_mid(k1,2)
                        Bglob_mid(4,k1*3-1) = dRdx_mid(k1,1)
                        Bglob_mid(5,k1*3-2) = dRdx_mid(k1,3)
                        Bglob_mid(5,k1*3  ) = dRdx_mid(k1,1)
                        Bglob_mid(6,k1*3-1) = dRdx_mid(k1,3)
                        Bglob_mid(6,k1*3  ) = dRdx_mid(k1,2)
                    end do
                    
                    Bloc_mid=0.0d0
                    Bloc_mid=matmul(TGL_mid,Bglob_mid)
                    
                    !-------------------------------------------------------------------------
                    !End-Configuration
                    !-------------------------------------------------------------------------
                    !Deformation Gradient for mid-configuration
                    call DefGrad3D(inodes,nds,dRde,dRdn,dRdc,jacinv,nodes,tdisp,MatF)
                    
                    !Polar Decomposition for mid configuration
                    matR = 0.0d0
                    call PolarDecomp2D(nds,matF,matR)
                    
                    !Update local reference system to the mid configuration
                    raux = rconv
                    call axisupdate(rconv,matR,rconv)
                    
                    !Global to natural jacobian in mid configuration
                    call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tdisp)
                    
                    jacinv = 0.0d0
                    temp1 = 0.0d0
                    jacinv = jac
                    call gaussj(jacinv,nds,temp1,ilixo)
                    
                    !Transformation matrices mid configuration 
                    Temp33 = 0.0d0
                    Temp33 = matmul(transpose(rconv),jacinv)
                    Trans = 0.0d0
                    Trans = transpose(rconv)
                    
                    !Natural to local transformation matrix for mid configuration
                    call TransformationMat3D(Temp33,TCL)
                    
                    !Global to local transformation matrix for mid configuration
                    call TransformationMat3D(Trans,TGL)
                    
                    !Update Weight factor
                    gwt = we(i)*wn(j)*wc(k)*detj
                    
                    continue
                    
                end if
                !############
                
                call ShapeFunc3(nel,xi,eta,zeta,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,tzero)
                
                jac = transpose(jac)
                
                jacinv = 0.0d0
                temp1 = 0.0d0
                jacinv = jac
                call gaussj(jacinv,nds,temp1,ilixo)
                
                !###################################################################
                !Corrigir no caso de nlgeom
                call TransformationMat3D(jac,TGC)
                call TransformationMat3D(jacinv,TCL)
                !###################################################################
                
                continue
                
                Bglob = 0.0d0
                Bdil = 0.0d0
                do k1=1,(p+1)*(q+1)*(w+1)
                    
                    Bglob(1,k1*3-2) = dRdx(k1,1)
                    Bglob(2,k1*3-1) = dRdx(k1,2)
                    Bglob(3,k1*3  ) = dRdx(k1,3)
                    
                    Bglob(4,k1*3-2) = dRdx(k1,2)
                    Bglob(4,k1*3-1) = dRdx(k1,1)
                    Bglob(5,k1*3-2) = dRdx(k1,3)
                    Bglob(5,k1*3  ) = dRdx(k1,1)
                    Bglob(6,k1*3-1) = dRdx(k1,3)
                    Bglob(6,k1*3  ) = dRdx(k1,2)
                    
                    Bdil(1,k1*3-2) = dRdx(k1,1)/3.0d0
                    Bdil(2,k1*3-2) = dRdx(k1,1)/3.0d0
                    Bdil(3,k1*3-2) = dRdx(k1,1)/3.0d0
                    Bdil(1,k1*3-1) = dRdx(k1,2)/3.0d0
                    Bdil(2,k1*3-1) = dRdx(k1,2)/3.0d0
                    Bdil(3,k1*3-1) = dRdx(k1,2)/3.0d0
                    Bdil(1,k1*3  ) = dRdx(k1,3)/3.0d0
                    Bdil(2,k1*3  ) = dRdx(k1,3)/3.0d0
                    Bdil(3,k1*3  ) = dRdx(k1,3)/3.0d0
                    
                end do
                
                Bdev = 0.0d0
                Bdev  = Bglob - Bdil

                BNG = 0.0d0

                
                !Let's Try the ANS ideia ---------------------------------------------------------------
                CA = dsqrt(1.0d0/3.0d0)
                CB = dsqrt(3.0d0/5.0d0)
                
                !Tying Points Coordinates
                TyPt = zeta
                TyPt(1,1)  =    CA; TyPt(1,2)  =    CB;
                TyPt(2,1)  =   -CA; TyPt(2,2)  =    CB;
                TyPt(3,1)  =    CA; TyPt(3,2)  = 0.0d0;
                TyPt(4,1)  =   -CA; TyPt(4,2)  = 0.0d0;
                TyPt(5,1)  =    CA; TyPt(5,2)  =   -CB;
                TyPt(6,1)  =   -CA; TyPt(6,2)  =   -CB;
                
                TyPt(7,1)  =    CB; TyPt(7,2)  =    CA;
                TyPt(8,1)  = 0.0d0; TyPt(8,2)  =    CA;
                TyPt(9,1)  =   -CB; TyPt(9,2)  =    CA;
                TyPt(10,1) =    CB; TyPt(10,2) =   -CA;
                TyPt(11,1) = 0.0d0; TyPt(11,2) =   -CA;
                TyPt(12,1) =   -CB; TyPt(12,2) =   -CA;
                
                TyPt(13,1) =  CA; TyPt(13,2) =  CA; TyPt(13,3) =  zeta
                TyPt(14,1) =  CA; TyPt(14,2) = -CA; TyPt(14,3) =  zeta
                TyPt(15,1) = -CA; TyPt(15,2) =  CA; TyPt(15,3) =  zeta
                TyPt(16,1) = -CA; TyPt(16,2) = -CA; TyPt(16,3) =  zeta
                
                ni = INN(IEN(nel,1),1)
                nj = INN(IEN(nel,1),2)
                nk = INN(IEN(nel,1),3)
    
                xib   = ((u_knot(ni+1) - u_knot(ni))*xi   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
                etab  = ((v_knot(nj+1) - v_knot(nj))*eta  + (v_knot(nj+1) + v_knot(nj)))/2.0d0
                zetab = ((w_knot(nk+1) - w_knot(nk))*zeta + (w_knot(nk+1) + w_knot(nk)))/2.0d0
                !-------------------------------------------------------------------------------------------

                !Bezier knot vectors -----------------------------------------------
                ub = 1.0d0
                vb = 1.0d0
                ubr = 1.0d0
                vbr = 1.0d0
                
                do k1=1,p+1
                    ub(k1) = 0.0d0
                end do
                do k1=1,p
                    ubr(k1) = 0.0d0
                end do
                do k1=1,q+1
                    vb(k1) = 0.0d0
                end do
                do k1=1,q
                    vbr(k1) = 0.0d0
                end do
                
                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)
                   
                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
                
                !Third TP set cycle ------------------------------------------------
                !Mass Matrix -------------------------------------------------------
                MTPr = 0.0d0
                do k1=1,4
                    uk = (TyPt(k1+12,1) + 1.0d0)/2.0d0
                    vk = (TyPt(k1+12,2) + 1.0d0)/2.0d0
                    
                    Nbr = 0.0d0
                    call BSplineBasisAndDeriv(ncpx+p,p-1,uk,ubr, Nbr, dNbrdxi)
                    Mbr = 0.0d0
                    call BSplineBasisAndDeriv(ncpy+q,q-1,vk,vbr, Mbr,dMbrdeta)
                    
                    count = 0
                    do k2=0,q-1
                        do k3=0,p-1
                            count = count + 1
                            MTPr(k1,count) = Nbr(p-k3)*Mbr(q-k2)
                        end do
                    end do      
                end do
                
                temp4 = 0.0d0
                call gaussj(MTPr,4,temp4,ilixo)
                
                
                !IP Mass Matrix ----------------------------------------------------
                xib  = ( xi + 1.0d0)/2.0d0
                etab = (eta + 1.0d0)/2.0d0
                
                call BSplineBasisAndDeriv(ncpx+p,p-1, xib,ubr,Nbr, dNbrdxi)
                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
                    
                count = 0
                Mr = 0.0d0
                do k2=0,q-1
                    do k3=0,p-1
                        count = count + 1
                        Mr(1,count) = Nbr(p-k3)*Mbr(q-k2)
                    end do
                end do 
                
                
                MMultr = matmul(Mr,MTPr)
                continue
                
                do k1=1,4
                    call ShapeFunc3(nel,TyPt(k1+12,1),TyPt(k1+12,2),TyPt(k1+12,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tdisp)
                    jacA = transpose(jacA)
                    
                    do k2=1,(p+1)*(q+1)*(w+1)
                        BNG(1,k2*3-2) = BNG(1,k2*3-2) + dRAdx(k2,1)*MMultr(1,k1)
                        BNG(2,k2*3-2) = BNG(2,k2*3-2) + dRAdx(k2,1)*MMultr(1,k1)
                        BNG(3,k2*3-2) = BNG(3,k2*3-2) + dRAdx(k2,1)*MMultr(1,k1)
                        BNG(1,k2*3-1) = BNG(1,k2*3-1) + dRAdx(k2,2)*MMultr(1,k1)
                        BNG(2,k2*3-1) = BNG(2,k2*3-1) + dRAdx(k2,2)*MMultr(1,k1)
                        BNG(3,k2*3-1) = BNG(3,k2*3-1) + dRAdx(k2,2)*MMultr(1,k1)
                        BNG(1,k2*3  ) = BNG(1,k2*3  ) + dRAdx(k2,3)*MMultr(1,k1)
                        BNG(2,k2*3  ) = BNG(2,k2*3  ) + dRAdx(k2,3)*MMultr(1,k1)
                        BNG(3,k2*3  ) = BNG(3,k2*3  ) + dRAdx(k2,3)*MMultr(1,k1)
                    end do
                end do
                !-------------------------------------------------------------------------------------------

                ! Strain-Displacement Matrix in local frame
                Bglob = Bdev + BNG/3.0d0
                
                Bloc = BGlob

                if(nlgeom == .false.) Bloc_mid = Bloc
                
                call MatPlastic3D(iprops,nds,(p+1)*(q+1)*(w+1),props,Bloc_mid,6,el_ddisp,stress(nel,cpi,:),&
                & strain(nel,cpi,:),dstrain(nel,cpi,:),dstress(nel,cpi,:),hard(nel,cpi),matD)
                
                !Store Constitutive Matrix
                TmatD(nel,cpi,:,:) = MatD(:,:)
                
                tempstr = 0.0d0
                do k1=1,ntens
                    tempstr(k1,1) = stress(nel,cpi,k1)
                    deform(k1,1) = dstrain(nel,cpi,k1)
                end do
                
                if (nlgeom == .true.) then
                
                    BNL = 0.0d0
                    do k1=1,inodes
                        dRdenc=0.0d0
                        dRdenc(1,1)=dRdx(k1,1) !
                        dRdenc(2,1)=dRdx(k1,2) !
                        dRdenc(3,1)=dRdx(k1,3) !
                        
                        dRdxyz = 0.0d0
                        dRdxyz = matmul(trans,dRdenc)

                        BNL(1,k1*3-2)=dRdxyz(1,1)
                        BNL(2,k1*3-2)=dRdxyz(2,1)
                        BNL(3,k1*3-2)=dRdxyz(3,1)
                        BNL(4,k1*3-1)=dRdxyz(1,1)
                        BNL(5,k1*3-1)=dRdxyz(2,1)
                        BNL(6,k1*3-1)=dRdxyz(3,1)
                        BNL(7,k1*3  )=dRdxyz(1,1)
                        BNL(8,k1*3  )=dRdxyz(2,1)
                        BNL(9,k1*3  )=dRdxyz(3,1)
                    end do 
                 
                    MSTR=0.0d0
                    do k1=1,nds
                        MSTR(k1*3-2,k1*3-2) = tempstr(1,1)
                        MSTR(k1*3-2,k1*3-1) = tempstr(4,1)
                        MSTR(k1*3-2,k1*3  ) = tempstr(5,1)
                        MSTR(k1*3-1,k1*3-2) = tempstr(4,1)
                        MSTR(k1*3-1,k1*3-1) = tempstr(2,1)
                        MSTR(k1*3-1,k1*3  ) = tempstr(6,1)
                        MSTR(k1*3  ,k1*3-2) = tempstr(5,1)
                        MSTR(k1*3  ,k1*3-1) = tempstr(6,1)
                        MSTR(k1*3  ,k1*3  ) = tempstr(3,1)
                    end do
                    
                    KNLG = KNLG + matmul(matmul(transpose(BNL),MSTR),BNL)*gwt
                    
                    !Export local axis ---------------
                    laxis(nel,cpi,:,:) = rconv
 
                end if
                
                if(nlgeom == .false.) then
                    Ke = Ke + matmul(matmul(transpose(Bloc),matD),Bloc)*gwt
                else
                    Ke = Ke + matmul(matmul(transpose(Bloc),matD),Bloc)*gwt + matmul(matmul(transpose(BNL),MSTR),BNL)*gwt
                end if
                
                
                Kua = Kua + matmul(matmul(transpose(Bloc),matD),BEAS)*gwt
                Kaa  =Kaa + matmul(matmul(transpose(BEAS),matD),BEAS)*gwt
                
                
                
                
                Finte = Finte + matmul(transpose(Bloc),tempstr)*gwt
                
                
                !Gravitic loads -------------------------------------------------------------
                if (gravity==.true.)then
                    call grvt(nshpl,nds,R,gwt,gconst,gravdir,density,Finte)
                end if
                
                !Pressure loads -------------------------------------------------------------
!                if(pressure==.true.) then
!                    do k1=1,nelems
!                        if(k1==nel .and. abs(pressvec(k1)) .gt. 0.0d0) then
!                            call prssr(nel,npi_xi,npi_eta,npi_zeta,jac,we(i),wn(j),wc(k),R,Finte)
!                        end if
!                    end do
!                end if
                
                !Strain Energy --------------------------------------------------------------
                do k1=1,ntens
                    SEnergy = SEnergy + tempstr(k1,1)*deform(k1,1)*gwt/2.0d0
                end do
                
                nj = INN(IEN(nel,1),2)
    
                yc  = ((v_knot(nj+1) - v_knot(nj))*eta + (v_knot(nj+1) + v_knot(nj)))/2.0d0*10.0
                zc = ((w_knot(nk+1) - w_knot(nk))*zeta + (w_knot(nk+1) + w_knot(nk)))/2.0d0
                
                continue
                
            end do  
        end do
    end do
    
    !tempa=0.0d0
    !call gaussj(kaa,ialpha,tempa,ilixo)
    
    !Ke = Ke  - matmul(matmul(Kua,Kaa),transpose(Kua))
    
    
    
    
    
    continue

end subroutine ElHex27EAS
