!------------------------------------------------------------------------------------------------------------
!
! General hexaedral element with a 2x2x2 gaussian quadrature
!
! Input: idx1, idx2 - NURBS coordinates (element indexes in the index space) (currently not being used)
!        nel - Element number
!        el_ddisp - Element displacement
!
! Output: Ke - Elementar stiffness matrix
!         Finte - Elementar internal forces
!------------------------------------------------------------------------------------------------------------

subroutine ElHex8DSG(nel,Ke,el_ddisp,Finte,inc)
    
    use Mod_Variables
    implicit none
    
    integer(4)::nel,cpi,count,inc
    
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::el_ddisp
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(OUT)::Ke
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(OUT)::Finte
    
    integer(4),parameter::npi_xi= 2, npi_eta = 2, npi_zeta = 2
    
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
    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BNat,BnatLoc
    
    !DSG
    real(8)::ge1,ge2,ge3,gn1,gn2,gn3,gc1,gc2,gc3
    
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
    
    !#####
    if(nlgeom==.true.)then
        tdisp = el_ddisp
        mdisp = el_ddisp/2.0d0
    else
        tdisp = 0.0d0
        mdisp= 0.0d0
        
        TGL = 0.0d0
        TCL = 0.0d0
        do i=1,6
            TGL(i,i) = 1.0d0
            TCL(i,i) = 1.0d0
        end do
    end if
    !#####
    
    !Gauss points cycle
    Ke = 0.0d0
    KNLG = 0.0d0
    BNL = 0.0d0
    MSTR = 0.0d0
    
    Finte = 0.0d0
    do i=1,npi_xi
        do j=1,npi_eta
            do k=1,npi_zeta
            
                cpi = cpi + 1
                MatD(:,:) = TmatD(nel,cpi,:,:)
                
                xi = e(i)
                eta = n(j)
                zeta = c(k)
                
                call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tzero)
                
                !Weight factor
                gwt = we(i)*wn(j)*wc(k)*detj
                
                
                !#########
                !if(nlgeom==.true.)then
                    
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
                        nodes(k1,1) = Points(IEN(nel,count),1) !+ u(IEN(nel,count)*3-2,1)
                        nodes(k1,2) = Points(IEN(nel,count),2) !+ u(IEN(nel,count)*3-1,1)
                        nodes(k1,3) = Points(IEN(nel,count),3) !+ u(IEN(nel,count)*3  ,1)
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
!                    matF = 0.0d0
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
                    
                !end if
                !############
                
                
                Bglob = 0.0d0
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
                    
                end do
                
                
                ge1 = jac(1,1)
                ge2 = jac(1,2)
                ge3 = jac(1,3)
                gn1 = jac(2,1)
                gn2 = jac(2,2)
                gn3 = jac(2,3)
                gc1 = jac(3,1)
                gc2 = jac(3,2)
                gc3 = jac(3,3)
                
                Bnat = 0.0d0
                
                Bnat(1,1)  =   ge1*(n(j)/8.0d0 + c(k)/8.0d0 - (n(j)*c(k))/8.0d0 - 1.0d0/8.0d0) 
                Bnat(1,2)  =   ge2*(n(j)/8.0d0 + c(k)/8.0d0 - (n(j)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(1,3)  =   ge3*(n(j)/8.0d0 + c(k)/8.0d0 - (n(j)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(1,4)  = - ge1*(n(j)/8.0d0 + c(k)/8.0d0 - (n(j)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(1,5)  = - ge2*(n(j)/8.0d0 + c(k)/8.0d0 - (n(j)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(1,6)  = - ge3*(n(j)/8.0d0 + c(k)/8.0d0 - (n(j)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(1,7)  =   ge1*(n(j)/8.0d0 - c(k)/8.0d0 - (n(j)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(1,8)  =   ge2*(n(j)/8.0d0 - c(k)/8.0d0 - (n(j)*c(k))/8.0d0 + 1.0d0/8.0d0) 
                Bnat(1,9)  =   ge3*(n(j)/8.0d0 - c(k)/8.0d0 - (n(j)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(1,10) = - ge1*(n(j)/8.0d0 - c(k)/8.0d0 - (n(j)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(1,11) = - ge2*(n(j)/8.0d0 - c(k)/8.0d0 - (n(j)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(1,12) = - ge3*(n(j)/8.0d0 - c(k)/8.0d0 - (n(j)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(1,13) =   ge1*(n(j)/8.0d0 - c(k)/8.0d0 + (n(j)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(1,14) =   ge2*(n(j)/8.0d0 - c(k)/8.0d0 + (n(j)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(1,15) =   ge3*(n(j)/8.0d0 - c(k)/8.0d0 + (n(j)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(1,16) = - ge1*(n(j)/8.0d0 - c(k)/8.0d0 + (n(j)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(1,17) = - ge2*(n(j)/8.0d0 - c(k)/8.0d0 + (n(j)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(1,18) = - ge3*(n(j)/8.0d0 - c(k)/8.0d0 + (n(j)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(1,19) =   ge1*(n(j)/8.0d0 + c(k)/8.0d0 + (n(j)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(1,20) =   ge2*(n(j)/8.0d0 + c(k)/8.0d0 + (n(j)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(1,21) =   ge3*(n(j)/8.0d0 + c(k)/8.0d0 + (n(j)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(1,22) = - ge1*(n(j)/8.0d0 + c(k)/8.0d0 + (n(j)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(1,23) = - ge2*(n(j)/8.0d0 + c(k)/8.0d0 + (n(j)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(1,24) = - ge3*(n(j)/8.0d0 + c(k)/8.0d0 + (n(j)*c(k))/8.0d0 + 1.0d0/8.0d0) 

                Bnat(2,1)  =   gn1*(e(i)/8.0d0 + c(k)/8.0d0 - (e(i)*c(k))/8.0d0 - 1.0d0/8.0d0) 
                Bnat(2,2)  =   gn2*(e(i)/8.0d0 + c(k)/8.0d0 - (e(i)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(2,3)  =   gn3*(e(i)/8.0d0 + c(k)/8.0d0 - (e(i)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(2,4)  = - gn1*(e(i)/8.0d0 - c(k)/8.0d0 - (e(i)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(2,5)  = - gn2*(e(i)/8.0d0 - c(k)/8.0d0 - (e(i)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(2,6)  = - gn3*(e(i)/8.0d0 - c(k)/8.0d0 - (e(i)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(2,7)  =   gn1*(e(i)/8.0d0 - c(k)/8.0d0 - (e(i)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(2,8)  =   gn2*(e(i)/8.0d0 - c(k)/8.0d0 - (e(i)*c(k))/8.0d0 + 1.0d0/8.0d0) 
                Bnat(2,9)  =   gn3*(e(i)/8.0d0 - c(k)/8.0d0 - (e(i)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(2,10) = - gn1*(e(i)/8.0d0 + c(k)/8.0d0 - (e(i)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(2,11) = - gn2*(e(i)/8.0d0 + c(k)/8.0d0 - (e(i)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(2,12) = - gn3*(e(i)/8.0d0 + c(k)/8.0d0 - (e(i)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(2,13) =   gn1*(e(i)/8.0d0 - c(k)/8.0d0 + (e(i)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(2,14) =   gn2*(e(i)/8.0d0 - c(k)/8.0d0 + (e(i)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(2,15) =   gn3*(e(i)/8.0d0 - c(k)/8.0d0 + (e(i)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(2,16) =   gn1*(e(i)/8.0d0 + c(k)/8.0d0 + (e(i)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(2,17) = - gn2*(e(i)/8.0d0 + c(k)/8.0d0 + (e(i)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(2,18) = - gn3*(e(i)/8.0d0 + c(k)/8.0d0 + (e(i)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(2,19) =   gn1*(e(i)/8.0d0 + c(k)/8.0d0 + (e(i)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(2,20) =   gn2*(e(i)/8.0d0 + c(k)/8.0d0 + (e(i)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(2,21) =   gn3*(e(i)/8.0d0 + c(k)/8.0d0 + (e(i)*c(k))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(2,22) = - gn1*(e(i)/8.0d0 - c(k)/8.0d0 + (e(i)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(2,23) = - gn2*(e(i)/8.0d0 - c(k)/8.0d0 + (e(i)*c(k))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(2,24) = - gn3*(e(i)/8.0d0 - c(k)/8.0d0 + (e(i)*c(k))/8.0d0 - 1.0d0/8.0d0) 

                Bnat(3,1)  =   gc1*(n(j)/8.0d0 + e(i)/8.0d0 - (n(j)*e(i))/8.0d0 - 1.0d0/8.0d0) 
                Bnat(3,2)  =   gc2*(n(j)/8.0d0 + e(i)/8.0d0 - (n(j)*e(i))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(3,3)  =   gc3*(n(j)/8.0d0 + e(i)/8.0d0 - (n(j)*e(i))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(3,4)  =   gc1*(n(j)/8.0d0 - e(i)/8.0d0 + (n(j)*e(i))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(3,5)  =   gc2*(n(j)/8.0d0 - e(i)/8.0d0 + (n(j)*e(i))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(3,6)  =   gc3*(n(j)/8.0d0 - e(i)/8.0d0 + (n(j)*e(i))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(3,7)  = - gc1*(n(j)/8.0d0 + e(i)/8.0d0 + (n(j)*e(i))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(3,8)  = - gc2*(n(j)/8.0d0 + e(i)/8.0d0 + (n(j)*e(i))/8.0d0 + 1.0d0/8.0d0) 
                Bnat(3,9)  = - gc3*(n(j)/8.0d0 + e(i)/8.0d0 + (n(j)*e(i))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(3,10) = - gc1*(n(j)/8.0d0 - e(i)/8.0d0 - (n(j)*e(i))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(3,11) = - gc2*(n(j)/8.0d0 - e(i)/8.0d0 - (n(j)*e(i))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(3,12) = - gc3*(n(j)/8.0d0 - e(i)/8.0d0 - (n(j)*e(i))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(3,13) = - gc1*(n(j)/8.0d0 + e(i)/8.0d0 - (n(j)*e(i))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(3,14) = - gc2*(n(j)/8.0d0 + e(i)/8.0d0 - (n(j)*e(i))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(3,15) = - gc3*(n(j)/8.0d0 + e(i)/8.0d0 - (n(j)*e(i))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(3,16) = - gc1*(n(j)/8.0d0 - e(i)/8.0d0 + (n(j)*e(i))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(3,17) = - gc2*(n(j)/8.0d0 - e(i)/8.0d0 + (n(j)*e(i))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(3,18) = - gc3*(n(j)/8.0d0 - e(i)/8.0d0 + (n(j)*e(i))/8.0d0 - 1.0d0/8.0d0)  
                Bnat(3,19) =   gc1*(n(j)/8.0d0 + e(i)/8.0d0 + (n(j)*e(i))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(3,20) =   gc2*(n(j)/8.0d0 + e(i)/8.0d0 + (n(j)*e(i))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(3,21) =   gc3*(n(j)/8.0d0 + e(i)/8.0d0 + (n(j)*e(i))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(3,22) =   gc1*(n(j)/8.0d0 - e(i)/8.0d0 - (n(j)*e(i))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(3,23) =   gc2*(n(j)/8.0d0 - e(i)/8.0d0 - (n(j)*e(i))/8.0d0 + 1.0d0/8.0d0)  
                Bnat(3,24) =   gc3*(n(j)/8.0d0 - e(i)/8.0d0 - (n(j)*e(i))/8.0d0 + 1.0d0/8.0d0) 

                Bnat(4,1)  =   ((ge1*(c(k)/8.0d0 - 1.0d0/8.0d0)) + (gn1*(c(k)/8.0d0 - 1.0d0/8.0d0))) 
                Bnat(4,2)  =   ((ge2*(c(k)/8.0d0 - 1.0d0/8.0d0)) + (gn2*(c(k)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(4,3)  =   ((ge3*(c(k)/8.0d0 - 1.0d0/8.0d0)) + (gn3*(c(k)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(4,4)  =   ((ge1*(c(k)/8.0d0 - 1.0d0/8.0d0)) - (gn1*(c(k)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(4,5)  =   ((ge2*(c(k)/8.0d0 - 1.0d0/8.0d0)) - (gn2*(c(k)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(4,6)  =   ((ge3*(c(k)/8.0d0 - 1.0d0/8.0d0)) - (gn3*(c(k)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(4,7)  = - ((ge1*(c(k)/8.0d0 - 1.0d0/8.0d0)) + (gn1*(c(k)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(4,8)  = - ((ge2*(c(k)/8.0d0 - 1.0d0/8.0d0)) + (gn2*(c(k)/8.0d0 - 1.0d0/8.0d0))) 
                Bnat(4,9)  = - ((ge3*(c(k)/8.0d0 - 1.0d0/8.0d0)) + (gn3*(c(k)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(4,10) = - ((ge1*(c(k)/8.0d0 - 1.0d0/8.0d0)) - (gn1*(c(k)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(4,11) = - ((ge2*(c(k)/8.0d0 - 1.0d0/8.0d0)) - (gn2*(c(k)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(4,12) = - ((ge3*(c(k)/8.0d0 - 1.0d0/8.0d0)) - (gn3*(c(k)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(4,13) = - ((ge1*(c(k)/8.0d0 + 1.0d0/8.0d0)) + (gn1*(c(k)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(4,14) = - ((ge2*(c(k)/8.0d0 + 1.0d0/8.0d0)) + (gn2*(c(k)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(4,15) = - ((ge3*(c(k)/8.0d0 + 1.0d0/8.0d0)) + (gn3*(c(k)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(4,16) = - ((ge1*(c(k)/8.0d0 + 1.0d0/8.0d0)) - (gn1*(c(k)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(4,17) = - ((ge2*(c(k)/8.0d0 + 1.0d0/8.0d0)) - (gn2*(c(k)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(4,18) = - ((ge3*(c(k)/8.0d0 + 1.0d0/8.0d0)) - (gn3*(c(k)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(4,19) =   ((ge1*(c(k)/8.0d0 + 1.0d0/8.0d0)) + (gn1*(c(k)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(4,20) =   ((ge2*(c(k)/8.0d0 + 1.0d0/8.0d0)) + (gn2*(c(k)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(4,21) =   ((ge3*(c(k)/8.0d0 + 1.0d0/8.0d0)) + (gn3*(c(k)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(4,22) =   ((ge1*(c(k)/8.0d0 + 1.0d0/8.0d0)) - (gn1*(c(k)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(4,23) =   ((ge2*(c(k)/8.0d0 + 1.0d0/8.0d0)) - (gn2*(c(k)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(4,24) =   ((ge3*(c(k)/8.0d0 + 1.0d0/8.0d0)) - (gn3*(c(k)/8.0d0 + 1.0d0/8.0d0))) 

                Bnat(5,1)  =   ((gc1*(n(j)/8.0d0 - 1.0d0/8.0d0)) + (ge1*(n(j)/8.0d0 - 1.0d0/8.0d0))) 
                Bnat(5,2)  =   ((gc2*(n(j)/8.0d0 - 1.0d0/8.0d0)) + (ge2*(n(j)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(5,3)  =   ((gc3*(n(j)/8.0d0 - 1.0d0/8.0d0)) + (ge3*(n(j)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(5,4)  = - ((gc1*(n(j)/8.0d0 - 1.0d0/8.0d0)) - (ge1*(n(j)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(5,5)  = - ((gc2*(n(j)/8.0d0 - 1.0d0/8.0d0)) - (ge2*(n(j)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(5,6)  = - ((gc3*(n(j)/8.0d0 - 1.0d0/8.0d0)) - (ge3*(n(j)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(5,7)  =   ((gc1*(n(j)/8.0d0 + 1.0d0/8.0d0)) - (ge1*(n(j)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(5,8)  =   ((gc2*(n(j)/8.0d0 + 1.0d0/8.0d0)) - (ge2*(n(j)/8.0d0 + 1.0d0/8.0d0))) 
                Bnat(5,9)  =   ((gc3*(n(j)/8.0d0 + 1.0d0/8.0d0)) - (ge3*(n(j)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(5,10) = - ((gc1*(n(j)/8.0d0 + 1.0d0/8.0d0)) + (ge1*(n(j)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(5,11) = - ((gc2*(n(j)/8.0d0 + 1.0d0/8.0d0)) + (ge2*(n(j)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(5,12) = - ((gc3*(n(j)/8.0d0 + 1.0d0/8.0d0)) + (ge3*(n(j)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(5,13) =   ((gc1*(n(j)/8.0d0 - 1.0d0/8.0d0)) - (ge1*(n(j)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(5,14) =   ((gc2*(n(j)/8.0d0 - 1.0d0/8.0d0)) - (ge2*(n(j)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(5,15) =   ((gc3*(n(j)/8.0d0 - 1.0d0/8.0d0)) - (ge3*(n(j)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(5,16) = - ((gc1*(n(j)/8.0d0 - 1.0d0/8.0d0)) + (ge1*(n(j)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(5,17) = - ((gc2*(n(j)/8.0d0 - 1.0d0/8.0d0)) + (ge2*(n(j)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(5,18) = - ((gc3*(n(j)/8.0d0 - 1.0d0/8.0d0)) + (ge3*(n(j)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(5,19) =   ((gc1*(n(j)/8.0d0 + 1.0d0/8.0d0)) + (ge1*(n(j)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(5,20) =   ((gc2*(n(j)/8.0d0 + 1.0d0/8.0d0)) + (ge2*(n(j)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(5,21) =   ((gc3*(n(j)/8.0d0 + 1.0d0/8.0d0)) + (ge3*(n(j)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(5,22) = - ((gc1*(n(j)/8.0d0 + 1.0d0/8.0d0)) - (ge1*(n(j)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(5,23) = - ((gc2*(n(j)/8.0d0 + 1.0d0/8.0d0)) - (ge2*(n(j)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(5,24) = - ((gc3*(n(j)/8.0d0 + 1.0d0/8.0d0)) - (ge3*(n(j)/8.0d0 + 1.0d0/8.0d0))) 

                Bnat(6,1)  =   ((gc1*(e(i)/8.0d0 - 1.0d0/8.0d0)) + (gn1*(e(i)/8.0d0 - 1.0d0/8.0d0))) 
                Bnat(6,2)  =   ((gc2*(e(i)/8.0d0 - 1.0d0/8.0d0)) + (gn2*(e(i)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(6,3)  =   ((gc3*(e(i)/8.0d0 - 1.0d0/8.0d0)) + (gn3*(e(i)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(6,4)  = - ((gc1*(e(i)/8.0d0 + 1.0d0/8.0d0)) + (gn1*(e(i)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(6,5)  = - ((gc2*(e(i)/8.0d0 + 1.0d0/8.0d0)) + (gn2*(e(i)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(6,6)  = - ((gc3*(e(i)/8.0d0 + 1.0d0/8.0d0)) + (gn3*(e(i)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(6,7)  =   ((gc1*(e(i)/8.0d0 + 1.0d0/8.0d0)) - (gn1*(e(i)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(6,8)  =   ((gc2*(e(i)/8.0d0 + 1.0d0/8.0d0)) - (gn2*(e(i)/8.0d0 + 1.0d0/8.0d0))) 
                Bnat(6,9)  =   ((gc3*(e(i)/8.0d0 + 1.0d0/8.0d0)) - (gn3*(e(i)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(6,10) = - ((gc1*(e(i)/8.0d0 - 1.0d0/8.0d0)) - (gn1*(e(i)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(6,11) = - ((gc2*(e(i)/8.0d0 - 1.0d0/8.0d0)) - (gn2*(e(i)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(6,12) = - ((gc3*(e(i)/8.0d0 - 1.0d0/8.0d0)) - (gn3*(e(i)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(6,13) =   ((gc1*(e(i)/8.0d0 - 1.0d0/8.0d0)) - (gn1*(e(i)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(6,14) =   ((gc2*(e(i)/8.0d0 - 1.0d0/8.0d0)) - (gn2*(e(i)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(6,15) =   ((gc3*(e(i)/8.0d0 - 1.0d0/8.0d0)) - (gn3*(e(i)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(6,16) = - ((gc1*(e(i)/8.0d0 + 1.0d0/8.0d0)) - (gn1*(e(i)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(6,17) = - ((gc2*(e(i)/8.0d0 + 1.0d0/8.0d0)) - (gn2*(e(i)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(6,18) = - ((gc3*(e(i)/8.0d0 + 1.0d0/8.0d0)) - (gn3*(e(i)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(6,19) =   ((gc1*(e(i)/8.0d0 + 1.0d0/8.0d0)) + (gn1*(e(i)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(6,20) =   ((gc2*(e(i)/8.0d0 + 1.0d0/8.0d0)) + (gn2*(e(i)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(6,21) =   ((gc3*(e(i)/8.0d0 + 1.0d0/8.0d0)) + (gn3*(e(i)/8.0d0 + 1.0d0/8.0d0)))  
                Bnat(6,22) = - ((gc1*(e(i)/8.0d0 - 1.0d0/8.0d0)) + (gn1*(e(i)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(6,23) = - ((gc2*(e(i)/8.0d0 - 1.0d0/8.0d0)) + (gn2*(e(i)/8.0d0 - 1.0d0/8.0d0)))  
                Bnat(6,24) = - ((gc3*(e(i)/8.0d0 - 1.0d0/8.0d0)) + (gn3*(e(i)/8.0d0 - 1.0d0/8.0d0)))
                
                BnatLoc = matmul(TCL,BNat)
                
                !Local strain-displacement matrix at end configuration
                Bloc=0.0d0
                Bloc=matmul(TGL,Bglob)
                
                Bloc = BnatLoc
                
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
                
                !#########
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
                !###########
                
                if(nlgeom == .false.) then
                    Ke = Ke + matmul(matmul(transpose(Bloc),matD),Bloc)*gwt
                else
                    Ke = Ke + matmul(matmul(transpose(Bloc),matD),Bloc)*gwt + matmul(matmul(transpose(BNL),MSTR),BNL)*gwt
                end if
                
                Finte = Finte + matmul(transpose(Bloc),tempstr)*gwt
                
                
                !Gravitic loads -------------------------------------------------------------
                if (gravity==.true.)then
                    call grvt(nshpl,nds,R,gwt,gconst,gravdir,density,Finte)
                end if
                
                !Pressure loads -------------------------------------------------------------
                if(pressure==.true.) then
                    do k1=1,nelems
                        if(k1==nel .and. abs(pressvec(k1)) .gt. 0.0d0) then
                            call prssr(nel,npi_xi,npi_eta,npi_zeta,jac,we(i),wn(j),wc(k),R,Finte)
                        end if
                    end do
                end if
                
                !Strain Energy --------------------------------------------------------------
                do k1=1,ntens
                    SEnergy = SEnergy + tempstr(k1,1)*deform(k1,1)*gwt/2.0d0
                end do
                
                continue
                
            end do  
        end do
    end do

    continue

end subroutine ElHex8DSG
