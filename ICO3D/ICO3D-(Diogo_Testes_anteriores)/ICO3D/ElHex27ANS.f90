!------------------------------------------------------------------------------------------------------------
!
! General hexaedral element with a 3x3x3 gaussian quadrature
!
! Input: idx1, idx2 - NURBS coordinates (element indexes in the index space) (currently not being used)
!        nel - Element number
!        el_ddisp - Element displacement
!
! Output: Ke - Elementar stiffness matrix
!         Finte - Elementar internal forces
!------------------------------------------------------------------------------------------------------------


!subroutine ElHex27ANS(nel,Ke,el_ddisp,Finte,inc)
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
!    integer(4)::i,j,k,ni,nj,nk,nel_nza,k1,k2
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
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob,Bloc,Bloc2,BCN
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
!    real(8),dimension(nds*2,nds*2)::TCL,TGL,TLG,TLC
!    real(8),dimension(9,(p+1)*(q+1)*(w+1)*nds)::BNL
!    real(8),dimension(9,9)::MSTR
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds)::KNLG
!    real(8),dimension(nds,1)::dRdenc,dRdxyz
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob_mid,Bloc_mid
!    real(8),dimension(nds*2)::temp6
!    real(8),dimension(nds*2,nds*2)::TCL0,TLC0
!    
!    !ANS
!    real(8)::detj0,detjA,A,B
!    real(8),dimension(nds,nds)::dxdxiA,dxdxi,jac0,jacA
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BNat,BNat2,BCN2
!    real(8),dimension((p+1)*(q+1)*(w+1))::RA
!    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdxi,dRdxii
!    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdxA,dRdxiA,dRdxiiA
!    real(8),dimension(16)::shpANS
!    real(8),dimension(16,3)::coordANS
!    real(8),dimension(3,1)::dRdx_teste,temp31
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BNatural,BN2
!    real(8),dimension(nds,nds)::Mat33
!    
!    real(8),dimension(nds*2,nds*2)::TGC, TGC0
!    
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
!    
!    
!    call ShapeFunc3(nel,0.0d0,0.0d0,0.0d0,R,dRdx,dRdxi,dRdxii,detj0,jac0,dxdxi,tzero)
!    
!    call TransformationMat3D(jac0,TCL0)
!                
!    temp6 = 0.0d0
!    TLC0 = TCL0
!    call gaussj(TLC0,6,temp6,ilixo)
!    
!    Mat33 = jac0 !matmul(transpose(rconv),jac)
!    call TransformationMat3D(mat33,TGC0)
!
!    !Gauss points cycle
!    Ke = 0.0d0
!    KNLG = 0.0d0
!    BNL = 0.0d0
!    MSTR = 0.0d0
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
!                call ShapeFunc3(nel,xi,eta,0.0d0,R,dRdx,dRdxi,dRdxii,detj0,jac0,dxdxi,tzero)
!    
!                call TransformationMat3D(jac0,TCL0)
!                
!                temp6 = 0.0d0
!                TLC0 = TCL0
!                call gaussj(TLC0,6,temp6,ilixo)
!    
!                Mat33 = jac0 !matmul(transpose(rconv),jac)
!                call TransformationMat3D(mat33,TGC0)
!                
!                call ShapeFunc3(nel,xi,eta,zeta,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,tzero)
!                
!                !Weight factor
!                gwt = we(i)*wn(j)*wc(k)*detj
!                
!                    
!                !Local Axis ---------------------------
!                if(inc==1)then
!                    call local_axis(nds,jac,rconv)
!                else
!                    rconv = laxis(nel,cpi,:,:)
!                end if
!                
!                !-------------------------------------------------------------------------
!                !Mid-Configuration
!                !-------------------------------------------------------------------------
!                matF_mid = 0.0d0
!                do k1=1,inodes
!                    dRde(k1) = dRdx(k1,1)
!                    dRdn(k1) = dRdx(k1,2)
!                    dRdc(k1) = dRdx(k1,3)
!                end do
!                
!                jacinv = 0.0d0
!                jacinv = jac
!                call gaussj(jacinv,nds,temp1,ilixo)
!                
!                count = 0
!                do k1=1,inodes
!                    count = count + 1
!                    nodes(k1,1) = Points(IEN(nel,count),1) !+ u(IEN(nel,count)*3-2,1)
!                    nodes(k1,2) = Points(IEN(nel,count),2) !+ u(IEN(nel,count)*3-1,1)
!                    nodes(k1,3) = Points(IEN(nel,count),3) !+ u(IEN(nel,count)*3  ,1)
!                end do
!                
!                !Deformation Gradient for mid-configuration
!                call DefGrad3D(inodes,nds,dRde,dRdn,dRdc,jacinv,nodes,mdisp,MatF_mid)
!                
!                !Polar Decomposition for mid configuration
!                matR_mid=0.0d0
!                call PolarDecomp2D(nds,matF_mid,matR_mid)
!                
!                !Update local reference system to the mid configuration
!                call axisupdate(rconv,matR_mid,rconv_mid)
!                
!                !Global to natural jacobian in mid configuration
!                call ShapeFunc2(nel,xi,eta,zeta,R_mid,dRdx_mid,detj_mid,jac_mid,mdisp)
!                
!                jacinv_mid = 0.0d0
!                temp1 = 0.0d0
!                jacinv_mid = jac_mid
!                call gaussj(jacinv_mid,nds,temp1,ilixo)
!                
!                !Transformation matrices mid configuration 
!                Temp33_mid=0.0d0
!                Temp33_mid= matmul(transpose(rconv_mid),jacinv_mid)
!                Trans_mid=0.0d0
!                Trans_mid=transpose(rconv_mid)
!                
!                !Natural to local transformation matrix for mid configuration
!                call TransformationMat3D(Temp33_mid,TCL_mid)
!                
!                !Global to local transformation matrix for mid configuration
!                call TransformationMat3D(Trans_mid,TGL_mid)
!                
!                Bglob_mid = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
!                    Bglob_mid(1,k1*3-2) = dRdx_mid(k1,1)
!                    Bglob_mid(2,k1*3-1) = dRdx_mid(k1,2)
!                    Bglob_mid(3,k1*3  ) = dRdx_mid(k1,3)
!                    Bglob_mid(4,k1*3-2) = dRdx_mid(k1,2)
!                    Bglob_mid(4,k1*3-1) = dRdx_mid(k1,1)
!                    Bglob_mid(5,k1*3-2) = dRdx_mid(k1,3)
!                    Bglob_mid(5,k1*3  ) = dRdx_mid(k1,1)
!                    Bglob_mid(6,k1*3-1) = dRdx_mid(k1,3)
!                    Bglob_mid(6,k1*3  ) = dRdx_mid(k1,2)
!                end do
!                
!                Bloc_mid=0.0d0
!                Bloc_mid=matmul(TGL_mid,Bglob_mid)
!                !-------------------------------------------------------------------------
!                !End-Configuration
!                !-------------------------------------------------------------------------
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
!                !Deformation Gradient for end-configuration
!                call DefGrad3D(inodes,nds,dRde,dRdn,dRdc,jacinv,nodes,tdisp,MatF)
!                
!                !Polar Decomposition for end configuration
!                matR = 0.0d0
!                call PolarDecomp2D(nds,matF,matR)
!                
!                !Update local reference system to the end configuration
!                raux = rconv
!                call axisupdate(rconv,matR,rconv)
!                
!                !Global to natural jacobian in end configuration
!                !call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tdisp)
!                call ShapeFunc3(nel,xi,eta,zeta,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,tzero)
!                
!                jacinv = 0.0d0
!                temp1 = 0.0d0
!                jacinv = jac
!                call gaussj(jacinv,nds,temp1,ilixo)
!                
!                
!                !Transformation matrices end configuration 
!                Temp33 = 0.0d0
!                Temp33 = matmul(transpose(rconv),jacinv)
!                Trans = 0.0d0
!                Trans = transpose(rconv)
!                
!                !Natural to local transformation matrix for end configuration
!                call TransformationMat3D(Temp33,TCL)
!                
!                temp6 = 0.0d0
!                TLC = TCL
!                call gaussj(TLC,6,temp6,ilixo)
!                
!                !Global to local transformation matrix for end configuration
!                call TransformationMat3D(Trans,TGL)
!                
!                temp6 = 0.0d0
!                TLG = TGL
!                call gaussj(TLG,6,temp6,ilixo)
!                
!                !Global to natural transformation matrix for end configuration
!                Mat33 = jac !matmul(transpose(rconv),jac)
!                call TransformationMat3D(mat33,TGC)
!
!                !Update Weight factor
!                gwt = we(i)*wn(j)*wc(k)*detj
!                
!                continue
!                
!                
!                
!                
!                
!                
!                
!                
!                !Compatible strain-displacement matrix in global frame
!                Bglob = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
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
!                end do
!                
!                !BNatural = matmul(TLC,Bglob)
!                BNatural = matmul(TGC,Bglob)
!                
!                BCN = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
!                    BCN(1,k1*3-2) = dRdx(k1,1)
!                    BCN(2,k1*3-1) = dRdx(k1,2)
!                    BCN(3,k1*3  ) = dRdx(k1,3)
!                    BCN(4,k1*3-2) = dRdx(k1,2)
!                    BCN(4,k1*3-1) = dRdx(k1,1)
!                    BCN(5,k1*3-2) = dRdx(k1,3)
!                    BCN(5,k1*3  ) = dRdx(k1,1)
!                    BCN(6,k1*3-1) = dRdx(k1,3)
!                    BCN(6,k1*3  ) = dRdx(k1,2)
!                end do
!                
!                BCN = BNatural
!                
!                !ANS Implementation I
!                A = sqrt(1.0d0/3.0d0)
!                B = sqrt(3.0d0/5.0d0)
!                
!                !-------------------------------------------------------------------
!                !-------------------------------------------------------------------
!                !Coordinates
!                coordANS = 0.0d0
!                coordANS(1,1)  =     A; coordANS(1,2)  =     B
!                coordANS(2,1)  =    -A; coordANS(2,2)  =     B
!                coordANS(3,1)  =     A; coordANS(3,2)  = 0.0d0
!                coordANS(4,1)  =    -A; coordANS(4,2)  = 0.0d0
!                coordANS(5,1)  =     A; coordANS(5,2)  =    -B
!                coordANS(6,1)  =    -A; coordANS(6,2)  =    -B
!                coordANS(7,1)  =     B; coordANS(7,2)  =     A
!                coordANS(8,1)  =     B; coordANS(8,2)  =    -A
!                coordANS(9,1)  = 0.0d0; coordANS(9,2)  =     A
!                coordANS(10,1) = 0.0d0; coordANS(10,2) =    -A
!                coordANS(11,1) =    -B; coordANS(11,2) =     A
!                coordANS(12,1) =    -B; coordANS(12,2) =    -A
!                coordANS(13,1) =     A; coordANS(13,2) =     A
!                coordANS(14,1) =     A; coordANS(14,2) =    -A
!                coordANS(15,1) =    -A; coordANS(15,2) =     A
!                coordANS(16,1) =    -A; coordANS(16,2) =    -A
!                
!                
!                coordANS = 0.0d0
!                coordANS(1,1)  =     A; coordANS(1,2)  =     B
!                coordANS(2,1)  =    -A; coordANS(2,2)  =     B
!                coordANS(3,1)  =     A; coordANS(3,2)  = 0.0d0
!                coordANS(4,1)  =    -A; coordANS(4,2)  = 0.0d0
!                coordANS(5,1)  =     A; coordANS(5,2)  =    -B
!                coordANS(6,1)  =    -A; coordANS(6,2)  =    -B
!                coordANS(7,1)  =     B; coordANS(7,2)  =     A
!                coordANS(8,1)  = 0.0d0; coordANS(8,2)  =     A
!                coordANS(9,1)  =    -B; coordANS(9,2)  =     A
!                coordANS(10,1) =     B; coordANS(10,2) =    -A
!                coordANS(11,1) = 0.0d0; coordANS(11,2) =    -A
!                coordANS(12,1) =    -B; coordANS(12,2) =    -A
!                coordANS(13,1) =     A; coordANS(13,2) =     A
!                coordANS(14,1) =    -A; coordANS(14,2) =     A
!                coordANS(15,1) =     A; coordANS(15,2) =    -A
!                coordANS(16,1) =    -A; coordANS(16,2) =    -A
!                !-------------------------------------------------------------------
!                !-------------------------------------------------------------------
!                
!
!                !Shape Functions
!                shpANS = 0.0d0
!                
!                !-------------------------------------------------------------------
!                shpANS(1)  = eta/(4.0d0*B)*(eta/B+1.0d0)*(1.0d0+xi/A) 
!                shpANS(2)  = eta/(4.0d0*B)*(eta/B+1.0d0)*(1.0d0-xi/A)
!                shpANS(3)  = 1.0d0/2.0d0*(1.0d0-eta*eta/(B*B))*(1.0d0+xi/A)
!                shpANS(4)  = 1.0d0/2.0d0*(1.0d0-eta*eta/(B*B))*(1.0d0-xi/A)
!                shpANS(5)  = eta/(4.0d0*B)*(eta/B-1.0d0)*(1.0d0+xi/A) 
!                shpANS(6)  = eta/(4.0d0*B)*(eta/B-1.0d0)*(1.0d0-xi/A)
!                
!!                shpANS(6)  = eta/(4.0d0*B)*(eta/B+1.0d0)*(1.0d0+xi/A) 
!!                shpANS(5)  = eta/(4.0d0*B)*(eta/B+1.0d0)*(1.0d0-xi/A)
!!                shpANS(4)  = 1.0d0/2.0d0*(1.0d0-eta*eta/(B*B))*(1.0d0+xi/A)
!!                shpANS(3)  = 1.0d0/2.0d0*(1.0d0-eta*eta/(B*B))*(1.0d0-xi/A)
!!                shpANS(2)  = eta/(4.0d0*B)*(eta/B-1.0d0)*(1.0d0+xi/A) 
!!                shpANS(1)  = eta/(4.0d0*B)*(eta/B-1.0d0)*(1.0d0-xi/A)
!                !-------------------------------------------------------------------
!                
!                
!                
!                !-------------------------------------------------------------------
!                shpANS(7)  = xi/(4.0d0*B)*(xi/B+1.0d0)*(1.0d0+eta/A)
!                shpANS(8)  = xi/(4.0d0*B)*(xi/B+1.0d0)*(1.0d0-eta/A)
!                shpANS(9)  = 1.0d0/2.0d0*(1.0d0-xi*xi/(B*B))*(1.0d0+eta/A)
!                shpANS(10) = 1.0d0/2.0d0*(1.0d0-xi*xi/(B*B))*(1.0d0-eta/A)
!                shpANS(11) = xi/(4.0d0*B)*(xi/B-1.0d0)*(1.0d0+eta/A)
!                shpANS(12) = xi/(4.0d0*B)*(xi/B-1.0d0)*(1.0d0-eta/A)
!                
!                shpANS(7)  = xi/(4.0d0*B)*(xi/B+1.0d0)*(1.0d0+eta/A)
!                shpANS(8)  = 1.0d0/2.0d0*(1.0d0-xi*xi/(B*B))*(1.0d0+eta/A)
!                shpANS(9)  = xi/(4.0d0*B)*(xi/B-1.0d0)*(1.0d0+eta/A)
!                shpANS(10) = xi/(4.0d0*B)*(xi/B+1.0d0)*(1.0d0-eta/A)
!                shpANS(11) = 1.0d0/2.0d0*(1.0d0-xi*xi/(B*B))*(1.0d0-eta/A)
!                shpANS(12) = xi/(4.0d0*B)*(xi/B-1.0d0)*(1.0d0-eta/A)
!                !-------------------------------------------------------------------
!                
!                
!                
!                
!                !-------------------------------------------------------------------
!                shpANS(13) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0+eta/A)
!                shpANS(14) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0-eta/A)
!                shpANS(15) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0+eta/A)
!                shpANS(16) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0-eta/A)
!                
!                shpANS(13) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0+eta/A)
!                shpANS(14) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0+eta/A)
!                shpANS(15) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0-eta/A)
!                shpANS(16) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0-eta/A)
!                !-------------------------------------------------------------------
!                
!                BNat = 0.0d0
!                
!                !xi-xi and eta-eta components
!                do k1=1,6
!                    Bnat(1,:) = Bnat(1,:) +   shpANS(k1)*BCN(1,:)
!                    Bnat(2,:) = Bnat(2,:) + shpANS(k1+6)*BCN(2,:)
!                end do
!                
!                !zeta-zeta components
!                BNat(3,:) = BCN(3,:)
!                
!                !xi-eta components
!                do k1=13,16
!                    Bnat(4,:) = Bnat(4,:) + shpANS(k1)*BCN(4,:)
!                end do
!                
!                !xi-zeta components
!                do k1=1,6
!                    call ShapeFunc3(nel,coordANS(k1,1),coordANS(k1,2),coordANS(k1,3),RA,dRdxA,dRdxiA,dRdxiiA,detjA,jacA,dxdxiA,tzero)
!                    
!                    BN2 = 0.0d0
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        BN2(1,k2*3-2) = dRdxA(k2,1)
!                        BN2(2,k2*3-1) = dRdxA(k2,2)
!                        BN2(3,k2*3  ) = dRdxA(k2,3)
!                    
!                        BN2(4,k2*3-2) = dRdxA(k2,2)
!                        BN2(4,k2*3-1) = dRdxA(k2,1)
!                        BN2(5,k2*3-2) = dRdxA(k2,3)
!                        BN2(5,k2*3  ) = dRdxA(k2,1)
!                        BN2(6,k2*3-1) = dRdxA(k2,3)
!                        BN2(6,k2*3  ) = dRdxA(k2,2)
!                    end do
!                
!                    !BN2 = matmul(TLC,BN2)
!                    BN2 = matmul(TGC,BN2)  !*detj0/detj
!                    
!                    Bnat(5,:) = Bnat(5,:) + shpANS(k1)*BN2(5,:) !*detjA/detj
!                    
!!                    do k2=1,(p+1)*(q+1)*(w+1)
!!                        BNat(5,k2*3-2) = BNat(5,k2*3-2) + shpANS(k1)*dRdxA(k2,3)
!!                        BNat(5,k2*3  ) = BNat(5,k2*3  ) + shpANS(k1)*dRdxA(k2,1)
!!                    end do
!
!                end do
!                
!                !eta-zeta components
!                do k1=7,12
!                    call ShapeFunc3(nel,coordANS(k1,1),coordANS(k1,2),coordANS(k1,3),RA,dRdxA,dRdxiA,dRdxiiA,detjA,jacA,dxdxiA,tzero)
!                    
!                    BN2 = 0.0d0
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        BN2(1,k2*3-2) = dRdxA(k2,1)
!                        BN2(2,k2*3-1) = dRdxA(k2,2)
!                        BN2(3,k2*3  ) = dRdxA(k2,3)
!                    
!                        BN2(4,k2*3-2) = dRdxA(k2,2)
!                        BN2(4,k2*3-1) = dRdxA(k2,1)
!                        BN2(5,k2*3-2) = dRdxA(k2,3)
!                        BN2(5,k2*3  ) = dRdxA(k2,1)
!                        BN2(6,k2*3-1) = dRdxA(k2,3)
!                        BN2(6,k2*3  ) = dRdxA(k2,2)
!                    end do
!                
!                    !BN2 = matmul(TLC,BN2)
!                    BN2 = matmul(TGC,BN2) !*detj0/detj
!                    
!                    Bnat(6,:) = Bnat(6,:) + shpANS(k1)*BN2(6,:) !*detjA/detj
!                    
!                    
!!                    do k2=1,(p+1)*(q+1)*(w+1)
!!                        BNat(6,k2*3-1) = BNat(6,k2*3-1) + shpANS(k1)*dRdxA(k2,3)
!!                        BNat(6,k2*3  ) = BNat(6,k2*3  ) + shpANS(k1)*dRdxA(k2,2)
!!                    end do
!
!                end do
!                
!                !Local strain-displacement matrix at end configuration
!                Bloc=0.0d0
!                Bloc=matmul(TGL,Bglob)
!                
!                Bloc=matmul(TCL,BNat) !*detj0/detj
!                
!                continue
!                
!                
!!                jacinv = 0.0d0
!!                temp1 = 0.0d0
!!                jacinv = dxdxi
!!                call gaussj(jacinv,nds,temp1,ilixo)
!!                
!!                !Transformation matrices mid configuration 
!!                Temp33 = 0.0d0
!!                Temp33 = matmul(transpose(rconv),jacinv)
!!                call TransformationMat3D(Temp33,TCL)
!!                
!!                do k1=1,(p+1)*(q+1)*(w+1)
!!                    BNat(1,k1*3-2) = dRdxii(k1,1)
!!                    BNat(2,k1*3-1) = dRdxii(k1,2)
!!                    BNat(3,k1*3  ) = dRdxii(k1,3)
!!                    
!!                    BNat(4,k1*3-2) = dRdxii(k1,2)
!!                    BNat(4,k1*3-1) = dRdxii(k1,1)
!!                end do
!!                
!!                !Bloc2=matmul(TGL,Bloc2)
!!                
!!                
!!                Bloc(5,:) = Bloc2(5,:)
!!                Bloc(6,:) = Bloc2(6,:)
!
!                !Bloc = Bloc2
!
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
!    continue
!
!end subroutine ElHex27ANS

!
!
!subroutine ElHex27ANS(nel,Ke,el_ddisp,Finte,inc)
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
!    integer(4)::i,j,k,ni,nj,nk,nel_nza,k1,k2
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
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob,Bloc,Bloc2,BCN
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
!    real(8),dimension(nds*2,nds*2)::TCL,TGL,TLG,TLC
!    real(8),dimension(9,(p+1)*(q+1)*(w+1)*nds)::BNL
!    real(8),dimension(9,9)::MSTR
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds)::KNLG
!    real(8),dimension(nds,1)::dRdenc,dRdxyz
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob_mid,Bloc_mid
!    real(8),dimension(nds*2)::temp6
!    real(8),dimension(nds*2,nds*2)::TCL0,TLC0
!    
!    !ANS
!    real(8)::detj0,detjA,A,B
!    real(8),dimension(nds,nds)::dxdxiA,dxdxi,jac0,jacA
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BNat,BNat2,BCN2
!    real(8),dimension((p+1)*(q+1)*(w+1))::RA
!    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdxi,dRdxii
!    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdxA,dRdxiA,dRdxiiA
!    real(8),dimension(50)::shpANS
!    real(8),dimension(50,3)::coordANS
!    real(8),dimension(3,1)::dRdx_teste,temp31
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BNatural,BN2
!    real(8),dimension(nds,nds)::Mat33
!    real(8),dimension(nds*2,nds*2)::TGC, TGC0
!    
!    integer(4),dimension(50)::logi
!    
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
!    
!    
!    call ShapeFunc3(nel,0.0d0,0.0d0,0.0d0,R,dRdx,dRdxi,dRdxii,detj0,jac0,dxdxi,tzero)
!    
!    call TransformationMat3D(jac0,TCL0)
!                
!    temp6 = 0.0d0
!    TLC0 = TCL0
!    call gaussj(TLC0,6,temp6,ilixo)
!    
!    Mat33 = jac0 !matmul(transpose(rconv),jac)
!    call TransformationMat3D(mat33,TGC0)
!
!    !Gauss points cycle
!    Ke = 0.0d0
!    KNLG = 0.0d0
!    BNL = 0.0d0
!    MSTR = 0.0d0
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
!                call ShapeFunc3(nel,xi,eta,0.0d0,R,dRdx,dRdxi,dRdxii,detj0,jac0,dxdxi,tzero)
!    
!                call TransformationMat3D(jac0,TCL0)
!                
!                temp6 = 0.0d0
!                TLC0 = TCL0
!                call gaussj(TLC0,6,temp6,ilixo)
!    
!                Mat33 = jac0 !matmul(transpose(rconv),jac)
!                call TransformationMat3D(mat33,TGC0)
!                
!                call ShapeFunc3(nel,xi,eta,zeta,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,tzero)
!                
!                !Weight factor
!                gwt = we(i)*wn(j)*wc(k)*detj
!                
!                    
!                !Local Axis ---------------------------
!                if(inc==1)then
!                    call local_axis(nds,jac,rconv)
!                else
!                    rconv = laxis(nel,cpi,:,:)
!                end if
!                
!                !-------------------------------------------------------------------------
!                !Mid-Configuration
!                !-------------------------------------------------------------------------
!                matF_mid = 0.0d0
!                do k1=1,inodes
!                    dRde(k1) = dRdx(k1,1)
!                    dRdn(k1) = dRdx(k1,2)
!                    dRdc(k1) = dRdx(k1,3)
!                end do
!                
!                jacinv = 0.0d0
!                jacinv = jac
!                call gaussj(jacinv,nds,temp1,ilixo)
!                
!                count = 0
!                do k1=1,inodes
!                    count = count + 1
!                    nodes(k1,1) = Points(IEN(nel,count),1) !+ u(IEN(nel,count)*3-2,1)
!                    nodes(k1,2) = Points(IEN(nel,count),2) !+ u(IEN(nel,count)*3-1,1)
!                    nodes(k1,3) = Points(IEN(nel,count),3) !+ u(IEN(nel,count)*3  ,1)
!                end do
!                
!                !Deformation Gradient for mid-configuration
!                call DefGrad3D(inodes,nds,dRde,dRdn,dRdc,jacinv,nodes,mdisp,MatF_mid)
!                
!                !Polar Decomposition for mid configuration
!                matR_mid=0.0d0
!                call PolarDecomp2D(nds,matF_mid,matR_mid)
!                
!                !Update local reference system to the mid configuration
!                call axisupdate(rconv,matR_mid,rconv_mid)
!                
!                !Global to natural jacobian in mid configuration
!                call ShapeFunc2(nel,xi,eta,zeta,R_mid,dRdx_mid,detj_mid,jac_mid,mdisp)
!                
!                jacinv_mid = 0.0d0
!                temp1 = 0.0d0
!                jacinv_mid = jac_mid
!                call gaussj(jacinv_mid,nds,temp1,ilixo)
!                
!                !Transformation matrices mid configuration 
!                Temp33_mid=0.0d0
!                Temp33_mid= matmul(transpose(rconv_mid),jacinv_mid)
!                Trans_mid=0.0d0
!                Trans_mid=transpose(rconv_mid)
!                
!                !Natural to local transformation matrix for mid configuration
!                call TransformationMat3D(Temp33_mid,TCL_mid)
!                
!                !Global to local transformation matrix for mid configuration
!                call TransformationMat3D(Trans_mid,TGL_mid)
!                
!                Bglob_mid = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
!                    Bglob_mid(1,k1*3-2) = dRdx_mid(k1,1)
!                    Bglob_mid(2,k1*3-1) = dRdx_mid(k1,2)
!                    Bglob_mid(3,k1*3  ) = dRdx_mid(k1,3)
!                    Bglob_mid(4,k1*3-2) = dRdx_mid(k1,2)
!                    Bglob_mid(4,k1*3-1) = dRdx_mid(k1,1)
!                    Bglob_mid(5,k1*3-2) = dRdx_mid(k1,3)
!                    Bglob_mid(5,k1*3  ) = dRdx_mid(k1,1)
!                    Bglob_mid(6,k1*3-1) = dRdx_mid(k1,3)
!                    Bglob_mid(6,k1*3  ) = dRdx_mid(k1,2)
!                end do
!                
!                Bloc_mid=0.0d0
!                Bloc_mid=matmul(TGL_mid,Bglob_mid)
!                !-------------------------------------------------------------------------
!                !End-Configuration
!                !-------------------------------------------------------------------------
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
!                !Deformation Gradient for end-configuration
!                call DefGrad3D(inodes,nds,dRde,dRdn,dRdc,jacinv,nodes,tdisp,MatF)
!                
!                !Polar Decomposition for end configuration
!                matR = 0.0d0
!                call PolarDecomp2D(nds,matF,matR)
!                
!                !Update local reference system to the end configuration
!                raux = rconv
!                call axisupdate(rconv,matR,rconv)
!                
!                !Global to natural jacobian in end configuration
!                !call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tdisp)
!                call ShapeFunc3(nel,xi,eta,zeta,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,tzero)
!                
!                jacinv = 0.0d0
!                temp1 = 0.0d0
!                jacinv = jac
!                call gaussj(jacinv,nds,temp1,ilixo)
!                
!                
!                !Transformation matrices end configuration 
!                Temp33 = 0.0d0
!                Temp33 = matmul(transpose(rconv),jacinv)
!                Trans = 0.0d0
!                Trans = transpose(rconv)
!                
!                !Natural to local transformation matrix for end configuration
!                call TransformationMat3D(Temp33,TCL)
!                
!                temp6 = 0.0d0
!                TLC = TCL
!                call gaussj(TLC,6,temp6,ilixo)
!                
!                !Global to local transformation matrix for end configuration
!                call TransformationMat3D(Trans,TGL)
!                
!                temp6 = 0.0d0
!                TLG = TGL
!                call gaussj(TLG,6,temp6,ilixo)
!                
!                !Global to natural transformation matrix for end configuration
!                Mat33 = jac !matmul(transpose(rconv),jac)
!                call TransformationMat3D(mat33,TGC)
!
!                !Update Weight factor
!                gwt = we(i)*wn(j)*wc(k)*detj
!                
!                continue
!                
!                
!                
!                
!                
!                
!                
!                
!                !Compatible strain-displacement matrix in global frame
!                Bglob = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
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
!                end do
!                
!                !BNatural = matmul(TLC,Bglob)
!                BCN = matmul(TGC,Bglob)
!                
!                !ANS Implementation I
!                A = sqrt(1.0d0/3.0d0)
!                B = sqrt(3.0d0/5.0d0)
!                
!                !-------------------------------------------------------------------
!                !-------------------------------------------------------------------
!                !Coordinates
!                coordANS = 0.0d0
!                coordANS(7:12,3) = -1.0d0
!                coordANS(13:18,3) = 1.0d0
!                
!                coordANS(1,1)  =  1.0; coordANS(1,2)  =  1.0;
!                coordANS(2,1)  = -1.0; coordANS(2,2)  =  1.0; 
!                coordANS(3,1)  =  1.0; coordANS(3,2)  =  0.0; 
!                coordANS(4,1)  = -1.0; coordANS(4,2)  =  0.0;
!                coordANS(5,1)  =  1.0; coordANS(5,2)  = -1.0; 
!                coordANS(6,1)  = -1.0; coordANS(6,2)  = -1.0;
!                coordANS(7,1)  =  1.0; coordANS(7,2)  =  1.0;
!                coordANS(8,1)  = -1.0; coordANS(8,2)  =  1.0; 
!                coordANS(9,1)  =  1.0; coordANS(9,2)  =  0.0; 
!                coordANS(10,1)  = -1.0; coordANS(10,2)  =  0.0;
!                coordANS(11,1)  =  1.0; coordANS(11,2)  = -1.0; 
!                coordANS(12,1)  = -1.0; coordANS(12,2)  = -1.0;
!                coordANS(13,1)  =  1.0; coordANS(13,2)  =  1.0;
!                coordANS(14,1)  = -1.0; coordANS(14,2)  =  1.0; 
!                coordANS(15,1)  =  1.0; coordANS(15,2)  =  0.0; 
!                coordANS(16,1)  = -1.0; coordANS(16,2)  =  0.0;
!                coordANS(17,1)  =  1.0; coordANS(17,2)  = -1.0; 
!                coordANS(18,1)  = -1.0; coordANS(18,2)  = -1.0;
!                
!                coordANS(25:30,3) = -1.0d0
!                coordANS(31:36,3) = 1.0d0
!                
!                coordANS(19,1)  =  1.0; coordANS(19,2)  =  1.0;
!                coordANS(20,1)  =  0.0; coordANS(20,2)  =  1.0; 
!                coordANS(21,1)  = -1.0; coordANS(21,2)  =  1.0; 
!                coordANS(22,1)  =  1.0; coordANS(22,2)  = -1.0;
!                coordANS(23,1)  =  0.0; coordANS(23,2)  = -1.0; 
!                coordANS(24,1)  = -1.0; coordANS(24,2)  = -1.0;
!                
!                coordANS(25,1)  =  1.0; coordANS(25,2)  =  1.0;
!                coordANS(26,1)  =  0.0; coordANS(26,2)  =  1.0; 
!                coordANS(27,1)  = -1.0; coordANS(27,2)  =  1.0; 
!                coordANS(28,1)  =  1.0; coordANS(28,2)  = -1.0;
!                coordANS(29,1)  =  0.0; coordANS(29,2)  = -1.0; 
!                coordANS(30,1)  = -1.0; coordANS(30,2)  = -1.0;
!                
!                coordANS(31,1)  =  1.0; coordANS(31,2)  =  1.0;
!                coordANS(32,1)  =  0.0; coordANS(32,2)  =  1.0; 
!                coordANS(33,1)  = -1.0; coordANS(33,2)  =  1.0; 
!                coordANS(34,1)  =  1.0; coordANS(34,2)  = -1.0;
!                coordANS(35,1)  =  0.0; coordANS(35,2)  = -1.0; 
!                coordANS(36,1)  = -1.0; coordANS(36,2)  = -1.0;
!                !-------------------------------------------------------------------
!                !-------------------------------------------------------------------
!                
!
!                !Shape Functions
!                shpANS = 0.0d0
!                
!                !-------------------------------------------------------------------
!!                shpANS(1)  = eta/(4.0d0*B)*(eta/B+1.0d0)*(1.0d0+xi/A) 
!!                shpANS(2)  = eta/(4.0d0*B)*(eta/B+1.0d0)*(1.0d0-xi/A)
!!                shpANS(3)  = 1.0d0/2.0d0*(1.0d0-eta*eta/(B*B))*(1.0d0+xi/A)
!!                shpANS(4)  = 1.0d0/2.0d0*(1.0d0-eta*eta/(B*B))*(1.0d0-xi/A)
!!                shpANS(5)  = eta/(4.0d0*B)*(eta/B-1.0d0)*(1.0d0+xi/A) 
!!                shpANS(6)  = eta/(4.0d0*B)*(eta/B-1.0d0)*(1.0d0-xi/A)
!!                
!!!                shpANS(6)  = eta/(4.0d0*B)*(eta/B+1.0d0)*(1.0d0+xi/A) 
!!!                shpANS(5)  = eta/(4.0d0*B)*(eta/B+1.0d0)*(1.0d0-xi/A)
!!!                shpANS(4)  = 1.0d0/2.0d0*(1.0d0-eta*eta/(B*B))*(1.0d0+xi/A)
!!!                shpANS(3)  = 1.0d0/2.0d0*(1.0d0-eta*eta/(B*B))*(1.0d0-xi/A)
!!!                shpANS(2)  = eta/(4.0d0*B)*(eta/B-1.0d0)*(1.0d0+xi/A) 
!!!                shpANS(1)  = eta/(4.0d0*B)*(eta/B-1.0d0)*(1.0d0-xi/A)
!!                !-------------------------------------------------------------------
!!                
!!                !-------------------------------------------------------------------
!!                shpANS(7)  = xi/(4.0d0*B)*(xi/B+1.0d0)*(1.0d0+eta/A)
!!                shpANS(8)  = xi/(4.0d0*B)*(xi/B+1.0d0)*(1.0d0-eta/A)
!!                shpANS(9)  = 1.0d0/2.0d0*(1.0d0-xi*xi/(B*B))*(1.0d0+eta/A)
!!                shpANS(10) = 1.0d0/2.0d0*(1.0d0-xi*xi/(B*B))*(1.0d0-eta/A)
!!                shpANS(11) = xi/(4.0d0*B)*(xi/B-1.0d0)*(1.0d0+eta/A)
!!                shpANS(12) = xi/(4.0d0*B)*(xi/B-1.0d0)*(1.0d0-eta/A)
!!                
!!                shpANS(7)  = xi/(4.0d0*B)*(xi/B+1.0d0)*(1.0d0+eta/A)
!!                shpANS(8)  = 1.0d0/2.0d0*(1.0d0-xi*xi/(B*B))*(1.0d0+eta/A)
!!                shpANS(9)  = xi/(4.0d0*B)*(xi/B-1.0d0)*(1.0d0+eta/A)
!!                shpANS(10) = xi/(4.0d0*B)*(xi/B+1.0d0)*(1.0d0-eta/A)
!!                shpANS(11) = 1.0d0/2.0d0*(1.0d0-xi*xi/(B*B))*(1.0d0-eta/A)
!!                shpANS(12) = xi/(4.0d0*B)*(xi/B-1.0d0)*(1.0d0-eta/A)
!!                !-------------------------------------------------------------------
!!                
!!                !-------------------------------------------------------------------
!!                shpANS(13) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0+eta/A)
!!                shpANS(14) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0-eta/A)
!!                shpANS(15) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0+eta/A)
!!                shpANS(16) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0-eta/A)
!!                
!!                shpANS(13) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0+eta/A)
!!                shpANS(14) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0+eta/A)
!!                shpANS(15) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0-eta/A)
!!                shpANS(16) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0-eta/A)
!!                
!!                shpANS(17) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0+zeta/A);
!!                shpANS(18) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0+zeta/A);
!!                shpANS(19) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0-zeta/A);
!!                shpANS(20) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0-zeta/A);
!!                
!!                shpANS(21) = 1.0d0/4.0d0*(1.0d0+eta/A)*(1.0d0+zeta/A);
!!                shpANS(22) = 1.0d0/4.0d0*(1.0d0-eta/A)*(1.0d0+zeta/A);
!!                shpANS(23) = 1.0d0/4.0d0*(1.0d0+eta/A)*(1.0d0-zeta/A);
!!                shpANS(24) = 1.0d0/4.0d0*(1.0d0-eta/A)*(1.0d0-zeta/A);
!
!                shpANS = 0.0d0
!
!                shpANS(1)  = R(18);
!                shpANS(2)  = R(16);
!                shpANS(3)  = R(15);
!                shpANS(4)  = R(13);
!                shpANS(5)  = R(12);
!                shpANS(6)  = R(10);
!                shpANS(7)  = R(27);
!                shpANS(8)  = R(25);
!                shpANS(9)  = R(24);
!                shpANS(10)  = R(22);
!                shpANS(11)  = R(21);
!                shpANS(12)  = R(19);
!                shpANS(13)  = R(9);
!                shpANS(14)  = R(7);
!                shpANS(15)  = R(6);
!                shpANS(16)  = R(4);
!                shpANS(17)  = R(3);
!                shpANS(18)  = R(1);
!                
!                shpANS(19)  = R(10);
!                shpANS(20)  = R(11);
!                shpANS(21)  = R(12);
!                shpANS(22)  = R(16);
!                shpANS(23)  = R(17);
!                shpANS(24)  = R(18);
!                shpANS(25)  = R(19);
!                shpANS(26)  = R(20);
!                shpANS(27)  = R(21);
!                shpANS(28)  = R(25);
!                shpANS(29)  = R(26);
!                shpANS(30)  = R(27);
!                shpANS(31)  = R(9);
!                shpANS(32)  = R(8);
!                shpANS(33)  = R(7);
!                shpANS(34)  = R(3);
!                shpANS(35)  = R(2);
!                shpANS(36)  = R(1);
!                !-------------------------------------------------------------------
!                
!                BNat = 0.0d0
!                
!                !xi-xi and eta-eta components
!!                do k1=1,6
!!                    Bnat(1,:) = Bnat(1,:) +   shpANS(k1)*BCN(1,:)
!!                    Bnat(2,:) = Bnat(2,:) + shpANS(k1+6)*BCN(2,:)
!!                end do
!                
!                BNat(1,:) = BCN(1,:)
!                BNat(2,:) = BCN(2,:)
!                BNat(3,:) = BCN(3,:)
!                
!                !xi-eta components
!!                do k1=13,16
!!                    Bnat(4,:) = Bnat(4,:) + shpANS(k1)*BCN(4,:)
!!                end do
!
!                BNat(4,:) = BCN(4,:)
!                !BNat(5,:) = BCN(5,:)
!                !BNat(6,:) = BCN(6,:)
!                
!                
!                !xi-zeta components
!                do k1=1,6
!                    call ShapeFunc3(nel,coordANS(k1,1),coordANS(k1,2),coordANS(k1,3),RA,dRdxA,dRdxiA,dRdxiiA,detjA,jacA,dxdxiA,tzero)
!                    
!                    BN2 = 0.0d0
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        BN2(1,k2*3-2) = dRdxA(k2,1)
!                        BN2(2,k2*3-1) = dRdxA(k2,2)
!                        BN2(3,k2*3  ) = dRdxA(k2,3)
!                    
!                        BN2(4,k2*3-2) = dRdxA(k2,2)
!                        BN2(4,k2*3-1) = dRdxA(k2,1)
!                        BN2(5,k2*3-2) = dRdxA(k2,3)
!                        BN2(5,k2*3  ) = dRdxA(k2,1)
!                        BN2(6,k2*3-1) = dRdxA(k2,3)
!                        BN2(6,k2*3  ) = dRdxA(k2,2)
!                    end do
!                
!                    BN2 = matmul(TGC,BN2)
!                    
!                    Bnat(5,:) = Bnat(5,:) + shpANS(k1)*BN2(5,:);
!
!                end do
!                
!                !xi-zeta components
!                do k1=19,24
!                    call ShapeFunc3(nel,coordANS(k1,1),coordANS(k1,2),coordANS(k1,3),RA,dRdxA,dRdxiA,dRdxiiA,detjA,jacA,dxdxiA,tzero)
!                    
!                    BN2 = 0.0d0
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        BN2(1,k2*3-2) = dRdxA(k2,1)
!                        BN2(2,k2*3-1) = dRdxA(k2,2)
!                        BN2(3,k2*3  ) = dRdxA(k2,3)
!                    
!                        BN2(4,k2*3-2) = dRdxA(k2,2)
!                        BN2(4,k2*3-1) = dRdxA(k2,1)
!                        BN2(5,k2*3-2) = dRdxA(k2,3)
!                        BN2(5,k2*3  ) = dRdxA(k2,1)
!                        BN2(6,k2*3-1) = dRdxA(k2,3)
!                        BN2(6,k2*3  ) = dRdxA(k2,2)
!                    end do
!                
!                    BN2 = matmul(TGC,BN2)
!                    
!                    Bnat(6,:) = Bnat(6,:) + shpANS(k1)*BN2(6,:);
!                    
!                end do
!
!                !Local strain-displacement matrix at end configuration
!                Bloc=0.0d0
!                Bloc=matmul(TGL,Bglob)
!                
!                Bloc=matmul(TCL,BNat)
!                
!                continue
!                
!                
!!                jacinv = 0.0d0
!!                temp1 = 0.0d0
!!                jacinv = dxdxi
!!                call gaussj(jacinv,nds,temp1,ilixo)
!!                
!!                !Transformation matrices mid configuration 
!!                Temp33 = 0.0d0
!!                Temp33 = matmul(transpose(rconv),jacinv)
!!                call TransformationMat3D(Temp33,TCL)
!!                
!!                do k1=1,(p+1)*(q+1)*(w+1)
!!                    BNat(1,k1*3-2) = dRdxii(k1,1)
!!                    BNat(2,k1*3-1) = dRdxii(k1,2)
!!                    BNat(3,k1*3  ) = dRdxii(k1,3)
!!                    
!!                    BNat(4,k1*3-2) = dRdxii(k1,2)
!!                    BNat(4,k1*3-1) = dRdxii(k1,1)
!!                end do
!!                
!!                !Bloc2=matmul(TGL,Bloc2)
!!                
!!                
!!                Bloc(5,:) = Bloc2(5,:)
!!                Bloc(6,:) = Bloc2(6,:)
!
!                !Bloc = Bloc2
!
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
!    continue
!
!end subroutine ElHex27ANS

!subroutine ElHex27ANS(nel,Ke,el_ddisp,Finte,inc)
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
!    integer(4)::ii,jj,kk,k2,k3
!    real(8)::CA,CB,detjA
!    real(8),dimension(nds,nds)::jacA,dxdxiA,dxdxi
!    real(8),dimension(npi,nds)::NatPoints
!    real(8),dimension(6)::ShpNat
!    real(8),dimension(npi)::Ra
!    real(8),dimension(npi,3)::DerShpNat,dRadx,dRAdxi,dRAdxii
!    real(8),dimension(npi,3)::dRdxi,dRdxii
!    real(8),dimension(22,nds)::TyPt
!    real(8)::xib,etab,zetab
!    
!    real(8),dimension(p+1)::Nb
!    real(8),dimension(q+1)::Mb
!    real(8),dimension(w+1)::Lb
!    real(8),dimension(p+1)::dNbdxi
!    real(8),dimension(q+1)::dMbdeta
!    real(8),dimension(w+1)::dLbdzeta
!    
!    real(8),dimension(2*(p+1)+p)::ub
!    real(8),dimension(2*(q+1)+q)::vb
!    real(8),dimension(2*(w+1)+w)::wb
!    
!    real(8),dimension(2*(p+1))::ubr
!    real(8),dimension(2*(q+1))::vbr
!    
!    real(8),dimension(nds*2,nds*2)::TGC
!    
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BNat,BANS,BNG,BNGold
!    
!    real(8)::shp,suma
!    
!    ShpNat = 0.0d0
!    derShpNat = 0.0d0
!    
!    suma = 0.0d0
!    
!    do k1=1,p+1
!        ubr(k1) = 0.0d0
!        ubr(k1+p+1) = 1.0d0
!    end do    
!    
!    do k1=1,q+1
!        vbr(k1) = 0.0d0
!        vbr(k1+q+1) = 1.0d0
!    end do 
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
!    KNLG = 0.0d0
!    BNL = 0.0d0
!    MSTR = 0.0d0
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
!                call ShapeFunc3(nel,xi,eta,zeta,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,tzero)
!                
!                jacinv = 0.0d0
!                temp1 = 0.0d0
!                jacinv = jac
!                call gaussj(jacinv,nds,temp1,ilixo)
!                
!                call TransformationMat3D(jac,TGC)
!                call TransformationMat3D(jacinv,TCL)
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
!                !Bloc=matmul(TGL,Bglob)
!                
!                
!                BNG = matmul(TGC,Bglob)
!                
!                BNGold = BNG
!                
!                !Let's Try the ANS ideia ---------------------------------------------------------------
!                CA = dsqrt(1.0d0/3.0d0)
!                CB = dsqrt(3.0d0/5.0d0)
!                
!                !Tying Points Coordinates
!                TyPt = 0.0d0
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
!                TyPt(13,1)  =   CA; TyPt(13,2)  =   CA
!                TyPt(14,1) =    CA; TyPt(14,2) =   -CA
!                TyPt(15,1) =   -CA; TyPt(15,2) =    CA
!                TyPt(16,1) =   -CA; TyPt(16,2) =   -CA
!                
!                !BNG(1,:) = 0.0d0
!                !BNG(2,:) = 0.0d0
!                !BNG(3,:) = 0.0d0
!                
!                !BNG(4,:) = 0.0d0
!                BNG(5,:) = 0.0d0
!                BNG(6,:) = 0.0d0
!                
!                ni = INN(IEN(nel,1),1)
!                nj = INN(IEN(nel,1),2)
!                nk = INN(IEN(nel,1),3)
!    
!                !Calculate location of the integration point in the parametric space (xi,eta) based on the coordinates of the parent element (xii,etai)
!                xib   = ((u_knot(ni+1) - u_knot(ni))*xi   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                etab  = ((v_knot(nj+1) - v_knot(nj))*eta  + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                zetab = ((w_knot(nk+1) - w_knot(nk))*zeta + (w_knot(nk+1) + w_knot(nk)))/2.0d0
!                
!
!                !TP1 -------------------------------------------
!                ub = u_knot(ni+1)
!                vb = v_knot(nj+1)
!                
!                ub(1) = 0.0d0; vb(1) = 0.0d0
!                ub(2) = 0.0d0; vb(2) = 0.0d0
!                ub(3) = 0.0d0; vb(3) = 0.0d0
!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(1,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(1,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(1,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(1,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                call BSplineBasisAndDeriv(5,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(5,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)
!                
!                
!                call ShapeFunc3(nel,TyPt(1,1),TyPt(1,2),TyPt(1,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    BNat(1,k2*3-2) = dRAdx(k2,1)
!                    BNat(2,k2*3-1) = dRAdx(k2,2)
!                    BNat(3,k2*3  ) = dRAdx(k2,3)
!                    BNat(4,k2*3-2) = dRAdx(k2,2)
!                    BNat(4,k2*3-1) = dRAdx(k2,1)
!                    BNat(5,k2*3-2) = dRAdx(k2,3)
!                    BNat(5,k2*3  ) = dRAdx(k2,1)
!                    BNat(6,k2*3-1) = dRAdx(k2,3)
!                    BNat(6,k2*3  ) = dRAdx(k2,2)
!                end do
!
!                BNat = matmul(TGC,BNat)
!                
!                !BNG(1,:) = BNG(1,:) + shp*BNGold(1,:)
!                BNG(5,:) = BNG(5,:) + shp*BNat(5,:)
!                
!                
!                !TP2 -------------------------------------------
!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(2,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(2,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(2,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(2,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                call BSplineBasisAndDeriv(5,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(5,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)
!                
!                call ShapeFunc3(nel,TyPt(2,1),TyPt(2,2),TyPt(2,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    BNat(1,k2*3-2) = dRAdx(k2,1)
!                    BNat(2,k2*3-1) = dRAdx(k2,2)
!                    BNat(3,k2*3  ) = dRAdx(k2,3)
!                    BNat(4,k2*3-2) = dRAdx(k2,2)
!                    BNat(4,k2*3-1) = dRAdx(k2,1)
!                    BNat(5,k2*3-2) = dRAdx(k2,3)
!                    BNat(5,k2*3  ) = dRAdx(k2,1)
!                    BNat(6,k2*3-1) = dRAdx(k2,3)
!                    BNat(6,k2*3  ) = dRAdx(k2,2)
!                end do
!                
!                BNat = matmul(TGC,BNat)
!                
!                !BNG(1,:) = BNG(1,:) + shp*BNGold(1,:)
!                BNG(5,:) = BNG(5,:) + shp*BNat(5,:)
!                
!                
!                !TP3 -------------------------------------------
!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(3,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(3,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(3,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(3,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                call BSplineBasisAndDeriv(5,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(5,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(3,1),TyPt(3,2),TyPt(3,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    BNat(1,k2*3-2) = dRAdx(k2,1)
!                    BNat(2,k2*3-1) = dRAdx(k2,2)
!                    BNat(3,k2*3  ) = dRAdx(k2,3)
!                    BNat(4,k2*3-2) = dRAdx(k2,2)
!                    BNat(4,k2*3-1) = dRAdx(k2,1)
!                    BNat(5,k2*3-2) = dRAdx(k2,3)
!                    BNat(5,k2*3  ) = dRAdx(k2,1)
!                    BNat(6,k2*3-1) = dRAdx(k2,3)
!                    BNat(6,k2*3  ) = dRAdx(k2,2)
!                end do
!                
!                BNat = matmul(TGC,BNat)
!                
!                !BNG(1,:) = BNG(1,:) + shp*BNGold(1,:)
!                BNG(5,:) = BNG(5,:) + shp*BNat(5,:)
!                
!                
!                !TP4 -------------------------------------------
!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(4,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(4,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(4,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(4,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                call BSplineBasisAndDeriv(5,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(5,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(4,1),TyPt(4,2),TyPt(4,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    BNat(1,k2*3-2) = dRAdx(k2,1)
!                    BNat(2,k2*3-1) = dRAdx(k2,2)
!                    BNat(3,k2*3  ) = dRAdx(k2,3)
!                    BNat(4,k2*3-2) = dRAdx(k2,2)
!                    BNat(4,k2*3-1) = dRAdx(k2,1)
!                    BNat(5,k2*3-2) = dRAdx(k2,3)
!                    BNat(5,k2*3  ) = dRAdx(k2,1)
!                    BNat(6,k2*3-1) = dRAdx(k2,3)
!                    BNat(6,k2*3  ) = dRAdx(k2,2)
!                end do
!                
!                BNat = matmul(TGC,BNat)
!                
!                !BNG(1,:) = BNG(1,:) + shp*BNGold(1,:)
!                BNG(5,:) = BNG(5,:) + shp*BNat(5,:)
!                
!                
!                !TP5 -------------------------------------------
!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(5,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(5,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(5,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(5,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                call BSplineBasisAndDeriv(5,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(5,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(5,1),TyPt(5,2),TyPt(5,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    BNat(1,k2*3-2) = dRAdx(k2,1)
!                    BNat(2,k2*3-1) = dRAdx(k2,2)
!                    BNat(3,k2*3  ) = dRAdx(k2,3)
!                    BNat(4,k2*3-2) = dRAdx(k2,2)
!                    BNat(4,k2*3-1) = dRAdx(k2,1)
!                    BNat(5,k2*3-2) = dRAdx(k2,3)
!                    BNat(5,k2*3  ) = dRAdx(k2,1)
!                    BNat(6,k2*3-1) = dRAdx(k2,3)
!                    BNat(6,k2*3  ) = dRAdx(k2,2)
!                end do
!                
!                BNat = matmul(TGC,BNat)
!                
!                !BNG(1,:) = BNG(1,:) + shp*BNGold(1,:)
!                BNG(5,:) = BNG(5,:) + shp*BNat(5,:)
!                
!                
!                !TP6 -------------------------------------------
!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(6,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(6,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(6,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(6,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                call BSplineBasisAndDeriv(5,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(5,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(6,1),TyPt(6,2),TyPt(6,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    BNat(1,k2*3-2) = dRAdx(k2,1)
!                    BNat(2,k2*3-1) = dRAdx(k2,2)
!                    BNat(3,k2*3  ) = dRAdx(k2,3)
!                    BNat(4,k2*3-2) = dRAdx(k2,2)
!                    BNat(4,k2*3-1) = dRAdx(k2,1)
!                    BNat(5,k2*3-2) = dRAdx(k2,3)
!                    BNat(5,k2*3  ) = dRAdx(k2,1)
!                    BNat(6,k2*3-1) = dRAdx(k2,3)
!                    BNat(6,k2*3  ) = dRAdx(k2,2)
!                end do
!                
!                BNat = matmul(TGC,BNat)
!                
!                !BNG(1,:) = BNG(1,:) + shp*BNGold(1,:)
!                BNG(5,:) = BNG(5,:) + shp*BNat(5,:)
!                
!                
!                !TP7 -------------------------------------------
!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(7,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(7,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(7,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(7,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                call BSplineBasisAndDeriv(5,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(5,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(7,1),TyPt(7,2),TyPt(7,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    BNat(1,k2*3-2) = dRAdx(k2,1)
!                    BNat(2,k2*3-1) = dRAdx(k2,2)
!                    BNat(3,k2*3  ) = dRAdx(k2,3)
!                    BNat(4,k2*3-2) = dRAdx(k2,2)
!                    BNat(4,k2*3-1) = dRAdx(k2,1)
!                    BNat(5,k2*3-2) = dRAdx(k2,3)
!                    BNat(5,k2*3  ) = dRAdx(k2,1)
!                    BNat(6,k2*3-1) = dRAdx(k2,3)
!                    BNat(6,k2*3  ) = dRAdx(k2,2)
!                end do
!                
!                BNat = matmul(TGC,BNat)
!                
!                !BNG(2,:) = BNG(2,:) + shp*BNGold(2,:)
!                BNG(6,:) = BNG(6,:) + shp*BNat(6,:)
!                
!                
!                !TP8 -------------------------------------------
!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(8,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(8,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(8,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(8,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                call BSplineBasisAndDeriv(5,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(5,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Mb(3)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(8,1),TyPt(8,2),TyPt(8,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    BNat(1,k2*3-2) = dRAdx(k2,1)
!                    BNat(2,k2*3-1) = dRAdx(k2,2)
!                    BNat(3,k2*3  ) = dRAdx(k2,3)
!                    BNat(4,k2*3-2) = dRAdx(k2,2)
!                    BNat(4,k2*3-1) = dRAdx(k2,1)
!                    BNat(5,k2*3-2) = dRAdx(k2,3)
!                    BNat(5,k2*3  ) = dRAdx(k2,1)
!                    BNat(6,k2*3-1) = dRAdx(k2,3)
!                    BNat(6,k2*3  ) = dRAdx(k2,2)
!                end do
!                
!                BNat = matmul(TGC,BNat)
!                
!                !BNG(2,:) = BNG(2,:) + shp*BNGold(2,:)
!                BNG(6,:) = BNG(6,:) + shp*BNat(6,:)
!                
!                
!                !TP9 -------------------------------------------
!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(9,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(9,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(9,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(9,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                call BSplineBasisAndDeriv(5,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(5,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Mb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(9,1),TyPt(9,2),TyPt(9,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    BNat(1,k2*3-2) = dRAdx(k2,1)
!                    BNat(2,k2*3-1) = dRAdx(k2,2)
!                    BNat(3,k2*3  ) = dRAdx(k2,3)
!                    BNat(4,k2*3-2) = dRAdx(k2,2)
!                    BNat(4,k2*3-1) = dRAdx(k2,1)
!                    BNat(5,k2*3-2) = dRAdx(k2,3)
!                    BNat(5,k2*3  ) = dRAdx(k2,1)
!                    BNat(6,k2*3-1) = dRAdx(k2,3)
!                    BNat(6,k2*3  ) = dRAdx(k2,2)
!                end do
!                
!                BNat = matmul(TGC,BNat)
!                
!                !BNG(2,:) = BNG(2,:) + shp*BNGold(2,:)
!                BNG(6,:) = BNG(6,:) + shp*BNat(6,:)
!                
!                
!                !TP10 -------------------------------------------
!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(10,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(10,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(10,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(10,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                call BSplineBasisAndDeriv(5,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(5,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Mb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(10,1),TyPt(10,2),TyPt(10,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    BNat(1,k2*3-2) = dRAdx(k2,1)
!                    BNat(2,k2*3-1) = dRAdx(k2,2)
!                    BNat(3,k2*3  ) = dRAdx(k2,3)
!                    BNat(4,k2*3-2) = dRAdx(k2,2)
!                    BNat(4,k2*3-1) = dRAdx(k2,1)
!                    BNat(5,k2*3-2) = dRAdx(k2,3)
!                    BNat(5,k2*3  ) = dRAdx(k2,1)
!                    BNat(6,k2*3-1) = dRAdx(k2,3)
!                    BNat(6,k2*3  ) = dRAdx(k2,2)
!                end do
!                
!                BNat = matmul(TGC,BNat)
!                
!                !BNG(2,:) = BNG(2,:) + shp*BNGold(2,:)
!                BNG(6,:) = BNG(6,:) + shp*BNat(6,:)
!                
!                
!                !TP11 -------------------------------------------
!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(11,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(11,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(11,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(11,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                call BSplineBasisAndDeriv(5,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(5,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp =       Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp =       Mb(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp =       Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(11,1),TyPt(11,2),TyPt(11,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    BNat(1,k2*3-2) = dRAdx(k2,1)
!                    BNat(2,k2*3-1) = dRAdx(k2,2)
!                    BNat(3,k2*3  ) = dRAdx(k2,3)
!                    BNat(4,k2*3-2) = dRAdx(k2,2)
!                    BNat(4,k2*3-1) = dRAdx(k2,1)
!                    BNat(5,k2*3-2) = dRAdx(k2,3)
!                    BNat(5,k2*3  ) = dRAdx(k2,1)
!                    BNat(6,k2*3-1) = dRAdx(k2,3)
!                    BNat(6,k2*3  ) = dRAdx(k2,2)
!                end do
!                
!                BNat = matmul(TGC,BNat)
!                
!                !BNG(2,:) = BNG(2,:) + shp*BNGold(2,:)
!                BNG(6,:) = BNG(6,:) + shp*BNat(6,:)
!                
!                
!                !TP12 -------------------------------------------
!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(12,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(12,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(12,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(12,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                call BSplineBasisAndDeriv(5,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(5,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp =       Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp =       Mb(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp =       Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(12,1),TyPt(12,2),TyPt(12,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    BNat(1,k2*3-2) = dRAdx(k2,1)
!                    BNat(2,k2*3-1) = dRAdx(k2,2)
!                    BNat(3,k2*3  ) = dRAdx(k2,3)
!                    BNat(4,k2*3-2) = dRAdx(k2,2)
!                    BNat(4,k2*3-1) = dRAdx(k2,1)
!                    BNat(5,k2*3-2) = dRAdx(k2,3)
!                    BNat(5,k2*3  ) = dRAdx(k2,1)
!                    BNat(6,k2*3-1) = dRAdx(k2,3)
!                    BNat(6,k2*3  ) = dRAdx(k2,2)
!                end do
!                
!                BNat = matmul(TGC,BNat)
!                
!                !BNG(2,:) = BNG(2,:) + shp*BNGold(2,:)
!                BNG(6,:) = BNG(6,:) + shp*BNat(6,:)
!                
!                
!                
!!                !TP13 -------------------------------------------
!!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(13,1) + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(13,1) + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(13,2) + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(13,2) + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                
!!                call BSplineBasisAndDeriv(5,p,xib,ub,Nb,dNbdxi)
!!                call BSplineBasisAndDeriv(5,q,etab,vb,Mb,dMbdeta)    
!!                
!!                
!!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                call ShapeFunc3(nel,TyPt(13,1),TyPt(13,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!!                    
!!                BNat = 0.0d0
!!                do k2=1,(p+1)*(q+1)*(w+1)
!!                    BNat(1,k2*3-2) = dRAdx(k2,1)
!!                    BNat(2,k2*3-1) = dRAdx(k2,2)
!!                    BNat(3,k2*3  ) = dRAdx(k2,3)
!!                    BNat(4,k2*3-2) = dRAdx(k2,2)
!!                    BNat(4,k2*3-1) = dRAdx(k2,1)
!!                    BNat(5,k2*3-2) = dRAdx(k2,3)
!!                    BNat(5,k2*3  ) = dRAdx(k2,1)
!!                    BNat(6,k2*3-1) = dRAdx(k2,3)
!!                    BNat(6,k2*3  ) = dRAdx(k2,2)
!!                end do
!!                
!!                BNat = matmul(TGC,BNat)
!!                BNG(4,:) = BNG(4,:) + shp*BNat(4,:)
!!                !BNG(4,:) = BNG(4,:) + shp*BNGold(4,:)
!!                
!!                suma = suma + shp
!!                
!!                
!!                !TP14 -------------------------------------------
!!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(14,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(14,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(14,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(14,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                
!!                call BSplineBasisAndDeriv(5,p,xib,ub,Nb,dNbdxi)
!!                call BSplineBasisAndDeriv(5,q,etab,vb,Mb,dMbdeta)    
!!                
!!                
!!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                call ShapeFunc3(nel,TyPt(14,1),TyPt(14,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!!                    
!!                BNat = 0.0d0
!!                do k2=1,(p+1)*(q+1)*(w+1)
!!                    BNat(1,k2*3-2) = dRAdx(k2,1)
!!                    BNat(2,k2*3-1) = dRAdx(k2,2)
!!                    BNat(3,k2*3  ) = dRAdx(k2,3)
!!                    BNat(4,k2*3-2) = dRAdx(k2,2)
!!                    BNat(4,k2*3-1) = dRAdx(k2,1)
!!                    BNat(5,k2*3-2) = dRAdx(k2,3)
!!                    BNat(5,k2*3  ) = dRAdx(k2,1)
!!                    BNat(6,k2*3-1) = dRAdx(k2,3)
!!                    BNat(6,k2*3  ) = dRAdx(k2,2)
!!                end do
!!                
!!                BNat = matmul(TGC,BNat)
!!                BNG(4,:) = BNG(4,:) + shp*BNat(4,:)
!!                !BNG(4,:) = BNG(4,:) + shp*BNGold(4,:)
!!                
!!                suma = suma + shp
!!                
!!                !TP15 -------------------------------------------
!!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(15,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(15,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(15,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(15,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                
!!                call BSplineBasisAndDeriv(5,p,xib,ub,Nb,dNbdxi)
!!                call BSplineBasisAndDeriv(5,q,etab,vb,Mb,dMbdeta)    
!!                
!!                
!!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                call ShapeFunc3(nel,TyPt(15,1),TyPt(15,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!!                    
!!                BNat = 0.0d0
!!                do k2=1,(p+1)*(q+1)*(w+1)
!!                    BNat(1,k2*3-2) = dRAdx(k2,1)
!!                    BNat(2,k2*3-1) = dRAdx(k2,2)
!!                    BNat(3,k2*3  ) = dRAdx(k2,3)
!!                    BNat(4,k2*3-2) = dRAdx(k2,2)
!!                    BNat(4,k2*3-1) = dRAdx(k2,1)
!!                    BNat(5,k2*3-2) = dRAdx(k2,3)
!!                    BNat(5,k2*3  ) = dRAdx(k2,1)
!!                    BNat(6,k2*3-1) = dRAdx(k2,3)
!!                    BNat(6,k2*3  ) = dRAdx(k2,2)
!!                end do
!!                
!!                BNat = matmul(TGC,BNat)
!!                BNG(4,:) = BNG(4,:) + shp*BNat(4,:)
!!                !BNG(4,:) = BNG(4,:) + shp*BNGold(4,:)
!!                
!!                suma = suma + shp
!!                
!!                !TP16 -------------------------------------------
!!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(16,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(16,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(16,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(16,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                
!!                call BSplineBasisAndDeriv(5,p,xib,ub,Nb,dNbdxi)
!!                call BSplineBasisAndDeriv(5,q,etab,vb,Mb,dMbdeta)    
!!                
!!                
!!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                call ShapeFunc3(nel,TyPt(16,1),TyPt(16,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!!                    
!!                BNat = 0.0d0
!!                do k2=1,(p+1)*(q+1)*(w+1)
!!                    BNat(1,k2*3-2) = dRAdx(k2,1)
!!                    BNat(2,k2*3-1) = dRAdx(k2,2)
!!                    BNat(3,k2*3  ) = dRAdx(k2,3)
!!                    BNat(4,k2*3-2) = dRAdx(k2,2)
!!                    BNat(4,k2*3-1) = dRAdx(k2,1)
!!                    BNat(5,k2*3-2) = dRAdx(k2,3)
!!                    BNat(5,k2*3  ) = dRAdx(k2,1)
!!                    BNat(6,k2*3-1) = dRAdx(k2,3)
!!                    BNat(6,k2*3  ) = dRAdx(k2,2)
!!                end do
!!                
!!                BNat = matmul(TGC,BNat)
!!                BNG(4,:) = BNG(4,:) + shp*BNat(4,:)
!!                !BNG(4,:) = BNG(4,:) + shp*BNGold(4,:)
!!                
!!                suma = suma + shp
!!                
!!                continue
!                
!
!                
!                
!                Bloc=matmul(TGL,Bglob)
!                
!                Bloc=matmul(TCL,BNG)
!                
!                !BLOC = matmul(TGL,BNG)
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
!    continue
!
!end subroutine ElHex27ANS



!subroutine ElHex27ANS(nel,Ke,el_ddisp,Finte,inc)
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
!    real(8),dimension(22,nds)::TyPt
!    real(8)::xib,etab,zetab
!    
!    real(8),dimension(p+1)::Nb
!    real(8),dimension(q+1)::Mb
!    real(8),dimension(w+1)::Lb
!    real(8),dimension(p+1)::dNbdxi
!    real(8),dimension(q+1)::dMbdeta
!    real(8),dimension(w+1)::dLbdzeta
!    
!    real(8),dimension(ncpx+p+1+p)::ub
!    real(8),dimension(ncpy+q+1+q)::vb
!    real(8),dimension(ncpz+w+1+w)::wb
!    
!    real(8),dimension(2*(p+1))::ubr
!    real(8),dimension(2*(q+1))::vbr
!    
!    real(8),dimension(nds*2,nds*2)::TGC,TGC2
!    
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BNat,BANS,BNG,BNGold
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BConv
!    
!    real(8)::shp,suma,uk,vk,wk,wei,sumw
!    
!    real(8),dimension(3,3)::jact
!    
!    logical::lsp,endcycle
!    
!    ShpNat = 0.0d0
!    derShpNat = 0.0d0
!    
!    suma = 0.0d0
!    
!    do k1=1,p+1
!        ubr(k1) = 0.0d0
!        ubr(k1+p+1) = 1.0d0
!    end do    
!    
!    do k1=1,q+1
!        vbr(k1) = 0.0d0
!        vbr(k1+q+1) = 1.0d0
!    end do 
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
!!    call ShapeFunc3(nel,0.0d0,0.0d0,0.0d0,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,tzero)
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
!    KNLG = 0.0d0
!    BNL = 0.0d0
!    MSTR = 0.0d0
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
!                call ShapeFunc3(nel,xi,eta,zeta,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,tzero)
!                
!                jac = transpose(jac)
!                
!                jacinv = 0.0d0
!                temp1 = 0.0d0
!                jacinv = jac
!                call gaussj(jacinv,nds,temp1,ilixo)
!                
!                
!                call TransformationMat3D(jac,TGC)
!                !call TransformationMat3D_PB(jac,TGC)
!                
!                call TransformationMat3D(jacinv,TCL)
!                
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
!                BConv = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
!                    do k2 =1,nds
!                        BConv(1,(k1-1)*3+k2) = dRdxii(k1,1)*jac(1,k2)
!                        BConv(2,(k1-1)*3+k2) = dRdxii(k1,2)*jac(2,k2)
!                        BConv(3,(k1-1)*3+k2) = dRdxii(k1,3)*jac(3,k2)
!                        
!                        BConv(4,(k1-1)*3+k2) = dRdxii(k1,1)*jac(2,k2) + dRdxii(k1,2)*jac(1,k2)
!                        BConv(5,(k1-1)*3+k2) = dRdxii(k1,1)*jac(3,k2) + dRdxii(k1,3)*jac(1,k2)
!                        BConv(6,(k1-1)*3+k2) = dRdxii(k1,2)*jac(3,k2) + dRdxii(k1,3)*jac(2,k2)
!                    end do
!                end do
!                
!                !Local strain-displacement matrix at end configuration
!                Bloc=0.0d0
!                
!                !BConv=matmul(TGC,Bglob)
!                !BNG = matmul(TGC,Bglob)
!                
!                BNG = BConv
!                
!                BNGold = BNG
!                
!                !Let's Try the ANS ideia ---------------------------------------------------------------
!                CA = dsqrt(1.0d0/3.0d0)
!                CB = dsqrt(3.0d0/5.0d0)
!                
!                !Tying Points Coordinates
!                TyPt = 0.0d0
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
!!                TyPt(13,1) =  1.0d0; TyPt(13,2) =  1.0d0; TyPt(13,3) =  zeta
!!                TyPt(14,1) =  1.0d0; TyPt(14,2) =  0.0d0; TyPt(14,3) =  zeta
!!                TyPt(15,1) =  1.0d0; TyPt(15,2) = -1.0d0; TyPt(15,3) =  zeta
!!                TyPt(16,1) =  0.0d0; TyPt(16,2) =  1.0d0; TyPt(16,3) =  zeta
!!                TyPt(17,1) =  0.0d0; TyPt(17,2) = -1.0d0; TyPt(17,3) =  zeta
!!                TyPt(18,1) = -1.0d0; TyPt(18,2) =  1.0d0; TyPt(18,3) =  zeta
!!                TyPt(19,1) = -1.0d0; TyPt(19,2) =  0.0d0; TyPt(19,3) =  zeta
!!                TyPt(20,1) = -1.0d0; TyPt(20,2) = -1.0d0; TyPt(20,3) =  zeta
!                
!                BNG(1,:) = 0.0d0
!                BNG(2,:) = 0.0d0
!                !BNG(3,:) = 0.0d0
!                
!                !BNG(4,:) = 0.0d0
!                BNG(5,:) = 0.0d0
!                BNG(6,:) = 0.0d0
!                
!                ni = INN(IEN(nel,1),1)
!                nj = INN(IEN(nel,1),2)
!                nk = INN(IEN(nel,1),3)
!    
!                xib   = ((u_knot(ni+1) - u_knot(ni))*xi   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                etab  = ((v_knot(nj+1) - v_knot(nj))*eta  + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                zetab = ((w_knot(nk+1) - w_knot(nk))*zeta + (w_knot(nk+1) + w_knot(nk)))/2.0d0
!                
!
!                !TP1 -------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(1,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(1,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                wk = ((w_knot(nk+1) - w_knot(nk))*TyPt(1,3)   + (w_knot(nk+1) + w_knot(nk)))/2.0d0
!                
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!                call kvecANS(ncpz,w,w_knot,wk,wb)
!
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)
!                call BSplineBasisAndDeriv(ncpz+w,w,zetab,wb,Lb,dLbdzeta)    
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)
!                
!                
!                call ShapeFunc3(nel,TyPt(1,1),TyPt(1,2),TyPt(1,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!
!                BNG(5,:) = BNG(5,:) + shp*BNat(5,:)
!                
!                call ShapeFunc3(nel,TyPt(1,1),TyPt(1,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!
!                BNG(1,:) = BNG(1,:) + shp*BNat(1,:)
!                
!                
!                !TP2 -------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(2,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(2,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)
!                
!                call ShapeFunc3(nel,TyPt(2,1),TyPt(2,2),TyPt(2,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!                
!                BNG(5,:) = BNG(5,:) + shp*BNat(5,:)
!                
!                call ShapeFunc3(nel,TyPt(2,1),TyPt(2,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!
!                BNG(1,:) = BNG(1,:) + shp*BNat(1,:)
!                
!                !TP3 -------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(3,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(3,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(3,1),TyPt(3,2),TyPt(3,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!                
!                BNG(5,:) = BNG(5,:) + shp*BNat(5,:)
!                
!                call ShapeFunc3(nel,TyPt(3,1),TyPt(3,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!
!                BNG(1,:) = BNG(1,:) + shp*BNat(1,:)
!                
!                !TP4 -------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(4,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(4,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(4,1),TyPt(4,2),TyPt(4,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!               
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!                
!                BNG(5,:) = BNG(5,:) + shp*BNat(5,:)
!                
!                call ShapeFunc3(nel,TyPt(4,1),TyPt(4,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!
!                BNG(1,:) = BNG(1,:) + shp*BNat(1,:)
!                
!                !TP5 -------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(5,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(5,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(5,1),TyPt(5,2),TyPt(5,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!                
!                BNG(5,:) = BNG(5,:) + shp*BNat(5,:)
!                
!                call ShapeFunc3(nel,TyPt(5,1),TyPt(5,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!
!                BNG(1,:) = BNG(1,:) + shp*BNat(1,:)
!                
!                !TP6 -------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(6,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(6,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(6,1),TyPt(6,2),TyPt(6,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!                
!                BNG(5,:) = BNG(5,:) + shp*BNat(5,:)
!                
!                call ShapeFunc3(nel,TyPt(6,1),TyPt(6,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!
!                BNG(1,:) = BNG(1,:) + shp*BNat(1,:)
!                
!                !TP7 -------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(7,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(7,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(7,1),TyPt(7,2),TyPt(7,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!                
!                BNG(6,:) = BNG(6,:) + shp*BNat(6,:)
!                
!                call ShapeFunc3(nel,TyPt(7,1),TyPt(7,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!
!                BNG(2,:) = BNG(2,:) + shp*BNat(2,:)
!                
!                !TP8 -------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(8,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(8,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Mb(3)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(8,1),TyPt(8,2),TyPt(8,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!                
!                BNG(6,:) = BNG(6,:) + shp*BNat(6,:)
!                
!                call ShapeFunc3(nel,TyPt(8,1),TyPt(8,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!
!                BNG(2,:) = BNG(2,:) + shp*BNat(2,:)
!                
!                !TP9 -------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(9,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(9,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Mb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(9,1),TyPt(9,2),TyPt(9,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!                
!                BNG(6,:) = BNG(6,:) + shp*BNat(6,:)
!                
!                call ShapeFunc3(nel,TyPt(9,1),TyPt(9,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!
!                BNG(2,:) = BNG(2,:) + shp*BNat(2,:)
!                
!                !TP10 -------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(10,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(10,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Mb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(10,1),TyPt(10,2),TyPt(10,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!                
!                BNG(6,:) = BNG(6,:) + shp*BNat(6,:)
!                
!                call ShapeFunc3(nel,TyPt(10,1),TyPt(10,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!
!                BNG(2,:) = BNG(2,:) + shp*BNat(2,:)
!                
!                !TP11 -------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(11,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(11,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp =       Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp =       Mb(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp =       Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(11,1),TyPt(11,2),TyPt(11,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!                
!                BNG(6,:) = BNG(6,:) + shp*BNat(6,:)
!                
!                call ShapeFunc3(nel,TyPt(11,1),TyPt(11,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!
!                BNG(2,:) = BNG(2,:) + shp*BNat(2,:)
!                
!                !TP12 -------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(12,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(12,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)    
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp =       Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp =       Mb(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp =       Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(12,1),TyPt(12,2),TyPt(12,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!               
!                    
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!                
!                BNG(6,:) = BNG(6,:) + shp*BNat(6,:)
!                
!                call ShapeFunc3(nel,TyPt(12,1),TyPt(12,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                BNat = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    do k3 =1,nds
!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!                        
!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!                    end do
!                end do
!
!                BNG(2,:) = BNG(2,:) + shp*BNat(2,:)
!                
!
!                !TP13 -------------------------------------------
!!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(13,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(13,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                call kvecANS(ncpx,p,u_knot,uk,ub)
!!                call kvecANS(ncpy,q,v_knot,vk,vb)
!!
!!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta) 
!!                
!!                
!!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                call ShapeFunc3(nel,TyPt(13,1),TyPt(13,2),TyPt(13,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!!                jacA = transpose(jacA)
!!                
!!                    
!!                BNat = 0.0d0
!!                do k2=1,(p+1)*(q+1)*(w+1)
!!                    do k3 =1,nds
!!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!!                        
!!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!!                    end do
!!                end do
!!
!!                BNG(4,:) = BNG(4,:) + shp*BNat(4,:)
!!                
!!                
!!                !TP14 -------------------------------------------
!!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(14,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(14,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                call kvecANS(ncpx,p,u_knot,uk,ub)
!!                call kvecANS(ncpy,q,v_knot,vk,vb)
!!
!!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)
!!                
!!                
!!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                call ShapeFunc3(nel,TyPt(14,1),TyPt(14,2),TyPt(14,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!!                jacA = transpose(jacA)
!!                    
!!                BNat = 0.0d0
!!                do k2=1,(p+1)*(q+1)*(w+1)
!!                    do k3 =1,nds
!!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!!                        
!!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!!                    end do
!!                end do
!!                
!!                BNG(4,:) = BNG(4,:) + shp*BNat(4,:)
!!
!!                !TP15 -------------------------------------------
!!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(15,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(15,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                call kvecANS(ncpx,p,u_knot,uk,ub)
!!                call kvecANS(ncpy,q,v_knot,vk,vb)
!!
!!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)
!!                
!!                
!!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                call ShapeFunc3(nel,TyPt(15,1),TyPt(15,2),TyPt(15,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!!                jacA = transpose(jacA)
!!                    
!!                BNat = 0.0d0
!!                do k2=1,(p+1)*(q+1)*(w+1)
!!                    do k3 =1,nds
!!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!!                        
!!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!!                    end do
!!                end do
!!
!!                BNG(4,:) = BNG(4,:) + shp*BNat(4,:)
!!
!!                
!!                !TP16 -------------------------------------------
!!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(16,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(16,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                call kvecANS(ncpx,p,u_knot,uk,ub)
!!                call kvecANS(ncpy,q,v_knot,vk,vb)
!!
!!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)
!!                
!!                
!!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                call ShapeFunc3(nel,TyPt(16,1),TyPt(16,2),TyPt(16,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!!                jacA = transpose(jacA)
!!                    
!!                BNat = 0.0d0
!!                do k2=1,(p+1)*(q+1)*(w+1)
!!                    do k3 =1,nds
!!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!!                        
!!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!!                    end do
!!                end do
!!                
!!                BNG(4,:) = BNG(4,:) + shp*BNat(4,:)
!!                
!!                continue
!                
!                
!                
!                
!                
!                
!                
!                
!                
!                
!                
!                !TP13 -------------------------------------------
!!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(13,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(13,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                call kvecANS(ncpx,p,u_knot,uk,ub)
!!                call kvecANS(ncpy,q,v_knot,vk,vb)
!!
!!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta) 
!!                
!!                
!!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(3)
!!                
!!                call ShapeFunc3(nel,TyPt(13,1),TyPt(13,2),TyPt(13,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!!                jacA = transpose(jacA)
!!                
!!                    
!!                BNat = 0.0d0
!!                do k2=1,(p+1)*(q+1)*(w+1)
!!                    do k3 =1,nds
!!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!!                        
!!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!!                    end do
!!                end do
!!
!!                BNG(1,:) = BNG(1,:) + shp*BNat(1,:)
!!                BNG(2,:) = BNG(2,:) + shp*BNat(2,:)
!!                BNG(4,:) = BNG(4,:) + shp*BNat(4,:)
!!                
!!                
!!                !TP14 -------------------------------------------
!!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(14,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(14,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                call kvecANS(ncpx,p,u_knot,uk,ub)
!!                call kvecANS(ncpy,q,v_knot,vk,vb)
!!
!!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta) 
!!                
!!                
!!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)
!!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)
!!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                
!!                call ShapeFunc3(nel,TyPt(14,1),TyPt(14,2),TyPt(14,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!!                jacA = transpose(jacA)
!!                
!!                    
!!                BNat = 0.0d0
!!                do k2=1,(p+1)*(q+1)*(w+1)
!!                    do k3 =1,nds
!!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!!                        
!!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!!                    end do
!!                end do
!!
!!                BNG(1,:) = BNG(1,:) + shp*BNat(1,:)
!!                
!!                
!!                !TP15 -------------------------------------------
!!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(15,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(15,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                call kvecANS(ncpx,p,u_knot,uk,ub)
!!                call kvecANS(ncpy,q,v_knot,vk,vb)
!!
!!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta) 
!!                
!!                
!!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(1)
!!                
!!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!!                
!!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                
!!                call ShapeFunc3(nel,TyPt(15,1),TyPt(15,2),TyPt(15,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!!                jacA = transpose(jacA)
!!                
!!                    
!!                BNat = 0.0d0
!!                do k2=1,(p+1)*(q+1)*(w+1)
!!                    do k3 =1,nds
!!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!!                        
!!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!!                    end do
!!                end do
!!
!!                BNG(1,:) = BNG(1,:) + shp*BNat(1,:)
!!                BNG(2,:) = BNG(2,:) + shp*BNat(2,:)
!!                BNG(4,:) = BNG(4,:) + shp*BNat(4,:)
!!                
!!                
!!                !TP16 -------------------------------------------
!!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(16,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(16,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                call kvecANS(ncpx,p,u_knot,uk,ub)
!!                call kvecANS(ncpy,q,v_knot,vk,vb)
!!
!!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta) 
!!                
!!                
!!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp =       Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta == 0.0d0) shp =       Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(3)
!!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp =       Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                call ShapeFunc3(nel,TyPt(16,1),TyPt(16,2),TyPt(16,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!!                jacA = transpose(jacA)
!!                
!!                    
!!                BNat = 0.0d0
!!                do k2=1,(p+1)*(q+1)*(w+1)
!!                    do k3 =1,nds
!!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!!                        
!!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!!                    end do
!!                end do
!!
!!                BNG(2,:) = BNG(2,:) + shp*BNat(2,:)
!!                
!!                !TP17 -------------------------------------------
!!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(17,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(17,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                call kvecANS(ncpx,p,u_knot,uk,ub)
!!                call kvecANS(ncpy,q,v_knot,vk,vb)
!!
!!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta) 
!!                
!!                
!!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp =       Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta == 0.0d0) shp =       Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mb(1)
!!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp =       Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                call ShapeFunc3(nel,TyPt(17,1),TyPt(17,2),TyPt(17,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!!                jacA = transpose(jacA)
!!                
!!                    
!!                BNat = 0.0d0
!!                do k2=1,(p+1)*(q+1)*(w+1)
!!                    do k3 =1,nds
!!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!!                        
!!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!!                    end do
!!                end do
!!
!!                BNG(2,:) = BNG(2,:) + shp*BNat(2,:)
!!                
!!                !TP18 -------------------------------------------
!!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(18,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(18,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                call kvecANS(ncpx,p,u_knot,uk,ub)
!!                call kvecANS(ncpy,q,v_knot,vk,vb)
!!
!!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta) 
!!                
!!                
!!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(3)
!!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                call ShapeFunc3(nel,TyPt(18,1),TyPt(18,2),TyPt(18,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!!                jacA = transpose(jacA)
!!                
!!                    
!!                BNat = 0.0d0
!!                do k2=1,(p+1)*(q+1)*(w+1)
!!                    do k3 =1,nds
!!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!!                        
!!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!!                    end do
!!                end do
!!
!!                BNG(1,:) = BNG(1,:) + shp*BNat(1,:)
!!                BNG(2,:) = BNG(2,:) + shp*BNat(2,:)
!!                BNG(4,:) = BNG(4,:) + shp*BNat(4,:)
!!                
!!                !TP19 -------------------------------------------
!!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(19,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(19,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                call kvecANS(ncpx,p,u_knot,uk,ub)
!!                call kvecANS(ncpy,q,v_knot,vk,vb)
!!
!!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta) 
!!                
!!                
!!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(3)
!!                
!!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)
!!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)
!!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)
!!                
!!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                call ShapeFunc3(nel,TyPt(19,1),TyPt(19,2),TyPt(19,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!!                jacA = transpose(jacA)
!!                
!!                    
!!                BNat = 0.0d0
!!                do k2=1,(p+1)*(q+1)*(w+1)
!!                    do k3 =1,nds
!!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!!                        
!!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!!                    end do
!!                end do
!!
!!                BNG(1,:) = BNG(1,:) + shp*BNat(1,:)
!!                
!!                
!!                !TP20 -------------------------------------------
!!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(20,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(20,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!!                call kvecANS(ncpx,p,u_knot,uk,ub)
!!                call kvecANS(ncpy,q,v_knot,vk,vb)
!!
!!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta) 
!!                
!!                
!!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(1)
!!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mb(1)
!!                
!!                call ShapeFunc3(nel,TyPt(20,1),TyPt(20,2),TyPt(20,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!!                jacA = transpose(jacA)
!!                
!!                    
!!                BNat = 0.0d0
!!                do k2=1,(p+1)*(q+1)*(w+1)
!!                    do k3 =1,nds
!!                        BNat(1,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(1,k3)
!!                        BNat(2,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(2,k3)
!!                        BNat(3,(k2-1)*3+k3) = dRAdxii(k2,3)*jacA(3,k3)
!!                        
!!                        BNat(4,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3)
!!                        BNat(5,(k2-1)*3+k3) = dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3)
!!                        BNat(6,(k2-1)*3+k3) = dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3)
!!                    end do
!!                end do
!!
!!                BNG(1,:) = BNG(1,:) + shp*BNat(1,:)
!!                BNG(2,:) = BNG(2,:) + shp*BNat(2,:)
!!                BNG(4,:) = BNG(4,:) + shp*BNat(4,:)
!                
!                
!                
!                
!                !Bloc=matmul(TGL,Bglob)
!                
!                Bloc=matmul(TCL,BNG)
!                
!                
!                
!                !###########################################
!                !Bloc=matmul(TCL,BConv)
!                !###########################################
!                
!                
!                
!                
!                
!                !BLOC = matmul(TGL,BNG)
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
!    continue
!
!end subroutine ElHex27ANS




!!------------------------------------------------------------------------------------------------------------
!!
!! Second order hexahedral element with ANS
!!
!! Input: nel - Element number
!!        el_ddisp - Element displacement
!!        inc - Increment number
!!
!! Output: Ke - Elementar stiffness matrix
!!         Finte - Elementar internal forces
!!------------------------------------------------------------------------------------------------------------
!
!!H2ANS(z)
!
!subroutine ElHex27ANS(nel,Ke,el_ddisp,Finte,inc)
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
!    real(8),dimension(ncpx+p+1+p)::ub
!    real(8),dimension(ncpy+q+1+q)::vb
!    real(8),dimension(ncpz+w+1+w)::wb
!    
!    real(8),dimension(ncpx+p+1+p-2)::ubr
!    real(8),dimension(ncpy+q+1+q-2)::vbr
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
!    logical::lsp,endcycle
!    logical::com1,com2,com3,com4,com5,com6
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
!!    call ShapeFunc3(nel,0.0d0,0.0d0,0.0d0,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,tzero)
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
!                xi = e(i)
!                eta = n(j)
!                zeta = c(k)
!                
!                call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tzero)
!                
!                !Weight factor
!                vol = vol + we(i)*wn(j)*wc(k)*detj
!            end do
!        end do
!    end do
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
!!                Bglob = 0.0d0
!!                do k1=1,(p+1)*(q+1)*(w+1)
!!                    
!!                    Bglob(1,k1*3-2) = dRdx(k1,1)
!!                    Bglob(2,k1*3-1) = dRdx(k1,2)
!!                    Bglob(3,k1*3  ) = dRdx(k1,3)
!!                    
!!                    Bglob(4,k1*3-2) = dRdx(k1,2)
!!                    Bglob(4,k1*3-1) = dRdx(k1,1)
!!                    Bglob(5,k1*3-2) = dRdx(k1,3)
!!                    Bglob(5,k1*3  ) = dRdx(k1,1)
!!                    Bglob(6,k1*3-1) = dRdx(k1,3)
!!                    Bglob(6,k1*3  ) = dRdx(k1,2)
!!                    
!!                end do
!                
!                BConv = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
!                    do k2 =1,nds
!                        BConv(1,(k1-1)*3+k2) = dRdxii(k1,1)*jac(1,k2)
!                        BConv(2,(k1-1)*3+k2) = dRdxii(k1,2)*jac(2,k2)
!                        BConv(3,(k1-1)*3+k2) = dRdxii(k1,3)*jac(3,k2)
!                        
!                        BConv(4,(k1-1)*3+k2) = dRdxii(k1,1)*jac(2,k2) + dRdxii(k1,2)*jac(1,k2)
!                        BConv(5,(k1-1)*3+k2) = dRdxii(k1,1)*jac(3,k2) + dRdxii(k1,3)*jac(1,k2)
!                        BConv(6,(k1-1)*3+k2) = dRdxii(k1,2)*jac(3,k2) + dRdxii(k1,3)*jac(2,k2)
!                    end do
!                end do
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
!                if(com1 == .true.) BNG(1,:) = 0.0d0
!                if(com2 == .true.) BNG(2,:) = 0.0d0
!                !if(com3 == .true.) BNG(3,:) = 0.0d0
!                if(com4 == .true.) BNG(4,:) = 0.0d0
!                if(com5 == .true.) BNG(5,:) = 0.0d0
!                if(com6 == .true.) BNG(6,:) = 0.0d0
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
!                A = CA
!                B = CB
!                
!                shpANS = 0.0d0
!                !-------------------------------------------------------------------
!                shpANS(1)  = eta/(4.0d0*B)*(eta/B+1.0d0)*(1.0d0+xi/A) 
!                shpANS(2)  = eta/(4.0d0*B)*(eta/B+1.0d0)*(1.0d0-xi/A)
!                shpANS(3)  = 1.0d0/2.0d0*(1.0d0-eta*eta/(B*B))*(1.0d0+xi/A)
!                shpANS(4)  = 1.0d0/2.0d0*(1.0d0-eta*eta/(B*B))*(1.0d0-xi/A)
!                shpANS(5)  = eta/(4.0d0*B)*(eta/B-1.0d0)*(1.0d0+xi/A) 
!                shpANS(6)  = eta/(4.0d0*B)*(eta/B-1.0d0)*(1.0d0-xi/A)
!                !-------------------------------------------------------------------
!                !-------------------------------------------------------------------
!                shpANS(7)  = xi/(4.0d0*B)*(xi/B+1.0d0)*(1.0d0+eta/A)
!                shpANS(8)  = 1.0d0/2.0d0*(1.0d0-xi*xi/(B*B))*(1.0d0+eta/A)
!                shpANS(9)  = xi/(4.0d0*B)*(xi/B-1.0d0)*(1.0d0+eta/A)
!                shpANS(10) = xi/(4.0d0*B)*(xi/B+1.0d0)*(1.0d0-eta/A)
!                shpANS(11) = 1.0d0/2.0d0*(1.0d0-xi*xi/(B*B))*(1.0d0-eta/A)
!                shpANS(12) = xi/(4.0d0*B)*(xi/B-1.0d0)*(1.0d0-eta/A)
!                !-------------------------------------------------------------------
!                
!                !-------------------------------------------------------------------
!                shpANS(13) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0+eta/A)
!                shpANS(14) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0+eta/A)
!                shpANS(15) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0-eta/A)
!                shpANS(16) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0-eta/A)
!                
!                
!                !TP1 ---------------------------------------------------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(1,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(1,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!                
!                ubr = 0.0d0
!                do k2=2,ncpx+2*p
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(2)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)
!                
!                call ShapeFunc3(nel,TyPt(1,1),TyPt(1,2),TyPt(1,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + RA(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !shp = shp*shpANS(1)
!                sumw = sumw + shp
!                
!                if(com5==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(5,(k2-1)*3+k3) = BNG(5,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(1,1),TyPt(1,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com1==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(1,(k2-1)*3+k3) = BNG(1,(k2-1)*3+k3) + dRAdxii(k2,1)*jacA(1,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                
!                !TP2 ---------------------------------------------------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(2,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(2,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                ubr = 0.0d0
!                do k2=2,ncpx+2*p
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi) 
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)
!
!                call ShapeFunc3(nel,TyPt(2,1),TyPt(2,2),TyPt(2,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + RA(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !shp = shp*shpANS(2)
!                sumw = sumw + shp
!                    
!                if(com5==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(5,(k2-1)*3+k3) = BNG(5,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(2,1),TyPt(2,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com1==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(1,(k2-1)*3+k3) = BNG(1,(k2-1)*3+k3) + dRAdxii(k2,1)*jacA(1,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                
!                !TP3 ---------------------------------------------------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(3,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(3,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                ubr = 0.0d0
!                do k2=2,ncpx+2*p
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi) 
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(2)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(3,1),TyPt(3,2),TyPt(3,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + RA(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !shp = shp*shpANS(3)
!                sumw = sumw + shp
!                
!                if(com5==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(5,(k2-1)*3+k3) = BNG(5,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(3,1),TyPt(3,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com1==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(1,(k2-1)*3+k3) = BNG(1,(k2-1)*3+k3) + dRAdxii(k2,1)*jacA(1,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                
!                !TP4 ---------------------------------------------------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(4,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(4,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                ubr = 0.0d0
!                do k2=2,ncpx+2*p
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(4,1),TyPt(4,2),TyPt(4,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + RA(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !shp = shp*shpANS(4)
!                sumw = sumw + shp
!                
!                if(com5==.true.) then
!                    BNat = 0.0d0
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(5,(k2-1)*3+k3) = BNG(5,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(4,1),TyPt(4,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com1==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(1,(k2-1)*3+k3) = BNG(1,(k2-1)*3+k3) + dRAdxii(k2,1)*jacA(1,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                
!                !TP5 ---------------------------------------------------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(5,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(5,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                ubr = 0.0d0
!                do k2=2,ncpx+2*p
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)*Mb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(5,1),TyPt(5,2),TyPt(5,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + RA(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !shp = shp*shpANS(5)
!                sumw = sumw + shp
!                    
!                if(com5==.true.) then
!                    BNat = 0.0d0
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(5,(k2-1)*3+k3) = BNG(5,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(5,1),TyPt(5,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com1==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(1,(k2-1)*3+k3) = BNG(1,(k2-1)*3+k3) + dRAdxii(k2,1)*jacA(1,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                !TP6 ---------------------------------------------------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(6,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(6,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                ubr = 0.0d0
!                do k2=2,ncpx+2*p
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(1)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)*Mb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)*Mb(1)
!                
!                call ShapeFunc3(nel,TyPt(6,1),TyPt(6,2),TyPt(6,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + RA(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !shp = shp*shpANS(6)
!                sumw = sumw + shp
!                    
!                if(com5==.true.) then
!                    BNat = 0.0d0
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(5,(k2-1)*3+k3) = BNG(5,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(6,1),TyPt(6,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com1==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(1,(k2-1)*3+k3) = BNG(1,(k2-1)*3+k3) + dRAdxii(k2,1)*jacA(1,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                
!                !Normalize Strain-Displacement Matrix ------------------------------------------------------
!                if(com1==.true.) BNG(1,:) = BNG(1,:)/sumw
!                if(com5==.true.) BNG(5,:) = BNG(5,:)/sumw
!                sumw = 0.0d0
!                !-------------------------------------------------------------------------------------------
!                
!                
!                !TP7 ---------------------------------------------------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(7,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(7,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                vbr = 0.0d0
!                do k2=2,ncpy+2*q
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta) 
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp =       Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp =       Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp =       Mbr(1)
!                
!                call ShapeFunc3(nel,TyPt(7,1),TyPt(7,2),TyPt(7,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + RA(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !shp = shp*shpANS(7)
!                sumw = sumw + shp
!                
!                if(com6==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(6,(k2-1)*3+k3) = BNG(6,(k2-1)*3+k3) + (dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3))*shp
!                        end do
!                    end do
!                end if
!
!                !call ShapeFunc3(nel,TyPt(7,1),TyPt(7,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com2==.true.)then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(2,(k2-1)*3+k3) = BNG(2,(k2-1)*3+k3) + dRAdxii(k2,2)*jacA(2,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                !TP8 ---------------------------------------------------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(8,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(8,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                vbr = 0.0d0
!                do k2=2,ncpy+2*q
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp =       Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp =       Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp =       Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mbr(1)
!                
!                call ShapeFunc3(nel,TyPt(8,1),TyPt(8,2),TyPt(8,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + RA(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !shp = shp*shpANS(8)
!                sumw = sumw + shp 
!                    
!                if(com6==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(6,(k2-1)*3+k3) = BNG(6,(k2-1)*3+k3) + (dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3))*shp
!                        end do
!                    end do
!                end if
! 
!                !call ShapeFunc3(nel,TyPt(8,1),TyPt(8,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com2==.true.)then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(2,(k2-1)*3+k3) = BNG(2,(k2-1)*3+k3) + dRAdxii(k2,2)*jacA(2,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                !TP9 ---------------------------------------------------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(9,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(9,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                vbr = 0.0d0
!                do k2=2,ncpy+2*q
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp =       Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp =       Mbr(2)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp =       Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mbr(1)
!                
!                call ShapeFunc3(nel,TyPt(9,1),TyPt(9,2),TyPt(9,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + RA(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !shp = shp*shpANS(9)
!                sumw = sumw + shp
!                
!                if(com6==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(6,(k2-1)*3+k3) = BNG(6,(k2-1)*3+k3) + (dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3))*shp
!                        end do
!                    end do
!                end if
!
!                !call ShapeFunc3(nel,TyPt(9,1),TyPt(9,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com2==.true.)then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(2,(k2-1)*3+k3) = BNG(2,(k2-1)*3+k3) + dRAdxii(k2,2)*jacA(2,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                !TP10 --------------------------------------------------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(10,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(10,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                vbr = 0.0d0
!                do k2=2,ncpy+2*q
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)  
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp =       Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp =       Mbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp =       Mbr(1)
!                
!                call ShapeFunc3(nel,TyPt(10,1),TyPt(10,2),TyPt(10,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + RA(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !shp = shp*shpANS(10)
!                sumw = sumw + shp
!                
!                if(com6==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(6,(k2-1)*3+k3) = BNG(6,(k2-1)*3+k3) + (dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(10,1),TyPt(10,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com2==.true.)then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(2,(k2-1)*3+k3) = BNG(2,(k2-1)*3+k3) + dRAdxii(k2,2)*jacA(2,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                !TP11 --------------------------------------------------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(11,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(11,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                vbr = 0.0d0
!                do k2=2,ncpy+2*q
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)   
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp =       Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp =       Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp =       Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mbr(1)
!                
!                call ShapeFunc3(nel,TyPt(11,1),TyPt(11,2),TyPt(11,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + RA(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !shp = shp*shpANS(11)
!                sumw = sumw + shp
!                
!                if(com6==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(6,(k2-1)*3+k3) = BNG(6,(k2-1)*3+k3) + (dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(11,1),TyPt(11,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com2==.true.)then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(2,(k2-1)*3+k3) = BNG(2,(k2-1)*3+k3) + dRAdxii(k2,2)*jacA(2,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                !TP12 --------------------------------------------------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(12,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(12,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                vbr = 0.0d0
!                do k2=2,ncpy+2*q
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp =       Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp =       Mbr(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp =       Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mbr(1)
!                
!                call ShapeFunc3(nel,TyPt(12,1),TyPt(12,2),TyPt(12,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + RA(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !shp = shp*shpANS(12)
!                sumw = sumw + shp
!                
!                if(com6==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(6,(k2-1)*3+k3) = BNG(6,(k2-1)*3+k3) + (dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(12,1),TyPt(12,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com2==.true.)then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(2,(k2-1)*3+k3) = BNG(2,(k2-1)*3+k3) + dRAdxii(k2,2)*jacA(2,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                !Normalize Strain-Displacement Matrix ------------------------------------------------------
!                if(com2==.true.) BNG(2,:) = BNG(2,:)/sumw
!                if(com6==.true.) BNG(6,:) = BNG(6,:)/sumw
!                sumw = 0.0d0
!                !-------------------------------------------------------------------------------------------
!
!                !TP13 --------------------------------------------------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(13,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(13,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!                
!                ubr = 0.0d0
!                do k2=2,ncpx+2*p
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                vbr = 0.0d0
!                do k2=2,ncpy+2*q
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)*Mbr(1)
!                
!                call ShapeFunc3(nel,TyPt(13,1),TyPt(13,2),TyPt(13,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + RA(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !shp = shp*shpANS(13)
!                sumw = sumw + shp
!                
!                if(com4==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(4,(k2-1)*3+k3) = BNG(4,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                !-------------------------------------------------------------------------------------------
!                
!                !TP14 --------------------------------------------------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(14,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(14,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                ubr = 0.0d0
!                do k2=2,ncpx+2*p
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                vbr = 0.0d0
!                do k2=2,ncpy+2*q
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)*Mbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)*Mbr(1)
!                
!                call ShapeFunc3(nel,TyPt(14,1),TyPt(14,2),TyPt(14,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                 
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + RA(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !shp = shp*shpANS(14)
!                sumw = sumw + shp
!                    
!                if(com4==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(4,(k2-1)*3+k3) = BNG(4,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                !-------------------------------------------------------------------------------------------
!                
!                !TP15 --------------------------------------------------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(15,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(15,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                ubr = 0.0d0
!                do k2=2,ncpx+2*p
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                vbr = 0.0d0
!                do k2=2,ncpy+2*q
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(1)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(1)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(1)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)*Mbr(1)
!                
!                call ShapeFunc3(nel,TyPt(15,1),TyPt(15,2),TyPt(15,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + RA(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !shp = shp*shpANS(15)
!                sumw = sumw + shp
!                    
!                if(com4==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(4,(k2-1)*3+k3) = BNG(4,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                !-------------------------------------------------------------------------------------------
!                
!                !TP16 --------------------------------------------------------------------------------------
!                uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(16,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(16,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                call kvecANS(ncpx,p,u_knot,uk,ub)
!                call kvecANS(ncpy,q,v_knot,vk,vb)
!
!                ubr = 0.0d0
!                do k2=2,ncpx+2*p
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                vbr = 0.0d0
!                do k2=2,ncpy+2*q
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(1)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(1)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)*Mbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(1)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)*Mbr(1)
!                
!                call ShapeFunc3(nel,TyPt(16,1),TyPt(16,2),TyPt(16,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                wei = 0.0d0
!                do k2=1,(p+1)*(q+1)*(w+1)
!                    wei = wei + RA(k2)*weight(IEN(nel,k2))
!                end do
!                
!                !shp = shp*shpANS(16)
!                sumw = sumw + shp
!                    
!                if(com4==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(4,(k2-1)*3+k3) = BNG(4,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                !-------------------------------------------------------------------------------------------
!                
!                !Normalize Strain-Displacement Matrix ------------------------------------------------------
!                if(com4==.true.) BNG(4,:) = BNG(4,:)/sumw
!                !-------------------------------------------------------------------------------------------
!
!                ! Strain-Displacement Matrix in local frame
!                Bloc=matmul(TCL,BNG)
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
!    continue
!
!end subroutine ElHex27ANS


!!------------------------------------------------------------------------------------------------------------
!!
!! Second order hexahedral element with ANS
!!
!! Input: nel - Element number
!!        el_ddisp - Element displacement
!!        inc - Increment number
!!
!! Output: Ke - Elementar stiffness matrix
!!         Finte - Elementar internal forces
!!------------------------------------------------------------------------------------------------------------
!
!! H2ANSE1N
!
!subroutine ElHex27ANS(nel,Ke,el_ddisp,Finte,inc)
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
!    real(8),dimension(2*(p+1)+p)::ub
!    real(8),dimension(2*(q+1)+q)::vb
!    real(8),dimension(2*(w+1)+w)::wb
!    
!    real(8),dimension(2*(p+1))::ubr
!    real(8),dimension(2*(q+1))::vbr
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
!    real(8)::shp,suma,uk,vk,wk,wei,sumw
!    
!    real(8),dimension(3,3)::jact
!    
!    real(8),dimension(6,6)::NTP
!    
!    logical::lsp,endcycle
!    logical::com1,com2,com3,com4,com5,com6
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
!!    call ShapeFunc3(nel,0.0d0,0.0d0,0.0d0,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,tzero)
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
!                xi = e(i)
!                eta = n(j)
!                zeta = c(k)
!                
!                call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tzero)
!                
!                !Weight factor
!                vol = vol + we(i)*wn(j)*wc(k)*detj
!            end do
!        end do
!    end do
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
!!                Bglob = 0.0d0
!!                do k1=1,(p+1)*(q+1)*(w+1)
!!                    
!!                    Bglob(1,k1*3-2) = dRdx(k1,1)
!!                    Bglob(2,k1*3-1) = dRdx(k1,2)
!!                    Bglob(3,k1*3  ) = dRdx(k1,3)
!!                    
!!                    Bglob(4,k1*3-2) = dRdx(k1,2)
!!                    Bglob(4,k1*3-1) = dRdx(k1,1)
!!                    Bglob(5,k1*3-2) = dRdx(k1,3)
!!                    Bglob(5,k1*3  ) = dRdx(k1,1)
!!                    Bglob(6,k1*3-1) = dRdx(k1,3)
!!                    Bglob(6,k1*3  ) = dRdx(k1,2)
!!                    
!!                end do
!                
!                BConv = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
!                    do k2 =1,nds
!                        BConv(1,(k1-1)*3+k2) = dRdxii(k1,1)*jac(1,k2)
!                        BConv(2,(k1-1)*3+k2) = dRdxii(k1,2)*jac(2,k2)
!                        BConv(3,(k1-1)*3+k2) = dRdxii(k1,3)*jac(3,k2)
!                        
!                        BConv(4,(k1-1)*3+k2) = dRdxii(k1,1)*jac(2,k2) + dRdxii(k1,2)*jac(1,k2)
!                        BConv(5,(k1-1)*3+k2) = dRdxii(k1,1)*jac(3,k2) + dRdxii(k1,3)*jac(1,k2)
!                        BConv(6,(k1-1)*3+k2) = dRdxii(k1,2)*jac(3,k2) + dRdxii(k1,3)*jac(2,k2)
!                    end do
!                end do
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
!                if(com1 == .true.) BNG(1,:) = 0.0d0
!                if(com2 == .true.) BNG(2,:) = 0.0d0
!                !if(com3 == .true.) BNG(3,:) = 0.0d0
!                if(com4 == .true.) BNG(4,:) = 0.0d0
!                if(com5 == .true.) BNG(5,:) = 0.0d0
!                if(com6 == .true.) BNG(6,:) = 0.0d0
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
!                
!                !TP1 ---------------------------------------------------------------------------------------
!                ub = u_knot(ni+1)
!                ubr = u_knot(ni+1)
!                vb = v_knot(nj+1)
!                vbr = v_knot(nj+1)
!                
!                ub(1) = 0.0d0; vb(1) = 0.0d0
!                ub(2) = 0.0d0; vb(2) = 0.0d0
!                ub(3) = 0.0d0; vb(3) = 0.0d0
!                
!                ubr(1) = 0.0d0
!                ubr(2) = 0.0d0
!                vbr(1) = 0.0d0
!                vbr(2) = 0.0d0
!                
!                ub(4) = ((u_knot(ni+1) - u_knot(ni))*TyPt(1,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5) = ((u_knot(ni+1) - u_knot(ni))*TyPt(1,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4) = ((v_knot(nj+1) - v_knot(nj))*TyPt(1,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5) = ((v_knot(nj+1) - v_knot(nj))*TyPt(1,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                ubr(3) = ((u_knot(ni+1) - u_knot(ni))*TyPt(1,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                
!                ubr = 0.0d0
!                do k2=2,2*(p+1)+1
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(2)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)
!                
!                sumw = sumw + shp
!                
!                call ShapeFunc3(nel,TyPt(1,1),TyPt(1,2),TyPt(1,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                if(com5==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(5,(k2-1)*3+k3) = BNG(5,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(1,1),TyPt(1,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com1==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(1,(k2-1)*3+k3) = BNG(1,(k2-1)*3+k3) + dRAdxii(k2,1)*jacA(1,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                
!                !TP2 ---------------------------------------------------------------------------------------
!                ub(4)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(2,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(2,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(2,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(2,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                ubr(3) = ((u_knot(ni+1) - u_knot(ni))*TyPt(2,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                
!                ubr = 0.0d0
!                do k2=2,2*(p+1)+1
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi) 
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)
!                
!                sumw = sumw + shp
!                
!                call ShapeFunc3(nel,TyPt(2,1),TyPt(2,2),TyPt(2,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                    
!                if(com5==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(5,(k2-1)*3+k3) = BNG(5,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(2,1),TyPt(2,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com1==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(1,(k2-1)*3+k3) = BNG(1,(k2-1)*3+k3) + dRAdxii(k2,1)*jacA(1,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                
!                !TP3 ---------------------------------------------------------------------------------------
!                ub(4)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(3,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(3,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(3,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(3,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                ubr(3) = ((u_knot(ni+1) - u_knot(ni))*TyPt(3,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                
!                ubr = 0.0d0
!                do k2=2,2*(p+1)+1
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi) 
!                
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(2)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)*Mb(1)
!                
!                sumw = sumw + shp
!                
!                call ShapeFunc3(nel,TyPt(3,1),TyPt(3,2),TyPt(3,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                if(com5==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(5,(k2-1)*3+k3) = BNG(5,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(3,1),TyPt(3,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com1==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(1,(k2-1)*3+k3) = BNG(1,(k2-1)*3+k3) + dRAdxii(k2,1)*jacA(1,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                
!                !TP4 ---------------------------------------------------------------------------------------
!                ub(4)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(4,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(4,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(4,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(4,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                ubr(3) = ((u_knot(ni+1) - u_knot(ni))*TyPt(4,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                
!                ubr = 0.0d0
!                do k2=2,2*(p+1)+1
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)*Mb(3)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(1)*Mb(3)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)*Mb(3)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)*Mb(1)
!                
!                sumw = sumw + shp
!                
!                call ShapeFunc3(nel,TyPt(4,1),TyPt(4,2),TyPt(4,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                if(com5==.true.) then
!                    BNat = 0.0d0
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(5,(k2-1)*3+k3) = BNG(5,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(4,1),TyPt(4,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com1==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(1,(k2-1)*3+k3) = BNG(1,(k2-1)*3+k3) + dRAdxii(k2,1)*jacA(1,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                
!                !TP5 ---------------------------------------------------------------------------------------
!                ub(4)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(5,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(5,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(5,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(5,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                ubr(3) = ((u_knot(ni+1) - u_knot(ni))*TyPt(5,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                
!                ubr = 0.0d0
!                do k2=2,2*(p+1)+1
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)*Mb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)*Mb(1)
!                
!                sumw = sumw + shp
!                
!                call ShapeFunc3(nel,TyPt(5,1),TyPt(5,2),TyPt(5,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                    
!                if(com5==.true.) then
!                    BNat = 0.0d0
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(5,(k2-1)*3+k3) = BNG(5,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(5,1),TyPt(5,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com1==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(1,(k2-1)*3+k3) = BNG(1,(k2-1)*3+k3) + dRAdxii(k2,1)*jacA(1,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                !TP6 ---------------------------------------------------------------------------------------
!                ub(4)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(6,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(6,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(6,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(6,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                ubr(3) = ((u_knot(ni+1) - u_knot(ni))*TyPt(6,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                
!                ubr = 0.0d0
!                do k2=2,2*(p+1)+1
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(1)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)*Mb(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)*Mb(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(1)*Mb(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)*Mb(1)
!                
!                sumw = sumw + shp
!                
!                call ShapeFunc3(nel,TyPt(6,1),TyPt(6,2),TyPt(6,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                    
!                if(com5==.true.) then
!                    BNat = 0.0d0
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(5,(k2-1)*3+k3) = BNG(5,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(6,1),TyPt(6,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com1==.true.) then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(1,(k2-1)*3+k3) = BNG(1,(k2-1)*3+k3) + dRAdxii(k2,1)*jacA(1,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                
!                !Normalize Strain-Displacement Matrix ------------------------------------------------------
!                if(com1==.true.) BNG(1,:) = BNG(1,:)/sumw
!                if(com5==.true.) BNG(5,:) = BNG(5,:)/sumw
!                sumw = 0.0d0
!                !-------------------------------------------------------------------------------------------
!                
!                
!                !TP7 ---------------------------------------------------------------------------------------
!                ub(4)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(7,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(7,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(7,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(7,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vbr(3) = ((v_knot(nj+1) - v_knot(nj))*TyPt(7,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                vbr = 0.0d0
!                do k2=2,2*(q+1)+1
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta) 
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp =       Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp =       Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp =       Mbr(1)
!                
!                call ShapeFunc3(nel,TyPt(7,1),TyPt(7,2),TyPt(7,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                sumw = sumw+shp
!                
!                if(com6==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(6,(k2-1)*3+k3) = BNG(6,(k2-1)*3+k3) + (dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3))*shp
!                        end do
!                    end do
!                end if
!
!                !call ShapeFunc3(nel,TyPt(7,1),TyPt(7,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com2==.true.)then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(2,(k2-1)*3+k3) = BNG(2,(k2-1)*3+k3) + dRAdxii(k2,2)*jacA(2,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                !TP8 ---------------------------------------------------------------------------------------
!                ub(4)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(8,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(8,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(8,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(8,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vbr(3) = ((v_knot(nj+1) - v_knot(nj))*TyPt(8,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                vbr = 0.0d0
!                do k2=2,2*(q+1)+1
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp =       Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp =       Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp =       Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mbr(1)
!                
!                sumw = sumw+shp
!                
!                call ShapeFunc3(nel,TyPt(8,1),TyPt(8,2),TyPt(8,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                    
!                if(com6==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(6,(k2-1)*3+k3) = BNG(6,(k2-1)*3+k3) + (dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3))*shp
!                        end do
!                    end do
!                end if
! 
!                !call ShapeFunc3(nel,TyPt(8,1),TyPt(8,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com2==.true.)then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(2,(k2-1)*3+k3) = BNG(2,(k2-1)*3+k3) + dRAdxii(k2,2)*jacA(2,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                !TP9 ---------------------------------------------------------------------------------------
!                ub(4)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(9,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(9,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(9,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(9,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vbr(3) = ((v_knot(nj+1) - v_knot(nj))*TyPt(9,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                vbr = 0.0d0
!                do k2=2,2*(q+1)+1
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp =       Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp =       Mbr(2)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp =       Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mbr(1)
!                
!                sumw = sumw+shp
!                
!                call ShapeFunc3(nel,TyPt(9,1),TyPt(9,2),TyPt(9,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                if(com6==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(6,(k2-1)*3+k3) = BNG(6,(k2-1)*3+k3) + (dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3))*shp
!                        end do
!                    end do
!                end if
!
!                !call ShapeFunc3(nel,TyPt(9,1),TyPt(9,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com2==.true.)then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(2,(k2-1)*3+k3) = BNG(2,(k2-1)*3+k3) + dRAdxii(k2,2)*jacA(2,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                !TP10 --------------------------------------------------------------------------------------
!                ub(4)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(10,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(10,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(10,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(10,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vbr(3) = ((v_knot(nj+1) - v_knot(nj))*TyPt(10,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                vbr = 0.0d0
!                do k2=2,2*(q+1)+1
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)  
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp =       Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp =       Mbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp =       Mbr(1)
!                
!                sumw = sumw+shp
!                
!                call ShapeFunc3(nel,TyPt(10,1),TyPt(10,2),TyPt(10,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                if(com6==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(6,(k2-1)*3+k3) = BNG(6,(k2-1)*3+k3) + (dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(10,1),TyPt(10,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com2==.true.)then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(2,(k2-1)*3+k3) = BNG(2,(k2-1)*3+k3) + dRAdxii(k2,2)*jacA(2,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                !TP11 --------------------------------------------------------------------------------------
!                ub(4)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(11,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(11,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(11,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(11,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vbr(3) = ((v_knot(nj+1) - v_knot(nj))*TyPt(11,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                vbr = 0.0d0
!                do k2=2,2*(q+1)+1
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)   
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(3)*Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp =       Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp =       Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(3)*Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp =       Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mbr(1)
!                
!                sumw = sumw+shp
!                
!                call ShapeFunc3(nel,TyPt(11,1),TyPt(11,2),TyPt(11,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                if(com6==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(6,(k2-1)*3+k3) = BNG(6,(k2-1)*3+k3) + (dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(11,1),TyPt(11,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com2==.true.)then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(2,(k2-1)*3+k3) = BNG(2,(k2-1)*3+k3) + dRAdxii(k2,2)*jacA(2,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                !TP12 --------------------------------------------------------------------------------------
!                ub(4)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(12,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(12,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(12,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(12,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vbr(3) = ((v_knot(nj+1) - v_knot(nj))*TyPt(12,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                
!                vbr = 0.0d0
!                do k2=2,2*(q+1)+1
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)   
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp =       Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nb(1)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nb(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp =       Mbr(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nb(1)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nb(1)*Mbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp =       Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nb(1)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nb(1)*Mbr(1)
!                
!                sumw = sumw+shp
!                
!                call ShapeFunc3(nel,TyPt(12,1),TyPt(12,2),TyPt(12,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!               
!                if(com6==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(6,(k2-1)*3+k3) = BNG(6,(k2-1)*3+k3) + (dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3))*shp
!                        end do
!                    end do
!                end if
!                
!                !call ShapeFunc3(nel,TyPt(12,1),TyPt(12,2),zeta,RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                !jacA = transpose(jacA)
!                
!                if(com2==.true.)then
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(2,(k2-1)*3+k3) = BNG(2,(k2-1)*3+k3) + dRAdxii(k2,2)*jacA(2,k3)*shp
!                        end do
!                    end do
!                endif
!                !-------------------------------------------------------------------------------------------
!                
!                !Normalize Strain-Displacement Matrix ------------------------------------------------------
!                if(com2==.true.) BNG(2,:) = BNG(2,:)/sumw
!                if(com6==.true.) BNG(6,:) = BNG(6,:)/sumw
!                sumw = 0.0d0
!                !-------------------------------------------------------------------------------------------
!
!                !TP13 --------------------------------------------------------------------------------------
!                ub(4)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(13,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(13,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(13,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(13,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                ubr(3) = ((u_knot(ni+1) - u_knot(ni))*TyPt(13,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vbr(3) = ((v_knot(nj+1) - v_knot(nj))*TyPt(13,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!
!                ubr = 0.0d0
!                do k2=2,2*(p+1)+1
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                vbr = 0.0d0
!                do k2=2,2*(q+1)+1
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)*Mbr(1)
!                
!                sumw = sumw + shp
!                
!                call ShapeFunc3(nel,TyPt(13,1),TyPt(13,2),TyPt(13,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                
!                if(com4==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(4,(k2-1)*3+k3) = BNG(4,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                !-------------------------------------------------------------------------------------------
!                
!                !TP14 --------------------------------------------------------------------------------------
!                ub(4)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(14,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(14,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(14,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(14,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                ubr(3) = ((u_knot(ni+1) - u_knot(ni))*TyPt(14,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vbr(3) = ((v_knot(nj+1) - v_knot(nj))*TyPt(14,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!
!                ubr = 0.0d0
!                do k2=2,2*(p+1)+1
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                vbr = 0.0d0
!                do k2=2,2*(q+1)+1
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)*Mbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)*Mbr(1)
!                
!                sumw = sumw + shp
!                
!                call ShapeFunc3(nel,TyPt(14,1),TyPt(14,2),TyPt(14,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                    
!                if(com4==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(4,(k2-1)*3+k3) = BNG(4,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                !-------------------------------------------------------------------------------------------
!                
!                !TP15 --------------------------------------------------------------------------------------
!                ub(4)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(15,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(15,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(15,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(15,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                ubr(3) = ((u_knot(ni+1) - u_knot(ni))*TyPt(15,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vbr(3) = ((v_knot(nj+1) - v_knot(nj))*TyPt(15,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!
!                ubr = 0.0d0
!                do k2=2,2*(p+1)+1
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                vbr = 0.0d0
!                do k2=2,2*(q+1)+1
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(1)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(1)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(1)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)*Mbr(1)
!                
!                sumw = sumw + shp
!                
!                call ShapeFunc3(nel,TyPt(15,1),TyPt(15,2),TyPt(15,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                    
!                if(com4==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(4,(k2-1)*3+k3) = BNG(4,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                !-------------------------------------------------------------------------------------------
!                
!                !TP16 --------------------------------------------------------------------------------------
!                ub(4)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(16,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                ub(5)  = ((u_knot(ni+1) - u_knot(ni))*TyPt(16,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vb(4)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(16,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                vb(5)  = ((v_knot(nj+1) - v_knot(nj))*TyPt(16,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                ubr(3) = ((u_knot(ni+1) - u_knot(ni))*TyPt(16,1)   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                vbr(3) = ((v_knot(nj+1) - v_knot(nj))*TyPt(16,2)   + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!
!                ubr = 0.0d0
!                do k2=2,2*(p+1)+1
!                    ubr(k2-1) = ub(k2)
!                end do
!                
!                vbr = 0.0d0
!                do k2=2,2*(q+1)+1
!                    vbr(k2-1) = vb(k2)
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                
!                if(xi .lt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(2)*Mbr(2)
!                if(xi == 0.0d0   .and. eta .lt. 0.0d0) shp = Nbr(1)*Mbr(2)
!                if(xi .gt. 0.0d0 .and. eta .lt. 0.0d0) shp = Nbr(1)*Mbr(2)
!                
!                if(xi .lt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi == 0.0d0   .and. eta == 0.0d0) shp = Nbr(1)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta == 0.0d0) shp = Nbr(1)*Mbr(1)
!                
!                if(xi .lt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(2)*Mbr(1)
!                if(xi == 0.0d0   .and. eta .gt. 0.0d0) shp = Nbr(1)*Mbr(1)
!                if(xi .gt. 0.0d0 .and. eta .gt. 0.0d0) shp = Nbr(1)*Mbr(1)
!                
!                sumw = sumw + shp
!                
!                call ShapeFunc3(nel,TyPt(16,1),TyPt(16,2),TyPt(16,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                jacA = transpose(jacA)
!                    
!                if(com4==.true.) then    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(4,(k2-1)*3+k3) = BNG(4,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3))*shp
!                        end do
!                    end do
!                end if
!                !-------------------------------------------------------------------------------------------
!                
!                !Normalize Strain-Displacement Matrix ------------------------------------------------------
!                if(com4==.true.) BNG(4,:) = BNG(4,:)/sumw
!                !-------------------------------------------------------------------------------------------
!
!                ! Strain-Displacement Matrix in local frame
!                Bloc=matmul(TCL,BNG)
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
!    continue
!
!end subroutine ElHex27ANS


!------------------------------------------------------------------------------------------------------------
!
! Second order hexahedral element with ANS
!
! Input: nel - Element number
!        el_ddisp - Element displacement
!        inc - Increment number
!
! Output: Ke - Elementar stiffness matrix
!         Finte - Elementar internal forces
!------------------------------------------------------------------------------------------------------------

!H2ANSProj

subroutine ElHex27ANS(nel,Ke,el_ddisp,Finte,inc)
    
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
    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BConv
    
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
    
!    call ShapeFunc3(nel,0.0d0,0.0d0,0.0d0,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,tzero)
    
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
    KNLG = 0.0d0
    BNL = 0.0d0
    MSTR = 0.0d0
    
    vol = 0.0d0
    
    Finte = 0.0d0
    do k=1,npi_zeta
        do j=1,npi_eta
            do i=1,npi_xi
                
                xi = e(i)
                eta = n(j)
                zeta = c(k)
                
                call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tzero)
                
                !Weight factor
                vol = vol + we(i)*wn(j)*wc(k)*detj
            end do
        end do
    end do
    
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
                
                BConv = 0.0d0
                do k1=1,(p+1)*(q+1)*(w+1)
                    do k2 =1,nds
                        BConv(1,(k1-1)*3+k2) = dRdxii(k1,1)*jac(1,k2)
                        BConv(2,(k1-1)*3+k2) = dRdxii(k1,2)*jac(2,k2)
                        BConv(3,(k1-1)*3+k2) = dRdxii(k1,3)*jac(3,k2)
                        
                        BConv(4,(k1-1)*3+k2) = dRdxii(k1,1)*jac(2,k2) + dRdxii(k1,2)*jac(1,k2)
                        BConv(5,(k1-1)*3+k2) = dRdxii(k1,1)*jac(3,k2) + dRdxii(k1,3)*jac(1,k2)
                        BConv(6,(k1-1)*3+k2) = dRdxii(k1,2)*jac(3,k2) + dRdxii(k1,3)*jac(2,k2)
                    end do
                end do
                
                !Local strain-displacement matrix at end configuration
                Bloc=0.0d0
                
                !BConv=matmul(TGC,Bglob)
                !BNG = matmul(TGC,Bglob)
                
                BNG = BConv
                
                !BNGold = BNG
                
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
                
!                TyPt(17,1) =     CA; TyPt(17,2) =  1.0d0; TyPt(17,3) =  zeta
!                TyPt(18,1) =    -CA; TyPt(18,2) =  1.0d0; TyPt(18,3) =  zeta
!                TyPt(19,1) =     CA; TyPt(19,2) = -1.0d0; TyPt(19,3) =  zeta
!                TyPt(20,1) =    -CA; TyPt(20,2) = -1.0d0; TyPt(20,3) =  zeta
!                TyPt(21,1) =  1.0d0; TyPt(21,2) =     CA; TyPt(21,3) =  zeta
!                TyPt(22,1) =  1.0d0; TyPt(22,2) =    -CA; TyPt(22,3) =  zeta
!                TyPt(23,1) = -1.0d0; TyPt(23,2) =     CA; TyPt(23,3) =  zeta
!                TyPt(24,1) = -1.0d0; TyPt(24,2) =    -CA; TyPt(24,3) =  zeta
                
                com1 = .true.
                com2 = .true.
                !com3 = .false.
                com4 = .true.
                com5 = .true.
                com6 = .true.
                
                if(com1 == .true.) BNG(1,:) = 0.0d0
                if(com2 == .true.) BNG(2,:) = 0.0d0
                !if(com3 == .true.) BNG(3,:) = 0.0d0
                if(com4 == .true.) BNG(4,:) = 0.0d0
                if(com5 == .true.) BNG(5,:) = 0.0d0
                if(com6 == .true.) BNG(6,:) = 0.0d0
                
                ni = INN(IEN(nel,1),1)
                nj = INN(IEN(nel,1),2)
                nk = INN(IEN(nel,1),3)
    
                xib   = ((u_knot(ni+1) - u_knot(ni))*xi   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
                etab  = ((v_knot(nj+1) - v_knot(nj))*eta  + (v_knot(nj+1) + v_knot(nj)))/2.0d0
                zetab = ((w_knot(nk+1) - w_knot(nk))*zeta + (w_knot(nk+1) + w_knot(nk)))/2.0d0
                
                sumw = 0.0d0
                
                wei = 0.0d0
                do k2=1,(p+1)*(q+1)*(w+1)
                    wei = wei + R(k2)*weight(IEN(nel,k2))
                end do
                
!                !-------------------------------------------------------------------------------------------
!                
!                A = CA
!                B = CB
!                
!                shpANS = 0.0d0
!                !-------------------------------------------------------------------
!                shpANS(1)  = eta/(4.0d0*B)*(eta/B+1.0d0)*(1.0d0+xi/A) 
!                shpANS(2)  = eta/(4.0d0*B)*(eta/B+1.0d0)*(1.0d0-xi/A)
!                shpANS(3)  = 1.0d0/2.0d0*(1.0d0-eta*eta/(B*B))*(1.0d0+xi/A)
!                shpANS(4)  = 1.0d0/2.0d0*(1.0d0-eta*eta/(B*B))*(1.0d0-xi/A)
!                shpANS(5)  = eta/(4.0d0*B)*(eta/B-1.0d0)*(1.0d0+xi/A) 
!                shpANS(6)  = eta/(4.0d0*B)*(eta/B-1.0d0)*(1.0d0-xi/A)
!                !-------------------------------------------------------------------
!                
!                !-------------------------------------------------------------------
!                shpANS(7)  = xi/(4.0d0*B)*(xi/B+1.0d0)*(1.0d0+eta/A)
!                shpANS(8)  = 1.0d0/2.0d0*(1.0d0-xi*xi/(B*B))*(1.0d0+eta/A)
!                shpANS(9)  = xi/(4.0d0*B)*(xi/B-1.0d0)*(1.0d0+eta/A)
!                shpANS(10) = xi/(4.0d0*B)*(xi/B+1.0d0)*(1.0d0-eta/A)
!                shpANS(11) = 1.0d0/2.0d0*(1.0d0-xi*xi/(B*B))*(1.0d0-eta/A)
!                shpANS(12) = xi/(4.0d0*B)*(xi/B-1.0d0)*(1.0d0-eta/A)
!                !-------------------------------------------------------------------
!                
!                !-------------------------------------------------------------------
!                shpANS(13) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0+eta/A)
!                shpANS(14) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0+eta/A)
!                shpANS(15) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0-eta/A)
!                shpANS(16) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0-eta/A)
!                !-------------------------------------------------------------------
                
                
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
                
                
                !First TP set cycle ------------------------------------------------
                !Mass Matrix -------------------------------------------------------
                MTP = 0.0d0
                do k1=1,6
                    uk = (TyPt(k1,1) + 1.0d0)/2.0d0
                    vk = (TyPt(k1,2) + 1.0d0)/2.0d0
                    
                    Nbr = 0.0d0
                    call BSplineBasisAndDeriv(ncpx+p,p-1,uk,ubr,Nbr,dNbrdxi)
                    Mb = 0.0d0
                    call BSplineBasisAndDeriv(ncpy+q,q  ,vk, vb, Mb,dMbdeta)
                    
                    count = 0
                    do k2=0,q
                        do k3=0,p-1
                            count = count + 1
                            MTP(k1,count) = Nbr(p-k3)*Mb(q+1-k2)
                        end do
                    end do      
                end do
                
                !IP Mass Matrix ----------------------------------------------------
                xib  = ( xi + 1.0d0)/2.0d0
                etab = (eta + 1.0d0)/2.0d0
                
                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
                call BSplineBasisAndDeriv(ncpy+q,q  ,etab, vb, Mb,dMbdeta)
                    
                count = 0
                do k2=0,q
                    do k3=0,p-1
                        count = count + 1
                        M(1,count) = Nbr(p-k3)*Mb(q+1-k2)
                    end do
                end do 
                
                temp6 = 0.0d0
                call gaussj(MTP,6,temp6,ilixo)
                
                continue
                
                MMult = matmul(M,MTP)
                
                do k1=1,6
                    call ShapeFunc3(nel,TyPt(k1,1),TyPt(k1,2),TyPt(k1,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
                    jacA = transpose(jacA)
                    
                    do k2=1,(p+1)*(q+1)*(w+1)
                        do k3 =1,nds
                            BNG(1,(k2-1)*3+k3) = BNG(1,(k2-1)*3+k3) + dRAdxii(k2,1)*jacA(1,k3)*MMult(1,k1)
                            BNG(5,(k2-1)*3+k3) = BNG(5,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3))*MMult(1,k1)
                        end do
                    end do
                    
                    continue
                    
                end do
                
                !Second TP set cycle -----------------------------------------------
                !Mass Matrix -------------------------------------------------------
                MTP = 0.0d0
                do k1=1,6
                    uk = (TyPt(k1+6,1) + 1.0d0)/2.0d0
                    vk = (TyPt(k1+6,2) + 1.0d0)/2.0d0
                    
                    Nb = 0.0d0
                    call BSplineBasisAndDeriv(ncpx+p,p  ,uk, ub,  Nb,  dNbdxi)
                    Mbr = 0.0d0
                    call BSplineBasisAndDeriv(ncpy+q,q-1,vk,vbr, Mbr,dMbrdeta)
                    
                    count = 0
                    do k2=0,q-1
                        do k3=0,p
                            count = count + 1
                            MTP(k1,count) = Nb(p+1-k3)*Mbr(q-k2)
                        end do
                    end do      
                end do
                
                temp6 = 0.0d0
                call gaussj(MTP,6,temp6,ilixo)
                continue
                
                !IP Mass Matrix ----------------------------------------------------
                xib  = ( xi + 1.0d0)/2.0d0
                etab = (eta + 1.0d0)/2.0d0
                
                call BSplineBasisAndDeriv(ncpx+p,p  , xib, ub, Nb,  dNbdxi)
                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
                    
                count = 0
                M = 0.0d0
                do k2=0,q-1
                    do k3=0,p
                        count = count + 1
                        M(1,count) = Nb(p+1-k3)*Mbr(q-k2)
                    end do
                end do 
                
                MMult = 0.0d0
                MMult = matmul(M,MTP)
                continue
                
                do k1=1,6
                    call ShapeFunc3(nel,TyPt(k1+6,1),TyPt(k1+6,2),TyPt(k1+6,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
                    jacA = transpose(jacA)
                    
                    do k2=1,(p+1)*(q+1)*(w+1)
                        do k3 =1,nds
                            BNG(2,(k2-1)*3+k3) = BNG(2,(k2-1)*3+k3) + dRAdxii(k2,2)*jacA(2,k3)*MMult(1,k1)
                            BNG(6,(k2-1)*3+k3) = BNG(6,(k2-1)*3+k3) + (dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3))*MMult(1,k1)
                        end do
                    end do
                end do
                
                
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
                    
                    continue
                         
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
                    call ShapeFunc3(nel,TyPt(k1+12,1),TyPt(k1+12,2),TyPt(k1+12,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
                    jacA = transpose(jacA)
                    
                    do k2=1,(p+1)*(q+1)*(w+1)
                        do k3 =1,nds
                            !BNG(1,(k2-1)*3+k3) = BNG(1,(k2-1)*3+k3) + dRAdxii(k2,1)*jacA(1,k3)*MMultr(1,k1)
                            !BNG(2,(k2-1)*3+k3) = BNG(2,(k2-1)*3+k3) + dRAdxii(k2,2)*jacA(2,k3)*MMultr(1,k1)
                            BNG(4,(k2-1)*3+k3) = BNG(4,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3))*MMultr(1,k1)
                            !BNG(5,(k2-1)*3+k3) = BNG(5,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3))*MMultr(1,k1)
                            !BNG(6,(k2-1)*3+k3) = BNG(6,(k2-1)*3+k3) + (dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3))*MMultr(1,k1)
                        end do
                    end do
                end do
                !-------------------------------------------------------------------------------------------

                ! Strain-Displacement Matrix in local frame
                Bloc=matmul(TCL,BNG)

                if(nlgeom == .false.) Bloc_mid = Bloc
                
                call MatPlastic3D(iprops,nds,(p+1)*(q+1)*(w+1),props,Bloc_mid,6,el_ddisp,stress(nel,cpi,:),strain(nel,cpi,:),dstrain(nel,cpi,:),dstress(nel,cpi,:),hard(nel,cpi),matD)
                
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
                
                Finte = Finte + matmul(transpose(Bloc),tempstr)*gwt
                
                
                !Gravitic loads -------------------------------------------------------------
                if (gravity==.true.)then
                    call grvt(nshpl,nds,R,gwt,gconst,gravdir,density,Finte)
                end if
                
                !Pressure loads -------------------------------------------------------------
                if(pressure==.true.) then
                    do k1=1,nelems
                        if(k1==nel .and. abs(pressvec(k1)) .gt. 0.0d0) then
                            call prssr(nel,npi_xi,npi_eta,npi_zeta,jac,xi,eta,zeta,we(i),wn(j),wc(k),Finte)
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

end subroutine ElHex27ANS

!!------------------------------------------------------------------------------------------------------------
!!
!! Second order hexahedral element with ANS
!!
!! Input: nel - Element number
!!        el_ddisp - Element displacement
!!        inc - Increment number
!!
!! Output: Ke - Elementar stiffness matrix
!!         Finte - Elementar internal forces
!!------------------------------------------------------------------------------------------------------------
!
!!H2ANSProjGlob
!
!subroutine ElHex27ANS(nel,Ke,el_ddisp,Finte,inc)
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
!    real(8),dimension(ncpx+p+1)::ub
!    real(8),dimension(ncpy+q+1)::vb
!    
!    real(8),dimension(ncpx+p-1)::ubr
!    real(8),dimension(ncpy+q-1)::vbr
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
!!    call ShapeFunc3(nel,0.0d0,0.0d0,0.0d0,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,tzero)
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
!                xi = e(i)
!                eta = n(j)
!                zeta = c(k)
!                
!                call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tzero)
!                
!                !Weight factor
!                vol = vol + we(i)*wn(j)*wc(k)*detj
!            end do
!        end do
!    end do
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
!!                Bglob = 0.0d0
!!                do k1=1,(p+1)*(q+1)*(w+1)
!!                    
!!                    Bglob(1,k1*3-2) = dRdx(k1,1)
!!                    Bglob(2,k1*3-1) = dRdx(k1,2)
!!                    Bglob(3,k1*3  ) = dRdx(k1,3)
!!                    
!!                    Bglob(4,k1*3-2) = dRdx(k1,2)
!!                    Bglob(4,k1*3-1) = dRdx(k1,1)
!!                    Bglob(5,k1*3-2) = dRdx(k1,3)
!!                    Bglob(5,k1*3  ) = dRdx(k1,1)
!!                    Bglob(6,k1*3-1) = dRdx(k1,3)
!!                    Bglob(6,k1*3  ) = dRdx(k1,2)
!!                    
!!                end do
!                
!                BConv = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
!                    do k2 =1,nds
!                        BConv(1,(k1-1)*3+k2) = dRdxii(k1,1)*jac(1,k2)
!                        BConv(2,(k1-1)*3+k2) = dRdxii(k1,2)*jac(2,k2)
!                        BConv(3,(k1-1)*3+k2) = dRdxii(k1,3)*jac(3,k2)
!                        
!                        BConv(4,(k1-1)*3+k2) = dRdxii(k1,1)*jac(2,k2) + dRdxii(k1,2)*jac(1,k2)
!                        BConv(5,(k1-1)*3+k2) = dRdxii(k1,1)*jac(3,k2) + dRdxii(k1,3)*jac(1,k2)
!                        BConv(6,(k1-1)*3+k2) = dRdxii(k1,2)*jac(3,k2) + dRdxii(k1,3)*jac(2,k2)
!                    end do
!                end do
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
!                if(com1 == .true.) BNG(1,:) = 0.0d0
!                if(com2 == .true.) BNG(2,:) = 0.0d0
!                !if(com3 == .true.) BNG(3,:) = 0.0d0
!                if(com4 == .true.) BNG(4,:) = 0.0d0
!                if(com5 == .true.) BNG(5,:) = 0.0d0
!                if(com6 == .true.) BNG(6,:) = 0.0d0
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
!                A = CA
!                B = CB
!                
!                shpANS = 0.0d0
!                !-------------------------------------------------------------------
!                shpANS(1)  = eta/(4.0d0*B)*(eta/B+1.0d0)*(1.0d0+xi/A) 
!                shpANS(2)  = eta/(4.0d0*B)*(eta/B+1.0d0)*(1.0d0-xi/A)
!                shpANS(3)  = 1.0d0/2.0d0*(1.0d0-eta*eta/(B*B))*(1.0d0+xi/A)
!                shpANS(4)  = 1.0d0/2.0d0*(1.0d0-eta*eta/(B*B))*(1.0d0-xi/A)
!                shpANS(5)  = eta/(4.0d0*B)*(eta/B-1.0d0)*(1.0d0+xi/A) 
!                shpANS(6)  = eta/(4.0d0*B)*(eta/B-1.0d0)*(1.0d0-xi/A)
!                !-------------------------------------------------------------------
!                !-------------------------------------------------------------------
!                shpANS(7)  = xi/(4.0d0*B)*(xi/B+1.0d0)*(1.0d0+eta/A)
!                shpANS(8)  = 1.0d0/2.0d0*(1.0d0-xi*xi/(B*B))*(1.0d0+eta/A)
!                shpANS(9)  = xi/(4.0d0*B)*(xi/B-1.0d0)*(1.0d0+eta/A)
!                shpANS(10) = xi/(4.0d0*B)*(xi/B+1.0d0)*(1.0d0-eta/A)
!                shpANS(11) = 1.0d0/2.0d0*(1.0d0-xi*xi/(B*B))*(1.0d0-eta/A)
!                shpANS(12) = xi/(4.0d0*B)*(xi/B-1.0d0)*(1.0d0-eta/A)
!                !-------------------------------------------------------------------
!                
!                !-------------------------------------------------------------------
!                shpANS(13) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0+eta/A)
!                shpANS(14) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0+eta/A)
!                shpANS(15) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0-eta/A)
!                shpANS(16) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0-eta/A)
!                
!                
!                !New knot vectors --------------------------------------------------
!                ub = u_knot
!                vb = v_knot
!                
!                do k1=2,ncpx+p
!                    ubr(k1-1) = ub(k1)
!                end do
!                
!                do k1=2,ncpy+q
!                    vbr(k1-1) = vb(k1)
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
!                    uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(k1,1) + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                    vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(k1,2) + (v_knot(nj+1) + v_knot(nj)))/2.0d0
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
!                    call ShapeFunc3(nel,TyPt(k1,1),TyPt(k1,2),TyPt(k1,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                    jacA = transpose(jacA)
!                    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(1,(k2-1)*3+k3) = BNG(1,(k2-1)*3+k3) + dRAdxii(k2,1)*jacA(1,k3)*MMult(1,k1)
!                            BNG(5,(k2-1)*3+k3) = BNG(5,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(3,k3) + dRAdxii(k2,3)*jacA(1,k3))*MMult(1,k1)
!                        end do
!                    end do
!                end do
!                
!                
!                
!                !Second TP set cycle ------------------------------------------------
!                !Mass Matrix -------------------------------------------------------
!                MTP = 0.0d0
!                do k1=1,6
!                    uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(k1+6,1) + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                    vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(k1+6,2) + (v_knot(nj+1) + v_knot(nj)))/2.0d0
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
!                    call ShapeFunc3(nel,TyPt(k1+6,1),TyPt(k1+6,2),TyPt(k1+6,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                    jacA = transpose(jacA)
!                    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(2,(k2-1)*3+k3) = BNG(2,(k2-1)*3+k3) + dRAdxii(k2,2)*jacA(2,k3)*MMult(1,k1)
!                            BNG(6,(k2-1)*3+k3) = BNG(6,(k2-1)*3+k3) + (dRAdxii(k2,2)*jacA(3,k3) + dRAdxii(k2,3)*jacA(2,k3))*MMult(1,k1)
!                        end do
!                    end do
!                end do
!                
!                
!                
!                !Third TP set cycle ------------------------------------------------
!                !Mass Matrix -------------------------------------------------------
!                MTPr = 0.0d0
!                do k1=1,4
!                    uk = ((u_knot(ni+1) - u_knot(ni))*TyPt(k1+12,1) + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                    vk = ((v_knot(nj+1) - v_knot(nj))*TyPt(k1+12,2) + (v_knot(nj+1) + v_knot(nj)))/2.0d0
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
!                    call ShapeFunc3(nel,TyPt(k1+12,1),TyPt(k1+12,2),TyPt(k1+12,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tzero)
!                    jacA = transpose(jacA)
!                    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        do k3 =1,nds
!                            BNG(4,(k2-1)*3+k3) = BNG(4,(k2-1)*3+k3) + (dRAdxii(k2,1)*jacA(2,k3) + dRAdxii(k2,2)*jacA(1,k3))*MMultr(1,k1)
!                        end do
!                    end do
!                end do
!                !-------------------------------------------------------------------------------------------
!
!                ! Strain-Displacement Matrix in local frame
!                Bloc=matmul(TCL,BNG)
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
!    continue
!
!end subroutine ElHex27ANS