subroutine ElHex27_ProjVol(nel,nelp,Ke,el_ddisp,Finte,inc)
    
    use Mod_Variables
      
    implicit none
    
    integer(4)::cpi,count
    integer(4),parameter::npi_xi= 3, npi_eta = 3, npi_zeta = 3
    integer(4),parameter::totaldof = 81
      
    integer(4),intent(IN)::nel,inc,nelp
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::el_ddisp
    real(8),dimension(totaldof,totaldof),intent(OUT)::Ke
    real(8),dimension(totaldof,1),intent(OUT)::Finte
      
    integer(4)::i,j,k,ni,nj,nk,nel_nza,k1
    real(8)::xi,eta,zeta
    real(8),dimension(npi_xi)::e,we
    real(8),dimension(npi_eta)::n,wn
    real(8),dimension(npi_zeta)::c,wc
    
    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob,Bloc
    
      real(8),dimension(6,6)::matD
      real(8),dimension(ntens,1)::tempstr,deform
      real(8)::density
    
      real(8),dimension((p+1)*(q+1)*(w+1))::R,R_mid
      real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdx,dRdx_mid

      real(8)::detJ,detj_mid,gwt
   
      real(8),dimension(nds,nds)::jac
    
      !Geometric Nonlinearity
      integer(4)::inodes,ilixo
      real(8),dimension(totaldof,1)::tdisp,mdisp,tzero
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
      real(8),dimension(9,totaldof)::BNL
      real(8),dimension(9,9)::MSTR
      real(8),dimension(totaldof,totaldof)::KNLG
      real(8),dimension(nds,1)::dRdenc,dRdxyz
      real(8),dimension(6,totaldof)::Bglob_mid,Bloc_mid
      
      integer(4)::ii,jj,kk,k2,k3,k4,loc_num
      real(8)::CA,CB,detjA
      real(8),dimension(nds,nds)::jacA,dxdxiA,dxdxi,dxdxi_mid
      real(8),dimension(npi,nds)::NatPoints
      real(8),dimension(6)::ShpNat
      real(8),dimension(npi)::Ra
      real(8),dimension(npi,3)::DerShpNat,dRadx,dRAdxi,dRAdxii
      real(8),dimension(npi,3)::dRdxi,dRdxii,dRdxi_mid,dRdxii_mid
      real(8),dimension(8,nds)::TyPt
      real(8)::xib,etab,zetab

      real(8),dimension(p+1)::Nb
      real(8),dimension(q+1)::Mb
      real(8),dimension(w+1)::Lb
      real(8),dimension(p+1)::dNbdxi
      real(8),dimension(q+1)::dMbdeta
      real(8),dimension(w+1)::dLbdzeta

      real(8),dimension(2*(p+1))::ub
      real(8),dimension(2*(q+1))::vb
      real(8),dimension(2*(w+1))::wb
      real(8),dimension(2*p)::ubr
      real(8),dimension(2*q)::vbr
      real(8),dimension(2*w)::wbr

      real(8),dimension(p)::Nbr
      real(8),dimension(q)::Mbr
      real(8),dimension(w)::Lbr
      real(8),dimension(p)::dNbrdxi
      real(8),dimension(q)::dMbrdeta
      real(8),dimension(w)::dLbrdzeta

      real(8),dimension(nds*2,nds*2)::TGC,TGC2
      real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BNat,BANS,BNG,BNGold
      real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BConv,BConv_mid

      real(8)::shp,suma,uk,vk,wk,wei,sumw,vol,A,B
      
      real(8),dimension(3,3)::jact
      real(8),dimension(6,6)::MTP
      
      real(8),dimension(1,6)::M,MMult
      
      real(8),dimension(8,8)::MTPr
      real(8),dimension(1,8)::Mr,MMultr
      
      real(8),dimension(4)::temp4
      real(8),dimension(6)::temp6
      real(8),dimension(8)::temp8
      
      real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bdil,Bdev
      real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bdil_mid,Bdev_mid
      
      
      logical::lsp,endcycle
      logical::com1,com2,com3,com4,com5,com6

      real(8)::yc,zc
      
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
    
      tdisp = 0.0d0
      mdisp= 0.0d0
      TCL = 0.0d0  
      TGL = 0.0d0
      
      if(nlgeom==.true.)then
        tdisp = el_ddisp
        mdisp = el_ddisp/2.0d0
      else
        do i=1,6
            TGL(i,i) = 1.0d0
            TCL(i,i) = 1.0d0
        end do
      end if
    
      !Gauss points cycle
      Ke = 0.0d0
      KNLG = 0.0d0
      BNL = 0.0d0
      MSTR = 0.0d0
    
      Finte = 0.0d0
      do k=1,npi_zeta
        do j=1,npi_eta
            do i=1,npi_xi
            
                cpi = cpi + 1
                
                !MatD(:,:) = TmatD(nel,cpi,:,:)
                
                xi = e(i)
                eta = n(j)
                zeta = c(k)
                
                !ANS Parameters --------------------------------------------------------
                CA = dsqrt(1.0d0/3.0d0)
                CB = dsqrt(3.0d0/5.0d0)
                
                !Tying Points Coordinates
                TyPt = 0.0d0
                
                TyPt(1,1) =  CA; TyPt(1,2) =  CA; TyPt(1,3) =  CA
                TyPt(2,1) =  CA; TyPt(2,2) = -CA; TyPt(2,3) =  CA
                TyPt(3,1) = -CA; TyPt(3,2) =  CA; TyPt(3,3) =  CA
                TyPt(4,1) = -CA; TyPt(4,2) = -CA; TyPt(4,3) =  CA
                
                TyPt(5,1) =  CA; TyPt(5,2) =  CA; TyPt(5,3) =  -CA
                TyPt(6,1) =  CA; TyPt(6,2) = -CA; TyPt(6,3) =  -CA
                TyPt(7,1) = -CA; TyPt(7,2) =  CA; TyPt(7,3) =  -CA
                TyPt(8,1) = -CA; TyPt(8,2) = -CA; TyPt(8,3) =  -CA
                !-----------------------------------------------------------------------
                
                
                call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tzero)
                
                !Weight factor
                gwt = we(i)*wn(j)*wc(k)*detj
                
                !Local Axis ---------------------------
                if(inc==1)then
                    call local_axis(nds,jac,rconv)
                else
                    rconv = laxis(nelp,cpi,:,:)
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
                
                jac_mid = transpose(jac_mid)
                
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
                
                BGlob_mid = 0.0d0
                Bdil_mid = 0.0d0
                do k1=1,(p+1)*(q+1)*(w+1)
                    do k2 =1,nds
                        Bglob_mid(1,k1*3-2) = dRdx_mid(k1,1)
                        Bglob_mid(2,k1*3-1) = dRdx_mid(k1,2)
                        Bglob_mid(3,k1*3  ) = dRdx_mid(k1,3)
                        
                        Bglob_mid(4,k1*3-2) = dRdx_mid(k1,2)
                        Bglob_mid(4,k1*3-1) = dRdx_mid(k1,1)
                        Bglob_mid(5,k1*3-2) = dRdx_mid(k1,3)
                        Bglob_mid(5,k1*3  ) = dRdx_mid(k1,1)
                        Bglob_mid(6,k1*3-1) = dRdx_mid(k1,3)
                        Bglob_mid(6,k1*3  ) = dRdx(k1,2)
                        
                        Bdil_mid(1,k1*3-2) = dRdx_mid(k1,1)/3.0d0
                        Bdil_mid(2,k1*3-2) = dRdx_mid(k1,1)/3.0d0
                        Bdil_mid(3,k1*3-2) = dRdx_mid(k1,1)/3.0d0
                        Bdil_mid(1,k1*3-1) = dRdx_mid(k1,2)/3.0d0
                        Bdil_mid(2,k1*3-1) = dRdx_mid(k1,2)/3.0d0
                        Bdil_mid(3,k1*3-1) = dRdx_mid(k1,2)/3.0d0
                        Bdil_mid(1,k1*3  ) = dRdx_mid(k1,3)/3.0d0
                        Bdil_mid(2,k1*3  ) = dRdx_mid(k1,3)/3.0d0
                        Bdil_mid(3,k1*3  ) = dRdx_mid(k1,3)/3.0d0
                    end do
                end do
                
                Bdev_mid = 0.0d0
                Bdev_mid  = Bglob_mid - Bdil_mid
                
                !-----------------------------------------------------------------------------------------------
                !Mid-configuration -----------------------------------------------------------------------------
                BNG = 0.0d0
                
                ni = INN(IEN(nel,1),1)
                nj = INN(IEN(nel,1),2)
                nk = INN(IEN(nel,1),3)
    
                xib   = ((u_knot(ni+1) - u_knot(ni))*xi   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
                etab  = ((v_knot(nj+1) - v_knot(nj))*eta  + (v_knot(nj+1) + v_knot(nj)))/2.0d0
                zetab = ((w_knot(nk+1) - w_knot(nk))*zeta + (w_knot(nk+1) + w_knot(nk)))/2.0d0
                
                !Bezier knot vectors -----------------------------------------------
                ub = 1.0d0
                vb = 1.0d0
                wb = 1.0d0
                ubr = 1.0d0
                vbr = 1.0d0
                wbr = 1.0d0
                
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
                do k1=1,w+1
                    wb(k1) = 0.0d0
                end do
                do k1=1,w
                    wbr(k1) = 0.0d0
                end do
                
                call BSplineBasisAndDeriv(ncpx,p,xib,ub,Nb,dNbdxi)
                call BSplineBasisAndDeriv(ncpy,q,etab,vb,Mb,dMbdeta)
                call BSplineBasisAndDeriv(ncpz,w,zetab,wb,Lb,dLbdzeta)
                   
                call BSplineBasisAndDeriv(ncpx-1,p-1,xib,ubr,Nbr,dNbrdxi)
                call BSplineBasisAndDeriv(ncpy-1,q-1,etab,vbr,Mbr,dMbrdeta)
                call BSplineBasisAndDeriv(ncpz-1,w-1,zetab,wbr,Lbr,dLbrdzeta)
                
                
                !-------------------------------------------------------------------------------------------------------------------------------
                MTPr = 0.0d0
                do k1=1,8
                    uk = (TyPt(k1,1) + 1.0d0)/2.0d0
                    vk = (TyPt(k1,2) + 1.0d0)/2.0d0
                    wk = (TyPt(k1,3) + 1.0d0)/2.0d0
                    
                    Nbr = 0.0d0
                    call BSplineBasisAndDeriv(ncpx-1,p-1,uk,ubr, Nbr,dNbrdxi)
                    Mbr = 0.0d0
                    call BSplineBasisAndDeriv(ncpy-1,q-1,vk,vbr, Mbr,dMbrdeta)
                    Lbr = 0.0d0
                    call BSplineBasisAndDeriv(ncpz-1,w-1,wk,wbr, Lbr,dLbrdzeta)
                    
                    count = 0
                    do k2=0,w-1
                        do k3=0,q-1
                            do k4=0,p-1
                                count = count + 1
                                MTPr(k1,count) = Nbr(p-k4)*Mbr(q-k3)*Lbr(w-k2)
                            end do
                        end do      
                    end do
                end do
                
                temp8 = 0.0d0
                call gaussj(MTPr,8,temp8,ilixo)
                
                
                !IP Mass Matrix ----------------------------------------------------
                xib   = (  xi + 1.0d0)/2.0d0
                etab  = ( eta + 1.0d0)/2.0d0
                zetab = (zeta + 1.0d0)/2.0d0
                
                call BSplineBasisAndDeriv(ncpx-1,p-1, xib,ubr,Nbr,dNbrdxi)
                call BSplineBasisAndDeriv(ncpy-1,q-1,etab,vbr,Mbr,dMbrdeta)
                call BSplineBasisAndDeriv(ncpz-1,w-1,zetab,wbr,Lbr,dLbrdzeta)
                    
                count = 0
                Mr = 0.0d0
                do k2=0,w-1
                    do k3=0,q-1
                        do k4=0,p-1
                            count = count + 1
                            Mr(1,count) = Nbr(p-k4)*Mbr(q-k3)*Lbr(w-k2)
                        end do
                    end do
                end do 
                
                
                MMultr = matmul(Mr,MTPr)
                continue
                
                do k1=1,8
                    call ShapeFunc3(nel,TyPt(k1,1),TyPt(k1,2),TyPt(k1,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tdisp)
                    jacA = transpose(jacA)
                    
                    do k2=1,(p+1)*(q+1)*(w+1)
                        BNG(1,k2*3-2) = BNG(1,k2*3-2) + dRAdx(k2,1)*MMultr(1,k1)
                        BNG(1,k2*3-1) = BNG(1,k2*3-1) + dRAdx(k2,2)*MMultr(1,k1)
                        BNG(1,k2*3  ) = BNG(1,k2*3  ) + dRAdx(k2,3)*MMultr(1,k1)
                    
                        BNG(2,k2*3-2) = BNG(2,k2*3-2) + dRAdx(k2,1)*MMultr(1,k1)
                        BNG(2,k2*3-1) = BNG(2,k2*3-1) + dRAdx(k2,2)*MMultr(1,k1)
                        BNG(2,k2*3  ) = BNG(2,k2*3  ) + dRAdx(k2,3)*MMultr(1,k1)
                    
                        BNG(3,k2*3-2) = BNG(3,k2*3-2) + dRAdx(k2,1)*MMultr(1,k1)
                        BNG(3,k2*3-1) = BNG(3,k2*3-1) + dRAdx(k2,2)*MMultr(1,k1)
                        BNG(3,k2*3  ) = BNG(3,k2*3  ) + dRAdx(k2,3)*MMultr(1,k1)
                    end do
                end do
                !-------------------------------------------------------------------------------------------------------------------------------
                
                !Update global matrix ---------
                Bglob_mid = Bdev_mid + BNG/3.0d0
                
                !Local strain displacement matrix -------
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
                !call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tdisp)
                call ShapeFunc3(nel,xi,eta,zeta,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,tdisp)
                jac = transpose(jac)
                
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
                
                !Local strain-displacement matrix at end configuration
                Bloc=0.0d0
                
                ni = INN(IEN(nel,1),1)
                nj = INN(IEN(nel,1),2)
                nk = INN(IEN(nel,1),3)
    
                xib   = ((u_knot(ni+1) - u_knot(ni))*xi   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
                etab  = ((v_knot(nj+1) - v_knot(nj))*eta  + (v_knot(nj+1) + v_knot(nj)))/2.0d0
                zetab = ((w_knot(nk+1) - w_knot(nk))*zeta + (w_knot(nk+1) + w_knot(nk)))/2.0d0
                
                !Bezier knot vectors -----------------------------------------------
                BNG = 0.0d0
                
                ni = INN(IEN(nel,1),1)
                nj = INN(IEN(nel,1),2)
                nk = INN(IEN(nel,1),3)
    
                xib   = ((u_knot(ni+1) - u_knot(ni))*xi   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
                etab  = ((v_knot(nj+1) - v_knot(nj))*eta  + (v_knot(nj+1) + v_knot(nj)))/2.0d0
                zetab = ((w_knot(nk+1) - w_knot(nk))*zeta + (w_knot(nk+1) + w_knot(nk)))/2.0d0
                
                !Bezier knot vectors -----------------------------------------------
                ub = 1.0d0
                vb = 1.0d0
                wb = 1.0d0
                ubr = 1.0d0
                vbr = 1.0d0
                wbr = 1.0d0
                
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
                do k1=1,w+1
                    wb(k1) = 0.0d0
                end do
                do k1=1,w
                    wbr(k1) = 0.0d0
                end do
                
                call BSplineBasisAndDeriv(ncpx,p,xib,ub,Nb,dNbdxi)
                call BSplineBasisAndDeriv(ncpy,q,etab,vb,Mb,dMbdeta)
                call BSplineBasisAndDeriv(ncpz,w,zetab,wb,Lb,dLbdzeta)
                   
                call BSplineBasisAndDeriv(ncpx-1,p-1,xib,ubr,Nbr,dNbrdxi)
                call BSplineBasisAndDeriv(ncpy-1,q-1,etab,vbr,Mbr,dMbrdeta)
                call BSplineBasisAndDeriv(ncpz-1,w-1,zetab,wbr,Lbr,dLbrdzeta)
                
                
                !-------------------------------------------------------------------------------------------------------------------------------
                MTPr = 0.0d0
                do k1=1,8
                    uk = (TyPt(k1,1) + 1.0d0)/2.0d0
                    vk = (TyPt(k1,2) + 1.0d0)/2.0d0
                    wk = (TyPt(k1,3) + 1.0d0)/2.0d0
                    
                    Nbr = 0.0d0
                    call BSplineBasisAndDeriv(ncpx-1,p-1,uk,ubr, Nbr,dNbrdxi)
                    Mbr = 0.0d0
                    call BSplineBasisAndDeriv(ncpy-1,q-1,vk,vbr, Mbr,dMbrdeta)
                    Lbr = 0.0d0
                    call BSplineBasisAndDeriv(ncpz-1,w-1,wk,wbr, Lbr,dLbrdzeta)
                    
                    count = 0
                    do k2=0,w-1
                        do k3=0,q-1
                            do k4=0,p-1
                                count = count + 1
                                MTPr(k1,count) = Nbr(p-k4)*Mbr(q-k3)*Lbr(w-k2)
                            end do
                        end do      
                    end do
                end do
                
                temp8 = 0.0d0
                call gaussj(MTPr,8,temp8,ilixo)
                
                
                !IP Mass Matrix ----------------------------------------------------
                xib   = (  xi + 1.0d0)/2.0d0
                etab  = ( eta + 1.0d0)/2.0d0
                zetab = (zeta + 1.0d0)/2.0d0
                
                call BSplineBasisAndDeriv(ncpx-1,p-1, xib,ubr,Nbr,dNbrdxi)
                call BSplineBasisAndDeriv(ncpy-1,q-1,etab,vbr,Mbr,dMbrdeta)
                call BSplineBasisAndDeriv(ncpz-1,w-1,zetab,wbr,Lbr,dLbrdzeta)
                    
                count = 0
                Mr = 0.0d0
                do k2=0,w-1
                    do k3=0,q-1
                        do k4=0,p-1
                            count = count + 1
                            Mr(1,count) = Nbr(p-k4)*Mbr(q-k3)*Lbr(w-k2)
                        end do
                    end do
                end do 
                
                
                MMultr = matmul(Mr,MTPr)
                continue
                
                do k1=1,8
                    call ShapeFunc3(nel,TyPt(k1,1),TyPt(k1,2),TyPt(k1,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tdisp)
                    jacA = transpose(jacA)
                    
                    do k2=1,(p+1)*(q+1)*(w+1)
                        BNG(1,k2*3-2) = BNG(1,k2*3-2) + dRAdx(k2,1)*MMultr(1,k1)
                        BNG(1,k2*3-1) = BNG(1,k2*3-1) + dRAdx(k2,2)*MMultr(1,k1)
                        BNG(1,k2*3  ) = BNG(1,k2*3  ) + dRAdx(k2,3)*MMultr(1,k1)
                    
                        BNG(2,k2*3-2) = BNG(2,k2*3-2) + dRAdx(k2,1)*MMultr(1,k1)
                        BNG(2,k2*3-1) = BNG(2,k2*3-1) + dRAdx(k2,2)*MMultr(1,k1)
                        BNG(2,k2*3  ) = BNG(2,k2*3  ) + dRAdx(k2,3)*MMultr(1,k1)
                    
                        BNG(3,k2*3-2) = BNG(3,k2*3-2) + dRAdx(k2,1)*MMultr(1,k1)
                        BNG(3,k2*3-1) = BNG(3,k2*3-1) + dRAdx(k2,2)*MMultr(1,k1)
                        BNG(3,k2*3  ) = BNG(3,k2*3  ) + dRAdx(k2,3)*MMultr(1,k1)
                    end do
                end do
                !-------------------------------------------------------------------------------------------------------------------------------
                
                !Update global matrix ---------
                Bglob = Bdev + BNG/3.0d0
                
                ! Strain-Displacement Matrix in local frame
                Bloc=0.0d0
                Bloc=matmul(TGL,Bglob)

                if(nlgeom == .false.) Bloc_mid = Bloc
                
                call MatPlastic3D(iprops,nds,(p+1)*(q+1)*(w+1),props,Bloc_mid,6,el_ddisp,stress(nelp,cpi,:),&
                & strain(nelp,cpi,:),dstrain(nelp,cpi,:),dstress(nelp,cpi,:),hard(nelp,cpi),matD)
                
                !Store Constitutive Matrix
                !TmatD(nel,cpi,:,:) = MatD(:,:)
                
                tempstr = 0.0d0
                do k1=1,ntens
                    tempstr(k1,1) = stress(nelp,cpi,k1)
                    deform(k1,1) = dstrain(nelp,cpi,k1)
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
                    laxis(nelp,cpi,:,:) = rconv
 
                end if
                
                if(nlgeom == .false.) then
                    Ke = Ke + matmul(matmul(transpose(Bloc),matD),Bloc)*gwt
                else
                    Ke = Ke + matmul(matmul(transpose(Bloc),matD),Bloc)*gwt + matmul(matmul(transpose(BNL),MSTR),BNL)*gwt
                end if
                
                Finte = Finte + matmul(transpose(Bloc),tempstr)*gwt
                
                !Strain Energy --------------------------------------------------------------
                do k1=1,ntens
                    SEnergy = SEnergy + tempstr(k1,1)*deform(k1,1)*gwt/2.0d0
                end do
                
                continue
                
            end do  
        end do
      end do

      continue

end subroutine ElHex27_ProjVol

!subroutine ElHex27_ProjVol(nel,Ke,el_ddisp,Finte,inc)
!    
!    use Mod_Variables
!      
!    implicit none
!    
!    integer(4)::cpi,count
!    integer(4),parameter::npi_xi= 3, npi_eta = 3, npi_zeta = 3
!    integer(4),parameter::totaldof = 81
!      
!    integer(4),intent(IN)::nel,inc
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::el_ddisp
!    real(8),dimension(totaldof,totaldof),intent(OUT)::Ke
!    real(8),dimension(totaldof,1),intent(OUT)::Finte
!      
!    integer(4)::i,j,k,ni,nj,nk,nel_nza,k1,l
!    real(8)::xi,eta,zeta
!    real(8),dimension(npi_xi)::e,we
!    real(8),dimension(npi_eta)::n,wn
!    real(8),dimension(npi_zeta)::c,wc
!    
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob,Bloc
!    
!      real(8),dimension(6,6)::matD
!      real(8),dimension(ntens,1)::tempstr,deform,STDstrain,ENHstrain
!      real(8)::density
!    
!      real(8),dimension((p+1)*(q+1)*(w+1))::R,R_mid
!      real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdx,dRdx_mid
!
!      real(8)::detJ,detj_mid,gwt
!   
!      real(8),dimension(nds,nds)::jac
!    
!      !Geometric Nonlinearity
!      integer(4)::inodes,ilixo
!      real(8),dimension(totaldof,1)::tdisp,mdisp,tzero
!      real(8),dimension(nds,nds)::rconv,rconv_mid,jacinv,raux
!      real(8),dimension(nds,nds)::jac_mid,jacinv_mid
!      real(8),dimension(nds,nds)::matF,matF_mid,matR,matR_mid
!      real(8),dimension((p+1)*(q+1)*(w+1))::dRde,dRdn,dRdc
!      real(8),dimension(nds)::temp1
!      real(8),dimension((p+1)*(q+1)*(w+1),nds)::nodes
!      real(8),dimension(nds,nds)::Temp33_mid,Temp33
!      real(8),dimension(nds,nds)::Trans_mid,Trans
!      real(8),dimension(nds*2,nds*2)::TCL_mid,TGL_mid
!      real(8),dimension(nds*2,nds*2)::TCL,TGL
!      real(8),dimension(9,totaldof)::BNL
!      real(8),dimension(9,9)::MSTR
!      real(8),dimension(totaldof,totaldof)::KNLG
!      real(8),dimension(nds,1)::dRdenc,dRdxyz
!      real(8),dimension(6,totaldof)::Bglob_mid,Bloc_mid
!      
!      integer(4)::ii,jj,kk,k2,k3,k4,loc_num
!      real(8)::CA,CB,detjA
!      real(8),dimension(nds,nds)::jacA,dxdxiA,dxdxi,dxdxi_mid
!      real(8),dimension(npi,nds)::NatPoints
!      real(8),dimension(6)::ShpNat
!      real(8),dimension(npi)::Ra
!      real(8),dimension(npi,3)::DerShpNat,dRadx,dRAdxi,dRAdxii
!      real(8),dimension(npi,3)::dRdxi,dRdxii,dRdxi_mid,dRdxii_mid
!      real(8),dimension(8,nds)::TyPt
!      real(8)::xib,etab,zetab
!
!      real(8),dimension(p+1)::Nb
!      real(8),dimension(q+1)::Mb
!      real(8),dimension(w+1)::Lb
!      real(8),dimension(p+1)::dNbdxi
!      real(8),dimension(q+1)::dMbdeta
!      real(8),dimension(w+1)::dLbdzeta
!
!      real(8),dimension(2*(p+1))::ub
!      real(8),dimension(2*(q+1))::vb
!      real(8),dimension(2*(w+1))::wb
!      real(8),dimension(2*p)::ubr
!      real(8),dimension(2*q)::vbr
!      real(8),dimension(2*w)::wbr
!
!      real(8),dimension(p)::Nbr
!      real(8),dimension(q)::Mbr
!      real(8),dimension(w)::Lbr
!      real(8),dimension(p)::dNbrdxi
!      real(8),dimension(q)::dMbrdeta
!      real(8),dimension(w)::dLbrdzeta
!
!      real(8),dimension(nds*2,nds*2)::TGC,TGC2
!      real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BNat,BANS,BNG,BNGold
!      real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BConv,BConv_mid
!
!      real(8)::shp,suma,uk,vk,wk,wei,sumw,vol,A,B
!      
!      real(8),dimension(3,3)::jact
!      real(8),dimension(6,6)::MTP
!      
!      real(8),dimension(1,6)::M,MMult
!      
!      real(8),dimension(8,8)::MTPr
!      real(8),dimension(1,8)::Mr,MMultr
!      
!      real(8),dimension(4)::temp4
!      real(8),dimension(6)::temp6
!      real(8),dimension(8)::temp8
!      
!      real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bdil,Bdev
!      real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bdil_mid,Bdev_mid
!      
!      integer(4),parameter:: nenh = 6
!      real(8)::detj0
!      real(8),dimension(nenh)::tempa,diagon
!      real(8),dimension(nds*2,nds*2)::MatT0
!      real(8),dimension(6,nenh)::Malpha,Balpha,MAT
!      real(8),dimension((p+1)*(q+1)*(w+1)*nds,nenh)::Kua
!      real(8),dimension(nenh,nenh)::Kaa,Kaainv
!      real(8),dimension(nenh,1):: dalpha
!      
!      
!      logical::lsp,endcycle
!      logical::com1,com2,com3,com4,com5,com6
!
!      real(8)::yc,zc
!      
!      inodes = (p+1)*(q+1)*(w+1)
!      density = props(3)
!    
!      tzero = 0.0d0
!    
!      cpi = 0
!
!      we = 1.0d0
!      wn = 1.0d0
!      wc = 1.0d0
!      call gauleg(npi_xi, e, we)
!      call gauleg(npi_eta, n, wn)
!      call gauleg(npi_zeta, c, wc)
!    
!      tdisp = 0.0d0
!      mdisp= 0.0d0
!      TCL = 0.0d0  
!      TGL = 0.0d0
!      
!      if(nlgeom==.true.)then
!        tdisp = el_ddisp
!        mdisp = el_ddisp/2.0d0
!      else
!        do i=1,6
!            TGL(i,i) = 1.0d0
!            TCL(i,i) = 1.0d0
!        end do
!      end if
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
!    jacinv = transpose(jac)
!    call gaussj(jacinv,nds,temp1,ilixo)
!    
!    call TransformationMat3D(JacInv,MatT0)  
!    
!    dalpha = 0.0d0    
!    do l=1,2
!        
!      !Gauss points cycle
!      Ke = 0.0d0
!      KNLG = 0.0d0
!      BNL = 0.0d0
!      MSTR = 0.0d0
!      Kua = 0.0d0
!      Kaa = 0.0d0
!      
!      cpi = 0
!    
!      Finte = 0.0d0
!      do k=1,npi_zeta
!        do j=1,npi_eta
!            do i=1,npi_xi
!            
!                cpi = cpi + 1
!                
!                MatD(:,:) = TmatD(nel,cpi,:,:)
!                
!                xi = e(i)
!                eta = n(j)
!                zeta = c(k)
!                
!                !ANS Parameters --------------------------------------------------------
!                CA = dsqrt(1.0d0/3.0d0)
!                CB = dsqrt(3.0d0/5.0d0)
!                
!                !Tying Points Coordinates
!                TyPt = 0.0d0
!                
!                TyPt(1,1) =  CA; TyPt(1,2) =  CA; TyPt(1,3) =  CA
!                TyPt(2,1) =  CA; TyPt(2,2) = -CA; TyPt(2,3) =  CA
!                TyPt(3,1) = -CA; TyPt(3,2) =  CA; TyPt(3,3) =  CA
!                TyPt(4,1) = -CA; TyPt(4,2) = -CA; TyPt(4,3) =  CA
!                
!                TyPt(5,1) =  CA; TyPt(5,2) =  CA; TyPt(5,3) =  -CA
!                TyPt(6,1) =  CA; TyPt(6,2) = -CA; TyPt(6,3) =  -CA
!                TyPt(7,1) = -CA; TyPt(7,2) =  CA; TyPt(7,3) =  -CA
!                TyPt(8,1) = -CA; TyPt(8,2) = -CA; TyPt(8,3) =  -CA
!                !-----------------------------------------------------------------------
!                
!                
!                call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tzero)
!                
!                !Weight factor
!                gwt = we(i)*wn(j)*wc(k)*detj
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
!                    nodes(k1,1) = Points(IEN(nel,count),1)
!                    nodes(k1,2) = Points(IEN(nel,count),2)
!                    nodes(k1,3) = Points(IEN(nel,count),3)
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
!                jac_mid = transpose(jac_mid)
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
!                BGlob_mid = 0.0d0
!                Bdil_mid = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
!                    do k2 =1,nds
!                        Bglob_mid(1,k1*3-2) = dRdx_mid(k1,1)
!                        Bglob_mid(2,k1*3-1) = dRdx_mid(k1,2)
!                        Bglob_mid(3,k1*3  ) = dRdx_mid(k1,3)
!                        
!                        Bglob_mid(4,k1*3-2) = dRdx_mid(k1,2)
!                        Bglob_mid(4,k1*3-1) = dRdx_mid(k1,1)
!                        Bglob_mid(5,k1*3-2) = dRdx_mid(k1,3)
!                        Bglob_mid(5,k1*3  ) = dRdx_mid(k1,1)
!                        Bglob_mid(6,k1*3-1) = dRdx_mid(k1,3)
!                        Bglob_mid(6,k1*3  ) = dRdx(k1,2)
!                        
!                        Bdil_mid(1,k1*3-2) = dRdx_mid(k1,1)/3.0d0
!                        Bdil_mid(2,k1*3-2) = dRdx_mid(k1,1)/3.0d0
!                        Bdil_mid(3,k1*3-2) = dRdx_mid(k1,1)/3.0d0
!                        Bdil_mid(1,k1*3-1) = dRdx_mid(k1,2)/3.0d0
!                        Bdil_mid(2,k1*3-1) = dRdx_mid(k1,2)/3.0d0
!                        Bdil_mid(3,k1*3-1) = dRdx_mid(k1,2)/3.0d0
!                        Bdil_mid(1,k1*3  ) = dRdx_mid(k1,3)/3.0d0
!                        Bdil_mid(2,k1*3  ) = dRdx_mid(k1,3)/3.0d0
!                        Bdil_mid(3,k1*3  ) = dRdx_mid(k1,3)/3.0d0
!                    end do
!                end do
!                
!                Bdev_mid = 0.0d0
!                Bdev_mid  = Bglob_mid - Bdil_mid
!                
!                !-----------------------------------------------------------------------------------------------
!                !Mid-configuration -----------------------------------------------------------------------------
!                BNG = 0.0d0
!                
!                ni = INN(IEN(nel,1),1)
!                nj = INN(IEN(nel,1),2)
!                nk = INN(IEN(nel,1),3)
!    
!                xib   = ((u_knot(ni+1) - u_knot(ni))*xi   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                etab  = ((v_knot(nj+1) - v_knot(nj))*eta  + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                zetab = ((w_knot(nk+1) - w_knot(nk))*zeta + (w_knot(nk+1) + w_knot(nk)))/2.0d0
!                
!                !Bezier knot vectors -----------------------------------------------
!                ub = 1.0d0
!                vb = 1.0d0
!                wb = 1.0d0
!                ubr = 1.0d0
!                vbr = 1.0d0
!                wbr = 1.0d0
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
!                do k1=1,w+1
!                    wb(k1) = 0.0d0
!                end do
!                do k1=1,w
!                    wbr(k1) = 0.0d0
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)
!                call BSplineBasisAndDeriv(ncpz+w,w,zetab,wb,Lb,dLbdzeta)
!                   
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                call BSplineBasisAndDeriv(ncpz+w,w-1,zetab,wbr,Lbr,dLbrdzeta)
!                
!                
!                !-------------------------------------------------------------------------------------------------------------------------------
!                MTPr = 0.0d0
!                do k1=1,8
!                    uk = (TyPt(k1,1) + 1.0d0)/2.0d0
!                    vk = (TyPt(k1,2) + 1.0d0)/2.0d0
!                    wk = (TyPt(k1,3) + 1.0d0)/2.0d0
!                    
!                    Nbr = 0.0d0
!                    call BSplineBasisAndDeriv(ncpx+2*p+1,p-1,uk,ubr, Nbr,dNbrdxi)
!                    Mbr = 0.0d0
!                    call BSplineBasisAndDeriv(ncpy+2*q+1,q-1,vk,vbr, Mbr,dMbrdeta)
!                    Lbr = 0.0d0
!                    call BSplineBasisAndDeriv(ncpz+2*w+1,w-1,wk,wbr, Lbr,dLbrdzeta)
!                    
!                    count = 0
!                    do k2=0,w-1
!                        do k3=0,q-1
!                            do k4=0,p-1
!                                count = count + 1
!                                MTPr(k1,count) = Nbr(p-k4)*Mbr(q-k3)*Lbr(w-k2)
!                            end do
!                        end do      
!                    end do
!                end do
!                
!                temp8 = 0.0d0
!                call gaussj(MTPr,8,temp8,ilixo)
!                
!                
!                !IP Mass Matrix ----------------------------------------------------
!                xib   = (  xi + 1.0d0)/2.0d0
!                etab  = ( eta + 1.0d0)/2.0d0
!                zetab = (zeta + 1.0d0)/2.0d0
!                
!                call BSplineBasisAndDeriv(ncpx+2*p,p-1, xib,ubr,Nbr,dNbrdxi)
!                call BSplineBasisAndDeriv(ncpy+2*q,q-1,etab,vbr,Mbr,dMbrdeta)
!                call BSplineBasisAndDeriv(ncpz+2*w,w-1,zetab,wbr,Lbr,dLbrdzeta)
!                    
!                count = 0
!                Mr = 0.0d0
!                do k2=0,w-1
!                    do k3=0,q-1
!                        do k4=0,p-1
!                            count = count + 1
!                            Mr(1,count) = Nbr(p-k4)*Mbr(q-k3)*Lbr(w-k2)
!                        end do
!                    end do
!                end do 
!                
!                
!                MMultr = matmul(Mr,MTPr)
!                continue
!                
!                do k1=1,8
!                    call ShapeFunc3(nel,TyPt(k1,1),TyPt(k1,2),TyPt(k1,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tdisp)
!                    jacA = transpose(jacA)
!                    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        BNG(1,k2*3-2) = BNG(1,k2*3-2) + dRAdx(k2,1)*MMultr(1,k1)
!                        BNG(1,k2*3-1) = BNG(1,k2*3-1) + dRAdx(k2,2)*MMultr(1,k1)
!                        BNG(1,k2*3  ) = BNG(1,k2*3  ) + dRAdx(k2,3)*MMultr(1,k1)
!                    
!                        BNG(2,k2*3-2) = BNG(2,k2*3-2) + dRAdx(k2,1)*MMultr(1,k1)
!                        BNG(2,k2*3-1) = BNG(2,k2*3-1) + dRAdx(k2,2)*MMultr(1,k1)
!                        BNG(2,k2*3  ) = BNG(2,k2*3  ) + dRAdx(k2,3)*MMultr(1,k1)
!                    
!                        BNG(3,k2*3-2) = BNG(3,k2*3-2) + dRAdx(k2,1)*MMultr(1,k1)
!                        BNG(3,k2*3-1) = BNG(3,k2*3-1) + dRAdx(k2,2)*MMultr(1,k1)
!                        BNG(3,k2*3  ) = BNG(3,k2*3  ) + dRAdx(k2,3)*MMultr(1,k1)
!                    end do
!                end do
!                !-------------------------------------------------------------------------------------------------------------------------------
!                
!                !Update global matrix ---------
!                Bglob_mid = Bdev_mid + BNG/3.0d0
!                
!                !Local strain displacement matrix -------
!                Bloc_mid=0.0d0
!                Bloc_mid=matmul(TGL_mid,Bglob_mid)
!                
!                
!                
!                
!                !-------------------------------------------------------------------------
!                !End-Configuration
!                !-------------------------------------------------------------------------
!                
!                !Deformation Gradient for mid-configuration
!                call DefGrad3D(inodes,nds,dRde,dRdn,dRdc,jacinv,nodes,tdisp,MatF)
!                
!                !Polar Decomposition for mid configuration
!                matR = 0.0d0
!                call PolarDecomp2D(nds,matF,matR)
!                
!                !Update local reference system to the mid configuration
!                raux = rconv
!                call axisupdate(rconv,matR,rconv)
!                
!                !Global to natural jacobian in mid configuration
!                !call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac,tdisp)
!                call ShapeFunc3(nel,xi,eta,zeta,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,tdisp)
!                jac = transpose(jac)
!                
!                jacinv = 0.0d0
!                temp1 = 0.0d0
!                jacinv = jac
!                call gaussj(jacinv,nds,temp1,ilixo)
!                
!                !Transformation matrices mid configuration 
!                Temp33 = 0.0d0
!                Temp33 = matmul(transpose(rconv),jacinv)
!                Trans = 0.0d0
!                Trans = transpose(rconv)
!                
!                !Natural to local transformation matrix for mid configuration
!                call TransformationMat3D(Temp33,TCL)
!                
!                !Global to local transformation matrix for mid configuration
!                call TransformationMat3D(Trans,TGL)
!                
!                !Update Weight factor
!                gwt = we(i)*wn(j)*wc(k)*detj
!                
!                continue
!                
!                Bglob = 0.0d0
!                Bdil = 0.0d0
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
!                    Bdil(1,k1*3-2) = dRdx(k1,1)/3.0d0
!                    Bdil(2,k1*3-2) = dRdx(k1,1)/3.0d0
!                    Bdil(3,k1*3-2) = dRdx(k1,1)/3.0d0
!                    Bdil(1,k1*3-1) = dRdx(k1,2)/3.0d0
!                    Bdil(2,k1*3-1) = dRdx(k1,2)/3.0d0
!                    Bdil(3,k1*3-1) = dRdx(k1,2)/3.0d0
!                    Bdil(1,k1*3  ) = dRdx(k1,3)/3.0d0
!                    Bdil(2,k1*3  ) = dRdx(k1,3)/3.0d0
!                    Bdil(3,k1*3  ) = dRdx(k1,3)/3.0d0
!                    
!                end do
!                
!                Bdev = 0.0d0
!                Bdev  = Bglob - Bdil
!
!                BNG = 0.0d0
!                
!                !Local strain-displacement matrix at end configuration
!                Bloc=0.0d0
!                
!                ni = INN(IEN(nel,1),1)
!                nj = INN(IEN(nel,1),2)
!                nk = INN(IEN(nel,1),3)
!    
!                xib   = ((u_knot(ni+1) - u_knot(ni))*xi   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                etab  = ((v_knot(nj+1) - v_knot(nj))*eta  + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                zetab = ((w_knot(nk+1) - w_knot(nk))*zeta + (w_knot(nk+1) + w_knot(nk)))/2.0d0
!                
!                !Bezier knot vectors -----------------------------------------------
!                BNG = 0.0d0
!                
!                ni = INN(IEN(nel,1),1)
!                nj = INN(IEN(nel,1),2)
!                nk = INN(IEN(nel,1),3)
!    
!                xib   = ((u_knot(ni+1) - u_knot(ni))*xi   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
!                etab  = ((v_knot(nj+1) - v_knot(nj))*eta  + (v_knot(nj+1) + v_knot(nj)))/2.0d0
!                zetab = ((w_knot(nk+1) - w_knot(nk))*zeta + (w_knot(nk+1) + w_knot(nk)))/2.0d0
!                
!                !Bezier knot vectors -----------------------------------------------
!                ub = 1.0d0
!                vb = 1.0d0
!                wb = 1.0d0
!                ubr = 1.0d0
!                vbr = 1.0d0
!                wbr = 1.0d0
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
!                do k1=1,w+1
!                    wb(k1) = 0.0d0
!                end do
!                do k1=1,w
!                    wbr(k1) = 0.0d0
!                end do
!                
!                call BSplineBasisAndDeriv(ncpx+p,p,xib,ub,Nb,dNbdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q,etab,vb,Mb,dMbdeta)
!                call BSplineBasisAndDeriv(ncpz+w,w,zetab,wb,Lb,dLbdzeta)
!                   
!                call BSplineBasisAndDeriv(ncpx+p,p-1,xib,ubr,Nbr,dNbrdxi)
!                call BSplineBasisAndDeriv(ncpy+q,q-1,etab,vbr,Mbr,dMbrdeta)
!                call BSplineBasisAndDeriv(ncpz+w,w-1,zetab,wbr,Lbr,dLbrdzeta)
!                
!                
!                !-------------------------------------------------------------------------------------------------------------------------------
!                MTPr = 0.0d0
!                do k1=1,8
!                    uk = (TyPt(k1,1) + 1.0d0)/2.0d0
!                    vk = (TyPt(k1,2) + 1.0d0)/2.0d0
!                    wk = (TyPt(k1,3) + 1.0d0)/2.0d0
!                    
!                    Nbr = 0.0d0
!                    call BSplineBasisAndDeriv(ncpx+2*p+1,p-1,uk,ubr, Nbr,dNbrdxi)
!                    Mbr = 0.0d0
!                    call BSplineBasisAndDeriv(ncpy+2*q+1,q-1,vk,vbr, Mbr,dMbrdeta)
!                    Lbr = 0.0d0
!                    call BSplineBasisAndDeriv(ncpz+2*w+1,w-1,wk,wbr, Lbr,dLbrdzeta)
!                    
!                    count = 0
!                    do k2=0,w-1
!                        do k3=0,q-1
!                            do k4=0,p-1
!                                count = count + 1
!                                MTPr(k1,count) = Nbr(p-k4)*Mbr(q-k3)*Lbr(w-k2)
!                            end do
!                        end do      
!                    end do
!                end do
!                
!                temp8 = 0.0d0
!                call gaussj(MTPr,8,temp8,ilixo)
!                
!                
!                !IP Mass Matrix ----------------------------------------------------
!                xib   = (  xi + 1.0d0)/2.0d0
!                etab  = ( eta + 1.0d0)/2.0d0
!                zetab = (zeta + 1.0d0)/2.0d0
!                
!                call BSplineBasisAndDeriv(ncpx+2*p+1,p-1, xib,ubr,Nbr,dNbrdxi)
!                call BSplineBasisAndDeriv(ncpy+2*q+1,q-1,etab,vbr,Mbr,dMbrdeta)
!                call BSplineBasisAndDeriv(ncpz+2*w+1,w-1,zetab,wbr,Lbr,dLbrdzeta)
!                    
!                count = 0
!                Mr = 0.0d0
!                do k2=0,w-1
!                    do k3=0,q-1
!                        do k4=0,p-1
!                            count = count + 1
!                            Mr(1,count) = Nbr(p-k4)*Mbr(q-k3)*Lbr(w-k2)
!                        end do
!                    end do
!                end do 
!                
!                
!                MMultr = matmul(Mr,MTPr)
!                continue
!                
!                do k1=1,8
!                    call ShapeFunc3(nel,TyPt(k1,1),TyPt(k1,2),TyPt(k1,3),RA,dRAdx,dRAdxi,dRAdxii,detjA,jacA,dxdxiA,tdisp)
!                    jacA = transpose(jacA)
!                    
!                    do k2=1,(p+1)*(q+1)*(w+1)
!                        BNG(1,k2*3-2) = BNG(1,k2*3-2) + dRAdx(k2,1)*MMultr(1,k1)
!                        BNG(1,k2*3-1) = BNG(1,k2*3-1) + dRAdx(k2,2)*MMultr(1,k1)
!                        BNG(1,k2*3  ) = BNG(1,k2*3  ) + dRAdx(k2,3)*MMultr(1,k1)
!                    
!                        BNG(2,k2*3-2) = BNG(2,k2*3-2) + dRAdx(k2,1)*MMultr(1,k1)
!                        BNG(2,k2*3-1) = BNG(2,k2*3-1) + dRAdx(k2,2)*MMultr(1,k1)
!                        BNG(2,k2*3  ) = BNG(2,k2*3  ) + dRAdx(k2,3)*MMultr(1,k1)
!                    
!                        BNG(3,k2*3-2) = BNG(3,k2*3-2) + dRAdx(k2,1)*MMultr(1,k1)
!                        BNG(3,k2*3-1) = BNG(3,k2*3-1) + dRAdx(k2,2)*MMultr(1,k1)
!                        BNG(3,k2*3  ) = BNG(3,k2*3  ) + dRAdx(k2,3)*MMultr(1,k1)
!                    end do
!                end do
!                !-------------------------------------------------------------------------------------------------------------------------------
!                
!                !Update global matrix ---------
!                Bglob = Bdev + BNG/3.0d0
!                
!                ! Strain-Displacement Matrix in local frame
!                Bloc=0.0d0
!                Bloc=matmul(TGL,Bglob)
!
!                
!                !-------------------------------------------------------
!                ! EAS
!                !-------------------------------------------------------
!                
!                !First Test ---------------------------------
!                Malpha = 0.0d0
!                Malpha(1,1) = 3.0d0*xi**2.0d0 - 1.0d0
!                Malpha(2,2) = 3.0d0*eta**2.0d0 - 1.0d0
!                Malpha(3,3) = 3.0d0*zeta**2.0d0 - 1.0d0
!                
!                Malpha(1,4) = xi
!                Malpha(2,5) = eta
!                Malpha(3,6) = zeta
!                
!                !Second Test ---------------------------------
!!                Malpha = 0.0d0
!!                Malpha(1,1) = xi
!!                Malpha(2,2) = eta
!!                Malpha(3,3) = zeta
!!                Malpha(1,4) = xi*eta
!!                Malpha(1,5) = xi*zeta
!!                Malpha(2,6) = xi*eta
!!                Malpha(2,7) = eta*zeta
!!                Malpha(3,8) = xi*zeta
!!                Malpha(3,9) = eta*zeta
!!!                
!!!                Malpha(1,13) = xi*eta*zeta
!!!                Malpha(2,14) = xi*eta*zeta
!!!                Malpha(3,15) = xi*eta*zeta
!!                
!                Balpha = 0.0d0
!                Balpha = matmul(MatT0,Malpha)*detj0/detj
!
!                if(nlgeom == .false.) Bloc_mid = Bloc
!                
!                !call MatPlastic3D(iprops,nds,(p+1)*(q+1)*(w+1),props,Bloc_mid,6,el_ddisp,stress(nel,cpi,:),strain(nel,cpi,:),dstrain(nel,cpi,:),dstress(nel,cpi,:),hard(nel,cpi),matD)
!                
!                !Calculate standard strain --------------------------
!	            STDStrain=matmul(Bloc,el_ddisp)
!				
!		        !Calculate enhanced strain --------------------------
!		        ENHStrain=matmul(Balpha,dalpha)
!
!		        deform = STDstrain + ENHstrain
!        
!                tempstr=matmul(Matd,deform)
!        
!                !Store Constitutive Matrix
!                !TmatD(nel,cpi,:,:) = MatD(:,:)
!                
!                !tempstr = 0.0d0
!                do k1=1,ntens
!                    stress(nel,cpi,k1) = tempstr(k1,1)
!                    dstrain(nel,cpi,k1) = deform(k1,1)
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
!                Kua = Kua + matmul(matmul(transpose(Bloc),MatD),Balpha)*gwt
!                Kaa = Kaa + matmul(matmul(transpose(Balpha),MatD),Balpha)*gwt
!                
!                if(nlgeom == .false.) then
!                    Ke = Ke + matmul(matmul(transpose(Bloc),matD),Bloc)*gwt
!                else
!                    Ke = Ke + matmul(matmul(transpose(Bloc),matD),Bloc)*gwt + matmul(matmul(transpose(BNL),MSTR),BNL)*gwt
!                end if
!                
!                Finte = Finte + matmul(transpose(Bloc),tempstr)*gwt
!                
!                !Strain Energy --------------------------------------------------------------
!                if(l==2)then
!                    do k1=1,ntens
!                        SEnergy = SEnergy + tempstr(k1,1)*deform(k1,1)*gwt/2.0d0
!                    end do
!                end if
!                
!                continue
!                
!            end do  
!        end do
!      end do
!
!      tempa = 0.0d0
!      Kaainv = 0.0d0
!      Kaainv = Kaa
!      call gaussj(Kaainv,nenh,tempa,ilixo)
!    
!      Ke = Ke - matmul(matmul(Kua,Kaainv),transpose(Kua))
!      
!      dalpha=0.0d0
!      dalpha=matmul(matmul(Kaainv,transpose(Kua)),el_ddisp)
!	  dalpha=dalpha*(-1.0d0)
!      
!    end do
!    
!    continue
!
!end subroutine ElHex27_ProjVol
