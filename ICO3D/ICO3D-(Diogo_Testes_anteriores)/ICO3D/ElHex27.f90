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

!subroutine ElHex27(nel,Ke,el_ddisp,Finte)
!    
!    use Mod_Variables
!    implicit none
!    
!    integer(4)::nel,count
!    integer(4)::cpi
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
!    real(8),dimension((p+1)*(q+1)*(w+1))::R
!    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdx
!
!    real(8)::detJ,gwt
!    
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob,Bloc
!    
!    real(8),dimension(6,6)::matD
!    real(8),dimension(ntens,1)::tempstr,deform
!    
!    real(8)::density
!    real(8),dimension(nds,nds)::jac
!    
!    density = props(3)
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
!    !Gauss points cycle
!    Ke = 0.0d0
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
!                call ShapeFunc(nel,xi,eta,zeta,R,dRdx,detj)
!                !call ShapeFunc2(nel,xi,eta,zeta,R,dRdx,detj,jac)
!                
!                !Weight factor
!                gwt = we(i)*wn(j)*wc(k)*detj
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
!                Bloc = Bglob
!                
!                call MatPlastic3D(iprops,nds,(p+1)*(q+1)*(w+1),props,Bloc,6,el_ddisp,stress(nel,cpi,:),strain(nel,cpi,:),dstrain(nel,cpi,:),dstress(nel,cpi,:),hard(nel,cpi),matD)
!                
!                tempstr = 0.0d0
!                do k1=1,ntens
!                    tempstr(k1,1) = stress(nel,cpi,k1)
!                    deform(k1,1) = dstrain(nel,cpi,k1)
!                end do
!                
!                Ke = Ke + matmul(matmul(transpose(Bglob),matD),Bglob)*gwt
!                
!                Finte = Finte + matmul(transpose(Bglob),tempstr)*gwt
!                
!                !Gravitic loads -------------------------------------------------------------
!                if (gravity==.true.)then
!                    call grvt((p+1)*(q+1)*(w+1),nds,R,gwt,gconst,gravdir,density,Finte)
!                end if
!                
!                !Pressure loads -------------------------------------------------------------
!                !Need to set the nlgeom first
!!                if(pressure==.true.) then
!!                    do k1=1,nelems
!!                        if(k1==nel .and. abs(pressvec(k1)) .gt. 0.0d0) then
!!                            call prssr(nel,npi_xi,npi_eta,npi_zeta,jac,we(i),wn(j),wc(k),R,Finte)
!!                        end if
!!                    end do
!!                end if
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
!end subroutine ElHex27

!------------------------------------------------------------------------------------------------------------
!
! General second-order hexaedral element with a 3x3x3 gaussian quadrature
!
! Input: nel - Element number
!        el_ddisp - Element displacement
!        inc - current increment
!
! Output: Ke - Elementar stiffness matrix
!         Finte - Elementar internal forces
!------------------------------------------------------------------------------------------------------------

subroutine ElHex27(nel,Ke,el_ddisp,Finte,inc)
    
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
    
    real(8)::yc,zc,xc
    
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
                
                !Local strain-displacement matrix at end configuration
                Bloc=0.0d0
                Bloc=matmul(TGL,Bglob)
                
                if(nlgeom == .false.) Bloc_mid = Bloc
                
                call MatPlastic3D(iprops,nds,(p+1)*(q+1)*(w+1),props,Bloc_mid,6,el_ddisp,stress(nel,cpi,:),strain(nel,cpi,:),dstrain(nel,cpi,:),dstress(nel,cpi,:),hard(nel,cpi),matD)
                
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
                            call prssr(nel,npi_xi,npi_eta,npi_zeta,jac,xi,eta,zeta,we(i),wn(j),wc(k),Finte)
                        end if
                    end do
                end if
                
                !Strain Energy --------------------------------------------------------------
                do k1=1,ntens
                    SEnergy = SEnergy + tempstr(k1,1)*deform(k1,1)*gwt/2.0d0
                end do
                
                ni = INN(IEN(nel,1),1)
                nj = INN(IEN(nel,1),2)
                nk = INN(IEN(nel,1),3)
                
                xc = 0.0d0
                yc = 0.0d0
                zc = 0.0d0
                do k1=1,(p+1)*(q+1)*(w+1)
                    xc = xc + R(k1)*Points(IEN(nel,k1),1)
                    yc = yc + R(k1)*Points(IEN(nel,k1),2)
                    zc = zc + R(k1)*Points(IEN(nel,k1),3)
                end do
                
                continue
                
            end do  
        end do
    end do
    
    continue

end subroutine ElHex27