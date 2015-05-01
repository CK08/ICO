!------------------------------------------------------------------------------------------------------------
!
! General quadrilateral element for plain strain problems with a 3x3 gaussian quadrature
!
! Input: inc - Increment Number
!        nel - Patch Element number
!        nelp - Global Element number
!        el_ddisp - Element displacement
!
! Output: Ke - Elementar stiffness matrix
!         Finte - Elementar internal forces
!------------------------------------------------------------------------------------------------------------

subroutine ElQuad9E(nel,nelp,Ke,el_ddisp,Finte,inc)
    
    use Mod_Variables
    implicit none
    
    integer(4)::nel,nelp,inc
    integer(4)::cpi,ipos
    
    real(8),dimension((p+1)*(q+1)*nds,1),intent(IN)::el_ddisp
    real(8),dimension((p+1)*(q+1)*nds,(p+1)*(q+1)*nds),intent(OUT)::Ke
    real(8),dimension((p+1)*(q+1)*nds,1),intent(OUT)::Finte
    
    integer(4),parameter::npi_xi= 3, npi_eta = 3
    
    integer(4)::i,j,ni,nj,nk,nel_nza,k1,k2,k
    real(8)::xi,eta
    real(8),dimension(npi_xi)::e,we
    real(8),dimension(npi_eta)::n,wn
    
    real(8),dimension((p+1)*(q+1))::R, R_mid
    real(8),dimension((p+1)*(q+1),nds)::dRdx, dRdx_mid
    real(8)::detJ,gwt
    
    real(8),dimension(3,(p+1)*(q+1)*nds)::Bglob,Bloc
    real(8),dimension(3,(p+1)*(q+1)*nds)::Bglob_mid,Bloc_mid
    
    real(8),dimension(3,3)::matD
    real(8),dimension(ntens,1)::tempstr,deform
    
    real(8)::dl,dx,dy
    real(8),dimension(nds,(p+1)*(q+1)*nds)::MatShp
    real(8),dimension(nds,1)::vec
    real(8),dimension((p+1)*(q+1)*nds,1)::Fdist
    
    !NLGeom
    integer(4)::ilixo,count
    real(8)::detj_mid
    real(8),dimension(nds)::temp1
    real(8),dimension((p+1)*(q+1),nds)::nodes
    real(8),dimension((p+1)*(q+1)*nds,1)::tdisp,mdisp,tzero
    real(8),dimension((p+1)*(q+1)*nds,(p+1)*(q+1)*nds)::KNLG
    real(8),dimension(4,4)::MSTR
    real(8),dimension(4,(p+1)*(q+1)*nds)::BNL
    real(8),dimension(2,2)::rconv_mid,rconv
    real(8),dimension(2,2)::jac_mid,jac,jacinv_mid,jacinv
    real(8),dimension(3,3)::matF_mid,matF,MatR,MatR_mid
    real(8),dimension((p+1)*(q+1))::dRde,dRdn
    real(8),dimension(2,2)::Trans,Trans_mid
    real(8),dimension(2,2)::Temp33,Temp33_mid
    real(8),dimension(3,3)::TCL_mid,TCL
    real(8),dimension(3,3)::TGL_mid,TGL
    real(8),dimension(nds,1)::dRden,dRdxy
    
    real(8)::cx,cy
    
    cpi = 0

    we=1.0d0
    wn=1.0d0
    call gauleg(npi_xi, e, we)
    call gauleg(npi_eta, n, wn)
    
    !Displacement for mid- and end-configuration
    if(nlgeom==.true.)then
        tzero = 0.0d0
        tdisp = el_ddisp
        mdisp = el_ddisp/2.0d0
    else
        tzero = 0.0d0
        tdisp = 0.0d0
        mdisp= 0.0d0
        
        TGL = 0.0d0
        do i=1,3
            TGL(i,i) = 1.0d0
        end do
    end if
    
    !Gauss points cycle
    Ke = 0.0d0
    Finte = 0.0d0
    
    KNLG = 0.0d0
    BNL = 0.0d0
    MSTR = 0.0d0
    
!    open(unit=10,file='StressI.txt')
    
    do i=1,npi_xi
        do j=1,npi_eta
            
            cpi = cpi + 1
            !MatD(:,:) = TmatD(nel,cpi,:,:)
            
            xi = e(i)
            eta = n(j)
            
            call ShapeFunc2(nel,xi,eta,R,dRdx,jac,detj,tzero)
            
            !Weight factor
            gwt = we(i)*wn(j)*detj
            
            !Coordinates of the initial configuration
            count = 0
            do k1=1,(p+1)*(q+1)
                count = count + 1
                nodes(k1,1) = Points(IEN(nel,count),1)
                nodes(k1,2) = Points(IEN(nel,count),2)
            end do
            
            if(nlgeom==.true.)then
                
                !Jacobian inverse of the initial configuration
                jacinv = 0.0d0
                jacinv = jac
                temp1 = 0.0d0
                call gaussj(jacinv,nds,temp1,ilixo)
                
                !Local Axis ---------------------------
                if(inc==1)then
                    call local_axis_2D(nds,jac,rconv)
                else
                    rconv = laxis(nelp,cpi,:,:)
                end if
                
                !-------------------------------------------------------------------------
                !Mid-Configuration
                !-------------------------------------------------------------------------
                matF_mid = 0.0d0
                do k1=1,(p+1)*(q+1)
                    dRde(k1) = dRdx(k1,1)
                    dRdn(k1) = dRdx(k1,2)
                end do
                
                !Deformation Gradient for mid-configuration
                call DefGrad2D((p+1)*(q+1),nds,dRde,dRdn,jacinv,nodes,mdisp,MatF_mid)
                
                !Polar Decomposition for mid configuration
                matR_mid=0.0d0
                call PolarDecomp2D(3,matF_mid,matR_mid)
                
                !Update local reference system to the mid configuration
                call axisupdate2D(rconv,matR_mid,rconv_mid)
                
                !Global to natural jacobian in mid configuration
                call ShapeFunc2(nel,xi,eta,R_mid,dRdx_mid,jac_mid,detj_mid,mdisp)
                
                !Jacobian inverse in mid-configuration
                jacinv_mid = 0.0d0
                temp1 = 0.0d0
                jacinv_mid = jac_mid
                temp1 = 0.0d0
                call gaussj(jacinv_mid,nds,temp1,ilixo)
                
                !Transformation matrices mid configuration 
                Temp33_mid=0.0d0
                Temp33_mid= matmul(transpose(rconv_mid),jacinv_mid)
                Trans_mid=0.0d0
                Trans_mid=transpose(rconv_mid)
                
                !Natural to local transformation matrix for mid configuration
                call TransformationMat2D(Temp33_mid,TCL_mid)
                    
                !Global to local transformation matrix for mid configuration
                call TransformationMat2D(Trans_mid,TGL_mid)
                
                Bglob_mid = 0.0d0
                do k1=1,(p+1)*(q+1)
                    Bglob_mid(1,k1*2-1) = dRdx_mid(k1,1)
                    Bglob_mid(2,k1*2  ) = dRdx_mid(k1,2)
                    Bglob_mid(3,k1*2-1) = dRdx_mid(k1,2)
                    Bglob_mid(3,k1*2  ) = dRdx_mid(k1,1)
                end do
                
                Bloc_mid=0.0d0
                Bloc_mid=matmul(TGL_mid,Bglob_mid)
                
                !-------------------------------------------------------------------------
                !End-Configuration
                !-------------------------------------------------------------------------
                
                !Deformation Gradient for end-configuration
                call DefGrad2D((p+1)*(q+1),nds,dRde,dRdn,jacinv,nodes,tdisp,MatF)
                
                !Polar Decomposition for mid configuration
                matR_mid=0.0d0
                call PolarDecomp2D(3,matF,matR)
                
                !Update local reference system to the mid configuration
                call axisupdate2D(rconv,matR,rconv)
                
                !Global to natural jacobian in mid configuration
                call ShapeFunc2(nel,xi,eta,R,dRdx,jac,detj,tdisp)
                
                !Jacobian inverse in mid-configuration
                jacinv = 0.0d0
                temp1 = 0.0d0
                jacinv = jac
                call gaussj(jacinv,nds,temp1,ilixo)
                
                !Transformation matrices mid configuration 
                Temp33=0.0d0
                Temp33= matmul(transpose(rconv),jacinv)
                Trans=0.0d0
                Trans=transpose(rconv)
                
                !Natural to local transformation matrix for mid configuration
                call TransformationMat2D(Temp33,TCL)
                    
                !Global to local transformation matrix for mid configuration
                call TransformationMat2D(Trans,TGL)
                
                !Update Weight factor
                gwt = we(i)*wn(j)*detj
                
                continue
                
            end if
            
            Bglob = 0.0d0
            do k1=1,(p+1)*(q+1)
                Bglob(1,k1*2-1) = dRdx(k1,1)
                Bglob(2,k1*2  ) = dRdx(k1,2)
                Bglob(3,k1*2-1) = dRdx(k1,2)
                Bglob(3,k1*2  ) = dRdx(k1,1)
            end do
            
            !Local strain-displacement matrix at end configuration
            Bloc=0.0d0
            Bloc=matmul(TGL,Bglob)
            
            if(nlgeom == .false.) Bloc_mid = Bloc
            
            call MaterialPStrain(iprops,nds,(p+1)*(q+1),props,Bloc_mid,3,el_ddisp,stress(nelp,cpi,:),strain(nelp,cpi,:),dstrain(nelp,cpi,:),dstress(nelp,cpi,:),hard(nelp,cpi),matD)
            
            tempstr = 0.0d0
            do k1=1,ntens
                tempstr(k1,1) = stress(nelp,cpi,k1)
                deform(k1,1) =  dstrain(nelp,cpi,k1)
            end do
            
            if(NLGeom == .true.)then
                
                BNL = 0.0d0
                do k1=1,(p+1)*(q+1)
                    dRden = 0.0d0
                    dRden(1,1) = dRdx(k1,1)
                    dRden(2,1) = dRdx(k1,2)
                    
                    dRdxy = 0.0d0
                    dRdxy = matmul(Trans,dRden)

                    BNL(1,k1*2-1) = dRdxy(1,1)
                    BNL(2,k1*2-1) = dRdxy(2,1)
                    BNL(3,k1*2  ) = dRdxy(1,1)
                    BNL(4,k1*2  ) = dRdxy(2,1)

                end do 
             
                MSTR=0.0d0
                do k1=1,nds
                    MSTR(k1*2-1,k1*2-1) = tempstr(1,1)
                    MSTR(k1*2-1,k1*2  ) = tempstr(3,1)
                    MSTR(k1*2  ,k1*2-1) = tempstr(3,1)
                    MSTR(k1*2  ,k1*2  ) = tempstr(2,1)
                end do
                
                KNLG = KNLG + matmul(matmul(transpose(BNL),MSTR),BNL)*gwt
                
                Ke = Ke + matmul(matmul(transpose(Bloc),matD),Bloc)*gwt + matmul(matmul(transpose(BNL),MSTR),BNL)*gwt
                    
                !Export local axis ---------------
                laxis(nelp,cpi,:,:) = rconv
                
            else
                Ke = Ke + matmul(matmul(transpose(Bloc),matD),Bloc)*gwt
            end if
            
            !Ke = Ke + matmul(matmul(transpose(Bloc),matD),Bloc)*gwt
            
            Finte = Finte + matmul(transpose(Bloc),tempstr)*gwt
            
            do k1=1,ntens
                SEnergy = SEnergy + tempstr(k1,1)*deform(k1,1)*gwt/2.0d0
            end do
            

            
            ipos = npi*(nelp-1)+cpi
            do k1=1,(p+1)*(q+1)
                GP_coords(ipos,1) = GP_coords(ipos,1) + R(k1)*(nodes(k1,1)+el_ddisp(k1*nds-1,1))
                GP_coords(ipos,2) = GP_coords(ipos,2) + R(k1)*(nodes(k1,2)+el_ddisp(k1*nds  ,1))
            end do
            
            continue
          
        end do
    end do
    
!    close(10)
    
    continue

end subroutine ElQuad9E