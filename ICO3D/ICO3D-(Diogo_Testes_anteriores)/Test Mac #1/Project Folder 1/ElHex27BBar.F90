!------------------------------------------------------------------------------------------------------------
!
! Hexaedral element with a 3x3x3 gaussian quadrature and B-bar method, Q2\Q1 from the work of Elguedj et.al.
! "B-bar and F-bar projection methods for nearly incompressible linear and nonlinear elsaticity and 
! plasticity using higher-order NURBS elements" (2008)
!
! Input: nel - Element number
!        el_ddisp - Element displacement
!        
! Output: Ke - Elementar stiffness matrix (Not needed)
!         Finte - Elementar internal forces (Still not working)
!         MatA - Stiffness for disp-disp component
!         MatC - Stiffness for disp-proj component
!         MatV - Stiffness for proj-proj component
!
! ATTENTION: NOT FULLY WORKING. Only possible to perform linear elastic analysis
!
!------------------------------------------------------------------------------------------------------------
subroutine ElHex27BBar(nel,Ke,el_ddisp,Finte,MatA,MatC,MatM)
    
    use Mod_Variables
    implicit none
    
    integer(4)::nel
    integer(4)::cpi
    
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::el_ddisp
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(OUT)::Ke
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(OUT)::Finte
    
    integer(4),parameter::npi_xi= 3, npi_eta = 3, npi_zeta = 3
    
    integer(4)::i,j,jj,k,k2,k3,k4,k5,k6,k7,k8,ni,nj,nk,nel_nza,k1,error,count,count2
    real(8)::xi,eta,zeta
    real(8),dimension(npi_xi)::e,we
    real(8),dimension(npi_eta)::n,wn
    real(8),dimension(npi_zeta)::c,wc
    
    real(8),dimension((p+1)*(q+1)*(w+1))::R
    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdx
    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob,Bloc
    
    real(8)::detj,detjnew,gwt
    
    real(8),dimension((p)*(q)*(w))::Rbar
    real(8),dimension((p)*(q)*(w),nds)::dRdxbar
    real(8),dimension(6,6)::matD,DDev,DDil
    real(8),dimension(ntens,1)::tempstr,deform
    real(8),dimension(p*q*w,p*q*w)::matM,Minv
    real(8),dimension(p*q*w,1)::aux1
    real(8),dimension(p*q*w)::temp1
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds)::MatA
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,p*q*w)::MatC
    real(8),dimension(6,1)::mm
    real(8)::dil,eg,lambda
    real(8),dimension(npi_xi-1)::ebar,webar
    real(8),dimension(npi_eta-1)::nbar,wnbar
    real(8),dimension(npi_zeta-1)::cbar,wcbar
    real(8),dimension(ncpx+p-1)::ubar
    real(8),dimension(ncpy+q-1)::vbar
    real(8),dimension(ncpz+w-1)::wbar
    
    cpi = 0
    
    we = 1.0d0
    wn = 1.0d0
    wc = 1.0d0
    call gauleg(npi_xi, e, we)
    call gauleg(npi_eta, n, wn)
    call gauleg(npi_zeta, c, wc)
    
    call gauleg(npi_xi-1, ebar, webar)
    call gauleg(npi_eta-1, nbar, wnbar)
    call gauleg(npi_zeta-1, cbar, wcbar)
    
    Ke = 0.0d0
    MatA = 0.0d0
    MatM = 0.0d0
    MatC = 0.0d0

    !Gauss points cycle -------------------------------------------------------------------------------
    cpi = 0
    Ke = 0.0d0
    Finte = 0.0d0
    do i=1,npi_xi
        do j=1,npi_eta
            do k=1,npi_zeta
            
                cpi = cpi + 1
                MatD(:,:) = TmatD(nel,cpi,:,:)
                
                xi = e(i)
                eta = n(j)
                zeta = c(k)
                
                call ShapeFunc(nel,xi,eta,zeta,R,dRdx,detj)
                
                !Weight factor
                gwt = we(i)*wn(j)*wc(k)*detj
                
                !Strain-displacement matrix -------------
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
                
                !Projected space knot vectors --------------------
                do k1=2,ncpx+p
                    ubar(k1-1) = u_knot(k1)
                end do
                
                do k1=2,ncpy+q
                    vbar(k1-1) = v_knot(k1)
                end do 
                
                do k1=2,ncpz+w
                    wbar(k1-1) = w_knot(k1)
                end do     
                
                !Projected Space basis functions -----------------
                call ShapeFuncBar(nel,xi,eta,zeta,INN,IEN,B_net,ubar,vbar,wbar,p-1,q-1,w-1,ncpx,ncpy,ncpz,nds,nnodes,nelems,RBar,dRdxBar,detjnew)
                
                Bloc = Bglob
                call MatPlastic3D(iprops,nds,(p+1)*(q+1)*(w+1),props,Bloc,6,el_ddisp,stress(nel,cpi,:),strain(nel,cpi,:),dstrain(nel,cpi,:),dstress(nel,cpi,:),hard(nel,cpi),matD)
                
                tempstr = 0.0d0
                do k1=1,ntens
                    tempstr(k1,1) = stress(nel,cpi,k1)
                    deform(k1,1) = dstrain(nel,cpi,k1)
                end do
                
                !Deviatoric Component of the constitutive matrix --------------
                DDev = 0.0d0
                do k1=1,3
                    DDev(k1,k1) = 1.0d0
                    DDev(k1+3,k1+3) = 1.0d0/2.0d0
                end do
                
                do k1=1,3
                    do k2=1,3
                        DDev(k1,k2) = DDev(k1,k2) - 1.0d0/3.0d0
                    end do
                end do
                
                Ddev = Ddev*2.0d0*props(1)/(2.0d0*(1.0d0+props(2)))
                
                !Dilatational Component of the constitutive matrix --------------
                Ddil = MatD - DDev
                eg = props(1)/(2.0d0*(1.0d0+props(2)))
                lambda = props(1)*props(2)/((1.0d0+props(2))*(1.0d0-2.0d0*props(2)));
                dil = 1.0d0/3.0d0*(3.0d0*lambda+2.0d0*eg)
                

                !Matrix A ------------------------------------------------------------
                MatA = MatA + matmul(matmul(transpose(Bglob),DDev),Bglob)*gwt
                !---------------------------------------------------------------------
                
                !Matrix V ------------------------------------------------------------
                do k1=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1)      !index A
                    do k2=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1)  !index B
                        MatM(k1,k2) = MatM(k1,k2) + rbar(k1)*rbar(k2)*gwt*dil
                    end do
                end do
                !---------------------------------------------------------------------
                
                !Matrix C ------------------------------------------------------------
                mm = 0.0d0
                do k1=1,3
                    mm(k1,1) = 1.0d0
                end do
                
                aux1 = 0.0d0
                do k1=1,p*q*w
                    aux1(k1,1) = rbar(k1)
                end do
                
                MatC = MatC + matmul(matmul(transpose(BGlob),mm),transpose(aux1))*gwt*dil
                !---------------------------------------------------------------------
                

                Ke = Ke + matmul(matmul(transpose(Bglob),matD),Bglob)*gwt
                
                Finte = Finte + matmul(transpose(Bglob),tempstr)*gwt
                
                do k1=1,ntens
                    SEnergy = SEnergy + tempstr(k1,1)*deform(k1,1)*gwt/2.0d0
                end do
                
                continue
                
            end do  
        end do
    end do

    continue

end subroutine ElHex27BBar






!subroutine ElHex27BBar(nel,Ke,el_ddisp,Finte)
!    
!    use Mod_Variables
!    implicit none
!    
!    integer(4)::nel
!    integer(4)::cpi
!    
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::el_ddisp
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(OUT)::Ke
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(OUT)::Finte
!    
!    integer(4),parameter::npi_xi= 3, npi_eta = 3, npi_zeta = 3
!    
!    integer(4)::i,j,jj,k,k2,k3,k4,k5,k6,k7,k8,ni,nj,nk,nel_nza,k1,error,count,count2
!    real(8)::xi,eta,zeta
!    real(8),dimension(npi_xi)::e,we
!    real(8),dimension(npi_eta)::n,wn
!    real(8),dimension(npi_zeta)::c,wc
!    
!    real(8),dimension(ncpx+p-1)::ubar
!    real(8),dimension(ncpy+q-1)::vbar
!    real(8),dimension(ncpz+w-1)::wbar
!    
!    real(8),dimension(npi_xi-1)::ebar,webar,eext
!    real(8),dimension(npi_eta-1)::nbar,wnbar,next
!    real(8),dimension(npi_zeta-1)::cbar,wcbar,cext
!    
!    real(8),dimension((p+1)*(q+1)*(w+1))::R
!    real(8),dimension((p)*(q)*(w))::Rbar
!    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdx
!    real(8),dimension((p)*(q)*(w),nds)::dRdxbar
!    real(8),dimension(8,3)::coef
!
!    real(8)::detJ,gwt,detjnew,lixo1,lixo2,lixo3,volnew
!    
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob,Bloc
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bdil,Bdev,BBar
!    
!    real(8),dimension(6,6)::matD,DDev,DDil
!    real(8),dimension(ntens,1)::tempstr,deform
!    
!    real(8),dimension(8,8)::matM,Minv
!    real(8),dimension(8,1)::aux1
!    real(8),dimension(8)::temp1,lixo
!    real(8),dimension(1,1)::temp11
!    real(8),dimension(27)::deter
!    real(8),dimension(3)::temp3
!    real(8)::NB
!    real(8),dimension(3,(p+1)*(q+1)*(w+1))::Bteste
!    
!    real(8),dimension(8,(p+1)*(q+1)*(w+1)*nds)::SecInt
!    
!    real(8),dimension(2,2,2,4)::NewCoord
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds)::KDil,Kdev
!    
!!    real(8),dimension((p+1)*(q+1)*(w+1)*nds)::SecInt
!    
!    real(8)::vol
!    
!    cpi = 0
!    
!    vol = 0.0d0
!    
!    we = 1.0d0
!    wn = 1.0d0
!    wc = 1.0d0
!    call gauleg(npi_xi, e, we)
!    call gauleg(npi_eta, n, wn)
!    call gauleg(npi_zeta, c, wc)
!    
!    call gauleg(npi_xi-1, ebar, webar)
!    call gauleg(npi_eta-1, nbar, wnbar)
!    call gauleg(npi_zeta-1, cbar, wcbar)
!    
!    eext(1) = -1.0d0
!    eext(2) =  1.0d0
!    
!    next(1) = -1.0d0
!    next(2) =  1.0d0
!    
!    cext(1) = -1.0d0
!    cext(2) =  1.0d0
!
!    !Calculate volume of the element
!    Ke = 0.0d0
!    KDil = 0.0d0
!    Kdev = 0.0d0
!    
!    Finte = 0.0d0
!    count = 0
!    do i=1,npi_xi
!        do j=1,npi_eta
!            do k=1,npi_zeta
!                
!                count = count + 1
!                
!                xi = e(i)
!                eta = n(j)
!                zeta = c(k)
!                
!                call ShapeFunc(nel,xi,eta,zeta,R,dRdx,detj)
!                
!                vol = vol + detj*we(i)*wn(j)*wc(k)
!                deter(count) = detj
!                
!                continue
!                
!            end do  
!        end do
!    end do
!    
!    
!    
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    !Projected coordinates
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!!    count2 = 0
!!    do i=1,npi_xi-1
!!        do j=1,npi_eta-1
!!            do k=1,npi_zeta-1
!!                
!!                !count = count + 1
!!                
!!                xi = eext(i)
!!                eta = next(j)
!!                zeta = cext(k)
!!                
!!                ni = INN(IEN(nel,1),1)
!!                nj = INN(IEN(nel,1),2)
!!                nk = INN(IEN(nel,1),3)
!!                
!!                call ShapeFunc(nel,xi,eta,zeta,R,dRdx,detj)
!!                
!!                count2 = count2 + 1
!!                
!!                count = 0
!!                do k1=0,w
!!                    do k2=0,q
!!                        do k3=0,p
!!                        
!!                            count = count + 1
!!                        
!!                            NewCoord(i,j,k,:) = NewCoord(i,j,k,:) + R(count)*B_net(ni-k3,nj-k2,nk-k1,:)
!!                        end do
!!                    end do
!!                end do
!!                
!!                
!!                continue
!!                
!!            end do  
!!        end do
!!    end do
!    
!    
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    ! M-MATRIX
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    matM = 0.0d0
!    NB = 0.0d0
!    do i=1,npi_xi
!        do j=1,npi_eta
!            do k=1,npi_zeta
!                
!                cpi = cpi + 1
!                
!                xi = e(i)
!                eta = n(j)
!                zeta = c(k)
!                
!                call ShapeFunc(nel,xi,eta,zeta,R,dRdx,detj)
!                
!                !#####################################################################################################################################
!                !#####################################################################################################################################
!                do k1=2,ncpx+p
!                    ubar(k1-1) = u_knot(k1)
!                end do
!                
!                do k1=2,ncpy+q
!                    vbar(k1-1) = v_knot(k1)
!                end do 
!                
!                do k1=2,ncpz+w
!                    wbar(k1-1) = w_knot(k1)
!                end do     
!                
!                call ShapeFuncBar(nel,xi,eta,zeta,INN,IEN,B_net,ubar,vbar,wbar,p-1,q-1,w-1,ncpx,ncpy,ncpz,nds,nnodes,nelems,RBar,dRdxBar,detjnew)
!                
!                aux1 = 0.0d0
!                do k1=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1)
!                    aux1(k1,1) = Rbar(k1)
!                end do
!                !#####################################################################################################################################
!                !#####################################################################################################################################
!                
!                do k1=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1)      !index A
!                    do k2=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1)  !index B
!                        
!                        count = 0
!                    
!                        do k3=1,npi_xi
!                            do k4=1,npi_eta
!                                do k5=1,npi_zeta
!                                    count = count + 1
!                                    MatM(k1,k2) = matM(k1,k2) + Rbar(k1)*Rbar(k2)*we(k3)*wn(k4)*wc(k5)*deter(count)
!                                end do
!                            end do
!                        end do
!                    enddo
!                end do
!                
!!                count = 0
!!                do k3=1,npi_xi
!!                    do k4=1,npi_eta
!!                        do k5=1,npi_zeta
!!                            count = count + 1
!!                            matM = matM + matmul(aux1,transpose(aux1))*we(k3)*wn(k4)*wc(k5)*deter(count)
!!                        end do
!!                    end do
!!                end do
!                
!                !matM = matM + matmul(aux1,transpose(aux1))*we(i)*wn(j)*wc(k)*detj
!                
!                !NB = NB + sum(Rbar)*we(i)*wn(j)*wc(k)*detj
!                
!            end do
!        end do
!    end do
!    
!    
!!    lixo=0.0d0
!!    do k1=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1)
!!        lixo(k1) = sum(matM(:,k1))
!!    end do
!!    
!!    MatM = 0.0d0
!!    do k1=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1)
!!        matM(k1,k1) = lixo(k1)
!!    end do
!    
!    Minv = 0.0d0
!    Minv = MatM
!    call gaussj(Minv,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1),temp1,error)
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!  
!  
!  
!  
!  
!  
!
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    ! SECOND INTEGRAL
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    cpi = 0
!    Ke = 0.0d0
!    Finte = 0.0d0
!    SecInt = 0.0d0
!    do i=1,npi_xi
!        do j=1,npi_eta
!            do k=1,npi_zeta
!            
!                cpi = cpi + 1
!                
!                xi = e(i)
!                eta = n(j)
!                zeta = c(k)
!                
!                call ShapeFunc(nel,xi,eta,zeta,R,dRdx,detj)
!                
!                !Weight factor
!                gwt = we(i)*wn(j)*wc(k)*detj
!                
!                Bglob = 0.0d0
!                Bdil = 0.0d0
!                Bdev = 0.0d0
!                BBar = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
!                    Bglob(1,k1*3-2) = dRdx(k1,1)
!                    Bglob(2,k1*3-1) = dRdx(k1,2)
!                    Bglob(3,k1*3  ) = dRdx(k1,3)
!                    Bglob(4,k1*3-2) = dRdx(k1,2)
!                    Bglob(4,k1*3-1) = dRdx(k1,1)
!                    Bglob(5,k1*3-2) = dRdx(k1,3)
!                    Bglob(5,k1*3  ) = dRdx(k1,1)
!                    Bglob(6,k1*3-1) = dRdx(k1,3)
!                    Bglob(6,k1*3  ) = dRdx(k1,2)
!                end do
!                
!                
!                !#####################################################################################################################################
!                !#####################################################################################################################################
!                do k1=2,ncpx+p
!                    ubar(k1-1) = u_knot(k1)
!                end do
!                
!                do k1=2,ncpy+q
!                    vbar(k1-1) = v_knot(k1)
!                end do 
!                
!                do k1=2,ncpz+w
!                    wbar(k1-1) = w_knot(k1)
!                end do     
!                
!                call ShapeFuncBar(nel,xi,eta,zeta,INN,IEN,B_net,ubar,vbar,wbar,p-1,q-1,w-1,ncpx,ncpy,ncpz,nds,nnodes,nelems,RBar,dRdxBar,detjnew)
!                
!                aux1 = 0.0d0
!                do k1=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1)
!                    aux1(k1,1) = Rbar(k1)
!                end do
!                !#####################################################################################################################################
!                !#####################################################################################################################################
!                
!                
!                do k1=1,(p+1)*(q+1)*(w+1) !index i
!                    do k2=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1) !index B
!                        count = 0
!                        do k3=1,npi_xi
!                            do k4=1,npi_eta
!                                do k5=1,npi_zeta
!                                    count = count + 1
!!                                    SecInt(k1*3-2) = SecInt(k1*3-2) + Rbar(k2)*Bglob(1,k1*3-2)*we(k3)*wn(k4)*wc(k5)*deter(count)
!!                                    SecInt(k1*3-1) = SecInt(k1*3-1) + Rbar(k2)*Bglob(2,k1*3-1)*we(k3)*wn(k4)*wc(k5)*deter(count)
!!                                    SecInt(k1*3  ) = SecInt(k1*3  ) + Rbar(k2)*Bglob(3,k1*3  )*we(k3)*wn(k4)*wc(k5)*deter(count)
!                                    SecInt(k2,k1*3-2) = SecInt(k2,k1*3-2) + Rbar(k2)*Bglob(1,k1*3-2)*we(k3)*wn(k4)*wc(k5)*deter(count)
!                                    SecInt(k2,k1*3-1) = SecInt(k2,k1*3-1) + Rbar(k2)*Bglob(2,k1*3-1)*we(k3)*wn(k4)*wc(k5)*deter(count)
!                                    SecInt(k2,k1*3  ) = SecInt(k2,k1*3  ) + Rbar(k2)*Bglob(3,k1*3  )*we(k3)*wn(k4)*wc(k5)*deter(count)
!                                end do
!                            end do
!                        end do
!                    end do
!                end do
!                
!                continue
!
!            end do  
!        end do
!    end do
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    
!    
!    
!    
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    !Gauss points cycle
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    cpi = 0
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
!                
!                !Weight factor
!                gwt = we(i)*wn(j)*wc(k)*detj
!                
!                Bglob = 0.0d0
!                Bdil = 0.0d0
!                Bdev = 0.0d0
!                BBar = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
!                    Bglob(1,k1*3-2) = dRdx(k1,1)
!                    Bglob(2,k1*3-1) = dRdx(k1,2)
!                    Bglob(3,k1*3  ) = dRdx(k1,3)
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
!                end do
!                
!                
!                
!                
!                !#####################################################################################################################################
!                !#####################################################################################################################################
!                do k1=2,ncpx+p
!                    ubar(k1-1) = u_knot(k1)
!                end do
!                
!                do k1=2,ncpy+q
!                    vbar(k1-1) = v_knot(k1)
!                end do 
!                
!                do k1=2,ncpz+w
!                    wbar(k1-1) = w_knot(k1)
!                end do     
!                
!                call ShapeFuncBar(nel,xi,eta,zeta,INN,IEN,B_net,ubar,vbar,wbar,p-1,q-1,w-1,ncpx,ncpy,ncpz,nds,nnodes,nelems,RBar,dRdxBar,detjnew)
!                
!                aux1 = 0.0d0
!                do k1=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1)
!                    aux1(k1,1) = Rbar(k1)
!                end do
!                !#####################################################################################################################################
!                !#####################################################################################################################################
!                
!                
!
!                
!!                SecInt = 0.0d0
!!                do k1=1,(p+1)*(q+1)*(w+1)
!!                    do k2=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1)
!!                        count = 0
!!                        do k3=1,npi_xi; do k4=1,npi_eta; do k5=1,npi_zeta
!!                            count = count + 1
!!                            SecInt(k2,k1*3-2) = SecInt(k2,k1*3-2) + Rbar(k2)*Bglob(1,k1*3-2)*we(k3)*wn(k4)*wc(k5)*deter(count)
!!                            SecInt(k2,k1*3-1) = SecInt(k2,k1*3-1) + Rbar(k2)*Bglob(2,k1*3-1)*we(k3)*wn(k4)*wc(k5)*deter(count)
!!                            SecInt(k2,k1*3  ) = SecInt(k2,k1*3  ) + Rbar(k2)*Bglob(3,k1*3  )*we(k3)*wn(k4)*wc(k5)*deter(count)
!!                        end do; end do; end do
!!                    end do
!!                end do
!                
!                
!                !matM = matmul(aux1,transpose(aux1))*we(i)*wn(j)*wc(k)*detj
!                
!
!                
!                BBar = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
!                    count = 0
!                    do k2=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1) !index A
!                        do k3=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1) !index B
!                            
!!                            BBar(1,k1*3-2) = BBar(1,k1*3-2) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3-2)
!!                            BBar(2,k1*3-2) = BBar(2,k1*3-2) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3-2)
!!                            BBar(3,k1*3-2) = BBar(3,k1*3-2) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3-2)
!!                            
!!                            BBar(1,k1*3-1) = BBar(1,k1*3-1) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3-1)
!!                            BBar(2,k1*3-1) = BBar(2,k1*3-1) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3-1)
!!                            BBar(3,k1*3-1) = BBar(3,k1*3-1) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3-1)
!!                            
!!                            BBar(1,k1*3  ) = BBar(1,k1*3  ) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3  )
!!                            BBar(2,k1*3  ) = BBar(2,k1*3  ) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3  )
!!                            BBar(3,k1*3  ) = BBar(3,k1*3  ) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3  )
!                            
!!                            BBar(1,k1*3-2) = BBar(1,k1*3-2) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3-2)
!!                            BBar(2,k1*3-2) = BBar(2,k1*3-2) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3-2)
!!                            BBar(3,k1*3-2) = BBar(3,k1*3-2) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3-2)
!!                            
!!                            BBar(1,k1*3-1) = BBar(1,k1*3-1) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3-1)
!!                            BBar(2,k1*3-1) = BBar(2,k1*3-1) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3-1)
!!                            BBar(3,k1*3-1) = BBar(3,k1*3-1) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3-1)
!!                            
!!                            BBar(1,k1*3  ) = BBar(1,k1*3  ) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3  )
!!                            BBar(2,k1*3  ) = BBar(2,k1*3  ) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3  )
!!                            BBar(3,k1*3  ) = BBar(3,k1*3  ) + Rbar(k2)*MInv(k2,k3)*SecInt(k1*3  )
!                            
!                            BBar(1,k1*3-2) = BBar(1,k1*3-2) + Rbar(k2)*MInv(k2,k3)*SecInt(k3,k1*3-2)
!                            BBar(2,k1*3-2) = BBar(2,k1*3-2) + Rbar(k2)*MInv(k2,k3)*SecInt(k3,k1*3-2)
!                            BBar(3,k1*3-2) = BBar(3,k1*3-2) + Rbar(k2)*MInv(k2,k3)*SecInt(k3,k1*3-2)
!                            BBar(1,k1*3-1) = BBar(1,k1*3-1) + Rbar(k2)*MInv(k2,k3)*SecInt(k3,k1*3-1)
!                            BBar(2,k1*3-1) = BBar(2,k1*3-1) + Rbar(k2)*MInv(k2,k3)*SecInt(k3,k1*3-1)
!                            BBar(3,k1*3-1) = BBar(3,k1*3-1) + Rbar(k2)*MInv(k2,k3)*SecInt(k3,k1*3-1)
!                            BBar(1,k1*3  ) = BBar(1,k1*3  ) + Rbar(k2)*MInv(k2,k3)*SecInt(k3,k1*3  )
!                            BBar(2,k1*3  ) = BBar(2,k1*3  ) + Rbar(k2)*MInv(k2,k3)*SecInt(k3,k1*3  )
!                            BBar(3,k1*3  ) = BBar(3,k1*3  ) + Rbar(k2)*MInv(k2,k3)*SecInt(k3,k1*3  )
!                            
!                        end do
!                    end do
!                end do
!                
!                
!                BBar  = Bbar/3.0d0
!                Bdev  = Bglob - Bdil
!                Bglob = BBar + Bdev
!                
!                continue
!                
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
!                
!                DDev = 0.0d0
!                do k1=1,3
!                    DDev(k1,k1) = 1.0d0
!                    DDev(k1+3,k1+3) = 1.0d0/2.0d0
!                end do
!                
!                do k1=1,3
!                    do k2=1,3
!                        DDev(k1,k2) = DDev(k1,k2) - 1.0d0/3.0d0
!                    end do
!                end do
!                
!                Ddev = Ddev*props(1)/(1.0d0+props(2))
!                
!                Ddil = MatD - DDev
!                
!                KDil = KDil + matmul(matmul(transpose(BBar),DDil),BBar)*gwt
!                KDev = KDev + matmul(matmul(transpose(BDev),DDev),BDev)*gwt
!                
!                !KDev = KDev + KDil
!                
!                Ke = Ke + matmul(matmul(transpose(Bglob),matD),Bglob)*gwt
!                
!                Finte = Finte + matmul(transpose(Bglob),tempstr)*gwt
!                
!!                do k1=1,ntens
!!                    SEnergy = SEnergy + tempstr(k1,1)*deform(k1,1)*gwt/2.0d0
!!                end do
!                
!                continue
!                
!            end do  
!        end do
!    end do
!      
!    !Ke = Kdev
!
!    continue
!
!end subroutine ElHex27BBar



                           
                            
                            
                            
!                            count = 0
!                            
!                            do k4=1,3; do k5=1,3; do k6=1,3
!                            count = count + 1
!                            BBar(1,k1*3-2) = BBar(1,k1*3-2) + Bglob(1,k1*3-2)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k4)*wn(k5)*wc(k6)*deter(count)
!                            BBar(2,k1*3-2) = BBar(2,k1*3-2) + Bglob(1,k1*3-2)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k4)*wn(k5)*wc(k6)*deter(count)
!                            BBar(3,k1*3-2) = BBar(3,k1*3-2) + Bglob(1,k1*3-2)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k4)*wn(k5)*wc(k6)*deter(count)
!                            
!                            BBar(1,k1*3-1) = BBar(1,k1*3-1) + Bglob(2,k1*3-1)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k4)*wn(k5)*wc(k6)*deter(count)
!                            BBar(2,k1*3-1) = BBar(2,k1*3-1) + Bglob(2,k1*3-1)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k4)*wn(k5)*wc(k6)*deter(count)
!                            BBar(3,k1*3-1) = BBar(3,k1*3-1) + Bglob(2,k1*3-1)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k4)*wn(k5)*wc(k6)*deter(count)
!                            
!                            BBar(1,k1*3  ) = BBar(1,k1*3  ) + Bglob(3,k1*3  )*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k4)*wn(k5)*wc(k6)*deter(count)
!                            BBar(2,k1*3  ) = BBar(2,k1*3  ) + Bglob(3,k1*3  )*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k4)*wn(k5)*wc(k6)*deter(count)
!                            BBar(3,k1*3  ) = BBar(3,k1*3  ) + Bglob(3,k1*3  )*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k4)*wn(k5)*wc(k6)*deter(count)
!                                
!                            end do;end do; end do    
!                                BBar(1,k1*3-2) = BBar(1,k1*3-2) + aux1(k2,1)*Bteste(1,k1)*Bglob(1,k1*3-2)*MInv(k2,k3)
!                                BBar(2,k1*3-2) = BBar(2,k1*3-2) + aux1(k2,1)*Bteste(1,k1)*Bglob(1,k1*3-2)*MInv(k2,k3)
!                                BBar(3,k1*3-2) = BBar(3,k1*3-2) + aux1(k2,1)*Bteste(1,k1)*Bglob(1,k1*3-2)*MInv(k2,k3)
!                                BBar(1,k1*3-1) = BBar(1,k1*3-1) + aux1(k2,1)*Bteste(2,k1)*Bglob(2,k1*3-1)*MInv(k2,k3)
!                                BBar(2,k1*3-1) = BBar(2,k1*3-1) + aux1(k2,1)*Bteste(2,k1)*Bglob(2,k1*3-1)*MInv(k2,k3)
!                                BBar(3,k1*3-1) = BBar(3,k1*3-1) + aux1(k2,1)*Bteste(2,k1)*Bglob(2,k1*3-1)*MInv(k2,k3)
!                                BBar(1,k1*3  ) = BBar(1,k1*3  ) + aux1(k2,1)*Bteste(3,k1)*Bglob(3,k1*3  )*MInv(k2,k3)
!                                BBar(2,k1*3  ) = BBar(2,k1*3  ) + aux1(k2,1)*Bteste(3,k1)*Bglob(3,k1*3  )*MInv(k2,k3)
!                                BBar(3,k1*3  ) = BBar(3,k1*3  ) + aux1(k2,1)*Bteste(3,k1)*Bglob(3,k1*3  )*MInv(k2,k3)
                                
!                                count = 0
!                                do k5=1,npi_xi
!                                    do k6=1,npi_eta
!                                        do k7=1,npi_zeta
!                                            count = count + 1
!                                            BBar(1,k1*3-2) = BBar(1,k1*3-2) + Bglob(1,k1*3-2)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!                                            BBar(2,k1*3-2) = BBar(2,k1*3-2) + Bglob(1,k1*3-2)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!                                            BBar(3,k1*3-2) = BBar(3,k1*3-2) + Bglob(1,k1*3-2)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!                                            
!                                            BBar(1,k1*3-1) = BBar(1,k1*3-1) + Bglob(2,k1*3-1)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!                                            BBar(2,k1*3-1) = BBar(2,k1*3-1) + Bglob(2,k1*3-1)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!                                            BBar(3,k1*3-1) = BBar(3,k1*3-1) + Bglob(2,k1*3-1)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!                                            
!                                            BBar(1,k1*3  ) = BBar(1,k1*3  ) + Bglob(3,k1*3  )*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!                                            BBar(2,k1*3  ) = BBar(2,k1*3  ) + Bglob(3,k1*3  )*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!                                            BBar(3,k1*3  ) = BBar(3,k1*3  ) + Bglob(3,k1*3  )*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!                                        end do
!                                    end do
!                                end do




!subroutine ElHex27BBar(nel,Ke,el_ddisp,Finte)
!    
!    use Mod_Variables
!    implicit none
!    
!    integer(4)::nel
!    integer(4)::cpi
!    
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::el_ddisp
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(OUT)::Ke
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(OUT)::Finte
!    
!    integer(4),parameter::npi_xi= 3, npi_eta = 3, npi_zeta = 3
!    
!    integer(4)::i,j,jj,k,k2,k3,k4,k5,k6,k7,k8,ni,nj,nk,nel_nza,k1,error,count
!    real(8)::xi,eta,zeta
!    real(8),dimension(npi_xi)::e,we
!    real(8),dimension(npi_eta)::n,wn
!    real(8),dimension(npi_zeta)::c,wc
!    
!    real(8),dimension(ncpx+p-1)::ubar
!    real(8),dimension(ncpy+q-1)::vbar
!    real(8),dimension(ncpz+w-1)::wbar
!    
!    real(8),dimension(npi_xi-1)::ebar,webar
!    real(8),dimension(npi_eta-1)::nbar,wnbar
!    real(8),dimension(npi_zeta-1)::cbar,wcbar
!    
!    real(8),dimension((p+1)*(q+1)*(w+1))::R
!    real(8),dimension((p)*(q)*(w))::Rbar
!    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdx
!    real(8),dimension((p)*(q)*(w),nds)::dRdxbar
!    real(8),dimension(8,3)::coef
!
!    real(8)::detJ,gwt,detjnew,lixo1,lixo2,lixo3,volnew
!    
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob,Bloc
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bdil,Bdev,BBar
!    
!    real(8),dimension(6,6)::matD
!    real(8),dimension(ntens,1)::tempstr,deform
!    
!    real(8),dimension(8,8)::matM,Minv
!    real(8),dimension(8,1)::aux1
!    real(8),dimension(8)::temp1,lixo
!    real(8),dimension(1,1)::temp11
!    real(8),dimension(27)::deter
!    
!    real(8)::vol
!    
!    cpi = 0
!    
!    vol = 0.0d0
!    
!    we = 1.0d0
!    wn = 1.0d0
!    wc = 1.0d0
!    call gauleg(npi_xi, e, we)
!    call gauleg(npi_eta, n, wn)
!    call gauleg(npi_zeta, c, wc)
!    
!    call gauleg(npi_xi-1, ebar, webar)
!    call gauleg(npi_eta-1, nbar, wnbar)
!    call gauleg(npi_zeta-1, cbar, wcbar)
!    
!    !Calculate volume of the element
!    Ke = 0.0d0
!    Finte = 0.0d0
!    count = 0
!    do i=1,npi_xi
!        do j=1,npi_eta
!            do k=1,npi_zeta
!                
!                count = count + 1
!                
!                xi = e(i)
!                eta = n(j)
!                zeta = c(k)
!                
!                call ShapeFunc(nel,xi,eta,zeta,R,dRdx,detj)
!                
!                vol = vol + detj*we(i)*wn(j)*wc(k)
!                deter(count) = detj
!                
!                continue
!                
!            end do  
!        end do
!    end do
!    
!    !Gauss points cycle
!    cpi = 0
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
!                
!                !Weight factor
!                gwt = we(i)*wn(j)*wc(k)*detj
!                
!                Bglob = 0.0d0
!                Bdil = 0.0d0
!                Bdev = 0.0d0
!                BBar = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
!                    Bglob(1,k1*3-2) = dRdx(k1,1)
!                    Bglob(2,k1*3-1) = dRdx(k1,2)
!                    Bglob(3,k1*3  ) = dRdx(k1,3)
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
!                end do
!                
!                
!                
!                !#####################################################################################################################################
!                !#####################################################################################################################################
!                do k1=2,ncpx+p
!                    ubar(k1-1) = u_knot(k1)
!                end do
!                
!                do k1=2,ncpy+q
!                    vbar(k1-1) = v_knot(k1)
!                end do 
!                
!                do k1=2,ncpz+w
!                    wbar(k1-1) = w_knot(k1)
!                end do     
!                
!                call ShapeFuncBar(nel,xi,eta,zeta,INN,IEN,B_net,ubar,vbar,wbar,p-1,q-1,w-1,ncpx,ncpy,ncpz,nds,nnodes,nelems,RBar,dRdxBar,detjnew)
!                
!                aux1 = 0.0d0
!                do k5=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1)
!                    aux1(k5,1) = Rbar(k5)
!                end do
!                !#####################################################################################################################################
!                !#####################################################################################################################################
!                
!                
!                
!                
!                count = 0
!                BBar = 0.0d0
!                matM = 0.0d0
!                volnew = 0.0d0
!                do k2=1,npi_xi-1
!                    do k3=1,npi_eta-1
!                        do k4=1,npi_zeta-1
!                            
!                            
!                            
!!                            xi   = e(i)
!!                            eta  = n(j)
!!                            zeta = c(k)
!!                            
!!                            call ShapeFunc(nel,xi,eta,zeta,R,dRdx,detj)
!!                            
!!                            do k1=2,ncpx+p
!!                                ubar(k1-1) = u_knot(k1)
!!                            end do
!!                            
!!                            do k1=2,ncpy+q
!!                                vbar(k1-1) = v_knot(k1)
!!                            end do 
!!                            
!!                            do k1=2,ncpz+w
!!                                wbar(k1-1) = w_knot(k1)
!!                            end do     
!!                            
!!                            call ShapeFuncBar(nel,xi,eta,zeta,INN,IEN,B_net,ubar,vbar,wbar,p-1,q-1,w-1,ncpx,ncpy,ncpz,nds,nnodes,nelems,RBar,dRdxBar,detjnew)
!!                            count = 0
!!                            do k5=1,npi_xi
!!                            do k6=1,npi_eta
!!                            do k7=1,npi_zeta
!!                            count = count + 1
!!                            matM = matM + matmul(aux1,transpose(aux1))*we(k5)*wn(k6)*wc(k7)*deter(count)
!!                            end do
!!                            end do
!!                            end do
!                            
!                            matM = matM + matmul(aux1,transpose(aux1))*we(i)*wn(j)*wc(k)*detj
!                            !matM = matM + matmul(aux1,transpose(aux1))*webar(k2)*wnbar(k3)*wcbar(k4)*detjnew
!                            
!                            volnew = volnew + webar(k2)*wnbar(k3)*wcbar(k4)*detjnew
!                            
!!                            do k1=1,(p+1)*(q+1)*(w+1)
!!                                do k5=1,8
!!                                    do k6=1,8
!!                                    
!!                                        Bbar(1,k1*3-2) =  Bbar(1,k1*3-2) + dRdx(k1,1)*Rbar(k5)/lixo1*Rbar(k6)*webar(k2)*wnbar(k3)*wcbar(k4)*detj/gwt
!!                                        Bbar(2,k1*3-2) =  Bbar(2,k1*3-2) + dRdx(k1,1)*Rbar(k5)/lixo1*Rbar(k6)*webar(k2)*wnbar(k3)*wcbar(k4)*detj/gwt
!!                                        Bbar(3,k1*3-2) =  Bbar(3,k1*3-2) + dRdx(k1,1)*Rbar(k5)/lixo1*Rbar(k6)*webar(k2)*wnbar(k3)*wcbar(k4)*detj/gwt
!!                                        Bbar(1,k1*3-1) =  Bbar(1,k1*3-1) + dRdx(k1,2)*Rbar(k5)/lixo2*Rbar(k6)*webar(k2)*wnbar(k3)*wcbar(k4)*detj/gwt
!!                                        Bbar(2,k1*3-1) =  Bbar(2,k1*3-1) + dRdx(k1,2)*Rbar(k5)/lixo2*Rbar(k6)*webar(k2)*wnbar(k3)*wcbar(k4)*detj/gwt
!!                                        Bbar(3,k1*3-1) =  Bbar(3,k1*3-1) + dRdx(k1,2)*Rbar(k5)/lixo2*Rbar(k6)*webar(k2)*wnbar(k3)*wcbar(k4)*detj/gwt
!!                                        Bbar(1,k1*3  ) =  Bbar(1,k1*3  ) + dRdx(k1,3)*Rbar(k5)/lixo3*Rbar(k6)*webar(k2)*wnbar(k3)*wcbar(k4)*detj/gwt
!!                                        Bbar(2,k1*3  ) =  Bbar(2,k1*3  ) + dRdx(k1,3)*Rbar(k5)/lixo3*Rbar(k6)*webar(k2)*wnbar(k3)*wcbar(k4)*detj/gwt
!!                                        Bbar(3,k1*3  ) =  Bbar(3,k1*3  ) + dRdx(k1,3)*Rbar(k5)/lixo3*Rbar(k6)*webar(k2)*wnbar(k3)*wcbar(k4)*detj/gwt
!!                                        
!!                                    end do
!!                                end do
!!                            end do
!                            
!                            continue
!                            
!                        enddo
!                    enddo
!                enddo
!                
!                lixo=0.0d0
!                do k1=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1)
!                    lixo(k1) = sum(matM(:,k1))
!                end do
!                
!                MatM = 0.0d0
!                do k1=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1)
!                    matM(k1,k1) = lixo(k1)
!                end do
!                
!                Minv = 0.0d0
!                Minv = MatM
!                call gaussj(Minv,8,temp1,error)
!                
!                !temp11 = matmul(matmul(transpose(aux1),Minv),aux1)
!
!!                temp11 = 0.0d0
!!                k1 = 0
!!                do k2=1,npi_xi-1
!!                    do k3=1,npi_eta-1
!!                        do k4=1,npi_zeta-1
!!                            k1 = k1 + 1
!!                            temp11(1,1) = temp11(1,1) + aux1(k1,1)*Minv(k1,k1)*aux1(k1,1)*we(i)*wn(j)*wc(k)*detj !*webar(k2)*wnbar(k3)*wcbar(k4)*detjnew !*we(i)*wn(j)*wc(k)*detj
!!                        end do
!!                    end do
!!                end do
!                
!                
!                
!!                BBar = 0.0d0
!!                do k1=1,(p+1)*(q+1)*(w+1)
!!                    count = 0
!!                    do k2=1,npi_xi-1
!!                        do k3=1,npi_eta-1
!!                            do k4=1,npi_zeta-1
!!                                
!!                                count = count + 1
!!                                
!!                                BBar(1,k1*3-2) = BBar(1,k1*3-2) + Bglob(1,k1*3-2)*aux1(count,1)*Minv(count,count)*aux1(count,1)*we(i)*wn(j)*wc(k)*detj
!!                                BBar(2,k1*3-2) = BBar(2,k1*3-2) + Bglob(1,k1*3-2)*aux1(count,1)*Minv(count,count)*aux1(count,1)*we(i)*wn(j)*wc(k)*detj
!!                                BBar(3,k1*3-2) = BBar(3,k1*3-2) + Bglob(1,k1*3-2)*aux1(count,1)*Minv(count,count)*aux1(count,1)*we(i)*wn(j)*wc(k)*detj
!!                                
!!                                BBar(1,k1*3-1) = BBar(1,k1*3-1) + Bglob(2,k1*3-1)*aux1(count,1)*Minv(count,count)*aux1(count,1)*we(i)*wn(j)*wc(k)*detj
!!                                BBar(2,k1*3-1) = BBar(2,k1*3-1) + Bglob(2,k1*3-1)*aux1(count,1)*Minv(count,count)*aux1(count,1)*we(i)*wn(j)*wc(k)*detj
!!                                BBar(3,k1*3-1) = BBar(3,k1*3-1) + Bglob(2,k1*3-1)*aux1(count,1)*Minv(count,count)*aux1(count,1)*we(i)*wn(j)*wc(k)*detj
!!                                
!!                                BBar(1,k1*3  ) = BBar(1,k1*3  ) + Bglob(3,k1*3  )*aux1(count,1)*Minv(count,count)*aux1(count,1)*we(i)*wn(j)*wc(k)*detj
!!                                BBar(2,k1*3  ) = BBar(2,k1*3  ) + Bglob(3,k1*3  )*aux1(count,1)*Minv(count,count)*aux1(count,1)*we(i)*wn(j)*wc(k)*detj
!!                                BBar(3,k1*3  ) = BBar(3,k1*3  ) + Bglob(3,k1*3  )*aux1(count,1)*Minv(count,count)*aux1(count,1)*we(i)*wn(j)*wc(k)*detj
!!                            end do
!!                        end do
!!                    end do
!!                end do
!                
!                BBar = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
!                    count = 0
!                    do k2=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1)
!                        do k3=1,(npi_xi-1)*(npi_eta-1)*(npi_zeta-1)
!                                
!                                BBar(1,k1*3-2) = BBar(1,k1*3-2) + Bglob(1,k1*3-2)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(i)*wn(j)*wc(k)*detj
!                                BBar(2,k1*3-2) = BBar(2,k1*3-2) + Bglob(1,k1*3-2)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(i)*wn(j)*wc(k)*detj
!                                BBar(3,k1*3-2) = BBar(3,k1*3-2) + Bglob(1,k1*3-2)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(i)*wn(j)*wc(k)*detj
!                                
!                                BBar(1,k1*3-1) = BBar(1,k1*3-1) + Bglob(2,k1*3-1)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(i)*wn(j)*wc(k)*detj
!                                BBar(2,k1*3-1) = BBar(2,k1*3-1) + Bglob(2,k1*3-1)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(i)*wn(j)*wc(k)*detj
!                                BBar(3,k1*3-1) = BBar(3,k1*3-1) + Bglob(2,k1*3-1)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(i)*wn(j)*wc(k)*detj
!                                
!                                BBar(1,k1*3  ) = BBar(1,k1*3  ) + Bglob(3,k1*3  )*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(i)*wn(j)*wc(k)*detj
!                                BBar(2,k1*3  ) = BBar(2,k1*3  ) + Bglob(3,k1*3  )*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(i)*wn(j)*wc(k)*detj
!                                BBar(3,k1*3  ) = BBar(3,k1*3  ) + Bglob(3,k1*3  )*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(i)*wn(j)*wc(k)*detj
!                                
!                                
!!                                count = 0
!!                                do k5=1,npi_xi
!!                                    do k6=1,npi_eta
!!                                        do k7=1,npi_zeta
!!                                            count = count + 1
!!                                            BBar(1,k1*3-2) = BBar(1,k1*3-2) + Bglob(1,k1*3-2)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!!                                            BBar(2,k1*3-2) = BBar(2,k1*3-2) + Bglob(1,k1*3-2)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!!                                            BBar(3,k1*3-2) = BBar(3,k1*3-2) + Bglob(1,k1*3-2)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!!                                            
!!                                            BBar(1,k1*3-1) = BBar(1,k1*3-1) + Bglob(2,k1*3-1)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!!                                            BBar(2,k1*3-1) = BBar(2,k1*3-1) + Bglob(2,k1*3-1)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!!                                            BBar(3,k1*3-1) = BBar(3,k1*3-1) + Bglob(2,k1*3-1)*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!!                                            
!!                                            BBar(1,k1*3  ) = BBar(1,k1*3  ) + Bglob(3,k1*3  )*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!!                                            BBar(2,k1*3  ) = BBar(2,k1*3  ) + Bglob(3,k1*3  )*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!!                                            BBar(3,k1*3  ) = BBar(3,k1*3  ) + Bglob(3,k1*3  )*aux1(k2,1)*MInv(k2,k3)*aux1(k3,1)*we(k5)*wn(k6)*wc(k7)*deter(count)
!!                                        end do
!!                                    end do
!!                                end do
!                                
!                        end do
!                    end do
!                end do
!                
!!                BBar = 0.0d0
!!                do k1=1,(p+1)*(q+1)*(w+1)
!!                    
!!                    BBar(1,k1*3-2) = Bglob(1,k1*3-2)*temp11(1,1)/3.0d0
!!                    BBar(2,k1*3-2) = Bglob(1,k1*3-2)*temp11(1,1)/3.0d0
!!                    BBar(3,k1*3-2) = Bglob(1,k1*3-2)*temp11(1,1)/3.0d0
!!                    
!!                    BBar(1,k1*3-1) = Bglob(2,k1*3-1)*temp11(1,1)/3.0d0
!!                    BBar(2,k1*3-1) = Bglob(2,k1*3-1)*temp11(1,1)/3.0d0
!!                    BBar(3,k1*3-1) = Bglob(2,k1*3-1)*temp11(1,1)/3.0d0
!!                    
!!                    BBar(1,k1*3  ) = Bglob(3,k1*3  )*temp11(1,1)/3.0d0
!!                    BBar(2,k1*3  ) = Bglob(3,k1*3  )*temp11(1,1)/3.0d0
!!                    BBar(3,k1*3  ) = Bglob(3,k1*3  )*temp11(1,1)/3.0d0
!!
!!                end do
!                
!                
!                BBar  = Bbar/3.0d0
!                Bdev  = Bglob - Bdil
!                Bglob = Bdev + Bbar
!                
!                continue
!                
!                
!                
!                
!                
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
!!                do k1=1,ntens
!!                    SEnergy = SEnergy + tempstr(k1,1)*deform(k1,1)*gwt/2.0d0
!!                end do
!                
!                continue
!                
!            end do  
!        end do
!    end do
!    
!
!    continue
!
!end subroutine ElHex27BBar


!subroutine ElHex27BBar(nel,Ke,el_ddisp,Finte)
!    
!    use Mod_Variables
!    implicit none
!    
!    integer(4)::nel
!    integer(4)::cpi
!    
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::el_ddisp
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(OUT)::Ke
!    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(OUT)::Finte
!    
!    integer(4),parameter::npi_xi= 3, npi_eta = 3, npi_zeta = 3
!    
!    integer(4)::i,j,jj,k,k2,k3,k4,ni,nj,nk,nel_nza,k1,error
!    real(8)::xi,eta,zeta
!    real(8),dimension(npi_xi)::e,we
!    real(8),dimension(npi_xi-1)::ebar,webar
!    real(8),dimension(npi_eta)::n,wn
!    real(8),dimension(npi_eta-1)::nbar,wnbar
!    real(8),dimension(npi_zeta)::c,wc
!    real(8),dimension(npi_zeta-1)::cbar,wcbar
!    
!    real(8),dimension((p+1)*(q+1)*(w+1))::R
!    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdx
!    
!    !B-Bar Method -----------------------------------------
!    real(8),dimension(ncpx+p-1)::ubar
!    real(8),dimension(ncpy+q-1)::vbar
!    real(8),dimension(ncpz+w-1)::wbar
!    real(8),dimension((p)*(q)*(w))::Rbar
!    real(8),dimension((p)*(q)*(w),nds)::dRdxbar
!
!    real(8)::detJ,gwt,detjnew
!    
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob,Bloc
!    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bdil,Bdev,BBar
!    
!    real(8),dimension(6,6)::matD
!    real(8),dimension(ntens,1)::tempstr,deform
!    
!    real(8),dimension(8,8)::matM
!    real(8),dimension(8)::temp1
!    
!    real(8)::vol
!    
!    cpi = 0
!    
!    vol = 0.0d0
!    
!    we = 1.0d0
!    wn = 1.0d0
!    wc = 1.0d0
!    call gauleg(npi_xi, e, we)
!    call gauleg(npi_eta, n, wn)
!    call gauleg(npi_zeta, c, wc)
!    
!    call gauleg(npi_xi-1, ebar, webar)
!    call gauleg(npi_eta-1, nbar, wnbar)
!    call gauleg(npi_zeta-1, cbar, wcbar)
!    
!    !Calculate volume of the element
!    Ke = 0.0d0
!    Finte = 0.0d0
!    do i=1,npi_xi
!        do j=1,npi_eta
!            do k=1,npi_zeta
!
!                xi = e(i)
!                eta = n(j)
!                zeta = c(k)
!                
!                call ShapeFunc(nel,xi,eta,zeta,R,dRdx,detj)
!                
!                vol = vol + detj
!                
!                continue
!                
!            end do  
!        end do
!    end do
!    
!    !Gauss points cycle
!    cpi = 0
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
!                
!                !Weight factor
!                gwt = we(i)*wn(j)*wc(k)*detj
!                
!                Bglob = 0.0d0
!                Bdil = 0.0d0
!                Bdev = 0.0d0
!                BBar = 0.0d0
!                do k1=1,(p+1)*(q+1)*(w+1)
!                    
!                    Bglob(1,k1*3-2) = dRdx(k1,1)
!                    Bglob(2,k1*3-1) = dRdx(k1,2)
!                    Bglob(3,k1*3  ) = dRdx(k1,3)
!                    Bglob(4,k1*3-2) = dRdx(k1,2)
!                    Bglob(4,k1*3-1) = dRdx(k1,1)
!                    Bglob(5,k1*3-2) = dRdx(k1,3)
!                    Bglob(5,k1*3  ) = dRdx(k1,1)
!                    Bglob(6,k1*3-1) = dRdx(k1,3)
!                    Bglob(6,k1*3  ) = dRdx(k1,2)
!                    
!                    Bdil(1,k1*3-2) = dRdx(k1,1)
!                    Bdil(2,k1*3-2) = dRdx(k1,1)
!                    Bdil(3,k1*3-2) = dRdx(k1,1)
!                    Bdil(1,k1*3-1) = dRdx(k1,2)
!                    Bdil(2,k1*3-1) = dRdx(k1,2)
!                    Bdil(3,k1*3-1) = dRdx(k1,2)
!                    Bdil(1,k1*3  ) = dRdx(k1,3)
!                    Bdil(2,k1*3  ) = dRdx(k1,3)
!                    Bdil(3,k1*3  ) = dRdx(k1,3)
!                    
!                end do
!                
!                do k2=1,2
!                    do k3=1,2
!                        do k4=1,2
!                            xi = ebar(k2)
!                            eta = nbar(k3)
!                            zeta = cbar(k4)
!                            
!                            call ShapeFunc(nel,xi,eta,zeta,R,dRdx,detjnew)
!                            
!                            do k1=2,ncpx+p
!                                ubar(k1-1) = u_knot(k1)
!                            end do
!                            
!                            do k1=2,ncpy+q
!                                vbar(k1-1) = v_knot(k1)
!                            end do 
!                            
!                            do k1=2,ncpz+w
!                                wbar(k1-1) = w_knot(k1)
!                            end do 
!                            
!                            call ShapeFuncBar(nel,xi,eta,zeta,INN,IEN,B_net,ubar,vbar,wbar,p-1,q-1,w-1,ncpx,ncpy,ncpz,nds,nnodes,nelems,RBar,dRdxBar,detjnew)
!                            
!                            do k1=1,(p+1)*(q+1)*(w+1)
!                               
!                                Bbar(1,k1*3-2) =  Bbar(1,k1*3-2) + dRdx(k1,1)
!                                Bbar(2,k1*3-2) =  Bbar(2,k1*3-2) + dRdx(k1,1)
!                                Bbar(3,k1*3-2) =  Bbar(3,k1*3-2) + dRdx(k1,1)
!                                Bbar(1,k1*3-1) =  Bbar(1,k1*3-1) + dRdx(k1,2)
!                                Bbar(2,k1*3-1) =  Bbar(2,k1*3-1) + dRdx(k1,2)
!                                Bbar(3,k1*3-1) =  Bbar(3,k1*3-1) + dRdx(k1,2)
!                                Bbar(1,k1*3  ) =  Bbar(1,k1*3  ) + dRdx(k1,3)
!                                Bbar(2,k1*3  ) =  Bbar(2,k1*3  ) + dRdx(k1,3)
!                                Bbar(3,k1*3  ) =  Bbar(3,k1*3  ) + dRdx(k1,3)
!                    
!                            end do
!                        end do
!                    end do
!                end do
!                
!                Bdil = Bdil/3.0d0
!                BBar = Bbar/3.0d0
!                
!                Bdev = Bglob - Bdil
!                
!                Bglob = Bdev + Bbar
!                
!!                Bglob = Bdev
!!                do k1=1,8
!!                    Bglob(1,k1*3-2) = Bdev(1,k1*3-2)  + coef(k1,1)
!!                    Bglob(2,k1*3-2) = Bdev(2,k1*3-2)  + coef(k1,1)
!!                    Bglob(3,k1*3-2) = Bdev(3,k1*3-2)  + coef(k1,1)
!!                    Bglob(1,k1*3-1) = Bdev(1,k1*3-1)  + coef(k1,2)
!!                    Bglob(2,k1*3-1) = Bdev(2,k1*3-1)  + coef(k1,2)
!!                    Bglob(3,k1*3-1) = Bdev(3,k1*3-1)  + coef(k1,2)
!!                    Bglob(1,k1*3  ) = Bdev(1,k1*3  )  + coef(k1,3)
!!                    Bglob(2,k1*3  ) = Bdev(2,k1*3  )  + coef(k1,3)
!!                    Bglob(3,k1*3  ) = Bdev(3,k1*3  )  + coef(k1,3)
!!                end do
!                
!                
!                continue
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
!!                do k1=1,ntens
!!                    SEnergy = SEnergy + tempstr(k1,1)*deform(k1,1)*gwt/2.0d0
!!                end do
!                
!                continue
!                
!            end do  
!        end do
!    end do
!    
!
!    continue
!
!end subroutine ElHex27BBar