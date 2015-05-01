!------------------------------------------------------------------------------------------------------------
!
! Hexaedral element with a 2x2x2 gaussian quadrature and B-bar method, Q1\Q0 from the work of Elguedj et.al.
! "B-bar and F-bar projection methods for nearly incompressible linear and nonlinear elsaticity and 
! plasticity using higher-order NURBS elements" (2008)
!
! Input: idx1, idx2 - NURBS coordinates (element indexes in the index space) (currently not being used)
!        nel - Element number
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

subroutine ElHex8BBar(nel,Ke,el_ddisp,Finte,MatA,MatC,MatM)
    
    use Mod_Variables
    implicit none
    
    integer(4)::nel
    integer(4)::cpi
    
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::el_ddisp
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(OUT)::Ke
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(OUT)::Finte
    
    integer(4),parameter::npi_xi= 2, npi_eta = 2, npi_zeta = 2
    
    integer(4)::i,j,jj,k,k2,k3,k4,k5,k6,k7,k8,ni,nj,nk,nel_nza,k1,error,count,count2
    real(8)::xi,eta,zeta
    real(8),dimension(npi_xi)::e,we
    real(8),dimension(npi_eta)::n,wn
    real(8),dimension(npi_zeta)::c,wc
    
    real(8),dimension((p+1)*(q+1)*(w+1))::R
    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdx
    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob,Bloc
    
    real(8)::detj,detjnew,gwt
    
    !real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bdil,Bdev,BBar
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
                !Bdil = 0.0d0
                !Bdev = 0.0d0
                !BBar = 0.0d0
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
                    
!                    Bdil(1,k1*3-2) = dRdx(k1,1)/3.0d0
!                    Bdil(2,k1*3-2) = dRdx(k1,1)/3.0d0
!                    Bdil(3,k1*3-2) = dRdx(k1,1)/3.0d0
!                    Bdil(1,k1*3-1) = dRdx(k1,2)/3.0d0
!                    Bdil(2,k1*3-1) = dRdx(k1,2)/3.0d0
!                    Bdil(3,k1*3-1) = dRdx(k1,2)/3.0d0
!                    Bdil(1,k1*3  ) = dRdx(k1,3)/3.0d0
!                    Bdil(2,k1*3  ) = dRdx(k1,3)/3.0d0
!                    Bdil(3,k1*3  ) = dRdx(k1,3)/3.0d0
                end do
                
!                Bdev  = Bglob - Bdil
                
                
                
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
                call ShapeFuncBar(nel,xi,eta,zeta,INN,IEN,B_net,ubar,vbar,wbar,p-1,q-1,w-1,ncpx,ncpy,&
                                & ncpz,nds,nnodes,nelems,RBar,dRdxBar,detjnew)
                
                Bloc = Bglob
                call MatPlastic3D(iprops,nds,(p+1)*(q+1)*(w+1),props,Bloc,6,el_ddisp,stress(nel,cpi,:),&
                            & strain(nel,cpi,:),dstrain(nel,cpi,:),dstress(nel,cpi,:),hard(nel,cpi),matD)
                
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
                
!                do k1=1,ntens
!                    SEnergy = SEnergy + tempstr(k1,1)*deform(k1,1)*gwt/2.0d0
!                end do
                
                continue
                
            end do  
        end do
    end do

    continue

end subroutine ElHex8BBar
