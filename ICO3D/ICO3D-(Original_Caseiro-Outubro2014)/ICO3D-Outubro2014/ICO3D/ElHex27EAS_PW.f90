

subroutine ElHex27EAS_PW(nel,Ke,el_ddisp,Finte,MatA,MatC,MatV)
    
    use Mod_Variables
    implicit none
    
    integer(4)::nel
    integer(4)::cpi
    
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::el_ddisp
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(OUT)::Ke
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(OUT)::Finte
    
    real(8),dimension(nalpha,nalpha),intent(OUT)::MatV
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(OUT)::MatA
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,nalpha),intent(OUT)::MatC
    
    integer(4),parameter::npi_xi= 3, npi_eta = 3, npi_zeta = 3
    
    integer(4)::i,j,jj,k,k2,k3,k4,k5,k6,k7,k8,ni,nj,nk,nel_nza,k1,error,count,count2
    real(8)::xi,eta,zeta
    real(8),dimension(npi_xi)::e,we
    real(8),dimension(npi_eta)::n,wn
    real(8),dimension(npi_zeta)::c,wc
    
    real(8),dimension((p+1)*(q+1)*(w+1))::R
    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdx
    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob,Bloc
    real(8),dimension(ntens,1)::tempstr,deform
    
    real(8)::detj,detjnew,gwt
    
    real(8),dimension(6,6)::matD,MatT
    real(8),dimension(6,nalpha)::Malpha,Balpha
    real(8),dimension(nalpha,1)::aux1
    real(8),dimension(nalpha)::temp1
    
    integer(4)::ilixo
    real(8)::detj0
    real(8),dimension((p+1)*(q+1)*(w+1))::R0
    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdx0
    real(8),dimension(3,3)::jac0,jac0inv
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1)::updtdisp
    real(8),dimension(3)::temp3
    real(8)::dBdxi,dBdeta,dBdzeta
    real(8)::a1,a2,a3
    
    cpi = 0
    
    we = 1.0d0
    wn = 1.0d0
    wc = 1.0d0
    call gauleg(npi_xi, e, we)
    call gauleg(npi_eta, n, wn)
    call gauleg(npi_zeta, c, wc)
    
    Ke = 0.0d0
    MatA = 0.0d0
    MatC = 0.0d0
    MatV = 0.0d0
    
    !Patch centre -------------------
    updtdisp = 0.0d0
    call ShapeFunc2(nel,0.0d0,0.0d0,0.0d0,R,dRdx,detj0,jac0,updtdisp)

    temp3 = 0.0d0
    ilixo = 0
    jac0inv = jac0
    call gaussj (jac0inv,nds,temp3,ilixo)
    call TransformationMat3D(jac0inv,MatT)
    
    
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
                
                Malpha = 0.0d0
                
                Malpha(1,1) = -xi*(1.0-eta*eta)*(1.0-zeta*zeta);
                Malpha(2,2) = -eta*(1.0-xi*xi)*(1.0-zeta*zeta);
                Malpha(3,3) = -zeta*(1.0-xi*xi)*(1.0-eta*eta);
                
                Malpha(1,4) = 2.0d0*xi*eta*(1.0-zeta*zeta);
                Malpha(2,4) = 2.0d0*xi*eta*(1.0-zeta*zeta);
                Malpha(3,4) = 2.0d0*xi*eta*(1.0-zeta*zeta);
                
                Malpha(1,5) = 2.0d0*xi*zeta*(1.0-eta*eta);
                Malpha(2,5) = 2.0d0*xi*zeta*(1.0-eta*eta);
                Malpha(3,5) = 2.0d0*xi*zeta*(1.0-eta*eta);
                
                Malpha(1,6) = 2.0d0*eta*zeta*(1.0-xi*xi);
                Malpha(2,6) = 2.0d0*eta*zeta*(1.0-xi*xi);
                Malpha(3,6) = 2.0d0*eta*zeta*(1.0-xi*xi);
            
                Malpha(1,7) = 4.0*xi*eta*zeta;
                Malpha(2,8) = 4.0*xi*eta*zeta;
                Malpha(3,9) = 4.0*xi*eta*zeta;
           
                Malpha(1,10) = -(1.0-eta*eta)*(1.0-zeta*zeta);
                Malpha(2,11) = -(1.0-xi*xi)*(1.0-zeta*zeta);
                Malpha(3,12) = -(1.0-xi*xi)*(1.0-eta*eta);

                
                
                
!                Malpha(4,1) = -xi*(1.0-eta*eta)*(1.0-zeta*zeta);
!                Malpha(4,2) = -eta*(1.0-xi*xi)*(1.0-zeta*zeta);
!                Malpha(5,3) = -xi*(1.0-eta*eta)*(1.0-zeta*zeta);
!                Malpha(5,4) = -zeta*(1.0-xi*xi)*(1.0-eta*eta);
!                Malpha(6,5) = -eta*(1.0-xi*xi)*(1.0-zeta*zeta);
!                Malpha(6,6) = -zeta*(1.0-xi*xi)*(1.0-eta*eta);
!                
!                Malpha(4,7) = 2.0d0*xi*zeta*(1.0-eta*eta);
!                Malpha(4,8) = 2.0d0*eta*zeta*(1.0-xi*xi);
!                Malpha(5,9) = 2.0d0*xi*eta*(1.0-zeta*zeta);
!                Malpha(5,10) = 2.0d0*eta*zeta*(1.0-xi*xi);
!                Malpha(6,11) = 2.0d0*xi*eta*(1.0-zeta*zeta);
!                Malpha(6,12) = 2.0d0*xi*zeta*(1.0-eta*eta);
                
                Balpha = 0.0d0
                Balpha = matmul(MatT,Malpha)*detj0/detj
                
                
                Bloc = Bglob
                call MatPlastic3D(iprops,nds,(p+1)*(q+1)*(w+1),props,Bloc,6,el_ddisp,stress(nel,cpi,:),&
                        & strain(nel,cpi,:),dstrain(nel,cpi,:),dstress(nel,cpi,:),hard(nel,cpi),matD)
                
                tempstr = 0.0d0
                do k1=1,ntens
                    tempstr(k1,1) = stress(nel,cpi,k1)
                    deform(k1,1) = dstrain(nel,cpi,k1)
                end do
                
                !Matrix A ------------------------------------------------------------
                MatA = MatA + matmul(matmul(transpose(Bglob),MatD),Bglob)*gwt
                !---------------------------------------------------------------------
                
                !Matrix C ------------------------------------------------------------
                MatC = MatC + matmul(matmul(transpose(Bglob),MatD),Balpha)*gwt
                !---------------------------------------------------------------------
                
                !Matrix V ------------------------------------------------------------
                MatV = MatV + matmul(matmul(transpose(Balpha),MatD),Balpha)*gwt
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

end subroutine
