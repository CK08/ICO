

subroutine ElHex27ANS_PW(nel,Ke,el_ddisp,Finte,MatA,MatC,MatV)
    
    use Mod_Variables
    implicit none
    
    integer(4)::nel
    integer(4)::cpi
    
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::el_ddisp
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(OUT)::Ke
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(OUT)::Finte
    
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(OUT)::MatV
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(OUT)::MatA
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(OUT)::MatC
    
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
    real(8),dimension(nds)::temp1
    
    integer(4)::ilixo
    real(8)::detj0
    real(8),dimension((p+1)*(q+1)*(w+1))::R0
    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdx0
    real(8),dimension(3,3)::jac0,jac0inv
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1)::updtdisp
    real(8),dimension(3)::temp3
    real(8)::dBdxi,dBdeta,dBdzeta
    real(8)::a1,a2,a3
    
    real(8)::detjA,A,B
    real(8),dimension(nds,nds)::dxdxiA,dxdxi,jacA,jac,jacinv
    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BNat,BNat2,BCN2,BCN
    real(8),dimension((p+1)*(q+1)*(w+1))::RA
    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdxi,dRdxii
    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdxA,dRdxiA,dRdxiiA
    real(8),dimension(16)::shpANS
    real(8),dimension(16,3)::coordANS
    real(8),dimension(3,1)::dRdx_teste,temp31
    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::BNatural,BN2
    real(8),dimension(nds,nds)::Mat33
    real(8),dimension(nds*2,nds*2)::TGC,TCL
    
    cpi = 0
    
    we = 1.0d0
    wn = 1.0d0
    wc = 1.0d0
    call gauleg(npi_xi, e, we)
    call gauleg(npi_eta, n, wn)
    call gauleg(npi_zeta, c, wc)
    
    updtdisp = 0.0d0
    
    Ke = 0.0d0
    MatA = 0.0d0
    MatC = 0.0d0
    MatV = 0.0d0
    
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
                
                !call ShapeFunc(nel,xi,eta,zeta,R,dRdx,detj)
                call ShapeFunc3(nel,xi,eta,zeta,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,updtdisp)
                
                jacinv = 0.0d0
                jacinv = jac
                call gaussj(jacinv,nds,temp1,ilixo)
                
                !Global to natural transformation matrix for end configuration
                Mat33 = jac
                call TransformationMat3D(mat33,TGC)
                
                Mat33 = jacinv
                call TransformationMat3D(mat33,TCL)
                
                !Weight factor
                gwt = we(i)*wn(j)*wc(k)*detj
                
                !Strain-displacement matrix -------------
                Bglob = 0.0d0
                do k1=1,(p+1)*(q+1)*(w+1)
                    Bglob(1,k1*3-2) = dRdx(k1,1)
                    Bglob(2,k1*3-1) = dRdx(k1,2)
                    Bglob(3,k1*3  ) = dRdx(k1,3)
                    !Bglob(4,k1*3-2) = dRdx(k1,2)
                    !Bglob(4,k1*3-1) = dRdx(k1,1)
                    Bglob(5,k1*3-2) = dRdx(k1,3)
                    Bglob(5,k1*3  ) = dRdx(k1,1)
                    Bglob(6,k1*3-1) = dRdx(k1,3)
                    Bglob(6,k1*3  ) = dRdx(k1,2)
                end do
                
                BNatural = matmul(TGC,Bglob)
                
                !ANS Stuff
                BCN = BNatural
                
                !ANS Implementation I
                A = sqrt(1.0d0/3.0d0)
                B = sqrt(3.0d0/5.0d0)
                
                !Coordinates
                coordANS = 0.0d0
                coordANS(1,1)  =     A; coordANS(1,2)  =     B
                coordANS(2,1)  =    -A; coordANS(2,2)  =     B
                coordANS(3,1)  =     A; coordANS(3,2)  = 0.0d0
                coordANS(4,1)  =    -A; coordANS(4,2)  = 0.0d0
                coordANS(5,1)  =     A; coordANS(5,2)  =    -B
                coordANS(6,1)  =    -A; coordANS(6,2)  =    -B
                coordANS(7,1)  =     B; coordANS(7,2)  =     A
                coordANS(8,1)  = 0.0d0; coordANS(8,2)  =     A
                coordANS(9,1)  =    -B; coordANS(9,2)  =     A
                coordANS(10,1) =     B; coordANS(10,2) =    -A
                coordANS(11,1) = 0.0d0; coordANS(11,2) =    -A
                coordANS(12,1) =    -B; coordANS(12,2) =    -A
                coordANS(13,1) =     A; coordANS(13,2) =     A
                coordANS(14,1) =    -A; coordANS(14,2) =     A
                coordANS(15,1) =     A; coordANS(15,2) =    -A
                coordANS(16,1) =    -A; coordANS(16,2) =    -A
                !-------------------------------------------------------------------
                !-------------------------------------------------------------------
                
                !Shape Functions
                shpANS = 0.0d0
                !-------------------------------------------------------------------
                shpANS(1)  = eta/(4.0d0*B)*(eta/B+1.0d0)*(1.0d0+xi/A) 
                shpANS(2)  = eta/(4.0d0*B)*(eta/B+1.0d0)*(1.0d0-xi/A)
                shpANS(3)  = 1.0d0/2.0d0*(1.0d0-eta*eta/(B*B))*(1.0d0+xi/A)
                shpANS(4)  = 1.0d0/2.0d0*(1.0d0-eta*eta/(B*B))*(1.0d0-xi/A)
                shpANS(5)  = eta/(4.0d0*B)*(eta/B-1.0d0)*(1.0d0+xi/A) 
                shpANS(6)  = eta/(4.0d0*B)*(eta/B-1.0d0)*(1.0d0-xi/A)
                !-------------------------------------------------------------------
                
                !-------------------------------------------------------------------
                shpANS(7)  = xi/(4.0d0*B)*(xi/B+1.0d0)*(1.0d0+eta/A)
                shpANS(8)  = 1.0d0/2.0d0*(1.0d0-xi*xi/(B*B))*(1.0d0+eta/A)
                shpANS(9)  = xi/(4.0d0*B)*(xi/B-1.0d0)*(1.0d0+eta/A)
                shpANS(10) = xi/(4.0d0*B)*(xi/B+1.0d0)*(1.0d0-eta/A)
                shpANS(11) = 1.0d0/2.0d0*(1.0d0-xi*xi/(B*B))*(1.0d0-eta/A)
                shpANS(12) = xi/(4.0d0*B)*(xi/B-1.0d0)*(1.0d0-eta/A)
                !-------------------------------------------------------------------

                !-------------------------------------------------------------------
                shpANS(13) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0+eta/A)
                shpANS(14) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0+eta/A)
                shpANS(15) = 1.0d0/4.0d0*(1.0d0+xi/A)*(1.0d0-eta/A)
                shpANS(16) = 1.0d0/4.0d0*(1.0d0-xi/A)*(1.0d0-eta/A)
                !-------------------------------------------------------------------
                
                BNat = 0.0d0
                
                !BNat(1,:) = BCN(1,:)
                !BNat(2,:) = BCN(2,:)
                !BNat(3,:) = BCN(3,:)
                !BNat(4,:) = BCN(4,:)
                !BNat(5,:) = BCN(5,:)
                !BNat(6,:) = BCN(6,:)
                
                !xi-eta components
                do k1=13,16
                    Bnat(4,:) = Bnat(4,:) + shpANS(k1)*BCN(4,:)
                end do
                
                Bloc=matmul(TCL,BNat)
                

                
                Bloc = Bglob
                call MatPlastic3D(iprops,nds,(p+1)*(q+1)*(w+1),props,Bloc,6,el_ddisp,stress(nel,cpi,:),strain(nel,cpi,:),dstrain(nel,cpi,:),dstress(nel,cpi,:),hard(nel,cpi),matD)
                
                tempstr = 0.0d0
                do k1=1,ntens
                    tempstr(k1,1) = stress(nel,cpi,k1)
                    deform(k1,1) = dstrain(nel,cpi,k1)
                end do
                
                !Matrix A ------------------------------------------------------------
                MatA = MatA + matmul(matmul(transpose(Bglob),MatD),Bglob)*gwt
                !---------------------------------------------------------------------
                
                !Matrix C ------------------------------------------------------------
                MatC = MatC + matmul(matmul(transpose(Bglob),MatD),Bnat)*gwt
                !---------------------------------------------------------------------
                
                !Matrix V ------------------------------------------------------------
                MatV = MatV + matmul(matmul(transpose(Bnat),MatD),Bnat)*gwt
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