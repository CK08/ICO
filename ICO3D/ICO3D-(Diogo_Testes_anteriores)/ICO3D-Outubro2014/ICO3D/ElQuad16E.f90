!------------------------------------------------------------------------------------------------------------
!
! General quadrilateral element for plain strain problems with a 4x4 gaussian quadrature
!
! Input: idx1, idx2 - NURBS coordinates (element indexes in the index space) (currently not being used)
!        nel - Element number
!
! Output: Ke - Elemental stiffness matrix
!------------------------------------------------------------------------------------------------------------


subroutine ElQuad16E(nel,Ke,el_ddisp,Finte)
    
    use Mod_Variables
    implicit none
    
    integer(4)::nel
    integer(4)::cpi
    
    real(8),dimension((p+1)*(q+1)*nds,1),intent(IN)::el_ddisp
    real(8),dimension((p+1)*(q+1)*nds,(p+1)*(q+1)*nds),intent(OUT)::Ke
    real(8),dimension((p+1)*(q+1)*nds,1),intent(OUT)::Finte
    
    
    integer(4),parameter::npi_xi= 4, npi_eta = 4
    
    integer(4)::i,j,ni,nj,nk,nel_nza,k1
    real(8)::xi,eta
    real(8),dimension(npi_xi)::e,we
    real(8),dimension(npi_eta)::n,wn
    
    real(8),dimension((p+1)*(q+1))::R
    real(8),dimension((p+1)*(q+1),nds)::dRdx
    real(8)::detJ,gwt
    
    real(8),dimension(3,(p+1)*(q+1)*nds)::Bglob,Bloc
    
    real(8),dimension(3,3)::matD
    real(8),dimension(ntens,1)::tempstr
    
    cpi = 0

    we=1.0d0
    wn=1.0d0
    call gauleg(npi_xi, e, we)
    call gauleg(npi_eta, n, wn)
    
    !Gauss points cycle
    Ke = 0.0d0
    Finte = 0.0d0
    do i=1,npi_xi
        do j=1,npi_eta
            
            cpi = cpi + 1
            MatD(:,:) = TmatD(nel,cpi,:,:)
            
            xi = e(i)
            eta = n(j)
            
            call ShapeFunc(nel,xi,eta,R,dRdx,detj)
            
            !Weight factor
            gwt = we(i)*wn(j)*detj
            
            Bglob = 0.0d0
            do k1=1,(p+1)*(q+1)
                Bglob(1,k1*2-1) = dRdx(k1,1)
                Bglob(2,k1*2  ) = dRdx(k1,2)
                Bglob(3,k1*2-1) = dRdx(k1,2)
                Bglob(3,k1*2  ) = dRdx(k1,1)
            end do
            
            Bloc = Bglob
            
            !call MaterialPStrain(ndi,nshr,ntens,iprops,nodedof,nnodes,props,B,Blines,ddisp,stress,strain,dstrain,dstress,hard,matD,matf0,matf,dtempo)
            call MaterialPStrain(iprops,nds,(p+1)*(q+1),props,Bloc,3,el_ddisp,stress(nel,cpi,:),&
                     & strain(nel,cpi,:),dstrain(nel,cpi,:),dstress(nel,cpi,:),hard(nel,cpi),matD)
            
            tempstr = 0.0d0
            do k1=1,ntens
                tempstr(k1,1) = stress(nel,cpi,k1)
            end do
            
            Ke = Ke + matmul(matmul(transpose(Bglob),matD),Bglob)*gwt
            
            Finte = Finte + matmul(transpose(Bloc),tempstr)*gwt
            
            continue
          
        end do
    end do

    continue

end subroutine ElQuad16E
