subroutine ElHex27SRI(nel,Ke,el_ddisp,Finte)
    
    use Mod_Variables
    implicit none
    
    integer(4)::nel,count
    integer(4)::cpi
    
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::el_ddisp
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,(p+1)*(q+1)*(w+1)*nds),intent(OUT)::Ke
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(OUT)::Finte
    
    integer(4),parameter::npi_xi= 3, npi_eta = 3, npi_zeta = 3
    
    integer(4)::i,j,k,ni,nj,nk,nel_nza,k1
    real(8)::xi,eta,zeta
    real(8),dimension(npi_xi)::e,we
    real(8),dimension(npi_eta)::n,wn
    real(8),dimension(npi_zeta)::c,wc
    
    real(8),dimension(npi_xi-1)::er,wer
    real(8),dimension(npi_eta-1)::nr,wnr
    real(8),dimension(npi_zeta-1)::cr,wcr
    
    real(8),dimension((p+1)*(q+1)*(w+1))::R
    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdx

    real(8)::detJ,gwt
    
    real(8),dimension(6,(p+1)*(q+1)*(w+1)*nds)::Bglob,Bloc
    
    real(8),dimension(6,6)::matD
    real(8),dimension(ntens,1)::tempstr,deform
    
    real(8)::density
    real(8),dimension(nds,nds)::jac
    
    density = props(3)
    
    cpi = 0

    we = 1.0d0
    wn = 1.0d0
    wc = 1.0d0
    call gauleg(npi_xi, e, we)
    call gauleg(npi_eta, n, wn)
    call gauleg(npi_zeta, c, wc)
    
    !Gauss points cycle
    Ke = 0.0d0
    Finte = 0.0d0
    
    wer = 1.0d0
    wnr = 1.0d0
    wcr = 1.0d0
    call gauleg(npi_xi-1, er, wer)
    call gauleg(npi_eta-1, nr, wnr)
    call gauleg(npi_zeta-1, cr, wcr)
    
    cpi = 0
    do i=1,npi_xi-1
        do j=1,npi_eta-1
            do k=1,npi_zeta-1
            
                cpi = cpi + 1
                MatD(:,:) = TmatD(nel,cpi,:,:)
                
                xi = er(i)
                eta = nr(j)
                zeta = cr(k)
                
                call ShapeFunc(nel,xi,eta,zeta,R,dRdx,detj)
                
                !Weight factor
                gwt = wer(i)*wnr(j)*wcr(k)*detj
                
                Bglob = 0.0d0
                do k1=1,(p+1)*(q+1)*(w+1)
                    !Bglob(1,k1*3-2) = dRdx(k1,1)
                    !Bglob(2,k1*3-1) = dRdx(k1,2)
                    !Bglob(3,k1*3  ) = dRdx(k1,3)
                    
                    !Bglob(4,k1*3-2) = dRdx(k1,2)
                    !Bglob(4,k1*3-1) = dRdx(k1,1)
                    Bglob(5,k1*3-2) = dRdx(k1,3)
                    Bglob(5,k1*3  ) = dRdx(k1,1)
                    Bglob(6,k1*3-1) = dRdx(k1,3)
                    Bglob(6,k1*3  ) = dRdx(k1,2)
                end do
                
                Bloc = Bglob
                
!                call MatPlastic3D(iprops,nds,(p+1)*(q+1)*(w+1),props,Bloc,6,el_ddisp,stress(nel,cpi,:),strain(nel,cpi,:),dstrain(nel,cpi,:),dstress(nel,cpi,:),hard(nel,cpi),matD)
!                
!                tempstr = 0.0d0
!                do k1=1,ntens
!                    tempstr(k1,1) = stress(nel,cpi,k1)
!                    deform(k1,1) = dstrain(nel,cpi,k1)
!                end do
                
                Ke = Ke + matmul(matmul(transpose(Bglob),matD),Bglob)*gwt
                
!                Finte = Finte + matmul(transpose(Bglob),tempstr)*gwt
!                
!                !Gravitic loads -------------------------------------------------------------
!                if (gravity==.true.)then
!                    call grvt((p+1)*(q+1)*(w+1),nds,R,gwt,gconst,gravdir,density,Finte)
!                end if
!
!                !Strain Energy --------------------------------------------------------------
!                do k1=1,ntens
!                    SEnergy = SEnergy + tempstr(k1,1)*deform(k1,1)*gwt/2.0d0
!                end do
                
                continue
                
            end do  
        end do
    end do
    
    !Ke = 0.0d0
    !Finte = 0.0d0
    
    cpi = 0
    
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
                
                Bglob = 0.0d0
                do k1=1,(p+1)*(q+1)*(w+1)
                    Bglob(1,k1*3-2) = dRdx(k1,1)
                    Bglob(2,k1*3-1) = dRdx(k1,2)
                    Bglob(3,k1*3  ) = dRdx(k1,3)
                    
                    Bglob(4,k1*3-2) = dRdx(k1,2)
                    Bglob(4,k1*3-1) = dRdx(k1,1)
                    !Bglob(5,k1*3-2) = dRdx(k1,3)
                    !Bglob(5,k1*3  ) = dRdx(k1,1)
                    !Bglob(6,k1*3-1) = dRdx(k1,3)
                    !Bglob(6,k1*3  ) = dRdx(k1,2)
                end do
                
                Bloc = Bglob
                
                call MatPlastic3D(iprops,nds,(p+1)*(q+1)*(w+1),props,Bloc,6,el_ddisp,stress(nel,cpi,:),strain(nel,cpi,:),dstrain(nel,cpi,:),dstress(nel,cpi,:),hard(nel,cpi),matD)
                
                tempstr = 0.0d0
                do k1=1,ntens
                    tempstr(k1,1) = stress(nel,cpi,k1)
                    deform(k1,1) = dstrain(nel,cpi,k1)
                end do
                
                Ke = Ke + matmul(matmul(transpose(Bglob),matD),Bglob)*gwt
                
                Finte = Finte + matmul(transpose(Bglob),tempstr)*gwt
                
                !Gravitic loads -------------------------------------------------------------
                if (gravity==.true.)then
                    call grvt((p+1)*(q+1)*(w+1),nds,R,gwt,gconst,gravdir,density,Finte)
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

end subroutine