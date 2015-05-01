!------------------------------------------------------------------------------------------------------------
!
! Q2/Q1 B-bar element
!
!------------------------------------------------------------------------------------------------------------

subroutine ElQuad9EBBar(nel,Ke,KBar,el_ddisp,Finte)
    
    use Mod_Variables
    implicit none
    
    integer(4)::nel
    integer(4)::cpi
    
    real(8),dimension((p+1)*(q+1)*nds,1),intent(IN)::el_ddisp
    real(8),dimension((p+1)*(q+1)*nds,(p+1)*(q+1)*nds),intent(OUT)::Ke
    real(8),dimension(p*q*nds,p*q*nds),intent(OUT)::KBar
    real(8),dimension((p+1)*(q+1)*nds,1),intent(OUT)::Finte
    
    
    integer(4),parameter::npi_xi= 3, npi_eta = 3
    
    integer(4)::i,j,ni,nj,nk,nel_nza,k1
    real(8)::xi,eta
    real(8),dimension(npi_xi)::e,we
    real(8),dimension(npi_eta)::n,wn
    
    real(8),dimension((p+1)*(q+1))::R
    real(8),dimension((p+1)*(q+1),nds)::dRdx
    real(8)::detJ,gwt
    
    real(8),dimension(3,(p+1)*(q+1)*nds)::Bglob,Bloc,BBar,BDev,Bdil
    
    real(8),dimension(3,3)::matD,MatDDev,MatDDil
    real(8),dimension(ntens,1)::tempstr,deform
    
    !B-Bar variables
    integer(4)::count,k2,k3,k4,k5,k6,error
    real(8)::vol
    real(8),dimension(npi_xi*npi_eta)::deter
    real(8),dimension(4,4)::MatM,Minv
    real(8),dimension(4)::temp1
    real(8),dimension(ncpx+p-1)::ubar
    real(8),dimension(ncpy+q-1)::vbar
    real(8),dimension((npi_xi-1)*(npi_eta-1))::Rbar,lixo
    real(8),dimension((npi_xi-1)*(npi_eta-1),1)::aux1
    real(8),dimension(4,(p+1)*(q+1)*nds)::SecInt
    
    real(8),dimension(3,1)::mm
    real(8),dimension(3,3)::I0,aux33
    real(8)::EG
    
    
!    IENbar(1,1) = 11
!    IENbar(1,2) = 9
!    IENbar(1,3) = 3
!    IENbar(1,4) = 1
!    
!    IENbar(2,1) = 12
!    IENbar(2,2) = 10
!    IENbar(2,3) = 4
!    IENbar(2,4) = 2
!    
!    IENbar(3,1) = 15
!    IENbar(3,2) = 13
!    IENbar(3,3) = 7
!    IENbar(3,4) = 5
!    
!    IENbar(4,1) = 16
!    IENbar(4,2) = 14
!    IENbar(4,3) = 8
!    IENbar(4,4) = 6
    
    
    cpi = 0
    
    MatDDev = 0.0d0
    MatDDil = 0.0d0
    
    mm = 1.0d0
    mm(3,1) = 0.0d0
    
    I0 = 0.0d0
    I0(1,1) = 1.0d0
    I0(2,2) = 1.0d0
    I0(3,3) = 0.5d0
    
    EG = props(1)/(2.0d0*(1.0d0+props(2)))
    
    aux33 = matmul(mm,transpose(mm))/3.0d0
    MatDDev = 2.0d0*EG*(I0-aux33) 
    
    we=1.0d0
    wn=1.0d0
    call gauleg(npi_xi, e, we)
    call gauleg(npi_eta, n, wn)
    
    do k1=2,ncpx+p
        ubar(k1-1) = u_knot(k1)
    end do
    
    do k1=2,ncpy+q
        vbar(k1-1) = v_knot(k1)
    end do 
    
    SecInt = 0.0d0
    
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Volume
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    count = 0
    vol = 0.0d0
    deter = 0.0d0
    do i=1,npi_xi
        do j=1,npi_eta
                
            count = count + 1
            
            xi = e(i)
            eta = n(j)
            
            call ShapeFunc(nel,xi,eta,R,dRdx,detj)
            
            vol = vol + detj*we(i)*wn(j)
            deter(count) = detj
                
        end do
    end do
    
    
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! M-MATRIX
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    matM = 0.0d0
    do i=1,npi_xi
        do j=1,npi_eta
                
            cpi = cpi + 1
            
            xi  = e(i)
            eta = n(j)
            
            
            !#####################################################################################################################################
            !#####################################################################################################################################
            call ShapeFunc(nel,xi,eta,R,dRdx,detj)
            call ShapeFuncBar(nel,xi,eta,ubar,vbar,RBar)   
            
!            aux1 = 0.0d0
!            do k1=1,(npi_xi-1)*(npi_eta-1)
!                aux1(k1,1) = Rbar(k1)
!            end do
            
            do k1=1,(npi_xi-1)*(npi_eta-1)      !index A
                do k2=1,(npi_xi-1)*(npi_eta-1)  !index B
                    
                    count = 0
                
                    do k3=1,npi_xi
                        do k4=1,npi_eta
                            count = count + 1
                            MatM(k1,k2) = MatM(k1,k2) + Rbar(k1)*Rbar(k2)*we(k3)*wn(k4)*deter(count)
                        end do
                    end do
                enddo
            end do
            
!            aux1 = 0.0d0
!            do k1=1,(npi_xi-1)*(npi_eta-1)
!                aux1(k1,1) = Rbar(k1)
!            end do
            
!            Minv = 0.0d0
!            Minv = MatM
!            call gaussj(Minv,(npi_xi-1)*(npi_eta-1),temp1,error)
            !#####################################################################################################################################
            !#####################################################################################################################################
            
            
            Bglob = 0.0d0
            do k1=1,(p+1)*(q+1)
                Bglob(1,k1*2-1) = dRdx(k1,1)
                Bglob(2,k1*2  ) = dRdx(k1,2)
                Bglob(3,k1*2-1) = dRdx(k1,2)
                Bglob(3,k1*2  ) = dRdx(k1,1)
            end do
            
            count = 0
            
            
            do k1=1,(p+1)*(q+1) !index i
                do k2=1,(npi_xi-1)*(npi_eta-1) !index B
                    count = 0
                    do k3=1,npi_xi
                        do k4=1,npi_eta
                            count = count + 1
                            SecInt(k2,k1*2-1) = SecInt(k2,k1*2-1) + Rbar(k2)*Bglob(1,k1*2-1)*we(k3)*wn(k4)*deter(count)
                            SecInt(k2,k1*2  ) = SecInt(k2,k1*2  ) + Rbar(k2)*Bglob(2,k1*2  )*we(k3)*wn(k4)*deter(count)
                        end do
                    end do
                end do
            end do
            
            
            
                
        end do
    end do
    
    
!    lixo=0.0d0
!    do k1=1,(npi_xi-1)*(npi_eta-1)
!        lixo(k1) = sum(matM(:,k1))
!    end do
!    
!    MatM = 0.0d0
!    do k1=1,(npi_xi-1)*(npi_eta-1)
!        matM(k1,k1) = lixo(k1)
!    end do
    
    Minv = 0.0d0
    Minv = MatM
    call gaussj(Minv,(npi_xi-1)*(npi_eta-1),temp1,error)
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    
    
    
    
    
    
    !Gauss points cycle
    Ke = 0.0d0
    KBar = 0.0d0
    Finte = 0.0d0
    cpi = 0
    do i=1,npi_xi
        do j=1,npi_eta
            
            cpi = cpi + 1
            MatD(:,:) = TmatD(nel,cpi,:,:)
            
            xi = e(i)
            eta = n(j)
            
            call ShapeFunc(nel,xi,eta,R,dRdx,detj)
            call ShapeFuncBar(nel,xi,eta,ubar,vbar,RBar)
            
            !Weight factor
            gwt = we(i)*wn(j)*detj
            
            Bglob = 0.0d0
            BDil  = 0.0d0
            BDev = 0.0d0
            do k1=1,(p+1)*(q+1)
                Bglob(1,k1*2-1) = dRdx(k1,1)
                Bglob(2,k1*2  ) = dRdx(k1,2)
                Bglob(3,k1*2-1) = dRdx(k1,2)
                Bglob(3,k1*2  ) = dRdx(k1,1)
                
                Bdil(1,k1*2-1) = dRdx(k1,1)/3.0d0
                Bdil(2,k1*2-1) = dRdx(k1,1)/3.0d0
                Bdil(1,k1*2  ) = dRdx(k1,2)/3.0d0
                Bdil(2,k1*2  ) = dRdx(k1,2)/3.0d0
            end do
            
            BBar = 0.0d0
            do k1=1,(p+1)*(q+1)
                count = 0
                do k2=1,(npi_xi-1)*(npi_eta-1) !index A
                    do k3=1,(npi_xi-1)*(npi_eta-1) !index B
                        
                        BBar(1,k1*2-1) = BBar(1,k1*2-1) + Rbar(k2)*MInv(k2,k3)*SecInt(k3,k1*2-1)
                        BBar(2,k1*2-1) = BBar(2,k1*2-1) + Rbar(k2)*MInv(k2,k3)*SecInt(k3,k1*2-1)

                        BBar(1,k1*2  ) = BBar(1,k1*2  ) + Rbar(k2)*MInv(k2,k3)*SecInt(k3,k1*2  )
                        BBar(2,k1*2  ) = BBar(2,k1*2  ) + Rbar(k2)*MInv(k2,k3)*SecInt(k3,k1*2  )

                    end do
                end do
            end do
            
!            Bbar = 0.0d0
!            do k1=1,(p+1)*(q+1)
!                do k2=1,(npi_xi-1)*(npi_eta-1)           !A
!                    do k3 = 1,(npi_xi-1)*(npi_eta-1)     !B
!                        
!                        count = 0
!                        
!                        do k4=1,npi_xi                   !intx
!                            do k5=1,npi_eta              !inty
!                                
!                                count = count + 1
!                                
!                                BBar(1,k1*2-1) = BBar(1,k1*2-1) + Rbar(k2)*Minv(k2,k3)*Rbar(k3)*we(k4)*wn(k5)*deter(count)*BGlob(1,k1*2-1)
!                                BBar(2,k1*2-1) = BBar(2,k1*2-1) + Rbar(k2)*Minv(k2,k3)*Rbar(k3)*we(k4)*wn(k5)*deter(count)*BGlob(1,k1*2-1)
!
!                                BBar(1,k1*2  ) = BBar(1,k1*2  ) + Rbar(k2)*Minv(k2,k3)*Rbar(k3)*we(k4)*wn(k5)*deter(count)*BGlob(2,k1*2  )
!                                BBar(2,k1*2  ) = BBar(2,k1*2  ) + Rbar(k2)*Minv(k2,k3)*Rbar(k3)*we(k4)*wn(k5)*deter(count)*BGlob(2,k1*2  )
!                            end do
!                        end do
!                    end do
!                end do
!            end do
            
            
            BBar  = Bbar/3.0d0
            Bdev  = Bglob - Bdil
            Bglob = Bdev + BBar
            
            
            
            Bloc = Bglob
            
            !call MaterialPStrain(ndi,nshr,ntens,iprops,nodedof,nnodes,props,B,Blines,ddisp,stress,strain,dstrain,dstress,hard,matD,matf0,matf,dtempo)
            call MaterialPStrain(iprops,nds,(p+1)*(q+1),props,Bloc,3,el_ddisp,stress(nel,cpi,:),&
                                strain(nel,cpi,:),dstrain(nel,cpi,:),dstress(nel,cpi,:),hard(nel,cpi),matD)
            !call MaterialPStress(iprops,nds,(p+1)*(q+1),props,Bloc,3,el_ddisp,stress(nel,cpi,:),strain(nel,cpi,:),dstrain(nel,cpi,:),dstress(nel,cpi,:),hard(nel,cpi),matD)
            
            tempstr = 0.0d0
            do k1=1,ntens
                tempstr(k1,1) = stress(nel,cpi,k1)
                deform(k1,1) = dstrain(nel,cpi,k1)
            end do
            
            Ke = Ke + matmul(matmul(transpose(Bglob),matD),Bglob)*gwt
            
            !Ke = Ke + matmul(matmul(transpose(Bglob),matDDev),Bglob)*gwt
            
            !Kbar = 
            
            Finte = Finte + matmul(transpose(Bloc),tempstr)*gwt
            
            do k1=1,ntens
                SEnergy = SEnergy + tempstr(k1,1)*deform(k1,1)*gwt/2.0d0
            end do
            
            continue
          
        end do
    end do

    continue

end subroutine ElQuad9EBBar
