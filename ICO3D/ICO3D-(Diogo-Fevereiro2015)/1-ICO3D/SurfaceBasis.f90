!----------------------------------------------------------------------------------------------
!
! Subroutine to compute the surface NURBS basis functions for the contact analysis
!
!----------------------------------------------------------------------------------------------

subroutine SurfaceBasis(ncpx,p,xi,u_knot,ncpy,q,eta,v_knot,Gcoords,conn,tnodes,nds,R,dR,ddR)
    
    implicit none
    
    integer(4),intent(IN)::p,q
    integer(4),intent(IN)::ncpx,ncpy
    integer(4),intent(IN)::nds,tnodes
    integer(4),dimension(ncpx*ncpy),intent(IN)::conn
    real(8),intent(IN)::xi,eta
    real(8),dimension(ncpx+p+1),intent(IN)::u_knot
    real(8),dimension(ncpy+q+1),intent(IN)::v_knot
    real(8),dimension(tnodes,nds+1),intent(IN)::GCoords
    
    real(8),dimension(ncpx*ncpy),intent(OUT)::R
    real(8),dimension(ncpx*ncpy,2),intent(OUT)::dR
    real(8),dimension(ncpx*ncpy,4),intent(OUT)::ddR
    
    integer(4)::ispan,jspan,k1,k2,k3,k4,count,count2,k
    real(8)::sumt,sumxi,sumeta,sumxixi,sumetaeta,sumxieta,sumetaxi
    
    integer(4),dimension(ncpx)::iN
    integer(4),dimension(ncpy)::iM
    integer(4),dimension(ncpx*ncpy)::iNM
    
    
    real(8),dimension(p+1)::N,dN,ddN
    real(8),dimension(q+1)::M,dM,ddM
    real(8),dimension((p+1)*(q+1))::R1
    real(8),dimension((p+1)*(q+1),2)::dR1
    real(8),dimension((p+1)*(q+1),4)::ddR1
    
    !------------------------------------------------------------------------------
    ! B-Spline basis functions...
    !------------------------------------------------------------------------------
    call BSplineBasisAndDeriv2(ncpx,p,xi,u_knot,N,dN,ddN,ispan)
    call BSplineBasisAndDeriv2(ncpy,q,eta,v_knot,M,dM,ddM,jspan)
    
    !------------------------------------------------------------------------------
    ! Auxiliary arrays
    !------------------------------------------------------------------------------
    iN = 0
    count = 0
    do k=ispan-p,ispan
        count = count + 1
        iN(k) = 1
    end do
    
    iM = 0
    do k=jspan-q,jspan
        count = count + 1
        iM(k) = 1
    end do
    
    count = 0
    do k2 = 1,ncpy
        do k1 = 1,ncpx
            count = count + 1
            iNM(count) = iN(k1)*iM(k2)
        end do
    end do
    
    !------------------------------------------------------------------------------
    ! Convert to NURBS basis
    !------------------------------------------------------------------------------
    R1 = 0.0d0
    dR1 = 0.0d0
    ddR1 = 0.0d0
    sumt = 0.0d0
    sumxi = 0.0d0
    sumxixi = 0.0d0
    sumeta = 0.0d0
    sumetaeta = 0.0d0
    sumxieta = 0.0d0
    sumetaxi = 0.0d0
    count = 0
    
    do k2 = 1,q+1
        do k1 = 1,p+1
            count = count + 1
            
            R1(count) = N(k1)*M(k2)*GCoords(conn(count),nds+1)
            sumt = sumt + R1(count)
            
            dR1(count,1) = dN(k1)*M(k2)*GCoords(conn(count),nds+1)
            sumxi = sumxi + dR1(count,1)
            
            dR1(count,2) = N(k1)*dM(k2)*GCoords(conn(count),nds+1)
            sumeta = sumeta + dR1(count,2)
            
            if(p .gt. 1) then
                !dR_dxi_dxi
                ddR1(count,1) = dN(k1)*M(k2)*GCoords(conn(count),nds+1)
                sumxixi = sumxixi + ddR1(count,1)
            end if
            if(p .gt. 1 .and. q .gt. 1) then    
                !dR_dxi_deta
                ddR1(count,2) = dN(k1)*dM(k2)*GCoords(conn(count),nds+1)
                sumxieta = sumxieta + ddR1(count,2)
                
                !dR_deta_dxi
                ddR1(count,3) = ddR1(count,2)
                sumetaxi = sumxieta
            end if
            if(q .gt. 1) then    
                !dR_deta_deta
                ddR1(count,4) = N(k1)*dM(k2)*GCoords(conn(count),nds+1)
                sumetaeta = sumetaeta + ddR1(count,4)
            end if
            
        end do
    end do
    
    do k1=1,(p+1)*(q+1)
        
        R1(k1) = R1(k1)/sumt
        
        dR1(k1,1) = (dR1(k1,1) - R1(k1)*sumxi)/sumt
        dR1(k1,2) = (dR1(k1,2) - R1(k1)*sumeta)/sumt
        
        !NURBS book pg. 137 ---
        if(p .gt. 1) then
            ddR1(k1,1) = (ddR1(k1,1) - 2.0d0*sumxi*dR1(k1,1) - sumxixi*R1(k1))/sumt
        end if
        if(p .gt. 1 .and. q .gt. 1) then
            ddR1(k1,2) = (ddR1(k1,2) - sumeta*dR1(k1,1) - sumxi*dR1(k1,2) - sumxieta*R1(k1))/sumt
            ddR1(k1,3) = ddR1(k1,2)
        end if
        if(q .gt. 1) then  
            ddR1(k1,4) = (ddR1(k1,4) - 2.0d0*sumeta*dR1(k1,2) - sumetaeta*R1(k1))/sumt
        end if
    end do
    
    !------------------------------------------------------------------------------
    ! Include vanishing terms in the basis
    !------------------------------------------------------------------------------
    R = 0.0d0
    dR = 0.0d0
    ddR = 0.0d0
    count = 0
    count2 = 0
    
    do k2=1,ncpy
        do k1=1,ncpx
            
            count = count + 1
            
            if(iNM(count) == 1)then
                
                count2 = count2 + 1
                
                R(count) = R1(count2)
                
                dR(count,1) = dR1(count2,1)
                dR(count,2) = dR1(count2,2)
                
                if(p .gt. 1)then
                    ddR(count,1) = ddR1(count2,1)
                end if
                if(p .gt. 1 .and. q .gt. 1) then
                    ddR(count,2) = ddR1(count2,2)
                    ddR(count,3) = ddR1(count2,3)
                end if
                if(q .gt. 1) then  
                    ddR(count,4) = ddR1(count2,4)
                end if
            end if  
        end do   
    end do

    
    10 continue

end subroutine

!subroutine SurfaceBasis(ncpx,p,xi,u_knot,ncpy,q,eta,v_knot,Gcoords,conn,tnodes,nds,R,dR,ddR)
!    
!    implicit none
!    
!    integer(4),intent(IN)::p,q
!    integer(4),intent(IN)::ncpx,ncpy
!    integer(4),intent(IN)::nds,tnodes
!    integer(4),dimension(ncpx*ncpy),intent(IN)::conn
!    real(8),intent(IN)::xi,eta
!    real(8),dimension(ncpx+p+1),intent(IN)::u_knot
!    real(8),dimension(ncpy+q+1),intent(IN)::v_knot
!    real(8),dimension(tnodes,nds+1),intent(IN)::GCoords
!    
!    real(8),dimension(ncpx*ncpy),intent(OUT)::R
!    real(8),dimension(ncpx*ncpy,2),intent(OUT)::dR
!    real(8),dimension(ncpx*ncpy,4),intent(OUT)::ddR
!    
!    integer(4)::ispan,jspan,k1,k2,k3,k4,count,count2
!    real(8)::sumt,sumxi,sumeta,sumxixi,sumetaeta,sumxieta,sumetaxi
!    
!    real(8),dimension(p+1)::N,dN,ddN
!    real(8),dimension(q+1)::M,dM,ddM
!    real(8),dimension((p+1)*(q+1))::R1
!    real(8),dimension((p+1)*(q+1),2)::dR1
!    real(8),dimension((p+1)*(q+1),4)::ddR1
!    
!    !------------------------------------------------------------------------------
!    ! B-Spline basis functions...
!    !------------------------------------------------------------------------------
!    call BSplineBasisAndDeriv2(ncpx,p,xi,u_knot,N,dN,ddN,ispan)
!    call BSplineBasisAndDeriv2(ncpy,q,eta,v_knot,M,dM,ddM,jspan)
!    
!    !------------------------------------------------------------------------------
!    ! ... convert to NURBS basis...
!    !------------------------------------------------------------------------------
!    R1 = 0.0d0
!    dR1 = 0.0d0
!    ddR1 = 0.0d0
!    sumt = 0.0d0
!    sumxi = 0.0d0
!    sumxixi = 0.0d0
!    sumeta = 0.0d0
!    sumetaeta = 0.0d0
!    sumxieta = 0.0d0
!    sumetaxi = 0.0d0
!    count = 0
!    
!    do k2 = 1,q+1
!        do k1 = 1,p+1
!            count = count + 1
!            
!            R1(count) = N(k1)*M(k2)*GCoords(conn(count),nds+1)
!            sumt = sumt + R1(count)
!            
!            dR1(count,1) = dN(k1)*M(k2)*GCoords(conn(count),nds+1)
!            sumxi = sumxi + dR1(count,1)
!            
!            dR1(count,2) = N(k1)*dM(k2)*GCoords(conn(count),nds+1)
!            sumeta = sumeta + dR1(count,2)
!            
!            if(p .gt. 1) then
!                !dR_dxi_dxi
!                ddR1(count,1) = dN(k1)*M(k2)*GCoords(conn(count),nds+1)
!                sumxixi = sumxixi + ddR1(count,1)
!            end if
!            if(p .gt. 1 .and. q .gt. 1) then    
!                !dR_dxi_deta
!                ddR1(count,2) = dN(k1)*dM(k2)*GCoords(conn(count),nds+1)
!                sumxieta = sumxieta + ddR1(count,2)
!                
!                !dR_deta_dxi
!                ddR1(count,3) = ddR1(count,2)
!                sumetaxi = sumxieta
!            end if
!            if(q .gt. 1) then    
!                !dR_deta_deta
!                ddR1(count,4) = N(k1)*dM(k2)*GCoords(conn(count),nds+1)
!                sumetaeta = sumetaeta + ddR1(count,4)
!            end if
!            
!        end do
!    end do
!    
!    do k1=1,(p+1)*(q+1)
!        
!        R1(k1) = R1(k1)/sumt
!        
!        dR1(k1,1) = (dR1(k1,1) - R1(k1)*sumxi)/sumt
!        dR1(k1,2) = (dR1(k1,2) - R1(k1)*sumeta)/sumt
!        
!        !NURBS book pg. 137 ---
!        if(p .gt. 1) then
!            ddR1(k1,1) = (ddR1(k1,1) - 2.0d0*sumxi*dR1(k1,1) - sumxixi*R1(k1))/sumt
!        end if
!        if(p .gt. 1 .and. q .gt. 1) then
!            ddR1(k1,2) = (ddR1(k1,2) - sumeta*dR1(k1,1) - sumxi*dR1(k1,2) - sumxieta*R1(k1))/sumt
!            ddR1(k1,3) = ddR1(k1,2)
!        end if
!        if(q .gt. 1) then  
!            ddR1(k1,4) = (ddR1(k1,4) - 2.0d0*sumeta*dR1(k1,2) - sumetaeta*R1(k1))/sumt
!        end if
!    end do
!    
!    !------------------------------------------------------------------------------
!    ! ... and include vanishing terms
!    !------------------------------------------------------------------------------
!    R = 0.0d0
!    dR = 0.0d0
!    ddR = 0.0d0
!    count = 0
!    count2 = 0
!    
!    do k2=1,ncpy
!        do k1=1,ncpx
!        
!            count = count + 1
!        
!            if(k1 == ispan-p .and. k2 == jspan-q)then
!                do k4=k2,jspan
!                    do k3=k1,ispan
!                    
!                        count2 = count2 + 1
!                            
!                        R(count) = R1(count2)
!                
!                        dR(count,1) = dR1(count2,1)
!                        dR(count,2) = dR1(count2,2)
!                        
!                        if(p .gt. 1)then
!                            ddR(count,1) = ddR1(count2,1)
!                        end if
!                        if(p .gt. 1 .and. q .gt. 1) then
!                            ddR(count,2) = ddR1(count2,2)
!                            ddR(count,3) = ddR1(count2,3)
!                        end if
!                        if(q .gt. 1) then  
!                            ddR(count,4) = ddR1(count2,4)
!                        end if
!                        
!                        count = count + 1
!                        
!                    end do
!                end do
!                
!                goto 10
!                
!            end if
!            
!            
!            
!        end do   
!    end do
!
!    
!    10 continue
!
!end subroutine

!subroutine SurfaceBasis(ncpx,p,xi,u_knot,ncpy,q,eta,v_knot,Gcoords,conn,tnodes,nds,R,dR,ddR)
!    
!    implicit none
!    
!    integer(4),intent(IN)::p,q
!    integer(4),intent(IN)::ncpx,ncpy
!    integer(4),intent(IN)::nds,tnodes
!    integer(4),dimension(ncpx*ncpy),intent(IN)::conn
!    real(8),intent(IN)::xi,eta
!    real(8),dimension(ncpx+p+1),intent(IN)::u_knot
!    real(8),dimension(ncpy+q+1),intent(IN)::v_knot
!    real(8),dimension(tnodes,nds+1),intent(IN)::GCoords
!    
!    real(8),dimension(ncpx*ncpy),intent(OUT)::R
!    real(8),dimension(ncpx*ncpy,2),intent(OUT)::dR
!    real(8),dimension(ncpx*ncpy,4),intent(OUT)::ddR
!    
!    integer(4)::ispan,jspan,k,k1,k2,k3,k4,count,count2
!    real(8)::sumt,sumxi,sumeta,sumxixi,sumetaeta,sumxieta,sumetaxi
!    
!    real(8),dimension(p+1)::N1,dN1,ddN1
!    real(8),dimension(q+1)::M1,dM1,ddM1
!    real(8),dimension(ncpx)::N,dN,ddN
!    real(8),dimension(ncpy)::M,dM,ddM
!    real(8),dimension((p+1)*(q+1))::R1
!    real(8),dimension((p+1)*(q+1),2)::dR1
!    real(8),dimension((p+1)*(q+1),4)::ddR1
!    
!    !------------------------------------------------------------------------------
!    ! B-Spline basis functions...
!    !------------------------------------------------------------------------------
!    call BSplineBasisAndDeriv2(ncpx,p,xi,u_knot,N1,dN1,ddN1,ispan)
!    call BSplineBasisAndDeriv2(ncpy,q,eta,v_knot,M1,dM1,ddM1,jspan)
!    
!    
!    !------------------------------------------------------------------------------
!    ! ... include vanishing terms...
!    !------------------------------------------------------------------------------
!    N = 0.0d0
!    dN = 0.0d0
!    ddN = 0.0d0
!    count = 0
!    do k=ispan-p,ispan
!        count = count + 1
!        if(p .gt. 1) then
!            ddN(k) = ddN1(count)
!        end if
!        
!        dN(k) = dN1(count)
!        N(k) = N1(count)
!        
!    end do
!    
!    M = 0.0d0
!    dM = 0.0d0
!    ddM = 0.0d0
!    count = 0
!    do k=jspan-q,jspan
!        count = count + 1
!        if(q .gt. 1) then
!            ddM(k) = ddM1(count)
!        end if
!        
!        dM(k) = dM1(count)
!        M(k) = M1(count)
!        
!    end do
!    
!    !------------------------------------------------------------------------------
!    ! ... convert to NURBS basis...
!    !------------------------------------------------------------------------------
!    R1 = 0.0d0
!    dR1 = 0.0d0
!    ddR1 = 0.0d0
!    sumt = 0.0d0
!    sumxi = 0.0d0
!    sumxixi = 0.0d0
!    sumeta = 0.0d0
!    sumetaeta = 0.0d0
!    sumxieta = 0.0d0
!    sumetaxi = 0.0d0
!    count = 0
!    
!    do k2 = jspan-q,jspan
!        do k1 = ispan-p,ispan
!            count = count + 1
!            
!            R(count) = N(k1)*M(k2)*GCoords(conn(count),nds+1)
!            sumt = sumt + R1(count)
!            
!            dR(count,1) = dN(k1)*M(k2)*GCoords(conn(count),nds+1)
!            sumxi = sumxi + dR1(count,1)
!            
!            dR(count,2) = N(k1)*dM(k2)*GCoords(conn(count),nds+1)
!            sumeta = sumeta + dR1(count,2)
!            
!            if(p .gt. 1) then
!                !dR_dxi_dxi
!                ddR(count,1) = dN(k1)*M(k2)*GCoords(conn(count),nds+1)
!                sumxixi = sumxixi + ddR1(count,1)
!            end if
!            if(p .gt. 1 .and. q .gt. 1) then    
!                !dR_dxi_deta
!                ddR(count,2) = dN(k1)*dM(k2)*GCoords(conn(count),nds+1)
!                sumxieta = sumxieta + ddR1(count,2)
!                
!                !dR_deta_dxi
!                ddR(count,3) = ddR1(count,2)
!                sumetaxi = sumxieta
!            end if
!            if(q .gt. 1) then    
!                !dR_deta_deta
!                ddR(count,4) = N(k1)*dM(k2)*GCoords(conn(count),nds+1)
!                sumetaeta = sumetaeta + ddR1(count,4)
!            end if
!            
!        end do
!    end do
!    
!    do k1=1,ncpx*ncpy
!        
!        R(k1) = R(k1)/sumt
!        
!        dR(k1,1) = (dR(k1,1) - R(k1)*sumxi)/sumt
!        dR(k1,2) = (dR(k1,2) - R(k1)*sumeta)/sumt
!        
!        !NURBS book pg. 137 ---
!        if(p .gt. 1) then
!            ddR(k1,1) = (ddR(k1,1) - 2.0d0*sumxi*dR(k1,1) - sumxixi*R(k1))/sumt
!        end if
!        if(p .gt. 1 .and. q .gt. 1) then
!            ddR(k1,2) = (ddR(k1,2) - sumeta*dR(k1,1) - sumxi*dR(k1,2) - sumxieta*R(k1))/sumt
!            ddR(k1,3) = ddR(k1,2)
!        end if
!        if(q .gt. 1) then  
!            ddR(k1,4) = (ddR(k1,4) - 2.0d0*sumeta*dR(k1,2) - sumetaeta*R(k1))/sumt
!        end if
!    end do
!    
!!    !------------------------------------------------------------------------------
!!    ! ... and include vanishing terms
!!    !------------------------------------------------------------------------------
!!    R = 0.0d0
!!    dR = 0.0d0
!!    ddR = 0.0d0
!!    count = 0
!!    count2 = 0
!!    
!!    do k2=1,ncpy
!!        do k1=1,ncpx
!!        
!!            count = count + 1
!!        
!!            if(k1 == ispan-p .and. k2 == jspan-q)then
!!                do k4=k2,jspan
!!                    do k3=k1,ispan
!!                    
!!                        count2 = count2 + 1
!!                            
!!                        R(count) = R1(count2)
!!                
!!                        dR(count,1) = dR1(count2,1)
!!                        dR(count,2) = dR1(count2,2)
!!                        
!!                        if(p .gt. 1)then
!!                            ddR(count,1) = ddR1(count2,1)
!!                        end if
!!                        if(p .gt. 1 .and. q .gt. 1) then
!!                            ddR(count,2) = ddR1(count2,2)
!!                            ddR(count,3) = ddR1(count2,3)
!!                        end if
!!                        if(q .gt. 1) then  
!!                            ddR(count,4) = ddR1(count2,4)
!!                        end if
!!                        
!!                        count = count + 1
!!                        
!!                    end do
!!                end do
!!                
!!                goto 10
!!                
!!            end if
!!            
!!            
!!            
!!        end do   
!!    end do
!
!    
!    10 continue
!
!end subroutine
