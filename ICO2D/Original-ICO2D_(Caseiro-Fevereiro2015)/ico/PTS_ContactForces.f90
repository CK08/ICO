subroutine PTS_ContactForces(iter,inc)
    
    use Mod_Variables
    implicit none
    integer(4)::i,j,k,k1,k2,kspan,ii
    integer(4),intent(IN)::iter,inc
    integer(4),dimension(ncp_m+1)::ctconn
    real(8),dimension(ncp_s)::norm,Lagrangex,Lagrangey
    real(8),dimension(p_s+1)::Rss,dRss,ddRss
    real(8),dimension(ncp_s)::Rs,dRsdxi,d2Rsdxi2
    real(8),dimension(ncp_s,ncp_s)::MCC
    real(8),dimension(ncp_s)::temp1
    real(8),dimension(ncp_s,1)::array
    
!    Lagrangex = 0.0d0
!    Lagrangey = 0.0d0
    
    MCC = 0.0d0
    array = 0.0d0
    
    FintLM = 0.0d0
    do i=1,ngrv
        
        !Compute B-Splines basis functions...
         call BSplineBasisAndDeriv2(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan)
         !... and convert to NURBS basis, including vanishing terms 
         call BSplinesToNURBS(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
        
        do k1=1,ncp_s
            MCC(k1,i) = Rs(k1)
        end do
        
        !if (PTS_active(i)==1 .and. iter .gt. 1)then
        !if (PTS_active(i)==1)then
            ctconn = 0
            ctconn(1) = conn_slave(i)
            do j=1,ncp_m
                ctconn(j+1)=conn_master(j)
            end do
            
            do k1=1,ncp_s
                FintLM(conn_slave(k1)*nds-1,1) = Fini(conn_slave(k1)*nds-1,1)/(1.0d0*incmax)*(1.0d0*inc) - Fint(conn_slave(k1)*nds-1,1)
                FintLM(conn_slave(k1)*nds  ,1) = Fini(conn_slave(k1)*nds  ,1)/(1.0d0*incmax)*(1.0d0*inc) - Fint(conn_slave(k1)*nds  ,1)
            end do
            
            do k1=2,ncp_m+1
                FintLM(ctconn(k1)*nds-1,1) = Fini(ctconn(k1)*nds-1,1)/(1.0d0*incmax)*(1.0d0*inc) - Fint(ctconn(k1)*nds-1,1)
                FintLM(ctconn(k1)*nds  ,1) = Fini(ctconn(k1)*nds  ,1)/(1.0d0*incmax)*(1.0d0*inc) - Fint(ctconn(k1)*nds  ,1)
            end do
            
!            do k1=1,ncp_s
!                Lagrangex(i) = Lagrangex(i) + Rs(k1)*FintLM(conn_slave(k1)*nds-1,1)
!                Lagrangey(i) = Lagrangey(i) + Rs(k1)*FintLM(conn_slave(k1)*nds  ,1)
!            end do
!            
!            Lagrange(i) = -1.0d0*dsqrt(Lagrangex(i)*Lagrangex(i)+Lagrangey(i)*Lagrangey(i))
!            
!            array(i,1) = Lagrange(i)
            
        !end if
    end do
    
!    PTS_Stress = 0.0d0
!    do i=1,ngrv
!        !PTS_Stress(i) = PTS_Stress(i) + dLagrange(i)/CP_area(i)
!        PTS_Stress(i) = Lagrange(i)/CP_area(i)
!    end do
    
!    FintLM = 0.0d0
!    if (iter .gt. 1)then
!        FintLM = Fini/(1.0d0*incmax)*(1.0d0*inc) - Fint
!    end if
    
!    call gaussj(MCC,ncp_s,temp1,ii)
!    
!    array = matmul(MCC,array)
!    
!    do k1=1,ncp_s
!        if (PTS_active(k1)==1 .and. iter .gt. 1) Lagrange(k1) = array(k1,1) 
!    end do
    
    continue
    
end subroutine