subroutine KLM_NoContact(i,iLM,KLM,FLM)
    
    use Mod_Variables
    implicit none
    
    real(8),dimension(tnodes*nds+ngrv,tnodes*nds+ngrv),intent(INOUT)::KLM
    real(8),dimension(tnodes*nds,ngrv),intent(INOUT)::FLM
    integer(4),intent(INOUT)::iLM 
    integer(4),intent(IN)::i
    integer(4),dimension(ncp_s+ncp_m)::smconn
    integer(4)::j,k1,count
    
    iLM = iLM + 1
    KLM(tnodes*nds+i,:) = 0.0d0
    KLM(:,tnodes*nds+i) = 0.0d0
    KLM(tnodes*nds+i,tnodes*nds+i) = 1.0d0
!    FLM(:,i) = 0.0d0
    
    do j=1,ncp_s
        smconn(j)=conn_slave(j)
    end do
    
    count = 0
    
    do j=ncp_s+1,ncp_m+ncp_s
        count = count + 1
        smconn(j)=conn_master(count)
    end do
            
    do k1=1,ncp_s+ncp_m
        FintLM(smconn(k1)*nds-1,1) = 0.0d0
        FintLM(smconn(k1)*nds  ,1) = 0.0d0
    end do
    
end subroutine
