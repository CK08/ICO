subroutine STPAlloc(istp)
    
    use Mod_Variables
    implicit none
    integer(4),intent(IN)::istp
    integer(4)::i,j,k,count
    
    incmax = STP_incmax(istp)
    itermax = STP_itermax(istp)
    nbc = STP_nbc(istp)
    
    if(istp==1) deallocate(KfEq,KfEqInv,FextEq,dispeqv)
    
    allocate(KfEq(tnodes*nds-nbc,tnodes*nds-nbc))
    KfEq = 0.0d0
    
    allocate(KfEqInv(tnodes*nds-nbc,tnodes*nds-nbc))
    KfEqInv = 0.0d0
    
    allocate(FextEq(tnodes*nds-nbc,1))
    FextEq = 0.0d0
    
    allocate(dispeqv(tnodes*nds-nbc,1))
    dispeqv = 0.0d0
    
    dispdofold = dispdof
    
    do j=1,tnodes*nds
        dispdof(j) = STP_dispdof(istp,j) !- dispdofold(j)
    end do
    
    
    loaddof(:) = STP_loaddof(istp,:)
    BCdof(:) = STP_BCdof(istp,:)
    distload(:,:) = STP_distload(istp,:,:)
    
    nbc   = STP_nbc(istp)
    !nload = STP_nload(istp)
    ndisp = STP_ndisp(istp)
    
    sumImpDisp = 0.0d0
    do j=1,tnodes*nds
        sumImpDisp = sumImpDisp + STP_dispdof(istp,j)
    end do
    
end subroutine

subroutine STPDalloc(istp)
    
    use Mod_Variables
    implicit none
    integer(4),intent(IN)::istp
    integer(4)::i,j,k,count
    
    
    deallocate(KfEq)
    deallocate(KfEqInv)
    deallocate(FextEq)
    deallocate(dispeqv)
    

end subroutine 