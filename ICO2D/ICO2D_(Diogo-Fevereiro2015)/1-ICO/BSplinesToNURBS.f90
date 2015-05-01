!-------------------------------------------------------------------------------------------------------
!
! Convert B-Spline Basis to NURBS
! (Used in contact mechanics)
!
!-------------------------------------------------------------------------------------------------------

subroutine BSplinesToNURBS(ncp,ps,xib,knot_vec,Rmm,dRmm,ddRmm,kspan,connect,Rm,dRmdxi,d2Rmdxi2)

    use mod_variables
    implicit none
    
    integer(4)::ncp,ps,kspan
    real(8)::xib
    real(8),dimension(ncp+ps+1)::knot_vec
    real(8),dimension(ps+1)::Rmm,dRmm,ddRmm
    real(8),dimension(ncp)::Rm,dRmdxi,d2Rmdxi2
    integer(4),dimension(ncp)::connect
    integer(4)::k,count
    real(8)::sumt,sumd,sumdd
    
    Rm       = 0.0d0
    dRmdxi   = 0.0d0
    d2Rmdxi2 = 0.0d0
    count = 0
    sumt = 0.0d0
    sumd = 0.0d0
    sumdd = 0.0d0
    do k=kspan-ps,kspan-ps+ps
        count = count + 1
        Rm(k)       =   Rmm(count)*GCoords(connect(k),3)
        sumt = sumt + Rm(k)
        
        dRmdxi(k)   =  dRmm(count)*GCoords(connect(k),3)
        sumd = sumd + dRmdxi(k)
        
        if(ps .gt. 1) then
            d2Rmdxi2(k) = ddRmm(count)*GCoords(connect(k),3)
            sumdd = sumdd + d2Rmdxi2(k)
        end if
    end do
    
    do k=kspan-ps,kspan-ps+ps
        count = count + 1
        if(ps .gt. 1) then
            d2Rmdxi2(k) = (d2Rmdxi2(k)*sumt - dRmdxi(k)*sumdd)/(sumt*sumt)
        end if
        
        dRmdxi(k) = (dRmdxi(k)*sumt - Rm(k)*sumd)/(sumt*sumt)
        
        Rm(k) = Rm(k)/sumt
    end do
    
    continue

end






subroutine BSplinesToNURBS2(ncps,ps,xib,knot_vec,Rmm,dRmm,ddRmm,kspan,connect,R,dRdxi,ddRdxi2)

    use mod_variables
    implicit none
    
    integer(4)::ncps,ps,kspan
    real(8)::xib
    real(8),dimension(ncps+ps+1)::knot_vec
    real(8),dimension(ps+1)::Rmm,dRmm,ddRmm
    real(8),dimension(ncps)::R,dRdxi,ddRdxi2
    integer(4),dimension(ncps)::connect
    integer(4)::i,j,k,count
    real(8)::sumtot,sumxi,sumxi2
    
    R = 0.0d0
    dRdxi = 0.0d0
    ddRdxi2 = 0.0d0
    
    sumtot = 0.0d0
    sumxi = 0.0d0
    sumxi2 = 0.0d0
    
    count = 0
    do i=kspan-ps,kspan
        count = count + 1

        R(i) = Rmm(count)*GCoords(Connect(i),3)
        sumtot = sumtot + R(i)
        
        dRdxi(i) = dRmm(count)*GCoords(Connect(i),3)
        sumxi = sumxi + dRdxi(i)
        
        if(ps .gt. 1)then
            ddRdxi2(i) = ddRmm(count)*GCoords(Connect(i),3)
            sumxi2 = sumxi2 + dRdxi(i)
        end if
            
            
    end do
    
    count = 0
    do k=kspan-ps,kspan
        count = count + 1
        if(ps .gt. 1) then
            ddRdxi2(k) = (ddRdxi2(k)*sumtot - dRdxi(k)*sumxi2)/(sumtot*sumtot)
        end if
        
        dRdxi(k) = (dRdxi(k)*sumtot - R(k)*sumxi)/(sumtot*sumtot)
        
        R(k) = R(k)/sumtot
    end do
    
    continue
    
    
    

end