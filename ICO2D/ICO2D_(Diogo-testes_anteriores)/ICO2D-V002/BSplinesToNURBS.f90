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
        
        dRmdxi(k) = (dRmdxi(k)*sumt - Rm(k)*sumd)/(sumt*sumt) ! Normalizacao das variais ???
        
        Rm(k) = Rm(k)/sumt
    end do
    
    continue

end
