!----------------------------------------------------------------------------------------------
!
! Subroutine determine the residual and force error
!
!----------------------------------------------------------------------------------------------

subroutine ConvCheck (inc,AbsRes,AbsDisp)

    use Mod_Variables
    implicit none
    
    integer(4),intent(IN)::inc
    real(8),intent(OUT)::AbsRes,AbsDisp
    
    integer(4)::i,j,k
    real(8)::sumRes, sumFext, sumd, sumdt
    
    absRes=0.0d0
    sumRes=0.0d0
    sumFext=0.0d0
    sumd=0.0d0
    sumdt=0.0d0
    if(nptch == 1)then
        do i=1,ncpx*ncpy*nds
            sumRes = sumRes + Res(i,1)**2.0d0
            sumFext = sumFext + (Fini(i,1)/(1.0d0*incmax)*(1.0d0*inc))**2.0d0

            sumd= sumd + dddisp(i,1)**2.0d0
            sumdt= sumdt + ddisp(i,1)**2.0d0
        end do
    else
        do i=1,tnodes*nds
            sumRes = sumRes + Res(i,1)**2.0d0
            sumFext = sumFext + (Fini(i,1)/(1.0d0*incmax)*(1.0d0*inc))**2.0d0

            sumd= sumd + dddisp(i,1)**2.0d0
            sumdt= sumdt + ddisp(i,1)**2.0d0
        end do
    end if
    
    if(sumFext==0.0d0)then
        !absres=0.0d0
        sumFext =1.0d0
        if(sumdt .lt. 1.0d-10 .and. sumImpDisp == 0.0d0)then
            absdisp = 999.0d0
        elseif(sumdt .lt. 1.0d-10) then
            absdisp=0.0d0
!                    elseif(sumdt==0.0d0 .and. sumd == 0.0d0) then
!                        absdisp=0.0d0
        else 
            absRes=sqrt(sumRes)/sqrt(sumFext)*100.0d0
            absdisp=sqrt(sumd)/sqrt(sumdt)*100.0d0
        end if
    else
        absRes=sqrt(sumRes)/sqrt(sumFext)*100.0d0
        absdisp=sqrt(sumd)/sqrt(sumdt)*100.0d0
    end if

end subroutine