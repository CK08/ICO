SUBROUTINE gaussCoord(R,disp,cpi,nelp,npi_xi,npi_eta,npi_zeta,nodes)

use mod_variables
IMPLICIT NONE

integer(4), intent(in)::cpi, nelp, npi_xi, npi_eta, npi_zeta
real(8), dimension((p+1)*(q+1)*(w+1)), intent(in):: R
real(8), dimension((p+1)*(q+1)*(w+1)*nds,1), intent(in):: disp

real(8),dimension((p+1)*(q+1)*(w+1),nds), intent(in)::nodes
integer(4)::k1,ipos


open(unit=10,file='GPCoords.txt', access = 'APPEND')

ipos = npi*(nelp-1)+cpi
do k1=1,(p+1)*(q+1)*(w+1)
        GP_coords(ipos,1) = GP_coords(ipos,1) + R(k1)*(nodes(k1,1)+disp(k1*nds-2,1))
        GP_coords(ipos,2) = GP_coords(ipos,2) + R(k1)*(nodes(k1,2)+disp(k1*nds-1,1))
        GP_coords(ipos,3) = GP_coords(ipos,3) + R(k1)*(nodes(k1,3)+disp(k1*nds  ,1))   
end do
 
!write(10,FMT=10) nel, cpi, GP_coords(1,1), GP_coords(1,2), GP_coords(1,3)
           
10 format(2(I4,','),3(E,','))

close(10)


end SUBROUTINE
