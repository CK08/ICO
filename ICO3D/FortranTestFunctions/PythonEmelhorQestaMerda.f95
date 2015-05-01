!PythonEmelhorQestaMerda.f95
program teste
use loading

implicit none

character*256::FileName
integer(4)::i,j,k,l,nze,k1,k2,kt,errorflag
integer(4)::count,iter,inc
integer(4)::isolver
integer(4)::nred
!    
!!Iterative Solver Variables ----
real(8),dimension(:),allocatable::r1,r2,vs,ws,ys
!external::aprod,msolve
integer(4)::istop,itn
real(8)::Anorm,Rnorm,ynorm,Acond,STol
integer(4)::nout,itnlim
real(8)::shift
logical::checkA,goodb,precon
    
!Variables for solver3
real(8),dimension(:,:),allocatable::coef
real(8),dimension(:),allocatable::AR
integer(4),dimension(:),allocatable::IA,JA
real(8),dimension(:),allocatable::wksp,rparm
integer(4),dimension(:),allocatable::colnz,iwksp,iparm
integer(4)::mnz,countx,county
    
    
real(8)::da,tst
    
real(8), dimension(:,:), allocatable::Ke,dispin,finte
real(8)::absRes,absDisp,sumRes,sumFext,sumd,sumdt
    
real(8), dimension(:,:), allocatable::MatA,MatC,MatV
real(8), dimension(:), allocatable::temp1
    
real(8),dimension(:,:), allocatable::coordi
real(8),dimension(:,:), allocatable::MMult
    
real(4)::CPUt
real(4),dimension(2)::atime
    
    
write(*,*)'!------------------------------------------------------------------------------'
write(*,*)'!                                                                             !'
write(*,*)'!      %%%%%%%%%%%%%%%%   %%%%%%%%%%%%   %%%%%%%%%%%%%%%  %%%%%%%%%  %%%%     !'
write(*,*)'!            %%          %%             %%           %%         %%  %%  %%    !'
write(*,*)'!           %%          %%             %%           %%         %%  %%    %%   !'
write(*,*)'!          %%          %%             %%           %%   %%%%%%%%  %%      %%  !'
write(*,*)'!         %%          %%             %%           %%         %%  %%      %%   !'
write(*,*)'!        %%          %%             %%           %%         %%  %%    %%      !'
write(*,*)'!       %%          %%             %%           %%         %%  %%  %%         !'
write(*,*)'! %%%%%%%%%%%%%%   %%%%%%%%%%%%   %%%%%%%%%%%%%%%  %%%%%%%%%  %%%%            !'
write(*,*)'!                                                                             !'
write(*,*)'!------------------------------------------------------------------------------'
write(*,*)'!                                                                             !'
write(*,*)'!                   Isogeometric COde (ICO) for 3D applications               !'
write(*,*)'!                                                                             !'
write(*,*)'!------------------------------------------------------------------------------' 
    
write(*,*)''
write(*,*)'Input file name (including file extension): '


end program
