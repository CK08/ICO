!-----------------------------------------------------------------------------------------!
! Subroutine to calculate the Kdelta and KLM has defined by Matzen et al. 
! "A point to segment contact formulation for isogeometric, NURBS based finite elements."
! but with Ns, N0s and Ts formulated by Dimitri et al. "Isogeometric large deformation 
! frictionless contact using T-splines"
!-----------------------------------------------------------------------------------------!
subroutine KLM_Matzen(iLM,lm,gn,nor,tang,d2xdxi2,Rs,Rm,dRmdxi,a11,b11,Ns,N0s,Ts,T0s,Kdelta,FLM,KLM,GLM)

use Mod_Variables
implicit none

real(8),intent(IN)::Rs,Rm,dRmdxi !,dimension(ncp_m) 
real(8),intent(IN)::tang,nor,d2xdxi2 !,dimension(nds,1)
real(8),intent(IN)::lm,gn
real(8),intent(INOUT)::a11,b11
real(8),intent(INOUT)::Ns,N0s,Ts,T0s !dimension((ncp_m+ncp_s)*nds,1),
real(8),intent(INOUT)::Kdelta !dimension((ncp_m+ncp_s)*nds,(ncp_m+ncp_s)*nds),
real(8),intent(INOUT)::FLM,KLM,GLM
integer(4),intent(INOUT)::iLM
! Variaveis para calculos
integer(4)::i,j
real(8)::m,c1,c2,c3,c4
!------------------------------------------!
! Subroutine para deformavel-rigido apenas !
!------------------------------------------!

! Determinar a metrica e a curvatura da fronteira (a11, b11)
a11 = lm*lm
b11 = d2xdxi2(1,1)*nor(1,1) + d2xdxi2(2,1)*nor(2,1)                            
! Calcular os vectores Ns, Ts e N0s 
Ns = 0.0d0
Ts = 0.0d0
N0s = 0.0d0
! Parte Slave dos vectores 
do i=1,ncp_s
        Ns(i*nds-1,1) = Rs(i)*nor(1,1)
        Ns(i*nds  ,1) = Rs(i)*nor(2,1)
        Ts(i*nds-1,1) = Rs(i)*tang(1,1)
        Ts(i*nds  ,1) = Rs(i)*tang(2,1)
end do
! Parte Master dos vectores
do j=ncp_s+1,ncp_s+ncp_m 
        N0s(j*nds-1,1) = dRmdxi(j-ncp_s)*nor(1,1)
        N0s(j*nds  ,1) = dRmdxi(j-ncp_s)*nor(2,1)
        Ns(j*nds-1,1) = -Rm(j-ncp_s)*nor(1,1)
        Ns(j*nds  ,1) = -Rm(j-ncp_s)*nor(2,1)
        Ts(j*nds-1,1) = -Rm(j-ncp_s)*tang(1,1)
        Ts(j*nds  ,1) = -Rm(j-ncp_s)*tang(2,1)
end do                             

! Calcular os coeficientes m, c1, c2, c3, c4
m = a11 - gn*b11
c1 = (b11*b11*gn*gn)/(lm*m*m) + (b11*gn)/(lm*m) - (b11*lm*gn)/(m*m) - (lm)/(m) 
c2 = (b11*b11*gn*gn)/(lm*m*m) + (b11*gn)/(lm*m) - (b11*lm*gn)/(m*m) - (lm)/(m)  
c3 = (b11*b11*gn*gn*gn)/(lm*lm*m*m) + (2*b11*gn*gn)/(lm*lm*m) + (gn)/(lm*lm) - (b11*gn*gn)/(m*m) - (2*gn)/(m) 
c4 = (b11*b11*gn)/(m*m) - (b11*lm*lm)/(m*m) 
                            
! Calcular o Kdelta (multiplicando ja pelo multiplicador de lagrange do ponto de Greville)
!
! Matzen formulation 
Kdelta = c1*matmul(N0s,transpose(Ts)) + c2*matmul(Ts,transpose(N0s)) + c3*matmul(N0s,transpose(N0s)) + c4*matmul(Ts,transpose(Ts))                          
 
! Incrementar o contador iLM 
iLM = iLM + 1
                            
! Alocar a matrix KLM (deformavel-rigido/deformavel-deformavel)







end subroutine
