subroutine ShapeFuncBar(nel,xii,etai,ubar,vbar,R)
    
    use Mod_Variables
    
    implicit none
    
    !Element number
    integer(4)::nel,idx1,idx2
    real(8),intent(IN)::xii,etai
    
    real(8),dimension(ncpx+p-1),intent(IN)::ubar
    real(8),dimension(ncpy+q-1),intent(IN)::vbar
    
    !Array of bivariate NURBS basis functions
    real(8),dimension(p*q),intent(OUT)::R
    !Array of bivariate NURBS basis functions derivatives wrt physical coordinates
    real(8),dimension(p*q,nds)::dRdx
    !Jacobian Determinant
    real(8)::detJ
    
    real(8)::detj1,detj2
    
    !NURBS coordinates
    integer(4)::ni,nj,i,j,aa,bb,cc,loc_num,errorflag
    real(8)::sumxi,sumeta,sumzeta,sumtot    !,nk
    real(8),dimension(nds)::temp1
    
    !Parametric coordinates
    real(8)::xi,eta  !,zeta
    
    !Arrays of univariate B-Spline basis functions
    real(8),dimension(p)::N
    real(8),dimension(q)::M
    
    !Basis functions derivatives wrt parametric coordinates
    real(8),dimension(p)::dNdxi
    real(8),dimension(q)::dMdeta
    
    !Array of bivariate NURBS basis functions derivatives wrt parametric coordinates
    real(8),dimension((p)*(q),nds)::dRdxi
    
    !Derivative of the physical coordinates wrt parametric coordinates and its inverse
    real(8),dimension(nds,nds)::dxdxi,dxidx
    
    !Derivatives of parametric coordinates wrt parent element coordinates
    real(8),dimension(nds,nds)::dxi_dtildexi
    
    !Jacobian matrix
    real(8),dimension(nds,nds)::jac
    
    real(8)::tmp,tol
    
    tol = 1.0d-3
    
    !Initializations
    R = 0.0d0
    dRdx = 0.0d0
    dRdxi = 0.0d0
    detJ = 0.0d0
    
    ni = 0
    nj = 0
    
    xi   = 0.0d0
    eta  = 0.0d0
    
    N  = 0.0d0
    M  = 0.0d0
    dNdxi = 0.0d0
    dMdeta = 0.0d0
    
    dxdxi = 0.0d0
    dxidx = 0.0d0
    
    dxi_dtildexi = 0.0d0
    
    i = 0
    j = 0
    aa = 0
    bb = 0
    
    loc_num = 0
    
    jac     = 0.0d0
    sumxi   = 0.0d0
    sumeta  = 0.0d0
    sumzeta = 0.0d0
    sumtot  = 0.0d0
    
    !NURBS coordinates for local node 1
    ni = INN(IEN(nel,1),1)
    nj = INN(IEN(nel,1),2)
    
    !Calculate location of the integration point in the parametric space (xi,eta) based on the coordinates of the parent element (xii,etai)
    !xi  = ((u_knot(ni+1) - u_knot(ni))*xii  + (u_knot(ni+1) + u_knot(ni)))/2.0d0
    !eta = ((v_knot(nj+1) - v_knot(nj))*etai + (v_knot(nj+1) + v_knot(nj)))/2.0d0
    
    xi  = ((ubar(ni) - ubar(ni-1))*xii  + (ubar(ni) + ubar(ni-1)))/2.0d0
    eta = ((vbar(nj) - vbar(nj-1))*etai + (vbar(nj) + vbar(nj-1)))/2.0d0
    
    !Determinant of the jacobian parametric-parent
    detj1 = (U_knot(ni+1)-U_knot(ni))*(V_knot(nj+1)-V_knot(nj))/4.0d0
    
!    open(unit=1,file='shp.txt')
!    xi = 0.0d0
!    do j=0,50
    
    !Calculate Basis functions and derivatives
    call BSplineBasisAndDeriv(ncpx-1,p-1,xi,ubar,N,dNdxi)
    call BSplineBasisAndDeriv(ncpy-1,q-1,eta,vbar,M,dMdeta)
    
!    write(1,*)N(1),N(2)
!    xi = xi + 1.0d0/50.0d0
!    end do
!    
!    close(1)
    
    !Build numerators and denominators
    do j=0,q-1
        do i=0,p-1
            loc_num = loc_num + 1     !Local basis function number
            R(loc_num) = N(p-i)*M(q-j) !*B_net(ni-i,nj-j,nds+1)
            sumtot = sumtot + R(loc_num)
            
            dRdxi(loc_num,1) = dNdxi(p-i)*M(q-j) !*B_net(ni-i,nj-j,nds+1)
            sumxi = sumxi + dRdxi(loc_num,1)
            
            dRdxi(loc_num,2) = N(p-i)*dMdeta(q-j) !*B_net(ni-i,nj-j,nds+1)
            sumeta = sumeta + dRdxi(loc_num,2)
        end do
    end do
    
    continue
    
!    do loc_num=1,p*q
!        R(loc_num) = R(loc_num)/sumtot
!        dRdxi(loc_num,1) = (dRdxi(loc_num,1)*sumtot - R(loc_num)*sumxi )/(sumtot*sumtot)
!        dRdxi(loc_num,2) = (dRdxi(loc_num,2)*sumtot - R(loc_num)*sumeta)/(sumtot*sumtot) 
!    end do
    
    !Gradient of mapping from parameter space to physical space
!    loc_num = 0
!    do j=0,q-1
!        do i=0,p-1
!            
!            loc_num = loc_num + 1
!            
!            dxdxi(1,1) = dxdxi(1,1) + Points(IENbar(nel,loc_num),1)*dRdxi(loc_num,1)
!            dxdxi(2,1) = dxdxi(2,1) + Points(IENbar(nel,loc_num),2)*dRdxi(loc_num,1)
!            dxdxi(1,2) = dxdxi(1,2) + Points(IENbar(nel,loc_num),1)*dRdxi(loc_num,2)
!            dxdxi(2,2) = dxdxi(2,2) + Points(IENbar(nel,loc_num),2)*dRdxi(loc_num,2)
!            
!        end do
!    end do
!    
!    detj2 = dxdxi(1,1)*dxdxi(2,2) - dxdxi(1,2)*dxdxi(2,1)
!    
!    !Inverse of the gradient
!    temp1 = 0.0d0
!    dxidx = dxdxi
!    call gaussj(dxidx,nds,temp1,errorflag)
!    
!    
!    dRdx = matmul(dRdxi,dxidx)
!    
!    
!    detj = detj1*detj2 
    
    continue
    
end subroutine ShapeFuncBar