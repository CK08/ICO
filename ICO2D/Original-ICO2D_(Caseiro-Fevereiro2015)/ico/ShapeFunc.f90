!--------------------------------------------------------------------------------------------------------
! Subroutine to calculate the shape functions of a 2D element
! similar to the one presented in the book
! Isogeometric Analysis -Toward Integration of CAD and FEA
! but with a few changes
!
! Input: nel - Element number
!        xii, etai - Integration point coordinates in the parent element space
!        
! Output: R - NURBS basis functions
!         dRdx - NURBS basis functions derivatives wrt physical coordinates
!         detj - determinant of the jacobian matrix from the parent element 
!                to the physical domain for integration purpouses 
!--------------------------------------------------------------------------------------------------------

subroutine ShapeFunc(nel,xii,etai,R,dRdx,detj)
    
    use Mod_Variables
    
    implicit none
    
    !Element number
    integer(4)::nel,idx1,idx2
    real(8),intent(IN)::xii,etai
    !Array of bivariate NURBS basis functions
    real(8),dimension((p+1)*(q+1)),intent(OUT)::R
    !Array of bivariate NURBS basis functions derivatives wrt physical coordinates
    real(8),dimension((p+1)*(q+1),nds),intent(OUT)::dRdx
    !Jacobian Determinant
    real(8),intent(OUT)::detJ
    
    real(8)::detj1,detj2
    
    !NURBS coordinates
    integer(4)::ni,nj,i,j,aa,bb,cc,loc_num,errorflag
    real(8)::sumxi,sumeta,sumzeta,sumtot    !,nk
    real(8),dimension(nds)::temp1
    
    !Parametric coordinates
    real(8)::xi,eta  !,zeta
    
    !Arrays of univariate B-Spline basis functions
    real(8),dimension(p+1)::N
    real(8),dimension(q+1)::M
    
    !Basis functions derivatives wrt parametric coordinates
    real(8),dimension(p+1)::dNdxi
    real(8),dimension(q+1)::dMdeta
    
    !Array of bivariate NURBS basis functions derivatives wrt parametric coordinates
    real(8),dimension((p+1)*(q+1),nds)::dRdxi
    
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
    xi  = ((u_knot(ni+1) - u_knot(ni))*xii  + (u_knot(ni+1) + u_knot(ni)))/2.0d0
    eta = ((v_knot(nj+1) - v_knot(nj))*etai + (v_knot(nj+1) + v_knot(nj)))/2.0d0
    
    !Determinant of the jacobian parametric-parent
    detj1 = (U_knot(ni+1)-U_knot(ni))*(V_knot(nj+1)-V_knot(nj))/4.0d0
    
    !Calculate Basis functions and derivatives
    call BSplineBasisAndDeriv(ncpx,p,xi,u_knot,N,dNdxi)
    call BSplineBasisAndDeriv(ncpy,q,eta,v_knot,M,dMdeta)
    
    !Build numerators and denominators
    do j=0,q
        do i=0,p
            loc_num = loc_num + 1     !Local basis function number
            R(loc_num) = N(p+1-i)*M(q+1-j)*Weights(IEN(nel,loc_num))
            sumtot = sumtot + R(loc_num)
            
            dRdxi(loc_num,1) = dNdxi(p+1-i)*M(q+1-j)*Weights(IEN(nel,loc_num))
            sumxi = sumxi + dRdxi(loc_num,1)
            
            dRdxi(loc_num,2) = N(p+1-i)*dMdeta(q+1-j)*Weights(IEN(nel,loc_num))
            sumeta = sumeta + dRdxi(loc_num,2)
        end do
    end do
    
    
    do loc_num=1,(p+1)*(q+1)
        dRdxi(loc_num,1) = (dRdxi(loc_num,1)*sumtot - R(loc_num)*sumxi )/(sumtot*sumtot)
        dRdxi(loc_num,2) = (dRdxi(loc_num,2)*sumtot - R(loc_num)*sumeta)/(sumtot*sumtot) 
        R(loc_num) = R(loc_num)/sumtot
    end do
    
    
    !Gradient of mapping from parameter space to physical space
    loc_num = 0
    do j=0,q
        do i=0,p
            
            loc_num = loc_num + 1
            
            dxdxi(1,1) = dxdxi(1,1) + dRdxi(loc_num,1)*Points(IEN(nel,loc_num),1)
            dxdxi(2,1) = dxdxi(2,1) + dRdxi(loc_num,1)*Points(IEN(nel,loc_num),2)
            dxdxi(1,2) = dxdxi(1,2) + dRdxi(loc_num,2)*Points(IEN(nel,loc_num),1)
            dxdxi(2,2) = dxdxi(2,2) + dRdxi(loc_num,2)*Points(IEN(nel,loc_num),2)
            
!            do aa = 1, nds
!                do bb = 1, nds
!                    dxdxi(aa,bb) = dxdxi(aa,bb) + B_net(ni-i,nj-j,aa)*dRdxi(loc_num,bb)
!                end do
!            end do
            
            !loc_num = loc_num - 1
            
        end do
    end do
    
    detj2 = dxdxi(1,1)*dxdxi(2,2) - dxdxi(1,2)*dxdxi(2,1)
    
    !Inverse of the gradient
    temp1 = 0.0d0
    dxidx = dxdxi
    call gaussj(dxidx,nds,temp1,errorflag)
    
!    do loc_num=1,(p+1)*(q+1)
!        do aa=1,nds
!            do bb=1,nds
!                dRdx(loc_num,aa) = dRdx(loc_num,aa) + dRdxi(loc_num,bb)*dxidx(bb,aa)
!            end do
!        end do
!    end do
    
    dRdx = matmul(dRdxi,dxidx)
    
    !Gradient mapping from parent element to parameter space
!    dxi_dtildexi(1,1) = (U_knot(ni+1)-U_knot(ni))/2.0d0
!    dxi_dtildexi(2,2) = (V_knot(nj+1)-V_knot(nj))/2.0d0
    
    
!    do aa=1,nds
!        do bb=1,nds
!            do cc=1,nds
!                Jac(aa,bb) = Jac(aa,bb) + dxdxi(aa,cc)*dxi_dtildexi(cc,bb)
!            end do
!        end do
!    end do
    
    !Gradient mapping from parent element to physical space 
!    jac = 0.0d0
!    jac = matmul(dxdxi,dxi_dtildexi)
    
    !Jacobian determinant
    !detj = 0.0d0
    !detj = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)
    
    detj = detj1*detj2 
    
    continue
    
end subroutine ShapeFunc





subroutine ShapeFunc2(nel,xii,etai,R,dRdx,jac,detj,udisp)
    
    use Mod_Variables
    
    implicit none
    
    !Element number
    integer(4)::nel,idx1,idx2
    real(8),intent(IN)::xii,etai
    !Array of bivariate NURBS basis functions
    real(8),dimension((p+1)*(q+1)),intent(OUT)::R
    !Array of bivariate NURBS basis functions derivatives wrt physical coordinates
    real(8),dimension((p+1)*(q+1),nds),intent(OUT)::dRdx
    !Jacobian Determinant
    real(8),intent(OUT)::detJ
    real(8),dimension((p+1)*(q+1)*nds,1),intent(IN)::udisp
    !Jacobian matrix
    real(8),dimension(nds,nds),intent(OUT)::jac
    
    real(8)::detj1,detj2
    
    !NURBS coordinates
    integer(4)::ni,nj,i,j,aa,bb,cc,loc_num,errorflag
    real(8)::sumxi,sumeta,sumzeta,sumtot    !,nk
    real(8),dimension(nds)::temp1
    
    !Parametric coordinates
    real(8)::xi,eta  !,zeta
    
    !Arrays of univariate B-Spline basis functions
    real(8),dimension(p+1)::N
    real(8),dimension(q+1)::M
    
    !Basis functions derivatives wrt parametric coordinates
    real(8),dimension(p+1)::dNdxi
    real(8),dimension(q+1)::dMdeta
    
    !Array of bivariate NURBS basis functions derivatives wrt parametric coordinates
    real(8),dimension((p+1)*(q+1),nds)::dRdxi
    
    !Derivative of the physical coordinates wrt parametric coordinates and its inverse
    real(8),dimension(nds,nds)::dxdxi,dxidx
    
    !Derivatives of parametric coordinates wrt parent element coordinates
    real(8),dimension(nds,nds)::dxi_dtildexi
    
    
    
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
    
    if(u_knot(ni+1) == u_knot(ni)) ni = ni + 1
    if(v_knot(nj+1) == v_knot(nj)) nj = nj + 1
        
    
    !Calculate location of the integration point in the parametric space (xi,eta) based on the coordinates of the parent element (xii,etai)
    xi  = ((u_knot(ni+1) - u_knot(ni))*xii  + (u_knot(ni+1) + u_knot(ni)))/2.0d0
    eta = ((v_knot(nj+1) - v_knot(nj))*etai + (v_knot(nj+1) + v_knot(nj)))/2.0d0
    
    !Determinant of the jacobian parametric-parent
    detj1 = (U_knot(ni+1)-U_knot(ni))*(V_knot(nj+1)-V_knot(nj))/4.0d0
    
    !Calculate Basis functions and derivatives
    call BSplineBasisAndDeriv(ncpx,p,xi,u_knot,N,dNdxi)
    call BSplineBasisAndDeriv(ncpy,q,eta,v_knot,M,dMdeta)
    
    !Build numerators and denominators
    do j=0,q
        do i=0,p
            loc_num = loc_num + 1     !Local basis function number
            R(loc_num) = N(p+1-i)*M(q+1-j)*Weights(IEN(nel,loc_num))
            sumtot = sumtot + R(loc_num)
            
            dRdxi(loc_num,1) = dNdxi(p+1-i)*M(q+1-j)*Weights(IEN(nel,loc_num))
            sumxi = sumxi + dRdxi(loc_num,1)
            
            dRdxi(loc_num,2) = N(p+1-i)*dMdeta(q+1-j)*Weights(IEN(nel,loc_num))
            sumeta = sumeta + dRdxi(loc_num,2)
        end do
    end do
    
    
    do loc_num=1,(p+1)*(q+1)
        dRdxi(loc_num,1) = (dRdxi(loc_num,1)*sumtot - R(loc_num)*sumxi )/(sumtot*sumtot)
        dRdxi(loc_num,2) = (dRdxi(loc_num,2)*sumtot - R(loc_num)*sumeta)/(sumtot*sumtot) 
        R(loc_num) = R(loc_num)/sumtot
    end do
    
    
    !Gradient of mapping from parameter space to physical space
    loc_num = 0
    do j=0,q
        do i=0,p
            
            loc_num = loc_num + 1
            
            dxdxi(1,1) = dxdxi(1,1) + dRdxi(loc_num,1)*(Points(IEN(nel,loc_num),1) + udisp(loc_num*nds-1,1))
            dxdxi(2,1) = dxdxi(2,1) + dRdxi(loc_num,1)*(Points(IEN(nel,loc_num),2) + udisp(loc_num*nds  ,1))
            dxdxi(1,2) = dxdxi(1,2) + dRdxi(loc_num,2)*(Points(IEN(nel,loc_num),1) + udisp(loc_num*nds-1,1))
            dxdxi(2,2) = dxdxi(2,2) + dRdxi(loc_num,2)*(Points(IEN(nel,loc_num),2) + udisp(loc_num*nds  ,1))

        end do
    end do
    
    detj2 = dxdxi(1,1)*dxdxi(2,2) - dxdxi(1,2)*dxdxi(2,1)
    
    !Inverse of the gradient
    temp1 = 0.0d0
    dxidx = dxdxi
    call gaussj(dxidx,nds,temp1,errorflag)
    
    detj = detj1*detj2 
    detj = abs(detj)
    
    
    dRdx = matmul(dRdxi,dxidx)
    
    !Gradient mapping from parent element to parameter space
    dxi_dtildexi(1,1) = (U_knot(ni+1)-U_knot(ni))/2.0d0
    dxi_dtildexi(2,2) = (V_knot(nj+1)-V_knot(nj))/2.0d0
    
!    do aa=1,nds
!        do bb=1,nds
!            do cc=1,nds
!                Jac(aa,bb) = Jac(aa,bb) + dxdxi(aa,cc)*dxi_dtildexi(cc,bb)
!            end do
!        end do
!    end do
    
    !Gradient mapping from parent element to physical space 
    jac = 0.0d0
    jac = matmul(dxdxi,dxi_dtildexi)
    
    !Jacobian determinant
    !detj = 0.0d0
    !detj = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)
    continue
    
end subroutine ShapeFunc2
