!--------------------------------------------------------------------------------------------------------
! Subroutine to calculate the shape functions of a 3D element
! similar to the one presented in the book
! Isogeometric Analysis -Toward Integration of CAD and FEA
! but with a few changes
!
! Input: nel - Element number
!        xii, etai, zetai - Integration point coordinates in the parent element space
!        
! Output: R - NURBS basis functions
!         dRdx - NURBS basis functions derivatives wrt physical coordinates
!         detj - determinant of the jacobian matrix from the parent element 
!                to the physical domain for integration purpouses 
!--------------------------------------------------------------------------------------------------------

subroutine ShapeFunc(nel,xii,etai,zetai,R,dRdx,detj)
    
    use Mod_Variables
    
    implicit none
    
    !Element number
    integer(4)::nel,idx1,idx2,idx3
    real(8),intent(IN)::xii,etai,zetai
    !Array of bivariate NURBS basis functions
    real(8),dimension((p+1)*(q+1)*(w+1)),intent(OUT)::R
    !Array of bivariate NURBS basis functions derivatives wrt physical coordinates
    real(8),dimension((p+1)*(q+1)*(w+1),nds),intent(OUT)::dRdx
    !Jacobian Determinant
    real(8),intent(OUT)::detJ
    
    real(8)::detj1,detj2
    
    !NURBS coordinates
    integer(4)::ni,nj,nk,i,j,k,aa,bb,cc,loc_num,errorflag
    real(8)::sumxi,sumeta,sumzeta,sumtot    !,nk
    real(8),dimension(nds)::temp1
    
    !Parametric coordinates
    real(8)::xi,eta,zeta
    
    !Arrays of univariate B-Spline basis functions
    real(8),dimension(p+1)::N
    real(8),dimension(q+1)::M
    real(8),dimension(w+1)::L
    
    !Basis functions derivatives wrt parametric coordinates
    real(8),dimension(p+1)::dNdxi
    real(8),dimension(q+1)::dMdeta
    real(8),dimension(w+1)::dLdzeta
    
    !Array of bivariate NURBS basis functions derivatives wrt parametric coordinates
    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdxi
    
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
    nk = 0
    
    xi   = 0.0d0
    eta  = 0.0d0
    zeta = 0.0d0
    
    N = 0.0d0
    M = 0.0d0
    L = 0.0d0
    
    dNdxi = 0.0d0
    dMdeta = 0.0d0
    dLdzeta = 0.0d0
    
    dxdxi = 0.0d0
    dxidx = 0.0d0
    
    dxi_dtildexi = 0.0d0
    
    i = 0
    j = 0
    k = 0
    
    aa = 0
    bb = 0
    cc = 0
    
    loc_num = 0
    
    jac     = 0.0d0
    sumxi   = 0.0d0
    sumeta  = 0.0d0
    sumzeta = 0.0d0
    sumtot  = 0.0d0
    
    !NURBS coordinates for local node 1
    ni = INN(IEN(nel,1),1)
    nj = INN(IEN(nel,1),2)
    nk = INN(IEN(nel,1),3)
    
    !Calculate location of the integration point in the parametric space (xi,eta) based on the coordinates of the parent element (xii,etai)
    xi   = ((u_knot(ni+1) - u_knot(ni))*xii   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
    eta  = ((v_knot(nj+1) - v_knot(nj))*etai  + (v_knot(nj+1) + v_knot(nj)))/2.0d0
    zeta = ((w_knot(nk+1) - w_knot(nk))*zetai + (w_knot(nk+1) + w_knot(nk)))/2.0d0
    
    !Determinant of the jacobian parametric-parent
    detj1 = (U_knot(ni+1)-U_knot(ni))*(V_knot(nj+1)-V_knot(nj))*(W_knot(nk+1)-W_knot(nk))/8.0d0
    
    
!    xi = 0.0d0
!    open(unit=999,file='Basis.txt',access='append')
!    do k=1,500
!    
!    xi = xi + 1.0d0/500.0d0
    
    !Calculate Basis functions and derivatives
    call BSplineBasisAndDeriv(ncpx,p,xi,u_knot,N,dNdxi)
    call BSplineBasisAndDeriv(ncpy,q,eta,v_knot,M,dMdeta)
    call BSplineBasisAndDeriv(ncpz,w,zeta,w_knot,L,dLdzeta)
    
!    
!    write(999,FMT=11)N(1),N(2),N(3),dNdxi(1),dNdxi(2),dNdxi(3)
!    
!    11 format(6(E,1x))
!    
!    end do
!    close(999)
    
    
    
    !Build numerators and denominators
    
    do k=0,w
        do j=0,q
            do i=0,p
                loc_num = loc_num + 1     !Local basis function number
                R(loc_num) = N(p+1-i)*M(q+1-j)*L(w+1-k)*weight(IEN(nel,loc_num)) !B_net(ni-i,nj-j,nk-k,nds+1)
                sumtot = sumtot + R(loc_num)
                
                dRdxi(loc_num,1) = dNdxi(p+1-i)*M(q+1-j)*L(w+1-k)*weight(IEN(nel,loc_num)) !*B_net(ni-i,nj-j,nk-k,nds+1)
                sumxi = sumxi + dRdxi(loc_num,1)
                
                dRdxi(loc_num,2) = N(p+1-i)*dMdeta(q+1-j)*L(w+1-k)*weight(IEN(nel,loc_num)) !*B_net(ni-i,nj-j,nk-k,nds+1)
                sumeta = sumeta + dRdxi(loc_num,2)
                
                dRdxi(loc_num,3) = N(p+1-i)*M(q+1-j)*dLdzeta(w+1-k)*weight(IEN(nel,loc_num)) !*B_net(ni-i,nj-j,nk-k,nds+1)
                sumzeta = sumzeta + dRdxi(loc_num,3)
                
!                loc_num = loc_num + 1     !Local basis function number
!                R(loc_num) = N(p+1-i)*M(q+1-j)*L(w+1-k)*B_net(ni-i,nj-j,nk-k,nds+1)
!                sumtot = sumtot + R(loc_num)
!                
!                dRdxi(loc_num,1) = dNdxi(p+1-i)*M(q+1-j)*L(w+1-k)*B_net(ni-i,nj-j,nk-k,nds+1)
!                sumxi = sumxi + dRdxi(loc_num,1)
!                
!                dRdxi(loc_num,2) = N(p+1-i)*dMdeta(q+1-j)*L(w+1-k)*B_net(ni-i,nj-j,nk-k,nds+1)
!                sumeta = sumeta + dRdxi(loc_num,2)
!                
!                dRdxi(loc_num,3) = N(p+1-i)*M(q+1-j)*dLdzeta(w+1-k)*B_net(ni-i,nj-j,nk-k,nds+1)
!                sumzeta = sumzeta + dRdxi(loc_num,3)
            end do
        end do
    end do
    
    
    do i=1,(p+1)*(q+1)*(w+1)
        dRdxi(i,1) = dRdxi(i,1)/sumtot - (R(i)*sumxi  )/(sumtot**2)
        dRdxi(i,2) = dRdxi(i,2)/sumtot - (R(i)*sumeta )/(sumtot**2)
        dRdxi(i,3) = dRdxi(i,3)/sumtot - (R(i)*sumzeta)/(sumtot**2) 
        R(i) = R(i)/sumtot
    end do
    
  
    !Gradient of mapping from parameter space to physical space
    loc_num = 0
    do k=0,w
        do j=0,q
            do i=0,p
                
                loc_num = loc_num + 1
                
!                dxdxi(1,1) = dxdxi(1,1) + B_net(ni-i,nj-j,nk-k,1)*dRdxi(loc_num,1)
!                dxdxi(2,1) = dxdxi(2,1) + B_net(ni-i,nj-j,nk-k,2)*dRdxi(loc_num,1)
!                dxdxi(3,1) = dxdxi(3,1) + B_net(ni-i,nj-j,nk-k,3)*dRdxi(loc_num,1)
!                
!                dxdxi(1,2) = dxdxi(1,2) + B_net(ni-i,nj-j,nk-k,1)*dRdxi(loc_num,2)
!                dxdxi(2,2) = dxdxi(2,2) + B_net(ni-i,nj-j,nk-k,2)*dRdxi(loc_num,2)
!                dxdxi(3,2) = dxdxi(3,2) + B_net(ni-i,nj-j,nk-k,3)*dRdxi(loc_num,2)
!                
!                dxdxi(1,3) = dxdxi(1,3) + B_net(ni-i,nj-j,nk-k,1)*dRdxi(loc_num,3)
!                dxdxi(2,3) = dxdxi(2,3) + B_net(ni-i,nj-j,nk-k,2)*dRdxi(loc_num,3)
!                dxdxi(3,3) = dxdxi(3,3) + B_net(ni-i,nj-j,nk-k,3)*dRdxi(loc_num,3)
                
                dxdxi(1,1) = dxdxi(1,1) + dRdxi(loc_num,1)*Points(IEN(nel,loc_num),1)
                dxdxi(2,1) = dxdxi(2,1) + dRdxi(loc_num,1)*Points(IEN(nel,loc_num),2)
                dxdxi(3,1) = dxdxi(3,1) + dRdxi(loc_num,1)*Points(IEN(nel,loc_num),3)
                
                dxdxi(1,2) = dxdxi(1,2) + dRdxi(loc_num,2)*Points(IEN(nel,loc_num),1)
                dxdxi(2,2) = dxdxi(2,2) + dRdxi(loc_num,2)*Points(IEN(nel,loc_num),2)
                dxdxi(3,2) = dxdxi(3,2) + dRdxi(loc_num,2)*Points(IEN(nel,loc_num),3)
                
                dxdxi(1,3) = dxdxi(1,3) + dRdxi(loc_num,3)*Points(IEN(nel,loc_num),1)
                dxdxi(2,3) = dxdxi(2,3) + dRdxi(loc_num,3)*Points(IEN(nel,loc_num),2)
                dxdxi(3,3) = dxdxi(3,3) + dRdxi(loc_num,3)*Points(IEN(nel,loc_num),3)
                
            end do
        end do
    end do
    

    !detj2 = dxdxi(1,1)*dxdxi(2,2) - dxdxi(1,2)*dxdxi(2,1)
    call fdeterm(dxdxi,detj2)
    
    !Inverse of the gradient
    temp1 = 0.0d0
    dxidx = dxdxi
    call gaussj(dxidx,nds,temp1,errorflag)
    
    dRdx = matmul(dRdxi,dxidx)
    
    detj = detj1*detj2 
    
    detj = abs(detj)
    
    continue
    
end subroutine ShapeFunc

!--------------------------------------------------------------------------------------------------------
! Subroutine to calculate the shape functions of a 3D element
! similar to the one presented in the book
! Isogeometric Analysis -Toward Integration of CAD and FEA
! but with a few more changes
!
! Input: nel - Element number
!        xii, etai, zetai - Integration point coordinates in the parent element space
!        
! Output: R - NURBS basis functions
!         dRdx - NURBS basis functions derivatives wrt physical coordinates
!         detj - determinant of the jacobian matrix from the parent element 
!                to the physical domain for integration purpouses
!         jac - jacobian matrix from the parent element to the physical domain 
!         udisp - array with the same size as the coordenates but with zeros
!--------------------------------------------------------------------------------------------------------
subroutine ShapeFunc2(nel,xii,etai,zetai,R,dRdx,detj,jac,udisp)
    
    use Mod_Variables
    
    implicit none
    
    !Element number
    integer(4)::nel,idx1,idx2,idx3
    real(8),intent(IN)::xii,etai,zetai
    !Array of bivariate NURBS basis functions
    real(8),dimension((p+1)*(q+1)*(w+1)),intent(OUT)::R
    !Array of bivariate NURBS basis functions derivatives wrt physical coordinates
    real(8),dimension((p+1)*(q+1)*(w+1),nds),intent(OUT)::dRdx
    !Jacobian Determinant
    real(8),intent(OUT)::detJ
    real(8),dimension(nds,nds),intent(OUT)::jac
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::udisp
    
    real(8)::detj1,detj2
    
    !NURBS coordinates
    integer(4)::ni,nj,nk,i,j,k,aa,bb,cc,loc_num,errorflag
    real(8)::sumxi,sumeta,sumzeta,sumtot    !,nk
    real(8),dimension(nds)::temp1
    
    !Parametric coordinates
    real(8)::xi,eta,zeta
    
    !Arrays of univariate B-Spline basis functions
    real(8),dimension(p+1)::N
    real(8),dimension(q+1)::M
    real(8),dimension(w+1)::L
    
    !Basis functions derivatives wrt parametric coordinates
    real(8),dimension(p+1)::dNdxi
    real(8),dimension(q+1)::dMdeta
    real(8),dimension(w+1)::dLdzeta
    
    !Array of bivariate NURBS basis functions derivatives wrt parametric coordinates
    real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdxi
    
    !Derivative of the physical coordinates wrt parametric coordinates and its inverse
    real(8),dimension(nds,nds)::dxidx
    
    !Derivatives of parametric coordinates wrt parent element coordinates
    real(8),dimension(nds,nds)::dxi_dtildexi
    
    !Jacobian matrix
    real(8),dimension(nds,nds)::dxdxi
    
    integer(4)::count
    
    real(8)::tmp,tol
    
    
    
    tol = 1.0d-3
    
    !Initializations
    R = 0.0d0
    dRdx = 0.0d0
    dRdxi = 0.0d0
    detJ = 0.0d0
    
    ni = 0
    nj = 0
    nk = 0
    
    xi   = 0.0d0
    eta  = 0.0d0
    zeta = 0.0d0
    
    N = 0.0d0
    M = 0.0d0
    L = 0.0d0
    
    dNdxi = 0.0d0
    dMdeta = 0.0d0
    dLdzeta = 0.0d0
    
    dxdxi = 0.0d0
    dxidx = 0.0d0
    
    dxi_dtildexi = 0.0d0
    
    i = 0
    j = 0
    k = 0
    
    aa = 0
    bb = 0
    cc = 0
    
    loc_num = 0
    
    jac     = 0.0d0
    sumxi   = 0.0d0
    sumeta  = 0.0d0
    sumzeta = 0.0d0
    sumtot  = 0.0d0
    
    !NURBS coordinates for local node 1
    ni = INN(IEN(nel,1),1)
    nj = INN(IEN(nel,1),2)
    nk = INN(IEN(nel,1),3)
    
    !Calculate location of the integration point in the parametric space (xi,eta) based on the coordinates of the parent element (xii,etai)
    xi   = ((u_knot(ni+1) - u_knot(ni))*xii   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
    eta  = ((v_knot(nj+1) - v_knot(nj))*etai  + (v_knot(nj+1) + v_knot(nj)))/2.0d0
    zeta = ((w_knot(nk+1) - w_knot(nk))*zetai + (w_knot(nk+1) + w_knot(nk)))/2.0d0
    
    !Determinant of the jacobian parametric-parent
    detj1 = (U_knot(ni+1)-U_knot(ni))*(V_knot(nj+1)-V_knot(nj))*(W_knot(nk+1)-W_knot(nk))/8.0d0
    
    !Calculate Basis functions and derivatives
    call BSplineBasisAndDeriv(ncpx,p,xi,u_knot,N,dNdxi)
    call BSplineBasisAndDeriv(ncpy,q,eta,v_knot,M,dMdeta)
    call BSplineBasisAndDeriv(ncpz,w,zeta,w_knot,L,dLdzeta)
    
    
    !Build numerators and denominators
    do k=0,w
        do j=0,q
            do i=0,p
                loc_num = loc_num + 1     !Local basis function number
                R(loc_num) = N(p+1-i)*M(q+1-j)*L(w+1-k)*weight(IEN(nel,loc_num)) !B_net(ni-i,nj-j,nk-k,nds+1)
                sumtot = sumtot + R(loc_num)
                
                dRdxi(loc_num,1) = dNdxi(p+1-i)*M(q+1-j)*L(w+1-k)*weight(IEN(nel,loc_num)) !*B_net(ni-i,nj-j,nk-k,nds+1)
                sumxi = sumxi + dRdxi(loc_num,1)
                
                dRdxi(loc_num,2) = N(p+1-i)*dMdeta(q+1-j)*L(w+1-k)*weight(IEN(nel,loc_num)) !*B_net(ni-i,nj-j,nk-k,nds+1)
                sumeta = sumeta + dRdxi(loc_num,2)
                
                dRdxi(loc_num,3) = N(p+1-i)*M(q+1-j)*dLdzeta(w+1-k)*weight(IEN(nel,loc_num)) !*B_net(ni-i,nj-j,nk-k,nds+1)
                sumzeta = sumzeta + dRdxi(loc_num,3)
            end do
        end do
    end do
    
    
    do i=1,(p+1)*(q+1)*(w+1)
        dRdxi(i,1) = dRdxi(i,1)/sumtot - (R(i)*sumxi  )/(sumtot**2)
        dRdxi(i,2) = dRdxi(i,2)/sumtot - (R(i)*sumeta )/(sumtot**2)
        dRdxi(i,3) = dRdxi(i,3)/sumtot - (R(i)*sumzeta)/(sumtot**2) 
        R(i) = R(i)/sumtot
    end do
    
  
    !Gradient of mapping from parameter space to physical space
    loc_num = 0
    count = (p+1)*(q+1)*(w+1)+1
    do k=0,w
        do j=0,q
            do i=0,p
                
                loc_num = loc_num + 1
                count = count - 1
                
                dxdxi(1,1) = dxdxi(1,1) + dRdxi(loc_num,1)*(Points(IEN(nel,loc_num),1) + udisp(loc_num*3-2,1))
                dxdxi(2,1) = dxdxi(2,1) + dRdxi(loc_num,1)*(Points(IEN(nel,loc_num),2) + udisp(loc_num*3-1,1))
                dxdxi(3,1) = dxdxi(3,1) + dRdxi(loc_num,1)*(Points(IEN(nel,loc_num),3) + udisp(loc_num*3  ,1))
                
                dxdxi(1,2) = dxdxi(1,2) + dRdxi(loc_num,2)*(Points(IEN(nel,loc_num),1) + udisp(loc_num*3-2,1))
                dxdxi(2,2) = dxdxi(2,2) + dRdxi(loc_num,2)*(Points(IEN(nel,loc_num),2) + udisp(loc_num*3-1,1))
                dxdxi(3,2) = dxdxi(3,2) + dRdxi(loc_num,2)*(Points(IEN(nel,loc_num),3) + udisp(loc_num*3  ,1))
                
                dxdxi(1,3) = dxdxi(1,3) + dRdxi(loc_num,3)*(Points(IEN(nel,loc_num),1) + udisp(loc_num*3-2,1))
                dxdxi(2,3) = dxdxi(2,3) + dRdxi(loc_num,3)*(Points(IEN(nel,loc_num),2) + udisp(loc_num*3-1,1))
                dxdxi(3,3) = dxdxi(3,3) + dRdxi(loc_num,3)*(Points(IEN(nel,loc_num),3) + udisp(loc_num*3  ,1))
                
            end do
        end do
    end do
    

    call fdeterm(dxdxi,detj2)
    
    !Inverse of the gradient
    temp1 = 0.0d0
    dxidx = dxdxi
    call gaussj(dxidx,nds,temp1,errorflag)
    
    dRdx = matmul(dRdxi,dxidx)
    
    detj = detj1*detj2 
    
    detj = abs(detj)
    
    !Gradient mapping from parent element to parameter space
    dxi_dtildexi(1,1) = (U_knot(ni+1)-U_knot(ni))/2.0d0
    dxi_dtildexi(2,2) = (V_knot(nj+1)-V_knot(nj))/2.0d0
    dxi_dtildexi(3,3) = (W_knot(nk+1)-W_knot(nk))/2.0d0
    
    
    Jac = 0.0d0
    do aa=1,nds
        do bb=1,nds
            do cc=1,nds
                Jac(aa,bb) = Jac(aa,bb) + dxdxi(aa,cc)*dxi_dtildexi(cc,bb)
            end do
        end do
    end do
    
    !Gradient mapping from parent element to physical space 
!    jac = 0.0d0
    jac = matmul(dxdxi,dxi_dtildexi)
    
    !jac = dxdxi
    
    continue
    
end subroutine ShapeFunc2










subroutine ShapeFuncEAS(nel,nnodes,nds,nelems,nshpl,p,q,w,xii,etai,zetai,R,dRdx,detj,jac,udisp,Points,IEN)
    
      !use ModVariables
    
      implicit none
    
      !Element number
      integer(4)::nel,idx1,idx2,idx3,p,q,w,nnodes,nds,nelems,nshpl
      real(8),intent(IN)::xii,etai,zetai
      !Array of bivariate NURBS basis functions
      real(8),dimension((p+1)*(q+1)*(w+1)),intent(OUT)::R
      !Array of bivariate NURBS basis functions derivatives wrt physical coordinates
      real(8),dimension((p+1)*(q+1)*(w+1),nds),intent(OUT)::dRdx
      !Jacobian Determinant
      real(8),intent(OUT)::detJ
      real(8),dimension(nds,nds),intent(OUT)::jac
      real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::udisp
      
      integer(4),dimension(nelems,nshpl)::IEN
      real(8),dimension(nnodes,nds)::Points
      
      
      real(8)::detj1,detj2
    
      !NURBS coordinates
      integer(4)::ni,nj,nk,i,j,k,aa,bb,cc,loc_num,errorflag
      real(8)::sumxi,sumeta,sumzeta,sumtot    !,nk
      real(8),dimension(nds)::temp1
    
      !Parametric coordinates
      real(8)::xi,eta,zeta
    
      !Arrays of univariate B-Spline basis functions
      real(8),dimension(p+1)::N
      real(8),dimension(q+1)::M
      real(8),dimension(w+1)::L
    
      !Basis functions derivatives wrt parametric coordinates
      real(8),dimension(p+1)::dNdxi
      real(8),dimension(q+1)::dMdeta
      real(8),dimension(w+1)::dLdzeta
    
      !Array of bivariate NURBS basis functions derivatives wrt parametric coordinates
      real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdxi
    
      !Derivative of the physical coordinates wrt parametric coordinates and its inverse
      real(8),dimension(nds,nds)::dxidx
    
      !Derivatives of parametric coordinates wrt parent element coordinates
      real(8),dimension(nds,nds)::dxi_dtildexi
    
      !Jacobian matrix
      real(8),dimension(nds,nds)::dxdxi
    
      integer(4)::count
    
      real(8)::tmp,tol
      
      real(8),dimension((p+1)*2)::u_knot
      real(8),dimension((q+1)*2)::v_knot
      real(8),dimension((w+1)*2)::w_knot
      
      tol = 1.0d-3
    
      !Initializations
      R = 0.0d0
      dRdx = 0.0d0
      dRdxi = 0.0d0
      detJ = 0.0d0
    
      ni = 0
      nj = 0
      nk = 0
    
      xi   = 0.0d0
      eta  = 0.0d0
      zeta = 0.0d0
    
      N = 0.0d0
      M = 0.0d0
      L = 0.0d0
    
      dNdxi = 0.0d0
      dMdeta = 0.0d0
      dLdzeta = 0.0d0
    
      dxdxi = 0.0d0
      dxidx = 0.0d0
    
      dxi_dtildexi = 0.0d0
    
      i = 0
      j = 0
      k = 0
    
      aa = 0
      bb = 0
      cc = 0
    
      loc_num = 0
    
      jac     = 0.0d0
      sumxi   = 0.0d0
      sumeta  = 0.0d0
      sumzeta = 0.0d0
      sumtot  = 0.0d0
    
      !NURBS coordinates for local node 1
      !ni = INN(IEN(nel,1),1)
      !nj = INN(IEN(nel,1),2)
      !nk = INN(IEN(nel,1),3)
    
    !Calculate location of the integration point in the parametric space (xi,eta) based on the coordinates of the parent element (xii,etai)
      xi   = (xii   + 1.0d0)/2.0d0
      eta  = (etai  + 1.0d0)/2.0d0
      zeta = (zetai + 1.0d0)/2.0d0
    
      !Determinant of the jacobian parametric-parent
      detj1 = 1.0d0/8.0d0
    
      u_knot = 1.0d0
      v_knot = 1.0d0
      w_knot = 1.0d0
      
      do i=1,p+1
        u_knot(i)=0.0d0
      end do 
      
      do i=1,q+1
        v_knot(i)=0.0d0
      end do 
      
      do i=1,w+1
        w_knot(i)=0.0d0
      end do
    
    
      !Calculate Basis functions and derivatives
!      open(unit=1,file='C:\JCaseiro\Doutoramento\Isogeometric\
!     #Validation Examples\Hexahedral\Cook Membrane\EAS\Elastic\oi.txt')
!      xi = 0.0d0
!      do k=1,49
!      xi = xi + 1.0d0/50.0d0
!      
!      write(1,*)N(1),N(2),N(3)
!      end do
!      close(1)
!      
!      pause
      
      call BSplineBasisAndDeriv(p+1,p,xi,u_knot,N,dNdxi)
      call BSplineBasisAndDeriv(q+1,q,eta,v_knot,M,dMdeta)
      call BSplineBasisAndDeriv(w+1,w,zeta,w_knot,L,dLdzeta)
    

      !Build numerators and denominators
      do k=0,w
        do j=0,q
            do i=0,p
                loc_num = loc_num + 1     !Local basis function number
                R(loc_num) = N(p+1-i)*M(q+1-j)*L(w+1-k)       
                sumtot = sumtot + R(loc_num)
                
                dRdxi(loc_num,1) = dNdxi(p+1-i)*M(q+1-j)*L(w+1-k)
                sumxi = sumxi + dRdxi(loc_num,1)
                
                dRdxi(loc_num,2) = N(p+1-i)*dMdeta(q+1-j)*L(w+1-k)
                sumeta = sumeta + dRdxi(loc_num,2)
                
                dRdxi(loc_num,3) = N(p+1-i)*M(q+1-j)*dLdzeta(w+1-k)
                sumzeta = sumzeta + dRdxi(loc_num,3)
            end do
        end do
      end do
    
    
      do i=1,(p+1)*(q+1)*(w+1)
        dRdxi(i,1) = dRdxi(i,1)/sumtot - (R(i)*sumxi  )/(sumtot**2)
        dRdxi(i,2) = dRdxi(i,2)/sumtot - (R(i)*sumeta )/(sumtot**2)
        dRdxi(i,3) = dRdxi(i,3)/sumtot - (R(i)*sumzeta)/(sumtot**2) 
        R(i) = R(i)/sumtot
      end do
    
  
      !Gradient of mapping from parameter space to physical space
      loc_num = 0
      count = (p+1)*(q+1)*(w+1)+1
      do k=0,w
        do j=0,q
            do i=0,p
                
                loc_num = loc_num + 1
                count = count - 1
                
                dxdxi(1,1) = dxdxi(1,1) + dRdxi(loc_num,1)*(Points(IEN(nel,loc_num),1) + udisp(loc_num*3-2,1))
                dxdxi(2,1) = dxdxi(2,1) + dRdxi(loc_num,1)*(Points(IEN(nel,loc_num),2) + udisp(loc_num*3-1,1))
                dxdxi(3,1) = dxdxi(3,1) + dRdxi(loc_num,1)*(Points(IEN(nel,loc_num),3) + udisp(loc_num*3  ,1))
                
                dxdxi(1,2) = dxdxi(1,2) + dRdxi(loc_num,2)*(Points(IEN(nel,loc_num),1) + udisp(loc_num*3-2,1))
                dxdxi(2,2) = dxdxi(2,2) + dRdxi(loc_num,2)*(Points(IEN(nel,loc_num),2) + udisp(loc_num*3-1,1))
                dxdxi(3,2) = dxdxi(3,2) + dRdxi(loc_num,2)*(Points(IEN(nel,loc_num),3) + udisp(loc_num*3  ,1))
                
                dxdxi(1,3) = dxdxi(1,3) + dRdxi(loc_num,3)*(Points(IEN(nel,loc_num),1) + udisp(loc_num*3-2,1))
                dxdxi(2,3) = dxdxi(2,3) + dRdxi(loc_num,3)*(Points(IEN(nel,loc_num),2) + udisp(loc_num*3-1,1))
                dxdxi(3,3) = dxdxi(3,3) + dRdxi(loc_num,3)*(Points(IEN(nel,loc_num),3) + udisp(loc_num*3  ,1))
                
            end do
        end do
      end do
    

      call fdeterm(dxdxi,detj2)
    
      !Inverse of the gradient
      temp1 = 0.0d0
      dxidx = dxdxi
      call gaussj(dxidx,nds,temp1,errorflag)
    
      dRdx = matmul(dRdxi,dxidx)
    
      detj = detj1*detj2 
    
      detj = abs(detj)
    
      !Gradient mapping from parent element to parameter space
      dxi_dtildexi(1,1) = 1.0d0/2.0d0
      dxi_dtildexi(2,2) = 1.0d0/2.0d0
      dxi_dtildexi(3,3) = 1.0d0/2.0d0
    
    
      Jac = 0.0d0
      do aa=1,nds
        do bb=1,nds
            do cc=1,nds
                Jac(aa,bb) = Jac(aa,bb)+dxdxi(aa,cc)*dxi_dtildexi(cc,bb)
            end do
        end do
      end do
    

      continue
    
      end subroutine ShapeFuncEAS
      
      
    
    
      
      
      
subroutine ShapeFunc3(nel,xii,etai,zetai,R,dRdx,dRdxi,dRdxii,detj,jac,dxdxi,udisp)
    
    use Mod_Variables
    
    implicit none
    
    !Element number
    integer(4)::nel,idx1,idx2,idx3
    real(8),intent(IN)::xii,etai,zetai
    !Array of bivariate NURBS basis functions
    real(8),dimension((p+1)*(q+1)*(w+1)),intent(OUT)::R
    !Array of bivariate NURBS basis functions derivatives wrt physical coordinates
    real(8),dimension((p+1)*(q+1)*(w+1),nds),intent(OUT)::dRdx,dRdxi,dRdxii
    !Jacobian Determinant
    real(8),intent(OUT)::detJ
    real(8),dimension(nds,nds),intent(OUT)::jac,dxdxi
    real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::udisp
    
    real(8)::detj1,detj2
    
    !NURBS coordinates
    integer(4)::ni,nj,nk,i,j,k,aa,bb,cc,loc_num,errorflag
    real(8)::sumxi,sumeta,sumzeta,sumtot    !,nk
    real(8),dimension(nds)::temp1
    
    !Parametric coordinates
    real(8)::xi,eta,zeta
    
    !Arrays of univariate B-Spline basis functions
    real(8),dimension(p+1)::N
    real(8),dimension(q+1)::M
    real(8),dimension(w+1)::L
    
    !Basis functions derivatives wrt parametric coordinates
    real(8),dimension(p+1)::dNdxi
    real(8),dimension(q+1)::dMdeta
    real(8),dimension(w+1)::dLdzeta
    
    !Derivative of the physical coordinates wrt parametric coordinates and its inverse
    real(8),dimension(nds,nds)::dxidx
    
    !Derivatives of parametric coordinates wrt parent element coordinates
    real(8),dimension(nds,nds)::dxi_dtildexi

    integer(4)::count
    real(8)::tmp,tol
    
    tol = 1.0d-3
    
    !Initializations
    R = 0.0d0
    dRdx = 0.0d0
    dRdxi = 0.0d0
    detJ = 0.0d0
    
    ni = 0
    nj = 0
    nk = 0
    
    xi   = 0.0d0
    eta  = 0.0d0
    zeta = 0.0d0
    
    N = 0.0d0
    M = 0.0d0
    L = 0.0d0
    
    dNdxi = 0.0d0
    dMdeta = 0.0d0
    dLdzeta = 0.0d0
    
    dxdxi = 0.0d0
    dxidx = 0.0d0
    
    dxi_dtildexi = 0.0d0
    
    i = 0
    j = 0
    k = 0
    
    aa = 0
    bb = 0
    cc = 0
    
    loc_num = 0
    
    jac     = 0.0d0
    sumxi   = 0.0d0
    sumeta  = 0.0d0
    sumzeta = 0.0d0
    sumtot  = 0.0d0
    
    !NURBS coordinates for local node 1
    ni = INN(IEN(nel,1),1)
    nj = INN(IEN(nel,1),2)
    nk = INN(IEN(nel,1),3)
    
    !Calculate location of the integration point in the parametric space (xi,eta) based on the coordinates of the parent element (xii,etai)
    xi   = ((u_knot(ni+1) - u_knot(ni))*xii   + (u_knot(ni+1) + u_knot(ni)))/2.0d0
    eta  = ((v_knot(nj+1) - v_knot(nj))*etai  + (v_knot(nj+1) + v_knot(nj)))/2.0d0
    zeta = ((w_knot(nk+1) - w_knot(nk))*zetai + (w_knot(nk+1) + w_knot(nk)))/2.0d0
    
    !Determinant of the jacobian parametric-parent
    detj1 = (U_knot(ni+1)-U_knot(ni))*(V_knot(nj+1)-V_knot(nj))*(W_knot(nk+1)-W_knot(nk))/8.0d0
    
    !Calculate Basis functions and derivatives
    call BSplineBasisAndDeriv(ncpx,p,xi,u_knot,N,dNdxi)
    call BSplineBasisAndDeriv(ncpy,q,eta,v_knot,M,dMdeta)
    call BSplineBasisAndDeriv(ncpz,w,zeta,w_knot,L,dLdzeta)
    
    
    !Build numerators and denominators
    
    do k=0,w
        do j=0,q
            do i=0,p
                loc_num = loc_num + 1     !Local basis function number
                R(loc_num) = N(p+1-i)*M(q+1-j)*L(w+1-k)*weight(IEN(nel,loc_num)) !B_net(ni-i,nj-j,nk-k,nds+1)
                sumtot = sumtot + R(loc_num)
                
                dRdxi(loc_num,1) = dNdxi(p+1-i)*M(q+1-j)*L(w+1-k)*weight(IEN(nel,loc_num)) !*B_net(ni-i,nj-j,nk-k,nds+1)
                sumxi = sumxi + dRdxi(loc_num,1)
                
                dRdxi(loc_num,2) = N(p+1-i)*dMdeta(q+1-j)*L(w+1-k)*weight(IEN(nel,loc_num)) !*B_net(ni-i,nj-j,nk-k,nds+1)
                sumeta = sumeta + dRdxi(loc_num,2)
                
                dRdxi(loc_num,3) = N(p+1-i)*M(q+1-j)*dLdzeta(w+1-k)*weight(IEN(nel,loc_num)) !*B_net(ni-i,nj-j,nk-k,nds+1)
                sumzeta = sumzeta + dRdxi(loc_num,3)
            end do
        end do
    end do
    
    
    do i=1,(p+1)*(q+1)*(w+1)
        dRdxi(i,1) = dRdxi(i,1)/sumtot - (R(i)*sumxi  )/(sumtot**2)
        dRdxi(i,2) = dRdxi(i,2)/sumtot - (R(i)*sumeta )/(sumtot**2)
        dRdxi(i,3) = dRdxi(i,3)/sumtot - (R(i)*sumzeta)/(sumtot**2) 
        R(i) = R(i)/sumtot
    end do
    
  
    !Gradient of mapping from parameter space to physical space
    loc_num = 0
    count = (p+1)*(q+1)*(w+1)+1
    do k=0,w
        do j=0,q
            do i=0,p
                
                loc_num = loc_num + 1
                count = count - 1
                
                dxdxi(1,1) = dxdxi(1,1) + dRdxi(loc_num,1)*(Points(IEN(nel,loc_num),1) + udisp(loc_num*3-2,1))
                dxdxi(2,1) = dxdxi(2,1) + dRdxi(loc_num,1)*(Points(IEN(nel,loc_num),2) + udisp(loc_num*3-1,1))
                dxdxi(3,1) = dxdxi(3,1) + dRdxi(loc_num,1)*(Points(IEN(nel,loc_num),3) + udisp(loc_num*3  ,1))
                
                dxdxi(1,2) = dxdxi(1,2) + dRdxi(loc_num,2)*(Points(IEN(nel,loc_num),1) + udisp(loc_num*3-2,1))
                dxdxi(2,2) = dxdxi(2,2) + dRdxi(loc_num,2)*(Points(IEN(nel,loc_num),2) + udisp(loc_num*3-1,1))
                dxdxi(3,2) = dxdxi(3,2) + dRdxi(loc_num,2)*(Points(IEN(nel,loc_num),3) + udisp(loc_num*3  ,1))
                
                dxdxi(1,3) = dxdxi(1,3) + dRdxi(loc_num,3)*(Points(IEN(nel,loc_num),1) + udisp(loc_num*3-2,1))
                dxdxi(2,3) = dxdxi(2,3) + dRdxi(loc_num,3)*(Points(IEN(nel,loc_num),2) + udisp(loc_num*3-1,1))
                dxdxi(3,3) = dxdxi(3,3) + dRdxi(loc_num,3)*(Points(IEN(nel,loc_num),3) + udisp(loc_num*3  ,1))
                
            end do
        end do
    end do
    

    call fdeterm(dxdxi,detj2)
    
    !Inverse of the gradient
    temp1 = 0.0d0
    dxidx = dxdxi
    call gaussj(dxidx,nds,temp1,errorflag)
    
    dRdx = matmul(dRdxi,dxidx)
    
    detj = detj1*detj2 
    
    detj = abs(detj)
    
    !Gradient mapping from parent element to parameter space
    dxi_dtildexi(1,1) = (U_knot(ni+1)-U_knot(ni))/2.0d0
    dxi_dtildexi(2,2) = (V_knot(nj+1)-V_knot(nj))/2.0d0
    dxi_dtildexi(3,3) = (W_knot(nk+1)-W_knot(nk))/2.0d0
    
    dRdxii = matmul(dRdxi,dxi_dtildexi)
    
    Jac = 0.0d0
    do aa=1,nds
        do bb=1,nds
            do cc=1,nds
                Jac(aa,bb) = Jac(aa,bb) + dxdxi(aa,cc)*dxi_dtildexi(cc,bb)
            end do
        end do
    end do
    
    !Gradient mapping from parent element to physical space 
!    jac = 0.0d0
    jac = matmul(dxdxi,dxi_dtildexi)
    
    !jac = dxdxi
    
    continue
    
end subroutine ShapeFunc3

!subroutine ShapeFuncEAS(nel,nnodes,nds,nelems,nshpl,p,q,w,xii,etai,zetai,R,dRdx,detj,jac,udisp,Points,IEN)
!    
!      !use ModVariables
!    
!      implicit none
!    
!      !Element number
!      integer(4)::nel,idx1,idx2,idx3,p,q,w,nnodes,nds,nelems,nshpl
!      real(8),intent(IN)::xii,etai,zetai
!      !Array of bivariate NURBS basis functions
!      real(8),dimension((p+1)*(q+1)*(w+1)),intent(OUT)::R
!      !Array of bivariate NURBS basis functions derivatives wrt physical coordinates
!      real(8),dimension((p+1)*(q+1)*(w+1),nds),intent(OUT)::dRdx
!      !Jacobian Determinant
!      real(8),intent(OUT)::detJ
!      real(8),dimension(nds,nds),intent(OUT)::jac
!      real(8),dimension((p+1)*(q+1)*(w+1)*nds,1),intent(IN)::udisp
!      
!      integer(4),dimension(nelems,nshpl)::IEN
!      real(8),dimension(nnodes,nds)::Points
!      
!      
!      real(8)::detj1,detj2
!    
!      !NURBS coordinates
!      integer(4)::ni,nj,nk,i,j,k,aa,bb,cc,loc_num,errorflag
!      real(8)::sumxi,sumeta,sumzeta,sumtot    !,nk
!      real(8),dimension(nds)::temp1
!    
!      !Parametric coordinates
!      real(8)::xi,eta,zeta
!    
!      !Arrays of univariate B-Spline basis functions
!      real(8),dimension(p+1)::N
!      real(8),dimension(q+1)::M
!      real(8),dimension(w+1)::L
!    
!      !Basis functions derivatives wrt parametric coordinates
!      real(8),dimension(p+1)::dNdxi
!      real(8),dimension(q+1)::dMdeta
!      real(8),dimension(w+1)::dLdzeta
!    
!      !Array of bivariate NURBS basis functions derivatives wrt parametric coordinates
!      real(8),dimension((p+1)*(q+1)*(w+1),nds)::dRdxi
!    
!      !Derivative of the physical coordinates wrt parametric coordinates and its inverse
!      real(8),dimension(nds,nds)::dxidx
!    
!      !Derivatives of parametric coordinates wrt parent element coordinates
!      real(8),dimension(nds,nds)::dxi_dtildexi
!    
!      !Jacobian matrix
!      real(8),dimension(nds,nds)::dxdxi
!    
!      integer(4)::count
!    
!      real(8)::tmp,tol
!      
!      real(8),dimension((p+1)*2)::u_knot
!      real(8),dimension((q+1)*2)::v_knot
!      real(8),dimension((w+1)*2)::w_knot
!      
!      tol = 1.0d-3
!    
!      !Initializations
!      R = 0.0d0
!      dRdx = 0.0d0
!      dRdxi = 0.0d0
!      detJ = 0.0d0
!    
!      ni = 0
!      nj = 0
!      nk = 0
!    
!      xi   = 0.0d0
!      eta  = 0.0d0
!      zeta = 0.0d0
!    
!      N = 0.0d0
!      M = 0.0d0
!      L = 0.0d0
!    
!      dNdxi = 0.0d0
!      dMdeta = 0.0d0
!      dLdzeta = 0.0d0
!    
!      dxdxi = 0.0d0
!      dxidx = 0.0d0
!    
!      dxi_dtildexi = 0.0d0
!    
!      i = 0
!      j = 0
!      k = 0
!    
!      aa = 0
!      bb = 0
!      cc = 0
!    
!      loc_num = 0
!    
!      jac     = 0.0d0
!      sumxi   = 0.0d0
!      sumeta  = 0.0d0
!      sumzeta = 0.0d0
!      sumtot  = 0.0d0
!    
!      !NURBS coordinates for local node 1
!      !ni = INN(IEN(nel,1),1)
!      !nj = INN(IEN(nel,1),2)
!      !nk = INN(IEN(nel,1),3)
!    
!      !Calculate location of the integration point in the parametric space (xi,eta) based on the coordinates of the parent element (xii,etai)
!      xi   = (xii   + 1.0d0)/2.0d0
!      eta  = (etai  + 1.0d0)/2.0d0
!      zeta = (zetai + 1.0d0)/2.0d0
!    
!      !Determinant of the jacobian parametric-parent
!      detj1 = 1.0d0/8.0d0
!    
!      u_knot = 1.0d0
!      v_knot = 1.0d0
!      w_knot = 1.0d0
!      
!      do i=1,p+1
!        u_knot(i)=0.0d0
!      end do 
!      
!      do i=1,q+1
!        v_knot(i)=0.0d0
!      end do 
!      
!      do i=1,w+1
!        w_knot(i)=0.0d0
!      end do
!    
!      call BSplineBasisAndDeriv(p+1,p,xi,u_knot,N,dNdxi)
!      call BSplineBasisAndDeriv(q+1,q,eta,v_knot,M,dMdeta)
!      call BSplineBasisAndDeriv(w+1,w,zeta,w_knot,L,dLdzeta)
!    
!      !Build numerators and denominators
!      do k=0,w
!        do j=0,q
!            do i=0,p
!                loc_num = loc_num + 1     !Local basis function number
!                R(loc_num) = N(p+1-i)*M(q+1-j)*L(w+1-k)       
!                sumtot = sumtot + R(loc_num)
!                
!                dRdxi(loc_num,1) = dNdxi(p+1-i)*M(q+1-j)*L(w+1-k)
!                sumxi = sumxi + dRdxi(loc_num,1)
!                
!                dRdxi(loc_num,2) = N(p+1-i)*dMdeta(q+1-j)*L(w+1-k)
!                sumeta = sumeta + dRdxi(loc_num,2)
!                
!                dRdxi(loc_num,3) = N(p+1-i)*M(q+1-j)*dLdzeta(w+1-k)
!                sumzeta = sumzeta + dRdxi(loc_num,3)
!            end do
!        end do
!      end do
!    
!    
!      do i=1,(p+1)*(q+1)*(w+1)
!        dRdxi(i,1) = dRdxi(i,1)/sumtot - (R(i)*sumxi  )/(sumtot**2)
!        dRdxi(i,2) = dRdxi(i,2)/sumtot - (R(i)*sumeta )/(sumtot**2)
!        dRdxi(i,3) = dRdxi(i,3)/sumtot - (R(i)*sumzeta)/(sumtot**2) 
!        R(i) = R(i)/sumtot
!      end do
!    
!  
!      !Gradient of mapping from parameter space to physical space
!      loc_num = 0
!      count = (p+1)*(q+1)*(w+1)+1
!      do k=0,w
!        do j=0,q
!            do i=0,p
!                
!                loc_num = loc_num + 1
!                count = count - 1
!                
!                dxdxi(1,1) = dxdxi(1,1) + dRdxi(loc_num,1)*(Points(IEN(nel,loc_num),1) + udisp(loc_num*3-2,1))
!                dxdxi(2,1) = dxdxi(2,1) + dRdxi(loc_num,1)*(Points(IEN(nel,loc_num),2) + udisp(loc_num*3-1,1))
!                dxdxi(3,1) = dxdxi(3,1) + dRdxi(loc_num,1)*(Points(IEN(nel,loc_num),3) + udisp(loc_num*3  ,1))
!                
!                dxdxi(1,2) = dxdxi(1,2) + dRdxi(loc_num,2)*(Points(IEN(nel,loc_num),1) + udisp(loc_num*3-2,1))
!                dxdxi(2,2) = dxdxi(2,2) + dRdxi(loc_num,2)*(Points(IEN(nel,loc_num),2) + udisp(loc_num*3-1,1))
!                dxdxi(3,2) = dxdxi(3,2) + dRdxi(loc_num,2)*(Points(IEN(nel,loc_num),3) + udisp(loc_num*3  ,1))
!                
!                dxdxi(1,3) = dxdxi(1,3) + dRdxi(loc_num,3)*(Points(IEN(nel,loc_num),1) + udisp(loc_num*3-2,1))
!                dxdxi(2,3) = dxdxi(2,3) + dRdxi(loc_num,3)*(Points(IEN(nel,loc_num),2) + udisp(loc_num*3-1,1))
!                dxdxi(3,3) = dxdxi(3,3) + dRdxi(loc_num,3)*(Points(IEN(nel,loc_num),3) + udisp(loc_num*3  ,1))
!                
!            end do
!        end do
!      end do
!    
!
!      call fdeterm(dxdxi,detj2)
!    
!      !Inverse of the gradient
!      temp1 = 0.0d0
!      dxidx = dxdxi
!      call gaussj(dxidx,nds,temp1,errorflag)
!    
!      dRdx = matmul(dRdxi,dxidx)
!    
!      detj = detj1*detj2 
!    
!      !detj = abs(detj)
!    
!      !Gradient mapping from parent element to parameter space
!      dxi_dtildexi(1,1) = 1.0d0/2.0d0
!      dxi_dtildexi(2,2) = 1.0d0/2.0d0
!      dxi_dtildexi(3,3) = 1.0d0/2.0d0
!    
!    
!      Jac = 0.0d0
!      do aa=1,nds
!        do bb=1,nds
!            do cc=1,nds
!                Jac(aa,bb) = Jac(aa,bb)+dxdxi(aa,cc)*dxi_dtildexi(cc,bb)
!            end do
!        end do
!      end do
!    
!
!      continue
!    
!      end subroutine ShapeFuncEAS
