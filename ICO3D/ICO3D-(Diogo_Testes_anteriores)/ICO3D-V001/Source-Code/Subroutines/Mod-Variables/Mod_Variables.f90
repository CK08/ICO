!---------------------------------------------------------------
!
! Global Variables Module
!
!---------------------------------------------------------------

Module Mod_Variables

    !-----------------------------------------------------------------------------------------------------------------------
    ! English::
    ! The following global variables have different functions in the main program, such as:         
    !        - RToler (Convergence Residual Tolerance)
    !        - nds (Number of dimensions in the problem)
    !        - incmax (Maximum number of increments in the analysis)
    !        - itermax (Maximum number of iterations per increment in the analysis)
    !        - p, q, w (degree in each direction: x, y e z, respectively)        
    !        - ncpx, ncpy, ncpz (Number of control points in each direction - x, y e z, respectively)
    !        - nnodes (Number of nodes)
    !        - nelems (Number of elements)
    !        - nshpl (Number of shape functions)
    !        - closed_u, closed_v, closed_w (Closed Knot Vectors - Defined as in IGA book (Hughes) but currently NOT IN USE!)
    !        - elcode (String with element type, for example: 'Hex8')
    !        - nbc (Number of boundary conditions)
    !        - ndisp (Number of degrees of freedom affected by prescribied displacements)
    !
    !        - npi (Number of integration points per element)
    !        - SEnergy (Strain Energy)
    !        - SEnergyConv (Strain Energy Converged)
    !        - props (Material properties vector. The order of the vector is given in table 1 of ICO User Manual)
    !        - iprops (Dimension of "props" vector)
    !        - npenal (Penalty Method applied to DupNode [Duplicated Nodes] - in Assembly3D.f90)
    !        - ipenal (Penalty Method applied to DupNode [Duplicated Nodes] - in Assembly3D.f90)
    !        - bc (Boundary Conditions, not for prescribed displacements)
    !        - iper ("Nao me lembro" - J.F. Caseiro)
    !        - dir_bc (not used)
    !        - ld_el (not used)
    !        - load (Matrix with point loads)
    !        - dispBC (Prescribed Displacement Boundary Conditions)
    !        - BW (Weighted control points)
    !        - EY (not used - possibly assign one material for each element [EY - Young modulus])
    !        - NU (not used - possibly assign one material for each element [NU - Poisson's coefficient])
    !        - Fext (External Forces)
    !        - FextEq (External Forces Reduced System)
    !        - Fini (Initial Forces)
    !        - FInc (Incremental external force)
    !        - Fint (Internal Forces)
    !        - Res (External Forces - Internal Forces = Residual)
    !        - ImpDisp (Implicit/Prescribed Displacement - Array)
    !        - Kf (Global Stiffeness Matrix)
    !        - KfEq (Reduced Global Stiffeness Matrix)
    !        - KfEqInv (Inverse of the Reduced Global Stiffeness Matrix - currently not used [slowwww])
    !        - dispeqv (Increment of the Increment of the Reduced Displacement)
    !        - RedCount (Auxiliar array to construct the global displacement array from the reduced one)
    !        - dDisp (Displacement Increment)
    !        - ddDisp (Incremente of the displacement increment)
    !        - u (Displacement)
    !        - stress (STRESS:: element number, integration point, stress component)
    !        - dstress (Stress Increment)
    !        - stress_conv (Converged Stress)
    !        - strain (Strain:: equal to stress)
    !        - dstrain (Strain Increment)
    !        - strain_conv (Converged Strain)
    !        - hard (Hardening:: element number, integration point)
    !        - hard_conv (Converged Hardening)
    !        - TmatD (Constitutive (D) matrix backup) - Allocated in ReadInFile
    !        - laxis (local axis:: element, integration point, axis)
    !        - laxisconv (Converged local axis)
    !        - PMesh (Control Points Coordinates in the global numbering system [i.e. GiD, contact])
    !        - stpfrac (Fraction of applied prescribed displacement in the increment)
    !        - Weight (Control Point Weights)
    !        - Points (Control Points Coordinates)
    !        - penalmeth (Penalty Method - it does something cool!!!)
    !        - nlgeom (bool: if True - accounts for GNL [Geometric Non-Linearity])
    !        - loaddof (Load with global numbering)
    !        - dispdof (Zero displacement BC with global numbering)
    !        - Eigen (Used to compute Eigen values of a 3x3 matrix)
    !        - vcp (Eigen Vectors)
    !        - vlp (Eigen Values) 
    !        
    !        MATERIALS - material properties relative dimensions 
    !        - ndi (Number of direct components) 
    !        - nshr (Number of Shear components)
    !        - ntens (ndi + nshr)
    !        
    !        KNOT VECTORS - 
    !        - u_knot, v_knot, w_knot (knot vectors in each direction: x, y e z, respectively)
    !        
    !        CONTROL POINTS - 
    !        - b_net (Control Points original coordenates matrix)
    !        - b_net_final (Control Points deformed coordenates matrix) 
    !        
    !        NURBS COORDENATES VECTORES 
    !        - INN (Array with the indices of the knots (NURBS coordenates).They are used to identify the knots at which the support of a function begins)
    !
    !        NODES PER ELEMENT VECTORS - 
    !        - IEN (Matrix with the connectivities for IG [IsoGeometric] elements - are given by the element number and knot vectors, first in the x direction, followed by y and z)
    !        - conn ("Nao serve para nada" - J.F. Caseiro)
    !
    !        GRAVITY -
    !        - gravity (Bool - if TRUE, gravity is considered in the analysis)
    !        - gravdir (Integer that takes values 1, 2 and 3 for gravity direction: x, y ou z, respectively)
    !        - gconst (Gravitacional Constant "g")
    !
    !        PRESSURE - 
    !        - ipress (Integer that defines the number of surfaces under pressure load)
    !        - ierror (error index)
    !        - tempc (read file - not important or is it?)
    !        - pressvec (Pressure Magnitude)
    !        - ftag (Element Face - read ICO User Manual)
    !        - pressure (Bool - True if pressure exists)
    !                
    !        B-BAR METHOD - 
    !        - IENb (Connectivities for B-Bar)
    !        - MatAf (B-Bar mixed method - zienkiewicz book [good luck ;)])
    !        - MatVf (B-Bar mixed method - zienkiewicz book [good luck ;)])
    !        - MatCf, MatCfT (B-Bar mixed method - zienkiewicz book [good luck ;)])
    !        - Kbar (Stiffeness after static condensation of B-Bar)
    !        - fbar (Load after static condensation of B-Bar)
    !        - nalpha (number of EAS enhanced parameters)
    !        - diag ()
    !        - nlin, ncol (Number of line and columns of the MatVf)
    !        - indx (index number)    
    !    
    ! Added by Diogo Cardoso, March 2014.
    !-----------------------------------------------------------------------------------------------------------------------
    real(8),parameter::RToler=0.0001d0
    
    integer(4)::nds,incmax,itermax
    integer(4)::p,q,w
    integer(4)::ncpx,ncpy,ncpz
    integer(4)::nnodes,nelems,nshpl
    integer(4)::closed_u,closed_v,closed_w
    
    character(16)::elcode
    
    integer(4)::nbc,ndisp
    
    !Material dimensions
    integer(4)::ndi    !Direct Components
    integer(4)::nshr   !Shear Components
    integer(4)::ntens  !Total
    
    integer(4)::npi    !Integration points per element
    
    real(8)::SEnergy,SEnergyConv
    
    integer(4)::iprops,npenal,ipenal
    real(8),dimension(:),allocatable::props
    
    integer(4),allocatable::bc(:,:,:,:),iper(:)
    
    real(8),allocatable::dir_bc(:,:),ld_el(:,:)
    
    !Knot vectors
    real(8),allocatable::u_knot(:),v_knot(:),w_knot(:)
    
    !Control points
    real(8),allocatable::b_net(:,:,:,:),b_net_final(:,:,:,:)
    
    real(8), dimension(:,:,:,:), allocatable::load,dispBC
    real(8), dimension(:,:,:,:),allocatable :: BW
    real(8), dimension(:), allocatable::EY
    real(8), dimension(:), allocatable:: NU
    
    !NURBS coordinate array
    integer(4),dimension(:,:),allocatable::INN
    
    !Element nodes array
    integer(4),dimension(:,:),allocatable::IEN,conn
    
    real(8), dimension(:,:), allocatable:: Fext,Fini,FInc,Res,FextEq,Fint,ImpDisp
    real(8), dimension(:,:), allocatable::Kf, KfEq, KfEqInv, dispeqv
    
    integer(4), dimension(:), allocatable::RedCount
    
    real(8), dimension(:,:), allocatable::dDisp, ddDisp,u
    
    real(8),dimension(:,:,:), allocatable::stress,dstress,strain,dstrain
    real(8),dimension(:,:,:), allocatable::stress_conv,strain_conv
    real(8),dimension(:,:), allocatable::hard,hard_conv
    real(8),dimension(:,:,:,:), allocatable:: TmatD,laxis,laxisconv
    
    real(8),dimension(:,:), allocatable::PMesh
    
    real(8),allocatable,dimension(:)::stpfrac,Weight
    real(8),allocatable,dimension(:,:)::Points
    
    integer(4),allocatable,dimension(:,:)::penalmeth
    
    !Gravity paramenters
    logical::gravity
    integer(4)::gravdir
    real(8)::gconst
    
    !Pressure paramenters
    integer(4)::ipress,ierror
    character(3)::tempc
    
    real(8),allocatable,dimension(:)::pressvec
    character(3),allocatable,dimension(:)::ftag
    logical::pressure
    
    logical::nlgeom
    
    real(8),allocatable,dimension(:)::loaddof,dispdof
    
    logical::Eigen
    real(8),allocatable,dimension(:,:)::vcp
    real(8),allocatable,dimension(:)::vlp
    
    !Allocations for B-bar method
    integer(4),dimension(:,:),allocatable::IENb
    real(8),dimension(:,:),allocatable::MatAf
    real(8),dimension(:,:),allocatable::MatVf
    real(8),dimension(:,:),allocatable::MatCf,MatCfT
    real(8),dimension(:,:),allocatable::Kbar, fbar
    
    integer(4)::nalpha
    
    real(8),dimension(:),allocatable::diag
    
    integer(4)::runner,nlin,ncol
    integer(4),dimension(:,:),allocatable::indx
    
end module
