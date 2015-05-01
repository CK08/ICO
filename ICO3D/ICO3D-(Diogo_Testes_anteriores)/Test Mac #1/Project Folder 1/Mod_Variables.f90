!---------------------------------------------------------------
!
! Global Variables Module
!
!---------------------------------------------------------------

Module Mod_Variables

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