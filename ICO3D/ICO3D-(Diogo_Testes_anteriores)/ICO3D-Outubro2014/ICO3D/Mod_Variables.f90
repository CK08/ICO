!---------------------------------------------------------------
!
! Global Variables Module
!
!---------------------------------------------------------------

Module Mod_Variables

    real(8),parameter::RToler=5.0d-3
    
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
    
    integer(4)::nptch !Number of patches
    
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
    
    real(8),dimension(:),allocatable::Reac
    real(8),dimension(:,:),allocatable::KTot
    
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
    
    !ANS
    real(8),dimension(:,:),allocatable::MatL1,MatL2
    real(8),dimension(:,:),allocatable::MatL3
    real(8),dimension(:,:),allocatable::TyPtANS
    
    !Allocations for Multi-Patch
    integer(4)::coin,tnodes
    integer(4)::pmax,qmax,wmax,npimax
    integer(4)::ncpxmax,ncpymax,ncpzmax
    integer(4),dimension(:),allocatable::MP_nds,MP_npi,MP_iprops
    integer(4),dimension(:),allocatable::MP_p,MP_q,MP_w
    integer(4),dimension(:),allocatable::MP_ncpx,MP_ncpy,MP_ncpz
    integer(4),dimension(:),allocatable::MP_closed_u,MP_closed_v,MP_closed_w
    
    real(8),dimension(:,:),allocatable::MP_uknot,MP_vknot,MP_wknot,MP_props
    real(8),dimension(:,:,:,:,:),allocatable::MP_b_net,MP_b_net_final
!    real(8),dimension(:,:),allocatable::AllPoints
    integer(4),dimension(:,:),allocatable::CoinCP
    integer(4),dimension(:,:),allocatable::MP_INN
    integer(4),dimension(:,:),allocatable::MP_IEN,MP_IEN_Temp
    
    character(16),dimension(:),allocatable::MP_elcode
    integer(4),allocatable,dimension(:)::bcdof
    real(8),dimension(:,:),allocatable::GCoords
    
    !Contact
    logical::PTS, MRigid
    
    !Master data
    integer(4)::p_mst, q_mst, ncpx_mst, ncpy_mst
    integer(4),dimension(:),allocatable::Conn_mst
    real(8),dimension(:),allocatable::u_knot_mst, v_knot_mst
    
    !Slave data
    integer(4)::p_slv, q_slv, ncpx_slv, ncpy_slv, ngrv_xi, ngrv_eta, islave
    integer(4),dimension(:),allocatable::Conn_slv
    real(8),dimension(:),allocatable::u_knot_slv, v_knot_slv, grv_xi, grv_eta
    
    real(8),dimension(:),allocatable::Lagrange, dLagrange
    real(8),dimension(:,:),allocatable::FintLM
    
    !Multibody contact
    integer(4)::npair
    integer(4),dimension(:),allocatable::MB_islave
    logical,dimension(:),allocatable::MB_Mrigid
    
    integer(4),dimension(:),allocatable::MB_p_mst, MB_q_mst
    integer(4),dimension(:),allocatable::MB_ncpx_mst,MB_ncpy_mst
    real(8),dimension(:,:),allocatable::MB_u_knot_mst, MB_v_knot_mst
    integer(4),dimension(:,:),allocatable::MB_conn_mst
    
    integer(4),dimension(:),allocatable::MB_p_slv, MB_q_slv
    integer(4),dimension(:),allocatable::MB_ncpx_slv,MB_ncpy_slv
    real(8),dimension(:,:),allocatable::MB_u_knot_slv, MB_v_knot_slv
    integer(4),dimension(:,:),allocatable::MB_conn_slv
    real(8), dimension(:,:),allocatable::GP_coords
    
end module
