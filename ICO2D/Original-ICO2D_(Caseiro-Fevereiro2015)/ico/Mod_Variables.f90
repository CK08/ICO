!---------------------------------------------------------------
!
! Global Variables Module
!
!---------------------------------------------------------------

Module Mod_Variables

    real(8),parameter::RToler=1.0d-7
    
    integer(4)::nds,p,q,ncpx,ncpy
    integer(4)::nnodes,nelems,nshpl,telems
    integer(4)::closed_u,closed_v
    
    character(16)::elcode
    
    integer(4)::nbc,ndisp,tnodes,nstp,itermax,incmax
    real(8)::sumImpDisp
    
    integer(4),dimension(:),allocatable::STP_nbc,STP_ndisp,STP_nload,STP_dload
    integer(4),dimension(:),allocatable::STP_incmax,STP_itermax
    
    integer(4),dimension(:,:),allocatable::STP_bcdof
    real(8),dimension(:,:),allocatable::STP_dispdof,STP_loaddof
    real(8),dimension(:,:,:),allocatable::STP_distload
    
    !Material dimensions
    integer(4)::ndi    !Direct Components
    integer(4)::nshr   !Shear Components
    integer(4)::ntens  !Total
    
    integer(4)::npi    !Integration points per element
    
    real(8)::SEnergy,SEnergyConv
    
    integer(4)::iprops,npenal,ipenal
    real(8),dimension(:),allocatable::props
    
    integer(4),allocatable::bc(:,:,:),iper(:)
    
    real(8),allocatable::dir_bc(:,:),ld_el(:,:)
    real(8),allocatable::u_knot(:),v_knot(:),b_net(:,:,:),b_net_final(:,:,:)
    
    real(8), dimension(:,:,:), allocatable::load,dispBC
    real(8), dimension(:,:,:),allocatable :: BW
    real(8), dimension(:), allocatable::EY
    real(8), dimension(:), allocatable:: NU
    
    !NURBS coordinate array
    integer(4),dimension(:,:),allocatable::INN
    
    !Element nodes array
    integer(4),dimension(:,:),allocatable::IEN,conn
    
    real(8), dimension(:,:), allocatable:: Fext,Fini,FInc,Res,FextEq,Fint,ImpDisp,Fct
    real(8), dimension(:,:), allocatable::Kf, KfEq, KfEqInv, dispeqv,KT
    
    integer(4), dimension(:), allocatable::RedCount
    
    real(8), dimension(:,:), allocatable::dDisp, ddDisp,u,Reac
    
    real(8),dimension(:,:,:), allocatable::stress,dstress,strain,dstrain
    real(8),dimension(:,:,:), allocatable::stress_conv,strain_conv
    real(8),dimension(:,:), allocatable::hard,hard_conv
    real(8),dimension(:,:,:,:), allocatable::laxis,laxis_conv
    real(8),dimension(:,:,:,:), allocatable:: TmatD
    
    real(8),dimension(:,:), allocatable::PMesh,PMesh0
    
    real(8),allocatable,dimension(:)::stpfrac
    
    integer(4),allocatable,dimension(:,:)::penalmeth
    
    real(8),allocatable::Points(:,:),Weights(:)
    real(8),allocatable::Points0(:,:)
    
    logical::nlgeom,firstinc,Multistep
    
    real(8),allocatable,dimension(:,:)::GP_Coords
    
    !Multipatch parameters
    integer(4)::nptch
    integer(4),allocatable,dimension(:,:)::order,ncpptch,statptch
    real(8),allocatable,dimension(:,:)::MP_u_knot,MP_v_knot
    real(8),allocatable,dimension(:,:,:,:)::MP_B_net,MP_B_net_final
    character(16),allocatable,dimension(:)::MP_elcode
    integer(4),allocatable,dimension(:)::MP_npi,MP_nnodes, MP_nelems
    real(8),allocatable,dimension(:,:)::MP_props
    integer(4),allocatable,dimension(:,:)::MP_conn
    
    integer(4),dimension(:),allocatable::bcdof
    real(8),dimension(:),allocatable::loaddof,dispdof,dispdofold
    real(8),dimension(:,:),allocatable::distload
    
    !Contact Mechanics
    logical::ContactPTS,contactGPTS
    character(16)::ContactMethod
    integer(4)::p_s,ncp_s,p_m,ncp_m, islv, SurfGP
    integer(4),dimension(:),allocatable::conn_slave,conn_master
    real(8),dimension(:),allocatable::u_knot_slave,u_knot_master
    real(8),dimension(:,:),allocatable::b_slave,b_master
    
    integer(4),dimension(:,:),allocatable:: INN_s,IEN_s
    real(8),dimension(:,:),allocatable::Points_slv
    real(8),dimension(:),allocatable::Weights_slv,Elem_area,Elem_CF,Elem_str
    real(8),dimension(:),allocatable::Elem_str_plot,Elem_xi_plot
    real(8),dimension(:),allocatable::CFs,CollPoint_CF
    
    real(8),dimension(:,:),allocatable::GCoords,plot_xy
    
    real(8),dimension(:,:),allocatable::FintLM,FintLMConv,FintLMi
    real(8),dimension(:),allocatable::xibs
    
    real(8),dimension(:),allocatable::Lagrange,LagrangeConv,dLagrange,ddLagrange,Lagmx,lagmy
    real(8), dimension(:,:), allocatable::KCont
    
    integer(4)::ngrv
    real(8)::sumLag
    real(8),dimension(:),allocatable::grev, CP_area, PTS_Stress, Press
    real(8),dimension(:),allocatable::initial_gap,final_gap
    integer(4),dimension(:),allocatable::PTS_active,Conn_El
    
    real(8),allocatable,dimension(:,:)::Ns_Conv
    
end module