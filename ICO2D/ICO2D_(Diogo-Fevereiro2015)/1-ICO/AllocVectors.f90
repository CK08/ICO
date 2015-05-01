!-----------------------------------------------------------------------
! 
! Subroutine to initialise some necessary vectors and matrices
!
!-----------------------------------------------------------------------

subroutine AllocVectors()
    
    use Mod_Variables
    implicit none
    
    !Displacement increment
    allocate(dDisp(tnodes*nds,1))
    dDisp=0.0d0

    !Increment of the displacement increment
    allocate(dddisp(tnodes*nds,1))
    ddDisp=0.0d0
    
    !Total converged displacement
    allocate(u(tnodes*nds,1))
    u=0.0d0
    
    !Reaction Forces
    allocate(Reac(tnodes*nds,1))
    Reac=0.0d0
    
    !Another Stiffness matrix for calculating the reaction forces
    allocate(KT(tnodes*nds,tnodes*nds))
    KT=0.0d0
    
    firstinc = .true.
    
    !Some Arrays Initializations
    !dispdof = 0.0d0
    SEnergyConv = 0.0d0
    FintLM = 0.0d0
    FintLMConv = 0.0d0
    
end subroutine