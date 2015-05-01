!------------------------------------------------------------------------------------------------------------
!
! Subroutine for plasticity with isotropic hardening for plane strain analysis based on 'Non-linear 
! Finite Element Analysis of Solids and Structures Vol 1' by M.A. Crisfield
!
! Input: iprops - number of material properties
!        nodedof - number of dof's per control point
!        nnodes - number of control points
!        props - properties array
!        B - Strain-displacement matrix
!        Blines - number of lines of the B matrix
!        ddisp - displacement increment
! Input/Output: stress - stress array of the current integration point
!               strain - strain array of the current integration point
!               dstress - stress increment of the current integration point
!               dstrain - strain increment of the current integration point
!               hard - hardening undegone by the current integration point
!               matD - consitent constitutive matrix
!
!------------------------------------------------------------------------------------------------------------
subroutine MaterialPStrain(iprops,nodedof,nnodes,props,B,Blines,ddisp,stress,strain,dstrain,dstress,hard,matD)

    implicit none
    
    real(8),parameter::YieldTol=1.0d-5
    
    integer(4),parameter::ndi=2,nshr=1,ntens=3
    
    integer(4),intent(IN)::iprops,nodedof,nnodes,Blines
    real(8),dimension(iprops),intent(IN)::props
    real(8),dimension(Blines,nodedof*nnodes),intent(IN)::B
    real(8),dimension(nodedof*nnodes,1),intent(IN)::ddisp
    real(8),dimension(1,1,ntens),intent(INOUT)::stress,strain
    real(8),dimension(1,1,ntens),intent(INOUT)::dstress,dstrain
    real(8),dimension(1,1),intent(INOUT)::hard
    real(8),dimension(ntens,ntens),intent(INOUT)::matD
    
    real(8),dimension(ntens,1)::temps,temps2
    
    integer(4)::i,j,k,knewton,k1,k2,errorflag,counter
    
    real(8),dimension(ntens,1)::stressold,stressin,destrain
    real(8),dimension(ntens)::flow,dpstrain,store
    
    real(8)::EMod,ve,esp,SY0,H,Diff
    real(8)::eg,ek,SVM,YieldFunc,SYield,SHydro,hard0,res,dp
    real(8)::effg,effg2,effg3,efflam,effhrd
    
    real(8),dimension(ntens,1)::DEP,aux1,rr,stressrate
    real(8),dimension(ntens,ntens)::dadS,matQ,Qinv,matR
    real(8),dimension(ntens,ntens)::matDInc,MatDCons
    
    real(8)::param,dpr,YFS,deltadp,YieldOUT
    real(8),dimension(1,1)::num,den
    real(8),dimension(1,3)::aux13,aux13b
    real(8),dimension(3,1)::aux31,aux31b
    real(8),dimension(3,3)::aux33,aux33b
    
    real(8)::bbar,norma,J2
    real(8),dimension(1,1)::bDa,HardPar
    real(8),dimension(ntens,1)::a,Db,dFdep
    real(8),dimension(1,ntens)::bD
    real(8),dimension(ntens,ntens)::Dba,DbaD,matDt,mdiff,matdelas
    
    real(8),dimension(ntens)::a1,a2,a3,crl
    real(8)::JJ2,normcrl
    
    real(8)::strzz,flowz,dpstrainz,depz
    
    !For consistent tangent stiffness
    real(8),dimension(ntens,ntens)::Bsigma,matDct,matDalg,mataux,matinv,Ident
    
    real(8)::dgama,dilat,rnorm,rdumi
    real(8),dimension(3,1)::devsten,devdef,ted,def
    real(8)::rteta,rtetb,cons1,cons2,conp5,co113
    real(8),dimension(3)::ti,nor
    
!    real(8),dimension(3,3)::TCL
    real(8),dimension(3)::ten,defin
    real(8),dimension(2,2)::fm0,fm1
    

    !Initialise variables -----
    Stressold=0.0d0
    
    ident=0.0d0
    do k1=1,ntens
        ident(k1,k1)=1.0d0
    end do
    
    !Material Properties ------------------------
    EMod=props(1)
    ve=props(2)
    SY0=props(5)
    H = 0.0d0

    !Shear Modulus ---------------------------------
    EG=Emod/(2.0d0*(1.0d0+ve))
    
    !Bulk Modulus for Plane Strain------------------
    ek=EMod/(2.0d0*(1.0d0+ve)*(1.0d0-2.0d0*ve))

    !Elastic Constitutive matrix
    call MatLinearElasticDP(Emod,ve,MatDelas)
    matD=matdelas
    
    !Save stress at the begining -------------------
    do j=1,ntens
        stressold(j,1)=stress(1,1,j)
    end do

    !Increment in strain --------------------------
    temps=matmul(B,ddisp)
    !dstrain=matmul(B,ddisp)
    
    do k1=1,ntens
        dstrain(1,1,k1) = temps(k1,1)
    end do
    
   
    !Elastic increment in stress ------------------
    !dstress=matmul(MatDelas,dstrain)
    temps2 = matmul(MatDelas,temps)
    matd=matdelas
    
    do k1=1,ntens
        dstress(1,1,k1) = temps2(k1,1)
    end do

    !Trial Stress ---------------------------------
    do j=1,ntens
        stress(1,1,j)=stress(1,1,j)+dstress(1,1,j)
    end do
    
    !Store stress before returning to yield surface
    do j=1,ntens
        stressold(j,1)=stress(1,1,j)
    end do
    
    !Stress along the out-of-plane axis -----------
    strzz=0.0d0
    strzz=ve*(stress(1,1,1) + stress(1,1,2))
    
    !Equivalent von Mises Stress ------------------
    do j=1,ntens
        stressin(j,1)=stress(1,1,j)
    end do
    
    SVM=(stress(1,1,1)-stress(1,1,2))**2.0d0+(stress(1,1,2)-strzz)**2.0d0+(stress(1,1,1)-strzz)**2.0d0
    SVM=SVM+6.0d0*(stress(1,1,3)**2.0d0)
    SVM=SVM/2.0d0
    SVM=sqrt(SVM)

    !Determine if Actively Yielding ---------------
    SYield=SY0+hard(1,1)
    YieldFunc=SVM-SYield
    
    continue
    
    !Yielding Condition NOT satisfied -------------
    store=0.0d0
    counter=0
    
        if (YieldFunc .gt. 0.0d0)then
    
        do while(YieldFunc .gt. YieldTol)
            
            !write(*,*)'In plasticity'
            
            !Hydrostatic Pressure ------------------------------------------------
            SHydro=0.0d0
            do j=1,ndi
                SHydro=SHydro + stress(1,1,j)
            end do
            
            SHydro=SHydro + strzz
            
            SHydro=SHydro/3.0d0
            !---------------------------------------------------------------------
                    
            !Flow direction from the trial surface -------------------------------
            flow=0.0d0
            do j=1,ndi
                flow(j)=(stress(1,1,j)-SHydro)/SVM
            end do
            do j=ndi+1,ntens
                flow(j)=stress(1,1,j)/SVM !*2.0d0
            end do
            
            flowz=0.0d0
            flowz=(strzz-SHydro)/SVM
            !---------------------------------------------------------------------
            
            !Determination of plastic multiplier and new yield stress ------------            
            hard0=hard(1,1)
            dp=0.0d0
            res=0.0d0
            SYield=SY0+hard0

            do knewton=1,15
                res=SVM-3.0d0*EG*dp-SYield
                dp=dp+res/(3.0d0*EG+H)
                !call MatIsoHardening(iprops,props,6,SYield,H)
                hard(1,1)=hard0+H*dp
                SYield=SY0+hard(1,1)
                if(abs(res) .lt. 1.D-5) goto 10
            end do

            write(*,*)'Plasticity algorithm failed to converge'
                
    10      continue
         
            !Increment in plastic strain -----------------------------------------
            dpstrain=0.0d0
            do j=1,ndi
                !stress(1,1,j)=flow(j)*syield + shydro
                dpstrain(j)= 3.0d0/2.0d0*flow(j)*dp
            end do
            do j=ndi+1,ntens
                !stress(1,1,j)=flow(j)*syield
                dpstrain(j)= 3.0d0*flow(j)*dp
            end do
            
            dpstrainz=3.0d0/2.0d0*flowz*dp
            store=store+dpstrain
            !---------------------------------------------------------------------
        
            !Stress update -------------------------------------------------------
            DEP=0.0d0
            aux1=0.0d0
            do j=1,ntens
                aux1(j,1)=dpstrain(j)
            end do
        
            !DEP=matmul(matD,aux1)
            DEP=matmul(matDelas,aux1)
            
            DEP(1,1) = DEP(1,1) + Emod*ve/((1.0d0+ve)*(1.0d0-2.0d0*ve))*dpstrainz
            DEP(2,1) = DEP(2,1) + Emod*ve/((1.0d0+ve)*(1.0d0-2.0d0*ve))*dpstrainz
            
            depz=Emod*ve/((1.0d0+ve)*(1.0d0-2.0d0*ve))*(aux1(1,1) + aux1(2,1)) + Emod*(1.0d0-ve)/((1.0d0+ve)*(1.0d0-2.0d0*ve))*dpstrainz
            
            do j=1,ntens
                stress(1,1,j)=stress(1,1,j) - DEP(j,1)
            end do
            
            strzz=strzz - depz
            !---------------------------------------------------------------------
        
            !Yield surface verification ------------------------------------------
            SVM=(stress(1,1,1)-stress(1,1,2))**2.0d0+(stress(1,1,2)-strzz)**2.0d0+(stress(1,1,1)-strzz)**2.0d0
            SVM=SVM+6.0d0*(stress(1,1,3)**2.0d0)
            SVM=SVM/2.0d0
            SVM=sqrt(SVM)
       
            YieldFunc=SVM-SYield
            continue
            !---------------------------------------------------------------------
         
        end do


        !Volumetric dilatation
        do k1=1, ntens
            def(k1,1) = strain(1,1,k1) + dstrain(1,1,k1)
        end do
        
        dilat=0.0d0
        do k1=1,ndi
            dilat = dilat + def(k1,1) 
        end do
        dilat=dilat/3.0d0
        
        devdef=0.0d0
        do k1=1,ndi
            devdef(k1,1)=def(k1,1) - dilat
        end do
        
        do k1=ndi+1,ntens
            devdef(k1,1) = def(k1,1)/2.0d0
        end do
        
        rnorm=0.0d0
        do k1=1,ndi
            rdumi = 2.0d0*EG*devdef(k1,1)
            rnorm = rnorm + rdumi*rdumi
        end do
        do k1=ndi+1,ntens
            rdumi = 2.0d0*EG*devdef(k1,1)
            rnorm = rnorm + rdumi*rdumi*2.0d0
        end do
        
        rnorm = dsqrt(rnorm)
        
        rdumi=1.0d0/rnorm

        ti=0.0d0
        do k1=1,ndi
            ti(k1)=1.0d0
        end do
        dgama=dp/sqrt(2.0d0/3.0d0)
        nor = flow/(dsqrt(2.0d0/3.0d0))
        rteta = 1.0d0 - (2.0d0*EG*rdumi*dgama)
        rtetb=(1.0d0/(1.0d0+(H/(3.0d0*EG))))-(1.0d0-rteta)
        
        cons1=2.0d0*EG*rteta
        cons2=2.0d0*EG*rtetb
        conp5=cons1*0.5d00
        
        matd=0.0d0
        do k1=1,ndi
            matD(k1,k1)=cons1
        enddo

        do k1=ndi+1,ntens
            matD(k1,k1)=conp5
        enddo
        
        co113=EK-cons1/3.0d0
        do k1=1,ntens
            do k2=k1,ntens
                matD(k1,k2)=matD(k1,k2)+co113*ti(k1)*ti(k2)-cons2*nor(k1)*nor(k2)
            enddo
        enddo

        do k1=1,ntens
            do k2=1,k1-1
                matD(k1,k2)=matD(k2,k1)
            enddo
        enddo
        
        continue

    endif
    
    continue
    
end subroutine