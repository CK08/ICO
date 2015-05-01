!----------------------------------------------------------------------------------------------
!
! Material Subroutine for Isotropic Plasticity
!
! Input: iprops - number of properties
!        nodedof - degrees-of-freedom of each node
!        nnodes - nuber of nodes per element
!        props - property array os size iprops
!        B - Compatible strain-displacement matrix
!        Balpha - Enhanced strain-displacement matrix
!        BLines - number of lines of the strain-displacement operator
!        ddisp - displacement vector
!        nalpha - number of enhanced parameters
!        dalpha - increment in internal enhanced variables
!        
! Output: dstress - stress increment vector at the current ip
!         dstrain - strain increment vector at the current ip
!         matD - 3D consistent constitutive matrix
!         
! Updates: stress - total stress vector at the current ip
!          strain - total strain vector at the current ip
!          hard - plstic hardening
!
!----------------------------------------------------------------------------------------------

subroutine MatPlastic3DEAS(iprops,nodedof,nnodes,props,B,Balpha,Blines,nalpha,ddisp,dalpha,stress,strain,dstrain,dstress,hard,matD)
  
    implicit none

    real(8),parameter::YieldTol=1.0d-10

    integer(4),parameter::ndi=3,nshr=3,ntens=6
    
    integer(4),intent(IN)::iprops,nodedof,nnodes,Blines,nalpha
    real(8),dimension(iprops),intent(IN)::props
    real(8),dimension(Blines,nodedof*nnodes),intent(IN)::B
    real(8),dimension(nodedof*nnodes,1),intent(IN)::ddisp
    real(8),dimension(1,1,ntens),intent(INOUT)::stress,strain
    real(8),dimension(1,1,ntens),intent(INOUT)::dstress,dstrain
    real(8),dimension(1,1),intent(INOUT)::hard
    real(8),dimension(ntens,ntens),intent(INOUT)::matD
    
    real(8),dimension(nalpha,1),intent(IN)::dalpha
    real(8),dimension(Blines,nalpha),intent(IN)::Balpha
    
    integer(4)::i,j,k,knewton,k1,k2,errorflag

    real(8),dimension(ntens,1)::stressold,stressin,destrain
    real(8),dimension(ntens)::flow,dpstrain,store

    real(8)::EMod,ve,esp,SY0,H,Diff
    real(8)::eg,ek,SVM,YieldFunc,SYield,SHydro,hard0,res,dp
    real(8)::effg,effg2,effg3,efflam,effhrd

    real(8)::teste
    integer(4)::count,ilixo

    real(8),dimension(ntens,1)::DEP,aux1

    real(8)::bbar,norma,J2
    real(8),dimension(1,1)::bDa,HardPar
    real(8),dimension(ntens,1)::a,Db,dFdep
    real(8),dimension(1,ntens)::bD
    real(8),dimension(ntens,ntens)::Dba,DbaD,matDt,mdiff,matdelas

    real(8),dimension(ntens)::a1,a2,a3,crl
    real(8)::JJ2

    real(8),dimension(1,1)::den
    real(8),dimension(ntens,ntens)::matQ,Qinv,dadS,matR,aux66,aux66b
    real(8),dimension(ntens,ntens)::matdcons,matDinc

    !For consistent tangent stiffness
    real(8),dimension(ntens,ntens)::Bsigma,matDct,matDalg,mataux
    real(8),dimension(ntens,ntens)::matinv,Ident
    !    real(8),dimension(3,3)::TCL
    real(8),dimension(ntens,1)::dstrainenh
    real(8),dimension(ntens)::temp1

    real(8)::dgama,dilat,rnorm,rdumi
    real(8),dimension(6,1)::devsten,devdef,ted,def
    real(8)::rteta,rtetb,cons1,cons2,conp5,co113
    real(8),dimension(6)::ti,nor
    
    real(8),dimension(ntens,1)::temps1,temps2,temps1a

    !Initialise variables -----
    Stressold=0.0d0

    !Material Properties ------------------------
    EMod=props(1)
    ve=props(2)
    SY0=props(5)

    !Shear Modulus ---------------------------------
    EG=Emod/(2.0d0*(1.0d0+ve))
      
    !Bulk Modulus for 3D Applications --------------
    ek=EMod/(3.0d0*(1.0d0-2.0d0*ve))

    call MatLinearElastic3D(Emod,ve,MatDelas)

    !Save stress at the begining -------------------
    do j=1,ntens
        stressold(j,1)=stress(1,1,j)
    end do

    !Increment in strain --------------------------
    temps1 = 0.0d0
    temps1  = matmul(B,ddisp) 
    
    temps1a = 0.0d0
    temps1a = matmul(Balpha,dalpha)
    
    temps1 = temps1 + temps1a
    
    !Elastic increment in stress ------------------
    temps2=matmul(MatDelas,temps1)

    !Trial Stress ---------------------------------
    do j=1,ntens
        dstrain(1,1,j) = temps1(j,1)
        dstress(1,1,j) = temps2(j,1)
        stress(1,1,j)=stress(1,1,j)+dstress(1,1,j)
        strain(1,1,j)=strain(1,1,j)+dstrain(1,1,j)
    end do
       
    !Equivalent von Mises Stress ------------------
    do j=1,ntens
        stressin(j,1)=stress(1,1,j)
    end do

    call vonMisesStress(stressin,ntens,SVM)

    !Determine if Actively Yielding ---------------
    SYield=SY0+hard(1,1)
    YieldFunc=SVM-SYield

    continue

    !Yielding Condition NOT satisfied -------------
    count=0
    store=0.0d0
    if (YieldFunc .gt. 0.0d0)then

    do while(YieldFunc .gt. YieldTol)

        !Hydrostatic Pressure ------------------------------------------------
        SHydro=0.0d0
        do j=1,ndi
            SHydro=SHydro + stress(1,1,j)
        end do

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
        !---------------------------------------------------------------------
        
        !Determination of plastic multiplier and new yield stress ------------            
        hard0=hard(1,1)
        dp=0.0d0
        res=0.0d0
        SYield=SY0+hard0

        do knewton=1,15
            teste=3.0d0*EG*dp
            res=SVM-3.0d0*EG*dp-SYield
            dp=dp+res/(3.0d0*EG+H)
            call MatIsoHardening(iprops,props,6,SYield,H)
            hard(1,1)=hard0+H*dp
            SYield=SY0+hard(1,1)
            if(abs(res) .lt. 1.D-5) goto 10
        end do

        write(*,*)'Plasticity algorithm failed to converge'
            
    10       continue
     
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

        do j=1,ntens
            stress(1,1,j)=stress(1,1,j) - DEP(j,1)
        end do
        !---------------------------------------------------------------------

        !Yield surface verification ------------------------------------------
        do j=1,ntens
            stressin(j,1)=stress(1,1,j)
        end do
    !        
        call vonMisesStress(stressin,ntens,SVM)

        YieldFunc=SVM-SYield
        count=count+1
        continue
        !---------------------------------------------------------------------
     
    end do

    continue

    !Inconsistant tangent matrix --------
    a=0.0d0
    ident=0.0d0
    do k1=1,ntens
        a(k1,1)=flow(k1)
        ident(k1,k1)=1.0d0
    end do

    call MatIsoHardening(iprops,props,6,SYield,H)

    aux66=matmul(matmul(a,transpose(a)),matdelas)
    den=matmul(matmul(transpose(a),matdelas),a)

    aux66=aux66/(den(1,1) + H)
    aux66b=Ident-aux66

    matDinc=0.0d0
    matDinc=matmul(matdelas,aux66b)

    matd=matDinc
    continue


    !Consistent tangent matrix ------
    dadS=0.0d0
    Ident=0.0d0
    do k1=1,ndi
        do k2=1,ndi
            dadS(k1,k2)=-1.0d0
        end do
        dadS(k1,k1)=2.0d0
        Ident(k1,k1)=1.0d0
    end do
    do k1=ndi+1, ntens
        dadS(k1,k1)=6.0d0
        Ident(k1,k1)=1.0d0
    end do
            
    dadS=dadS/(2.0d0*SVM)
                
    a=0.0d0
    do k1=1,ntens
        do k2=1,ntens
            dadS(k1,k2)=dadS(k1,k2)-flow(k1)*flow(k2)/svm
        end do
        a(k1,1)=flow(k1)
    end do
                
    matQ=(Ident + dp*matmul(matdelas,dadS))
                
    Qinv=0.0d0
    Qinv=matQ
    call gaussj(Qinv,6,temp1,ilixo)
    !call FINDInv(matQ, Qinv, ntens, errorflag)

    matR=0.0d0
    matR=matmul(Qinv,matdelas)

    aux66 = matmul(a,transpose(a))
    aux66b = matmul(aux66,matR)

    den=0.0d0
    den=matmul(matmul(transpose(a),matR),a)

    call MatIsoHardening(iprops,props,6,SYield,H)

    aux66=0.0d0
    aux66=Ident - aux66b/(den(1,1)+H)

    matDCons=matmul(matR,aux66)

    matd=matDCons
    continue

    !Volumetric dilatation
    do k1=1, ntens
        def(k1,1) = strain(1,1,k1) !+ dstrain(k1,1)
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
!----------------------------------------------------------------------------------------------