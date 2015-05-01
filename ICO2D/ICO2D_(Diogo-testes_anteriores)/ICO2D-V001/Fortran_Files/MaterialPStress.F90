!------------------------------------------------------------------------------------------------------------
!
! Subroutine for plasticity with isotropic hardening for plane stress analysis based on 'Non-linear 
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

!subroutine EPMaterialPStress(ndi,nshr,ntens,iprops,nodedof,nnodes,props,B,Blines,ddisp,stress,strain,dstrain,dstress,hard,matD)
subroutine MaterialPStress(iprops,nodedof,nnodes,props,B,Blines,ddisp,stress,strain,dstrain,dstress,hard,matD)

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
    real(8),dimension(3)::temp1
    
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
    
    !For consistent tangent stiffness
    real(8),dimension(ntens,ntens)::Bsigma,matDct,matDalg,mataux,matinv,Ident
    
    real(8)::dgama,dilat,rnorm,rdumi,strainzz
    real(8),dimension(3,1)::devsten,devdef,ted,def
    real(8)::rteta,rtetb,cons1,cons2,conp5,co113
    real(8),dimension(3)::ti,nor
    
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
    
    !Shear Modulus ---------------------------------
    EG=Emod/(2.0d0*(1.0d0+ve))
    
    !Bulk Modulus for Plane Stress------------------
    ek=EMod/(2.0d0*(1.0d0-ve))

    !Elastic Constitutive matrix
    call MatLinearElasticTP(Emod,ve,MatDelas)
    
    !Save stress at the begining -------------------
    do j=1,ntens
        stressold(j,1)=stress(1,1,j)
    end do

    !Increment in strain --------------------------
    temps=matmul(B,ddisp)
    do k1=1,ntens
        dstrain(1,1,k1) = temps(k1,1)
    end do
      
    !Elastic increment in stress ------------------
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
    store=0.0d0
    counter=0
    
    if (YieldFunc .gt. 0.0d0)then
        
        deltadp=0.0d0
        dpr=0.0d0
        YFS=YieldFunc
        
        hard0=hard(1,1)
        
        YieldOUT=100.0d0
        
        do while(YieldOUT .gt. YieldTol)
            
            counter=counter+1
            
            !Flow direction (from Crisfield)
            flow(1)=2.0d0*stress(1,1,1)-stress(1,1,2)
            flow(2)=2.0d0*stress(1,1,2)-stress(1,1,1)
            flow(3)=6.0d0*stress(1,1,3)
            
            flow=flow/(2.0d0*SVM)
            
            continue

            !Determine initial plastic multiplier and Yield stress
            if (counter==1)then

                SYield=SY0+hard0
            
                param=0.0d0
                do k1=1,3
                    do k2=1,3
                        param=param + flow(k1)*matdelas(k1,k2)*flow(k2)
                    end do
                end do
                
                call MatIsoHardening(iprops,props,6,SYield,H)
                
                dp=YieldFunc/(param + H)
                hard(1,1)=hard0+H*dp
                SYield=SY0+hard(1,1)

            end if
            
            !Increment in plastic strain -----------------------------------------
            dpstrain=0.0d0
            do j=1,ndi
                !stress(1,1,j)=flow(j)*syield + shydro
                dpstrain(j)= 3.0d0/2.0d0*flow(j)*dp
                !dpstrain(j)= flow(j)*dp !Hinton and Owen
            end do
            do j=ndi+1,ntens
                !stress(1,1,j)=flow(j)*syield
                dpstrain(j)= 3.0d0*flow(j)*dp
                !dpstrain(j)= flow(j)*dp !Hinton and Owen
            end do

            store=store+dpstrain
            !---------------------------------------------------------------------
        
            !Determine backward stress -------------------------------------------
            DEP=0.0d0
            aux1=0.0d0
            do k1=1,3
                do k2=1,3
                    DEP(k1,1)=DEP(k1,1) + matdelas(k1,k2)*flow(k2)
                end do
            end do
            
            DEP=DEP*dp
            
            !Update stress -------------------------------------------------------
            if (counter==1)then
                do j=1,ntens
                    stress(1,1,j)=stress(1,1,j) - DEP(j,1)
                end do
            else
                
                !Difference between current stress and bakward-Euler stress -----
                do j=1,ntens
                    rr(j,1)=stress(1,1,j) - (stressold(j,1) - DEP(j,1))
                end do
                    
                dadS=0.0d0
                dadS(1,1)=2.0d0
                dadS(1,2)=-1.0d0
                dadS(2,1)=-1.0d0
                dadS(2,2)=2.0d0
                dadS(3,3)=6.0d0
                dadS=dadS/(2.0d0*SVM)
                    
                a=0.0d0
                do k1=1,3
                    do k2=1,3
                        dadS(k1,k2)=dadS(k1,k2)-flow(k1)*flow(k2)/svm
                    end do
                    a(k1,1)=flow(k1)
                end do
                    
                !Change in stress -----
                matQ=(ident + dp*matmul(matdelas,dadS))
                
                temp1 = 0.0d0    
                Qinv=0.0d0
                Qinv = matQ
                call gaussj(Qinv,ntens,temp1,errorflag)
                    
                num=0.0d0
                den=0.0d0
                aux13=matmul(transpose(a),Qinv)
                num=matmul(aux13,rr)
                
                aux13b=matmul(aux13,matdelas)
                den=matmul(aux13b,a)
                 
                !Plastic multiplier rate ------------------------------------
                dpr=0.0d0
                dpr=(YieldFunc-num(1,1))/(den(1,1) + H)
                    
                aux31 = matmul(Qinv,rr)
                aux31b = matmul(matmul(Qinv,matdelas),a)
                
                !Stress rate ------------------------------------------------
                stressrate=-1.0d0*aux31 - dpr*aux31b
                
                do k1=1,3
                    stress(1,1,k1)=stress(1,1,k1) + stressrate(k1,1)
                end do

            end if
            !---------------------------------------------------------------------
        
            !Yield surface verification ------------------------------------------
            do j=1,ntens
                stressin(j,1)=stress(1,1,j)
            end do
!        
            call vonMisesStress(stressin,ntens,SVM)
        
            call MatIsoHardening(iprops,props,6,SYield,H)
            
            !New plastic multiplier ----
            deltadp=deltadp + YieldFunc/(param + H)
            dp=deltadp
            
            !Update yield surface ----
            hard(1,1)=hard0+H*dp
            SYield=SY0+hard(1,1)
            
            !Yield Function ----
            YieldFunc=SVM-SYield 
            YieldOUT=(SVM-SYield)/SVM*100.0d0
            continue
            !---------------------------------------------------------------------
         
        end do
        
        !Inconsistent tangent matrix taken from Crisfield -------
        aux33=matmul(matmul(a,transpose(a)),matdelas)
        den=matmul(matmul(transpose(a),matdelas),a)
        
        aux33b=0.0d0
        aux33b=Ident - aux33/(den(1,1) + H)
        
        matDInc=matmul(matdelas,aux33b)
        
        !Consistent tangent matrix taken from Crisfield -------
        dadS=0.0d0
        dadS(1,1)=2.0d0
        dadS(1,2)=-1.0d0
        dadS(2,1)=-1.0d0
        dadS(2,2)=2.0d0
        dadS(3,3)=6.0d0
        dadS=dadS/(2.0d0*SVM)
                    
        a=0.0d0
        do k1=1,3
            do k2=1,3
                dadS(k1,k2)=dadS(k1,k2)-flow(k1)*flow(k2)/svm
            end do
            a(k1,1)=flow(k1)
        end do
                    
        matQ=(ident + dp*matmul(matdelas,dadS))
                    
        temp1 = 0.0d0    
        Qinv=0.0d0
        Qinv = matQ
        call gaussj(Qinv,ntens,temp1,errorflag)
        
        matR=0.0d0
        matR=matmul(Qinv,matdelas)
        
        aux33 = matmul(a,transpose(a))
        aux33b = matmul(aux33,matR)
        
        den=0.0d0
        den=matmul(matmul(transpose(a),matR),a)
        
        aux33=0.0d0
        aux33=Ident - aux33b/(den(1,1) + H)
        
        matDCons=matmul(matR,aux33)
        
        matd=matDCons
        continue
        
!        !Volumetric dilatation
!        do k1=1, ntens
!            def(k1,1) = strain(1,1,k1) + dstrain(k1,1)
!        end do
!        
!        dilat=0.0d0
!        do k1=1,ndi
!            dilat = dilat + def(k1,1) 
!        end do
!        
!        !strainzz = -1.0d0*(ve/Emod)*(stress(1,1,1)+stress(1,1,2))
!        !dilat = dilat + strainzz
!        
!        dilat=dilat/3.0d0
!        
!        devdef=0.0d0
!        do k1=1,ndi
!            devdef(k1,1)=def(k1,1) - dilat
!        end do
!        
!        do k1=ndi+1,ntens
!            devdef(k1,1) = def(k1,1)/2.0d0
!        end do
!        
!        rnorm=0.0d0
!        do k1=1,ndi
!            rdumi = 2.0d0*EG*devdef(k1,1)
!            rnorm = rnorm + rdumi*rdumi
!        end do
!        do k1=ndi+1,ntens
!            rdumi = 2.0d0*EG*devdef(k1,1)
!            rnorm = rnorm + rdumi*rdumi*2.0d0
!        end do
!        
!        rnorm = dsqrt(rnorm)
!        
!        rdumi=1.0d0/rnorm
!
!        ti=0.0d0
!        do k1=1,ndi
!            ti(k1)=1.0d0
!        end do
!        dgama=dp/sqrt(2.0d0/3.0d0)
!        nor = flow/(dsqrt(2.0d0/3.0d0))
!        rteta = 1.0d0 - (2.0d0*EG*rdumi*dgama)
!        rtetb=(1.0d0/(1.0d0+(H/(3.0d0*EG))))-(1.0d0-rteta)
!        
!        cons1=2.0d0*EG*rteta
!        cons2=2.0d0*EG*rtetb
!        conp5=cons1*0.5d00
!        
!        matd=0.0d0
!        do k1=1,ndi
!            matD(k1,k1)=cons1
!        enddo
!
!        do k1=ndi+1,ntens
!            matD(k1,k1)=conp5
!        enddo
!        
!        co113=EK-cons1/3.0d0
!        do k1=1,ntens
!            do k2=k1,ntens
!                matD(k1,k2)=matD(k1,k2)+co113*ti(k1)*ti(k2)-cons2*nor(k1)*nor(k2)
!            enddo
!        enddo
!
!        do k1=1,ntens
!            do k2=1,k1-1
!                matD(k1,k2)=matD(k2,k1)
!            enddo
!        enddo
!        
!        continue

    endif
    
    continue
    
end subroutine