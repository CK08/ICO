subroutine PTS_SlaveSegment(iter,inc,istp,iptch)
    
    use Mod_variables
    implicit none
    
    integer(4)::loc_num, i, j, k, k1, k2, mnel, ni, count,iter,inc,istp,iptch,kspan,ilixo
    real(8)::sumtot, sumxi
    real(8)::det1, det2, det, xi 
    real(8),dimension(2)::dxdxi
    real(8),dimension(p_s+1)::xii,we
    real(8),dimension(p_s+1)::N
    real(8),dimension(p_s+1)::dNdxi
    real(8),dimension(p_s+1)::R,dRdxi
    
    real(8)::length_slv,real1
    real(8),dimension(ngrv)::suma, bl, br, cntrl_pt, tst,temp1,ratio
    real(8),dimension(ngrv-1)::slv_mid
    real(8),dimension(p_s+1)::Rss,dRss,ddRss
    real(8),dimension(ncp_s)::Rs,dRsdxi,d2Rsdxi2
    real(8),dimension(ngrv,ngrv)::MatCol
    
    mnel = 0
    length_slv = 0.0d0
    
    !Compute length of the NURBS master segment
    do k1=1, ncp_s+p_s
        
        !Master Element
        if((u_knot_slave(k1) .ne. u_knot_slave(k1+1))) then
            
            mnel = mnel + 1
            
            xi = 0.0d0
            we = 0.0d0
            call gauleg(p_s+1, xii, we)
            
            
            !Integration points cycle
            do k2=1, p_s+1
                
                ni = INN_s(IEN_s(mnel,1),1)
                xi  = ((u_knot_slave(ni+1) - u_knot_slave(ni))*xii(k2)  + (u_knot_slave(ni+1) + u_knot_slave(ni)))/2.0d0
                
                det1 = (U_knot_slave(ni+1)-U_knot_slave(ni))/2.0d0
                
                call BSplineBasisAndDeriv(ncp_s,p_s,xi,u_knot_slave,N,dNdxi)
                
                loc_num = 0
                R = 0.0d0
                dRdxi = 0.0d0
                sumtot = 0.0d0
                sumxi = 0.0d0
                do j=0,p_s
                    loc_num = loc_num + 1
                    R(loc_num) = N(p_s+1-j)*Weights_slv(IEN_s(mnel,loc_num))
                    sumtot = sumtot + R(loc_num)
                    
                    dRdxi(loc_num) = dNdxi(p_s+1-j)*Weights_slv(IEN_s(mnel,loc_num))
                    sumxi = sumxi + dRdxi(loc_num)
                end do
                
                do loc_num=1,(p+1)
                    dRdxi(loc_num) = (dRdxi(loc_num)*sumtot - R(loc_num)*sumxi )/(sumtot*sumtot)
                    R(loc_num) = R(loc_num)/sumtot
                end do
                
                loc_num = 0
                
                dxdxi = 0.0d0
                do i=0,p
                    loc_num = loc_num + 1
                    dxdxi(1) = dxdxi(1) + dRdxi(loc_num)*Points_slv(IEN_s(mnel,loc_num),1)
                    dxdxi(2) = dxdxi(2) + dRdxi(loc_num)*Points_slv(IEN_s(mnel,loc_num),2)
                end do
                
                det2 = dsqrt(dxdxi(1)*dxdxi(1)+dxdxi(2)*dxdxi(2))
                
                det = det1*det2
                
                continue
                
                !Length of the slave segment---------------
                length_slv = length_slv + det*we(k2)
                
            end do !Integration points cycle
            
        end if !Master Element
        
    end do
    
    
    !Compute mid points between consecutive collocation points
    slv_mid = 0.0d0
    do i=1,ngrv-1
       slv_mid(i) = (grev(i)+grev(i+1))/2.0d0 
    end do
    
    !Compute area associated with each collocation point
    CP_area = 0.0d0
    do i=1,ngrv
        if(i==1)then
            CP_area(i) = slv_mid(i)*length_slv
        elseif(i==ngrv)then
            CP_area(i) = (1.0d0 - slv_mid(i-1))*length_slv
        else
            CP_area(i) = (slv_mid(i) - slv_mid(i-1))*length_slv
        end if
    end do
    
    !Compute Contact Stress based on the area associated with the collocation point
    PTS_Stress = 0.0d0
    do i=1,ngrv
        PTS_Stress(i) = Lagrange(i)/CP_area(i)
    end do 
    
    do i=1,ngrv
        ratio(i) = CP_area(i)
    end do
    
    continue
    
!    bl = 0.0d0
!    br = 0.0d0
!    !Second Test
!    do i=1,ngrv
!        do j=1,(i-1)
!            bl(i) = bl(i) + (u_knot_slave(i+p+1)-u_knot_slave(i))
!        end do
!        
!        do j=1,i
!            br(i) = br(i) + (u_knot_slave(i+p+1)-u_knot_slave(i))
!        end do
!        
!        CP_area(i) = (br(i)-bl(i))*length_slv
!    end do
    
    bl = 0.0d0
    do i=1,ncp_s
        bl(i) = u_knot_slave(i+p+1)-u_knot_slave(i)
        CP_area(i) = bl(i)*length_slv
    end do
    
    continue
    
    cntrl_pt = 0.0d0
    do i=1,ncp_s
        do j=1,ngrv

            call BSplineBasisAndDeriv2(ncp_s,p_s,grev(j),u_knot_slave,Rss,dRss,ddRss,kspan)
            call BSplinesToNURBS(ncp_s,p_s,grev(j),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
        
            cntrl_pt(i) = cntrl_pt(i) + Rs(i)*Lagrange(j)
        end do
    end do
    
!    MatCol = 0.0d0
!    do i=1,ngrv
!        
!        call BSplineBasisAndDeriv2(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan)
!        call BSplinesToNURBS(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
!    
!        do j=1,ncp_s
!            MatCol(i,j) = Rs(j)
!        end do
!    end do
!    
!    temp1 = 0.0d0
!    ilixo = 0
!    call gaussj(MatCol,ngrv,temp1,ilixo)
!    
!    cntrl_pt = 0.0d0
!    do i=1,ncp_s
!        do j=1,ngrv
!            cntrl_pt(i) = cntrl_pt(i) + MatCol(i,j)*Lagrange(j)
!        end do
!    end do
    
    
    
    do i=1,ngrv
        ratio(i) = ratio(i)/CP_area(i)
    end do
    
    real1 = sum(ratio(:))/(ngrv*1.0d0)
    
    real1 = 0.0d0
    real1 = sum(bl(:))
    
    PTS_Stress = 0.0d0
    do i=1,ngrv
        !PTS_Stress(i) = cntrl_pt(i)/CP_area(i)*(1.0d0*p_s+1.0d0) !/ratio(i)
        if(PTS_active(i) == 1)then
            PTS_Stress(i) = cntrl_pt(i)/CP_area(i)*real1
        else
            PTS_Stress(i) = 0.0d0
        end if
    end do 
    
    continue
      
!    tst = 0.0d0
!    do i=1,ngrv
!        call BSplineBasisAndDeriv2(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan)
!        call BSplinesToNURBS(ncp_s,p_s,grev(i),u_knot_slave,Rss,dRss,ddRss,kspan,conn_slave,Rs,dRsdxi,d2Rsdxi2)
!        
!        do j=1,ncp_s
!            tst(i) = tst(i) + Rs(j)*PTS_Stress(j)
!        end do
!        
!    end do
    
    !Compute contact forces by element area ----------------------------------------- 
    count = 0
    do k1=1, ncp_s+p_s
        if((u_knot_slave(k1) .ne. u_knot_slave(k1+1))) then
            count = count + 1
        end if
    end do 
    
    if(iter==1 .and. inc==1 .and. istp==1)then
        allocate(Elem_Area(count))
        Elem_Area = 0.0d0
        allocate(Elem_CF(count))
        Elem_CF = 0.0d0
        allocate(Elem_Str(count))
        Elem_Str = 0.0d0
        
        allocate(Elem_Str_plot(2*count))
        Elem_Str_plot = 0.0d0
        allocate(Elem_xi_plot(2*count))
        Elem_xi_plot = 0.0d0
        
        allocate(CFs(ncp_s))
        
        allocate(CollPoint_CF(ngrv))
        CollPoint_CF = 0.0d0
    end if
    
    count = 0
    Elem_Area = 0.0d0
    Elem_CF = 0.0d0
    do k1=1, ncp_s+p_s
        if(u_knot_slave(k1) .ne. u_knot_slave(k1+1)) then
            count = count + 1
            
            Elem_Area(count) = (u_knot_slave(k1+1)-u_knot_slave(k1))*length_slv
            
            do k2=1,ngrv
                if(grev(k2) >= u_knot_slave(k1) .and. grev(k2) <= u_knot_slave(k1+1))then
                    Elem_CF(count) = Elem_CF(count) + Lagrange(k2)
                    continue
                end if
            end do
        end if
    end do 
    
    !Compute contact forces by element area (aux arrays) ---------------------------- 
    count = 0
    Elem_Str = 0.0d0
    Elem_Str_plot = 0.0d0
    Elem_xi_plot = 0.0d0
    do k1=1, ncp_s+p_s
        if(u_knot_slave(k1) .ne. u_knot_slave(k1+1)) then
        
            count = count + 1
            
            Elem_Str(count) = Elem_CF(count)/Elem_Area(count)
            
            Elem_Str_plot(count*2-1) = Elem_CF(count)/Elem_Area(count)
            Elem_Str_plot(count*2  ) = Elem_CF(count)/Elem_Area(count)
            
            Elem_xi_plot(count*2-1) = u_knot_slave(k1)
            Elem_xi_plot(count*2  ) = u_knot_slave(k1+1)
        end if
    end do
    
    count = 0
    CollPoint_CF = 0.0d0
    do k1=1, ncp_s+p_s
        if(u_knot_slave(k1) .ne. u_knot_slave(k1+1)) then
            count = count + 1
            do k2=1,ngrv
                if(grev(k2) >= u_knot_slave(k1) .and. grev(k2) <= u_knot_slave(k1+1))then
                    CollPoint_CF(k2) = CollPoint_CF(k2) + Elem_Str(count)
                end if
            end do
        end if
    end do 
    
    
    continue
    
end subroutine