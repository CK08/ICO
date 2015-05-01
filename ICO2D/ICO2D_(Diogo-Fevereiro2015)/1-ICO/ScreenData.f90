!----------------------------------------------------------------------------------------------
!
! Subroutine write the analysis information to the screen and to a file
!
!----------------------------------------------------------------------------------------------
subroutine ScreenData(itag,inc,iter,AbsRes,absdisp,cput)

    implicit none
    integer(4)::itag,inc,iter
    real(8)::Absres,absdisp
    real(4)::cput

    if(itag==1)then
        open(unit=10,file='Analysis_Summary.txt', status = 'REPLACE')
        open(unit=11,file='GPCoords.txt', status = 'REPLACE')
        close(11)
    else
        open(unit=10,file='Analysis_Summary.txt', access = 'APPEND')
    end if
    
    if(itag==1)then
        write(*,*)'!---------------------------------------------------------!'
        write(*,*)'!                                                         !'
        write(*,*)'!      %%%%%%%%%%%%%%%%   %%%%%%%%%%%%   %%%%%%%%%%%%%%%  !'
        write(*,*)'!            %%          %%             %%           %%   !'
        write(*,*)'!           %%          %%             %%           %%    !'
        write(*,*)'!          %%          %%             %%           %%     !'
        write(*,*)'!         %%          %%             %%           %%      !'
        write(*,*)'!        %%          %%             %%           %%       !'
        write(*,*)'!       %%          %%             %%           %%        !'
        write(*,*)'! %%%%%%%%%%%%%%   %%%%%%%%%%%%   %%%%%%%%%%%%%%%         !'
        write(*,*)'!                                                         !'
        write(*,*)'!---------------------------------------------------------!'
        write(*,*)'!                                                         !'
        write(*,*)'!      Isogeometric COde (ICO) for 2D applications        !'
        write(*,*)'!                                                         !'
        write(*,*)'!                            GRIDS - University of Aveiro !' 
        write(*,*)'!---------------------------------------------------------!'
        write(*,*)''
        write(*,*)'Input file name (including file extension): '
        
        
        write(10,*)'!---------------------------------------------------------!'
        write(10,*)'!                                                         !'
        write(10,*)'!      %%%%%%%%%%%%%%%%   %%%%%%%%%%%%   %%%%%%%%%%%%%%%  !'
        write(10,*)'!            %%          %%             %%           %%   !'
        write(10,*)'!           %%          %%             %%           %%    !'
        write(10,*)'!          %%          %%             %%           %%     !'
        write(10,*)'!         %%          %%             %%           %%      !'
        write(10,*)'!        %%          %%             %%           %%       !'
        write(10,*)'!       %%          %%             %%           %%        !'
        write(10,*)'! %%%%%%%%%%%%%%   %%%%%%%%%%%%   %%%%%%%%%%%%%%%         !'
        write(10,*)'!                                                         !'
        write(10,*)'!---------------------------------------------------------!'
        write(10,*)'!                                                         !'
        write(10,*)'!      Isogeometric COde (ICO) for 2D applications        !'
        write(10,*)'!                                                         !'
        write(10,*)'!                            GRIDS - University of Aveiro !' 
        write(10,*)'!---------------------------------------------------------!'
        write(10,*)''
        write(10,*)'Input file name (including file extension): '
        
    elseif(itag==2)then
        
        write(*,*)''
        write(*,*)'Reading Input File'
        
        write(10,*)''
        write(10,*)'Reading Input File'
        
    elseif(itag==3)then
        
        write(*,*)'...'
        write(*,*)'Done'
        
        write(10,*)'...'
        write(10,*)'Done'
        
    elseif(itag==4)then
        
        write(*,*)''
        write(*,*)'Generating Data Structure'
        
        write(10,*)''
        write(10,*)'Generating Data Structure'
    
    elseif(itag==5)then
    
        write(*,*)'...'
        write(*,*)'Done'
        write(*,*)''
        write(*,*)'|-----------------------------------------------------------------------------|'
        write(*,*)'|                               Begin Analysis                                |'
        write(*,*)'|-----------------------------------------------------------------------------|'
        
        write(10,*)'...'
        write(10,*)'Done'
        write(10,*)''
        write(10,*)'|-----------------------------------------------------------------------------|'
        write(10,*)'|                               Begin Analysis                                |'
        write(10,*)'|-----------------------------------------------------------------------------|'  
    
    elseif(itag==6)then
        
        write(*,*)''
        write(*,*)'|-----------------------------------------------------------------------------|'
        write(*,FMT=60)inc
        write(*,*)'|-----------------------------------------------------------------------------|'
        write(*,*)'|  Iteration number   |       Force Residual      |   Displacement Residual   |'
        write(*,*)'|---------------------|---------------------------|---------------------------|'
        
        write(10,*)''
        write(10,*)'|-----------------------------------------------------------------------------|'
        write(10,FMT=60)inc
        write(10,*)'|-----------------------------------------------------------------------------|'
        write(10,*)'|  Iteration number   |       Force Residual      |   Displacement Residual   |'
        write(10,*)'|---------------------|---------------------------|---------------------------|'
        
        60 format(' |                               INCREMENT ', I3,'                                  |')
        
    elseif(itag==7)then
        write(*,*)'Warning - No load or displacement applied to the model'
        write(10,*)'Warning - No load or displacement applied to the model'
    
    elseif(itag==8)then
        write(*,*)'|-----------------------------------------------------------------------------|'
        write(*,*)''
        write(*,*)''
        write(*,*)'The solution has converged at iteration', iter
        write(*,FMT=65)AbsRes
        write(*,FMT=66)AbsDisp
        
        write(10,*)'|-----------------------------------------------------------------------------|'
        write(10,*)''
        write(10,*)''
        write(10,*)'The solution has converged at iteration', iter
        write(10,FMT=65)AbsRes
        write(10,FMT=66)AbsDisp
        
        65 format('  -> The absolute residual of the forces is R=        ', E)
        66 format('  -> The absolute residual of the displacement is Rd= ', E)
        write(*,*)''
    
    elseif(itag==9)then
        
        write(*,*)'|-----------------------------------------------------------------------------|'
        write(*,*)''
        write(*,*)' | Solution has not converged after iteration', iter
        write(*,*)' | Teminating analysis'
        write(*,*)'|-----------------------------------------------------------------------------|'
        
        write(10,*)'|-----------------------------------------------------------------------------|'
        write(10,*)''
        write(10,*)' | Solution has not converged after iteration', iter
        write(10,*)' | Teminating analysis'
        write(10,*)'|-----------------------------------------------------------------------------|'
    
    elseif(itag==10)then
    
        write(*,*)''
        write(*,*)'|-----------------------------------------------------------------------------|'
        write(*,*)'|                              Analysis Complete                              |'
        write(*,*)'|-----------------------------------------------------------------------------|'
        
        write(10,*)''
        write(10,*)'|-----------------------------------------------------------------------------|'
        write(10,*)'|                              Analysis Complete                              |'
        write(10,*)'|-----------------------------------------------------------------------------|'
        
    elseif(itag==11)then
        
        write(*,*)''
        write(*,*)'|-----------------------------------------------------------------------------|'
        write(*,*)'|Total CPU time  ',CPUt 
        write(*,*)'|-----------------------------------------------------------------------------|'
        
        write(10,*)''
        write(10,*)'|-----------------------------------------------------------------------------|'
        write(10,*)'|Total CPU time  ',CPUt 
        write(10,*)'|-----------------------------------------------------------------------------|' 
        
    elseif(itag==12)then
        
        write(*,FMT=78)iter, Absres,AbsDisp  
        write(10,FMT=78)iter, Absres,AbsDisp  
        78 format(' |        ',I3 , '          |', E, '  |'E , '  |')
    end if
    
    close(10)
end subroutine
