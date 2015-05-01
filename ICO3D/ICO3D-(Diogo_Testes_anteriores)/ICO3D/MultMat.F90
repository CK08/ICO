!--------------------------------------------------------
!
! Subroutine for matrix multiplication
!
!--------------------------------------------------------

subroutine multmat(amatr,bmatr,rowa,cola,colb,cmatr)

      integer rowa, cola, colb,i,j,k
      real*8 amatr(rowa,cola),cmatr(rowa,colb),bmatr(cola,colb)



      do i = 1,rowa
          do j = 1,colb
              cmatr(i,j) = 0.d0
              do k = 1,cola
              cmatr(i,j) = amatr(i,k)*bmatr(k,j) + cmatr(i,j)
              enddo
          enddo
      enddo

      end subroutine multmat