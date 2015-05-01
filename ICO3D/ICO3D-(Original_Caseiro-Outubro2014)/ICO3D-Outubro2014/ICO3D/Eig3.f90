!      * * F E A P * * A Finite Element Analysis Program
!....  Copyright (c) 1984-2002: Regents of the University of California
!                               All rights reserved
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Compute eigenvalues/vectors for 3x3 symmetric matrix
!
!      Inputs:
!         v(3,3) - matrix with initial values (only upper half used)
!
!      Outputs:
!         v(3,3) - matrix of eigenvectors (by column)
!         d(3)   - eigenvalues associated with columns of v
!         rot    - number of rotations to diagonalize
!-----[--.----+----.----+----.-----------------------------------------]
!     Storage done as follows:
!
!       | v(1,1) v(1,2) v(1,3) |     |  d(1)  a(1)  a(3)  |
!       | v(2,1) v(2,2) v(2,3) |  =  |  a(1)  d(2)  a(2)  |
!       | v(3,1) v(3,2) v(3,3) |     |  a(3)  a(2)  d(3)  |
!
!        Transformations performed on d(i) and a(i) and v(i,j) become
!        eigenvectors.  Thus, original array is destroyed.
!-----[--.----+----.----+----.-----------------------------------------]     
      
      subroutine eig3(v,d,rot)

      implicit  none

      integer   rot, its, i,j,k
      real*8    g,h, aij, sm,thresh, t, c,s,tau
      real*8    v(3,3), d(3), a(3), b(3), z(3)

!     Move array into one-d arrays

      a(1) = v(1,2)
      a(2) = v(2,3)
      a(3) = v(1,3)

      do i = 1,3
        d(i) = v(i,i)
        b(i) = d(i)
        z(i) = 0.0d0
        do j = 1,3
          v(i,j) = 0.0d0
        end do ! j
        v(i,i) = 1.d0
      end do ! i

      rot = 0

      do its = 1,50

!       Set convergence test and threshold

        sm = abs(a(1)) + abs(a(2)) + abs(a(3))
        if (sm.eq.0.0d0) return

        if(its.lt.4) then
          thresh = 0.011d0*sm
        else
          thresh = 0.0d0
        end if

!       Perform sweeps for rotations

        do i = 1,3
          j = mod(i,3) + 1
          k = mod(j,3) + 1

          aij  = a(i)
          g    = 100.d0*abs(aij)
          if(abs(d(i))+g.ne.abs(d(i)) .or. abs(d(j))+g.ne.abs(d(j))) then

            if(abs(aij).gt.thresh) then
              a(i) = 0.0d0
              h    = d(j) - d(i)
              if(abs(h)+g.eq.abs(h)) then
                t = aij/h
              else
                t = sign(2.d0,h/aij)/(abs(h/aij)+sqrt(4.d0+(h/aij)**2))
              endif
!
!ccccccccccccccccccccccccccccccccccccccccc
              if(t.eq.1.d10) write(*,*)t
!ccccccccccccccccccccccccccccccccccccccccc	      
!
!             Set rotation parameters

              c    = 1.d0/sqrt(1.d0+t*t)
              s    = t*c
              tau  = s/(1.d0+c)

!             Rotate diagonal terms

              h    = t*aij
              z(i) = z(i) - h
              z(j) = z(j) + h
              d(i) = d(i) - h
              d(j) = d(j) + h

!             Rotate off-diagonal terms

              h    = a(j)
              g    = a(k)
              a(j) = h + s*(g - h*tau)
              a(k) = g - s*(h + g*tau)

!             Rotate eigenvectors

              do k = 1,3
                g      = v(k,i)
                h      = v(k,j)
                v(k,i) = g - s*(h + g*tau)
                v(k,j) = h + s*(g - h*tau)
              end do ! k

              rot = rot + 1

            end if
          else
            a(i) = 0.0d0
          end if
        end do ! i

!       Update diagonal terms

        do i = 1,3
          b(i) = b(i) + z(i)
          d(i) = b(i)
          z(i) = 0.0d0
        end do ! i

      end do ! its

      end
      
      
      !*** DETERMINA OS VALORES PRÓPRIOS E
!*** VECTORES PRÓPRIOS DE UMA
!*** MATRIZ SIMÉTRICA 2x2 OU 3x3
!*** NDI{I}:DIMENSÃO DO PROBLEMA
!*** A(NDI,NDI){O}:MATRIZ SIMÉTRICA
!*** VLP(NDI){O}:VALORES PRÓPRIOS DECRESCENTES
!*** VCP(NDI,NDI){O}:VECTORES PRÓPRIOS POR COLUNA
!*** IER{O}:INDICADOR DE ERRO
!*** [0]:NÃO OCORREU QUALQUER ERRO
!*** [1]:DIMENSÃO DO PROBLEMA É ILEGAL
!*** [2]:ERRO NA DECOMPOSIÇÃO ESPECTRAL
!
	SUBROUTINE VECP23(NDI,A,VLP,VCP,IER)
!	USE NUMERICAL_LIBRARIES
	IMPLICIT REAL(8) (A-H,O-Z)
	REAL(8),DIMENSION(NDI,NDI)::A,VCP,B
	REAL(8),DIMENSION(NDI)::VLP
	IER=0
!*** VERIFICA AS DIMENSÕES
	IF(NDI.NE.2.AND.NDI.NE.3)THEN
	IER=1
	RETURN
	ENDIF
!*** MÉTODO DE JACOBI
	!CALL DEVCSF (NDI, A, NDI, VLP,VCP,NDI)
	CALL JACOB(A,VLP,VCP,NDI,IER) 
	IF(IER.NE.0)THEN
	IER=2
	ENDIF
	END SUBROUTINE VECP23
!
!*** DECOMPOSIÇÃO ESPECTRAL DE UMA
!*** MATRIZ SIMÉTRICA
!*** A(N,N)...........MATRIZ A DECOMPOR 
!*** D(N).............VALORES PRÓPRIOS  
!*** V(N,N)...........VECTORES PRÓPRIOS 
!*** IERR.............INDICADOR DE ERRO 
!
	SUBROUTINE JACOB(A,D,V,N,IERR) 
	IMPLICIT REAL(8) (A-H,O-Z)
	REAL(8),DIMENSION(N,N)::A,V
	REAL(8),DIMENSION(N)::D
	REAL(8),DIMENSION(N)::B,Z
	REAL(8),PARAMETER::RZE=0.0D00,RUM=1.0D00,RCE=100.0D00
	INTEGER,PARAMETER::MITER=900
	IERR=0
	DO IP=1,N
	DO IQ=1,N
	V(IP,IQ)=RZE
	ENDDO
	V(IP,IP)=RUM
	ENDDO
	DO IP=1,N
	B(IP)=A(IP,IP)
	D(IP)=B(IP)
	Z(IP)=RZE
	ENDDO
	DO I=1,MITER
	SM=RZE
	DO IP=1,N-1
	DO IQ=IP+1,N
	SM=SM+ABS(A(IP,IQ))
	ENDDO
	ENDDO
	IF(SM.EQ.RZE)RETURN
	IF(I.LT.4)THEN
	TRESH=0.2D00*SM/N**2
	ELSE
	TRESH=RZE
	ENDIF
	DO IP=1,N-1
	DO IQ=IP+1,N
	G=RCE*ABS(A(IP,IQ))
	IF((I.GT.4).AND.(ABS(D(IP))+ G.EQ.ABS(D(IP))).AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
	A(IP,IQ)=RZE
	ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN
	H=D(IQ)-D(IP)
	IF(ABS(H)+G.EQ.ABS(H))THEN
	T=A(IP,IQ)/H
	ELSE
	THETA=0.5D00*H/A(IP,IQ)
	T=1.0D00/(ABS(THETA)+SQRT(RUM+THETA**2))
	IF(THETA.LT.RZE)T=-T
	ENDIF
	C=RUM/SQRT(RUM+T**2)
	S=T*C
	TAU=S/(RUM+C)
	H=T*A(IP,IQ)
	Z(IP)=Z(IP)-H
	Z(IQ)=Z(IQ)+H
	D(IP)=D(IP)-H
	D(IQ)=D(IQ)+H
	A(IP,IQ)=RZE
	DO J=1,IP-1
	G=A(J,IP)
	H=A(J,IQ)
	A(J,IP)=G-S*(H+G*TAU)
	A(J,IQ)=H+S*(G-H*TAU)
	ENDDO
	DO J=IP+1,IQ-1
	G=A(IP,J)
	H=A(J,IQ)
	A(IP,J)=G-S*(H+G*TAU)
	A(J,IQ)=H+S*(G-H*TAU)
	ENDDO
	DO J=IQ+1,N
	G=A(IP,J)
	H=A(IQ,J)
	A(IP,J)=G-S*(H+G*TAU)
	A(IQ,J)=H+S*(G-H*TAU)
	ENDDO
	DO J=1,N
	G=V(J,IP)
	H=V(J,IQ)
	V(J,IP)=G-S*(H+G*TAU)
	V(J,IQ)=H+S*(G-H*TAU)
	ENDDO
	ENDIF
	ENDDO
	ENDDO
	DO IP=1,N
	B(IP)=B(IP)+Z(IP)
	D(IP)=B(IP)
	Z(IP)=RZE
	ENDDO
	ENDDO
	IERR=1
	RETURN
	END SUBROUTINE JACOB