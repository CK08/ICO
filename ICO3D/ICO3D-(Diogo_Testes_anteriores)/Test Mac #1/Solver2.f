*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File symmlq.f
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine SYMMLQ( n, b, r1, r2, v, w, x, y,
     $                   Aprod, Msolve, checkA, goodb, precon, shift,
     $                   nout , itnlim, rtol,
     $                   istop, itn, Anorm, Acond, rnorm, ynorm )

      external           Aprod, Msolve
      integer            n, nout, itnlim, istop, itn
      logical            checkA, goodb, precon
      double precision   shift, rtol, Anorm, Acond, rnorm, ynorm,
     $                   b(n), r1(n), r2(n), v(n), w(n), x(n), y(n)
*     ------------------------------------------------------------------
*
*     SYMMLQ  is designed to solve the system of linear equations
*
*                Ax = b
*
*     where A is an n by n symmetric matrix and b is a given vector.
*     The matrix A is not required to be positive definite.
*     (If A is known to be definite, the method of conjugate gradients
*     might be preferred, since it will require about the same number of
*     iterations as SYMMLQ but slightly less work per iteration.)
*
*
*     The matrix A is intended to be large and sparse.  It is accessed
*     by means of a subroutine call of the form
*
*                call Aprod ( n, x, y )
*
*     which must return the product y = Ax for any given vector x.
*
*
*     More generally, SYMMLQ is designed to solve the system
*
*                (A - shift*I) x = b
*
*     where  shift  is a specified scalar value.  If  shift  and  b
*     are suitably chosen, the computed vector x may approximate an
*     (unnormalized) eigenvector of A, as in the methods of
*     inverse iteration and/or Rayleigh-quotient iteration.
*     Again, the matrix (A - shift*I) need not be positive definite.
*     The work per iteration is very slightly less if  shift = 0.
*
*
*     A further option is that of preconditioning, which may reduce
*     the number of iterations required.  If M = C C' is a positive
*     definite matrix that is known to approximate  (A - shift*I)
*     in some sense, and if systems of the form  My = x  can be
*     solved efficiently, the parameters precon and Msolve may be
*     used (see below).  When  precon = .true., SYMMLQ will
*     implicitly solve the system of equations
*
*             P (A - shift*I) P' xbar  =  P b,
*
*     i.e.                  Abar xbar  =  bbar
*     where                         P  =  C**(-1),
*                                Abar  =  P (A - shift*I) P',
*                                bbar  =  P b,
*
*     and return the solution       x  =  P' xbar.
*     The associated residual is rbar  =  bbar - Abar xbar
*                                      =  P (b - (A - shift*I)x)
*                                      =  P r.
*
*     In the discussion below, eps refers to the machine precision.
*     eps is computed by SYMMLQ.  A typical value is eps = 2.22e-16
*     for IBM mainframes and IEEE double-precision arithmetic.
*
*     Parameters
*     ----------
*
*     n       input      The dimension of the matrix A.
*
*     b(n)    input      The rhs vector b.
*
*     r1(n)   workspace
*     r2(n)   workspace
*     v(n)    workspace
*     w(n)    workspace
*
*     x(n)    output     Returns the computed solution  x.
*
*     y(n)    workspace
*
*     Aprod   external   A subroutine defining the matrix A.
*                        For a given vector x, the statement
*
*                              call Aprod ( n, x, y )
*
*                        must return the product y = Ax
*                        without altering the vector x.
*
*     Msolve  external   An optional subroutine defining a
*                        preconditioning matrix M, which should
*                        approximate (A - shift*I) in some sense.
*                        M must be positive definite.
*                        For a given vector x, the statement
*
*                              call Msolve( n, x, y )
*
*                        must solve the linear system My = x
*                        without altering the vector x.
*
*                        In general, M should be chosen so that Abar has
*                        clustered eigenvalues.  For example,
*                        if A is positive definite, Abar would ideally
*                        be close to a multiple of I.
*                        If A or A - shift*I is indefinite, Abar might
*                        be close to a multiple of diag( I  -I ).
*
*                        NOTE.  The program calling SYMMLQ must declare
*                        Aprod and Msolve to be external.
*
*     checkA  input      If checkA = .true., an extra call of Aprod will
*                        be used to check if A is symmetric.  Also,
*                        if precon = .true., an extra call of Msolve
*                        will be used to check if M is symmetric.
*
*     goodb   input      Usually, goodb should be .false.
*                        If x is expected to contain a large multiple of
*                        b (as in Rayleigh-quotient iteration),
*                        better precision may result if goodb = .true.
*                        See Lewis (1977) below.
*                        When goodb = .true., an extra call to Msolve
*                        is required.
*
*     precon  input      If precon = .true., preconditioning will
*                        be invoked.  Otherwise, subroutine Msolve
*                        will not be referenced; in this case the
*                        actual parameter corresponding to Msolve may
*                        be the same as that corresponding to Aprod.
*
*     shift   input      Should be zero if the system Ax = b is to be
*                        solved.  Otherwise, it could be an
*                        approximation to an eigenvalue of A, such as
*                        the Rayleigh quotient b'Ab / (b'b)
*                        corresponding to the vector b.
*                        If b is sufficiently like an eigenvector
*                        corresponding to an eigenvalue near shift,
*                        then the computed x may have very large
*                        components.  When normalized, x may be
*                        closer to an eigenvector than b.
*
*     nout    input      A file number.
*                        If nout .gt. 0, a summary of the iterations
*                        will be printed on unit nout.
*
*     itnlim  input      An upper limit on the number of iterations.
*
*     rtol    input      A user-specified tolerance.  SYMMLQ terminates
*                        if it appears that norm(rbar) is smaller than
*                              rtol * norm(Abar) * norm(xbar),
*                        where rbar is the transformed residual vector,
*                              rbar = bbar - Abar xbar.
*
*                        If shift = 0 and precon = .false., SYMMLQ
*                        terminates if norm(b - A*x) is smaller than
*                              rtol * norm(A) * norm(x).
*
*     istop   output     An integer giving the reason for termination...
*
*              -1        beta2 = 0 in the Lanczos iteration; i.e. the
*                        second Lanczos vector is zero.  This means the
*                        rhs is very special.
*                        If there is no preconditioner, b is an
*                        eigenvector of A.
*                        Otherwise (if precon is true), let My = b.
*                        If shift is zero, y is a solution of the
*                        generalized eigenvalue problem Ay = lambda My,
*                        with lambda = alpha1 from the Lanczos vectors.
*
*                        In general, (A - shift*I)x = b
*                        has the solution         x = (1/alpha1) y
*                        where My = b.
*                        
*               0        b = 0, so the exact solution is x = 0.
*                        No iterations were performed.
*
*               1        Norm(rbar) appears to be less than
*                        the value  rtol * norm(Abar) * norm(xbar).
*                        The solution in  x  should be acceptable.
*
*               2        Norm(rbar) appears to be less than
*                        the value  eps * norm(Abar) * norm(xbar).
*                        This means that the residual is as small as
*                        seems reasonable on this machine.
*
*               3        Norm(Abar) * norm(xbar) exceeds norm(b)/eps,
*                        which should indicate that x has essentially
*                        converged to an eigenvector of A
*                        corresponding to the eigenvalue shift.
*
*               4        Acond (see below) has exceeded 0.1/eps, so
*                        the matrix Abar must be very ill-conditioned.
*                        x may not contain an acceptable solution.
*
*               5        The iteration limit was reached before any of
*                        the previous criteria were satisfied.
*
*               6        The matrix defined by Aprod does not appear
*                        to be symmetric.
*                        For certain vectors y = Av and r = Ay, the
*                        products y'y and r'v differ significantly.
*
*               7        The matrix defined by Msolve does not appear
*                        to be symmetric.
*                        For vectors satisfying My = v and Mr = y, the
*                        products y'y and r'v differ significantly.
*
*               8        An inner product of the form  x' M**(-1) x
*                        was not positive, so the preconditioning matrix
*                        M does not appear to be positive definite.
*
*                        If istop .ge. 5, the final x may not be an
*                        acceptable solution.
*
*     itn     output     The number of iterations performed.
*
*     Anorm   output     An estimate of the norm of the matrix operator
*                        Abar = P (A - shift*I) P',   where P = C**(-1).
*
*     Acond   output     An estimate of the condition of Abar above.
*                        This will usually be a substantial
*                        under-estimate of the true condition.
*
*     rnorm   output     An estimate of the norm of the final
*                        transformed residual vector,
*                           P (b  -  (A - shift*I) x).
*
*     ynorm   output     An estimate of the norm of xbar.
*                        This is sqrt( x'Mx ).  If precon is false,
*                        ynorm is an estimate of norm(x).
*
*
*
*     To change precision
*     -------------------
*
*     Alter the words
*            double precision,
*            daxpy, dcopy, ddot, dnrm2
*     to their single or double equivalents.
*     ------------------------------------------------------------------
*
*
*     This routine is an implementation of the algorithm described in
*     the following references:
*
*     C.C. Paige and M.A. Saunders,  Solution of Sparse Indefinite
*          Systems of Linear Equations,
*          SIAM J. Numer. Anal. 12, 4, September 1975, pp. 617-629.
*
*     J.G. Lewis,  Algorithms for Sparse Matrix Eigenvalue Problems,
*          Report STAN-CS-77-595, Computer Science Department,
*          Stanford University, Stanford, California, March 1977.
*
*     Applications of SYMMLQ and the theory of preconditioning
*     are described in the following references:
*
*     D.B. Szyld and O.B. Widlund,  Applications of Conjugate Gradient
*          Type Methods to Eigenvalue Calculations,
*          in R. Vichnevetsky and R.S. Steplman (editors),
*          Advances in Computer Methods for Partial Differential
*          Equations -- III, IMACS, 1979, 167-173.
*
*     D.B. Szyld,  A Two-level Iterative Method for Large Sparse
*          Generalized Eigenvalue Calculations,
*          Ph. D. dissertation, Department of Mathematics,
*          New York University, New York, October 1983.
*
*     P.E. Gill, W. Murray, D.B. Ponceleon and M.A. Saunders,
*          Preconditioners for indefinite systems arising in
*          optimization, SIMAX 13, 1, 292--311, January 1992.
*          (SIAM J. on Matrix Analysis and Applications)
*     ------------------------------------------------------------------
*
*
*     SYMMLQ development:
*            1972: First version.
*            1975: John Lewis recommended modifications to help with
*                  inverse iteration:
*                  1. Reorthogonalize v1 and v2.
*                  2. Regard the solution as x = x1  +  bstep * b,
*                     with x1 and bstep accumulated separately
*                     and bstep * b added at the end.
*                     (In inverse iteration, b might be close to the
*                     required x already, so x1 may be a lot smaller
*                     than the multiple of b.)
*            1978: Daniel Szyld and Olof Widlund implemented the first
*                  form of preconditioning.
*                  This required both a solve and a multiply with M.
*            1979: Implemented present method for preconditioning.
*                  This requires only a solve with M.
*            1984: Sven Hammarling noted corrections to tnorm and x1lq.
*                  SYMMLQ added to NAG Fortran Library.
*     15 Sep 1985: Final F66 version.  SYMMLQ sent to "misc" in netlib.
*     16 Feb 1989: First F77 version.
*
*     22 Feb 1989: Hans Mittelmann observed beta2 = 0 (hence failure)
*                  if Abar = const*I.  istop = -1 added for this case.
*
*     01 Mar 1989: Hans Mittelmann observed premature termination on
*                  ( 1  1  1 )     (   )                   ( 1  1    )
*                  ( 1  1    ) x = ( 1 ),  for which  T3 = ( 1  1  1 ).
*                  ( 1     1 )     (   )                   (    1  1 )
*                  T2 is exactly singular, so estimating cond(A) from
*                  the diagonals of Lbar is unsafe.  We now use
*                  L       or  Lbar         depending on whether
*                  lqnorm  or  cgnorm       is least.
*
*     03 Mar 1989: eps computed internally instead of coming in as a
*                  parameter.
*     07 Jun 1989: ncheck added as a parameter to say if A and M
*                  should be checked for symmetry.
*                  Later changed to checkA (see below).
*     20 Nov 1990: goodb added as a parameter to make Lewis's changes
*                  an option.  Usually b is NOT much like x.  Setting
*                  goodb = .false. saves a call to Msolve at the end.
*     20 Nov 1990: Residual not computed exactly at end, to save time
*                  when only one or two iterations are required
*                  (e.g. if the preconditioner is very good).
*                  Beware, if precon is true, rnorm estimates the
*                  residual of the preconditioned system, not Ax = b.
*     04 Sep 1991: Parameter list changed and reordered.
*                  integer ncheck is now logical checkA.
*     22 Jul 1992: Example from Lothar Reichel and Daniela Calvetti
*                  showed that beta2 = 0 (istop = -1) means that
*                  b is an eigenvector when M = I.
*                  More complicated if there is a preconditioner;
*                  not clear yet how to describe it.
*     20 Oct 1999: Bug.  alfa1 = 0 caused Anorm = 0, divide by zero.
*                  Need to estimate Anorm from column of Tk.
*
*     Michael A. Saunders                   na.msaunders@na-net.ornl.gov
*     Department of EESOR                  mike@sol-michael.stanford.edu
*     Stanford University
*     Stanford, CA 94305-4023                             (650) 723-1875
*     ------------------------------------------------------------------
*
*
*     Subroutines and functions
*
*     USER       Aprod, Msolve
*     BLAS       daxpy, dcopy, ddot , dnrm2
*
*
*     Intrinsics and local variables

      intrinsic          abs, max, min, mod, sqrt
      double precision   ddot, dnrm2
      double precision   alfa, b1, beta, beta1, bstep, cs,
     $                   cgnorm, dbar, delta, denom, diag,
     $                   eps, epsa, epsln, epsr, epsx,
     $                   gamma, gbar, gmax, gmin, gpert,
     $                   lqnorm, oldb, qrnorm, rhs1, rhs2,
     $                   s, sn, snprod, t, tnorm,
     $                   x1cg, x1lq, ynorm2, zbar, z
      integer            i

      double precision   zero         ,  one         ,  two
      parameter        ( zero = 0.0d+0,  one = 1.0d+0,  two = 2.0d+0 )

      character*16       enter, exit
      character*52       msg(-1:8)

      data               enter /' Enter SYMMLQ.  '/,
     $                   exit  /' Exit  SYMMLQ.  '/

      data               msg
     $ / 'beta2 = 0.  If M = I, b and x are eigenvectors of A',
     $   'beta1 = 0.  The exact solution is  x = 0',
     $   'Requested accuracy achieved, as determined by rtol',
     $   'Reasonable accuracy achieved, given eps',
     $   'x has converged to an eigenvector',
     $   'Acond has exceeded 0.1/eps',
     $   'The iteration limit was reached',
     $   'Aprod  does not define a symmetric matrix',
     $   'Msolve does not define a symmetric matrix',
     $   'Msolve does not define a pos-def preconditioner' /
*     ------------------------------------------------------------------


*     Compute eps, the machine precision.  The call to daxpy is
*     intended to fool compilers that use extra-length registers.

      eps    = one / 16.0d+0

   10 eps    = eps / two
      x(1)   = eps
      y(1)   = one
      call daxpy ( 1, one, x, 1, y, 1 )
      if (y(1) .gt. one) go to 10

      eps    = eps * two

*     Print heading and initialize.

      if (nout .gt. 0) then
         write(nout, 1000) enter, n, checkA, goodb, precon,
     $                     itnlim, rtol, shift
      end if
      istop  = 0
      itn    = 0
      Anorm  = zero
      Acond  = zero
      rnorm  = zero
      ynorm  = zero

      do 50 i = 1, n
         x(i) = zero
   50 continue

*     Set up y for the first Lanczos vector v1.
*     y is really beta1 * P * v1  where  P = C**(-1).
*     y and beta1 will be zero if b = 0.

      call dcopy ( n, b, 1, y , 1 )
      call dcopy ( n, b, 1, r1, 1 )
      if ( precon ) call Msolve( n, r1, y )
      if ( goodb  ) then
         b1  = y(1)
      else
         b1  = zero
      end if
      beta1  = ddot  ( n, r1, 1, y, 1 )

*     See if Msolve is symmetric.

      if (checkA  .and.  precon) then
         call Msolve( n, y, r2 )
         s      = ddot  ( n, y, 1, y, 1 )
         t      = ddot  ( n,r1, 1,r2, 1 )
         z      = abs( s - t )
         epsa   = (s + eps) * eps**0.33333
         if (z .gt. epsa) then
            istop = 7
            go to 900
         end if
      end if

*     Test for an indefinite preconditioner.

      if (beta1 .lt. zero) then
         istop = 8
         go to 900
      end if

*     If b = 0 exactly, stop with x = 0.

      if (beta1 .eq. zero) then
         go to 900
      end if

*     Here and later, v is really P * (the Lanczos v).

      beta1  = sqrt( beta1 )
      s      = one / beta1
      do 100 i = 1, n
         v(i)  = s * y(i)
  100 continue

*     See if Aprod  is symmetric.

      call Aprod ( n, v, y )
      if (checkA) then
         call Aprod ( n, y, r2 )
         s      = ddot  ( n, y, 1, y, 1 )
         t      = ddot  ( n, v, 1,r2, 1 )
         z      = abs( s - t )
         epsa   = (s + eps) * eps**0.33333
         if (z .gt. epsa) then
            istop = 6
            go to 900
         end if
      end if

*     Set up y for the second Lanczos vector.
*     Again, y is beta * P * v2  where  P = C**(-1).
*     y and beta will be zero or very small if b is an eigenvector.

      call daxpy ( n, (- shift), v, 1, y, 1 )
      alfa   = ddot  ( n, v, 1, y, 1 )
      call daxpy ( n, (- alfa / beta1), r1, 1, y, 1 )

*     Make sure  r2  will be orthogonal to the first  v.

      z      = ddot  ( n, v, 1, y, 1 )
      s      = ddot  ( n, v, 1, v, 1 )
      call daxpy ( n, (- z / s), v, 1, y, 1 )

      call dcopy ( n, y, 1, r2, 1 )
      if ( precon ) call Msolve( n, r2, y )
      oldb   = beta1
      beta   = ddot  ( n, r2, 1, y, 1 )

      if (beta .lt. zero) then
         istop = 8
         go to 900
      end if

*     Cause termination (later) if beta is essentially zero.

      beta   = sqrt( beta )
      if (beta .le. eps) then
         istop = -1
      end if

*     See if the local reorthogonalization achieved anything.

      denom  = sqrt( s ) * dnrm2( n, r2, 1 )  +  eps
      s      = z / denom
      t      = ddot  ( n, v, 1, r2, 1 ) / denom
      if (nout .gt. 0  .and.  goodb) then
         write(nout, 1100) beta1, alfa, s, t
      end if

*     Initialize other quantities.

      cgnorm = beta1
      gbar   = alfa
      dbar   = beta
      rhs1   = beta1
      rhs2   = zero
      bstep  = zero
      snprod = one
      tnorm  = alfa**2 + beta**2
      ynorm2 = zero
      gmax   = abs( alfa ) + eps
      gmin   = gmax

      if ( goodb ) then
         do 200 i = 1, n
            w(i)  = zero
  200    continue
      else
         call dcopy ( n, v, 1, w, 1 )
      end if

*     ------------------------------------------------------------------
*     Main iteration loop.
*     ------------------------------------------------------------------

*     Estimate various norms and test for convergence.

  300 Anorm  = sqrt( tnorm  )
      ynorm  = sqrt( ynorm2 )
      epsa   = Anorm * eps
      epsx   = Anorm * ynorm * eps
      epsr   = Anorm * ynorm * rtol
      diag   = gbar
      if (diag .eq. zero) diag = epsa

      lqnorm = sqrt( rhs1**2 + rhs2**2 )
      qrnorm = snprod * beta1
      cgnorm = qrnorm * beta / abs( diag )

*     Estimate  cond(A).
*     In this version we look at the diagonals of  L  in the
*     factorization of the tridiagonal matrix,  T = L*Q.
*     Sometimes, T(k) can be misleadingly ill-conditioned when
*     T(k+1) is not, so we must be careful not to overestimate Acond.

      if (lqnorm .le. cgnorm) then
         Acond  = gmax / gmin
      else
         denom  = min( gmin, abs( diag ) )
         Acond  = gmax / denom
      end if

*     See if any of the stopping criteria are satisfied.
*     In rare cases, istop is already -1 from above (Abar = const * I).

      if (istop .eq. 0) then
         if (itn    .ge. itnlim ) istop = 5
         if (Acond  .ge. 0.1/eps) istop = 4
         if (epsx   .ge. beta1  ) istop = 3
         if (cgnorm .le. epsx   ) istop = 2
         if (cgnorm .le. epsr   ) istop = 1
      end if
*     ==================================================================

*     See if it is time to print something.

      if (nout .le.  0)          go to 600
      if (n    .le. 40)          go to 400
      if (itn  .le. 10)          go to 400
      if (itn  .ge. itnlim - 10) go to 400
      if (mod(itn,10)  .eq.   0) go to 400
      if (cgnorm .le. 10.0*epsx) go to 400
      if (cgnorm .le. 10.0*epsr) go to 400
      if (Acond  .ge. 0.01/eps ) go to 400
      if (istop  .ne. 0)         go to 400
      go to 600

*     Print a line for this iteration.

  400 zbar   = rhs1 / diag
      z      = (snprod * zbar  +  bstep) / beta1
      x1lq   = x(1)  +  b1 * bstep / beta1
      x1cg   = x(1)  +  w(1) * zbar  +  b1 * z

      if (    itn     .eq. 0) write(nout, 1200)
      write(nout, 1300) itn, x1cg, cgnorm, bstep/beta1, Anorm, Acond
      if (mod(itn,10) .eq. 0) write(nout, 1500)
*     ==================================================================


*     Obtain the current Lanczos vector  v = (1 / beta)*y
*     and set up  y  for the next iteration.

  600 if (istop .ne. 0) go to 800
      s      = one / beta

      do 620 i = 1, n
         v(i)  = s * y(i)
  620 continue

      call Aprod ( n, v, y )
      call daxpy ( n, (- shift), v, 1, y, 1 )
      call daxpy ( n, (- beta / oldb), r1, 1, y, 1 )
      alfa   = ddot( n, v, 1, y, 1 )
      call daxpy ( n, (- alfa / beta), r2, 1, y, 1 )
      call dcopy ( n, r2, 1, r1, 1 )
      call dcopy ( n, y, 1, r2, 1 )
      if ( precon ) call Msolve( n, r2, y )
      oldb   = beta
      beta   = ddot  ( n, r2, 1, y, 1 )
      if (beta .lt. zero) then
         istop = 6
         go to 800
      end if
      beta   = sqrt( beta )
      tnorm  = tnorm  +  alfa**2  +  oldb**2  +  beta**2

*     Compute the next plane rotation for  Q.

      gamma  = sqrt( gbar**2 + oldb**2 )
      cs     = gbar / gamma
      sn     = oldb / gamma
      delta  = cs * dbar  +  sn * alfa
      gbar   = sn * dbar  -  cs * alfa
      epsln  = sn * beta
      dbar   =            -  cs * beta

*     Update  x.

      z      = rhs1 / gamma
      s      = z * cs
      t      = z * sn

      do 700 i = 1, n
         x(i)  = (w(i) * s   +   v(i) * t)  +  x(i)
         w(i)  =  w(i) * sn  -   v(i) * cs
  700 continue

*     Accumulate the step along the direction  b,
*     and go round again.

      bstep  = snprod * cs * z  +  bstep
      snprod = snprod * sn
      gmax   = max( gmax, gamma )
      gmin   = min( gmin, gamma )
      ynorm2 = z**2  +  ynorm2
      rhs1   = rhs2  -  delta * z
      rhs2   =       -  epsln * z
      itn    = itn   +  1
      go to 300

*     ------------------------------------------------------------------
*     End of main iteration loop.
*     ------------------------------------------------------------------

*     Move to the CG point if it seems better.
*     In this version of SYMMLQ, the convergence tests involve
*     only cgnorm, so we're unlikely to stop at an LQ point,
*     EXCEPT if the iteration limit interferes.

  800 if (cgnorm .le. lqnorm) then
         zbar   = rhs1 / diag
         bstep  = snprod * zbar  +  bstep
         ynorm  = sqrt( ynorm2  +  zbar**2 )
         rnorm  = cgnorm
         call daxpy ( n, zbar, w, 1, x, 1 )
      else
         rnorm  = lqnorm
      end if

      if ( goodb ) then

*        Add the step along  b.

         bstep  = bstep / beta1
         call dcopy ( n, b, 1, y, 1 )
         if ( precon ) call Msolve( n, b, y )
         call daxpy ( n, bstep, y, 1, x, 1 )
      end if

*     ==================================================================
*     Display final status.
*     ==================================================================
  900 if (nout  .gt. 0) then
         write(nout, 2000) exit, istop, itn,
     $                     exit, Anorm, Acond,
     $                     exit, rnorm, ynorm
         write(nout, 3000) exit, msg(istop)
      end if

      return

*     ------------------------------------------------------------------
 1000 format(// 1p,    a, 5x, 'Solution of symmetric   Ax = b'
     $       / ' n      =', i7, 5x, 'checkA =', l4, 12x,
     $          'goodb  =', l4, 7x, 'precon =', l4
     $       / ' itnlim =', i7, 5x, 'rtol   =', e11.2, 5x,
     $          'shift  =', e23.14)
 1100 format(/ 1p, ' beta1  =', e10.2, 3x, 'alpha1 =', e10.2
     $       / ' (v1,v2) before and after ', e14.2
     $       / ' local reorthogonalization', e14.2)
 1200 format(// 5x, 'itn', 7x, 'x1(cg)', 10x,
     $         'norm(r)', 5x, 'bstep', 7x, 'norm(A)', 3X, 'cond(A)')
 1300 format(1p, i8, e19.10, e11.2, e14.5, 2e10.2)
 1500 format(1x)
 2000 format(/ 1p, a, 6x, 'istop =', i3,   15x, 'itn   =', i8
     $       /     a, 6x, 'Anorm =', e12.4, 6x, 'Acond =', e12.4
     $       /     a, 6x, 'rnorm =', e12.4, 6x, 'ynorm =', e12.4)
 3000 format(      a, 6x, a )
*     ------------------------------------------------------------------
*     end of SYMMLQ
      end



*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File symmlqblas.f
*
*     daxpy    dcopy    ddot     dnrm2
*
*** from netlib, Thu May 16 21:00:13 EDT 1991 ***
*** Declarations of the form dx(1) changed to dx(*)
*
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine  dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function dnrm2 ( n, dx, incx)
      integer          next
      double precision   dx(*), cutlo, cuthi, hitest, sum, xmax,zero,one
      data   zero, one /0.0d0, 1.0d0/
c
c     euclidean norm of the n-vector stored in dx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson, 1978 jan 08
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
c         cuthi = minimum of  dsqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
c
      if(n .gt. 0) go to 10
         dnrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      nn = n * incx
c                                                 begin main loop
      i = 1
   20    go to next,(30, 50, 70, 110)
   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
c
c                        phase 1.  sum is zero
c
   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 i = j
      assign 110 to next
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( dabs(dx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         go to 200
c
  115 sum = sum + (dx(i)/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
   85 hitest = cuthi/float( n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
      do 95 j =i,nn,incx
      if(dabs(dx(j)) .ge. hitest) go to 100
   95    sum = sum + dx(j)**2
      dnrm2 = dsqrt( sum )
      go to 300
c
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
c
c              end of main loop.
c
c              compute square root and adjust for scaling.
c
      dnrm2 = xmax * dsqrt(sum)
  300 continue
      return
      end
