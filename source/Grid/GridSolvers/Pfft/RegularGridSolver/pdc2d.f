c************************************************************************
c
c Subroutine: pdc2d
c
c Purpose:
c   A fast direct method for solving the block tridiagonal
c   linear system
c 
c   Au = f,
c
c   where the matrix A is separable, that is, the system
c   can be expressed in the form:
c
c   (A1(x)M2 + M1(x)A2 + ch*M1(x)M2)u = f.
c
c   The notation "(x)" denotes the tensor product of the matrices,
c   that is, A(x)B = {a_ij*B}.
c   A1 and A2 are symmetric tridiagonal matrices of dimension
c   n1 and n2, respectively. M1 and M2 are diagonal matrices of the
c   same dimension. Restriction: the matrix M1 must be positive
c   definite and the matrix A must be nonsingular.
c
c   The above system can be written in a block form
c
c   C_i u_i-1 + D_i u_i + C_i+1 u_i+1 = f_i,                      (1)
c
c   where u_i and f_i are the i:th blocks of length n2 of the vectors
c   u and f, respectively. Here C_i = a1(i)*M2, i=2,...,n1, and
c   D_i = (b1(i) + ch*d1(i))*M2 + d1(i)*A2, i=1,...,n1.
c
c
c Version: 1.02
c
c Date: Apr 28 2009
c
c Parameters:
c
c   Input:
c          n1     - The dimension of the matrices A1 and M1
c          n2     - The dimension of the matrices A2 and M2
c          ldf    - The leading dimension of the two-dimensional 
c                   array f; ldf should be at least n2
c          a1     - The codiagonal of the matrix A1;
c                   the first components is in position a1(2)
c          b1     - The diagonal of the matrix A1
c          d1     - The diagonal of the matrix M1
c          a2     - The codiagonal of the matrix A2;
c                   the first components is in position a2(2)
c          b2     - The diagonal of the matrix A2
c          d2     - The diagonal of the matrix M2
c          ch     - The coefficient in the equation
c          ldw    - The length of double precision workspace;
c                   the minimum value is
c                   6*nl*min((n1 + 2*size - 1)/size, n1)
c                   + max(9*n1, 11*n2),
c                   where nl = 1 + max(int(log4(n1)), 0) and
c                   size is the number of processes
c          liw    - The length of integer workspace;
c                   the minimum value is
c                   6*n1 + (4**nl - 1)/3 + 2*nl + int(log4(2*size)) + 7,
c                   where nl and size are the same as previously
c          init   - Flag which indicates whether the purpose of
c                   the call is to initialize the data structures
c                   (init(1) = .true. or init(2) = .true. or
c                    init(3) = .true.) or to otherwise solve the problem;
c                   The three initialization stages are:
c                     init(1) = .true.: the communicators created earlier
c                       by init(2) = .true. are freed; can be used to
c                       alter the number of processes used in the solution;
c                     init(2) = .true.: create new internal communicators;
c                       this should be performed at least once and every
c                       time after (or within) a call with init(1) = .true.
c                     init(3) = .true.: initialize other internal workspace
c                       of the subroutine; this should be performed at least
c                       once and every time after (or within) a call with
c                       init(2) = .true.; can be used to change some of the
c                       parameters n1, a1, b1 or d1. Note: the values of
c                       the parameters n2, a2, b2, d2, ldf, f and ch can be
c                       changed without reinitialization
c          inicom - The MPI communicator to be used in the parallel
c                   execution
c
c   Input/output:
c          f      - On entry f contains the vector blocks of
c                   the right hand side vector stored by this process;
c                   on successful exit f contains the correponding
c                   blocks of the solution;
c                   The components are stored according to the
c                   block representation (1), that is, the first
c                   ldf components of f contain the block f_ilf etc;
c                   The solution is returned in the same order;
c                   The user should partition the vector according
c                   the values of ilf and iuf which are computed
c                   in the initialization when init(3) = .true.
c          ilf    - The number of the first vector block of f
c                   of size ldf stored by this process;
c                   Computed in the initialization if init(3) = .true.
c                   This should not be changed by user;
c          iuf    - The number of the last vector block of f
c                   of size ldf stored by this process;
c                   Computed in the initialization if init(3) = .true.
c                   This should not be changed by user
c
c   Output:
c          ierr   - Error flag indicating failure.
c                   Possible return values:
c                   0  no error
c                   1  n1 < 4*size-1 or n2 < 1 or ldf < n2,
c                      where size is the number of processes
c                   2  real workspace too short;
c                      the required amount of workspace can
c                      be found in iw(1) if liw > 4
c                   3  integer workspace too short;
c                      the required amount of workspace can
c                      be found in iw(2) if liw > 4
c                   4  failure in the LAPACK subroutine dstebz
c                      while solving eigenvalues;
c                      possibly one of the arrays a1, b1 or d1
c                      is incorrect
c                   5  failure in the LAPACK subroutine dstein
c                      while solving eigenvectors;
c                      possibly one of the arrays a1, b1 or d1
c                      is incorrect
c                   6  The number of processes associated to the
c                      communicator inicom is not of the form
c                      2**k, k >= 0.
c                   7  Error in a MPI subroutine
c
c   Workspace:
c          dw     - Double precision workspace, length at least ldw
c          iw     - Integer workspace, length at least liw
c
c
c Subroutines called:
c   dstebz and dstein from the LAPACK library.
c   dcopy, dscal, daxpy and dnrm2 from the BLAS1 library.
c   These subroutines can be obtained from the NETLIB archive.
c
c
c Language: FORTRAN
c
c Portability: FORTRAN-77 with do-enddo extension
c
c
c Algorithm:
c   The partial solution variant of the cyclic reduction method
c   (PSCR method) for linear systems with separable block
c   tridiagonal matrices.
c
c Complexity estimate: about 44*n1*n2*log4(n1+1) - 41*n1*n2 flops.
c
c
c References:
c
c   T. Rossi and J. Toivanen:
c   A Parallel Fast Direct Solver for Block Tridiagonal Systems with
c   Separable Matrices of Arbitrary Dimension,
c   SIAM Journal on Scientific Computing, Volume 20, Number 5, 1999,
c   pp. 1778-1793.
c   http://dx.doi.org/10.1137/S1064827597317016
c
c   T. Rossi, J. Toivanen:
c   DC2D and DC3D Version 1.0 User Manual,
c   Reports on Applied Mathematics and Computing, No. 1,
c   University of Jyvaskyla, 1997.
c
c
c Authors: Tuomo Rossi   (tro@mit.jyu.fi),
c          Jari Toivanen (tene@mit.jyu.fi)
c
c Address: Department of Mathematical Information Technology
c          FI-40014 University of Jyvaskyla
c          Finland
c
c Copyright: Tuomo Rossi and Jari Toivanen, 2009
c
c************************************************************************
      subroutine pdc2d(n1,n2,f,ldf,ilf,iuf,a1,b1,d1,a2,b2,d2,ch,
     &                 dw,ldw,iw,liw,inicom,init,ierr)
      integer n1, n2, ldf, ilf, iuf, ldw, liw, iw(liw), inicom, ierr
      double precision f(ldf*n1), a1(n1), b1(n1), d1(n1)
      double precision a2(n2), b2(n2), d2(n2), ch, dw(ldw)
      logical init(3)
c
      integer icomm, ieig, iwev, iv1, iv3, ivr, ig, ir, ix, itri
      integer ip4, iiwev, isplit
      integer iep, ife, ilb, ile, ilen, ilenf, iloc, ilocf, ilocl, isp
      integer iub, j, k, level, ll, locs, nl, np, ns
      integer comlvl, rank, size
      double precision c
c
      ierr = 0
c
      if (n1.lt.1.or.n2.lt.1.or.ldf.lt.n2) then
         ierr = 1
         return
      end if
      if (liw.lt.5) then
         ierr = 3
         return
      end if
c
      icomm  = 6
c
      if (init(1)) then
         comlvl = iw(5)
         call frcoms(iw(icomm),comlvl,ierr)         
         if (ierr.ne.0) return
      end if
c
      if (init(2)) then
c
c Initialize the communicators
c
         call mkcoms(iw(icomm),size,rank,comlvl,inicom,ierr)
         if (ierr.ne.0) return
         iw(3) = size
         iw(4) = rank
         iw(5) = comlvl
      else
         size   = iw(3)
         rank   = iw(4)
         comlvl = iw(5)
      end if
c
c nl = 1 + max(int(log(dble(n1))/log(4.d0)),0)
c
      nl = 1
      k = n1
 100  if (k.ge.4) then
         nl = nl + 1
         k = k/4
         goto 100
      end if
c
c Pointers to the real work space
c
      k = (n1 + 2*size - 1)/size
      k = min(k, n1)
      ieig  = 1
      iwev  = ieig + 6*k*nl
      iv1   = ieig + 6*k*nl
      iv3   = iv1  + n2
      ig    = iv3  + n2
      ir    = ig   + 3*n2
      ix    = ir   + 3*n2
      itri  = ix   + n2
      ivr   = itri + n2
      iw(1) = max(iwev+9*n1,ivr+n2) - 1
c
      if (iw(1).gt.ldw) then
         ierr = 2
         return
      end if
c
c Pointers to the integer work space
c
      ip4    = icomm  + comlvl
      iiwev  = ip4    + nl + 1
      isplit = iiwev  + 6*n1
      iw(2)  = isplit + nl + (4**nl - 1)/3 - 1
c
      if (iw(2).gt.liw) then
         ierr = 3
         return
      end if
c
      if (init(3)) then
         iw(ip4) = 1
         do k=1,nl
            iw(ip4+k) = 4*iw(ip4+k-1)
         end do
c
c Make the division into strips
c
         call inispl(n1,iw(isplit),iw(ip4))
c
         level = nl
         call getloc(locs,ilocf,ilocl,np,level,size,rank,iw(ip4))
         call getbnd(level,ilocf,ilb,iub,iw(isplit),iw(ip4))
         ilf = max(ilb,1)
         call getbnd(level,ilocl,ilb,iub,iw(isplit),iw(ip4))
         iuf = iub - 1
c
         ilenf = iuf - ilf + 1
c
c Compute the eigenvalues and eigenvectors for the problems
c
         do level=1,nl
            call getloc(locs,ilocf,ilocl,np,level,size,rank,iw(ip4))
c
            do iloc=ilocf,ilocl
               call getbnd(level,iloc,ilb,iub,iw(isplit),iw(ip4))
               ilen = iub - ilb - 1
               if (ilen.gt.0) then
                  isp = isplit + level + (iw(ip4+level) - 1)/3
     &                + 4*(iloc - 1)
                  if (level.eq.nl) isp = isplit
                  if (np.le.1) then
                     ife = 1
                     ile = ilen
                     iep = ieig + 6*(ilb - ilf + 1 + ilenf*(level - 1))
                  else
                     ife = max(ilf-ilb,1)
                     ile = iuf - ilb
                     iep = ieig + 6*ilenf*(level - 1)
                  end if
                  call eigval(ilen,ife,ile,a1(ilb+1),b1(ilb+1),
     &                        d1(ilb+1),dw(iep),dw(iwev),iw(iiwev),
     &                        rank,np,size,ierr)
                  if (ierr.ne.0) return
                  call eigvec(ilen,ife,ile,a1(ilb+1),b1(ilb+1),
     &                        d1(ilb+1),dw(iep),iw(isp),dw(iwev),
     &                        iw(iiwev),ierr)
                  if (ierr.ne.0) return
               end if
            end do
         end do
      end if
c
      if (init(1).or.init(2).or.init(3)) return
c
      ilenf = iuf - ilf + 1
c
c First recursion, bottom level
c
      level = nl
      call getloc(locs,ilocf,ilocl,np,level,size,rank,iw(ip4))
c
      do iloc=ilocf,ilocl
c
c Find the bounds for the strip
c
         call getbnd(level,iloc,ilb,iub,iw(isplit),iw(ip4))
         ilen = iub - ilb - 1
         if (ilen.gt.0) then
            k = ilb - ilf + 1
            iep = ieig + 6*(k + ilenf*(level - 1))
            ll = k*ldf + 1
c
            if (ilen.eq.1) then
c
c Problem with one grid column
c
               c = dw(iep+1)**2
               do k=0,n2-1
                  dw(ix+k) = c*f(ll+k)
               end do
               c = dw(iep) + ch
               call soltri(n2,a2,b2,c,d2,dw(ix),dw(itri))
               call upforc(n1,n2,ilb,iub,dw(iv1),dw(iv3),
     &                     a1,d2,dw(ix),dw(ix),.false.)
            else if (ilen.eq.2) then
c
c Problem with two grid columns
c
               call soldbl(n2,dw(iep),a2,b2,d2,ch,
     &                     f(ll),ldf,dw(ir),n2,dw(itri))
               call upforc(n1,n2,ilb,iub,dw(iv1),dw(iv3),
     &                     a1,d2,dw(ir),dw(ir+n2),.false.)
            else
c
c Problem with three grid columns
c
               call soltrb(n2,dw(iep),a2,b2,d2,ch,
     &                     f(ll),ldf,dw(ir),n2,dw(itri),.false.)
               call upforc(n1,n2,ilb,iub,dw(iv1),dw(iv3),
     &                     a1,d2,dw(ir),dw(ir+2*n2),.false.)
            end if
c
            call upsol1(f,ldf,ilf,iuf,dw(ivr),dw(iv1),dw(iv3),
     &                  n1,n2,ilb,iub,.false.)
         end if
      end do
c
c First recursion, levels through bottom - 1 to top + 1
c
      do level=nl-1,2,-1
         call getloc(locs,ilocf,ilocl,np,level,size,rank,iw(ip4))
c
         do iloc=ilocf,ilocl
c
c Find the bounds for the strip
c
            call getbnd(level,iloc,ilb,iub,iw(isplit),iw(ip4))
            ilen = iub - ilb - 1
            ns = min(ilen/4,3)
            if (ilen.gt.3.and.ns.gt.0) then
c
c Problem with 'ns' grid columns
c
               isp = isplit + level + (iw(ip4+level) - 1)/3
     &             + 4*(iloc - 1)
c
               call getfor(dw(ig),n2,ns,f,ldf,ilf,iw(isp+1),
     &                     np,rank,iw(icomm+level-1),ierr)
               if (ierr.ne.0) return
c
               if (np.le.1) then
                  ife = ilb + 1
                  ile = iub - 1
               else
                  ife = ilf
                  ile = iuf
                  if (ife.eq.ilb) ile = ile - 1
               end if
c
               do k=iv1,iv1+2*n2-1
                  dw(k) = 0.d0
               end do
c
c Go through eigenvalues one by one
c
               do j=ife-ilf,ile-ilf
                  iep = ieig + 6*(j + ilenf*(level - 1))
c
                  call ftrans(ns,n2,dw(ix),dw(ig),dw(iep+2))
c
                  c = dw(iep) + ch
                  call soltri(n2,a2,b2,c,d2,dw(ix),dw(itri))
c
                  if (ilb.ne.0)
     &               call daxpy(n2,dw(iep+1),dw(ix),1,dw(iv1),1)
                  if (iub.ne.n1+1)
     &               call daxpy(n2,dw(iep+5),dw(ix),1,dw(iv3),1)
               end do
c
               call upforc(n1,n2,ilb,iub,dw(iv1),dw(iv3),
     &                     a1,d2,dw(iv1),dw(iv3),.false.)
c
               if (np.le.1)
     &            call upsol1(f,ldf,ilf,iuf,dw(ivr),dw(iv1),dw(iv3),
     &                        n1,n2,ilb,iub,.true.)
            end if
         end do
c
         call upsol2(f,ldf,ilf,dw(iv1),dw(iv3),dw(ivr),dw(ix),n1,n2,np,
     &               level,size,rank,iw(icomm),iw(isplit),iw(ip4),ierr)
         if (ierr.ne.0) return
      end do
c
c Second recursion, levels through top to bottom - 1
c      
      do level=1,nl-1
         call getloc(locs,ilocf,ilocl,np,level,size,rank,iw(ip4))
c
         do iloc=ilocf,ilocl
c
c Find the bounds for the strip
c
            call getbnd(level,iloc,ilb,iub,iw(isplit),iw(ip4))
            ilen = iub - ilb - 1
            ns = min(ilen/4,3)
            if (ilen.gt.3.and.ns.gt.0) then
c
c Problem with 'ns' grid columns
c
               isp = isplit + level + (iw(ip4+level) - 1)/3
     &             + 4*(iloc - 1)
c
               call getfor(dw(ig),n2,ns,f,ldf,ilf,iw(isp+1),
     &                     np,rank,iw(icomm+level-1),ierr)
               if (ierr.ne.0) return
c
               if (np.le.1) then
                  ife = ilb + 1
                  ile = iub - 1
               else
                  ife = ilf
                  ile = iuf
                  if (ife.eq.ilb) ile = ile - 1
               end if
c
c Set the nonhomogenous boundary conditions
c
               call getbv(dw(iv1),dw(iv3),f,ldf,ilf,iuf,dw(ir),dw(ivr),
     &                    ilb,iub,n1,n2,np,rank)
               call upforc(n1,n2,ilb,iub,dw(iv1),dw(iv3),
     &                     a1,d2,dw(iv1),dw(iv3),.false.)
c
               do k=ir,ir+ns*n2-1
                  dw(k) = 0.d0
               end do
c
c Go through eigenvalues one by one
c
               do j=ife-ilf,ile-ilf
                  iep = ieig + 6*(j + ilenf*(level - 1))
c
                  call ftrans(ns,n2,dw(ix),dw(ig),dw(iep+2))
c
                  if (ilb.ne.0)
     &               call daxpy(n2,dw(iep+1),dw(iv1),1,dw(ix),1)
                  if (iub.ne.n1+1)
     &               call daxpy(n2,dw(iep+5),dw(iv3),1,dw(ix),1)
c
                  c = dw(iep) + ch
                  call soltri(n2,a2,b2,c,d2,dw(ix),dw(itri))
c     
                  do k=0,ns-1
                     call daxpy(n2,dw(iep+k+2),dw(ix),1,dw(ir+k*n2),1)
                  end do
               end do
c
               call upsol3(f,ldf,ilf,dw(iv1),dw(ivr),dw(ir),dw(ig),ns,
     &                     n2,np,rank,iw(icomm+level-1),iw(isp+1),ierr)
               if (ierr.ne.0) return
            end if
         end do
      end do
c
c Second recursion, bottom level
c
      level = nl
      call getloc(locs,ilocf,ilocl,np,level,size,rank,iw(ip4))
c
      do iloc=ilocf,ilocl
c
c Find the bounds for the strip
c
         call getbnd(level,iloc,ilb,iub,iw(isplit),iw(ip4))
         ilen = iub - ilb - 1
c
         if (ilen.gt.0) then
            call getbv(dw(iv1),dw(iv3),f,ldf,ilf,iuf,dw(ir),dw(ivr),
     &                 ilb,iub,n1,n2,np,rank)
c
            k = ilb - ilf + 1
            iep = ieig + 6*(k + ilenf*(level - 1))
            ll = k*ldf + 1
c
            if (ilen.eq.1) then
c
c Problem with one grid column
c
               call upforc(n1,n2,ilb,iub,f(ll),f(ll),
     &                     a1,d2,dw(iv1),dw(iv3),.true.)
               c = dw(iep+1)**2
               call dscal(n2,c,f(ll),1)
               c = dw(iep) + ch
               call soltri(n2,a2,b2,c,d2,f(ll),dw(itri))
            else if (ilen.eq.2) then
c
c Problem with two grid columns
c
               call upforc(n1,n2,ilb,iub,f(ll),f(ll+ldf),
     &                     a1,d2,dw(iv1),dw(iv3),.true.)
               call soldbl(n2,dw(iep),a2,b2,d2,ch,
     &                     f(ll),ldf,f(ll),ldf,dw(itri))
            else
c
c Problem with three grid columns
c
               call upforc(n1,n2,ilb,iub,f(ll),f(ll+2*ldf),
     &                     a1,d2,dw(iv1),dw(iv3),.true.)
               call soltrb(n2,dw(iep),a2,b2,d2,ch,
     &                     f(ll),ldf,f(ll),ldf,dw(itri),.true.)
            end if
         end if
      end do
c
      return
      end
c
c************************************************************************
c
c Initialization of the data structure containing
c the division into strips
c
      subroutine inispl(n,split,p4)
      integer n, split(*), p4(*)
c
      integer i, ipp, icp, iend, ilen, id, im, j, k, level
c
      split(1) = 0
      split(2) = n + 1
      ipp   = 1
      icp   = 3
      level = 1
c
 100  split(icp) = split(ipp)
      icp = icp + 1
      do i=1,p4(level)
         ipp  = ipp + 1
         iend = split(ipp)
         ilen = iend - split(ipp-1)
         k = min((ilen - 1)/4 + 1, 4)
         id = ilen/k
         im = mod(ilen,k)
         do j=1,4
            k = split(icp-1) + id + min(im,1)
            split(icp) = min(k,iend)
            im = max(im-1,0)
            icp  = icp + 1
         end do
      end do
      ipp   = ipp + 1
      level = level + 1
      if (split(ipp+1)-split(ipp).gt.4) goto 100
c
      return
      end
c
c************************************************************************
c
c Find the bounds for a given strip from the data structure
c
      subroutine getbnd(level,loc,ilb,iub,split,p4)
      integer level, loc, ilb, iub, split(*), p4(*)
c
      integer i
c
      i = level + (p4(level) - 1)/3 + loc
      ilb = split(i-1)
      iub = split(i)
c
      return
      end
c
c************************************************************************
c
c Compute the eigenvalues for the generalized eigensystem of length n
c
      subroutine eigval(n,jf,jl,a,b,d,eigen,dw,iw,rank,np,size,ierr)
      integer n, jf, jl, iw(6*n), rank, np, size, ierr
      double precision a(n), b(n), d(n), eigen(6,jl-jf+1), dw(8*n)
c
c      integer i, id, ic, ie, iu, ib, is, iv, kt, m
      integer i, id, ic, ie, iu, ib, is, iv, m
      double precision c
c
c Pointers to the workspace
c
      id = 1
      ic = id + n
      ie = ic + n
      iu = ie + n
c
      ib = 1
      is = ib + n
      iv = is + n
c
c Eliminate the mass matrix
c
      dw(id) = b(1)/d(1)
      do i=2,n
         dw(id+i-1) = b(i)/d(i)
         dw(ic+i-2) = a(i)/sqrt(d(i-1)*d(i))
      end do
c
c      if (np.gt.1) then
c         call blacs_get(-1,0,kt)
c         m = np*(rank/np)
c         do i=0,np-1
c            iw(iv+i) = m + i
c         end do
c         call blacs_gridmap(kt,iw(iv),1,1,np)
c         call pdstebz(kt,'A','E',n,c,c,i,i,0.d0,dw(id),dw(ic),m,i,
c     &                dw(ie),iw(ib),iw(is),dw(iu),5*n,iw(iv),4*n,ierr)
c         call blacs_gridexit(kt)
c      else
         call dstebz('A','E',n,c,c,i,i,0.d0,dw(id),dw(ic),m,i,dw(ie),
     &               iw(ib),iw(is),dw(iu),iw(iv),ierr)
c      end if
      if (ierr.ne.0) then
         ierr = 4
         return
      end if
c
      do i=1,jl-jf+1
         eigen(1,i) = dw(ie+i+jf-2)
      end do
c
      return
      end
c
c************************************************************************
c
c Compute the required components of eigenvectors
c for the generalized eigensystem of length n
c
      subroutine eigvec(n,jf,jl,a,b,d,eigen,isp,dw,iw,ierr)
      integer n, jf, jl, isp(*), iw(3*n), ierr
      double precision a(n), b(n), d(n), eigen(6,jl-jf+1), dw(9*n)
c
      double precision s, dnrm2
      integer i, j, k, ipos, id, ic, ie, iz, iu, ib, is, iv
c
c Pointers to the workspace
c
      id  = 1
      ic = id + n
      ie = ic + n
      iz = ie + n
      iu = iz + n
c
      ib = 1
      is = ib + n
      iv = is + n
c
c Eliminate the mass matrix
c
      dw(id) = b(1)/d(1)
      do i=2,n
         dw(id+i-1) = b(i)/d(i)
         dw(ic+i-2) = a(i)/sqrt(d(i-1)*d(i))
      end do
c
      iw(ib) = 1
      iw(is) = n
c
      do j=1,jl-jf+1
         dw(ie) = eigen(1,j)
         call dstein(n,dw(id),dw(ic),1,dw(ie),iw(ib),iw(is),dw(iz),n,
     &               dw(iu),iw(iv),i,ierr)
         if (ierr.ne.0) then
            ierr = 5
            return
         end if
c
c Normalize the eigenvector
c
         s = 1.d0/dnrm2(n,dw(iz),1)
         do i=1,n
            dw(iz+i-1) = s*dw(iz+i-1)/sqrt(d(i))
         end do
c
c Copy the required components
c
         if (n.le.3) then
            do k=1,n
               eigen(k+1,j) = dw(iz+k-1)
            end do
         else
            eigen(2,j) = dw(iz)
            eigen(6,j) = dw(iz+n-1)
            do k=1,3
               ipos = isp(k+1) - isp(1)
               if (ipos.lt.n) then
                  eigen(k+2,j) = dw(iz+ipos-1)
               else
                  eigen(k+2,j) = 0.d0
               end if
            end do
         end if
      end do
c
      return
      end
c
c************************************************************************
c
c Solve a coupled problem with 3 columns and n rows
c using separation technique
c
      subroutine soltrb(n,eigen,a,b,d,ch,f,ldf,u,ldu,w,midcol)
      integer n, ldf, ldu
      double precision eigen(6,3), a(n), b(n), d(n)
      double precision ch, f(ldf,3), u(ldu,3), w(n)
      logical midcol
c
      integer i
      double precision ev11, ev12, ev13, ev21, ev22, ev23
      double precision ev31, ev32, ev33, u1, u2, u3, c
c
      ev11 = eigen(2,1)
      ev12 = eigen(3,1)
      ev13 = eigen(4,1)
      ev21 = eigen(2,2)
      ev22 = eigen(3,2)
      ev23 = eigen(4,2)
      ev31 = eigen(2,3)
      ev32 = eigen(3,3)
      ev33 = eigen(4,3)
c
c First Fourier transform
c
      do i=1,n
         u1     = f(i,1)
         u2     = f(i,2)
         u3     = f(i,3)
         u(i,1) = ev11*u1 + ev12*u2 + ev13*u3
         u(i,2) = ev21*u1 + ev22*u2 + ev23*u3
         u(i,3) = ev31*u1 + ev32*u2 + ev33*u3
      end do
c
c Solve the tridiagonal systems
c
      do i=1,3
         c = eigen(1,i) + ch
         call soltri(n,a,b,c,d,u(1,i),w)
      end do
c
c Second Fourier transform
c
      if (midcol) then
         do i=1,n
            u1     = u(i,1)
            u2     = u(i,2)
            u3     = u(i,3)
            u(i,1) = ev11*u1 + ev21*u2 + ev31*u3
            u(i,2) = ev12*u1 + ev22*u2 + ev32*u3
            u(i,3) = ev13*u1 + ev23*u2 + ev33*u3
         end do
      else
         do i=1,n
            u1     = u(i,1)
            u2     = u(i,2)
            u3     = u(i,3)
            u(i,1) = ev11*u1 + ev21*u2 + ev31*u3
            u(i,3) = ev13*u1 + ev23*u2 + ev33*u3
         end do
      end if
c
      return
      end
c
c************************************************************************
c
c Solve a coupled problem with 2 columns and n rows
c using separation technique
c
      subroutine soldbl(n,eigen,a,b,d,ch,f,ldf,u,ldu,w)
      integer n, ldf, ldu
      double precision eigen(6,2), a(n), b(n), d(n)
      double precision ch, f(ldf,2), u(ldu,2), w(n)
c
      integer i
      double precision ev11, ev12, ev21, ev22, u1, u2, c
c
      ev11 = eigen(2,1)
      ev12 = eigen(3,1)
      ev21 = eigen(2,2)
      ev22 = eigen(3,2)
c
c First Fourier transform
c
      do i=1,n
         u1     = f(i,1)
         u2     = f(i,2)
         u(i,1) = ev11*u1 + ev12*u2
         u(i,2) = ev21*u1 + ev22*u2
      end do
c
c Solve the tridiagonal systems
c
      do i=1,2
         c = eigen(1,i) + ch
         call soltri(n,a,b,c,d,u(1,i),w)
      end do
c
c Second Fourier transform
c
      do i=1,n
         u1     = u(i,1)
         u2     = u(i,2)
         u(i,1) = ev11*u1 + ev21*u2
         u(i,2) = ev12*u1 + ev22*u2
      end do
c
      return
      end
c
c************************************************************************
c
c Solve a tridiagonal linear system
c
      subroutine soltri(n,a,b,s,c,x,w)
      integer n
      double precision a(n), b(n), s, c(n), x(n), w(n)
c
      integer i
      double precision d, ai, an, xp, wp
c
      d = 1.d0/(b(1) + s*c(1))
      xp = x(1)*d
      x(1) = xp
      if (n.eq.1) return
      an = a(2)
      wp = -an*d
      w(1) = wp
      do i=2,n-1
         ai = an
         d = -1.d0/(b(i) + s*c(i) + ai*wp)
         xp = (ai*xp - x(i))*d
         x(i) = xp
         an = a(i+1)
         wp = an*d
         w(i) = wp
      end do
      d = -1.d0/(b(n) + s*c(n) + an*wp)
      xp = (an*xp - x(n))*d
      x(n) = xp
      do i=n-1,1,-1
         xp = x(i) + w(i)*xp
         x(i) = xp
      end do
c
      return
      end
c
c************************************************************************
c
c Do Fourier transform for 1, 2 or 3 columns
c
      subroutine ftrans(ns,n,x,g,evec)
      integer ns, n
      double precision x(n), g(n,ns), evec(ns)
c
      integer k
      double precision e1, e2, e3
c
      if (ns.eq.3) then
         e1 = evec(1)
         e2 = evec(2)
         e3 = evec(3)
         do k=1,n
            x(k) = e1*g(k,1) + e2*g(k,2) + e3*g(k,3)
         end do
      else if (ns.eq.2) then
         e1 = evec(1)
         e2 = evec(2)
         do k=1,n
            x(k) = e1*g(k,1) + e2*g(k,2)
         end do
      else
         e1 = evec(1)
         do k=1,n
            x(k) = e1*g(k,1)
         end do
      end if
c
      return
      end
c
c************************************************************************
c
c Update a force vector
c
      subroutine upforc(n1,n2,ilb,iub,fl,fu,a1,d2,vl,vu,add)
      integer n1, n2, ilb, iub
      double precision fl(n2), fu(n2), a1(n1), d2(n2), vl(n2), vu(n2)
      logical add
c
      integer k
      double precision c
c
      if (ilb.ne.0) then
         c = -a1(ilb+1)
         if (add) then
            do k=1,n2
               fl(k) = fl(k) + c*d2(k)*vl(k)
            end do
         else
            do k=1,n2
               fl(k) = c*d2(k)*vl(k)
            end do
         end if
      end if
      if (iub.ne.n1+1) then
         c = -a1(iub)
         if (add) then
            do k=1,n2
               fu(k) = fu(k) + c*d2(k)*vu(k)
            end do
         else
            do k=1,n2
               fu(k) = c*d2(k)*vu(k)
            end do
         end if
      end if
c
      return
      end
c
c
c
      subroutine mkcoms(comm,size,rank,levels,inicom,ierr)
      integer comm(*), size, rank, levels, inicom, ierr
c
      include 'mpif.h'
c
      integer level, n, color, mpierr
c
      ierr = 7
c
      call MPI_COMM_SIZE(inicom,size,mpierr)
      if (mpierr.ne.MPI_SUCCESS) return
c
      if (size.eq.0) then
         ierr = 6
         return
      end if
      levels = int(log(dble(size))/log(2.d0) + 0.5d0)
      if (2**levels.ne.size) then
         ierr = 6
         return
      end if
c
      call MPI_COMM_RANK(inicom,rank,mpierr)
      if (mpierr.ne.MPI_SUCCESS) return
c
      levels = (levels + 1)/2
      comm(1) = inicom
      n = 1
c
      do level=2,levels
         n = 4*n
         color = n*rank/size
         call MPI_COMM_SPLIT(comm(level-1),color,rank,comm(level),
     &                       mpierr)
         if (mpierr.ne.MPI_SUCCESS) return
      end do
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine frcoms(comm,levels,ierr)
      integer levels, comm(levels), ierr
c
      include 'mpif.h'
c
      integer level, mpierr
c
      ierr = 7
c
      do level=2,levels
         call MPI_COMM_FREE(comm(level),mpierr)
         if (mpierr.ne.MPI_SUCCESS) return
      end do
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine getloc(locs,ilocf,ilocl,np,level,size,rank,p4)
      integer locs, ilocf, ilocl, np, level, size, rank, p4(level)
c
      integer m
c
      m = p4(level)
      locs = m/size
      ilocf = m*rank/size + 1
      ilocl = ilocf + max(locs-1,0)
      np = size/m
c
      return
      end
c
c
c
      subroutine getfor(g,n2,ns,f,ldf,ilf,split,np,rank,comm,ierr)
      integer n2, ns, ldf, ilf, split(ns), np, rank, comm, ierr
      double precision g(n2,3), f(*)
c
      include 'mpif.h'
c
      integer k, ll
      integer status(MPI_STATUS_SIZE), mpierr
c
      ierr = 7
c
      if (np.le.1) then
         do k=1,ns
            ll = (split(k) - ilf)*ldf + 1
            call dcopy(n2,f(ll),1,g(1,k),1)
         end do
      else if (np.eq.2) then
         k = mod(rank,2)
         if (k.eq.0) then
            ll = (split(1) - ilf)*ldf + 1
            call MPI_SEND(f(ll),n2,MPI_DOUBLE_PRECISION,1,0,comm,mpierr)
            if (mpierr.ne.MPI_SUCCESS) return
         else 
            call MPI_RECV(g(1,1),n2,MPI_DOUBLE_PRECISION,0,0,comm,
     &                    status,mpierr)
            if (mpierr.ne.MPI_SUCCESS) return
            ll = (split(2) - ilf)*ldf + 1
            call dcopy(n2,f(ll),1,g(1,2),1)
            ll = (split(3) - ilf)*ldf + 1
            call dcopy(n2,f(ll),1,g(1,3),1)
         end if
         call MPI_BCAST(g,3*n2,MPI_DOUBLE_PRECISION,1,comm,mpierr)
         if (mpierr.ne.MPI_SUCCESS) return
      else
         k = 4*mod(rank,np)
         if (k.eq.np.or.k.eq.3*np) then
            ll = (split(k/np) - ilf)*ldf + 1
            call MPI_SEND(f(ll),n2,MPI_DOUBLE_PRECISION,
     &                    np/2,0,comm,mpierr)
            if (mpierr.ne.MPI_SUCCESS) return
         else if (k.eq.2*np) then
            call MPI_RECV(g(1,1),n2,MPI_DOUBLE_PRECISION,np/4,0,comm,
     &                    status,mpierr)
            if (mpierr.ne.MPI_SUCCESS) return
            ll = (split(2) - ilf)*ldf + 1
            call dcopy(n2,f(ll),1,g(1,2),1)
            call MPI_RECV(g(1,3),n2,MPI_DOUBLE_PRECISION,3*np/4,0,comm,
     &                    status,mpierr)
            if (mpierr.ne.MPI_SUCCESS) return
         end if
         call MPI_BCAST(g,3*n2,MPI_DOUBLE_PRECISION,np/2,comm,mpierr)
         if (mpierr.ne.MPI_SUCCESS) return
      end if
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine upsol1(f,ldf,ilf,iuf,vr,v1,v3,n1,n2,ilb,iub,add)
      integer ldf, ilf, iuf, n1, n2, ilb, iub
      double precision f(*), vr(n2), v1(n2), v3(n2)
      logical add
c
      integer ll
c
      if (ilb.ne.0) then
         ll = (ilb - ilf)*ldf + 1
         call daxpy(n2,1.d0,v1,1,f(ll),1)
      end if
      if (iub.ne.iuf+1) then
         ll = (iub - ilf)*ldf + 1
         call daxpy(n2,1.d0,v3,1,f(ll),1)
      else if (iub.ne.n1+1) then
         if (add) then
            call daxpy(n2,1.d0,v3,1,vr,1)
         else
            call dcopy(n2,v3,1,vr,1)
         end if
      end if
c
      return
      end
c
c
c
      subroutine upsol2(f,ldf,ilf,v1,v3,vr,x,n1,n2,np,level,
     &                  size,rank,comm,split,p4,ierr)
      integer ldf, ilf, n1, n2, np, level, size, rank, comm(level)
      integer split(*), p4(level), ierr
      double precision f(*), v1(n2), v3(n2), vr(n2), x(n2)
c
      include 'mpif.h'
c
      integer ilb, ilocf, iub, k, ll, m
      integer dst, src, status(MPI_STATUS_SIZE), mpierr
c
      ierr = 7
c
      m = p4(level)
c
      k = 2*size/m
      if (k.eq.1.or.k.eq.2) then
         if (rank.ne.0) then
            do k=1,n2
               v1(k) = 0.d0
            end do
         end if
         if (rank.ne.size-1)
     &      call dcopy(n2,vr,1,v3,1)
      else if (np.gt.1) then
         call dcopy(2*n2,v1,1,x,1)
         call MPI_REDUCE(x,v1,2*n2,MPI_DOUBLE_PRECISION,
     &                   MPI_SUM,0,comm(level),mpierr)
         if (mpierr.ne.MPI_SUCCESS) return
      end if
c
      k = 2*size/m
      if (k.ge.1) then
         if (mod(2*rank,k).eq.0) then
            k = max(np,1)
            dst = rank + k
            if (dst.ge.size) dst = MPI_PROC_NULL
            src = rank - k
            if (src.lt.0) src = MPI_PROC_NULL
            call MPI_SENDRECV(v3,n2,MPI_DOUBLE_PRECISION,dst,0,
     &                        x,n2,MPI_DOUBLE_PRECISION,src,0,
     &                        comm(1),status,mpierr)
            if (mpierr.ne.MPI_SUCCESS) return
            if (src.ne.MPI_PROC_NULL) then
               ilocf = m*rank/size + 1
               call getbnd(level,ilocf,ilb,iub,split,p4)
               ll = (ilb - ilf)*ldf
               do k=1,n2
                  f(ll+k) = f(ll+k) + x(k) + v1(k)
               end do
            end if
         end if
      end if
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine upsol3(f,ldf,ilf,v1,vr,r,g,ns,n2,np,rank,comm,
     &                  split,ierr)
      integer ldf, ilf, ns, n2, np, rank, comm, split(3), ierr
      double precision f(*), v1(n2), vr(n2), r(n2,3), g(n2,3)
c
      include 'mpif.h'
c
      integer k, ll
      integer mpierr
c
      ierr = 7
c
      if (np.le.1) then
         do k=1,ns
            ll = (split(k) - ilf)*ldf + 1
            call dcopy(n2,r(1,k),1,f(ll),1)
         end do
      else
         call dcopy(3*n2,r,1,g,1)
         call MPI_ALLREDUCE(g,r,3*n2,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                      comm,mpierr)
         if (mpierr.ne.MPI_SUCCESS) return
c
         if (np.eq.2) then
            if (mod(rank,2).eq.0) then
               ll = (split(1) - ilf)*ldf + 1
               call dcopy(n2,r(1,1),1,f(ll),1)
               call dcopy(n2,r(1,2),1,vr,1)
            else
               ll = (split(2) - ilf)*ldf + 1
               call dcopy(n2,r(1,2),1,f(ll),1)
               ll = (split(3) - ilf)*ldf + 1
               call dcopy(n2,r(1,3),1,f(ll),1)
            end if
         else
            k = 4*mod(rank,np)
            if (k.eq.np.or.k.eq.2*np.or.k.eq.3*np) then
               ll = (split(k/np) - ilf)*ldf + 1
               call dcopy(n2,r(1,k/np),1,f(ll),1)
            end if
            k = k/np
            if (k.ne.0)
     &         call dcopy(n2,r(1,k),1,v1,1)
            if (k.ne.3)
     &         call dcopy(n2,r(1,k+1),1,vr,1)
         end if
      end if
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine getbv(v1,v3,f,ldf,ilf,iuf,r,vr,ilb,iub,n1,n2,np,rank)
      integer ldf, ilf, iuf, ilb, iub, n1, n2, np, rank
      double precision v1(n2), v3(n2), f(*), r(n2,3), vr(n2)
c
      integer k, ll
c
      if (np.le.1) then
         if (ilb.ne.0) then
            ll = (ilb - ilf)*ldf + 1
            call dcopy(n2,f(ll),1,v1,1)
         end if
         if (iub.ne.iuf+1) then
            ll = (iub - ilf)*ldf + 1
            call dcopy(n2,f(ll),1,v3,1)
         else if (iub.ne.n1+1) then
            call dcopy(n2,vr,1,v3,1)
         end if
      else
         k = mod(rank,4*np)/np
         if (k.ne.0)
     &      call dcopy(n2,r(1,k),1,v1,1)
         if (iub.ne.n1+1) then
            if (k.ne.3) then
               call dcopy(n2,r(1,k+1),1,v3,1)
            else
               call dcopy(n2,vr,1,v3,1)
            end if
         end if
      end if
c
      return
      end
