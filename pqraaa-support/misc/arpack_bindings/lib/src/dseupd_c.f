c
c Copyright (c) 2015
c Karl Meerbergen
c
c Distributed under the Boost Software License, Version 1.0.
c (See accompanying file LICENSE_1_0.txt
c

      subroutine dseupd_c (rvec_c, howmny, select_c, select, d, z,
     &                   ldz, sigma, bmat,
     &                   n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &                   ipntr, workd, workl, lworkl, info )
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat, howmny, which*2
      integer    rvec_c
      logical    select(ncv)
      integer    select_c(ncv)
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Double precision     
     &           sigma, tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(7), ipntr(11)
      Double precision
     &           d(nev), resid(n), v(ldv,ncv), z(ldz, nev), 
     &           workd(2*n), workl(lworkl)
c
c     %------------%
c     | Parameters |
c     %------------%
c
c     %----------------%
c     | Local varables |
c     %----------------%
c
      logical    rvec
      integer i
c
      rvec = rvec_c .ne. 0
c
      if (howmny.eq.'S') then
        do 1 i=1,ncv
          select(i) = select_c(i) .ne. 0
 1      continue
      end if
c
      call dseupd (rvec, howmny, select, d, z, ldz, sigma, bmat,
     &               n, which, nev, tol, resid, ncv, v, ldv, iparam, 
     &               ipntr, workd, workl, lworkl, info )
      end
