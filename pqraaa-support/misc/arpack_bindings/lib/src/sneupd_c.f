c
c Copyright (c) 2015
c Karl Meerbergen
c
c Distributed under the Boost Software License, Version 1.0.
c (See accompanying file LICENSE_1_0.txt
c

c
c-----------------------------------------------------------------------
      subroutine sneupd_c (rvecc, howmny, select_c, select, dr, di, z,
     &                   ldz, sigmar, 
     &                   sigmai, workev, bmat, n, which, nev, tol, 
     &                   resid, ncv, v, ldv, iparam, ipntr, workd, 
     &                   workl, lworkl, info)
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat, howmny, which*2
      integer    rvecc
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Real
     &           sigmar, sigmai, tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(14)
      integer    select_c(ncv)
      logical    select(ncv)
      Real
     &           dr(nev+1), di(nev+1), resid(n), v(ldv,ncv), z(ldz,*), 
     &           workd(3*n), workl(lworkl), workev(3*ncv)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer i
      logical rvec
c
      rvec = rvecc .ne. 0
c
      if (howmny.eq.'S') then
        do 1 i=1,ncv
          select(i) = select_c(i) .ne. 0
 1      continue
      end if
c
      call sneupd (rvec, howmny, select, dr, di, z, ldz, sigmar, sigmai,
     &               workev, bmat,
     &               n, which, nev, tol, resid, ncv, v, ldv, iparam, 
     &               ipntr, workd, workl, lworkl, info )
      end
