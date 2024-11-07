c-----------------------------------------------------------------------
c     some routines extracted/ modified from SPARSKIT2 + one from blas
c
c     The original subroutines have been modified to support 
c     the complex case.
c-----------------------------------------------------------------------
      subroutine zreadmtc (nmax,nzmax,job,fname,a,ja,ia,rhs,nrhs,
     *     guesol,nrow,ncol,nnz,title,key,type,ierr)
c-----------------------------------------------------------------------
c     this  subroutine reads  a boeing/harwell matrix,  given the
c     corresponding file. handles right hand sides in full format
c     only (no sparse right hand sides). Also the matrix must  be
c     in assembled forms.
c     It differs from readmt, in that  the name of the file needs
c     to be passed, and then the file is opened and closed within
c     this routine. 
c     Author: Youcef Saad - Date: Oct 31, 1989
c     updated Jul 20, 1998 by Irene Moulitsas
c     
c-----------------------------------------------------------------------
c on entry:
c---------
c nmax   =  max column dimension  allowed for matrix. The array ia should
c           be of length at least ncol+1 (see below) if job.gt.0
c nzmax  = max number of nonzeros elements allowed. the arrays a,
c          and ja should be of length equal to nnz (see below) if these
c          arrays are to be read (see job).
c
c job    = integer to indicate what is to be read. (note: job is an
c          input and output parameter, it can be modified on return)
c          job = 0    read the values of ncol, nrow, nnz, title, key,
c                     type and return. matrix is not read and arrays
c                     a, ja, ia, rhs are not touched.
c          job = 1    read srtucture only, i.e., the arrays ja and ia.
c          job = 2    read matrix including values, i.e., a, ja, ia
c          job = 3    read matrix and right hand sides: a,ja,ia,rhs.
c                     rhs may contain initial guesses and exact
c                     solutions appended to the actual right hand sides.
c                     this will be indicated by the output parameter
c                     guesol [see below].
c
c fname = name of the file where to read the matrix from.
c
c nrhs   = integer. nrhs is an input as well as ouput parameter.
c          at input nrhs contains the total length of the array rhs.
c          See also ierr and nrhs in output parameters.
c
c on return:
c----------
c job    = on return job may be modified to the highest job it could
c          do: if job=2 on entry but no matrix values are available it
c          is reset to job=1 on return. Similarly of job=3 but no rhs
c          is provided then it is rest to job=2 or job=1 depending on
c          whether or not matrix values are provided.
c          Note that no error message is triggered (i.e. ierr = 0
c          on return in these cases. It is therefore important to
c          compare the values of job on entry and return ).
c
c a      = the a matrix in the a, ia, ja (column) storage format
c ja     = column number of element a(i,j) in array a.
c ia     = pointer  array. ia(i) points to the beginning of column i.
c
c rhs    = real array of size nrow + 1 if available (see job)
c
c nrhs   = integer containing the number of right-hand sides found
c          each right hand side may be accompanied with an intial guess
c          and also the exact solution.
c
c guesol = a 2-character string indicating whether an initial guess
c          (1-st character) and / or the exact solution (2-nd
c          character) is provided with the right hand side.
c          if the first character of guesol is 'G' it means that an
c          an intial guess is provided for each right-hand side.
c          These are appended to the right hand-sides in the array rhs.
c          if the second character of guesol is 'X' it means that an
c          exact solution is provided for each right-hand side.
c          These are  appended to the right hand-sides
c          and the initial guesses (if any) in the array rhs.
c
c nrow   = number of rows in matrix
c ncol   = number of columns in matrix
c nnz    = number of nonzero elements in A. This info is returned
c          even if there is not enough space in a, ja, ia, in order
c          to determine the minimum storage needed.
c
c title  = character*72 = title of matrix test ( character a*72).
c key    = character*8  = key of matrix
c type   = charatcer*3  = type of matrix.
c          for meaning of title, key and type refer to documentation
c          Harwell/Boeing matrices.
c
c ierr   = integer used for error messages
c         * ierr  =  0 means that  the matrix has been read normally.
c         * ierr  =  1 means that  the array matrix could not be read
c         because ncol+1 .gt. nmax
c         * ierr  =  2 means that  the array matrix could not be read
c         because nnz .gt. nzmax
c         * ierr  =  3 means that  the array matrix could not be read
c         because both (ncol+1 .gt. nmax) and  (nnz .gt. nzmax )
c         * ierr  =  4 means that  the right hand side (s) initial
c         guesse (s) and exact solution (s)   could  not be
c         read because they are stored in sparse format (not handled
c         by this routine ...)
c         * ierr  =  5 means that the right-hand-sides, initial guesses
c         and exact solutions could not be read because the length of
c         rhs as specified by the input value of nrhs is not
c         sufficient to store them. The rest of the matrix may have
c         been read normally.
c
c Notes:
c-------
c 1) This routine can be interfaced with the C language, since only
c    the name of the  file needs to be passed and no iounti number.
c
c 2) Refer to the  documentation on  the Harwell-Boeing formats for
c    details on the format assumed by readmt.
c    We summarize the format here for convenience.
c
c    a) all  lines in  inout are  assumed to be 80  character long.
c    b) the file  consists of a header followed by the block of the
c       column start  pointers followed  by  the  block  of the row
c       indices,  followed  by  the  block  of the  real values and
c       finally the  numerical  values of  the right-hand-side if a
c       right hand side is supplied.
c    c) the file starts by a header which contains four lines if no
c       right hand side is supplied and five lines otherwise.
c       * first  line  contains  the  title  (72  characters  long)
c         followed  by  the  8-character  identifier (name  of  the
c         matrix, called key) [ A72,A8 ]
c       * second line  contains the number of lines for each of the
c         following data blocks (4 of them) and the total number of
c         lines excluding the header.  [5i4]
c       * the   third  line  contains  a  three   character  string
c         identifying the type of  matrices as they  are referenced
c         in  the Harwell  Boeing documentation [e.g., rua, rsa,..]
c         and the number of rows, columns, nonzero entries.
c         [A3,11X,4I14]
c       * The fourth  line contains the variable fortran format for
c         the following data blocks. [2A16,2A20]
c       * The fifth  line is  only present if  right-hand-sides are
c         supplied. It  consists  of  three  one  character-strings
c         containing the  storage  format for the  right-hand-sides
c         ('F'= full,'M'=sparse=same as matrix), an  initial  guess
c         indicator  ('G' for yes),  an  exact  solution  indicator
c         ('X' for yes), followed by the number of right-hand-sides
c         and then the number of row indices.  [A3,11X,2I14]
c     d) The three  following blocks follow the header as described
c        above.
c     e) In case the right hand-side are in sparse formats then the
c        fourth  block  uses the  same  storage  format as  for the
c        matrix to  describe  the NRHS right  hand  sides provided,
c        with a column being replaced by a right hand side.
c-----------------------------------------------------------------------
      character title*72, key*8, type*3, ptrfmt*16, indfmt*16,
     &     valfmt*20, rhsfmt*20, rhstyp*3, guesol*2
      integer totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow, ncol,
     &     nnz, neltvl, nrhs, nmax, nzmax, nrwindx
      integer ia (nmax+1), ja (nzmax)
c      real*8 a(nzmax), rhs(*)

      double complex a(nzmax), rhs(*)

      character fname*100
c-----------------------------------------------------------------------
      ierr = 0
      lenrhs = nrhs
c
      iounit=15
      open(iounit,file = fname)
      read (iounit,10) title, key, totcrd, ptrcrd, indcrd, valcrd,
     &     rhscrd, type, nrow, ncol, nnz, neltvl, ptrfmt, indfmt,
     &     valfmt, rhsfmt
 10   format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
c
      if (rhscrd .gt. 0) read (iounit,11) rhstyp, nrhs, nrwindx
 11   format (a3,11x,i14,i14)
c
c     anything else to read ?
c
      if (job .le. 0) goto 12
c     ---- check whether matrix is readable ------
      n = ncol
      if (ncol .gt. nmax) ierr = 1
      if (nnz .gt. nzmax) ierr = ierr + 2
      if (ierr .ne. 0) goto 12
c     ---- read pointer and row numbers ----------
      read (iounit,ptrfmt) (ia (i), i = 1, n+1)
      read (iounit,indfmt) (ja (i), i = 1, nnz)
c     --- reading values of matrix if required....
      if (job .le. 1)  goto 12
c     --- and if available -----------------------
      if (valcrd .le. 0) then
         job = 1
         goto 12
      endif
      read (iounit,valfmt) (a(i), i = 1, nnz)
c     --- reading rhs if required ----------------
      if (job .le. 2)  goto 12
c     --- and if available -----------------------
      if ( rhscrd .le. 0) then
         job = 2
         nrhs = 0
         goto 12
      endif
c
c     --- read right-hand-side.--------------------
c
      if (rhstyp(1:1) .eq. 'M') then
         ierr = 4
         goto 12
      endif
c
      guesol = rhstyp(2:3)
c
      nvec = 1
      if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g') nvec=nvec+1
      if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x') nvec=nvec+1
c
      len = nrhs*nrow
c
      if (len*nvec .gt. lenrhs) then
        ierr = 5
        goto 12
      endif
c
c     read right-hand-sides
c
      next = 1
      iend = len
      read(iounit,rhsfmt) (rhs(i), i = next, iend)
c
c     read initial guesses if available
c
      if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g') then
        next = next+len
        iend = iend+ len
        read(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c
c     read exact solutions if available
c
      if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x') then
        next = next+len
        iend = iend+ len
        read(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c
 12   close(iounit)
      return
c---------end-of-readmt_c-----------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine  zcsrcsc (n,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n+1),ja(*),jao(*)
c      real*8  a(*),ao(*)
      double complex a(*), ao(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row     to      Compressed Sparse Column
c
c (transposition operation)   Not in place.
c-----------------------------------------------------------------------
c -- not in place --
c this subroutine transposes a matrix stored in a, ja, ia format.
c ---------------
c on entry:
c----------
c n     = dimension of A.
c job   = integer to indicate whether to fill the values (job.eq.1) of the
c         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
c
c ipos  = starting position in ao, jao of the transposed matrix.
c         the iao array takes this into account (thus iao(1) is set to ipos.)
c         Note: this may be useful if one needs to append the data structure
c         of the transpose to that of A. In this case use for example
c                call csrcsc (n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2))
c         for any other normal usage, enter ipos=1.
c a     = real array of length nnz (nnz=number of nonzero elements in input
c         matrix) containing the nonzero elements.
c ja    = integer array of length nnz containing the column positions
c         of the corresponding elements in a.
c ia    = integer of size n+1. ia(k) contains the position in a, ja of
c         the beginning of the k-th row.
c
c on return:
c ----------
c output arguments:
c ao    = real array of size nzz containing the "a" part of the transpose
c jao   = integer array of size nnz containing the column indices.
c iao   = integer array of size n+1 containing the "ia" index array of
c         the transpose.
c
c-----------------------------------------------------------------------
      call zcsrcsc2 (n,n,job,ipos,a,ja,ia,ao,jao,iao)
      end
c-----------------------------------------------------------------------
      subroutine zcsrcsc2 (n,n2,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n2+1),ja(*),jao(*)
c      real*8  a(*),ao(*)
      double complex a(*), ao(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row     to      Compressed Sparse Column
c
c (transposition operation)   Not in place.
c-----------------------------------------------------------------------
c Rectangular version.  n is number of rows of CSR matrix,
c                       n2 (input) is number of columns of CSC matrix.
c-----------------------------------------------------------------------
c -- not in place --
c this subroutine transposes a matrix stored in a, ja, ia format.
c ---------------
c on entry:
c----------
c n     = number of rows of CSR matrix.
c n2    = number of columns of CSC matrix.
c job   = integer to indicate whether to fill the values (job.eq.1) of the
c         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
c
c ipos  = starting position in ao, jao of the transposed matrix.
c         the iao array takes this into account (thus iao(1) is set to ipos.)
c         Note: this may be useful if one needs to append the data structure
c         of the transpose to that of A. In this case use for example
c                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2))
c         for any other normal usage, enter ipos=1.
c a     = real array of length nnz (nnz=number of nonzero elements in input
c         matrix) containing the nonzero elements.
c ja    = integer array of length nnz containing the column positions
c         of the corresponding elements in a.
c ia    = integer of size n+1. ia(k) contains the position in a, ja of
c         the beginning of the k-th row.
c
c on return:
c ----------
c output arguments:
c ao    = real array of size nzz containing the "a" part of the transpose
c jao   = integer array of size nnz containing the column indices.
c iao   = integer array of size n+1 containing the "ia" index array of
c         the transpose.
c
c-----------------------------------------------------------------------
c----------------- compute lengths of rows of transp(A) ----------------
      do 1 i=1,n2+1
         iao(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i), ia(i+1)-1
            j = ja(k)+1
            iao(j) = iao(j)+1
 2       continue
 3    continue
c---------- compute pointers from lengths ------------------------------
      iao(1) = ipos
      do 4 i=1,n2
         iao(i+1) = iao(i) + iao(i+1)
 4    continue
c--------------- now do the actual copying -----------------------------
      do 6 i=1,n
         do 62 k=ia(i),ia(i+1)-1
            j = ja(k)
            next = iao(j)
            if (job .eq. 1)  ao(next) = a(k)
            jao(next) = i
            iao(j) = next+1
 62      continue
 6    continue
c-------------------------- reshift iao and leave ----------------------
      do 7 i=n2,1,-1
         iao(i+1) = iao(i)
 7    continue
      iao(1) = ipos
c--------------- end of csrcsc2 ----------------------------------------
c-----------------------------------------------------------------------
      end

        subroutine qsplit(a,ind,n,ncut)
        real*8 a(n)
        integer ind(n), n, ncut
c-----------------------------------------------------------------------
c     does a quick-sort split of a real array.
c     on input a(1:n). is a real array
c     on output a(1:n) is permuted such that its elements satisfy:
c
c     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
c     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
c
c     ind(1:n) is an integer array which permuted in the same way as a(*).
c-----------------------------------------------------------------------
        real*8 tmp, abskey
        integer itmp, first, last
c-----
        first = 1
        last = n
        if (ncut .lt. first .or. ncut .gt. last) return
c
c     outer loop -- while mid .ne. ncut do
c
 1      mid = first
        abskey = abs(a(mid))
        do 2 j=first+1, last
           if (abs(a(j)) .gt. abskey) then
              mid = mid+1
c     interchange
              tmp = a(mid)
              itmp = ind(mid)
              a(mid) = a(j)
              ind(mid) = ind(j)
              a(j)  = tmp
              ind(j) = itmp
           endif
 2      continue
c
c     interchange
c
        tmp = a(mid)
        a(mid) = a(first)
        a(first)  = tmp
c
        itmp = ind(mid)
        ind(mid) = ind(first)
        ind(first) = itmp
c
c     test for while loop
c
        if (mid .eq. ncut) return
        if (mid .gt. ncut) then
           last = mid-1
        else
           first = mid+1
        endif
        goto 1
c----------------end-of-qsplit------------------------------------------
c-----------------------------------------------------------------------
        end


c-----------------------------------------------------------------------
      function zdistdot(n,x,ix,y,iy)
      integer n, ix, iy
c      real*8 distdot, x(*), y(*), ddot
       double complex zdistdot, x(*), y(*), zdotu
c      external ddot
       external zdotu
c      distdot = ddot(n,x,ix,y,iy)
       zdistdot = zdotu(n,x,ix,y,iy)
      return
      end
c-----end-of-zdistdot
c-----------------------------------------------------------------------
      subroutine zprtmtc(nrow,ncol,a,ja,ia,rhs,guesol,title,key,type,
     1     ifmt,job,fname)
c-----------------------------------------------------------------------
c writes a matrix in Harwell-Boeing format into a file.
c assumes that the matrix is stored in COMPRESSED SPARSE COLUMN FORMAT.
c some limited functionality for right hand sides. 
c Author: Youcef Saad - Date: Sept., 1989 - updated Oct. 31, 1989 to
c cope with new format. 
c-----------------------------------------------------------------------
c on entry:
c---------
c nrow   = number of rows in matrix
c ncol	 = number of columns in matrix 
c a	 = real*8 array containing the values of the matrix stored 
c          columnwise
c ja 	 = integer array of the same length as a containing the column
c          indices of the corresponding matrix elements of array a.
c ia     = integer array of containing the pointers to the beginning of 
c	   the row in arrays a and ja.
c rhs    = real array  containing the right-hand-side (s) and optionally
c          the associated initial guesses and/or exact solutions
c          in this order. See also guesol for details. the vector rhs will
c          be used only if job .gt. 2 (see below). Only full storage for
c          the right hand sides is supported. 
c
c guesol = a 2-character string indicating whether an initial guess 
c          (1-st character) and / or the exact solution (2-nd)
c          character) is provided with the right hand side.
c	   if the first character of guesol is 'G' it means that an
c          an intial guess is provided for each right-hand sides. 
c          These are assumed to be appended to the right hand-sides in 
c          the array rhs.
c	   if the second character of guesol is 'X' it means that an
c          exact solution is provided for each right-hand side.
c          These are assumed to be appended to the right hand-sides 
c          and the initial guesses (if any) in the array rhs.
c
c title  = character*72 = title of matrix test ( character a*72 ).
c key    = character*8  = key of matrix 
c type   = charatcer*3  = type of matrix.
c
c ifmt	 = integer specifying the format chosen for the real values
c	   to be output (i.e., for a, and for rhs-guess-sol if 
c          applicable). The meaning of ifmt is as follows.
c	  * if (ifmt .lt. 100) then the D descriptor is used,
c           format Dd.m, in which the length (m) of the mantissa is 
c           precisely the integer ifmt (and d = ifmt+6)
c	  * if (ifmt .gt. 100) then prtmt will use the 
c           F- descriptor (format Fd.m) in which the length of the 
c           mantissa (m) is the integer mod(ifmt,100) and the length 
c           of the integer part is k=ifmt/100 (and d = k+m+2)
c	    Thus  ifmt= 4   means  D10.4  +.xxxxD+ee    while
c	          ifmt=104  means  F7.4   +x.xxxx
c	          ifmt=205  means  F9.5   +xx.xxxxx
c	    Note: formats for ja, and ia are internally computed.
c
c job	 = integer to indicate whether matrix values and
c	   a right-hand-side is available to be written
c          job = 1   write srtucture only, i.e., the arrays ja and ia.
c          job = 2   write matrix including values, i.e., a, ja, ia
c          job = 3   write matrix and one right hand side: a,ja,ia,rhs.
c	   job = nrhs+2 write matrix and nrhs successive right hand sides
c	   Note that there cannot be any right-hand-side if the matrix
c	   has no values. Also the initial guess and exact solutions when 
c          provided are for each right hand side. For example if nrhs=2 
c          and guesol='GX' there are 6 vectors to write.
c          
c 
c fname = name of the file where to write the matrix to.
c
c on return:
c---------- 
c the matrix a, ja, ia will be written in output unit iounit
c in the Harwell-Boeing format. None of the inputs is modofied.
c  
c Notes: 1) This code attempts to pack as many elements as possible per
c        80-character line. 
c        2) this code attempts to avoid as much as possible to put
c        blanks in the formats that are written in the 4-line header
c	 (This is done for purely esthetical reasons since blanks
c        are ignored in format descriptors.)
c        3) sparse formats for right hand sides and guesses are not
c        supported.
c-----------------------------------------------------------------------
      character title*72,key*8,type*3,ptrfmt*16,indfmt*16,valfmt*20,
     *        guesol*2, rhstyp*3
      integer totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow, ncol,
     1     nnz, nrhs, len, nperli, nrwindx
      integer ja(*), ia(*)
c      real*8 a(*),rhs(*)
      double complex a(*), rhs(*)
      character fname*100
c
      iounit = 15
      open(iounit,file = fname)
c--------------
c     compute pointer format
c--------------
      nnz    = ia(ncol+1) -1
      if (nnz .eq. 0) then
         return
      endif
      len    = int ( alog10(0.1+real(nnz+1))) + 1
      nperli = 80/len
      ptrcrd = ncol/nperli + 1
      if (len .gt. 9) then
         assign 101 to ix
      else
         assign 100 to ix
      endif
      write (ptrfmt,ix) nperli,len
 100  format(1h(,i2,1HI,i1,1h) )
 101  format(1h(,i2,1HI,i2,1h) )
c----------------------------
c compute ROW index format
c----------------------------
      len    = int ( alog10(0.1+real(nrow) )) + 1
      nperli = min0(80/len,nnz)
      indcrd = (nnz-1)/nperli+1
      write (indfmt,100) nperli,len
c---------------
c compute values and rhs format (using the same for both)
c--------------- 
      valcrd = 0
      rhscrd  = 0
c quit this part if no values provided.
      if (job .le. 1) goto 20
c     
      if (ifmt .ge. 100) then
         ihead = ifmt/100
         ifmt = ifmt-100*ihead
         len = ihead+ifmt+2
         nperli = 80/len
c     
         if (len .le. 9 ) then
            assign 102 to ix
         elseif (ifmt .le. 9) then
            assign 103 to ix
         else 
            assign 104 to ix
         endif
c     
         write(valfmt,ix) nperli,len,ifmt
 102     format(1h(,i2,1hF,i1,1h.,i1,1h) )
 103     format(1h(,i2,1hF,i2,1h.,i1,1h) )
 104     format(1h(,i2,1hF,i2,1h.,i2,1h) )
C
      else
         len = ifmt + 6
         nperli = 80/len
c     try to minimize the blanks in the format strings.
         if (nperli .le. 9) then
           if (len .le. 9 ) then
              assign 105 to ix
           elseif (ifmt .le. 9) then
              assign 106 to ix
           else 
              assign 107 to ix
           endif
        else 
           if (len .le. 9 ) then
              assign 108 to ix
           elseif (ifmt .le. 9) then
              assign 109 to ix
           else 
               assign 110 to ix
            endif
         endif
c-----------
         write(valfmt,ix) nperli,len,ifmt
 105     format(1h(,i1,1hD,i1,1h.,i1,1h) )
 106     format(1h(,i1,1hD,i2,1h.,i1,1h) )
 107     format(1h(,i1,1hD,i2,1h.,i2,1h) )
 108     format(1h(,i2,1hD,i1,1h.,i1,1h) )
 109     format(1h(,i2,1hD,i2,1h.,i1,1h) )
 110     format(1h(,i2,1hD,i2,1h.,i2,1h) )
c     
      endif
      valcrd = (nnz-1)/nperli+1
      nrhs   = job -2
      if (nrhs .ge. 1) then
         i = (nrhs*nrow-1)/nperli+1
         rhscrd = i
         if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g')
     +      rhscrd = rhscrd+i
         if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x')
     +      rhscrd = rhscrd+i
         rhstyp = 'F'//guesol
      endif 
 20   continue
c     
      totcrd = ptrcrd+indcrd+valcrd+rhscrd
c     write 4-line or five line header
      write(iounit,10) title,key,totcrd,ptrcrd,indcrd,valcrd,
     1     rhscrd,type,nrow,ncol,nnz,nrhs,ptrfmt,indfmt,valfmt,valfmt
c-----------------------------------------------------------------------
      nrwindx = 0
      if (nrhs .ge. 1) write (iounit,11) rhstyp, nrhs, nrwindx
 10   format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
 11   format(A3,11x,i14,i14)
c     
      write(iounit,ptrfmt) (ia (i), i = 1, ncol+1)
      write(iounit,indfmt) (ja (i), i = 1, nnz)
      if (job .le. 1) return
      write(iounit,valfmt) (a(i), i = 1, nnz)
      if (job .le. 2) return 
      len = nrow*nrhs 
      next = 1
      iend = len
      write(iounit,valfmt) (rhs(i), i = next, iend)
c     
c     write initial guesses if available
c     
      if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g') then
         next = next+len
         iend = iend+ len
         write(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c     
c     write exact solutions if available
c     
      if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x') then
         next = next+len
         iend = iend+ len
         write(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c     
      close(iounit)
      return
c----------end of prtmt ------------------------------------------------ 
c-----------------------------------------------------------------------
      end
      
      subroutine zrnrms   (nrow, nrm, a, ia, diag) 
      double complex a(*) 
      real*8 scal, diag(nrow)
      integer ia(nrow+1) 
c-----------------------------------------------------------------------
c gets the norms of each row of A. (choice of three norms)
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = real vector of length nrow containing the norms
c
c-----------------------------------------------------------------
      do 1 ii=1,nrow
c
c     compute the norm if each element.
c     
         scal = 0.0d0
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         if (nrm .eq. 0) then
            do 2 k=k1, k2
               scal = max(scal,abs(a(k) ) ) 
 2          continue
         elseif (nrm .eq. 1) then
            do 3 k=k1, k2
               scal = scal + abs(a(k) ) 
 3          continue
         else
            do 4 k=k1, k2
               scal = scal+abs(a(k))**2
 4          continue
         endif 
         if (nrm .eq. 2) scal = sqrt(scal) 
         diag(ii) = scal
 1    continue
      return
c-----------------------------------------------------------------------
c-------------end-of-rnrms----------------------------------------------
      end 
c----------------------------------------------------------------------- 
      subroutine zcnrms   (nrow, nrm, a, ja, ia, diag) 
      double complex a(*)
      real*8 diag(nrow) 
      integer ja(*), ia(nrow+1) 
c-----------------------------------------------------------------------
c gets the norms of each column of A. (choice of three norms)
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = real vector of length nrow containing the norms
c
c-----------------------------------------------------------------
      do 10 k=1, nrow 
         diag(k) = 0.0d0
 10   continue
      do 1 ii=1,nrow
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         do 2 k=k1, k2
            j = ja(k) 
c     update the norm of each column
            if (nrm .eq. 0) then
               diag(j) = max(diag(j),abs(a(k) ) ) 
            elseif (nrm .eq. 1) then
               diag(j) = diag(j) + abs(a(k) ) 
            else
               diag(j) = diag(j)+abs(a(k))**2
            endif 
 2       continue
 1    continue
      if (nrm .ne. 2) return
      do 3 k=1, nrow
         diag(k) = sqrt(diag(k))
 3    continue
      return
c-----------------------------------------------------------------------
c------------end-of-cnrms-----------------------------------------------
      end 
c----------------------------------------------------------------------- 
      subroutine zroscal(nrow,job,nrm,a,ja,ia,diag,b,jb,ib,ierr) 
      double complex a(*), b(*)
      real*8 diag(nrow) 
      integer nrow,job,nrm,ja(*),jb(*),ia(nrow+1),ib(nrow+1),ierr 
c-----------------------------------------------------------------------
c scales the rows of A such that their norms are one on return
c 3 choices of norms: 1-norm, 2-norm, max-norm.
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = diagonal matrix stored as a vector containing the matrix
c        by which the rows have been scaled, i.e., on return 
c        we have B = Diag*A.
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c	    
c ierr  = error message. ierr=0     : Normal return 
c                        ierr=i > 0 : Row number i is a zero row.
c Notes:
c-------
c 1)        The column dimension of A is not needed. 
c 2)        algorithm in place (B can take the place of A).
c-----------------------------------------------------------------
      call zrnrms (nrow,nrm,a,ia,diag)
      ierr = 0
      do 1 j=1, nrow
         if (diag(j) .eq. 0.0d0) then
            ierr = j 
            print *, 'Row scaling with a zero diagonal: ierr = ', ierr
            return
         else
            diag(j) = 1.0d0/diag(j)
         endif
 1    continue
      call zdiamua(nrow,job,a,ja,ia,diag,b,jb,ib)
      return
c-------end-of-roscal---------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine zcoscal(nrow,job,nrm,a,ja,ia,diag,b,jb,ib,ierr) 
c----------------------------------------------------------------------- 
      double complex a(*),b(*)
      real*8 diag(nrow) 
      integer nrow,job,ja(*),jb(*),ia(nrow+1),ib(nrow+1),ierr 
c-----------------------------------------------------------------------
c scales the columns of A such that their norms are one on return
c result matrix written on b, or overwritten on A.
c 3 choices of norms: 1-norm, 2-norm, max-norm. in place.
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = diagonal matrix stored as a vector containing the matrix
c        by which the columns have been scaled, i.e., on return 
c        we have B = A * Diag
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c
c ierr  = error message. ierr=0     : Normal return 
c                        ierr=i > 0 : Column number i is a zero row.
c Notes:
c-------
c 1)     The column dimension of A is not needed. 
c 2)     algorithm in place (B can take the place of A).
c-----------------------------------------------------------------
      call zcnrms (nrow,nrm,a,ja,ia,diag)
      ierr = 0
      do 1 j=1, nrow
         if (diag(j) .eq. 0.0) then
            ierr = j 
            print *, 'Column scaling with a zero diagonal: ierr = ', 
     &          ierr
            return
         else
            diag(j) = 1.0d0/diag(j)
         endif
 1    continue
      call zamudia (nrow,job,a,ja,ia,diag,b,jb,ib)
      return
c--------end-of-coscal-------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine zdiamua (nrow,job, a, ja, ia, diag, b, jb, ib)
      double complex a(*), b(*)
      real*8 diag(nrow), scal
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1) 
c-----------------------------------------------------------------------
c performs the matrix by matrix product B = Diag * A  (in place) 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c diag = diagonal matrix stored as a vector dig(1:n)
c
c on return:
c----------
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c	    
c Notes:
c-------
c 1)        The column dimension of A is not needed. 
c 2)        algorithm in place (B can take the place of A).
c           in this case use job=0.
c-----------------------------------------------------------------
      do 1 ii=1,nrow
c     
c     normalize each row 
c     
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         scal = diag(ii) 
         do 2 k=k1, k2
            b(k) = a(k)*scal
 2       continue
 1    continue
c     
      if (job .eq. 0) return
c     
      do 3 ii=1, nrow+1
         ib(ii) = ia(ii)
 3    continue
      do 31 k=ia(1), ia(nrow+1) -1 
         jb(k) = ja(k)
 31   continue
      return
c----------end-of-diamua------------------------------------------------
c-----------------------------------------------------------------------
      end 
c----------------------------------------------------------------------- 
      subroutine zamudia (nrow,job, a, ja, ia, diag, b, jb, ib)
      double complex a(*), b(*)
      real*8 diag(nrow) 
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1) 
c-----------------------------------------------------------------------
c performs the matrix by matrix product B = A * Diag  (in place) 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c diag = diagonal matrix stored as a vector dig(1:n)
c
c on return:
c----------
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c	    
c Notes:
c-------
c 1)        The column dimension of A is not needed. 
c 2)        algorithm in place (B can take the place of A).
c-----------------------------------------------------------------
      do 1 ii=1,nrow
c     
c     scale each element 
c     
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         do 2 k=k1, k2
            b(k) = a(k)*diag(ja(k)) 
 2       continue
 1    continue
c     
      if (job .eq. 0) return
c     
      do 3 ii=1, nrow+1
         ib(ii) = ia(ii)
 3    continue
      do 31 k=ia(1), ia(nrow+1) -1 
         jb(k) = ja(k)
 31   continue
      return
c-----------------------------------------------------------------------
c-----------end-of-amudiag----------------------------------------------
      end 
c----------------------------------------------------------------------- 
