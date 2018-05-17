            PROGRAM  POLRUN
c------                  MAINLINE FOR N(H) ANALYSIS USING 'POLAN'.    Feb 1985.
c 28.1.93  Basic heads to 80 chr; Data,Peak.. 98/100; debug 132.  Del init '0'.
c 28.1.93  Read input file name;
c     Read (1):   field  (and station heading);
c     Read (2):   data,  data,  data,   . . .
c
c     Use 1 blank line to reread a station,field line (1);   2 blanks to exit.
c
c     An initial data line with FH = 9.  gives quick-check without data lists.
c                          then FH =-9.  reverts to normal output.
C-------------------------------------------------------------------------------
      DIMENSION FV(399), HT(399), QQ(50)
      data ndim / 399 /
      character head*25, dat*8
      character(256) :: buf
      character(:), allocatable :: datin
      integer nargin


      common /test/ test1, test2, test3              !put in POLAN... for debugs
    8 format (/, a,/ ( f8.2, 9f10.2) )
    9 format (/, a,/ (f6.2,f6.1, 6(f8.2, f6.1)) )
      N = ndim                                                   !set array size
      call date_and_time (dat)

      nargin = command_argument_count()
      if (nargin.eq.0) stop('must specify input file')

      call get_command_argument(1,buf)
      datin = trim(buf)

      OPEN (UNIT=1, FILE= datin, STATUS='OLD', action='read')
      OPEN (UNIT=2, FILE='out.dat',status='replace', action='write')
      write (2,40) datin, dat
40    format (' P O L A N    of  FEBRUARY 1993.     ',
     &        'Data file: ',a, '   Run: ',A8,'.')
cc      QQ(1) = -1.          ! no coefficients back
10       quik = 0.         ! quick-look off
         GO TO 100
20          quik = 1.      ! quick-look on
c                                      read field and mode
100   READ (1,120) HEAD, FH,DIP, AMODE, VALLEY, LIST                    read (1)
120   format (a25, 4f5.0 , i5, f5.0)
              if (fh.eq. 0.)  stop
              if (fh.eq. 9.)  go to 20
              if (fh.eq.-9.)  go to 10
c
150   WRITE (2,160)
      WRITE (2,160) HEAD, FH,DIP, AMODE, VALLEY, LIST
      write (*,160) head, fh,dip, amode,valley,list                         show
160   format (a25,'   FH',f5.2, '  Dip',f5.1,
     1       '    Amode',f5.1,'  Valy',f6.2,'  List',i2)
c
c-------------------------------------- read data - loop for new ionogram
200   READ (1,220) HEAD, START, (FV(I),HT(I), I=1,5)                     read(2)-
220   format (a25, f5.3, 5(f5.3,f5.2) )
         if (fv(1).eq.0.) go to 100                    !read new field constants
      WRITE (2,240)
      WRITE (2,240) HEAD, START
      write (*,240) head, start                                             show
240   format (a25,'   Start =',F7.3)
      NH = 5
300      IF (FV(NH)+HT(NH).eq.0..and.HT(NH-1).eq.0.) GO TO 350               end
         IF (FV(NH).eq.-1. .OR. NH.gt.NDIM-40)  GO TO 350                    end
         READ (1,320) (FV(I),HT(I), I= NH+1, NH+8)                      -read(2)
320      format (8(f5.3,f5.2))
         NH = NH+8
         GO TO 300
350         nh = nh-1
            if (fv(nh)+ht(nh-1).eq.0.)  go to 350
      if (quik<1.) write (2,9) 'Input data', (fv(i),ht(i), i= 1,nh)         list
c-------------------------------------- analysis
c
            if (quik.ge.1.)  list = 0
      CALL POLAN (N,FV,HT, QQ, FH,DIP, START, AMODE, VALLEY, LIST)
c
c-------------------------------------- output
               if (quik.ge.1.) go to 200
      write (2,9) 'Real Heights', (fv(i),ht(i), i= 1,n+2)
c
      if (qq(1).gt.2.) write(2,8) 'Coefficients QQ',
     &                            (qq(i), i=1, ifix(qq(1))+1)
      write (2,*)'===================================================='
      GO TO 200

      END PROGRAM
