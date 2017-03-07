            PROGRAM  POLRUNEW
c21/10/96 =POLRUN, but read free-format data, with known number of data points.
c  (e.g. Bill Wrights data lists); and O/P=> 80 cols.  Initial Header line.
c
c  1st line (field, station) as before;  2nd line = Name, start, nh (a25,f5,i5).

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
      character head*25, dat*8, datin*15
      common /test/ test1, test2, test3              !put in POLAN... for debugs
    8 format (/, a,/ ( f8.2, 7f10.2) )
    9 format (/, a,/ (f6.2,f6.1, 4(f8.2, f6.1)) )
      N = ndim                                                   !set array size
      call date_and_time (dat)
      call get_command(datin)
      if (datin.eq.' ') print*, 'Name of data file : '
      if (datin.eq.' ') read *, datin
      OPEN (UNIT=1, FILE= datin, STATUS='OLD')
      read (1,*)
      OPEN (UNIT=2, FILE='polout.t')
      write (2,40) datin, dat
40    format (' P O L A N    of  FEBRUARY 1993.     ',
     &        'Data file: ',a15, '   Run: ',A8,'.')
cc      QQ(1) = -1.          ! no coefficients back
10       quik = 0.         ! quick-look off
         GO TO 100
20          quik = 1.      ! quick-look on
c                                      read field and mode
100   READ (1,120) HEAD, FH,DIP, AMODE, VALLEY, LIST                   !read (1)
120   format (a25, 4f5.0 , i5, f5.0)
              if (fh.eq. 0.)  stop
              if (fh.eq. 9.)  go to 20
              if (fh.eq.-9.)  go to 10
c
150   WRITE (2,160)
      WRITE (2,160) HEAD, FH,DIP, AMODE, VALLEY, LIST
      write (*,160) head, fh,dip, amode,valley,list                     !show
160   format (a25,'   FH',f5.2, '  Dip',f5.1,
     1       '  Amode',f5.1,'  Valy',f6.2,'  List',i2)
c
c-------------------------------------- read data - loop for new ionogram
c200   READ (1,220) HEAD,START, (FV(I),HT(I), I=1,5)                    !read(2)-
c220   format (a25, f5.3, 5(f5.3,f5.2) )
200   READ (1,210) HEAD,START, NH
210   format (a25, f5.0, i5)
         if (NH.eq.0) go to 100                    !read new field constants
      WRITE (2,240)
      WRITE (2,240) HEAD, START,  NH
      write (*,240) head, start,  nh                                     ! show
240   format (a25,'   Start =',F7.3, '    NH =',i4)
      READ (1,*) (FV(I),HT(I), I= 1, NH)                                !-read(2)
      if (quik<1.) write (2,9) 'Input data', (fv(i),ht(i), i= 1,nh)        !list
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
      END

