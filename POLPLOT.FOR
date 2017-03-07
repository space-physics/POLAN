      PROGRAM  POLPLOT
c  Read data (as polrun);  Plot h', f;  Analyse;  Plot h', hr.    Aug'85/Feb'93.
c                          Plots calculated heights, & curve from QQ when given.
c  1'93: to new graf library.   Show Chapman peak to ñ1.5SH.
c  2'93: amode=-99. to read mode array  domode(12).
      common /test/ test1, test2, test3              !put in POLAN... for debugs
      dimension fv(380), ht(380), qq(90),  fs(120), hs(120)
      dimension ph(120), phx(80), pf(120), pfx(80)                  !plot arrays
      character dat*8, datin*15, head1*25, head2*25, text*40
      data ndim, qq(1), qq(90) / 120, 1., -1. /
      data tfmin,tfmax, thmin,thmax / 1, 8, 80, 350 /                   !default
      real domode(12) / 12*0. /
      n = ndim
    1 format (i9)
    9 format (/, 1x,a,/ (f6.2,f6.1, 6(f8.2, f6.1)) )
      test1 = 0.
      test2 = 0.
      call getcl(datin)
      if (datin==' ') then
         print*, 'Name of input data file : '
         read *, datin
         endif
      open (unit=1, file=datin, status='old')
      print *, '0 = no output,  1 -> polplot.out,  2 -> printer, '
      print*,'    3,4 add plot points,  >4 debug.  -ve set list  ? '
      read "(9f9.0)", quik, test1, test2, test3
      lout= iabs(quik)
      pfmin= 0.
      pfmax= 0.
         print*, 'fmin, fmax, hmin, hmax [0def,1auto] ? '
         read"(9f9.0)", pfmin, pfmax, phmin, phmax
      if (lout==2.or.lout==4) then
            open (unit=2, file='lpt1')
         else if (lout/=0) then
            open (unit=2, file='polplot.out')
         endif
      call date (dat)
      if (lout>0) write (2,40) datin, dat
40    format ('POLPLOT:-  P O L A N  of  FEB. 1988.',
     &        '   Data file: ',a15, ' Run: ',A8)
         if (quik<0.) write (2,9) '// QUIK-LOOK RUN //'
c                                      read field and mode
100   read (1,120) head1, fh,dip, amode, valley, list                   READ (1)
120   format (a25, 4f5.0 , i5, f5.0)
         if (fh==0.) STOP 'RUN ENDS ##'
         if (quik<0.) list = lout
      if (lout>0) write (2,*)' ----------------------------------------'
      if (lout>0) write (2,160) head1, fh,dip, amode, valley, list
160   format (/,a25,3x, 'fh',f5.2, '  Dip',f5.1,
     1       '   Amode',F6.1,'  Valy',f6.2,'  List',I2,f8.1)
      if (lout==0) list = 0
      if (n==ndim.and.list>1) write(2,9)
     &               char(27)//'x'//char(0)//char(27)//'!'//char(4)       !17cpi
      if (amode==-99.) read(1,"(12f9.0)") domode
c
c======================================    READ DATA - loop for new ionogram
200   nh = -3
      read (1,220) head2,start, (fv(i),ht(i), i=1,5)                    READ(2)-
220   format (a25, f5.3, 5(f5.3,f5.2) )
            if (fv(1)==0.) go to 100                   !read new field constants
      if (lout>0) write (2,240) head2, start
240   format (/,a25,4x,'Start ='f7.3)
300      nh = nh+8
         if (fv(nh)+ht(nh)==0..AND.ht(nh-1)==0.)  go to 350                  END
         if (fv(nh)==-1. .or. nh>ndim-40) go to 350                          END
         read (1,320) (fv(i),ht(i), i=nh+1,nh+8)                        -read(2)
320      format (8(f5.3,f5.2))
         go to 300
350         nh = nh-1
            if (fv(nh)==0.)  go to 350
c                                         LIST DATA  (and save to re-run)
      if (quik>0.) write(2,9) 'Input data:', (fv(i),ht(i), i= 1,nh)
         nsav = nh + 2
         call amove (nsav, fs, fv)
         call amove (nsav, hs, ht)
c                                         FORM PLOT ARRAYS
      nx = 0
      nv = 0
      hvmax= 0.
      ffh = fh*.92                                               !fh near 160 km
      if (fh<0.)  ffh = -fh                                         !fh constant
      do 480 i = 1, nh
         if (fv(i)<0.) then
            nx = nx+1                            !plot x data at reflection freq
            pfx(nx)= sqrt( (fv(i) + ffh) *fv(i) )
            phx(nx)= ht(i)
         else if (fv(i).ne.0.) then
            nv = nv+1                                                !o ray data
            pf(nv)= fv(i)
            ph(nv)= ht(i)
         endif
480      hvmax= amax1(hvmax, ht(i))
c========================================       ANALYSIS
      if (lout==0) list = -10
      modedo= 0
      if (amode==-99.)  modedo= 1
      if (amode==-99.)  amode = domode(modedo)
c
500   CALL POLAN (N,FV,HT, QQ, FH,DIP, START, AMODE, VALLEY, LIST)
c
      fc = fv(n-3)
      hmax= ht(n-3)
      sh = abs(ht(n+2))                                    !to show Chapman peak
c                                               LIST OUTPUT
      if (quik>0.) write(2,9) 'Calcd Heights:', (fv(i),ht(i), i= 1,n)
c
      if (qq(1)>2..and.list>=0.and.lout>0)
     &            write (2,620) 'Array QQ', (qq(i), i=1, ifix(qq(1))+1)
  620 format ( a / (f8.3, 9f10.3))
      if (fv(n).le.0.)  n = n-1
c========================================= PLOT VIRTUAL HEIGHTS, as o or x.
      if (pfmax<pfmin)  then
         pfmin= int((fv(1) -.09)/.20) *.20                      !get plot scales
         pfmax= int((pf(nv)+.75)/.50) *.50
         phmin= int((ht(1) -10.)/20.) * 20.
         phmax= int((ht(n) +30.)/20.) * 20.
         if (pfmax>8.*pfmin) pfmin= 0.
       else if (pfmax==0.) then
         pfmin= tfmin
         pfmax= tfmax     !or min(, fc*1.1)
         phmin= thmin
         phmax= thmax     !, hmax+sh*.1)
       endif
      call plstart (pfmin, pfmax,  phmin, phmax, 0.)
      call plabel ('POLPLOT/'//datin, 'f')
      call plabel (head1, 'h')
      call plabel ('Frequency, MHz.', 'x')
      call plabel ('Height, km.',     'y')
      call plabel (head2, 'tl')
      write(text,"(a,f6.1,a,i3,a,i3,a)") 'Mode=',amode,';',
     &               int(qq(4)),' poly terms.'
      call plabel (text, 'tl')
      write(text,"(a,f5.2,f6.1,f5.1)") 'fc,hm,H=', fc, hmax, sh
      call plabel (text, 'tl')
      call setpen (2.)                                                    !green
      call plsymb (nv, pf, ph,  11)                                !vh data= o,x
      if (nx>0)  call plsymb (nx, pfx,phx, 15)
      call plabel ('Virtual Heights = o, x.', 'tl')
c                                         PLOT PROFILE, as +
      call plabline ('Calculated Heights', 'br', 1., 1., 24)
      call plotn (n, fv, ht, 241.)                             !calcd hts= --+--
c====================================== END
950   if (modedo>0.and.modedo<12) then
         modedo= modedo + 1                                    !cycle read modes
         amode= domode(modedo)
         if (amode/=0.) goto 960
         endif
      modedo= 0
      amode= 0.
cc      call plabel ('NEXT AMODE ?', 'tl')
cc      read"(f9.0)", amode
960   call plend                                    !end plot, close output file
      if (amode==0.) go to 200                                 !read more data
         call amove (nsav, fv, fs)                             !or restore data,
         call amove (nsav, ht, hs)
        if(lout>4) write (2,9) 'Input data', (fv(i),ht(i), i= 1,nsav)
         goto 500                                              ! & run new MODE.
      end

