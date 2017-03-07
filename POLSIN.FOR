c  POLSIN.FOR  =  SETUP,  SELDAT,  STAVAL.
c****************************************  -  input processing for Polan.
c
      SUBROUTINE   SETUP  (FV, HT, NDIM, START)
c ---                                (called from section c1.2 of polan)  feb84.
c    NDIM is the dimension of the data arrays FV and HT.
c
c 1. Move the virtual height data up 44 places, to start at  fv(45), ht(45).
c 2. Set  fmin =  the lowest plasma frequency in the data.
c         hmin =  the lowest virtual height in the data.
c         (these are used for start extrapolation in the absence of x-ray data).
c 3. Identify the start method, and set the following parameters:
c      js = the number of points added below fmin, = 2/1/0 for oray/xray/direct.
c      fa, ha  define the starting point for the first polynomial.
c      hs =  the height to evaluate the gyrofrequency for an x-ray start.
c      lk = 1 / 0 / -1   for   o-ray / x-ray poly / x-ray slab   starts.
c The real and virtual origins are at  kr = 1  and  kv = 45-js.
c-------------------------------------------------------------------------------
Constants which SETUP sets in /pol/ are:-  lk,kr,kv, fa,ha, hs,vdepth.
c
      DIMENSION  FV(*), HT(*)
      COMMON /POL/ B(99,20),Q(20), FH,ADIP, MODE,MOD, FA,HA, tcont,lbug
      COMMON /POL/ HS, FC,FCC, HMAX,SH, PARHT, HVAL,VWIDTH,VDEPTH, XWAT
      COMMON /POL/ MAXB,NF, NR,NL, NX, MS,MT,JM, LK, KR,KRM, KV,MF, NC  counters
      VDEPTH = 0.
      VWIDTH = 0.
      TCONT = 0.
      SH = 0.
      FC = 0.
c1.a --                        Move virtual data to start at  fv(78), ht(78)
c1.a --                        Find index j of the 1st o ray; fmin= lowest fn
      DO 10 K = NDIM-44, 1, -1
         FV(K+44)= FV(K)
         HT(K+44)= HT(K)
            IF (FV(K).GT.0.)  J= K
   10    CONTINUE
12    FMIN = FV(J)                                                      1st oray
      IF (J.GT.1)  FMIN = AMIN1(FMIN,SQRT(FV(1)*(FV(1)+.9*FH)))         1st fn
c1.b --
c1.b --                     Set starting conditions (direct, extrap or model).
c                              Calculate initial fa, ha;- start from fmin:
      HMIN = HT(J)
      DO 20 I = 1, 5
         IF (HT(J+I).LT.45.)  GO TO 22                                  endlayer
20       HMIN = AMIN1(HMIN, HT(J+I))                                       h'min
22    HA = HMIN
      FA = FMIN
      HS = HA
      JS = 0
      LK = 1
      IF (START.LT.0..AND.START.GE.-1..AND.J.EQ.1)  GO TO 40              direct
c
c#12'95    O ray Start extrapn:  Omit initial retarded points from ho calcn.
           h1max= min( hmin+40., ht(2)+200*(fv(2)-fv(1)) )
           if (j.eq.1.and.ht(1).gt.h1max.and.I.gt.4) then
              do 24 i = 1, 43
                 fv(i) = fv(i+1)
   24            ht(i) = ht(i+1)
              goto 12
            endif
c                                                       - start extrapn:
      JS = 2                                                            add 2pts
      if(fv(j+2).gt.0.)DH = ABS(HT(J+2)-HT(J)) *FMIN /(FV(J+2)-FMIN)    to f = 0
         IF (FV(J+2).LE.0..or.HT(J+2).LE.0.)  DH = HT(J)/6.
      HS = AMIN1(HMIN-DH, HMIN/2.+50.)                                  start fh
      HS = AMAX1(HS, HMIN/4.+55.)                                       lowlimit
      FA = AMIN1(.5, .6*FMIN)
c                                                       - model  height:
      IF (ABS(START).GE.45.) HS = AMIN1( ABS(START), HS*.4+HMIN*.6 )       model
      if (j.le.1)  then                                                  o start
         HA = HS                                                        start ht
         IF (MODE.EQ.10 .OR. MODE.GE.20)  FA = FMIN*.6                  snglpoly
         IF (START.GE.45..OR.START.EQ.0.) GO TO 40                       o start
c                                                       - model density:
           H = INT(START*.1)*10
           HA = 90.+H+H
           FA = START-H
       else
c                                                       - x-ray calculation:
         JS = 1                                                         add 1pt.
         HA = HS*.4 + HMIN*.6                                            x start
           LK = 0
           IF (START.GT.-3.)  FA = -START - 1.                          polystrt
           IF (START.LT.-1.)  GO TO 40
         LK = -1
         FA = FMIN*.6                                                   slabstrt
         VDEPTH = FMIN*.3
        endif
c1.c --
c1.c        Store initial points. Virtual data starts at index kv (43 to 45).
40    KV = 45-JS
      FV(44) = (FA + FMIN)/2.0                                          o start,
      HT(44) = HS + (HMIN-HS)*FV(44)/FMIN                               add virt
               if (j.eq.1)  hs= ht(44)                                  (& show)
      FV(KV) = FA
      HT(KV) = HA
      FV(1) = FA
      HT(1) = HA
      KR = 1                                                            realstrt
      RETURN
      END
cc******************************************************************************
      SUBROUTINE  SELDAT  (NV, FV, HT)
c                                   called from section c2.2 of polan.    feb84.
c     Select the virtual-height data to be used in the next step of polan.
c     Define the following parameters:-
c  nf = number of o-ray data points to be used
c       (equal to nv, if sufficient good points exist before a peak/restart);
c  nx = number of x-rays points;
c  fm = fv(mf) = highest frequency included in this step;
c  fcc is set equal to zero for a normal (non-peak) step;
c             equal to  -.1 for  a cusp at fm;
c             equal to a scaled critical frequency (fc),  or to 0.1, for a peak.
c-------------------------------------------------------------------------------
      DIMENSION  FV(*), HT(*)
      COMMON /POL/ B(99,20),Q(20), FH,ADIP, MODE,MOD, FA,HA, tcont,lbug
      COMMON /POL/ HS, FC,FCC, HMAX,SH, PARHT, HVAL,VWIDTH,VDEPTH, XWAT
      COMMON /POL/ MAXB,NF, NR,NL, NX, MS,MT,JM, LK, KR,KRM, KV,MF, NC  counters
c
c   ffit gives the minimum desirable frequency range, for fitting x-ray data.
c   gfit gives the maximum slope  dh'/df  which should be used in an x-start.
      data   ffit,  gfit    /  0.411, 40.1 /
c-------------------------------------------------------------------------------
c2.a                      Count initial x rays. Check frequency sequencing.
c2.a                      Check for cusp, peak, or end of data.
2     FCC = 0.                                                           no peak
      FSX = 0.
      FRX = 0.
      MF = KV
      FS = FV(KV)                                                        start f
      F1 = 0.
      FH = GIND(0., HT(KR))
      NF = 0
c.....                                 Frequency loop  (nf = 1 to nv):
10    NF = NF+1
12       MF = MF+1                                                      toppoint
         FM = FV(MF)                                                    top freq
         FN = FV(MF+1)                                                  nextfreq
         HV = HT(MF)                                                    top virt
                      if(lbug.eq.6.and.kr.eq.1)write(2,*)nv,nf, f1,fm   ##----->
         IF (FM.GT.FS)  GO TO 20
            FRX = SQRT(FM*(FM+FH))                                      top xray
            IF (MF.EQ.KV+1)  FSX=FRX                                    1st xray
            IF (FM.LT.0..AND.FRX.GE.FSX.AND.NF.EQ.1)  GO TO 12            x ray.
15                 nv = -mf
                   return                                               error>>>
20         IF (F1.EQ.0.)  F1 = FM
           IF (HV.LT.0..OR.FN.LT.0.)  GO TO 40                          cusp,end
           IF (ABS(HT(MF+1)).LE.50.)  GO TO 50                              peak
         FS = FM
                   if (fn.le.fm)  go to 15                              error>>>
      IF (NF+1.LT.NV.AND.MF-KV.LT.MAXB-8)  GO TO 10
c.....
c                                      Check final point
           GRAD = (HT(MF+1)-HV) / (FN-FM)
      IF (NF+1.EQ.NV.and.(NF.LE.2.or.KR.GT.1.or.GRAD.LT.100.)) GO TO 10    add 1
c
c                           Leave loop, with nf = no of o rays (.le.nv).
c2.b
c2.b                                   Check need for additional points.
      IF (MOD.GT.10.OR.MOD.LT.4)  GO TO 70                              end freq
         IF (HT(MF+1)-HV.GT.(FN-FM)*GFIT)  GO TO 70                     retarded
         IF (FRX.GT.FSX  .AND. FM.GT.FRX)  GO TO 70                     all x in
         IF (FM-F1+.04*(nf-nv).LT.FFIT .OR.FM.LT.FSX)  GO TO 10          incr fm
      GO TO 70
c                                      Cusp  (fcc =-0.1) or end data
40          FCC = AMIN1(FN, -.1)                                        discont.
            HT(MF) = ABS(HV)                                             (=cusp)
            GO TO 70
c                                      Peak  (fcc = 0.1  or scaled fc)
50          FCC = AMAX1(FN,.1)                                           f crit.
c
c2.c       Shift data arrays to delete out-of-range x rays (fnr > fm+0.1).
c2.c
cc70         if (nf.lt.1.or.(nf.lt.2.and.fcc.eq.0.1))  go to 15          end>>>>>
70         if (nf.lt.1) goto 15     !!!2'93 allow peakfit, seems OK??  !error>>>
      MX = MF - NF                                                      lastxray 
      KV1= KV
80    NX = MX - KV                                                      # x rays
      IF (NX.EQ.0)  GO TO 100
         F = FV(MX)
         IF (SQRT((F+FH)*F).GT.FM+0.1) GO TO 90                          delete
         IF (NX.EQ.1)  GO TO 100
         IF (HT(MX)-HT(MX-1).LT.GFIT*(FV(MX-1)-F))  GO TO 100           grad ok.
90          J = MX
95             FV(J)= FV(J-1)                                           move up1
               HT(J)= HT(J-1)
               J = J-1
               IF (J.GE.KV)  GO TO 95
            KV = KV+1
            IF (KV.GT.KV1+NV)  GO TO 2                                  redo all
            GO TO 80
100   IF (NX.EQ.0.AND.KR.EQ.1)  LK = 1                                  no xstrt
      RETURN
      END
cc******************************************************************************
      SUBROUTINE  STAVAL  (FV,HT, DIP,VALLEY, LOOP)
c                                                                         feb84.
c      Staval handles all start and valley calculations, in which the 
c      real-height polynomial has a constant term  [at q(jm), where jm=mt+1]  
c      so that the height at the origin  [fa,ha]  is recalculated.
c-------------------------------------------------------------------------------
c--> The initialising call enters and exits with  loop = 0.
c       If this is not a valley calculation (i.e.  hval = 0.)  then staval just
c           sets the initial weight  xwat  for any x-ray data, and returns.
c       For a valley calculation (shown by hval.ne.0.)
c           the iteration-flags  devl, devll  are initialised,  and
c           the standard values for valley depth and width are calculated
c                           (based on the 'standard scale height'  sha).
c
c--> A call with  loop = 1  is to adjust the valley depth for a single o-ray
c                           iteration.  this exits with  loop = 4.
c-------------------------------------------------------------------------------
c--> The main call enters with  loop = 2.            ( -2 if fb is negative. )
c      Exit with  loop = 2  shows that the start/valley calculation is completed
c                           and the required real heights have been added to ht.
c      If  loop > 2  on exit,  polan recycles the valley calculation:-
c          loop = 3  requests a new calculation with a changed gyrofrequency.
c          loop = 4  is to recalculate with a new value of valley depth  vdepth.
c-------------------------------------------------------------------------------
      DIMENSION  FV(*), HT(*)
      COMMON /POL/ B(99,20),Q(20), FH,ADIP, MODE,MOD, FA,HA, tcont,lbug
      COMMON /POL/ HS, FC,FCC, HMAX,SH, PARHT, HVAL,VWIDTH,VDEPTH, XWAT
      COMMON /POL/ MAXB,NF, NR,NL, NX, MS,MT,JM, LK, KR,KRM, KV,MF, NC  counters
constants:    base,    depthfac,   depthconst
      DATA    VBASE,    VDEEP,       VCONST     / 0.6, .008, 20. /
constants:  1st depth,  scale,  peak-top increase           for valley iteration
      DATA    VAL1,     DVAL1,       VPEAK      / .1001, 6., 1.4 /
constants:   xweight,    fhht itern limit
      DATA     WVX,          HXERR              / 1.0,  2.0 /
c
      SAVE  VAL, DVAL, DEVL, DEVLL, HDEC, FHHT
      sq(x) = sqrt((1.-x)*(1.+x))
c-------------------------------------------------------------------------------
          NFB  = LOOP                                                   fix fb ?
          LOOP = IABS(LOOP)
      IF (LOOP.GT.1)  GO TO 200
      IF (LOOP.EQ.1)  LOOP = 4                                          o ray
      IF (LOOP.EQ.4)  GO TO 40
c###############################################################################
cx1---                       First call (with loop = 0) from polan section c3.
          tras = -3.1
      XWAT = WVX
      PARHT= 0.
cx1.a                                  Check for valley, set flag hval
      HVAL = HT(KV)                                                       valley
          IF (HVAL.EQ. 0.) HVAL = VALLEY
          IF (HVAL.EQ. 0.) HVAL = 1.
          IF (HVAL.GE.10.) HVAL = 0.                                    novalley
      IF (HVAL.EQ.0.)  RETURN                                           novalley
c-------------------------------------------------------------------------------
c                                      Initialise valley width and depth
      HDEC = (INT(HVAL) - HVAL) *2.
      DEVL = 1.E6
      DEVLL= 1.E7
      DVAL = DVAL1                                                      valscale
      HMAX = HT(KR)
      SHA = HMAX*.25 - 20.                                              scale ht
      HS = HMAX + SHA                                                   vally fh
      SH = SH *VPEAK                                                     topside
cx1.b                        Set standard width
20     VWIDTH = (SHA + 7) * ABS(HVAL)                                   standard
            IF (HVAL.LT.-2.)  VWIDTH = INT(-HVAL) *2.                   specifid
c-----------------------------------------------------------------------
c                   Set initial depth,  or depth adjustment for valley iteration
cx1.c                        Standard depth (width calculations iterate here).
c
40     VAL = VWIDTH**2 *VDEEP/(VCONST+VWIDTH)                           v. depth
           if (hmax.gt.140.)  val = 0.5*val                             F1 vally
        IF (HDEC.GT.0.)  VAL = HDEC                                     specifid
            IF (HVAL.EQ.-1..AND.NX.GT.0)  VAL = VAL1                    1stdepth
            IF (INT(HVAL).EQ.-1.)  HVAL = -1.                            iterate
c
cx1.d                        Peak section  (depth calculations iterate here).
60     vdepth = val
       if(hdec.le.0..or.val.gt..8*fa) VDEPTH = VAL*FA/(VAL+FA)          depth<fa
        PARHT = 2.0*SH* SQ(1.-VDEPTH/FA)                                par.vall
        HA = HMAX + PARHT                                               valy bot
        IF (NC.EQ.0)  FHHT = HA+20.                                     field ht
      RETURN
c###############################################################################
c###############################################################################
c                   Main call  with  loop = 2,  from polan section c4.3.
c                              Start / valley  entry,  after  "solve".
200   tras = -4.3
      DEVN = Q(20)
      FNL = 0.
      DHL = Q(JM)
c
cx2--                                  Calculate and check real heights
220   IF (NX.LE.0..OR.(DEVN.EQ.0..AND.HVAL.EQ.-1.)) DVAL = -1.          list,end
      QM = Q(JM)
      HR = HA + QM
      HN = HR
      FN = FA
      MC = NC
      MV = NF + NX
      MQ = MT + MIN0(LK,0)
c.....                                     (real hts stored for coefic recycle)
      DO 300  I = 1, MV
         HL = HN
         FL = FN
         FN = FV(KV+I)
         IF (FN.LT.0.)  FN = SQRT((FN+FH)*FN)
         DELTF = FN - FA
         HN = HR + SUMVAL(MQ,Q,DELTF,1)
c
            if (nc+mc.ge.25.or.dip.lt.0.)  go to 290                    omittest
            if ( (hn-hl)/(fn-fl) .gt. 2.0)  go to 290                   grad ok.
c                                               dh / dfn  .lt.  2.0
            mc = mc + 3
            if (nx.gt.0.and.xwat.eq.wvx.and.nc.gt.1) go to 380              redo
               qmm = qm-amax1(5.,abs(qm)/4.)*amax0(1,3-kr)               next qm
               if (fn.ne.fnl.or.hn-hl.gt.dhl) go to 350                 lower qm
                  if (lbug.ge.-9) write(2,250) fv(kv+i)                 **----->
250               format('Data/Gyrofrequency incompatible at f ='f6.2)
                  mc = 30
290      IF (KR+I.LT.KV)  HT(KR+I) = HN
300      CONTINUE
c.....
      KRM = MIN0(KR+MV, KV-1)                                           topcalcd
cx3---                                 Set field ht  (and recycle fb).
      IF (NC.LE.2)  XWAT = WVX
      KHX = KR + 1 + NX/3
      DHX = HT(KHX) - FHHT
      FHHT = HT(KHX)
         IF (NC.EQ.1 .OR. ABS(DHX).GT.HXERR)   LOOP = 3                 reset fb
         IF (NX.LE.0.OR.NFB.LT.0.OR.NC.GE.20)  LOOP = 2                 no reset
      IF (KR.GT.1)  GO TO 500                                             valley
      GO TO 400                                                          x start
c-------------------------------------------------------------------------------
c               (jump here to get new solution with reduced start offset)
350              call solve (-900, -jm, b,q, qmm, devn,lbug)            decr. qm
               if(lbug.gt.0) write(2,360) nc,qm,q(jm), fn,hn,q(19),devn **----->
360      format(i3,' STAVAL: Qm reduced from'f6.1,' to'f6.1,', to avo',
     &  'id -ve slope at f,h ='2f7.2,'   Devn increases'f6.2,' to'f6.2)
                 fnl = fn
                 dhl = hn-hl
                 go to 220
c-------------------------------------------------------------------------------
c                                           Re-calculate a bad result
380                       xwat = wvx/2.                                 bad data
                          if(lbug.ge.-9) write(2,390)                   **----->
390                       format (' X ray Weights Reduced to 1/4.')
                          loop = 4
                          return
c*******************************************************************************
cx4---*******************************  X ray Start: slab (.3-.6f1) +poly********
cx4---
400   SLAB = Q(MT)
         IF (LK.EQ.0)  SLAB = 0.                                        polystrt
          if(lbug.ge.-9)write(2,420) nc,qm, slab,devn, jm,nf,nx,ms,fhht ------->
420       format (i2,' Start Ofset',f6.1,' km,  Slab',f6.1,' km.  Devn',
     1 f5.2,' km.',i3,' terms fit',i3,' O +',i2,' X +',i2,'.  hx=',f6.1)
      IF (LOOP.EQ.3)  RETURN                                            reset fb
      KD = 1-LK                                                         addedpts
      GO TO 700
c*******************************************************************************
cx5---*******************************  Valley:  Iterate and print       ********
cx5---
500   VWIDTH = QM + PARHT
            if (hdec.gt.0.)  dval= -hdec                                no itern
            if ((lbug.gt.0.or.dval.lt.0..or.nc.ge.20).and.lbug.ge.-9)
     1         write(2,520) nc, vwidth,vdepth, devn,jm,nf,nx,ms, fhht   ------->
520         format ( i2,' Valley',f6.1,' km wide,',f6.2,' MHz deep.',
     1              '  Devn',f5.2,' km.', i3,' terms fit',
     2              i3,' O +',i2,' X +'i2,'.  hx=',f6.1)
      IF (LOOP.EQ.3) RETURN                                             reset fb
c
      KD = 4
      LOOP = 4
      VBOT = QM*VBASE
      SLAB = QM-VBOT                                                    lin.slab
         IF (DVAL.LT.0..OR.NC.GT.25)  GO TO 700                         end valy
         IF (HVAL.EQ.-1.)  GO TO 600                                     iterate
c
c - - - - - - - - - - - - Normal cycle: Scale depth as (width)**2  - - - - - - -
      IF (DEVN.GT.DEVL/.9.and.nc.lt.4)  GO TO 20                          revert
      IF (DEVL.LT.1.E6)  GO TO 560                                      end valy
          DEVL = DEVN
          GO TO 40                                                      loopvaly
560   DVAL = -1.                                                        
          if (lbug.eq.0)  go to 500                                     listvaly
      GO TO 700                                                         end valy
cx6--- - - - - - - - Vdepth Iteration: choose depth val1 or val1*dval1. - - - -
cx6---
600     IF (DVAL.NE.DVAL1) GO TO 620                                     iterate
        IF (VAL.EQ.VAL1)   GO TO 675                                    1stround
        DVAL = 1.+.5*COS(.015*ADIP)                                     2ndround
        IF (DEVN.LT.DEVL)  GO TO 675                                    continue
          VAL = VAL1                                                      revert
          GO TO 680
c                                      Iterate from chosen depth
620       IF (DEVN.GT.DEVL*.97-.003.AND.DEVLL.LT.1.E6)  GO TO 660            end
          IF (DEVN.LT.DEVL)   GO TO 670                                 continue
            DVAL = 1./DVAL                                              reverse,
            VAL = VAL*DVAL                                              2ndentry
            DEVLL = DEVN
            GO TO 680
c - - - - - - - - - - - -              Interpolate to minimum deviation
660       DMIN = (DEVLL-DEVN)/ABS(DEVLL+DEVN-DEVL-DEVL)*.5  - 1.
          DVAL = -DVAL**AMIN1(DMIN, .5)
          GO TO 680                                                         last
c                                      Recalculate valley width
670     DEVLL= DEVL
675     DEVL = DEVN
680   VAL = VAL*ABS(DVAL)                                                   next
      GO TO 60
c------------------------------------------------------------------------------
cx7---******************************   Tidy up.    Delete x rays.   ************
cx7---             Add initial points;   kd = 1,2,4  for poly, slab, valy
700   LOOP = 2                                                          all done
                       if (lbug.eq.6) call trace (fv,ht, tras)          ##----->
      KR = KR + KD                                                      stepreal
      KV = KV + IABS(NX)                                                delete x
      IF (KD.EQ.0)  RETURN
      IF (KD.EQ.1)  GO TO 750                                           polystrt
c                                      Store valley,start slab at kd = 4,2
        FBOT = FA - VDEPTH
        HT(KR-2) = HA                                                   v.bottom
        FV(KR-2) = FBOT                                                 v.bottom
        HT(KR-1) = HA+QM*VBASE                                          v.bottom
        FV(KR-1) = FV(KR-2)                                             v.bottom
        IF (KD.LT.4)  GO TO 750
            HT(KR-3) = HA - 0.5*PARHT                                   peak top
            FV(KR-3) = 0.5*SQRT(3.*FA**2 + FBOT**2)                     peak top
            FV(KV) = FA
c
750   HA = HA + QM
      IF (KD.LE.2)  KR = KD
      HT(KR) = HA
      FV(KR) = FA
      IF (KR.GT.3)  GO TO 800                                           notstart
c                                      Add points on start polynomial
          HT(1) = HA - SLAB
          FV(KV-1) = FA                                                  calc ht
          F1 = FV(KV+1)                                                 1st freq
            IF (LK.EQ.0)  FV(KV-1) = (FA+F1)*.5                          (poly)
          FV(KV) = (FV(KV-1)+F1)*.5                                      calc ht
          KV = KV-2-LK
          VBOT = 0.
800   IF(KR.GT.1) TCONT=TCONT+FV(KR-1)**2*(VBOT+.5*SLAB)+FA**2*.5*SLAB  lin.vall
     1                  + PARHT*(FA*FA-VDEPTH*(FA+FV(KR-1))/3.)         par.vall
                      if (lbug.eq.6) call trace (fv,ht,tras)            ##----->
                      if (lbug.eq.6) lbug = 0                           end strt
                      if (lbug.eq.9) lbug = -9                          end valy
      RETURN
      END

