      SUBROUTINE POLAN (N,FV,HT, QQ, FB,DIP, START, AMODE, VALLEY, LIST)
c N.B.--> 2'93  delete ndim from call;  set by N on 1st call.
c
c  A generalised POLynomial real-height ANalysis, for ionograms.  May 1986.
c
c  BEFORE USING POLAN STUDY THE DESCRIPTION OF THE PARAMETERS, AND THE
c  NECESSARY ARRAY DIMENSIONS, AS GIVEN IN THE TEXT FILE  POLAN.TXT.
c          (Please list this file, and keep it by you for reference.)
c-CHANGES-----------------------------------------------------------------------
c 12'95:    Move data up 44 places (was 30) to ensure room for complex calcns.
c      Mode 9 to NT= 6 (was 7).  Max NT= 14 (was 15, which sometimes distorted).
c    POLSIN: Omit initial retarded points from ho calcn, for better O starts.
c            Valley Width to SHA +7 = hm/4-13.  => F1 Width*2/3, Depth*0.5.
c            Valley depth halved for F1 (hmax > 140 km).  => typical 0.07 MHz.
c  2'93  delete ndim;  1st call needs n = ndim.  (& ndim > #data+45; 12'95)
c  3'88  MT +1 for physical equns, but slower increase with NX.
c  8'87  to fully FORTRAN77.   9'86 correct valley calcn at mode 1.
c  5'86  allow dip= 0.         4'86 add NDIM and QQ,  add modes xyz.
c  PEAK revised aug87,nov87,feb88;  2'93 merge chap/snglpoly, set curve(fm)
c-OUTPUTS-----------------------------------------------------------------------
c  Error checking, debug and output statements are in lower case.
c               Normal printed outputs from Polan are indicated by    ------->
c               Trace outputs and abnormal conditions are shown by    **----->
c               Fuller debug outputs, obtained with  list > 0,  are   ##----->
c  Loops are delimited by  c.....
c-------------------------------------------------------------------------------
 !     integer, value :: N
      real, intent(inout) :: FV(*),HT(*)

      real, intent(out) ::  QQ(*)

      integer numq, mode

      DIMENSION  IT(20),IV(20),IR(20),IH(20)
      COMMON /POL/ B(99,20),Q(20), FH,ADIP, MODE,MOD, FA,HA, tcont,lbug
      COMMON /POL/ HS, FC,FCC, HMAX,SH, PARHT, HVAL,VWIDTH,VDEPTH, XWAT
      COMMON /POL/ MAXB,NF, NR,NL, NX, MS,MT,JM, LK, KR,KRM, KV,MF, NC  counters
      common /test/ test1, test2, test3
      integer :: ndim = 0
cBottom set below is original
c             |------ first step -------|  |---- following steps -----|
c  At MODE =  1 2  3  4  5  6  7  8  9 10  1  2  3  4  5  6  7  8  9 10
      DATA IT/1,2, 3, 4, 4, 5, 6, 6, 6,73, 1, 2, 3, 4, 5, 6, 6, 6, 6,73/  nterms
      DATA IV/1,2, 3, 4, 5, 7, 8,10,12,90, 1, 1, 2, 3, 4, 5, 7, 8,13,90/  nvirts
      DATA IR/0,0, 0, 1, 1, 2, 2, 3, 5, 2, 0,-1,-1, 1,-2,-3,-3,-4,-6,-3/ realhts
      DATA IH/1,1, 2, 3, 3, 4, 5, 6, 8,28, 1, 1, 1, 1, 1, 1, 2, 2, 3,28/ calchts
c
         tcontf(x) = (fa*(fa/j+2.*x/(j+1)) + x*x/(j+2)) *x**j           polycont
         maxb= 98          !=dimension of array B, -1.  Set IV(10)=IV(20)=maxb-9
         if (ndim.eq.0) then                                            1st call
            ndim= max(n, 100)
            if (n.lt.100) print*,'SET INITIAL N = NDIM, > #DATA + 46'
            endif
c                       Set number of coefs.  Don't store any if qq(1) = -1.
                        if (qq(1).ne.-1.) QQ(1) = 1.
          FV(NDIM-45) = -1.                                             end data
          HT(NDIM-45) =  0.
c1 ---                       (1) Initialisation
c
      if(list.gt.1) write(2,10)
     &                    n,fv(1),ht(1),fb,dip,start,amode,valley,list
10           format (/'#Args: N,f1,h1=',I4,F5.2,F6.1,'   fb,dip,start',
     1                 F6.2,2F6.1,'   mode,valy,list',F5.1,F6.2, I3 /)
c1 ---
c1.1                        Set trace, gyrofreq  and  mode  constants.
      ADIP= MAX( ABS(DIP),.001)                                         fordip=0
      FH  = GIND(FB,-ADIP)                                              groundfh
         lbug= list
         if (lbug.gt.7) lbug = -lbug
         if (lbug.gt.1) call trace (fv,ht, 0.)                          ##----->
      MODE = int(ABS(AMODE))
         if (mode.eq.0)  mode = 6 + int(adip/60.)*10
      MODB = 0                                                          normal
         if (mode.ge.30)  modb = mode/10                                one poly
         moda = mode - modb*10                                          1stlayer
      MODC = 0
         if (modb.ge.30)  modc = modb/10                                toplayer
         modb = modb - modc*10                                          midlayer
         if (moda.eq.0)  moda = modb
         if (modb.eq.0)  modb = moda
         if (modc.eq.0)  modc = modb
      mod1 = moda - (moda/10)*10
      if (mode.eq.10.OR.mode.ge.30)  mod1 = 10                          snglpoly
c-------------------------------------------------------------------------------
c1.2--  Move the virtual height data to start at  fv(31), ht(31).
c    The real, virtual height origins are at  kr, kv  =  1, 31-js
c                          (where js is the number of freqs added below fmin).
c Set fa,ha = the origin (starting point) of the first real-height polynomial;
c Set    lk = 1 / 0 / -1   for   o-ray / poly / slab start.
c       Add an interpolated point at fv(30),ht(30) to control o-ray starts.
c
      CALL  SETUP (FV, HT, NDIM, START)
c
      DHS = MAX(HA-HS, 0.)                                              for c3.3
c===============================================================================
c////////////   R E A L   H E I G H T   A N A L Y S I S   L O O P   \\\\\\\\\\\\
c
c Calculate one polynomial, with NT terms, from the point fa=fv(k), ha=ht(k)
c   to fit the next NV virtual and NR real hts.
c                     (Include one real height below ha if NR is negative).
c   Calculate a further NH real heights,  and set  k = k+NH.
c
c   If a critical frequency is found in k+1 to k+NV+1, calculate up to
c         the preceding freq, and determine a least-squares Chapman peak.
c With x-ray data (-ve freqs) at start or after critical, recalculate ha.
c
c Real height origin (fa,ha) is at k = kr, virtual at k = kv.  krm = top real.
c-----------------------------------------------------------------------------
c2 ---                             (2) Prepare  Data
c
c2.1                                Start of analysis (or restart after a peak)
200   MOD = MOD1                                                        =1 to 10
      KRM = KR                                                          top real
      NNR = 0
      GO TO 230
c2.1a                                   Later segments (using mod+10)
210   IF (MOD.LE.10)  then
         KR = KR-IR(MOD)                                                step bak
         KV = KV-IR(MOD)
         MOD = MOD+10
         endif
      NNR = IR(MOD)                                                     real hts
c2.1b                                   All segments
230   NT = IT(MOD)                                                         terms
      NV = IV(MOD)                                                      virtuals
      NH = IH(MOD)                                                       to calc
      NR = IABS(NNR)
      NL = MIN(1,NR-NNR)                                                 prev ht
      MS = 0
         MSA = 0
         IF (LK.EQ.1.and.KV.EQ.29.and.AMODE.GE.0.) MSA = 1              GradEqun
         IF (NV.GT.1.and.NV.LT.8)  NV = NV - MSA                        at start
c-------------------------------------------------------------------------------
c2.2                   Count initial x-rays. Check frequency sequencing.
c                      Check for cusp, peak, or end of data.
c Set NF= no of o-rays (=NV, if sufficient points exist before a peak/restart).
c     NX= x-rays;   MV= NF+NX;   FM= fv(mf) = top frequency in this step.
c     FCC = FC or 0.1 for a peak,  = -.1 for cusp at fm,  = 0. otherwise.
c
      CALL  SELDAT  (NV, FV, HT)
                                  tras= -2.2
                                  if (nv.lt.0)  go to 730               exit>>>>
      FM = FV(MF)                                                       top freq
      MV = MF-KV                                                        #  freqs
      IF (NF.GT.NV)  NH = NH + (NF-NV)/2                                to calc.
            if (fcc.gt.0..and.lbug.eq.-8) lbug = 8
         if(lbug.gt.2.OR.(lbug.eq.2.AND.mod.lt.10))calltrace(fv,ht,2.2) ##----->
      IF (KR.EQ.1)  GO TO 300                                              start
c-------------------------------------------------------------------------------
c2.3    Subtract the group retardation due to the last calculated real-height
c  section.  This modifies all h' at f > fa (where fa = fv(KR)),  and increases
c  the index LK (giving the height to which group retardation is removed) to KR.
c                                  Will set KV= KV+1, if an h' error is removed.
      CALL  REDUCE (FV, HT)
      MV = MF - KV                                    !in case kv was increased
                                  tras = -2.3
                                  if (jm.le.0) go to 740                exit>>>>
                                  if (lbug.gt.2) call trace (fv,ht,2.3) ##----->
c===============================================================================
c3 ---                     (3) Set Up Equations for next profile step
c
c3.1                                   Set number of terms, and origin.
300   MT = NT
cc 2'93 to below      IF (NT.GT.20)  MT = (NF+2) *NT/100                !mode 10
      IF (NT.GT.20)  MT = NINT(2.*SQRT(REAL(NF)))                        mode 10
      IF (LK.LT.0.OR.HVAL.NE.0.) MT = MT + 1                            physequn
      IF (NX.GT.0)   MT = MT + NX/3                                     + x-rays
      IF (NF.GT.NV)  MT = MT + (NF-NV)/3
      IF (MODE.GE.30)  then
         MT = MODA                                                      1stlayer
         IF (KR.GT.1)  MT = MODB                                        2ndlayer
         IF (FV(MF+2).EQ.0. .OR. NF.GE.NV)  MT = MODC                   last
         endif
      MT = MIN(14, MT, MV+NR+MSA)                                       maxterms
      JM = MT
      FA = FV(KR)                                                       newstart
      HA = HT(KR)                                                       newstart
      NC = 0                                                            iteratns
c                                      Check for valley, set valley flag hval
c                                      Initialise valley width & depth.
c                                   Increase MT by 0-2, for phys equns (c3.3)
      LOOP = 0
      CALL STAVAL (FV,HT, DIP,VALLEY, LOOP)
      if (hval.ne.0..and.amode.ge.0.)  MT = MT + MAX( (9-MT)/3, 0 )
c                                      Start or valley require q(mt+1), = offset
      IF (LK.LE.0.or.HVAL.NE.0.)  JM = MT + 1
                                  if (lbug.gt.3) call trace (fv,ht,3.1) ##----->
c-------------------------------------------------------------------------------
c3.2                                   Set Up Equations in B
c                                             (start, valley iterate here)
330   CALL COEFIC (MV,FV,HT)                                            setcoefs
                   tras = -3.2
                   if (jm.le.0)  go to 740                              exit>>>>
      IF (AMODE.LT.0..or.MV+4.GT.MAXB)  GO TO 380                       no phys
      IF (MSA.EQ.1) then
            G = 1.0 + 1.8/FV(30)                                        o start,
            B(MV+1, 1)  = .5                                              set Q1
            B(MV+1,JM+1)= .5*G*(HT(30)-HA)                              LK=1only
            MS = 1
            endif
c3.3                                    ADD PHYSICAL RELATIONS
      IF (LK.LT.0)  then               !Add ms physical start relations; x only
            WS = 0.1                                                      weight
            IF (START.GT.44.)  WS = 0.3
         B(MV+1,JM)   = - WS                                            -q(jm)
         B(MV+1,MT)   = .25 *WS                                         +q(mt)/4
         B(MV+1,JM+1) = DHS *.25/FA**2 *WS                                =ha-hs
         B(MV+2,MT)   = WS*.3                                           q(mt)
         B(MV+2,JM+1) =(HS/2.-40.) /FA**2 *WS*.3                           =slab
         B(MV+3,MT-1) = WS*.15                                          loworder
         MS = 3
      else IF (HVAL.NE.0.) then        !Add ms physical valley relations; x or o
            WV = 1.0                                                      weight
            IF (HVAL.LT.-2.)  WV = 10.                                  fix valy
            IF (NX.GT.0)  WV = 0.2
         B(MV+1,JM)   = WV
         B(MV+1,JM+1) = (VWIDTH - PARHT) *WV                            stddvaly
         B(MV+2,MT) = 0.5                                               loworder
         B(MV+3,1)  = 0.4                                                grad-
         B(MV+3,JM) = -.1/VDEPTH                                          -cont.
         B(MV+4,MT-1) = 0.15                                            loworder
         MS = MIN(4, MT)
      endif
c
380   NC = NC+1
          if(lbug.gt.2.OR.(ms.gt.0.AND.lbug.ne.0)) call trace(fv,ht,3.3)##----->
c===============================================================================
c4 ---                            (4) Real-Height Solution
c4.1                                  Solve equations in B
      NS = MIN(MV+NR+MS, MAXB)
      CALL  SOLVE (NS, JM, B, Q,  0., devn,lbug)                         calc. q
                                 tras = 4.1
                                 if (ns.lt.jm)  go to 720               error>>>
c
c4.2                                  Check & adjust physical constraints
c                                                         (incl peak curve 2'93)
                          if (msa.eq. 1) lbug= lbug+200                 leave q1
                          if (dip.lt.0.) lbug= lbug+400                 noadjust
      CALL ADJUST (HA,FA,FM,FCC, LK,JM,MT, B,Q, lbug)
c
c4.3                                  Iterate a valley or x-start calcn.
      LOOP = 2
      IF (JM.NE.MT.and.NX.LE.0)  then
         VWIDTH = Q(JM) + PARHT         !O-ray Valley; loop once to adjust depth
         IF (NC.EQ.1.AND.HVAL.NE.0.)  LOOP = 1                          adj.vdep
         endif
         if(lbug.gt.1.AND.hval.ne.0..OR.lbug.gt.2)calltrace(fv,ht,4.3)  ##----->
c
      IF (FB.LT.0.)  LOOP = -IABS(LOOP)                                 fixed fh
      IF (JM.GT.MT)  CALL STAVAL (FV,HT, START,VALLEY, LOOP)   !strt/val; adj HA
      IF (LOOP.EQ.4)  GO TO 330                                         valyloop
      IF (LOOP.EQ.3.and.FB.GT.0.)  GO TO 330                            xrayloop
c===============================================================================
c5 ---                             (5) Store Real Heights
c
      KRM= KR +NR -NL                                                   lastreal
      KVM= KV +NR -NL
      FL = FV(KVM)                                                      fortcont
      NH = MIN(NH,NF-MOD1/10*2)                                         to calc.
      MQ = MT + MIN(LK,0)                                               polterms
      KA = KVM + 1                     !(peak itern loops here)         next pt.
c.....
      DO 520  K = KA, MF
         FR = FV(K)
         KRM = KRM+1
         DELTF = FR-FA
         FV(KRM) = FR
520      HT(KRM) = HA + SUMVAL(MQ,Q,DELTF,1)
c.....
      FN = FV(KVM+NH)
      IF (FCC.NE.0.)  then
         NH = KRM-KR                                                    to peak.
         FN = FM
         endif
      DO 560  J = 1, MQ
560      tcont = tcont+j*q(j)*(tcontf(fn-fa) -tcontf(fl-fa))            polycont
c
      KR = KR+NH                                                         step on
      KV = KV+NH                                                          origin
      IF (FCC.EQ. 0.)  GO TO 210                                         do next
      IF (FCC.EQ.-.1)  GO TO 200                                            cusp
      FC = FCC
c///////////////////////////////////// End  of  Normal  Steps //////////////////
c6.0 ---                          store coefficients of last polynomial in qq
      if (qq(1).lt.1.)  go to 620                                        no save
c
            numq = int(qq(1)) + 4
            do 610 j = numq-3, numq+mq
               if (qq(j).eq.-1.)  go to 620                            !no store
               if (j.ge.numq) qq(j)= q(j-numq+1)      !store poly coefs q, in qq
610            continue
            qq(numq-3) = fa
            qq(numq-2) = ha
            qq(numq-1) = mq
            numq = numq + mq
            qq(numq)= devn
            if (fcc.lt.0.)  numq= numq + 1
            if (fcc.lt.0..and.qq(numq).ne.-1.)  qq(numq)= fm
            qq(1) = numq
620   continue
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ (6) Calculate and List Peak data $$$$$$$$$$$$
c6 ---                    Least-squares fitting of a Chapman-layer peak.
c
            tras= 6.
            if (fcc.lt.0.)  goto 760                                      end>>>
            if (ht(krm).lt.70.)  go to 740                              abort>>>
      CALL  PEAK  (FV, HT,  QQ)
c-------------------------------------------------------------------------------
c7 ---                             (7)  Restart for a New Layer,  or end.
c7.1
             if (lbug.ne.-10) lbug = iabs(lbug)                         trace on
             if (lbug.eq. 8)  lbug = -8
      HS = (HMAX + HT(KRM))/2.                                            for fh
      IF (HT(KV+1).GT.0.)  GO TO 200                                    newlayer
c
c7.2                All finished.  Add points at z= 0.35, 0.85, 1.5  above peak,
c                                          using a scale height gradient of 0.1.

      N = KR + 3
      FV(N-2) = FC *.98642                   !0.3=.98984, 0.4=.98257, 0.5=.97372
      FV(N-1) = FC *.93300                   !0.8=.93957, 0.9=.92622, 1.0=.91213
      FV(N  ) = FC *.83462                   !z = 1.5
      HT(N-2) = HMAX + SH *0.3565            !z=.30+0044,.35+0065,.4+0081,5+0128
      HT(N-1) = HMAX + SH *0.8878            !z=.8+0332, .85+0378,.9+0424,1+0526
      HT(N  ) = HMAX + SH *1.6216            !z= 1.5 +.1216
      FV(N+2) = tcont*.124                                                T.E.C.
      HT(N+2) = SH                                                      scale ht
      FV(N+3) = 0.                              !( f,h(n+1)=devf,devh from PEAK)
      HT(N+3) = 0.
      NUMQ = int(QQ(1)) + 1
      if (numq.gt.1)  QQ(NUMQ) = -99.
      RETURN
c-------------------------------------------------------------------------------
c7.3                                   E x i t   o n   E r r o r
720           if (fc.eq.-1.)  go to 760
730           if (nv.lt.0)  kv = -nv -1                                 >>seldat
740           if (lbug.eq.-10)  go to 760                               no lists
                     write(2,750) (fv(i),ht(i), i= kv-1,kv+4), tras, kv **----->
750                        format ('>>>>>error at f,h ='12f7.2,
     $                             ' section',f4.1,i5,'  >>>> end', /)
c
                                    call trace (fv,ht,tras)             **----->
760           n = kr + 1
              fv(n) = fv(kv+1)
              ht(n) = ht(kv+1)
              fv(n+1) = 0.
              ht(n+1) = 0.
                 lq = int(qq(1)) + 1
                 if (lq.gt.1.and.qq(lq).ne.-1.)  qq(lq) = -98.          err flag
              return
      END Subroutine Polan

