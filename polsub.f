c   POLSUB.FOR  =  COEFIC,  ADJUST,  REDUCE.
      SUBROUTINE   COEFIC  (MV, FV, HT)
c                                                                jan77/may86
c   Calculate coefficients  b(i,j),  for the real-height polynomial
c            h-ha = q(j)*(f-fa)**j;     MT terms,  fa=fv(k),  ha=ht(k).
c   The first mv rows of b give virtual ht coefs at freqs k+1 to k+mv;
c       following nr rows give real height coefs at freqs k+1 to k+nr.
c   Adds terms for a linear section at x ray start or above peak.
c   Subtracts a parabolic segment (with s.h.=1.40*sh) above peak.
c
c   Called with MV -ve (from "reduce", in section 2.3) to reduce the next
c   |MV| virtual heights by the delay in the previous section fa to fv(k).
c-------------------------------------------------------------------------------
      COMMON /POL/ B(99,20),Q(20), FH,ADIP, MODE,MOD, FA,HA, TCONT,LBUG
      COMMON /POL/ HS, FC,FCC, HMAX,SH, PARHT, HVAL,VWIDTH,VDEPTH, XWAT
      COMMON /POL/ MAXB,NF, NR,NL, NX, MS,MT,JM, LK, KR,KRM, KV,MF, NC  counters
c
      DIMENSION TR(17),W(17), GAUSS(17),FNR(17), FV(*),HT(*)
      DATA TR  /  .046910077,.23076534, .50     ,.76923466,.95308992,
     1 .009219683,.047941372,.11504866,.20634102,.31608425,.43738330,
     2  .56261670,.68391575,.79365898, .88495134,.95205863,.99078032 /
      DATA  W  /  .11846344,.23931434, .28444444,.23931434,.11846344,
     a  .02358767,.05346966,.08003917, .10158371,.11674627,.12457352,
     b  .12457352,.11674627,.10158371, .08003917,.05346966,.02358767 /
      DATA  VB, FTC, GTC  /  0.6,  0., 180. /
         sq(x)  =  sqrt((1.-x)*(1.+x))
    1 format (a,15f8.2)
c-------------------------------------------------------------------------------
cc1 ---
cc1 ---                          (1).  Initialise
      F1 = FV(KV+1)
      NFF = IABS(MV)                                                      freqs.
c                                      end points for tapered weights
      KM = KV + NFF
      FW = FV(KM+1)                                                       for wv

      IF (FW.LT.1.)  FW = FV(KM)*2. - FV(KM-1)                           at fw=0
      I = KV + MIN0(NFF, MAX0(MOD-15,2)+NR)
      DW = FW - FV(I)                                                   for wv=1
      DA = FV(KV+NR-NL) - FA                                            for wv=1
c
c                                      gyrofrequency height
      HR = HA + 10.*AMIN0(1,KR-1)
      IF (KRM.GT.KR)  HR = HT(KR+1)
      IF (KRM.EQ.KR)  FTC = FA*1.2
      IF (HVAL.NE.0.)  HS = HMAX + VWIDTH/3.
      IF (MV.LT.0)  GO TO 10
      IF (LK.GE.0.OR.KRM.EQ.KR)  GO TO 6
c                                      fh level, for slab start
      DELTF = 0.
      GRAD2 = Q(1) - GTC
      IF (GRAD2.ge.0.)  then                                            abnormal
2           grad1 = grad2
            deltf = deltf + .1*fa
            grad2 = sumval(mt-1,q,deltf,2) -gtc
              if (grad2.gt.0..and.deltf.lt.fa*.45)  go to 2
            deltf = deltf + .1*fa* amin1(grad2,0.) /(grad1-grad2)
         endif
      FTC = (FTC + FA + DELTF)*.5                                        average

c                                      Clear matrix B
6         DO 8  I = 1, MAXB
          DO 8  J = 1, 19
8         B(I,J) = 0.
10    CONTINUE
c===============================================================================
c###########################  Loop frequencies: i= 1 to nff normally,
c                                                = 0 to nff to fit 1 back height
      I = 1 - NL
      IF (I.EQ.0.AND.MOD.NE.12)  I = -1                                  prev ht
      IF (MV.LT.0)  I = 1                                               reducing
c==============================================  Loop frequencies:
20    IREAL = NFF + NL + MAX0(I,0)
      KVI = KV + I
      KRI = KR + I
cc2 ---
cc2 ---                          (2).  Frequency and weight
      F  = FV(KVI)
      HV = HT(KVI)
      FR2 = F*F
         IF (KRI.LE.KRM)  HR = HT(KRI)*.8 + HT(KR+1)*.2                 refln ht
         HFH = HR
         IF (F1.LT.0.)  HFH = HS
      FH = GIND(0.,HFH)                                                 set fh
         IF (F.LT.0.)  FR2 = F*(F+FH)                                    x ray
      FR = SQRT(FR2)
c                                      Weighting for least squares analysis
      WREAL = 10.
      IF (KVI.LT.KV)  WREAL = 3.                                        back ht.
          IF (I.LE.0)   GO TO 60                                        realonly
                    if (fr.le.fv(kr).or.fr.le.fc) jm = -jm
                    if (parht.gt.2.*sh.or.nff+nr.gt.maxb)  jm = -jm     fatal>>>
                    if(jm.le.0)write(2,*)'coef2:',i,kri,kvi,nr,mf,mv,        >>>
     &                                 hv,f,fr,fc,fw,fv(kr),fv(km),parht     >>>
                                       if (jm.le.0)  return               end>>>
      IF (MV.LT.0)  GO TO 30
      WVIRT = 1.
      IF (F.LT.0..or.FW.LE.FR)  then
            WVIRT = XWAT
         ELSE IF (NR.GT.NL.and.FCC.EQ.0.) then
            WVIRT = SQRT( AMIN1((FR-FA)/DA, (FW-FR)/DW) )               !taper w
         END IF
      B(I,18) = F
      B(I,19) = WVIRT
c                                      vary weight with freq interval (o ray)
      FVL = FV(KVI-1)
      IF (nx.gt.0.or.F.LT.0..or.FVL.LT.0.)  GO TO 30
              DELTF = (FV(KVI+1) - FVL) /2.
              IF (DELTF.LT.0.)  DELTF = F - FVL
              WVIRT = WVIRT *SQRT( 8.*ABS(DELTF)/SQRT(F) )
ccc              wvirt = wvirt * 150./(ht(kvi) - 50.)  !!!!2'93 lower h'->better fit
c-------------------------------------------------------------------------------
c                      Integrate from TB (at f=fa) to TA (=0 at reflection)
c          17-point integrals are 5-pt to mid point TC, plus 12-pt from TC to TA
cc3 ---
cc3 ---                          (3).  Set integration range
30    TB = SQ(FA/FR)                                                    lowerlim
      TA = 0.                                                           reflectn
      IF (MV.lt.0) then             !REDUCTION STEP, set upper limit FP
         FP = FV(KR)                       !2'93 added point below peak, so:
         IF (FP.EQ.FC) FP = FV(KR-2)       !(kr-1 ->-2 with added ht)   end poly
          if(fp.ge.fr) write(2,1)'coef3',fc,fp,fr                       !##--->>
         TA = SQ(FP/FR)                                                 upperlim
         inum = 5
         if (adip.gt.73.and.ta.lt..3)  inum = 12
        else
         FTC = AMIN1(FTC, FA*.2 + FR*.8)                                for xray
         TC = 0.39-.05/COS(.016*ADIP)                                   for 17pt
c                                      Set integration order inum
         inum = 5
         if (adip.ge.60. .or. mode.ge.8)   inum = 12
         if (adip.ge.80..and.tb.gt.1.2*TC)  inum = 17
               if (f1.lt.0.) inum = 17                                  x calcn.
               if (f1.lt.0.) TC = SQ(FTC/FR)
        endif
c
      ira = 1
      if (inum.eq.12)  ira = 6
      irb = ira + inum-1
      if (inum.eq.17)  TA  = TC                                         17-pt(a)
38    TD = TB - TA                                                      range
cc4 ---
cc4 ---                          (4).  Retardation in start/peak section
      DEPAR = 0.                                                        peak sum
      DELIN = 0.                                                        slab sum
      IF (LK.GE.0.AND.PARHT.EQ.0.)  GO TO 48
         IF (SH.NE.0.)  DZ = .5*PARHT/SH                                 z at fm
            DO 45 IR = 6, 17
            if (lk.ge.0)  then                                          lin only
               TP = SQ(FC/FR*SQ(TR(IR)*DZ))
               DEPAR = DEPAR + GIND(F,TP) *W(IR)                        pardelay
               endif
            GINDV = GIND(F,SQ((FA-TR(IR)*VDEPTH)/FR))
   45       DELIN = DELIN + GINDV*W(IR)                                 lindelay
                 IF(LBUG.EQ.7)WRITE(2,1)'COEF4',FC,F, DEPAR,DELIN,GINDV ##===>>>
             if(lbug.eq.9)write(2,1)'fh,fr,ftc,hfh,hr=',fh,fr,ftc,hfh,hr
            if (parht.gt.0.)  fh = gind (0., (hmax+vwidth)*.8 + hr*.2)
   48 CONTINUE
cc5 ---
cc5 ---                          (5).  Group index for polynomial terms
      DO 58 IR = IRA, IRB
         IF (inum.ne.17.or.IR.NE.6)  GO TO 55
            TD = TA                                                     17-pt(b)
            TA = 0.
            FH = GIND(0.,HR)                                            refln fh
            IF (F.GT.0.)  GO TO 55                                      not xray
               FR2 = F*(F+FH)
               FR = SQRT(FR2)                                           refln fr
               TD = SQ(FTC/FR)                                          corr. tc
55       T = TR(IR)*TD+TA
         FN = SQ(T)*FR
         GAUSS(IR) = GIND(F,T)*T*W(IR)/FN *TD*FR2
58       FNR(IR) = FN - FA
c-------------------------------------------------------------------------------
cc6 ---                          (6).  Store coefficients in array B
c
60    DO 90 J = 1, JM
         RH = (FR-FA)**J                                                real ht.
         DPEAK = 0.
         SUM = 0.
            IF (KVI.EQ.KV.AND. J.EQ.1)  RH = 1.                           mode 2
            IF (I.LE.0)  GO TO 88                                       realonly
            IF (J.EQ.JM)  DPEAK = DEPAR*PARHT                           par.peak
         if (J.ge.MT.and.(J.ne.MT.or.LK.lt.0))  then
c                                                  Start or valley terms
               RH = J-MT
               IF (J.EQ.MT)  SUM = DELIN                                lin.slab
               IF (LK.GT.0)  SUM = GINDV*VB+DELIN*(1.-VB)               lin.valy
          else
c                                                  Polynomial terms
            DO 70 IR = IRA,IRB
               A = GAUSS(IR)
               GAUSS(IR) = A*FNR(IR)
70             SUM = SUM+A
            SUM = J*SUM                                                    delay
          endif
c                                                  Coefficients B(i,j).
80       B(I,J) = (SUM+RH)*WVIRT                                         virtual
         IF (MV.LT.0) HT(KVI) = HT(KVI)-Q(J)*SIGN(SUM,HV)-SIGN(DPEAK,HV)  reduce
88       IF (IREAL.LE.NFF+NR)  B(IREAL,J) = RH*WREAL                        real
90       CONTINUE
c
cd      IF(LBUG.EQ.5)WRITE(2,*)'COEF6', NFF,I,IREAL, F,FR,HV,dpeak,WVIRT ##-->>>
      IF (MV.LT.0)  GO TO 100
c-----------------------------
cc7 ---                          (7).  Store Virt/Real heights as r.h.s.
      IF (I.GT.0)  B(I,JM+1) = (HV-HA-DPEAK)*WVIRT                      V r.h.s.
c
      IF (IREAL.GT.NFF+NR)  GO TO 100
         REAL = HT(KRI) - HA
            IF(MOD.EQ.12) REAL = Q(1)+2.*Q(2)*(FA-FV(KR-1))             gradient
         B(IREAL,JM+1) = REAL*WREAL                                     R r.h.s.
100   I = MAX0(I,0)+1
      IF (I.LE.NFF)  GO TO 20                                           loopfreq
c############################# End of frequency cycle ==========================
      RETURN
      END
c===============================================================================
      SUBROUTINE  ADJUST (HA,FA, FM,FC, LK,JM,MT, BB, Q, LBUG)
c                                                                 jan82/jun83.
c Called from section c4.2 of Polan, after the initial call to Solve,
c Adjust is used to condition the real-height polynomial  (in  df = f-fa)
c        h = ha + df*q(1) + df**2 *q(2) + ... + df**mq *q(mq)
c     by ensuring that
c 1-  The gradient  dh/df  is positive at  f = fa.
c     For an O-ray start calculation, an initial gradient  Q1  of
c                        greater than 80 km/MHz is reduced to  160.Q1/(80+Q1).
c     For an X-ray slab start (lk =-1) the gradient is limited above 160 km.
c 2-  The solution is well-defined:
c               -  the last 3 coefficients do not alternate in sign, and
c               -  these coefficients have absolute values less than 999.
c     These conditions are omitted for the peak and x-ray start, when a large
c                          gradient may be required at one end of the fit range.
c 3-  For x start calculations:  the underlying slab is positive,
c                                the initial height offset is negative,
c                                the initial height (+offset) is > 60 km.
c 4-  For a valley calculation:  the initial gradient  dh/df > sha,
c                                the initial curvature  q(2) < -1.5,
c                                the initial height offset is positive.
c 5-  For the last polynomial in a layer (with a peak to be fitted):
c                       the curvature d\FDh/df\FD is large at the top frequency fm.
c
c Adjustments alter the least-squares coefficients  q(1)  to  q(jm),
c     and the rms deviation of the fit  devn  (stored in q(20)).
c-----------------------------------------------------------------------
c     When  jm = mt,  the polynomial ends at  mq = mt.
c-----
c     jm = mt+1  and  lk.le.0  indicates an x start calculation.
c        The poly then ends at  mq = mt-1,  and
c        q(mt) = underlying slab thickness,  preferably positive;
c        q(jm) = starting height offset,  which must be negative.
c-----
c     jm = mt+1  and  lk.gt.0  indicates a valley calculation.
c        The poly then ends at  mq = mt,  
c        and  q(jm) = the valley width,   which must be positive.
c-----------------------------------------------------------------------
c  Adjustments are made by Call Solve (-mm, -n, bb, q, qset, devn,lbug).
c  This modifies the least-squares solution by adding the constraint
c  q(n) = qset,  with a weight  mm/100.   Any constraints previously
c  applied at this or larger values of n are first removed by  Solve.
c-----------------------------------------------------------------------
      DIMENSION BB(99,20), Q(20)
c
            lflag= (lbug+50)/100                         !2=leave q1; 4=omit adj
            lbug = lbug- lflag*100
            ladj = 0
c                    Set mq = no of poly terms,   and head
      DEVN = Q(20)
      MQ = MT
      IF (LK.LT.0)  MQ = MT-1
         if (lbug.gt.2.or.(lbug.gt.0.and.lk.lt.2)) write(2,7)(m,m=2,mt) ------->
      if(lbug.gt.0)write(2,9)'--- ',lk,jm,mt,ha,fa,fm,(q(i),i=1,jm),devn------->
7           format ('*adjust----  lk jm mt    ha     fa    fm',6x,
     1              'q1',9i9/40x,10i9)
9           format ('*adjust', a5,3i3, f8.2,2f6.2, 10f9.2 /41x,10f9.2)
         if (lflag.gt.3.or.jm.lt.4) return                              noadjust
c-----------------------------------------------------------------------
c(a)                          Initial gradient;  must exceed 1.5
      QMIN = 2.
      if (lflag.gt.1) goto 100                      !initial O start; leave q(1)
      if (lk.le.1)  then   !====  O,X Start: limit gradient to 160,320.
         QMAX = 80.
         IF (LK.EQ.-1)  QMAX = 160.
         IF (Q(1).GT.QMAX)  QMIN = 2.*QMAX*Q(1) / (QMAX+Q(1))
        else               !====  Valley:  gradient > sha
         IF (JM.GT.MT)  QMIN = HA/4. - 20.
        endif
c
      IF (Q(1).lt.QMIN.or.QMIN.ge.85.)  then
         CALL SOLVE (-900, -1, BB, Q, QMIN, DEVN,lbug)                  q(1) adj
         ladj= ladj + 1
      if(lbug.gt.0)write(2,9)'q(1)',lk,jm,mt,ha,fa,fm,(q(i),i=1,jm),devn------->
        endif
c                    Valley:  initial curvature < -1.5
      IF (JM.gt.MT.and.LK.gt.0.and.Q(2).ge.-1.5) then
         CALL SOLVE (-900, -2, BB, Q, -1.5, DEVN,lbug)                  q(2) adj
         ladj= ladj + 2
      if(lbug.gt.0)write(2,9)'q(2)',lk,jm,mt,ha,fa,fm,(q(i),i=1,jm),devn------->
       endif
c-----------------------------------------------------------------------
c(b)                    Check for unnecessary high-order coefficients.
c                Reduce order if last 3 signs alternate, and sizes increasing,
c                                       or if any size exceeds 999.
100   QL = Q(MQ)
      IF (MQ.LT.5.OR.ABS(QL).LT.150. .OR.FC.GT.0..OR.LK.LE.0) GO TO 200     o.k.
          IF (QL/Q(MQ-2).GT.2. .AND. QL/Q(MQ-1).LT.-1.)  GO TO 150        fix it
          IF (ABS(QL).LT.999..AND.ABS(Q(MQ-1)).LT.999.)  GO TO 200          o.k.
c                       Add constraint  q(mq) = 0.,  and reduce mq.
150           CALL SOLVE (900, -MQ, BB, Q, 0., DEVN,lbug)
      if(lbug.gt.0)write(2,9)'mq ',lk,jm,mt,ha,fa,fm,(q(i),i=1,jm),devn ------->
              MQ = MQ-1
              GO TO 100
200   CONTINUE
c-----------------------------------------------------------------------
c(c)       S T A R T  and  V A L L E Y:   check sign, and limit offset term
c
      IF (JM.EQ.MT)  GO TO 500                                   !no offset term
      IF (LK.GT.0 )  GO TO 400                                            valley
c--------               Start:-  slab
      QM = Q(MT)
      QMT= 0.1
      IF (QM.lt.0.)  then
         CALL SOLVE (-900, -MT, BB, Q, QMT, DEVN,lbug)                  slab= .1
      if(lbug.gt.0)write(2,9)'Slab',lk,jm,mt,ha,fa,fm,(q(i),i=1,jm),devn------->
         IF(DEVN.GT.Q(19)*1.25)CALL SOLVE(0,-MT,BB,Q, QMT,DEVN,lbug)      cancel
        endif
c                       Start:-  offset
         OFST = -0.1
         IF (Q(JM).GT.0.)  GO TO 450
         IF (LK.GE.0 .OR. HA+Q(JM).GT.60.)  GO TO 500
            OFST = 60. - HA
            GO TO 450
c--------               Valley
400   IF (Q(JM).GE.0.)  GO TO 500
c                 Correct initial real-height offset
         OFST = 0.1
450      CALL  SOLVE (-900, -JM, BB,Q, OFST, DEVN,lbug)
            if(lbug.gt.0) write(2,9) 'Ofset',
     &                           lk,jm,mt, ha,fa,fm, (q(i),i=1,jm),devn ------->
c-----------------------------------------------------------------------
c(d)       P O L Y N O M I A L   to   P E A K:   check curvature at fm.   (2'93)
c
500   if (fc.le.0..or.mq.lt.4) goto 600                  !no peak, or cubic only
      hfm = ha + sumval (mq, q, fm-fa, 1)
      shm = max(.27*hfm - 21., 5.)                               !model sh at fm
      d = fc - fm
      if (d.lt..03) d= (0.1*sqrt(fm) +.03*fm)     !fm=3,7,11,16->.26,.47,.66,.88
      curve= 0.28*shm /(d*sqrt(fm*d))             !parab=0.35*..,  80% ->0.28*..
      if (sumval(mq,q,fm-fa,3).gt.curve) goto 600                          !O.K.
         devn = fm - fa
         CALL SOLVE (100,-300, BB,Q, curve, devn, lbug)       !set curv(fm), W=1
                     if(lbug.gt.0) write(2,9) 'Curvm',
     &                           lk,jm,mt, ha,fa,fm, (q(i),i=1,jm),devn ------->
         if (q(1).gt.qmin .or.mod(ladj,2).eq.1) goto 600
                  CALL SOLVE (-900, -1, BB, Q, QMIN, DEVN,lbug)         q(1) adj
      if(lbug.gt.0)write(2,9)'q(1)',lk,jm,mt,ha,fa,fm,(q(i),i=1,jm),devn------->
c-----------------------------------------------------------------------
600   RETURN
      END
c===============================================================================
      SUBROUTINE   REDUCE  (FV, HT)
c                                                                         feb84.
c   To reduce the remaining virtual heights by the group retardation in
c       the last calculated real-height section  (from  fa  to  fv(kr)).
c
c   Modifies all h' at  f > fv(kr),  by subtracting the group retardation.
c   Increases lk (the height to which group retardation is removed) to kr.
c       Use  LIST = 7  to show details of the reduction steps.
c-------------------------------------------------------------------------------
      DIMENSION  FV(*), HT(*)
      COMMON /POL/ B(99,20),Q(20), FH,ADIP, MODE,MOD, FA,HA, TCONT,LBUG
      COMMON /POL/ HS, FC,FCC, HMAX,SH, PARHT, HVAL,VWIDTH,VDEPTH, XWAT
      COMMON /POL/ MAXB,NF, NR,NL, NX, MS,MT,JM, LK, KR,KRM, KV,MF, NC  counters
   7  format (' reduce:',15f8.3)
c-------------------------------------------------------------------------------
c2.a--       Reduction using full polynomial expression (at f<fm+mode/25).
      K = MF
      FRED = FV(MF) + .04*FLOAT(MODE)
20      K = K+1                                                         nextvirt
        IF (ABS(HT(K)).GT.40..AND.FV(K).LE.FRED)  GO TO 20              reduce k
c---
        KV1 = KV+1                                                      minreduc
        KVM = MIN0(K -1, KV+MAXB-NR)                                    maxreduc
                  if(lbug.eq.7)write(2,7) (fv(j),j=kv1,kvm)             ##----->
                  if(lbug.eq.7)write(2,7) (ht(j),j=kv1,kvm)             ##----->
      CALL COEFIC (KV-KVM, FV, HT)                                      prevpoly
                  if(lbug.eq.7)write(2,7) (ht(j),j=kv1,kvm)             ##----->
                  if (jm.le.0)  return                                  exit>>>>
c----------                            Error check on reduced virt hts.
         do 40  k = kv1, kvm
            if (abs(ht(k)).lt.ha)  go to 50                              low ht.
40          continue
         GO TO 60
c                                      list bad point (h'<ha) and delete it
50                 if(lbug.ge.-9) write(2,52) (fv(j),ht(j), j=k-1,k+1)  **----->
52                   format ('*****reduce: data error at f, h =',
     $                                     6f8.3,  i5,' end >>>>')
55                fv(k)= fv(k-1)                                         move up
                  ht(k)= ht(k-1)
                  k= k-1
                  if (k.gt.kv-4)  go to 55
                  kv = kv1
c-------------------------------------------------------------------------------
c2.b--       Approx reduction at f>fred for delay in fv(lk) (=fa) to fv(kr)
60    MQ = MT + MIN0(LK,0)                                              polterms
      LK = MAX0(LK,1)+1
      FAR = FV(KR)                                                      top freq
      G2 = Q(1)
c.....                            For section  fv(l-1) to fv(l)
      DO 100 L = LK,KR                                                  loopslab
c                                 :-  Get mean density, thickness
        FL = FV(L-1)
        FN = FV(L)
        AVN = (FL**2+FN*(FL+FN))/3.                                     mean f*f
        HA = HT(L)
        DH = HA - HT(L-1)
c                                 :-  Correct for curvature
        DELTF = FN - FA
        G1 = G2
        IF (FN.NE.FC.AND.FN.GT.FA)  G2 = SUMVAL(MQ,Q,DELTF,2)             dh/df.
        IF (FN.EQ.FC .OR.FN.LT.FA)  AVN = AVN+FN*(FN-FL)/3.                 peak
        DH = DH+(G2-G1)*(FL+FN)*(FN-FL)**2/(12.*AVN)                    correctn
        FH = GIND(0.,HA)                                                scale fh
c                                 :-  Loop frequency to end
        K = KVM
c---
80      K = K+1                                                         loopfreq
          F = FV(K)
          IF (F.EQ.-1.)  GO TO 100                                      end freq
          FR = F
          HV = HT(K)
      IF (HV.EQ.0. .AND.HT(K+1).EQ.0.)  GO TO 100                       end freq
              if (abs(hv).le.40.)  go to 80                             skippeak
              if (f.lt.0..and.f.gt.-fh)  go to 90                       f error.
              if (f.lt.0.)  fr = sqrt((f+fh)*f)
              if (fr.lt.far)  go to 90                                  f error.
          TAV = SQRT(1.-AVN/FR**2)
          HT(K) = HV - GIND(F,TAV)*DH*SIGN(1.,HV)                         reduce
                  if(lbug.eq.7)write(2,7) fl,fn,tav,dh, f,hv,ht(k)      ##----->
        GO TO 80                                                        loopfreq
c---                              End freq loop.   Delete bad data and end
90                if(lbug.ge.-9)write(2,52) (fv(j),ht(j), j=k-1,k+1),k  *******>
                  fv(k) = -1.                                           end data
c---
100   CONTINUE
c.....
      LK = KR

      END Subroutine REDUCE
