c  POLMIS.FOR  =  PEAK,  TRACE,  SOLVE,  SUMVAL,  GIND. - less peak
      SUBROUTINE  PEAK  (FV,HT, QQ)
c*************************************************************jan78/dec84/nov87.
c  Get FC, SH  by least-squares fitting of a Chapman-layer expression 
c     to the profile gradients at the last NK calculated points.
c  Include any scaled  FC, FCX,   and iterate from a model scale height.
c
c  Called with KR, KV pointing at the last calculated real height before peak.
c  Exit with KR, KV at peak.      LIST = 8 adds detailed trace outputs.
c-------------------------------------------------------------------------------
      DIMENSION   FV(*), HT(*), QQ(*), PQ(20)
      COMMON /POL/ B(99,20),Q(20), FH,ADIP, MODE,MOD, FA,HA, tcont,lbug
      COMMON /POL/ HS, FC,FCC, HMAX,SH, PARHT, HVAL,VWIDTH,VDEPTH, XWAT
      COMMON /POL/ MAXB,NF, NR,NL, NX, MS,MT,JM, LK, KR,KRM, KV,MF, NC  counters
!          sq(x) = sqrt((1.-x)*(1.+x))
    9 format (a, 10f10.3/(15x,10f10.3))
c-------------------------------------------------------------------------------
c                                      Set frequencies
      KR = KR + 1                                                        to peak
      KV = KV + 1                                                        to peak
      FM = FV(KRM)                                                      top freq
      HN = HT(KRM)                                                      ht at fm
      HT(KRM) = HN
      NK = MIN(NF*9/(NF+9), KR/3) + 2                     !2'93 10->9.  f to fit
         if (fv(kv-nk).lt.fm*.7)  nk = nk-1               !2'93 .5-> .7  AND .85
         if (fv(kv-nk).lt.fm*.85.and.nk.gt.5) nk = nk-1
         if (mode.eq.10 .or. mode.ge.30)  nk = nk-1       !& -1 for single poly
      KA = KV-1 -NK
      KB = KRM - NK
      SH = MAX(.27*HN - 21., 5.)                                        model SH
      HMAX = HN + .3*SH                                                   for fh
      FH = GIND(0., HMAX)
      MK = NK
        IF (FCC.GT..1)  MK = MK + 1                                     incl. fc
        FCX = ABS(FV(KV+1))
        HCX = HT(KV+1)
        IF (ABS(HCX).le.30..and.FCX.ge.FM+.5*FH)  then                  x data:
           MK = MK + 1                                                  incl fcx
           FV(KV+1)= FV(KV)                                             swap fc
           FV(KV) = FCX                                                  and fcx
           HT(KV+1)= HT(KV)
           HT(KV) = HCX
           KV = KV + 1
         endif
      fv(kr) = fcc
      ht(kr) = ht(kv)
        fnext= fv(kv+1)
        if (fnext.eq.-1.) fnext = 0.
        if (fnext.lt.0.)  then
           fnext = sqrt(fnext*(fnext+.99*fh))                           1st xray
              do 600  i = 2, 60
              if (fv(kv+i).ge.0.)  go to 610                            1st oray
600        continue
610        fnext = min(fnext, fv(kv+i))
        endif
c-------------------------------------------------------------------------------
c.....                                 Set up and solve peak equations
      DO 400 LOOP = 1, 3                                                 iterate
         pgrad = 1.0
c...
      DO 300 I = 1, MK
         F = FV(KA+I)
            IF (F.EQ.FCX)  F = SQRT(FCX*(FCX - FH))                      x crit.
         W = 2.**(14.*(F-FM)/F)  / .7                                     weight
         GRAD = 0.                                                          peak
         IF (F.le.FM)  then                                              f crit.
            W = W  * .7
            MQ = MT + MIN(LK, 0)                                        polterms
            GRAD = 4./F/ MAX( SUMVAL(MQ,Q,F-FA,2), 3.)                  2dN/dh/N
            if (i.eq.nk) grad= grad - grad*.2*(1.-float(mt)/nf)**3      HN corrn
          endif
         T = SH*GRAD
         pgrad = pgrad * 20.*(max(grad,.002))
         B(I,1) = W                                                     *log(FC)
         B(I,5) = GRAD                                                  save 5,6
         B(I,6) = (ALOG(1.+T)-T)*.25                                    for errs
         B(I,2) = B(I,6) *W                                               * SH^2
         B(I,3) = ALOG(F) *W                                              r.h.s.
300   CONTINUE
c                    add condition  SH^2 = SHA^2
            w = 0.1*sqrt(pgrad) + .02
            b(mk+1,1) = 0.
            b(mk+1,2) = w
            b(mk+1,3) = w
C...
               if(lbug.eq.8) write(2,320) nk,mk,sh,w
               if(lbug.eq.8.and.loop.eq.1)
     &         write (2,322) (fv(kb+k),ht(kb+k),b(k,5),b(k,1), k=1,mk)  !#----->
320            format (' peak: nk,mk,sh,w=', 2i3, f6.1, f7.3 )
322            format (7x,'F, H, G, W=', 16f7.2 / (18x,16f7.2))
C
      CALL SOLVE (MK+1, 2, B, PQ, qset, devn,lbug)
C                                   Set  sh = avge (sh,sha)  if ill-defined
      GMAX = B(NK, 5)
         erfac= 400./(1+mk-nk)/(nk-.8) *devn
         qfac = .1/gmax - .1/b((nk+1)/2,5)
         if (qfac.lt.erfac)  then
            avsh = max(PQ(2), 0.) *.5 + .5                              use mean
            if (PQ(2).gt.1.)  avsh = PQ(2)/avsh                         avg 1/SH
            call  solve (10, -2, B, PQ, avsh, devn,lbug)                set avsh
            endif
c                                   If no scaled fc, limit calcd fc to 1.5fm
      fc = exp(PQ(1))
      if (mk.eq.nk.and.fc.gt.1.5*fm)
     $            call solve (400, -1, B, PQ, alog(1.5*fm), devn,lbug)  limit fc
      SH = SH*SQRT( ABS(PQ(2)) )
      PARHT= SH*MIN( ALOG(1.+SH*GMAX), 1.8)
      HMAX = HN + PARHT
         FH = GIND(0., HMAX)
           if(lbug.eq.8) write(2,9) 'qfac,erfac,prht', qfac,erfac,parht
           if (parht.gt.sh/loop)  go to 500                             no itern
           if (qfac .lt. erfac *2.*loop**2) go to 500                   no itern
400   CONTINUE
C.....
C---------------------       Fit peak to last nk real heights.
500     fc = exp(pq(1))
        fcn = 0.
        fcm = fnext*.8 + fm*.2 -.015                        !upper limit (nov87)
        if (fc.gt.fcm.and.fcm.gt.fm)  fcn = fcm
        if (fc.lt.fm*1.002)  fcn = fm*1.002         !lowerlim(8'87; ->.002 2'88)
        if (fcn.ne.0.)  call solve (300,-1, B,PQ, alog(fcn), devn,lbug)
      FC = EXP(PQ(1))
      SHR = PQ(2)
      SUMHW= 0.                 !get hmax as weighted mean of fits to each point
      SUMW = 0.                             !'85 fit hn,hl; '88-91 fit G's only.
      fw = fv(max(kb,1))
      DO 520 I = 1, NK          !2'93: use mostly zf(fits f/fc), & zg(fits grad)
         W = (fv(kb+i)-fw)**2
         r = 4.*alog((fv(kb+i)/fc))       !at f,=1+z-e^z for Chap, -z^2/2 parab
         zf = sqrt(-2.*r)                              !1st (parab) approx to zf
         do 505 ii = 1, 3                                      !iterate to Chap:
            rf = 1. + zf - exp(zf)                                 !(z +ve down)
            gf = rf - zf                                       !Chap dr/dz at zf
  505       zf = zf + (r-rf)/gf                                     !next approx
         zg = alog(1. + sh*b(i,5))                    !Chap z to match poly grad
         ZM = zf  !!!*.8 + zg*.2          !fit grad is no help   !fit mean
         hlx = hmx
         HMX = HT(KB+I) + SH*ZM                               !peak ht for this z
            if (lbug.eq.8) write(2,9) ' f,h, hm,w, zf,g=',
     &                     fv(kb+i),ht(kb+i), hmx,w, zf,zg              !#----->
cc         if(lbug.eq.8)write(2,9)' r,z1;r2,g2;zf,zg', r,z1, r2,g2, zf,zg
         SUMHW= SUMHW + HMX*W
520      SUMW = SUMW + W
      HMAX = SUMHW/SUMW                                     !final mean estimate
C                             Adjust last point towards chap layer (ht [& grad])
      if (nf.ge.3) then
         HADJ = HMAX - HMX                                  ! (adds to hmax err)
         HT(KRM) = HN + HADJ*.80                     !2'93: move 80% to fit hmax
         ht(krm-1)= ht(krm-1) + (hmax-hlx)*.40       !2'93: move 40% to fit hmax
           fchap= fc*exp( .25*(zm - sh*gmax) )               !freq to match grad
           fadj = fchap - fm                                 !(added to fc devn)
cc2'93 ??best no shift?  fv(krm) = fm + fadj *0.1    ! ??move sl toward same grad
          if(lbug.eq.8) write(2,9)' fc,fad,hmax,had=',fc,fadj,hmax,hadj !#----->
        endif
      PARHT = HMAX - HT(KRM)                       !parabolic sectn for 'reduce'
         ht(kr)= hmax - 0.5*parht                  !2'93 add point -way to peak
         z = 0.5*parht/sh
         fv(kr)= FC*exp(0.25*(1. + z - exp(z)))    !on Chapman layer.
         kr = kr + 1                               !(use fp at (kr-2) in coefic)
      FV(KR)= FC
      HT(KR)= HMAX
      DEVN = DEVN*2/(1+MK-NK)
      IF (FNEXT.EQ.0.)  GO TO 700
c-------------------------------------------------------------------------------
c Check overlap - this section can be omitted with manual scalings (F monotonic)
      fc = min(fc, fnext-.03)
      dfc= fv(kr) - fc                                                  REDUCE f
      fm = fm - dfc
      devn= min( devn, .5*(fnext-fm)/fc)
      if (dfc.le.0.01)  go to 700
c                                           Layer overlap;  reduce frequencies
             if(lbug.ge.-9)write(2,640) fv(kr), fc, fnext               **----->
640             format (' * FC  reduced from',f7.3,'   to',F7.3,
     $                  '   since next layer starts at',f7.3,' MHz.' /)
            df = dfc/fv(kr)**3
               do 660  k = 1, kr
660            fv(k) = fv(k) - df*fv(k)**3
            fa = fa - df*fa**3
c-------------------------------------------------------------------------------
c                                           List peak parameters
700   devf = FC*devn + abs(fadj)                                        freq err
      bavge= -b(2,6) - b(nk,6)
      devs = SH/shr *devn/bavge                                         sh err.
      devh = devn *parht/bavge + devs*sh*gmax + abs(hadj)               hmax err
      TCONT= TCONT + PARHT*(FC**2+.5*FM**2) /1.5                        parabola
      SSH = SH
      IF (LOOP.EQ.1)  SSH = -SH                                           -model
      if(lbug.ge.-9)write(2,750) fc,devf,hmax,devh,ssh,devs,tcont/fc**2 **----->
750    format('PEAK',F7.3,' (',F5.3,') MHz,  Height',F6.1,' (',F4.1,
     1     ') km.    ScaleHt',f5.1,' (',F4.1,') km  SlabT',F6.1,' km')
      FV(KR+4)= devf
      HT(KR+4)= devh
      FV(KV) = FC
C                                   TO SAVE PARAMETERS IN ARRAY  QQ
      if (qq(1).le.1.)  go to 800                                !poly not saved
         numq = INT(QQ(1))
         do 770 j = numq+1, numq+4
            if (qq(j).eq.-1.)  go to 800
770         continue
         wf = 0.1 + 0.9/nf                    !nf=1,2,3,5,9..->1,.55,.4,.28,.2,.
         fp = fv(krm)*wf + fv(krm-1)*(1.-wf)          !2'93 best poly/chap match
         fp = min(fp, 0.98*fc)
         QQ(numq+1) = max(fp, 0.90*fm)                 !but keep in 0.9fm-0.97fc
         QQ(numq+2) = FC
         QQ(numq+3) = HMAX
         QQ(numq+4) = SSH
         QQ(1) = numq + 4                                            != new numq
800   RETURN
      END
c===============================================================================
      SUBROUTINE  TRACE ( FV,HT, TRAS )
c************************************************************************jan'79.
c          Debug listing,  produced when  0 < lbug < 10.
c  (lbug = 1 gives a single trace line for each real-height step, in ADJUST.)
c  lbug > 3 adds intermediate listings of freq,height  (& Q,devn, in ADJUST).
c
      DIMENSION  FV(*), HT(*)
      COMMON /POL/ B(99,20),Q(20), FH,ADIP, MODE,MOD, FA,HA, tcont,lbug
      COMMON /POL/ HS, FC,FCC, HMAX,SH, PARHT, HVAL,VWIDTH,VDEPTH, XWAT
      COMMON /POL/ MAXB,NF, NR,NL, NX, MS,MT,JM, LK, KR,KRM, KV,MF, NC  counters
c-----------------------------------------------------------------------
      IF (IABS(LBUG).GT.9)  RETURN                                       no list
      IF (TRAS.NE.0.)  GO TO 10
              WRITE(2,9) 'FRQ', (FV(I), I=1,15)
              WRITE(2,9) 'HTS', (HT(I), I=1,15)
9             FORMAT (' #',A3, 16F8.2/ (5X,16F8.2))
              RETURN
c-----------------------------------------------------------------------
10    IF (TRAS.EQ.2.2) WRITE(2,15)
15    FORMAT('---------------------------------------------------------'
     1     / '#TRACE:  kr  lk jm mt    ha     fa   frm    fm',
     2       '   kv  nf  nr nl nx ms  mode mod     hs' ,
     3       '    fc   fcc    sh   parht   hval vwidth vdepth' )
c
      FRM = FV(KRM)
      FM = FV(MF)
      WRITE(2,25) TRAS, KR,LK, JM,MT, HA,FA,FRM,FM, KV,NF,NR, NL,NX,MS,
     $               MODE,MOD, HS, FC,FCC,SH, PARHT,HVAL,VWIDTH,VDEPTH
25    FORMAT('#=',F4.1,':',2I4, 2I3, F8.2,3F6.2, 3I4, 3I3,
     $                 I5,I4, 3X,4F6.2, 4F7.2)
c
      IF (LBUG.LT.4.OR.ABS(TRAS-3.3).NE.1.)  RETURN
c
c                           list data, at lbug > 3 and tras =2.3 OR 4.3
      KM = KV+NF+NX
      WRITE(2,9) 'FRQ', (FV(K), K=KV,KM)
      WRITE(2,9) 'VHT', (HT(K), K=KV,KM+3)
      WRITE(2,9) 'RHT', (HT(K), K=KR,KRM)
60    RETURN
      END
c===============================================================================
      SUBROUTINE  SOLVE (M, N, B, Q,  QSET,devn, lbug)
c******************************************************************jan77.
c Solve m simultaneous equns in n unknowns, given in the array  B(m,n+1).
c The result is returned in  Q(1)  to  Q(n).     qset,lbug added sept'86.
c The rms fit error is returned in q(20), and the previous value in q(19)
c-----------------------------------------------------------------------
c  Adjustments to a previously obtained solution may be made by:-
c         Call  SOLVE (mm, -n,  B, Q, Qset, devn, lbug).
c  This modifies the least-squares solution by removing any previous 
c  constraint on Q(n),  and adding the constraint
c         Q(n) = Qset,  with a weight  |mm|/100.
c  Any number of the coefficients  Q  can be modified in this way.
c      If mm is negative, any constraints previously applied 
c      at larger values of n are also removed by Solve.
c      (Use mm = 0 to just remove a constraint on q(n).)
c      Addition of constraints is numerically stable.
c      Removal of constraints is less accurate but still o.k.
c 2'93 To add profile constraint on ht/grad/curv at  n = (-)100/200/300,
c      Set:  devn= freq-FA;  & Value= Qset,  Weight= |mm|/100 as before.
c-----------------------------------------------------------------------
      DIMENSION  B(99,20), Q(20), WSETQ(17), VSETQ(17)
cc      double precision s, r, a, c, d   !!results worse??
      SAVE  NQ, NP, WSETQ, VSETQ              !Delete initial CF77 for fortran77
      DATA  NQ, NP, WSETQ, VSETQ  / 2*0, 34*0. /
   98 format (a,i5,9f8.3)
         IF (N.LT.0)  GO TO 3                                    !add constraint
      NQ = N
      NP = N+1
      if (lbug.ne.5.and.m.ge.n)  go to 99
c                        LIST MATRIX B,  at lbug = 5 or error in parameters.
88                    write(2,92) (j, j = 1, np)
                      do 90 i = 1, iabs(m)
90                       write(2,94) i, (b(i,j), j = 1, np)
92                    format ('>>Solve:  J  = ',I6,12I9/ 12x,12i9)
94                    format ('Matrix B row',I2, 13F9.4/ 14x,12f9.4)
99                    if (m.lt.n) return                                error
c-----------------------------------------------------------------------
c                                 Householder transformations
      K = 0
1     K = K+1
         WSETQ(K) = 0.                                 !zero constraints on q(n)
         S = 0.
         MK = M-K+1
            IF (MK.LE.0)  GO TO 7
         S = SUMVAL (MK, B(K,K),B(K,K), 4)                  !sums B(k to m, K)^2
            IF (K.EQ.NP)  GO TO 7
         A = B(K,K)
         D = SIGN(SQRT(S), A)
         B(K,K) = A + D
         C = A*D + S
               if(c.eq.0.) m= -m                            !terminates on error
               if(c.eq.0.) goto 88
            DO 2  J = K+1, NP
            S = SUMVAL (MK, B(K,K),B(K,J), 4) /C
               DO 2  I = K, M
2              B(I,J) = B(I,J) - B(I,K)*S
      B(K,K) = -D
      GO TO 1
c=====------------------------------------------------------------------
c----------       Adjust previous solution.  All entries at top(3), exit at bot.
3     ww  = IABS(M) /100.                                                !weight
      nad= iabs(n)/100                               !1-3 adds ht,grad,curv equn
               if (lbug.ge.4.or.(iabs(n).gt.nq.and.nad.eq.0))
     $         write(2,98) '>>>reSolve with n,qset,w=',n,qset,.01*m
         if (nad.gt.0) then                                 !set up equn, to add
            df = devn                                       !freq (-FA) for equn
            fj = 1.
            if (nad.gt.1) fj = 1./df
            nn = max(1, nad-1)                              !h,g,c from Q(1,1,2)
            do 35 j = nn, nq                                !store equatn in Q()
               fj = fj * df                                 !height = f, f^2,..
               qj = fj                                      !grad,curv=1,f,f^2,.
               if (nad.gt.1) qj = qj * j                    !grad: sum(j.Qj)
               if (nad.gt.2) qj = qj * (j - 1)              !curv: sum(j.j-1.Qj)
   35          q(j) = qj * ww
               q(np)= qset * ww                             !r.h.s.
            remov = 1.
            goto 52                                         !add in new equation
            endif
c                       m -ve; remove any previous constraints on q(n) to q(nq)
      REMOV = -1.
      NN = IABS(N)-1
4        NN = NN+1
         IF (WSETQ(NN).NE.0.)  GO TO 5          !remove constraint, return to 51
51       IF (NN.LT.NQ.and.M.LT.0) GO TO 4       !find any more constraints
c----------             Add constraint (m/100)*Q(n) = (m/100)*Qset
      REMOV = 1.                                         !to add equation
      NN = IABS(N)
      WSETQ(NN) = ww                                     !record weight & value
         IF (M.EQ.0)  GO TO 8                            !nothing to add
      VSETQ(NN) = QSET                                   !Record value set
5     Q(NN) = WSETQ(NN)                                  !
      Q(NP) = WSETQ(NN)*VSETQ(NN)                        !r.h.s.
      IF (REMOV.LT.0.)  WSETQ(NN) = 0.                   !record it's cancelled
c
52    DO 6  I = NN, NP                 !Givens transform, to add in new equation
         BII = B(I,I)
         R = SQRT( MAX(BII**2+Q(I)**2*REMOV, .1E-9) )
         C = BII/R
         S = Q(I)/R
            DO 6  J = I, NP
            QJ = Q(J)
            IF (I.EQ.NN.AND.J.NE.I.AND.J.NE.NP .and.nad.eq.0)  QJ = 0.
            Q(J) = C*QJ - S*B(I,J)
6     B(I,J) = C*B(I,J) + S*QJ*REMOV
c
      IF (REMOV.LT.0.)  GO TO 51                 !check for any more constraints
      GO TO 8
c=====------------------------------------------------------------------
c                                 Back substitution
7     B(NP,1) = M
      B(NP,NP)= SQRT(S)
8     devn = ABS(B(NP,NP))/SQRT(B(NP,1))
      DO 10  II = 2, NP
         I = NP-II+1
         S = B(I,NP)
            DO 9  J = I+1, NQ
9           IF (II.GT.2)  S = S - Q(J)*B(I,J)
10       Q(I) = S/B(I,I)
cc               if (lbug.ge.4) write(2,96) (q(i),i=1,np-1), devn
      q(19) = q(20)
      q(20) = devn
      RETURN
      END
c===============================================================================
      FUNCTION  SUMVAL (N, A, B, L)
c***********************************************************************jan77.
c  Returns the sum for  j = 1 to N  of:-
c     a(j)*b**j  at  L = 1 ;     j*a(j)*b**(j-1)  at  L = 2  (using b = B(1)).
c 2'93                     (j-1)*j*a(j)*b**(j-2)  at  L = 3  curvature, (2-n)
c     a(j) *b(j) at  L = 4 ;     a(j)*b(1-L+j*L)  at  L > 4  (giving  LL = 4).
c-----------------------------------------------------------------------------
      DIMENSION  A(*), B(*)
      double precision s                                                   !2'93
      J = N
      IF (L.EQ.3)  J = N - 1
      S = 0.
      JD = 1
      IF (L.GT.4)  JD = L
      JB = (N-1)*JD + 1
      LL = MIN(L, 4)
c.....................           loop from  j = n  to  j = 1.
  9   GO TO (1,2,3,4), LL
c                       L = 1;  evaluate polynomial  a(j)*b1**j
   1   S = (S + A(J)) * B(1)
      GO TO 8
c                       L = 2;  calculate gradient,   a(j)*j*b1**(j-1)
   2  S = S*B(1) + A(J)*J
      GO TO 8
c                       L = 3;  calculate curvature   a(j)*j*(j-1)*b1**(j-2)
   3  S = S*B(1) + A(J+1)*(J+1)*J                         !this starts at j= n-1
      GO TO 8
c                       L = 4;  calculate sum of products  a(j)*b(j)   (JD = 1)
c                       L > 4;  sum A(j)*B(1,j), with L= 1st dim of B  (JD = L)
   4  S = S + A(J)*B(JB)
      JB= JB - JD
c
  8   J = J-1
      IF (J.GT.0)  GO TO 9
c.....................
10    SUMVAL = S
      RETURN
      END
c===============================================================================
      FUNCTION  GIND (F, T)
c*****************************************************************jan77.
c  Group index subroutine, for o or x rays.      (f negative for x ray).
c  Gives (group index -1) to full machine accuracy,  for any value of t.
c Initialise by  call gind(GH,-dip) to set gyrofreq (MHz) and dip (deg).
c GH is ground value;  scaled to height h by  FH = gind(0.,h), with h>2.
c       An initial negative value of GH suppresses the scaling.
c-----------------------------------------------------------------------
      SAVE  GH,GHSN,GCSCT, FH,FHSN,FCSCT,HFH, C,C2
      data ibug / 0 /
c
         IF (F.EQ.0.)  GO TO 2
         IF (T.LT.0.)  GO TO 1
         T2 = MAX(T, .1E-9)**2
         goto 5
c                              store ground constants
1     GH = F
      FH = ABS(GH)
      DIP = MAX(-.01745329*T, .001)                !so can use DIP = 0
      GHSN = GH*SIN(DIP)
      GCSCT= (GH*COS(DIP))**2*.5/GHSN
      GO TO 3
c                              scale constants to height  h  =  t.
2     IF (GH.LT.0. .OR.T.LT.2.)  GO TO 4
               if(ibug.eq.1) write (2,*) 'gind fh to ht=',t
      FH = GH/(1.+T/6371.2)**3
3     FHSN = GHSN*FH/GH
      FCSCT= GCSCT*FH/GH
      HFH = .5*FHSN
         C = FCSCT+FHSN
         C2= C*C
4     GIND = FH
      RETURN
c                              calculate group index  gind  =  mu' - 1.
c                              f < 0. sections increase accuracy for x ray.
5     if (f.lt.0.) then
         X  = (F + FH)*T2
         G1 = X - FH
        else
         G1 = F*T2
        endif
      G2 = FCSCT/G1
      G3 = SQRT(G2*G2 + 1.)
      if (f.lt.0.) then
         G4 = FHSN*(G3 - G2)
         X  = X*(G1 - FH)
         G5 = -G4*X /(FHSN*(SQRT(X+C2) + C))
        else
         G4 = FHSN/(G2 + G3)
         G5 = G1 + G4
       endif
9     G6 = F + G4
      G7 = (F*G2*G4/G1 - HFH) *(F - G1)/(G3*G6)
      GIND = ABS(G7 + G6)/SQRT(G5*G6) -1.

      END Function GIND
