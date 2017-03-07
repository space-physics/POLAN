=====
POLAN
=====

POLAN is a classic Fortran program used to calculate real-height profiles from chirp ionosonde data from the ionosphere.

:revised:  Nov'87/ Mar'88/ Feb'93/ (Nov'95)
:author: J. E. Titheridge

.. contents::

Compile POLAN
=============
::

    gfortran polrun.f polan.f POLMIS.FOR POLSIN.FOR POLSUB.FOR -o polan

Basic POLAN Example
===================
::

    ./polan examples/in.dat

this creates a big output text file ``out.dat``

POLynomial ANalysis subroutine
================================
::

    POLAN  (N, FV,HT, QQ, FB,DIP,   START, AMODE, VALLEY, LIST).

for the calculation of real-height profiles from sweep-frequency ionograms.

     This outline is for Versions after September 1986.  It supplements the
information contained in the 194-page report "Ionogram analysis with the
generalised program POLAN",  obtainable as report UAG-93 from:-
     NOAA/NGDC,  325 Broadway,   E/GC2, Dept. 884
     Boulder,    Colorado 80303, USA.   Phone (303)497-6761.

     If problems arise run one data set with  LIST = 3, and mail all output
to:  J.E. Titheridge,  Physics Dept.,  University of Auckland,  New Zealand.
     or (better) use EMAIL, to:  J.TITHERIDGE@AUCKLAND.AC.NZ

NOTE: Profile summaries, and any debug information, are written to unit 2.
    Changes made in recent versions of POLAN are summarised in section F. below.

A. DATA ARRAYS
--------------

  Polan is called with frequency, height data in the arrays  FV, HT.
        The dimension of these arrays (NDIM) must be greater than
        30 + the number of data points in the arrays.
    --> For the first call to POLAN, you must set  N = NDIM.
  QQ is an output array, used primarily with single-polynomial calculations.
        For 3 layers and 8-term polynomials, the dimension of QQ must exceed 50.
        POLAN will not write over a value of -1. in QQ, so setting the last 
        element of QQ equal to -1.0 will prevent any overflow.  
        The data returned in QQ is described in D.2 below.
    --> If this data is not required, use Dimension QQ(2) and set QQ(1)= -1.

  In the data arrays, intermediate layers are terminated by a scaled (or 
        zero) value for the critical frequency, with:-
            h' = 0.0  for a Chapman peak and normal valley,
            h' = 10.0 for a peak  with no following valley,
            h' = negative and less than 50 is used to set the valley
                     constants, for this valley only, as described below.
            h' = negative (equal to minus the scaled virtual height)
                     for a cusp-type discontinuity only.
        Note that profiles are normally continuous across a cusp point,
        so h' is scaled normally.  (Or preferably scale points either side
        of a cusp, and not at the cusp itself;  see JATP 44,657,1982.)

  The o-ray FC (scaled or zero) may be followed by an x-ray value (-FCX).
 
  The final layer is terminated by at least 2 null points, with  h = f = 0.
  Data can be terminated without a peak by using a final frequency of -1.0.

     Data for the extraordinary ray, if any, precedes the o-ray data for
  each layer.  This is because x-ray data is used only to calculate the
  (start or valley) corrections to be made at the beginning of the
  calculation for that layer.  x-ray data are distinguished by using -f.

     The format for input data is best seen by study of the examples in the 
 test file ``examples/in.dat``.

B.  INPUT PARAMETERS
-------------------- 
Input parameter in the call to POLAN are here described.
 
   FB  gives the gyrofrequency at the ground in MHz, for an inverse cube 
variation.   If you have only the gyrofrequency  FH  at a height  h km,
the ground value is obtained from    FB = FH * (1. + h/6371.2)**3.
   To use a gyrofrequency (FH, say) which is independent of height, 
set  FB = - FH.
------------------------

   DIP  is the magnetic dip angle  IN  DEGREES.   Use of a negative value
for  DIP suppresses the physical checks which are normally applied to the
calculated profile,  so that the result obtained is the best mathematical
(but possibly non-physical) fit to the virtual-height data.  
  [Some physically based equations are still included in start and valley 
calculations, unless AMODE is negative.]
------------------------

   START normally gives a model height at 0.5 MHz.  Typical values are:    
noon   sunset-2/rise+2hr   set/rise    set+1hr   set+2   set+4 to rise-1    
85km    88km(E layer)    90(E)/80(F)   100 km    130 km     150 km. 
 
   A preferred procedure is to calculate model values of START from the 
equations (10) to (13) given in J. Atmosph. Terr. Phys. 48, 435-446, 1986.

   Use of START = 0.0 makes some allowance for underlying ionisation based 
on a limited extrapolation of the first few virtual heights.

   With initial x-ray data, START is taken to give the gyrofrequency height
for underlying ionisation calculations; the values listed above are still
suitable for this purpose.  The x-ray data is used to calculate a slab start
correction from 0.3*fmin  (adding points at 0.3, 0.6 and 0.8 *fmin). 

[Alternative procedures can be obtained using non-standard values of START:-
   START between 0. and 44.  defines the plasma frequency for a model start.
   Start = -1.0   uses a direct start, from the first scaled point.
   Start < -1.0   for x-starts to use a polynomial from (-Start -1.0) MHz. ]
------------------------

THE final three parameters - AMODE, VALLEY and LIST, are zero for most work.

   AMODE  sets the type of analysis, as listed below.   Zero uses mode 6.
     Use Amode+10. for 12-point integrals, for high accuracy at large dip
     angles (this is done automatically, at  DIP > 60, when Amode=0).
   For denser (e.g. digital) data, with more than 30 points in one layer,
     use a higher-order mode.  Thus AMODE = 9. gives maximum detail,  or
     AMODE= 95. gives single-polynomials with 5, 9 terms for the E, F2 layers.

  Values of Amode greater than 29.0 are used to specify the number of
     polynomial constants to be used to describe each ionospheric layer.
     e.g. 80.  uses an 8-term real height polynomial for each separate layer.
          85.  uses 8 terms for the final layer and 5 terms for lower layers.
          853. uses 8 terms for the last, 3 terms for the first, and 5 terms
               for any intermediate layer.

     Setting AMODE negative causes physical relations to be omitted from the
start and valley calculations. 
------------------------

   VALLEY= 0.0 or 1.0  uses a valley width equal to the initial default
value of twice the local scale height.  The initial default depth is 0.05
MHz.  The calculated depth is scaled according to (calculated width)**2. 

     Alternative solutions may be obtained as follows:

  VALLEY = 10.0  gives a monotonic (no valley) analysis.
  Valley =  5.0  gives a maximum valley (upper reasonable limit) analysis.
  Valley =  0.1 to 5.0  multiplies the standard valley width by this factor.
  Valley = -.01 to -.99 uses  -2.0 * Valley  as the initial depth
                             (instead of the default value of 0.05 MHz).
  Valley = -1.0  iterates both valley depth and width for best fit, with 
              x-ray data.  (-1.D iterates from an initial depth of 0.D MHz).
  Valley = -2.0 to -50. specifies a fixed valley width of 2*int(-Valley) km.
                        Any decimal part D specifies a depth of 2*D in MHz.
------------------------

  LIST = 0   prints results for the start, peak and valley regions only.
         1   adds one line of output showing the frequency range and the
             polynomial coefficients calculated at each step.
         2, 3   add additional output.
         4 to 9 show the data used at each step, and the calculated
                polynomial coefficients:
            5   shows each set of simult equations, in the call to SOLVE;
            6/7/8/9 give detail in the start/reduction/peak/valley steps.

         LIST negative  suppresses most trace output below the first peak.
         LIST= -10 suppresses all output, even the normal layer summaries.

C.  OUTPUT PARAMETERS,  returned by POLAN.
------------------------------------------
 
  The arrays  FV, HT contain the calculated frequencies and real heights.

  N  gives the number of calculated real-height data points.

  The peak of the last layer is at  FC = fv(N-3),  Hmax = ht(N-3).
  A point at (N-4) is added, on the fitted Chapman-layer peak; this and the
         points above the peak permit accurate 2nd-difference interpolation.
  Points at  N-2, N-1 and  N  in the output arrays are extrapolated heights
         at  0.35, 0.85 and 1.5  scale heights above the peak (calculated from
             the Chapman expression with a scale height gradient of 0.1).

  fv(N+1)  gives the standard error of the last critical frequency, in MHz.
  ht(N+1)  gives the standard error of the last peak height  Hmax,  in km.
  fv(N+2)  gives the slab thickness, in km.   This is equal to the 
             sub-peak electron content divided by the peak density.
  ht(N+2)  gives the scale height SH of the last peak, in km.
             A negative value of SH shows that a model value was used for
             the scale height, to limit an unreasonable peak extrapolation.

  QQ returns the real-height coefficients, for single-polynomial calculations,
             as described under D.2 below.  For overlapping polynomial modes,
             coefficients are returned for the last polynomial in each layer.

D.  MODES OF ANALYSIS
---------------------

D.1 THE TEN STANDARD MODES
~~~~~~~~~~~~~~~~~~~~~~~~~~

    MODE is obtained from the input parameter AMODE, modified to the range 
    1 to 10, and is used to select the type of analysis as summarised below.
    All Modes include an estimated start correction,  a Chapman-layer peak,
    and a model valley between layers.

MODE=1.- The Linear-Lamination analysis.
     2.- A Parabolic-Lamination analysis, matching end gradients  ( = Paul).
     3.- Overlapping Cubics, with no spurious oscillations (JATP 1982 p657).
     4.- Fourth Order Overlapping Polynomials   (Radio Science 1967, p1169).
     5.- Fifth Order Least-Squares fit to 6 points  (4 virtual + 2 real).
     6.- Sixth Order Least-Squares fit to 8 points  (5 virtual + 3 real).
     7.- Sixth Order fit to 7 virtual +3 real heights; calculates 2 new hts.
     8.- Sixth Order fit to 8 virtual +4 real heights; calculates 2 new hts.
     9.- Seventh Order fit to 13 virtual + 6 real hts; calculates 3 new hts.
     10. A Single Polynomial,  fitting  2*sqrt(NV)  terms to  NV heights.
         A maximum of 90 (=MAXB-9) points can be included in one polynomial.

   The basic parameters which define the type of analysis depend on the
parameter MODE, and are obtained from the arrays given below.  
   NT is the number of terms used in the polynomial representation of each 
real-height segment.
   NV is the number of virtual heights which are fitted in this step.
   NR is the number of previously-calculated real heights which are fitted
(in addition to the origin FA, HA).  A negative value of NR indicates that
one of the fitted real heights is below the origin.   If  NT = NV + NR  we
get an exact fit to the data, and if  NT < NV + NR  the calculated profile
segment is a least-squares fit. 

   NH is the number of new real heights to be calculated.  
   'First step' values are used at the beginning of an analysis, or when
starting on a new layer, when no real heights are known above the starting
point.  In this case the number of known real heights is zero, and the
tabulated values of NR define the position of the origin (counting backwards
from the last calculated real height) for the following step. 

       |-------- First step --------|    |------- Following steps --------|
MODE=  1, 2  3  4  5  6  7  8   9  10    1  2  3  4   5   6   7   8   9  10 
 NT =  1, 2, 3, 4, 4, 5, 6, 6,  7, 73,   1, 2, 3, 4,  5,  6,  6,  6,  7, 73
 NV =  1, 2, 3, 4, 5, 7, 8,10, 12, 90,   1, 1, 2, 3,  4,  5,  7,  8, 13, 90
 NR =  0, 0, 0, 1, 1, 2, 2, 3,  5,  2,   0,-1,-1, 1, -2, -3, -3, -4, -6, -3
 NH =  1, 1, 2, 3, 3, 4, 5, 6,  8, 28,   1, 1, 1, 1,  1,  1,  2,  2,  3, 28


D.2 SINGLE-POLYNOMIAL MODES
~~~~~~~~~~~~~~~~~~~~~~~~~~~

  These use a defined number of real-height coefficients for each layer, 
and return all profile parameters in the array QQ.  The order of the 
analysis is set by the parameter AMODE, as follows.

AMODE = 10L,  where L is an integer in the range 3 to 14, uses a single
              polynomial with L terms to describe each ionospheric layer.
AMODE = 10L+M   uses  L terms for the final layer, and M for earlier layers.
AMODE = 100L+10M+F is L terms for Last, M for Middle and  F for First layer
                                             (M and F must be less than 10).

QQ
++

returns the real-height parameters which describe the profile, for
single-polynomial modes of analysis (unless QQ(1) was set equal to -1.0 by
the calling program).  (For normal [overlapping polynomial] runs, QQ returns
the coefficients for the last polynomial, and the peak, in each layer.)

The returned value of QQ(1) gives the total number of stored values (numq).
Starting at QQ(2), the parameters returned for each layer are:
     FA, HA,  nq,  q1, q2, .. qn,  devn,   FP, FC, Hmax, and SH.

nq
++

is the number of polynomial coefficients (q1 to qn) used for this layer.
This is normally equal to the number of coefficients requested in AMODE.
   
HA is the true height at FA, after any start or valley adjustments, so the 
real-height profile is 
              h  =  HA + q1.(f-FA) + q2.(f-FA)^2 + ... qn.(f-FA)^nq.

devn is the rms deviation (in km) of the fit to the virtual height data.

FC, Hmax and SH
+++++++++++++++

 are the constants which define the Chapman-layer peak;
this joins the polynomial section at the frequency FP (close to the second to
highest scaled frequency for the layer, but limited to 0.9FM < FP < 0.97FC).

   For a 2nd (or 3rd) layer,  FA, HA give the new real-height origin at the 
top of the valley region.   Thus FA is equal to the previous FC,  and the
valley width is   W = HA - Hmax  in km.   The valley depth (D, in MHz) can be
obtained from the width using equations (14) of the report UAG-93, which give
     D = 0.008 W**2/(20 + W) MHz,  followed by   D = D.FC/(D + FC).

   The end point of the data in QQ is verified by a value  QQ(numq+1) = -99.
for a normal exit, and  -98. for an error (or no-peak) exit.

E.   PROCESSING 
---------------
 Outline of the REAL-HEIGHT ANALYSIS LOOP within POLAN.

E.1  THE OVERALL PROCEDURE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FOR ONE CYCLE OF THE CALCULATION

Analysis can proceed with any number of scaled virtual heights (even
1 height and no critical frequency) for each layer.  If the number of data
points NV is less than the number of polynomial terms NT (as specified by 
AMODE), NT is automatically decreased.

-    Calculate one polynomial, with NT terms, from the point  FA = fv(K),
HA = ht(K)  to fit the next NV virtual and NR real heights.  (The fitted 
real heights include one point below HA, if NR is negative.)   
The real-height origin (FA,HA) is at K = KR, in the data arrays FV, HT;
the corresponding virtual height is at K = KV. 

-    With x-ray data (-ve frequencies), at the start or after a peak,
recalculate HA to include the correction for underlying or valley ionisation. 

-    Calculate a further NH real heights, and set KR = KR + NH; KV = KV + NH.
                                                                           
-    Repeat this loop, calculating successive overlapping real-height
sections, until a critical frequency (or end-of-layer) is found in the range
KV +1  to  KV +NV +1.   Then calculate real heights at the remaining scaled 
frequencies and determine a least-squares Chapman-layer peak. 

E.2  INDIVIDUAL STEPS WITHIN EACH CYCLE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
numbered according to the corresponding section in the program POLAN.

SECTION 2.2  Count initial x-rays.  Check frequency sequencing.
             Check for cusp, peak, or end of data.
 Set NF = number of o-rays 
       (= NV, if sufficient points exist before a following peak);
     NX = number of x-rays;      MV = NF+NX.
     FM = fv(mf) = the top frequency used in this step.
     FCC= FC or 0.1 for a peak,  = -.1 for a cusp (gradient discontinuity)
                        at FM,   = 0.0 otherwise.

SECTION 2.3  Subtract the group retardation due to the last calculated
                real-height section.
     This modifies all the virtual heights at f > FA  (where FA = fv(KR)),
     and increases the index LK (which gives the point up to which the
     group retardation has been removed) to KR.

SECTION 3.  Set up equations for the next profile step.

          Check for the occurrence of a valley; if this is required, set
     the valley flag HVAL and set initial values for the width and depth.

          Set up equations in the matrix B.   For start calculations using 
     x-ray data, or for any valley calculations, add suitably weighted
     equations specifying desired physical properties of the solution.

SECTION 4.  Solve the set of simultaneous equations in the array B.

          Check that the solution satisfies basic physical constraints.
     If it does not, obtain a new least-squares solution with the limiting 
     constraints imposed (in the subroutine ADJUST).

          For an x-start or valley calculation, iterate the solution as
     required to ensure the use of a correct gyrofrequency height, and 
     the correct relation between depth and width of the valley.
          For an o-ray valley, loop once to adjust the valley depth.

SECTION 5.  Calculate and store the real heights.

          Set KRM as the index for the highest calculated real height.

SECTION 6.  Least-squares fitting of a Chapman layer peak.

          Calculate the critical frequency and the scale height of a
     layer peak, by an iterative fit to the real-height gradients at the 
     last few calculated points  (as in Radio Science 20, 247, 1985).
          Determine the height of the peak by fitting the peak shape to a 
     weighted mean of the last few calculated real heights.  Adjust the
     last real height to agree closely with the Chapman peak (Sept'86).
     Add an interpolated point between the 'last' height and the peak(2'93).

SECTION 7.  Go to section 2, to restart for a new layer.

     If there are no further data:-   add one point half-way to the peak;
extrapolate 3 points for the topside ionosphere (assuming a Chapman layer
with a scale height gradient of 0.1 km/km);  store constants relating to
the last layer peak;  and return.

F.  POLAN changes made in recent versions.
------------------------------------------

F.2 CHANGES  February 1993 (marked 2'93)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Deleted NDIM from call.  First call must have N = NDIM (or ndim is set to 100).
Added extra output point below peak, and spaced those above for best interpoln.
Reduced # points over which Chapman peak is fitted, for single polynomials.
Imposed lower limit on profile curvature at top point, before peak fit.
Mode 10 to NT= 2.*sqrt(NV), so 20/40/60 data -> Nt= 9/13/14 (prev NV>18->NT=15)

NOTE: I now use ! for comments; you may need to change this for your compiler.

F.1 CHANGES  September 1986
~~~~~~~~~~~~~~~~~~~~~~~~~~~

(a)  Addition of the parameters  NDIM  and  QQ  in the call to POLAN.
     Use of NDIM makes it unnecessary to reset N (to the dimension of the
input arrays) on each call.

     QQ returns the coefficients for single-polynomial representations.  
It is now a required parameter in the call to POLAN,  but is not used if
(initially) QQ(1) = -1.   (Previous use of QQ returned 1 less coefficient 
than described in section D.2, since the count nq was taken to include
the constant HA).  For normal (overlapping polynomial) runs, QQ returns the
coefficients for the last polynomial, and the peak, in each layer.

(b)  Use of a negative scale height, to indicate use of a model value rather
than one derived from the data, is restricted to the output listing (and the
output array QQ).  In some previous versions, -SH was accidentally carried
over to later stages creating numerous problems. 

(c)  The default analysis (obtained at AMODE = 0.0) has been changed from
Mode 5 to Mode 6.  Experience has shown some benefits and no problems with
the higher modes, particularly since the change (d) below which gives good
results even when the scaled frequency interval varies considerably. 

(d)  Weighting of different points in the least-squares calculation has
been made proportional to the scaled frequency interval.  This stops smooth
sections of the profile, where fewer points may have been scaled, from
getting too low a weight.  It reduces spurious fluctuations in high order
modes to well below the levels described in J. Atmosph. Terr. Phys. 44,
657-669, 1982. 

(e)  The START model has been revised to the procedure described in J.
Atmosph. Terr. Phys. 48, 435-446, 1986. 

(f)  Minor improvements have been made in several steps of the calculation. 
Programs will now run at DIP = 0.  Calculations proceed normally with 2 or
more data points for each layer;  even a layer with only one point (with
or without FC) is handled.

(g)  Descriptive comments have been extracted from the listing of POLAN.FOR (polan.f),
into this file.

