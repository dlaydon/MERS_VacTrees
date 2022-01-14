/* randlib_par.cpp, part of the Global Epidemic Simulation v2.0 BETA
   Parallel random number generator and distributions
 
   Copyright 2013, MRC Centre for Outbreak Analysis and Modelling
 
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include "randlib_par.h"
#include "math.h"
#include "omp.h"

#include <stdio.h>
#include <stdlib.h>
#include "DataStructure.h"



#define CACHE_LINE_SIZE 128
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define minF(a,b) ((a) <= (b) ? (a) : (b))
#define maxF(a,b) ((a) >= (b) ? (a) : (b))
#define ERR_CRITICAL(x) {fprintf(stderr,x);exit(1);}

/* RANDLIB static variables */
long Xm1,Xm2,Xa1,Xa2,*Xcg1,*Xcg2,Xa1vw,Xa2vw;

#include "Libraries.h"

vector<int> SampleFromContainer(vector<int> ContainerToSample, int NumElementsToSample)
{
	vector<int> OutputVector;
	std::sample(ContainerToSample.begin(), ContainerToSample.end(), std::back_inserter(OutputVector), NumElementsToSample, mt19937{ std::random_device{}() });
	//std::cout << NumElementsToSample << " random letters out of Container " << " : " << out << '\n';

	return (OutputVector); 

	///* ORIG - from https://en.cppreference.com/w/cpp/algorithm/sample*/
	//std::string in = "abcdefgh", out;
	//std::sample(in.begin(), in.end(), std::back_inserter(out),
	//	5, std::mt19937{ std::random_device{}() });
	//std::cout << "five random letters out of " << in << " : " << out << '\n';
}
double RandomNumber_Mt(int thread_no)
{
	//// wrapper of ranf_mt function - used to clean ifdefs etc out of main code 
	double RandomNo = 0;
	RandomNo = ranf_mt(thread_no); /// this will change depending on thread_no (obviously), and thus answer will change with max_threads.
	//while (RandomNo == 0) RandomNo = ranf_mt(thread_no);
	return RandomNo;
}
double RandomNumber_Mt() //// non-threaded version defaults to using "thread" zero.
{
	return RandomNumber_Mt(0);
}
bool BernoulliTrial(double ProbSuccess, int thread_no)
{
	double randomNo		= RandomNumber_Mt(thread_no); /// draws random number between 0 and 1, keeping track of threading. 
	bool EventSuccess	= randomNo < ProbSuccess;

	return EventSuccess;
}
bool BernoulliTrial(double ProbSuccess)
{
	return BernoulliTrial(ProbSuccess, 0); 
}

double ranf(void)
{
long k,s1,s2,z;
unsigned int curntg;
double dev;

#ifdef _USE_OPENMP
	curntg=CACHE_LINE_SIZE*omp_get_thread_num();
#else
	curntg=0;
#endif
    s1 = Xcg1[curntg];
    s2 = Xcg2[curntg];
   
    k = s1/53668L;
    s1 = Xa1*(s1-k*53668L)-k*12211;
    if(s1 < 0) s1 += Xm1;
    k = s2/52774L;
    s2 = Xa2*(s2-k*52774L)-k*3791;
    if(s2 < 0) s2 += Xm2;
    Xcg1[curntg] = s1;
    Xcg2[curntg] = s2;
    z = s1-s2;
    if(z < 1) z += (Xm1-1);
    dev=((double) z)/((double) Xm1);
	return dev;
}

double ranf_mt(int tn)
{
  long k,s1,s2,z;
	int curntg;
	double dev;
	curntg=CACHE_LINE_SIZE*tn;
	s1 = Xcg1[curntg];
	s2 = Xcg2[curntg];
	k = s1/53668L;
	s1 = Xa1*(s1-k*53668L)-k*12211;
	if(s1 < 0) s1 += Xm1;
	k = s2/52774L;
	s2 = Xa2*(s2-k*52774L)-k*3791;
	if(s2 < 0) s2 += Xm2;
	Xcg1[curntg] = s1;
	Xcg2[curntg] = s2;
	z = s1-s2;
	if(z < 1) z += (Xm1-1);
	dev=((double) z)/((double) Xm1);
	return dev;
}

void initSeeds(long iseed1,long iseed2, int max_threads)
/*
**********************************************************************
     void setall(long iseed1,long iseed2)
               SET ALL random number generators
     Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
     initial seeds of the other generators are set accordingly, and
     all generators states are set to these seeds.
     This is a transcription from Pascal to Fortran of routine
     Set_Initial_Seed from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     iseed1 -> First of two integer seeds
     iseed2 -> Second of two integer seeds
**********************************************************************
*/
{
	int g;
    // Max 32 threads?
  	Xcg1=(long *) calloc(max_threads*CACHE_LINE_SIZE,sizeof(long));
    Xcg2=(long *) calloc(max_threads*CACHE_LINE_SIZE,sizeof(long));
   
  
    Xm1 = 2147483563L;
    Xm2 = 2147483399L;
    Xa1 = 40014L;
    Xa2 = 40692L;
    Xa1vw = 2082007225L;
    Xa2vw = 784306273L;
    
    
    *Xcg1 = iseed1;
    *Xcg2 = iseed2;

    for (g=1; g<max_threads; g++) {
        *(Xcg1+(g*CACHE_LINE_SIZE)) = mltmod(Xa1vw,*(Xcg1+((g-1)*CACHE_LINE_SIZE)),Xm1);
        *(Xcg2+(g*CACHE_LINE_SIZE)) = mltmod(Xa2vw,*(Xcg2+((g-1)*CACHE_LINE_SIZE)),Xm2);
    }
    
    //for (g=0; g<7; g++) {
//      printf("Random seeds[%d] = (%d,%d)\n",g,Xcg1[CACHE_LINE_SIZE*g],Xcg2[CACHE_LINE_SIZE*g]);
    //}
}

long mltmod(long a,long s,long m)
/*
**********************************************************************
     long mltmod(long a,long s,long m)
                    Returns (A*S) MOD M
     This is a transcription from Pascal to Fortran of routine
     MULtMod_Decompos from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     a, s, m  -->
**********************************************************************
*/
{
#define h 32768L
long mltmod,a0,a1,k,p,q,qh,rh;
/*
     H = 2**((b-2)/2) where b = 32 because we are using a 32 bit
      machine. On a different machine recompute H
*/
    if(!(a <= 0 || a >= m || s <= 0 || s >= m)) goto S10;
    fputs(" a, m, s out of order in mltmod - ABORT!\n",stderr);
    fprintf(stderr," a = %12ld s = %12ld m = %12ld\n",a,s,m);
    fputs(" mltmod requires: 0 < a < m; 0 < s < m\n",stderr);
    exit(1);
S10:
    if(!(a < h)) goto S20;
    a0 = a;
    p = 0;
    goto S120;
S20:
    a1 = a/h;
    a0 = a-h*a1;
    qh = m/h;
    rh = m-h*qh;
    if(!(a1 >= h)) goto S50;
    a1 -= h;
    k = s/qh;
    p = h*(s-k*qh)-k*rh;
S30:
    if(!(p < 0)) goto S40;
    p += m;
    goto S30;
S40:
    goto S60;
S50:
    p = 0;
S60:
/*
     P = (A2*S*H)MOD M
*/
    if(!(a1 != 0)) goto S90;
    q = m/a1;
    k = s/q;
    p -= (k*(m-a1*q));
    if(p > 0) p -= m;
    p += (a1*(s-k*q));
S70:
    if(!(p < 0)) goto S80;
    p += m;
    goto S70;
S90:
S80:
    k = p/qh;
/*
     P = ((A2*H + A1)*S)MOD M
*/
    p = h*(p-k*qh)-k*rh;
S100:
    if(!(p < 0)) goto S110;
    p += m;
    goto S100;
S120:
S110:
    if(!(a0 != 0)) goto S150;
/*
     P = ((A2*H + A1)*H*S)MOD M
*/
    q = m/a0;
    k = s/q;
    p -= (k*(m-a0*q));
    if(p > 0) p -= m;
    p += (a0*(s-k*q));
S130:
    if(!(p < 0)) goto S140;
    p += m;
    goto S130;
S150:
S140:
    mltmod = p;
    return mltmod;
#undef h
}

long ignbin(long n,double pp)
/*
**********************************************************************
     long ignbin(long n,double pp)
                    GENerate BINomial random deviate
                              Function
     Generates a single random deviate from a binomial
     distribution whose number of trials is NPat and whose
     probability of an event in each trial is P.
                              Arguments
     n  --> The number of trials in the binomial distribution
            from which a random deviate is to be generated.
	    JJV (NPat >= 0)
     pp --> The probability of an event in each trial of the
            binomial distribution from which a random deviate
            is to be generated.
	    JJV (0.0 <= PP <= 1.0)
     ignbin <-- A random deviate yielding the number of events
                from NPat independent trials, each of which has
                a probability of event P.
                              Method
     This is algorithm BTPE from:
         Kachitvichyanukul, V. and Schmeiser, B. W.
         Binomial Random Variate Generation.
         Communications of the ACM, 31, 2
         (February, 1988) 216.
**********************************************************************
     SUBROUTINE BTPEC(NPat,PP,ISEED,JX)
     BINOMIAL RANDOM VARIATE GENERATOR
     MEAN .LT. 30 -- INVERSE CDF
       MEAN .GE. 30 -- ALGORITHM BTPE:  ACCEPTANCE-REJECTION VIA
       FOUR REGION COMPOSITION.  THE FOUR REGIONS ARE A TRIANGLE
       (SYMMETRIC IN THE CENTER), A PAIR OF PARALLELOGRAMS (ABOVE
       THE TRIANGLE), AND EXPONENTIAL LEFT AND RIGHT TAILS.
     BTPE REFERS TO BINOMIAL-TRIANGLE-PARALLELOGRAM-EXPONENTIAL.
     BTPEC REFERS TO BTPE AND "COMBINED."  THUS BTPE IS THE
       RESEARCH AND BTPEC IS THE IMPLEMENTATION OF A COMPLETE
       USABLE ALGORITHM.
     REFERENCE:  VORATAS KACHITVICHYANUKUL AND BRUCE SCHMEISER,
       "BINOMIAL RANDOM VARIATE GENERATION,"
       COMMUNICATIONS OF THE ACM, FORTHCOMING
     WRITTEN:  SEPTEMBER 1980.
       LAST REVISED:  MAY 1985, JULY 1987
     REQUIRED SUBPROGRAM:  RAND() -- A UNIFORM (0,1) RANDOM NUMBER
                           GENERATOR
     ARGUMENTS
       NPat : NUMBER OF BERNOULLI TRIALS            (INPUT)
       PP : PROBABILITY OF SUCCESS IN EACH TRIAL (INPUT)
       ISEED:  RANDOM NUMBER SEED                (INPUT AND OUTPUT)
       JX:  RANDOMLY GENERATED OBSERVATION       (OUTPUT)
     VARIABLES
       PSAVE: VALUE OF PP FROM THE LAST CALL TO BTPEC
       NSAVE: VALUE OF NPat FROM THE LAST CALL TO BTPEC
       XNP:  VALUE OF THE MEAN FROM THE LAST CALL TO BTPEC
       P: PROBABILITY USED IN THE GENERATION PHASE OF BTPEC
       FFM: TEMPORARY VARIABLE EQUAL TO XNP + P
       M:  INTEGER VALUE OF THE CURRENT MODE
       FM:  FLOATING POINT VALUE OF THE CURRENT MODE
       XNPQ: TEMPORARY VARIABLE USED IN SETUP AND SQUEEZING STEPS
       P1:  AREA OF THE TRIANGLE
       NumC:  HEIGHT OF THE PARALLELOGRAMS
       XM:  CENTER OF THE TRIANGLE
       XL:  LEFT END OF THE TRIANGLE
       XR:  RIGHT END OF THE TRIANGLE
       AL:  TEMPORARY VARIABLE
       XLL:  RATE FOR THE LEFT EXPONENTIAL TAIL
       XLR:  RATE FOR THE RIGHT EXPONENTIAL TAIL
       P2:  AREA OF THE PARALLELOGRAMS
       P3:  AREA OF THE LEFT EXPONENTIAL TAIL
       P4:  AREA OF THE RIGHT EXPONENTIAL TAIL
       U:  A U(0,P4) RANDOM VARIATE USED FIRST TO SELECT ONE OF THE
           FOUR REGIONS AND THEN CONDITIONALLY TO GENERATE A VALUE
           FROM THE REGION
       V:  A U(0,1) RANDOM NUMBER USED TO GENERATE THE RANDOM VALUE
           (REGION 1) OR TRANSFORMED INTO THE VARIATE TO ACCEPT OR
           REJECT THE CANDIDATE VALUE
       IX:  INTEGER CANDIDATE VALUE
       X:  PRELIMINARY CONTINUOUS CANDIDATE VALUE IN REGION 2 LOGIC
           AND A FLOATING POINT IX IN THE ACCEPT/REJECT LOGIC
       K:  ABSOLUTE VALUE OF (IX-M)
       F:  THE HEIGHT OF THE SCALED DENSITY FUNCTION USED IN THE
           ACCEPT/REJECT DECISION WHEN BOTH M AND IX ARE SMALL
           ALSO USED IN THE INVERSE TRANSFORMATION
       R: THE RATIO P/Q
       G: CONSTANT USED IN CALCULATION OF PROBABILITY
       MP:  MODE PLUS ONE, THE LOWER INDEX FOR EXPLICIT CALCULATION
            OF F WHEN IX IS GREATER THAN M
       IX1:  CANDIDATE VALUE PLUS ONE, THE LOWER INDEX FOR EXPLICIT
             CALCULATION OF F WHEN IX IS LESS THAN M
       I:  INDEX FOR EXPLICIT CALCULATION OF F FOR BTPE
       AMAXP: MAXIMUM ERROR OF THE LOGARITHM OF NORMAL BOUND
       YNORM: LOGARITHM OF NORMAL BOUND
       ALV:  NATURAL LOGARITHM OF THE ACCEPT/REJECT VARIATE V
       X1,F1,Z,W,Z2,X2,F2, AND W2 ARE TEMPORARY VARIABLES TO BE
       USED IN THE FINAL ACCEPT/REJECT TEST
       QN: PROBABILITY OF NO SUCCESS IN NPat TRIALS
     REMARK
       IX AND JX COULD LOGICALLY BE THE SAME VARIABLE, WHICH WOULD
       SAVE A MEMORY POSITION AND A LINE OF CODE.  HOWEVER, SOME
       COMPILERS (E.G.,CDC MNF) OPTIMIZE BETTER WHEN THE ARGUMENTS
       ARE NOT INVOLVED.
     ISEED NEEDS TO BE DOUBLE PRECISION IF THE IMSL ROUTINE
     GGUBFS IS USED TO GENERATE UNIFORM RANDOM NUMBER, OTHERWISE
     TYPE OF ISEED SHOULD BE DICTATED BY THE UNIFORM GENERATOR
**********************************************************************
*****DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY
*/
{
/* JJV changed initial values to ridiculous values */
double psave = -1.0E37;
long nsave = -214748365;
long ignbin,i,ix,ix1,k,m,mp,T1;
double al,alv,amaxp,c,f,f1,f2,ffm,fm,g,p,p1,p2,p3,p4,q,qn,r,u,v,w,w2,x,x1,
    x2,xl,xll,xlr,xm,xnp,xnpq,xr,ynorm,z,z2;

//S10:
/*
*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE
JJV added checks to ensure 0.0 <= PP <= 1.0
*/
    if(pp < 0.0F) ERR_CRITICAL("PP < 0.0 in IGNBIN");
    if(pp > 1.0F) ERR_CRITICAL("PP > 1.0 in IGNBIN");
    psave = pp;
    p = minF(psave,1.0-psave);
    q = 1.0-p;
//S20:
/*
JJV added check to ensure NPat >= 0
*/
    if(n < 0L) ERR_CRITICAL("NPat < 0 in IGNBIN");
    xnp = n*p;
    nsave = n;
    if(xnp < 30.0) goto S140;
    ffm = xnp+p;
    m = (long) ffm;
    fm = m;
    xnpq = xnp*q;
    p1 = (long) (2.195*sqrt(xnpq)-4.6*q)+0.5;
    xm = fm+0.5;
    xl = xm-p1;
    xr = xm+p1;
    c = 0.134+20.5/(15.3+fm);
    al = (ffm-xl)/(ffm-xl*p);
    xll = al*(1.0+0.5*al);
    al = (xr-ffm)/(xr*q);
    xlr = al*(1.0+0.5*al);
    p2 = p1*(1.0+c+c);
    p3 = p2+c/xll;
    p4 = p3+c/xlr;
S30:
/*
*****GENERATE VARIATE
*/
    u = ranf()*p4;
    v = ranf();
/*
     TRIANGULAR REGION
*/
    if(u > p1) goto S40;
    ix = (long) (xm-p1*v+u);
    goto S170;
S40:
/*
     PARALLELOGRAM REGION
*/
    if(u > p2) goto S50;
    x = xl+(u-p1)/c;
    v = v*c+1.0-ABS(xm-x)/p1;
    if(v > 1.0 || v <= 0.0) goto S30;
    ix = (long) x;
    goto S70;
S50:
/*
     LEFT TAIL
*/
    if(u > p3) goto S60;
    ix = (long) (xl+log(v)/xll);
    if(ix < 0) goto S30;
    v *= ((u-p2)*xll);
    goto S70;
S60:
/*
     RIGHT TAIL
*/
    ix = (long) (xr-log(v)/xlr);
    if(ix > n) goto S30;
    v *= ((u-p3)*xlr);
S70:
/*
*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST
*/
    k = ABS(ix-m);
    if(k > 20 && k < xnpq/2-1) goto S130;
/*
     EXPLICIT EVALUATION
*/
    f = 1.0;
    r = p/q;
    g = (n+1)*r;
    T1 = m-ix;
    if(T1 < 0) goto S80;
    else if(T1 == 0) goto S120;
    else  goto S100;
S80:
    mp = m+1;
    for(i=mp; i<=ix; i++) f *= (g/i-r);
    goto S120;
S100:
    ix1 = ix+1;
    for(i=ix1; i<=m; i++) f /= (g/i-r);
S120:
    if(v <= f) goto S170;
    goto S30;
S130:
/*
     SQUEEZING USING UPPER AND LOWER BOUNDS ON ALOG(F(X))
*/
    amaxp = k/xnpq*((k*(k/3.0+0.625)+0.1666666666666)/xnpq+0.5);
    ynorm = -(k*k/(2.0*xnpq));
    alv = log(v);
    if(alv < ynorm-amaxp) goto S170;
    if(alv > ynorm+amaxp) goto S30;
/*
     STIRLING'S FORMULA TO MACHINE ACCURACY FOR
     THE FINAL ACCEPTANCE/REJECTION TEST
*/
    x1 = ix+1.0;
    f1 = fm+1.0;
    z = n+1.0-fm;
    w = n-ix+1.0;
    z2 = z*z;
    x2 = x1*x1;
    f2 = f1*f1;
    w2 = w*w;
    if(alv <= xm*log(f1/x1)+(n-m+0.5)*log(z/w)+(ix-m)*log(w*p/(x1*q))+(13860.0-
      (462.0-(132.0-(99.0-140.0/f2)/f2)/f2)/f2)/f1/166320.0+(13860.0-(462.0-
      (132.0-(99.0-140.0/z2)/z2)/z2)/z2)/z/166320.0+(13860.0-(462.0-(132.0-
      (99.0-140.0/x2)/x2)/x2)/x2)/x1/166320.0+(13860.0-(462.0-(132.0-(99.0
      -140.0/w2)/w2)/w2)/w2)/w/166320.0) goto S170;
    goto S30;
S140:
/*
     INVERSE CDF LOGIC FOR MEAN LESS THAN 30
*/
    qn = pow(q,(double)n);
    r = p/q;
    g = r*(n+1);
S150:
    ix = 0;
    f = qn;
    u = ranf();
S160:
    if(u < f) goto S170;
    if(ix > 110) goto S150;
    u -= f;
    ix += 1;
    f *= (g/ix-r);
    goto S160;
S170:
    if(psave > 0.5) ix = n-ix;
    ignbin = ix;
    return ignbin;
}


long ignpoi_mt(double mu,int tn)
/*
**********************************************************************
long ignpoi_mt(double mu)
GENerate POIsson random deviate
Function
Generates a single random deviate from a Poisson
distribution with mean MU.
Arguments
mu --> The mean of the Poisson distribution from which
a random deviate is to be generated.
(mu >= 0.0)
ignpoi_mt <-- The random deviate.
Method
Renames KPOIS from TOMS as slightly modified by BWB to use RANF
instead of SUNIF.
For details see:
Ahrens, J.H. and Dieter, U.
Computer Generation of Poisson Deviates
From Modified Normal Distributions.
ACM Trans. Math. Software, 8, 2
(June 1982),163-179
**********************************************************************
**********************************************************************


P O I S S O NPat  DISTRIBUTION                                      


**********************************************************************
**********************************************************************

FOR DETAILS SEE:                                                 

AHRENS, J.H. AND DIETER, U.                            
COMPUTER GENERATION OF POISSON DEVIATES                
FROM MODIFIED NORMAL DISTRIBUTIONS.                    
ACM TRANS. MATH. SOFTWARE, 8,2 (JUNE 1982), 163 - 179. 

(SLIGHTLY MODIFIED VERSION OF THE PROGRAM IN THE ABOVE ARTICLE)  

**********************************************************************
INTEGER FUNCTION IGNPOI(IR,MU)
INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR
MU=MEAN MU OF THE POISSON DISTRIBUTION
OUTPUT: IGNPOI=SAMPLE FROM THE POISSON-(MU)-DISTRIBUTION
MUPREV=PREVIOUS MU, MUOLD=MU AT LAST EXECUTION OF STEP P OR B.
TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL
SEPARATION OF CASES A AND B
*/
	{
	extern double fsign( double num, double sign );
	static double a0 = -0.5;
	static double a1 = 0.3333333;
	static double a2 = -0.2500068;
	static double a3 = 0.2000118;
	static double a4 = -0.1661269;
	static double a5 = 0.1421878;
	static double a6 = -0.1384794;
	static double a7 = 0.125006;
	/* JJV changed the initial values of MUPREV and MUOLD */
	static double fact[10] = {
		1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0
		};
	/* JJV added ll to the list, for Case A */
	long ignpoi_mt,j,k,kflag,l,ll,m;
	double b1,b2,c,c0,c1,c2,c3,d,del,difmuk,e,fk,fx,fy,g,omega,p,p0,px,py,q,s,
		t,u,v,x,xx,pp[35];

	if(mu < 10.0) goto S120;
	/*
	NumC A S E  A. (RECALCULATION OF S,D,LL IF MU HAS CHANGED)
	JJV changed l in Case A to ll
	*/
	s = sqrt(mu);
	d = 6.0*mu*mu;
	/*
	THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
	PROBABILITIES FK WHENEVER K >= M(MU). LL=IFIX(MU-1.1484)
	IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .
	*/
	ll = (long) (mu-1.1484);
//S10:
	/*
	STEP NPat. NORMAL SAMPLE - SNORM(IR) FOR STANDARD NORMAL DEVIATE
	*/
	g = mu+s*snorm_mt(tn);
	if(g < 0.0) goto S20;
	ignpoi_mt = (long) (g);
  /*
	STEP I. IMMEDIATE ACCEPTANCE IF IGNPOI IS LARGE ENOUGH
	*/
	if(ignpoi_mt >= ll) return ignpoi_mt;
	/*
	STEP S. SQUEEZE ACCEPTANCE - SUNIF(IR) FOR (0,1)-SAMPLE U
	*/
	fk = (double)ignpoi_mt;
	difmuk = mu-fk;
	u = ranf_mt(tn);
	if(d*u >= difmuk*difmuk*difmuk) return ignpoi_mt;
S20:
	/*
	STEP P. PREPARATIONS FOR STEPS Q AND H.
	(RECALCULATIONS OF PARAMETERS IF NECESSARY)
	.3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
	THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
	APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
	NumC=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.
	*/
	omega = 0.3989423/s;
	b1 = 4.166667E-2/mu;
	b2 = 0.3*b1*b1;
	c3 = 0.1428571*b1*b2;
	c2 = b2-15.0*c3;
	c1 = b1-6.0*b2+45.0*c3;
	c0 = 1.0-b1+3.0*b2-15.0*c3;
	c = 0.1069/mu;
//S30:
	if(g < 0.0) goto S50;
	/*
	'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)
	*/
	kflag = 0;
	goto S70;
S40:
	/*
	STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)
	*/
	if(fy-u*fy <= py*exp(px-fx)) return ignpoi_mt;
S50:
	/*
	STEP E. EXPONENTIAL SAMPLE - SEXPO(IR) FOR STANDARD EXPONENTIAL
	DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
	(IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)
	*/
	e = sexpo_mt(tn);


	u = ranf_mt(tn);
	u += (u-1.0);
	t = 1.8+fsign(e,u);
	if(t <= -0.6744) goto S50;
	ignpoi_mt = (long) (mu+s*t);
  fk = (double)ignpoi_mt;
	difmuk = mu-fk;
	/*
	'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)
	*/
	kflag = 1;
	goto S70;
S60:
	/*
	STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)
	*/
	if(c*fabs(u) > py*exp(px+e)-fy*exp(fx+e)) goto S50;
	return ignpoi_mt;
S70:
	/*
	STEP F. 'SUBROUTINE' F. CALCULATION OF PX,PY,FX,FY.
	CASE IGNPOI .LT. 10 USES FACTORIALS FROM TABLE FACT
	*/
	if(ignpoi_mt >= 10) goto S80;
	px = -mu;
	py = pow(mu,(double)ignpoi_mt)/ *(fact+ignpoi_mt);
	goto S110;
S80:
	/*
	CASE IGNPOI .GE. 10 USES POLYNOMIAL APPROXIMATION
	A0-A7 FOR ACCURACY WHEN ADVISABLE
	.8333333E-1=1./12.  .3989423=(2*PI)**(-.5)
	*/
	del = 8.333333E-2/fk;
	del -= (4.8*del*del*del);
	v = difmuk/fk;
	if(fabs(v) <= 0.25) goto S90;
	px = fk*log(1.0+v)-difmuk-del;
	goto S100;
S90:
	px = fk*v*v*(((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0)-del;
S100:
	py = 0.3989423/sqrt(fk);
S110:
	x = (0.5-difmuk)/s;
	xx = x*x;
	fx = -0.5*xx;
	fy = omega*(((c3*xx+c2)*xx+c1)*xx+c0);
	if(kflag <= 0) goto S40;
	goto S60;
S120:
	/*
	NumC A S E  B. (START NEW TABLE AND CALCULATE P0 IF NECESSARY)
	JJV changed MUPREV assignment to initial value
	*/
	m = maxF(1L,(long) (mu));
	l = 0;
	p = exp(-mu);
	q = p0 = p;
S130:
	/*
	STEP U. UNIFORM SAMPLE FOR INVERSION METHOD
	*/
	u = ranf_mt(tn);
	ignpoi_mt = 0;
	if(u <= p0) return ignpoi_mt;
	/*
	STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
	PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
	(0.458=PP(9) FOR MU=10)
	*/
	if(l == 0) goto S150;
	j = 1;
	if(u > 0.458) j = minF(l,m);
	for(k=j; k<=l; k++) {
		if(u <= *(pp+k-1)) goto S180;
		}
	if(l == 35) goto S130;
S150:
/*
STEP NumC. CREATION OF NEW POISSON PROBABILITIES P
AND THEIR CUMULATIVES Q=PP(K)
*/
	l += 1;
	for(k=l; k<=35; k++) {
		p = p*mu/(double)k;
		q += p;
		*(pp+k-1) = q;
		if(u <= q) goto S170;
		}
	l = 35;
	goto S130;
	S170:
	l = k;
	S180:
	ignpoi_mt = k;
	return ignpoi_mt;
}


long ignbin_mt(long n,double pp,int tn)
/*
**********************************************************************
long ignbin_mt(long n,double pp)
GENerate BINomial random deviate
Function
Generates a single random deviate from a binomial
distribution whose number of trials is NPat and whose
probability of an event in each trial is P.
Arguments
n  --> The number of trials in the binomial distribution
from which a random deviate is to be generated.
JJV (NPat >= 0)
pp --> The probability of an event in each trial of the
binomial distribution from which a random deviate
is to be generated.
JJV (0.0 <= PP <= 1.0)
ignbin <-- A random deviate yielding the number of events
from NPat independent trials, each of which has
a probability of event P.
Method
This is algorithm BTPE from:
Kachitvichyanukul, V. and Schmeiser, B. W.
Binomial Random Variate Generation.
Communications of the ACM, 31, 2
(February, 1988) 216.
**********************************************************************
SUBROUTINE BTPEC(NPat,PP,ISEED,JX)
BINOMIAL RANDOM VARIATE GENERATOR
MEAN .LT. 30 -- INVERSE CDF
MEAN .GE. 30 -- ALGORITHM BTPE:  ACCEPTANCE-REJECTION VIA
FOUR REGION COMPOSITION.  THE FOUR REGIONS ARE A TRIANGLE
(SYMMETRIC IN THE CENTER), A PAIR OF PARALLELOGRAMS (ABOVE
THE TRIANGLE), AND EXPONENTIAL LEFT AND RIGHT TAILS.
BTPE REFERS TO BINOMIAL-TRIANGLE-PARALLELOGRAM-EXPONENTIAL.
BTPEC REFERS TO BTPE AND "COMBINED."  THUS BTPE IS THE
RESEARCH AND BTPEC IS THE IMPLEMENTATION OF A COMPLETE
USABLE ALGORITHM.
REFERENCE:  VORATAS KACHITVICHYANUKUL AND BRUCE SCHMEISER,
"BINOMIAL RANDOM VARIATE GENERATION,"
COMMUNICATIONS OF THE ACM, FORTHCOMING
WRITTEN:  SEPTEMBER 1980.
LAST REVISED:  MAY 1985, JULY 1987
REQUIRED SUBPROGRAM:  RAND() -- A UNIFORM (0,1) RANDOM NUMBER
GENERATOR
ARGUMENTS
NPat : NUMBER OF BERNOULLI TRIALS            (INPUT)
PP : PROBABILITY OF SUCCESS IN EACH TRIAL (INPUT)
ISEED:  RANDOM NUMBER SEED                (INPUT AND OUTPUT)
JX:  RANDOMLY GENERATED OBSERVATION       (OUTPUT)
VARIABLES
PSAVE: VALUE OF PP FROM THE LAST CALL TO BTPEC
NSAVE: VALUE OF NPat FROM THE LAST CALL TO BTPEC
XNP:  VALUE OF THE MEAN FROM THE LAST CALL TO BTPEC
P: PROBABILITY USED IN THE GENERATION PHASE OF BTPEC
FFM: TEMPORARY VARIABLE EQUAL TO XNP + P
M:  INTEGER VALUE OF THE CURRENT MODE
FM:  FLOATING POINT VALUE OF THE CURRENT MODE
XNPQ: TEMPORARY VARIABLE USED IN SETUP AND SQUEEZING STEPS
P1:  AREA OF THE TRIANGLE
NumC:  HEIGHT OF THE PARALLELOGRAMS
XM:  CENTER OF THE TRIANGLE
XL:  LEFT END OF THE TRIANGLE
XR:  RIGHT END OF THE TRIANGLE
AL:  TEMPORARY VARIABLE
XLL:  RATE FOR THE LEFT EXPONENTIAL TAIL
XLR:  RATE FOR THE RIGHT EXPONENTIAL TAIL
P2:  AREA OF THE PARALLELOGRAMS
P3:  AREA OF THE LEFT EXPONENTIAL TAIL
P4:  AREA OF THE RIGHT EXPONENTIAL TAIL
U:  A U(0,P4) RANDOM VARIATE USED FIRST TO SELECT ONE OF THE
FOUR REGIONS AND THEN CONDITIONALLY TO GENERATE A VALUE
FROM THE REGION
V:  A U(0,1) RANDOM NUMBER USED TO GENERATE THE RANDOM VALUE
(REGION 1) OR TRANSFORMED INTO THE VARIATE TO ACCEPT OR
REJECT THE CANDIDATE VALUE
IX:  INTEGER CANDIDATE VALUE
X:  PRELIMINARY CONTINUOUS CANDIDATE VALUE IN REGION 2 LOGIC
AND A FLOATING POINT IX IN THE ACCEPT/REJECT LOGIC
K:  ABSOLUTE VALUE OF (IX-M)
F:  THE HEIGHT OF THE SCALED DENSITY FUNCTION USED IN THE
ACCEPT/REJECT DECISION WHEN BOTH M AND IX ARE SMALL
ALSO USED IN THE INVERSE TRANSFORMATION
R: THE RATIO P/Q
G: CONSTANT USED IN CALCULATION OF PROBABILITY
MP:  MODE PLUS ONE, THE LOWER INDEX FOR EXPLICIT CALCULATION
OF F WHEN IX IS GREATER THAN M
IX1:  CANDIDATE VALUE PLUS ONE, THE LOWER INDEX FOR EXPLICIT
CALCULATION OF F WHEN IX IS LESS THAN M
I:  INDEX FOR EXPLICIT CALCULATION OF F FOR BTPE
AMAXP: MAXIMUM ERROR OF THE LOGARITHM OF NORMAL BOUND
YNORM: LOGARITHM OF NORMAL BOUND
ALV:  NATURAL LOGARITHM OF THE ACCEPT/REJECT VARIATE V
X1,F1,Z,W,Z2,X2,F2, AND W2 ARE TEMPORARY VARIABLES TO BE
USED IN THE FINAL ACCEPT/REJECT TEST
QN: PROBABILITY OF NO SUCCESS IN NPat TRIALS
REMARK
IX AND JX COULD LOGICALLY BE THE SAME VARIABLE, WHICH WOULD
SAVE A MEMORY POSITION AND A LINE OF CODE.  HOWEVER, SOME
COMPILERS (E.G.,CDC MNF) OPTIMIZE BETTER WHEN THE ARGUMENTS
ARE NOT INVOLVED.
ISEED NEEDS TO BE DOUBLE PRECISION IF THE IMSL ROUTINE
GGUBFS IS USED TO GENERATE UNIFORM RANDOM NUMBER, OTHERWISE
TYPE OF ISEED SHOULD BE DICTATED BY THE UNIFORM GENERATOR
**********************************************************************
*****DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY
*/
	{
	/* JJV changed initial values to ridiculous values */
	double psave = -1.0E37;
	long nsave = -214748365;
	long ignbin_mt,i,ix,ix1,k,m,mp,T1;
	double al,alv,amaxp,c,f,f1,f2,ffm,fm,g,p,p1,p2,p3,p4,q,qn,r,u,v,w,w2,x,x1,
		x2,xl,xll,xlr,xm,xnp,xnpq,xr,ynorm,z,z2;

//S10:
	/*
	*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE
	JJV added checks to ensure 0.0 <= PP <= 1.0
	*/
	if(pp < 0.0F) ERR_CRITICAL("PP < 0.0 in IGNBIN");
	if(pp > 1.0F) ERR_CRITICAL("PP > 1.0 in IGNBIN");
	psave = pp;
	p = minF(psave,1.0-psave);
	q = 1.0-p;
//S20:
	/*
	JJV added check to ensure NPat >= 0
	*/
	if(n < 0L) 
    ERR_CRITICAL("NPat < 0 in IGNBIN");
	xnp = n*p;
	nsave = n;
	if(xnp < 30.0) goto S140;
	ffm = xnp+p;
	m = (long) ffm;
	fm = m;
	xnpq = xnp*q;
	p1 = (long) (2.195*sqrt(xnpq)-4.6*q)+0.5;
	xm = fm+0.5;
	xl = xm-p1;
	xr = xm+p1;
	c = 0.134+20.5/(15.3+fm);
	al = (ffm-xl)/(ffm-xl*p);
	xll = al*(1.0+0.5*al);
	al = (xr-ffm)/(xr*q);
	xlr = al*(1.0+0.5*al);
	p2 = p1*(1.0+c+c);
	p3 = p2+c/xll;
	p4 = p3+c/xlr;
S30:
	/*
	*****GENERATE VARIATE
	*/
	u = ranf_mt(tn)*p4;
	v = ranf_mt(tn);
	/*
	TRIANGULAR REGION
	*/
	if(u > p1) goto S40;
	ix = (long) (xm-p1*v+u);
	goto S170;
S40:
	/*
	PARALLELOGRAM REGION
	*/
	if(u > p2) goto S50;
	x = xl+(u-p1)/c;
	v = v*c+1.0-ABS(xm-x)/p1;
	if(v > 1.0 || v <= 0.0) goto S30;
	ix = (long) x;
	goto S70;
S50:
	/*
	LEFT TAIL
	*/
	if(u > p3) goto S60;
	ix = (long) (xl+log(v)/xll);
	if(ix < 0) goto S30;
	v *= ((u-p2)*xll);
	goto S70;
S60:
	/*
	RIGHT TAIL
	*/
	ix = (long) (xr-log(v)/xlr);
	if(ix > n) goto S30;
	v *= ((u-p3)*xlr);
S70:
	/*
	*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST
	*/
	k = ABS(ix-m);
	if(k > 20 && k < xnpq/2-1) goto S130;
	/*
	EXPLICIT EVALUATION
	*/
	f = 1.0;
	r = p/q;
	g = (n+1)*r;
	T1 = m-ix;
	if(T1 < 0) goto S80;
	else if(T1 == 0) goto S120;
	else  goto S100;
S80:
	mp = m+1;
	for(i=mp; i<=ix; i++) f *= (g/i-r);
	goto S120;
S100:
	ix1 = ix+1;
	for(i=ix1; i<=m; i++) f /= (g/i-r);
S120:
	if(v <= f) goto S170;
	goto S30;
S130:
	/*
	SQUEEZING USING UPPER AND LOWER BOUNDS ON ALOG(F(X))
	*/
	amaxp = k/xnpq*((k*(k/3.0+0.625)+0.1666666666666)/xnpq+0.5);
	ynorm = -(k*k/(2.0*xnpq));
	alv = log(v);
	if(alv < ynorm-amaxp) goto S170;
	if(alv > ynorm+amaxp) goto S30;
	/*
	STIRLING'S FORMULA TO MACHINE ACCURACY FOR
	THE FINAL ACCEPTANCE/REJECTION TEST
	*/
	x1 = ix+1.0;
	f1 = fm+1.0;
	z = n+1.0-fm;
	w = n-ix+1.0;
	z2 = z*z;
	x2 = x1*x1;
	f2 = f1*f1;
	w2 = w*w;
	if(alv <= xm*log(f1/x1)+(n-m+0.5)*log(z/w)+(ix-m)*log(w*p/(x1*q))+(13860.0-
		(462.0-(132.0-(99.0-140.0/f2)/f2)/f2)/f2)/f1/166320.0+(13860.0-(462.0-
		(132.0-(99.0-140.0/z2)/z2)/z2)/z2)/z/166320.0+(13860.0-(462.0-(132.0-
		(99.0-140.0/x2)/x2)/x2)/x2)/x1/166320.0+(13860.0-(462.0-(132.0-(99.0
		-140.0/w2)/w2)/w2)/w2)/w/166320.0) goto S170;
	goto S30;
S140:
	/*
	INVERSE CDF LOGIC FOR MEAN LESS THAN 30
	*/
	qn = pow(q,(double)n);
	r = p/q;
	g = r*(n+1);
S150:
	ix = 0;
	f = qn;
	u = ranf_mt(tn);
S160:
	if(u < f) goto S170;
	if(ix > 110) goto S150;
	u -= f;
	ix += 1;
	f *= (g/ix-r);
	goto S160;
S170:
	if(psave > 0.5) ix = n-ix;
	ignbin_mt = ix;
	return ignbin_mt;
}



double sexpo(void)
/*
**********************************************************************
                                                                      
                                                                      
     (STANDARD-)  E X P O NPat E NPat T I A L   DISTRIBUTION                
                                                                      
                                                                      
**********************************************************************
**********************************************************************
                                                                      
     FOR DETAILS SEE:                                                 
                                                                      
               AHRENS, J.H. AND DIETER, U.                            
               COMPUTER METHODS FOR SAMPLING FROM THE                 
               EXPONENTIAL AND NORMAL DISTRIBUTIONS.                  
               COMM. ACM, 15,10 (OCT. 1972), 873 - 882.               
                                                                      
     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       
     'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       
                                                                      
     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
     SUNIF.  The argument IR thus goes away.                          
                                                                      
**********************************************************************
     Q(NPat) = SUM(ALOG(2.0)**K/K!)    K=1,..,NPat ,      THE HIGHEST NPat
     (HERE 8) IS DETERMINED BY Q(NPat)=1.0 WITHIN STANDARD PRECISION
*/
{
static double q[8] = {
    0.6931472,0.9333737,0.9888778,0.9984959,0.9998293,0.9999833,0.9999986,
    .9999999
};
long i;
double sexpo,a,u,ustar,umin;

	a = 0.0;
    u = ranf();
    goto S30;
S20:
    a += q[0];
S30:
    u += u;
/*
 * JJV changed the following to reflect the true algorithm and prevent
 * JJV unpredictable behavior if U is initially 0.5.
 *  if(u <= 1.0) goto S20;
 */
    if(u < 1.0) goto S20;
    u -= 1.0;
    if(u > q[0]) goto S60;
    sexpo = a+u;
    return sexpo;
S60:
    i = 1;
    ustar = ranf();
    umin = ustar;
S70:
    ustar = ranf();
    if(ustar < umin) umin = ustar;
    i += 1;
    if(u > q[i-1]) goto S70;
    return  a+umin*q[0];
}

double sexpo_mt(int tn)
/*
**********************************************************************


(STANDARD-)  E X P O NPat E NPat T I A L   DISTRIBUTION                


**********************************************************************
**********************************************************************

FOR DETAILS SEE:                                                 

AHRENS, J.H. AND DIETER, U.                            
COMPUTER METHODS FOR SAMPLING FROM THE                 
EXPONENTIAL AND NORMAL DISTRIBUTIONS.                  
COMM. ACM, 15,10 (OCT. 1972), 873 - 882.               

ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       
'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       

Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
SUNIF.  The argument IR thus goes away.                          

**********************************************************************
Q(NPat) = SUM(ALOG(2.0)**K/K!)    K=1,..,NPat ,      THE HIGHEST NPat
(HERE 8) IS DETERMINED BY Q(NPat)=1.0 WITHIN STANDARD PRECISION
*/
	{
	static double q[8] = {
		0.6931472,0.9333737,0.9888778,0.9984959,0.9998293,0.9999833,0.9999986,
			.9999999
		};
	long i;
	double sexpo_mt,a,u,ustar,umin;

	a = 0.0;
	u = ranf_mt(tn);
	goto S30;
S20:
	a += q[0];
S30:
	u += u;
	/*
	* JJV changed the following to reflect the true algorithm and prevent
	* JJV unpredictable behavior if U is initially 0.5.
	*  if(u <= 1.0) goto S20;
	*/
	if(u < 1.0) goto S20;
	u -= 1.0;
	if(u > q[0]) goto S60;
	sexpo_mt = a+u;

  if (sexpo_mt<0) {
    printf("sexpo_mt<0!");
  }

	return sexpo_mt;
S60:
	i = 1;
	ustar = ranf_mt(tn);
	umin = ustar;
S70:
	ustar = ranf_mt(tn);
	if(ustar < umin) umin = ustar;
	i += 1;
	if(u > q[i-1]) goto S70;
  
  if (a+umin*q[0]<0) {
    printf("a+umin*q<0!");
  }

	return  a+umin*q[0];
}

double genexp(double av)
/*
**********************************************************************
     float genexp(float av)
                    GENerate EXPonential random deviate
                              Function
     Generates a single random deviate from an exponential
     distribution with mean AV.
                              Arguments
     av --> The mean of the exponential distribution from which
            a random deviate is to be generated.
        JJV (av >= 0)
                              Method
     Renames SEXPO from TOMS as slightly modified by BWB to use RANF
     instead of SUNIF.
     For details see:
               Ahrens, J.H. and Dieter, U.
               Computer Methods for Sampling From the
               Exponential and Normal Distributions.
               Comm. ACM, 15,10 (Oct. 1972), 873 - 882.
**********************************************************************
*/
{
double genexp;

/* JJV added check that av >= 0 */
    if(av >= 0.0) goto S10;
    fputs(" AV < 0 in GENEXP - ABORT\n",stderr);
    fprintf(stderr," Value of AV: %16.6E\n",av);
    exit(1);
S10:
    genexp = sexpo()*av;
    return genexp;
}

double genexp_mt(int thread_no, double av)
/* Multi threaded version (Wes)
**********************************************************************
     float genexp(float av)
                    GENerate EXPonential random deviate
                              Function
     Generates a single random deviate from an exponential
     distribution with mean AV.
                              Arguments
     av --> The mean of the exponential distribution from which
            a random deviate is to be generated.
        JJV (av >= 0)
                              Method
     Renames SEXPO from TOMS as slightly modified by BWB to use RANF
     instead of SUNIF.
     For details see:
               Ahrens, J.H. and Dieter, U.
               Computer Methods for Sampling From the
               Exponential and Normal Distributions.
               Comm. ACM, 15,10 (Oct. 1972), 873 - 882.
**********************************************************************
*/
{
double genexp;

/* JJV added check that av >= 0 */
    if(av >= 0.0) goto S10;
    fputs(" AV < 0 in GENEXP - ABORT\n",stderr);
    fprintf(stderr," Value of AV: %16.6E\n",av);
    exit(1);
S10:
    genexp = sexpo_mt(thread_no)*av;
    return genexp;
}


double snorm(void)
/*
**********************************************************************
                                                                      
                                                                      
     (STANDARD-)  NPat O R M A L  DISTRIBUTION                           
                                                                      
                                                                      
**********************************************************************
**********************************************************************
                                                                      
     FOR DETAILS SEE:                                                 
                                                                      
               AHRENS, J.H. AND DIETER, U.                            
               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             
               SAMPLING FROM THE NORMAL DISTRIBUTION.                 
               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          
                                                                      
     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  
     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  
                                                                      
     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
     SUNIF.  The argument IR thus goes away.                          
                                                                      
**********************************************************************
     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
*/
{
static double a[32] = {
    0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
    0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
    0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
    1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
    1.862732,2.153875
};
static double d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
    0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
    0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
    0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
};
static double t[31] = {
    7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
    1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
    2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
    4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
    9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
};
static double h[31] = {
    3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
    4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
    4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
    5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
    8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
};
static long i;
static double snorm,u,s,ustar,aa,w,y,tt;
    u = ranf();
    s = 0.0;
    if(u > 0.5) s = 1.0;
    u += (u-s);
    u = 32.0*u;
    i = (long) (u);
    if(i == 32) i = 31;
    if(i == 0) goto S100;
/*
                                START CENTER
*/
    ustar = u-(double)i;
    aa = *(a+i-1);
S40:
    if(ustar <= *(t+i-1)) goto S60;
    w = (ustar-*(t+i-1))**(h+i-1);
S50:
/*
                                EXIT   (BOTH CASES)
*/
    y = aa+w;
    snorm = y;
    if(s == 1.0) snorm = -y;
    return snorm;
S60:
/*
                                CENTER CONTINUED
*/
    u = ranf();
    w = u*(*(a+i)-aa);
    tt = (0.5*w+aa)*w;
    goto S80;
S70:
    tt = u;
    ustar = ranf();
S80:
    if(ustar > tt) goto S50;
    u = ranf();
    if(ustar >= u) goto S70;
    ustar = ranf();
    goto S40;
S100:
/*
                                START TAIL
*/
    i = 6;
    aa = *(a+31);
    goto S120;
S110:
    aa += *(d+i-1);
    i += 1;
S120:
    u += u;
    if(u < 1.0) goto S110;
    u -= 1.0;
S140:
    w = u**(d+i-1);
    tt = (0.5*w+aa)*w;
    goto S160;
S150:
    tt = u;
S160:
    ustar = ranf();
    if(ustar > tt) goto S50;
    u = ranf();
    if(ustar >= u) goto S150;
    u = ranf();
    goto S140;
}

double snorm_mt(int tn)
/*
**********************************************************************


(STANDARD-)  NPat O R M A L  DISTRIBUTION                           


**********************************************************************
**********************************************************************

FOR DETAILS SEE:                                                 

AHRENS, J.H. AND DIETER, U.                            
EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             
SAMPLING FROM THE NORMAL DISTRIBUTION.                 
MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          

ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  
(M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  

Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
SUNIF.  The argument IR thus goes away.                          

**********************************************************************
THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
*/
	{
	static double a[32] = {
		0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
			0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
			0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
			1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
			1.862732,2.153875
		};
	static double d[31] = {
		0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
			0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
			0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
			0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
		};
	static double t[31] = {
		7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
			1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
			2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
			4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
			9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
		};
	static double h[31] = {
		3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
			4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
			4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
			5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
			8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
		};
	long i;
	double snorm_mt,u,s,ustar,aa,w,y,tt;
	
  u = ranf_mt(tn);
	s = 0.0;
	if(u > 0.5) s = 1.0;
	u += (u-s);
	u = 32.0*u;
	i = (long) (u);
	if(i == 32) i = 31;
	if(i == 0) goto S100;
	/*
	START CENTER
	*/
	ustar = u-(double)i;
	aa = *(a+i-1);
S40:
	if(ustar <= *(t+i-1)) goto S60;
	w = (ustar-*(t+i-1))**(h+i-1);
S50:
	/*
	EXIT   (BOTH CASES)
	*/
	y = aa+w;
	snorm_mt = y;
	if(s == 1.0) snorm_mt = -y;
  return snorm_mt;
S60:
	/*
	CENTER CONTINUED
	*/
	u = ranf_mt(tn);
	w = u*(*(a+i)-aa);
	tt = (0.5*w+aa)*w;
	goto S80;
S70:
	tt = u;
	ustar = ranf_mt(tn);
S80:
	if(ustar > tt) goto S50;
	u = ranf_mt(tn);
	if(ustar >= u) goto S70;
	ustar = ranf_mt(tn);
	goto S40;
S100:
	/*
	START TAIL
	*/
	i = 6;
	aa = *(a+31);
	goto S120;
S110:
	aa += *(d+i-1);
	i += 1;
S120:
	u += u;
	if(u < 1.0) goto S110;
	u -= 1.0;
S140:
	w = u**(d+i-1);
	tt = (0.5*w+aa)*w;
	goto S160;
S150:
	tt = u;
S160:
	ustar = ranf_mt(tn);
	if(ustar > tt) goto S50;
	u = ranf_mt(tn);
	if(ustar >= u) goto S150;
	u = ranf_mt(tn);
	goto S140;
}
double fsign( double num, double sign )
/* Transfers sign of argument sign to argument num */
{
if ( ( sign>0.0f && num<0.0f ) || ( sign<0.0f && num>0.0f ) )
    return -num;
else return num;
}

double gengam(double a,double r)
/*
**********************************************************************
double gengam(double a,double r)
GENerates random deviates from GAMma distribution
Function
Generates random deviates from the gamma distribution whose
density is
(A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)
Arguments
a --> Location parameter of Gamma distribution
JJV   (a > 0)
r --> Shape parameter of Gamma distribution
JJV   (r > 0)
Method
Renames SGAMMA from TOMS as slightly modified by BWB to use RANF
instead of SUNIF.
For details see:
(Case R >= 1.0)
Ahrens, J.H. and Dieter, U.
Generating Gamma Variates by a
Modified Rejection Technique.
Comm. ACM, 25,1 (Jan. 1982), 47 - 54.
Algorithm GD
JJV altered following to reflect argument ranges
(Case 0.0 < R < 1.0)
Ahrens, J.H. and Dieter, U.
Computer Methods for Sampling from Gamma,
Beta, Poisson and Binomial Distributions.
Computing, 12 (1974), 223-246/
Adapted algorithm GS.
**********************************************************************
*/
{
	static double gengam;
	/* JJV added argument checker */
	if(a > 0.0 && r > 0.0) goto S10;
	fputs(" A or R nonpositive in GENGAM - abort!\n",stderr);
	fprintf(stderr," A value: %16.6E R value: %16.6E\n",a,r);
	exit(1);
S10:
	gengam = sgamma(r);
	gengam /= a;
	return gengam;
}

//double sgamma(double a)
///*
//**********************************************************************
//
//
//(STANDARD-)  G A M M A  DISTRIBUTION                             
//
//
//**********************************************************************
//**********************************************************************
//
//PARAMETER  A >= 1.0  !                                 
//
//**********************************************************************
//
//FOR DETAILS SEE:                                                 
//
//AHRENS, J.H. AND DIETER, U.                            
//GENERATING GAMMA VARIATES BY A                         
//MODIFIED REJECTION TECHNIQUE.                          
//COMM. ACM, 25,1 (JAN. 1982), 47 - 54.                  
//
//STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER     
//(STRAIGHTFORWARD IMPLEMENTATION)     
//
//Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
//SUNIF.  The argument IR thus goes away.                          
//
//**********************************************************************
//
//PARAMETER  0.0 < A < 1.0  !                            
//
//**********************************************************************
//
//FOR DETAILS SEE:                                                 
//
//AHRENS, J.H. AND DIETER, U.                            
//COMPUTER METHODS FOR SAMPLING FROM GAMMA,              
//BETA, POISSON AND BINOMIAL DISTRIBUTIONS.              
//COMPUTING, 12 (1974), 223 - 246.                       
//
//(ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)    
//
//**********************************************************************
//INPUT: A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
//OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
//COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
//COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
//COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
//PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
//SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
//*/
//{
//	extern double fsign( double num, double sign );
//	static double q1 = 4.166669E-2;
//	static double q2 = 2.083148E-2;
//	static double q3 = 8.01191E-3;
//	static double q4 = 1.44121E-3;
//	static double q5 = -7.388E-5;
//	static double q6 = 2.4511E-4;
//	static double q7 = 2.424E-4;
//	static double a1 = 0.3333333;
//	static double a2 = -0.250003;
//	static double a3 = 0.2000062;
//	static double a4 = -0.1662921;
//	static double a5 = 0.1423657;
//	static double a6 = -0.1367177;
//	static double a7 = 0.1233795;
//	static double e1 = 1.0;
//	static double e2 = 0.4999897;
//	static double e3 = 0.166829;
//	static double e4 = 4.07753E-2;
//	static double e5 = 1.0293E-2;
//	static double aa = 0.0;
//	static double aaa = 0.0;
//	static double sqrt32 = 5.656854;
//	/* JJV added b0 to fix rare and subtle bug */
//	static double sgamma,s2,s,d,t,x,u,r,q0,b,b0,si,c,v,q,e,w,p;
//	if(a == aa) goto S10;
//	if(a < 1.0) goto S120;
//	/*
//	STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
//	*/
//	aa = a;
//	s2 = a-0.5;
//	s = sqrt(s2);
//	d = sqrt32-12.0*s;
//S10:
//	/*
//	STEP  2:  T=STANDARD NORMAL DEVIATE,
//	X=(S,1/2)-NORMAL DEVIATE.
//	IMMEDIATE ACCEPTANCE (I)
//	*/
//	t = snorm();
//	x = s+0.5*t;
//	sgamma = x*x;
//	if(t >= 0.0) return sgamma;
//	/*
//	STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
//	*/
//	u = ranf();
//	if(d*u <= t*t*t) return sgamma;
//	/*
//	STEP  4:  RECALCULATIONS OF Q0,B,SI,NumC IF NECESSARY
//	*/
//	if(a == aaa) goto S40;
//	aaa = a;
//	r = 1.0/ a;
//	q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r;
//	/*
//	APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
//	THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
//	NumC WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
//	*/
//	if(a <= 3.686) goto S30;
//	if(a <= 13.022) goto S20;
//	/*
//	CASE 3:  A .GT. 13.022
//	*/
//	b = 1.77;
//	si = 0.75;
//	c = 0.1515/s;
//	goto S40;
//S20:
//	/*
//	CASE 2:  3.686 .LT. A .LE. 13.022
//	*/
//	b = 1.654+7.6E-3*s2;
//	si = 1.68/s+0.275;
//	c = 6.2E-2/s+2.4E-2;
//	goto S40;
//S30:
//	/*
//	CASE 1:  A .LE. 3.686
//	*/
//	b = 0.463+s+0.178*s2;
//	si = 1.235;
//	c = 0.195/s-7.9E-2+1.6E-1*s;
//S40:
//	/*
//	STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
//	*/
//	if(x <= 0.0) goto S70;
//	/*
//	STEP  6:  CALCULATION OF V AND QUOTIENT Q
//	*/
//	v = t/(s+s);
//	if(fabs(v) <= 0.25) goto S50;
//	q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
//	goto S60;
//S50:
//	q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
//S60:
//	/*
//	STEP  7:  QUOTIENT ACCEPTANCE (Q)
//	*/
//	if(log(1.0-u) <= q) return sgamma;
//S70:
//	/*
//	STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
//	U= 0,1 -UNIFORM DEVIATE
//	T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
//	*/
//	e = sexpo();
//	u = ranf();
//	u += (u-1.0);
//	t = b+fsign(si*e,u);
//	/*
//	STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
//	*/
//	if(t < -0.7187449) goto S70;
//	/*
//	STEP 10:  CALCULATION OF V AND QUOTIENT Q
//	*/
//	v = t/(s+s);
//	if(fabs(v) <= 0.25) goto S80;
//	q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
//	goto S90;
//S80:
//	q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
//S90:
//	/*
//	STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
//	*/
//	if(q <= 0.0) goto S70;
//	if(q <= 0.5) goto S100;
//	/*
//	* JJV modified the code through line 115 to handle large Q case
//	*/
//	if(q < 15.0) goto S95;
//	/*
//	* JJV Here Q is large enough that Q = log(exp(Q) - 1.0) (for real Q)
//	* JJV so reformulate test at 110 in terms of one EXP, if not too big
//	* JJV 87.49823 is close to the largest real which can be
//	* JJV exponentiated (87.49823 = log(1.0E38))
//	*/
//	if((q+e-0.5*t*t) > 87.49823) goto S115;
//	if(c*fabs(u) > exp(q+e-0.5*t*t)) goto S70;
//	goto S115;
//S95:
//	w = exp(q)-1.0;
//	goto S110;
//S100:
//	w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q;
//S110:
//	/*
//	IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
//	*/
//	if(c*fabs(u) > w*exp(e-0.5*t*t)) goto S70;
//S115:
//	x = s+0.5*t;
//	sgamma = x*x;
//	return sgamma;
//S120:
//	/*
//	ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))
//
//	JJV changed B to B0 (which was added to declarations for this)
//	JJV in 120 to END to fix rare and subtle bug.
//	JJV Line: 'aa = 0.0' was removed (unnecessary, wasteful).
//	JJV Reasons: the state of AA only serves to tell the A >= 1.0
//	JJV case if certain A-dependent constants need to be recalculated.
//	JJV The A < 1.0 case (here) no longer changes any of these, and
//	JJV the recalculation of B (which used to change with an
//	JJV A < 1.0 call) is governed by the state of AAA anyway.
//	aa = 0.0;
//	*/
//	b0 = 1.0+0.3678794*a;
//S130:
//	p = b0*ranf();
//	if(p >= 1.0) goto S140;
//	sgamma = exp(log(p)/ a);
//	if(sexpo() < sgamma) goto S130;
//	return sgamma;
//S140:
//	sgamma = -log((b0-p)/ a);
//	if(sexpo() < (1.0-a)*log(sgamma)) goto S130;
//	return sgamma;
//}

double gengam_mt(double a,double r,int tn)
/*
**********************************************************************
double gengam(double a,double r)
GENerates random deviates from GAMma distribution
Function
Generates random deviates from the gamma distribution whose
density is
(A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)
Arguments
a --> Location parameter of Gamma distribution
JJV   (a > 0)
r --> Shape parameter of Gamma distribution
JJV   (r > 0)
Method
Renames SGAMMA from TOMS as slightly modified by BWB to use RANF
instead of SUNIF.
For details see:
(Case R >= 1.0)
Ahrens, J.H. and Dieter, U.
Generating Gamma Variates by a
Modified Rejection Technique.
Comm. ACM, 25,1 (Jan. 1982), 47 - 54.
Algorithm GD
JJV altered following to reflect argument ranges
(Case 0.0 < R < 1.0)
Ahrens, J.H. and Dieter, U.
Computer Methods for Sampling from Gamma,
Beta, Poisson and Binomial Distributions.
Computing, 12 (1974), 223-246/
Adapted algorithm GS.
**********************************************************************
*/
	{
	double gengam_mt;
	/* JJV added argument checker */
	if(a > 0.0 && r > 0.0) goto S10;
	fputs(" A or R nonpositive in GENGAM - abort!\n",stderr);
	fprintf(stderr," A value: %16.6E R value: %16.6E\n",a,r);
	exit(1);
S10:
	gengam_mt = sgamma_mt(r,tn);
	gengam_mt /= a;
	return gengam_mt;
}

double sgamma_mt(double a,int tn)
/*
**********************************************************************


(STANDARD-)  G A M M A  DISTRIBUTION                             


**********************************************************************
**********************************************************************

PARAMETER  A >= 1.0  !                                 

**********************************************************************

FOR DETAILS SEE:                                                 

AHRENS, J.H. AND DIETER, U.                            
GENERATING GAMMA VARIATES BY A                         
MODIFIED REJECTION TECHNIQUE.                          
COMM. ACM, 25,1 (JAN. 1982), 47 - 54.                  

STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER     
(STRAIGHTFORWARD IMPLEMENTATION)     

Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
SUNIF.  The argument IR thus goes away.                          

**********************************************************************

PARAMETER  0.0 < A < 1.0  !                            

**********************************************************************

FOR DETAILS SEE:                                                 

AHRENS, J.H. AND DIETER, U.                            
COMPUTER METHODS FOR SAMPLING FROM GAMMA,              
BETA, POISSON AND BINOMIAL DISTRIBUTIONS.              
COMPUTING, 12 (1974), 223 - 246.                       

(ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)    

**********************************************************************
INPUT: A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
*/
	{
	extern double fsign( double num, double sign );
	static double q1 = 4.166669E-2;
	static double q2 = 2.083148E-2;
	static double q3 = 8.01191E-3;
	static double q4 = 1.44121E-3;
	static double q5 = -7.388E-5;
	static double q6 = 2.4511E-4;
	static double q7 = 2.424E-4;
	static double a1 = 0.3333333;
	static double a2 = -0.250003;
	static double a3 = 0.2000062;
	static double a4 = -0.1662921;
	static double a5 = 0.1423657;
	static double a6 = -0.1367177;
	static double a7 = 0.1233795;
	static double e1 = 1.0;
	static double e2 = 0.4999897;
	static double e3 = 0.166829;
	static double e4 = 4.07753E-2;
	static double e5 = 1.0293E-2;
	static double aa = 0.0;
	static double aaa = 0.0;
	static double sqrt32 = 5.656854;
	/* JJV added b0 to fix rare and subtle bug */
	static double sgamma_mt,s2,s,d,t,x,u,r,q0,b,b0,si,c,v,q,e,w,p;
	if(a == aa) goto S10;
	if(a < 1.0) goto S120;
	/*
	STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
	*/
	aa = a;
	s2 = a-0.5;
	s = sqrt(s2);
	d = sqrt32-12.0*s;
S10:
	/*
	STEP  2:  T=STANDARD NORMAL DEVIATE,
	X=(S,1/2)-NORMAL DEVIATE.
	IMMEDIATE ACCEPTANCE (I)
	*/
	t = snorm();
	x = s+0.5*t;
	sgamma_mt = x*x;
	if(t >= 0.0) return sgamma_mt;
	/*
	STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
	*/
	u = ranf_mt(tn);
	if(d*u <= t*t*t) return sgamma_mt;
	/*
	STEP  4:  RECALCULATIONS OF Q0,B,SI,NumC IF NECESSARY
	*/
	if(a == aaa) goto S40;
	aaa = a;
	r = 1.0/ a;
	q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r;
	/*
	APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
	THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
	NumC WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
	*/
	if(a <= 3.686) goto S30;
	if(a <= 13.022) goto S20;
	/*
	CASE 3:  A .GT. 13.022
	*/
	b = 1.77;
	si = 0.75;
	c = 0.1515/s;
	goto S40;
S20:
	/*
	CASE 2:  3.686 .LT. A .LE. 13.022
	*/
	b = 1.654+7.6E-3*s2;
	si = 1.68/s+0.275;
	c = 6.2E-2/s+2.4E-2;
	goto S40;
S30:
	/*
	CASE 1:  A .LE. 3.686
	*/
	b = 0.463+s+0.178*s2;
	si = 1.235;
	c = 0.195/s-7.9E-2+1.6E-1*s;
S40:
	/*
	STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
	*/
	if(x <= 0.0) goto S70;
	/*
	STEP  6:  CALCULATION OF V AND QUOTIENT Q
	*/
	v = t/(s+s);
	if(fabs(v) <= 0.25) goto S50;
	q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
	goto S60;
S50:
	q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
S60:
	/*
	STEP  7:  QUOTIENT ACCEPTANCE (Q)
	*/
	if(log(1.0-u) <= q) return sgamma_mt;
S70:
	/*
	STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
	U= 0,1 -UNIFORM DEVIATE
	T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
	*/
	e = sexpo_mt(tn);
	u = ranf_mt(tn);
	u += (u-1.0);
	t = b+fsign(si*e,u);
	/*
	STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
	*/
	if(t < -0.7187449) goto S70;
	/*
	STEP 10:  CALCULATION OF V AND QUOTIENT Q
	*/
	v = t/(s+s);
	if(fabs(v) <= 0.25) goto S80;
	q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
	goto S90;
S80:
	q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
S90:
	/*
	STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
	*/
	if(q <= 0.0) goto S70;
	if(q <= 0.5) goto S100;
	/*
	* JJV modified the code through line 115 to handle large Q case
	*/
	if(q < 15.0) goto S95;
	/*
	* JJV Here Q is large enough that Q = log(exp(Q) - 1.0) (for real Q)
	* JJV so reformulate test at 110 in terms of one EXP, if not too big
	* JJV 87.49823 is close to the largest real which can be
	* JJV exponentiated (87.49823 = log(1.0E38))
	*/
	if((q+e-0.5*t*t) > 87.49823) goto S115;
	if(c*fabs(u) > exp(q+e-0.5*t*t)) goto S70;
	goto S115;
S95:
	w = exp(q)-1.0;
	goto S110;
S100:
	w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q;
S110:
	/*
	IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
	*/
	if(c*fabs(u) > w*exp(e-0.5*t*t)) goto S70;
S115:
	x = s+0.5*t;
	sgamma_mt = x*x;
	return sgamma_mt;
S120:
	/*
	ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))

	JJV changed B to B0 (which was added to declarations for this)
	JJV in 120 to END to fix rare and subtle bug.
	JJV Line: 'aa = 0.0' was removed (unnecessary, wasteful).
	JJV Reasons: the state of AA only serves to tell the A >= 1.0
	JJV case if certain A-dependent constants need to be recalculated.
	JJV The A < 1.0 case (here) no longer changes any of these, and
	JJV the recalculation of B (which used to change with an
	JJV A < 1.0 call) is governed by the state of AAA anyway.
	aa = 0.0;
	*/
	b0 = 1.0+0.3678794*a;
S130:
	p = b0*ranf_mt(tn);
	if(p >= 1.0) goto S140;
	sgamma_mt = exp(log(p)/ a);
	if(sexpo_mt(tn) < sgamma_mt) goto S130;
	return sgamma_mt;
S140:
	sgamma_mt = -log((b0-p)/ a);
	if(sexpo_mt(tn) < (1.0-a)*log(sgamma_mt)) goto S130;
	return sgamma_mt;
	}
