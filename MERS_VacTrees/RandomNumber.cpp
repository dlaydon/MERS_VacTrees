#include "DataStructure.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


double gasdev(long *idum)
{
	double ran1(long *idum);
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if (*idum < 0) iset=0;
	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}


#include <math.h>
/* (C) Copr. 1986-92 Numerical Recipes Software Vn2-5. */
#define PI 3.1415926535897932384626433832795029

double bnldev(double pp, int n, long *idum)
{
	double gammln(double xx);
	double ran1(long *idum);
	int j;
	static int nold=(-1);
	double am,em,g,angle,p,bnl,sq,t,y;
	static double pold=(-1.0),pc,plog,pclog,en,oldg;

	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=1;j<=n;j++)
			if (ran1(idum) < p) ++bnl;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= ran1(idum);
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*ran1(idum);
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
				-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (ran1(idum) > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}
#undef PI

//------------------
#include <time.h>
double runif()
{
	return ran1(&SEED);
}

double rnorm()
{
	return gasdev(&SEED);
}

int rbin(int n, double p)
{
	return (int)bnldev(p,n,&SEED);
}

void rmultinomial(int n, double *proba,int Size,int* deviate)
{
	int i;
	double cumulativeSumProba = 0;
	int cumulativeSumN = 0;
	for (i = 0; i<Size; i++) deviate[i] = 0;
	for (i = 0; i<Size; i++)
	{
		deviate[i] = rbin(n - cumulativeSumN, proba[i] / (1 - cumulativeSumProba));
		cumulativeSumProba += proba[i];
		cumulativeSumN += deviate[i];
		if (cumulativeSumN >= n) break;
	}
}

int rmultinomialOneValue(std::vector<double> &proba, int size)
{
	int i;
	//int n=1;
	int deviate = 0;
	double cumulativeSumProba = 0.;
	//int cumulativeSumN=0;
	for (i = 0; i<size; i++)
	{
		deviate = rbin(1, proba[i] / (1 - cumulativeSumProba));
		cumulativeSumProba += proba[i];
		if (deviate == 1) break;
	}
	if (i == size) i = size - 1;
	return i;
}


int rmultinomialOneValue(std::vector<double> &proba, std::set<int> &possibleInfectors) 
{
    int deviate					= 0;
	double cumulativeSumProba	= 0.;
	set<int>::iterator pInfector;
	for (pInfector = possibleInfectors.begin(); pInfector != possibleInfectors.end(); ++pInfector) 
	{
		deviate				= rbin(1, proba[*pInfector] / (1. - cumulativeSumProba));
		cumulativeSumProba	+= proba[*pInfector];
		if (deviate == 1) break;
	}
	return *pInfector;
}

//------------------fonctions
#include <math.h>

double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
/* (C) Copr. 1986-92 Numerical Recipes Software Vn2-5. */

double gammq(double a, double x)
{// j'ai supprim� nerror
	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	double gamser,gammcf,gln;
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software Vn2-5. */

/* (C) Copr. 1986-92 Numerical Recipes Software Vn2-5. */

/* fonction incomplete gamma function (de 0 a x), tir�e de Numerical Recipies*/
double gammp(double a, double x)
{
	// j'ai supprim� nerror
	void gcf(double* gammcf, double a, double x, double* gln);
	void gser(double* gamser, double a, double x, double* gln);
	double gamser, gammcf, gln;
	if (x < (a + 1.0)) {
		gser(&gamser, a, x, &gln);
		return gamser;
	}
	else {
		gcf(&gammcf, a, x, &gln);
		return 1.0 - gammcf;
	}
}

double pgamma(double x,double a, double b)
{
	return gammp(a, x * b);
}

double dgamma(double x,double a,double b)
{
	return exp(a * log(b) - gammln(a) + (a - 1) * log(x) - b * x);
}

double ldgamma(double x,double a,double b)
{
	return a * log(b) - gammln(a) + (a - 1) * log(x) - b * x;
}

#include <math.h>
#define ITMAX 100
#define EPS 3.0e-7

void gser(double *gamser, double a, double x, double *gln)
{
	double gammln(double xx);
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		return;
	}
}
#undef ITMAX
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software Vn2-5. */

#include <math.h>
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(double *gammcf, double a, double x, double *gln)
{
	double gammln(double xx);
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN
/* (C) Copr. 1986-92 Numerical Recipes Software Vn2-5. */

double choose(int n,int k)
{return exp(gammln(n+1)-gammln(n-k+1)-gammln(k+1));
}


// POUR NEGATIVE BINOMIAL: fsign, sgamma, poidev

/* Transfers sign of argument sign to argument num */
static double fsign( double num, double sign )
{
  if ( ( sign>0.0f && num<0.0f ) || ( sign<0.0f && num>0.0f ) )
    return -num;
  else return num;
}

double sgamma(double a)
/*  SOURCE : http://www.math.iastate.edu/burkardt/c_src/ranlib/ranlib.c
**********************************************************************
     (STANDARD-)  G A M M A  DISTRIBUTION
**********************************************************************
*/
{static double q1 = 4.166669E-2;
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
static double sgamma,s2,s,d,t,x,u,r,q0,b,si,c,v,q,e,w,p;
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
    t = rnorm();  //*********modif perso: snorm() devient rnorm()**************
    x = s+0.5*t;
    sgamma = x*x;
    if(t >= 0.0) return sgamma;
/*
     STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
*/
    u = runif();     //************modif perso ranf() devient runif()*****************
    if(d*u <= t*t*t) return sgamma;
/*
     STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
*/
    if(a == aaa) goto S40;
    aaa = a;
    r = 1.0/ a;
    q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r;
/*
               APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
               THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
               C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
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
    if(log(1.0-u) <= q) return sgamma;
S70:
/*
     STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
               U= 0,1 -UNIFORM DEVIATE
               T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
*/
    e = -log(runif()); //modif perso: sexpo est generer directement
    u = runif();    //************modif perso ranf() devient runif()*****************
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
    w = exp(q)-1.0;
    goto S110;
S100:
    w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q;
S110:
/*
               IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
*/
    if(c*fabs(u) > w*exp(e-0.5*t*t)) goto S70;
    x = s+0.5*t;
    sgamma = x*x;
    return sgamma;
S120:
/*
     ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))
*/
    aa = 0.0;
    b = 1.0+0.3678794*a;
S130:
    p = b*runif();
    if(p >= 1.0) goto S140;
    sgamma = exp(log(p)/ a);
    if(-log(runif()) < sgamma) goto S130;  //modif perso: sexpo est generer directement
    return sgamma;
S140:
    sgamma = -log((b-p)/ a);
    if(-log(runif()) < (1.0-a)*log(sgamma)) goto S130; //modif perso: sexpo est generer directement
    return sgamma;
}


double rgamma(double a,double b)
{return sgamma(a)/b;}

double rgammaBigger0(double a,double b)
{double x=rgamma(a,b);
while(x==0) x=rgamma(a,b);
return x;
}

#include <math.h>
#define PI 3.141592654

// poisson deviate
int poidev(double xm)
{
	if(xm==0) return 0;
	if(xm>100000)	//normal approximation for xm v large
		{return int(xm+pow(xm,0.5)*rnorm());}

	static double sq,alxm,g,oldm=(-1.0);
	double em,t,y;

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= runif();
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*runif());
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (runif() > t);
	}
	return (long)em;
}
#undef PI

int rnegbin(double n,double p) //April 2016 - potential pb
//Modified on Jan 2016 to match definition used for dnegbin: replacing p by 1-p
// has density f(x)=Gamma(x+n)/x!Gamma(n)*p^n*(1-p)^x -- checked with following code
// f(x): proba it takes x successes before one get n failues. where p is the probability of failure
/*//check function rnegbin
double R=2;
double k=0.14;
double p=1./(1.+k/R);
double r=R*(1-p)/p;
int j;
cout<<"r="<<r<<";p="<<p<<"\n";
int nbSimul=100000;
double *Rval=new double[nbSimul];
double meanR=0;
double varR=0;
for(j=0;j<nbSimul;j++)
	{Rval[j]=rnegbin(r,1-p);
	meanR+=Rval[j];
	}
meanR=meanR/nbSimul;
for(j=0;j<nbSimul;j++)
	{varR+=pow(Rval[j]-meanR,2);
	}
varR=varR/nbSimul;
cout<<"meanR="<<meanR<<","<<"varR="<<varR<<"\n";
cout<<"mRexact="<<p*r/(1-p)<<" / "<<R<<"\n";
cout<<"vRexact="<<p*r/(1-p)/(1-p)<<" / "<<R*(1+R/k)<<"\n";
char q;
cin>>q;*/
{
static int ignnbn;
static double y,a,r;
/*
     Generate Y, a random gamma (n,(1-p)/p) variable
     JJV Note: the above parametrization is consistent with Devroye,
     JJV       but gamma (p/(1-p),n) is the equivalent in our code
*/
    r = (double)n;
    a = (1.0F-p)/p;//before change in Jan 2016: a = p/(1.0F-p);
/*
 * JJV changed this to call SGAMMA directly
 *  y = gengam(a,r); <- OLD
 */
    y = sgamma(r)/a;
/*
     Generate a random Poisson(y) variable
*/
    ignnbn = poidev(y);
    return ignnbn;
}



// digamma function trouvee sur le site: http://www.strw.leidenuniv.nl/~mathar/progs/digamma.c
/*************************************
 An ANSI-C implementation of the digamma-function for real arguments based
 on the Chebyshev expansion proposed in appendix E of
 http://arXiv.org/abs/math.CA/0403344 . For other implementations see
 the GSL implementation for Psi(Digamma) in
 http://www.gnu.org/software/gsl/manual/gsl-ref_toc.html

Richard J. Mathar, 2005-11-24
**************************************/
#include <math.h>

/* constants: pi, Euler's constant, and log(2) */
#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif
#ifndef M_GAMMAl
#define M_GAMMAl 0.5772156649015328606065120900824024L
#endif
#ifndef M_LN2l
#define M_LN2l 0.6931471805599453094172321214581766L
#endif

double digammal(double x)
{
	/* force into the interval 1..3 */
	if( x < 0.0L )
		return digammal(1.0L-x)+M_PIl/tanl(M_PIl*(1.0L-x)) ;	/* reflection formula */
	else if( x < 1.0L )
		return digammal(1.0L+x)-1.0L/x ;
	else if ( x == 1.0L)
		return -M_GAMMAl ;
	else if ( x == 2.0L)
		return 1.0L-M_GAMMAl ;
	else if ( x == 3.0L)
		return 1.5L-M_GAMMAl ;
	else if ( x > 3.0L)
		/* duplication formula */
		return 0.5L*(digammal(x/2.0L)+digammal((x+1.0L)/2.0L))+M_LN2l ;
	else
	{
		/* Just for your information, the following lines contain
		* the Maple source code to re-generate the table that is
		* eventually becoming the Kncoe[] array below
		* interface(prettyprint=0) :
		* Digits := 63 :
		* r := 0 :
		*
		* for l from 1 to 60 do
		* 	d := binomial(-1/2,l) :
		* 	r := r+d*(-1)^l*(Zeta(2*l+1) -1) ;
		* 	evalf(r) ;
		* 	print(%,evalf(1+Psi(1)-r)) ;
		*o d :
		*
		* for N from 1 to 28 do
		* 	r := 0 :
		* 	n := N-1 :
		*
 		*	for l from iquo(n+3,2) to 70 do
		*		d := 0 :
 		*		for s from 0 to n+1 do
 		*		 d := d+(-1)^s*binomial(n+1,s)*binomial((s-1)/2,l) :
 		*		od :
 		*		if 2*l-n > 1 then
 		*		r := r+d*(-1)^l*(Zeta(2*l-n) -1) :
 		*		fi :
 		*	od :
 		*	print(evalf((-1)^n*2*r)) ;
 		*od :
 		*quit :
		*/
		static double Kncoe[] = { .30459198558715155634315638246624251L,
		.72037977439182833573548891941219706L, -.12454959243861367729528855995001087L,
		.27769457331927827002810119567456810e-1L, -.67762371439822456447373550186163070e-2L,
		.17238755142247705209823876688592170e-2L, -.44817699064252933515310345718960928e-3L,
		.11793660000155572716272710617753373e-3L, -.31253894280980134452125172274246963e-4L,
		.83173997012173283398932708991137488e-5L, -.22191427643780045431149221890172210e-5L,
		.59302266729329346291029599913617915e-6L, -.15863051191470655433559920279603632e-6L,
		.42459203983193603241777510648681429e-7L, -.11369129616951114238848106591780146e-7L,
		.304502217295931698401459168423403510e-8L, -.81568455080753152802915013641723686e-9L,
		.21852324749975455125936715817306383e-9L, -.58546491441689515680751900276454407e-10L,
		.15686348450871204869813586459513648e-10L, -.42029496273143231373796179302482033e-11L,
		.11261435719264907097227520956710754e-11L, -.30174353636860279765375177200637590e-12L,
		.80850955256389526647406571868193768e-13L, -.21663779809421233144009565199997351e-13L,
		.58047634271339391495076374966835526e-14L, -.15553767189204733561108869588173845e-14L,
		.41676108598040807753707828039353330e-15L, -.11167065064221317094734023242188463e-15L } ;

		double Tn_1 = 1.0L ;	/* T_{n-1}(x), started at n=1 */
		double Tn = x-2.0L ;	/* T_{n}(x) , started at n=1 */
		double resul = Kncoe[0] + Kncoe[1]*Tn ;

		x -= 2.0L ;

		for(unsigned int n = 2 ; n < sizeof(Kncoe)/sizeof(double) ;n++)
		{
			const double Tn1 = 2.0L * x * Tn - Tn_1 ;	/* Chebyshev recursion, Eq. 22.7.4 Abramowitz-Stegun */
			resul += Kncoe[n]*Tn1 ;
			Tn_1 = Tn ;
			Tn = Tn1 ;
		}
		return resul ;
	}
}

#ifdef TEST

/* an alternate implementation for test purposes, using formula 6.3.16 of Abramowitz/Stegun with the
   first n terms */
#include <stdio.h>

double digammalAlt(double x, int n)
{
	/* force into the interval 1..3 */
	if( x < 0.0L )
		return digammalAlt(1.0L-x,n)+M_PIl/tanl(M_PIl*(1.0L-x)) ;	/* reflection formula */
	else if( x < 1.0L )
		return digammalAlt(1.0L+x,n)-1.0L/x ;
	else if ( x == 1.0L)
		return -M_GAMMAl ;
	else if ( x == 2.0L)
		return 1.0L-M_GAMMAl ;
	else if ( x == 3.0L)
		return 1.5L-M_GAMMAl ;
	else if ( x > 3.0L)
		return digammalAlt(x-1.0L,n)+1.0L/(x-1.0L) ;
	else
	{
		x -= 1.0L ;
		register double resul = -M_GAMMAl ;

		for( ; n >= 1 ;n--)
			resul += x/(n*(n+x)) ;
		return resul ;
	}
}


#endif /* TEST */

// trigamma trouvee sur: http://www.codecogs.com/index.php
// CodeCogs GNU General Public License Agreement
// Copyright (C) 2004-2005 CodeCogs, Zyba Ltd, Broadwood, Holford, TA5 1DU, England.
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by CodeCogs.
// You must retain a copy of this licence in all copies.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the Adapted GNU General Public License for more details.
//
// *** THIS SOFTWARE CAN NOT BE USED FOR COMMERCIAL GAIN. ***
// ---------------------------------------------------------------------------------

double trigamma(double x)
{
  if (x <= 1e-4)
    return 1.0 / x / x;

  double value = 0.0;
  while (x < 5.0)
    value += 1.0 / x / x, x += 1.0;
  double y = 1.0 / x / x;
  return value + 0.5 * y + (1.0 + y * (1/6.0 + y * (y * (1/42.0 - y /30.0) -
1/30.0))) / x;
}

//-------distribution non-central chi2
// algorithme de Ding et al., ne marche pas pour non-central par grand; il faut faire approx normale!
/*
#define ERRMAX 1.0e-6
#define ITRMAX 1000
#define ZERO 0.0
#define ONE 1.0
#define TWO 2.0

double pchisqNC(double X, double THETA,double F)
{bool FLAG;
double BOUND;
double CHI2NC=X;
double LAM=THETA/2;
double N=1;
double U=exp(-LAM);
double V=U;
double X2=X/2;
double F2=F/2;
double T=pow(X2,F2)*exp(-X2-gammln(F2+1));
double TERM=V*T;
cout<<pow(X2,F2)*exp(-X2-gammln(F2+1)-LAM)<<"/";
CHI2NC=TERM;
FLAG=0;
label10:
         if(F+2*N-X<=0) goto label30;
FLAG=1;
label20: BOUND=T*X/(F+2*N-X);
    if((BOUND>ERRMAX)&((int)N<=ITRMAX)) goto label30;
    if(BOUND>ERRMAX) {char pause;
                     cout<<"La fonction pchisqNC n'a pas converge \n Appuyer sur n'importe quelle touche pour continuer...";
                     cin>>pause;}
return CHI2NC;
label30: U=U*LAM/N;
V=V+U;
T=T*X/(F+2*N);
TERM=V*T;
CHI2NC=CHI2NC+TERM;
N=N+1;
if(FLAG) goto label20;
goto label10;
}
*/

//-------distribution non-central chi2
// algorithme de Ding et al., avec approximation normale; il faut faire approx normale!

double erfcc(double x)
{
	double t,z,ans;

	z=fabs(x);
	t=1.0/(1.0+0.5*z);
	ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
	return x >= 0.0 ? ans : 2.0-ans;
}

#define max(a,b)  ((a)<(b)?(b):(a))
#define BIG 300.

#define ERRMAX 1.0e-6
#define ITRMAX 1000

double pnorm(double x){return 0.5*erfcc(-x/1.41421356237309504880168872421);}

double logPLnormDiscrete(double x,double mean,double sd)
{return log(pnorm((log(x+1)-log(mean))/sd)-pnorm((log(x)-log(mean))/sd));}

double pchisqNCnull(double X, double THETA)
{if(X<0) return 0;
if ( THETA > BIG ) {return 0.5*erfcc(-(X-THETA)/(2.*sqrt(2.*THETA)));}
bool FLAG;
double BOUND;
double CHI2NC=X;
double LAM=THETA/2;
double N=1;
double U=exp(-LAM);
double V=U;
double X2=X/2;
double T=exp(-X2);
double TERM=V*T;
CHI2NC=TERM;
FLAG=0;
label10:
         if(2*N-X<=0) goto label30;
FLAG=1;
label20: BOUND=T*X/(2*N-X);
    if((BOUND>ERRMAX)&((int)N<=ITRMAX)) goto label30;
    if(BOUND>ERRMAX) {char pause;
                     cout<<"La fonction pchisqNC n'a pas converge \n Appuyer sur n'importe quelle touche pour continuer...";
                     cin>>pause;}
return CHI2NC;
label30: U=U*LAM/N;
V=V+U;
T=T*X/(2*N);
TERM=V*T;
CHI2NC=CHI2NC+TERM;
N=N+1;
if(FLAG) goto label20;
goto label10;
}

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Define function for computing non-central       */
/* chi-squared distribution function (NCCSDF).     */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*  Non-central chi-square tail probabilities are computed
    using a recursive algorithm proposed by Posten (1989),
    exploiting the fact that such probabilities may be
    written as an infinite weighted sum of central chi-square
    tail probabilities, where the weights are Poisson
    probability masses.  */
/* Note: Need 'gammp.c' and 'gammln.c'.  */
double pchisqNCPosten(double x,double d,double ff)
/* x := critical value, d := non-centrality par, ff := degrees of freedom.                    */
{
/* Initialize and calculate initial values.  */
  double lam, p, w, lp, lg, ls, s, psum;

/* if d large, then use the Gaussian approximation; otherwise do it	*/
/* the hard way.
							*/
  if ( d > BIG ) {
    psum = 0.5*erfcc(-(x-ff-d)/(2.*sqrt(ff+2.*d)));
  } else {
    int i, n;
    lam = d/2.;
    lp = -lam;
    p = exp(lp);
    w = gammp(ff/2., x/2.);
    lg = gammln(ff/2.);
    ls = ((ff-2)/2.);
    ls = ls*log(x) - (x/2.) - ((ff-2.)/2.)*log(2.) - lg;
    s = exp(ls);
    psum = p*w;
    /*psum = exp(lp+log(w)); DDD */

    /* Calculate Poisson weighted sum in a recursive fashion, */
    /* up to n terms, where n is calculated to be 10 std      */
    /* above the mean corresponding to the Poisson weights    */
    /* (i.e. lam + 10*sqrt(lam)).  This seems to yield pvals  */
    /* which can be expected to be accurate up to 10^(-20).   */

//MODIF FAITE PAR SC: AJOUT DE (int)
    n =(int)floor(lam + 10.*pow(lam,0.5));
    for(i=1;i<n+1;++i) {
      lp = log(lam) - log((double) i) + log(p);
      p = exp(lp);
      s = (x/(ff + 2.*i - 2.))*s;
      w = -s + w;
      psum = psum + p*w;
      /*psum = psum + exp(lp+log(w)); DDD*/
    }
  }
  return(psum);
}

double bessi1(double x)
{
	double ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
			+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
	} else {
		y=3.75/ax;
		ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
			-y*0.420059e-2));
		ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
			+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
		ans *= (exp(ax)/sqrt(ax));
	}
	return x < 0.0 ? -ans : ans;
}

//----------LES DENSITES

double pnorm(double x,double mean,double var)    //fonction cumulee de la loi normale
{if (x<mean) return 1-0.5*(1+gammp(0.5,0.5*pow((x-mean),2)/var));
else return 0.5*(1+gammp(0.5,0.5*pow((x-mean),2)/var));}

double plnorm(double x,double mean,double var)
{return pnorm(log(x),mean,var);}

double dlnorm(double x,double mean,double var)
{return exp(-pow(log(x)-mean,2)/2./var) /x/pow(var,0.5)/2.50662827;
//return -0.918938533204672741780329736406-0.5*log(var)-pow(x-mean,2)/2/var;
}

double logdlnorm(double x,double mean,double var)
{return -pow(log(x)-mean,2)/2./var -log(x)-0.5*log(var)-0.918938533;
//return -0.918938533204672741780329736406-0.5*log(var)-pow(x-mean,2)/2/var;
}

double ldnorm(double x,double mean,double var)
{return -pow(x-mean,2)/2./var -0.5*log(var)-0.918938533;
}



double pweibull(double x,double lambda,double k)
{return 1-exp(-pow(x/lambda,k));}

double dweibull(double x,double lambda,double k)
{return (k/lambda)*pow(x/lambda,k-1)*exp( - pow(x/lambda,k) );}



double dlchisqNCnull(double x,double theta)
{if(theta>200) {return dlnorm(x,theta,4*theta);}
return -0.6931471805599453+0.5*(log(theta)-log(x))-0.5*(theta+x)+log(bessi1(sqrt(theta*x)));//on supprime le conditionnement le 06.05.11:-log(1-exp(-0.5*theta));
}
double dlnegbin(int x,double r,double p) //Modified in July 2015 p proba success, r number of failures, based on wikipedia definition
{return gammln(r+x)-gammln(r)-gammln(x+1)+x*log(p)+r*log(1-p);}

double dlpois(int x,double r)
{return -r+x*log(r)-gammln(x+1);}

double ldbin(int n, int k, double p)
{if(n==0) return 0;
//if((n*p>5)&(n*(1-p)>5)) return dlnorm(k,n*p,n*p*(1-p));
//double X=gammln(n+1)-gammln(k+1)-gammln(n-k+1)+k*log(p)+(n-k)*log(1-p);
/*if(((X<100000)&(X>-100000))==0)
    {cout<<"Z="<<X<<","<<n<<","<<k<<","<<p<<"\n";}*/
return gammln(n+1)-gammln(k+1)-gammln(n-k+1)+k*log(p)+(n-k)*log(1-p);}

double ldmultinomial(int n, double *proba,int size,int* deviate)
{int i;
double ldensity=gammln(n+1);
for(i=0;i<size;i++) ldensity+=-gammln(deviate[i]+1)+deviate[i]*log(proba[i]);
return ldensity;
}

//-------------incomplete beta function from numerical recipes
/*
double betai(double a, double b, double x)
{
	double betacf(double a, double b, double x);
	double gammln(double xx);
	void nrerror(char error_text[]);
	double bt;

	if (x < 0.0 || x > 1.0) cout<<"Bad x in routine betai\n";
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*betacf(a,b,x)/a;
	else
		return 1.0-bt*betacf(b,a,1.0-x)/b;
}*/

double betai(double a, double b, double x) //modification de numerical recipes
{
	double betacf(double a, double b, double x);
	double gammln(double xx);
	void nrerror(char error_text[]);
	double lbt;

	if (x < 0.0 || x > 1.0)
		{//cout<<"Bad x in routine betai\n";
		return -1;}//-1 pour erreur
	else
		lbt=gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x);
	if (x < (a+1.0)/(a+b+2.0))
		{double K=betacf(a,b,x);
		if(K==-1) return -1;
		else return exp(lbt+log(K)-log(a));
		}
	else{double KK=betacf(b,a,1.0-x);
		if(KK==-1) return -1;
		else return 1.0-exp(lbt+log(KK)-log(b));
		}
}


double pbeta(double x,double a,double b)
{return betai(a,b,x);}

double dbeta(double x,double a, double b)
{return exp( gammln(a+b)-gammln(a)-gammln(b) + (a-1)*log(x) + (b-1)*log(1-x) );}

#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

double betacf(double a, double b, double x)
{
	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;

	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (m > MAXIT) {//cout<<"a or b too big, or MAXIT too small in betacf \n";
					return -1;}
	return h;
}
#undef MAXIT
#undef EPS
#undef FPMIN

double pnegbin(double x, double n, double p)
{double K=betai(n,x+1,p);
if(K==-1) return -1;
else return K;}
