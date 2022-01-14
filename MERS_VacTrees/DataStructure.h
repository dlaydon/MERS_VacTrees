#ifndef DATA_STRUCTURE_H
#define DATA_STRUCTURE_H

#include <vector>
#include <set>

//RANDOM NUMBERS
extern double ran1(long *);
extern double runif();
extern double rnorm();
extern int rbin(int n, double p);
extern int poidev(double x);
extern double ldbin(int n, int k, double p);
extern double ldmultinomial(int n, double *proba,int size,int* deviate);
extern double gammln(double x);
extern double gammq(double a, double x);
extern void rmultinomial(int n, double *proba,int size,int* deviate);
extern int rmultinomialOneValue(std::vector<double> &proba, int size);
extern int rmultinomialOneValue(std::vector<double> &proba, std::set<int> &possibleInfectors);
extern int rnegbin(double n,double p);
extern double digammal(double x);
extern double trigamma(double x);
extern double pgamma(double x,double a, double b);
extern double dgamma(double x,double a,double b);
extern double sgamma(double a); 

extern double ldgamma(double x,double a,double b);
extern double rgamma(double a,double b); // careful - can general 0 for small a ->rgammaBigger0
extern double pnorm(double x,double mean,double var);
extern double plnorm(double x,double mean,double var);
extern double dlnorm(double x,double mean,double var);
extern double ldnorm(double x,double mean,double var);
extern double logdlnorm(double x,double mean,double var);
extern double pweibull(double x,double lambda,double k);
extern double dweibull(double x,double lambda,double k);
extern double pchisqNCnull(double X, double nonCentralPar);
extern double pnorm(double x);
extern double logPLnormDiscrete(double x,double mean,double sd);
extern double bessi1(double x);
extern double dlchisqNCnull(double x,double theta);
extern double dlnegbin(int x,double n,double p);
extern double dlpois(int x,double r);
extern double pnegbin(double k, double a, double b);
extern double erff(double x);
extern double rTruncNormalRight(double mu,double var,double muMax);
extern double choose(int n,int k);
extern double pbeta(double x,double a,double b);
extern double dbeta(double x,double a, double b);

//============
extern long SEED; // root random number

#endif
