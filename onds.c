#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <complex.h>
#include "onds.h"
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
double gammafunc(double x, double a, double b, double U, double v);
double ReMBSigma(struct params *pms, char c, double x);
double ImMBSigma(struct params *pms, char c, double x);
double SigmaMBLesser(struct params *pms, char c, double x);
double ReMBSigmaDeriv(struct params *pms, char c, double x);
double ImMBSigmaDeriv(struct params *pms, char c, double x);
double SigmaMBLesserDeriv(struct params *pms, char c, double x);
double SolveFunc(double n, void* p);
double SolveFuncDeriv(double n, void* p);
double IntFunc(double x, void* p);
double SpectrumFunc(double x, void* p);
double CurrentFunc(double x, void* p);
double IntFuncDeriv(double x, void* p);
double IntFuncDerivVg(double x, void* p);
double fd(double en, double beta, double mu);
double func2L(double x, void* p);
double func2R(double x, void* p);
double gammafunc(double x, double a, double b, double U, double v);
void SolveFuncFd(double n, void* p, double *y, double *dy);
void print_approx(char c);
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
gsl_integration_workspace *ws1;
gsl_integration_workspace *ws;
gsl_integration_workspace *ws2;
gsl_integration_workspace *ws3;
gsl_integration_workspace *ws4;
gsl_integration_workspace *ws5;
gsl_integration_workspace *ws6;
gsl_function F;
gsl_function F1;
gsl_function F2;
gsl_function F3;
gsl_function F4;
gsl_function F5;
gsl_function F6;
gsl_function FD;
struct params pm1;
struct params pm2;
struct params pm3;
struct params pm3;
struct params pm4;
struct params pm5;
struct params pm6;
const gsl_root_fdfsolver_type *T;
gsl_root_fdfsolver *s;
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
// Fermi function
//-------------------------------------------------------------------------------
double fd(double en, double beta, double mu)
{
    return ( 1.0/( exp(beta*(en-mu)) + 1.0) );
}
//===============================================================================
// lammbda2 and gamma 2 determine the embedding self-energy according to the
// equation (16) on J. Phys.: Conf. Ser. 220 012018 (2010).
//===============================================================================
//-------------------------------------------------------------------------------
//  Real part of the embedding self-energy
//  x = energy
//  a = on-site energy of the lead-sites
//  b = hopping between the lead sites
//  U = Bias voltage
//  v = coupling coeffient
//-------------------------------------------------------------------------------
double lambda2(double x, double a, double b, double U, double v){
    double em;
    if( x-a-U > 2.0*fabs(b) ){
        //em = (2.0*v*v*sqrt( 1 - ( (x-a-U)*(x-a-U) ) / ( 4.0*b*b ) ) / b);
        em = v*v/(2.0*b*b)*((x-a-U)-sqrt((x-a-U)*(x-a-U)-4.0*b*b));
    }
    else if( x-a-U < -2.0*fabs(b) ){
        //em = (2.0*v*v*sqrt( 1 + ( (x-a-U)*(x-a-U) ) / ( 4.0*b*b ) ) / b);
        em = v*v/(2.0*b*b)*((x-a-U)+sqrt((x-a-U)*(x-a-U)-4.0*b*b));
    }
    else if ( fabs(x-a-U) < 2*fabs(b)){
        em = v*v*(x-a-U)/(2.0*b*b);
    }
    else em = 0.0;

    return em;

}
//-------------------------------------------------------------------------------
//  Imaginary part of the embedding self-energy
//  x = energy
//  a = on-site energy of the lead-sites
//  b = hopping between the lead sites
//  U = Bias voltage
//  v = coupling coeffient
//-------------------------------------------------------------------------------
double gamma2(double x, double a, double b, double U, double v){
    double em
    if( fabs(x-a-U)<=2.0*b ){
        em = v*v/(4.0*b*b)*(sqrt((x-a-U)*(x-a-U)-4.0*b*b));
    }
    else em = 0.0;

    return em;
}
//===============================================================================
// Otherway of calcualting the embedding self-energy with gammafunc and lambdaEM
// routines.
//===============================================================================
//-------------------------------------------------------------------------------
// Imaginary part of the embedding self-energy
//-------------------------------------------------------------------------------
double gammafunc(double x, double a, double b, double U, double v)
{
    if( fabs(x-a-U)<=2.0*b ){
        return(2.0*v*v*sqrt( 1 - ( (x-a-U)*(x-a-U) ) / ( 4.0*b*b ) ) / b);
    }
    else{
        return 0.0;
    }
}
//-------------------------------------------------------------------------------
//  Calculates the integrand in the principal value integral
//  evalated to obtain the real part of the embedding self-energy
//-------------------------------------------------------------------------------
double lambdaEM(struct params *pms, double x){

    double lambdaL = 0.0;
    double lambdaR = 0.0;
    double llimit, ulimit, result;
    double epsabs, epsrel, abserr;
    int limit;
    double ul = pms->ul;
    double ur = pms->ur;
    double a = pms->a;
    double b = pms->b;

    struct params pm1 = *pms;

    epsabs = 0.00000001;       // absolute error limit of the result
    epsrel = 0.00000001;       // relative error limit of the result
    abserr = 0.0;              // absolute error of the result
    limit  = 2000;             // maximum number of subintervals in integration routine (must be < size of workspace)

    // calculate LambdaL
    llimit = a - 2.0*b + ul ;
    ulimit = a + 2.0*b + ul ;

    F1.params = &pm1;
    F1.function = &func2L;

    // cauchy principal value integration I = \int_a^b dx' f(x') / (x' - x)
    gsl_integration_qawc(&F1,llimit,ulimit,x,epsabs,epsrel,limit,ws1,&result,&abserr);
    result*=( -1.0/(2.0*M_PI) );
    lambdaL = result;

    // calculate LambdaR
    llimit = a - 2.0*b + ur;
    ulimit = a + 2.0*b + ur;

    F1.params = &pm1;
    F1.function = &func2R;

    // cauchy principal value integration I = \int_a^b dx' f(x') / (x' - x)
    gsl_integration_qawc(&F1,llimit,ulimit,x,epsabs,epsrel,limit,ws1,&result,&abserr);
    result*=(-1.0/(2.0*M_PI));
    lambdaR = result;

    return lambdaL + lambdaR;
}
//-------------------------------------------------------------------------------
// Left component of the embedding self-energy
//-------------------------------------------------------------------------------
double func2L(double x, void* p){
    struct params* pms = (struct params*)p;
    double a  = pms->a;
    double b  = pms->b;
    double ul = pms->ul;
    double vl = pms->V_l;
    double f  = ( gammafunc(x,a,b,ul,vl) );
    return f;
}

//-------------------------------------------------------------------------------
//  Right component of the embedding self-energy
//-------------------------------------------------------------------------------
double func2R(double x, void* p){
    struct params* pms = (struct params*)p;
    double a  = pms->a;
    double b  = pms->b;
    double ur = pms->ur;
    double vr = pms->V_r;
    double f  = ( gammafunc(x,a,b,ur,vr) );
    return f;
}
//-------------------------------------------------------------------------------
// Returns the real part of the retarded component of the many-body
// self-energy
//-------------------------------------------------------------------------------
double ReMBSigma(struct params *pms, char c, double x){

    double U_h  = pms->U_h;
    double ds   = pms->ds;
    double damp = pms->damp;
    double e0   = pms->e0;

    if(c == 'B'){  	// 2B
    	return creal((U_h*U_h*ds*(1.0-ds)) / (x-e0-U_h*ds+I*3.0*damp) );
    }
    else if (c == 'G'){	// GW
    	return creal((2.0*U_h*U_h*ds*(1.0-ds)) / (x-e0-U_h*ds+I*3.0*damp));
    }
    else if( c =='T'){	//T-Matrix
    	return 0.0;
    }
    else if (c == 'H'){ //Hartree-Fock
    	return 0.0;
    }
}
//-------------------------------------------------------------------------------
// Returns the imaginary part of the retarded component of the many-body
// self-energy
//-------------------------------------------------------------------------------
double ImMBSigma(struct params *pms, char c, double x){

    double U_h  = pms->U_h;
    double ds   = pms->ds;
    double damp = pms->damp;
    double e0   = pms->e0;

    if (c == 'B'){ 	//2B
        return cimag((U_h*U_h*ds*(1.0-ds)) / (x-e0-U_h*ds+I*3.0*damp) );
    }
    else if (c == 'G'){ //GW
        return cimag((2.0*U_h*U_h*ds*(1.0-ds)) / (x-e0-U_h*ds+I*3.0*damp) );
    }
    else if (c == 'T'){ // T-matrix
    	return 0.0;
    }
    else if( c == 'H'){ //Hartree-Fock
        return  0.0;
    }
}
//-------------------------------------------------------------------------------
// Returns the lesser component of the many-body
// self-energy
//-------------------------------------------------------------------------------
double SigmaMBLesser(struct params *pms, char c, double x){

    double U_h  = pms->U_h;
    double ds   = pms->ds;
    double damp = pms->damp;
    double e0   = pms->e0;

    if(c == 'B'){ 	//2B
        return (-I*6.0*U_h*U_h*ds*ds*(ds-1.0)*damp) / ( (x-e0-U_h*ds)*(x-e0-U_h*ds)+9.0*damp*damp);
    }
    else if (c == 'G') { //GW
        return (-I*12.0*U_h*U_h*ds*ds*(ds-1.0)*damp) / ( (x-e0-U_h*ds)*(x-e0-U_h*ds)+9.0*damp*damp);
    }
    else if (c == 'T'){ // T-marix
    	return 0.0;
    }
    else{		//Hartree-Fock
        return 0.0;
    }
}
//-------------------------------------------------------------------------------
// Returns the real part of the retarded component of the many-body
// self-energy
//-------------------------------------------------------------------------------
double ReMBSigmaDeriv(struct params *pms, char c, double x){

    double U_h  = pms->U_h;
    double ds   = pms->ds;
    double damp = pms->damp;
    double e0   = pms->e0;
    double numer, denum;

    if(c == 'B'){  	// 2B
    	numer = ((U_h*U_h*(1-ds)-U_h*U_h*ds)*(x-e0-U_h*ds+I*3.0*damp) + U_h*U_h*U_h*ds*(1-ds));
    	denum = (x-e0-U_h*ds+I*3.0*damp)*(x-e0-U_h*ds+I*3.0*damp);
    	return creal(numer/denum);
    }
    else if (c == 'G'){	// GW
       	numer = ((U_h*U_h*(1-ds)-U_h*U_h*ds)*(x-e0-U_h*ds+I*3.0*damp) + U_h*U_h*U_h*ds*(1-ds));
    	denum = (x-e0-U_h*ds+I*3.0*damp)*(x-e0-U_h*ds+I*3.0*damp);
    	return 2.0*creal(numer/denum);
    }
    else if( c =='T'){	//T-Matrix
    	return 0.0;
    }
    else if (c == 'H'){ //Hartree-Fock
    	return 0.0;
    }
}
//-------------------------------------------------------------------------------
// Returns the imaginary part of the retarded component of the many-body
// self-energy
//-------------------------------------------------------------------------------
double ImMBSigmaDeriv(struct params *pms, char c, double x){

    double U_h  = pms->U_h;
    double ds   = pms->ds;
    double damp = pms->damp;
    double e0   = pms->e0;
    double numer, denum;

    if(c == 'B'){  // 2B
       	numer = ((U_h*U_h*(1-ds)-U_h*U_h*ds)*(x-e0-U_h*ds+I*3.0*damp) + U_h*U_h*U_h*ds*(1-ds));
    	denum = (x-e0-U_h*ds+I*3.0*damp)*(x-e0-U_h*ds+I*3.0*damp);
    	return cimag(numer/denum);
    }
    else if (c == 'G'){	// GW
       	numer = ((U_h*U_h*(1-ds)-U_h*U_h*ds)*(x-e0-U_h*ds+I*3.0*damp) + U_h*U_h*U_h*ds*(1-ds));
    	denum = (x-e0-U_h*ds+I*3.0*damp)*(x-e0-U_h*ds+I*3.0*damp);
    	return 2.0*cimag(numer/denum);
    }
    else if( c =='T'){	//T-Matrix
    	return 0.0;
    }
    else if (c == 'H'){ //Hartree-Fock
    	return 0.0;
    }
}
//-------------------------------------------------------------------------------
//  Returns the lesser component of the many-body
//  self-energy
//-------------------------------------------------------------------------------
double SigmaMBLesserDeriv(struct params *pms, char c, double x){

    double U_h  = pms->U_h;
    double ds   = pms->ds;
    double damp = pms->damp;
    double e0   = pms->e0;
    double denum, numer1, numer2;

    if(c == 'B'){ 	//2B
        denum = ((x-e0-U_h*ds)*(x-e0-U_h*ds)+9.0*damp*damp)*((x-e0-U_h*ds)*(x-e0-U_h*ds)+9.0*damp*damp);
        numer1 = (12*U_h*U_h*ds*(ds-1)*damp+6.0*U_h*U_h*ds*ds*damp)*((x-e0-U_h*ds)*(x-e0-U_h*ds)+9.0*damp*damp);
        numer2 = 12.0*U_h*U_h*U_h*ds*ds*(ds-1)*damp*(x-e0-U_h*ds);
        return ( -I*(numer1+numer2) / denum);
    }
    else if (c == 'G') { //GW
        denum = ((x-e0-U_h*ds)*(x-e0-U_h*ds)+9.0*damp*damp)*((x-e0-U_h*ds)*(x-e0-U_h*ds)+9.0*damp*damp);
        numer1 = (12*U_h*U_h*ds*(ds-1)*damp+6.0*U_h*U_h*ds*ds*damp)*((x-e0-U_h*ds)*(x-e0-U_h*ds)+9.0*damp*damp);
        numer2 = 12.0*U_h*U_h*U_h*ds*ds*(ds-1)*damp*(x-e0-U_h*ds);
        return ( -I*2.0*(numer1+numer2) / denum);
    }
    else if (c == 'T'){ // T-marix
    	return 0.0;
    }
    else{		//Hartree-Fock
        return 0.0;
    }
}
//-------------------------------------------------------------------------------
//  Calculates the integrand in the density integral
//-------------------------------------------------------------------------------
double IntFunc(double x, void* p){
    struct params* pms = (struct params*)p;

    double a  = pms->a;
    double b  = pms->b;
    double vl = pms->V_l;
    double vr = pms->V_r;
    double mu = pms->mu;
    double e0 = pms->e0;
    double ul = pms->ul;
    double ur = pms->ur;
    double Vg = pms->Vg;
    double beta = pms->beta;
    double ds = pms->ds;
    double U_h = pms->U_h;
    char c = pms->c;
    double v_xc = 0.0;

    double denum = 0.0;
    double numer;

    double gamL = 0.0;
    double gamR = 0.0;

    double lambdaL = 0.0;
    double lambdaR = 0.0;

    double ReMB, ImMB, SigmaLS;
    double lambda, gamma;
    double f_L,f_R;


    double hub = U_h;
    double tih = ds;

    struct params params = *pms;


    // calculate GammaL and GammaR
    gamL = gammafunc(x,a,b,ul,vl);
    gamR = gammafunc(x,a,b,ur,vr);
    //gamL = gamma2(x,a,b,ul,vl);
    //gamR = gamma2(x,a,b,ur,vr);
    gamma = gamL + gamR;

    // calculate fermi functions.
    f_L = fd(x-ul,beta,mu);
    f_R = fd(x-ur,beta,mu);

    // calculate LambdaL
    //lambda = lambdaEM(&params,x);
    //lambda = lambda2(&params,x);
    lambdaL = lambda2(x,a,b,ul,vl);
    lambdaR = lambda2(x,a,b,ur,vr);
    lambda = lambdaL + lambdaR;


    ReMB = ReMBSigma(&params, c, x);
    ImMB = ImMBSigma(&params, c, x);
    SigmaLS = SigmaMBLesser(&params,c,x);

    numer = f_L*gamL + f_R*gamR - I*SigmaLS;
    denum = (x-e0-U_h*ds-v_xc-ReMB-Vg-lambda)*(x-e0-U_h*ds-v_xc-ReMB-Vg-lambda) + (gamma/2.0-ImMB)*(gamma/2.0-ImMB);

    return (numer/denum);
}

//-------------------------------------------------------------------------------
//   Calculates the integrand in the current integral
//-------------------------------------------------------------------------------
double CurrentFunc(double x, void* p){
    struct params* pms = (struct params*)p;

    double a  = pms->a;
    double b  = pms->b;
    double vl = pms->V_l;
    double vr = pms->V_r;
    double mu = pms->mu;
    double e0 = pms->e0;
    double ul = pms->ul;
    double ur = pms->ur;
    double Vg = pms->Vg;
    double beta = pms->beta;
    double ds = pms->ds;
    double U_h = pms->U_h;
    char c = pms->c;

    double hub = U_h;
    double tih = ds;
    double denum = 0.0;
    double numer1, numer2, numer3;

    double gamL = 0.0;
    double gamR = 0.0;

    double ReMB, ImMB, SigmaLS;
    double lambda, gamma;
    double f_L,f_R;

    double v_xc = 0.0;
    struct params params = *pms;

    // calculate GammaL and GammaR
    gamL = gammafunc(x,a,b,ul,vl);
    gamR = gammafunc(x,a,b,ur,vr);
    gamma = gamL + gamR;

    // calculate fermi functions.
    f_L = fd(x-ul,beta,mu);
    f_R = fd(x-ur,beta,mu);

    // calculate LambdaL
    lambda = lambdaEM(&params,x);

    //MB-sigma
    ReMB = ReMBSigma(&params, c, x);
    ImMB = ImMBSigma(&params, c, x);
    SigmaLS = SigmaMBLesser(&params,c,x);

    numer1 = (f_L - f_R)*gamL*gamR;
    numer2 = I*0.5*(gamL-gamR)*SigmaLS;
    numer3 = 0.5*(f_L*gamL-f_R*gamR)*ImMB;
    denum = (x-e0-U_h*ds-v_xc-Vg-lambda-ReMB)*(x-e0-U_h*ds-v_xc-Vg-lambda-ReMB) + (ImMB-gamma/2.0)*(ImMB-gamma/2.0);

    return ((numer1+numer2-numer3)/denum);
}

//-------------------------------------------------------------------------------
//   for the spectral function
//-------------------------------------------------------------------------------
double SpectrumFunc(double x, void* p){
    struct params* pms = (struct params*)p;

    double a  = pms->a;
    double b  = pms->b;
    double vl = pms->V_l;
    double vr = pms->V_r;
    double mu = pms->mu;
    double e0 = pms->e0;
    double ul = pms->ul;
    double ur = pms->ur;
    double Vg = pms->Vg;
    double beta = pms->beta;
    double ds = pms->ds;
    double U_h = pms->U_h;
    char c = pms->c;

    double v_xc = 0.0;

    double denum = 0.0;
    double numer;
    double res;

    double gamL = 0.0;
    double gamR = 0.0;

    double ReMB, ImMB, SigmaLS;
    double lambda, gamma;
    double f_L,f_R;

    struct params params = *pms;

    FILE* fd2;

    // calculate GammaL and GammaR
    gamL = gammafunc(x,a,b,ul,vl);
    gamR = gammafunc(x,a,b,ur,vr);
    gamma = gamL + gamR;

    // calculate fermi functions.
    f_L = fd(x-ul,beta,mu);
    f_R = fd(x-ur,beta,mu);

    // calculate LambdaL
    lambda = lambdaEM(&params,x);

    ReMB = ReMBSigma(&params, c, x);
    ImMB = ImMBSigma(&params, c, x);

    numer = gamma - 2.0*ImMB;
    denum = (x-e0-U_h*ds-v_xc-ReMB-Vg-lambda)*(x-e0-v_xc-U_h*ds-ReMB-Vg-lambda) + (gamma/2.0-ImMB)*(gamma/2.0-ImMB);
    res = (numer/denum)/(2.0*M_PI);

    fd2 = fopen("spectrum","a");
    fprintf(fd2,"%f %f\n",x,res);
    fclose(fd2);

    return (numer/denum);
}


//-------------------------------------------------------------------------------
//   Prints the Green function
//-------------------------------------------------------------------------------
void GreenFunc(double x, void* p){
    struct params* pms = (struct params*)p;

    double a  = pms->a;
    double b  = pms->b;
    double vl = pms->V_l;
    double vr = pms->V_r;
    double mu = pms->mu;
    double e0 = pms->e0;
    double ul = pms->ul;
    double ur = pms->ur;
    double Vg = pms->Vg;
    double beta = pms->beta;
    double ds = pms->ds;
    double U_h = pms->U_h;
    char c = pms->c;

    double denum = 0.0;
    double numer;
    double res;

    double gamL = 0.0;
    double gamR = 0.0;

    double ReMB, ImMB, SigmaLS;
    double lambda, gamma;
    double f_L,f_R;

    struct params params = *pms;

    FILE* fd2;

    // calculate GammaL and GammaR
    gamL = gammafunc(x,a,b,ul,vl);
    gamR = gammafunc(x,a,b,ur,vr);
    gamma = gamL + gamR;


    // calculate LambdaL
    lambda = lambdaEM(&params,x);

    ReMB = ReMBSigma(&params, c, x);
    ImMB = ImMBSigma(&params, c, x);

    numer = 1.0;
    denum = (x-e0-U_h*ds-ReMB-Vg-lambda)*(x-e0-U_h*ds-ReMB-Vg-lambda) + (gamma/2.0-ImMB)*(gamma/2.0-ImMB);
    res = (1/denum);

    fd2 = fopen("green_function.dat","a");
    fprintf(fd2,"%f %f\n",x,res);
    fclose(fd2);
}
//-------------------------------------------------------------------------------
//  Calculates the derivate of the integrand of the density integral with respect
//  to the density
//-------------------------------------------------------------------------------
double IntFuncDeriv(double x, void* p){
    struct params* pms = (struct params*)p;

    double a = pms->a;
    double b = pms->b;
    double vl = pms->V_l;
    double vr = pms->V_r;
    double mu = pms->mu;
    double e0 = pms->e0;
    double ul = pms->ul;
    double ur = pms->ur;
    double Vg = pms->Vg;
    double beta = pms->beta;
    double ds = pms->ds;
    double U_h = pms->U_h;
    double damp = pms->damp;
    char c = pms->c;

    double denum = 0.0;
    double numer1, numer2, numer;
    double f_L, f_R;
    double gamL = 0.0;
    double gamR = 0.0;
    double ReMB, ImMB, SigmaLS;
    double ReMBDeriv, ImMBDeriv, SigmaLSDeriv;
    double lambda, gamma;

    struct params params =*pms;

    // calculate GammaL and GammaR
    gamL = gammafunc(x,a,b,ul,vl);
    gamR = gammafunc(x,a,b,ur,vr);
    gamma = gamL + gamR;

    // calculate fermi functions.
    f_L = fd(x-ul,beta,mu);
    f_R = fd(x-ur,beta,mu);

    // calculate LambdaL
    lambda = lambdaEM(&params, x);

    //MB-sigma
    ReMB = ReMBSigma(&params, c, x);
    ImMB = ImMBSigma(&params, c, x);
    SigmaLS = SigmaMBLesser(&params,c,x);

    //Derivatives of MB-sigma
    ReMBDeriv = ReMBSigmaDeriv(&params, c, x);
    ImMBDeriv = ImMBSigmaDeriv(&params, c, x);
    SigmaLSDeriv = SigmaMBLesserDeriv(&params,c,x);

    numer = f_L*gamL+f_R*gamR-I*SigmaLS;
    denum = (x-e0-U_h*ds-Vg-ReMB-lambda)*(x-e0-U_h*ds-Vg-ReMB-lambda) + (ImMB-gamma/2.0)*(ImMB-gamma/2.0);

    numer1 = -I*0*SigmaLSDeriv*denum;
    numer2 = numer*(2.0*(x-e0-U_h*ds-Vg-ReMB-lambda)*(U_h+ReMBDeriv)-2.0*(ImMB-gamma/2.0)*ImMBDeriv);
    denum*=denum;

    return ((numer1+numer2)/denum);

    /* numer = (f_L*gamL+f_R*gamR)*2.0*U_h*(x-e0-U_h*ds-Vg-lambda);
     denum = (x-e0-U_h*ds-Vg-lambda)*(x-e0-U_h*ds-Vg-lambda) + (gamma/2.0)*(gamma/2.0);
     denum*=denum;
     return (numer/denum);
     */
}
//-------------------------------------------------------------------------------
//  Calculates the derivate of the integrand of the density integral with respect
//  to the gate voltage
//-------------------------------------------------------------------------------
double IntFuncDerivVg(double x, void* p){
    struct params* pms = (struct params*)p;

    double a = pms->a;
    double b = pms->b;
    double vl = pms->V_l;
    double vr = pms->V_r;
    double mu = pms->mu;
    double e0 = pms->e0;
    double ul = pms->ul;
    double ur = pms->ur;
    double Vg = pms->Vg;
    double beta = pms->beta;
    double ds = pms->ds;
    double U_h = pms->U_h;
    double damp = pms->damp;
    char c = pms->c;

    double denum = 0.0;
    double numer;
    double f_L, f_R;
    double gamL = 0.0;
    double gamR = 0.0;
    double ReMB, ImMB, SigmaLS, SigmaLSDE;
    double lambda, gamma;

    struct params params =*pms;

    // calculate GammaL and GammaR
    gamL = gammafunc(x,a,b,ul,vl);
    gamR = gammafunc(x,a,b,ur,vr);
    gamma = gamL + gamR;

    // calculate fermi functions.
    f_L = fd(x-ul,beta,mu);
    f_R = fd(x-ur,beta,mu);

    // calculate LambdaL
    lambda = lambdaEM(&params, x);

    numer = 2.0*(f_L*gamL + f_R*gamR)*(x-e0-U_h*ds-Vg-lambda);
    denum = (x-e0-U_h*ds-Vg-lambda)*(x-e0-U_h*ds-Vg-lambda) + (gamma/2.0)*(gamma/2.0);
    denum*=denum;

    return (numer/denum);
}
//-------------------------------------------------------------------------------
//  Function given to the numerical root finding routine
//-------------------------------------------------------------------------------
double SolveFunc(double n, void* p){

    struct params* pms = (struct params*)p;

    double llimit, ulimit;
    double result;
    double a = pms->a;
    double b = pms->b;
    double ur = pms->ur;
    double ul = pms->ul;

    struct params pm2 =*pms;
    pm2.ds = n;

    double epsabs = 0.000001; // absolute error limit of the result
    double epsrel = 0.000001; // relative error limit of the result
    double abserr = 0.0;      // absolute error of the result
    int limit = 1500;         // maximum number of subintervals in integration routine (must be < size of workspace)

    llimit = a - 2.0*b + ur;   // lower limit of the integration
    ulimit = a + 2.0*b + ul;   // upper limit if the integration

    F2.function = &IntFunc;
    F2.params = &pm2;
    gsl_integration_qags(&F2,llimit,ulimit,epsabs,epsrel,limit,ws2,&result,&abserr);
    result*=(1.0/(2.0*M_PI));

    return n-result;
}
//-------------------------------------------------------------------------------
//  Derivative of the function whose roots are to be determined
//-------------------------------------------------------------------------------
double SolveFuncDeriv(double n, void* p){

    struct params* pms = (struct params*)p;

    double llimit, ulimit;
    double result;
    double a = pms->a;
    double b = pms->b;
    double ur = pms->ur;
    double ul = pms->ul;

    struct params pm3 =*pms;
    pm3.ds = n;

    double epsabs = 0.000001; // absolute error limit of the result
    double epsrel = 0.000001; // relative error limit of the result
    double abserr = 0.0;      // absolute error of the result
    int limit = 1500;         // maximum number of subintervals in integration routine (must be < size of workspace)

    llimit = a - 2.0*b + ur;   // lower limit of the integration
    ulimit = a + 2.0*b + ul;   // upper limit if the integration

    F3.function = &IntFuncDeriv;
    F3.params = &pm3;
    gsl_integration_qags(&F3,llimit,ulimit,epsabs,epsrel,limit,ws3,&result,&abserr);
    result*=(1.0/(2.0*M_PI));

    return 1.0-result;
}
//-------------------------------------------------------------------------------

void SolveFuncFd(double n, void* p, double *y, double *dy){
    struct params *pms = (struct params*)p;

    *y  = SolveFunc(n,pms);
    *dy = SolveFuncDeriv(n,pms);
}
//-------------------------------------------------------------------------------
// MAIN
//-------------------------------------------------------------------------------
int main(void){

    int limit, status, val, sweep, i = 0;
    double epsabs,epsrel,result,abserr, root, resultdn, resultdvg;
    int iter = 0, max_iter = 100;
    double z, z0, y;
    double a;
    double b;
    double ll,ul, Vl, Vr;
    double ds;
    double U_h;
    double Vg;
    double a0;
    double damp;
    double eta;
    char c;
    double density [5];

    struct params pm;

    FILE *log;
    FILE *fd4;

    log = fopen("output.dat","w");
    fclose(log);

    epsabs = 0.000001;   // absolute error limit of the result
    epsrel = 0.000001;   // relative error limit of the result
    abserr = 0.0;        // absolute error of the result
    limit  = 1500;       // maximum number of subintervals in integration routine (must be < size of workspace)

    ws  = gsl_integration_workspace_alloc(2000);
    ws1 = gsl_integration_workspace_alloc(2000);
    ws2 = gsl_integration_workspace_alloc(2000);
    ws3 = gsl_integration_workspace_alloc(2000);
    ws4 = gsl_integration_workspace_alloc(2000);
    ws5 = gsl_integration_workspace_alloc(2000);
    ws6 = gsl_integration_workspace_alloc(2000);
    printf("\nIntegration workspace initialized..\n");

    T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc(T);
    gsl_function_fdf FDF;

    gsl_set_error_handler_off ();

    read_params(&pm);
    print_approx(pm.c);

    if(pm.sweep == 0) pm.ustep = 100.0;
    if(pm.damptype == 1) pm.damp = 2.0*pm.V_l*pm.V_r/pm.b;
    if(pm.aval == 0) pm.a = pm.mu;

    FDF.f = &SolveFunc;
    FDF.df = &SolveFuncDeriv;
    FDF.fdf = &SolveFuncFd;
    FDF.params = &pm;

    log = fopen("output.dat","a");

    switch(pm.val){

      //--------------- Graphical solution-------------------------
        case 1:
            ds = 0.0;
            printf("\nGraphical solution\n\n");
            fprintf(log,"\nGraphical solution\n\n");

            fd4 = fopen("densities.dat","w");
            fprintf(fd4,"#   n    n _integ.   ul        ur      Vg\n");
            fclose(fd4);

            while( pm.ul <= pm.ulmax){
                while (pm.ur >= pm.urmax){

                    ll = pm.a - 2.0*pm.b + pm.ur;   // lower limit of the integration
                    ul = pm.a + 2.0*pm.b + pm.ul;   // upper limit if the integration

                    while( ds <= 0.99 ){
                        pm.ds = ds;
                        F.params = &pm;
                        F.function = &IntFunc;

                        gsl_integration_qags(&F,ll,ul,epsabs,epsrel,limit,ws,&result,&abserr);
                        result*=(1.0/(2.0*M_PI));

                        fd4 = fopen("densities.dat","a");
                        fprintf(fd4,"%f %f %f %f %f\n",ds ,result, pm.ul, pm.ur, pm.Vg);
                        fclose(fd4);

                        printf("%f %f\n",ds,result);

                        ds+=0.009;
                    }
                    ds = 0.0;
                    pm.ur -= pm.ustep;
                }
                pm.ur = pm.urmin;
                pm.ul += pm.ustep;
            }
            break;

        //----------------Numerical solution-----------------------------
        case 2:
            printf ("\nNumerical solution Using %s method.\n\n", gsl_root_fdfsolver_name (s));
            fprintf (log,"\nNumerical solution Using %s method.\n\n", gsl_root_fdfsolver_name (s));

            fd4 = fopen("roots.dat", "w");
            fprintf(fd4, "#   n        err        ul       ur      Vg\n");
            fclose(fd4);

            root=0.0;
            y=0.1;

            while ( pm.ul <= pm.ulmax){
                while ( pm.ur >= pm.urmax){
                    printf("ul = %f, ur =%f ...\n", pm.ul, pm.ur);

    			   	while(y<0.9){
                        z=y;
                        gsl_root_fdfsolver_set (s, &FDF, z);

                        do{
                            iter++;
                            status = gsl_root_fdfsolver_iterate (s);
                            z0 = z;
                            z = gsl_root_fdfsolver_root (s);
                            status = gsl_root_test_delta (z, z0, 0, 1e-3);
                            if (status == GSL_SUCCESS ){
                                if( z > 0.0 && z < 1.0){
                                    fd4 = fopen("roots.dat", "a");
                                    fprintf(fd4,"%f %f %f %f %f\n",z, z - z0, pm.ul, pm.ur, pm.Vg);
                                    fclose(fd4);
                                    printf("%f %f %f %f %f\n",z, z - z0, pm.ul, pm.ur, pm.Vg);
                                    root=z;
                                }
       			    	          }
                            //printf("%d %f %f\n", iter, z, z - z0);
                        }while (status == GSL_CONTINUE && iter < max_iter);

                        y+=0.1;
                        iter = 0;
                    }
                    y=0.1;
                    pm.ur -=pm.ustep;
                }
                y=0.1;
                pm.ur = pm.urmin;
                pm.ul += pm.ustep;
    	    }
          break;
       //----------------Current calculation-----------------------------------
        case 3:

            printf("\nCalculating the current\n");
            fprintf(log,"\nCalculating the current\n");

            fd4 = fopen("current.dat","w");
            fprintf(fd4, "#  n      current      ul       ur     vg\n");
            fclose(fd4);

            ll = pm.a - 2.0*pm.b + pm.ur;   // lower limit of the integration
            ul = pm.a + 2.0*pm.b + pm.ul;   // upper limit if the integration

            printf("ds     current \n");

            for(i = 0; i < pm.nsol; i++){
                pm.ds = pm.sol[i];
                F.params = &pm;
                F.function = &CurrentFunc;

                gsl_integration_qags(&F,ll,ul,epsabs,epsrel,limit,ws,&result,&abserr);
                result*=(1.0/(M_PI));

                printf("%f %f\n", pm.ds, result);
                fd4 = fopen("current.dat", "a");
                fprintf(fd4,"%f %f %f %f %f\n", pm.ds, result, pm.ul, pm.ur, pm.Vg);
 	            fclose(fd4);

    	    }
          break;
        //----------------Spectral function-----------------------------------
        case 4:

            printf("\nCalculating the spectral funcion\n");
            fprintf(log,"\nCalculating the spectral funcion\n");

            fd4 = fopen("spectrum.dat", "w");
            fprintf(fd4, "#   x    A\n");
            fclose(fd4);

            ll = pm.a - 2.0*pm.b + pm.ur;   // lower limit of the integration
            ul = pm.a + 2.0*pm.b + pm.ul;   // upper limit if the integration

            printf("ds     dst \n");

            for(i = 0; i < pm.nsol; i++){
                pm.ds = pm.sol[i];
                F.params = &pm;
                F.function = &SpectrumFunc;

                gsl_integration_qags(&F,ll,ul,epsabs,epsrel,limit,ws,&result,&abserr);
                result*=(1.0/(2.0*M_PI));

                printf("%f %f\n", pm.ds, result);

    	    }
            break;

        //----------------Checking the derivatives-----------------------------------
        // WORKS AT  HARTREE LEVEL
        case 5:

            fd4 = fopen("derivatives.dat","w");
            fprintf(fd4, "#  n      df/dVg      df/dn      dn/dVg      ul       ur     vg\n");
            fclose(fd4);

    	      ll = pm.a - 2.0*pm.b + pm.ur;   // lower limit of the integration
            ul = pm.a + 2.0*pm.b + pm.ul;   // upper limit if the integration

            printf("  ds      df/dVg     df/dn   (df/dVg)/(1-df/dn) \n");

            for(i = 0; i < pm.nsol; i++){
                pm.ds = pm.sol[i];

                F.params = &pm;
                F.function = &IntFuncDeriv;
                gsl_integration_qags(&F,ll,ul,epsabs,epsrel,limit,ws,&resultdn,&abserr);
                resultdn*=(1.0/(2.0*M_PI));

                F6.params = &pm;
                F6.function = &IntFuncDerivVg;
                gsl_integration_qags(&F6,ll,ul,epsabs,epsrel,limit,ws6,&resultdvg,&abserr);
                resultdvg*=(1.0/(2.0*M_PI));

                result = resultdvg / (1-resultdn);

                printf("%f %f %f %f\n", pm.ds, resultdvg, resultdn, result);
                fd4 = fopen("derivatives", "a");
                fprintf(fd4,"%f %f %f %f %f %f %f\n", pm.ds, resultdvg, resultdn, result, pm.ul, pm.ur, pm.Vg);
 	            fclose(fd4);

    	    }
        //----------------Green function-----------------------------------
        case 6:

            printf("\n Printing the spectral function \n");

            fd4 = fopen("green_function.dat", "w");
            fprintf(fd4, "#   x    G\n");
            fclose(fd4);

            ll = pm.a - 2.0*pm.b + pm.ur;   // lower limit of the integration
            ul = pm.a + 2.0*pm.b + pm.ul;   // upper limit if the integration

            printf("ds     dst \n");

            double x=5*ll;
            while(x<5*ul)
            {
                GreenFunc(x, &pm);
                x+=0.02;

            }

            /* for(i = 0; i < pm.nsol; i++){
             pm.ds = pm.sol[i];
             F.params = &pm;
             F.function = &SpectrumFunc;

             gsl_integration_qags(&F,ll,ul,epsabs,epsrel,limit,ws,&result,&abserr);
             result*=(1.0/(2.0*M_PI));

             printf("%f %f\n", pm.ds, result);

             }*/
            break;

            break;
    }

    fclose(fd4);
    fclose(log);
    gsl_integration_workspace_free(ws);
    gsl_integration_workspace_free(ws1);
    gsl_integration_workspace_free(ws2);
    gsl_integration_workspace_free(ws3);
    gsl_integration_workspace_free(ws4);
    gsl_integration_workspace_free(ws5);
    gsl_integration_workspace_free(ws6);
    gsl_root_fdfsolver_free(s);

    return status;
}


//----------------------------------------------------------------------------------
//  Prints the used approximation
//----------------------------------------------------------------------------------
void print_approx(char c) {
    FILE *fd;
    fd = fopen("output.dat","a");
    fprintf(fd,"\n");
    if (c == 'H')
    {
        printf("Hartree-Fock\n");
        fprintf(fd,"Hartree-Fock\n");
    }
    else if (c == 'B')
    {
        printf("2B\n");
        fprintf(fd,"2B\n");
    }
    else if (c == 'G')
    {
        printf("GW\n");
        fprintf(fd,"GW\n");
    }
    else if (c == 'T')
    {
        printf("T-matrix\n");
        fprintf(fd,"T-matrix\n");
    }
    else
    {
        printf("CHECK your approximation!!!!\n");
        fprintf(fd,"CHECK your approximation!!!!\n");
    }
    fclose(fd);
}
//----------------------------------------------------------------------------------


//EOF
