#include "num.h"

struct params{
    double a;
    double b;
    double beta;
    double V_l;
    double V_r;
    double mu;
    double e0;
    double ul;
    double ur;
    double ulmin;
    double ulmax;
    double urmin;
    double urmax;
    double ustep;
    double ds;
    double U_h;
    double Vg;
    double damp;
    int damptype;
    double eta;
    char c;
    double xx;
    int sweep;
    int val;
    int aval;
    int nsol;
    double *sol;
};

void read_params(struct params *pm);
double *vector(long nl, long nh);
void free_vector(double *v, long nl, long nh);
