#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "onds.h"

void read_params(struct params *pms){
    
	double eta, mu, beta, gamma;
	double V_l, V_r, ul, ur;
	double a,b,a0, U_h, Vg;
	char c;
	int sweep, val,aval, gammatype, nsol, i;
	double ulmin, ulmax, urmin, urmax, ustep;
	
        FILE *out;
	FILE *fname;

	out = fopen("output.dat","a");
	fprintf(out, "Parameters:\n\n");
	fname = fopen("input.in", "r");
    
	// -----------------------------------------
    
	c = getc(fname); while(c!='\n') c = getc(fname);
	c = getc(fname); while(c!='\n') c = getc(fname);
    
	fscanf(fname,"%lf", &eta);
	pms->eta = eta;
	printf("read eta..... eta = %f\n", pms->eta);
	fprintf(out,"eta = %f\n", pms->eta);
    
	c = getc(fname); while(c!='\n') c = getc(fname);
    
	fscanf(fname, "%lf", &mu); 
	pms->mu = mu;
	printf("read mu..... mu = %f\n", pms->mu);
      	fprintf(out, "mu = %f\n", pms->mu);    

	c=getc(fname); while(c!='\n') c=getc(fname);
    
	fscanf(fname, "%lf", &beta);
	pms->beta = beta;
	printf("read beta..... beta = %f\n", pms->beta);
	fprintf(out,"beta = %f\n",pms->beta);    

	c = getc(fname); while(c!='\n') c = getc(fname);
    
	fscanf(fname, "%lf", &V_l);
	pms->V_l = V_l;
	printf("read V_l..... V_l = %f\n", pms->V_l);
	fprintf(out, "V_l = %f\n", pms->V_l);
    
	c=getc(fname); while(c!='\n') c=getc(fname);
    
	fscanf(fname, "%lf", &V_r);
	pms->V_r = V_r;
	printf("read V_r..... V_r = %f\n", pms->V_r);
	fprintf(out,"V_r = %f\n", pms->V_r);
    
	// -----------------------------------------
	// starting to read central region parameters
	// -----------------------------------------
    	fprintf(out, "\nCentral region parameters: \n");

	c = getc(fname); while(c!='\n') c = getc(fname);
	c = getc(fname); while(c!='\n') c = getc(fname);
	c = getc(fname); while(c!='\n') c = getc(fname);
    
	fscanf(fname, "%lf", &a0);
	pms->e0 = a0;
	printf("read e0..... e0 = %f\n", pms->e0);
	fprintf(out, "e0 = %f\n", pms->e0);
    
	c = getc(fname); while(c!='\n') c = getc(fname);
    
	fscanf(fname, "%lf", &U_h);
	pms->U_h = U_h;
	printf("read U_h..... U_h = %f\n", pms->U_h);
	fprintf(out," U_h = %f\n", pms->U_h);
    
	c = getc(fname); while(c!='\n') c = getc(fname);
    
	fscanf(fname, "%lf", &Vg);
	pms->Vg = Vg;
	printf("read Vg..... V_g = %f\n", pms->Vg);
	fprintf(out,"V_g = %f\n", pms->Vg);

	// -----------------------------------------    
	// starting to read lead parameters
	// -----------------------------------------
    	fprintf(out, "\nLead parameters:\n");

	c = getc(fname); while(c!='\n') c = getc(fname);
	c = getc(fname); while(c!='\n') c = getc(fname);
	c = getc(fname); while(c!='\n') c = getc(fname);
    
	fscanf(fname, "%d", &aval);
	pms->aval = aval;
	printf("read aval..... aval = %d\n", pms->aval);
	fprintf(out, "aval = %d\n", pms->aval);
    
	c = getc(fname); while(c!='\n') c = getc(fname);
    
	fscanf(fname, "%lf", &a);
	pms->a = a;
	printf("read a.....a = %f\n", pms->a);
	fprintf(out,"a = %f\n", pms->a);
    
	c = getc(fname); while(c!='\n') c = getc(fname);
    
	fscanf(fname, "%lf", &b);
	pms->b = b;
	printf("read b..... b = %f\n", pms->b);
	fprintf(out, "b = %f", pms->b);    

	c=getc(fname); while(c!='\n') c=getc(fname);
    
	fscanf(fname, "%lf", &ul);
	pms->ul = ul;
	printf("read U_l..... U_l = %f\n", pms->ul);
	fprintf(out, "U_l = %f", pms->ul);    

	c = getc(fname); while(c!='\n') c = getc(fname);
    
	fscanf(fname, "%lf", &ur);
	pms->ur = ur;
	printf("read U_r..... U_r = %f\n", pms->ur);
	fprintf(out, "U_r = %f\n", pms->ur);
    
	// -----------------------------------------
	fprintf(out, "\n");
    
	c = getc(fname); while(c!='\n') c = getc(fname);
	c = getc(fname); while(c!='\n') c = getc(fname);
	c = getc(fname); while(c!='\n') c = getc(fname);
    
	fscanf(fname, "%d", &val);
	pms->val = val;
	printf("read val..... val = %d\n", pms->val);
	fprintf(out, "val = %d\n", pms->val);  

	c = getc(fname); while(c!='\n') c = getc(fname);
    
	fscanf(fname, "%c", &c);
	pms->c = c;
	printf("read c..... c = %c\n", pms->c);
	fprintf(out, "c = %c\n", pms->c);
    
	c = getc(fname); while(c!='\n') c = getc(fname);
    
	fscanf(fname, "%d", &gammatype);
	pms->damptype = gammatype;
	printf("read gammatype..... gammatype = %d\n", pms->damptype);
	fprintf(out, "gammatype = %d\n", pms->damptype);
    
	c = getc(fname); while(c!='\n') c = getc(fname);
    
	fscanf(fname, "%lf", &gamma);
	pms->damp = gamma;
	printf("read gamma..... gamma = %f\n", pms->damp);
	fprintf(out,"gamma = %f\n", pms->damp);
    
	// -----------------------------------------
    
	fprintf(out, "\n");

	c=getc(fname); while(c!='\n') c=getc(fname);
	c=getc(fname); while(c!='\n') c=getc(fname);
	c=getc(fname); while(c!='\n') c=getc(fname);
    
	fscanf(fname, "%d", &sweep);
	pms->sweep = sweep;
	printf("read sweep..... sweep = %d\n", pms->sweep);
	fprintf(out, "sweep = %d\n", pms->sweep);
    
	c=getc(fname); while(c!='\n') c=getc(fname);
    
	fscanf(fname,"%lf %lf", &ulmin, &ulmax);
	pms->ulmin = ulmin;
    	pms->ulmax = ulmax;
    	printf("read ulmin, ulmax... ulmin = %f, ulmax = %f\n",pms->ulmin, pms->ulmax);
    	fprintf(out,"ulmin = %f, ulmax = %f\n",pms->ulmin, pms->ulmax);
    
    	c=getc(fname); while(c!='\n') c=getc(fname);
    
    	fscanf(fname,"%lf %lf", &urmin, &urmax);
    	pms->urmin = urmin;
    	pms->urmax = urmax;
    	printf("read urmin, urmax... urmin = %f, urmax = %f\n",pms->urmin, pms->urmax);
    	fprintf(out,"urmin = %f, urmax = %f\n",pms->urmin, pms->urmax);
    
   	c = getc(fname); while(c!='\n') c = getc(fname);
    
	fscanf(fname, "%lf", &ustep);
	pms->ustep = ustep;
	printf("read ustep..... ustep = %f\n", pms->ustep);
	fprintf(out,"ustep = %f\n", pms->ustep);
    
	// -----------------------------------------
    
    	c = getc(fname); while(c!='\n') c = getc(fname);
    	c = getc(fname); while(c!='\n') c = getc(fname);
    	c = getc(fname); while(c!='\n') c = getc(fname);
    
	fscanf(fname, "%d", &nsol);
	pms->nsol = nsol;
	printf("read nsol..... nsol = %d\n", pms->nsol);
	fprintf(out,"nsol = %d\n", pms->nsol);
	
  	 c = getc(fname); while(c!='\n') c = getc(fname);
    
	double sol[nsol];	
	pms->sol = vector(0,nsol);
	
	for (i = 0; i < nsol; i++)
	{
		fscanf(fname, "%lf", &sol[i]);
        	pms->sol[i] = sol[i];
		printf("read %f\n", pms->sol[i] );
		c=getc(fname);
	}
	
	fclose(out);
	fclose(fname);
    
}



