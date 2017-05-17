//Compiled in Matlab09:
//mex -lgsl -lgslcblas eiModel.cpp
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <time.h>
#include <stdio.h>
#include </usr/include/gsl/gsl_rng.h>
#include </usr/include/gsl/gsl_randist.h>
using namespace std;
gsl_rng *rng;

// model parameters
struct modelParams {
  double Jee;
  double Jei;
  double Jie;
  double Jii;
  double thetaE;
  double thetaI;
  double Eslope;
  double Edesp;
  double Islope;
  double Idesp;
  double deltaGE;
  double tauE;
  double tauI;
  double tauAdapt;
  double tauN;
  double sigmaN;
} ;


void allderivs(double* y, double* n, modelParams* fParams, double* dy)
{
    double aux; 
    //RE
    aux=fParams->Jee*y[0]-fParams->Jei*y[1]+fParams->thetaE-y[2]+n[0];
    if (aux<=fParams->Edesp)
        dy[0]=-y[0]/fParams->tauE;
    else
        dy[0]=(-y[0]+fParams->Eslope*(aux-fParams->Edesp))/fParams->tauE;
    //RI
    aux=fParams->Jie*y[0]-fParams->Jii*y[1]+fParams->thetaI+n[1];
    if (aux<=fParams->Idesp)
        dy[1]=-y[1]/fParams->tauI;
    else
        dy[1]=(-y[1]+fParams->Islope*(aux-fParams->Idesp))/fParams->tauI;
    //A
    dy[2]=(-y[2]+fParams->deltaGE*y[0])/fParams->tauAdapt;
}

void integrador(double val[], double noiseVal[], double dt, int nrVar, int nrSteps, modelParams* fParams)
{
    double nAux1=exp(-dt/fParams->tauN);
	double nAux2=sqrt(((2*pow(fParams->sigmaN,2)/fParams->tauN)*fParams->tauN*0.5)*(1-pow(exp(-dt/fParams->tauN),2)));
    double rk4Aux1=dt*0.500000000;
    double rk4Aux2=dt*0.166666666;

    int i,j,k=0,l=0,Tdt=(int)(1/dt);
    double aux1[3],aux2[3],aux3[3],aux4[3];
    
    double yAct[3];
    double noiseAct[2];
    
    yAct[0]=val[0];
	yAct[1]=val[1];
    yAct[2]=val[2];
    noiseAct[0]=noiseVal[0];
    noiseAct[1]=noiseVal[1];

	for (j=0;j<nrSteps-Tdt;j++){

        //RK4
		allderivs(yAct,noiseAct,fParams,aux1);
        for (i=0;i<nrVar;i++) 
            aux2[i]=yAct[i]+rk4Aux1*aux1[i];
        
        allderivs(aux2,noiseAct,fParams,aux3);
        for (i=0;i<nrVar;i++) 
            aux2[i]=yAct[i]+rk4Aux1*aux3[i];
        
        allderivs(aux2,noiseAct,fParams,aux4);
        for (i=0;i<nrVar;i++) {
            aux2[i]=yAct[i]+dt*aux4[i];
            aux4[i] += aux3[i];
        }
        
        allderivs(aux2,noiseAct,fParams,aux3);
        for (i=0;i<nrVar;i++)
            yAct[i] += rk4Aux2*(aux1[i]+aux3[i]+2.0*aux4[i]);
        
        //noise update
		noiseAct[0]=noiseAct[0]*nAux1+nAux2*gsl_ran_ugaussian(rng);
		noiseAct[1]=noiseAct[1]*nAux1+nAux2*gsl_ran_ugaussian(rng);
        
        k++;
        //save result
        if (k==Tdt)
        {
            l++;
            val[l*nrVar+0]=yAct[0];
            val[l*nrVar+1]=yAct[1];
            val[l*nrVar+2]=yAct[2];
            noiseVal[l*2+0]=noiseAct[0];
            noiseVal[l*2+1]=noiseAct[1];
            k=0;
        }
	}
}

// MAIN
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{  
    const gsl_rng_type * T;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    rng = gsl_rng_alloc (T);
    gsl_rng_set( rng, time(NULL) );
    
    //associate inputs
    double *fParamsIn,*simParamsIn,*initCondIn,*simOut,*noiseOut;
    simParamsIn = mxGetPr(mxDuplicateArray(prhs[0]));
    fParamsIn = mxGetPr(mxDuplicateArray(prhs[1]));
    initCondIn = mxGetPr(mxDuplicateArray(prhs[2]));
    
    //associate outputs
    int nrVars=3;//re,ri,ge
    int nrSteps=(int)((simParamsIn[1]-simParamsIn[0])/simParamsIn[2]);
    plhs[0] = mxCreateDoubleMatrix(nrVars,(int)(nrSteps*simParamsIn[2]),mxREAL); 
    simOut= mxGetPr(plhs[0]);
    //noise values
    plhs[1] = mxCreateDoubleMatrix(2,(int)(nrSteps*simParamsIn[2]),mxREAL); 
    noiseOut = mxGetPr(plhs[1]);
    
    modelParams fParams;
    fParams.Jee=fParamsIn[0];
    fParams.Jei=fParamsIn[1];
    fParams.Jie=fParamsIn[2];
    fParams.Jii=fParamsIn[3];
    fParams.thetaE=fParamsIn[4];
    fParams.thetaI=fParamsIn[5];
    fParams.Eslope=fParamsIn[6];
    fParams.Edesp=fParamsIn[7];
    fParams.Islope=fParamsIn[8];
    fParams.Idesp=fParamsIn[9];
    fParams.deltaGE=fParamsIn[10];
    fParams.tauE=fParamsIn[11];
    fParams.tauI=fParamsIn[12];
    fParams.tauAdapt=fParamsIn[13];
    fParams.tauN=fParamsIn[14];
    fParams.sigmaN=fParamsIn[15];
    
    simOut[0]=initCondIn[0];                //re0
    simOut[1]=initCondIn[1];                //ri0
    simOut[2]=initCondIn[2];                //ge0
    noiseOut[0]=0;                          //nE
    noiseOut[1]=0;                          //nI

    //integration
    integrador(simOut, noiseOut, simParamsIn[2], nrVars, nrSteps, &fParams);

    return;
}
