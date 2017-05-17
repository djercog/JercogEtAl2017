//==================================================================================================
// Daniel Jercog
//
//   To compile:
// 	    $g++ *.cpp -O3 -lm -o execSim
//
//   To execute the code model, 2 additional input parameters are required:
//      - To use Fig 7 parameters set value to 1 (set to 0 for Fig 7 Supp 2)
//      - To save rasters, second parameter in the call must be set to 1 (0 will not generate raster files)
//
//   Example call:
//	    $./execSim 1 1
//
//==================================================================================================
#include <iostream>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include "functions.h"
using namespace std;

int main(int argc, char *argv[]){

	int isKicked=atoi(argv[1]);         // Parameters used flag: 1 = Fig7 , 0 = Fig7Supp2.
	int saveRasters= atoi(argv[2]);  	// Save raster flag: 1 = save, 0 = do not save.

	// MODEL PARAMETERS FIG 7 ============================================================================
	double Ii=6.500000;					// Mean current to inhibitory neurons (I)
	double Ie=7.600000;					// Mean current to excitatory neurons (E)
	double sigma_e=2.5000000;			// Std of E membrane potential
	double sigma_i=2.5000000;			// Std of I membrane potential
	double commonNoise=220;				// Size of the external fluctuations		
	double lambdaKicks=400; 		 	// Mean isi of kicks (ms)
	double kickLength=2.0;              // Duration of external kicks (ms)
	double percKicked=0.10000;			// % of cells kicked by external fluctuation
	double upKickFactor=1.0/2.0;		// kicks during UP states are reduced in amplitude by this factor
	double iniSeed=time(NULL);			// initial seed
	double timeTot=100000;				// simulation length (ms)
	int coutMaxTm=100000;				// maximum time showing cout (ms)		
	bool saveOutputRaster=saveRasters;	// save rasters for E and I populations (set nr of cells saved below)
	bool coutV=0;		        		// print in console several values (Vms, Adapt, etc)		

	// model parameters
	int Ne = 4000;						// nr of excitatory cells
	int Ni = 1000;						// nr of inhibitory cells
	double tau_e = 20.000000;			// time constant of excitatory membrane potential
	double tau_i = 10.000000;			// time constant of inhibitory membrane potential
	double Vr = 14.000000;				// reset membrane potential (after spike; no refractory period)
	double Vth = 20.000000;				// membrane potential spike-threshold 
	double Jee = 280.0/Ne; 				// E to E connectivity
	double Jei = 70.0/Ni; 				// I to E connectivity	
	double Jie = 500.0/Ne; 				// E to I connectivity
	double Jii = 100.0/Ni;				// I to I connectivity	
	double Ja = 15.0;					//adaptation strenght: \beta=Ja*tau_e

	//time constants
	double decay_e = 23;				// EPSC decay time constant 
	double rise_e = 8;					// EPSC rise time constant 
	double decay_i = 1;					// IPSC decay time constant 
	double rise_i = 1;					// IPSC rise time constant 
	double tau_a = 500;					// adaptation decay time constant 

    // Parameters for Fig7Supp2
	if (!isKicked) {
        Ii=12.500;							// Mean current to inhibitory neurons (I)
    	Ie=13.2000;							// Mean current to excitatory neurons (E)
    	commonNoise=0;						// Size of the external fluctuations		
    	Ja = 10.0;							// adaptation strenght: \beta=Ja*tau_e
    	decay_e = 3;						// EPSC decay time constant 
    	rise_e = 2;							// EPSC rise time constant 
     }

	// integration times
	double dt = 0.05, delay = 1;		//step for integration, maximum delay for synapses
	int countBin=10;                    //bin size in ms for spike counts.


	// Internal variables ============================================================================
	int seed = iniSeed; 				// seed for the random number generator //int seed = 987654; //(unsigned int) time(NULL);  
	int nDelay = round(delay/dt), nTime = round(timeTot/dt);
	double ge=(sigma_e/sqrt(tau_e))*sqrt(dt)*(1.0-(dt/(2.0*tau_e)));	// stochastic terms of dV/dt
	double gi=(sigma_i/sqrt(tau_i))*sqrt(dt)*(1.0-(dt/(2.0*tau_i)));	
	double parXi=exp(-dt/rise_i), parXe=exp(-dt/rise_e);
	double parSe=exp(-dt/decay_e), parSi=exp(-dt/decay_i);
	double parIa = exp(-dt/tau_a);
	double fe,fi,commNoise,IaCumsum,VEAverage,VIAverage,kickIntensity;
	int i,j,spikeTrackTime,auxCount;
	int kickStepNr=round(kickLength/dt),kickStepCount=0;	//kickStepNr: length of kick in nr of steps
	
	//heterogeneous delays: delay matrices for E and I cells:
	double xi = 0, xe = 0, se = 0, si = 0;
	bool spikeTrackE[Ne], spikeTrackI[Ni];		// tells if cell is firing at current timestamp (I use bool as int later)
	bool spksENDelayed[nDelay][Ne];				// this guys below save the history of spiking activity in the last delay ms
	for (i=0; i<nDelay; i++) {
		for (j=0; j<Ne; j++) {
			spksENDelayed[i][j]=0;
		}
	}
	bool spksINDelayed[nDelay][Ni];
	for (i=0; i<nDelay; i++) {
		for (j=0; j<Ni; j++) {
			spksINDelayed[i][j]=0;
		}
	}

	int matDelaySize=1000;
	int matDelayE[matDelaySize];
	int matDelayI[matDelaySize];
	for (i=0; i<matDelaySize; i++) {
		matDelayE[i]=(int)((double) ((1.0+nDelay) *r8_uniform_01 (&seed)));
		matDelayI[i]=(int)((double) ((1.0+(nDelay/2.0)) *r8_uniform_01 (&seed)));	//delayI U(0,0.5) ms
	}

	// initial conditions for voltage and adaptation
	double Ve[Ne], Vi[Ni], Ia[Ne];//, IaCumsum;
	for (i = 0; i < Ne; i++) {
		Ve[i] = Vr + (Vth-Vr)*r8_uniform_01 ( &seed );
		Ia[i] = 0;
	}
	for (i = 0; i < Ni; i++) {
		Vi[i] = Vr + (Vth-Vr)*r8_uniform_01 ( &seed );
	}

	//kicks: size fixed, arrival as poisson process
	double commNoiseTime[50+(int)(timeTot/lambdaKicks)][2];
	int tmKick=round(r8_exp_01(&seed,lambdaKicks))/dt;
	commNoiseTime[0][0]=tmKick;
	commNoiseTime[0][1]=commonNoise;	//*r8_exp_01(&seed,1.0);
	int nrKicks=0,itKick=0;				//iterator for main LOOP
	while (tmKick<timeTot/dt) {
		nrKicks++;		
		tmKick=round(r8_exp_01(&seed,lambdaKicks))/dt + commNoiseTime[nrKicks-1][0];
		commNoiseTime[nrKicks][0]= tmKick;
		commNoiseTime[nrKicks][1]= commonNoise;
	}
	
	double currRateE=0;		// rate of the network in the previous count bin
	int spksE,spksI;		// spkes for each timesteps (used in COUT)

	// Saving output  ============================================================================
	//rate
	ostringstream strBaseName;
	strBaseName.str("");
	strBaseName << "Ie" << Ie*1000 << "Ii" << Ii*1000 << "sig" << sigma_e*10 << "coN" << commonNoise*10 << "Ja" << (int)Ja*10;
	int spikebinE=0,spikebinI=0; //accumulator for pop rate.
	int sizeBuffRate = 5000;
	sizeBuffRate=sizeBuffRate/countBin;		// <- find sweet spot for buff size
	countBin=round(countBin/dt);			//change variable to step number
	coutMaxTm=round(coutMaxTm/dt);			//same...
	ostringstream tmpStr;
	char fileRateE[50];tmpStr.str("");
	tmpStr << strBaseName.str().c_str() << "countE.bin";
	strncpy(fileRateE, tmpStr.str().c_str(), sizeof(fileRateE));
	fileRateE[sizeof(fileRateE) - 1] = 0;
	char fileRateI[50];tmpStr.str("");
	tmpStr << strBaseName.str().c_str() << "countI.bin";
	strncpy(fileRateI, tmpStr.str().c_str(), sizeof(fileRateI));
	fileRateI[sizeof(fileRateI) - 1] = 0;
	char fileAda[50];tmpStr.str("");
	tmpStr << strBaseName.str().c_str() << "Ada.bin";
	strncpy(fileAda, tmpStr.str().c_str(), sizeof(fileAda));
	fileAda[sizeof(fileAda) - 1] = 0;

	int buffRateIt=0;
	FILE* pFile;
	int outRateEBuff[sizeBuffRate],outRateIBuff[sizeBuffRate],outAdaBuff[sizeBuffRate];
	//kicks (NOTE: this are NOT the kicks given, since now I ve added the rule of min E rate to kick the network)
	tmpStr.str("");
	tmpStr << strBaseName.str().c_str() << "kicks.bin";
	pFile = fopen(tmpStr.str().c_str(), "w+b");
	for (i=0;i<nrKicks;i++){
		fwrite(&commNoiseTime[i][0],sizeof(double),1,pFile);
		fwrite(&commNoiseTime[i][1],sizeof(double),1,pFile);
	}
	fclose(pFile);
	//rasters 
	int nrCellsE=800, nrCellsI=200;	//numero de cells para salvar el raster
	vector<int> rastE,rastI;
	vector<int> rastEId,rastIId;

	//for UD detection
	int spikebin2=0, buffCountIt2=0;
	int countBin2=round(10.0/dt);	
	int sizeBuffRate2 = 500; sizeBuffRate2=sizeBuffRate2/countBin;
	int outCount[sizeBuffRate2];
	char fileCount[50];tmpStr.str("");
	tmpStr << strBaseName.str().c_str() << "countdet.bin";
	strncpy(fileCount, tmpStr.str().c_str(), sizeof(fileCount));
	fileCount[sizeof(fileCount)-1] = 0;


	// MAIN LOOP ============================================================================
	for(int timestep=0; timestep < nTime; timestep++) {

		//donde estoy actualmente en el vector de delays
		spikeTrackTime = timestep % nDelay;

		// reset spikeTracks of current time to zero
		for (i = 0; i < Ne; i++) {
			spikeTrackE[i] = 0;
		}
		for (i = 0; i < Ni; i++) {
			spikeTrackI[i] = 0;
		}	

		//new kick arrived!
		if (timestep>commNoiseTime[itKick][0]){
			if (currRateE<0.5){	//only kick if rateE is less than 0.5Hz
				kickIntensity=1.0;				
			} else {
				kickIntensity=upKickFactor;
			}
			kickStepCount=kickStepNr;		
			itKick++;						
		}
		if (kickStepCount>0){	//if there are kick steps left...
 			commNoise=commNoiseTime[itKick-1][1]*(1-exp(-(1.0+kickStepNr-kickStepCount)/10.0))*kickIntensity;
			kickStepCount--;
		} else {
			commNoise=0.0;
		}
		// calculate Ve and Ia for exc neurons
		spksE=0;spksI=0;
		IaCumsum=0;
		VEAverage=0;
		for (i = 0; i < Ne; i++) {
		  if (i<=round(percKicked*Ne)){
		     fe = (-Ve[i] + tau_e*Jee*se - tau_e*Jei*si - tau_e*Ja*Ia[i]  + Ie + commNoise)/tau_e;	
		  } else {
		     fe = (-Ve[i] + tau_e*Jee*se - tau_e*Jei*si - tau_e*Ja*Ia[i]  + Ie)/tau_e;	
		  }
		  //}Integration RK2, receives commNoise in case it has to be kicked
		  Ve[i] = Ve[i] + fe * dt * (1. - dt/(2.0*tau_e)) + ge * r8_normal_01(&seed);
		  
		  // calculate new Ia, using exact formula
		  Ia[i] = Ia[i]*parIa;

		  IaCumsum+=Ia[i];
		  VEAverage+=Ve[i];

		  // when reaching theshold
		  if (Ve[i] >= Vth) {
			spksE++;
			Ve[i] = Vr; 			// reset Ve
			Ia[i] += 1./tau_a; 		// adaptive current increases for particular neuron
			spikeTrackE[i] = 1;		//spikeTrackE[spikeTrackTime]++; // add spike to spikeTrack
			spikebinE++; 			// add spike to spikebinE for rate calc
			if (saveOutputRaster) {		//antes salvaba el raster hasta coutMaxTm, ahora lo salvo entero
				if (i<nrCellsE){
					rastE.push_back(timestep);
					rastEId.push_back(i);
				}
			}
			//para la deteccion de UD
			// E contribuye 90 cells, 9 are kicked
			if (i>=991 && i<=1080) {
			  spikebin2++;
			}
		  }
		}

		// calculate Vi for inh neurons
		VIAverage=0;
		for (i = 0; i < Ni; i++) {
		  if (i<=round(percKicked*Ni)){
			fi = (-Vi[i] + tau_i*Jie*se - tau_i*Jii*si + Ii + 0.4*commNoise)/tau_i;
		  } else {		   
			 fi = (-Vi[i] + tau_i*Jie*se - tau_i*Jii*si + Ii)/tau_i;
		  }
		  Vi[i] = Vi[i] + fi * dt * (1. - dt/(2.0*tau_i)) + gi * r8_normal_01(&seed);
		
		  VIAverage+=Vi[i];

		  // when reaching theshold
		  if (Vi[i] >= Vth) {
			  spksI++;
			  Vi[i] = Vr; 			// reset Vi
			  spikeTrackI[i] = 1;	//spikeTrackI[spikeTrackTime]++; // add spike to spikeTrack
			  spikebinI++; 			// add spike to spikebinI
			  if (saveOutputRaster) {
				  if (i<nrCellsI){
					  rastI.push_back(timestep);
					  rastIId.push_back(i);
				  }
			  }
			  //para la deteccion de UD
			  // I contribuye 10 cells, 1 is kicked
			  if (i>=99 && i<=108) {
			    spikebin2++;
			  }
		  }
		} // end of processing spike of inh neuron

		// console:
		if (coutV && (timestep<=coutMaxTm)){
			cout <<Ve[389]<<", "<<Vi[89]<<", "<<se<<", "<<si<<", "<<Ie-tau_e*Ja*IaCumsum/Ne <<", "<<VEAverage/(1.0*Ne)<<", "<<VIAverage/(1.0*Ni)<<", "<<Ve[409]<<", "<<Vi[109]<<", "<<spksE<<", "<<spksI<<", "<<commNoise<<'\n';				//<<<<<< COUT RES
		}	

		// calculate se and si according to exact formula for E cells, FOR EACH CELL
		se = xe + (se - xe) * parSe;
		si = xi + (si - xi) * parSi;	
		// calculate xe and xi according to exact formula
		// i.e., exact solution to ODE of variable xe of exc neurons: (rise_e * dxe/dt = - xe) => xe = xe *exp(-dt/rise_e);
		xe = xe*parXe;
		xi = xi*parXi;
		// xe and xi increase if there were spikes nDelay timesteps ago
		auxCount=0;
		for (j=0;j<Ne;j++){
			auxCount+=spksENDelayed[(spikeTrackTime-matDelayE[j%matDelaySize]+nDelay)%nDelay][j];
		}
		xe += auxCount/rise_e;
		auxCount=0;
		for (j=0;j<Ni;j++){
			auxCount+=spksINDelayed[(spikeTrackTime-matDelayI[j%matDelaySize]+nDelay)%nDelay][j];
		}
		xi += auxCount/rise_i;

		//update spkcount for current step
		for (i=0; i<Ne; i++) {
			spksENDelayed[spikeTrackTime][i]=spikeTrackE[i];
		}
		for (i=0; i<Ni; i++) {
			spksINDelayed[spikeTrackTime][i]=spikeTrackI[i];
		}

		//log the events if necesary
		if ((timestep+1) % countBin == 0){								//end of bin
			currRateE=(spikebinE/(1.0*Ne))*(1000.0/(countBin*dt));		//update current rate 				
			outRateEBuff[buffRateIt]=spikebinE;							//update output arrays
			outRateIBuff[buffRateIt]=spikebinI;
			//outAdaBuff[buffRateIt]=round(1000.0*(Ie-tau_e*Ja*IaCumsum/Ne));
            outAdaBuff[buffRateIt]=round(1000.0*(tau_e*Ja*IaCumsum/Ne));
			spikebinE=0,spikebinI=0,buffRateIt++;						//reset counters and move forward iterator
			if (buffRateIt==sizeBuffRate){								//if buffer is full, write and reset buffer iterator
				pFile = fopen(fileRateE, "ab");
				fwrite(outRateEBuff, 1, sizeBuffRate*sizeof(int), pFile);
				fclose(pFile);
				pFile = fopen(fileRateI, "ab");
				fwrite(outRateIBuff, 1, sizeBuffRate*sizeof(int), pFile);
				fclose(pFile);
				pFile = fopen(fileAda, "ab");
				fwrite(outAdaBuff, 1, sizeBuffRate*sizeof(int), pFile);
				fclose(pFile);
				buffRateIt=0;
			}
		}

		//log the events if necesary
		if ((timestep+1) % countBin2 == 0){		
			outCount[buffCountIt2]=spikebin2;
			spikebin2=0,buffCountIt2++;						//reset counters and move forward iterator
			if (buffCountIt2==sizeBuffRate2){								//if buffer is full, write and reset buffer iterator
				pFile = fopen(fileCount, "ab");
				fwrite(outCount, 1, sizeBuffRate2*sizeof(int), pFile);
				fclose(pFile);
				buffCountIt2=0;
			}
		}


	} // end of main routine

	//if saveRasters, write files...
	if (saveOutputRaster) {
		char fileRasterE[100];tmpStr.str("");
		tmpStr << strBaseName.str().c_str() << "rasterE.bin";
		strncpy(fileRasterE, tmpStr.str().c_str(), sizeof(fileRasterE));
		fileRasterE[sizeof(fileRasterE) - 1] = 0;
		char fileRasterI[100];tmpStr.str("");
		tmpStr << strBaseName.str().c_str() << "rasterI.bin";
		strncpy(fileRasterI, tmpStr.str().c_str(), sizeof(fileRasterI));
		fileRasterI[sizeof(fileRasterI) - 1] = 0;
		pFile = fopen(fileRasterE, "w+b");
		for (i=0;i<(int)rastE.size();i++){
			fwrite(&rastE[i],sizeof(int),1,pFile);
			fwrite(&rastEId[i],sizeof(int),1,pFile);
		}
		pFile = fopen(fileRasterI, "w+b");
		for (i=0;i<(int)rastI.size();i++){
			fwrite(&rastI[i],sizeof(int),1,pFile);
			fwrite(&rastIId[i],sizeof(int),1,pFile);
		}
		fclose(pFile);
	}

	return 0;

}
