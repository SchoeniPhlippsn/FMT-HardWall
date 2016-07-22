#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <malloc.h>
#include <complex.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coupling.h>
#include <fftw3.h>
#include <vector>
#include <cstring>
#include "Headers/Headers.h"


int main (){

	n_bins = 1024;
	n_bins_2 = n_bins/2;
	inv_n = 1./n_bins;

	H = 20;
	H2 = H/2;
	
	RII = 0.5*sqrt(0.5);

	DeltaR = H*inv_n;
	DeltaK = M_PI/100;

    inv_nDeltaR = inv_n*DeltaR;

	bulk = 0.3*H/((H-2*RII)*4*M_PI*M_PI*RII*RII*RII);

	rho_sum_o = (H-2*RII)*bulk*n_bins/H;
    
    alpha = 0.001;

    Setup(); 
    
    CalcRho();

    CalcW0();
        
    w0.print(n_bins,DeltaR); 

    CalcW1();
    
    w1[0][0].print(n_bins,DeltaR); 
    w1[0][20].print(n_bins,DeltaR); 
    
    CalcW2();
    
    w2[0][0].print(n_bins,DeltaR); 
    w2[0][20].print(n_bins,DeltaR); 
    
    CalcW3();
    
    w3.print(n_bins,DeltaR); 
          
    for( st = 0 ; st <= 10000; st++){
        for(j=0; j<= 100; j++) rho[j].fft();


        CalcN();
        
        
        CalcPHI();

        
        CalcSC();
    
        
        CalcC();
		
        if(st==0){
            rho[0].print(n_bins,DeltaR); 
            n0.print(n_bins,DeltaR); 
            n1[0].print(n_bins,DeltaR); 
            n2[0].print(n_bins,DeltaR); 
            n3.print(n_bins,DeltaR); 
            PHI0.print(n_bins,DeltaR); 
            PHI1[0].print(n_bins,DeltaR); 
            PHI2[0].print(n_bins,DeltaR); 
            PHI3.print(n_bins,DeltaR); 
            sc0.print(n_bins,DeltaR); 
            sc1[0][0].print(n_bins,DeltaR); 
            sc1[0][20].print(n_bins,DeltaR); 
            sc2[0][0].print(n_bins,DeltaR); 
            sc2[0][20].print(n_bins,DeltaR); 
            sc3.print(n_bins,DeltaR); 
            c1[0].print(n_bins,DeltaR); 
            c1[20].print(n_bins,DeltaR); 
        } 

        for(j=0; j<= 100; j++){
            
            rho_sum = 0;
            
            theta = DeltaK*j;

            b2 = WallDistance(theta);
            
            theta = DeltaK*(100-j);
            
            b3 = WallDistance(theta);

            for ( iz = 0; iz < n_bins; iz++){
                z = DeltaR*(iz+0.5);
                if(z <= H2){ 
                    if( z <= b2 ) rhon[j].real[iz] = 0.;
                    else rhon[j].real[iz] = exp(creal(c1[j].real[iz]));// + rho_sum[(ix]));
                }else{
                    z = H - z;
                    if( z <= b3 ) rhon[j].real[iz] = 0.;
                    else rhon[j].real[iz] = exp(creal(c1[j].real[iz]));// + rho_sum[(ix]));
                }
                rho_sum += creal(rhon[j].real[iz]);
            }
		
            sclf = rho_sum_o/rho_sum;
            if(st%1==0 && j== 100 ) printf("%i %f %f\n",st, rho_sum, sclf);

            for( iz =0; iz < n_bins; iz++){
                rhon[j].real[iz] *= sclf;
                rho[j].real[iz] = rho[j].real[iz]*(1-alpha) + rhon[j].real[iz]*alpha;     
            }
        }
    }    
    
    FILE * oFile;
    FILE * oFile2;

    oFile = fopen ("Results/rho_fin80.dat","w");
    oFile2 = fopen ("Results/rho_fin20.dat","w");


    for( iz =0; iz < n_bins; iz++){
        fprintf (oFile, "%f %f %f\n", DeltaR*(iz+0.5), creal(rho[80].real[iz]), cimag(rho[80].real[iz]));
        fprintf (oFile2, "%f %f %f\n", DeltaR*(iz+0.5), creal(rho[20].real[iz]), cimag(rho[20].real[iz]));
    }
    
}
