#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <malloc.h>
#include <complex.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
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
	DeltaK = inv_n*2*M_PI/DeltaR;
	DeltaT = M_PI/100;
	DeltaP = b1/50;

    inv_nDeltaR = inv_n*DeltaR;

	bulk = 0.3*H/((H-2*RII)*4*M_PI*RII*RII*RII *M_PI);

    alpha = 0.1;

    Setup(); 
    
    CalcRho();

    CalcW0();
       
    w0[0].print(n_bins,DeltaR); 
    
    CalcW1();
    
    w1[0][0].print(n_bins,DeltaR); 
    
    CalcW2();
    
    w2[0][0].print(n_bins,DeltaR); 
    
    CalcW3();
    
    w3[0].print(n_bins,DeltaR); 

    st =0;
    
    now = 0;
    prev = 1; 
    while( fabs(now-prev)>1e-8  ){
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
            sc0[0].print(n_bins,DeltaR); 
            sc1[0][0].print(n_bins,DeltaR); 
            sc1[0][20].print(n_bins,DeltaR); 
            sc2[0][0].print(n_bins,DeltaR); 
            sc2[0][20].print(n_bins,DeltaR); 
            sc3[0].print(n_bins,DeltaR); 
            c1[0].print(n_bins,DeltaR); 
            c1[20].print(n_bins,DeltaR); 
        }

        CalcRhoN();
        st++;
    }    
    
    for( j=0; j <= 100; j++){
        std::string js = toString(j); 
        rho[j].Filename = "rho_fin" + js + ".dat"; 
        rho[j].printReal(n_bins, DeltaR);
    }
   
    for( iz = 0; iz < n_bins; iz++){
        rho[0].real[iz] = 0.5*(rho[0].real[iz]+rho[100].real[iz]);

        for ( j = 1; j < 100; j++){
            rho[0].real[iz] += rho[j].real[iz];
            rho[0].real[iz] += rho[j].real[iz];
        }
        rho[0].real[iz] *= DeltaT;
    }
    rho[0].Filename = "rho_fin.dat";
    rho[0].printReal(n_bins, DeltaR);
}
