#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <malloc.h>
//#include <boost/math/special_functions/bessel.hpp>
#include <complex.h>
//#include <gsl/gsl_integration.h>
//#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coupling.h>
#include <fftw3.h>
#include <vector>
#include <cstring>
#include "Headers/Headers.h"


int main (){
  /*  
    z = 0.5; 
    std::cout << "0 " << gsl_sf_bessel_J0(z) << std::endl;
    std::cout << "1 " << gsl_sf_bessel_J1(z) << std::endl;
    std::cout << "-1 " << gsl_sf_bessel_Jn(-1,z) << std::endl;
    for( i = 2; i <= 8; i++){
        std::cout << i << " " << gsl_sf_bessel_Jn(i,z) << std::endl;
        std::cout << -i << " " << gsl_sf_bessel_Jn(-i,z) << std::endl;
    }
*/	
    n_bins = 512;
	n_bins_2 = n_bins/2;
	inv_n = 1./n_bins;

	H = 20;
	H2 = H/2;
	
	RII = 0.5*sqrt(0.5);
    
	DeltaR = H*inv_n;
	DeltaK = inv_n*2*M_PI/DeltaR;
	DeltaT = M_PI/100;
	DeltaP = b1/100;

    inv_nDeltaR = inv_n*DeltaR;

    double V_pear = 0.596981045;//(4*M_PI*RII*RII*RII);
	rho_sum_o = 0.3*n_bins/V_pear;

    alpha = 0.01;

    Setup(); 
    
    CalcRho();

/*
    CalcW0();

    CalcW1();
    CalcW2();
    CalcW3();
   */ 
    CalcW();

    for( j=0; j<= 10; j++) w0[10*j].print(n_bins,DeltaR); 
    for ( i = 0; i <= lmax; i++)
        for(j=0; j<= 10; j++) w1[i][10*j].print(n_bins,DeltaR); 
    
    for ( i = 0; i <= lmax; i++)
        for(j=0; j<= 10; j++) w2[i][10*j].print(n_bins,DeltaR); 
   
    for(j=0; j<= 10; j++) w3[10*j].print(n_bins,DeltaR); 

    st =0;
    
    now = 0;
    prev = 1; 
    while( fabs(now-prev) > 1e-9 && st < 5000 ){
        for(j=0; j<= 100; j++) rho[j].fft();

        CalcN();
        
        
        CalcPHI();

        
        CalcSC();
    
        
        CalcC();
		/*
        if(st==0){
            for( j=0; j<= 100; j++) rho[j].print(n_bins,DeltaR); 
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
        */

        CalcRhoN();

        st++;
    }    
    
    for( j=0; j <= 100; j++){
        std::string js = toString(j); 
        rho[j].Filename = "rho_fin" + js + ".dat"; 
        rho[j].printReal(n_bins, DeltaR);
    }
   
    for( iz = 0; iz < n_bins; iz++){
        rho_fin.real[iz] = 0.5*(rho[0].real[iz]+rho[100].real[iz]);

        for ( j = 1; j < 100; j++){
            rho_fin.real[iz] += rho[j].real[iz];
        }
        rho_fin.real[iz] *= 0.02*V_pear;
    }
    rho_fin.print(n_bins, DeltaR);
}
