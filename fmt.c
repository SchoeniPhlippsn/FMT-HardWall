//g++ -o fmt -lgsl -lgslcblas -lm -lfftw3 fmt.c -std=gnu++98

#include <iostream>
#include <fstream>
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

    n_bins = 200;
	n_bins_2 = n_bins/2;
	inv_n = 1./n_bins;

	H = 20;
	H2 = H/2;
	
	RII = 0.5;
    
	DeltaR = H*inv_n;
	DeltaK = inv_n*2*M_PI/DeltaR;
	DeltaT = M_PI/100;
	DeltaP = b1/100;

    inv_nDeltaR = inv_n*DeltaR;

    V_pear = 2.65072*a1*a1*b1+0.687223*a1*a2*b1+0.147262*a2*a2*b1+0.687223*a1*a3*b1+0.147262*a3*a3*b1;// Volume Pear
    //V_pear = 0.596981045;// Volume Pear 3.8
    //V_pear = asp*(4*M_PI*RII*RII*RII)/3.0;//Volume Ellipsoid
    //V_pear = (4*M_PI*RII*RII*RII)/3;// Volume Sphere
	
    rho_sum_o = 0.3*n_bins/V_pear;

    alpha = 0.01;

    Setup(); 
    
    CalcRho();

    CalcW();
/*    for(j=0; j<= 1; j++){
        w0[j].read(n_bins, "fftw0_0.dat");
        w3[j].read(n_bins, "fftw3_0.dat");
        for ( i = 0; i <= lmax; i++){
            w1[i][j].read(n_bins,"fftw1-" + toString(i) + "_0.dat");
            w2[i][j].read(n_bins,"fftw2-" + toString(i) + "_0.dat");
        }
        std::cout << "             W[" << j << "]" << std::endl;
    }*/


    st =0;
    
    now = 0;
    prev = 1; 
    while( fabs(now-prev) > 1e-9 || st < 1000 ){
        if(st == 10000) break;
        for(j=0; j<= 100; j++) rho[j].fft();

        CalcN();
        
        CalcPHI();

        
        CalcSC();
    
        
        CalcC();
		
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
        
        //P1.real[iz]= 0.5*(-creal(rho[0].real[iz])+creal(rho[100].real[iz]));
        P1.real[iz]= 1*(-creal(rho[0].real[iz])+creal(rho[100].real[iz]));
        P2.real[iz]= 1*(creal(rho[0].real[iz])+creal(rho[100].real[iz]));

        for ( j = 1; j < 100; j++){
            //atheta = j*M_PI/100;
            //theta = cos(atheta);
            theta = j*0.02-1;
            rho_fin.real[iz] += rho[j].real[iz];
            P1.real[iz] += theta*creal(rho[j].real[iz]);
            P2.real[iz] += 0.5*(3*theta*theta-1)*creal(rho[j].real[iz]);
        }
        if( creal(rho_fin.real[iz]) < 1e-9){
            P1.real[iz]=0;
            P2.real[iz]=0;
        }else{
            P1.real[iz]=P1.real[iz]/creal(rho_fin.real[iz]);
            P2.real[iz]=P2.real[iz]/creal(rho_fin.real[iz]);
        }
        rho_fin.real[iz] *= 0.02*V_pear;
        //rho_fin.real[iz] *= 0.01*M_PI*V_pear;
    }


    rho_fin.print(n_bins, DeltaR);
    P1.print(n_bins, DeltaR);
    P2.print(n_bins, DeltaR);
}
