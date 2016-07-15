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

	n_bins = 64;
	n_bins_2 = n_bins/2;
	inv_n = 1./n_bins;

	H = 10;
	H2 = H/2;
	
	RII = 0.5*sqrt(0.5);
	RT = 1.5*sqrt(0.5);

	DeltaR = H*inv_n;
	DeltaK = inv_n*2*M_PI/RII;

	bulk = 0.3*H*H*H/((H-2*RII)*(H-2*RII)*(H-2*RT)*4*M_PI*RII*RII*RII);

	rho_sum_o = (H-2*RII)*(H-2*RII)*(H-2*RT)*bulk*n_bins*n_bins*n_bins/(H*H*H);
    
    alpha = 0.1;

    Setup(); 
    
    CalcRho();
	
    CalcW0();
    
    CalcW1();
    CalcW2();
    CalcW3();
    
    for( st = 0 ; st <= 1000; st++){

        rho.fft();
       // if(st==0) rho.print(n_bins,DeltaR);

        CalcN();
    
   /*     if(st==0){		
            n0.print(n_bins,DeltaR); 
            n12.print(n_bins,DeltaR); 
            n222.print(n_bins,DeltaR); 
            n3.print(n_bins,DeltaR); 
        }
    */
        CalcPHI();

    
      /*  if(st==0){		
            PHI0.print(n_bins,DeltaR); 
            PHI12.print(n_bins,DeltaR); 
            PHI222.print(n_bins,DeltaR); 
            PHI3.print(n_bins,DeltaR); 
        }*/

        CalcSC();


        /*if(st==0){		
            sc0.print(n_bins,DeltaR); 
            sc12.print(n_bins,DeltaR); 
            sc222.print(n_bins,DeltaR); 
            sc3.print(n_bins,DeltaR); 
        }*/

        CalcC();

       // if(st==0) c1.print(n_bins,DeltaR); 
		
        rho_sum = 0;

        for ( ix = 0; ix < n_bins; ix++){
            x = DeltaR*(ix+0.5);
            if(x > H2) x = H-x;
            for ( iy = 0; iy < n_bins; iy++){
                y = DeltaR*(iy+0.5);
                if(y > H2) y = H-y;
                for ( iz = 0; iz < n_bins; iz++){
                    z = DeltaR*(iz+0.5);
                    if(z > H2) z = H - z;
                    if( x <= RII || y <= RII || z <= RT ) rho_n.real[(ix*n_bins+iy)*n_bins+iz] = 0.;
                    else rho_n.real[(ix*n_bins+iy)*n_bins+iz] = exp(creal(c1.real[(ix*n_bins+iy)*n_bins+iz]));// + rho_sum[(ix]));
                    rho_sum += creal(rho_n.real[(ix*n_bins+iy)*n_bins+iz]);
                }
            }
		}
		
        sclf = rho_sum_o/rho_sum;

        if(st==0) rho_n.print(n_bins,DeltaR); 
		
        if(st%1==0)printf("%i %f %f\n",st,rho_sum, sclf);
		for( ix =0; ix < n_bins; ix++){
			for( iy =0; iy < n_bins; iy++){
                for( iz =0; iz < n_bins; iz++){
                    rho_n.real[(ix*n_bins+iy)*n_bins+iz] *= sclf;
                    rho.real[(ix*n_bins+iy)*n_bins+iz] = rho.real[(ix*n_bins+iy)*n_bins+iz]*(1-alpha) + rho_n.real[(ix*n_bins+iy)*n_bins+iz]*alpha;     
                }
            }
		}
    }    
   
	FILE * oFile;
	FILE * oFile1;

	oFile = fopen ("rho_fin.dat","w");
	oFile1 = fopen ("rho_fin0.dat","w");

	for( ix =0; ix < n_bins; ix++){
		for( iy =0; iy < n_bins; iy++){
            fprintf (oFile1, "%f %f %f %f\n", DeltaR*(ix+0.5), DeltaR*(iy+0.5), creal(rho.real[(ix*n_bins+iy)*n_bins+n_bins_2]), cimag(rho.real[(ix*n_bins+iy)*n_bins+n_bins_2]));
            for( iz =0; iz < n_bins; iz++){
                fprintf (oFile, "%f %f %f %f %f\n", DeltaR*(ix+0.5), DeltaR*(iy+0.5), DeltaR*(iz+0.5), creal(rho.real[(ix*n_bins+iy)*n_bins+iz]), cimag(rho.real[(ix*n_bins+iy)*n_bins+iz]));
            }
        }
	}
}
