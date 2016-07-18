
void CalcRho(){

	for ( ix = 0; ix < n_bins; ix++){
		for ( iy = 0; iy < n_bins; iy++){
            for ( iz = 0; iz < n_bins; iz++){
                z = DeltaR*(iz+0.5);
                if(z > H2) z = H - z;
                if( z <= RT ) rho.real[(ix*n_bins+iy)*n_bins+iz] = 0;
                else rho.real[(ix*n_bins+iy)*n_bins+iz] = bulk;
            }
		}
	}
    std::cout << "Calculated:  Rho" << std::endl;
}


void CalcW0(){

	for ( ix = 0; ix < n_bins; ix++){
        if(ix < n_bins_2) kx = DeltaK*(ix+0.5);
		else kx = DeltaK*(ix-n_bins+0.5);

		for ( iy = 0; iy < n_bins; iy++){
            if(iy < n_bins_2) ky = DeltaK*(iy+0.5);
            else ky = DeltaK*(iy-n_bins+0.5);
           
            params.kr = sqrt(kx*kx+ky*ky); 

            for ( iz = 0; iz < n_bins; iz++){
                
                if(iz < n_bins_2) kz = DeltaR*(iz+0.5);
                else kz = DeltaR*(iz-n_bins+0.5);

                if(fabs(kz) < b1) wreal[iz]= w0func(kz,params);
                else wreal[iz]= 0;
            }

            fftw_execute(wr2f); 
         
            for ( iz = 0; iz < n_bins; iz++){
                w0.fourier[(ix*n_bins+iy)*n_bins+iz] = inv_nDeltaR*wfourier[iz];
            }
        }
    }

    std::cout << "Calculated:  W0" << std::endl;
}

void CalcW1(){
            
    for ( ix = 0; ix < n_bins; ix++){
        for ( iy = 0; iy < n_bins; iy++){
            for ( iz = 0; iz < n_bins; iz++){
                for ( i = 0; i <= lmax; i++){
                    w1[i].fourier[(ix*n_bins+iy)*n_bins+iz] = 0;
                }
            }
        }
    }
    
    for ( i = 0; i <= lmax; i++){
        params.l = i;
        for( j = 0; j<=2*i; j++){
            params.m = j-i;
            for ( ix = 0; ix < n_bins; ix++){
                if(ix < n_bins_2) kx = DeltaK*(ix+0.5);
                else kx = DeltaK*(ix-n_bins+0.5);

                for ( iy = 0; iy < n_bins; iy++){
                    if(iy < n_bins_2) ky = DeltaK*(iy+0.5);
                    else ky = DeltaK*(iy-n_bins+0.5);
                   
                    params.kr = sqrt(kx*kx+ky*ky); 

                    for ( iz = 0; iz < n_bins; iz++){
                        
                        if(iz < n_bins_2) kz = DeltaR*(iz+0.5);
                        else kz = DeltaR*(iz-n_bins+0.5);

                        if(fabs(kz) < b1) wreal[iz]= w1func(kz,params);
                        else wreal[iz]= 0;
                    }

                    fftw_execute(wr2f); 
                 
                    for ( iz = 0; iz < n_bins; iz++){
                        w1[i].fourier[(ix*n_bins+iy)*n_bins+iz] += WignerdInt(l,0,params.m)*inv_nDeltaR*wfourier[iz];
                    }
                }
            }
            std::cout << "             W1^" << params.l << "_" << params.m << std::endl;
        }
    }
}

void CalcW2(){

    for ( ix = 0; ix < n_bins; ix++){
        for ( iy = 0; iy < n_bins; iy++){
            for ( iz = 0; iz < n_bins; iz++){
                for ( i = 0; i <= lmax; i++){
                    w2[i].fourier[(ix*n_bins+iy)*n_bins+iz] = 0;
                }
            }
        }
    }

    for ( i = 0; i <= lmax; i++){
        params.l = i;
        for( j = 0; j<=2*i; j++){
            params.m = j-i;
            for ( ix = 0; ix < n_bins; ix++){
                if(ix < n_bins_2) kx = DeltaK*(ix+0.5);
                else kx = DeltaK*(ix-n_bins+0.5);

                for ( iy = 0; iy < n_bins; iy++){
                    if(iy < n_bins_2) ky = DeltaK*(iy+0.5);
                    else ky = DeltaK*(iy-n_bins+0.5);
                   
                    params.kr = sqrt(kx*kx+ky*ky); 

                    for ( iz = 0; iz < n_bins; iz++){
                        
                        if(iz < n_bins_2) kz = DeltaR*(iz+0.5);
                        else kz = DeltaR*(iz-n_bins+0.5);

                        if(fabs(kz) < b1) wreal[iz]= w2func(kz,params);
                        else wreal[iz]= 0;
                        
                    }

                    fftw_execute(wr2f); 
                 
                    for ( iz = 0; iz < n_bins; iz++){
                        w2[i].fourier[(ix*n_bins+iy)*n_bins+iz] += WignerdInt(l,0,params.m)*inv_nDeltaR*wfourier[iz];
                    }
                }
            }
            std::cout << "             W2^" << params.l << "_" << params.m << std::endl;
        }
    }
}

void CalcW3(){

   // FILE * o1File;
   // o1File = fopen ("Realfftw3.dat","w");

	for ( ix = 0; ix < n_bins; ix++){
        if(ix < n_bins_2) kx = DeltaK*(ix+0.5);
		else kx = DeltaK*(ix-n_bins+0.5);

		for ( iy = 0; iy < n_bins; iy++){
            if(iy < n_bins_2) ky = DeltaK*(iy+0.5);
            else ky = DeltaK*(iy-n_bins+0.5);
           
            params.kr = sqrt(kx*kx+ky*ky); 

            for ( iz = 0; iz < n_bins; iz++){
                
                if(iz < n_bins_2) kz = DeltaR*(iz+0.5);
                else kz = DeltaR*(iz-n_bins+0.5);

                if(fabs(kz) < b1) wreal[iz]= w3func(kz,params);
                else wreal[iz]= 0;
            }

            fftw_execute(wr2f); 
         
            for ( iz = 0; iz < n_bins; iz++){
                w3.fourier[(ix*n_bins+iy)*n_bins+iz] = inv_nDeltaR*wfourier[iz];
       //         fprintf (o1File, "%f %i %.10f\n", sqrt(ix*ix+iy*iy), iz, wfourier[iz]);
                //std::cout << DeltaR*(ix+0.5) << " " << DeltaR*(iy+0.5) << " " << DeltaR*(iz+0.5) << " " << creal(w3.fourier[(ix*n_bins+iy)*n_bins+iz]) << " " << cimag(w3.fourier[(ix*n_bins+iy)*n_bins+iz]) << std::endl;
            }
        }
    }
    
    std::cout << "             W3" << std::endl;
    
    gsl_integration_workspace_free (w);
    gsl_integration_qawo_table_free (table);
}

void CalcN(){

    for( ix = 0; ix < n_bins; ix++){
        for( iy = 0; iy < n_bins; iy++){
            for( iz = 0; iz < n_bins; iz++){
                n0.fourier[(ix*n_bins+iy)*n_bins+iz] = rho.fourier[(ix*n_bins+iy)*n_bins+iz]*w0.fourier[(ix*n_bins+iy)*n_bins+iz];
                for ( i = 0; i <= lmax; i++){
                    n1[i].fourier[(ix*n_bins+iy)*n_bins+iz] = rho.fourier[(ix*n_bins+iy)*n_bins+iz]*w1[i].fourier[(ix*n_bins+iy)*n_bins+iz];
                    n2[i].fourier[(ix*n_bins+iy)*n_bins+iz] = rho.fourier[(ix*n_bins+iy)*n_bins+iz]*w2[i].fourier[(ix*n_bins+iy)*n_bins+iz];
                }
                n3.fourier[(ix*n_bins+iy)*n_bins+iz] = rho.fourier[(ix*n_bins+iy)*n_bins+iz]*w3.fourier[(ix*n_bins+iy)*n_bins+iz];
            }
        }
    }
    
    n0.invfft();
    for ( i = 0; i <= lmax; i++){
        n1[i].invfft();
        n2[i].invfft();
    }
    n3.invfft();
}


void CalcPHI(){

    for( ix = 0; ix < n_bins; ix++){
        for( iy = 0; iy < n_bins; iy++){
            for( iz = 0; iz < n_bins; iz++){
                on3 = 1-creal(n3.real[(ix*n_bins+iy)*n_bins+iz]);
                inv_on3 = 1/on3;
                
                Phi0();
                Phi3();
               
                for ( i = 0; i <= lmax; i++){
                    Phi1();
                    Phi2();
                }
            }
        }
    }

    PHI0.fft();
    for ( i = 0; i <= lmax; i++){
        PHI1[i].fft();
        PHI2[i].fft();
    }
    PHI3.fft();
}

void CalcSC(){

    for( ix = 0; ix < n_bins; ix++){
        for( iy = 0; iy < n_bins; iy++){
            for( iz = 0; iz < n_bins; iz++){
                sc0.fourier[(ix*n_bins+iy)*n_bins+iz] =  PHI0.fourier[(ix*n_bins+iy)*n_bins+iz]*w0.fourier[(ix*n_bins+iy)*n_bins+iz];
                for ( i = 0; i <= lmax; i++){
                    sc1[i].fourier[(ix*n_bins+iy)*n_bins+iz] =  -PHI1[i].fourier[(ix*n_bins+iy)*n_bins+iz]*w1[i].fourier[(ix*n_bins+iy)*n_bins+iz];
                    sc2[i].fourier[(ix*n_bins+iy)*n_bins+iz] =  -PHI2[i].fourier[(ix*n_bins+iy)*n_bins+iz]*w2[i].fourier[(ix*n_bins+iy)*n_bins+iz];
                }
                sc3.fourier[(ix*n_bins+iy)*n_bins+iz] =  PHI3.fourier[(ix*n_bins+iy)*n_bins+iz]*w3.fourier[(ix*n_bins+iy)*n_bins+iz];
            }
        }
    }

    sc0.invfft();
    for ( i = 0; i <= lmax; i++){
        sc1[i].invfft();
        sc2[i].invfft();
    }
    sc3.invfft();
}

void CalcC(){

    for( ix = 0; ix < n_bins; ix++){ 
        for( iy = 0; iy < n_bins; iy++){ 
            for( iz = 0; iz < n_bins; iz++){ 
                c1.real[(ix*n_bins+iy)*n_bins+iz] = 0;
                c1.real[(ix*n_bins+iy)*n_bins+iz] -= sc0.real[(ix*n_bins+iy)*n_bins+iz];
                for ( i = 0; i <= lmax; i++){
                    c1.real[(ix*n_bins+iy)*n_bins+iz] -= sc1[i].real[(ix*n_bins+iy)*n_bins+iz];
                    c1.real[(ix*n_bins+iy)*n_bins+iz] -= sc2[i].real[(ix*n_bins+iy)*n_bins+iz];
                }
                c1.real[(ix*n_bins+iy)*n_bins+iz] -= sc3.real[(ix*n_bins+iy)*n_bins+iz];
            }
        }
    }

}
