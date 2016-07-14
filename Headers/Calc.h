
void CalcRho(){

	for ( ix = 0; ix < n_bins; ix++){
		x = DeltaR*(ix+0.5);
		if(x > H2) x = H -x;
		for ( iy = 0; iy < n_bins; iy++){
			y = DeltaR*(iy+0.5);
            if(y > H2) y = H -y;
            for ( iz = 0; iz < n_bins; iz++){
                z = DeltaR*(iz+0.5);
                if(z > H2) z = H - z;
                if( x <= RII || y <= RII || z <= RT ) rho.real[(ix*n_bins+iy)*n_bins+iz] = 0;
                else rho.real[(ix*n_bins+iy)*n_bins+iz] = bulk;
            }
		}
	}
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
                w0.fourier[(ix*n_bins+iy)*n_bins+iz] = DeltaR*DeltaR*DeltaR*inv_n*inv_n*inv_n*wfourier[iz];
            }
        }
    }

    std::cout << "Calculated:  W0" << std::endl;
}

void CalcW1(){

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
                        w1[i][j].fourier[(ix*n_bins+iy)*n_bins+iz] = DeltaR*DeltaR*DeltaR*inv_n*inv_n*inv_n*wfourier[iz];
                    }
                }
            }
            std::cout << "             W1^" << params.l << "_" << params.m << std::endl;
        }
    }
}

void CalcW2(){

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
                        w2[i][j].fourier[(ix*n_bins+iy)*n_bins+iz] = DeltaR*DeltaR*DeltaR*inv_n*inv_n*inv_n*wfourier[iz];
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
                w3.fourier[(ix*n_bins+iy)*n_bins+iz] = DeltaR*DeltaR*DeltaR*inv_n*inv_n*inv_n*wfourier[iz];
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
                    for( j = 0; j<=2*i; j++){
                        n1[i][j].fourier[(ix*n_bins+iy)*n_bins+iz] = rho.fourier[(ix*n_bins+iy)*n_bins+iz]*w1[i][j].fourier[(ix*n_bins+iy)*n_bins+iz];
                        n2[i][j].fourier[(ix*n_bins+iy)*n_bins+iz] = rho.fourier[(ix*n_bins+iy)*n_bins+iz]*w2[i][j].fourier[(ix*n_bins+iy)*n_bins+iz];
                    }
                }
                n3.fourier[(ix*n_bins+iy)*n_bins+iz] = rho.fourier[(ix*n_bins+iy)*n_bins+iz]*w3.fourier[(ix*n_bins+iy)*n_bins+iz];
            }
        }
    }
    
    n0.invfft();
    for ( i = 0; i <= lmax; i++){
        for( j = 0; j<=2*i; j++){
        n1[i][j].invfft();
        n2[i][j].invfft();
        }
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
                    for( j = -i; j<=i; j++){
                        Phi1();
                        Phi2();
                    }
                }
            }
        }
    }

    PHI0.fft();
    for ( i = 0; i <= lmax; i++){
        for( j = 0; j<=2*i; j++){
            PHI1[i][j].fft();
            PHI2[i][j].fft();
        }
    }
    PHI3.fft();
}

void CalcSC(){

    for( ix = 0; ix < n_bins; ix++){
        for( iy = 0; iy < n_bins; iy++){
            for( iz = 0; iz < n_bins; iz++){
                sc0.fourier[(ix*n_bins+iy)*n_bins+iz] =  PHI0.fourier[(ix*n_bins+iy)*n_bins+iz]*w0.fourier[(ix*n_bins+iy)*n_bins+iz];
                for ( i = 0; i <= lmax; i++){
                    for( j = 0; j<=2*i; j++){
                        sc1[i][j].fourier[(ix*n_bins+iy)*n_bins+iz] =  -PHI1[i][j].fourier[(ix*n_bins+iy)*n_bins+iz]*w1[i][j].fourier[(ix*n_bins+iy)*n_bins+iz];
                        sc2[i][j].fourier[(ix*n_bins+iy)*n_bins+iz] =  -PHI2[i][j].fourier[(ix*n_bins+iy)*n_bins+iz]*w2[i][j].fourier[(ix*n_bins+iy)*n_bins+iz];
                    }
                }
                sc3.fourier[(ix*n_bins+iy)*n_bins+iz] =  PHI3.fourier[(ix*n_bins+iy)*n_bins+iz]*w3.fourier[(ix*n_bins+iy)*n_bins+iz];
            }
        }
    }

    sc0.invfft();
    for ( i = 0; i <= lmax; i++){
        for( j = 0; j<=2*i; j++){
            sc1[i][j].invfft();
            sc2[i][j].invfft();
        }
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
                    for( j = 0; j<=2*i; j++){
                        c1.real[(ix*n_bins+iy)*n_bins+iz] -= sc1[i][j].real[(ix*n_bins+iy)*n_bins+iz];
                        c1.real[(ix*n_bins+iy)*n_bins+iz] -= sc2[i][j].real[(ix*n_bins+iy)*n_bins+iz];
                    }
                }
                c1.real[(ix*n_bins+iy)*n_bins+iz] -= sc3.real[(ix*n_bins+iy)*n_bins+iz];
            }
        }
    }

}
