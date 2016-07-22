
void CalcRho(){
    
    for(j=0; j<= 100; j++){
        for ( iz = 0; iz < n_bins; iz++){
            z = DeltaR*(iz+0.5);
            if(z > H2) z = H - z;
            if( z <= RII ) rho[j].real[iz] = 0;
            else rho[j].real[iz] = bulk;
        }
    }
    std::cout << "Calculated:  Rho" << std::endl;
}


void CalcW0(){

    for ( iz = 0; iz < n_bins; iz++){
                
        if(iz < n_bins_2) kz = DeltaR*(iz+0.5);
        else kz = DeltaR*(iz-n_bins+0.5);

        if(fabs(kz) < b1) w0.real[iz]= w0func(kz,params);
        else w0.real[iz]= 0;
    }

    w0.fft();
 
    for ( iz = 0; iz < n_bins; iz++){
        w0.fourier[iz] *= inv_nDeltaR;
    }

    std::cout << "             W0" << std::endl;
}


void CalcW1(){
            
    for ( j = 0; j <= 100; j++){
        theta = DeltaK*j;
        for ( i = 0; i <= lmax; i++){
            params.l = i;
            for ( iz = 0; iz < n_bins; iz++){
                
                if(iz < n_bins_2) kz = DeltaR*(iz+0.5);
                else kz = DeltaR*(iz-n_bins+0.5);
                w1[i][j].real[iz] = 0;
                if(fabs(kz) < b1){
                    for( v = 0; v<=2*i; v++){
                        params.m = v-i;
                        
                        w1[i][j].real[iz] += Wignerd(theta, i ,0, params.m)*w1func(kz,params);
                        
                    }
                }
            }

            w1[i][j].fft(); 
         
            for ( iz = 0; iz < n_bins; iz++){
                w1[i][j].fourier[iz] *= inv_nDeltaR;
            }
        }
    }
    std::cout << "             W1" << std::endl;
}    


void CalcW2(){
    
    for ( j = 0; j <= 100; j++){
        theta = DeltaK*j;
        for ( i = 0; i <= lmax; i++){
            params.l = i;
            for ( iz = 0; iz < n_bins; iz++){
                
                if(iz < n_bins_2) kz = DeltaR*(iz+0.5);
                else kz = DeltaR*(iz-n_bins+0.5);
                w2[i][j].real[iz] = 0;
                if(fabs(kz) < b1){
                    for( v = 0; v<=2*i; v++){
                        params.m = v-i;
                        
                        w2[i][j].real[iz] += Wignerd(theta, params.l ,0, params.m)*w2func(kz,params);
                    }
                }
            }

            w2[i][j].fft(); 
         
            for ( iz = 0; iz < n_bins; iz++){
                w2[i][j].fourier[iz] *= inv_nDeltaR;
            }
        }
    }
    std::cout << "             W2" << std::endl;
}

void CalcW3(){

    for ( iz = 0; iz < n_bins; iz++){
        
        if(iz < n_bins_2) kz = DeltaR*(iz+0.5);
        else kz = DeltaR*(iz-n_bins+0.5);

        if(fabs(kz) < b1) w3.real[iz]= w3func(kz,params);
        else w3.real[iz]= 0;
    }

    w3.fft(); 
 
    for ( iz = 0; iz < n_bins; iz++){
        w3.fourier[iz] *= inv_nDeltaR;
    }
    std::cout << "             W3" << std::endl;
}


void CalcN(){

    for( iz = 0; iz < n_bins; iz++){
        n0.fourier[iz] = 0;
        n3.fourier[iz] = 0;
        for ( j = 1; j < 100; j++){
            n0.fourier[iz] += rho[j].fourier[iz]*w0.fourier[iz];
            n3.fourier[iz] += rho[j].fourier[iz]*w3.fourier[iz];
        }
        n0.fourier[iz] += 0.5*rho[0].fourier[iz]*w0.fourier[iz];
        n3.fourier[iz] += 0.5*rho[0].fourier[iz]*w3.fourier[iz];
        n0.fourier[iz] += 0.5*rho[100].fourier[iz]*w0.fourier[iz];
        n3.fourier[iz] += 0.5*rho[100].fourier[iz]*w3.fourier[iz];

        n0.fourier[iz] *= DeltaK;
        n3.fourier[iz] *= DeltaK;
        
        for ( i = 0; i <= lmax; i++){
            n1[i].fourier[iz] = 0;
            n2[i].fourier[iz] = 0;
            for ( j = 1; j < 100; j++){
                n1[i].fourier[iz] += rho[j].fourier[iz]*w1[i][j].fourier[iz];
                n2[i].fourier[iz] += rho[j].fourier[iz]*w2[i][j].fourier[iz];
            }
            n1[i].fourier[iz] += 0.5*rho[0].fourier[iz]*w1[i][0].fourier[iz];
            n2[i].fourier[iz] += 0.5*rho[0].fourier[iz]*w2[i][0].fourier[iz];
            n1[i].fourier[iz] += 0.5*rho[100].fourier[iz]*w1[i][100].fourier[iz];
            n2[i].fourier[iz] += 0.5*rho[100].fourier[iz]*w2[i][100].fourier[iz];

            n1[i].fourier[iz] *= DeltaK;
            n2[i].fourier[iz] *= DeltaK;
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

    for( iz = 0; iz < n_bins; iz++){
        on3 = 1-creal(n3.real[iz]);
        inv_on3 = 1/on3;
        
        Phi0();
        Phi3();
       
        for ( i = 0; i <= lmax; i++){
            Phi1();
            Phi2();
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

    for( iz = 0; iz < n_bins; iz++){
        sc0.fourier[iz] =  PHI0.fourier[iz]*w0.fourier[iz];
        for ( i = 0; i <= lmax; i++){
            for ( j = 0; j <= 100; j++){
                sc1[i][j].fourier[iz] =  PHI1[i].fourier[iz]*w1[i][j].fourier[iz];
                sc2[i][j].fourier[iz] =  PHI2[i].fourier[iz]*w2[i][j].fourier[iz];
            }
        }
        sc3.fourier[iz] =  PHI3.fourier[iz]*w3.fourier[iz];
    }

    sc0.invfft();
    for ( i = 0; i <= lmax; i++){
        for ( j = 0; j <= 100; j++){
            sc1[i][j].invfft();
            sc2[i][j].invfft();
        }
    }
    sc3.invfft();
}


void CalcC(){

    for( iz = 0; iz < n_bins; iz++){ 
        for ( j = 0; j <= 100; j++){
            c1[j].real[iz] = 0;
            c1[j].real[iz] -= sc0.real[iz];
            for ( i = 0; i <= lmax; i++){
                c1[j].real[iz] -= sc1[i][j].real[iz];
                c1[j].real[iz] -= sc2[i][j].real[iz];
            }
            c1[j].real[iz] -= sc3.real[iz];
        }
    }
}
