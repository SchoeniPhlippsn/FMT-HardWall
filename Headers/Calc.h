
void CalcRho(){
    
    for(j=0; j<= 100; j++){
        
        rho_sum_o[j] = 0;

        theta = DeltaT*j;

        b2 = WallDistance(theta);
        
        theta = DeltaT*(100-j);
        
        b3 = WallDistance(theta);
        
        for ( iz = 0; iz < n_bins; iz++){
            z = DeltaR*(iz+0.5);
            if(z <= H2){ 
                if( z <= b2 ) rho[j].real[iz] = 0.;
                else rho[j].real[iz] = bulk;
            }else{
                z = H - z;
                if( z <= b3 ) rho[j].real[iz] = 0.;
                else rho[j].real[iz] = bulk;
            }
            rho_sum_o[j] += creal(rho[j].real[iz]);
        }


    }
    std::cout << "Calculated:  Rho" << std::endl;
}


void CalcW0(){

    for( j = 0; j <= 100; j++){
        theta = DeltaT*j;
        costheta = DeltaK*cos(theta);
        sintheta = DeltaK*sin(theta);
        for ( iz = 0; iz < n_bins; iz++){
            if( iz < n_bins/2){
                params.kr = sintheta*iz;
                kz = costheta*iz;
            }else{
                params.kr = sintheta*(n_bins-iz);
                kz = costheta*(iz-n_bins);
            }
            
            w0[j].fourier[iz] = 0.5*(w0func(b1,params)*(cos(kz*b1)-I*sin(kz*b1))+w0func(-b1,params)*(cos(kz*b1)+I*sin(kz*b1))) ;
            z= -b1+DeltaP;
            while(z < b1){ 
                w0[j].fourier[iz] += w0func(z,params)*(cos(kz*z)-I*sin(kz*z));
                z += DeltaP;
            }
            w0[j].fourier[iz] *= inv_n*b1/50;
        }
        std::cout << "             W0[" << j << "]" << std::endl;
    }
}

void CalcW1(){
            
    for ( j = 0; j <= 100; j++){
        theta = DeltaT*j;
        costheta = DeltaK*cos(theta);
        sintheta = DeltaK*sin(theta);
        for ( i = 0; i <= lmax; i++){
            params.l = i;
            for ( iz = 0; iz < n_bins; iz++){
                if( iz < n_bins/2){
                    params.kr = sintheta*iz;
                    kz = costheta*iz;
                }else{
                    params.kr = sintheta*(n_bins-iz);
                    kz = costheta*(iz-n_bins);
                }

                w1[i][j].real[iz] = 0;
                for( v = 0; v<=2*i; v++){
                    params.m = v-i;

                    w1[i][j].real[iz] = 0.5*(w1func(b1,params)*(cos(kz*b1)-I*sin(kz*b1))+w1func(-b1,params)*(cos(kz*b1)+I*sin(kz*b1)));
                    z= -b1+DeltaP;
                    while(z < b1){ 
                        w1[i][j].real[iz] += w1func(z,params)*(cos(kz*z)-I*sin(kz*z));
                        z += DeltaP;
                    }
                    w1[i][j].real[iz] *= inv_n*b1/50;
                    
                    w1[i][j].fourier[iz] += Wignerd(-theta, i ,0, params.m)*w1[i][j].real[iz];
                } 
            }
            std::cout << "             W1_" << i << "[" << j << "]" << std::endl;
        }
    }
}    


void CalcW2(){
    
    for ( j = 0; j <= 100; j++){
        theta = DeltaT*j;
        costheta = DeltaK*cos(theta);
        sintheta = DeltaK*sin(theta);
        for ( i = 0; i <= lmax; i++){
            params.l = i;
            for ( iz = 0; iz < n_bins; iz++){
                if( iz < n_bins/2){
                    params.kr = sintheta*iz;
                    kz = costheta*iz;
                }else{
                    params.kr = sintheta*(n_bins-iz);
                    kz = costheta*(iz-n_bins);
                }

                w2[i][j].real[iz] = 0;
                for( v = 0; v<=2*i; v++){
                    params.m = v-i;

                    w2[i][j].real[iz] = 0.5*(w2func(b1,params)*(cos(kz*b1)-I*sin(kz*b1))+w2func(-b1,params)*(cos(kz*b1)+I*sin(kz*b1)));
                    z= -b1+DeltaP;
                    while(z < b1){ 
                        w2[i][j].real[iz] += w2func(z,params)*(cos(kz*z)-I*sin(kz*z));
                        z += DeltaP;
                    }
                    w2[i][j].real[iz] *= inv_n*b1/50;
                    
                    w2[i][j].fourier[iz] += Wignerd(-theta, i ,0, params.m)*w1[i][j].real[iz];
                } 
            }
            std::cout << "             W2_" << i << "[" << j << "]" << std::endl;
        }
    }
}

void CalcW3(){

    for( j = 0; j <= 100; j++){
        theta = DeltaT*j;
        costheta = DeltaK*cos(theta);
        sintheta = DeltaK*sin(theta);
        for ( iz = 0; iz < n_bins; iz++){
            if( iz < n_bins/2){
                params.kr = sintheta*iz;
                kz = costheta*iz;
            }else{
                params.kr = sintheta*(n_bins-iz);
                kz = costheta*(iz-n_bins);
            }
            
            w3[j].fourier[iz] = 0.5*(w3func(b1,params)*(cos(kz*b1)-I*sin(kz*b1))+w3func(-b1,params)*(cos(kz*b1)+I*sin(kz*b1))) ;
            z= -b1+DeltaP;
            while(z < b1){ 
                w3[j].fourier[iz] += w3func(z,params)*(cos(kz*z)-I*sin(kz*z));
                z += DeltaP;
            }
            w3[j].fourier[iz] *= inv_n*b1/50;
        }
        w3[j].invfft();
        std::cout << "             W3[" << j << "]" << std::endl;
    }
}



void CalcN(){

    for( iz = 0; iz < n_bins; iz++){
        n0.fourier[iz] = 0;
        n3.fourier[iz] = 0;
        for ( j = 1; j < 100; j++){
            n0.fourier[iz] += rho[j].fourier[iz]*w0[j].fourier[iz];
            n3.fourier[iz] += rho[j].fourier[iz]*w3[j].fourier[iz];
        }
        n0.fourier[iz] += 0.5*rho[0].fourier[iz]*w0[j].fourier[iz];
        n3.fourier[iz] += 0.5*rho[0].fourier[iz]*w3[j].fourier[iz];
        n0.fourier[iz] += 0.5*rho[100].fourier[iz]*w0[j].fourier[iz];
        n3.fourier[iz] += 0.5*rho[100].fourier[iz]*w3[j].fourier[iz];

        n0.fourier[iz] *= DeltaT;
        n3.fourier[iz] *= DeltaT;

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

            n1[i].fourier[iz] *= DeltaT;
            n2[i].fourier[iz] *= DeltaT;
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
        for ( j = 0; j <= 100; j++){
            sc0[j].fourier[iz] =  PHI0.fourier[iz]*conj(w0[j].fourier[iz]);
            for ( i = 0; i <= lmax; i++){
                sc1[i][j].fourier[iz] =  PHI1[i].fourier[iz]*conj(w1[i][j].fourier[iz]);
                sc2[i][j].fourier[iz] =  PHI2[i].fourier[iz]*conj(w2[i][j].fourier[iz]);
            }
            sc3[j].fourier[iz] =  PHI3.fourier[iz]*conj(w3[j].fourier[iz]);
        }
    }
}


void CalcC(){

    for( iz = 0; iz < n_bins; iz++){ 
        for ( j = 0; j <= 100; j++){
            c1[j].fourier[iz] = sc0[j].fourier[iz];
            for ( i = 0; i <= lmax; i++){
                c1[j].fourier[iz] += sc1[i][j].fourier[iz];
                c1[j].fourier[iz] += sc2[i][j].fourier[iz];
            }
            c1[j].fourier[iz] += sc3[j].fourier[iz];
        }
    }
    
    for ( j = 0; j <= 100; j++){
        c1[j].invfft();
    }
}

void CalcRhoN(){

    for(j=0; j<= 100; j++){
        
        rho_sum = 0;
        
        theta = DeltaT*j;

        b2 = WallDistance(theta);
        
        theta = DeltaT*(100-j);
        
        b3 = WallDistance(theta);

        for ( iz = 0; iz < n_bins; iz++){
            z = DeltaR*(iz+0.5);
            if(z <= H2){ 
                if( z <= b2 ) rhon[j].real[iz] = 0.;
                else rhon[j].real[iz] = exp(creal(-c1[j].real[iz]));// + rho_sum[(ix]));
            }else{
                z = H - z;
                if( z <= b3 ) rhon[j].real[iz] = 0.;
                else rhon[j].real[iz] = exp(creal(-c1[j].real[iz]));// + rho_sum[(ix]));
            }
            rho_sum += creal(rhon[j].real[iz]);
        }
    
        sclf = rho_sum_o[j]/rho_sum;
        if(j== 100 ){ 
            printf("%i %f %f\n",st, rho_sum, sclf);
            prev = now;
            now = rho_sum;
        }

        for( iz =0; iz < n_bins; iz++){
            rhon[j].real[iz] *= sclf;
            rho[j].real[iz] = rho[j].real[iz]*(1-alpha) + rhon[j].real[iz]*alpha;     
        }
    }
}
