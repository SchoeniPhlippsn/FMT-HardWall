
void CalcRho(){
    
    for(j=0; j<= 100; j++){
        
        theta = j*0.02-1;
        b3 = WallDistance(theta);
        
        theta = -theta;
        
        b2 = WallDistance(theta);
        
        for ( iz = 0; iz < n_bins; iz++){
            z = DeltaR*(iz+0.5);
            if(z <= H2){ 
                if( z <= b2 ) rho[j].real[iz] = 0.;
                else rho[j].real[iz] = 1;
            }else{
                z = H - z;
                if( z <= b3 ) rho[j].real[iz] = 0.;
                else rho[j].real[iz] = 1;
            }
        }
    }
    
    rho_sum = 0;
    for( iz = 0; iz < n_bins; iz++){

        rho_fin.real[iz] = 0.5*(rho[0].real[iz]+rho[100].real[iz]);

        for ( j = 1; j < 100; j++){
            rho_fin.real[iz] += rho[j].real[iz];
        }
        rho_fin.real[iz] *= 0.02;
        
        rho_sum += creal(rho_fin.real[iz]);
    }
    sclf = rho_sum_o/rho_sum;
    
    for( iz =0; iz < n_bins; iz++){ 
        rho_fin.real[iz] *= sclf; 
        for(j=0; j<= 100; j++) rho[j].real[iz] *= sclf;
    }
     
    std::cout << "Calculated:  Rho" << std::endl;
}


void CalcW(){

    for( j = 0; j <= 100; j++){
        theta = j*0.02-1;
        atheta = -acos(theta);
        costheta = DeltaK*theta;
        sintheta = DeltaK*sqrt(1-theta*theta);
        
        params.tt = 0.5;
        
        params.gauss = gaussC(0.5,a2);
        params.mean = meanC(0.5,a2);
        params.diff = diffC(0.5,a2);
        params.RR = RR(0.5,a2);
        params.AblRR = RRAbl2(b1,a2);
        params.theta = getTheta(0.5,a2);
        
        for ( iz = 0; iz < n_bins; iz++){
            if( iz < n_bins_2){
                params.kr = sintheta*iz;
                kz = costheta*iz;
            }else{
                params.kr = sintheta*(n_bins-iz);
                kz = costheta*(iz-n_bins);
            }
            
            
            exponent = cexp(-I*kz*b1); 

            w0[j].fourier[iz] = 0.5*w0func(params)*exponent;
            w3[j].fourier[iz] = 0.5*w3func(params)*exponent;
            
            for ( i = 0; i <= lmax; i++){
                params.l = i;
                w1[i][j].fourier[iz] = 0;
                w2[i][j].fourier[iz] = 0;
                for( v = 0; v<=2*i; v++){
                    params.m = v-i;
                    Wig = Wignerd(atheta, i ,0, params.m);

                    w1[i][j].real[iz] = 0.5*w1func(params)*exponent;
                    w1[i][j].fourier[iz] += Wig*w1[i][j].real[iz];
                    w2[i][j].real[iz] = 0.5*w1func(params)*exponent;
                    w2[i][j].fourier[iz] += Wig*w1[i][j].real[iz];
                }
            }
        }
        
        params.gauss = gaussC(0.5,a3);
        params.mean = meanC(0.5,a3);
        params.diff = diffC(0.5,a3);
        params.RR = RR(0.5,a3);
        params.AblRR = RRAbl2(b1,a3);
        params.theta = getTheta(0.5,a3);
        for ( iz = 0; iz < n_bins; iz++){
            if( iz < n_bins_2){
                params.kr = sintheta*iz;
                kz = costheta*iz;
            }else{
                params.kr = sintheta*(n_bins-iz);
                kz = costheta*(iz-n_bins);
            }
            
            exponent = cexp(I*kz*b1); 
            
            w0[j].fourier[iz] += 0.5*w0func(params)*exponent;
            w3[j].fourier[iz] += 0.5*w3func(params)*exponent ;
            
            for ( i = 0; i <= lmax; i++){
                params.l = i;
                for( v = 0; v<=2*i; v++){
                    params.m = v-i;
                    Wig = Wignerd(atheta, i ,0, params.m);

                    w1[i][j].real[iz] = 0.5*w1func(params)*exponent;
                    w1[i][j].fourier[iz] += Wig*w1[i][j].real[iz];
                    w2[i][j].real[iz] = 0.5*w1func(params)*exponent;
                    w2[i][j].fourier[iz] += Wig*w1[i][j].real[iz];
                }
            }
        } 
            
        z= -b1+DeltaP;

        while(z < b1){
            if( z < 0){
                z2 = -z;
                params.aa = a3;
            }else{
                z2 = z;
                params.aa = a2;
            }
            params.tt = 0.5*(1-sqrt(1-z2*b1_inv));
            params.gauss = gaussC(params.tt,params.aa);
            params.mean = meanC(params.tt,params.aa);
            params.diff = diffC(params.tt,params.aa);
            params.RR = RR(params.tt,params.aa);
            params.AblRR = RRAbl2(z2,params.aa);
            params.theta = getTheta(params.tt,params.aa);
            
            for ( iz = 0; iz < n_bins; iz++){
                if( iz < n_bins_2){
                    params.kr = sintheta*iz;
                    kz = costheta*iz;
                }else{
                    params.kr = sintheta*(n_bins-iz);
                    kz = costheta*(iz-n_bins);
                }

               // double bessel[2] = { gsl_sf_bessel_J0(params.kr*params.RR), gsl_sf_bessel_J1(params.kr*params.RR)};
                exponent = cexp(-I*kz*z);

                w0[j].fourier[iz] += w0func(params)*exponent;
                w3[j].fourier[iz] += w3func(params)*exponent;
                for ( i = 0; i <= lmax; i++){
                    params.l = i;
                    for( v = 0; v<=2*i; v++){
                        params.m = v-i;
                        Wig = Wignerd(atheta, i ,0, params.m);
                        w1[i][j].real[iz] = w1func(params)*exponent;
                        w1[i][j].fourier[iz] += Wig*w1[i][j].real[iz];
                        w2[i][j].real[iz] = w2func(params)*exponent;
                        w2[i][j].fourier[iz] += Wig*w2[i][j].real[iz];
                    }
                }
            }
            z += DeltaP;
        }
        for ( iz = 0; iz < n_bins; iz++){
            w0[j].fourier[iz] *= DeltaP*inv_n;
            w3[j].fourier[iz] *= DeltaP*inv_n;
            for ( i = 0; i <= lmax; i++){
                w1[i][j].fourier[iz] *= DeltaP*inv_n;
                w2[i][j].fourier[iz] *= DeltaP*inv_n;
            }
        }
        std::cout << "             W[" << j << "]" << std::endl;
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

        n0.fourier[iz] *= 0.02;
        n3.fourier[iz] *= 0.02;

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

            n1[i].fourier[iz] *= 0.02;
            n2[i].fourier[iz] *= 0.02;
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

    prev = now;
    now = 0;
    
    for(j=0; j<= 100; j++){
    
        theta = j*0.02-1;

        b3 = WallDistance(theta);
        
        theta = -theta;
        
        b2 = WallDistance(theta);

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
        }
    } 
    
    
    rho_sum =0;
    for( iz = 0; iz < n_bins; iz++){
        rho_fin.real[iz] = 0.5*(rhon[0].real[iz]+rhon[100].real[iz]);

        for ( j = 1; j < 100; j++){ 
            rho_fin.real[iz] += rhon[j].real[iz];
        }
        rho_fin.real[iz] *= 0.02;
        
        rho_sum += creal(rho_fin.real[iz]);
    }
    now = rho_sum; 
    sclf = rho_sum_o/rho_sum;
    printf("%i %f %f\n", st, now, sclf);
    
    for( iz =0; iz < n_bins; iz++){ 
        rho_fin.real[iz] *= sclf; 
        for(j=0; j<= 100; j++){ 
            rhon[j].real[iz] *= sclf;
            rho[j].real[iz] = rho[j].real[iz]*(1-alpha) + rhon[j].real[iz]*alpha;     
        }
    }
}
