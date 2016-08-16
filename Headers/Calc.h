
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
        atheta = acos(theta);
        costheta = DeltaK*theta;
        sintheta = DeltaK*sqrt(1-theta*theta);

        for ( iz = 0; iz < n_bins; iz++){
            w0[j].fourier[iz] = 0;
            w3[j].fourier[iz] = 0; 
            for ( i = 0; i <= lmax; i++){
                w1[i][j].fourier[iz] = 0;
                w2[i][j].fourier[iz] = 0;
            }
        }

        for( st = 0; st < 50; st ++){
            
            z = b1*Int_step100[st];
            params.tt = 0.5*(1-sqrt(1-z*b1_inv));
            params.gauss = gaussC(params.tt,a2);
            params.mean = meanC(params.tt,a2);
            params.diff = diffC(params.tt,a2);
            params.RR = RR(params.tt,a2);
            params.AblRR = RRAbl2(z,a2);
            params.theta = getTheta(params.tt,a2);
            
            PARAMS.tt = params.tt;

            PARAMS.gauss = gaussC(PARAMS.tt,a3);
            PARAMS.mean = meanC(PARAMS.tt,a3);
            PARAMS.diff = diffC(PARAMS.tt,a3);
            PARAMS.RR = RR(PARAMS.tt,a3);
            PARAMS.AblRR = RRAbl2(z,a3);
            PARAMS.theta = getTheta(PARAMS.tt,a3);
            
            for ( iz = 0; iz < n_bins; iz++){
                if( iz < n_bins_2){
                    params.kr = sintheta*iz;
                    PARAMS.kr = sintheta*iz;
                    kz = costheta*iz;
                }else{
                    params.kr = sintheta*(n_bins-iz);
                    PARAMS.kr = sintheta*(n_bins-iz);
                    kz = costheta*(iz-n_bins);
                }

                exponent = cexp(-I*kz*z);
            
                bessel[lmax] = gsl_sf_bessel_J0(params.kr*params.RR);  
                bessel[lmax+1] = gsl_sf_bessel_J1(params.kr*params.RR);  
                bessel[lmax-1] = -bessel[lmax+1];  

                for( i=2; i <= lmax; i++){
                    bessel[lmax+i] = gsl_sf_bessel_Jn(i,params.kr*params.RR);
                    if( i % 2 == 0) bessel[lmax-i] = bessel[lmax+i];
                    else bessel[lmax-i] = -bessel[lmax+i];
                }

                w0[j].fourier[iz] += Int_weight100[st]*w0func(bessel,params)*exponent;
                w3[j].fourier[iz] += Int_weight100[st]*w3func(bessel,params)*exponent;
                for ( i = 0; i <= lmax; i++){
                    params.l = i;
                    for( v = 0; v<=2*i; v++){
                        params.m = v-i;
                        Wig = Wignerd(atheta, i ,0, params.m);
                        w1[i][j].real[iz] = Int_weight100[st]*w1func(bessel,params)*exponent;
                        w1[i][j].fourier[iz] += Wig*w1[i][j].real[iz];
                        w2[i][j].real[iz] = Int_weight100[st]*w2func(bessel,params)*exponent;
                        w2[i][j].fourier[iz] += Wig*w2[i][j].real[iz];
                    }
                }
                exponent = cexp(I*kz*z);
            
                bessel[lmax] = gsl_sf_bessel_J0(PARAMS.kr*PARAMS.RR);  
                bessel[lmax+1] = gsl_sf_bessel_J1(PARAMS.kr*PARAMS.RR);  
                bessel[lmax-1] = -bessel[lmax+1];  

                for( i=2; i <= lmax; i++){
                    bessel[lmax+i] = gsl_sf_bessel_Jn(i,PARAMS.kr*PARAMS.RR);
                    if( i % 2 == 0) bessel[lmax-i] = bessel[lmax+i];
                    else bessel[lmax-i] = -bessel[lmax+i];
                }

                w0[j].fourier[iz] += Int_weight100[st]*w0func(bessel,PARAMS)*exponent;
                w3[j].fourier[iz] += Int_weight100[st]*w3func(bessel,PARAMS)*exponent;
                for ( i = 0; i <= lmax; i++){
                    PARAMS.l = i;
                    for( v = 0; v<=2*i; v++){
                        PARAMS.m = v-i;
                        Wig = Wignerd(atheta, i ,0, PARAMS.m);
                        w1[i][j].real[iz] = Int_weight100[st]*w1func(bessel,PARAMS)*exponent;
                        w1[i][j].fourier[iz] += Wig*w1[i][j].real[iz];
                        w2[i][j].real[iz] = Int_weight100[st]*w2func(bessel,PARAMS)*exponent;
                        w2[i][j].fourier[iz] += Wig*w2[i][j].real[iz];
                    }
                }
            }
        }
        for ( iz = 0; iz < n_bins; iz++){
            w0[j].fourier[iz] *= b1*inv_n;
            w3[j].fourier[iz] *= b1*inv_n;
            for ( i = 0; i <= lmax; i++){
                w1[i][j].fourier[iz] *= b1*inv_n;
                w2[i][j].fourier[iz] *= b1*inv_n;
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
    sclf = rho_sum_o/rho_sum;
    now = sclf; 
    printf("%i %f %f\n", st, now, sclf);
    
    for( iz =0; iz < n_bins; iz++){ 
        rho_fin.real[iz] *= sclf; 
        for(j=0; j<= 100; j++){ 
            rhon[j].real[iz] *= sclf;
            rho[j].real[iz] = rho[j].real[iz]*(1-alpha) + rhon[j].real[iz]*alpha;     
        }
    }
}
