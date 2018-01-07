
void CalcRho(){

/*    for(j=0; j<= 100; j++){

        std::string Name = "Results/rho_fin" + toString(j) + ".dat";
        char file[30];
        strcpy(file,Name.c_str());
	    std::ifstream oFile(file);
        for( int iz =0; iz < n_bins; iz++){
            double z, freal,fimag;
            oFile >> z >> freal >> fimag;
            rho[j].real[iz] = freal+I*fimag;
        }
        oFile.close();
    }
*/    
    for(j = 0; j <= 100; j++){
        //atheta = j*M_PI*0.01; 
        //theta = cos(atheta); 
        theta = j*0.02-1;
        b3 = WallDistance(theta); //Pear
        //b3 = RII; //Sphere
        theta = -theta;
        
        b2 = WallDistance(theta); //Pear
        //b2 =RII; //Sphere
        
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
        for(j = 0; j <= 100; j++){
            rho[j].real[iz] *= sclf;
            rho[j].fft();
        }
    }

    std::cout << "Calculated:  Rho" << std::endl;
}

void CalcW(){

    for( j = 0; j <= 100; j++){
        /*if(j==0) theta=-0.9999999;
        else{ 
            if(j==100) theta=0.9999999;
            else theta = j*0.02-1;
        }
        atheta = asin(theta);
        sintheta = DeltaK*theta;
        costheta = DeltaK*sqrt(1-theta*theta);
        */
        //atheta = j*M_PI*0.01; 
        //theta = cos(atheta); 
        theta = j*0.02-1;
        if( theta < -0.9999) theta=-0.9999;
        if( theta > 0.9999) theta=0.9999;
        atheta = acos(theta);
/*        for ( i = 0; i <= lmax; i++){
            int ii = i;
            for( v = 0; v<=2*i; v++){
                int vv = v-i;
                Wig = Wignerd(atheta, ii ,0, vv);
                while(fabs(Wig) < 1e-9){
                   // std::cout << Wig << " " << theta << std::endl;
                    if(j==100) theta-=0.0001;
                    else theta += 0.0001;
                    atheta = acos(theta);
                    Wig = Wignerd(atheta, ii ,0, vv);
                    //std::cout << Wig << " " << theta << std::endl;
                }
            }
        }*/
        costheta = DeltaK*theta;
        sintheta = DeltaK*sqrt(1-theta*theta);
        if(fabs(sintheta)<1e-4) sintheta = 1e-3;

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
            params.RR = RR(params.tt,a2);//Pear
            //params.RR = RR(z,a2);
            params.AblRR = RRAbl2(z,a2);
            params.theta = getTheta(params.tt,a2); //Pear
            //params.theta = getTheta(z,a2);
            
            PARAMS.tt = params.tt;

            PARAMS.gauss = gaussC(PARAMS.tt,a3);
            PARAMS.mean = meanC(PARAMS.tt,a3);
            PARAMS.diff = diffC(PARAMS.tt,a3);
            PARAMS.RR = RR(PARAMS.tt,a3); //Pear
            //PARAMS.RR = RR(z,a3); 
            PARAMS.AblRR = RRAbl2(z,a3);
            PARAMS.theta = getTheta(-PARAMS.tt,a3);//Pear
            //PARAMS.theta = getTheta(-z,a3);
            
            for ( iz = 0; iz < n_bins; iz++){
//                if(iz < n_bins_2){
                    params.kr = sintheta*(iz+0.5)*params.RR;
                    PARAMS.kr = sintheta*(iz+0.5)*PARAMS.RR;
                    kz = costheta*(iz+0.5);
/*                }else{
                    params.kr = sintheta*(iz+0.5-n_bins);
                    PARAMS.kr = sintheta*(iz+0.5-n_bins);
                    kz = costheta*(iz+0.5-n_bins);
                }// */
                exponent = cexp(-I*kz*z);
            
                bessel[lmax] = gsl_sf_bessel_J0(params.kr);  
                bessel[lmax+1] = -gsl_sf_bessel_J1(params.kr);  
                bessel[lmax-1] = -bessel[lmax+1];  

                for( i=2; i <= lmax; i++){
                     if( i % 2 == 0){
                        bessel[lmax+i] = gsl_sf_bessel_Jn(i,params.kr);
                        bessel[lmax-i] = bessel[lmax+i];
                    }else{ 
                        bessel[lmax+i] = -gsl_sf_bessel_Jn(i,params.kr);
                        bessel[lmax-i] = -bessel[lmax+i];
                    }
                }
                params.kr = -sintheta*(iz+0.5);
                

                w0[j].fourier[iz] += Int_weight100[st]*w0func(bessel,params)*exponent;
                w3[j].fourier[iz] += Int_weight100[st]*w3func(bessel,params)*exponent;
                for ( i = 0; i <= lmax; i++){
                    params.l = i;
                    for( v = 0; v<=2*i; v++){
                        params.m = v-i;
                        //Wig = 1/Wignerd(atheta, i ,0, params.m);
                        Wig = Wignerd(atheta, i ,0, params.m);
                        w1[i][j].fourier[iz] += Wig*Int_weight100[st]*w1func(bessel,params)*exponent;
                        w2[i][j].fourier[iz] += Wig*Int_weight100[st]*w2func(bessel,params)*exponent;
                      //  if(j==0 || j==2) if(iz==0 && st==0)
                      //  std::cout << i << " " << v-i << ": " << Wig << " " << creal(w1func(bessel,params)) << " " << std::endl;
                    }
                   // if(j==0 || j==2) if(iz==0 && st==0)
                   //     std::cout << i << ": " << creal(w2[i][j].fourier[iz]) << std::endl;
                }
               // if(j==2) exit(0);
                exponent = cexp(I*kz*z);
            
                bessel[lmax] = gsl_sf_bessel_J0(PARAMS.kr);  
                bessel[lmax+1] = -gsl_sf_bessel_J1(PARAMS.kr);  
                bessel[lmax-1] = -bessel[lmax+1];  

                for( i=2; i <= lmax; i++){
                     if( i % 2 == 0){
                        bessel[lmax+i] = gsl_sf_bessel_Jn(i,PARAMS.kr);
                        bessel[lmax-i] = bessel[lmax+i];
                    }else{ 
                        bessel[lmax+i] = -gsl_sf_bessel_Jn(i,PARAMS.kr);
                        bessel[lmax-i] = -bessel[lmax+i];
                    }
                }
                PARAMS.kr = -sintheta*(iz+0.5);

                w0[j].fourier[iz] += Int_weight100[st]*w0func(bessel,PARAMS)*exponent;
                w3[j].fourier[iz] += Int_weight100[st]*w3func(bessel,PARAMS)*exponent;
                for ( i = 0; i <= lmax; i++){
                    PARAMS.l = i;
                    for( v = 0; v<=2*i; v++){
                        PARAMS.m = v-i;
                        //Wig = 1/Wignerd(atheta, i ,0, PARAMS.m);
                        Wig = Wignerd(atheta, i ,0, PARAMS.m);
                        w1[i][j].fourier[iz] += Wig*Int_weight100[st]*w1func(bessel,PARAMS)*exponent;
                        w2[i][j].fourier[iz] += Wig*Int_weight100[st]*w2func(bessel,PARAMS)*exponent;
                    }
                }

                params.kr = sintheta*(iz+0.5-n_bins)*params.RR;
                PARAMS.kr = sintheta*(iz+0.5-n_bins)*PARAMS.RR;
                kz = costheta*(iz+0.5-n_bins);

                exponent = cexp(-I*kz*z);
            
                bessel[lmax] = gsl_sf_bessel_J0(params.kr);  
                bessel[lmax+1] = -gsl_sf_bessel_J1(params.kr);  
                bessel[lmax-1] = -bessel[lmax+1];  

                for( i=2; i <= lmax; i++){
                     if( i % 2 == 0){
                        bessel[lmax+i] = gsl_sf_bessel_Jn(i,params.kr);
                        bessel[lmax-i] = bessel[lmax+i];
                    }else{ 
                        bessel[lmax+i] = -gsl_sf_bessel_Jn(i,params.kr);
                        bessel[lmax-i] = -bessel[lmax+i];
                    }
                }
                params.kr = -sintheta*(iz+0.5-n_bins);

                w0[j].fourier[iz] += Int_weight100[st]*w0func(bessel,params)*exponent;
                w3[j].fourier[iz] += Int_weight100[st]*w3func(bessel,params)*exponent;
                for ( i = 0; i <= lmax; i++){
                    params.l = i;
                    for( v = 0; v<=2*i; v++){
                        params.m = v-i;
                        //Wig = 1/Wignerd(atheta, i ,0, params.m);
                        Wig = Wignerd(atheta, i ,0, params.m);
                        w1[i][j].fourier[iz] += Wig*Int_weight100[st]*w1func(bessel,params)*exponent;
                        w2[i][j].fourier[iz] += Wig*Int_weight100[st]*w2func(bessel,params)*exponent;
                    }
                }
                exponent = cexp(I*kz*z);
            
                bessel[lmax] = gsl_sf_bessel_J0(PARAMS.kr);  
                bessel[lmax+1] = -gsl_sf_bessel_J1(PARAMS.kr);  
                bessel[lmax-1] = -bessel[lmax+1];  

                for( i=2; i <= lmax; i++){
                     if( i % 2 == 0){
                        bessel[lmax+i] = gsl_sf_bessel_Jn(i,PARAMS.kr);
                        bessel[lmax-i] = bessel[lmax+i];
                    }else{ 
                        bessel[lmax+i] = -gsl_sf_bessel_Jn(i,PARAMS.kr);
                        bessel[lmax-i] = -bessel[lmax+i];
                    }
                }
                PARAMS.kr = -sintheta*(iz+0.5-n_bins);

                w0[j].fourier[iz] += Int_weight100[st]*w0func(bessel,PARAMS)*exponent;
                w3[j].fourier[iz] += Int_weight100[st]*w3func(bessel,PARAMS)*exponent;
                for ( i = 0; i <= lmax; i++){
                    PARAMS.l = i;
                    for( v = 0; v<=2*i; v++){
                        PARAMS.m = v-i;
                        //Wig = 1/Wignerd(atheta, i ,0, PARAMS.m);
                        Wig = Wignerd(atheta, i ,0, PARAMS.m);
                        w1[i][j].fourier[iz] += Wig*Int_weight100[st]*w1func(bessel,PARAMS)*exponent;
                        w2[i][j].fourier[iz] += Wig*Int_weight100[st]*w2func(bessel,PARAMS)*exponent;
                    }
                }
               /* if(iz==0){
                    for ( i = 1; i <= lmax; i++){
                        w1[i][j].fourier[iz] = 0;
                        w2[i][j].fourier[iz] = 0;
                    }
                }*/
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

        n0.fourier[iz] += 0.5*rho[0].fourier[iz]*w0[0].fourier[iz];
        n3.fourier[iz] += 0.5*rho[0].fourier[iz]*w3[0].fourier[iz];
        for ( j = 1; j < 100; j++){
            n0.fourier[iz] += rho[j].fourier[iz]*w0[j].fourier[iz];
            n3.fourier[iz] += rho[j].fourier[iz]*w3[j].fourier[iz];
        }
        n0.fourier[iz] += 0.5*rho[100].fourier[iz]*w0[100].fourier[iz];
        n3.fourier[iz] += 0.5*rho[100].fourier[iz]*w3[100].fourier[iz];

        
        n0.fourier[iz] *= 0.02;
        n3.fourier[iz] *= 0.02;
        //n0.fourier[iz] *= 0.01*M_PI;
        //n3.fourier[iz] *= 0.01*M_PI;

        for ( i = 0; i <= lmax; i++){
            n1[i].fourier[iz] = 0;
            n2[i].fourier[iz] = 0;

            n1[i].fourier[iz] += 0.5*rho[0].fourier[iz]*w1[i][0].fourier[iz];
            n2[i].fourier[iz] += 0.5*rho[0].fourier[iz]*w2[i][0].fourier[iz];
//            if(i==4 && iz == 2) std::cout << 0 << ": " << creal(rho[0].fourier[iz]) << " " << creal(w1[i][0].fourier[iz]) << " " << creal(rho[0].fourier[iz]*w1[i][0].fourier[iz]) << std::endl;
            for ( j = 1; j < 100; j++){
                n1[i].fourier[iz] += rho[j].fourier[iz]*w1[i][j].fourier[iz];
                n2[i].fourier[iz] += rho[j].fourier[iz]*w2[i][j].fourier[iz];
 //               if(i==4 && iz == 2) std::cout << j << ": " << creal(rho[j].fourier[iz]) << " " << creal(w1[i][j].fourier[iz]) << " " << creal(rho[j].fourier[iz]*w1[i][j].fourier[iz]) << std::endl;
            }
            n1[i].fourier[iz] += 0.5*rho[100].fourier[iz]*w1[i][100].fourier[iz];
            n2[i].fourier[iz] += 0.5*rho[100].fourier[iz]*w2[i][100].fourier[iz];
//            if(i==4 && iz == 2) std::cout << 100 << ": " << creal(rho[100].fourier[iz]) << " " << creal(w1[i][100].fourier[iz]) << " " << creal(rho[100].fourier[iz]*w1[i][100].fourier[iz]) << std::endl;

            n1[i].fourier[iz] *= 0.02;
            n2[i].fourier[iz] *= 0.02;
            //n1[i].fourier[iz] *= 0.01*M_PI;
            //n2[i].fourier[iz] *= 0.01*M_PI;
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

    for ( j = 0; j <= 100; j++){
        sc0[j].invfft();
        for ( i = 0; i <= lmax; i++){
            sc1[i][j].invfft();
            sc2[i][j].invfft();
        }
        sc3[j].invfft();
    }
}


void CalcC(){

    for( iz = 0; iz < n_bins; iz++){ 
        for ( j = 0; j <= 100; j++){
            c1[j].real[iz] = -sc0[j].real[iz];
            for ( i = 0; i <= lmax; i++){
                c1[j].real[iz] -= sc1[i][j].real[iz];
                c1[j].real[iz] -= sc2[i][j].real[iz];
            }
            c1[j].real[iz] -= sc3[j].real[iz];
        }
    }
    
}

void CalcRhoN(){

    prev = now;
    now = 0;
    
    for(j = 0; j <= 100; j++){
    
        //atheta = j*M_PI*0.01; 
        //theta = cos(atheta); 
        theta = j*0.02-1;

        b3 = WallDistance(theta);//Pear
        //b3=RII;

        theta = -theta;
        
        b2 = WallDistance(theta);//Pear
        //b2=RII;

        for ( iz = 0; iz < n_bins; iz++){
            z = DeltaR*(iz+0.5);
            if(z <= H2){ 
                if( z <= b2 ) rhon[j].real[iz] = 0.;
                else rhon[j].real[iz] = exp(creal(c1[j].real[iz]));// + rho_sum[(ix]));
            }else{
                z = H - z;
                if( z <= b3 ) rhon[j].real[iz] = 0.;
                else rhon[j].real[iz] = exp(creal(c1[j].real[iz]));// + rho_sum[(ix]));
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
        //rho_fin.real[iz] *= 0.01*M_PI;
        
        rho_sum += creal(rho_fin.real[iz]);
    }
    sclf = rho_sum_o/rho_sum;
    now = sclf; 
    printf("%i %f %f\n", st, now, sclf);
    
    for( iz =0; iz < n_bins; iz++){ 
        rho_fin.real[iz] *= sclf; 
        for(j = 0; j <= 100; j++){ 
            rhon[j].real[iz] *= sclf;
            rho[j].real[iz] = rho[j].real[iz]*(1-alpha) + rhon[j].real[iz]*alpha;     
        }
    }
}
