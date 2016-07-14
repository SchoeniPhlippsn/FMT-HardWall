
void Phi0(){
    
    PHI0.real[(ix*n_bins+iy)*n_bins+iz] = -log(on3);

}

void Phi1(){

    PHI1[i][i+j].real[(ix*n_bins+iy)*n_bins+iz] = n2[i][i+j].real[(ix*n_bins+iy)*n_bins+iz]*inv_on3;

}

void Phi2(){

    PHI2[i][i+j].real[(ix*n_bins+iy)*n_bins+iz] = n1[i][i+j].real[(ix*n_bins+iy)*n_bins+iz]*inv_on3;

    if(i==0){
        PHI2[i][j].real[(ix*n_bins+iy)*n_bins+iz] += 3*9*ThirdTermTR(0,0,0,0,0,0)*n2[0][0].real[(ix*n_bins+iy)*n_bins+iz]*n2[0][0].real[(ix*n_bins+iy)*n_bins+iz]*inv_on3*inv_on3;
        l2=2;
        l3=2;  
        
        for ( m2 = -l2; m2 <= l2; m2++){
            for ( m3 = -l3; m3 <= l3; m3++){
                Coupling3 = ThirdTermTR(0,l2,l3,0,m2,m3);
                if( Coupling3 > 1e-6 || Coupling3 < -1e-6){
                    PHI2[i][i+j].real[(ix*n_bins+iy)*n_bins+iz] += 9*ThirdTermTR(0,l2,l3,0,m2,m3)*n2[l2][l2+m2].real[(ix*n_bins+iy)*n_bins+iz]*n2[l3][l3+m3].real[(ix*n_bins+iy)*n_bins+iz]*inv_on3*inv_on3;   
                }
            }
        } 
    }else{
        if(i==2){
            l2=2;
            l3=0;
            while(l2<4){
                for ( m2 = -l2; m2 <= l2; m2++){
                    for ( m3 = -l3; m3 <= l3; m3++){
                        Coupling3 = ThirdTermTR(i,l2,l3,j,m2,m3);
                        if( Coupling3 > 1e-6 || Coupling3 < -1e-6){
                            factor = 1;
                            if(i==l2 && j==m2) factor++; 
                            if(i==l3 && j==m3) factor++; 
                            
                            PHI2[i][i+j].real[(ix*n_bins+iy)*n_bins+iz] += 9*factor*Coupling3*n2[l2][l2+m2].real[(ix*n_bins+iy)*n_bins+iz]*n2[l3][l3+m3].real[(ix*n_bins+iy)*n_bins+iz]*inv_on3*inv_on3;
                        }
                    }
                }

                if(l3 == 0){
                    l3=2;
                    l2=0;
                }else l2+=2;
            }
        }
    }
}

void Phi3(){
    
    PHI3.real[(ix*n_bins+iy)*n_bins+iz] = n0.real[(ix*n_bins+iy)*n_bins+iz]*inv_on3;
    l2=0;
    l3=0;
    for ( l1 = 0; l1 <= lmax; l1++){
        for( m1 = -l1; m1<=l1; m1++){
            PHI3.real[(ix*n_bins+iy)*n_bins+iz] += n1[l1][l1+m1].real[(ix*n_bins+iy)*n_bins+iz]*n2[l1][l1+m1].real[(ix*n_bins+iy)*n_bins+iz]*inv_on3*inv_on3;
            if(l1==1 || l1>2) continue;
            if(l1==0){
                l2=0;
                l3=0;
            }else{
                l2=2;
                l3=2;
            }
            while(l2 < 3){
                for( m2 = -l2; m2<=l2; m2++){
                    for( m3 = -l3; m3<=l3; m3++){
                        Coupling3 = ThirdTermTR(l1,l2,l3,m1,m2,m3);
                        if( Coupling3 > 1e-6 || Coupling3 < -1e-6){
                            PHI3.real[(ix*n_bins+iy)*n_bins+iz] += 2*9*Coupling3*n2[l1][l1+m1].real[(ix*n_bins+iy)*n_bins+iz]*n2[l2][l2+m2].real[(ix*n_bins+iy)*n_bins+iz]*n2[l3][l3+m3].real[(ix*n_bins+iy)*n_bins+iz]*inv_on3*inv_on3*inv_on3;
                        }
                    }
                }
                l2+=2;
                l3+=2;
            }
        }
    }
}


