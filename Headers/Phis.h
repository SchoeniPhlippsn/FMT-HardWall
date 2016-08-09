
void Phi0(){
    
    PHI0.real[iz] = -log(on3);

}

void Phi1(){

    PHI1[i].real[iz] = n2[i].real[iz]*inv_on3;

}

void Phi2(){

    PHI2[i].real[iz] = n1[i].real[iz]*inv_on3;

    if(i==0){
        PHI2[i].real[iz] += 3*9*ThirdTermTR(0,0,0,0,0,0)*n2[0].real[iz]*n2[0].real[iz]*inv_on3*inv_on3;
        PHI2[i].real[iz] += 9*ThirdTermTR(0,2,2,0,0,0)*n2[2].real[iz]*n2[2].real[iz]*inv_on3*inv_on3;   
    
    }else{
        if(i==2){
            PHI2[i].real[iz] += 2*9*ThirdTermTR(2,2,0,0,0,0)*n2[2].real[iz]*n2[0].real[iz]*inv_on3*inv_on3;   
            PHI2[i].real[iz] += 2*9*ThirdTermTR(2,0,2,0,0,0)*n2[0].real[iz]*n2[2].real[iz]*inv_on3*inv_on3;   
            PHI2[i].real[iz] += 3*9*ThirdTermTR(2,2,2,0,0,0)*n2[2].real[iz]*n2[2].real[iz]*inv_on3*inv_on3;   
        }
    }
}

void Phi3(){
    
    PHI3.real[iz] = n0.real[iz]*inv_on3;
        
    l2=0;
    l3=0;
    for ( l1 = 0; l1 <= lmax; l1++){
        PHI3.real[iz] += n1[l1].real[iz]*n2[l1].real[iz]*inv_on3*inv_on3;
    }
    
    PHI3.real[iz] += 2*9*ThirdTermTR(0,0,0,0,0,0)*n2[0].real[iz]*n2[0].real[iz]*n2[0].real[iz]*inv_on3*inv_on3*inv_on3;
    PHI3.real[iz] += 2*9*ThirdTermTR(2,2,0,0,0,0)*n2[2].real[iz]*n2[2].real[iz]*n2[0].real[iz]*inv_on3*inv_on3*inv_on3;
    PHI3.real[iz] += 2*9*ThirdTermTR(2,0,2,0,0,0)*n2[2].real[iz]*n2[0].real[iz]*n2[2].real[iz]*inv_on3*inv_on3*inv_on3;
    PHI3.real[iz] += 2*9*ThirdTermTR(0,2,2,0,0,0)*n2[0].real[iz]*n2[2].real[iz]*n2[2].real[iz]*inv_on3*inv_on3*inv_on3;
    PHI3.real[iz] += 2*9*ThirdTermTR(2,2,2,0,0,0)*n2[2].real[iz]*n2[2].real[iz]*n2[2].real[iz]*inv_on3*inv_on3*inv_on3;

}


