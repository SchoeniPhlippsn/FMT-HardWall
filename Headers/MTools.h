
int min(int a, int b){
    if(a < b) return a;
    else return b;
}


double factorial ( int n ){
    int nn=1;
    int k;
    for(k = 1; k <= n; k++) nn *= k;
    return nn;
}


int getNumber(int l, int m, int n){
    int k;
    int numb=0;
    for( k=l; k!=0; k--) numb += (2*k-1)*(2*k-1);
    numb += (2*l+1)*(m+l);
    numb += n+l;
    return numb;
}


double Wignerd(double theta, int l, int m, int n){
    double d=0;
    int i,j,k;
    int sgn = -1;
    if(abs(m+n) % 2 == 0){
        i = getNumber(l,m,0);
        j = getNumber(l,0,n);

        d += WList[i]*WList[j];

        for(k=1; k <= l; k++){
           i = getNumber(l,m,k);
           j = getNumber(l,k,n);

           d += 2*sgn*WList[i]*WList[j]*cos(k*theta);
           sgn *= -1;
        } 
        
        if(abs(m+n) % 4 == 2) d *= -1;
    
    }else{

        for(k=1; k <= l; k++){
               i = getNumber(l,m,k);
               j = getNumber(l,k,n);

               d += 2*sgn*WList[i]*WList[j]*sin(k*theta);
               sgn *= -1;
        } 

        if(abs(m+n+1) % 4 == 2) d *= -1;
    
    }
    
    return d;   
} 

double SHN(int l, int m){
   
    int l_m = factorial(l-m);
    int lm = 1/factorial(l+m);

    double d = 0.25*InvPi*(2*l+1);
    if(m >0) return sqrt(l_m*d*lm);
    else return pow(-1,m)*sqrt(lm*d*l_m);
}

double RR (double t, double aa) {   
    double rr = a1*(1-t)*(1-t)*(1-t) + 3*aa*(1-t)*(1-t)*t - 3*aa*(1-t)*t*t - a1*t*t*t;
    return rr;
}

double RRAbl2 (double z, double aa) {   
    double term1 = (a1*b1-0.25*a1*z+0.75*aa*z)*b1_inv*b1_inv*b1_inv;
    double term2 = b1*b1*b1*b1 - b1*b1*b1*z;
    double term3 = -0.75*a1*b1+0.75*aa*b1+0.375*a1*z-1.125*aa*z;
    term2 += term3*term3;
    return term1*sqrt(term2);
}

double getTheta( double t, double aa) {
    double nz = 3.0*(a1-aa)*(1-t)*(1-t) + 12.0*aa*t*(1-t) + 3.0*(a1-aa)*t*t;
    if ( a1 < aa ) nz *= -1; 
   
    double nrho = 4*b1*(1-2*t);

    double nn = sqrt(nz*nz + nrho*nrho);
    nz = nz/nn;

    return acos(nz);
}
