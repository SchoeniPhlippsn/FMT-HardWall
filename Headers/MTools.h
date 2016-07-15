
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
        
        if(abs(m+n)/2 % 2 == 1) d *= -1;
    
    }else{

        for(k=1; k <= l; k++){
               i = getNumber(l,m,k);
               j = getNumber(l,k,n);

               d += 2*sgn*WList[i]*WList[j]*sin(k*theta);
               sgn *= -1;
        } 

        if(abs(m+n+1)/2 % 2 == 1) d *= -1;
    
    }
    
    return d;   
} 

double WignerdInt(int l, int m, int n){
    double d=0;
    int i,j,k;
    if(abs(m+n) % 2 == 0){
        i = getNumber(l,m,0);
        j = getNumber(l,0,n);

        d += WList[i]*WList[j]*M_PI;
        
        if(abs(m+n)/2 % 2 == 1) d *= -1;

    }else{
        for(k=1; k <= l; k+=2){
               i = getNumber(l,m,k);
               j = getNumber(l,k,n);

               d -= 4*WList[i]*WList[j]/k;
        } 

        if(abs(m+n+1)/2 % 2 == 1) d *= -1;
    }
    
    return d;   
} 

double SHN(int l, int m){
   
    int l_m = factorial(l-m);
    int lm = factorial(l+m);

    double d = (2*l+1)/(4*M_PI);
    if(m >0) return sqrt(l_m*d/lm);
    else return pow(-1,m)*sqrt(lm*d/l_m);
}

double t (double z, double b){
    double tt = 0.5*(1-sqrt(1-z/b));
    return tt;
}

double RR (double z, double a1, double a2, double b) {   
    double tt = t(z,b);
    double rr = a1*(1-tt)*(1-tt)*(1-tt) + 3*a2*(1-tt)*(1-tt)*tt - 3*a2*(1-tt)*tt*tt - a1*tt*tt*tt;
    return rr;
}

double getTheta( double z, double a1, double a2, double b) {
    double tt = t(z,b);
    double nz = 3.0*(a1-a2)*(1-tt)*(1-tt) + 12.0*a2*tt*(1-tt) + 3.0*(a1-a2)*tt*tt;;
    if ( a1 < a2 ) nz *= -1; 
   
    double nrho = 4*b*(1-2*tt);

    double nn = sqrt(nz*nz + nrho*nrho);
    nz = nz/nn;

    return acos(nz);
}
