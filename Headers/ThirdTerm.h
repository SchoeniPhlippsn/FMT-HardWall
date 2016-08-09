
double Ofunction(int l, int m){
    if(l==0 & m==0){
        return 2;
    }else{
        if( abs(l+m) % 2 == 1 || m > l || m == 0){
            return 0;
        }
        else{
            //if(m==0) return l;
            //else{
                double Of=pow(2,m-1)*m;
                Of = Of*gsl_sf_gamma(l/2.0)*gsl_sf_gamma((m+l+1)/2.0);
                Of = Of*gsl_sf_gammainv((l+3)/2.0)*gsl_sf_gammainv((l-m+2)/2.0);
                
                if(m % 2 == 0 && m > 0) return Of;
                return -Of;
            //}
        }
    }
}

double IInt(int l, int m)
{
    double N = sqrt(factorial(l-m)/factorial(l+m));
    N = N/(4*l+2);

    double If = N*(Ofunction(l-1,m+1)-Ofunction(l+1,m+1));
    
    return If;
}

double IIntTR(int l, int m)
{
    if(l== 0 && m==0) return 2.0/3;
    else{
        if(l==2 && m==0) return -2.0/15;
        else{
            if( l==2 && abs(m) == 2) return 4.0/(5*sqrt(6.0));
            return 0;
        }
    }
}

double JInt(int m)
{
    if(abs(m)==1) return -0.25;
    else{
        int sgn;
        
        if(abs(m)%2==0) sgn=1;
        else sgn = -1;

        double Jf = (1.0-sgn*2.0)/(3.0*(m*m-1.0));
        return Jf;
    }
}

double JIntTR(int m)
{
    if(m==0) return 0.5;
    else{
        if(abs(m)==2) return -0.25;
        else return 0;
    }
}

bool SecondTerm(int l1, int l2, int m1, int m2){
    if(l1==l2 && m1==m2) return 1;
    else return 0;
}

double ThirdTerm(int l1, int l2, int l3, int m1, int m2, int m3){
    double C0 = (2*l1+1)*(2*l2+1)*(2*l3+1)/(8.0*M_PI)*gsl_sf_coupling_3j(2*l1,2*l2,2*l3,2*m1,2*m2,2*m3);

    int minl = min(l1,l3); 

    int p = -minl;
    double sum=0;
    for(; p <= minl; p++){
        sum += gsl_sf_coupling_3j(2*l1,2*l2,2*l3,-2*p,0,2*p)*IInt(l1,-p)*IInt(l3,p)*JInt(p); 
    }
    return sum*C0;
}

double ThirdTermTR(int l1, int l2, int l3, int m1, int m2, int m3){
    //double C0 = (2*l1+1)*(2*l2+1)*(2*l3+1)/(48.0*M_PI)*gsl_sf_coupling_3j(2*l1,2*l2,2*l3,2*m1,2*m2,2*m3);
    
    //double C0 = gsl_sf_coupling_3j(2*l1,2*l2,2*l3,2*m1,2*m2,2*m3);
    
    double sum=0;
    //if( l1==0 && l2==0 && l3==0 ) sum = C0/(3*48*M_PI); 
    if( l1==0 && l2==0 && l3==0 ) sum = 1/(3*48*M_PI); 
    //if(l1==2 && l2 == 2 && l3 == 2) sum = -2.788866755*C0/(48*M_PI); 
    if(l1==2 && l2 == 2 && l3 == 2) sum = 2/(3*48*M_PI); 
    //if(l1==2 && l2 == 2 && l3 == 0) sum = -0.745355992*C0/(48*M_PI); 
    //if(l1==2 && l2 == 0 && l3 == 2) sum = -0.745355992*C0/(48*M_PI); 
    //if(l1==0 && l2 == 2 && l3 == 2) sum = -0.745355992*C0/(48*M_PI);  
    if(l1==2 && l2 == 2 && l3 == 0) sum = -1/(3*48*M_PI); 
    if(l1==2 && l2 == 0 && l3 == 2) sum = -1/(3*48*M_PI); 
    if(l1==0 && l2 == 2 && l3 == 2) sum = -1/(3*48*M_PI);  
    return sum;
}
