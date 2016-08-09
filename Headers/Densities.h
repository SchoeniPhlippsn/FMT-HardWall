double w0func (double zz, struct paras p) {
    double kr = p.kr;
    double aa = a2;

    if(zz < 0){
        zz *= -1;
        aa=a3;
    }
     
    double tt = 0.5*(1-sqrt(1-zz*b1_inv));
    double f  = 0.5*gaussC(tt,aa)*RRAbl2(zz,aa)*gsl_sf_bessel_J0(kr*RR(tt,aa));
    return f;
}

double w1func (double zz, struct paras p) {
    double kr = p.kr;
    int l = p.l;
    int m = p.m;
    double aa = a2;

    if(zz < 0){
        zz *= -1;
        aa=a3;
    }
     
    double tt = 0.5*(1-sqrt(1-zz*b1_inv));
    
    double theta = getTheta(tt,aa);
    double f;
    if(l==0) f = 0.5*RRAbl2(zz,aa)*Wignerd(theta,l,m,0)*gsl_sf_bessel_J0(kr*RR(tt,aa))*meanC(tt,aa);
    if(l==1) f = -0.5*RRAbl2(zz,aa)*Wignerd(theta,l,m,0)*gsl_sf_bessel_Jn(m,kr*RR(tt,aa))*meanC(tt,aa);
    if(l>1){
        double Wig1 = Wignerd(theta,l,m,2)+Wignerd(theta,l,m,-2);
        f = 2*M_PI*RRAbl2(zz,aa)*Wig1*gsl_sf_bessel_Jn(m,kr*RR(tt,aa))*diffC(tt,aa)*SHN(l,2)*SHN(l,0);
        if( l % 2 == 1) f *= -1;
    }
    return f;
}

double w2func (double zz, struct paras p) {
    double kr = p.kr;
    int l = p.l;
    int m = p.m;
    double aa = a2;

    if(zz < 0){
        zz *= -1;
        aa=a3;
    }
     
    double tt = 0.5*(1-sqrt(1-zz*b1_inv));
    
    double theta = getTheta(tt,aa);
    double f = 2*M_PI*RRAbl2(zz,aa)*Wignerd(theta,l,m,0)*gsl_sf_bessel_Jn(m,kr*RR(tt,aa));
    return f;
}

double w3func (double zz, struct paras p) {
    double kr = p.kr;
    double aa = a2;

    if(zz < 0){
        zz *= -1;
        aa=a3;
    }
    
    double tt = 0.5*(1-sqrt(1-zz*b1_inv));
    double f;
    if(kr>1e-7){
        f=2*M_PI*RR(tt,aa)*gsl_sf_bessel_J1(kr*RR(tt,aa));
        return f/(b1*kr);
    }else{
        f=M_PI*RR(tt,aa)*RR(tt,aa);
        return f;
    }
}
