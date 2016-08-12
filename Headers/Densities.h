double w0func (double bessel[2*lmax+1], struct paras p) {

    double f  = 0.5*p.gauss*p.AblRR*bessel[lmax];//*gsl_sf_bessel_J0(p.kr*p.RR);
    return f;
}


double w1func (double bessel[2*lmax+1], struct paras p) {
    double f;
    if(p.l==0) f = 0.5*p.AblRR*Wignerd(p.theta,p.l,p.m,0)*p.mean*bessel[lmax];//*gsl_sf_bessel_J0(p.kr*p.RR);
    if(p.l==1) f = -0.5*p.AblRR*Wignerd(p.theta,p.l,p.m,0)*p.mean*bessel[lmax+p.m];//*gsl_sf_bessel_Jn(p.m,p.kr*p.RR);
    if(p.l>1){
        double Wig1 = Wignerd(p.theta,p.l,p.m,2)+Wignerd(p.theta,p.l,p.m,-2);
        f = 2*M_PI*p.AblRR*Wig1*p.diff*SHN(p.l,2)*SHN(p.l,0)*bessel[lmax+p.m];//*gsl_sf_bessel_Jn(p.m,p.kr*p.RR);
        if( p.l % 2 == 1) f *= -1;
    }
    return f;
}

double w2func (double bessel[2*lmax+1], struct paras p) {
    
    double f = 2*M_PI*p.AblRR*Wignerd(p.theta,p.l,p.m,0)*bessel[lmax+p.m];//*gsl_sf_bessel_Jn(p.m,p.kr*p.RR);
    return f;
}

double w3func (double bessel[2*lmax+1], struct paras p) {
    double f;
    if(p.kr>1e-7){
        f=2*M_PI*p.RR*bessel[lmax+1];//*gsl_sf_bessel_J1(p.kr*p.RR);
        return f/(b1*p.kr);
    }else{
        f=M_PI*p.RR*p.RR;
        return f;
    }
}
