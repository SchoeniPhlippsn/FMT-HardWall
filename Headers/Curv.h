
double kappa1 (double z, double a1, double a2, double b) {
    double tt = t(z,b);
    double nenner = 16*b*b*(1-2*tt)*(1-2*tt)+9*a1*a1*(1-2*tt+2*tt*tt)*(1-2*tt+2*tt*tt)+9*a2*a2*(1-6*tt+6*tt*tt)*(1-6*tt+6*tt*tt)-18*a1*a2*(1-8*tt+20*tt*tt-24*tt*tt*tt+12*tt*tt*tt*tt);
    nenner = sqrt(nenner);
    nenner = nenner*nenner*nenner;
    double zaehler = 48*b*(a1*(1 - tt)*tt + a2*(1 - 3*tt + 3*tt*tt));
    return zaehler/nenner;
}

double kappa2 (double z, double a1, double a2, double b) {
    double tt = t(z,b);
    double nenner = a1*(1-2*tt+2*tt*tt)-a2*(1-6*tt+6*tt*tt);
    nenner = (3*a2*(1-tt)*tt+a1*(1-tt+tt*tt))*sqrt(16*b*b*(1-2*tt)*(1-2*tt)+9*nenner*nenner);
    double zaehler = 4*b;
    return zaehler/nenner;
}

double gaussC (double z, double a1, double a2, double b){
    return kappa1(z,a1,a2,b)*kappa2(z,a1,a2,b);
}

double meanC (double z, double a1, double a2, double b){
    return 0.5*(kappa1(z,a1,a2,b)+kappa2(z,a1,a2,b));
}

double diffC (double z, double a1, double a2, double b){
    return 0.5*(kappa1(z,a1,a2,b)-kappa2(z,a1,a2,b));
}
