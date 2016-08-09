
double kappa1 (double t, double aa) {
    double nenner = 16*b1*b1*(1-2*t)*(1-2*t)+9*a1*a1*(1-2*t+2*t*t)*(1-2*t+2*t*t)+9*aa*aa*(1-6*t+6*t*t)*(1-6*t+6*t*t)-18*a1*aa*(1-8*t+20*t*t-24*t*t*t+12*t*t*t*t);
    nenner = sqrt(nenner);
    nenner = nenner*nenner*nenner;
    double zaehler = 48*b1*(a1*(1 - t)*t + aa*(1 - 3*t + 3*t*t));
    return zaehler/nenner;
}

double kappa2 (double t, double aa) {
    double nenner = a1*(1-2*t+2*t*t)-aa*(1-6*t+6*t*t);
    nenner = (3*aa*(1-t)*t+a1*(1-t+t*t))*sqrt(16*b1*b1*(1-2*t)*(1-2*t)+9*nenner*nenner);
    double zaehler = 4*b1;
    return zaehler/nenner;
}

double gaussC (double t, double aa){
    return kappa1(t,aa)*kappa2(t,aa);
}

double meanC (double t, double aa){
    return 0.5*(kappa1(t,aa)+kappa2(t,aa));
}

double diffC (double t, double aa){
    return 0.5*(kappa1(t,aa)-kappa2(t,aa));
}
