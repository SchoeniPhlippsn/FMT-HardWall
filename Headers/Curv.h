
double kappa1 (double t, double aa) {  //Pear
    double nenner = a1*(1-2*t+2*t*t)-aa*(1-6*t+6*t*t);
    nenner = sqrt(16*b1*b1*(1-2*t)*(1-2*t)+9*nenner*nenner);
    nenner = nenner*nenner*nenner;
    double zaehler = 48*b1*(a1*(1 - t)*t + aa*(1 - 3*t + 3*t*t));
    double zn = zaehler/nenner;
    if (zn < 0) return -zn;
    else return zn;
}

/*double kappa1 (double t, double aa){
    return 1/RII;
}*/

double kappa2 (double t, double aa) { //Pear
    double nenner = a1*(1-2*t+2*t*t)-aa*(1-6*t+6*t*t);
    nenner = (3*aa*(1-t)*t+a1*(1-t+t*t))*sqrt(16*b1*b1*(1-2*t)*(1-2*t)+9*nenner*nenner);
    double zaehler = 4*b1;
    double zn = zaehler/nenner;
    if (zn < 0) return -zn;
    else return zn;
}

/*double kappa2 (double t, double aa){
    return 1/RII;
}*/

double gaussC (double t, double aa){
    return kappa1(t,aa)*kappa2(t,aa);
}

double meanC (double t, double aa){
    return 0.5*(kappa1(t,aa)+kappa2(t,aa));
}

double diffC (double t, double aa){
    if(kappa1(t,aa)>kappa2(t,aa)) return 0.5*(kappa1(t,aa)-kappa2(t,aa));
    else return 0.5*(kappa2(t,aa)-kappa1(t,aa));
    //return 0.5*(kappa1(t,aa)-kappa2(t,aa));
}
