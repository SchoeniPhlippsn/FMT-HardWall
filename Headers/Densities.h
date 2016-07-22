
double w0func (double z, struct paras p) {
    double kr = p.kr;
    double a1 = p.x_para1;
    double a2 = p.x_para2;
    double a3 = p.x_para3;
    double b  = p.y_para;


    double f;
    if(z>0){ 
        f = 0.5*gaussC(z,a1,a2,b)*RRAbl2(z,a1,a2,b);
    }else{
        z *= -1;
        f = 0.5*gaussC(z,a1,a3,b)*RRAbl2(z,a1,a3,b);
    }
    
    return f;
}

double w1func (double z, struct paras p) {
    double kr = p.kr;
    double a1 = p.x_para1;
    double a2 = p.x_para2;
    double a3 = p.x_para3;
    double b = p.y_para;
    int l = p.l;
    int m = p.m;
    
    double theta;
    double f;
    if(l==0){ 
        if(z>0){ 
            theta = getTheta(z,a1,a2,b);
            f = 0.5*Wignerd(theta,l,m,0)*meanC(z,a1,a2,b)*RRAbl2(z,a1,a2,b);
        }else{
            z *= -1;
            theta = getTheta(z,a1,a3,b);
            f = 0.5*Wignerd(theta,l,m,0)*meanC(z,a1,a3,b)*RRAbl2(z,a1,a3,b);
        }
    }
    if(l==1){
        if(z>0){
            theta = getTheta(z,a1,a2,b);
            f = -0.5*Wignerd(theta,l,m,0)*meanC(z,a1,a2,b)*RRAbl2(z,a1,a2,b);
        }else{ 
            z *= -1;
            theta = getTheta(z,a1,a3,b);
            f = -0.5*Wignerd(theta,l,m,0)*meanC(z,a1,a3,b)*RRAbl2(z,a1,a3,b);
        }
    }
    if(l>1){
        if(z>0){ 
            theta = getTheta(z,a1,a2,b);
            double Wig = Wignerd(theta,l,m,2)+Wignerd(theta,l,m,-2);
            f = 2*M_PI*Wig*diffC(z,a1,a2,b)*SHN(l,2)*SHN(l,0)*RRAbl2(z,a1,a2,b);
        }else{
            z *= -1;
            theta = getTheta(z,a1,a3,b);
            double Wig = Wignerd(theta,l,m,2)+Wignerd(theta,l,m,-2);
            f = 2*M_PI*Wig*diffC(z,a1,a3,b)*SHN(l,2)*SHN(l,0)*RRAbl2(z,a1,a3,b);
        }
        if( l % 2 == 1) f *= -1;
    }
    return f;
}

double w2func (double z, struct paras p) {
    double kr = p.kr;
    double a1 = p.x_para1;
    double a2 = p.x_para2;
    double a3 = p.x_para3;
    double b = p.y_para;
    int l = p.l;
    int m = p.m;
    
    double f;
    double theta;
    if(z>0){ 
        theta = getTheta(z,a1,a2,b);
        f = M_PI*Wignerd(theta,l,m,0)*RRAbl2(z,a1,a2,b);
    }else{ 
        z *= -1;
        theta = getTheta(z,a1,a3,b);
        f = M_PI*Wignerd(theta,l,m,0)*RRAbl2(z,a1,a3,b);
    }

    return f;
}

double w3func (double z, struct paras p) {
    double kr = p.kr;
    double a1 = p.x_para1;
    double a2 = p.x_para2;
    double a3 = p.x_para3;
    double b = p.y_para;
    
    double f;
    if(z>0) f = M_PI*RR(z,a1,a2,b)*RR(z,a1,a2,b);
    else{ 
        z *= -1;
        f = M_PI*RR(z,a1,a3,b)*RR(z,a1,a3,b);
    }
     
    return f;
}
