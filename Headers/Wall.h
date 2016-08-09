
double WallDistance(double theta){
    
    if( theta == 1 || theta == -1) return b1;
    else{
        
        double co = theta;
        double si = sqrt(1-co*co);
         
        double T21 = (18*a3 - 6*a1)*si;
        double T22 = (6*a1-18*a3)*si + 8*b1*co;
        double T23 = -4*b1*co + 3*(a3 - a1)*si;

        double t2 = (-T22-sqrt(T22*T22-4*T21*T23))/(2*T21);
        double x,y;
        if( t2 >= 0){
            x = a1*(1-t2)*(1-t2)*(1-t2) + 3*a3*(1-t2)*(1-t2)*t2 - 3*a3*(1-t2)*t2*t2 - a1*t2*t2*t2;
            y = -4*b1*t2*(1-t2);

            return si*x + co*y;
        }else{
            double T11 = (18*a2-6*a1)*si;
            double T12 = (6*a1-18*a2)*si - 8*b1*co;
            double T13 = 4*b1*co + 3*(a2 - a1)*si;

            double t1 = (-T12-sqrt(T12*T12-4*T11*T13))/(2*T11);
            
            x = a1*(1-t1)*(1-t1)*(1-t1) + 3*a2*(1-t1)*(1-t1)*t1 - 3*a2*(1-t1)*t1*t1 - a1*t1*t1*t1;
            y = 4*b1*t1*(1-t1);

            return si*x + co*y;
        }
    }
}
