
double WallDistance(double theta){
    
    if( j == 0 || j == 100) return b1;
    else{ 
        double T21 = 18*a3*sin(theta) - 6*a1*sin(theta);
        double T22 = 6*a1*sin(theta) + 8*b1*cos(theta) - 18*a3*sin(theta);
        double T23 = -4*b1*cos(theta) + 3*(a3 - a1)*sin(theta);

        double t2 = (-T22-sqrt(T22*T22-4*T21*T23))/(2*T21);
        double x,y;
        if( t2 >= 0){
            x = a1*(1-t2)*(1-t2)*(1-t2) + 3*a3*(1-t2)*(1-t2)*t2 - 3*a3*(1-t2)*t2*t2 - a1*t2*t2*t2;
            y = -4*b1*t2*(1-t2);

            return sin(theta)*x + cos(theta)*y;
        }else{
            double T11 = 18*a2*sin(theta) - 6*a1*sin(theta);
            double T12 = 6*a1*sin(theta) - 8*b1*cos(theta) - 18*a2*sin(theta);
            double T13 = 4*b1*cos(theta) + 3*(a2 - a1)*sin(theta);

            double t1 = (-T12-sqrt(T12*T12-4*T11*T13))/(2*T11);
            
            x = a1*(1-t1)*(1-t1)*(1-t1) + 3*a2*(1-t1)*(1-t1)*t1 - 3*a2*(1-t1)*t1*t1 - a1*t1*t1*t1;
            y = 4*b1*t1*(1-t1);

            return sin(theta)*x + cos(theta)*y;
        }
    }
}
