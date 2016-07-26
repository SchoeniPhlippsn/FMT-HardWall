//Constants

    const int lmax = 8;
    const double kth = 3.8;
    const double asp = 3;
    const double a1 = sqrt(0.5)*0.5;
    const double a2 = a1*(kth-asp)/kth;
    const double a3 = a1*(kth+asp)/kth;
    const double b1 = 1.5*sqrt(0.5);

//WeightDensities
    class weight w0[101];
    class weight w1[lmax+1][101];
    class weight w2[lmax+1][101];
	class weight w12[101];
	class weight w222[101];
	class weight w3[101];

	class weight rho[101];
	class weight rhon[101];
	class weight n0;
    class weight n1[lmax+1];
    class weight n2[lmax+1];
    class weight n3;
	class weight PHI0;
    class weight PHI1[lmax+1];
    class weight PHI2[lmax+1];
	class weight PHI3;
	class weight sc0[101];
    class weight sc1[lmax+1][101];
    class weight sc2[lmax+1][101];
	class weight sc3[101];
	class weight c1[101];

//Universal Variables
	int n_bins, n_bins_2;
    double z, kz, theta;
	double inv_n, H, H2, RII, RT;
    double DeltaR, DeltaK, DeltaT, DeltaP, bulk, inv_nDeltaR;
	int i, j, iz, v, st;
    int l,m, l1, l2, l3, m1, m2, m3;
    double prev, now;
     
//Varaibles Calc.h    
    double real, imag, costheta, sintheta;
    double errorReal, errorImag;
    struct paras params = {0, a1, a2, a3, b1, 0, 0};
    
//Varaiables PHI

    int factor;
    double on3, inv_on3, sclf;
    double alpha;
    double Coupling3;

//Varaiables rhon

    double b2, b3, rho_sum_o[101], rho_sum;

