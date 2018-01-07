//Constants

    const int lmax =8;
    const double kth=3.8;
    const double asp=3;
    const double a1=0.5;
    const double a2=a1*(kth-2.0/3.0*asp)/kth;
    const double a3=a1*(kth+2.0/3.0*asp)/kth;
    const double b1 = asp*a1;
    const double b1_inv = 1/b1;
    const double InvPi=1/M_PI;

//WeightDensities
    class weight w0[101];
    class weight w1[lmax+1][101];
    class weight w2[lmax+1][101];
	class weight w12[101];
	class weight w222[101];
	class weight w3[101];

	class weight rho[101];
	class weight rhon[101];
	class weight rho_fin;
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

	class weight P1;
	class weight P2;

//Universal Variables
	int n_bins, n_bins_2;
    double z, z2, kz, theta, atheta;
	double inv_n, H, H2, RII, RT;
    double DeltaR, DeltaK, DeltaT, DeltaP, inv_nDeltaR;
	int i, j, iz, v, st;
    int l,m, l1, l2, l3, m1, m2, m3;
    double prev, now, tt, aa;
    double V_pear;
     
//Varaibles Calc.h    
    int numb;
    double real, imag, costheta, sintheta;
    complex double coskz, sincosb1, cossinb1, exponent;
    double errorReal, errorImag, Wig;
    struct paras params = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    struct paras PARAMS = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double bessel[2*lmax+1];
      
//Varaiables PHI

    int factor;
    double on3, inv_on3, sclf;
    double alpha;
    double Coupling3;

//Varaiables rhon

    double b2, b3, rho_sum_o, rho_sum;

