//Constants

    const int lmax = 4;
    const double kth = 3.8;
    const double asp = 3;
    const double a1 = sqrt(2)*0.5;
    const double a2 = a1*(kth-asp)/kth;
    const double a3 = a1*(kth+asp)/kth;
    const double b1 = 1.5*sqrt(0.5);

//WeightDensities
    class weight w0;
    class weight w1[lmax+1][2*lmax+1];
    class weight w2[lmax+1][2*lmax+1];
	class weight w12;
	class weight w222;
	class weight w3;

	class weight rho;
	class weight rho_n;
	class weight n0;
    class weight n1[lmax+1][2*lmax+1];
    class weight n2[lmax+1][2*lmax+1];
	//class weight n12;
	//class weight n222;
    class weight n3;
	class weight PHI0;
    class weight PHI1[lmax+1][2*lmax+1];
    class weight PHI2[lmax+1][2*lmax+1];
	//class weight PHI12;
	//class weight PHI222;
	class weight PHI3;
	class weight sc0;
    class weight sc1[lmax+1][2*lmax+1];
    class weight sc2[lmax+1][2*lmax+1];
	//class weight sc12;
	//class weight sc222;
	class weight sc3;
	class weight c1;

    fftw_complex *wreal;
    fftw_complex *wfourier;
	fftw_plan wr2f;
	fftw_plan wf2r;

//Universal Variables
	int n_bins, n_bins_2;
    double x, y, z, kx, ky, kz;
	double inv_n, H, H2, RII, RT;
    double DeltaR, DeltaK, bulk;
	int i, j, ix, iy, iz, v, st;
    int l,m, l1, l2, l3, m1, m2, m3;
     
//Varaibles Calc.h    
    double real, imag, omega;
    double errorReal, errorImag;
    struct paras params = {0, a1, a2, a3, b1, 0, 0};
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_integration_qawo_table * table = gsl_integration_qawo_table_alloc(omega, b1, GSL_INTEG_COSINE, 1000);
    
    gsl_function F;

//Varaiables PHI

    int factor;
    double on3, inv_on3, rho_sum_o, rho_sum, sclf;
    double alpha;
    double Coupling3;
