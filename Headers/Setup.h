
template<class T>
inline std::string toString(const T& t){
	std::ostringstream os;
	os << t;
  	return os.str();
}


void Setup(){
    rho.setup(n_bins,"rho.dat");
    rho_n.setup(n_bins,"rho_n.dat");
    w0.setup(n_bins,"w0.dat");
    w3.setup(n_bins,"w3.dat");
    n0.setup(n_bins,"n0.dat");
    //n12.setup(n_bins,"n12.dat");
    //n222.setup(n_bins,"n222.dat");
    n3.setup(n_bins,"n3.dat");
    PHI0.setup(n_bins,"PHI0.dat");
    //PHI12.setup(n_bins,"PHI12.dat");
    //PHI222.setup(n_bins,"PHI222.dat");
    PHI3.setup(n_bins,"PHI3.dat");
    sc0.setup(n_bins,"sc0.dat");
    //sc12.setup(n_bins,"sc12.dat");
    //sc222.setup(n_bins,"sc222.dat");
    sc3.setup(n_bins,"sc3.dat");
    c1.setup(n_bins,"c1.dat");

    for( i = 0; i <= lmax ; i++){
        std::string is = toString(i); 
        for( j = 0; j<=2*i; j++){
            std::string js = toString(j-i); 
            w1[i][j].setup(n_bins,"w1_" + is + "_" + js + ".dat");
            w2[i][j].setup(n_bins,"w2_" + is + "_" + js + ".dat");
            n1[i][j].setup(n_bins,"n1_" + is + "_" + js + ".dat");
            n2[i][j].setup(n_bins,"n2_" + is + "_" + js + ".dat");
            PHI1[i][j].setup(n_bins,"PHI1_" + is + "_" + js + ".dat");
            PHI2[i][j].setup(n_bins,"PHI2_" + is + "_" + js + ".dat");
            sc1[i][j].setup(n_bins,"sc1_" + is + "_" + js + ".dat");
            sc2[i][j].setup(n_bins,"sc2_" + is + "_" + js + ".dat");
        }
    }
    
    wreal = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_bins);
    wfourier = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_bins);
	wr2f = fftw_plan_dft_1d(n_bins, wreal, wfourier, -1, FFTW_ESTIMATE);
	wf2r = fftw_plan_dft_1d(n_bins, wfourier, wreal, 1, FFTW_ESTIMATE);
}
