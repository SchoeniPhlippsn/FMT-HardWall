
template<class T>
inline std::string toString(const T& t){
	std::ostringstream os;
	os << t;
  	return os.str();
}


void Setup(){

    n0.setup(n_bins,"n0.dat");
    n3.setup(n_bins,"n3.dat");
    PHI0.setup(n_bins,"PHI0.dat");
    PHI3.setup(n_bins,"PHI3.dat");
    rho_fin.setup(n_bins,"rho_fin.dat");
    
    for( i = 0; i <= lmax ; i++){
        std::string is = toString(i); 
        n1[i].setup(n_bins,"n1-" + is + ".dat");
        n2[i].setup(n_bins,"n2-" + is + ".dat");
        PHI1[i].setup(n_bins,"PHI1-" + is + ".dat");
        PHI2[i].setup(n_bins,"PHI2-" + is + ".dat");
    }

    for( j=0; j <= 100; j++){
        std::string js = toString(j); 

        w0[j].setup(n_bins,"w0_" + js + ".dat");
        w3[j].setup(n_bins,"w3_" + js + ".dat");
        rho[j].setup(n_bins,"rho_" + js + ".dat");
        rhon[j].setup(n_bins,"rhon_" + js + ".dat");
        sc0[j].setup(n_bins,"sc0_" + js + ".dat");
        sc3[j].setup(n_bins,"sc3_" + js + ".dat");
        c1[j].setup(n_bins,"c1_" + js + ".dat");

        for( i = 0; i <= lmax ; i++){
            std::string is = toString(i); 
            w1[i][j].setup(n_bins,"w1-" + is + "_" + js + ".dat");
            w2[i][j].setup(n_bins,"w2-" + is + "_" + js + ".dat");
            sc1[i][j].setup(n_bins,"sc1-" + is + "_" + js + ".dat");
            sc2[i][j].setup(n_bins,"sc2-" + is + "_" + js + ".dat");
        }
    }
}
