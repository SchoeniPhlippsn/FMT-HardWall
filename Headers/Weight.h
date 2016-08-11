
struct paras { double kr; double tt; double aa; double theta; double gauss; double mean; double diff; double RR; double AblRR; int l; int m; };

class weight { 
	public:
		fftw_complex *real;
		fftw_complex *fourier;
		std::string Filename;

		void setup (int,std::string);
		void kill_plan ();
		void print (int, double);
		void printReal (int, double);
		void fft (){fftw_execute(r2f);};
		void invfft (){fftw_execute(f2r);};
		void copy_real(const weight *,int);

	private:
		fftw_plan r2f;
		fftw_plan f2r;
};

void weight::setup (int n_bins, std::string name){
    
    real = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_bins);
    fourier = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_bins);
	r2f = fftw_plan_dft_1d(n_bins, real, fourier, -1, FFTW_ESTIMATE);
	f2r = fftw_plan_dft_1d(n_bins, fourier, real, 1, FFTW_ESTIMATE);
    Filename = name;
}

void weight::kill_plan (){
    fftw_destroy_plan(r2f);
    fftw_destroy_plan(f2r);
}

void weight::copy_real(const weight *w,int n_bins){
    for(int iz=0;iz<n_bins;iz++){
        real[iz]=w->real[iz];
	}
}

void weight::print (int n_bins, double Delta){
	FILE * reFile;
	FILE * ftFile;
    std::string Name;
	char file[30];
    Name = "Results/" + Filename;
	strcpy(file,Name.c_str());
	reFile = fopen (file,"w");
    
    Name = "Results/fft" + Filename;
	strcpy(file,Name.c_str());
	ftFile = fopen (file,"w");
	
    for( int iz =0; iz < n_bins; iz++){
        double z = (iz+0.5)*Delta;
        fprintf (reFile, "%.18f %.18f %.18f\n", z, creal(real[iz]), cimag(real[iz]));
        fprintf (ftFile, "%.18f %.18f %.18f\n", z, creal(fourier[iz]), cimag(fourier[iz]));
    }
}

void weight::printReal (int n_bins, double Delta){
	FILE * reFile;
	FILE * ftFile;
    std::string Name;
	char file[30];
    Name = "Results/" + Filename;
	strcpy(file,Name.c_str());
	reFile = fopen (file,"w");
    
    for( int iz =0; iz < n_bins; iz++){
        double z = iz*Delta;
        fprintf (reFile, "%.18f %.18f %.18f\n", z, creal(real[iz]), cimag(real[iz]));
    }
}
