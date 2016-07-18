
struct paras { double kr; double x_para1; double x_para2; double x_para3; double y_para; int l; int m; };

class weight { 
	public:
		fftw_complex *real;
		fftw_complex *fourier;
		std::string Filename;

		void setup (int,std::string);
		void kill_plan ();
		void print (int, double);
		void fft (){fftw_execute(r2f);};
		void invfft (){fftw_execute(f2r);};
		void copy_real(const weight *,int);

	private:
		fftw_plan r2f;
		fftw_plan f2r;
};

void weight::setup (int n_bins, std::string name){
    
    real = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_bins * n_bins * n_bins);
    fourier = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_bins * n_bins * n_bins);
	r2f = fftw_plan_dft_3d(n_bins, n_bins, n_bins, real, fourier, -1, FFTW_ESTIMATE);
	f2r = fftw_plan_dft_3d(n_bins, n_bins, n_bins, fourier, real, 1, FFTW_ESTIMATE);
    Filename = name;
}

void weight::kill_plan (){
    fftw_destroy_plan(r2f);
    fftw_destroy_plan(f2r);
}

void weight::copy_real(const weight *w,int n_bins){
	for(int ix=0;ix<n_bins;ix++){
		for(int iy=0;iy<n_bins;iy++){
            for(int iz=0;iz<n_bins;iz++){
                real[(ix*n_bins + iy)*n_bins+iz]=w->real[(ix*n_bins + iy)*n_bins+iz];
            }
		}
	}
}

void weight::print (int n_bins, double Delta){
	FILE * reFile;
	FILE * ftFile;
	FILE * TFile;
	FILE * File;
    std::string Name;
	char file[20];
	strcpy(file,Filename.c_str());
	reFile = fopen (file,"w");
    
    Name = "T" + Filename;
    strcpy(file,Name.c_str());
	File = fopen (file,"w");
	
    Filename = "fft" + Filename;
	strcpy(file,Filename.c_str());
	ftFile = fopen (file,"w");
    
    Name = "T" + Filename;
    strcpy(file,Name.c_str());
	TFile = fopen (file,"w");

	for( int ix =0; ix < n_bins; ix++){
		double x = ix*Delta;
        double xx;
        if(ix < n_bins/2) xx = Delta*(ix+n_bins/2);
        else xx = Delta*(ix-n_bins/2);

		for( int iy =0; iy < n_bins; iy++){
			double y = iy*Delta;
            double yy;
            if(iy < n_bins/2) yy = Delta*(iy+n_bins/2);
            else yy = Delta*(iy-n_bins/2);
            
            fprintf (File, "%.18f %.18f %.18f %.18f\n", x+0.5*Delta, y+0.5*Delta, creal(real[(ix*n_bins + iy)*n_bins + n_bins/2]), cimag(real[(ix*n_bins + iy)*n_bins + n_bins/2]));
            fprintf (TFile, "%.18f %.18f %.18f %.18f\n", x, y, creal(fourier[(ix*n_bins + iy)*n_bins + n_bins/2]), cimag(fourier[(ix*n_bins + iy)*n_bins + n_bins/2]));
            
            for( int iz =0; iz < n_bins; iz++){
                double z = iz*Delta;
                double zz;
                if(iz < n_bins/2) zz = Delta*(iz+n_bins/2);
                else zz = Delta*(iz-n_bins/2);

                fprintf (reFile, "%.18f %.18f %.18f %.18f %.18f\n", xx, yy, zz, creal(real[(ix*n_bins + iy)*n_bins+iz]), cimag(real[(ix*n_bins + iy)*n_bins+iz]));
                fprintf (ftFile, "%.18f %.18f %.18f %.18f %.18f\n", x, y, z, creal(fourier[(ix*n_bins + iy)*n_bins + iz]), cimag(fourier[(ix*n_bins + iy)*n_bins + iz]));
            }
        }
	}
}
