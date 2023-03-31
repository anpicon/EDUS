class Laser//use methods E(t) and A(t) to calculate the pulse at a particular time
{
    public:
        double lambdanm;   //wavelength in nm
        double wl;          //frequency in au
        double Period;      //Period in au
        double Intensity;   //intensity in au
        double A0;          //A0 in au
        double E0;          //E0 in au
        double t0;          //CARE! for sin2-> initial time, for gaussian -> peak of the pulse 
        double tf;          //useful for sin2. tf = t0+ncycles*period
        double w1, w2;      //the window of pseudoPW. w1-> start time; w2-> end time
        double ncycle;      //in gaussian pulses is sigma, in sin2 is the number of cycles
        double sigma;
        bool envelope;      //true if you dont want the oscillations in the envelope
        bool gaussian;      //true = gaussian
        bool PW;            //true = pw, false = sin2
        bool sin2;          //true = sin2
        double phase;
        Coord_B pol1, pol2; //vectors spanning the plane of polarization
        double t_track;     // time to integrate the electric field and get the vector potential
        multivec1D<double>  A_track_cart, A_track_crys;// vector potential at t_track
        double resolution;  //resolution in au to calculate the vector potential via integral
        //constructor        
        Laser() : phase(0.),lambdanm(0.),wl(0.),Period(0.),Intensity(0.),A0(0.),E0(0.),t0(0.),tf(0.),ncycle(0.),sigma(0.),t_track(-100),resolution(0.001), w1(0.), w2(0.)
        {A_track_cart.resize(3);A_track_crys.resize(3);};

        //methods to set variables
        void set_freq (double& a) {wl=a; lambdanm=2.*pi*c*0.0529/wl; Period=2*pi/wl; if(Intensity!=0){A0=sqrt(Intensity/(3.51e+16))/wl; E0=A0*wl;} if(ncycle != 0) tf=t0+ncycle*Period; }
        void set_lambdanm (double& a){lambdanm=a; wl=2.*pi*c*0.0529/lambdanm; Period=2*pi/wl; if(Intensity!=0){A0=sqrt(Intensity/(3.51e+16))/wl; E0=A0*wl;} if(ncycle != 0) tf=t0+ncycle*Period; }
        void set_IntensityWcm2 (double& a){Intensity=a; if(wl != 0) {A0=sqrt(Intensity/(3.51e+16))/wl; E0=A0*wl;}}
        void set_ncycle_sigma (double& a){ncycle=a; sigma=a; if(Period != 0) tf=t0+ncycle*Period; } 
        void set_t0 (double& a) {t0=a; if(tf!=0.) tf+=t0; }
        void set_boolean(bool& env, bool& gaus, bool& pw, bool& sin){gaussian = gaus; envelope = env; PW = pw; sin2 = sin;}
        void set_pol(double&,double&,double&,double&,double&,double&);
        void set_window(double& window1, double& window2){w1=window1; w2=window2;}
        void set_phase(double& phi){phase=phi;}
        //methods to recap variables
        void print_par(string&);

        //functions for laser
        double gauss(double&);
        double Sin2(double&);
		double trapezoid(double& t);

        //methods to get electric field and vector potential
        vec1d E(double& t);
        vec1d& A(double& t);
        vec1d E_crys(double& t);
        vec1d& A_crys(double& t);
};




double Laser::gauss(double& t)
{
    if(pow(t-t0,2)<40*sigma*sigma)return exp(-pow(t-t0,2)/(2*sigma*sigma));
    else return 0.;
}


double Laser::Sin2(double& t)
{
    if(t>=t0 && t<=tf)return pow(sin(pi*(t-t0)/(tf-t0)),2);
    else return 0.;
}


double Laser::trapezoid(double& t)
{
	if(t>=t0 && t<w1)return (t-t0)/(w1-t0);
	if(t>=w1 && t<w2)return 1.;
	if(t>=w2 && t<=tf)return (t-tf)/(w2-tf);
	else return 0.;
}





void Laser::set_pol(double& a0, double& a1, double& a2, double& b0, double& b1, double& b2)
{
    if( (b0 != 0. && b1 != 0. && b2 != 0.) && ( abs(a0/(b0+1e-11) - a1/(b1+1e-11)) <1e-9 && abs(a1/(b1+1e-11) - a2/(b2+1e-11)) <1e-9))
    {
        printf("care! you defined a plane with two parallel vectors!");
        exit(1);
    }
    else if( abs(a0*b0 + a1*b1 + a2*b2) >= 1e-7 )//gramd-schmidt orthogonalization
    {
        double fac_proj_pol1 = (a0*b0 + a1*b1 + a2*b2)/(a0*a0+a1*a1+a2*a2); 
        b0 -= fac_proj_pol1*a0;
        b1 -= fac_proj_pol1*a1;
        b2 -= fac_proj_pol1*a2;
        double value = sqrt(b0*b0+b1*b1+b2*b2);
        b0 /= value; b1/= value; b2/=value;

        if(abs(a0*b0 + a1*b1 + a2*b2) >= 1e-7 )
        {
            printf("WARNING! unsuccessful orthogonalization!! care!!\n");
        }
    }
    pol1.setcart(a0,a1,a2);
    pol2.setcart(b0,b1,b2);

}//new





void Laser::print_par(string& title)
{
    cout << title << "\n";
    printf("*  Energy                                %11.5f au   %11.5f eV                    *\n", wl, wl*energy_au_eV); 
    printf("*  Wavelength                            %11.5f au   %11.5f nm                    *\n", 10.*lambdanm*space_A_au,lambdanm);
    printf("*  Period                                %11.5f au   %11.5f fs                    *\n", Period, Period*time_au_fs);
    printf("*  A0                                    %11.5e au   %11.5e V s/m                 *\n", A0, A0);
    printf("*  E0                                    %11.5e au   %11.5e V/ang                 *\n", E0, E0/(Efield_au_Vm*1.e10));
    printf("*  Intensity                             %11.5e au   %11.5e W/cm^2                *\n", Intensity*intensity_Wcm2_au,Intensity);
    printf("*  t0=                                   %11.5e au   %11.5e fs                    *\n", t0, t0*time_au_fs);
    printf("*  tf=                                   %11.5e au   %11.5e fs                    *\n", tf, tf*time_au_fs);
    printf("*  polarization vectors ->         (%2.4f,%2.4f,%2.4f) (%2.4f,%2.4f,%2.4f)            *\n", pol1.cart[0],  pol1.cart[1],  pol1.cart[2], pol2.cart[0],  pol2.cart[1], pol2.cart[2]);
}



multivec1D<double> Laser::E(double& t)//return electric field in cartesian coordinates
{
    double E_envelope;
    //choose envelope
    if(gaussian)   E_envelope = E0*Laser::gauss(t);
    else if(sin2)  E_envelope = E0*Laser::Sin2(t);
    else if(PW)    E_envelope = E0*Laser::trapezoid(t);
    else E_envelope = 0.;	
    double plane_1, plane_2;
    //choose varying part
    plane_1 = (envelope ? 1. : sin(wl*(t-t0)+phase));
    plane_2 = (envelope ? 1. : cos(wl*(t-t0)+phase));
    vec1d E_cart(3);
    for(int i=0; i<3; i++)
         E_cart[i] = E_envelope*(plane_1*pol1.cart[i]+plane_2*pol2.cart[i]);
    return E_cart;
}


multivec1D<double> Laser::E_crys(double& t)//return electric field in cartesian coordinates
{
    double E_envelope;
    //choose envelope
    if(gaussian)   E_envelope = E0*Laser::gauss(t);
    else if(sin2)  E_envelope = E0*Laser::Sin2(t);
    else if(PW)    E_envelope = E0*Laser::trapezoid(t);
    else E_envelope = 0.;	
    
    double plane_1, plane_2;
    //choose varying part
    plane_1 = (envelope ? 1. : sin(wl*(t-t0)+phase));
    plane_2 = (envelope ? 1. : cos(wl*(t-t0)+phase));
    vec1d E_crys(3);
    for(int i=0; i<3; i++)
         E_crys[i] = E_envelope*(plane_1*pol1.crys[i]+plane_2*pol2.crys[i]);
    return E_crys;
}

multivec1D <double>& Laser::A(double& t)//return vector potential via numerical integration in cartesian coordinates
{
    if(t < t_track) 
    {
	   t_track = t-100;
       A_track_cart[0] = 0.;
       A_track_cart[1] = 0.;
       A_track_cart[2] = 0.;
    }
    int nstep = round((t-t_track)/resolution);
    if(nstep < 1) nstep = 1;
    double res_double = (t-t_track)/double(nstep);

    A_track_cart -= res_double*this->E(t_track)*0.5;
    for(int istep=1; istep<nstep; istep++)
    { 
        double tt = t_track + istep*res_double; 
        A_track_cart -= res_double*this->E(tt);
    }
    A_track_cart = A_track_cart - this->E(t)*0.5*res_double;
    t_track = t;
    return A_track_cart;
}


multivec1D <double>& Laser::A_crys(double& t)//return vector potential via numerical integration in cartesian coordinates
{
    if(t < t_track) 
    {
 	    t_track = t-100;
    	A_track_crys[0] = 0.;
        A_track_crys[1] = 0.;
        A_track_crys[2] = 0.;
    }
    int nstep = round((t-t_track)/resolution);
    if(nstep < 1) nstep = 1;
    double res_double = (t-t_track)/double(nstep);

    A_track_crys -= res_double*this->E_crys(t_track)*0.5;
    for(int istep=1; istep<nstep; istep++)
    { 
        double tt = t_track + istep*res_double; 
        A_track_crys -= res_double*this->E_crys(tt);
    }
    A_track_crys = A_track_crys - this->E_crys(t)*0.5*res_double;
    t_track = t;
    return A_track_crys;
}
