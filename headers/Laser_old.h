/*
 *  Laser.h
 *  
 *
 *  Created by Antonio Picon on 4/15/11.
 *  Copyright 2011 JILA, University of Colorado. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

//const double pi=4.*atan(1.); Define it in Constants.h
//const double c=137.036; Define it in Constants.h

class Laser{
public:
	Laser() : lambdanm(0.),wl(0.),Period(0.),Intensity(0.),A0(0.),E0(0.) {};
	double lambdanm,wl,Period,Intensity,A0,E0;
	void set_freq (double a) {wl=a; lambdanm=2.*pi*c*0.0529/wl; Period=2*pi/wl;          if(Intensity!=0){A0=sqrt(Intensity/(3.51e+16))/wl; E0=A0*wl;}}
	void set_lambdanm (double a){lambdanm=a; wl=2.*pi*c*0.0529/lambdanm; Period=2*pi/wl; if(Intensity!=0){A0=sqrt(Intensity/(3.51e+16))/wl; E0=A0*wl;}}
	void set_IntensityWcm2 (double a){Intensity=a; if(wl != 0) {A0=sqrt(Intensity/(3.51e+16))/wl; E0=A0*wl;}}
	void set_A0 (double a){A0=a; if(wl!=0){Intensity=A0*A0*(3.51e+16)*wl*wl; E0=A0*wl;}}
	void set_E0 (double a){E0=a; if(wl!=0){A0=E0/wl; Intensity=A0*A0*(3.51e+16)*wl*wl;}}
	void print_par(string);
	double get_Period(void){return Period;}
	double get_Intensity(void){return Intensity;}
	double get_E0(void){return E0;}
	double get_wl(void){return wl;}
	double PlaneWave (double,double,double,double);
	double PlaneWave_A (double,double,double,double);
	double PseudoPW (double,double,double,double,double,double);
	double PseudoPW_A (double,double,double,double,double,double);
	double PseudoPW_Envelope (double,double,double,double,double);
	double Sin2 (double,double,double,double);
	double Sin2_A (double,double,double,double);
	double Sin2_Envelope (double,double,double);
	double Gaussian(double,double,double,double);
	double Gaussian_A(double,double,double,double);
	double Gaussian_Envelope(double,double,double);
	double Gaussian_Chirp(double,double,double,double);
	double NumericalGaussian_A(double, double,double, double, double);
	double NumericalSin2_A(double, double , double, double, double);

	friend class Laser_FT;
};

void Laser::print_par(string title){
	cout << title << "\n";
	printf("*  wl(au)=                                         %10.2f                               *\n", wl); 
	printf("*  lambda(nm)=                                     %10.2f                               *\n", lambdanm);
	printf("*  Period(au)=                                     %10.2f                               *\n", Period);
	printf("*  A0(au)=                                         %10.2e                               *\n", A0);
	printf("*  E0(au)=                                         %10.2e                               *\n", E0);
	printf("*  Intensity(W/cm2)=                               %10.2e                               *\n", Intensity);
}
double Laser::PlaneWave(double tinitial, double tfinal,double time, double phase){
	if(time>=tinitial && time<=tfinal)return (E0*sin(wl*(time-tinitial)+phase));
	else return 0.;
}
double Laser::PlaneWave_A(double tinitial, double tfinal,double time, double phase){
	if(time>=tinitial && time<=tfinal)return (A0*sin(wl*(time-tinitial)+phase));
	else return 0.;
}
double Laser::PseudoPW(double tinitial,double t1,double t2,double tfinal,double time,double phase){
	if(time>=tinitial && time<t1)return (E0*(time-tinitial)*sin(wl*(time-tinitial)+phase)/(t1-tinitial));
	if(time>=t1 && time<t2)return (E0*sin(wl*(time-tinitial)+phase));
	if(time>=t2 && time<=tfinal)return (E0*(time-tfinal)*sin(wl*(time-tinitial)+phase)/(t2-tfinal));
	else return 0.;
}
double Laser::PseudoPW_A(double tinitial,double t1,double t2,double tfinal,double time,double phase){
	if(time>=tinitial && time<t1)return (A0*(time-tinitial)*sin(wl*(time-tinitial)+phase)/(t1-tinitial));
	if(time>=t1 && time<t2)return (A0*sin(wl*(time-tinitial)+phase));
	if(time>=t2 && time<=tfinal)return (A0*(time-tfinal)*sin(wl*(time-tinitial)+phase)/(t2-tfinal));
	else return 0.;
}
double Laser::PseudoPW_Envelope(double tinitial,double t1,double t2,double tfinal,double time){
	if(time>=tinitial && time<t1)return (E0*(time-tinitial)/(t1-tinitial));
	if(time>=t1 && time<t2)return (E0);
	if(time>=t2 && time<=tfinal)return (E0*(time-tfinal)/(t2-tfinal));
	else return 0.;
}
double Laser::Sin2(double tinitial, double tfinal,double time, double phase){
	if(time>=tinitial && time<=tfinal)return E0*sin(wl*(time-tinitial)+phase)*pow(sin(pi*(time-tinitial)/(tfinal-tinitial)),2);
	else return 0.;
}
double Laser::Sin2_A(double tinitial, double tfinal,double time, double phase){
	if(time>=tinitial && time<=tfinal)return A0*sin(wl*(time-tinitial)+phase)*pow(sin(pi*(time-tinitial)/(tfinal-tinitial)),2);
	else return 0.;
}
double Laser::Sin2_Envelope(double tinitial, double tfinal,double time){
	if(time>=tinitial && time<=tfinal)return E0*pow(sin(pi*(time-tinitial)/(tfinal-tinitial)),2);
	else return 0.;
}
double Laser::Gaussian(double b, double sigma,double time, double phase){
	if(pow(time-b,2)<40*sigma*sigma)return E0*sin(wl*(time-b)+phase)*exp(-pow(time-b,2)/(2*sigma*sigma));
	else return 0.;
}
double Laser::Gaussian_A(double b, double sigma,double time, double phase){
	if(pow(time-b,2)<40*sigma*sigma)return A0*cos(wl*(time-b)+phase)*exp(-pow(time-b,2)/(2*sigma*sigma));
	else return 0.;
}
double Laser::Gaussian_Envelope(double b, double sigma,double time){
	if(pow(time-b,2)<40*sigma*sigma)return E0*exp(-pow(time-b,2)/(2*sigma*sigma));
	else return 0.;
}
double Laser::Gaussian_Chirp(double b, double sigma,double time, double chirp){
	if(pow(time-b,2)<16.*sigma*sigma)return E0*sin((wl+chirp)*(time-b))*exp(-pow(time-b,2)/(2*sigma*sigma));
	else return 0.;
}



double Laser::NumericalGaussian_A(double b, double sigma,double time, double t1, double phase)
{
	int resolution = 10;
	double A=0.;
	double dt1 = (time - t1)/resolution;
	A -= Laser::Gaussian(b, sigma, t1, phase)/2.;
	for(int i=1; i<resolution; i++)
	{ 
    	A -= Laser::Gaussian(b, sigma, t1 + i*dt1, phase);
	}
	A = dt1*(A - Laser::Gaussian(b, sigma, time, phase))/2.;
	return (   A   );

}


double Laser::NumericalSin2_A(double tinitial, double tfinal, double time, double t1, double phase)
{
	int resolution = 10;
	double A=0.;
double dt1 = (time - t1)/double(resolution); 
	A -= Laser::Sin2(tinitial, tfinal,t1, phase)/2.;
	for(int i=1; i<resolution; i++)
	{ 
    	A -= Laser::Sin2(tinitial, tfinal, t1 + i*dt1, phase);
	}
	A = dt1*(A - Laser::Sin2(tinitial, tfinal, time, phase)/2.);
	return (   A   );

}
