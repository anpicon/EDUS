/*
 *  Constants.h
 *  
 *  Created by Antonio Picon on 5/31/11.
 *  Copyright 2011 JILA, University of Colorado. All rights reserved.
 *
 */

#include <iostream>
#include <cmath>
#include <complex>
using namespace std;

const double pi = 4.0*atan(1.0);
const double epsthd = 1.0e-8;
const double c  = 137.036;
const complex<double> c1(0.0,1.0);

double time_au_fs       =   0.02418884326505;
double time_fs_au       =   1./(0.02418884326505);
double energy_au_eV     =   27.211396;
double energy_eV_au     =   1.0/27.211396;
double energy_au_J      =   4.35974417e-18;
double energy_au_mJ     =   4.35974417e-15;
double energy_J_au      =   1./energy_au_J;
double energy_mJ_au     =   1./energy_au_mJ;
double energy_au_cm     =   219474.63;
double energy_cm_au     =   1./energy_au_cm;
double space_au_m       =   5.2917721092e-11;
double space_m_au       =   1./space_au_m;
double space_au_A       =   0.52917721092;
double space_A_au       =   1./space_au_A;
double XSections_au_Mb  =   28.0029;
double intensity_au_Wcm2=   3.50944758e+16; 
double intensity_Wcm2_au=   1./intensity_au_Wcm2;
double Efield_Vm_au     =   5.14221e11;
double Efield_au_Vm     =   1./Efield_Vm_au;