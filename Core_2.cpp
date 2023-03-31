#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <complex>
#include <ctime>
#include <utility>
#include <memory>
#include <armadillo>

using namespace std;
using namespace arma;

#include "BoostArrays.h"
#include "typedef.h"
#include "Constants.h"
#include "Interpolation.h"
#include "Coordinate.h"
#include "crystal.h"
#include "Laser.h"
#include "ReadInput.h"
#include "stdarg.h"

/*********************************************************************************************************************************
******************************                   Class for core orbitals                  ****************************************
*********************************************************************************************************************************/
class Core
{
    public:
        Core(vec1d& Origin, vec1d& at_coord, string& type);
        Core(){};
        void init(vec1d& Origin, vec1d& at_coord, string& type);
        double operator()(double x, double y, double z);
        double energy(void);   

    private:
        int _N;
        vec1d _Origin;
        vec1d _at_coord;
        vec2d _R;
        vec1d _n; //= {1, 1, 3, 2, 2, 2, 2} ;
        vec1d _Cij; //= {0.352872, 0.473621, -0.001199, 0.210887, 0.000886, 0.000465, -0.000119} ;
        vec1d _Z; //= {8.4936, 4.8788, 15.4660, 7.0500, 2.2640, 1.4747, 1.1639} ;
        string _type;
};

Core::Core(vec1d& Origin, vec1d& at_coord, string& type):  _type(type)
{
    _Origin.resize(3);
    for(int i=0; i<3; i++) _Origin[i]= Origin[i];

    _at_coord.resize(3);
    for(int i=0; i<3; i++) _at_coord[i]= at_coord[i];
    

    if(_type == "C1s")
    {
        _N = 7;
        _n.resize(7);
        _Cij.resize(7);
        _Z.resize(7);

        _Z[0]   = 8.4936   ;   _Z[1]   = 4.8788  ;     _Z[2]   = 15.4660  ;   _Z[3]   = 7.0500  ;    _Z[4]   = 2.2640  ;      _Z[5]   = 1.4747   ;    _Z[6]   = 1.1639;  
        _Cij[0] = 0.352872 ;   _Cij[1] = 0.473621;     _Cij[2] = -0.001199;   _Cij[3] = 0.210887;    _Cij[4] = 0.000886;      _Cij[5] = 0.000465 ;    _Cij[6] = -0.000119 ;  
        _n[0]   = 1.        ;  _n[1]   = 1.       ;    _n[2]   = 3.       ;   _n[3]   = 2.      ;    _n[4]   = 2.      ;      _n[5]   = 2.       ;    _n[6]   = 2. ; 
    }
    else if (_type == "C2p")
    {
        _N = 5;
        _n.resize(5);
        _Cij.resize(5);
        _Z.resize(5);

        _Z[0]   = 7.0500   ;   _Z[1]   = 3.2275  ;     _Z[2]   = 2.1908   ;   _Z[3]   = 1.4413  ;    _Z[4]   = 1.0242  ;
        _Cij[0] = 0.006977 ;   _Cij[1] = 0.070877;     _Cij[2] = 0.230802 ;   _Cij[3] = 0.411931;    _Cij[4] = 0.350701;
        _n[0]   = 2.       ;   _n[1]   = 2.       ;    _n[2]   = 2.       ;   _n[3]   = 2.      ;    _n[4]   = 2.      ;
    }
    else
    {
        printf("Error in the Core constructor. Possible choices for the orbitals: C1s, C2p.");
    }

}


void Core::init(vec1d& Origin, vec1d& at_coord, string& type)
{
    _type = type;
    _Origin.resize(3);
    for(int i=0; i<3; i++) _Origin[i]= Origin[i];

    _at_coord.resize(3);
    for(int i=0; i<3; i++) _at_coord[i]= at_coord[i];
    

    if(_type == "C1s")
    {
        _N = 7;
        _n.resize(7);
        _Cij.resize(7);
        _Z.resize(7);

        _Z[0]   = 8.4936   ;   _Z[1]   = 4.8788  ;     _Z[2]   = 15.4660  ;   _Z[3]   = 7.0500  ;    _Z[4]   = 2.2640  ;      _Z[5]   = 1.4747   ;    _Z[6]   = 1.1639;  
        _Cij[0] = 0.352872 ;   _Cij[1] = 0.473621;     _Cij[2] = -0.001199;   _Cij[3] = 0.210887;    _Cij[4] = 0.000886;      _Cij[5] = 0.000465 ;    _Cij[6] = -0.000119 ;  
        _n[0]   = 1.        ;  _n[1]   = 1.       ;    _n[2]   = 3.       ;   _n[3]   = 2.      ;    _n[4]   = 2.      ;      _n[5]   = 2.       ;    _n[6]   = 2. ; 
    }
    else if (_type == "C2p")
    {
        _N = 5;
        _n.resize(5);
        _Cij.resize(5);
        _Z.resize(5);

        _Z[0]   = 7.0500   ;   _Z[1]   = 3.2275  ;     _Z[2]   = 2.1908   ;   _Z[3]   = 1.4413  ;    _Z[4]   = 1.0242  ;
        _Cij[0] = 0.006977 ;   _Cij[1] = 0.070877;     _Cij[2] = 0.230802 ;   _Cij[3] = 0.411931;    _Cij[4] = 0.350701;
        _n[0]   = 2.       ;   _n[1]   = 2.       ;    _n[2]   = 2.       ;   _n[3]   = 2.      ;    _n[4]   = 2.      ;
    }
    else
    {
        printf("Error in the Core constructor. Possible choices for the orbitals: C1s, C2p.");
    }

}



double Core::energy()
{
    if(_type == "C1s")
        return(-11.325519);
    else if(_type == "C2p")
        return(-0.433341);
    else return 0.;
}

//Core::~Core()   {    }

double Core::operator()(double x, double y, double z) 
{       
    double xtot = x - _at_coord[0] + _Origin[0];    //here everything passed to the class is in a.u.
    double ytot = y - _at_coord[1] + _Origin[1];    // the final + sign is for change of reference system
    double ztot = z - _at_coord[2] + _Origin[2];    // if we have O and O' and x_O' are the Coord_Rs of O' with respect to O, then r=x_O'+r'
    double r = sqrt(xtot*xtot + ytot*ytot + ztot*ztot); // r in a.u.                   //this way we center the atom in the origin O.
    double r2 = r*r;
    double value=0.;

    for(int ii=0; ii< _N; ii++)
    {   
        double factorial=1.;
        double nfact=2*_n[ii];
        while (nfact > 1) //calculating denominator of NN
        {
            factorial *= nfact;
            nfact -= 1;
        }
        factorial = pow(factorial, 0.5);
        nfact = Core::_Cij[ii]*pow( 2.*Core::_Z[ii], Core::_n[ii]+0.5 )/factorial; //part of R(r) that does not depend on r
    
        value += nfact * pow(r, Core::_n[ii]-1) * exp ( -Core::_Z[ii]*r );              
        nfact = 0.; 
    }


    //adding angular part 
    if(_type == "C1s")  value /= 2.*sqrt(pi);
    else if (_type == "C2p")
    {
        if(r!=0)        
        {
            double Yx = .5*sqrt(3./pi)*xtot/r;
            value *= Yx;
        }
    } 
    return value;
    
}



/*********************************************************************************************************************************
******************************                   Functions of the program                 ****************************************
*********************************************************************************************************************************/
void test_integral();
void read_xsf(string&);
void set_nat_ntype();
void set_Npoints();
void set_unitcell();
void set_R();
void set_origin();
void set_r_at();
void setvariables();
void print_info();
void setwf();
void interpolate();
void Plotxsf_X0_Y_Zmiddle(int& );
void PlotFitted_X0_Y_Zmiddle(FittedData<double>& PsiValence, int& at);
void Plotxsf_X_Y0_Zmiddle(int& at);
void PlotFitted_X_Y0_Zmiddle(FittedData<double>& PsiValence, int& at);
double Integrate3D(vec3d&, double&, vec1d&);
void Normalize();
/*********************************************************************************************************************************
******************************                   Variables of the program                 ****************************************
*********************************************************************************************************************************/

string seedname = "graphene";                      //name of the files from w90
int nat;                                           //number of atoms in the unit cell
int ntype;                                         //number of different types of atoms 
multivec1D<vec3d> Psi;                                         //wannier functions. first indexx-> which wannier from second to last(static)->R3, R2, R1
vector<string> lines;                              //dynamic array of strings to save the content of the xsf files
stringstream sname;                                //to add numbers and variables to the string
Coord_R origin;                                    //origin of axis in xsf file 
vector<Coord_R> r_at;                              //coordinates for each atom in the system of reference
vec1i Npoints(3);                                  //points in the grid of the xsf files for each direction
vec2d R(3,3);                                      //supercell vectors -> to use to construct Coord_R
vector<Coord_R> a;                              //unit cell vectors
vector<string> line;                               //vector with substrings of the line
string file;    
size_t indexx;                                     //dummy variable to fix lines when needed
size_t pos;                                        //is the final position in the string. Needed for function Separatestring

vec1d spacing(3);                                  //spacing in crystal coordinates to interpolate

int density = 10;                                  //number of points between two points after interpolation
multivec1D<FittedData<double>> Psi_all(nat);
vec1i resolution(3);

Core PsiC1s;




int main()
{   


    //recover the name of the wannier files with the indexx and the extension
    string name = seedname+"_00001.xsf";
    read_xsf(name);    //read the xsf file 
    printf("read information from %s file.\n", name.c_str());
    setvariables();
    printf("set global variables from the file.\n");
    print_info();
    setwf();
    printf("Wavefunctions read in all xsf.\n");
    Normalize();
    interpolate();
    printf("Interpolation completed.\n"); 
    ofstream fp_out;
    fp_out.open("dipole.txt");

    resolution[0]=300; resolution[1]=300; resolution[2]=300;
    vec3d PsirPsi1(resolution[0], resolution[1], resolution[2]);
    vec3d PsirPsi2(resolution[0], resolution[1], resolution[2]);
    vec3d PsirPsi3(resolution[0], resolution[1], resolution[2]);

    vec1d finer_spacing(3);
    for(int i=0; i<3; i++) finer_spacing[i] = 1./double(resolution[i]-1);
        
    for(int iR1=-2; iR1<=2; iR1++)
    {
        for(int iR2=-2; iR2<=2; iR2++)
        {
            for(int iR3=0; iR3<=0; iR3++)
            {
                fp_out << iR1 << " " << iR2 << " " << iR3 << endl;
                for(int iat=0; iat<nat; iat++)
                {
                    for(int iatR=0; iatR<nat; iatR++)
                    {
                        cout << "iR1 " << iR1 <<"  iR2 " << iR2 <<"  iR3 " << iR3 << endl;
                        vec1d O(3); O.fill(0);
                        string typeat = "C1s";
                        PsiC1s.init(O, r_at[iat].cart, typeat);
                        

                        for(int ix=0; ix<resolution[0]; ix++)
                        {
                            for(int iy=0; iy<resolution[1]; iy++)
                            {
                                for(int iz=0; iz<resolution[2]; iz++)
                                {
                                    Coord_R r, r_shift; 
                                    r.setcrys(double(ix)*finer_spacing[0]+ origin.crys[0], double(iy)*finer_spacing[1]+ origin.crys[1], double(iz)*finer_spacing[2]+origin.crys[2]);
                                    //r_core.setcrys(double(ix)*finer_spacing[0], double(iy)*finer_spacing[1], double(iz)*finer_spacing[2]);
                                    vec1d shift_crys(3); 
                                    for(int i=0; i<3; i++)
                                        shift_crys[i] = iR1*a[0].crys[i] + iR2*a[1].crys[i] + iR3*a[2].crys[i];
                                    
                                    //PsirPsi1[ix][iy][iz] = r.crys[0]*Psi_all[iat](r.crys[0], r.crys[1], r.crys[2])*Psi_all[iatR](r.crys[0]+shift_crys[0], r.crys[1]+shift_crys[1], r.crys[2]+shift_crys[2]);
                                    //PsirPsi2[ix][iy][iz] = r.crys[1]*Psi_all[iat](r.crys[0], r.crys[1], r.crys[2])*Psi_all[iatR](r.crys[0]+shift_crys[0], r.crys[1]+shift_crys[1], r.crys[2]+shift_crys[2]);
                                    //PsirPsi3[ix][iy][iz] = r.crys[2]*Psi_all[iat](r.crys[0], r.crys[1], r.crys[2])*Psi_all[iatR](r.crys[0]+shift_crys[0], r.crys[1]+shift_crys[1], r.crys[2]+shift_crys[2]);
                                    PsirPsi1[ix][iy][iz] = r.crys[0]*PsiC1s(r.cart[0], r.cart[1], r.cart[2])*Psi_all[iatR](r.crys[0]-shift_crys[0], r.crys[1]-shift_crys[1], r.crys[2]-shift_crys[2]);
                                    PsirPsi2[ix][iy][iz] = r.crys[1]*PsiC1s(r.cart[0], r.cart[1], r.cart[2])*Psi_all[iatR](r.crys[0]-shift_crys[0], r.crys[1]-shift_crys[1], r.crys[2]-shift_crys[2]);
                                    PsirPsi3[ix][iy][iz] = r.crys[2]*PsiC1s(r.cart[0], r.cart[1], r.cart[2])*Psi_all[iatR](r.crys[0]-shift_crys[0], r.crys[1]-shift_crys[1], r.crys[2]-shift_crys[2]);

                                }
                            }
                        }//end saving functions to integrate

                        double Jacobian = Coord_R::getJ();
                        double r0 = Integrate3D(PsirPsi1,Jacobian,finer_spacing);
                        double r1 = Integrate3D(PsirPsi2,Jacobian,finer_spacing);
                        double r2 = Integrate3D(PsirPsi3,Jacobian,finer_spacing);
                        Coord_R rtot;
                        rtot.setcrys(r0,r1,r2);
                        fp_out << iat << " " << iatR << " " << rtot.cart[0]*space_au_A << " " << rtot.cart[1]*space_au_A << " " << rtot.cart[2]*space_au_A << endl; 
                    }//end iatR
                }//end iat
                fp_out << endl;
            }//end iR3
        }//end iR2
    }//end iR1



}






//function to read xsf files from wannier90
//arguments: filename -> string with the name of the file to open
//           line     -> dynamic array of strings to save the content of the file
void read_xsf(string& file)
{
    ifstream fp_input;
    string counter;
    fp_input.open(file);
    int i=0;
    lines.push_back("");
    while(getline(fp_input,lines[i]))   {lines.push_back(""); i++;}  
    fp_input.close();
    sname.str(std::string());  //clear the sname variable to use it again
}



void set_nat_ntype()
{
    indexx = 14; //nat and ntype are always at line #15 
    pos = lines[indexx].length();
    Separate_string(lines[indexx], line, pos);
    nat   = atoi(line[0].c_str()); 
    ntype = atoi(line[1].c_str()); 
    line.clear();  //erase memory
}

void set_Npoints()
{
    indexx = 22; //nat and ntype are always at line #23 
    pos = lines[indexx].length();
    Separate_string(lines[indexx], line, pos);
    for(int i=0; i<3; i++) Npoints[i] = atoi(line[i].c_str());
    for(int i=0; i<3; i++) spacing[i]=1./double(Npoints[i]-1); 
    //for(int i=0; i<3; i++) cout << Npoints[i] << endl;
    line.clear();  //erase memory
}
void set_R()
{
    indexx = 24;                        //the required info start from line #25-28
    for(int iR=0; iR<3; iR++)          //through rows of R
    {
        size_t pos = lines[indexx+iR].length();
        Separate_string(lines[indexx+iR], line, pos);//row

        for(int ix=0; ix<3; ix++)  //through columns
        {//cout << line[ix] << endl;
            R[iR][ix] = atof(line[ix].c_str())*space_A_au;
        }
        line.clear();
    }    
    Coord_R::set_crys_to_cart(R);//the transpose is in this function

}
void set_origin()
{
    indexx = 23; //nat and ntype are always at line #15 
    pos = lines[indexx].length();
    Separate_string(lines[indexx], line, pos);
    origin.setcart( (atof(line[0].c_str())*space_A_au),  (atof(line[1].c_str())*space_A_au),  (atof(line[2].c_str())*space_A_au) );  
    //cout << origin.cart[0]*space_au_A << " "<<origin.cart[1]*space_au_A <<" "<<origin.cart[2]*space_au_A <<endl;

    line.clear();  //erase memory

}

void set_r_at()
{
    indexx = 15;                        //the required info start from line #16
    r_at.resize(nat);                  //resize r_at with the number of atoms
    
    for(int iat=0; iat<r_at.size(); iat++)
    {
        size_t pos = lines[indexx+iat].length();
        Separate_string(lines[indexx+iat], line, pos);

        r_at[iat].setcart( (atof(line[1].c_str())*space_A_au),  (atof(line[2].c_str())*space_A_au),  (atof(line[3].c_str())*space_A_au) );  
        //cout << r_at[iat].cart[0]*space_au_A << " "<<r_at[iat].cart[1]*space_au_A <<" "<<r_at[iat].cart[2]*space_au_A <<endl;
        line.clear();
    }

}


void setvariables()
{
    set_nat_ntype();
    set_Npoints();
    set_R();
    set_unitcell();
    set_origin();
    set_r_at();
}


void print_info()
{
      printf("Calculation ready to start.\n\n");
          printf("\n*********************************************************************************************\n");
          printf(  "*                              %20s                                         *\n",seedname.c_str());
          printf("\n*********************************************************************************************\n");
          printf(  "*                              Summary of parameters                                        *\n");
          printf(  "*********************************************************************************************\n");  
          printf(  "*  nat       %20s        %3d        *\n", " ", nat);
          printf(  "*  Npoints   %20s   %3d  %3d   %3d  *\n", " ", Npoints[0],Npoints[1], Npoints[2]);
          printf(  "*  origin    %20s   (%8.2f, %8.2f, %8.2f) angstrom     (%8.2f, %8.2f, %8.2f) au  *\n", " ", origin.cart[0]*space_au_A, origin.cart[1]*space_au_A, origin.cart[2]*space_au_A, origin.cart[0], origin.cart[1], origin.cart[2]);
          printf(  "*  a1        %20s   (%8.2f, %8.2f, %8.2f) angstrom     (%8.2f, %8.2f, %8.2f) au  *\n", " ", a[0].cart[0]*space_au_A, a[0].cart[1]*space_au_A, a[0].cart[2]*space_au_A, a[0].cart[0], a[0].cart[1], a[0].cart[2]);
          printf(  "*  a2        %20s   (%8.2f, %8.2f, %8.2f) angstrom     (%8.2f, %8.2f, %8.2f) au  *\n", " ", a[1].cart[0]*space_au_A, a[1].cart[1]*space_au_A, a[1].cart[2]*space_au_A, a[1].cart[0], a[1].cart[1], a[1].cart[2]);
          printf(  "*  a3        %20s   (%8.2f, %8.2f, %8.2f) angstrom     (%8.2f, %8.2f, %8.2f) au  *\n", " ", a[2].cart[0]*space_au_A, a[2].cart[1]*space_au_A, a[2].cart[2]*space_au_A, a[2].cart[0], a[2].cart[1], a[2].cart[2]);
          printf(  "*  R1        %20s   (%8.2f, %8.2f, %8.2f) angstrom     (%8.2f, %8.2f, %8.2f) au  *\n", " ", R[0][0]*space_au_A, R[0][1]*space_au_A, R[0][2]*space_au_A, R[0][0], R[0][1], R[0][2]);
          printf(  "*  R2        %20s   (%8.2f, %8.2f, %8.2f) angstrom     (%8.2f, %8.2f, %8.2f) au  *\n", " ", R[1][0]*space_au_A, R[1][1]*space_au_A, R[1][2]*space_au_A, R[1][0], R[1][1], R[1][2]);
          printf(  "*  R3        %20s   (%8.2f, %8.2f, %8.2f) angstrom     (%8.2f, %8.2f, %8.2f) au  *\n", " ", R[2][0]*space_au_A, R[2][1]*space_au_A, R[2][2]*space_au_A, R[2][0], R[2][1], R[2][2]);
          printf(  "*  Jacobian  %20s         %8.2f A^3                            %8.2f au^3        *\n", " ",Coord_R::getJ()*space_au_A*space_au_A*space_au_A, Coord_R::getJ());
          for(int i=0; i<nat; i++)
          printf(  "*  r_at_%2d  %20s   (%8.2f, %8.2f, %8.2f) angstrom     (%8.2f, %8.2f, %8.2f) au  *\n", i, " ", r_at[i].cart[0]*space_au_A,r_at[i].cart[1]*space_au_A, r_at[i].cart[2]*space_au_A, r_at[i].cart[0], r_at[i].cart[1], r_at[i].cart[2] );
          printf(  "*********************************************************************************************\n");  
}




void set_unitcell()
{
    line.clear();
    indexx = 6;                        //the required info start from line #7
    a.resize(3);
    for(int i=0; i<3; i++)
    {
        size_t pos = lines[indexx+i].length();
        Separate_string(lines[indexx+i], line, pos);
        a[i].setcart( (atof(line[0].c_str())*space_A_au),  (atof(line[1].c_str())*space_A_au),  (atof(line[2].c_str())*space_A_au) );  
        line.clear();
    } 
}

void setwf()
{
    lines.clear();
    Psi.resize(nat);

    for(int i=0; i<nat; i++)
    {
        Psi[i].resize(Npoints[0], Npoints[1], Npoints[2]);

        sname.seekp(0,ios::beg);
        sname << seedname << "_0000" << i+1 << ".xsf";
        file = sname.str(); //cout << sname.str() << endl;//to avoid error in calling read_xsf
        read_xsf(file);    //read the xsf file 
        
        string wavefunction;
        for(indexx = 27; indexx<lines.size(); indexx++)
        {
            wavefunction += lines[indexx];
        }
        pos = wavefunction.length();
        Separate_string(wavefunction, line, pos);
        int counter = 0;
        for(int i2=0; i2<Npoints[2]; i2++)
        {
            for(int i1=0; i1<Npoints[1]; i1++)
            {
                for(int i0=0; i0<Npoints[0]; i0++)
                {
                    Psi[i][i0][i1][i2] = atof(line[counter].c_str());
                    counter++;
                }
            }

        }  
        line.clear(); 
        //cout << i << " " << Psi[i][Npoints[0]-1][Npoints[1]-1][Npoints[2]-1] << endl;
    }
}



void Plotxsf_X0_Y_Zmiddle(int& at)
{
    ofstream fp_out;
    fp_out.open("Psivalence"+to_string(at)+"_X0_Y_Zmiddle.txt");

    for(int ix=0; ix<Npoints[0]; ix++)
    {
        for(int iy=0; iy<Npoints[1]; iy++)
        {
            for(int iz=0; iz<Npoints[2]; iz++)
            {
                if(ix==Npoints[1]-iy-1 && iz== int(Npoints[2]/2.))
                {
                    Coord_R r;
                    r.setcrys(ix*spacing[0] + origin.crys[0], iy*spacing[1] + origin.crys[1], iz*spacing[2] + origin.crys[2]);
                    //if(abs(abs(Psi[at][ix][iy][iz]) - 3) < 1.e-01)
                    fp_out << r.cart[0]*space_au_A << " " << r.cart[1]*space_au_A << " "<< r.cart[2]*space_au_A << " " << Psi[at][ix][iy][iz] << endl;
                }
            }
        }

    }
}


void PlotFitted_X0_Y_Zmiddle(FittedData<double>& PsiValence, int& at)
{
    ofstream fp_out;
    fp_out.open("Fitted"+to_string(at)+"_X0_Y_Zmiddle.txt");
    int density = 5;
    ifstream fp_Psi;
    double x, y, z;
    fp_Psi.open("Psivalence"+to_string(at)+"_X0_Y_Zmiddle.txt");
    fp_Psi >> x >> y >> z; //cout << x << " " << y << " "<< z << endl;
    for(int iy=0; iy<density*Npoints[1]; iy++)
    {
        Coord_R r;
        r.setcart(x*space_A_au, y*space_A_au*(1-2*double(iy)/double(density*(Npoints[1]-1))),z*space_A_au);

        fp_out << r.cart[0]*space_au_A << " " << r.cart[1]*space_au_A << " "<< r.cart[2]*space_au_A << " " << PsiValence(r.crys[0], r.crys[1], r.crys[2]) << endl;

    }

}


void Plotxsf_X_Y0_Zmiddle(int& at)
{
    ofstream fp_out;
    fp_out.open("Psivalence"+to_string(at)+"_X_Y0_Zmiddle.txt");

    for(int ix=0; ix<Npoints[0]; ix++)
    {
        for(int iy=0; iy<Npoints[1]; iy++)
        {
            for(int iz=0; iz<Npoints[2]; iz++)
            {
                if(ix==iy && iz== int(Npoints[2]/2.))
                {
                    Coord_R r;
                    r.setcrys(ix*spacing[0] + origin.crys[0], iy*spacing[1] + origin.crys[1], iz*spacing[2] + origin.crys[2]);
                    //if(abs(abs(Psi[at][ix][iy][iz]) - 3) < 1.e-01)
                    fp_out << r.cart[0]*space_au_A << " " << r.cart[1]*space_au_A << " "<< r.cart[2]*space_au_A << " " << Psi[at][ix][iy][iz] << endl;
                }
            }
        }

    }
}



void PlotFitted_X_Y0_Zmiddle(FittedData<double>& PsiValence, int& at)
{
    ofstream fp_out;
    fp_out.open("Fitted"+to_string(at)+"_X_Y0_Zmiddle.txt");
    ifstream fp_Psi;
    double x, y, z;
    fp_Psi.open("Psivalence"+to_string(at)+"_X_Y0_Zmiddle.txt");
    fp_Psi >> x >> y >> z; //cout << x << " " << y << " "<< z << endl;
    for(int iy=0; iy<density*Npoints[1]; iy++)
    {
        Coord_R r;
        r.setcart(x*space_A_au*(1-2*double(iy)/double(density*(Npoints[1]-1))), y*space_A_au, z*space_A_au);

        fp_out << r.cart[0]*space_au_A << " " << r.cart[1]*space_au_A << " "<< r.cart[2]*space_au_A << " " << PsiValence(r.crys[0], r.crys[1], r.crys[2]) << endl;

    }

}




void interpolate()
{
    Psi_all.resize(nat);

    for(int at =0; at<nat; at++)
    {
        Psi_all[at].construct(Psi[at],spacing,origin.crys);
        Plotxsf_X0_Y_Zmiddle(at);//for test
        PlotFitted_X0_Y_Zmiddle(Psi_all[at],at);//for test
        Plotxsf_X_Y0_Zmiddle(at);//for test
        PlotFitted_X_Y0_Zmiddle(Psi_all[at], at);//for test
    }
}



double Integrate3D(vec3d& f, double& Jacobian, vec1d& deltax)
{//omp is working 4 times faster if you use 12 threads instead of 1
    double I = 0.;
    double factor;
    omp_set_num_threads(1);

#pragma omp parallel for reduction(+: I)
    for(int ix=0; ix<f.n1(); ix++ )
    {
        double I_x = 0., I_xy = 0.;
        for(int iy=0; iy<f.n2(); iy++ )
        {
            for(int iz=0; iz<f.n3(); iz++ )
            {
                factor = (iz==0 || iz==f.n3()-1 ? .5 : 1);
                I_xy += factor*f[ix][iy][iz];
            }
            factor = (iy==0 || iy==f.n2()-1 ? .5 : 1);
            I_x += factor*I_xy;
            I_xy = 0.;
        }
        factor = (ix==0 || ix==f.n1()-1 ? .5 : 1);
        I += I_x; //cout << I << endl;
        I_x = 0.;
    }
    I *= Jacobian*deltax[0]*deltax[1]*deltax[2]; //cout << I <<" " << deltax[0] << " " << deltax[1] <<" "<< deltax[2] <<  endl;
    return I;
}


void Normalize()
{
    for(int at=0; at<nat; at++)
    {
        vec3d Psi2(Psi[at].n1(), Psi[at].n2(), Psi[at].n3());
        double normsquared = 0.;
        for(int ix = 0; ix < Psi[at].n1(); ix++)
        {
            for(int iy = 0; iy < Psi[at].n2(); iy++)
            {
                for(int iz = 0; iz < Psi[at].n3(); iz++)
                {
                    Psi2[ix][iy][iz] = Psi[at][ix][iy][iz]*Psi[at][ix][iy][iz]; 
                }
            }
        }
        printf("normalization of %d-th wannier function...\n", at);
        double Jacobian = abs(Coord_R::getJ()); 
        normsquared = Integrate3D(Psi2, Jacobian, spacing); 
        printf("norm = %f\n", sqrt(normsquared));
        for(int ix = 0; ix < Psi[at].n1(); ix++)
        {
            for(int iy = 0; iy < Psi[at].n2(); iy++)
            {
                for(int iz = 0; iz < Psi[at].n3(); iz++)
                {
                    Psi[at][ix][iy][iz] /= sqrt(normsquared);
                }
            }
        }
    }
}
















void test_integral() // test function: e^(-x^2)*e^(-y^2)*e^(-z^2)
{   //result of the test :: 5.56163 while the analytical solution is = sqrt(pi^3) = 5,568081662
    int n=300;
    vec3d f(n,n,n); f.fill(0.);
    ofstream f_f;
    f_f.open("f.txt");
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            for(int k=0; k<n; k++)
            {
                //f[i][j][k] = exp(-(20.*(-.5 +i/(n-1))))*exp(-(20.*(-.5 +j/(n-1))));//*exp(-(20.*(-.5 +k/(n-1))));
                f[i][j][k] = exp(-5*5*(-.5 +double(i)/double(n-1))*(-.5 +double(i)/double(n-1)))*exp(-5*5*(-.5 +double(j)/double(n-1))*(-.5 +double(j)/double(n-1)))*exp(-5*5*(-.5 +double(k)/double(n-1))*(-.5 +double(k)/double(n-1)));
                //cout << i << " " << j << endl;
                //if(k==0) f_f << i << " "<< " " << j << " " << f[i][j][k] << endl;
            }
        }
    }
    double Jacobian = 1.;
    spacing[0] = 5./(n-1);
    spacing[1] = 5./(n-1);
    spacing[2] = 5./(n-1);

    double start = omp_get_wtime();    
    for(int i=0; i<20; i++) cout << Integrate3D(f, Jacobian, spacing);
    cout << "Elapsed time in seconds: " << omp_get_wtime()-start << endl;
}