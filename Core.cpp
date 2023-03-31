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


/*******************************************************************************************************
********************************************************************************************************
**********************          This function can produce atomic          ******************************
**********************     orbitals as the ones described in the paper    ******************************
**********************     and replication in the neighborhood             ***************************** 
********************************************************************************************************
********************************************************************************************************/



class Core
{
	public:
		Core(vec1d& Origin, vec1d& at_coord, string& type);
		//virtual ~Core()  {  };
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
	double r = sqrt(xtot*xtot + ytot*ytot + ztot*ztot); // r in a.u.				   //this way we center the atom in the origin O.
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
	if(_type == "C1s")	value /= 2.*sqrt(pi);
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


void loading_corecoordinates(string& filename, multivec1D<string>& elements, int& ncore, vec1d& Origin, vec2d& a, vec2d& R)
{
	ifstream fp_input;
	double value;
    string seedname=filename+"_00001" + ".xsf"; //to_string(valence) + ".xsf";
    cout << seedname<< endl;
    fp_input.open(seedname.c_str());
	ncore = 0;
		
  		string trash;
   	    for(int i=0; i<6; i++) getline(fp_input, trash); //cout << trash << endl;} 
   		for(int i=0; i<3; i++) for(int j=0; j<3; j++)   { fp_input >> a[i][j];  a[i][j] *= space_A_au;}
   	    for(int i=0; i<6; i++) getline(fp_input, trash); //cout << trash << endl;} 
   	    int nat;
    	fp_input >> nat; cout << nat; //cout << nat  << endl;
    	fp_input >> trash;
    	string symb;
    	for (int i=0; i < nat; i++)
		{
			//x_at.push_back(Origin);
			fp_input >> symb;
			for(int k=0; k<elements.n1(); k++)
			{	
				if(symb == elements[k])
					ncore++;
			}
			for (int j=0; j<3; j++)
			{
				double atom;
				fp_input >> atom;
			}
		}
		for(int i=0; i<7; i++) {getline(fp_input, trash); }


	    for(int i=0; i<3; i++) 
		{
			fp_input >> Origin[i];	
			Origin[i] *= space_A_au;
		}

		for(int i=0; i<3; i++) 
		{
			for(int j=0; j<3; j++)
			{
				fp_input >> R[i][j];
			    cout << R[i][j] << endl;
				R[i][j] *= space_A_au;
			}
		}
fp_input.close();
}



void loading2(string& filename, vec2d& x_at, vec1i& ix_at, int& ncore, multivec1D<string>& elements)
{


	ifstream fp_input;
	double value;
    string seedname=filename+"_00001" + ".xsf"; //to_string(valence) + ".xsf";
    cout << seedname<< endl;
    fp_input.open(seedname.c_str());
    //x_at.resize(ncore, 3);
  	string trash;
    for(int i=0; i<14; i++) {getline(fp_input, trash); cout << trash << endl; }

    int nat;
    fp_input >> nat;  //cout << nat  << endl;
    fp_input >> trash;
    string symb;
    int count = 0;
    for (int i=0; i < nat; i++)
	{    

		fp_input >> symb;

			for(int k=0; k<elements.n1(); k++)
			{	
				if(symb == elements[k])
				{
					ix_at[i] = k;
					for(int coor=0; coor<3; coor++) 
					{	    
						fp_input >> x_at[count][coor];
						x_at[count][coor] *= space_A_au;
					}
				}
				else getline(fp_input, trash);
				count++;
			}

	}




    fp_input.close();

}


void loading_valence1 (string filename, int& Nx, int& Ny, int& Nz, int valence) 
{
	stringstream ss;
    ss << valence;
	ifstream fp_input;
	double value;
    string seedname=filename+"_0000" + ss.str() + ".xsf"; //to_string(valence) + ".xsf";
    cout << seedname << endl;
    fp_input.open(seedname.c_str());
  		string trash;
   	    for(int i=0; i<14; i++) {getline(fp_input, trash); }
   		
   	    int nat;
    	fp_input >> nat;  //cout << nat  << endl;
    	fp_input >> trash;
    	string symb;
		
    	for (int i=0; i < nat; i++)
    		getline(fp_input, trash);

		//for (int i=0; i < nat; i++)
		//{
		//	x_at.push_back(Origin);
		//	fp_input >> symb;
		//	for (int j=0; j<3; j++)
		//	{
		//		double num;
		//		fp_input >> num;
		//		x_at[i].push_back(num);
		//	}
		//}

		for(int i=0; i<6; i++) {getline(fp_input, trash); }


		fp_input >> Nx >> Ny >> Nz;  // cout << Nx <<"   " << Ny << "    "  << Nz;
fp_input.close();
}



void loading_valence2(string& filename, vec3d& PsiValence,int& valence)
{
	stringstream ss;
    ss << valence;
	ifstream fp_input;
	double value;
    string seedname=filename+"_0000" + ss.str() + ".xsf"; //to_string(valence) + ".xsf";
    cout << seedname << endl;
  		string trash;
    fp_input.open(seedname.c_str());
   	    for(int i=0; i<14; i++) {getline(fp_input, trash); cout << trash << endl;}
   		
   	    int nat;
    	fp_input >> nat;  //cout << nat  << endl;
    	fp_input >> trash;
    	string symb;
		
    	for (int i=0; i < nat; i++)
    		getline(fp_input, trash);

		for(int i=0; i<(7+4); i++) {getline(fp_input, trash);
cout << trash << endl;}
		for(int iz=0; iz<PsiValence.n3(); iz++)
		{
			for(int iy=0; iy<PsiValence.n2(); iy++)
			{
				for(int ix=0; ix<PsiValence.n1(); ix++)
				{
					double num;
					fp_input >> num;
					PsiValence[ix][iy][iz] = num;
				}//end ix
			}//end iy
		}//end iz
		fp_input.close();
}



template<class T>
double CavalieriIntegral_norm(int nx, int ny, int nz, T& PsiValence, double Detj)
{
  double Ixyz    =       0.;
  double Ixy     =       0.;   
  double Ix      =       0.;  
  double I       =       0.;

    //cout << Detj << endl;
    //cout << "nx " << nx << "  ny " << ny << "    nz " << nz << endl;
    //cout << "Origin"  <<Origin.cart[0] << "     "<< Origin.cart[1] << "    " <<Origin.cart[2]<<endl;


    int xin = 0; 
    int xf = nx;
    int yin = 0; 
    int yf = ny;
    int zin = 0; 
    int zf = nz;



      for( int ix = xin; ix < xf; ix++ )
      {
       // cout << ix; double a;
        //cout << "iz      " << iz << "      ";
        double x=double(ix)/double(nx);
        for( int iy = yin; iy < yf; iy++ )
        {
          //cout << "iy       "<<  iy << "      " << endl;
                double y=double(iy)/double(ny);
            
            for( int iz  = zin; iz < zf; iz++ )
            {
              //cout << ix << "       ";
              double z=double(iz)/double(nz);

              Coord_R k;
              k.setcrys(x,y,z); 

              Ixyz += PsiValence(k.crys[0],k.crys[1],k.crys[2])*PsiValence(k.crys[0],k.crys[1],k.crys[2]);
              
                if      ( iz == zin || iz == zf - 1 )   Ixy +=    Ixyz; 
                else if ( iz%2 == 0 )               Ixy += 2.*Ixyz;
                else                                Ixy += 4.*Ixyz;
                
                Ixyz = 0.;

            }//end ix

            if      ( iy == yin || iy == yf - 1 )   Ix +=    Ixy; 
            else if ( iy%2 == 0 )               Ix += 2.*Ixy;
            else                                Ix += 4.*Ixy;

            Ixy = 0.;

        }//end iy

        if      ( ix == xin || ix == xf-1 )   I +=    Ix; 
        else if ( ix%2 == 0 )               I += 2.*Ix;
        else                                I += 4.*Ix;
        
        Ix = 0.;
    }//end iz
      
    double factor=1.;
    if (nx > 1) factor/=3.;
    if (ny > 1) factor/=3.;
    if (nz > 1) factor/=3.;

    I *= Detj/(nx*ny*nz)*factor;
    
    printf ("norm = %3.3f \n" , I);
    return I;
}





template<class T>
double CavalieriIntegral_X(int nx, int ny, int nz, vec1d& Origin, Core& PsiCore, T& PsiValence, double& Detj, double& norm)
{
	double Ixyz    =       0.;
	double Ixy     =       0.;   
	double Ix      =       0.;	
	double I       =       0.;


      for( int ix = 0; ix < nx; ix++ )
      {
      	double x=double(ix)/double(nx-1)-0.5;
        for( int iy = 0; iy < ny; iy++ )
        {
              	double y=double(iy)/double(ny-1)-0.5;
            
            for( int iz  = 0; iz < nz; iz++ )
            {
            	double z=double(iz)/double(nz-1)-0.5;

            	Coord_R k;
            	k.setcrys(x,y,z); 

				Ixyz += (k.crys[0]+Origin[0])*PsiCore(k.cart[0],k.cart[1],k.cart[2])*PsiValence(k.crys[0],k.crys[1],k.crys[2]);

                if      ( iz == 0 || iz == nz-1 )   Ixy +=    Ixyz; 
                else if ( iz%2 == 0 )               Ixy += 2.*Ixyz;
                else                                Ixy += 4.*Ixyz;
                
                Ixyz = 0.;
            }//end ix

            if      ( iy == 0 || iy == ny-1 )   Ix +=    Ixy; 
            else if ( iy%2 == 0 )               Ix += 2.*Ixy;
            else                                Ix += 4.*Ixy;

            Ixy = 0.;

        }//end iy

        if      ( ix == 0 || ix == nx-1 )   I +=    Ix; 
        else if ( ix%2 == 0 )               I += 2.*Ix;
        else                                I += 4.*Ix;
        

        Ix = 0.;
    }//end iz

      
    double factor=1.;
    if (nx > 1) factor/=3.;
    if (ny > 1) factor/=3.;
    if (nz > 1) factor/=3.;
    
    I *= Detj/(nx*ny*nz*sqrt(norm))*factor;
    return I;
}




template<class T>
double CavalieriIntegral_Y(int nx, int ny, int nz, vec1d& Origin, Core& PsiCore, T& PsiValence, double Detj, double&  norm)
{
	double Ixyz    =       0.;
	double Ixy     =       0.;   
	double Ix      =       0.;	
	double I       =       0.;


      for( int ix = 0; ix < nx; ix++ )
      {
      	double x=double(ix)/double(nx-1)-0.5;
        for( int iy = 0; iy < ny; iy++ )
        {
              	double y=double(iy)/double(ny-1)-0.5;
            
            for( int iz  = 0; iz < nz; iz++ )
            {
            	double z=double(iz)/double(nz-1)-0.5;

            	Coord_R k;
            	k.setcrys(x,y,z); 

				Ixyz = (y+Origin[1])*PsiCore(k.getcart1(),k.getcart2(),k.getcart3())*PsiValence(x,y,z);

                if      ( iz == 0 || iz == nz-1 )   Ixy +=    Ixyz; 
                else if ( iz%2 == 0 )               Ixy += 2.*Ixyz;
                else                                Ixy += 4.*Ixyz;
                
                Ixyz = 0.;

            }//end ix

            if      ( iy == 0 || iy == ny-1 )   Ix +=    Ixy; 
            else if ( iy%2 == 0 )               Ix += 2.*Ixy;
            else                                Ix += 4.*Ixy;

            Ixy = 0.;

        }//end iy

        if      ( ix == 0 || ix == nx-1 )   I +=    Ix; 
        else if ( ix%2 == 0 )               I += 2.*Ix;
        else                                I += 4.*Ix;
        

        Ix = 0.;
    }//end iz
      
    double factor=1.;
    if (nx > 1) factor/=3.;
    if (ny > 1) factor/=3.;
    if (nz > 1) factor/=3.;
    
    I *= Detj/(nx*ny*nz*sqrt(norm))*factor;
    return I;
}




template<class T>
double CavalieriIntegral_Z(int nx, int ny, int nz, vec1d& Origin, Core& PsiCore, T& PsiValence, double& Detj, double& norm)
//double CavalieriIntegral_Z(int nx, int ny, int nz, vec1d& Origin, T& PsiCore, T& PsiValence, double& Detj, double& norm)
{
	double Ixyz    =       0.;
	double Ixy     =       0.;   
	double Ix      =       0.;	
	double I       =       0.;



      for( int ix = 0; ix < nx; ix++ )
      {
          double x=double(ix)/double(nx-1)-0.5;
        for( int iy = 0; iy < ny; iy++ )
        {
                  double y=double(iy)/double(ny-1)-0.5;
            
            for( int iz  = 0; iz < nz; iz++ )
            {
                double z=double(iz)/double(nz-1)-0.5;

            	Coord_R k;
            	k.setcrys(x,y,z); 

				Ixyz = (z+Origin[2])*PsiCore(k.getcart1(),k.getcart2(),k.getcart3())*PsiValence(x,y,z);

                if      ( iz == 0 || iz == nz-1 )   Ixy +=    Ixyz; 
                else if ( iz%2 == 0 )               Ixy += 2.*Ixyz;
                else                                Ixy += 4.*Ixyz;
                
                Ixyz = 0.;

            }//end ix

            if      ( iy == 0 || iy == ny-1 )   Ix +=    Ixy; 
            else if ( iy%2 == 0 )               Ix += 2.*Ixy;
            else                                Ix += 4.*Ixy;

            Ixy = 0.;

        }//end iy

        if      ( ix == 0 || ix == nx-1 )   I +=    Ix; 
        else if ( ix%2 == 0 )               I += 2.*Ix;
        else                                I += 4.*Ix;
        

        Ix = 0.;
    }//end iz

      
    double factor=1.;
    if (nx > 1) factor/=3.;
    if (ny > 1) factor/=3.;
    if (nz > 1) factor/=3.;
    
    I *= Detj/(nx*ny*nz*sqrt(norm))*factor;
    return I;
}
      


double ScalarProductW(int nx, int ny, int nz, FittedData<double>& PsiValence1, FittedData<double>& PsiValence2, double& Detj, double& norm1,double& norm2)
{
    double Ixyz    =       0.;
    double Ixy     =       0.;
    double Ix      =       0.;
    double I       =       0.;

    double dx=1./double(nx-1);
    double dy=1./double(ny-1);
    double dz=1./double(nz-1);
    double dx2=dx/2.;
    double dy2=dy/2.;
    double dz2=dz/2.;

      for( int ix = 0; ix < nx-1; ix++ )
      {
          double x=double(ix)/double(nx-1)-0.5;
        for( int iy = 0; iy < ny-1; iy++ )
        {
                  double y=double(iy)/double(ny-1)-0.5;
            
            for( int iz  = 0; iz < nz-1; iz++ )
            {
                double z=double(iz)/double(nz-1)-0.5;

                Coord_R k;
                k.setcrys(x,y,z);

                Ixyz = PsiValence1(x+dx2,y+dy2,z+dz2)*PsiValence2(x+dx2,y+dy2,z+dz2);
                
                Ixy +=    Ixyz;
                
                Ixyz = 0.;

            }//end ix

            Ix +=    Ixy;

            Ixy = 0.;

        }//end iy

          I +=    Ix;
        

        Ix = 0.;
    }//end iz

      
    double factor=1.;
    
    I *= Detj*dx*dy*dz/(sqrt(norm1)*sqrt(norm2))*factor;
    return I;
}







int main (int argc, char* argv[])
{
    bool iTest=false;
    ifstream fp_input;
    if ( argc != 2 )
    {
        cout<<"usage: "<< argv[0] <<" <filename>" << endl;
    }
    else
    {
        fp_input.open(argv[1]);
        if (!fp_input.is_open())
        {
            cout << "error opening file " << argv[1] << endl;
        }
    }
    //%%%%%%%%%%%%%%% DEFINING ARRAYS AND OUTPUTS %%%%%%%%%%%%%%%//
	vec2d a; a.resize(3,3); a.fill(0.); //vectors of the unit cell
    vec2d R; R.resize(3,3); R.fill(0.); //vectors of the super cell
    vec1d Origin; Origin.resize(3); Origin.fill(0.); //vector to know where is the origin of the Wannier cube files
    
	int nwannier; //number of Wannier functions
	int nelements; //element type of core holes
	string seedname; //the string filename to call the Wannier files
	multivec1D<string> elements; //array to store the different types of core-hole elements
	multivec1D<string> orbital; //array to store the different core-hole orbitals
    vec2d x_at; //position of atoms with a core hole
    vec1i ix_at; //index for array "elements"
    int ncore = 0; //total number of core orbitals
    vec1i Ncells(3); //number of unit cells in the supercell, 0-> direction of a1, 1-> a2, and 2-> a3
    
    ofstream fp_outputCW; fp_outputCW.open("Core_Wannier_tb.dat"); //REMINDER: check the format of Wannier .dat files
    ofstream fp_outputCC; fp_outputCC.open("Core_Core_tb.dat"); //REMINDER: check the format of Wannier .dat files
    ofstream fp_outputHcc; fp_outputHcc.open("Hcc.dat");

    //%%%%%%%%%%%%%%% READING INPUT %%%%%%%%%%%%%%%//
	//cout << "number of wannier" << endl;
	fp_input >> nwannier; cout << nwannier << endl;
	//cout << "filename of wannier " << endl;
	fp_input >> seedname; cout << seedname << endl;
	//cout << "number of core elements" << endl;
	fp_input >> nelements; cout << nelements << endl;
	elements.resize(nelements);
	orbital.resize(nelements);

	//cout << "type of core elements" << endl;
	for(int i=0; i<nelements; i++)   { fp_input >> elements[i]; cout << elements[i] << endl;}
	//cout << "type of core orbitals" << endl;
	for(int i=0; i<nelements; i++)    {fp_input >> orbital[i]; cout << orbital[i] << endl;}
	fp_input.close();
    
    clock_t clocktime = clock();
    
    //we read ncore (total number of core orbitals), Origin (Wannier cube files), a(vectors of unit cell), and R(supercell)
	loading_corecoordinates(seedname, elements, ncore, Origin, a, R);
    
    //Finding the number of unit cells that composes the super cell
    for (int k=0; k<3; k++) {
        for (int i=0; i<3; i++) {
            if(abs(a[k][i])>1.e-8)
            {
                Ncells[k]=int(round(R[k][i]/a[k][i]));
                //cout << "Ncells " << k << " " << i << " " << Ncells[k] << endl;
            }
        }
    }
    
    //We redefing the supercell R, previously the last point was not included
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            R[i][j]=Ncells[i]*a[i][j];
        }
    }
    
    //we define a supercell vector to change crystal <-> cartesian coordinates
    Coord_R Origine;
    Coord_R::set_crys_to_cart(R);
    Origine.setcart( Origin[0], Origin[1], Origin[2]);
    
    printf("Origin        : (%12.5f, %12.5f, %12.5f) a.u.\n",Origin[0],Origin[1],Origin[2]);
    printf("Unit vector a1: (%12.5f, %12.5f, %12.5f) a.u.\n",a[0][0],a[0][1],a[0][2]);
    printf("Unit vector a2: (%12.5f, %12.5f, %12.5f) a.u.\n",a[1][0],a[1][1],a[1][2]);
    printf("Unit vector a3: (%12.5f, %12.5f, %12.5f) a.u.\n",a[2][0],a[2][1],a[2][2]);
    printf("Supercell   R1: (%12.5f, %12.5f, %12.5f) a.u.\n",R[0][0],R[0][1],R[0][2]);
    printf("Supercell   R2: (%12.5f, %12.5f, %12.5f) a.u.\n",R[1][0],R[1][1],R[1][2]);
    printf("Supercell   R3: (%12.5f, %12.5f, %12.5f) a.u.\n",R[2][0],R[2][1],R[2][2]);
    
    //we read x_at and defining ix_at
    x_at.resize(ncore, 3);
    ix_at.resize(ncore);
	loading2(seedname, x_at, ix_at, ncore, elements);
    
   
   //%%%%%%%%%%%%%%% WANNIER-CORE DIPOLES %%%%%%%%%%%%%%%//
    cout << "Calculating Dipole off diagonal terms between core and valence orbitals...." << endl;
    
    //creating the output file Wannier-like format for the dipoles
    fp_outputCW << "number of core orbitals: " <<  ncore << endl;
    fp_outputCW << "electric dipole transitions:" << endl;
    


   for(int ic=0; ic<nwannier; ic++)
   {
       cout << "Valence # " << ic << endl;
       int nx=300; //resolution of interpolation
       int ny=300; //resolution of interpolation
       int nz=300; //resolution of interpolation
       int N1=0;                                                     //sampling in first direction
       int N2=0;                                                     //sampling in second direction
       int N3=0;                                                     //sampling in third direction
       double Determinant = Coord_R::getJ();
       
       int ic1=ic+1;
       loading_valence1 (seedname, N1, N2, N3, ic1 );
       vec3d PsiValence1;   PsiValence1.resize(N1,N2,N3);                                         //array for the valence band: PsiValence is after spline, PsiValence1 before
       loading_valence2 (seedname, PsiValence1, ic1 );
       cout << "Valence wavefunction loaded..." << endl;
       if(iTest==true) cout << "Psi: " << PsiValence1[N1-2][N2-1][N3-1] << " " << PsiValence1[N1-1][N2-1][N3-1] << endl;
       cout << endl << "Determinant: " << Determinant << endl;
       vec1d spacing; spacing.resize(3);
       spacing[0]=1./double(N1); spacing[1]=1./double(N2); spacing[2]=1./double(N3);
       printf("Spacing       : (%12.5f, %12.5f, %12.5f) a.u.\n",spacing[0],spacing[1],spacing[2]);
       
       //vec1d shift; shift.resize(3);
       //shift[0]=Origine.crys[0]; shift[1]=Origine.crys[1];  shift[2]=Origine.crys[2];

       //interpolating the Wannier function
       FittedData<double> PsiValence(PsiValence1,spacing,Origine.crys);
       
       double norm = CavalieriIntegral_norm(nx, ny, nz, PsiValence, Determinant);
       
       if(iTest==true)
       {
           //checking wave function
           ofstream fp_wan;
           int NNN=100;
           double delta=0.01;
           fp_wan.open("Wannierxz.txt");
           for (int ix=0; ix<NNN; ix++) {
               for (int iz=0; iz<NNN; iz++) {
                   Coord_R xyz;
                   double x=delta*(ix-NNN/2); //delta*(i-NNN/2);
                   double y=0.;
                   double z=delta*(iz-NNN/2);
                   xyz.setcart(x,y,z);
                   fp_wan << x << " "  << z << " " << PsiValence(xyz.crys[0],xyz.crys[1],xyz.crys[2]) << endl;
               }
               fp_wan << endl;
           }
           fp_wan.close();
           fp_wan.open("Wannieryz.txt");
           for (int ix=0; ix<NNN; ix++) {
               for (int iz=0; iz<NNN; iz++) {
                   Coord_R xyz;
                   double x=0.; //delta*(i-NNN/2);
                   double y=delta*(ix-NNN/2);
                   double z=delta*(iz-NNN/2);
                   xyz.setcart(x,y,z);
                   fp_wan << y << " "  << z << " " << PsiValence(xyz.crys[0],xyz.crys[1],xyz.crys[2]) << endl;
               }
               fp_wan << endl;
           }
           fp_wan.close();
           
           fp_wan.open("Wannierx.txt");
           for (int i=0; i<NNN; i++) {
               Coord_R xyz;
               double x=delta*(i-NNN/2); //delta*(i-NNN/2);
               double y=0.;
               double z=0.;
               xyz.setcrys(x,y,z);
               fp_wan << x << " " << PsiValence(xyz.crys[0],xyz.crys[1],xyz.crys[2]) << endl;
           }
           fp_wan.close();
           
           fp_wan.open("Wannierx_nointerp.txt");
           for (int ix=0; ix<N1; ix++) {
               double x=spacing[0]*(ix-N1/2); //delta*(i-NNN/2);
               fp_wan << x << " " << PsiValence1[ix][N2/2+6][N3/2+6] << endl;
           }
           fp_wan.close();
           
           fp_wan.open("Wanniery.txt");
           for (int i=0; i<NNN; i++) {
               Coord_R xyz;
               double x=0.; //delta*(i-NNN/2);
               double y=delta*(i-NNN/2);
               double z=0.;
               xyz.setcart(x,y,z);
               fp_wan << y << " " << PsiValence(xyz.crys[0],xyz.crys[1],xyz.crys[2]) << endl;
           }
           fp_wan.close();
           
           fp_wan.open("Wannierz.txt");
           for (int i=0; i<NNN; i++) {
               Coord_R xyz;
               double x=0.; //delta*(i-NNN/2);
               double y=0.;
               double z=delta*(i-NNN/2);
               xyz.setcart(x,y,z);
               fp_wan << z << " " << PsiValence(xyz.crys[0],xyz.crys[1],xyz.crys[2]) << endl;
           }
           fp_wan.close();
       }
       for(int jc=0; jc<ncore; jc++)
       {
           string atomicorbital = elements[ix_at[jc]] + orbital[ix_at[jc]];//"C1s";

           for (int i=-1; i<2; i++)
           {
               for (int j=-1; j<2; j++)
               {
                   for (int k=-1; k<2; k++)
                   {
                       printf("Calculating dipole matrix core-Wannier elements.... for %1i %1i %1i \n",i, j, k);                             //vec3d PsiCore(nz, vec2d(ny, vec1d(nx,0)));
                       vec1d at_position;
                       at_position.resize(3);
                       for (int l=0; l<3; l++)
                       {
                           at_position[l] =  x_at[jc][l] + i*a[0][l] + j*a[1][l] + k*a[2][l];
                       }
                       Coord_R Orishift;
                       Orishift.setcart(0., 0., 0.);
                       Core PsiCore( Orishift.cart, at_position, atomicorbital);
                       Coord_R r;
                       
                       if (iTest==true) {
                           ofstream fp_wan;
                           int NNN=140;
                           double delta=0.1;
                           cout << "Position core: x " << at_position[0] << " y " << at_position[1] << " z " << at_position[2] << endl;
                           fp_wan.open("Corexz.txt");
                           for (int ix=0; ix<NNN; ix++) {
                               for (int iz=0; iz<NNN; iz++) {
                                   Coord_R xyz;
                                   double x=delta*(ix-NNN/2); //delta*(i-NNN/2);
                                   double y=0.;
                                   double z=delta*(iz-NNN/2);
                                   xyz.setcart(x,y,z);
                                   fp_wan << x << " "  << z << " " << PsiCore(xyz.cart[0],xyz.cart[1],xyz.cart[2]) << endl;
                               }
                               fp_wan << endl;
                           }
                           fp_wan.close();
                           fp_wan.open("Coreyz.txt");
                           for (int ix=0; ix<NNN; ix++) {
                               for (int iz=0; iz<NNN; iz++) {
                                   Coord_R xyz;
                                   double x=0.; //delta*(i-NNN/2);
                                   double y=delta*(ix-NNN/2);
                                   double z=delta*(iz-NNN/2);
                                   xyz.setcart(x,y,z);
                                   fp_wan << y << " "  << z << " " << PsiCore(xyz.cart[0],xyz.cart[1],xyz.cart[2]) << endl;
                               }
                               fp_wan << endl;
                           }
                           fp_wan.close();
                       }

                       double x1=CavalieriIntegral_X(nx, ny, nz, Orishift.crys, PsiCore, PsiValence, Determinant, norm);
                       double x2=CavalieriIntegral_Y(nx, ny, nz, Orishift.crys, PsiCore, PsiValence, Determinant, norm);
                       double x3=CavalieriIntegral_Z(nx, ny, nz, Orishift.crys, PsiCore, PsiValence, Determinant, norm);

                       r.setcrys(x1,x2,x3);        //dipole moment <0val|r|core R> in crystal coord
                       
                       fp_outputCW << ic+ncore << "    " << jc << "    " << i << "    " << j << "     " << k << "    " << r.cart[0]<< "      "   << r.cart[1] << "    " << r.cart[2] << endl;

                   }//end k
               }//end j
           }//end i
       }//end jc
   }//end ic , out of the loop of reading and interpolating Wannier functions
    
    
    //%%%%%%%%%%%%%%% CORE-CORE DIPOLES %%%%%%%%%%%%%%%//
    cout << "Calculating Dipole off diagonal terms between core and core orbitals...." << endl;
    //creating the output file Wannier-like format for the dipoles
    fp_outputCC << "number of core orbitals: " <<  ncore << endl;
    fp_outputCC << "electric dipole transitions:" << endl;
    
    for(int ic=0; ic< ncore; ic++)
    {
        int nx = 300;
        int ny = 300;
        int nz = 300;
        string atomicorbitali = elements[ix_at[ic]] + orbital[ix_at[ic]];
        double Determinant = Coord_R::getJ();
        for(int jc=0; jc<ncore; jc++)
        {
            string atomicorbital = elements[ix_at[jc]] + orbital[ix_at[jc]];//"C1s";
            for (int i=-1; i<2; i++)
            {
                for (int j=-1; j<2; j++)
                {
                    for (int k=-1; k<2; k++)
                    {
                        vec1d at_positioni; at_positioni.resize(3);
                        vec1d at_positionj; at_positionj.resize(3);
                        for (int l=0; l<3; l++)
                        {
                            at_positionj[l] =  x_at[jc][l] + i*a[0][l] + j*a[1][l] + k*a[2][l];
                            at_positioni[l] =  x_at[ic][l];
                        }

                        Coord_R Orishift;
                        Orishift.setcart(0., 0., 0.);
                        Core PsiCorei( Orishift.cart, at_positioni, atomicorbitali);
                        Core PsiCorej( Orishift.cart, at_positionj, atomicorbitali);
                        Coord_R r;
                        double normone = 1.;
                        
                        
                        double x1=CavalieriIntegral_X(nx, ny, nz, Orishift.crys, PsiCorei, PsiCorej, Determinant, normone);
                        double x2=CavalieriIntegral_Y(nx, ny, nz, Orishift.crys, PsiCorei, PsiCorej, Determinant, normone);
                        double x3=CavalieriIntegral_Z(nx, ny, nz, Orishift.crys, PsiCorei, PsiCorej, Determinant, normone);
                                            
                        r.setcrys(x1,x2,x3);        //dipole moment <0val|r|core R> in crystal coord
                        cout << "\n ---> integral done.    Crystal Values = "<<r.crys[0] << "      "  << r.crys[1] << "       " << r.crys[2] << endl;
                        cout << "Cartesian values:    " << r.cart[0] << "      "   << r.cart[1] << "       " << r.cart[2] << endl;
                        
                        fp_outputCC << ic << "    " << jc << "    " << i << "    " << j << "     " << k << "    " << r.cart[0] << "      "   << r.cart[1] << "      " << r.cart[2] << endl;
                    }
                }
            }
        }
    }
    
    
    //%%%%%%%%%%%%%%% CORE-CORE ENERGIES %%%%%%%%%%%%%%%//
    for (int jc=0; jc<ncore; jc++) {
        string atomicorbital = elements[ix_at[jc]] + orbital[ix_at[jc]];
        Core PsiCore( Origin, Origin, atomicorbital);
        fp_outputHcc << PsiCore.energy() << endl;
    }
    
    clock_t clock2 = clock();
    cout << "time required:     " <<  (double) (clock2-clocktime)/CLOCKS_PER_SEC << endl ;

} //end main
