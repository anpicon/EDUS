#include "Coordinate.h"
//#include <exception>
class Wannier
{
  public:
    Wannier(){};
    void init(string name, int Nch, int Ncv);
    //complexd operator()(int& ic, int&jc, double& kx, double& ky, double& kz);
    complexd energy(int& ic, int& jc, Coord_B& k);
    void energy(vec2x&, Coord_B&);
    void energy_U(vec2x&, vec2x&, Coord_B&);
    complexd Xgradientenergy(int&ic, int& jc, Coord_B& k);
    complexd Ygradientenergy(int&ic, int& jc, Coord_B& k);
    complexd Zgradientenergy(int&ic, int& jc, Coord_B& k);
    double energy_Bloch(int& ic, int& jc, Coord_B& k);
    complexd dipolex(int& ic, int& jc, Coord_B& k);
    complexd dipoley(int& ic, int& jc, Coord_B& k);
    complexd dipolez(int& ic, int& jc, Coord_B& k);
    void dipole(vec3x&, Coord_B&);
    //complexd unitary(int& ic, int& jc, Coord_B& k);
    void GradientEnergy(vec3x& GE, Coord_B& k);
    string _name;
    int nrpts;

  private:
    int _Ncv;
    int _Nch;
    int nwannier;
    vec3x Hvv;
    vec1d Hcc;
    vec2d Rvv;
    vec2d Rcc;
    vec2d Rcv;
    vec1d ndegen;
    vec3x xvv;
    vec3x yvv;
    vec3x zvv;
    vec3d xcc;
    vec3d ycc;
    vec3d zcc;
    vec3d xcv;
    vec3d ycv;
    vec3d zcv;
};


void Wannier::init(string name, int Nch, int Ncv)
{
    _name = name;
    _Ncv = Ncv;
    _Nch = Nch;
   ifstream fp_input;
   string sname=_name + "_tb.dat";
   fp_input.open(sname.c_str());
   nwannier=0;


            /****************************************************************************
            ******************       Reading the _tb file       *************************
            *************          from the w90 calculation        **********************
            *****************************************************************************
            ****************************************************************************/
/* reading wannier output */
   if (fp_input)
   {
       string trash;
       for(int ii=0; ii<4; ii++)   getline(fp_input, trash);
       fp_input >> nwannier;   //number of wannier functions
       fp_input >> nrpts;      //number of points in the wigner seitz grid
       Rvv.resize(nrpts, 3);
       Hvv.resize(nwannier, nwannier, nrpts); //matrix HR[nrpts][nwannier][nwannier]
       ndegen.resize( nrpts);
       for (int i=0; i<nrpts; i++)
       {
         fp_input >> ndegen[i]; //degeneracy of the i-th R point
       }
       for(int iR=0; iR<nrpts; iR++)
       {
             getline(fp_input, trash);
             double Hr=0.;
             double Hi=0.;
             fp_input >> Rvv[iR][0] >> Rvv[iR][1] >> Rvv[iR][2];
    
             for(int in=0; in<nwannier; in++)
             {    for(int im=0; im<nwannier; im++)
                  {
                       fp_input >>trash;
                       fp_input >> trash;
                       fp_input >> Hr;
                       Hvv[im][in][iR].real(Hr/ndegen[iR]*energy_eV_au);
                       fp_input >> Hi;
                       Hvv[im][in][iR].imag(Hi/ndegen[iR]*energy_eV_au);
                      // cout <<iR<<"   " << im << "    " << in << endl;
 //cout << iR << " " << Rvv[iR][0] << " " << Rvv[iR][1] << " " << Rvv[iR][2] << " " << Hr << " " << Hi <<  endl;    

                  }//end im
              }//end in
         }//end iR

        string trashD;
        //getline(fp_input, trashD);
        //cout << trashD;
        //fp_input >> trashD;
        double Dr, Di;
        xvv.resize(nwannier, nwannier, nrpts);          //components of the r matrix (along x, y, z)
        yvv.resize(nwannier, nwannier, nrpts);           // The matrix elements are
        zvv.resize(nwannier, nwannier, nrpts);           // <m0|r|Rn>
        for(int iR=0; iR<nrpts; iR++)
        {
           getline(fp_input, trash);
           getline(fp_input, trash); 
           getline(fp_input, trash); 
           for(int in=0; in<nwannier; in++)
           {
              for(int im=0; im<nwannier; im++)
              {
 //cout << iR << " " << Rvv[iR][0] << " " << Rvv[iR][1] << " " << Rvv[iR][2] << " ";
                    fp_input >> trash; 
                    fp_input >> trash; 
                    //------read the matrix elements of x, y, z.
                    fp_input >> Dr;
                    fp_input >> Di;
                    xvv[im][in][iR].real(Dr/ndegen[iR]*space_A_au);
                    xvv[im][in][iR].imag(Di/ndegen[iR]*space_A_au);
                    //cout  << Dr << " " << Di << " ";
                    fp_input >> Dr;
                    fp_input >> Di;
                    yvv[im][in][iR].real(Dr/ndegen[iR]*space_A_au);
                    yvv[im][in][iR].imag(Di/ndegen[iR]*space_A_au);
                    //cout  << Dr << " " << Di << " ";
                    fp_input >> Dr;
                    fp_input >> Di;
                    zvv[im][in][iR].real(Dr/ndegen[iR]*space_A_au);
                    zvv[im][in][iR].imag(Di/ndegen[iR]*space_A_au);
                    //cout  << Dr << " " << Di << " ";
                    //cout << endl;     
 
              }
           }//end im
        }//end in

}
 else
 {
     printf(  "No path %20s file found. \n",sname.c_str());
     //cerr << "No path graphite_tb.dat file found.\n";
     exit(1); 
 }
 fp_input.close();

/* reading core calculation output */

 if(_Nch!=0)
 {
    int neighbors = 27;
    Hcc.resize(_Nch);
    Rcc.resize(neighbors,3);
    xcc.resize(_Nch, _Nch, neighbors);
    ycc.resize(_Nch, _Nch, neighbors);
    zcc.resize(_Nch, _Nch, neighbors);
    int neighborscv=5*5*5;
    Rcv.resize(neighborscv,3);
    xcv.resize(_Nch,_Ncv-_Nch, neighborscv);
    ycv.resize(_Nch,_Ncv-_Nch, neighborscv);
    zcv.resize(_Nch,_Ncv-_Nch, neighborscv);

    fp_input.open("dipole_cc.dat");
    if(fp_input)
    {
    	string trash;
        for(int iR=0; iR<neighbors; iR++)
        {
            fp_input >> Rcc[iR][0] >> Rcc[iR][1] >> Rcc[iR][2];
            for(int icjc=0; icjc<nwannier*nwannier; icjc++)
            {    
            	int ic,jc;
                fp_input >> ic >> jc;
                fp_input >> xcc[ic][jc][iR] >> ycc[ic][jc][iR] >> zcc[ic][jc][iR];
                //cout <<setw(20) << "iR = " << iR;
                //cout <<setw(20) << Rcc[iR][0];
                //cout <<setw(20) << Rcc[iR][1]; 
                //cout <<setw(20)<< Rcc[iR][2];
                //cout <<setw(20) << setprecision(4) << xcc[ic][jc][iR];
                //cout <<setw(20) << setprecision(4) << ycc[ic][jc][iR];
                //cout <<setw(20) << setprecision(4) << zcc[ic][jc][iR]<<endl;
                getline(fp_input, trash);
            }//end icjc
        }//end IR
    }
    else
    {
    printf(  "No path dipole_cc.dat file found. \n");
    //cerr << "No path graphite_tb.dat file found.\n";
    exit(1); 
    } 
    fp_input.close();
    fp_input.open("dipole_cv.dat");
    if(fp_input)
    {
    	string trash;
        for(int iR=0; iR<neighborscv; iR++)
        {
            fp_input >> Rcv[iR][0] >> Rcv[iR][1] >> Rcv[iR][2];
            for(int icjc=0; icjc<_Nch*(_Ncv-_Nch); icjc++)
            {    
            	int ic,jc;
                fp_input >> ic >> jc;
                //cout << ic << jc << endl;
                bool invert=true;
                //if(ic>=_Nch) 
                //	{ic -=_Nch; invert=true;}
                //else if(jc>=_Nch) 
                //	{jc -=_Nch; invert=false;}
                //else 
                //	cout <<"not required value in dipole_cv.dat\n"; 
            
                double x,y,z;
                fp_input >> x>>y>>z;
                int i,j; 
                i=ic; j=jc;
                //if(invert) {i=jc; j=ic;}
                //else {i=ic; j=jc;}
                xcv[i][j][iR] +=x;
                ycv[i][j][iR] +=y;
                zcv[i][j][iR] +=z;
            }//end icjc
            getline(fp_input, trash);
            getline(fp_input, trash);
            //for(int ic=0; ic<_Nch; ic++)
            //{
            //	for(int jc=0; jc<_Ncv-_Nch; jc++)
            //	{
            //		xcv[ic][jc][iR]*=.5; ycv[ic][jc][iR]*=.5; zcv[ic][jc][iR]*=.5; 
            //    cout <<setw(20) << "iR = " << iR;
            //    cout <<setw(20) << Rcv[iR][0];
            //    cout <<setw(20) << Rcv[iR][1]; 
            //    cout <<setw(20) << Rcv[iR][2];
            //		cout <<setw(20) << setprecision(4) << xcv[ic][jc][iR];
            //    	cout <<setw(20) << setprecision(4) << ycv[ic][jc][iR];
            //    	cout <<setw(20) << setprecision(4) << zcv[ic][jc][iR]<<endl;

            //	}
            //}

        }//end IR
    }
    else
    {
    printf(  "No path dipole_cv.dat file found. \n");
    //cerr << "No path graphite_tb.dat file found.\n";
    exit(1); 
    } 
    fp_input.close();

    fp_input.open("Hcc.dat");
    if(fp_input)
    {
        Hcc.resize(_Nch);
        for(int i=0; i<_Nch; i++)
          fp_input >> Hcc[i];
    }
    else
    {
    printf(  "No path Hcc.dat file found. \n");
    //cerr << "No path graphite_tb.dat file found.\n";
    exit(1); 
    } 
    fp_input.close();


 }
 
     
} //end of the constructor Wannier::Wannier



complexd Wannier::energy(int& ic, int& jc, Coord_B& k)
{
  if(ic >= _Nch && jc >= _Nch)
  {
    int i = ic - _Nch;
    int j = jc - _Nch; 
    complexd value=0.;
    for(int iR=0; iR<nrpts; iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*Rvv[iR][l];
        kR *= 2*pi;
        value+= exp(c1*kR)* Hvv[i][j][iR];
    }//loop R
    return value;
  }
  else if (ic<_Nch && jc < _Nch)
  {
    complexd value=0.;
    if(ic==jc)
    {
        return Hcc[ic];
    }
    else return 0;
  }
  else return 0.;
}


void Wannier::GradientEnergy(vec3x& GE, Coord_B& k)
{
    GE.fill(0.);
    for(int ic = 0; ic < _Ncv; ic++)
    {
        for(int jc = 0; jc<_Ncv; jc++)
        {
            if (ic<_Nch || jc<_Nch) {
                for(int l=0; l<3; l++) GE[ic][jc][l]=0.;
            }
            else{
                int i = ic - _Nch;
                int j = jc - _Nch;
                for(int iR=0; iR<nrpts; iR++)
                {
                    double kR=0.;
                    for(int l=0; l<3; l++) kR += k.crys[l]*Rvv[iR][l];
                    kR *= 2.*pi;
                    for(int l=0; l<3; l++) GE[ic][jc][l]+= c1*2.*pi*Rvv[iR][l]*exp(c1*kR)*Hvv[i][j][iR];
                }//loop R
            }
        }
    }
}





complexd Wannier::dipolex(int& ic, int& jc, Coord_B& k)
{
  if(ic >= _Nch && jc >= _Nch)
  {
    int i = ic - _Nch;
    int j = jc - _Nch; 
    complexd value=0.;
    for(int iR=0; iR<nrpts; iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*Rvv[iR][l];
        kR *= 2*pi;
        value += exp(c1*kR)* xvv[i][j][iR];
    }//loop R
    return value;
  }
  else if (ic<_Nch && jc < _Nch)
  {
    complexd value=0.;
    for(int iR=0; iR<xcc.n3(); iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*Rcc[iR][l];
        kR *= 2*pi;
        value+= exp(c1*kR)* xcc[ic][jc][iR];
    }//loop R
    return value;
  }
  else if (ic >= _Nch && jc < _Nch)//we can put it first to save time in the if statement when we need the conj
  {
    int i=ic-_Nch;
    complexd value=0.;
    for(int iR=0; iR<xcv.n3(); iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*Rcv[iR][l];
        kR *= 2*pi;
        value+= exp(c1*kR)* xcv[i][jc][iR];
    }//loop R
    return value;
  }//loop R
  else return conj(Wannier::dipolex(jc,ic,k));
}



complexd Wannier::dipoley(int& ic, int& jc, Coord_B& k)
{
  if(ic >= _Nch && jc >= _Nch)
  {
    int i = ic - _Nch;
    int j = jc - _Nch; 
    complexd value=0.;

    for(int iR=0; iR<nrpts; iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*Rvv[iR][l];
        kR *= 2*pi;
        value += exp(c1*kR)* yvv[i][j][iR];
    }//loop R
    return value;
  }
  else if (ic<_Nch && jc < _Nch)
  {
    complexd value=0.;
    for(int iR=0; iR<ycc.n3(); iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*Rcc[iR][l];
        kR *= 2*pi;
        value+= exp(c1*kR)* ycc[ic][jc][iR];
    }//loop R
    return value;
  }
  else if (ic >= _Nch && jc < _Nch)
  {
    int i=ic-_Nch;
    complexd value=0.;
    for(int iR=0; iR<ycv.n3(); iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*Rcv[iR][l];
        kR *= 2*pi;
        value+= exp(c1*kR)* ycv[i][jc][iR];
    }//loop R
    return value;
  }//loop R
  else return conj(Wannier::dipoley(jc,ic,k));
}



complexd Wannier::dipolez(int& ic, int& jc, Coord_B& k)
{
  if(ic >= _Nch && jc >= _Nch)
  {
    int i = ic - _Nch;
    int j = jc - _Nch; 
    complexd value=0.;
    for(int iR=0; iR<nrpts; iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*Rvv[iR][l];
        kR *= 2*pi;
        value += exp(c1*kR)* zvv[i][j][iR];
    }//loop R
    return value;
  }
  else if (ic<_Nch && jc < _Nch)
  {
    complexd value=0.;
    for(int iR=0; iR<zcc.n3(); iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*Rcc[iR][l];
        kR *= 2*pi;
        value+= exp(c1*kR)* zcc[ic][jc][iR];
    }//loop R
    return value;
  }
  else if (ic >= _Nch && jc < _Nch)
  {
    int i=ic-_Nch;
    complexd value=0.;
    for(int iR=0; iR<zcv.n3(); iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*Rcv[iR][l];
        kR *= 2*pi;
        value+= exp(c1*kR)* zcv[i][jc][iR];
    }//loop R
    return value;
  }//loop R
  else return conj(Wannier::dipolez(jc,ic,k));
}





void Wannier::energy_U(vec2x& H, vec2x& Uk, Coord_B& k)    
{

//k.crys[0] = 0.;
//k.crys[1] = 0.33333333;
//k.crys[2] =-0.33333333;
    H.fill(0.);
    for(int ic = 0; ic < _Ncv; ic++)
    {
        for(int jc = 0; jc<_Ncv; jc++)
        {
            H[ic][jc] = energy(ic, jc, k);
            if(ic==jc) 
              {
                H[ic][jc].imag(0.);
              }      
        }
    }

    //copy hamiltonian in U to get eigenvectors
    for(int ic = 0; ic < _Ncv; ic++)
    {
        for(int jc = 0; jc<_Ncv; jc++)
        {
          Uk[ic][jc] = H[ic][jc];
        }        
    }
    //V->compute eigenvalues and eigenvectors
    std::vector<double> Energy(_Ncv);
    MKL_INT info = LAPACKE_zheev( LAPACK_ROW_MAJOR, 'V', 'L', _Ncv, &(Uk[0][0]), _Ncv, &(Energy[0]) );
    //get conjugate transpose for Uk
    for(int ic=0; ic<_Ncv; ic++)
    {
         for(int jc = 0; jc<_Ncv; jc++)
        {
            complexd u = Uk[jc][ic];
            Uk[jc][ic] = conj(Uk[ic][jc]);
            Uk[ic][jc] = conj(u);
        }
    }
    //vec epsilon;
    //cx_mat U;
    //cx_mat Hw; 
    //epsilon.zeros(nwannier);
    //U.zeros(nwannier, nwannier);
    //Hw.zeros(nwannier, nwannier);
    ////copy the Ek matrix in an armadillo matrix
    //for (int ic=_Nch; ic<_Ncv; ic++)
    //{
    //    for (int jc=_Nch; jc < _Ncv; jc++)
    //    {
    //        Hw(ic-_Nch, jc-_Nch) = H[ic][jc];
    //    }
    //}
    //eig_sym(epsilon, U, Hw);




//if(k.crys[1] < 0.06 && k.crys[2] < 0.06)
//{
//cout << "k " << k.crys[0] << " " << k.crys[1] << " " << k.crys[2] << endl;
//cout << "H" << endl;
//cout << Hw << endl;
//cout << "U" << endl;
//cout << U.t() << endl;
//cout << "eigenvalues" << endl;
//cout << epsilon << endl;
//}
//cx_mat UU;
//UU.zeros(nwannier,nwannier);      

//  Uk.fill(0.);
//  for(int ic=0; ic<_Nch; ic++)
//  {
//      Uk[ic][ic] = 1.;
//  }
//
//  for (int ic=_Nch; ic<_Ncv; ic++)
//  {
//      for (int jc=_Nch; jc < _Ncv; jc++)
//      {
//          complexd el = U(jc-_Nch,ic-_Nch);
//          Uk[ic][jc]=conj(el);
//      }    
//  } 

//if(k.crys[1] < 0.06 && k.crys[2] < 0.06)
//cout << "UHUd" << endl;
//cout << UU*Hw*UU.t() << endl;
//cout << "UdHU " << endl;
//cout << UU.t()*Hw*UU << endl;
//cout << endl << endl;
//}
}

void  Wannier::energy(vec2x& H, Coord_B& k)
{
    for(int ic = 0; ic < _Ncv; ic++)
    {
        for(int jc = 0; jc<_Ncv; jc++)
        {
            H[ic][jc] = energy(ic, jc, k);

            if(ic==jc) 
                H[ic][jc].imag(0.);
 
        }
    }
}

void  Wannier::dipole(vec3x& Dx, Coord_B& k)
{
    for(int ic = 0; ic < _Ncv; ic++)
    {
        for(int jc = 0; jc<_Ncv; jc++)
        {
            Dx[ic][jc][0] = dipolex(ic, jc, k);
            Dx[ic][jc][1] = dipoley(ic, jc, k);
            Dx[ic][jc][2] = dipolez(ic, jc, k);
            if(ic==jc) 
              for(int i=0; i<3; i++) Dx[ic][jc][i].imag(0.);
        }
      
    }
/*
    for(int ic = 0; ic < _Ncv; ic++)
    {
        for(int jc = 0; jc<_Ncv; jc++)
        {
    cout << Dx[ic][jc][0]*Dx[ic][jc][0] + Dx[ic][jc][1]*Dx[ic][jc][1] + Dx[ic][jc][2]*Dx[ic][jc][2]<< " ";
    }
    cout << endl ; 
}cout << endl;*/
}

