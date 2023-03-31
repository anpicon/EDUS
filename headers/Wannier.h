#include "Coordinate.h"
//#include <exception>
class Wannier
{
  public:
    Wannier(){};
    void Wannier_start(string& name);

    void init(string name, int Nch, int Ncv, double FermiEnergy);
    //complexd operator()(int& ic, int&jc, double& kx, double& ky, double& kz);
    complexd energy(int& ic, int& jc, Coord_B& k);
    void energy(vec2x&, Coord_B&);
    void energy_U(vec2x&, vec2x&, Coord_B&);
    void energy_Bloch(vec1d& energy, Coord_B& k);

    complexd unitary(int& i, int& j, Coord_B& kk);
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
    void print_dipole();
    string _name;
    int nrpts;

  private:
    int _Ncv;
    int _Nch;
    int nwannier;
    vec3x Hvv;
    vec1d Hcc;
    //vec2d Rvv;
    vec2i RvvH;
    vec2i Rvvx;
    vec2i Rvvy;
    vec2i Rvvz;
    vec2i RvvU;
    vec2d Rvc;
    vec2d Rcc;
    vec1d ndegen;
    vec3x xvv;
    vec3x yvv;
    vec3x zvv;
    vec3d xcc;
    vec3d ycc;
    vec3d zcc;
    vec3d xvc;
    vec3d yvc;
    vec3d zvc;
    vec3x Uvv;

    multivec1D<bool> Singularity;

};

#include<Wannier_start.h>

void Wannier::init(string name, int Nch, int Ncv, double FermiEnergy)
{
    printf("starting wannier....\n");
      Wannier_start(name);//
    printf(  "*********************************************************************************************\n");  
    printf(  "*  Number of R points H:                       %3i       %3i                                *\n", int(RvvH.n1()), int(Hvv.n3()));
    printf(  "*  Number of R points XiX:                     %3i       %3i   %3i %3i                      *\n", int(Rvvx.n1()), int(xvv.n1()), int(xvv.n2()),int(xvv.n3()));
    printf(  "*  Number of R points XiY:                     %3i       %3i                                *\n", int(Rvvy.n1()), int(yvv.n3()));
    printf(  "*  Number of R points XiZ:                     %3i       %3i                                *\n", int(Rvvz.n1()), int(zvv.n3()));
    printf(  "*  Number of R points U:                       %3i       %3i                                *\n", int(RvvU.n1()), int(Uvv.n3()));
    printf(  "*********************************************************************************************\n");  

      _Nch = Nch; 
      _Ncv = Ncv;

/* reading core calculation output */
      ifstream fp_input;
 if(_Nch!=0)
 {

    int neighbors = 27;
    Hcc.resize(_Nch);
    Rcc.resize(neighbors,neighbors);
    xcc.resize(_Nch, _Nch, neighbors);
    ycc.resize(_Nch, _Nch, neighbors);
    zcc.resize(_Nch, _Nch, neighbors);


    printf("reading dipole_cc.dat....\n");
    fp_input.open("dipole_cc.dat");
    if(fp_input)
    {
        string trash;
            int ic, jc;
            for(int iR=0; iR<neighbors; iR++)
            {
                fp_input >>  Rcc[iR][0] >> Rcc[iR][1] >> Rcc[iR][2] ;
                //cout <<  " " << Rcc[iR][0] << " " <<  Rcc[iR][1] << " " << Rcc[iR][2] << endl;
                for(int i=0; i<_Nch*_Nch; i++)
                {
                  fp_input >> ic >> jc >> xcc[ic][jc][iR] >> ycc[ic][jc][iR] >> zcc[ic][jc][iR]>> trash;
                  if(abs(xcc[ic][jc][iR]>1.e-08))  xcc[ic][jc][iR]*=space_A_au;
                  else                        xcc[ic][jc][iR] = 0;
                  if(abs(ycc[ic][jc][iR]>1.e-08))  ycc[ic][jc][iR]*=space_A_au;
                  else                        ycc[ic][jc][iR] = 0;
                  if(abs(zcc[ic][jc][iR]>1.e-08))  zcc[ic][jc][iR]*=space_A_au;
                  else                        zcc[ic][jc][iR] = 0;
                }
            }
            printf("Last info in dipole_cc.dat: \n");
            printf("%3i, %3i, %6.3f  %6.3f  %6.3f\n", ic,jc,xcc[_Nch-1][_Nch-1][Rcc.n1()-1]*space_au_A,
                    ycc[_Nch-1][_Nch-1][Rcc.n1()-1]*space_au_A,
                    zcc[_Nch-1][_Nch-1][Rcc.n1()-1]*space_au_A);
    }
    else
    {
    printf(  "No path dipole_cc.dat file found. \n");
    //cerr << "No path graphite_tb.dat file found.\n";
    exit(1); 
    } 
    fp_input.close();

    neighbors = 5*5*5;

    Rvc.resize(neighbors, 3);
    xvc.resize(_Ncv, _Ncv, neighbors);
    yvc.resize(_Ncv, _Ncv, neighbors);
    zvc.resize(_Ncv, _Ncv, neighbors);

    fp_input.open("dipole_cv.dat");
    
    int mex = 2*(_Ncv-_Nch)*_Nch;
    //cout <<" mex = " << mex << endl;
    if(fp_input)
    {
        string trash;

            int ic, jc;
            for(int iR=0; iR<neighbors; iR++)
            {
                fp_input >> Rvc[iR][0] >> Rvc[iR][1] >> Rvc[iR][2];
                //cout << Rvc[iR][0] << " " << Rvc[iR][1] << " " <<  Rvc[iR][2] << endl;
                for(int i=0; i<mex; i++)
                {

                  fp_input >> ic >> jc >> xvc[ic][jc][iR] >> yvc[ic][jc][iR] >> zvc[ic][jc][iR]; //>> trash;   
                  if(abs(xvc[ic][jc][iR]>1.e-08))  xvc[ic][jc][iR]*=space_A_au;
                  else                        xvc[ic][jc][iR] = 0;
                  if(abs(yvc[ic][jc][iR]>1.e-08))  yvc[ic][jc][iR]*=space_A_au;
                  else                        yvc[ic][jc][iR] = 0;
                  if(abs(zvc[ic][jc][iR]>1.e-08))  zvc[ic][jc][iR]*=space_A_au;
                  else                        zvc[ic][jc][iR] = 0;
                }
            }

            printf("Last info in dipole_cv.dat: \n");
            printf("%3i, %3i, %6.3f  %6.3f  %6.3f\n", ic,jc,xvc[ic][jc][Rvc.n1()-1]*space_au_A,
                    yvc[ic][jc][Rvc.n1()-1]*space_au_A,
                    zvc[ic][jc][Rvc.n1()-1]*space_au_A);
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
        printf("Energy of core-hole: %4.3f\n", Hcc[Hcc.n1()-1]);
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
    for(int iR=0; iR<RvvH.n1(); iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*RvvH[iR][l];
        kR *= 2*pi;
        value+= exp(c1*kR)* Hvv[i][j][iR];
    }//loop R
    return value;
  }
  else if (ic<_Nch && jc < _Nch && ic==jc)
  {
    return Hcc[ic];
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
                for(int iR=0; iR<RvvH.n1(); iR++)
                {
                    double kR=0.;
                    for(int l=0; l<3; l++) kR += k.crys[l]*RvvH[iR][l];
                    kR *= 2.*pi;
                    for(int l=0; l<3; l++) GE[ic][jc][l]+= c1*2.*pi*RvvH[iR][l]*exp(c1*kR)*Hvv[i][j][iR];
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
    for(int iR=0; iR<Rvvx.n1(); iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*Rvvx[iR][l];
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
    complexd value=0.;
    for(int iR=0; iR<xvc.n3(); iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*Rvc[iR][l];
        kR *= 2*pi;
        value+= exp(c1*kR)* xvc[ic][jc][iR];
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

    for(int iR=0; iR<Rvvy.n1(); iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*Rvvy[iR][l];
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
    for(int iR=0; iR<yvc.n3(); iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*Rvc[iR][l];
        kR *= 2*pi;
        value+= exp(c1*kR)* yvc[ic][jc][iR];
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
    for(int iR=0; iR<Rvvz.n1(); iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*Rvvz[iR][l];
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
    for(int iR=0; iR<zvc.n3(); iR++)
    {
        double kR=0.;
        for(int l=0; l<3; l++) kR += k.crys[l]*Rvc[iR][l];
        kR *= 2*pi;
        value+= exp(c1*kR)* zvc[ic][jc][iR];
    }//loop R
    return value;
  }//loop R
  else return conj(Wannier::dipolez(jc,ic,k));
}





void Wannier::energy_U(vec2x& H, vec2x& Uk, Coord_B& k)    
{
    Coord_B ktemp;
    ktemp.set_crys(k.crys[2],k.crys[0],k.crys[1]);
    for(int ic = 0; ic < _Ncv; ic++)
    {
        for(int jc = 0; jc<_Ncv; jc++)
        {
            //H[ic][jc] = energy(ic, jc, k);
            H[ic][jc] = energy(ic, jc, ktemp);
            if(ic==jc) 
            {
		          H[ic][jc].imag(0.);
            }      
        }
    }
    //vector<int> ix(3);
    //for(int icoor=0; icoor<3; icoor++)
    //{
    //  ix[icoor] = (int)floor(fmod((k.crys[icoor])*Nk[icoor],Nk[icoor]));
    //  while(ix[icoor] < 0) ix[icoor]+=1;
    //}  
    Uk.fill(0.);
    //if(!Singularity[ix[2]+Nk[2]*ix[1]+Nk[2]*Nk[1]*ix[0]]) 
    ////if(Singularity[ix[2]+Nk[2]*ix[1]+Nk[2]*Nk[1]*ix[0]] || !Singularity[ix[2]+Nk[2]*ix[1]+Nk[2]*Nk[1]*ix[0]]) 
    //{
    //cout << " IK! " << ix[2]+Nk[2]*ix[1]+Nk[2]*Nk[1]*ix[0] << endl;
    //  for(int ic=0; ic<_Nch; ic++)
    //    Uk[ic][ic] = 1.;
    //  
    //  for(int ic=_Nch; ic<_Ncv; ic++)
    //  {
    //    for(int jc=_Nch; jc<_Ncv; jc++)
    //    {
    //      double kR = 0;
    //      for(int iR=0; iR<RvvU.n1();iR++)
    //      {
    //        for(int icoor=0; icoor<3; icoor++)
    //        //  kR += k.crys[icoor]*RvvU[iR][icoor];
    //        //Uk[ic][jc] += exp(2.*pi*c1*kR)*Uvv[ic-_Nch][jc-_Nch][iR];
    //        kR += k.crys[icoor]*RvvU[0][icoor];
    //        Uk[ic][jc] += exp(2.*pi*c1*kR)*Uvv[ic-_Nch][jc-_Nch][0];
    //      }
//
    //    }
    //  }
    //  return;
    //}

    vec epsilon;
    cx_mat U;
    cx_mat Hw; 
    epsilon.zeros(_Ncv-_Nch);
    U.zeros(_Ncv-_Nch, _Ncv-_Nch);
    Hw.zeros(_Ncv-_Nch, _Ncv-_Nch);
    //copy the Ek matrix in an armadillo matrix
    for (int ic=0; ic<_Ncv-_Nch; ic++)
    {
        for (int jc=0; jc < _Ncv-_Nch; jc++)
        {
            Hw(ic, jc) = H[ic+_Nch][jc+_Nch];
        }
    }
    eig_sym(epsilon, U, Hw);

  Uk.fill(0.);
  for(int ic=0; ic<_Nch; ic++)
  {
      Uk[ic][ic] = 1.;
  }

  for (int ic=_Nch; ic<_Ncv; ic++)
  {
      for (int jc=_Nch; jc < _Ncv; jc++)
      {
          complexd el = U(jc-_Nch,ic-_Nch);
          Uk[ic][jc]=conj(el);
      }    
  } 


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
    Coord_B ktemp;
    ktemp.set_crys(k.crys[2],k.crys[0],k.crys[1]);
    for(int ic = 0; ic < _Ncv; ic++)
    {
        for(int jc = 0; jc<_Ncv; jc++)
        {
            Dx[ic][jc][0] = dipolex(ic, jc, ktemp);
            Dx[ic][jc][1] = dipoley(ic, jc, ktemp);
            Dx[ic][jc][2] = dipolez(ic, jc, ktemp);
            //Dx[ic][jc][0] = dipolex(ic, jc, k);
            //Dx[ic][jc][1] = dipoley(ic, jc, k);
            //Dx[ic][jc][2] = dipolez(ic, jc, k);
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


void Wannier::print_dipole()
{
  for(int iR=0; iR<nrpts; iR++)
  {
    for(int ic=0; ic < _Ncv; ic++)
    {
      for(int jc=ic; jc< _Ncv; jc++)
      {
        cout << ic << " " << jc << " " << xvv[ic][jc][iR] << " " << xvv[jc][ic][nrpts-iR-1] << " " << (xvv[ic][jc][iR] - conj(xvv[jc][ic][nrpts-iR-1])) << " " << (xvv[ic][jc][iR] - conj(xvv[jc][ic][nrpts-iR-1]))/xvv[ic][jc][iR] << endl; 
        cout << ic << " " << jc << " " << yvv[ic][jc][iR] << " " << yvv[jc][ic][nrpts-iR-1] << " " << (yvv[ic][jc][iR] - conj(yvv[jc][ic][nrpts-iR-1])) << " " << (yvv[ic][jc][iR] - conj(yvv[jc][ic][nrpts-iR-1]))/yvv[ic][jc][iR] << endl; 
        cout << ic << " " << jc << " " << zvv[ic][jc][iR] << " " << zvv[jc][ic][nrpts-iR-1] << " " << (zvv[ic][jc][iR] - conj(zvv[jc][ic][nrpts-iR-1])) << " " << (zvv[ic][jc][iR] - conj(zvv[jc][ic][nrpts-iR-1]))/zvv[ic][jc][iR] << endl; 
      }
    }
  }
}


void Wannier::energy_Bloch(vec1d& energy, Coord_B& k)
{
    vec2x H(nwannier, nwannier);


    for(int ic = _Nch; ic < _Ncv; ic++)
    {
        for(int jc = _Nch; jc<_Ncv; jc++)
        {
            H[ic-_Nch][jc-_Nch] = this->energy(ic, jc, k);
            if(ic==jc) 
              {
                H[ic][jc].imag(0.);
              }      
        }
    }

    vec epsilon;
    cx_mat U;
    cx_mat Hw; 
    epsilon.zeros(nwannier);
    U.zeros(nwannier, nwannier);
    Hw.zeros(nwannier, nwannier);
    //copy the Ek matrix in an armadillo matrix
    for (int ic=0; ic<nwannier; ic++)
    {
        for (int jc=0; jc < nwannier; jc++)
        {
            Hw(ic, jc) = H[ic][jc];
        }
    }
    eig_sym(epsilon, U, Hw);


    for(int ic=0; ic<nwannier; ic++)
      energy[ic] = epsilon(ic);
}
