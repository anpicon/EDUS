#include "Coordinate.h"
#include <exception>
class TightBinding 
{
  public:
    TightBinding(){};
    void init(string name, int& Nch, int& Ncv);
    //complexd operator()(int& ic, int&jc, double& kx, double& ky, double& kz);
    void set_hopping(vec1d& hopping);
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
    complexd unitary(int& ic, int& jc, Coord_B& k);
    void GradientEnergy(vec3x& GE, Coord_B& k);
    string _name;

  private:
    int _Ncv;
    int _Nch;
    double a; // Lattice constant
    double b; // hopping parameter in plane
    double c; // distance between  lattices 
    double d; // hopping parameter between planes
    double g; // band gap
    double f;
    double t2;   // second-order hopping parameter
    double phi0; // phase of the second-order hopping parameter
    double core; // energy of the core band
    double rchv;  // dipole between core and valence
    double rchc;  // dipole between core and conduction
    double Romb_dist; // rombic distortion rate
    vec3d coeff;
};

void TightBinding::init(string name,int& Nch, int& Ncv)
{
    _Ncv = Ncv;
    _Nch = Nch;
    _name = name;
  //printf("name :  %s\n", _name.c_str());
    if (_name=="Graphene")
    {
        a=2.46*space_A_au;//Parameter in Angstroms
        b=2.97*energy_eV_au;//Parameter in eV
        c=0;
        d=0*energy_eV_au;
        g=0*energy_eV_au;//Bandgap in eV
    }
    else if(_name == "Graphene_costD" || _name == "Graphene_symmetricD")
    {
        a=2.46*space_A_au;//Parameter in Angstroms
        b=2.97*energy_eV_au;//Parameter in eV
        c=0;
        d=0*energy_eV_au;
        g=0*energy_eV_au;//Bandgap in eV
    }
    else if(_name == "CoreGraphene")
    {
        a=2.46*space_A_au;//Parameter in Angstroms
        b=2.97*energy_eV_au;//Parameter in eV
        c=0;
        d=0*energy_eV_au;
        g=0*energy_eV_au;//Bandgap in eV
    }
    else if (_name == "GrapheneZurron")
    {
        a=2.46*space_A_au;//Parameter in Angstroms
        b=2.97*energy_eV_au;//Parameter in eV
        c=0;
        d=0*energy_eV_au;
        g=0*energy_eV_au;//Bandgap in eV
    }
    else if (_name == "BoronNitride" || _name == "BN1b" || _name == "CoreBN_Nedge" || _name == "CoreBN_Bedge" || _name == "Haldane_CoreBN_Nedge" || _name == "Haldane_CoreBN_Bedge" || _name == "gen2d_hexagonal")
    {
        a=2.5*space_A_au;//Parameter in Angstroms
        b=3.16*energy_eV_au;//Parameter in eV
        c=0;
        d=0*energy_eV_au;
        g=6.2*energy_eV_au;//Bandgap in eV
        Romb_dist = 0.0; // 
    }
    else if (_name=="DLGraphene")
    {
        a=2.46*space_A_au;//Parameter in Angstroms
        b=2.97*energy_eV_au;//Parameter in eV
        c=3.35*space_A_au; // Parameter in Angstroms
        d=0.39*energy_eV_au; //Parameter in eV
        g=0*energy_eV_au;//Bandgap in eV
    }
    else if(_name == "Vampa2015")
    {
      coeff.resize(2,3,6); //first -> val, cond second -> x,y,z     third -> numbers

      coeff[0][0][0] = -0.0928;     coeff[1][0][0] =  0.0898;
      coeff[0][0][1] =  0.0705;     coeff[1][0][1] = -0.0814;
      coeff[0][0][2] =  0.0200;     coeff[1][0][2] = -0.0024;
      coeff[0][0][3] = -0.0012;     coeff[1][0][3] = -0.0048;
      coeff[0][0][4] =  0.0029;     coeff[1][0][4] = -0.0003;
      coeff[0][0][5] =  0.0006;     coeff[1][0][5] = -0.0009;
      coeff[0][1][0] = -0.0307;     coeff[1][1][0] =  0.1147;
      coeff[0][1][1] =  0.0307;     coeff[1][1][1] = -0.1147;
      coeff[0][1][2] =  0.0000;     coeff[1][1][2] =  0.0000;
      coeff[0][1][3] =  0.0000;     coeff[1][1][3] =  0.0000;
      coeff[0][1][4] =  0.0000;     coeff[1][1][4] =  0.0000;
      coeff[0][1][5] =  0.0000;     coeff[1][1][5] =  0.0000;

      coeff[0][2][0] = -0.0059;     coeff[1][2][0] =  0.0435;
      coeff[0][2][1] =  0.0059;     coeff[1][2][1] = -0.0435;    
      coeff[0][2][2] =  0.0000;     coeff[1][2][2] =  0.0000;    
      coeff[0][2][3] =  0.0000;     coeff[1][2][3] =  0.0000;    
      coeff[0][2][4] =  0.0000;     coeff[1][2][4] =  0.0000;    
      coeff[0][2][5] =  0.0000;     coeff[1][2][5] =  0.0000;    

      g = 0.1213;
    }
    else if(_name=="GaSe_kp2" || _name=="GaSe_kp2_cutoff")
    {
        a=-1.40*energy_eV_au;
        b=22.44*energy_eV_au;
        c=-54.41*energy_eV_au;
        d=74.14*energy_eV_au;
        f=2.49*energy_eV_au;
        g=67.56*energy_eV_au;
    }
    else if(_name=="GeS") //  A. M. Cook, B. M Fregoso, F. De Juan, S. Coh, and J. E. Moore, Nature communications 8, 1 (2017).
    {
        g=0.41*energy_eV_au;  //gap
        f=-2.33*energy_eV_au; //t1 hopping
        t2=0.61*energy_eV_au; //t2 hopping
        b=0.13*energy_eV_au;  //t3 hopping
        c=0.07*energy_eV_au;  //tt1 hopping
        d=-0.09*energy_eV_au; //tt2 hopping
    }
    else if(_name=="GeS_HSE06") // Alejandro Parameters with larger gap ((hybrid) exchange-correlation functionals)
    {
        g  =  0.62  *energy_eV_au;  //gap
        f  =  2.91  *energy_eV_au; //t1 hopping
        t2 = -0.797 *energy_eV_au; //t2 hopping
        b  = -0.05  *energy_eV_au;  //t3 hopping
        c  = -0.14  *energy_eV_au;  //tt1 hopping
        d  =  0.106 *energy_eV_au; //tt2 hopping
    }
}

void TightBinding::set_hopping(vec1d& hopping)
{
    if (_name=="Haldane_CoreBN_Nedge" || _name == "Haldane_CoreBN_Bedge")
    {
        g    = hopping[0];
        b    = hopping[1];
        t2   = hopping[2];
        phi0 = hopping[3];
    }
    if (_name=="gen2d_hexagonal")
    {
        a    = hopping[0];
        g    = hopping[1];
        b    = hopping[2];
        t2   = hopping[3];
        phi0 = hopping[4];
        core = hopping[5];
        rchv = hopping[6];
        rchc = hopping[7];
    }
}

complexd TightBinding::energy(int& ic, int& jc, Coord_B& k)
{
	if(_name == "CoreBN_Nedge" || _name == "CoreBN_Bedge")
    {
    	if(ic == 0)
    	{
    		if(jc == 0 && _name == "CoreBN_Nedge")
    			return -15.06; //409.9 eV, Nitrogen edge 
    		else if (jc == 0 && _name == "CoreBN_Bedge")
    			return -6.91; //188 eV, Boron k edge
    		else return 0;
    	}
    	else if (ic==1)
    	{
    		if (jc == 1)
    			return -g/2.;
    		else if (jc == 2)
    		{
    			complexd f = 1.+exp(2.*pi*c1*k.crys[1])+exp(2.*pi*c1*k.crys[2]);
    			return b*conj(f);
    		}
    		else return 0.;
    	}
    	else if (ic==2)
    	{
    		if(jc==1)
    		{
    			complexd f = 1.+exp(2.*pi*c1*k.crys[1])+exp(2.*pi*c1*k.crys[2]);
    			return b*f;
    		}
    		else if (jc==2)
    			return g/2.;
    		else return 0.;
    	}
    	else return 0.;
    }
    else if (_name=="Graphene") 
    {
        if(ic == jc)
              return( 0. );
        else
        {
          complexd f = exp(c1*2.*pi*k.crys[1])+1.+exp(c1*2.*pi*k.crys[2]);
          //complexd f=exp(-c1*a*k.cart[1]/sqrt(3.))*(1.+2.*exp(c1*sqrt(3.)*a*k.cart[1]/2.)*cos(a*k.cart[2]/2.)); //we change kx->ky, ky->kz
          if(ic==0)
              return(conj(f)*b);
          else
              return(f*b);
        }
    }//End Graphene Tight binding model
    else if (_name=="Graphene_costD" || _name == "Graphene_symmetricD") 
    {
        if(ic==jc && ic < 2)   
        {
            return(-10.892178);  
        } 
        else if ( ic >= 2 && jc >= 2 && ic == jc )
        {
            double modf = abs(exp(c1*2.*pi*k.crys[1])+1.+exp(c1*2.*pi*k.crys[2]));

            if ( ic == 2 ) return (-b*modf );
            else             return(b*modf );        
        }
        else return (0.);
    }

    else if(_name == "CoreGraphene")
    {
        if(ic==jc && ic < 2)   
        {
            return(-10.892178);  
        } 
        else if ( ic >= 2 && jc >= 2 && ic != jc )
        {
            complexd f = exp(c1*2.*pi*k.crys[1])+1.+exp(c1*2.*pi*k.crys[2]);

            if ( ic == 2 ) return (conj(f)*b);
            else             return(f*b);        
        }
        else return (0.);
        
    }
    else if (_name == "GrapheneZurron")
    {
        if(ic == jc)
            return( 0. );
        else
        {

            complexd f=exp(-c1*a*k.cart[1]/sqrt(3.))*(1.+2.*exp(c1*sqrt(3.)*a*k.cart[1]/2.)*cos(a*k.cart[2]/2.)); 
             if(ic==0)
                 return(f*b);
             else
                 return(conj(f)*b);    
        }
    }
    else if (_name == "BoronNitride")
    {
      if(ic == jc)
      {
        if(ic == 0)
            return(-g/2);
        else if(ic == 1)
            return(g/2);
        else
            return(0);
      }
      else
      {
        complexd f = exp(c1*2.*pi*k.crys[1])+1.+exp(c1*2.*pi*k.crys[2]);
        if(ic==0)
            return(conj(f)*b);
        else
            return(f*b);    
      }
    }
    else if (_name=="DLGraphene")
    {
        if(ic == jc)
              return( 0. );
        else
        {
          complexd f = exp(c1*2.*pi*k.crys[1])+1.+exp(c1*2.*pi*k.crys[2]);
          //complexd  f=exp(-c1*a*k.cart[1]/sqrt(3.))*(1.+2.*exp(c1*sqrt(3.)*a*k.cart[1]/2.)*cos(a*k.cart[2]/2.)); //we change kx->ky,  ky->kz
          if(ic==0 && jc==1)
              return(conj(f)*b);
          else if(ic==1 && jc==0)
              return(f*b);
          else if(ic==2 && jc==3)
              return(conj(f)*b);
          else if(ic==3 && jc==2)
              return(f*b);
          else if(ic==3 && jc==0)
              return(d);
          else if(ic==0 && jc==3)
              return(d);
          else
              return(0);
        }
    }
    else if(_name=="BN1b")
    {
      complexd f = exp(c1*2.*pi*k.crys[1])+1.+exp(c1*2.*pi*k.crys[2]);
      return(0.5*sqrt(g*g+4*b*b*f*conj(f)).real());
    }
    else if(_name == "Vampa2015")
    {
      if(ic!=jc) return (0.);
      else 
      {
        complexd ener=0.;
        for(int x=0; x<3; x++)
        {
          for(int j=0; j< coeff.n3(); j++)
          {
            ener += coeff[ic][x][j]*cos(j*k.crys[x]);
          }
        }
        if(ic==1) ener +=g;
        return(ener);
      }
    }
    else if(_name == "GaSe_kp2")
    {
        if(ic==0 && jc==0){
            return(a+b*(k.cart[1]*k.cart[1]+k.cart[2]*k.cart[2]));
        }
        else if((ic==0 && jc==2) || (ic==2 && jc==0)) {
            return(c*(k.cart[1]*k.cart[1]+k.cart[2]*k.cart[2]));
        }
        else if(ic==1 && jc==1) {
            return(-a + d*(k.cart[1]*k.cart[1]+k.cart[2]*k.cart[2]));
        }
        else if(ic==2 && jc==2) {
            return(f + g*(k.cart[1]*k.cart[1]+k.cart[2]*k.cart[2]));
        }
        else return (0.);
    }
    else if(_name == "GaSe_kp2_cutoff")
    {
        double cutoff=0.3*0.3;
        if (k.cart[1]*k.cart[1]+k.cart[2]*k.cart[2]<=cutoff) {
            if(ic==0 && jc==0){
                return(a+b*(k.cart[1]*k.cart[1]+k.cart[2]*k.cart[2]));
            }
            else if((ic==0 && jc==2) || (ic==2 && jc==0)) {
                return(c*(k.cart[1]*k.cart[1]+k.cart[2]*k.cart[2]));
            }
            else if(ic==1 && jc==1) {
                return(-a + d*(k.cart[1]*k.cart[1]+k.cart[2]*k.cart[2]));
            }
            else if(ic==2 && jc==2) {
                return(f + g*(k.cart[1]*k.cart[1]+k.cart[2]*k.cart[2]));
            }
            else return (0.);
        }
        else {
            if(ic==0 && jc==0){
                return(a+b*(cutoff));
            }
            else if((ic==0 && jc==2) || (ic==2 && jc==0)) {
                return(c*(cutoff));
            }
            else if(ic==1 && jc==1) {
                return(-a + d*(cutoff));
            }
            else if(ic==2 && jc==2) {
                return(f + g*(cutoff));
            }
            else return (0.);
        }
    }
    if(_name == "Haldane_CoreBN_Nedge" || _name == "Haldane_CoreBN_Bedge" || _name == "gen2d_hexagonal")
    {
        if(ic == 0)
        {
            if(jc == 0 && _name == "Haldane_CoreBN_Nedge")
                return -15.06; //409.9 eV, Nitrogen edge
            else if (jc == 0 && _name == "Haldane_CoreBN_Bedge")
                return -6.91; //188 eV, Boron k edge
            else if (jc == 0 && _name == "gen2d_hexagonal")
                return core; // specified in the input
            else return 0;
        }
        else if (ic==1)
        {
            if (jc == 1)
            {   // MINUS sign
                complexd f2 = exp(-2.*pi*c1*k.crys[2]+2.*pi*c1*k.crys[1])+exp(2.*pi*c1*k.crys[2])+exp(-2.*pi*c1*k.crys[1]);
                // Plus displacement lattice
                // complexd f2 = exp(2.*pi*c1*k.crys[2]-2.*pi*c1*k.crys[1])+exp(-2.*pi*c1*k.crys[2])+exp(2.*pi*c1*k.crys[1]);
                return -g/2.+2.*real(t2*exp(c1*phi0)*f2);
            }
            else if (jc == 2)
            {
                // MINUS DISPLACEMENT SIGN
                complexd f = 1.+exp(-2.*pi*c1*k.crys[1])+exp(-2.*pi*c1*k.crys[2]);
                // // PLUS DISPLACEMENT SIGN
                // complexd f = 1.+exp(2.*pi*c1*k.crys[1])+exp(2.*pi*c1*k.crys[2]);
                return b*(f);
            }
            else return 0.;
        }
        else if (ic==2)
        {
            if(jc==1)
            {
                // MINUS DISPLACEMENT SIGN
                complexd f = 1.+exp(-2.*pi*c1*k.crys[1])+exp(-2.*pi*c1*k.crys[2]);
                // // PLUS DISPLACEMENT SIGN
                // complexd f = 1.+exp(2.*pi*c1*k.crys[1])+exp(2.*pi*c1*k.crys[2]);
                return b*conj(f);
            }
            else if (jc==2)
            {
                // MINUS DISPLACEMENT SIGN
                complexd f2 = exp(-2.*pi*c1*k.crys[2]+2.*pi*c1*k.crys[1])+exp(2.*pi*c1*k.crys[2])+exp(-2.*pi*c1*k.crys[1]);

                // PLUS DISPLACEMENT SIGN
                // complexd f2 = exp(2.*pi*c1*k.crys[2]-2.*pi*c1*k.crys[1])+exp(-2.*pi*c1*k.crys[2])+exp(2.*pi*c1*k.crys[1]);
                return g/2.+2.*real(t2*exp(-c1*phi0)*f2);
            }
            else return 0.;
        }
        else return 0.;
    }
    else if (_name=="GeS" || _name=="GeS_HSE06")
    {
        if(ic == jc)
        {
            double tor = -2.*c*(cos(2.*pi*k.crys[1])+cos(2.*pi*k.crys[2])) -2.*d*(cos(2.*pi*(k.crys[1]-k.crys[2])));
            if(ic==0)
                return(g+tor);
            else
                return(-g+tor);
        }
        else
        {
            complexd phix = exp(c1*2.*pi*k.crys[1])+exp(c1*2.*pi*k.crys[2]);
            complexd tor  = -f-t2*phix-b*conj(phix);
          if(ic==0)
              return(conj(tor));
          else
              return(tor);
        }
    }//End GeS tight-binding
    
    else      return(0.);
}

complexd TightBinding::Xgradientenergy(int&ic, int&jc, Coord_B& k)
{
  if(_name == "Vampa2015")
  {
    if(ic!=jc) return (0.);
    else 
    {
      complexd gradener=0.;
      for(int j=0; j< coeff.n3(); j++)
      {
          gradener -= j*coeff[ic][0][j]*sin(j*k.crys[0]);
      }
      return(gradener);
    }
  }
  else                    return(0.);
}

complexd TightBinding::Ygradientenergy(int&ic, int&jc, Coord_B& k)
{
  if (_name=="Graphene") 
  {
    if(ic == jc)    return( 0. );
    else
    {
      complexd f = c1*2.*pi*exp(c1*2.*pi*k.crys[1]);
      if(ic==0)     return(conj(f)*b);
      else          return(f*b);
    }
  }
  else if(_name =="CoreGraphene")
  {
        if (ic == 2 && jc == 3)
        {
            complexd f = -c1*2.*pi*exp(-c1*2.*pi*k.crys[1]);
            return (f*b);
        }
        else if( ic == 3 && jc == 2)
        {
            complexd f = c1*2.*pi*exp(c1*2.*pi*k.crys[1]);
            return (f*b);
        }
        else return (0.);
  }
  else if (_name=="Graphene_costD") 
  {return 0.;
  }

  else if (_name == "GrapheneZurron")
  {
        if (ic==jc) return(0.);
        else
        {
            complexd f = -c1*2.*pi*exp(-c1*2.*pi*(k.crys[1] + k.crys[2])) + 2*pi*c1*exp(c1*pi*(k.crys[1]+k.crys[2]))*cos(3.*pi*(k.crys[1]-k.crys[2])) - 6.*pi*exp(c1*pi*(k.crys[1]+k.crys[2]))*sin(3.*pi*(k.crys[1]-k.crys[2]));//-c1*2.*pi/3.* exp(-c1*2.*pi/3.*(k.crys[1]+k.crys[2])) +c1*2.*pi/3.* exp(c1*pi/3.*(k.crys[1]+k.crys[2]))*(cos(pi*(k.crys[1]-k.crys[2])) -2*pi*exp(c1*pi/3.*(k.crys[1]+k.crys[2]))*sin(pi*(k.crys[1]-k.crys[2])));
            if (ic == 0) return (f*b);
            else         return(conj(f)*b);
        }
  }
  else if (_name=="BoronNitride") 
  {
    if(ic == jc)   return( 0. );
    else
    {
      complexd gf = c1*2.*pi*exp(c1*2.*pi*k.crys[1]);
      if(ic==0)    return(b*conj(gf));
      else         return(b*gf);
    }
  }
  else if (_name=="DLGraphene")
  {
    if(ic == jc)   return( 0. );
    else
    {
      complexd gf = c1*2.*pi*exp(c1*2.*pi*k.crys[1]);
      if(ic==0 && jc==1)
                    return(conj(gf)*b);
      else if(ic==1 && jc==0)
                    return(gf*b);
      else if(ic==2 && jc==3)
                    return(conj(gf)*b);
      else if(ic==3 && jc==2)
                    return(gf*b);
      else          return(0);
    }
  }
  else if(_name=="BN1b")
  {
    complexd f = exp(c1*2.*pi*k.crys[1])+1.+exp(c1*2.*pi*k.crys[2]);
    double ener = 0.5*sqrt(g*g+4*b*b*f*conj(f)).real();
    return 4*b*b*pi/ener*(exp(-c1*2.*pi*k.crys[1])*f).imag();
  }
  else if(_name == "Vampa2015")
  {
    if(ic!=jc) return (0.);
    else 
    {
      complexd gradener=0.;
      for(int j=0; j< coeff.n3(); j++)
      {
          gradener -= j*coeff[ic][1][j]*sin(j*k.crys[1]);
      }
      return(gradener);
    }
  }
  else if(_name == "GaSe_kp2")
  {
      if(ic==0 && jc==0){
          return(2.*b*(k.cart[1]));
      }
      else if((ic==0 && jc==2) || (ic==2 && jc==0)) {
          return(2.*c*(k.cart[1]));
      }
      else if(ic==1 && jc==1) {
          return(2.*d*(k.cart[1]));
      }
      else if(ic==2 && jc==2) {
          return(2.*g*(k.cart[1]));
      }
      else return (0.);
  }
  else if(_name == "GaSe_kp2_cutoff")
  {
      double cutoff=0.3*0.3;
      if(k.cart[1]*k.cart[1]+k.cart[2]*k.cart[2]<=cutoff){
          if(ic==0 && jc==0){
              return(2.*b*(k.cart[1]));
          }
          else if((ic==0 && jc==2) || (ic==2 && jc==0)) {
              return(2.*c*(k.cart[1]));
          }
          else if(ic==1 && jc==1) {
              return(2.*d*(k.cart[1]));
          }
          else if(ic==2 && jc==2) {
              return(2.*g*(k.cart[1]));
          }
          else return (0.);
      }
      else{
          return (0.);
      }
  }

  else              return(0.);
}

complexd TightBinding::Zgradientenergy(int&ic, int&jc, Coord_B& k)
{
  if (_name=="Graphene") 
  {
    if(ic == jc)    return( 0. );
    else
    {
      complexd f = c1*2.*pi*exp(c1*2.*pi*k.crys[2]);
      if(ic==0)     return(conj(f)*b);
      else          return(f*b);
    }
  }
  else if(_name =="CoreGraphene")
  {
        if (ic == 2 && jc == 3)
        {
            complexd f = -c1*2.*pi*exp(-c1*2.*pi*k.crys[2]);
            return (f*b);
        }
        else if( ic == 3 && jc == 2)
        {
            complexd f = c1*2.*pi*exp(c1*2.*pi*k.crys[2]);
            return (f*b);
        }
        else return (0.);
    }
      else if (_name=="Graphene_costD") 
    { return 0.;
    }

  else if (_name == "GrapheneZurron")
  {
          if (ic==jc) return(0.);
          else
          {
            complexd f = -c1*2.*pi*exp(-c1*2.*pi*(k.crys[1] + k.crys[2])) + 2*pi*c1*exp(c1*pi*(k.crys[1]+k.crys[2]))*cos(3.*pi*(k.crys[1]-k.crys[2])) + 6.*pi*exp(c1*pi*(k.crys[1]+k.crys[2]))*sin(3.*pi*(k.crys[1]-k.crys[2]));//        complexd f = -c1*2.*pi/3.* exp(-c1*2.*pi/3.*(k.crys[1]+k.crys[2])) +c1*2.*pi/3.* exp(c1*pi/3.*(k.crys[1]+k.crys[2]))*(cos(pi*(k.crys[1]-k.crys[2])) +2*pi*exp(c1*pi/3.*(k.crys[1]+k.crys[2]))*sin(pi*(k.crys[1]-k.crys[2])));
              if (ic == 0) return (f*b);
              else         return(conj(f)*b);
          }
  }
  else if (_name=="BoronNitride") 
  {
    if(ic == jc)    return( 0. );
    else
    {
      complexd gf = c1*2.*pi*exp(c1*2.*pi*k.crys[2]);
      if(ic==0)     return(b*conj(gf));
      else          return(b*gf);
    }
  }
  else if (_name=="DLGraphene")
  {
    if(ic == jc)   return( 0. );
    else
    {
      complexd gf = c1*2.*pi*exp(c1*2.*pi*k.crys[2]);
      if(ic==0 && jc==1)
                    return(conj(gf)*b);
      else if(ic==1 && jc==0)
                    return(gf*b);
      else if(ic==2 && jc==3)
                    return(conj(gf)*b);
      else if(ic==3 && jc==2)
                    return(gf*b);
      else          return(0);
    }
  }
  else if(_name=="BN1b")
  {
    complexd f = exp(c1*2.*pi*k.crys[1])+1.+exp(c1*2.*pi*k.crys[2]);
    double ener = 0.5*sqrt(g*g+4*b*b*f*conj(f)).real();
    return 4*b*b*pi/ener*(exp(-c1*2.*pi*k.crys[2])*f).imag();
  }
  else if(_name == "Vampa2015")
  {
    if(ic!=jc) return (0.);
    else 
    {
      complexd gradener=0.;
      for(int j=0; j< coeff.n3(); j++)
      {
          gradener -= j*coeff[ic][2][j]*sin(j*k.crys[2]);
      }
      return(gradener);
    }
  }
  else if(_name == "GaSe_kp2")
  {
      if(ic==0 && jc==0){
          return(2.*b*(k.cart[2]));
      }
      else if((ic==0 && jc==2) || (ic==2 && jc==0)) {
          return(2.*c*(k.cart[2]));
      }
      else if(ic==1 && jc==1) {
          return(2.*d*(k.cart[2]));
      }
      else if(ic==2 && jc==2) {
          return(2.*g*(k.cart[2]));
      }
      else return (0.);
  }
  else if(_name == "GaSe_kp2_cutoff")
  {
      double cutoff=0.3*0.3;
      if(k.cart[1]*k.cart[1]+k.cart[2]*k.cart[2]<=cutoff){
          if(ic==0 && jc==0){
              return(2.*b*(k.cart[2]));
          }
          else if((ic==0 && jc==2) || (ic==2 && jc==0)) {
              return(2.*c*(k.cart[2]));
          }
          else if(ic==1 && jc==1) {
              return(2.*d*(k.cart[2]));
          }
          else if(ic==2 && jc==2) {
              return(2.*g*(k.cart[2]));
          }
          else return (0.);
      }
      else return (0.);
  }

  else              return(0.);
}

double TightBinding::energy_Bloch(int& ic, int& jc, Coord_B& k)
{
	if(_name == "CoreBN_Nedge" || _name == "CoreBN_Bedge")
    {
    	if(ic==jc)
    	{
    		if(ic == 0 && _name == "CoreBN_Nedge")
    			return -15.06;
    		else if (ic == 0 && _name == "CoreBN_Bedge")
    			return -6.91;
    		else if (ic==1)
    		{
          complexd f = 1.+exp(2.*pi*c1*k.crys[1])+exp(2.*pi*c1*k.crys[2]);
    			return -sqrt(g*g+4.*b*b*f*conj(f)).real()/2.;
    		}
    		else if (ic==2)
    		{
          complexd f = 1.+exp(2.*pi*c1*k.crys[1])+exp(2.*pi*c1*k.crys[2]);
    			return sqrt(g*g+4.*b*b*f*conj(f)).real()/2.;
    		}
    		else return 0.;
    	}
    	else return 0.;
    }
    else if (_name == "BoronNitride")
    {
        if(ic != jc)    return( 0. );
        else
        {
          complexd f = exp(c1*2.*pi*k.crys[1])+1.+exp(c1*2.*pi*k.crys[2]);
          if(ic==0)
                        return(-0.5*sqrt(g*g+4*b*b*f*conj(f)).real());
          else
                        return(0.5*sqrt(g*g+4*b*b*f*conj(f)).real());
        }
    }
    else if(_name == "Haldane_CoreBN_Nedge" || _name == "Haldane_CoreBN_Bedge" || _name == "gen2d_hexagonal")
    {
        if(ic==jc)
        {
            if(ic == 0 && _name == "Haldane_CoreBN_Nedge")
                return -15.06;
            else if (ic == 0 && _name == "Haldane_CoreBN_Bedge")
                return -6.91;
            else if (ic == 0 && _name == "gen2d_hexagonal")
                return core;
            else if (ic==1 || ic==2)
            {
                double B0 = 2.*t2*cos(phi0)*(cos(2.*pi*k.crys[2]-2.*pi*k.crys[1])+cos(-2.*pi*k.crys[2])+cos(2.*pi*k.crys[1]));
                double Bx = b*(1.+cos(2.*pi*k.crys[1])+cos(2.*pi*k.crys[2]));
                double By = -b*(sin(2.*pi*k.crys[1])+sin(2.*pi*k.crys[2])); // MINUS! NEGATIVE DISPLACEMENT
                // MINUS DISPLACEMENT
                double Bz = -g/2. -2.*t2*sin(phi0)*(sin(-2.*pi*k.crys[2]+2.*pi*k.crys[1])+sin(2.*pi*k.crys[2])+sin(-2.*pi*k.crys[1]));

                // // PLUS DISPLACEMENT
                // double Bz = -g/2. -2.*t2*sin(phi0)*(sin(2.*pi*k.crys[2]-2.*pi*k.crys[1])+sin(-2.*pi*k.crys[2])+sin(2.*pi*k.crys[1]));
                double sign=(ic==1 ? -1.: 1.);
                return B0 + sign*sqrt(Bx*Bx+By*By+Bz*Bz);
            }
            else return 0.;
        }
        else return 0.;
    }

  	else if (_name=="Graphene") 
  	{
  	  if(ic != jc)    return( 0. );
  	  else
  	  {
  	    complexd f = exp(c1*2.*pi*k.crys[1])+1.+exp(c1*2.*pi*k.crys[2]);
  	    //complexd f=exp(-c1*a*k.cart[1]/sqrt(3.))*(1.+2.*exp(c1*sqrt(3.)*a*k.cart[1]/2.)*cos(a*k.cart[2]/2.)); //we change kx->ky, ky->kz
  	    if(ic==0)
  	                  return(-b*abs(f));
  	    else
  	                  return(b*abs(f));
  	    }
  	}//End Graphene Tight binding model
  	  else if (_name == "CoreGraphene")
  	  {
  	      if (ic < 2)
  	      {
  	          if(ic == jc) return (-10.892178);
  	          else return (0.);
  	      }
  	      else if(ic != jc)
  	          return( 0. );
  	      else
  	      {
  	          complexd f = exp(c1*2.*pi*k.crys[1])+1.+exp(c1*2.*pi*k.crys[2]);
  	          //complexd f=exp(-c1*a*k.cart[1]/sqrt(3.))*(1.+2.*exp(c1*sqrt(3.)*a*k.cart[1]/2.)*cos(a*k.cart[2]/2.)); //we change kx->ky, ky->kz
  	          if(ic==2)
  	              return(-b*abs(f));
  	          else
  	              return(b*abs(f));
  	      }
  	  }
  	  else if (_name=="Graphene_costD" || _name== "Graphene_symmetricD") 
  	  {
  	      return(real(TightBinding::energy(ic, jc, k)));
  	  }
	
  	  else if(_name == "GrapheneZurron")
  	  {
  	          if (ic!=jc)
  	              return 0.;
  	          else{
  	          complexd f=exp(-c1*a*k.cart[1]/sqrt(3.))*(1.+2.*exp(c1*sqrt(3.)*a*k.cart[1]/2.)*cos(a*k.cart[2]/2.));          
  	          if(ic==0)
  	               return(-1.*b*abs(f));
  	           else
  	               return(b*abs(f));
  	       }
	
  	  }
  	else if (_name=="DLGraphene")
  	{
  	  if(ic != jc)    return( 0. );
  	  else
  	  {
  	    complexd f = exp(c1*2.*pi*k.crys[1])+1.+exp(c1*2.*pi*k.crys[2]);
  	    //complexd  f=exp(-c1*a*k.cart[1]/sqrt(3.))*(1.+2.*exp(c1*sqrt(3.)*a*k.cart[1]/2.)*cos(a*k.cart[2]/2.)); //we change kx->ky,  ky->kz
  	    if(ic==0)
  	                  return( -0.5*(d-sqrt(d*d+4*b*b*f*conj(f))).real() );
  	    else if (ic==1)
  	                  return( -0.5*(d+sqrt(d*d+4*b*b*f*conj(f))).real() );
  	    else if (ic==2)
  	                  return(  0.5*(d-sqrt(d*d+4*b*b*f*conj(f))).real() );
  	    else if (ic==3)
  	                  return(  0.5*(d+sqrt(d*d+4*b*b*f*conj(f))).real() );
  	    else          return 0.;
  	    }
  	}
  	else if(_name=="BN1b")
  	{
  	  complexd f = exp(c1*2.*pi*k.crys[1])+1.+exp(c1*2.*pi*k.crys[2]);
  	  return(0.5*sqrt(g*g+4*b*b*f*conj(f)).real());
  	}
  	else if(_name == "Vampa2015")
  	{
  	  
  	    double ener=0.;
  	    for(int x=0; x<3; x++)
  	    {
  	      for(int j=0; j< coeff.n3(); j++)
  	      {
  	        ener += coeff[ic][x][j]*cos(j*k.crys[x]);
  	      }
  	    }
  	    if(ic==1) ener +=g;
  	    return(ener);
  	  
  	}
  else    return(0);
}

complexd TightBinding::unitary(int& ic, int& jc, Coord_B& k)
{
	if(_name == "CoreBN_Nedge" || _name == "CoreBN_Bedge")
    {
    	if(ic == 0)
    	{
    		if(jc == 0)
    			return 1.; //409.9 eV, Nitrogen edge 
    		else return 0;
    	}
    	else if(jc == 0) return 0.;
      else
    	{    	
          complexd f = 1.+exp(2.*pi*c1*k.crys[1])+exp(2.*pi*c1*k.crys[2]);
    		double E = sqrt(g*g+4.*b*b*f*conj(f)).real()/2.;
    		double a = sqrt((g+2.*E)/(4.*E));
    		double h = sqrt((2.*E-g)/(4.*E));
    		if (ic == 1)
    		{
    			if (jc == 1)
        			return a;
    			else return -h*conj(f)/abs(f);
    		}
    		else if ( ic == 2 )
    		{
    			if ( jc == 1 )
    				return h;
    			else return a*conj(f)/abs(f);
    		}
    		else return 0.;
    	}
    }
    else if (_name == "BoronNitride")
    {
        double Bx = b*(1.+cos(2.*pi*k.crys[1])+cos(2.*pi*k.crys[2]));
        double By = b*(sin(2.*pi*k.crys[1])+sin(2.*pi*k.crys[2]));
        double Bz = -g/2.;
        double B  = sqrt(Bx*Bx+By*By+Bz*Bz);
        double phi    = atan2(By,Bx)/2.;
        double costh2 = sqrt((B+Bz)/2./B);
        double sinth2 = sqrt((B-Bz)/2./B);
        
        if (ic == 0)
        {
            if (jc == 0)
                return exp(+c1*phi)*sinth2;
            else return -exp(-c1*phi)*costh2;
        }
        else if ( ic == 1 )
        {
            if ( jc == 0 )
                return exp(+c1*phi)*costh2;
            else return exp(-c1*phi)*sinth2;
        }
        else return 0.;
        /* Old version
         complexd f = exp(c1*2.*pi*k.crys[1])+1.+exp(c1*2.*pi*k.crys[2]);
         double aa = sqrt((4.*b*b*f*conj(f)+g*g).real());
         if(ic==0)
         {
           double den = sqrt(4+pow(abs((g+aa)/(b*f)),2));
           if(jc==0)
           {
             return (-(g+aa)/(b*den*conj(f)));
           }
           else
           {
             return (2./den);
           }
         }
         else
         {
           double den = sqrt(4+pow(abs((g-aa)/(b*f)),2));
           
       
           if (jc == 0)    return ((-g+aa)/(b*den*conj(f)));
           else            return (2./den);
         }*/
    }
    else if(_name == "Haldane_CoreBN_Nedge" || _name == "Haldane_CoreBN_Bedge" || _name == "gen2d_hexagonal")
    {
        if(ic == 0)
        {
            if(jc == 0)
                return 1.; //for the core band
            else return 0;
        }
        else if(jc == 0) return 0.;
        else
        {
            //double B0 = 2.*t2*cos(phi0)*(cos(2.*pi*k.crys[2]-2.*pi*k.crys[1])+cos(-2.*pi*k.crys[2])+cos(2.*pi*k.crys[1]));
            double Bx = b*(1.+cos(2.*pi*k.crys[1])+cos(2.*pi*k.crys[2]));
            double By = b*(sin(2.*pi*k.crys[1])+sin(2.*pi*k.crys[2])); // SIGN MINUS!!!
            double Bz = -g/2. -2.*t2*sin(phi0)*(sin(-2.*pi*k.crys[2]+2.*pi*k.crys[1])+sin(+2.*pi*k.crys[2])+sin(-2.*pi*k.crys[1]));
            double B  = sqrt(Bx*Bx+By*By+Bz*Bz);
            double phi    = atan2(By,Bx)/2.;
            double costh2 = sqrt((B+Bz)/2./B);
            double sinth2 = sqrt((B-Bz)/2./B);
            if (ic == 1)
            {
                if (jc == 1)
                    return exp(+c1*phi)*sinth2;
                else return -exp(-c1*phi)*costh2;
            }
            else if ( ic == 2 )
            {
                if ( jc == 1 )
                    return exp(+c1*phi)*costh2;
                else return exp(-c1*phi)*sinth2;
            }
            else return 0.;
        }
    }
    else if(_name == "GeS" || _name=="GeS_HSE06")
    {
        complexd phix = exp(c1*2.*pi*k.crys[1])+exp(c1*2.*pi*k.crys[2]);
        complexd tor  = -f-t2*phix-b*conj(phix);
        
        double Bx = real(tor);
        double By = imag(tor);
        double Bz = g;
        double B  = sqrt(Bx*Bx+By*By+Bz*Bz);
        double phi    = atan2(By,Bx)/2.;
        double costh2 = sqrt((B+Bz)/2./B);
        double sinth2 = sqrt((B-Bz)/2./B);
        if (ic == 0)
        {
            if (jc == 0)
                return exp(+c1*phi)*sinth2;
            else return -exp(-c1*phi)*costh2;
        }
        else if ( ic == 1 )
        {
            if ( jc == 0 )
                return exp(+c1*phi)*costh2;
            else return exp(-c1*phi)*sinth2;
        }
        else return 0.;
    }

  	else if (_name=="Graphene") 
  	{
  	  if(jc==1)
  	  {
  	    complexd f = exp(-c1*2.*pi*k.crys[1])+1.+exp(-c1*2.*pi*k.crys[2]);
  	    //complexd f=exp(-c1*a*k.cart[1]/sqrt(3.))*(1.+2.*exp(c1*sqrt(3.)*a*k.cart[1]/2.)*cos(a*k.cart[2]/2.)); //we change kx->ky, ky->kz
  	    if (ic == 0)  return(-f/(sqrt(2.)*abs(f)));
  	    else          return (f/(sqrt(2.)*abs(f)));
  	  }
  	  else            return (1./sqrt(2.));  
  	}//End Graphene Tight binding model
	
  	else if (_name == "CoreGraphene")
  	{
  	  if ( ic < 2 )
  	  {
  	    if(ic == jc)     return (1.);
  	    else             return (0.);
  	  }
  	  else
  	  {
  	    if (jc ==3)
  	    {
  	      complexd f = exp(-c1*2.*pi*k.crys[1])+1.+exp(-c1*2.*pi*k.crys[2]);
  	      if (ic == 2)    return(-f/(sqrt(2.)*abs(f)));
  	      else            return (+f/(sqrt(2.)*abs(f)));
  	    }
  	    else if (jc == 2) return (1./sqrt(2.));
  	    else              return (0.);
  	  } 
  	}
  	  else if (_name=="Graphene_costD" || _name == "Graphene_symmetricD") 
  	{
  	    if(ic==jc)   
  	    {
  	        return(1.);  
  	    }  
  	    else return (0.);
  	}
	
  	else if (_name =="GrapheneZurron")
  	{
  	      if(jc==1)
  	      {
  	          complexd f=exp(-c1*a*k.cart[1]/sqrt(3.))*(1.+2.*exp(c1*sqrt(3.)*a*k.cart[1]/2.)*cos(a*k.cart[2]/2.));          
	
  	          if (ic == 0) return(-f/(sqrt(2.)*abs(f)));
  	          else return (f/(sqrt(2.)*abs(f)));
  	      }
  	      else return (1./sqrt(2.));
  	}
 else if(_name == "DLGraphene" )
 {	
  	  complexd f = exp(c1*2.*pi*k.crys[1])+1.+exp(c1*2.*pi*k.crys[2]); 
  	  if(ic == 0)
  	  {
  	    if(jc == 0)
  	      return   
  	      -2/sqrt(8 + pow(abs((d - sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))/(b*conj(f))),2) + 16*pow(abs((b*f)/(d + sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))),2));
  	    else if (jc == 1)
  	      return 
  	      (-d + conj(sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f))))/(b*f*sqrt(8 + pow(abs((d - sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))/(b*conj(f))),2) + 16*pow(abs((b*f)/(d + sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))),2)));
  	    else if (jc == 2)
  	      return 
  	      (-4*b*f)/(sqrt(8 + pow(abs((d - sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))/(b*conj(f))),2) + 16*pow(abs((b*f)/(d + sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))),2))*(d + conj(sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))));
  	    else 
  	      return
  	      2/sqrt(8 + pow(abs((d - sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))/(b*conj(f))),2) + 16*pow(abs((b*f)/(d + sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))),2));    
  	  }
  	  else if (ic == 1)
  	  {
  	    if(jc == 0)
  	      return  
  	      2/sqrt(8 + 16*pow(abs((b*f)/(d - sqrt(pow(d,2)  + 4*pow(b,2)*f*conj(f)))),2) + pow(abs((d + sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))/(b*conj(f))),2)); 
  	    else if (jc == 1)
  	      return 
  	      -((d + conj(sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f))))/(b*f*sqrt(8 + 16*pow(abs((b*f)/(d - sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))),2) + pow(abs((d + sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))/(b*conj(f))),2))));
  	    else if (jc == 2)
  	      return 
  	      (4*b*f)/(sqrt(8 + 16*pow(abs((b*f)/(d - sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))),2) + pow(abs((d + sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))/(b*conj(f))),2))*(d - conj(sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))));
  	    else 
  	      return
  	      2/sqrt(8 + 16*pow(abs((b*f)/(d - sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))),2) + pow(abs((d + sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))/(b*conj(f))),2));
  	  }
  	  else if (ic == 2)
  	  {
  	     if(jc == 0)
  	      return 
  	      -2/sqrt(8 + 16*pow(abs((b*f)/(d - sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))),2) + pow(abs((d + sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))/(b*conj(f))),2));  
  	    else if (jc == 1)
  	      return 
  	      -((d + conj(sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f))))/(b*f*sqrt(8 + 16*pow(abs((b*f)/(d - sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))),2) + pow(abs((d + sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))/(b*conj(f))),2))));
  	    else if (jc == 2)
  	      return 
  	      (-4*b*f)/(sqrt(8 + 16*pow(abs((b*f)/(d - sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))),2) + pow(abs((d + sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))/(b*conj(f))),2))*(d - conj(sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))));
  	    else 
  	      return
  	      2/sqrt(8 + 16*pow(abs((b*f)/(d - sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))),2) + pow(abs((d + sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))/(b*conj(f))),2));
  	  }
  	  else
  	  {
  	    if(jc == 0)
  	      return   
  	      2/sqrt(8 + pow(abs((d - sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))/(b*conj(f))),2) + 16*pow(abs((b*f)/(d + sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))),2));
  	    else if (jc == 1)
  	      return 
  	      (-d + conj(sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f))))/(b*f*sqrt(8 + pow(abs((d - sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))/(b*conj(f))),2) + 16*pow(abs((b*f)/(d + sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))),2)));
  	    else if (jc == 2)
  	      return 
  	      (4*b*f)/(sqrt(8 + pow(abs((d - sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))/(b*conj(f))),2) + 16*pow(abs((b*f)/(d + sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))),2))*(d + conj(sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))));
  	    else 
  	      return
  	      2/sqrt(8 + pow(abs((d - sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))/(b*conj(f))),2) + 16*pow(abs((b*f)/(d + sqrt(pow(d,2) + 4*pow(b,2)*f*conj(f)))),2));
  	  }
 }	
  	else if(_name=="BN1b")
  	{
  	  return 1;
  	}
  	else if(_name == "Vampa2015")
  	{
  	  if(ic == jc) return (1.);
  	  else return (0.);
  	}
  	else if (_name == "GaSe_kp2")
  	{
  	    double kx=k.cart[1];
  	    double ky=k.cart[2];
  	    if(ic == 0 && jc == 0) {
  	        complexd f = (3.*(2252.*kx*kx + 2252.*ky*ky + 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(45.*kx*kx - sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321)/200. + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	        complexd norm = (3.*(2252.*kx*kx + 2252.*ky*ky + 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(45.*kx*kx - sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321.)/200. + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	        norm*=norm; norm+=1.;
  	        return f/sqrt(norm);
  	    }
  	    else if(ic == 0 && jc == 2) {
  	        complexd f = 1.;
  	        complexd norm = (3.*(2252.*kx*kx + 2252.*ky*ky + 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(45.*kx*kx - sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321.)/200. + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	        norm*=norm; norm+=1.;
  	        return f/sqrt(norm);
  	    }
  	    else if(ic == 2 && jc == 0) {
  	        complexd f = (3.*(2252.*kx*kx + 2252.*ky*ky+ 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321.)/200. + 45.*kx*kx + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	        complexd norm = (3.*(2252.*kx*kx + 2252.*ky*ky + 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321.)/200. + 45.*kx*kx + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	        norm*=norm; norm+=1.;
  	        return f/sqrt(norm);
  	    }
  	    else if(ic == 2 && jc == 2) {
  	        complexd f = 1.;
  	        complexd norm = (3.*(2252.*kx*kx + 2252.*ky*ky + 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321.)/200. + 45.*kx*kx + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	        norm*=norm; norm+=1.;
  	        return f/sqrt(norm);
  	    }
  	    else if(ic == 1 && jc == 1) return 1.;
  	    else return (0.);
  	}
  	else if (_name == "GaSe_kp2_cutoff")
  	{
  	    double cutoff=0.3*0.3;
  	    double kx=k.cart[1];
  	    double ky=k.cart[2];
  	    if (kx*kx+ky*ky<=cutoff) {
  	        if(ic == 0 && jc == 0) {
  	            complexd f = (3.*(2252.*kx*kx + 2252.*ky*ky + 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(45.*kx*kx - sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321)/200. + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	            complexd norm = (3.*(2252.*kx*kx + 2252.*ky*ky + 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(45.*kx*kx - sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321.)/200. + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	            norm*=norm; norm+=1.;
  	            return f/sqrt(norm);
  	        }
  	        else if(ic == 0 && jc == 2) {
  	            complexd f = 1.;
  	            complexd norm = (3.*(2252.*kx*kx + 2252.*ky*ky + 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(45.*kx*kx - sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321.)/200. + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	            norm*=norm; norm+=1.;
  	            return f/sqrt(norm);
  	        }
  	        else if(ic == 2 && jc == 0) {
  	            complexd f = (3.*(2252.*kx*kx + 2252.*ky*ky+ 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321.)/200. + 45.*kx*kx + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	            complexd norm = (3.*(2252.*kx*kx + 2252.*ky*ky + 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321.)/200. + 45.*kx*kx + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	            norm*=norm; norm+=1.;
  	            return f/sqrt(norm);
  	        }
  	        else if(ic == 2 && jc == 2) {
  	            complexd f = 1.;
  	            complexd norm = (3.*(2252.*kx*kx + 2252.*ky*ky + 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321.)/200. + 45.*kx*kx + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	            norm*=norm; norm+=1.;
  	            return f/sqrt(norm);
  	        }
  	        else if(ic == 1 && jc == 1) return 1.;
  	        else return (0.);
  	    }
  	    else{
  	        kx=sqrt(cutoff);
  	        ky=0.;
  	        if(ic == 0 && jc == 0) {
  	            complexd f = (3.*(2252.*kx*kx + 2252.*ky*ky + 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(45.*kx*kx - sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321)/200. + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	            complexd norm = (3.*(2252.*kx*kx + 2252.*ky*ky + 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(45.*kx*kx - sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321.)/200. + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	            norm*=norm; norm+=1.;
  	            return f/sqrt(norm);
  	        }
  	        else if(ic == 0 && jc == 2) {
  	            complexd f = 1.;
  	            complexd norm = (3.*(2252.*kx*kx + 2252.*ky*ky + 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(45.*kx*kx - sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321.)/200. + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	            norm*=norm; norm+=1.;
  	            return f/sqrt(norm);
  	        }
  	        else if(ic == 2 && jc == 0) {
  	            complexd f = (3.*(2252.*kx*kx + 2252.*ky*ky+ 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321.)/200. + 45.*kx*kx + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	            complexd norm = (3.*(2252.*kx*kx + 2252.*ky*ky + 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321.)/200. + 45.*kx*kx + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	            norm*=norm; norm+=1.;
  	            return f/sqrt(norm);
  	        }
  	        else if(ic == 2 && jc == 2) {
  	            complexd f = 1.;
  	            complexd norm = (3.*(2252.*kx*kx + 2252.*ky*ky + 83.))/(5441.*(kx*kx + ky*ky)) - (100.*(sqrt(138776068.*kx*kx*kx*kx + 277552136.*kx*kx*ky*ky + 3510336.*kx*kx + 138776068.*ky*ky*ky*ky + 3510336.*ky*ky + 151321.)/200. + 45.*kx*kx + 45.*ky*ky + 109./200.))/(5441.*(kx*kx + ky*ky));
  	            norm*=norm; norm+=1.;
  	            return f/sqrt(norm);
  	        }
  	        else if(ic == 1 && jc == 1) return 1.;
  	        else return (0.);
  	    }
  	}


else return 0;
}

complexd TightBinding::dipolex(int& ic, int& jc, Coord_B& k)
{
   	if(_name == "gen2d_hexagonal")
    {
           //we have always the freedom to put the atom with the core hole
        //in 0 and the other one in -delta_3 = -a/sqrt(3)
        if ((ic == 0 && jc == 1) || (ic == 1 && jc == 0))
            return rchv;
        else if ((ic == 0 && jc == 2) || (ic == 2 && jc == 0))
            return rchc;
        else return 0.;
    }
    else if(_name == "CoreBN_Nedge" || _name == "Haldane_CoreBN_Nedge")
    {
	   	//we have always the freedom to put the atom with the core hole 
    	//in 0 and the other one in -delta_3 = -a/sqrt(3)
    	if ((ic == 0 && jc == 1) || (ic == 1 && jc == 0))
    		return 0.0697595094;
    	else return 0.;
    }
    else if (_name == "CoreBN_Bedge" || _name == "Haldane_CoreBN_Bedge")
    {
    	if ((ic == 0 && jc == 2) || (ic == 2 && jc == 0))
    		return 0.0869151662;
    	else return 0.;
    }
    else if(_name == "CoreGraphene")
    {
        if( (ic == 0 && jc == 2) || (ic == 1 && jc == 3) || (ic == 2 && jc == 0) || (ic == 3 && jc == 1))
                                    return (0.0780088 );
        else                        return (0.);
    }
    else if (_name=="Graphene_costD" || _name == "Graphene_symmetricD") 
    {
        if(ic <= 1 && jc > 1)
        {
            complexd dx = 0.0780088/sqrt(2.);
            if(ic == 0 || (ic == 1 && jc == 2)) return dx;
            else return -dx;  
        }
        else if (ic > 1 && jc <= 1)
          return conj(TightBinding::dipolex(jc, ic, k));
        else return 0.;
    }
    else if(_name == "DLGraphene")
    {
        if(ic==jc)
        {
          if (ic==2 || ic==3)       return(c);
          else                      return(0.);
        }
        else                        return(0.);
    }
    else if(_name == "Vampa2015")
    {
      if (ic!=jc) return (3.46);
      else return(0.);
    }
    else                            return(0.);
}

complexd TightBinding::dipoley(int& ic, int& jc, Coord_B& k)
{
    if(_name == "CoreBN_Nedge" || _name == "Haldane_CoreBN_Nedge" )
    {
	   	//we have always the freedom to put the atom with the core hole 
    	//in 0 and the other one in -delta_3 = -a/sqrt(3)
    	if (ic == 2 && jc == 2)
    		return a/sqrt(3.);
    	else return 0.;
    }
    else if (_name == "CoreBN_Bedge" || _name == "Haldane_CoreBN_Bedge")
    {
        //we have always the freedom to put the atom with the core hole
        //in 0 and the other one in -delta_3 = -a/sqrt(3)
    	if ((ic == 2 && jc == 2) || (ic == 0 && jc == 0) )
    		return a/sqrt(3.); // CAREFULL WITH SIGN!!!
    	else return 0.;
    }
    else if (_name == "gen2d_hexagonal")
    {
        //we have always the freedom to put the atom with the core hole
        //in 0 and the other one in -delta_3 = -a/sqrt(3)
        if (ic == 2 && jc == 2)
            return -a/sqrt(3.);
        else if ((ic == 0 && jc == 0) && abs(rchc) > 1.e-6 )
            return -a/sqrt(3.);
        else return 0.;
    }
    else if (_name == "GeS" || _name=="GeS_HSE06")
    {
        if (ic == 1 && jc == 1)
            return 0.52*space_A_au;
        else return 0.;
    }
  	else if (_name=="Graphene")
  	{
  	  if(ic == 1 && jc == 1 )         return(-2.68341111);
  	  else                            return(0.);
  	}//End Graphene Tight binding model
  	else if (_name=="Graphene_costD" || _name == "Graphene_symmetricD") 
  	{
  	    if (ic == 1 && jc == 1) return -2.68341111;
  	    else if (ic > 1 && jc > 1)
  	    {
  	        if (ic == jc) return (-2.68341111/2);//-a*.5/sqrt(3.));
  	        else return (2.68341111/2);//+a*.5/sqrt(3.));
  	    }
  	    else return(0.);
  	}
  	else if(_name == "CoreGraphene")
  	{
  	  if (ic == jc && ic % 2 == 1)    return(-2.68341111);
  	  else                            return 0.;
  	}
  	else if (_name=="BoronNitride") 
  	{
  	  if(ic == 1 && jc == 1 )         return(a/sqrt(3));
  	  else                            return(0.);
  	}//End Graphene Tight binding model
  	else if(_name == "DLGraphene")
  	{
  	    if(ic==jc)
  	    {
  	      if (ic==1)                return(-a/sqrt(3));
  	      else if (ic==2 || ic==3)  return(a/(2*sqrt(3)));
  	      else                      return(0.);
  	    }
  	    else                        return(0.);
  	}
  	else if(_name == "Vampa2015")
  	{
  	  if (ic!=jc) return (3.46);
  	  else return(0.);
  	}
	
  	else                              return(0);
}

complexd TightBinding::dipolez(int& ic, int& jc, Coord_B& k)
{
  if (_name=="Graphene")            return 0.;
  else if(_name == "CoreGraphene")  return(0.);
  else if(_name == "Graphene_costD")  return(0.);
  else if ( _name == "Graphene_symmetricD") 
  {
      return (TightBinding::dipoley(ic, jc, k));
  }
  else if(_name == "DLGraphene")
  {
    if(ic==jc)
    {
      if (ic==2)                    return(a/2.);
      else if(ic==3)                return(-a/2.);
      else                          return(0.);
    }
    else                            return(0.);
  }
  else if(_name == "Vampa2015")
  {
    if (ic!=jc) return (3.94);
    else return(0.);
  }

  else                              return(0);
}



void TightBinding::energy_U(vec2x& H, vec2x& Uk,Coord_B& k)
{
    for(int ic = 0; ic < _Ncv; ic++)
    {
        for(int jc = 0; jc<_Ncv; jc++)
        {
            H[ic][jc] = energy(ic, jc, k);
            Uk[ic][jc] = unitary(ic, jc, k);
        }
    }
}

void TightBinding::GradientEnergy(vec3x& GE, Coord_B& k)
{
    GE.fill(0.);
    for(int ic = _Nch; ic < _Ncv; ic++)
    {
        for(int jc = _Nch; jc<_Ncv; jc++)
        {
            GE[ic][jc][0] = Xgradientenergy(ic, jc, k);
            GE[ic][jc][1] = Ygradientenergy(ic, jc, k);
            GE[ic][jc][2] = Zgradientenergy(ic, jc, k);
        }
    }
}

void  TightBinding::energy(vec2x& H, Coord_B& k)
{
    for(int ic = 0; ic < _Ncv; ic++)
    {
        for(int jc = 0; jc<_Ncv; jc++)
        {
            H[ic][jc] = energy(ic, jc, k);
        }
    }
}

void  TightBinding::dipole(vec3x& Dx, Coord_B& k)
{
    for(int ic = 0; ic < _Ncv; ic++)
    {
        for(int jc = 0; jc<_Ncv; jc++)
        {
            Dx[ic][jc][0] = dipolex(ic, jc, k);
            Dx[ic][jc][1] = dipoley(ic, jc, k);
            Dx[ic][jc][2] = dipolez(ic, jc, k);
        }
    }
}


