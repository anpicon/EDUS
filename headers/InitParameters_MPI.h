void BZ_parameters_MPI(vec1i& Nk,vec1d& dk, Laser& pulse1, Laser& pulse2, vec2d& b,vec2d& kpt,int& nTAk, vec2d& TAkpt, vec1i& tagTAk, int& nkk, int& rank, int& num_procs)
{
    Coord_B::set_crys_to_cart(b);
    pulse1.set_pol(pulse1.pol1.cart[0], pulse1.pol1.cart[1], pulse1.pol1.cart[2], pulse1.pol2.cart[0], pulse1.pol2.cart[1], pulse1.pol2.cart[2]);
    pulse2.set_pol(pulse2.pol1.cart[0], pulse2.pol1.cart[1], pulse2.pol1.cart[2], pulse2.pol2.cart[0], pulse2.pol2.cart[1], pulse2.pol2.cart[2]);
    //defining the k-grid
    int nk = Nk[0]*Nk[1]*Nk[2];
    kpt.resize(nkk,3);
    for(int coor=0; coor<3; coor++)
        dk[coor] = 1./Nk[coor];

    ofstream fp_k;
    stringstream sname;
    sname.seekp(0,ios::beg); sname << rank;
    string ssname= "kpt_" + sname.str() + ".txt";
    fp_k.open(ssname.c_str());
/*

        int Subset0=0;
        int Subset1=0;
    for (int i = 0; i < num_procs; i++)
    {
        if ( (num_procs-i-1) < nk%num_procs ) Subset1++;
        else Subset0++;
    }    cout << endl << "S0 " << Subset0 << "   S1 "<< Subset1 << endl;
    
    */for(int ik=0; ik<nkk; ik++)
    {
        int i;   
    /*    
        if(rank<Subset0){
          i = rank*nkk+ik;
        }
        else{
          i = (rank-Subset0)*nkk+ik+Subset0*(nkk-1);
        }

*/
        i = rank*nk/num_procs + ik + ( (num_procs-rank-1) < nk%num_procs )*(rank-num_procs+nk%num_procs);

        vec1i iik(3);
        iik[0] = (i/(Nk[2]*Nk[1]))%Nk[0];
        iik[1] = (i/Nk[2])%Nk[1];
        iik[2] = i%Nk[2];
        fp_k << i << " " << iik[0] << " " << iik[1] << " " << iik[2];
        for(int coor=0; coor<3; coor++)
        {
            kpt[ik][coor]=-0.5+dk[coor]/2.+dk[coor]*iik[coor];
            fp_k <<" "<< kpt[ik][coor];
        }
        fp_k << endl;
    }
    fp_k.close();
    
    //function to print the absorption at specific points of the k-grid
    /* //Difficult to implement in the MPI version, as we need to find the minimum distance for a given k point
    for(int I=0; I<nTAk; I++)
    {
        for(int coor = 0; coor < 3; coor++)
        {
            if(TAkpt[I][coor] < -0.5 || TAkpt[I][coor] > 0.5)
            {
                printf("Error. The coordinates of your k points have to stay in the range [-0.5 : 0.5]\n");
                exit(1);
            }
        }

        Coord_B k1;
        k1.setcrys(TAkpt[I][0], TAkpt[I][1], TAkpt[I][2]);
        double distance = 10.;
        double ik_right;
        for(int ik=0; ik<kpt.n1(); ik++)
        {
            Coord_B k;
            k.setcrys(kpt[ik][0], kpt[ik][1], kpt[ik][2]);
            double distance_ik = 0.;

            for(int ix=0; ix<3; ix++) 
            {
                double distance_ix = k1.cart[ix] - k.cart[ix];
                distance_ik += distance_ix*distance_ix;
            }

            if(distance_ik < distance)
            {
              ik_right = ik;
              distance = distance_ik;         
            }
        }//end ik

        for(int ix=0; ix<3; ix++)  TAkpt[I][ix] = kpt[ik_right][ix];
        tagTAk[I] = ik_right;

    }//end I
     */
}//end BZ_parameters








void Print_Input
   (string RefCoord,
    vec2d& a, vec2d& b, //unit cell
   vec1i& Nk, vec1d& dk,          //k space
   vec1i& Nb,  //bands
    //int& Nc, int& Nv, int& Nch,              //bands
   Laser& pulse1, Laser& pulse2, double& DELAY,//Coord_B& u1, Coord_B& u2, bool& gaussian1, bool& gaussian2, double& sigma1, double& sigma2, double& DELAY,//lasers
   string& iMode, string& TBtype, double& dt,  //tdse options
   bool& iWFDs, bool& iCurrent, bool& iTAbs, bool& iTAbsK,    //observables
   double& T1, double& T2, double& Tch,
   int& nTAk, vec2d& TAkpt, vec1i& tagTAk //x-ray absorption
   )
{

      printf("Calculation ready to start.\n\n");
          printf("\n*********************************************************************************************\n");
          printf(  "*                              %20s                                     *\n",RefCoord.c_str());
          printf(  "*********************************************************************************************\n");
          printf("\n*********************************************************************************************\n");
          printf(  "*                              Summary of parameters                                        *\n");
          printf(  "*********************************************************************************************\n");  
          printf(  "*  dt                      %2.4f fs   - %2.4f au%42s*\n", dt*time_au_fs, dt, " ");
          printf(  "*  iCurrent:               %23s%42s*\n",iCurrent ? "true" : "false", " ");
          printf(  "*  Model:                  %23s%42s*\n",  iMode.c_str(), " ");
          printf(  "*  TB type:                %23s%42s*\n", TBtype.c_str(), " ");
          printf(  "*  FermiEnergy:            %8.4f%15s%42s*\n", FermiE," ", " ");
          printf(  "*********************************************************************************************\n");  
          printf(  "*                                    Unit cell                                              *\n");
          printf(  "*********************************************************************************************\n");  
          printf(  "*                      a.u.                                   Angstrom                      *\n");                  
          printf(  "*  a1  ->     (%7.2f, %7.2f, %7.2f)%11s(%7.2f, %7.2f, %7.2f)%13s*\n", a[0][0], a[0][1], a[0][2], " ",a[0][0]*space_au_A, a[0][1]*space_au_A, a[0][2]*space_au_A, " ");
          printf(  "*  a2  ->     (%7.2f, %7.2f, %7.2f)%11s(%7.2f, %7.2f, %7.2f)%13s*\n", a[1][0], a[1][1], a[1][2], " ",a[1][0]*space_au_A, a[1][1]*space_au_A, a[1][2]*space_au_A, " ");
          printf(  "*  a3  ->     (%7.2f, %7.2f, %7.2f)%11s(%7.2f, %7.2f, %7.2f)%13s*\n", a[2][0], a[2][1], a[2][2], " ",a[2][0]*space_au_A, a[2][1]*space_au_A, a[2][2]*space_au_A, " ");
          printf(  "*                                                                                           *\n");
          printf(  "*                               Reciprocal lattice vectors:                                 *\n");
          printf(  "*                         a.u.                               Angstrom^-1                    *\n");                  
          printf(  "*  b1  ->     (%7.2f, %7.2f, %7.2f)%11s(%7.2f, %7.2f, %7.2f)%13s*\n", b[0][0], b[0][1], b[0][2]," ", b[0][0]*space_A_au, b[0][1]*space_A_au, b[0][2]*space_A_au," ");
          printf(  "*  b2  ->     (%7.2f, %7.2f, %7.2f)%11s(%7.2f, %7.2f, %7.2f)%13s*\n", b[1][0], b[1][1], b[1][2]," ", b[1][0]*space_A_au, b[1][1]*space_A_au, b[1][2]*space_A_au," ");
          printf(  "*  b3  ->     (%7.2f, %7.2f, %7.2f)%11s(%7.2f, %7.2f, %7.2f)%13s*\n", b[2][0], b[2][1], b[2][2]," ", b[2][0]*space_A_au, b[2][1]*space_A_au, b[2][2]*space_A_au," ");
          printf(  "*                                                                                           *\n");
          printf(  "*  Jacobian determinant: %2.4f                                                             *\n", Coord_B::getJ());  
          printf(  "*********************************************************************************************\n");      
          printf(  "*  k-space:      (#k1,#k2,#k3) -> (%3d, %3d, %3d)                                           *\n", Nk[0], Nk[1], Nk[2]);
          printf(  "*  k-space resolution:         -> (%2.4f, %2.4f, %2.4f)                                  *\n", dk[0], dk[1], dk[2]);
          printf(  "*********************************************************************************************\n");  
          printf(  "*                                   Band types:                                             *\n");
          printf(  "*********************************************************************************************\n");  
          printf(  "*   # conduction bands: %3d                                                                 *\n", Nb[2]);
          printf(  "*   # valence bands:    %3d                                                                 *\n", Nb[1]);
          printf(  "*   # core-hole bands:  %3d                                                                 *\n", Nb[0]);
          printf(  "*********************************************************************************************\n");  
          printf(  "*                               Dissipation terms                                           *\n");
          printf(  "*********************************************************************************************\n");  
          printf(  "*   Dissipation factor T1:      %3.3fa.u.     -     %3.3feV                                 *\n", T1, T1/energy_eV_au);
          printf(  "*   Dissipation factor T2:      %3.3fa.u.     -     %3.3feV                                 *\n", T2, T2/energy_eV_au);
          printf(  "*   Decay core-hole factor Tch: %3.3fa.u.     -     %3.3feV                                 *\n", Tch,Tch/energy_eV_au);
          printf(  "*********************************************************************************************\n");  
          string title = "*                   Optical/IR laser parameters - Pump:                                     *\n*********************************************************************************************";
pulse1.print_par( title );
           
      if(pulse1.gaussian == true)
      {
          printf(  "*  Gaussian profile: Sigma(fs):%22s%8.2f%31s*\n", " ",pulse1.sigma*time_au_fs, " "                       );
          printf(  "*  FWHM(fs):                   %22s%8.2f%31s*\n", " ",2.0*sqrt(2.0*log(2.0))*pulse1.sigma*time_au_fs, " ");
          printf(  "*  FWHM Intensity(fs):         %22s%8.2f%31s*\n", " ",2.0*sqrt(log(2.0))*pulse1.sigma*time_au_fs, " "    );
      }
      else
      {
          printf(  "*  sin2 profile:     number of cycles:               %8.2f                               *\n", pulse1.ncycle                       );
      }
     
          //printf(  "*  Polarization in cartesian coordinates:%9s(%5.2f, %5.2f, %5.2f)%21s*\n"," ", u1.cart[0],u1.cart[1],u1.cart[2], " ");
          //printf(  "*  Polarization in crystal coordinates:  %9s(%5.2f, %5.2f, %5.2f)%21s*\n"," ", u1.crys[0],u1.crys[1],u1.crys[2], " ");
          printf(  "*********************************************************************************************\n");  
          title = "*                         XUV laser parameters - Probe:                                     *\n*********************************************************************************************";
pulse2.print_par(title);
      if(pulse2.gaussian == true)
      {
          printf(  "*  Gaussian profile: Sigma(fs):%22s%8.2f%31s*\n", " ",pulse2.sigma*time_au_fs, " "                       );
          printf(  "*  FWHM(fs):                   %22s%8.2f%31s*\n", " ",2.0*sqrt(2.0*log(2.0))*pulse2.sigma*time_au_fs, " ");
          printf(  "*  FWHM Intensity(fs):         %22s%8.2f%31s*\n", " ",2.0*sqrt(log(2.0))*pulse2.sigma*time_au_fs, " "    );
      }
      else
      {
          printf(  "*  sin2 profile:     number of cycles:               %8.2f                               *\n", pulse2.ncycle                       );
      }
          //printf(  "*  Polarization in cartesian coordinates:         (%5.2f, %5.2f, %5.2f)%21s*\n",u2.cart[0],u2.cart[1],u2.cart[2], " ");
          //printf(  "*  Polarization in crystal coordinates:           (%5.2f, %5.2f, %5.2f)%21s*\n",u2.crys[0],u2.crys[1],u2.crys[2], " ");
          printf(  "*  Delay time:                                %8.2f au - %8.2f fs%21s*\n",DELAY, DELAY*time_au_fs, " ");


      if(iTAbsK)
      {
          printf(  "*********************************************************************************************\n");  
          printf(  "*                  Information about printing on specific k points                          *\n");  
          printf(  "*********************************************************************************************\n");  

          for(int J=0; J<nTAk; J++)
          {
          printf(  "*   %6i                    ---          %10.6f    %10.6f  %10.6f %12s*\n", tagTAk[J], TAkpt[J][0], TAkpt[J][1], TAkpt[J][2], " ");
          }
      }
      printf(  "*********************************************************************************************\n\n");  


}




void Initialize_T(vec2d& T, double& T1, double& T2, double& Tch, int& Nch)
{
    for (int ic=Nch; ic<T.n1(); ic++) {
        for (int jc=Nch; jc<T.n2(); jc++) {
            if(ic==jc) T[ic][jc]=T1;
            else T[ic][jc]=T2;
        }   
    }
    for (int ic=0; ic<Nch; ic++) {
        for (int jc=0; jc<Nch; jc++) {
            T[ic][jc]=Tch;
        }
        for (int jc=Nch; jc<T.n2(); jc++) {
            T[ic][jc]=Tch/2.;
        }
    }
    for (int ic=0; ic<T.n1(); ic++){
        for (int jc=ic+1; jc<T.n2(); jc++) {
            T[jc][ic]=T[ic][jc];
        }
    }
}

template <typename T>
void Initial_Population(vec3x& P0,vec2d& k, T& Model,double FermiE, int Nbands)
{
    //Initial population: Fermi level
    vec3x P_bloch(P0.n1(), P0.n2(), P0.n3());
    //density matrix at the beginning (Bloch)
    for (int ik=0; ik<P0.n1(); ik++) 
    {
         for (int ic=0; ic<Nbands; ic++) 
                 P_bloch[ik][ic][ic]=1.;
    }// End ik
    

    ofstream fp_output;
    
    fp_output.open("P0.txt");
    fp_output << "k1      k2      k3      n1      n2        P0_WANNIER                     P0_Bloch" << endl;
    
    for (int ik=0; ik<P0.n1(); ik++) 
    {   Coord_B kpt;
        kpt.setcrys(k[ik][0],k[ik][1],k[ik][2]);
        vec2x H; 
        vec2x U;
        H.resize(P0.n2(), P0.n2());
        U.resize(P0.n2(), P0.n2());
        Model.energy_U(H, U, kpt);    
        for (int ia=0; ia<P0.n2(); ia++) 
        {
            for(int ib=0; ib<P0.n2(); ib++)
            {
                for(int ii=0; ii<P0.n2(); ii++)
                {
                    for(int ij=0; ij<P0.n2(); ij++)
                    {
                        P0[ik][ia][ib]+=U[ii][ia]*P_bloch[ik][ii][ij]*conj(U[ij][ib]);   
                    }//end ij
                }//end ii
                fp_output << ik << "   "  <<ia << "     " << ib << "     " << P0[ik][ia][ib]<<"       "<< P_bloch[ik][ia][ib]<<endl;
            }//end ib
        }//end ia
        fp_output<<endl;
       
    }//end k
    fp_output.close();
}


