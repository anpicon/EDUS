// Functions to determine Fourier image of Coloumb repulsion
struct Coulomb_parameters {
  bool Coulomb_calc; // by default it false, if we don't initialize in in input file
  int Ncut; // number of terms in Fourier series
  vec1d A; // coefficients of Fourier transform
  vec1d D;
  int N_BZ_points_total;

  vec1x X_cc, X_sc, X_cs, X_ss; // coefficients for Fourier transformed potential
  vec1x X_cc0, X_sc0, X_cs0, X_ss0; // coefficients for Fourier transformed potential
//  vec2x cos_nx_crys, sin_nx_crys, cos_ny_crys, sin_ny_crys; // matrices of cos and sin in crystal basis
  vec3x P_static;
  vec3x Xk_storage;
  vec3x P0_dia, P_dkx, P_dky, P_d2kx, P_d2ky,  P_dky_dkx;
  vec3x Hk_renorm_shared;

  double qTF; // Thomas-Fermi screening
  double epsilon_static;
  double max;
  double min;
  double n_cond;
  double det_J_for_k_sum;

  vec1d Screen_const; // [0] - without derivatives
  // [1] - d^2 /dx^2 // [2] - d^2/dy^2 // [3] - d^2/(dxdy) 

  vec1x f_cc, f_sc, f_cs, f_ss;
  double test;
  bool Wannie_basis;
  bool Diagonal_basis;
  bool Rytova_Keldysh;
  bool Print_new_band_dispersion;

  vec2d Ekdia;
  vec1x Pk_shared; // shared array to print in file 
  vec1d N_shared; // population for shared memory
  vec1d N_exciton;
  ofstream Ek_file;
  ofstream P_stream;
  ofstream exciton_file;

};


// function to make borders smooth
double falpha(double alpha, double q_cut, double q){
  double f = 1.0 / (exp(- alpha * (q_cut - abs(q)) ) + 1 );
  return f;
}






// Thomas Fermi potential
double Coloumb_THh_F(double qx, double qy, double qTF, double eps_area, 
  double q_cut, double Vbord_x, double Vbord_y, double V_corner,
  bool Rytova_Keldysh){
  Coord_B k1; // !!!! THIS LINE SLOW
  k1.setcrys(0.0, qx, qy);
  // k1.setcrys(0.0, qy, qx); // rotate 90
  double Num, Den;

  double q_mod = k1.dot(k1); 
  double q_border = 0.2; // borderline where we use Thomas-Fermi screening
  double a = 2.5/space_au_A; // unit cell

  // function to make borders smooth
  double alpha = 250;
  double fx = falpha(alpha, q_cut, qx);
  double fy = falpha(alpha, q_cut, qy);
  double f_m_x = falpha(-alpha, q_cut, qx);
  double f_m_y = falpha(-alpha, q_cut, qy);

  if(Rytova_Keldysh){
    double r0 = 18.89726; //a.u. BN on quarz Henriques et al https://doi.org/10.1088/1361-648X/ab47b3
    // double r0 = 20; //a.u. BN on quarz Henriques et al https://doi.org/10.1088/1361-648X/ab47b3
    //double eps1 = 2.4; // quarz Henrique 2019 
    double eps1 = 1; // vacuum Galvani 2016, Ridolfi 2019
    Num = 2* M_PI;
    //if (q_mod < q_border)
    // {
    //   q_mod += qTF*qTF;

    // }
    Den = eps_area * sqrt(q_mod + qTF*qTF)*( sqrt(q_mod) * r0 +  eps1);

  } else{ // simpe 1/q potential


    Num= 1; // just model random numbers to debug
    //  double Den= pow((qx/2),2) + pow((sqrt(3.0)*qx/2 + qy),2) + qTF*qTF;
    Den= eps_area*sqrt(q_mod + qTF*qTF);
  }


  double Vq = Num/Den;
  // potential should be 0 in the borders as a periodic function
  Vq *= fx * fy;
  Vq += Vbord_x * f_m_x * fy + Vbord_y * f_m_y * fx +  V_corner * f_m_x * f_m_y;

  return Vq;
}








// integrations to get Fourier series coefficients. Input: m,n - order of term, ABCD result coefficient, N_sum Number of points, q - array of 1d grid, sp_arr - array of simpson coefficients
double Four_integr_Coloumb(int m, int n, vector<double>  &AD, int N_sum, 
  vec1d  &q, vec1d  &sp_arr, double qTF, double eps0,
  double q_cut, double Vbord_x, double Vbord_y, double V_corner,
   bool Rytova_Keldysh){
  AD[0]=0.0;
  AD[1]=0.0;
  double Coef=0.0;
  // integral in simpson method
  for(int x=0; x<N_sum; x++){
    for(int y=0; y<N_sum; y++){
      double trigm_cos,trigm_sin, trign_cos,trign_sin, spxy, mq, nq;
      spxy=sp_arr[x*N_sum + y]*Coloumb_THh_F(q[x], q[y],  qTF,  eps0, 
         q_cut, Vbord_x, Vbord_y, V_corner, Rytova_Keldysh);// 2D Simpson coefficient and Coloumb
      mq=2* M_PI * m * q[x];
      nq=2* M_PI * n * q[y];
      trigm_cos = cos(mq); // cos(2pi m qx)
      trigm_sin = sin(mq); // sin(2pi m qx)
      trign_cos = cos(nq); // cos(2pi m qy)
      trign_sin = sin(nq); // sin(2pi m qy)

      double dA=4*spxy*trigm_cos*trign_cos ;
      double dD=4*spxy*trigm_sin*trign_sin ;

      AD[0]=AD[0]+dA;
      AD[1]=AD[1]+dD;
    }
  }

  for(int i=0; i<2; i++){ //round off nearly zero coefficients
 //   if ((AD[i]*AD[i])<0.000001) {  AD[i]=0.0;    }
  }

  if (m==0 && n==0){ AD[0]=AD[0]/4;} // first coefficient is A/4, B,C,D=0
  else if (m==0 || n==0){ // second coefficient is A/2, B/2, C/2, D=0
    AD[0]=AD[0]/2;
    AD[1]=AD[1]/2;
    }


  return 0;
}








// calulates Fourie-Series approximation
double Fourier_series_Coloumb(double kx_temp,double ky_temp, int Ncut, vec1d &A, vec1d &D){
  double VF=0.0;
  double trigm_cos,trigm_sin, trign_cos,trign_sin, dA, dD;
  for (int m =0; m<Ncut; m++ ){ //column
    for (int n =0; n<Ncut; n++ ){ // row
      trigm_cos = cos(2* M_PI * m * kx_temp); // cos(2pi m qx)
      trigm_sin = sin(2* M_PI * m * kx_temp); // sin(2pi m qx)
      trign_cos = cos(2* M_PI * n * ky_temp); // cos(2pi m qy)
      trign_sin = sin(2* M_PI * n * ky_temp); // sin(2pi m qy)
      dA = A[n+Ncut*m]*trigm_cos*trign_cos;
      dD = D[n+Ncut*m]*trigm_sin*trign_sin;
    VF= VF + dA + dD;
    }
   }


  return VF;
}



// function to calculate coefficients for Coulomb band reconstruction X
// For MPI implementation
void Calculate_X_coefficients_MPI(
  vec3x & P,
  Private_omp_parameters& OMP_private,
  int Ncv,
  Coulomb_parameters& Coulomb_set, 
  int root_rank, int my_rank, 
  trig_coefficients & trig_k)
{

  double N_BZ = Coulomb_set.N_BZ_points_total;
  int Ncut = Coulomb_set.Ncut;
  int Ncut2Ncv2 = Ncut*Ncut*Ncv*Ncv;
  double Nk= (sqrt( N_BZ ));
  double dk = 1.0 / Nk;


  vec1x f_cc(Ncv*Ncv), f_sc(Ncv*Ncv), f_cs(Ncv*Ncv), f_ss(Ncv*Ncv);


  std::complex<double> P_current, P_d2ky_d2kx, P_dkx, P_dky;
  complex<double> corr_cc =0.0;
  complex<double> corr_cs =0.0;
  complex<double> corr_sc =0.0;
  complex<double> corr_ss =0.0; // corerections due to the derivatives
  double cos_cos_mn, cos_sin_mn, sin_cos_mn, sin_sin_mn;
  double integrWeight_ik;

  #pragma omp barrier

  for (int mm =0; mm< Ncut*Ncut; mm++ ){ //row summation over Fourier series
    int m = mm / Ncut; // integer division
    int n = mm % Ncut; // leftover   
      f_cc.fill(0.); // fill f with zeros for summation
      f_cs.fill(0.);
      f_sc.fill(0.);
      f_ss.fill(0.);
      if (OMP_private.id_masters[OMP_private.thr_total - 1]) Coulomb_set.f_ss.fill(0.); // shared memory
      if (OMP_private.id_masters[OMP_private.thr_total - 2]) Coulomb_set.f_cc.fill(0.);
      if (OMP_private.id_masters[OMP_private.thr_total - 3]) Coulomb_set.f_sc.fill(0.);
      if (OMP_private.id_masters[OMP_private.thr_total - 4]) Coulomb_set.f_cs.fill(0.);


      // if (id_masters[thr_total - 4] == true) Coulomb_set.test = 0.0;
      // #pragma omp barrier
      for (int ik=0; ik<OMP_private.lenght_k; ik++){ // summation over k-vectors
        
        integrWeight_ik = OMP_private.integrWeight[ik];
        cos_cos_mn = trig_k.cos_mkx[m*OMP_private.lenght_k + ik] * trig_k.cos_nky[n*OMP_private.lenght_k + ik];
        cos_sin_mn = trig_k.cos_mkx[m*OMP_private.lenght_k + ik] * trig_k.sin_nky[n*OMP_private.lenght_k + ik];
        sin_cos_mn = trig_k.sin_mkx[m*OMP_private.lenght_k + ik] * trig_k.cos_nky[n*OMP_private.lenght_k + ik];
        sin_sin_mn = trig_k.sin_mkx[m*OMP_private.lenght_k + ik] * trig_k.sin_nky[n*OMP_private.lenght_k + ik];
        



        //  #pragma omp simd !!!

        for (int ic=0; ic<Ncv; ic++){// summation over bands
          for (int jc=ic; jc<Ncv; jc++){ // summation over bands
            //Pi[ic][jc] = P[ik][ic][jc]; for simd
            P_d2ky_d2kx = Coulomb_set.P_d2ky[ik + OMP_private.begin_count][ic][jc]+ \
              Coulomb_set.P_d2kx[ik + OMP_private.begin_count][ic][jc];

              // see section in notes
            P_current =  P[ik + OMP_private.begin_count][ic][jc]; // not to call element of large array many times            
            P_dkx = Coulomb_set.P_dkx[ik + OMP_private.begin_count][ic][jc];
            P_dky = Coulomb_set.P_dky[ik + OMP_private.begin_count][ic][jc];

            // corr_cc =  cos_cos_mn*(- 4* M_PI * M_PI *(m*m + n*n )* P_current + P_d2ky_d2kx) \
            //   - 4*M_PI*m*sin_cos_mn * P_dkx  \
            //   - 4*M_PI*n*cos_sin_mn * P_dky;

            // corr_cs =  cos_sin_mn*(- 4* M_PI * M_PI *(m*m + n*n )* P_current + P_d2ky_d2kx) \
            //   - 4*M_PI*m*sin_sin_mn * P_dkx \
            //   + 4*M_PI*n*cos_cos_mn * P_dky;

            // corr_sc =  sin_cos_mn*(- 4* M_PI * M_PI *(m*m + n*n )* P_current + P_d2ky_d2kx) \
            //   + 4*M_PI*m*cos_cos_mn * P_dkx \
            //   - 4*M_PI*n*sin_sin_mn * P_dky;
  
            // corr_ss =  sin_sin_mn*(- 4* M_PI * M_PI *(m*m + n*n )* P_current + P_d2ky_d2kx) \
            //   + 4*M_PI*m*cos_sin_mn * P_dkx \
            //   + 4*M_PI*n*sin_cos_mn * P_dky;

            // corr_cc *= (dk*dk / 27);  
            // corr_cs *= (dk*dk / 27);  
            // corr_sc *= (dk*dk / 27);  
            // corr_ss *= (dk*dk / 27);  

            // f_cc[ic*Ncv + jc] += integrWeight_ik * (P_current * cos_cos_mn + corr_cc);
            // f_cs[ic*Ncv + jc] += integrWeight_ik * (P_current * cos_sin_mn + corr_cs);
            // f_sc[ic*Ncv + jc] += integrWeight_ik * (P_current * sin_cos_mn + corr_sc);
            // f_ss[ic*Ncv + jc] += integrWeight_ik * (P_current * sin_sin_mn + corr_ss);


            f_cc[ic*Ncv + jc] += (P_current * cos_cos_mn) / N_BZ;
            f_cs[ic*Ncv + jc] += (P_current * cos_sin_mn) / N_BZ;
            f_sc[ic*Ncv + jc] += (P_current * sin_cos_mn) / N_BZ;
            f_ss[ic*Ncv + jc] += (P_current * sin_sin_mn) / N_BZ;



          }// summation over bands
        }// summation over bands
      } // // summation over k-vectors
      //  created f coefficients for this particular (m,n) and band [ic][jc]
   
    #pragma omp barrier
    // summing local f_ into shared f_
    #pragma omp critical
    { // summation over threads

      for (int ic=0; ic<Ncv; ic++){// summation over bands
        for (int jc=ic; jc<Ncv; jc++){ // summation over bands OMP
          Coulomb_set.f_ss[ic*Ncv + jc] += Coulomb_set.det_J_for_k_sum * f_ss[ic*Ncv + jc];
          Coulomb_set.f_cc[ic*Ncv + jc] += Coulomb_set.det_J_for_k_sum * f_cc[ic*Ncv + jc];
          Coulomb_set.f_sc[ic*Ncv + jc] += Coulomb_set.det_J_for_k_sum * f_sc[ic*Ncv + jc];
          Coulomb_set.f_cs[ic*Ncv + jc] += Coulomb_set.det_J_for_k_sum * f_cs[ic*Ncv + jc];

          
        }
      }
    }
   #pragma omp barrier
   //if (id_masters[thr_total - 1] == true) std::cout << " Coulomb_set.test after " << Coulomb_set.test  << endl << endl;
    

  #pragma omp master
  {
  // cout << f_cc[0][1*Ncv + 1] << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    // this function sums all elements of vec2x objects
    // result is stored in root_rank memory
    MPI_reduce_vec1x(Coulomb_set.f_cc, root_rank, my_rank, Ncv*Ncv);
    MPI_reduce_vec1x(Coulomb_set.f_cs, root_rank, my_rank, Ncv*Ncv);
    MPI_reduce_vec1x(Coulomb_set.f_sc, root_rank, my_rank, Ncv*Ncv); 
    MPI_reduce_vec1x(Coulomb_set.f_ss, root_rank, my_rank, Ncv*Ncv);  
    // cout << f_cc[0][1*Ncv + 1] << endl;
    MPI_Barrier(MPI_COMM_WORLD);
  } // omp master
  
  

  #pragma omp barrier
  if (my_rank == root_rank) { 
    
    #pragma omp master
    {
     

      int index_dia, index_offdia;
      double Amn = Coulomb_set.A[mm];
      double Dmn = Coulomb_set.D[mm];
      for (int ic=0; ic<Ncv; ic++){// summation over bands
        for (int jc=ic; jc<Ncv; jc++){ // summation over bands
          index_dia = ic*Ncv + jc;
          index_offdia = jc*Ncv + ic;

          Coulomb_set.X_cc[mm*Ncv*Ncv + index_dia] =  (Amn* Coulomb_set.f_cc[index_dia] + Dmn* Coulomb_set.f_ss[index_dia]) - Coulomb_set.X_cc0[mm*Ncv*Ncv +index_dia];
          if(jc != ic){
            Coulomb_set.X_cc[mm*Ncv*Ncv +index_offdia] = conj(Coulomb_set.X_cc[mm*Ncv*Ncv + index_dia]);
          }

          Coulomb_set.X_cs[mm*Ncv*Ncv +index_dia] =  (Amn* Coulomb_set.f_cs[index_dia] - Dmn* Coulomb_set.f_sc[index_dia]) - Coulomb_set.X_cs0[mm*Ncv*Ncv +index_dia];
          if(jc != ic){
            Coulomb_set.X_cs[mm*Ncv*Ncv +index_offdia] = conj(Coulomb_set.X_cs[mm*Ncv*Ncv +index_dia]);
          }

          Coulomb_set.X_sc[mm*Ncv*Ncv +index_dia] =  (Amn* Coulomb_set.f_sc[index_dia] - Dmn* Coulomb_set.f_cs[index_dia]) - Coulomb_set.X_sc0[mm*Ncv*Ncv +index_dia];
          if(jc != ic){
            Coulomb_set.X_sc[mm*Ncv*Ncv +index_offdia] = conj(Coulomb_set.X_sc[mm*Ncv*Ncv +index_dia]);
          }

          Coulomb_set.X_ss[mm*Ncv*Ncv +index_dia] =  (Amn* Coulomb_set.f_ss[index_dia] + Dmn* Coulomb_set.f_cc[index_dia]) - Coulomb_set.X_ss0[mm*Ncv*Ncv +index_dia];
          if(jc != ic){
            Coulomb_set.X_ss[mm*Ncv*Ncv +index_offdia] = conj(Coulomb_set.X_ss[mm*Ncv*Ncv +index_dia]);
          }
        }// summation over bands
      }// summation over bands
      
     }// omp master

    } //if (my_rank == root_rank)

    #pragma omp barrier

  } //Fourier series

  


  #pragma omp barrier
  
  
  #pragma omp master
  { 
    MPI_Barrier(MPI_COMM_WORLD);
    

    MPI_Bcast(& Coulomb_set.X_cc[0], Ncut2Ncv2, MPI_C_DOUBLE_COMPLEX, root_rank, MPI_COMM_WORLD);
    
    MPI_Bcast(& Coulomb_set.X_cs[0], Ncut2Ncv2, MPI_C_DOUBLE_COMPLEX, root_rank, MPI_COMM_WORLD);

    MPI_Bcast(& Coulomb_set.X_sc[0], Ncut2Ncv2, MPI_C_DOUBLE_COMPLEX, root_rank, MPI_COMM_WORLD);

    MPI_Bcast(& Coulomb_set.X_ss[0], Ncut2Ncv2, MPI_C_DOUBLE_COMPLEX, root_rank, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
  } // omp master

  
  
  #pragma omp barrier
  
  
} // end of func





// calculate
void Calculate_X_MPI(int ik, vec2x& Xk, 
  int Ncv, Coulomb_parameters& Coulomb_set,
  trig_coefficients & trig_k, int nk_pr)
{
  int Ncut = Coulomb_set.Ncut;
  int index_dia;
  double cos_cos_mn, cos_sin_mn, sin_cos_mn, sin_sin_mn;
  // Matrix of Coulomb contribution to hamiltonian in Wannier
  Xk.fill(0.);
  vec1x alpha(Ncv); 
  alpha.fill(0.);
  int m, n;
  for (int mm =0; mm< Ncut*Ncut; mm++ ){ //row summation over Fourier series
    m = mm / Ncut; // integer division
    n = mm % Ncut; // leftover
    cos_cos_mn = trig_k.cos_mkx[m*nk_pr + ik] * trig_k.cos_nky[n*nk_pr + ik];
    cos_sin_mn = trig_k.cos_mkx[m*nk_pr + ik] * trig_k.sin_nky[n*nk_pr + ik];
    sin_cos_mn = trig_k.sin_mkx[m*nk_pr + ik] * trig_k.cos_nky[n*nk_pr + ik];
    sin_sin_mn = trig_k.sin_mkx[m*nk_pr + ik] * trig_k.sin_nky[n*nk_pr + ik];


    for (int ic=0; ic<Ncv; ic++){// summation over bands
      #pragma omp simd
      for (int jc=0; jc<Ncv; jc++){ // summation over bands
        index_dia = mm*Ncv*Ncv + ic*Ncv + jc;
        alpha[jc]  =  Coulomb_set.X_cc[index_dia] *cos_cos_mn +  Coulomb_set.X_cs[index_dia] *cos_sin_mn +
	  Coulomb_set.X_sc[index_dia] *sin_cos_mn + Coulomb_set.X_ss[index_dia] *sin_sin_mn;
        Xk[jc][ic] -= alpha[jc]; // INDEX ORDER AND SIGN see notes: X_{jc,ic}(k) = - sum_{k'} V(k-k') \rho_{ic,jc}(k')
      }
    }
  }
}


//








void get_dkx_dky_general( vec3x& P_dky, vec3x& P_dkx, vec3x& P0_dia,
 Coulomb_parameters& Coulomb_set, Private_omp_parameters& OMP_private, bool ifdky){
  int ik_max = P0_dia.n1();


  int Ncv = P0_dia.n3();
  

  int Nk= int(sqrt( Coulomb_set.N_BZ_points_total ));
  double Nk_doub = Nk;
  double dk = 1.0 / Nk_doub; 

  int Nk_y_top = (P0_dia.n1() / Nk) - 1;

  int ikx0, iky0; // coordinates of current position
  int ikx_s, iky_s; // coordinates of shifted position

  vec2i index_dik(3,3); // contain indices of neighbours 
  vec2x P_neighb(3, 3);

  int ik;
  for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++){
    ik = ik_pr + OMP_private.begin_count;
    ikx0 = ik % Nk;
    iky0 = (ik / Nk);
    for (int diky = -1; diky < 2; diky++){ // calculate index of all neigbours
      for (int dikx = -1; dikx < 2; dikx++){
        ikx_s = ikx0 + dikx;
        if (ikx_s < 0) {
          ikx_s += Nk;
        } else if (ikx_s > (Nk-1)) {
          ikx_s -= Nk;
        }

        iky_s = iky0 + diky;
        if (iky_s < 0) {
          iky_s = Nk_y_top;
        } else if (iky_s > Nk_y_top ) {
          iky_s = 0;
        }

        index_dik[diky+1][dikx+1] = iky_s * Nk + ikx_s;

        // if (ik == 0){
        //   cout << "ik = " << ik << "iky_s = " << iky_s << " ikx_s = " << ikx_s << endl;
        // }

      }
    }
    // the derivatives: 
    for (int ic=0; ic < Ncv; ic++){
      for (int jc = ic; jc<Ncv; jc++){
        for (int diky = -1; diky < 2; diky++){ // calculate index of all neigbours
          for (int dikx = -1; dikx < 2; dikx++){
            P_neighb[diky+1][dikx+1] = P0_dia[ index_dik[diky+1][dikx+1] ][ic][jc];
          }
        }

        if(ifdky){ // only if we need d/dy
          P_dky[ik][ic][jc] = 0.0;
          for (int dikx = -1; dikx < 2; dikx++){
            P_dky[ik][ic][jc] += (P_neighb[2][dikx+1] - P_neighb[0][dikx+1]);
          }
          P_dky[ik][ic][jc] /= (6 * dk);
          P_dky[ik][jc][ic] = conj(P_dky[ik][ic][jc]);
        }



        P_dkx[ik][ic][jc] = 0;
        for (int diky = -1; diky < 2; diky++){
          P_dkx[ik][ic][jc] += (P_neighb[diky+1][2] - P_neighb[diky+1][0]);
        }
         P_dkx[ik][ic][jc] /= (6 * dk);
         P_dkx[ik][jc][ic] = conj(P_dkx[ik][ic][jc]);

      }
    }
  } // ik

}





// this function and get_dkx are just a wrappers for
//void get_dkx_dky_general() but depending on name we can take both derivatives or only one
void get_dkx_dky( vec3x & P_dky, vec3x & P_dkx, vec3x & P0_dia,
 Coulomb_parameters & Coulomb_set, Private_omp_parameters & OMP_private ){
  
  bool ifdky = true; //need d/dky
  get_dkx_dky_general( P_dky, P_dkx, 
                P0_dia, Coulomb_set, OMP_private, ifdky);
}




// this function and get_dkx_dky are just a wrappers for
//void get_dkx_dky_general() but depending on name we can take both derivatives or only one
void get_dkx( vec3x& P_dkx, vec3x& P0,
 Coulomb_parameters& Coulomb_set, Private_omp_parameters& OMP_private){
  vec3x P_dky(2,2,2); // just emty vector, don't need it further
  bool ifdky = false; // dont need d/dky
  get_dkx_dky_general(P_dky, P_dkx, 
                P0, Coulomb_set, OMP_private, ifdky);
}


