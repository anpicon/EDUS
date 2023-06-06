// Functions to determine Fourier image of Coloumb repulsion
struct Coulomb_parameters {
  bool Coulomb_calc; // by default it false, if we don't initialize in in input file
  int Ncut; // number of terms in Fourier series
  vec2x A, B, C, D; // coefficients of Fourier transform

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
  double P_min;

  vec2x Screen_const; // [0] - without derivatives
  // [1] - d^2 /dx^2 // [2] - d^2/dy^2 // [3] - d^2/(dxdy) 

  vec1x f_cc, f_sc, f_cs, f_ss;
  double test;
  bool Wannie_basis;
  bool Diagonal_basis;
  bool Rytova_Keldysh;
  bool Print_new_band_dispersion;
  string Sample_orientation; // by default in yz plane and is 011
  double r0; // R-K polarizability parameter


  vec2d Ekdia;
  vec1x Pk_shared; // shared array to print in file 
  vec1d N_shared; // population for shared memory
  vec1d N_exciton;
  ofstream Ek_file;
  ofstream P_stream;
  ofstream exciton_file;
  vec1d Taylor_Delta_Pk_max, Taylor_alpha_max;
  double Taylor_Pv_max;

};


// function to make borders smooth
double falpha(double alpha, double q_cut, double q){
  double f = 1.0 / (exp(- alpha * (q_cut - abs(q)) ) + 1 );
  return f;
}






// Thomas Fermi potential
double Coloumb_THh_F(Coord_B& k1, double qx, double qy, double qTF, double eps_area, 
            double r0, double Vbord_x, double Vbord_y, double V_corner,
            bool Rytova_Keldysh){

    double Num, Den;

    double q_mod = k1.dot(k1); 
    // double q_border = 0.2; // borderline where we use Thomas-Fermi screening


    // function to make borders smooth
    // double alpha = 250;
    // double fx = falpha(alpha, q_cut, qx);
    // double fy = falpha(alpha, q_cut, qy);
    // double f_m_x = falpha(-alpha, q_cut, qx);
    // double f_m_y = falpha(-alpha, q_cut, qy);

    if(Rytova_Keldysh){
        // double r0 = 18.89726; //a.u. BN on quarz Henriques et al https://doi.org/10.1088/1361-648X/ab47b3
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
    // Vq *= fx * fy;
    // Vq += Vbord_x * f_m_x * fy + Vbord_y * f_m_y * fx +  V_corner * f_m_x * f_m_y;

    return Vq;
}


// in parallel region
void Create_Wq_with_G_sum(vec3x &Wq, int Ncv, vec1d  &q,
    int & iq_loc, int & qx_global,int &  qy_global,
    vec2d & Displacement_orb, int & G_num, 
    bool Rytova_Keldysh, double q_cut, double Vbord_x, double Vbord_y, double V_corner,
    double qTF, double eps0, int N_sum, 
    Coulomb_parameters & Coulomb_set, vec2d & G_vec){


    Coord_B k1; 


    double Vq, phi;

    complex<double> Vph_q, phase_q;  // Vq bare Coulomb, Wq ater exponent and G summation


    
        // sum over neighbouring unit cells to make potential symmetric
    for (int i_G = 0; i_G < G_num; i_G++){
        if (Coulomb_set.Sample_orientation == "011"){ // sample yz
            k1.setcrys(0.0, q[qx_global] + G_vec[i_G][0], q[qy_global] + G_vec[i_G][1]); //to create k-points in crystal coordinates for phase and Coulomb
        } else if (Coulomb_set.Sample_orientation == "110"){ // sample in xy
            k1.setcrys( q[qx_global] + G_vec[i_G][0], q[qy_global] + G_vec[i_G][1] , 0.0); //to create k-points in crystal coordinates for phase and Coulomb
        }

        
        
        Vq = Coloumb_THh_F(k1, q[qx_global], q[qy_global],  qTF,  eps0, 
                                        Coulomb_set.r0, Vbord_x, Vbord_y, V_corner, Rytova_Keldysh);// 2D Simpson coefficient and Coloumb
     

        int index_linear;
        for (int ic=0; ic<Ncv; ic++){// summation over bands
            for (int jc=ic; jc<Ncv; jc++){ // summation over bands
                index_linear = ic*Ncv + jc;
                if (ic == 0 and jc == 0){
                    Wq[iq_loc][ic][jc] += Vq;

                }
                if(jc != ic){
                    //phi = -k1.cart[1] * 2.72758; // Simple TB
                    phi = 0.0;
                    for(int i_coord=0; i_coord<3; i_coord++){ // phase depends on displacements
                        phi += k1.cart[i_coord] * (Displacement_orb[ic][i_coord] - Displacement_orb[jc][i_coord]);
                    }

                    phase_q =  exp(c1 * phi);
                    //real coefficients contain A and D Coefficients when imaginary contain B and C coefficients  
                    Vph_q = Vq * phase_q;
                    Wq[iq_loc][ic][jc] += Vph_q;

                }
            }
        }
    }


    // for (int ic=0; ic<Ncv; ic++){// summation over bands
    //     for (int jc=ic+1; jc<Ncv; jc++){ // summation over bands
    //         Wq[iq_loc][jc][ic] = conj(Wq[iq_loc][ic][jc]);
    //     }
    // }


}





// integrations to get Fourier series coefficients. Input: m,n - order of term, ABCD result coefficient, N_sum Number of points, q - array of 1d grid, sp_arr - array of simpson coefficients
double Four_integr_Coloumb(int m, int n, vec2x  &AD, int N_sum, 
    vec1d  &q, vec1d  &sp_arr, 
    int & begin_count, int & end_count, int & lenght_q,
     vec3x& Wq, int Ncv){

    AD.fill(0.0);
    double Coef=0.0;
    // integral in simpson method
    int iq_global, qx_global, qy_global;
    double trigm_cos,trigm_sin, trign_cos,trign_sin, spxy, mq, nq, dA, dD, dB, dC;
    complex<double> Vq;
    int index_linear, index_linear_tr;


    for (int iq_loc = 0; iq_loc < lenght_q; iq_loc++){// create 2d simpson array
        iq_global = begin_count + iq_loc;
        qx_global = iq_global / (N_sum); // integer division
        qy_global = iq_global % (N_sum); // leftover; doesn't matter which is x and y


        spxy=sp_arr[iq_loc];
        mq=2* M_PI * m * q[qx_global];
        nq=2* M_PI * n * q[qy_global];
        trigm_cos = cos(mq); // cos(2pi m qx)
        trigm_sin = sin(mq); // sin(2pi m qx)
        trign_cos = cos(nq); // cos(2pi m qy)
        trign_sin = sin(nq); // sin(2pi m qy)

        dA = 4 * spxy * trigm_cos * trign_cos;
        dD = 4 * spxy * trigm_sin * trign_sin;
        dB = 4 * spxy * trigm_cos * trign_sin;
        dC = 4 * spxy * trigm_sin * trign_cos;


        
        for (int ic=0; ic<Ncv; ic++){// summation over bands
            for (int jc=ic; jc<Ncv; jc++){ // summation over bands

                index_linear = ic*Ncv + jc;
                Vq = Wq[iq_loc][ic][jc];
                if (ic == 0 and jc == 0){
                    //diagonal only for 0,0
                    AD[0][0] += Vq * dA; // diagonal elements, all the same because no phase here
                    AD[3][0] += Vq * dD;
                }
                if(jc != ic){
                         
                    AD[0][index_linear] += dA * Vq;
                    AD[1][index_linear] += dB * Vq;
                    AD[2][index_linear] += dC * Vq;
                    AD[3][index_linear] += dD * Vq;

                }
            }
        }
    }
          





        
    
    for (int i1=0; i1<AD.n1(); i1++){
        for(int i2=0; i2<AD.n2(); i2++){ //
            if (m==0 && n==0){ 
                AD[i1][i2] *= 0.25;
            } // first coefficient is A/4, B,C,D=0
            else if (m==0 || n==0){ // second coefficient is A/2, B/2, C/2, D=0
                AD[i1][i2] *= 0.5;
            }
        }
    }

    // for (int i1=0; i1<AD.n1(); i1++){
    //     for (int ic=0; ic<Ncv; ic++){// summation over bands
    //         for (int jc=ic+1; jc<Ncv; jc++){ // summation over bands
    //             index_linear = ic*Ncv + jc;
    //             index_linear_tr = jc*Ncv + ic;
    //             AD[i1][index_linear_tr] = conj(AD[i1][index_linear]);
    //         }
    //     }
    // }

    return 0;
}








// calulates Fourie-Series approximation
complex<double> Fourier_series_Coloumb(double kx_temp,double ky_temp, int Ncut, 
                                vec2x &A, vec2x &B,vec2x &C, vec2x &D, int i_band){
    complex<double> VF=0.0;
    double trigm_cos,trigm_sin, trign_cos,trign_sin;
    complex<double> dA, dB, dC, dD;
    for (int m =0; m<Ncut; m++ ){ //column
        for (int n =0; n<Ncut; n++ ){ // row
            trigm_cos = cos(2* M_PI * m * kx_temp); // cos(2pi m qx)
            trigm_sin = sin(2* M_PI * m * kx_temp); // sin(2pi m qx)
            trign_cos = cos(2* M_PI * n * ky_temp); // cos(2pi m qy)
            trign_sin = sin(2* M_PI * n * ky_temp); // sin(2pi m qy)
            dA = A[n+Ncut*m][i_band] * trigm_cos * trign_cos;
            dB = B[n+Ncut*m][i_band] * trigm_cos * trign_sin;
            dC = C[n+Ncut*m][i_band] * trigm_sin * trign_cos;
            dD = D[n+Ncut*m][i_band] * trigm_sin * trign_sin;

            VF += dA + dB + dC + dD;
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


  std::complex<double> P_current;

  double cos_cos_mn, cos_sin_mn, sin_cos_mn, sin_sin_mn, corr_trig;
  double integrWeight_ik;
  integrWeight_ik = 1.0 / N_BZ;

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

    // correction terms come from fast oscillations of sin and cos
    // we have \int P * cos(pi* m* kx) * sin (pi* n* kx)
    corr_trig = - 4* M_PI * M_PI *(m*m + n*n )/(N_BZ* 27) ;// 
    // we also neglect everything linear in m and n or that doesn't contain it
    // since our population function is smooth and slow
    // we can analytically integrate trigonometrics between points assuming that population is constant there (neglecting P' and P'')
    // for example, Nk=200, m = 10, n = 10
    // corr_trig = - 4 *pi*pi*(100+100)/ (200*200*27) = -0.0073
    // not so small!
    // corr_trig= 0.;
      // if (id_masters[thr_total - 4] == true) Coulomb_set.test = 0.0;
    // #pragma omp barrier
    
      for (int ik=0; ik<OMP_private.lenght_k; ik++){ // summation over k-vectors
        
        // integrWeight_ik = OMP_private.integrWeight[ik];
        cos_cos_mn = trig_k.cos_mkx[m*OMP_private.lenght_k + ik] * trig_k.cos_nky[n*OMP_private.lenght_k + ik];
        cos_sin_mn = trig_k.cos_mkx[m*OMP_private.lenght_k + ik] * trig_k.sin_nky[n*OMP_private.lenght_k + ik];
        sin_cos_mn = trig_k.sin_mkx[m*OMP_private.lenght_k + ik] * trig_k.cos_nky[n*OMP_private.lenght_k + ik];
        sin_sin_mn = trig_k.sin_mkx[m*OMP_private.lenght_k + ik] * trig_k.sin_nky[n*OMP_private.lenght_k + ik];
        


        //  #pragma omp simd !!!

        for (int ic=0; ic<Ncv; ic++){// summation over bands
          for (int jc=ic; jc<Ncv; jc++){ // summation over bands

              // see section in notes
            P_current =  integrWeight_ik * P[ik + OMP_private.begin_count][ic][jc]; // not to call element of large array many times            

            f_cc[ic*Ncv + jc] += P_current * cos_cos_mn* (1. + corr_trig);
            f_cs[ic*Ncv + jc] += P_current * cos_sin_mn* (1. + corr_trig);
            f_sc[ic*Ncv + jc] += P_current * sin_cos_mn* (1. + corr_trig);
            f_ss[ic*Ncv + jc] += P_current * sin_sin_mn* (1. + corr_trig);


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
          Coulomb_set.f_ss[ic*Ncv + jc] += f_ss[ic*Ncv + jc];
          Coulomb_set.f_cc[ic*Ncv + jc] += f_cc[ic*Ncv + jc];
          Coulomb_set.f_sc[ic*Ncv + jc] += f_sc[ic*Ncv + jc];
          Coulomb_set.f_cs[ic*Ncv + jc] += f_cs[ic*Ncv + jc];

          
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
            complexd Amn = 0.0;// = Coulomb_set.A[mm];
            complexd Bmn = 0.0;// = Coulomb_set.D[mm];
            complexd Cmn = 0.0;// = Coulomb_set.A[mm];
            complexd Dmn = 0.0;// = Coulomb_set.D[mm];
            for (int ic=0; ic<Ncv; ic++){// summation over bands
                for (int jc=ic; jc<Ncv; jc++){ // summation over bands
                    
                    index_dia = ic*Ncv + jc;
                    index_offdia = jc*Ncv + ic;

                    if (ic == jc){
                        Amn = Coulomb_set.A[mm][0];
                        Dmn = Coulomb_set.D[mm][0];
                    } else{
                        Amn = Coulomb_set.A[mm][index_dia];
                        Bmn = Coulomb_set.B[mm][index_dia];
                        Cmn = Coulomb_set.C[mm][index_dia];
                        Dmn = Coulomb_set.D[mm][index_dia];
                    }




                    Coulomb_set.X_cc[mm*Ncv*Ncv + index_dia] =  (Amn* Coulomb_set.f_cc[index_dia] + Dmn* Coulomb_set.f_ss[index_dia]) - Coulomb_set.X_cc0[mm*Ncv*Ncv +index_dia];
                    if(jc != ic){ // for off diagonal part we have additional imaginary term which includes B and C coeff
                        Coulomb_set.X_cc[mm*Ncv*Ncv + index_dia] += (Bmn* Coulomb_set.f_cs[index_dia] + Cmn* Coulomb_set.f_sc[index_dia]);
                        Coulomb_set.X_cc[mm*Ncv*Ncv +index_offdia] = conj(Coulomb_set.X_cc[mm*Ncv*Ncv + index_dia]);
                    }

                    Coulomb_set.X_cs[mm*Ncv*Ncv +index_dia] =  (Amn* Coulomb_set.f_cs[index_dia] - Dmn* Coulomb_set.f_sc[index_dia]) - Coulomb_set.X_cs0[mm*Ncv*Ncv +index_dia];
                    if(jc != ic){
                        Coulomb_set.X_cs[mm*Ncv*Ncv + index_dia] += (-Bmn* Coulomb_set.f_cc[index_dia] + Cmn* Coulomb_set.f_ss[index_dia]);
                        Coulomb_set.X_cs[mm*Ncv*Ncv +index_offdia] = conj(Coulomb_set.X_cs[mm*Ncv*Ncv +index_dia]);
                    }

                    Coulomb_set.X_sc[mm*Ncv*Ncv +index_dia] =  (Amn* Coulomb_set.f_sc[index_dia] - Dmn* Coulomb_set.f_cs[index_dia]) - Coulomb_set.X_sc0[mm*Ncv*Ncv +index_dia];
                    if(jc != ic){
                        Coulomb_set.X_sc[mm*Ncv*Ncv + index_dia] += (Bmn* Coulomb_set.f_ss[index_dia] - Cmn* Coulomb_set.f_cc[index_dia]);
                        Coulomb_set.X_sc[mm*Ncv*Ncv +index_offdia] = conj(Coulomb_set.X_sc[mm*Ncv*Ncv +index_dia]);
                    }

                    Coulomb_set.X_ss[mm*Ncv*Ncv +index_dia] =  (Amn* Coulomb_set.f_ss[index_dia] + Dmn* Coulomb_set.f_cc[index_dia]) - Coulomb_set.X_ss0[mm*Ncv*Ncv +index_dia];
                    if(jc != ic){
                        Coulomb_set.X_ss[mm*Ncv*Ncv + index_dia] += (-Bmn* Coulomb_set.f_sc[index_dia] - Cmn* Coulomb_set.f_cs[index_dia]);
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


