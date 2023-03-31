// double det_sum_k = Coord_B::getJ(); // Jacobian determinant. we need it to multiply sums over k

Coulomb_set.det_J_for_k_sum = 1.0;// Coord_B::getJ();
if (rank_ == root_rank) { 



// double border_BZ = 1.0 ; // non periodic Coulomb
double border_BZ = 0.5 ; // periodic Coulomb
double  q_cut = border_BZ - 0.02;//0.48;
 
double Vbord_x = Coloumb_THh_F(-border_BZ, 0.0, Coulomb_set.qTF, 
  Coulomb_set.epsilon_static,  2*border_BZ, 0.0, 0.0, 0.0, 
  Coulomb_set.Rytova_Keldysh);
double Vbord_y = Coloumb_THh_F(0.0, -border_BZ, Coulomb_set.qTF, 
  Coulomb_set.epsilon_static,  2*border_BZ, 0.0, 0.0, 0.0, 
  Coulomb_set.Rytova_Keldysh);
double V_corner = Coloumb_THh_F(border_BZ, border_BZ, Coulomb_set.qTF, 
  Coulomb_set.epsilon_static,  2*border_BZ, 0.0, 0.0, 0.0, 
  Coulomb_set.Rytova_Keldysh);


  std::vector<std::vector<double>> x, y, z; // matrices of coordinates for surface
  std::vector<float> Vq(Nk[1]*Nk[1]); // vector of z coordinates for imshow (supports only float)
  double kx_temp, ky_temp;
  for (int i = 0; i < Nk[1];  i++) {
      std::vector<double> x_row, y_row, z_row;
      for (int j = 0; j < Nk[1];  j++) {
        kx_temp = kpt[i+ Nk[1]*j][1];
        ky_temp = kpt[i+ Nk[1]*j][2];
        x_row.push_back(kx_temp);
        y_row.push_back(ky_temp);
        Vq.at(Nk[1] * j + i) = Coloumb_THh_F(kx_temp,ky_temp, Coulomb_set.qTF, Coulomb_set.epsilon_static, 1.0, 0.0, 0.0, 0.0, Coulomb_set.Rytova_Keldysh); // calculate Thomas-Fermi Coloumb in every point of BZ
      }
      x.push_back(x_row);
      y.push_back(y_row);
  }




  // Foufier series
  //  Matrices of Fourier transform coefficients
  int N_sum=2999; // Number of points to integrate, should be odd for Simpson
  vec1d q(N_sum); // allocate grid
  vec1d integrWeightq(N_sum*N_sum);
  getSimpson2D(N_sum, integrWeightq, q,  border_BZ);


  #pragma omp parallel
  {
    std::vector<double> AD(2, 0); // vector where we will put calculated ABCD results
    int FC;
  #pragma omp  for schedule(dynamic)
  for (int mm =0; mm< Ncut*Ncut; mm++ ){ //row summation over Fourier series
      int m = mm / Ncut; // integer division
      int n = mm % Ncut; // leftover
      if ((n % 5) == 0){
         std::cout << endl;
      }
        FC= Four_integr_Coloumb(m, n, AD,N_sum, q, integrWeightq, 
          Coulomb_set.qTF, Coulomb_set.epsilon_static, q_cut, Vbord_x, Vbord_y, V_corner, Coulomb_set.Rytova_Keldysh);
        Coulomb_set.A[mm]=AD[0];
        Coulomb_set.D[mm]=AD[1];
        #pragma omp critical
        { // execute only by one threads, others wait
          // std::cout << " m: " << m << " n: " << n ;
          std::cout << " | A(" << m << ","<< n << ")=" <<  Coulomb_set.A[n+Ncut*m];
          std::cout << " | D(" << m<< "," << n << ")=" <<  Coulomb_set.D[n+Ncut*m];
        }
      }
    
  }

  //Creating Fourier transformed function
  double errorFourieTrans = 0.0;
  double maxCoulomb = -10.0;
  double k_mod;
  Coord_B k2; // !!!! THIS LINE SLOW
  std::vector<float> VqF(Nk[1]*Nk[1]); // vector of z coordinates for imshow (supports only float)
  std::vector<float> Delta_V(Nk[1]*Nk[1]); //Difference between approximation and Thomas Fermi to estimate error
  for (int i = 0; i < Nk[1];  i++) {
      for (int j = 0; j < Nk[1];  j++) {
         kx_temp = kpt[i+ Nk[1]*j][1];
         ky_temp = kpt[i+ Nk[1]*j][2];
         k2.setcrys(0.0, kx_temp, ky_temp);
         k_mod = k2.dot(k2);
          VqF.at(Nk[1] * j + i) = Fourier_series_Coloumb(kx_temp,ky_temp,Ncut,
                                                Coulomb_set.A, Coulomb_set.D);
          Delta_V.at(Nk[1] * j + i) = VqF.at(Nk[1] * j + i) - Vq.at(Nk[1] * j + i);
          double errCurij = abs(Delta_V[Nk[1] * j + i]);

          if ( (errCurij >  errorFourieTrans) and (k_mod > (Coulomb_set.qTF * Coulomb_set.qTF))){
            errorFourieTrans = errCurij;
          }

       }
   }
   std::cout << endl<< endl<< endl << " maximum error of Fourier: " << errorFourieTrans  << endl << endl<< endl;


system("mkdir Output");
string A_txt, D_txt; 
ofstream A_coeff, D_coeff;
A_txt= "Output/A_coeff.txt";
D_txt= "Output/D_coeff.txt";
A_coeff.open(A_txt.c_str());
D_coeff.open(D_txt.c_str());
// Print A, D
for (int mm =0; mm< Ncut*Ncut; mm++ ){ //row summation over Fourier series
  int m = mm / Ncut; // integer division
  int n = mm % Ncut; // leftover
  A_coeff  << setw(25) << setprecision(15) << scientific << Coulomb_set.A[n+Ncut*m] << endl;
  D_coeff  << setw(25) << setprecision(15) << scientific << Coulomb_set.D[n+Ncut*m] << endl;
}




double q_TF0 = 0; //pow(10, -10);
Coulomb_set.Screen_const.fill(0.0);

double border_BZ_correction = 0.5;
N_sum=29999;
q.resize(N_sum);
integrWeightq.resize(N_sum*N_sum);
getSimpson2D(N_sum, integrWeightq, q, border_BZ_correction);


#pragma omp parallel
{ 
  double Vq0;
  double Vq_screened;
  double qx, qy, q_mod, q_peak;
  Coord_B k1; // !!!! THIS LINE SLOW
  double sum_loc = 0;
  double sum_loc0 = 0;
  double sum_loc1 = 0;
  double sum_loc2 = 0;
  double sum_loc3 = 0;

  if (Coulomb_set.qTF < 0.02) q_peak = 0.02;


   #pragma omp  for schedule(dynamic) 
   for(int ik_x = N_sum/4 ; ik_x< 3* N_sum/4; ik_x++ ){
    for (int ik_y = N_sum/4; ik_y< 3* N_sum/4; ik_y++ ){
      qx = q[ik_x];
      qy = q[ik_y];
      k1.setcrys(0.0, qx, qy);
      q_mod = k1.dot(k1);
      if(q_mod < (q_peak * q_peak) and q_mod > pow(10, -15)){
        Vq0 = Coloumb_THh_F(qx, qy, q_TF0, 
                      Coulomb_set.epsilon_static, 1.0, 0.0, 0.0, 0.0, Coulomb_set.Rytova_Keldysh); // calculate Thomas-Fermi Coloumb in every point of BZ
      
        Vq_screened = Fourier_series_Coloumb(qx, qy ,Ncut, 
                                          Coulomb_set.A, Coulomb_set.D);

        // sum_loc += integrWeightq[ik_x + N_sum*ik_y] * (Vq0 - Vq_screened); 
        sum_loc = (Vq0 - Vq_screened)/(N_sum*N_sum); 
        if(sum_loc >0){ // only when it's peak higher than should be
          sum_loc0 += sum_loc;
          sum_loc1 += sum_loc * qx * qx/2;
          sum_loc2 += sum_loc * qy * qy/2;
          sum_loc3 += sum_loc * qx * qy;
        }

      }

    }

   } //#pragma omp  for schedule(dynamic) 

   #pragma omp critical 
   {
    Coulomb_set.Screen_const[0] += Coulomb_set.det_J_for_k_sum * sum_loc0;
    Coulomb_set.Screen_const[1] += Coulomb_set.det_J_for_k_sum * sum_loc1;
    Coulomb_set.Screen_const[2] += Coulomb_set.det_J_for_k_sum * sum_loc2;
    Coulomb_set.Screen_const[3] += Coulomb_set.det_J_for_k_sum * sum_loc3;
   }


 }
  cout << "Coulomb_set.Screen_const[0] = " << Coulomb_set.Screen_const[0] << endl;
  cout << "Coulomb_set.Screen_const[1] = " << Coulomb_set.Screen_const[1] << endl;
  cout << "Coulomb_set.Screen_const[2] = " << Coulomb_set.Screen_const[2] << endl;
  cout << "Coulomb_set.Screen_const[3] = " << Coulomb_set.Screen_const[3] << endl << endl;




  // matplotlib to draw results
  // C++ library from here https://github.com/lava/matplotlib-cpp

  //
  


  // bool show =true; // show image in interactive mode
  // string Title("Thomas-Fermi potential");
  // string File_name("ColoumbPictures/TF.png");
  // Imshow_matrix(Vq, Nk[1], Nk[1], Title, File_name, show);

  // string Title1("Fourier_series_Coloumb");
  // string File_name1("ColoumbPictures/Fourier_series_Coloumb.png");
  // Imshow_matrix(VqF, Nk[1], Nk[1], Title1, File_name1, show);


  // string Title2("Delta_V");
  // string File_name2("ColoumbPictures/Delta_V.png");
  // Imshow_matrix(Delta_V, Nk[1], Nk[1], Title2, File_name2, show);


} // end if root rank

MPI_Barrier(MPI_COMM_WORLD);



MPI_Bcast(& Coulomb_set.A[0], Ncut*Ncut, 
  MPI_DOUBLE, root_rank, MPI_COMM_WORLD);

MPI_Bcast(& Coulomb_set.D[0], Ncut*Ncut, 
  MPI_DOUBLE, root_rank, MPI_COMM_WORLD);

MPI_Bcast(& Coulomb_set.Screen_const[0], 4, 
  MPI_DOUBLE, root_rank, MPI_COMM_WORLD);



MPI_Barrier(MPI_COMM_WORLD);
for (int i = 0; i < num_procs; i++){
  MPI_Barrier(MPI_COMM_WORLD);
    // Cheking kpt borders
    if (rank_ == i){
      cout << " ******************** "  << endl;
        cout << " my_rank " <<  rank_ << endl;
        cout << " Coulomb_set.D[3] " << Coulomb_set.D[3] << " Coulomb_set.A[3] " << Coulomb_set.A[3]   << endl << endl << endl; 
    }

  MPI_Barrier(MPI_COMM_WORLD);
}
MPI_Barrier(MPI_COMM_WORLD);


