// integration wieght


// Function creates grid (-0.5, 0.5) and array of Simpson coefficients
void Grid_and_Simpson(int N_sum, vec1d  &q, vec1d  &sp_arr, double border_BZ){
	double h=  2*border_BZ /(N_sum-1); // step of integration
	for(int i=0; i<N_sum; i++){  q[i] = -border_BZ + i*h; }// create grid

	double h_simp = (3*(N_sum-1)); // simpson weight factor
	sp_arr[0]=1/h_simp;

	for(int i=1; i<(N_sum); i=i+2){ sp_arr[i] = 4.0/h_simp; }
	for(int i=2; i<(N_sum); i=i+2){ sp_arr[i] = 2.0/h_simp; } // Simson array is [1,4,2,4,2,4 ....2,4,1]/h_simp
	sp_arr[N_sum-1]=1/h_simp;
}




void getSimpson2D(int N_sum, vec1d& simpWeight, vec1d  &q,  double border_BZ){
// create 2D simpson grid on a N_sum*N_sum lattice
	vec1d sp_arr_k(N_sum);     // allocate auxilary 1d simpson array
	Grid_and_Simpson(N_sum, q, sp_arr_k, border_BZ); // fill grid and 1d simpson array
	for (int kx = 0; kx < N_sum; kx++){ // create 2d simpson array
		for (int ky = 0; ky < N_sum; ky++){
			simpWeight[ky +  N_sum*kx] = sp_arr_k[kx] * sp_arr_k[ky];
		}
	}
}



// simpson integration coefficients on a grid
void getSimpson2D_MPI(int N_sum, vec1d& simpWeight, vec1d  &q,  double border_BZ, 
	int & begin_count, int & end_count, int & lenght_q){
	// create 2D simpson grid on a N_sum*N_sum lattice
	int N_sum2D = N_sum*N_sum;
	lenght_q = int(N_sum2D/num_procs);
	int leftover = N_sum2D % num_procs;

	end_count = lenght_q *rank_ +  lenght_q; // MPI 
	begin_count = lenght_q *rank_;


	if (rank_ < leftover){
	  begin_count += rank_;
	  end_count += rank_ + 1;
	}
	else{
	  begin_count += leftover;
	  end_count += leftover;
	}

	lenght_q = end_count - begin_count;
	q.resize(N_sum);
	simpWeight.resize(lenght_q);


	vec1d sp_arr_k(N_sum);     // allocate auxilary 1d simpson array
	Grid_and_Simpson(N_sum, q, sp_arr_k, border_BZ); // fill grid and 1d simpson array


	
	int iq_global, qx_global, qy_global;
	for (int iq_loc = 0; iq_loc < lenght_q; iq_loc++){// create 2d simpson array
		iq_global = begin_count + iq_loc;
		qx_global = iq_global / N_sum; // integer division
		qy_global = iq_global % N_sum; // leftover; doesn't matter which is x and y

		simpWeight[iq_loc] = sp_arr_k[qx_global] * sp_arr_k[qy_global];

	}



}



void get2DPrismInterrationCoeff(int N_sum, vec1d& coeffMatrix){
	// create quasi-prism approximation coefficients, meaning that elementary
	// volume is 0.5dkx * dky * (Fij + F(i+1)j + Fi(j+1))/3
	// triangle 0.5*dkx * dky and medium point inside it
	double h= 1.0 /(N_sum-1); // step of integration
	h = h*h;
	for(int kx=0; kx<N_sum; kx++){
		for(int ky=0; ky<(N_sum); ky++){
			coeffMatrix[ky +  N_sum*kx] = h;
			if (kx ==0 || ky==0 || kx == (N_sum - 1) || ky == (N_sum - 1)){
				// borders
				coeffMatrix[ky +  N_sum*kx] = h/3;
			}
			if (kx ==0 && ky==0){ // corners
				coeffMatrix[ky +  N_sum*kx] = h/6;
			}
			if (kx ==0 && ky==(N_sum - 1)){ // corners
				coeffMatrix[ky +  N_sum*kx] = h/6;
			}
			if (kx ==(N_sum - 1) && ky==0){ // corners
				coeffMatrix[ky +  N_sum*kx] = h/6;
			}
			if (kx ==(N_sum - 1) && ky==(N_sum - 1)){ // corners
				coeffMatrix[ky +  N_sum*kx] = h/6;
			}

		}
	}
}


void get2DPrismInterrationCoeff_MPI(int N_sum,  
 vec2d& kpt, vec1d& coeffMatrix){
  // create quasi-prism approximation coefficients, meaning that elementary
  // volume is 0.5dkx * dky * (Fij + F(i+1)j + Fi(j+1))/3
  // triangle 0.5*dkx * dky and medium point inside it
  double h= 1./(N_sum-1); // step of integration
  h = h*h;
  double eps = h/1000;
  double dkx, dky;
  for(int k = 0; k< coeffMatrix.n1(); k++){
	dkx = abs(abs(kpt[k][1]) - 0.5);
	dky = abs(abs(kpt[k][2]) - 0.5);
	  coeffMatrix[k] = h;

	  if (dkx < eps || dky < eps){
		// borders
		// coef 2 besause we have periodic grid
		coeffMatrix[k] = 2 * h/3;
	  }
	  if (dkx < eps && dky < eps){ // corners
		// coef 4 besause we have periodic grid
		coeffMatrix[k] = 4 * h/6;
	  }
	  
	
  }
}

void get_trig_coef(trig_coefficients & trig, vec2d & k, 
  int Ncut, int n )
{
  trig.cos_mkx.resize(n*Ncut);
  trig.sin_mkx.resize(n*Ncut);
  trig.cos_nky.resize(n*Ncut);
  trig.sin_nky.resize(n*Ncut);


  #pragma omp parallel // num_threads(3)
  {

	double mkx, nky;
	#pragma omp for schedule(dynamic)
	for (int m =0; m< Ncut; m++){

	  for (int ik=0; ik<n; ik++){ //

        if (trig.Sample_orientation == "011"){ // sample yz
					mkx=2* M_PI * m * (k[ik][1]); // 
					nky=2* M_PI * m * (k[ik][2]);
        } else if (trig.Sample_orientation == "110"){ // sample in xy
 					mkx=2* M_PI * m * (k[ik][0]); // 
					nky=2* M_PI * m * (k[ik][1]);
        }


		trig.cos_mkx[m*n + ik] = cos(mkx);
		trig.sin_mkx[m*n + ik] = sin(mkx);
		trig.cos_nky[m*n + ik] = cos(nky);
		trig.sin_nky[m*n + ik] = sin(nky);

	  }
	}
  }

}




void print_max_min_conduct(vec3x & P, 
  Private_omp_parameters& OMP_private, 
  Coulomb_parameters& Coulomb_set, int root_rank){
  int Ncv = P.n2();
  int ik; // position in shared variable
  double maxP = -pow(10,17);
  double minP = pow(10,17);
  double P_curr;
  double max_MPI= 0;
  double min_MPI= 0;

  #pragma omp master
	{  
	  Coulomb_set.max = maxP; // shared variable
	  Coulomb_set.min = minP;
	}



  for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++)    {
	ik = ik_pr + OMP_private.begin_count;
	P_curr = real(P[ik][Ncv-1][Ncv-1]);
	if (P_curr >  maxP){
	  maxP = P_curr;
	}
	if (P_curr <  minP){
	  minP = P_curr;
	}

  }

  #pragma omp barrier
  #pragma omp critical
  {
	if (maxP > Coulomb_set.max){
	  Coulomb_set.max = maxP;
	}
	if (minP < Coulomb_set.min){
	  Coulomb_set.min = minP;
	}
  }


  #pragma omp barrier
 
  #pragma omp master
  {
	 MPI_Barrier(MPI_COMM_WORLD);

	MPI_Reduce(&Coulomb_set.max, & max_MPI, 1, MPI_DOUBLE, MPI_MAX, root_rank, MPI_COMM_WORLD);
	MPI_Reduce(&Coulomb_set.min, & min_MPI, 1, MPI_DOUBLE, MPI_MIN, root_rank, MPI_COMM_WORLD);

	if (rank_ == root_rank){

	  cout << " max: " << max_MPI << " min: " << min_MPI << endl;

	}
  }
  #pragma omp barrier
  

}


