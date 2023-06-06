Ncut = Coulomb_set.Ncut;
Coulomb_set.epsilon_static *= det_a; // \epsilon* A
// Coulomb_set.Sample_orientation = "110"; // by default in yz plane and is 011

Coulomb_set.X_cc.resize(Ncut*Ncut*Ncv*Ncv);
Coulomb_set.X_sc.resize(Ncut*Ncut*Ncv*Ncv);
Coulomb_set.X_cs.resize(Ncut*Ncut*Ncv*Ncv);
Coulomb_set.X_ss.resize(Ncut*Ncut*Ncv*Ncv);
Coulomb_set.X_cc0.resize(Ncut*Ncut*Ncv*Ncv);
Coulomb_set.X_sc0.resize(Ncut*Ncut*Ncv*Ncv);
Coulomb_set.X_cs0.resize(Ncut*Ncut*Ncv*Ncv);
Coulomb_set.X_ss0.resize(Ncut*Ncut*Ncv*Ncv);
Coulomb_set.X_ss0.fill(0.);
Coulomb_set.X_cc0.fill(0.);
Coulomb_set.X_sc0.fill(0.);
Coulomb_set.X_cs0.fill(0.);
Coulomb_set.f_cc.resize(Ncv*Ncv);
Coulomb_set.f_sc.resize(Ncv*Ncv);
Coulomb_set.f_cs.resize(Ncv*Ncv);
Coulomb_set.f_ss.resize(Ncv*Ncv);
Coulomb_set.f_ss.fill(0.);
Coulomb_set.f_cc.fill(0.);
Coulomb_set.f_sc.fill(0.);
Coulomb_set.f_cs.fill(0.);

// creating displacements:
MPI_Barrier(MPI_COMM_WORLD);
vec3x DD_C(Ncv,Ncv,3);//vec2x dx(Ncv,Ncv); //vec2x dy(Ncv,Ncv); //vec2x dz(Ncv,Ncv);
vec2x H_C(Ncv,Ncv);
vec2x U_C(Ncv, Ncv);
Coord_B kkk_C;

if (rank_ == root_rank){
	if (Coulomb_set.Sample_orientation == "011"){ // sample yz
		kkk_C.setcrys(0,0.3,0.1); // random point, we need diagonal dipole, it's constant
	} else if (Coulomb_set.Sample_orientation == "110"){ // sample in xy
		kkk_C.setcrys(0.3,0.1, 0); // random point, we need diagonal dipole, it's constant
	}
    gHUD(H_C, U_C, DD_C, kkk_C, iMode);

    cout << endl << "diagonal terms of dipole" << endl ;
    for (int ic=0; ic<Ncv; ic++){
        cout << " band number" << ic;
        for (int i_coord=0; i_coord<3; i_coord++){
            Displacement_orb[ic][i_coord] = real(DD_C[ic][ic][i_coord]);
            cout << " coord number: " << i_coord << " Displacement_orb[ic][i_coord]=" << Displacement_orb[ic][i_coord];
        }
        cout << endl;
        
    }
}

MPI_Bcast(& Displacement_orb[0][0], Ncv*3, 
    MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);





int n_offdia_terms = Ncv*Ncv; // don't want to think too much on coeff order now. Half of this points will never be used. Real dimension is int((Ncv*Ncv - Ncv) / 2) + 1; // number of off-dia terms
Coulomb_set.A.resize(Ncut*Ncut, n_offdia_terms) ; // !!!! allocate memory for coefficients
Coulomb_set.B.resize(Ncut*Ncut, n_offdia_terms); // !!!! allocate memory for coefficients
Coulomb_set.C.resize(Ncut*Ncut, n_offdia_terms) ; // !!!! allocate memory for coefficients
Coulomb_set.D.resize(Ncut*Ncut, n_offdia_terms); // !!!! allocate memory for coefficients
Coulomb_set.Screen_const.resize(4, n_offdia_terms);// [0] - without derivatives
  // [1] - d^2 /dx^2 // [2] - d^2/dy^2 // [3] - d^2/(dxdy)

// double det_sum_k = Coord_B::getJ(); // Jacobian determinant. we need it to multiply sums over k

// double border_BZ = 1.0 ; // non periodic Coulomb
double border_BZ = 0.5 ; // just crap, need to delete in next ver
double  q_cut = border_BZ - 0.02;//0.48;


Coord_B k0;
double Vbord_x = 0; //don't need it, should delete in nex ver
double Vbord_y = 0;
double V_corner = 0;


double kx_temp, ky_temp;
std::vector<float> Vq(Nk[1]*Nk[1]); // vector of z coordinates for imshow (supports only float)
std::vector<std::vector<double>> x, y, z; // matrices of coordinates for surface



MPI_Barrier(MPI_COMM_WORLD);


// Foufier series
//  Matrices of Fourier transform coefficients
int N_sum=4999; // Number of points to integrate, should be odd for Simpson
vec1d q; // allocate grid
vec1d integrWeightq;
int begin_count, end_count, lenght_q; // MPI counters for sums
getSimpson2D_MPI(N_sum, integrWeightq, q,  border_BZ,
					begin_count, end_count, lenght_q);

if (rank_ == root_rank) std::cout << "Coulomb coefficients calulation" << endl;
vec3x Wq;

int G_num = 25;
vec2d G_vec(G_num, 2); //vectors to sum over neibouring BZ
int i_G = 0;
int bord = int(sqrt(G_num)/2); // G_num = 9 -> bord = 1
for (int i_x = -bord; i_x < (bord + 1); i_x++){
    for(int i_y = -bord; i_y < (bord + 1); i_y++){
        G_vec[i_G][0] = i_x;
        G_vec[i_G][1] = i_y;
        i_G += 1;
    }
}
Wq.resize(lenght_q, Ncv, Ncv); // save Coulomb to plot it later
Wq.fill(0.0);
#pragma omp parallel
{
	int qx_global, qy_global, iq_global;
	#pragma omp  for schedule(dynamic)
	for (int iq_loc = 0; iq_loc < lenght_q; iq_loc++){// create 2d simpson array
	    iq_global = begin_count + iq_loc;
	    qx_global = iq_global / (N_sum); // integer division
	    qy_global = iq_global % (N_sum); // leftover; doesn't matter which is x and y
	    // if (qy_global == 0 and qx_global%10 == 0) cout << "qy_global " << qy_global << " qx_global " << qx_global << endl;

		Create_Wq_with_G_sum(Wq, Ncv, q, iq_loc, qx_global, qy_global,
		    Displacement_orb, G_num, 
		   Coulomb_set.Rytova_Keldysh, q_cut, Vbord_x, Vbord_y, V_corner,
		    Coulomb_set.qTF, Coulomb_set.epsilon_static, N_sum, Coulomb_set, G_vec );

	}

	#pragma omp barrier
}



MPI_Barrier(MPI_COMM_WORLD);
#pragma omp parallel
{
	vec2x AD(4, n_offdia_terms); // vector where we will put calculated ABCD results
	int FC;
	#pragma omp  for schedule(dynamic)
	for (int mm =0; mm< Ncut*Ncut; mm++ ){ //row summation over Fourier series
		int m = mm / Ncut; // integer division
		int n = mm % Ncut; // leftover
			if (((n % 5) == 0) and (rank_ == root_rank)){
				std::cout << endl;
			}
		FC= Four_integr_Coloumb(m, n, AD, N_sum, 
			q, integrWeightq, 
			begin_count, end_count, lenght_q, 
			Wq, Ncv);
		for(int i=0; i<AD.n2(); i++){ //
			Coulomb_set.A[mm][i]=real(AD[0][i]);
			Coulomb_set.B[mm][i]= c1 * imag(AD[1][i]);
			Coulomb_set.C[mm][i]= c1 * imag(AD[2][i]);
			Coulomb_set.D[mm][i]=real(AD[3][i]);
		}

		if (rank_ == root_rank){
			#pragma omp critical
			{ // execute only by one thread, others wait
				std::cout << " m: " << m << " n: " << n ;
				// std::cout << " | A(" << m << ","<< n << ")=" <<  Coulomb_set.A[n+Ncut*m];
				// std::cout << " | D(" << m<< "," << n << ")=" <<  Coulomb_set.D[n+Ncut*m];
			}
		}
	}
	#pragma omp barrier
}

MPI_Barrier(MPI_COMM_WORLD);

MPI_reduce_vec2x_safe(Coulomb_set.A, root_rank, rank_, Ncut*Ncut, n_offdia_terms);
MPI_reduce_vec2x_safe(Coulomb_set.B, root_rank, rank_, Ncut*Ncut, n_offdia_terms);
MPI_reduce_vec2x_safe(Coulomb_set.C, root_rank, rank_, Ncut*Ncut, n_offdia_terms);
MPI_reduce_vec2x_safe(Coulomb_set.D, root_rank, rank_, Ncut*Ncut, n_offdia_terms);

MPI_Barrier(MPI_COMM_WORLD);

if (rank_ == root_rank) { // calculate correction


	string A_txt, B_txt, C_txt, D_txt; 
	ofstream A_coeff, B_coeff, C_coeff, D_coeff;
	A_txt= "Output/A_coeff.txt";
	B_txt= "Output/B_coeff.txt";
	C_txt= "Output/C_coeff.txt";
	D_txt= "Output/D_coeff.txt";
	A_coeff.open(A_txt.c_str());
	B_coeff.open(B_txt.c_str());
	C_coeff.open(C_txt.c_str());
	D_coeff.open(D_txt.c_str());
	// Print A, D
	for (int mm =0; mm< Ncut*Ncut; mm++ ){ //row summation over Fourier series
		int m = mm / Ncut; // integer division
		int n = mm % Ncut; // leftover
		A_coeff  << setprecision(15) << scientific << real(Coulomb_set.A[n+Ncut*m][0]) << " " << real(Coulomb_set.A[n+Ncut*m][1]) << " " << imag(Coulomb_set.A[n+Ncut*m][1]) << endl;
		B_coeff  << setprecision(15) << scientific << real(Coulomb_set.B[n+Ncut*m][0]) << " " << real(Coulomb_set.B[n+Ncut*m][1]) << " " << imag(Coulomb_set.B[n+Ncut*m][1]) << endl;
		C_coeff  << setprecision(15) << scientific << real(Coulomb_set.C[n+Ncut*m][0]) << " " << real(Coulomb_set.C[n+Ncut*m][1]) << " " << imag(Coulomb_set.C[n+Ncut*m][1]) << endl;
		D_coeff  << setprecision(15) << scientific << real(Coulomb_set.D[n+Ncut*m][0]) << " " << real(Coulomb_set.D[n+Ncut*m][1]) << " " << imag(Coulomb_set.D[n+Ncut*m][1]) << endl;
	}
	A_coeff.close();
	B_coeff.close();
	C_coeff.close();
	D_coeff.close();



	double q_TF0 = 0; //pow(10, -10);
	Coulomb_set.Screen_const.fill(0.0);

	double border_BZ_correction = 0.5;
	N_sum=19999;
	double dk2 = 1. /(N_sum*N_sum);

	q.resize(N_sum);
	integrWeightq.resize(N_sum*N_sum);
	getSimpson2D(N_sum, integrWeightq, q, border_BZ_correction);


	#pragma omp parallel
	{ 
		double Vq0;
		complexd Vq_screened, dV_V0_Vq_screened;
		double qx, qy, q_mod, q_peak, phi;
		Coord_B k1; // !!!! THIS LINE SLOW
		vec2x sum_loc;
		sum_loc.resize(4, n_offdia_terms);
		sum_loc.fill(0.0);
		double delta_Vq0_ = 0;

		complex<double> phase_q = 1;
		int index_linear;
		q_peak = Coulomb_set.qTF;
		if (Coulomb_set.qTF < 0.02) q_peak = 0.02;

		vec3x Wq0, W_Screened_Fourier;
		Wq0.resize(1, Ncv, Ncv);// Potential without Thomas-Fermi
		W_Screened_Fourier.resize(1, Ncv, Ncv);// Potential without Thomas-Fermi
		int iq_loc0 = 0; // do not need all matrix

		#pragma omp  for schedule(dynamic) 
		for(int ik_x = N_sum/4 ; ik_x< 3* N_sum/4; ik_x++ ){
			for (int ik_y = N_sum/4; ik_y< 3* N_sum/4; ik_y++ ){

				qx = q[ik_x];
				qy = q[ik_y];

				if (Coulomb_set.Sample_orientation == "011"){ // sample yz
					k1.setcrys(0.0, qx, qy); //to create k-points in crystal coordinates for phase and Coulomb
				} else if (Coulomb_set.Sample_orientation == "110"){ // sample in xy
					k1.setcrys( qx, qy , 0.0); //to create k-points in crystal coordinates for phase and Coulomb
				}
				q_mod = k1.dot(k1);
				if(q_mod < (q_peak * q_peak) and q_mod > pow(10, -15)){
					Wq0.fill(0.0);
					
					// Potential without Thomas-Fermi:
					Create_Wq_with_G_sum(Wq0, Ncv, q, iq_loc0, ik_x, ik_y,
					    Displacement_orb, G_num, 
					   Coulomb_set.Rytova_Keldysh, q_cut, Vbord_x, Vbord_y, V_corner,
					    q_TF0, Coulomb_set.epsilon_static, N_sum, Coulomb_set, G_vec );


					for (int ic=0; ic<Ncv; ic++){// summation over bands
	            		for (int jc=ic; jc<Ncv; jc++){ // summation over bands
	            			index_linear = ic*Ncv + jc;
	            			if (ic == jc) {
	            				index_linear = 0;
	            			}
	            			Vq_screened = Fourier_series_Coloumb(qx, qy ,Ncut, 
											Coulomb_set.A, Coulomb_set.B, Coulomb_set.C, Coulomb_set.D, index_linear);
	            			dV_V0_Vq_screened = (Wq0[0][ic][jc] - Vq_screened) * dk2; 
	            			if(real(dV_V0_Vq_screened) > 0)
	            			{
		            			sum_loc[0][index_linear] += dV_V0_Vq_screened;
		            			sum_loc[1][index_linear] += dV_V0_Vq_screened* qx * (qx/2.0);
		            			sum_loc[2][index_linear] += dV_V0_Vq_screened* qx * (qx/2.0);
		            			sum_loc[3][index_linear] += dV_V0_Vq_screened* qx * qy;
		            		}

	            		}
	            	}
	            }

				

			}
		} //#pragma omp  for schedule(dynamic) 

		#pragma omp critical 
		{
			for (int i_order = 0; i_order < 4; i_order++){
				for (int i_band = 0; i_band < n_offdia_terms; i_band++){
					Coulomb_set.Screen_const[i_order][i_band] += sum_loc[i_order][i_band];
				}
			}
			
		}

		#pragma omp barrier
	}



} // end if root rank



MPI_Barrier(MPI_COMM_WORLD);

MPI_Bcast_vec2x_safe( Coulomb_set.A, root_rank, rank_, Ncut*Ncut, n_offdia_terms);
MPI_Bcast_vec2x_safe( Coulomb_set.B, root_rank, rank_, Ncut*Ncut, n_offdia_terms);
MPI_Bcast_vec2x_safe( Coulomb_set.C, root_rank, rank_, Ncut*Ncut, n_offdia_terms);
MPI_Bcast_vec2x_safe( Coulomb_set.D, root_rank, rank_, Ncut*Ncut, n_offdia_terms);

MPI_Bcast_vec2x_safe( Coulomb_set.Screen_const, root_rank, rank_, 4, n_offdia_terms);



if (rank_ == root_rank) {
	cout << "Coulomb_set.Screen_const[0][0] = " << Coulomb_set.Screen_const[0][0] << endl;
	cout << "Coulomb_set.Screen_const[1][0] = " << Coulomb_set.Screen_const[1][0] << endl;
	cout << "Coulomb_set.Screen_const[2][0] = " << Coulomb_set.Screen_const[2][0] << endl;
	cout << "Coulomb_set.Screen_const[3][0] = " << Coulomb_set.Screen_const[3][0] << endl << endl;

	cout << "Coulomb_set.Screen_const[0][1] = " << Coulomb_set.Screen_const[0][1] << endl;
	cout << "Coulomb_set.Screen_const[1][1] = " << Coulomb_set.Screen_const[1][1] << endl;
	cout << "Coulomb_set.Screen_const[2][1] = " << Coulomb_set.Screen_const[2][1] << endl;
	cout << "Coulomb_set.Screen_const[3][1] = " << Coulomb_set.Screen_const[3][1] << endl << endl;



	std::cout  << "   Coulomb_set.qTF   " << Coulomb_set.qTF << endl;
	std::cout  << "   Coulomb_set.epsilon_static   " << Coulomb_set.epsilon_static << endl;
	std::cout  << "   Coulomb_set.Coulomb_calc   " << Coulomb_set.Coulomb_calc << endl;
	if(Coulomb_set.Wannie_basis){
		std::cout  << "   Coulomb interaction defined in Wannier Basis   " << endl;
	}
	if(Coulomb_set.Diagonal_basis){
	 	std::cout  << "   Coulomb interaction defined in diagonal Basis   " << endl;
	}
	if(Coulomb_set.Rytova_Keldysh){
	 	std::cout  << "   Bare Coulomb in Rytova-Keldysh form   " << endl;

	}

}

MPI_Barrier(MPI_COMM_WORLD);
std::cout << " rank=" << rank_ << endl;
MPI_Barrier(MPI_COMM_WORLD);
