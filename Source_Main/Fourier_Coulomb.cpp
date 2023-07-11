Ncut = Coulomb_set.Ncut;
Coulomb_set.epsilon_static *= det_a; // permitivity with unit cell area \epsilon* A


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


Coulomb_set.V_Hartree.resize(Ncv, Ncv); // Hartree term correction \sum e^{G(ts - tm)} V_{RK}(G)
Coulomb_set.E_Hartree.resize(Ncv); // Hartree term correction to diagonal term \sum V_Hartree(m,s) f_{s}(t) 
Coulomb_set.E_Hartree0.resize(Ncv); // Hartree term correction to diagonal term \sum V_Hartree(m,s) f_{s}(t) 

Coulomb_set.f_ss.fill(0.);
Coulomb_set.f_cc.fill(0.);
Coulomb_set.f_sc.fill(0.);
Coulomb_set.f_cs.fill(0.);


Coulomb_set.V_Hartree.fill(0.);
Coulomb_set.E_Hartree.fill(0.);
Coulomb_set.E_Hartree0.fill(0.);

// creating displacements:
MPI_Barrier(MPI_COMM_WORLD);
vec3x DD_C(Ncv,Ncv,3);
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

MPI_Barrier(MPI_COMM_WORLD);
MPI_Bcast(& Displacement_orb[0][0], Ncv*3, 
    MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);





int n_offdia_terms = Ncv*Ncv; // don't want to think too much on coeff order now. Half of this points will never be used. Real dimension is int((Ncv*Ncv - Ncv) / 2) + 1; // number of off-dia terms
Coulomb_set.A.resize(Ncut*Ncut, n_offdia_terms) ; // !!!! allocate memory for coefficients
Coulomb_set.B.resize(Ncut*Ncut, n_offdia_terms); // !!!! allocate memory for coefficients
Coulomb_set.C.resize(Ncut*Ncut, n_offdia_terms) ; // !!!! allocate memory for coefficients
Coulomb_set.D.resize(Ncut*Ncut, n_offdia_terms); // !!!! allocate memory for coefficients
Coulomb_set.A.fill(0.0);
Coulomb_set.B.fill(0.0);
Coulomb_set.C.fill(0.0);
Coulomb_set.D.fill(0.0);

Coulomb_set.Screen_const.resize(4, n_offdia_terms);// [0] - without derivatives
  // [1] - d^2 /dx^2 // [2] - d^2/dy^2 // [3] - d^2/(dxdy)

// double border_BZ = 1.0 ; // non periodic Coulomb
double border_BZ = 0.5 ; // periodic


MPI_Barrier(MPI_COMM_WORLD);

if (Coulomb_set.Calculate){


	// Foufier series to compute
	//  Matrices of Fourier transform coefficients with fine precision
	int N_sum=4999; // Number of points to integrate, should be odd for Simpson
	// int N_sum=39999; // Number of points to integrate, should be odd for Simpson

	vec1d q; // allocate grid
	vec1d integrWeightq;
	int begin_count, end_count, lenght_q; // MPI counters for sums
	getSimpson2D_MPI(N_sum, integrWeightq, q,  border_BZ,
						begin_count, end_count, lenght_q);

	
	if (rank_ == root_rank) std::cout << "Coulomb coefficients calulation" << endl;
	cout << "rank_=" << rank_ << " begin_count=" << begin_count << " end_count=" << end_count << " lenght_q=" << lenght_q << endl;

	int G_num = 9999*9999;
	vec2d G_vec_temp(G_num, 2); //vectors to sum over neibouring BZ
	Coulomb_set.G_distance *= sqrt(det_b); // shell size in cartesian coordinates where we sum, here det_b -  is an area of unit cell; by default it was 50
	double norm_G;
	Coord_B k1;
	int i_G = 0;
	int bord = int(sqrt(G_num)/2); // G_num = 9 -> bord = 1
	for (int i_x = -bord; i_x < (bord + 1); i_x++){
		for(int i_y = -bord; i_y < (bord + 1); i_y++){
	        if (Coulomb_set.Sample_orientation == "011"){ // sample yz
	            k1.setcrys(0.0, i_x, i_y); //to create k-points in crystal coordinates for phase and Coulomb
	        } else if (Coulomb_set.Sample_orientation == "110"){ // sample in xy
	            k1.setcrys(i_x, i_y , 0.0); //to create k-points in crystal coordinates for phase and Coulomb
	        }
	        norm_G = sqrt( k1.cart[0]*k1.cart[0] + k1.cart[1]*k1.cart[1] + k1.cart[2]*k1.cart[2]);
	        if (norm_G < Coulomb_set.G_distance){
				G_vec_temp[i_G][0] = i_x;
				G_vec_temp[i_G][1] = i_y;
				i_G += 1;
	        }
		}
	}

	G_num = i_G;

	if (rank_ == root_rank){
		cout << "G_num = " << G_num << " G_distance = " << Coulomb_set.G_distance << endl;
	}


	double Re_V_OMP_max_eV = -1e10;
	double Im_V_OMP_max_eV = -1e10;
	double Re_V_OMP_min_eV = 1e10;
	double Im_V_OMP_min_eV = 1e10;

	string label_Coulomb; 
	ofstream V_Hartree_stream;
	if (rank_ == root_rank){
		system("mkdir Output/Coulomb");

		stringstream sname;
        sname.seekp(0,ios::beg); 
        sname << G_num << "_Ncut" << Ncut;
		label_Coulomb = "G_num" + sname.str() + Coulomb_set.labelInput;
	}

	#pragma omp parallel
	{
		vec2d G_vec(G_num, 2); //vectors to sum over neibouring BZ
		for (int i_G0 = 0; i_G0 < G_num; i_G0++){
			G_vec[i_G0][0] = G_vec_temp[i_G0][0];
			G_vec[i_G0][1] = G_vec_temp[i_G0][1];
		}

		Coord_B k_local;

		int qx_global, qy_global, iq_global;
		double qx, qy, spxy, mq, nq, qTF0;


		#pragma omp master // Hartree term correction \sum e^{G(ts - tm)} V_{RK}(G)
		{
			qx = 0.0;
			qy = 0.0;
			qTF0 = Coulomb_set.qTF;
			Coulomb_set.Hartree_calculation = true;
			Create_Wq_with_G_sum(Coulomb_set.V_Hartree, Ncv, qx, qy, 
			    Displacement_orb, G_num, 
			    qTF0, Coulomb_set, G_vec, k_local, Re_V_OMP_max_eV, Im_V_OMP_max_eV, Re_V_OMP_min_eV, Im_V_OMP_min_eV);
			for (int ic=0; ic<Ncv; ic++){// summation over bands
				if (ic > 0) Coulomb_set.V_Hartree[ic][ic] = Coulomb_set.V_Hartree[0][0];
				for (int jc=ic; jc<Ncv; jc++){ // summation over bands
					
					if (rank_ == root_rank){
						cout << "ic= " << ic << "; jc= " << jc << " real(V_Hartree[ic][jc])= " << real(Coulomb_set.V_Hartree[ic][jc])*energy_au_eV<< " eV; imag(V_Hartree[ic][jc])= " << imag(Coulomb_set.V_Hartree[ic][jc])* energy_au_eV << " eV" << endl;
					}
					Coulomb_set.V_Hartree[jc][ic] = real(Coulomb_set.V_Hartree[ic][jc]);
					Coulomb_set.V_Hartree[ic][jc] = Coulomb_set.V_Hartree[jc][ic];
				}
			}
			Coulomb_set.Hartree_calculation = false;
			if (rank_ == root_rank){

			
				string V_Hartree_txt; 
				ofstream V_Hartree_stream;
				
				V_Hartree_txt= "Output/Coulomb/V_Hartree_" + label_Coulomb + ".txt";
				V_Hartree_stream.open(V_Hartree_txt.c_str());
				for (int ic=0; ic<Ncv; ic++){// summation over bands
					for (int jc=0; jc<Ncv; jc++){
						V_Hartree_stream << setprecision(20) << real(Coulomb_set.V_Hartree[ic][jc]) << "    ";
					}
					V_Hartree_stream << endl;
				}
				V_Hartree_stream.close();
			}

		}



		vec2x Wq;
		Wq.resize(Ncv, Ncv); // Coulomb for each band in point (qx, qy) to add it to the sum
		vec2x AD(4, n_offdia_terms); // vector where we will put calculated ABCD results
		vec2x A_local, B_local, C_local, D_local;
		A_local.resize(Ncut*Ncut, n_offdia_terms) ; // !!!! allocate memory for coefficients
		B_local.resize(Ncut*Ncut, n_offdia_terms) ; // !!!! allocate memory for coefficients
		C_local.resize(Ncut*Ncut, n_offdia_terms) ; // !!!! allocate memory for coefficients
		D_local.resize(Ncut*Ncut, n_offdia_terms) ; // !!!! allocate memory for coefficients
		A_local.fill(0.0);
		B_local.fill(0.0);
		C_local.fill(0.0);
		D_local.fill(0.0);
		double Re_Vq_max_eV, Im_Vq_max_eV, Re_Vq_min_eV, Im_Vq_min_eV;
	    double Re_VBZ_max_eV = -1e10; // maximum in q vectors
	    double Im_VBZ_max_eV = -1e10;
	    double Re_VBZ_min_eV = 1e10;
	    double Im_VBZ_min_eV = 1e10;
		#pragma omp barrier
		#pragma omp  for schedule(dynamic)
		for (int iq_loc = 0; iq_loc < lenght_q; iq_loc++){// create 2d simpson array
		    iq_global = begin_count + iq_loc;
		    qx_global = iq_global / (N_sum); // integer division
		    qy_global = iq_global % (N_sum); // leftover; 
		   
		    
			spxy=integrWeightq[iq_loc];
			qx = q[qx_global];
			qy = q[qy_global];

			Create_Wq_with_G_sum(Wq, Ncv, qx, qy, 
			    Displacement_orb, G_num, 
			    Coulomb_set.qTF, Coulomb_set, G_vec, k_local, Re_Vq_max_eV, Im_Vq_max_eV, Re_Vq_min_eV, Im_Vq_min_eV);

	        if(Re_Vq_max_eV > Re_VBZ_max_eV) Re_VBZ_max_eV = Re_Vq_max_eV;
	        if(Im_Vq_max_eV > Im_VBZ_max_eV) Im_VBZ_max_eV = Im_Vq_max_eV;

	        if(Re_Vq_min_eV < Re_VBZ_min_eV) Re_VBZ_min_eV = Re_Vq_min_eV;
	        if(Im_Vq_min_eV < Im_VBZ_min_eV) Im_VBZ_min_eV = Im_Vq_min_eV;

			for (int mm =0; mm< Ncut*Ncut; mm++ ){ //row summation over Fourier series
				int m = mm / Ncut; // integer division
				int n = mm % Ncut; // leftover

			    mq=2* M_PI * m * qx;
	    		nq=2* M_PI * n * qy;
				Four_integr_Coloumb(mq, nq, AD, N_sum, 
					spxy, Wq, Ncv, m, n);

				for(int i=0; i<n_offdia_terms; i++){ //
					A_local[mm][i] +=real(AD[0][i]);
					B_local[mm][i] += c1 * imag(AD[1][i]);
					C_local[mm][i] += c1 * imag(AD[2][i]);
					D_local[mm][i] +=real(AD[3][i]);
				}
			}
			
			if ((iq_loc % (lenght_q/5) == 0)){
				cout << "rank_=" << rank_ << " done " << 100.0 * (iq_loc + 0.0)/(lenght_q + 0.0) << " %" << " iq_loc=" << iq_loc << endl;
			}
		} // end of omp  for schedule(dynamic)

	    #pragma omp barrier
	    // summing local ABCD into shared ABCD
	    #pragma omp critical
	    { // summation over threads
	    	for (int mm =0; mm< Ncut*Ncut; mm++ ){ //row summation over Fourier series
	    		for(int i=0; i<n_offdia_terms; i++){ //
					Coulomb_set.A[mm][i] += A_local[mm][i];
					Coulomb_set.B[mm][i] += B_local[mm][i];
					Coulomb_set.C[mm][i] += C_local[mm][i];
					Coulomb_set.D[mm][i] += D_local[mm][i];



				}
			}
			if(Re_VBZ_max_eV > Re_V_OMP_max_eV) Re_V_OMP_max_eV = Re_VBZ_max_eV;
			if(Im_VBZ_max_eV > Im_V_OMP_max_eV) Im_V_OMP_max_eV = Im_VBZ_max_eV;

			if(Re_VBZ_min_eV < Re_V_OMP_min_eV) Re_V_OMP_min_eV = Re_VBZ_min_eV;
			if(Im_VBZ_min_eV < Im_V_OMP_min_eV) Im_V_OMP_min_eV = Im_VBZ_min_eV;
	    }


		#pragma omp barrier
		// if (rank_ == root_rank) cout << "Re_V(q+G)_min_eV = " << Re_V_OMP_min_eV << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_reduce_vec2x_safe(Coulomb_set.A, root_rank, rank_, Ncut*Ncut, n_offdia_terms);
	MPI_reduce_vec2x_safe(Coulomb_set.B, root_rank, rank_, Ncut*Ncut, n_offdia_terms);
	MPI_reduce_vec2x_safe(Coulomb_set.C, root_rank, rank_, Ncut*Ncut, n_offdia_terms);
	MPI_reduce_vec2x_safe(Coulomb_set.D, root_rank, rank_, Ncut*Ncut, n_offdia_terms);

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast_vec2x_safe( Coulomb_set.A, root_rank, rank_, Ncut*Ncut, n_offdia_terms);
	MPI_Bcast_vec2x_safe( Coulomb_set.B, root_rank, rank_, Ncut*Ncut, n_offdia_terms);
	MPI_Bcast_vec2x_safe( Coulomb_set.C, root_rank, rank_, Ncut*Ncut, n_offdia_terms);
	MPI_Bcast_vec2x_safe( Coulomb_set.D, root_rank, rank_, Ncut*Ncut, n_offdia_terms);

	MPI_Barrier(MPI_COMM_WORLD);

	double V_buff;
	MPI_Reduce(&Re_V_OMP_max_eV, &V_buff, 1, MPI_DOUBLE, MPI_MAX, root_rank, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank_ == root_rank) cout << "|Re_V(q+G)_max|_eV = " << V_buff << endl;
	MPI_Reduce(&Im_V_OMP_max_eV, &V_buff, 1, MPI_DOUBLE, MPI_MAX, root_rank, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank_ == root_rank) cout << "|Im_V(q+G)_max|_eV = " << V_buff << endl;
	MPI_Reduce(&Re_V_OMP_min_eV, &V_buff, 1, MPI_DOUBLE, MPI_MIN, root_rank, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank_ == root_rank) cout << "|Re_V(q+G)_min|_eV = "  << V_buff << endl;
	MPI_Reduce(&Im_V_OMP_min_eV, &V_buff, 1, MPI_DOUBLE, MPI_MIN, root_rank, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank_ == root_rank) cout << "|Im_V(q+G)_min|_eV = "  << V_buff << endl;




	if (rank_ == root_rank) { 

		std::cout << " Coulomb coefficients A,B,C,D calculated " << endl;

		string A_txt, B_txt, C_txt, D_txt; 
		ofstream A_coeff, B_coeff, C_coeff, D_coeff;
		A_txt= "Output/Coulomb/A_coeff" + label_Coulomb + ".txt";
		B_txt= "Output/Coulomb/B_coeff" + label_Coulomb + ".txt";
		C_txt= "Output/Coulomb/C_coeff" + label_Coulomb + ".txt";
		D_txt= "Output/Coulomb/D_coeff" + label_Coulomb + ".txt";
		A_coeff.open(A_txt.c_str());
		B_coeff.open(B_txt.c_str());
		C_coeff.open(C_txt.c_str());
		D_coeff.open(D_txt.c_str());
		// Print A, D
		for (int mm =0; mm< Ncut*Ncut; mm++ ){ //row summation over Fourier series
			int m = mm / Ncut; // integer division
			int n = mm % Ncut; // leftover
			for (int i_band = 0; i_band < n_offdia_terms; i_band++){
				A_coeff  << setprecision(20) << scientific << real(Coulomb_set.A[mm][i_band]) << "    " << imag(Coulomb_set.A[mm][i_band]) << "    " ;
				B_coeff  << setprecision(20) << scientific << real(Coulomb_set.B[mm][i_band]) << "    " << imag(Coulomb_set.B[mm][i_band]) << "    " ;
				C_coeff  << setprecision(20) << scientific << real(Coulomb_set.C[mm][i_band]) << "    " << imag(Coulomb_set.C[mm][i_band]) << "    " ;
				D_coeff  << setprecision(20) << scientific << real(Coulomb_set.D[mm][i_band]) << "    " << imag(Coulomb_set.D[mm][i_band]) << "    " ;
			}
			A_coeff  << endl;
			B_coeff  << endl;
			C_coeff  << endl;
			D_coeff  << endl;
		}
		A_coeff.close();
		B_coeff.close();
		C_coeff.close();
		D_coeff.close();

	}

	double q_TF0 = 0; //pow(10, -10);
	Coulomb_set.Screen_const.fill(0.0);

	double border_BZ_correction = 0.5;
	N_sum=39999;
	double dk2 = 1. /(N_sum*N_sum);


	getSimpson2D_MPI(N_sum, integrWeightq, q,  border_BZ,
						begin_count, end_count, lenght_q);

	#pragma omp parallel
	{ // calculate correction

		vec2d G_vec(G_num, 2); //vectors to sum over neibouring BZ
		for (int i_G0 = 0; i_G0 < G_num; i_G0++){
			G_vec[i_G0][0] = G_vec_temp[i_G0][0];
			G_vec[i_G0][1] = G_vec_temp[i_G0][1];
		}
		double Re_Vq_max_eV, Im_Vq_max_eV, Re_Vq_min_eV, Im_Vq_min_eV;

		complexd Vq_screened, dV_V0_Vq_screened;
		double qx, qy, q_mod, q_peak, phi;
		int iq_global, qx_global, qy_global;
		Coord_B k_local; // 
		vec2x sum_loc;
		sum_loc.resize(4, n_offdia_terms);
		sum_loc.fill(0.0);
		double delta_Vq0_ = 0;

		complex<double> phase_q = 1;
		int index_linear;
		q_peak = Coulomb_set.qTF;
		if (Coulomb_set.qTF < 0.02) q_peak = 0.02; // area of integration around q=0

		vec2x Wq0;
		Wq0.resize(Ncv, Ncv);// Potential without Thomas-Fermi


		#pragma omp  for schedule(dynamic) 
		for (int iq_loc = 0; iq_loc < lenght_q; iq_loc++){// create 2d simpson array
		    iq_global = begin_count + iq_loc;
		    qx_global = iq_global / (N_sum); // integer division
		    qy_global = iq_global % (N_sum); // leftover; doesn't matter which is x and y

			qx = q[qx_global];
			qy = q[qy_global];

			if (Coulomb_set.Sample_orientation == "011"){ // sample yz
				k_local.setcrys(0.0, qx, qy); //to create k-points in crystal coordinates for phase and Coulomb
			} else if (Coulomb_set.Sample_orientation == "110"){ // sample in xy
				k_local.setcrys( qx, qy , 0.0); //to create k-points in crystal coordinates for phase and Coulomb
			}
			q_mod = k_local.dot(k_local);
			if(q_mod < (q_peak * q_peak) and q_mod > pow(10, -17)){
				Wq0.fill(0.0);
				
				// Potential without Thomas-Fermi:
				Create_Wq_with_G_sum(Wq0, Ncv, qx, qy, 
				    Displacement_orb, G_num, 
				    q_TF0, Coulomb_set, G_vec, k_local, Re_Vq_max_eV, Im_Vq_max_eV, Re_Vq_min_eV, Im_Vq_min_eV ); //q_TF0 =0


				for (int ic=0; ic<Ncv; ic++){// summation over bands
	        		for (int jc=ic; jc<Ncv; jc++){ // summation over bands
	        			index_linear = ic*Ncv + jc;

	        			Vq_screened = Fourier_series_Coloumb(qx, qy ,Ncut, 
										Coulomb_set.A, Coulomb_set.B, Coulomb_set.C, Coulomb_set.D, index_linear);
	        			dV_V0_Vq_screened = (Wq0[ic][jc] - Vq_screened) * dk2; 
	        			// if(real(dV_V0_Vq_screened) > 0)
	        			{
	            			sum_loc[0][index_linear] += dV_V0_Vq_screened;
	            			sum_loc[1][index_linear] += dV_V0_Vq_screened* qx * (qx/2.0); // need it to check higher order parameters
	            			sum_loc[2][index_linear] += dV_V0_Vq_screened* qy * (qy/2.0);
	            			sum_loc[3][index_linear] += dV_V0_Vq_screened* qx * qy;
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



	MPI_Barrier(MPI_COMM_WORLD);


	MPI_reduce_vec2x_safe(Coulomb_set.Screen_const, root_rank, rank_, 4, n_offdia_terms);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast_vec2x_safe(Coulomb_set.Screen_const, root_rank, rank_, 4, n_offdia_terms);
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank_ == root_rank){
		string Screen_const_txt; 
		ofstream Screen_const_stream;
		
		Screen_const_txt= "Output/Coulomb/Screen_const_" + label_Coulomb + ".txt";;
		Screen_const_stream.open(Screen_const_txt.c_str());
		for (int i_band = 0; i_band < n_offdia_terms; i_band++){// summation over bands
			Screen_const_stream << setprecision(20) << real(Coulomb_set.Screen_const[0][i_band]) << " " << imag(Coulomb_set.Screen_const[0][i_band]) << endl;
		}
		Screen_const_stream.close();

		cout << "Correction factor. First index is order, second is band ic*Ncv + jc" << endl;
		cout << "Coulomb_set.Screen_const[0][0] = " << Coulomb_set.Screen_const[0][0]*energy_au_eV << " eV" << endl;
		cout << "Coulomb_set.Screen_const[1][0] = " << Coulomb_set.Screen_const[1][0]*energy_au_eV << " eV"<< endl;
		cout << "Coulomb_set.Screen_const[2][0] = " << Coulomb_set.Screen_const[2][0]*energy_au_eV << " eV"<< endl;
		cout << "Coulomb_set.Screen_const[3][0] = " << Coulomb_set.Screen_const[3][0]*energy_au_eV << " eV"<< endl << endl;

		cout << "Coulomb_set.Screen_const[0][1] = " << Coulomb_set.Screen_const[0][1]*energy_au_eV << " eV"<< endl;
		cout << "Coulomb_set.Screen_const[1][1] = " << Coulomb_set.Screen_const[1][1]*energy_au_eV << " eV"<< endl;
		cout << "Coulomb_set.Screen_const[2][1] = " << Coulomb_set.Screen_const[2][1]*energy_au_eV << " eV"<< endl;
		cout << "Coulomb_set.Screen_const[3][1] = " << Coulomb_set.Screen_const[3][1]*energy_au_eV << " eV"<< endl << endl;
	}
} // end if (Coulomb_set.Calculate){




if (Coulomb_set.Read_from_files){

    cout << "reading Coulomb parameters..." << endl;
    ifstream V_Hartree_stream;

    string V_Hartree_txt = "Output/Coulomb/V_Hartree_" + Coulomb_set.labelInput + ".txt";
    V_Hartree_stream.open(V_Hartree_txt);
    string trash; 
    size_t pos;
    vector<string> str;
    Coulomb_set.V_Hartree.fill(0.0);



    for (int ic=0; ic<Ncv; ic++){// summation over bands
        getline(V_Hartree_stream, trash);   
        pos = trash.length();
        Separate_string(trash, str, pos);

        // if(rank_ == root_rank){
	    //     cout << "trash " << trash << endl;
	    //     cout << "pos " << pos << endl;
	    //     // cout << "str[jc] " << str[jc] << endl << endl;
	    //     cout << " atof(str[0].c_str()) " << atof(str[0].c_str()) << endl;
	    //     cout << " atof(str[1].c_str()) " << atof(str[1].c_str()) << endl;
	    //     // cout << " atof(str[1].c_str()) " << atof(str[1].c_str()) << endl;
    	// }
    	MPI_Barrier(MPI_COMM_WORLD);
        for (int jc=0; jc<Ncv; jc++){
	        if(rank_ == root_rank){
		        // cout << "trash " << trash << endl;
		        // cout << "pos " << pos << endl;
		        // // cout << "str[jc] " << str[jc] << endl << endl;
		        // cout << " atof(str[jc].c_str()) " << atof(str[jc].c_str()) << endl;
	    	}
            Coulomb_set.V_Hartree[ic][jc] = atof(str[jc].c_str()) + 0.0*c1;
			if (rank_ == root_rank){
				cout << "ic= " << ic << "; jc= " << jc << " real(V_Hartree[ic][jc])= " << real(Coulomb_set.V_Hartree[ic][jc])*energy_au_eV<< " eV; imag(V_Hartree[ic][jc])= " << imag(Coulomb_set.V_Hartree[ic][jc])* energy_au_eV << " eV" << endl;
			}
        }
    }
    V_Hartree_stream.close();

    string A_txt = "Output/Coulomb/A_coeff" + Coulomb_set.labelInput + ".txt";
    string B_txt = "Output/Coulomb/B_coeff" + Coulomb_set.labelInput + ".txt";
    string C_txt = "Output/Coulomb/C_coeff" + Coulomb_set.labelInput + ".txt";
    string D_txt = "Output/Coulomb/D_coeff" + Coulomb_set.labelInput + ".txt";
    ifstream A_stream, B_stream, C_stream, D_stream;
    A_stream.open(A_txt);
    B_stream.open(B_txt);
    C_stream.open(C_txt);
    D_stream.open(D_txt);
    string trash_A, trash_B, trash_C, trash_D;
    size_t pos_A, pos_B, pos_C, pos_D;
    vector<string> str_A, str_B, str_C, str_D;

    for (int mm =0; mm< Ncut*Ncut; mm++ ){ //row summation over Fourier series

        getline(A_stream, trash_A); 
        getline(B_stream, trash_B); 
        getline(C_stream, trash_C); 
        getline(D_stream, trash_D); 
        pos_A = trash_A.length();
        pos_B = trash_B.length();
        pos_C = trash_C.length();
        pos_D = trash_D.length();

        Separate_string(trash_A, str_A, pos_A); 
        Separate_string(trash_B, str_B, pos_B); 
        Separate_string(trash_C, str_C, pos_C); 
        Separate_string(trash_D, str_D, pos_D); 
        // if (rank_ == root_rank) cout << " trash_A " << trash_A << endl;
        // if (rank_ == root_rank) cout << " atof(str_A[0].c_str()) " << atof(str_A[0].c_str()) << endl;
        // if (rank_ == root_rank) cout << " atof(str_A[1].c_str()) " << atof(str_A[1].c_str()) << endl;
        MPI_Barrier(MPI_COMM_WORLD);

        for (int i_band = 0; i_band < n_offdia_terms; i_band ++){
        	// if (rank_ == root_rank) cout << "int(i_band/2)= " << int(i_band/2) << endl;
            Coulomb_set.A[mm][i_band] = atof(str_A[2*i_band].c_str()) + c1*atof(str_A[2*i_band + 1].c_str());
            Coulomb_set.B[mm][i_band] = atof(str_B[2*i_band].c_str()) + c1*atof(str_B[2*i_band + 1].c_str());
            Coulomb_set.C[mm][i_band] = atof(str_C[2*i_band].c_str()) + c1*atof(str_C[2*i_band + 1].c_str());
            Coulomb_set.D[mm][i_band] = atof(str_D[2*i_band].c_str()) + c1*atof(str_D[2*i_band + 1].c_str());
           
			if ((rank_ == root_rank) and (mm ==  1)){
				cout << i_band << " ReA[mm=1] = " << real(Coulomb_set.A[mm][i_band]) << " ImA[mm=1] = " << imag(Coulomb_set.A[mm][i_band]) << endl;
			}


        }

    }
    if (rank_ == root_rank) cout << "Fourier coefficients A,B,C,D read correctly " << endl;
    A_stream.close();
    B_stream.close();
    C_stream.close();
    D_stream.close();


    string Screen_const_txt = "Output/Coulomb/Screen_const_" + Coulomb_set.labelInput + ".txt";
    ifstream Screen_const_stream;
    Screen_const_stream.open(Screen_const_txt);
    Coulomb_set.Screen_const.fill(0.);

    for (int i_band = 0; i_band < n_offdia_terms; i_band ++){// summation over bands
        getline(Screen_const_stream, trash);   
        pos = trash.length();
        Separate_string(trash, str, pos);
        Coulomb_set.Screen_const[0][i_band] = atof(str[0].c_str()) + c1*atof(str[1].c_str());

    }
    if (rank_ == root_rank){
    	cout << "Coulomb_set.Screen_cons read correctly " << endl;
    	for (int i_band = 0; i_band < n_offdia_terms; i_band ++){
    		cout << "i_band= " << i_band << "; Coulomb_set.Screen_const = " << Coulomb_set.Screen_const[0][i_band]*energy_au_eV << " eV" << endl;
    	}
    } 
    Screen_const_stream.close();
}








if (rank_ == root_rank) {
	std::cout  << "   Coulomb_set.qTF   " << Coulomb_set.qTF << endl;
	std::cout  << "   Coulomb_set.epsilon_static * Area   " << Coulomb_set.epsilon_static << endl;
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