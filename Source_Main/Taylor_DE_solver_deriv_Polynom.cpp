// Taylor time evolution
// see Methods of DE propagations in Coulomb notes 

for (int iE_2 = 0; iE_2 < 3; iE_2++) EF_pr[0][iE_2] = 0.;

for(int i_pump = 0; i_pump < Laser_pumps.size(); i_pump++){ // sum over all pump pulses
    vector<vec1d> EF_pr_copy(2); 
    EF_pr_copy[0] = Laser_pumps[i_pump].E(time_loop);
    
    for (int iE_2 = 0; iE_2 < 3; iE_2++) EF_pr[0][iE_2] += EF_pr_copy[0][iE_2];
}
EF_pr[1] = .5*pulse2.E(time_loop);

// electric field modulus:
Emax = 0; // all components of field. We need it to know, if we need this field
for (int iE_1 = 0; iE_1 < 2; iE_1++){
    for (int iE_2 = 0; iE_2 < 3; iE_2++){
        Emax += EF_pr[iE_1][iE_2]*EF_pr[iE_1][iE_2];
    }
}
Emax = sqrt(Emax); 

E_non0 = false;
if (Emax > E_eps){
    E_non0 = true;
}

#pragma omp barrier
#pragma omp master
{ // shared variables, only master change it

    Coulomb_set.Taylor_Delta_Pk_max.fill(0.);
    Coulomb_set.Taylor_alpha_max.fill(0.);
    if (E_non0){
        Mpi_communication(P0,  message);
        MPI_Barrier(MPI_COMM_WORLD);
    }
}


/*
Index shift means:

previous time
moments          t0   t1  t2  t3  ...     tn

derivatIndex
step 0            0   1   2   3  ...      n    

derivatIndex
step 1            n   0   1   2  ...      n-1  

derivatIndex
step 2            n-1 n   0   1   2  ...  n-2 

and etc.

So for example if we want something in the current moment 0, 
we have to addres the element with index derivatIndex[0]
*/
Taylor_index_shift(derivatIndex, TaylorOrder);
#pragma omp barrier





get_derivative_Df(kpt, P0, OMP_private.Pv, T, Nb, 
    EF_pr, pulse2.wl, 
    Coulomb_set, trig_k_omp, OMP_private,
    GradientIndex, Weigths, Bvector, 
    root_rank, rank_, E_non0);

// store population before we did time step in OMP_private.P_W_prev
if (Diff_Eq.dynamical_dt_evolution){
    for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++){ // all local wave vectors
        for(int ic=0; ic<Ncv; ic++){
            OMP_private.P_W_prev[ik_pr][ic][ic] = P0[ik_pr + OMP_private.begin_count][ic][ic];
            // OMP_private.P_dia_prev[ik_pr][ic][jc] = OMP_private.P_diag[ik_pr][ic][jc];
            
        }
    }
}


/*
Here we fill the matrix of first derivatives with new values,
thats why we pass the (derivativesMatrices[1]) wich contain
all first derivatives and current index (derivatIndex[1])[0]
which means current time, 
The matrix
(derivativesMatrices[1])[(derivatIndex[1])[0]][lenght_k][Ncv*Ncv])
we rewrite with new values
*/
// copy OMP_private.Pv into (derivativesMatrices[1])
copy_Pv_derivMatrices(OMP_private.Pv, 
    (derivativesMatrices[1]), (derivatIndex[1])[0], OMP_private, Ncv);


vec1d t_i(TaylorOrder-1); // here we store all t_i
t_i[0] = dt_prev[0];
for (int i = 1; i< TaylorOrder-1; i++){
    t_i[i] = t_i[i-1] + dt_prev[i-1];
}

int scale_digit = 7;
int round_digit = 5;
double scale = pow(10,scale_digit); 
mat tau_ip(TaylorOrder-1, TaylorOrder-1);

for (int i = 0; i< (TaylorOrder-1); i++){ //
    for (int p = 0; p< (TaylorOrder-1); p++){ //
        tau_ip(i, p) = scale * pow(t_i[i] , (p+1)); // we count matrix from 0, that's why it's p+1
    }
}

mat tau_inv = inv(tau_ip); // creating tau^{-1} to calculate coefficients

for (int i = 0; i< (TaylorOrder-1); i++){ //
    for (int p = 0; p< (TaylorOrder-1); p++){ //
        tau_inv(i, p) *= scale; // we count matrix from 0, that's why it's p+1
    }
}



cx_vec dP_i(TaylorOrder-1);
cx_vec alpha_i(TaylorOrder-1);
std::complex<double> Delta_Pk;
std::complex<double> Delta_Pk_i;
double i_double;
vec1d Delta_Pk_max(TaylorOrder);
vec1d alpha_max(TaylorOrder);
Delta_Pk_max.fill(0.0);
alpha_max.fill(0.0);
double Delta_Pk_max_loc, alpha_max_loc, Pv_max;

vector<int>& P_1_k_index = (derivatIndex[1]);
vec3x& P_1_k = (derivativesMatrices[1]);
// calculate alpha coefficients for each k- point and calculate next step
int ik; // position in shared variable
for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++){
    ik = ik_pr + OMP_private.begin_count;
    for (int ic = 0; ic < Ncv; ic++){ // bands
        for (int jc = 0; jc < Ncv; jc++){ //

            for (int i = 0; i< TaylorOrder-1; i++){
                dP_i(i) = P_1_k[P_1_k_index[0] ][ik_pr][Ncv*ic + jc] - P_1_k[P_1_k_index[i+1]][ik_pr][Ncv*ic + jc];
            }
            alpha_i = tau_inv * dP_i;

            Delta_Pk = P_1_k[P_1_k_index[0] ][ik_pr][Ncv*ic + jc] * dt;
            Delta_Pk_max_loc = abs(Delta_Pk);
            
            if (Delta_Pk_max_loc > Delta_Pk_max[0] ){
                Delta_Pk_max[0] = Delta_Pk_max_loc;
                Pv_max = abs(P_1_k[P_1_k_index[0] ][ik_pr][Ncv*ic + jc]); 
            }
            for (int i = 1; i < TaylorOrder ; i++){ // here we messed up with indices because we have arrays from 0
                i_double = i+1;
                Delta_Pk_i = pow(dt, i+1) * alpha_i(i-1) /i_double;
                Delta_Pk_max_loc = abs(Delta_Pk_i);

                if (Delta_Pk_max_loc > Delta_Pk_max[i] ){
                    Delta_Pk_max[i] = Delta_Pk_max_loc;
                    alpha_max[i] =  abs(alpha_i(i-1));
                }
                // Delta_Pk += round_to_digits(real(Delta_Pk_i), round_digit) + round_to_digits(imag(Delta_Pk_i), round_digit) * 1i  ;
                Delta_Pk += Delta_Pk_i;
            }



            P0[ik][ic][jc] += Delta_Pk;
        }
    }
}
    

    #pragma omp barrier // sinchronise threads
    if (it % (100*it_resolution) == 0){
        #pragma omp critical
        {
            if (Coulomb_set.Taylor_Delta_Pk_max[0] < Delta_Pk_max[0] ){
                Coulomb_set.Taylor_Delta_Pk_max[0] = Delta_Pk_max[0];
                Coulomb_set.Taylor_Pv_max = Pv_max;
            }

            for (int i = 1; i < TaylorOrder ; i++){ // here we messed up with indices because we have arrays from 0
                if (Coulomb_set.Taylor_Delta_Pk_max[i] < Delta_Pk_max[i] ){
                    Coulomb_set.Taylor_Delta_Pk_max[i] = Delta_Pk_max[i];
                    Coulomb_set.Taylor_alpha_max[i] = alpha_max[i];
                }

            }


            
        }

        #pragma omp barrier // sinchronise threads
        #pragma omp master
        {
            vec1d Taylor_Delta_Pk_max_MPI(TaylorOrder);
            vec1d Taylor_alpha_max_MPI(TaylorOrder);
            double Taylor_Pv_max_MPI;
            
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Reduce(&Coulomb_set.Taylor_Delta_Pk_max[0], &Taylor_Delta_Pk_max_MPI[0], 
                TaylorOrder, MPI_DOUBLE, MPI_MAX, root_rank, MPI_COMM_WORLD);

            MPI_Reduce(&Coulomb_set.Taylor_alpha_max[0],    &Taylor_alpha_max_MPI[0], 
                TaylorOrder, MPI_DOUBLE, MPI_MAX, root_rank, MPI_COMM_WORLD);

            MPI_Reduce(&Coulomb_set.Taylor_Pv_max, &Taylor_Pv_max_MPI, 
                1, MPI_DOUBLE, MPI_MAX, root_rank, MPI_COMM_WORLD);

            MPI_Barrier(MPI_COMM_WORLD);

            if ((rank_ == root_rank) ){

                cout  << "  order i= " << 1  << " Delta_Pk_max = " << Taylor_Delta_Pk_max_MPI[0] << "Pv_max = " << Taylor_Pv_max_MPI <<  endl;
                for (int i = 1; i < TaylorOrder ; i++){

                    cout  << "  order i= " << i+1  << " Delta_Pk_max = " << Taylor_Delta_Pk_max_MPI[i] << "  alpha_max[i]= "<< Coulomb_set.Taylor_alpha_max[i] << endl;
                }
                 cout  << endl;
            }
        }
    }
    #pragma omp barrier // sinchronise threads