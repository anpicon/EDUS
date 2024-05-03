// Linear Multistep methods for ODE


// get electric field
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
Multistep_index_shift(First_derivatIndex, Diff_Eq.SolverOrder);
#pragma omp barrier


get_derivative_Df(P0, OMP_private.Pv, T, Nb, 
    EF_pr, pulse2.wl, 
    Coulomb_set, trig_k_omp, OMP_private,
    GradientIndex, Weigths, Bvector, 
    root_rank, rank_, E_non0);



/*
Here we fill the matrix of first derivatives with new values,
thats why we pass the derivativesMultisteps wich contain
all first derivatives and current index First_derivatIndex[0
The matrix
(derivativesMultisteps[First_derivatIndex[0]])[lenght_k][Ncv][Ncv])
we rewrite with new values
*/

copy_Pv_Multistep(OMP_private.Pv, (derivativesMultisteps[First_derivatIndex[0]]), OMP_private.lenght_k,  Ncv);

complex<double> Delta_Pk;
int ik; // position in shared variable
for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++){
    ik = ik_pr + OMP_private.begin_count;
    for (int ic = 0; ic < Ncv; ic++){ // bands
        for (int jc = ic; jc < Ncv; jc++){ //
        	Delta_Pk = 0.0;
        	for(int step = 0; step < Diff_Eq.SolverOrder; step++){ // create sum over Coef_Multistep[step] * f[step]
        		Delta_Pk += Coef_Multistep[step] * (derivativesMultisteps[First_derivatIndex[step]])[ik_pr][ic][jc];
        	}
        	P0[ik][ic][jc] += dt *Delta_Pk;
        	if(ic != jc) P0[ik][jc][ic] = conj(P0[ik][ic][jc]);
        }
    }
}