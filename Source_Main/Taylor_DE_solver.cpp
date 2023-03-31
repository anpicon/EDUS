// Taylor time evolution
#pragma omp barrier
#pragma omp master
{ // shared variables, only master change it
    Mpi_communication(P0,  message);
    MPI_Barrier(MPI_COMM_WORLD);
}

Taylor_index_shift(derivatIndex, TaylorOrder);
#pragma omp barrier

EF_pr[0] = pulse1.E(time_loop);
EF_pr[1] = .5*pulse2.E(time_loop);



get_derivative_Df(kpt, P0, OMP_private.Pv, T, Nb, 
    EF_pr, pulse2.wl, 
    Coulomb_set, trig_k_omp, OMP_private,
    GradientIndex, Weigths, Bvector, 
    root_rank, rank_, message);


/*
Here we fill the matrix of first derivatives with new values,
thats why we pass the (derivativesMatrices[1]) wich contain
all first derivatives and current index (derivatIndex[1])[0]
which means current time, which first index of
(derivativesMatrices[1])[(derivatIndex[1])[0]][lenght_k][Ncv*Ncv])
we rewrite with new values
*/
// copy OMP_private.Pv into (derivativesMatrices[1])
copy_Pv_derivMatrices(OMP_private.Pv, 
    (derivativesMatrices[1]), (derivatIndex[1])[0], OMP_private, Ncv);


// Now calculate all higher derivatives and shift their indeces
for (int derivOrder = 2; derivOrder < (TaylorOrder + 1); derivOrder++){
    vec3x& currDerivMatr = (derivativesMatrices[derivOrder]);
    vector<int>& currDerivIndex = (derivatIndex[derivOrder]);
    int lenghtDeriv = TaylorOrder - derivOrder + 1;
    int indexFill = currDerivIndex[lenghtDeriv-1];
    for (int numPoint = (lenghtDeriv-1); numPoint >  0; numPoint--){
        currDerivIndex[numPoint] = currDerivIndex[numPoint-1];
    }
    currDerivIndex[0] = indexFill; // see comments for first derivative

    //d^{n}F(t)= (d^{n-1}F(t) .- d^{n-1}F(t-dt))./dt
    getDerivative(currDerivMatr,
        (derivativesMatrices[derivOrder-1]),
        currDerivIndex, (derivatIndex[derivOrder - 1]), dt_prev);
}



// now we do the Taylor step:
double ht = 1;
double factorial = 1;
Delta_P.fill(0.0);
for (int derivOrder = 1; derivOrder < (TaylorOrder + 1); derivOrder++){
    ht *= dt;
    factorial *= derivOrder;
    int index0 = (derivatIndex[derivOrder])[0];
    vec3x& currDerivMatr = (derivativesMatrices[derivOrder]);


    for(int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++){ // wave vectors
        for (int ic = 0; ic < Ncv*Ncv; ic++){ // bands
            Delta_P[ik_pr][ic] += \
            currDerivMatr[index0][ik_pr][ic] * ht / factorial;
        }
    }
}



int ik; // position in shared variable
for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++){
    ik = ik_pr + OMP_private.begin_count;
    for (int ic = 0; ic < Ncv; ic++){ // bands
        for (int jc = 0; jc < Ncv; jc++){ //
            P0[ik][ic][jc] += Delta_P[ik_pr][Ncv*ic + jc];
        }
    }
}