// taylor solver

void init_Taylor_matrices(vector<vec3x> &derivativesMatrices, 
	vector<vector<int>> & derivatIndex, 
	int& TaylorOrder, int& lenght_k, int& Ncv ){
	/*
	vector<vec3x> derivativesMatrices(TaylorOrder + 1);
	here we need (TaylorOrder - order + 1) derivatives for each order.
	Wich means that for example that for TaylorOrder 8 we need to store
	8 derivatives of order 1 and only 1 derivatives of order 2.
	In the beginning all the derivatives are 0, because we start from equilibrium
	*/

	//vector<vector<int>> derivatIndex(TaylorOrder + 1);
	// we need to store many derivations, but we need also to store
	// the right position of this derivation in time scale
	// (derivatIndex[order])[0] means that we have the derivation at the current
	// moment t. (derivatIndex[order])[i] means the derivation at moment
	// t - i*dt


	// allocate all matrix of derivatives:
	for (int derivOrder = 1; derivOrder < (TaylorOrder + 1); derivOrder++){
		int lenghtDeriv = TaylorOrder - derivOrder + 1;
		(derivatIndex[derivOrder]).resize(lenghtDeriv);
		(derivativesMatrices[derivOrder]).resize(lenghtDeriv, lenght_k,  Ncv*Ncv);
		(derivativesMatrices[derivOrder]).fill(0.0);
		for(int numPoint = 0; numPoint< lenghtDeriv; numPoint++){
			(derivatIndex[derivOrder])[numPoint] = numPoint;
		}
	}

}


void Taylor_index_shift(vector<vector<int>> & derivatIndex, 
	int& TaylorOrder){

	int indexFill = (derivatIndex[1])[(TaylorOrder-1)];
	for (int numPoint = (TaylorOrder-1); numPoint >  0; numPoint--){
	(derivatIndex[1])[numPoint] =(derivatIndex[1])[numPoint-1];
	}
	(derivatIndex[1])[0] = indexFill; // shifting index
	// matrix with index (derivatIndex[1])[0] is always to be rewrited in func


}

void copy_Pv_derivMatrices(vec3x& Pv, vec3x& Pv_Taylor,
	int indexDeriv, Private_omp_parameters& OMP_private, int Ncv){
	
	int round_digit = 5;
	double scale = pow(10,round_digit);
	std::complex<double> Pv_k;
    for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++) {
		for (int ic=0; ic<Ncv; ic++){
			// #pragma omp simd
			for (int jc=0; jc<Ncv; jc++){
				Pv_k = Pv[ik_pr][ic][jc];
				Pv_Taylor[indexDeriv][ik_pr][Ncv*ic + jc] = Pv_k; // floor_to_digits(real(Pv_k), round_digit) + floor_to_digits(imag(Pv_k), round_digit) * 1i ;
			}
		}

    }
}

//d^{n}P(t)/dt^{n}= 2* (d^{n-1}P(t)/dt^{n-1} .- d^{n-1}P(t-dt)//dt^{n-1})./dt 
// -  d^{n}P(t-dt)/dt^{n} 
void getDerivative (vec3x& Pn, vec3x& Pn_minus_1, 
	vector<int>& Pn_index_arr,  vector<int>&  Pn_minus_1_index, vec1d & dt){
	int Pn_index0 = Pn_index_arr[0];
	int index_size = Pn_minus_1_index.size();
	double dt_n = 0;
	for (int i = 0; i < index_size; i++){
		dt_n += dt[i];
	}

	vec1x Pn_prev(Pn.n3());
	for(int ik = 0; ik < Pn.n2(); ik++){ // wave vectors
		// #pragma omp simd
		for (int ic = 0; ic < Pn.n3(); ic++){ // bands
			// Pn_prev[ic] = Pn[Pn_index1][ik][ic];
			Pn[Pn_index0][ik][ic] =  ((Pn_minus_1[Pn_minus_1_index[0]][ik][ic] - \
			 Pn_minus_1[Pn_minus_1_index[index_size-1]][ik][ic])/dt_n); //- Pn_prev[ic];
		}

	}


}


