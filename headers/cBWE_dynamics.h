// Functions that get dynamicd of the system 
template <typename T>
void DD(T& Model, vec1d& EF1, Coord_B& kk, vec2x& Dipole,int shifti, int shiftj)
{
    for (int i=0; i<Dipole.n1(); i++) {
        for (int j=0; j<Dipole.n2(); j++) {
            Dipole[i][j]=(Model.dipolex(i+shifti,j+shiftj,kk)*EF1[0]+Model.dipoley(i+shifti,j+shiftj,kk)*EF1[1]+Model.dipolez(i+shifti,j+shiftj,kk)*EF1[2]);
        }
    }
}

// template <typename T>
// void DDD_old(T& Model, vec1d& ue1,vec1d& ue2, Coord_B& kk, double& EF, double& EFX, vec2x& Dk,int Nch)
// {
//     for (int i=0; i<Nch; i++) {
//         for (int j=0; j<Nch; j++) {
//            Dk[i][j]=(Model.dipolex(i,j,kk)*EF1[0]+Model.dipoley(i,j,kk)*EF1[1]+Model.dipolez(i,j,kk)*EF1[2]);
// 	       //Dk[i][j]=Model.dipoley(i,j,kk)*EF;
//         }
//         for (int j=Nch; j<Dk.n2(); j++) {
//             Dk[i][j]=(Model.dipolex(i,j,kk)*EFX[0]+Model.dipoley(i,j,kk)*EFX[1]+Model.dipolez(i,j,kk)*EFX[2])*0.5; //we already include the RWA
//         }
//     }
//     for (int i=Nch; i<Dk.n1(); i++) {
//         for (int j=0; j<Nch; j++) {
//             Dk[i][j]=(Model.dipolex(i,j,kk)*EFX[0]+Model.dipoley(i,j,kk)*EFX[1]+Model.dipolez(i,j,kk)*EFX[2])*0.5; //we already include the RWA
//         }
//         for (int j=Nch; j<Dk.n2(); j++) {
//            Dk[i][j]=(Model.dipolex(i,j,kk)*EF1[0]+Model.dipoley(i,j,kk)*EF1[1]+Model.dipolez(i,j,kk)*EF1[2]);
//         	//Dk[i][j]=Model.dipoley(i,j,kk)*EF;
//         }
//
//     }
// }

template <typename T>
void DDD(T& Model, vec1d& EF1,vec1d& EFX, Coord_B& kk, vec2x& Dk,int Nch)
{
    for (int i=0; i<Nch; i++) {
        for (int j=0; j<Nch; j++) {
           Dk[i][j]=.5*( (Model.dipolex(i,j,kk)+conj(Model.dipolex(j,i,kk)))*EF1[0]+ (Model.dipoley(i,j,kk)+conj(Model.dipoley(j,i,kk)))*EF1[1]+ (Model.dipolez(i,j,kk)+conj(Model.dipolez(j,i,kk)) )*EF1[2]);
           //Dk[i][j]=Model.dipoley(i,j,kk)*EF;
        }
        for (int j=Nch; j<Dk.n2(); j++) {
            Dk[i][j]=.5*((Model.dipolex(i,j,kk)+conj(Model.dipolex(j,i,kk)))*EFX[0]+(Model.dipoley(i,j,kk)+conj(Model.dipoley(j,i,kk)))*EFX[1]+(Model.dipolez(i,j,kk)+conj(Model.dipolez(j,i,kk)))*EFX[2])*0.5; //we already include the RWA
        }
    }
    for (int i=Nch; i<Dk.n1(); i++) {
        for (int j=0; j<Nch; j++) {
            Dk[i][j]=.5*((Model.dipolex(i,j,kk)+conj(Model.dipolex(j,i,kk)))*EFX[0]+(Model.dipoley(i,j,kk)+conj(Model.dipoley(j,i,kk)))*EFX[1]+(Model.dipolez(i,j,kk)+conj(Model.dipolez(j,i,kk)))*EFX[2])*0.5; //we already include the RWA
        }
        for (int j=Nch; j<Dk.n2(); j++) {
           Dk[i][j]=.5*( (Model.dipolex(i,j,kk)+conj(Model.dipolex(j,i,kk)))*EF1[0]+ (Model.dipoley(i,j,kk)+conj(Model.dipoley(j,i,kk)))*EF1[1]+ (Model.dipolez(i,j,kk)+conj(Model.dipolez(j,i,kk)) )*EF1[2]);
            //Dk[i][j]=Model.dipoley(i,j,kk)*EF;
        }

    }
}






template <typename TT>
void get_cBWE_dynamic(TT& Model, vec2d& k, Coord_B& kk,  vec3x& P, 
	vec3x& Pv,  int indexDeriv,
	vec2d& T, vec1i& Nb,vec1d& EF1, vec1d& EF2,double wx2,
	vec1d& AF1, Coulomb_parameters& Coulomb_set, 
	vector<bool>  &id_masters,vec1d& integrWeight, vec3x& P_diag_local){
/*
Fucntion calculate the derivative of population
dP/dt from cBWE. The result is stored in variable Pv
*/
	int Ncv=Nb[0]+Nb[1]+Nb[2];
	int Nch=Nb[0];



	if (Coulomb_set.Coulomb_calc)
		{// calculate set of coefficients to obtain the Coulomb terms
		if(Coulomb_set.Wannie_basis){
			Calculate_X_coefficients(P, k, AF1, Ncv, Coulomb_set, 
			id_masters, integrWeight);
		} else if(Coulomb_set.Diagonal_basis){	
			Calculate_X_coefficients(P_diag_local, k, AF1, Ncv, Coulomb_set, 
			id_masters, integrWeight);
		}

		
	}
    complexd alpha;
    vec2x Hk(Ncv,Ncv);
    vec2x Uk(Ncv,Ncv);
    vec2x Dk(Ncv,Ncv);
    vec2x Xk_W(Ncv,Ncv);
    vec2x Xk_dia(Ncv,Ncv);
    P_diag_local.fill(0.0);

for (int ik=0; ik < P.n1(); ik++)    { // sum over all wave vectors, available at current node

    //int tid = omp_get_thread_num();

    kk.setcrys(k[ik][0]+AF1[0], k[ik][1]+AF1[1], k[ik][2]+AF1[2]); // calculate cartesian coord for this crystal
    //kk.setcrys(k[ik][0]+AF1[0],k[ik][1]+AF1[1],k[ik][2]+AF1[2]);
    DDD(Model, EF1, EF2, kk, Dk, Nb[0]);
    Model.energy_U(Hk, Uk, kk);



    // taking into account the Coulomb term - it renorms the energy matrix!
	if (Coulomb_set.Coulomb_calc)
	{ // calculating X(k, lambda, lambda') for every particular vave vector
		if(Coulomb_set.Wannie_basis){
			Calculate_X(kk, Xk_W, Ncv, Coulomb_set);
		}
		if(Coulomb_set.Diagonal_basis){
		    // get population in diagonal basis
		    // yes, we use population of previous step
		    // to calculate Xk_dia
		    // because I don't know how to implement it not affecting other code
			for (int ic=0; ic<Ncv; ic++){
				for (int jc=0; jc<Ncv; jc++){
					for (int ii=0; ii<Ncv; ii++){
						for (int jj=0; jj<Ncv; jj++){
							P_diag_local[ik][ic][jc] += \
								conj(Uk[ic][ii])*P[ik][ii][jj]*Uk[jc][jj];
						}//end jj
					}//end ii
				}
			}//End ic
			Calculate_X(kk, Xk_dia, Ncv, Coulomb_set);
		}
		
		if(Coulomb_set.Diagonal_basis){
			// from diagonal back to wannier basis:
			Xk_W.fill(0.0);
			for (int ic=0; ic<Ncv; ic++){
				for (int ii=0; ii<Ncv; ii++){
					for (int jj=0; jj<Ncv; jj++){
						for (int jc=0; jc<Ncv; jc++){
							 // !!!!CHECK TRANSFORM!!!!
							 // !!!!CHECK TRANSFORM!!!! 
							 // !!!!CHECK TRANSFORM!!!!
							Xk_W[ic][jc] += conj(Uk[ii][ic]) *Xk_dia[ii][jj] *Uk[jj][jc];
						}
					}//end jj
				}//end ii
			}//End ic
		} 
		for (int ic=0; ic<Ncv; ic++){// summation over bands
			for (int jc=0; jc<Ncv; jc++){ // summation over bands
				Hk[ic][jc] -= Xk_W[jc][ic]; // see theory why this order of jc ic
			}
		}
	} // end of Coulomb_calc
/*
*****************************************************************************
******************       Loop for core holes and valence        *************
******************        Population(rho_ii)        *************************
****************************************************************************/
	for (int ic=0; ic<Ncv; ic++) { //we shouldn't account for j==ic, but in principle the imag() function will deliver 0 for this case
		alpha=0.;
		for (int j=0; j<Ncv; j++) {
			alpha += (Hk[ic][j] + Dk[ic][j]) *P[ik][ic][j];
		}
		Pv[indexDeriv][ik][Ncv*ic + ic] = 2.*imag(alpha);
	} //end ic


//2nd loop ic!=jc
/****************************************************************************
******************             Loop for             *************************
******************        Coherences(rho_ij)        *************************
*****************************************************************************
****************************************************************************/
	for (int ic=0; ic<Ncv; ic++) {
		for (int jc= ic+1; jc<Ncv; jc++){
			alpha=0.;
			for (int j=0; j<Ncv; j++) {
				alpha += (Hk[jc][j] + Dk[jc][j]) *P[ik][ic][j]; 
			}
			for (int j=0; j<Ncv; j++){ 
				alpha -= conj((Hk[ic][j] + Dk[ic][j]) *P[ik][jc][j]);
			}
			Pv[indexDeriv][ik][Ncv*ic + jc] = -c1*alpha;
			if( (ic < Nch) && (jc >= Nch)) { //term due to RWA
				Pv[indexDeriv][ik][Ncv*ic + jc] += -c1 *(-wx2 *P[ik][ic][jc]);
			} 
		} //end jc
	} //end ic


/****************************************************************************
******************         Dissipation term         *************************
******************       in the Wannier gauge       *************************
****************************************************************************/
	for (int ic = 0; ic < Ncv; ic++){
		for (int jc = ic; jc < Ncv; jc++){
			complexd Tii = 0.; //calculation of T[ic][jc][ii][ij]
			for (int ii = 0; ii < Ncv; ii++){
				for (int jj = 0; jj < Ncv; jj++){
					complexd Pij = 0.;
					for (int i = 0; i < Ncv; i++){
						for (int j = 0; j < Ncv; j++){
							Pij += conj(Uk[ii][i])*P[ik][i][j]*Uk[jj][j]; //this term depends on ii,jj,i,j
						}//end sum over j (in the notes n')
					}//end sum over i (in the notes n)
					Pv[indexDeriv][ik][Ncv*ic + jc] -= conj(Uk[jj][jc]) *Uk[ii][ic] *T[ii][jj] *Pij; //part of dissipation coming from the electron equations
				}//end sum over jj (in the notes j)
			}//end sum over ii (in the notes i)
			for (int i=0; i<(Nb[0]+Nb[1]); i++) {
				Tii += conj(Uk[i][jc]) *Uk[i][ic] *T[i][i];
			}
			//for (int i=0; i<(Nch+Nc); i++) Tii+=conj(Model.unitary(i,jc,kk))*Model.unitary(i,ic,kk)*T[i][i];
			Pv[indexDeriv][ik][Ncv*ic + jc] += Tii; //part of dissipation coming only from the holes
			}//end jc
		}//end ic

	} //end ikz



	for (int ik=0; ik<P.n1(); ik++) {
		for (int ic=0; ic<Ncv; ic++) {
		  //#pragma omp simd
		    for (int jc=ic+1; jc<Ncv; jc++){
		      Pv[indexDeriv][ik][Ncv*jc + ic] = \
		      conj(Pv[indexDeriv][ik][Ncv*ic + jc]);
		    }
		}
	}



}






//d^{n}P(t)= (d^{n-1}P(t) .- d^{n-1}P(t-dt))./dt
void getDerivative (vec3x& Pn, vec3x& Pn_minus_1, 
	int Pn_index, vector<int>&  Pn_minus_1_index, double dt){
	double scaleCoef = pow(10, 10);
	for(int ik = 0; ik < Pn.n2(); ik++){ // wave vectors
		#pragma omp simd
		for (int ic = 0; ic < Pn.n3(); ic++){ // bands
			Pn[Pn_index][ik][ic] = \
			((scaleCoef* Pn_minus_1[Pn_minus_1_index[0]][ik][ic] - \
			 scaleCoef* Pn_minus_1[Pn_minus_1_index[1]][ik][ic])/dt) / scaleCoef;
		}

	}


}