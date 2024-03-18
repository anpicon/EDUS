#include "Coordinate.h"

// FORMER Runge_Kutta_Df FUNCTION
// MPI version
// changed the name of function
// to get more meaningful names
void get_derivative_Df(
    vec2d& kpt, vec3x& P,vec3x& Pv,vec2d& T,
    vec1i& Nb,vector<vec1d>&  EF, double wx2, 
    Coulomb_parameters& Coulomb_set, 
    trig_coefficients & trig_k, 
    Private_omp_parameters& OMP_private,
    vector<vector<vector<int>>>& GradientIndex, 
    vector<double>& Weigths, 
    vector<vector<vector<double>>>& Bvector, 
    int root_rank, int my_rank, bool E_non0)
{

    int Ncv=P.n2();
    int Nch=Nb[0];
    complex<double> corr; // correction of peak in Coulomb case
    
    OMP_private.P_grad.fill(0.0); // variable to store gradient



    vec2x Hk(Ncv, Ncv); // temporary container for Hamiltonian

    if (Coulomb_set.Coulomb_calc)
    {// calculate set of coefficients to obtain the Coulomb terms
        if(Coulomb_set.Wannie_basis){
            Calculate_X_coefficients_MPI(P, OMP_private,  
             Ncv, Coulomb_set, 
             root_rank, my_rank, trig_k);
        } else if(Coulomb_set.Diagonal_basis){ 
        
            // create P in diagonal basis:
            int ik_pr;
            for (int ik = OMP_private.begin_count; ik < OMP_private.end_count; ik++){ // all local wave vectors
                ik_pr = ik - OMP_private.begin_count;
                
                // function transforms P[ik] with OMP_private.Uk[ik_pr] and writes it to the 
                // OMP_private.P_diag[ik]
                get_P_in_dia_vect(P,  OMP_private.Uk, OMP_private.P_diag, Ncv, 
                    ik, ik_pr, OMP_private);

            } // k loop
          

 
            Calculate_X_coefficients_MPI(OMP_private.P_diag, OMP_private,  
            Ncv, Coulomb_set, 
            root_rank, my_rank, trig_k);

        }
    } //end if (Coulomb_set.Coulomb_calc)
    

        
    vec2x Xk_W(Ncv,Ncv); // Coulomb term in Wannier Basis
    vec2x Xk_d(Ncv,Ncv); // Coulomb term in eigenstate Basis

    int ik; // position in shared variable

    for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++) // loop over position in private variable
    {
        ik = ik_pr + OMP_private.begin_count; // coordinate in global array
        
        if (Coulomb_set.Coulomb_calc){
            // Coulomb coefficients into Hamiltonian

            if(Coulomb_set.Wannie_basis){
                Calculate_X_MPI(ik_pr, Xk_W, Ncv, Coulomb_set, trig_k, OMP_private.lenght_k);

                //corrections due to the peak:
                int index_linear_tr, index_linear; 
                if (Coulomb_set.n_cond > pow(10, -20)){
                    for (int ic=0; ic<Ncv; ic++){// summation over bands   // 
                        //Xk_W[ic][ic] -= Coulomb_set.Screen_const[0][0] * (P[ik][ic][ic] - OMP_private.P_Wannier_0[ik_pr][ic][ic]);           
                        for (int jc=ic; jc<Ncv; jc++){ // summation over bands
                            index_linear_tr = jc*Ncv + ic; // transposed index for transposed population
                            if (ic != jc){
                                index_linear = ic*Ncv + jc;
                            } else {
                                index_linear = 0; // all diagonal elements of Screen_const are equal. We calculate it only for band [0][0]
                            }
                            
                            Xk_W[ic][jc] -= conj(Coulomb_set.Screen_const[0][index_linear]) * (P[ik][jc][ic] - OMP_private.P_Wannier_0[ik_pr][jc][ic]);
                            // conj(Coulomb_set.Screen_const[0][1]) here because we apply operation for reverse order of indices in population
                            
                            if (ic != jc) Xk_W[jc][ic] = conj(Xk_W[ic][jc]);
                        }
                    }
                }
            }  
            if(Coulomb_set.Diagonal_basis){
                Calculate_X_MPI(ik_pr, Xk_d, Ncv, Coulomb_set, trig_k, OMP_private.lenght_k);
                for (int ic=0; ic<Ncv; ic++){// summation over bands
                    Xk_d[ic][ic] -= Coulomb_set.Screen_const[0][0] * (OMP_private.P_diag[ik][ic][ic] - OMP_private.P_Bloch_0[ik_pr][ic][ic]);
                    // #pragma omp simd
                    for (int jc=ic+1; jc<Ncv; jc++){ // summation over bands
                        //corrections due to the peak
                        //if (Coulomb_set.n_cond > pow(10, -18))
                        {
                            Xk_d[ic][jc] -= conj(Coulomb_set.Screen_const[0][1]) * OMP_private.P_diag[ik][jc][ic];// \
                                    + Coulomb_set.Screen_const[1] * Coulomb_set.P_d2kx[ik][jc][ic] \
                                    + Coulomb_set.Screen_const[2] * Coulomb_set.P_d2ky[ik][jc][ic] \
                                    + Coulomb_set.Screen_const[3] * Coulomb_set.P_dky_dkx[ik][jc][ic];
                            Coulomb_set.Xk_storage[ik][ic][jc] = Xk_d[ic][jc];
                            Xk_d[jc][ic] = conj(Xk_d[ic][jc]);
                            Coulomb_set.Xk_storage[ik][jc][ic] = conj(Xk_d[ic][jc]);
                        }

                    }
                }



                // from diagonal back to wannier basis:
                get_Xd_in_Wannier_vect(Xk_d, Xk_W, OMP_private.Uk, 
                    Ncv,  ik_pr, OMP_private);
                
            }

            

            // new energy dispersion renormed by excitons
            for (int ic=0; ic<Ncv; ic++){// summation over bands
                Xk_W[ic][ic] = real(Xk_W[ic][ic]);
                for (int jc=ic; jc<Ncv; jc++){ // summation over bands
                    Hk[ic][jc] = OMP_private.Hk[ik_pr][ic][jc] + Xk_W[ic][jc]; 
                    if (ic == jc) Hk[ic][ic] += (Coulomb_set.E_Hartree[ic]);
                    if (ic != jc) Hk[jc][ic] = conj(Hk[ic][jc]);

                    OMP_private.Hk_renorm[ik_pr][ic][jc] = Hk[ic][jc];
                    if (ic != jc) OMP_private.Hk_renorm[ik_pr][jc][ic] = Hk[jc][ic];
                }
            }
        } else { // No Coulomb coefficients
            // Hk = Hk0;
            for (int ic=0; ic<Ncv; ic++){// summation over bands
                // #pragma omp simd safelen(4)
                for (int jc=ic; jc<Ncv; jc++){ // summation over bands
                    Hk[ic][jc] = OMP_private.Hk[ik_pr][ic][jc];
                    if (ic != jc) Hk[jc][ic] = conj(Hk[ic][jc]);
                }
            }
        }

        
        int whichEF;
        /*
         *****************************************************************************
         ******************       Loop for core holes and valence        *************
         ******************        Population(rho_ii)        *************************
         *****************************************************************************
         ****************************************************************************/
        vec1x alpha(Ncv);
        for (int ic=0; ic<Ncv; ic++) 
        {
            alpha.fill(0.);
            // #pragma omp simd
            for (int j=0; j<Ncv; j++) 
            {
                whichEF = (ic>=Nch && j<Nch) || (ic < Nch && j >= Nch);
                alpha[j] += (Hk[ic][j] + EF[whichEF][0]* OMP_private.Dk[0][ik_pr][ic][j]
                     + EF[whichEF][1] *OMP_private.Dk[1][ik_pr][ic][j] +
                    EF[whichEF][2]*OMP_private.Dk[2][ik_pr][ic][j]) *P[ik][ic][j];
            }
            complexd beta=0.; for (int j=0; j<Ncv; j++) beta+=alpha[j];
            Pv[ik_pr][ic][ic]= 2.*imag(beta);

        } //end ic

                
        
        //2nd loop ic!=jc
        /****************************************************************************
        ******************             Loop for             *************************
        ******************        Coherences(rho_ij)        *************************
        *****************************************************************************
        ****************************************************************************/
        
        for (int ic=0; ic<Ncv; ic++)
        {
            for (int jc=ic+1; jc<Ncv; jc++)
            {
                alpha.fill(0.);
                // #pragma omp simd
                for (int j=0; j<Ncv; j++)
                {
                    whichEF = (jc>=Nch && j<Nch) || (jc < Nch && j >= Nch);
                    alpha[j] += (Hk[jc][j]+EF[whichEF][0] *OMP_private.Dk[0][ik_pr][jc][j] +  \
                                EF[whichEF][1] *OMP_private.Dk[1][ik_pr][jc][j]+ \
                                EF[whichEF][2] *OMP_private.Dk[2][ik_pr][jc][j])*P[ik][ic][j];

                }
                // #pragma omp simd
                for (int j=0; j<Ncv; j++) 
                {
                    whichEF = (ic>=Nch && j<Nch) || (ic < Nch && j >= Nch);
                    alpha[j] -= conj((Hk[ic][j] + \
                                EF[whichEF][0] * OMP_private.Dk[0][ik_pr][ic][j]+ \
                                EF[whichEF][1] * OMP_private.Dk[1][ik_pr][ic][j]+ \
                                EF[whichEF][2] * OMP_private.Dk[2][ik_pr][ic][j]) *P[ik][jc][j]);
                }
                complexd beta=0.; for (int j=0; j<Ncv; j++) beta+=alpha[j];
                Pv[ik_pr][ic][jc]=-c1*beta;
                if(ic<Nch && jc>=Nch) Pv[ik_pr][ic][jc]+=-c1*(-wx2*P[ik][ic][jc]); //term due to RWA 

            } //end jc
        } //end ic

        /****************************************************************************
        ******************         Dissipation term         *************************
        ******************       in the Wannier gauge       *************************
        *****************************************************************************
        ****************************************************************************/
        if (OMP_private.n_diss_terms > 0) { // perform this only if we have non zero dissipation terms
            OMP_private.Mk.fill(0.0); // !!!if
            int ii, jj;
            for (int i_0=0; i_0 <OMP_private.n_diss_terms; i_0++){
                ii = OMP_private.T_dissip_index_0[i_0]; // we need only terms with non zero T
                jj = OMP_private.T_dissip_index_1[i_0];
                for (int i=0; i<Ncv; i++){
                    alpha.fill(0.);
                    // #pragma omp simd
                    for (int j=0; j<Ncv; j++){
                        alpha[j]= conj(OMP_private.Uk[ik_pr][ii][i])* \
                                    OMP_private.Uk[ik_pr][jj][j] *P[ik][i][j]; //this term depends on ii,jj,i,j
                    }//end j
                    for (int j=0; j<Ncv; j++) OMP_private.Mk[ii][jj] += alpha[j];
                }//end i
            }//end ii
            
            for (int ic=0; ic<Ncv; ic++) {
                for (int jc=ic; jc<Ncv; jc++) {
                    for (int i_0=0; i_0 <OMP_private.n_diss_terms; i_0++){
                        ii = OMP_private.T_dissip_index_0[i_0]; // we need only terms with non zero T
                        jj = OMP_private.T_dissip_index_1[i_0];

                        Pv[ik_pr][ic][jc] -= OMP_private.hUk[ik_pr][jc][jj] * \
                                        conj(OMP_private.Uk[ik_pr][ic][ii]) * \
                                        T[ii][jj] *OMP_private.Mk[ii][jj];//part of dissipation coming from the electron equations
                    }//end sum over i_0 (in the notes i)
                    complexd Tii=0.;
                    for (int i=0; i<(Nb[0]+Nb[1]); i++) {
                        Tii += OMP_private.hUk[ik_pr][jc][i] *conj(OMP_private.hUk[ik_pr][ic][i])* T[i][i];
                    }
                    Pv[ik_pr][ic][jc] += Tii; //part of dissipation coming only from the holes
                    
                }//end jc
            }//end ic
        } // end if (OMP_private.n_diss_terms > 0)
        

        //-------- GRADIENT ----------//
        if (E_non0)
        { // we don't have gradient term when field EF is close to 0
            for (int ic=0; ic<Ncv; ic++)
            {
                for (int jc=ic; jc<Ncv; jc++) 
                {
                    for(int ishell=0; ishell<Bvector.size(); ishell++){
                        for(int ib=0; ib<Bvector[ishell].size(); ib++){
                           for(int i =0; i<3; i++){
                            OMP_private.P_grad[ik][ic][jc]+= Weigths[ishell]*EF[0][i]*Bvector[ishell][ib][i]*
                                (P[GradientIndex[ik][ishell][ib]][ic][jc]);//-P[ik][ic][jc]);
                           }
                       }
                   }
                   Pv[ik_pr][ic][jc] +=  OMP_private.P_grad[ik][ic][jc]; // SIGN!!!!!
                   OMP_private.P_grad[ik][jc][ic] = conj(OMP_private.P_grad[ik][ic][jc]);
                }//end jc
            }//end ic
        } //fi (Emax > E_eps)

    } //end ik




    for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++)
    {
        for (int ic=0; ic<Ncv; ic++)
        {
            Pv[ik_pr][ic][ic]=real(Pv[ik_pr][ic][ic]);
            //#pragma omp simd 
            for (int jc=ic+1; jc<Ncv; jc++)
            {
                Pv[ik_pr][jc][ic]=conj(Pv[ik_pr][ic][jc]);
            }
        }
    }
    #pragma omp barrier

} //End Runge_Kutta_Df_W



void Runge_Kutta_Ad(vec3x& P0,vec3x& P1,vec3x& Pv,
    double& dt, Private_omp_parameters& OMP_private)
{
    int Ncv=P0.n2();
    int ik_pr;
    for (int ik = OMP_private.begin_count; ik < OMP_private.end_count; ik++)
    {
        ik_pr = ik - OMP_private.begin_count;
        for (int ic=0; ic<Ncv; ic++)
        {
            if (OMP_private.Vectorization){
                #pragma omp simd  //#pragma ivdep
                for (int jc=0; jc<Ncv; jc++)
                {
                    P1[ik][ic][jc] = P0[ik][ic][jc] + Pv[ik_pr][ic][jc]*dt;
                }//End jc
            } else {
                for (int jc=0; jc<Ncv; jc++)
                {
                    P1[ik][ic][jc] = P0[ik][ic][jc] + Pv[ik_pr][ic][jc]*dt;
                }//End jc
            }

        }//End ic
    }// End iky
    #pragma omp barrier

}


void Runge_Kutta_Ac(vec3x& P0,vec3x& Pv,
    double& dt, Private_omp_parameters& OMP_private)
{
    int ik_pr;
    for (int ik = OMP_private.begin_count; ik < OMP_private.end_count; ik++)
    {
        ik_pr = ik - OMP_private.begin_count;
        for (int ic=0; ic<P0.n2(); ic++)
        {
            if (OMP_private.Vectorization){
                #pragma omp simd //#pragma ivdep
                for (int jc=0; jc<P0.n2(); jc++)
                {
                    P0[ik][ic][jc]+= Pv[ik_pr][ic][jc]*dt;
                }//End jc
            } else {
                for (int jc=0; jc<P0.n2(); jc++){
                    P0[ik][ic][jc]+= Pv[ik_pr][ic][jc]*dt;
                }//End jc
            }
        }//End ic
    }// End ikz
    #pragma omp barrier

}




/*
for (int ic=0; ic<Ncv; ic++){
    for (int ii=0; ii<Ncv; ii++){
        for (int jj=0; jj<Ncv; jj++){
            for (int jc=0; jc<Ncv; jc++){
                Xk_W[ic][jc] += \
                conj(OMP_private.Uk[ik_pr][ii][ic]) \
                                     *Xk_d[ii][jj]  \
                    *OMP_private.Uk[ik_pr][jj][jc];
            }
        }//end jj
    }//end ii
}//End ic
*/

/*
complexd Tii=0.;
//calculation of T[ic][jc][ii][ij]
for (int ii=0; ii<Ncv; ii++)
{
    for (int jj=0; jj<Ncv; jj++)
    {
        complexd Pij=0.;
        for (int i=0; i<Ncv; i++)
            for (int j=0; j<Ncv; j++)
                Pij+=conj(OMP_private.Uk[ik_pr][ii][i])*OMP_private.Uk[ik_pr][jj][j]*P[ik][i][j]; //this term depends on ii,jj,i,j
        Pv[ik][ic][jc]-=conj(OMP_private.Uk[ik_pr][jj][jc])*OMP_private.Uk[ik_pr][ii][ic]*T[ii][jj]*Pij; //part of dissipation coming from the electron equations
    }//end sum over jj (in the notes j)
}//end sum over ii (in the notes i)
for (int i=0; i<(Nb[0]+Nb[1]); i++) Tii+=conj(OMP_private.Uk[ik_pr][i][jc])*OMP_private.Uk[ik_pr][i][ic]*T[i][i];
Pv[ik][ic][jc] += Tii; //part of dissipation coming only from the holes
*/
