
void Initial_Population(vec3x& P0,int& nk, vec2d& kpt,  int nktotal, vec3x& Uk,double FermiE, int Nbands)
{
    //Initial population: Fermi level
    vec3x P_bloch(nk, P0.n2(), P0.n3()); P_bloch.fill(0.);
    //density matrix at the beginning (Bloch)
    #pragma omp parallel for schedule(dynamic)
    for (int ik=0; ik<nk; ik++) 
    {
        for (int ic=0; ic<Nbands; ic++) 
            P_bloch[ik][ic][ic]=1.;
        // for(int ic=Nbands; ic<P0.n2(); ic++)
        //     P_bloch[ik][ic][ic]=0.;
        for (int ia=0; ia<P0.n2(); ia++) 
            for(int ib=0; ib<P0.n2(); ib++)
                for(int ii=0; ii<P0.n2(); ii++)
                    for(int ij=0; ij<P0.n2(); ij++)
                        // this indexes are arranged in a weird way. pay attention to it!! rows and columns are swapped in P0!!
                        P0[ik][ia][ib]+=Uk[ik][ii][ia]*P_bloch[ik][ii][ij]*conj(Uk[ik][ij][ib]);
                        /* let Ak be the conjugate transpose of Uk. then here we do

                        P0[ik][ia][ib]+=Uk[ik][ii][ia]*P_bloch[ik][ii][ij]*conj(Uk[ik][ij][ib])
                                       =Uk[ik][ii][ia]*P_bloch[ik][ii][ij]*Ak[ik][ib][ij]
                                       =Ak[ik][ib][ij]*P_bloch[ik][ii][ij]*Uk[ik][ii][ia]

                        which in conventional matrix notation would be P0[ik][ib][ia].

                        another way to look at it would be

                        P0[ik][ia][ib]+=Uk[ik][ii][ia]*P_bloch[ik][ii][ij]*conj(Uk[ik][ij][ib])
                                       =conj(Ak[ik][ia][ii])*P_bloch[ik][ii][ij]*conj(Uk[ik][ij][ib])

                        which is conj(Ak P_bloch Uk)[ia][ib] 

                        therefore, when P0 is concerned in future operations, we must write following the conventional matrix notation
                        and swap P0s rows and columns at the end */
       
    }
}


// diagonalize vec3x
// for example Hamiltonian from Wannier to Eigenstate basis
// works in OMP multithread part
void Diagonalize_unitary_vec3x(
vec3x& Hamiltonian_nondia, // variable we want to diagonalize
vec3x& H_Eigen, // variable we store result
Private_omp_parameters& OMP_private){
    int Ncv =  Hamiltonian_nondia.n2();

    #pragma omp barrier
    #pragma omp master
    {
        H_Eigen.fill(0);
    }
    #pragma omp barrier

    vec epsilon; // temporary things for armadillo
    cx_mat Uk_arm; // temporary things for armadillo
    cx_mat Hk_arm; // temporary things for armadillo
    epsilon.zeros(Ncv);
    Uk_arm.zeros(Ncv, Ncv);
    Hk_arm.zeros(Ncv, Ncv);
    bool Private_index_Hamiltonian_nondia = false;
    if (Hamiltonian_nondia.n1() < H_Eigen.n1()){ // means that we 
        //use different k indexeces in case when we have private array of Hamiltonian_nondia
        // or shares in OMP threading
        Private_index_Hamiltonian_nondia = true;

    }
    int ik;
    int index_Hamiltonian_nondia;
    // loop in positions of private variable
    for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++)
    {
        ik = ik_pr + OMP_private.begin_count;
        index_Hamiltonian_nondia = ik; // global index for shared variable
        if(Private_index_Hamiltonian_nondia){ // private index for private variable
            index_Hamiltonian_nondia = ik_pr;
        }

        for (int ic=0; ic<Ncv; ic++){// summation over bands
            for (int jc=0; jc<Ncv; jc++){ // summation over bands
                Hk_arm(ic, jc) = Hamiltonian_nondia[index_Hamiltonian_nondia][ic][jc]; // copy to the temporary armadillo variable
            }
        }

        eig_sym(epsilon, Uk_arm, Hk_arm); // diagonalization

        // copy result to the matrix we can print
        for (int ic=0; ic<Ncv; ic++){// summation over bands
                H_Eigen[ik][ic][ic] = epsilon(ic);   
        }

    }


    // string label = "H_Eigen_cond";
    // Print_vec3x_MPI(H_Eigen, Ncv-1,Ncv-1, scaled_t,
    //  OMP_private,
    //  Coulomb_set,
    //  label, it);


    // label = "H_Eigen_val";
    // Print_vec3x_MPI(H_Eigen, Ncv-2,Ncv-2, scaled_t,
    //  OMP_private,
    //  Coulomb_set,
    //  label, it);

}

//works in OMP multithread region
void create_noneq_population(
    Private_omp_parameters& OMP_private, 
    vec3x& P0, // variable to store new population
    vec3x& H_Eigen// eigenstate basis energy
    ){
    
    Diagonalize_unitary_vec3x(OMP_private.Hk, H_Eigen, OMP_private);


    OMP_private.P_eigen.fill(0);
    int Ncv = OMP_private.P_eigen.n2();

    int ik;
    double D_P;// excited population
    double D_F; // excited phase
    double gap = 6.2 * energy_eV_au;
    double D_E_k; // energy gap in a particular k point
    // loop in positions of private variable
    for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++) {
        ik = ik_pr + OMP_private.begin_count;

        D_E_k = abs(H_Eigen[ik][Ncv-1][Ncv-1] - H_Eigen[ik][Ncv-2][Ncv-2]);
        D_P = 0.;

        if (((D_E_k - gap) < 0.5 * energy_eV_au)  and (OMP_private.k[ik_pr][1] > 0 )){
            D_P = 1e-6 * exp(-100000* pow(abs(D_E_k - gap),2));
            D_F = 2 * pi *D_P / 1e-6;
        }



        OMP_private.P_eigen[ik][Ncv-1][Ncv-1] = D_P;
        OMP_private.P_eigen[ik][Ncv-2][Ncv-2] = 1.0 - D_P;
        // for non dia terms:
        //<a+_v a_c> <a+_c a_v> = <a+_c a_c>
        OMP_private.P_eigen[ik][Ncv-1][Ncv-2] = sqrt(D_P)*exp(c1* D_F);
        OMP_private.P_eigen[ik][Ncv-2][Ncv-1] = conj(OMP_private.P_eigen[ik][Ncv-1][Ncv-2]);


    }

    #pragma omp barrier
    #pragma omp master
    {
        P0.fill(0);
    }
    #pragma omp barrier


    for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++) {
        ik = ik_pr + OMP_private.begin_count;

        for (int ia=0; ia<P0.n2(); ia++) 
            for(int ib=0; ib<P0.n2(); ib++)
                for(int ii=0; ii<P0.n2(); ii++)
                    for(int ij=0; ij<P0.n2(); ij++)
                        P0[ik][ia][ib]+=OMP_private.Uk[ik_pr][ii][ia]*OMP_private.P_eigen[ik][ii][ij]*conj(OMP_private.Uk[ik_pr][ij][ib]);

    }
}