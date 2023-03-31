#include <mpi.h>


struct Message
{
    bool send;
    int buf;
    int count;
    int partner_rank;
    int tag;
};



int rank_, num_procs;
#include "Source_Main/include_headers.cpp"
#include "Observables_MPI.h"
#include "Organize_kspaceMPI.h"




int main (int argc, char* argv[])
{
    
    #include "Source_Main/init_mpi.cpp"
       stringstream srank; srank<<rank_;

    if(rank_==0)
    {
        printf("Running the program with %3i cpu, ", num_procs);
        #pragma omp parallel sections
        {
            printf("%3i threads \n", omp_get_num_threads());
        }
    }
    #include "Source_Main/variables.cpp"
    #include "Source_Main/input.cpp" //we define the iMode=TB, W, or CY
    
    if(iMode == "W")
    {
        cout << "Care. you are using wannier input and the crystal directions are being changed 012->201\n";
        vec2d atemp(3,3);
        for(int i=0; i<3; i++)
            for(int j=0; j<3; j++)
                atemp[i][j] = a[i][j];
        for(int j=0; j<3; j++)
        {
            a[0][j] = atemp[2][j];
            a[1][j] = atemp[0][j];
            a[2][j] = atemp[1][j];
        }         
        VecProd(b[0],a[1],a[2]); VecProd(b[1],a[2],a[0]); VecProd(b[2],a[0],a[1]);
        double value=DotProd(a[0],b[0])/2./pi; b[0][0]/=value; b[0][1]/=value; b[0][2]/=value;
        value=DotProd(a[1],b[1])/2./pi; b[1][0]/=value; b[1][1]/=value; b[1][2]/=value;
        value=DotProd(a[2],b[2])/2./pi; b[2][0]/=value; b[2][1]/=value; b[2][2]/=value;
    }       
    int nktot = Nk[0]*Nk[1]*Nk[2];
    vector<Message> message;
    
    // here we calculate indices, grid kpt and nk on each CPU
    CalculateIndicesMPI(kpt, GradientIndex, Bvector, 
    dk, nk, message);

    #include "Source_Main/test_matrix_k.cpp"
    #include "Source_Main/allocate_arr.cpp"
    #include "Source_Main/open_files.cpp"
    MPI_Barrier(MPI_COMM_WORLD);
    #include "Source_Main/save_HUD.cpp"
    //#include "Source_Main/print_HUD.cpp"
    #include "Source_Main/runge_kutta_init.cpp"
    int cnt=0;
    if (t_fin > 0){
        itfi = t_fin / dt; 
    } else {
        t_fin = itfi*dt;
    }
    double dt0 = dt;
    iti /= 5;
    // itfi /= 2;
    double time_loop;

    // // START WITH NON-EQUILIBRIUM STATE
    // for (int ik = 0; ik < nk; ik++){ // all wave vectors
    //     double D_P = 1e-7;
    //     vec2x P_exiton(Ncv, Ncv);
    //     P_exiton.fill(0.);
    //     P_exiton[Ncv-1][Ncv-1] += D_P;
    //     P_exiton[Ncv-2][Ncv-2] -= D_P;
    //     P_exiton[Ncv-1][Ncv-2] += sqrt(D_P)*exp(c1* pi / 2);
    //     P_exiton[Ncv-2][Ncv-1]  = conj(P_exiton[Ncv-1][Ncv-2]);

    //     for ()

    //     // P0[ik][Ncv-1][Ncv-1] += D_P;
    //     // P0[ik][Ncv-2][Ncv-2] -= D_P;

    // }// START WITH NON-EQUILIBRIUM STATE

MPI_Barrier(MPI_COMM_WORLD);
#pragma omp parallel // num_threads(3)
{
    // struct with omp parameters which will be 
    //private for each omp thread
    Private_omp_parameters OMP_private;

    
    OMP_private.thr_id = omp_get_thread_num();
    OMP_private.thr_total = omp_get_num_threads();
    // in every thread define border, where we operate
    OMP_private.lenght_k = int(nk/OMP_private.thr_total);
    int leftover = nk % OMP_private.thr_total;
    OMP_private.end_count = ( OMP_private.lenght_k *OMP_private.thr_id +  OMP_private.lenght_k);
    OMP_private.begin_count = ( OMP_private.lenght_k *OMP_private.thr_id);
    if (OMP_private.thr_id < leftover){
        OMP_private.begin_count += OMP_private.thr_id;
        OMP_private.end_count += OMP_private.thr_id + 1;
    }
    else{
        OMP_private.begin_count += leftover;
        OMP_private.end_count += leftover;
    }

    OMP_private.lenght_k = OMP_private.end_count - OMP_private.begin_count;
    
    OMP_private.id_masters.resize(OMP_private.thr_total);
    std::fill(OMP_private.id_masters.begin(), OMP_private.id_masters.end(), false);
    
    // if (OMP_private.thr_total < 4){
    // std::cout << " you need more threads! "  << endl;
    // exit(1);
    // }
    OMP_private.id_masters[OMP_private.thr_id] = true;
    #pragma omp critical
    {
    std::cout << " leftover " << leftover << " lenght_k_omp: " << OMP_private.lenght_k<< endl;
    std::cout << " nk on node: " << nk << " begin_count: " << OMP_private.begin_count << " end_count: " << OMP_private.end_count << endl << endl;
    }
    


    OMP_private.integrWeight.resize(OMP_private.lenght_k);
    OMP_private.P0.resize(OMP_private.lenght_k,  Ncv,Ncv);
    OMP_private.P1.resize(OMP_private.lenght_k,  Ncv,Ncv);
    OMP_private.P2.resize(OMP_private.lenght_k,  Ncv,Ncv);
    OMP_private.P2_dia.resize(OMP_private.lenght_k,  Ncv,Ncv);
    OMP_private.Pv.resize(OMP_private.lenght_k,  Ncv,Ncv);

    // population at previous step
    OMP_private.P_W_prev.resize(OMP_private.lenght_k,  Ncv,Ncv);
    OMP_private.P_dia_prev.resize(OMP_private.lenght_k,  Ncv,Ncv);
    OMP_private.P_dia_prev.fill(0.0);
    
    OMP_private.P_static.resize(OMP_private.lenght_k,  Ncv,Ncv);
    
    // P_diag - population in basis where equilibrium energy is diagonal
    // P_eigen -  population in basis where non-equlibrium hamiltonian is diagonal, 
    // different from P_diag only if we have Coulomb term
    OMP_private.P_diag.resize(nk,  Ncv,Ncv);
    OMP_private.P_diag.fill(0.0);
    OMP_private.P_eigen.resize(nk,  Ncv,Ncv);
    OMP_private.P_eigen.fill(0.0);
    OMP_private.P_grad.resize(nk,  Ncv, Ncv);



    OMP_private.k.resize(OMP_private.lenght_k, 3);
    OMP_private.Hk.resize(OMP_private.lenght_k, Ncv,Ncv);
    OMP_private.Uk.resize(OMP_private.lenght_k, Ncv,Ncv);
    OMP_private.hUk.resize(OMP_private.lenght_k, Ncv,Ncv);
    OMP_private.Dk.resize(3,OMP_private.lenght_k, Ncv, Ncv);
    OMP_private.Mk.resize(Ncv, Ncv);
    OMP_private.Ak.resize(Ncv);
    OMP_private.Hk_renorm.resize(OMP_private.lenght_k, Ncv,Ncv);
    // OMP_private.Dk.fill(0.0);

    OMP_private.Dk0.resize(nk, Ncv, Ncv);
    OMP_private.Dk1.resize(nk, Ncv, Ncv);
    OMP_private.Dk2.resize(nk, Ncv, Ncv);
    

          

    // creating local variables

    for (int ik = OMP_private.begin_count; ik < OMP_private.end_count; ik++){ // all local wave vectors
        OMP_private.k[ik - OMP_private.begin_count][1] = kpt[ik][1];
        OMP_private.k[ik - OMP_private.begin_count][2] = kpt[ik][2];
        OMP_private.integrWeight[ik - OMP_private.begin_count] = integrWeight[ik];
        // P0[ik][Ncv-1][Ncv-1] += pow(10,-14);
        // P0[ik][Ncv-2][Ncv-2] -= pow(10,-14);


        for (int jc=0; jc<Ncv; jc++){ //  over bands
            if (jc < (Ncv-1)) OMP_private.P_dia_prev[ik - OMP_private.begin_count][jc][jc] = 1.0;
            for (int ic=0; ic<Ncv; ic++){ //

                OMP_private.P0[ik - OMP_private.begin_count][ic][jc] = P0[ik][ic][jc];
                OMP_private.Hk[ik - OMP_private.begin_count][ic][jc] = Hamiltonian[ik][ic][jc];
                OMP_private.Hk_renorm[ik - OMP_private.begin_count][ic][jc] = OMP_private.Hk[ik - OMP_private.begin_count][ic][jc];
                OMP_private.Uk[ik - OMP_private.begin_count][ic][jc] = Unitary[ik][ic][jc];
                OMP_private.hUk[ik - OMP_private.begin_count][ic][jc] = conj(Unitary[ik][jc][ic]);

                OMP_private.P_W_prev[ik - OMP_private.begin_count][ic][jc] = OMP_private.P0[ik - OMP_private.begin_count][ic][jc];

                OMP_private.Dk0[ik][ic][jc] = Dipole[ik][ic][jc][0];
                OMP_private.Dk1[ik][ic][jc] = Dipole[ik][ic][jc][1];
                OMP_private.Dk2[ik][ic][jc] = Dipole[ik][ic][jc][2];
                
                for (int id=0; id< 3; id++){
                    OMP_private.Dk[id][ik - OMP_private.begin_count][ic][jc] = Dipole[ik][ic][jc][id];
                }
            }
        }
    }

    // create_noneq_population(OMP_private, P0, H_Eigen); // creates initial population in conduction band
    //     //modifies P0
    // if(true){ // to make variables introduced inside local
    //     string label = "P_Eigen_initial_cond";
    //     int it=0;
    //     double scaled_t = 0.0;

    //     Print_vec3x_MPI(OMP_private.P_eigen, Ncv-1,Ncv-1, scaled_t,
    //      OMP_private,
    //      Coulomb_set,
    //      label, it);


    //     label = "P_Eigen_initial_val";
    //     Print_vec3x_MPI(OMP_private.P_eigen, Ncv-2,Ncv-2, scaled_t,
    //      OMP_private,
    //      Coulomb_set,
    //      label, it);
    // }

    // struct with trigonometric functions
    trig_coefficients trig_k_omp;
    get_trig_coef(trig_k_omp, OMP_private.k, Ncut, OMP_private.lenght_k);

    vector<vec1d> EF_pr(2); // private for each thread, 
    //so they don't have to acces one block on memoryP_static
    for(int i=0; i<2; i++) {EF_pr[i].resize(3); EF_pr[i].fill(0);}
    
    #pragma omp barrier // sinchronise threads
    if (Coulomb_set.Coulomb_calc)
    {// calculate set of coefficients 
        // in equilibrium to obtain the Coulomb terms
        #pragma omp master
        {
            Mpi_communication(P0,  message);
            MPI_Barrier(MPI_COMM_WORLD);
        }
        #pragma omp barrier
        
        if(Coulomb_set.Wannie_basis){

            Calculate_X_coefficients_MPI(P0, OMP_private,  
                Ncv, Coulomb_set, root_rank, rank_, trig_k_omp);

        } else if(Coulomb_set.Diagonal_basis){  
            
            // create equilibrium P in diagonal basis:
            int Nbands = Nb[0]+Nb[1]; // number of valence bands
            OMP_private.P_diag.fill(0.0);
            OMP_private.P_static.fill(0.0);
            for (int ik = OMP_private.begin_count; ik < OMP_private.end_count; ik++){ // all local wave vectors
                for (int ic=0; ic<Nbands; ic++){ // fill only valence bands
                    OMP_private.P_diag[ik][ic][ic] = 1.0;
                    OMP_private.P_static[ik - OMP_private.begin_count][ic][ic] = 1.0;
                }
            } // k loop

            Calculate_X_coefficients_MPI(OMP_private.P_diag, OMP_private,  
                Ncv, Coulomb_set, root_rank, rank_, trig_k_omp);
        }
        
        #pragma omp master
        {
            for (int m = 0; m < Ncut*Ncut*Ncv*Ncv; m++ ){ //row summation over Fourier series
                Coulomb_set.X_cc0[m] =Coulomb_set.X_cc[m];
                Coulomb_set.X_cs0[m] =Coulomb_set.X_cs[m];
                Coulomb_set.X_sc0[m] =Coulomb_set.X_sc[m];
                Coulomb_set.X_ss0[m] =Coulomb_set.X_ss[m];
            }    // end of creation equilibrium coulomb part
        }
        
    } // fi Coulomb
    
    #pragma omp master
    {   
        MPI_Barrier(MPI_COMM_WORLD);
        for (int i = 0; i < num_procs; i++){
          MPI_Barrier(MPI_COMM_WORLD);
            // Cheking kpt borders
            if (rank_ == i){
                cout << " ******************** "  << endl;
                cout << " my_rank " <<  rank_ << endl;
                // cout << " Coulomb_set.X_cc0[0] " << Coulomb_set.X_cc0[0] << endl; 
                // cout << " Coulomb_set.X_cc0[end] " << Coulomb_set.X_cc0[Ncut*Ncut*Ncv*Ncv -1] << endl << endl << endl; 
            }

          MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    vec2x Delta_P; // need it for Taylor Solver
    
    vector<vec3x> derivativesMatrices(TaylorOrder + 1); // need it for Taylor Solver
    vector<vector<int>> derivatIndex(TaylorOrder + 1);  // need it for Taylor Solver

    


    if (Diff_Eq.Taylor){ // if true init matrices for drivatives
        init_Taylor_matrices(derivativesMatrices, derivatIndex, 
            TaylorOrder, OMP_private.lenght_k, Ncv);

        Delta_P.resize(OMP_private.lenght_k,  Ncv*Ncv);
        #pragma omp master
        {
            dt_prev.fill(dt); // matrix with previous dt values
        }

    }

    if (Diff_Eq.const_dt_evolution){
        #include "Source_Main/const_dt_evolution.cpp"
    }
    if (Diff_Eq.dynamical_dt_evolution){
        #include "Source_Main/dyn_dt_evolution.cpp"
    }


} // parallel section

#include "Source_Main/end_program.cpp"

} //END MAIN



// compilation on desktop:
// mpicxx main_MPI.cpp  -fopenmp -I/usr/include/python2.7 -lpython2.7 -larmadillo -I./headers
// export OMP_NUM_THREADS=7
//  mpirun -np 1 ./a.out input.txt
