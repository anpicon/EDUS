// if the contdition for printing is true
 Emax = sqrt(EF_pr[0][0]*EF_pr[0][0] + EF_pr[0][1]*EF_pr[0][1] + EF_pr[0][2]*EF_pr[0][2]); // electric field modulus

#pragma omp master
{   
    
    if ((rank_ == root_rank) and ( it % (10*it_resolution) == 0))
    {
    cout << " it: " << it << " time: " << time_loop*time_au_fs << " fs, time fin: " << t_fin*time_au_fs << " fs" << endl; 
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

   if(rank_==0){
        PrintEF(fp_E, time_loop, EF_pr);
    }
}

print_matrices_in_kspace = false;
if (    
        (Coulomb_set.Print_new_band_dispersion or Diff_Eq.PrintPopulation) // if there are flag to print any type of matrices 
    and (time_loop > Diff_Eq.start_print_time + icont *Diff_Eq.step_print) 
    and (time_loop < Diff_Eq.end_print_time) 
    ){
        print_matrices_in_kspace = true;
        #pragma omp barrier
        #pragma omp master
        {
            icont++; // counter for snapshots
            if (rank_ == root_rank) {
                t_Ek <<   time_loop*time_au_fs << endl;
            }
        }      
    }


if (Coulomb_set.Coulomb_calc) {

    // creates a file that contains new bands energies 
    if (print_matrices_in_kspace){
        if (rank_ == root_rank){ // only root rank writes to file
            #pragma omp master
            {   
                
                Coulomb_set.Ek_file.close(); // close previously opened
                string name_file; 
                stringstream sname;
                sname.seekp(0,ios::beg); 
                sname << icont;
                name_file= "Output/Edia_" + sname.str() + ".txt";
                Coulomb_set.Ek_file.open(name_file.c_str());
                
                
            }
        }
    }

    PrintLossesMPI_exciton(fp_Loss, time_loop,
    Ncv, Nb[0], (Nb[0]+Nb[1]),
    P0, Coulomb_set, OMP_private, trig_k_omp, EF_pr,
    root_rank, rank_, print_matrices_in_kspace, n_cond, P_cond_max_loc_mpi);

    // #pragma omp master
    // {
    //     MPI_Bcast(& n_cond, 1, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
    //     Coulomb_set.n_cond = n_cond;
    // }
    
    // #pragma omp master
    //     {
    //     PrintLossesMPI(nk, nktot, fp_Loss, 
    //     Nb[0], (Nb[0]+Nb[1]), 
    //     P0, time_loop, Unitary, n_cond, P_cond_max_loc_mpi);
    //     }
    #pragma omp barrier

   




 } else{
    PrintLossesMPI_exciton(fp_Loss, time_loop,
    Ncv, Nb[0], (Nb[0]+Nb[1]),
    P0, Coulomb_set, OMP_private, trig_k_omp, EF_pr,
    root_rank, rank_, print_matrices_in_kspace, n_cond, P_cond_max_loc_mpi);
     // #pragma omp master
     //    {
     //    PrintLossesMPI(nk, nktot, fp_Loss, 
     //    Nb[0], (Nb[0]+Nb[1]), 
     //    P0, time_loop, Unitary, n_cond, P_cond_max_loc_mpi);
     //    }
    #pragma omp barrier
 }
   


#pragma omp master
{
    if (Print_and_sort_k){ // warning: function works only
    // inside master omp thread 
        print_k_MPI(kpt, nk); // print k- grid
    
        
        Create_sorted_indices_k_MPI(kpt, Diff_Eq); // create sorted k- grid

        // if (rank_ == root_rank){
        //     for (int ik = 0; ik < Diff_Eq.k_index.n1(); ik++){ // all wave vectors
        //         cout << "ik =" << ik << "  Diff_Eq.k_index= " << Diff_Eq.k_index[ik] << endl;
        //     }
        // }
        
    }    
   Print_and_sort_k = false;
                   
}



if(iCurrent)
{


    // calculate new energy gradioent because we renormed hamiltonian
    if(Coulomb_set.Coulomb_calc){
        int ik; // position in shared variable

        //gather all local arrays Hk_renorm into one shared array
        for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++) // loop over position in private variable
        {
            ik = ik_pr + OMP_private.begin_count; // coordinate in global array
            for (int ic=0; ic<Ncv; ic++){ //  over bands
                for (int jc=ic; jc<Ncv; jc++){ //
                    Coulomb_set.Hk_renorm_shared[ik][ic][jc] = OMP_private.Hk_renorm[ik_pr][ic][jc];
                }
            }
        }

        // exchange the borders of Coulomb_set.Hk_renorm_shared to calculate gradient:
        #pragma omp barrier
        #pragma omp master
        { // shared variables, only master change it
            MPI_Barrier(MPI_COMM_WORLD);
            Mpi_communication(Coulomb_set.Hk_renorm_shared,  message);
            MPI_Barrier(MPI_COMM_WORLD);
            GradientEnergy.fill(0);
        }
        #pragma omp barrier

        

        // calculate new energy gradient
        for(int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++){
            ik = ik_pr + OMP_private.begin_count; // coordinate in global array
            for(int i=0; i<3; i++){
                for(int ic=0; ic<Ncv; ic++){
                    for(int jc=ic; jc<Ncv; jc++){
                        for(int ishell=0; ishell<Bvector.size(); ishell++){
                            for(int ib=0; ib<Bvector[ishell].size(); ib++){
                                GradientEnergy[ik][ic][jc][i] += Weigths[ishell]*Bvector[ishell][ib][i]*
                                (Coulomb_set.Hk_renorm_shared[GradientIndex[ik][ishell][ib]][ic][jc]-Coulomb_set.Hk_renorm_shared[ik][ic][jc]);
                            }
                        }
                        GradientEnergy[ik][jc][ic][i] = conj(GradientEnergy[ik][ic][jc][i]);
                    }
                }
            }
        }
    }

    #pragma omp barrier

    Current1(fp_J1, 
        OMP_private, P0,
        GradientEnergy, time_loop, dk,
        nktot, det_a, J1);

    #pragma omp barrier

    if(Coulomb_set.Coulomb_calc){ // in case of Coulomb correlation, we insert renormed Hamiltonian
        Current2(fp_J2, 
        OMP_private, P0,
        Dipole, Coulomb_set.Hk_renorm_shared, time_loop, dk,
        nktot, det_a, J2);

    }else{
        Current2(fp_J2, 
        OMP_private, P0,
        Dipole, Hamiltonian, time_loop, dk,
        nktot, det_a, J2);
    }
    



}
        
#pragma omp master
 {
    if(iTAbs){
        PrintTransientAbsMPI(nk,nktot,fp_TAbs, P0, time_loop, EF_pr[1], dk, Nb[0], Coord_B::getJ(), Dipole);
    }
 }





// printing maximum values of population and it's derivatives in k
if (Coulomb_set.Coulomb_calc and ((it-iti) % (100*it_resolution) == 0)) {
    #pragma omp barrier
    #pragma omp master
    {   
    
        if (rank_ == root_rank){
            cout << " values of P0_dia[Ncv-1][Ncv-1] ";
        }
    }
    #pragma omp barrier
    print_max_min_conduct(Coulomb_set.P0_dia, 
    OMP_private, Coulomb_set, root_rank);



}



// plot populations
double scaled_t = round( time_loop*time_au_fs* pow(10,4)) / pow(10,4);
if ( Diff_Eq.PrintPopulation and print_matrices_in_kspace)// ((it-iti) % (500*it_resolution) == 0))
{


    string label = "P_cond_in_Bloch";
    int ik_pr; // transform into Bloch basis
    for (int ik = OMP_private.begin_count; ik < OMP_private.end_count; ik++){ // all local wave vectors
        ik_pr = ik - OMP_private.begin_count;
        
        // function transforms P[ik] with OMP_private.Uk[ik_pr] and writes it to the 
        // OMP_private.P_diag[ik]
        get_P_in_dia_vect(P0,  OMP_private.Uk, OMP_private.P_diag, Ncv, 
            ik, ik_pr, OMP_private);
    }
    double scaled_t = round( time_loop*time_au_fs* pow(10,4)) / pow(10,4);
    for (int ii = 0; ii<Ncv; ii++){
        Print_vec3x_MPI(OMP_private.P_diag, ii, ii, scaled_t,
         OMP_private,
         Coulomb_set, Diff_Eq,
         label, icont); // icont - the same number as in in Ec
    }

    label = "P_val_cond_in_Bloch";
    Print_vec3x_MPI(OMP_private.P_diag, Ncv-2, Ncv-1, scaled_t,
     OMP_private,
     Coulomb_set, Diff_Eq,
     label, icont); // icont - the same number as in in Ec


    label = "P_Wannier_";
    Print_vec3x_MPI(P0, Ncv-1, Ncv-1, scaled_t,
     OMP_private,
     Coulomb_set,Diff_Eq,
     label, icont); // icont - the same number as in in Ec

    label = "P_Wannier_";
    Print_vec3x_MPI(P0, Ncv-2, Ncv-1, scaled_t,
     OMP_private,
     Coulomb_set,Diff_Eq,
     label, icont); // icont - the same number as in in Ec

    label = "P_Wannier_";
    Print_vec3x_MPI(P0, Ncv-2, Ncv-2, scaled_t,
     OMP_private,
     Coulomb_set,Diff_Eq,
     label, icont); // icont - the same number as in in Ec
 


}


if (Print_initial_H_P){
    Diagonalize_unitary_vec3x(OMP_private.Hk, H_Eigen, OMP_private);

    string label = "H_Eigen_cond";
    Print_vec3x_MPI(H_Eigen, Ncv-1,Ncv-1, scaled_t,
     OMP_private,
     Coulomb_set, Diff_Eq,
     label, it);

    label = "H_Eigen_val";
    Print_vec3x_MPI(H_Eigen, Ncv-2,Ncv-2, scaled_t,
     OMP_private,
     Coulomb_set,Diff_Eq,
     label, it);

    label = "H_Eigen_val_cond";
    Print_vec3x_MPI(H_Eigen, Ncv-2,Ncv-1, scaled_t,
     OMP_private,
     Coulomb_set,Diff_Eq,
     label, it);


    label = "H_Wannier_cond";
    Print_vec3x_MPI(Hamiltonian, Ncv-1,Ncv-1, scaled_t,
     OMP_private,
     Coulomb_set,Diff_Eq,
     label, it);

    label = "H_Wannier_val_cond";
    Print_vec3x_MPI(Hamiltonian, Ncv-2,Ncv-1, scaled_t,
     OMP_private,
     Coulomb_set,Diff_Eq,
     label, it);

    label = "H_Wannier_val";
    Print_vec3x_MPI(Hamiltonian, Ncv-2,Ncv-2, scaled_t,
     OMP_private,
     Coulomb_set,Diff_Eq,
     label, it);


    label = "P_initial_Wannier_";
    Print_vec3x_MPI(P0, Ncv-1, Ncv-1, scaled_t,
     OMP_private,
     Coulomb_set,Diff_Eq,
     label, it); // icont - the same number as in in Ec

    label = "P_initial_Wannier_";
    Print_vec3x_MPI(P0, Ncv-2, Ncv-1, scaled_t,
     OMP_private,
     Coulomb_set,Diff_Eq,
     label, it); // icont - the same number as in in Ec

    label = "P_initial_Wannier_";
    Print_vec3x_MPI(P0, Ncv-2, Ncv-2, scaled_t,
     OMP_private,
     Coulomb_set,Diff_Eq,
     label, it); // icont - the same number as in in Ec


    label = "Uk_";
    Print_vec3x_MPI(Unitary, Ncv-1,Ncv-1, scaled_t,
     OMP_private,
     Coulomb_set,Diff_Eq,
     label, it);

    label = "Uk_";
    Print_vec3x_MPI(Unitary, Ncv-1,Ncv-2, scaled_t,
     OMP_private,
     Coulomb_set,Diff_Eq,
     label, it);
    
    #pragma omp master
    {   
        Print_initial_H_P = false;
    }

}