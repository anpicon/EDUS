//Gather all data from treads and print
//
void Print_vec3x_MPI(vec3x& P, int ic, int jc,  double& time,
     Private_omp_parameters& OMP_private,
     Coulomb_parameters& Coulomb_set, methods_Diff_Eq &Diff_Eq,
     string& label, int it){
    #pragma omp barrier
    #pragma omp master
    { 
        MPI_Barrier(MPI_COMM_WORLD);
        Coulomb_set.Pk_shared.fill(0.0);
        
        string name_file; 
        
        stringstream sname;
        sname.seekp(0,ios::beg); 
        if (rank_ == 0){ 
            // sname << time;
            sname << ic << jc <<"_" <<  it;

            name_file= "Output/" + label + "_" + sname.str() + ".txt";
            cout<< name_file << endl;
            Coulomb_set.P_stream.open(name_file.c_str());
        }
    }

    #pragma omp barrier
    int ik; // position in shared variable

    for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++){
        ik = ik_pr + OMP_private.begin_count;

        // shared array to print in file 
        Coulomb_set.Pk_shared[ik] = P[ik][ic][jc];
    }
    
    #pragma omp barrier
    #pragma omp master
    { // rearrange k- points for printing
        vec1x Array_temp(Coulomb_set.Pk_shared.n1());

        int nk_y = Diff_Eq.Nk_node_MPI / (Diff_Eq.Nk_1 * Diff_Eq.Nk_0);
        int nk_x = Diff_Eq.Nk_1 * Diff_Eq.Nk_0;
        if (Diff_Eq.Nk_1 > Diff_Eq.Nk_2){
            nk_y = Diff_Eq.Nk_node_MPI / (Diff_Eq.Nk_1 * Diff_Eq.Nk_2);
            nk_x = Diff_Eq.Nk_1 * Diff_Eq.Nk_2;
        }


        for (int ik = 0; ik < Diff_Eq.Nk_node_MPI; ik++){ // all wave vectors
            int i_kx = ik / nk_x;
            int i_ky = ik % nk_x;
            // if (rank_ == 0) cout << "i_kx=" << i_kx <<  " i_ky=" << i_ky << " Diff_Eq.k_index[i_kx]=" << Diff_Eq.k_index[i_kx] << "   Diff_Eq.k_index[i_kx]*nk_x + i_ky=" << Diff_Eq.k_index[i_kx]*nk_x + i_ky <<  endl;
            Array_temp[ik] = Coulomb_set.Pk_shared[Diff_Eq.k_index[i_kx]*nk_x + i_ky];
        }
        // if (rank_ == 0) cout << "Copy arr" << endl;
        for (int ik = 0; ik < Diff_Eq.Nk_node_MPI; ik++){ // all wave vectors
            Coulomb_set.Pk_shared[ik] = Array_temp[ik];
        }

    }


    #pragma omp barrier

    #pragma omp master
    {   
        MPI_Barrier(MPI_COMM_WORLD);
        int Nk_node = Diff_Eq.Nk_node_MPI; // every MPI process gets number of points on node
        vec1x * Pk_print = & Coulomb_set.Pk_shared; // pointer to object
        vec1x Pk_receive; // variable, where we store what we will receive
        Pk_receive.resize(Nk_node);
        vec1x Array_print(Diff_Eq.Nk_total_All_Grid);
        Array_print.fill(0.0);
        int ik_MPI = 0;
        for (int i_proc = 0; i_proc < num_procs; i_proc++){
            MPI_Barrier(MPI_COMM_WORLD);
            if (i_proc > 0){
                // sending array to 0 node and then print

                // at first process 0 finds out how much data it will receive
                if (rank_ == i_proc){ // send
                    // destination - 0
                    MPI_Send(&Nk_node, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                }
                if (rank_ == 0){ // receive
                    // source - process number i_proc
                    MPI_Recv(&Nk_node, 1, MPI_INT, i_proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    Pk_receive.resize(Nk_node);
                }

            }
            MPI_Barrier(MPI_COMM_WORLD);
            if (i_proc > 0){
                if (rank_ == i_proc){ // send array
                    // destination - 0
                    MPI_Send(& Coulomb_set.Pk_shared[0], Nk_node, MPI_C_DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD);
                }
                if (rank_ == 0){ // receive
                    // source - process number i_proc
                    Pk_receive.fill(0.0);
                    MPI_Recv(& Pk_receive[0], Nk_node, MPI_C_DOUBLE_COMPLEX, i_proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            } // if (i_proc > 0)
            MPI_Barrier(MPI_COMM_WORLD); 
            if (rank_ == 0){ // only rank 0 prints!
                for(int ik = 0; ik < Nk_node; ik++){
                    if (i_proc > 0) Array_print[ik_MPI] = Pk_receive[ik];
                    if (i_proc == 0) Array_print[ik_MPI] = Coulomb_set.Pk_shared[ik];
                    ik_MPI+=1;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // first we create array, than print it, to avoid printing errors 
        if (rank_ == 0){ // only rank 0 prints!
            for(int ik = 0; ik < Diff_Eq.Nk_total_All_Grid; ik++){ 
                Coulomb_set.P_stream << real(Array_print[ik]) << " ";
                    Coulomb_set.P_stream  << imag(Array_print[ik]) << endl;
            }

            Coulomb_set.P_stream.close();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    #pragma omp barrier

}








// prints population
void PrintLossesMPI(int& nk,int& nktot,ofstream& fp_Loss, 
    int& Nch, int Ncc,vec3x& P0,double& time,vec3x& U, double &  n_cond, 
    double &  P_cond_max_loc_mpi)
{
    vec1d sumNcv(P0.n2());
    CalculateLosses(nk,sumNcv,Nch,Ncc,P0,time,U, P_cond_max_loc_mpi);

    vec1d global_sum(P0.n2());

    MPI_Reduce(&sumNcv[0], &global_sum[0], (int)P0.n2(), MPI_DOUBLE, MPI_SUM, 0,
           MPI_COMM_WORLD);    


    if(rank_==0)
    {
        fp_Loss << setw(15) << setprecision(8) << time*time_au_fs;
        for(int ic=0; ic<P0.n2(); ic++) fp_Loss << setw(15) << setprecision(8) <<global_sum[ic]/nktot;
        fp_Loss << endl;
        n_cond = global_sum[P0.n2() -1]/nktot;    
    }
    MPI_Barrier(MPI_COMM_WORLD);
}










// makes diagonalisation into exciton basis prints population
void PrintLossesMPI_exciton(ofstream& fp_Loss, double& time, 
    int Ncv, int& Nch, int Ncc, vec3x& P0,
    Coulomb_parameters& Coulomb_set, Private_omp_parameters& OMP_private,
    trig_coefficients & trig_k, vector<vec1d>&  EF, 
    int root_rank, int my_rank, bool print_EK, double &  n_cond, 
    double &  P_cond_max_loc_mpi)
{

    vec1d N_dia(Ncv);
    vec1d P_k_dia(Ncv);
    vec1d P_k_offdia(Ncv);

    vec1d H_k_dia(Ncv);
    vec1d sumNcv(Ncv);
    vec2x Hk(Ncv, Ncv);
    vec3x Uk_exc(1, Ncv, Ncv); // 1 - for unification, so we can use function get_P_in_dia
    P_k_dia.fill(0);
    P_k_offdia.fill(0);


    // No need to integrate with good precision here
    double nktot = Coulomb_set.N_BZ_points_total;



    int ik; // position in shared variable
    double P_cond_max_loc_omp = 0.0;
    double P_cond_min_loc_omp = 1.0;
    #pragma omp barrier
    #pragma omp master
    {
        P_cond_max_loc_mpi = 0.0;
        Coulomb_set.P_min = 1.0;
    }
    #pragma omp barrier
    // Now we have Hk and diagonalize it with armadillo
    // we don't need much precision or control phases here

    vec epsilon; // temporary things for armadillo
    cx_mat Uk_arm; // temporary things for armadillo
    cx_mat Hk_arm; // temporary things for armadillo
    epsilon.zeros(Ncv);
    Uk_arm.zeros(Ncv, Ncv);
    Hk_arm.zeros(Ncv, Ncv);

    // loop in positions of private variable
    for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++)
    {
        ik = ik_pr + OMP_private.begin_count;
        
        if (Coulomb_set.Coulomb_calc){
            // Coulomb coefficients into Hamiltonian
            for (int ic=0; ic<Ncv; ic++){// summation over bands
                for (int jc=0; jc<Ncv; jc++){ // summation over bands
                    Hk_arm(ic, jc) = OMP_private.Hk_renorm[ik_pr][ic][jc];
                }
            }

            eig_sym(epsilon, Uk_arm, Hk_arm); // diagonalization

                    // print Hk(t) - need to put here

            for (int ic=0; ic < Ncv; ic++){
                for (int jc = 0; jc < Ncv; jc++){
                    Uk_exc[0][ic][jc] = conj(Uk_arm(jc, ic));
                }
                if (print_EK )  Coulomb_set.Ekdia[ik][ic] = epsilon(ic);
            }
     
        } else { // No Coulomb coefficients
                    
            for (int ic=0; ic < Ncv; ic++){
                for (int jc = 0; jc < Ncv; jc++){
                    Uk_exc[0][ic][jc] = OMP_private.Uk[ik_pr][ic][jc];
                }
            }
        }

        // here ik_U is to make function work both in RungeKutta.h and here
        // I also had to introduce Uk_exc so this function can work both in exciton and non-exciton case
        int ik_U = 0;
        get_P_in_dia_vect(P0, Uk_exc, OMP_private.P_eigen, Ncv,  
            ik, ik_U, OMP_private);

        // P_k_offdia[Ncv-1] += abs(OMP_private.P_eigen[ik][Ncv-1][Ncv-2]); // valence - conduction exciton
        // if (Ncv>2) P_k_offdia[Ncv-2] += abs(OMP_private.P_eigen[ik][Ncv-2][Ncv-3]); // core - valence
        // if (Ncv>2) P_k_offdia[Ncv-3] += abs(OMP_private.P_eigen[ik][Ncv-1][Ncv-3]); // core- conduction

        for (int ic=0; ic < Ncv; ic++){
            double Pik = real(OMP_private.P_eigen[ik][ic][ic]); 
            

            P_k_dia[ic] += Pik;


            if (ic < Ncc){
                if(abs(1.-Pik) > P_cond_max_loc_omp){
                    P_cond_max_loc_omp = abs(1.-Pik);
                }
                if((1.-Pik) < P_cond_min_loc_omp){
                    P_cond_min_loc_omp = 1.-Pik;
                }
            } else {
                if (abs(Pik) > P_cond_max_loc_omp){
                    P_cond_max_loc_omp = abs(Pik);
                }
                if(Pik < P_cond_min_loc_omp){
                    P_cond_min_loc_omp = Pik;
                }
            }
                
            // if(Pik < -pow(10, -15) or Pik > (1 + pow(10, -15)) ){
            //  //   cout << "P diag= "<< Pik <<" band" << ic << " ik="<< ik << endl;
            // }
        }
    }
    
    #pragma omp master
    {
        Coulomb_set.N_shared.fill(0);
        Coulomb_set.N_exciton.fill(0);
    }


    #pragma omp barrier

    #pragma omp critical
    { // summation over threads and finding minimum over OpenMP threads
        if (P_cond_max_loc_mpi < P_cond_max_loc_omp){
            P_cond_max_loc_mpi = P_cond_max_loc_omp;
        }
        if (Coulomb_set.P_min > P_cond_min_loc_omp){
            Coulomb_set.P_min = P_cond_min_loc_omp;
        }
        
        // Coulomb_set.N_exciton[Ncv-1] += P_k_offdia[Ncv-1]/nktot;
        // if (Ncv>2) Coulomb_set.N_exciton[Ncv-2] += P_k_offdia[Ncv-2]/nktot;
        // if (Ncv>2) Coulomb_set.N_exciton[Ncv-3] += P_k_offdia[Ncv-2]/nktot;
        for (int ic=0; ic < Ncv; ic++){
            Coulomb_set.N_shared[ic] += P_k_dia[ic]/nktot;
        }
    }


    #pragma omp barrier
    
    #pragma omp master
    {
        vec1d P_k_dia_global(Ncv);
        vec1d P_k_offdia_global(Ncv);

        double P_min_MPI;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(&Coulomb_set.P_min, &P_min_MPI, 1, MPI_DOUBLE, MPI_MIN, root_rank,
               MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if ((P_min_MPI < -1e-13) and (rank_ == root_rank))
        {
            cout << "Negative population, P_min = " << P_min_MPI << endl;
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(&Coulomb_set.N_shared[0], &P_k_dia_global[0], Ncv, MPI_DOUBLE, MPI_SUM, root_rank,
               MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        // MPI_Reduce(&Coulomb_set.N_exciton[0], &P_k_offdia_global[0], Ncv, MPI_DOUBLE, MPI_SUM, root_rank,
               // MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        if(rank_==root_rank)
        {
            double P_trace = -Ncc;
            
            fp_Loss << setprecision(16) << time*time_au_fs << "   ";
            for(int ic=0; ic<Ncv; ic++){
                P_trace += P_k_dia_global[ic];
                // Coulomb_set.exciton_file << setprecision(16) << P_k_offdia_global[ic]<< "   ";
                if (ic < Ncc) P_k_dia_global[ic] = 1.0 - P_k_dia_global[ic]; 

                fp_Loss  << setprecision(16) <<P_k_dia_global[ic] << "   ";
            } 
            fp_Loss << endl;
            Coulomb_set.exciton_file << endl;
            if (abs(P_trace) > 1e-13){ 
                cout << "Preservation of particle num fail P_trace = " << P_trace << endl;
            }
            

            n_cond = 0; // population of conduction we will call the maxim population on band
            
            for(int ic=0; ic<Ncv; ic++){
                if (P_k_dia_global[ic] > n_cond){
                    n_cond = P_k_dia_global[ic];
                }      
            } 
        }


        
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(& n_cond, 1, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
        Coulomb_set.n_cond = n_cond;
    }

    #pragma omp barrier

    
    // printing diagonal Ek
    if (print_EK and Coulomb_set.Coulomb_calc){
        #pragma omp master
        {   
            MPI_Barrier(MPI_COMM_WORLD);
            int Nk_node = Coulomb_set.Ekdia.n1();
            vec2d * Ek_print = & Coulomb_set.Ekdia; // pointer to object
            vec2d Ek_receive; // variable, where we store what we will receive

            for (int i_proc = 0; i_proc < num_procs; i_proc++){
                
                if (i_proc > 0){
                    // sending array to 0 node and then print

                    // at first process 0 finds out how much data it will receive
                    if (rank_ == i_proc){ // send
                        // destination - 0
                        MPI_Send(&Nk_node, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                        // cout << "0" << endl;
                    }
                    if (rank_ == 0){ // receive
                        // source - process number i_proc
                        // cout << "Nk_node old" << Nk_node << endl;
                        MPI_Recv(&Nk_node, 1, MPI_INT, i_proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        Ek_receive.resize(Nk_node, Ncv);
                        // cout << "1" << endl;
                        // cout << "Nk_node new" << Nk_node << endl;
                    }

                    

                    if (rank_ == i_proc){ // send array
                        // destination - 0
                        MPI_Send(& Coulomb_set.Ekdia[0][0], Nk_node*Ncv, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                        // cout << "2" << endl;
                    }
                    if (rank_ == 0){ // receive
                        // source - process number i_proc
                        MPI_Recv(& Ek_receive[0][0], Nk_node*Ncv, MPI_DOUBLE, i_proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        // cout << "3" << endl;
                        Ek_print = & Ek_receive;
                        // cout << "4" << endl;
                    }
                } // if (i_proc > 0)
                 
                if (rank_ == 0){ // only rank 0 prints!
                    for(int ik = 0; ik < Nk_node; ik++){
                        for (int ic = 0; ic < Ncv; ic++){
                            Coulomb_set.Ek_file << setprecision(20) << (*Ek_print)[ik][ic] << " " ;
                        }
                        Coulomb_set.Ek_file << endl;
                    }
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }
    }
    
    #pragma omp barrier
}







void print_k_MPI(vec2d& kpt, int Nk_node) // print k- grid
{ // warning: we are inside master omp thread 
    ofstream k_stream;
    if (rank_ == 0)
    { // only root rank writes to file
        string name_file; 
        name_file= "Output/k_grid.txt";
         
        k_stream.open(name_file.c_str());
    }



    
    vec2d * k_print = & kpt; // pointer to object
    vec2d k_receive; // variable, where we store what we will receive


    for (int i_proc = 0; i_proc < num_procs; i_proc++){

        if (i_proc > 0){
            // sending array to 0 node and then print

            // at first process 0 finds out how much data it will receive
            if (rank_ == i_proc){ // send
                // destination - 0
                MPI_Send(&Nk_node, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
            if (rank_ == 0){ // receive
                // source - process number i_proc
                MPI_Recv(&Nk_node, 1, MPI_INT, i_proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                k_receive.resize(Nk_node, 3);

            }

            

            if (rank_ == i_proc){ // send array
                // destination - 0
                MPI_Send(& kpt[0][0], Nk_node*3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
            if (rank_ == 0){ // receive
                // source - process number i_proc
                MPI_Recv(& k_receive[0][0], Nk_node*3, MPI_DOUBLE, i_proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                k_print = & k_receive;
            }
        } // if (i_proc > 0)
       
        if (rank_ == 0){ // only rank 0 prints!
            for(int ik = 0; ik < Nk_node; ik++){
                k_stream << (*k_print)[ik][0] << " " << (*k_print)[ik][1] << " " << (*k_print)[ik][2] << endl;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

    }

    if (rank_ == 0){
         k_stream.close();
    }
   

}






void PrintTransientAbsMPI(int&nk,int& nktot,ofstream& fp_J, vec3x& P, double& time, vec1d& EFx,  vec1d& dk, int& Nch,double detM,vec4x& Dk)
{
    vec1x a(3); 
    CalculateTransientAbs(nk,a,P,time,EFx,dk,Nch,detM,Dk);

    vec1x global_sum(3);
    MPI_Reduce(&a[0], &global_sum[0], 3, MPI_C_DOUBLE_COMPLEX, MPI_SUM, 0,
       MPI_COMM_WORLD);    
    
    if(rank_==0)
    {
        fp_J << setw(15) << setprecision(8) << time*time_au_fs;
        fp_J << setw(15) << setprecision(8) << EFx[0];
        fp_J << setw(15) << setprecision(8) << EFx[1];
        fp_J << setw(15) << setprecision(8) << EFx[2];
        for(int ix=0; ix<3; ix++) 
        {
            fp_J << setw(15) << setprecision(8) << global_sum[ix].real()/nktot;
            fp_J << setw(15) << setprecision(8) << global_sum[ix].imag()/nktot;
        }
        fp_J << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}


void find_index_min(vec1d& k_temp, int &ik_min){
    double min_k = 1e15;
    ik_min = 0;
    for (int ik = 0; ik < k_temp.n1(); ik++){ // all wave vectors
        // yes it's quite stupid way to do this, but we have to do it once only for small array
        if (k_temp[ik] < min_k){
            min_k = k_temp[ik];
            ik_min = ik;
        }
    }
    k_temp[ik_min] = 1e15;

}





void Create_sorted_indices_k_MPI(vec2d& kpt, methods_Diff_Eq &Diff_Eq){

    int i_kx, i_ky, ik_min;
    int pol = 1;

    int nk_y = Diff_Eq.Nk_node_MPI / (Diff_Eq.Nk_1 * Diff_Eq.Nk_0);
    int nk_x = Diff_Eq.Nk_1 * Diff_Eq.Nk_0;

    if (Diff_Eq.Nk_1 > Diff_Eq.Nk_2){
        pol = 0; // find polarization
        nk_y = Diff_Eq.Nk_node_MPI / (Diff_Eq.Nk_1 * Diff_Eq.Nk_2);
        nk_x = Diff_Eq.Nk_1 * Diff_Eq.Nk_2;
    }
    // if (rank_ == 0) cout << "nk_x= " << nk_x << " nk_y= " << nk_y << endl;
    Diff_Eq.k_index.resize(nk_y);
    
    vec1d k_temp(nk_y);
    for (int ik = 0; ik < Diff_Eq.Nk_node_MPI; ik++){ // all wave vectors

        i_kx = ik / nk_x;
        i_ky = ik % nk_x;
        if (rank_ == 0){
            // cout << "ik= " << ik << "      i_kx =" << i_kx << " kpt[ik][pol]=" << kpt[ik][pol]<< endl;
        }
        if (i_ky == 0){
            k_temp[i_kx] = kpt[ik][pol];
        }
        
    }

    for (int ik = 0; ik < nk_y; ik++){ // all wave vectors
        if (rank_ == 0){
            // cout << "ik= " << ik << "      k_temp[ik] =" << k_temp[ik] << endl;
        }
        find_index_min(k_temp, ik_min);
        Diff_Eq.k_index[ik] = ik_min;
    }


}








