// time loop with constant dt
#pragma omp master
{
    it= 0;//iti - 10;
    time_loop = -0.21 / time_au_fs;// it*dt;
    itfi = (t_fin - time_loop)/dt;
}
#pragma omp barrier // sinchronise threads

for(int it_pr=0; it_pr<=itfi; it_pr++)
{   
    #pragma omp master
    { // shared variable, only master change it
        it = it_pr;
    }

    #pragma omp barrier

    if (it % it_resolution == 0){   
        #include "if_resolution.cpp" // printig all we need
    } // end if it_resolution
    
    if(Diff_Eq.Taylor){
        #include "Taylor_DE_solver_deriv_Polynom.cpp"
    } else if(Diff_Eq.Adams_Bashforth){
        #include "Linear_Multistep_ODE_solvers.cpp"
    } else{
        // RK time evolution
        #include "RungeKutta.cpp" 
    }

    #pragma omp barrier // sinchronise threads
    #pragma omp master
    {
        if (it % (10*it_resolution) == 0){
            P_cond_max = P_cond_max_loc_mpi;
            MPI_Reduce(&P_cond_max_loc_mpi, &P_cond_max, 1, MPI_DOUBLE, MPI_MAX, root_rank, MPI_COMM_WORLD);    
            MPI_Bcast(& P_cond_max, 1, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
            if (rank_ == root_rank){
                cout <<  " P_cond_max " << P_cond_max << endl;
                cout <<  " n_cond " << n_cond << endl;
                if (OMP_private.Pk_min < -1e-13)  cout << "Negative population, P_min = " << OMP_private.Pk_min << endl;
            }
        }


        time_loop += dt;

        it += 1;

    }
    #pragma omp barrier // sinchronise threads
    

}//end time evolution