// time loop with dynamic dt


double epsStepAbs = Diff_Eq.epsStepAbs; // maximum accuracy of each step 
double epsStep = 1;
double epsA = pow(10, -16); // need it to avoid dividing by zero
double epsE1 = pow(10, -9);
double epsE2 = pow(10, -9);

// vec1x P_cond_prev(OMP_private.lenght_k);



double dP_ik_ic_jc;



#pragma omp master
{
    it=iti - 10;
    time_loop = -0.21 / time_au_fs;// it*dt;
    n_cond = epsA; 
}



#pragma omp barrier // sinchronise threads
while(time_loop < t_fin){


    #pragma omp barrier

    if ((it-iti) % it_resolution == 0)
    {   
        #include "if_resolution.cpp" // printig all we need

        #pragma omp master
        {
            P_cond_max = P_cond_max_loc_mpi;

            // cout <<  " P_cond_max_loc_mpi " << P_cond_max_loc_mpi << endl;
            
            MPI_Reduce(&P_cond_max_loc_mpi, &P_cond_max, 1, MPI_DOUBLE, MPI_MAX, root_rank, MPI_COMM_WORLD);    
            MPI_Bcast(& P_cond_max, 1, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
            cout <<  " P_cond_max " << P_cond_max << endl;
            cout <<  " n_cond " << n_cond << endl;
        }
    
    } // end if it_resolution



    if(Diff_Eq.Taylor){
        #include "Taylor_DE_solver_deriv_Polynom.cpp"

    } else{
        // RK time evolution
        #include "RungeKutta.cpp" 
    }





        // to estimate the step we calculate value of Pv*dt
    // such step shouldn't change any point more then a epsStepAbs
    double loc_max_dP = 0.0;
    #pragma omp master
    {
    max_dP = 0.0; 
    }
    
    for (int ik = OMP_private.begin_count; ik < OMP_private.end_count; ik++){ // all local wave vectors

        dP_ik_ic_jc = \
        3* abs(OMP_private.P_W_prev[ik - OMP_private.begin_count][Ncv-1][Ncv-1] - P0[ik][Ncv-1][Ncv-1]  )/(P_cond_max + n_cond + abs(P0[ik][Ncv-1][Ncv-1]));  // !!!!! avoid dividing by zero!
        
        if (dP_ik_ic_jc > loc_max_dP){
            loc_max_dP = dP_ik_ic_jc;
        }

    }


    #pragma omp barrier // sinchronise threads
    #pragma omp critical
    { // find maximum between different threads
        if(max_dP < loc_max_dP){
            max_dP = loc_max_dP;
        }
    }
    #pragma omp barrier // sinchronise threads

    // find maximum between nodes
    #pragma omp master
    {   
        MPI_Barrier(MPI_COMM_WORLD);
        double reduction_result;

        MPI_Reduce(&max_dP, &reduction_result, 1, MPI_DOUBLE, MPI_MAX, root_rank, MPI_COMM_WORLD);
        if (rank_ == root_rank){
            max_dP = reduction_result;

            if((max_dP >  1.2*epsStepAbs) and (rank_ == root_rank)){
                cout  << "!!!!!!!! dt too large "  << " max_dP = " << max_dP  << " dtDynamic: " << dt  << " a.u. or " << dt * time_au_fs*1000 << " as" << endl;

            }

            if ((max_dP >  epsStepAbs) and (dt > 0.001)){
                    dt *= epsStepAbs/ max_dP;
                    // cout  << "! New dt_scaled: " << dt << " epsStepAbs = " << epsStep <<"  max_dP = " << max_dP <<  " time fs: " << time_loop *time_au_fs << endl;

            } else {

                // if ((it % (it_resolution) == 0) and time_loop > 0 and dt < dt0 ){
                if ((dt < dt0) ){
                    dt *= 1.02;
                }
            }
        } // if (rank_ == root_rank)
        MPI_Bcast(& dt, 1, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
        if(Diff_Eq.Taylor){
            for (int idt = (TaylorOrder-1); idt > 0; idt--){
                dt_prev[idt] = dt_prev[idt-1];
            }
            dt_prev[0] = dt;
        } else{
            // RK time evolution
            dt6=dt/6.0;  // define RK steps in dynamic
            dt3=dt/3.0;  
            dt2=dt/2.0;
        }




        
    } //master


    
    #pragma omp barrier // sinchronise threads
    #pragma omp master
    {
        if ((rank_ == root_rank) and it % (it_resolution) == 0){
            cout  << " dtDynamic: " << dt  << " a.u. or ";
            cout  << dt * time_au_fs*1000 << " as  |||";
            cout  << " max_dP = " << max_dP << " epsStep =" << epsStepAbs  <<  " time_loop fs: " << time_loop *time_au_fs <<   endl;
        }


        time_loop += dt;

        it += 1;

    }
    #pragma omp barrier // sinchronise threads

    

}//end time evolution