//save current population in OMP_private.P0 before apply RungeKutta
if (Diff_Eq.dynamical_dt_evolution){
    for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++){ // all local wave vectors
        for(int ic=0; ic<Ncv; ic++){
            OMP_private.P0[ik_pr][ic][ic] = P0[ik_pr + OMP_private.begin_count][ic][ic];
            // OMP_private.P2_dia[ik_pr][ic][jc] = OMP_private.P_diag[ik_pr + OMP_private.begin_count][ic][jc];
        }
    }
}

//------ 1st step Runge-Kutta
#pragma omp barrier
#pragma omp master
{ // shared variables, only master change it

	Mpi_communication(P0,  message);

	MPI_Barrier(MPI_COMM_WORLD);
    time_RK = time_loop;
}
#pragma omp barrier
EF_pr[0] = pulse1.E(time_RK);
EF_pr[1] = .5*pulse2.E(time_RK);

get_derivative_Df(kpt, P0, OMP_private.Pv, T, Nb, 
    EF_pr, pulse2.wl, 
    Coulomb_set, trig_k_omp, OMP_private,
    GradientIndex, Weigths, Bvector, 
    root_rank, rank_, message);

Runge_Kutta_Ad(P0,P1,OMP_private.Pv ,dt6,OMP_private);
Runge_Kutta_Ad(P0,P2,OMP_private.Pv ,dt2,OMP_private);



//------ 2nd step Runge-Kutta
#pragma omp master
{  // shared variables, only master change it
	time_RK = time_loop + dt2;
	Mpi_communication(P2, message);

	MPI_Barrier(MPI_COMM_WORLD);
}

#pragma omp barrier
EF_pr[0] = pulse1.E(time_RK);
EF_pr[1] = .5*pulse2.E(time_RK);


get_derivative_Df(kpt, P2, OMP_private.Pv, T, Nb,
    EF_pr, pulse2.wl, 
    Coulomb_set, trig_k_omp, OMP_private,
    GradientIndex, Weigths, Bvector, 
    root_rank, rank_, message);
/*
    ofstream Ptest;
    stringstream sn; sn << cnt;
    Ptest.open("test_"+ sn.str()+srank.str()+".txt");
    for(int ik=0; ik<nk; ik++)
    {
        Ptest << kpt[ik][0] << " " << kpt[ik][1] << " " << kpt[ik][2] << " ";
        for(int ic=0; ic< P0.n2(); ic++)
            for(int jc=0; jc<P0.n2(); jc++)
                Ptest << P2[ik][ic][jc].real() << " " << P2[ik][ic][jc].imag()<<" ";
        Ptest << endl;
    }   
*/
Runge_Kutta_Ac(P1,OMP_private.Pv,dt3,OMP_private);
Runge_Kutta_Ad(P0,P2,OMP_private.Pv,dt2,OMP_private);
// Runge_Kutta_Ad(P0,P2,OMP_private.Pv,dt2,OMP_private);



//------ 3rd step Runge-Kutta

#pragma omp master
{  // shared variables, only master change it
	Mpi_communication(P2, message);
	MPI_Barrier(MPI_COMM_WORLD);
}
#pragma omp barrier

get_derivative_Df(kpt, P2, OMP_private.Pv, T, Nb,
    EF_pr, pulse2.wl, 
    Coulomb_set, trig_k_omp, OMP_private,
    GradientIndex, Weigths, Bvector, 
    root_rank, rank_, message);


Runge_Kutta_Ac(P1, OMP_private.Pv, dt3,OMP_private);
Runge_Kutta_Ad(P0,P2, OMP_private.Pv, dt,OMP_private);


#pragma omp barrier

//------ 4th step Runge-Kutta

#pragma omp master
{  // shared variables, only master change it
	time_RK = time_loop + dt;
	Mpi_communication(P2, message);
	MPI_Barrier(MPI_COMM_WORLD);
}

#pragma omp barrier
EF_pr[0] = pulse1.E(time_RK);
EF_pr[1] = .5*pulse2.E(time_RK);

get_derivative_Df(kpt, P2, OMP_private.Pv, T, Nb,
    EF_pr, pulse2.wl, 
    Coulomb_set, trig_k_omp, OMP_private,
    GradientIndex, Weigths, Bvector, 
    root_rank, rank_, message);

Runge_Kutta_Ad(P1,P0,OMP_private.Pv,dt6,OMP_private);


// store population before we did RK step in OMP_private.P_W_prev
if (Diff_Eq.dynamical_dt_evolution){
    for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++){ // all local wave vectors
        for(int ic=0; ic<Ncv; ic++){
                OMP_private.P_W_prev[ik_pr][ic][ic] = OMP_private.P0[ik_pr][ic][ic];
                // OMP_private.P_dia_prev[ik_pr][ic][jc] = OMP_private.P2_dia[ik_pr][ic][jc];
            
        }
    }
}