
    printf(" \n \n \n*****************************************************\n");
    printf("*    Allocating memory for the calculation...       *\n");
    printf("*****************************************************\n");
    
    /*********************************************************************************************************************************
     ****************                   Large arrays of the system: definition                 ***************************************
     *********************************************************************************************************************************
     *                                                                                                                               *
     *    DM array                         P0[k][i][j]      ->  will be the density matrix at previous step                          *
     *    DM array                         P1[k][i][j]      ->  with this we construct the density matrix at the n-th step           *
     *    k array                          P2[k][i][j]      ->  will contain information for the variables in the RK steps           *
     *    f array                          Pv[k][i][j]      ->  f in the differential equation at each stpe of Runge-Kutta           *
     *    Dissipation term array           T[i][j]                                                                                   *
     *    k point array                                                                                                              *
     *                                                                                                                               *
     *********************************************************************************************************************************/

    printf("->    Allocated k points  \n");
    printf("                Memory required ~ %4.4f Mb\n", kpt.n1()*kpt.n2()*8./(1024.*1024.));  


    //NOTE!!! IN mpi nk != kpt.n1()!!
    // nk -> k to propagate ,
    // kpt.n1() -> k to propagate + k to receive
    MPI_Barrier(MPI_COMM_WORLD);

    vec3x Hamiltonian(kpt.n1(), Ncv,Ncv);
    vec3x Unitary(nk, Ncv,Ncv);
    vec3x H_Eigen(nk, Ncv,Ncv);
    vec4x Dipole(nk, Ncv,Ncv,3);

    MPI_Barrier(MPI_COMM_WORLD);
    printf("->    Allocated H,U,D  \n");
    printf("                Memory required ~ %4.4f Mb\n", 5*nk*Ncv*Ncv*8.*2/(1024.*1024.));  

    //note: P0 and P2 also have external k in order to calculate the gradient
    vec3x Pb(nk, Ncv, Ncv);              Pb.fill(0.);
    vec3x P0(kpt.n1(), Ncv, Ncv);        P0.fill(0.);
    vec3x P1(nk, Ncv, Ncv);              P1.fill(0.);
    vec3x P2(kpt.n1(), Ncv, Ncv);        P2.fill(0.);
    vec3x Pv(nk, Ncv, Ncv);              Pv.fill(0.);

    // shared array to print in file 
    Coulomb_set.Pk_shared.resize(nk);

    if (Coulomb_set.Coulomb_calc){
        Coulomb_set.Hk_renorm_shared.resize(kpt.n1(), Ncv, Ncv);
        Coulomb_set.P_dkx.resize(kpt.n1(), Ncv, Ncv);
        Coulomb_set.P_dky.resize(kpt.n1(), Ncv, Ncv);
        Coulomb_set.P0_dia.resize(kpt.n1(), Ncv, Ncv);
        Coulomb_set.P_d2kx.resize(kpt.n1(), Ncv, Ncv);
        Coulomb_set.P_d2ky.resize(kpt.n1(), Ncv, Ncv);
        Coulomb_set.P_dky_dkx.resize(kpt.n1(), Ncv, Ncv);

        Coulomb_set.Xk_storage.resize(nk, Ncv, Ncv);

        Coulomb_set.P_dkx.fill(0.);
        Coulomb_set.P_dky.fill(0.);
        Coulomb_set.P0_dia.fill(0.);
        Coulomb_set.P_d2kx.fill(0.);
        Coulomb_set.P_d2ky.fill(0.);
        Coulomb_set.P_dky_dkx.fill(0.);
        Coulomb_set.Xk_storage.fill(0.);



    }

    printf("->    Allocated Pb, P0, P1, P2, Pv for permanent and temporary density matrix and for Runge-Kutta steps  \n");
    printf("                Memory required ~ %4.4f Mb\n", 2*(kpt.n1()+nk)*P0.n2()*P0.n3()*8.*2./(1024.*1024.)  );



    vec2d   T(Ncv, Ncv);           T.fill(0.);       


    printf("->    Allocated decoherent tensor T  \n");
    printf("                Memory required ~ %4.4f Mb\n", T.n1()*T.n2()*8./(1024.*1024.));

    printf(" \n \n \n*****************************************************\n");
    printf(         "*     Initializing vectors for the calculation      *\n");
    printf(         "*****************************************************\n");


    Initialize_T(T, T1, T2, Tch, Nb[0]);
    printf("->    Decoherence tensor T initialized. \n");
    