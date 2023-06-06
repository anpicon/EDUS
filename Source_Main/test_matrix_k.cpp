vec2x H(Ncv,Ncv);
vec2x U(Ncv, Ncv);
//vec2d E(Ncv, Ncv);
vec3x DD(Ncv,Ncv,3);//vec2x dx(Ncv,Ncv); //vec2x dy(Ncv,Ncv); //vec2x dz(Ncv,Ncv);
Coord_B kkk;
ofstream Test;

if (rank_ == root_rank){
    printf("Printing test matrices for check...\n");
    
    kkk.setcrys(0,0.3,0.1);
    gHUD(H,U,DD,kkk,iMode);
    
    Test.open("Test_matrices.txt");

    Test << "k = " << kkk.crys[0] << " " << kkk.crys[1] << " " << kkk.crys[2] << endl;
    Test << "Energy: \n";
    for(int ic=0; ic<Ncv; ic++)    
    {
        for(int jc=0; jc<Ncv; jc++)
            Test << H[ic][jc] << " ";
        Test << "\n";
    }
    Test << "\n";
    Test << "Unitary: \n";
    for(int ic=0; ic<Ncv; ic++)
    {
        for(int jc=0; jc<Ncv; jc++)
            Test << U[ic][jc] << " ";
        Test << "\n";
    }
    Test << "\n";
    /*
        Test << "Energy_bloch: \n";
        for(int ic=0; ic<Ncv; ic++)
        {
            for(int jc=0; jc<Ncv; jc++)
                Test << E[ic][jc] << " ";
            Test << "\n";
        }
        Test << "\n";
     */
    Test << "U E_W U_dagger (to compare with eq. 104 in notes).\n";
    vec2x E2(Ncv,Ncv);  
    for(int ic=0; ic<Ncv; ic++)
    {
        for(int jc=0; jc<Ncv; jc++)
        {
            E2[ic][jc] = 0.;
            
            for(int mc = 0; mc < Ncv; mc++)
                for(int nc = 0; nc<Ncv; nc++)
                    E2[ic][jc] += U[ic][mc]*H[mc][nc]*conj(U[jc][nc]);

            Test << E2[ic][jc] << " ";
        }
        Test << "\n";
    }
    Test << "U  U_dagger (should be 1).\n";
    for(int ic=0; ic<Ncv; ic++) {
        for(int jc=0; jc<Ncv; jc++) {
            E2[ic][jc] = 0.;
            
            for(int mc = 0; mc < Ncv; mc++)
                E2[ic][jc] += U[ic][mc]*conj(U[jc][mc]);

            Test << E2[ic][jc] << " ";
            if (abs(E2[ic][jc])> 1e-10 and ic != jc ){
                Test << "ERROR" << " ";
            } 
            if (abs(E2[ic][jc] - 1.)> 1e-10 and ic == jc ){
                Test << "ERROR" << " ";
            }
        }
        Test << "\n";
    }

    Test << "\n";
    Test << "dx: \n";
    for(int ic=0; ic<Ncv; ic++)
    {
        for(int jc=0; jc<Ncv; jc++)
            Test << DD[ic][jc][0] << " ";
        Test << "\n";
    }
    Test << "\n";
    Test << "dy: \n";
    for(int ic=0; ic<Ncv; ic++)
    {
        for(int jc=0; jc<Ncv; jc++)
            Test << DD[ic][jc][1] << " ";
        Test << "\n";
    }
    Test << "\n";    
    Test << "dz: \n";
    for(int ic=0; ic<Ncv; ic++)
    {
        for(int jc=0; jc<Ncv; jc++)
            Test << DD[ic][jc][2] << " ";
        Test << "\n";
    }
    Test << "\n";    
    Test.close();
}

MPI_Barrier(MPI_COMM_WORLD);





/*
cout << "fermi energy" <<  endl;
Coord_B kkk;
kkk.setcrys(0.,0.3333333333,-0.3333333333);
vec2x ene(Ncv,Ncv);
for(int ic=0; ic<Ncv; ic++) for(int jc=0; jc<Ncv; jc++) ene[ic][jc] = WModel.energy(ic,jc,kkk);


cout << "0 0 " <<  ene[0][0] << endl;
cout << "0 1 " <<  ene[0][1] << endl;
cout << "1 0 " <<  ene[1][0] << endl;
cout << "1 1 " <<  ene[1][1] << endl;
*/
int my_rank;
int p;
MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
MPI_Comm_size (MPI_COMM_WORLD, &p); // number of processes

for (int i = 0; i < p; i++){
    // Cheking kpt borders
    if (my_rank == i){
        MPI_Barrier(MPI_COMM_WORLD);
        cout << " my_rank " <<  my_rank << endl;
        cout << " kpt.n1 " <<  kpt.n1() << " kpt.n2 " <<  kpt.n2() << endl;
        cout << " kpt[0][0] " << kpt[0][0]  << " kpt[0][1] " << kpt[0][1]  << " kpt[0][2] " << kpt[0][2] << endl; 
    }

     MPI_Barrier(MPI_COMM_WORLD);
}
MPI_Barrier(MPI_COMM_WORLD); 
