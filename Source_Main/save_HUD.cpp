cout << "Start saving HUD" << endl;
#pragma omp parallel // num_threads(3)
{
    complex<double> H_aver;
    vec2d S(Ncv,Ncv);
    S.fill(0.0);
    if (change_gap_constant > 1e-16){ // Scissor transformation - change the gap manually
        for (int ic = 0; ic < Ncv; ic++){
            for (int jc = 0; jc < Ncv; jc++){
                if (ic == jc && ic >= Nb[1]){
                    //S[ic][jc] = change_gap_constant * energy_au_eV;
                    S[ic][jc] = change_gap_constant;
                }
            }
        }
    }

    int ik;
    #pragma omp for schedule(dynamic)
    for(ik=0; ik<nk; ik++){ // Here we create Hamiltonian, Unitary matrix and Dipoles
        Coord_B ktemp;
        vec2x Hk(Ncv,Ncv);
        vec2x Uk(Ncv,Ncv);
        vec3x Dk(Ncv,Ncv,3); //vec2x Dk(Ncv,Ncv);
        ktemp.setcrys(kpt[ik][0], kpt[ik][1], kpt[ik][2]);

        gHUD(Hk, Uk, Dk, ktemp, iMode); // get Hamiltonian, U matrix and dipoles for particular k point ktemp 
                                        // iMode is a model we choose in the input

        if(ik==nk/4 || ik==nk/2 || ik==3*nk/4) cout << " Saving " << 100*ik / nk << "%" << endl;

        for(int ic=0; ic<Ncv; ic++){
            for(int jc=0; jc<Ncv; jc++){
                complex<double> value = 0.;
                if (change_gap_constant > 1e-16){// change the gap manually
                    for (int ii = 0; ii < Ncv; ii++){
                        value += conj(Uk[ii][ic])*S[ii][ii]*Uk[ii][jc]; // add scissor matrix in wannier representation to wannier hamiltonian
                    }
                }

                Hamiltonian[ik][ic][jc] = Hk[ic][jc] + value;
                Unitary[ik][ic][jc] = Uk[ic][jc];
                Dipole[ik][ic][jc][0] = Dk[ic][jc][0];
                Dipole[ik][ic][jc][1] = Dk[ic][jc][1];
                Dipole[ik][ic][jc][2] = Dk[ic][jc][2];

                if(iMode=="W"){ // make dipoles hermitian
                    Dipole[ik][ic][jc][0] = 0.5*(Dk[ic][jc][0]+Dk[jc][ic][0]);
                    Dipole[ik][ic][jc][1] = 0.5*(Dk[ic][jc][1]+Dk[jc][ic][1]);
                    Dipole[ik][ic][jc][2] = 0.5*(Dk[ic][jc][2]+Dk[jc][ic][2]);
                }
            }
        }

        // Check if Hamiltonian is hermitian
        for (int ic=0; ic<Ncv; ic++){// summation over bands
            for (int jc=ic+1; jc<Ncv; jc++){ // summation over bands
                double Hermit_check;
                Hermit_check = abs(Hk[ic][jc] - conj(Hk[jc][ic]));
                if (Hermit_check > 1e-15){
                    cout << "Non-Hermitian Hamiltonian " << " ik" << ik << " ic=" << ic<< " jc=" << jc << " |H-H*|=" << Hermit_check << endl;
                    Hk[ic][jc] = 0.5* (Hk[ic][jc] + conj(Hk[jc][ic]) ); // symmetrization
                    Hk[jc][ic] = conj(Hk[ic][jc]);
                }
            }
        }

    }
    #pragma omp barrier
}

// I WROTE THIS vv //
MPI_Barrier(MPI_COMM_WORLD);

// if (rank_ == 0) cout << "PRINTING K-POINTS\n";
// #pragma omp barrier
// #pragma omp master
// {
//     orderly_printing_kpoints("kpts.txt", rank_, nk, kpt);
// }

// if (rank_ == 0) cout << "PRINTING ENERGIES\n";
// #pragma omp barrier
// #pragma omp master
// {
//     orderly_printing_energies("energies.txt", rank_, nk, Ncv, Hamiltonian, Unitary);
// }

// I WROTE THIS ^^ //

printf("end of saving matrices. \n");

Mpi_communication(Hamiltonian,  message);
MPI_Barrier(MPI_COMM_WORLD);
    
if(iCurrent)
{ // 
    //
    printf("Saving gradient of energy...\n");
    GradientEnergy.resize(nk,Ncv,Ncv,3); 
    GradientEnergy.fill(0);
    
    // should work in parallel but really no point to do it
    // #pragma omp parallel for schedule(dynamic)
    for(int ik=0; ik<nk; ik++){
        for(int i=0; i<3; i++)
            for(int ic=0; ic<Ncv; ic++)
                for(int jc=0; jc<Ncv; jc++)
                    for(int ishell=0; ishell<Bvector.size(); ishell++)
                        for(int ib=0; ib<Bvector[ishell].size(); ib++){
                            GradientEnergy[ik][ic][jc][i] += Weigths[ishell]*Bvector[ishell][ib][i]*
                            (Hamiltonian[GradientIndex[ik][ishell][ib]][ic][jc]-Hamiltonian[ik][ic][jc]);
                        }
    }
                            
    printf("done\n");

    // string name_file;     
    // stringstream sname;
    // sname.seekp(0,ios::beg);   
    // name_file= "Output/GradientEnergy.txt";
    // cout<< name_file << endl; 
    // Coulomb_set.P_stream.open(name_file.c_str());
    // for(int ik=0; ik<nk; ik++){
    //     for(int ic=0; ic<Ncv; ic++){
    //         for(int jc=ic; jc<Ncv; jc++){
    //             for(int i=0; i<3; i++){

    //                 Coulomb_set.P_stream << setprecision(16)  << real(GradientEnergy[ik][ic][jc][i]) << " ";
    //                 Coulomb_set.P_stream << setprecision(16)  << imag(GradientEnergy[ik][ic][jc][i]) << " ";
    //             }
    //         }
    //     }
    //     Coulomb_set.P_stream  << endl;
    // }

}




