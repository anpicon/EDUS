/* its not the most elegant solution, but ill have separate functions since i cannot make the printing for a general matrix. i dont think you can get the 
shape of an arbitrary array ¯\_(ツ)_/¯ */

bool folderExists(const std::string& folderPath) {
    // Declare a struct to store file status information
    struct stat info;

    // Use stat function to retrieve file status
    // Returns 0 on success, -1 on failure
    // S_ISDIR is a macro that checks if the file is a directory
    return stat(folderPath.c_str(), &info) == 0 && S_ISDIR(info.st_mode);
}

/*

    Parameters: The function takes a single parameter, const std::string& folderPath, which represents the path to the folder/directory you want to check.

    struct stat info;: This line declares a structure named info of type struct stat, which is used to store file status information.

    stat(folderPath.c_str(), &info) == 0: This line uses the stat function to obtain information about the file or directory specified by folderPath. The function returns 0 on success and -1 on failure. If it succeeds, the file status information is stored in the info structure.

    S_ISDIR(info.st_mode): This line uses the S_ISDIR macro to check if the file specified by folderPath is a directory. The st_mode field of the info structure contains information about the file type and permissions.

    Return Statement: The function returns true if the stat function is successful (returns 0) and if the file is a directory (as determined by S_ISDIR(info.st_mode)). It returns false otherwise.

In summary, the folderExists function checks if a folder exists by using the stat function to retrieve information about the specified path. It then verifies if the file is a directory using the S_ISDIR macro. If both conditions are met, the function returns true; otherwise, it returns false.

*/

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int mpi_largest_nk(int nk){ // for some reason i dont have to define num_procs. even more curious, idk where it is defined ¯\_(ツ)_/¯

	// just to be safe i stop every thing here
	MPI_Barrier(MPI_COMM_WORLD);

	/* because not all the mpi ranks have the same number of k-points, in order to use the mpi_barriers without getting stuck anywhere, i have to know what
	the largest number of k-points in a rank is */

	int largest_nk = 0;

	// it sends to rank 0 the memory address in the current rank of the number of kpts (&nk) with tag 10. description of function's arguments below
	// MPI_Send(&data_to_send, length of the data, MPI_TYPE_OF_DATA, destination rank, tag, MPI_COMM_WORLD);
	MPI_Send(&nk, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);

	/* in rank 0 i receive the message sent just above from the other ranks (mpi uses a combination of tag, source rank and destination rank to identify the
	messages) and save its number of k-points if its larger than the overall largest number. description of the function used to receive the data is below */
	// MPI_Recv(&data_to_receive, length of the data, MPI_TYPE_OF_DATA, source rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	if (rank_ == 0){
		// there are num_procs ranks labeled by rank_
	    for (int i = 0; i < num_procs; i++){
	    	int first_received_nk;
	    	// receives the message from the other ranks
	        MPI_Recv(&first_received_nk, 1, MPI_INT, i, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	        // if this ranks number of k-points is larger than the largest number so far it becomes the overall largest number
	        if (first_received_nk > largest_nk) largest_nk = first_received_nk;
	    }

	    for (int i = 0; i < num_procs; i++){
	    	// sends the largest number to the ith rank
	        MPI_Send(&largest_nk, 1, MPI_INT, i, 11, MPI_COMM_WORLD);
	    }
	}

	// each rank receives the largest number of k-points
	MPI_Recv(&largest_nk, 1, MPI_INT, 0, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	// mpi_barrier again just to be safe
	MPI_Barrier(MPI_COMM_WORLD);

	return largest_nk;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void orderly_printing_kpoints(std::string file_name, int rank_, int nk, vec2d& kpt){

	int largest_nk = mpi_largest_nk(nk);

	// here i send the information in all ranks to rank 0. once in rank 0, i print the information by order of rank
	ofstream kpts("Output/" + file_name,ios::app);

	// now that we know the largest number of kpts in each rank, largest_nk, we make a for loop for every rank from 0 to largest_nk
	for (int ik=0; ik<largest_nk; ik++){

		// if is not in rank 0, send to rank 0
		if (rank_ != 0){
			// maybe i could move this and where it is received to outside the loop in ik
            MPI_Send(&nk, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

            if (ik < nk){
            	// kpts			
				MPI_Send(&kpt[ik][0], 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
                MPI_Send(&kpt[ik][1], 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);

                // kpt index
                MPI_Send(&ik, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);

                // rank_
                MPI_Send(&rank_, 1, MPI_INT, 0, 4, MPI_COMM_WORLD);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // in rank 0 we print the information from rank 0 and the the ones that arrive from the other ranks
        if (rank_ == 0){
            if (kpts.is_open()){
            	// if we still have points that we didnt print in rank 0, we print them here 
                if (ik < nk){
                    kpts << kpt[ik][0] << " " << kpt[ik][1] << " " << ik << " " << rank_ << "\n";
                }
                // here we receive the messages from all the other ranks from 1 to num_procs
                for (int i = 1; i < num_procs; i++){

                    int received_nk;
                    MPI_Recv(&received_nk, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    // if there are still points to be written to the file in rank i, we write them here
                    if (ik < received_nk){
                    	//kpts
                        double received_kx, received_ky;
                        MPI_Recv(&received_kx, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&received_ky, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                        // kpt index
                        int received_ik;
                        MPI_Recv(&received_ik, 1, MPI_INT, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                        // rank_
                        int received_rank_;
                        MPI_Recv(&received_rank_, 1, MPI_INT, i, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                        kpts << received_kx << " " << received_ky << " " << received_ik << " " << received_rank_ << "\n";

                    }
                }
            }
        }

 		// i think i dont need this barrier
        MPI_Barrier(MPI_COMM_WORLD);
	}

kpts.close();

}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// now i do the same thing for the energies

void orderly_printing_energies(std::string file_name, int rank_, int nk, int Ncv, vec3x& Hamiltonian, vec3x& Unitary){

	int largest_nk = mpi_largest_nk(nk);

	/* for some reason, when i first tried this, i couldnt create an array of dimension Ncv x Ncv x nk. maybe there is a smart way to do this. for now,
	i have to do the basis transformation */

	ofstream energies("Output/" + file_name, ios::app);
	for(int ik=0; ik<largest_nk; ik++){

	    complex<double> hamiltonian_in_eigen_basis[Ncv][Ncv];
	    for (int ic = 0; ic < Ncv; ic ++){
	        if (ik < nk){
	            for (int jc = 0; jc < Ncv; jc++){
	                for (int ii = 0; ii < Ncv; ii++){
	                    for (int ij = 0; ij < Ncv; ij++){
	                        hamiltonian_in_eigen_basis[ic][jc] += Unitary[ik][ic][ii]*Hamiltonian[ik][ii][ij]*conj(Unitary[ik][jc][ij]);
	                    }
	                }
	            }
	        }

	        MPI_Barrier(MPI_COMM_WORLD);

	        if (rank_ != 0){

	            // number of kpts in this rank
	            MPI_Send(&nk, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

	            if (ik < nk){

	                // band number
	                MPI_Send(&ic, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);

	                // energy complex<double>
	                double sent_energy_real = hamiltonian_in_eigen_basis[ic][ic].real();
	                double sent_energy_imag = hamiltonian_in_eigen_basis[ic][ic].imag();
	                MPI_Send(&sent_energy_real, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
	                MPI_Send(&sent_energy_imag, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);

	                // kpt index
	                MPI_Send(&ik, 1, MPI_INT, 0, 4, MPI_COMM_WORLD);

	                // rank_
	                MPI_Send(&rank_, 1, MPI_INT, 0, 5, MPI_COMM_WORLD);

	            }
	        }

	        MPI_Barrier(MPI_COMM_WORLD);

	        if (rank_ == 0){
	            if (energies.is_open()){
	                if (ik < nk){
	                    energies << ic << " " << hamiltonian_in_eigen_basis[ic][ic].real() << " " << hamiltonian_in_eigen_basis[ic][ic].imag() << " " << ik << " " << rank_ << "\n";
	                }
	                for (int i = 1; i < num_procs; i++){

	                    int received_nk;
	                    MPI_Recv(&received_nk, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	                    if (ik < received_nk){

	                        int received_ic;
	                        MPI_Recv(&received_ic, 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	                        double received_energy_real, received_energy_imag;
	                        MPI_Recv(&received_energy_real, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	                        MPI_Recv(&received_energy_imag, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	                        int received_ik;
	                        MPI_Recv(&received_ik, 1, MPI_INT, i, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	                        int received_rank_;
	                        MPI_Recv(&received_rank_, 1, MPI_INT, i, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	                        energies << received_ic << " " << received_energy_real << " " << received_energy_imag << " " << received_ik << " " << received_rank_ << "\n";

	                    }
	                }
	            }
	        }
	    }

	    // i think i dont need this barrier
	    MPI_Barrier(MPI_COMM_WORLD);
	}

	energies.close();

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// finally, the population

void orderly_printing_population(std::string file_name, int rank_, int nk, int Ncv, vec3x& population_matrix, vec3x& change_of_basis_matrix){

	int largest_nk = mpi_largest_nk(nk);
	if (!folderExists("Output/population/") && rank_ == 0) system("mkdir Output/population/");

	ofstream population("Output/population/" + file_name, ios::app);
	for(int ik=0; ik<largest_nk; ik++){

		complex<double> population_in_bands[Ncv][Ncv];
	    for (int ic = 0; ic < Ncv; ic ++){
	    	if (ik < nk){
	            for (int jc = 0; jc < Ncv; jc++){
	                for (int ii = 0; ii < Ncv; ii++){
	                    for (int ij = 0; ij < Ncv; ij++){
	                    	// if the indexes look weird, see note in Initial_Population function, in InitialPopulation.h header
	                        population_in_bands[ic][jc] += change_of_basis_matrix[ik][ic][ii]*population_matrix[ik][ij][ii]*conj(change_of_basis_matrix[ik][jc][ij]);

	                    }
	                }
	            }
	        }

	        MPI_Barrier(MPI_COMM_WORLD);

	        if (rank_ != 0){

	            // number of kpts in this rank
	            MPI_Send(&nk, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

	            if (ik < nk){

	                // band number
	                MPI_Send(&ic, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);

	                // population complex<double>
	                double sent_population_real = population_in_bands[ic][ic].real();
	                double sent_population_imag = population_in_bands[ic][ic].imag();
	                MPI_Send(&sent_population_real, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
	                MPI_Send(&sent_population_imag, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);

	                // kpt index
	                MPI_Send(&ik, 1, MPI_INT, 0, 4, MPI_COMM_WORLD);

	                // rank_
	                MPI_Send(&rank_, 1, MPI_INT, 0, 5, MPI_COMM_WORLD);

	            }
	        }

	        MPI_Barrier(MPI_COMM_WORLD);

	        if (rank_ == 0){
	            if (population.is_open()){
	                if (ik < nk){
	                    population << ic << " " << population_in_bands[ic][ic].real() << " " << population_in_bands[ic][ic].imag() << " " << ik << " " << rank_ << "\n";
	                }
	                for (int i = 1; i < num_procs; i++){

	                    int received_nk;
	                    MPI_Recv(&received_nk, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	                    if (ik < received_nk){

	                        int received_ic;
	                        MPI_Recv(&received_ic, 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	                        double received_population_real, received_population_imag;
	                        MPI_Recv(&received_population_real, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	                        MPI_Recv(&received_population_imag, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	                        int received_ik;
	                        MPI_Recv(&received_ik, 1, MPI_INT, i, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	                        int received_rank_;
	                        MPI_Recv(&received_rank_, 1, MPI_INT, i, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	                        population << received_ic << " " << received_population_real << " " << received_population_imag << " " << received_ik << " " << received_rank_ << "\n";

	                    }
	                }
	            }
	        }
	    }

	    // i think i dont need this barrier
	    MPI_Barrier(MPI_COMM_WORLD);
	}

	population.close();

}