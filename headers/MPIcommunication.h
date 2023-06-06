void Mpi_communication(vec3x& P, vector<Message>& message)
{
	int Ncv = P.n2();
        if(num_procs > 1){ // many nodes version
    
		for(Message& m : message){
			if(m.send){
				MPI_Send(&P[m.buf][0][0], Ncv*Ncv*m.count, MPI_C_DOUBLE_COMPLEX, 
								m.partner_rank, m.tag, MPI_COMM_WORLD);
			} else {
				MPI_Recv(&P[m.buf][0][0], Ncv*Ncv*m.count, MPI_C_DOUBLE_COMPLEX, 
					m.partner_rank, m.tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		
	} 
	
}	




void MPI_reduce_vec1x(vec1x& V, int root_rank, int my_rank, int N){
// this function sums all elements of vec2x objects

	vec1x V_buff(N);
	MPI_Reduce(&V[0], &V_buff[0], 
	    N, MPI_C_DOUBLE_COMPLEX, MPI_SUM, root_rank, MPI_COMM_WORLD);

	if (my_rank == root_rank){
		V =  V_buff;
	}



}


// MUCH SLOWER THAN MPI_reduce_vec1x
void MPI_reduce_vec2x(vec2x& V, int root_rank, int my_rank, int N1, int N2 ){
// this function sums all elements of vec2x objects

	vec2x V_buff(N1, N2);
	MPI_Reduce(&V[0][0], &V_buff[0][0], 
	    N1*N2, MPI_C_DOUBLE_COMPLEX, MPI_SUM, root_rank, MPI_COMM_WORLD);

	if (my_rank == root_rank){
		for (int i1 = 0; i1 < N1; i1++){
			for (int i2 = 0; i2 < N2; i2++){
				V[i1][i2] =  V_buff[i1][i2];
			}
		}
		
	}



}



// MUCH SLOWER THAN MPI_reduce_vec1x
void MPI_reduce_vec2x_safe(vec2x& V, int root_rank, int my_rank, int N1, int N2 ){
// this function sums all elements of vec2x objects
	MPI_Barrier(MPI_COMM_WORLD);
	// vec2x V_buff(N1, N2);
	complex<double> V_buff, V_curr;
	for (int i1 = 0; i1 < N1; i1++){
		for (int i2 = 0; i2 < N2; i2++){
			V_curr = V[i1][i2];
			V_buff = V[i1][i2];
			MPI_Reduce(&V_curr, &V_buff, 
		    	1, MPI_C_DOUBLE_COMPLEX, MPI_SUM, root_rank, MPI_COMM_WORLD);

			V[i1][i2] =  V_buff;
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	// if (my_rank == root_rank){
	// 	for (int i1 = 0; i1 < N1; i1++){
	// 		for (int i2 = 0; i2 < N2; i2++){
	// 			V[i1][i2] =  V_buff[i1][i2];
	// 		}
	// 	}
		
	// }

}


// MUCH SLOWER 
void MPI_Bcast_vec2x_safe(vec2x& V, int root_rank, int my_rank, int N1, int N2 ){
// this function sums all elements of vec2x objects
	MPI_Barrier(MPI_COMM_WORLD);
	// vec2x V_buff(N1, N2);
	complex<double> V_buff;
	for (int i1 = 0; i1 < N1; i1++){
		for (int i2 = 0; i2 < N2; i2++){

			V_buff = V[i1][i2];
			MPI_Bcast(& V_buff, 1, 
				MPI_C_DOUBLE_COMPLEX, root_rank, MPI_COMM_WORLD);

			V[i1][i2] =  V_buff;
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);


}