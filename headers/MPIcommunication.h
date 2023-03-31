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


