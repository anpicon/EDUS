double norm(vec2d& A)
{
	double n=0;
	for(int ic=0; ic<A.n1(); ic++)
	{
		n+=A[ic][0]*A[ic][0];
	}
	return n;
}

void CalculateWeigths_new(vec2d& kpt, vector<double>& Weigths, vector<vector<vector<double>>>& Bvector, vec1i& Nk)
{
	//is the material 2D?
	int plane_direction = -1;
    for(int i=0; i<3; i++)
		if(Nk[i]==1) plane_direction=i; 

	vector<double> min(1,1.e-08); 
	bool ShellsAreEnough = false;
	Coord_B ktest; 
        int ShellCounter = 0;
        if (rank_ == 0){
		printf("\n\n\n\n*********************************************************************************************\n");
		printf("*                                   Shell calculation                                       *\n");
		printf("*********************************************************************************************\n");
		printf("*         | Shell |                Bvectors              |   Norm    | norm(q1-q) | Dependent? |         *");
	}
	while(!ShellsAreEnough)
	{
		//1. Look for next minimum
		min.push_back(1000.);
		for(int ik=0; ik<kpt.n1(); ik++)
		{
			bool IsInShells = false;
			ktest.setcrys(kpt[ik][0], kpt[ik][1], kpt[ik][2]);
			for(int i=0; i<min.size()-1; i++)
				IsInShells = IsInShells || (ktest.norm() <= min[i]+1.e-11);

			if(!IsInShells && ktest.norm() < min.back())
				min[min.size()-1] = ktest.norm();
		}


		//2. Find all k vectors with length = minimum
		Bvector.resize(Bvector.size()+1);
		int Shell = Bvector.size()-1;
		for(int ik=0; ik<kpt.n1(); ik++)
		{ 
			ktest.setcrys(kpt[ik][0], kpt[ik][1], kpt[ik][2]);
			if( abs(ktest.norm() - min[min.size()-1]) < 1.e-11 )
			{
				vector<double> kcart(3);
				for(int i=0; i<3; i++) kcart[i] = ktest.cart[i];
				Bvector[Shell].push_back(kcart);
			}
		}

		//3. Find weigth of the shell
		// Mostofi, A. A., J. R. Yates, Y.-S. Lee, I. Souza, D. Vanderbilt, and
        // N. Marzari, 2008, Comput. Phys. Commun. 178, 685
		// section 3.2
		//mat A(21,Bvector.size()), S(21,Bvector.size()), q(21,1), U, V, w;
		vec2d A(21,Bvector.size());
		vec2d q(21,1);
		vec2d w(Bvector.size(),1);
		//vec s;
		
		//A.zeros();		S.zeros();
		A.fill(0);
		for(int ishell=0; ishell<=Shell; ishell++)
		{
			int row =0; 
			for(int idir=0; idir<3; idir++)
				if(idir!= plane_direction)
				{
					for(size_t i=0; i<Bvector[ishell].size(); i++)
						A[row][ishell]+= Bvector[ishell][i][idir]*Bvector[ishell][i][idir];
					row++;
				}				
			for(int idir=0; idir<3; idir++)
				for(int jdir=idir+1; jdir<3; jdir++)
					if(idir!= plane_direction && jdir != plane_direction)
					{	
						for(size_t i=0; i<Bvector[ishell].size(); i++)
							A[row][ishell] += Bvector[ishell][i][idir]*Bvector[ishell][i][jdir];
						row++;
					}
			row = 5;
                        for(int idir=0; idir<3; idir++)
                                for(int jdir=idir; jdir<3; jdir++)
					for(int ldir=jdir; ldir<3; ldir++)
						for(int mdir=ldir; mdir<3; mdir++)
						{
 							for(size_t i=0; i<Bvector[ishell].size(); i++)											A[row][ishell] += Bvector[ishell][i][idir]*Bvector[ishell][i][jdir]*Bvector[ishell][i][ldir]*Bvector[ishell][i][mdir];

							row ++;
						}
		}					

		q[0][0] = 1;      q[1][0] = 1;      q[2][0] = !(plane_direction+1);
		for (int row=3; row<21; row++)
			q[row][0] = 0.;    

		//LU factorization of A
		//--> to change with pseudoinverse!!     vec2d LU(A.n1(), A.n2());
		//--> to change with pseudoinverse!!     for(int ic=0; ic<A.n1(); ic++)
		//--> to change with pseudoinverse!!     {
		//--> to change with pseudoinverse!!     	for(int jc=0; jc<A.n2(); jc++)
		//--> to change with pseudoinverse!!     	{
		//--> to change with pseudoinverse!!     		LU[ic][jc]=A[ic][jc];
		//--> to change with pseudoinverse!!     	}
		//--> to change with pseudoinverse!!     }
		//--> to change with pseudoinverse!!     std::vector<int> IPIV(A.n1());
		//--> to change with pseudoinverse!!     for(int ic=0; ic<IPIV.n1(); ic++)
		//--> to change with pseudoinverse!!     {
		//--> to change with pseudoinverse!!     	IPIV[ic]=ic;
		//--> to change with pseudoinverse!!     }
		//--> to change with pseudoinverse!!     int info = LAPACKE_dgetrf (LAPACK_ROW_MAJOR, A.n1(), A.n2(), &A[0][0], A.n2(), &IPIV[0]);
		//--> to change with pseudoinverse!!     //Linear system A*w = q
		//--> to change with pseudoinverse!!     info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', A.n1(), 1, &A[0][0], &q[0][0], A.n2(), )
	    //    w=pinv(A,1.e-16)*q;
		Weigths.resize(w.n1());
		for(int i=0; i<w.n1(); i++)
			Weigths[i] = w[i][0];
		//mat q1 = A*w;
		vec2d q1(21,1);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
           21, Bvector.size(), 1, 1., &A[0][0], 1, &w[0][0], Bvector.size(), 0., &q1[0][0], Bvector.size());
		vec2d deltaq(21,1);
		for(int ic=0; ic<21; ic++){
			deltaq[ic][0] = q1[ic][0]-q[ic][0];
		}
		if (rank_ == 0){
			printf("\n*         -------------------------------------------------------------------------         *");
        	}
        	for(int ib=0; ib<Bvector[Shell].size(); ib++)
        		if (rank_ == 0){
				printf(  "\n*         |  %3i  |   % 8.5f    % 8.5f    % 8.5f   | % 8.5f  |  %5.2e  |            |         *", Shell, 
        				Bvector[Shell][ib][0],Bvector[Shell][ib][1], Bvector[Shell][ib][2], 
        				sqrt(pow(Bvector[Shell][ib][0],2)+pow(Bvector[Shell][ib][1],2)+pow(Bvector[Shell][ib][2],2)),norm(deltaq) );
			}

		//4. Check if:
		//   - The new shell is linear dependent on previous ones. If so, one of the singular values is very close to 0
		//   - A*w = q is satisfied. If it isn't, we go forward adding other shells
		//		if(LinearDependent) 
		//		{
		//			Weigths.erase(Weigths.begin(), Weigths.end());
		//			Bvector.pop_back();
		//		}
		//                else
		ShellCounter++;
		ShellsAreEnough = ( norm(deltaq)<1.e-07 );

		//cout << "ShellsAreEnough? " << (int) ShellsAreEnough<< endl;
		if(ShellCounter==8 and (rank_ == 0)) 
		{
			printf("More than 8 shells needed. The program could be very slow!\n");
			exit(1);
		}
	}
	if(rank_ == 0){
		printf("\n*         -------------------------------------------------------------------------         *\n");
	}

}//end of function CalculateWeigths




