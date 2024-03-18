void CalculateWeigths(vec2d& kpt, vector<double>& Weigths, vector<vector<vector<double>>>& Bvector, vec1i& Nk)
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
		mat A(21,Bvector.size()), S(21,Bvector.size()), q(21,1), U, V, w;
		vec s;
		
		A.zeros();		S.zeros();

		for(int ishell=0; ishell<=Shell; ishell++)
		{
			int row =0; 
			for(int idir=0; idir<3; idir++)
				if(idir!= plane_direction)
				{
					for(size_t i=0; i<Bvector[ishell].size(); i++)
						A(row,ishell) += Bvector[ishell][i][idir]*Bvector[ishell][i][idir];
					row++;
				}				
			for(int idir=0; idir<3; idir++)
				for(int jdir=idir+1; jdir<3; jdir++)
					if(idir!= plane_direction && jdir != plane_direction)
					{	
						for(size_t i=0; i<Bvector[ishell].size(); i++)
							A(row,ishell) += Bvector[ishell][i][idir]*Bvector[ishell][i][jdir];
						row++;
					}
			row = 5;
                        for(int idir=0; idir<3; idir++)
                                for(int jdir=idir; jdir<3; jdir++)
					for(int ldir=jdir; ldir<3; ldir++)
						for(int mdir=ldir; mdir<3; mdir++)
						{
 							for(size_t i=0; i<Bvector[ishell].size(); i++)											A(row,ishell) += Bvector[ishell][i][idir]*Bvector[ishell][i][jdir]*Bvector[ishell][i][ldir]*Bvector[ishell][i][mdir];

							row ++;
						}
		}					

		q(0,0) = 1;      q(1,0) = 1;      q(2,0) = !(plane_direction+1);
		for (int row=3; row<21; row++)
			q(row,0) = 0.;    
		
		svd( U, s, V, A );
		for(int index=0; index<s.size(); index++)
			if(s(index)!=0) S(index,index)=s(index);
		mat Sinv(Bvector.size(),21);           Sinv.zeros();
		bool LinearDependent = false;
		for(int index=0; index<s.size(); index++)
			if(abs(s(index)) < 1.e-11) LinearDependent = true;
			else Sinv(index,index)=1./s(index);
		//w = V*Sinv*U.t()*q;
	        w=pinv(A,1.e-16)*q;
		Weigths.resize(w.size());
		for(int i=0; i<w.size(); i++)
			Weigths[i] = w(i,0);
		mat q1 = A*w;
		if (rank_ == 0){
			printf("\n*         -------------------------------------------------------------------------         *");
        	}
        	for(int ib=0; ib<Bvector[Shell].size(); ib++)
        		if (rank_ == 0){
				printf(  "\n*         |  %3i  |   % 8.5f    % 8.5f    % 8.5f   | % 8.5f  |  %5.2e  |    %3s     |         *", Shell, 
        				Bvector[Shell][ib][0],Bvector[Shell][ib][1], Bvector[Shell][ib][2], 
        				sqrt(pow(Bvector[Shell][ib][0],2)+pow(Bvector[Shell][ib][1],2)+pow(Bvector[Shell][ib][2],2)),norm(q1-q), LinearDependent?"YES":"NO");
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
		ShellsAreEnough = ( norm(q1-q)<1.e-07 );

		//cout << "ShellsAreEnough? " << (int) ShellsAreEnough<< endl;
		if(ShellCounter==32 and (rank_ == 0)) 
		{
			printf("More than 32 shells needed. The program could be very slow!\n");
			exit(1);
		}
	}
	if(rank_ == 0){
		printf("\n*         -------------------------------------------------------------------------         *\n");
	}

}//end of function CalculateWeigths





//nonmpi version
void CalculateIndices(vec1d& dk, vec2d& kpt_shared, vector<vector<vector<int>>>& GradientIndex, vector<vector<vector<double>>>& Bvector)
{
	int n_kpt = kpt_shared.n1();
	GradientIndex.resize(n_kpt);
	
	#pragma omp parallel
	{ 
		Coord_B kpb;
		vec2d kpt(n_kpt, kpt_shared.n2());
		
		// local copies
		#pragma omp critical
		{
			for(int ik=0; ik<n_kpt; ik++){
				kpt[ik][0] = kpt_shared[ik][0];
				kpt[ik][1] = kpt_shared[ik][1];
				kpt[ik][2] = kpt_shared[ik][2];
			}
		}

		#pragma omp  for schedule(dynamic)
		for(int ik=0; ik<n_kpt; ik++)
		{
			GradientIndex[ik].resize(Bvector.size());
			for(int ishell=0; ishell<Bvector.size(); ishell++)
			{
				GradientIndex[ik][ishell].resize(Bvector[ishell].size());
				for(int ib=0; ib<Bvector[ishell].size(); ib++)
				{
					int ikpb;
					Coord_B kk;
					kk.setcrys(kpt[ik][0], kpt[ik][1], kpt[ik][2]);
					kpb.setcart(kk.cart[0]+Bvector[ishell][ib][0],
								kk.cart[1]+Bvector[ishell][ib][1],
								kk.cart[2]+Bvector[ishell][ib][2]);	
			
					for(int icoor=0; icoor<3; icoor++)				
					{
						//while(kpb.crys[icoor] < -0.5-1.e-07) kpb.crys[icoor] += 1.;
						//while(kpb.crys[icoor] >= 0.5-1.e-07) kpb.crys[icoor] -= 1.;
						//double aus = kpb.crys[icoor];
						//cout << setw(20) << setprecision(10) << scientific << kpb.crys[icoor];
						kpb.crys[icoor] = fmod(kpb.crys[icoor]+0.5,1.) -0.5;				
						//cout << setw(20) << setprecision(10) << scientific << kpb.crys[icoor];
						//cout << setw(20) << setprecision(10) << scientific << kpb.crys[icoor]-aus;
					}
			
					bool found = false;
					for(int ik1 =0; ik1 < n_kpt; ik1++)
					{
						if(abs(kpb.crys[0]-kpt[ik1][0]) < 1.e-08 || abs(kpb.crys[0]-kpt[ik1][0]-1.) < 1.e-08 ||  abs(kpb.crys[0]-kpt[ik1][0]+1.) < 1.e-08)
							if(abs(kpb.crys[1]-kpt[ik1][1]) < 1.e-08 || abs(kpb.crys[1]-kpt[ik1][1]-1.) < 1.e-08 ||  abs(kpb.crys[1]-kpt[ik1][1]+1.) < 1.e-08)
								if(abs(kpb.crys[2]-kpt[ik1][2]) < 1.e-08 || abs(kpb.crys[2]-kpt[ik1][2]-1.) < 1.e-08 ||  abs(kpb.crys[2]-kpt[ik1][2]+1.) < 1.e-08)
								{
										ikpb = ik1;
										found = true;
										break;	
								}
					}
					GradientIndex[ik][ishell][ib] = ikpb;/*
					cout << (found ? "Found" : "Not found") << "    ";
					if(found)
					{

						cout << setw(20) << setprecision(10) << scientific << kpt[ikpb][0];
						cout << setw(20) << setprecision(10) << scientific << kpt[ikpb][1];
						cout << setw(20) << setprecision(10) << scientific << kpt[ikpb][2];

					}
					cout <<endl;*/
					if(!found) 
					{
						printf("Take care: I couldn't find k + b = (%5.3e,%5.3e,%5.3e) relative to k = (%5.3e,%5.3e,%5.3e), ", 
						kpb.crys[0], kpb.crys[1],kpb.crys[2], kk.crys[0], kk.crys[1],
						kk.crys[2]);
						printf("b=(%5.3f,%5.3f,%5.3f). (ik=%4i, ib=%4i)\n", Bvector[ishell][ib][0], Bvector[ishell][ib][1], Bvector[ishell][ib][2], ik, ib);
					}
				}
			}
		}
	}
}
    

