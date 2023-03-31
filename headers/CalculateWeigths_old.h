#include "quicksort.h"
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
    printf("\n\n\n\n*********************************************************************************************\n");
    printf("*                                   Shell calculation                                       *\n");
    printf("*********************************************************************************************\n");




	printf("*         | Shell |                Bvectors              |   Norm    | Dependent? |         *");

	while(!ShellsAreEnough)
	{
		//1. Look for next minimum
		min.push_back(1000.);
		for(int ik=0; ik<kpt.n1(); ik++)
		{
			bool IsInShells = false;
			ktest.setcrys(kpt[ik][0], kpt[ik][1], kpt[ik][2]);
			for(int i=0; i<min.size()-1; i++)
				IsInShells = IsInShells || (ktest.norm() <= min[i]+1.e-08);

			if(!IsInShells && ktest.norm() < min.back())
				min[min.size()-1] = ktest.norm();
		}

                ShellCounter++;
		//2. Find all k vectors with length = minimum
		Bvector.resize(Bvector.size()+1);
		int Shell = Bvector.size()-1;
		for(int ik=0; ik<kpt.n1(); ik++)
		{ 
			ktest.setcrys(kpt[ik][0], kpt[ik][1], kpt[ik][2]);
			if( abs(ktest.norm() - min[min.size()-1]) < 1.e-08 )
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
		mat A(6,Bvector.size()), S(6,Bvector.size()), q(6,1), U, V, w;
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
		}

		q(0,0) = 1;      q(1,0) = 1;      q(2,0) = !(plane_direction+1);
		q(3,0) = 0;      q(4,0) = 0;      q(5,0) = 0;
		
		svd( U, s, V, A );
		
		for(int index=0; index<s.size(); index++)
			if(s(index)!=0) S(index,index)=s(index);

		mat Sinv(Bvector.size(),6);           Sinv.zeros();
		bool LinearDependent = false;
		for(int index=0; index<s.size(); index++)
			if(abs(s(index)) < 1.e-08) LinearDependent = true;
			else Sinv(index,index)=1./s(index);
		w = V*Sinv*U.t()*q;
	
		Weigths.resize(w.size());
		for(int i=0; i<w.size(); i++)
			Weigths[i] = w(i,0);
		mat q1 = A*w;
		printf("\n*         -------------------------------------------------------------------------         *");
        for(int ib=0; ib<Bvector[Shell].size(); ib++)
        		printf(  "\n*         |  %3i  |   % 8.5f    % 8.5f    % 8.5f   | % 8.5f  |    %3s     |         *", Shell, 
        			Bvector[Shell][ib][0],Bvector[Shell][ib][1], Bvector[Shell][ib][2], 
        			pow(Bvector[Shell][ib][0],2)+pow(Bvector[Shell][ib][1],2)+pow(Bvector[Shell][ib][2],2),(LinearDependent?"YES":"NO"));

        //4. Check if:
        //   - The new shell is linear dependent on previous ones. If so, one of the singular values is very close to 0
        //   - A*w = q is satisfied. If it isn't, we go forward adding other shells
		if(LinearDependent) 
		{
			Weigths.erase(Weigths.begin(), Weigths.end());
			Bvector.pop_back();
		}
		
		ShellsAreEnough = (abs(q1(0,0)-q(0,0))<1.e-08);
		for(int irow=1; irow<6; irow++) ShellsAreEnough = ShellsAreEnough || (abs(q1(irow,0)-q(irow,0))<1.e-08);
		//cout << "ShellsAreEnough? " << (int) ShellsAreEnough<< endl;
		if(ShellCounter==8) 
		{
			printf("More than 8 shells needed. The program could be very slow!\n");
			exit(1);
		}
	}
	printf("\n*         -------------------------------------------------------------------------         *\n");
}//end of function CalculateWeigths




//nonmpi version
void CalculateIndices(vec2d& kpt, vector<vector<vector<int>>>& GradientIndex, vector<vector<vector<double>>>& Bvector)
{
	Coord_B kpb;
	GradientIndex.resize(kpt.n1());
	 
	for(int ik=0; ik<kpt.n1(); ik++)
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
					while(kpb.crys[icoor] < -0.5-1.e-07) kpb.crys[icoor] += 1.;
					while(kpb.crys[icoor] >= 0.5-1.e-07) kpb.crys[icoor] -= 1.;
				}
		
				bool found = false;
				for(int ik1 =0; ik1 < kpt.n1(); ik1++)
				{
					if(abs(kpb.crys[0]-kpt[ik1][0]) < 1.e-08 && 
					   abs(kpb.crys[1]-kpt[ik1][1]) < 1.e-08 && 
					   abs(kpb.crys[2]-kpt[ik1][2]) < 1.e-08)
					{
						ikpb = ik1;
						found = true;
						break;	
					}
				}
				GradientIndex[ik][ishell][ib] = ikpb;
				if(!found) 
				{
					printf("Take care: I couldn't find k + b = (%5.3f,%5.3f,%5.3f) relative to k = (%5.3f,%5.3f,%5.3f), ", 
					kpb.crys[0], kpb.crys[1],kpb.crys[2], kk.crys[0], kk.crys[1],
					kk.crys[2]);
					printf("b=(%5.3f,%5.3f,%5.3f). (ik=%4i, ib=%4i)\n", Bvector[ishell][ib][0], Bvector[ishell][ib][1], Bvector[ishell][ib][2], ik, ib);
				}
			}
		}
	}
}
    

int FindInside(int& ik, vector<int> const &cluster)
{
	for(int i=0; i<cluster.size(); i++)
		if( cluster[i] == ik ) 
			return i;
	return -1;
}


bool IsInside(int& ikpb, vector<int> const &cluster)
{
	for(int i : cluster)
		if( i == ikpb ) 
			return true;
	return false;
}






//this function reorganize the k space
void Organize_kspace(vector<vector<int>>& ik_central, vector<vector<vector<int>>>& ik_border, vector<vector<vector<int>>>& ik_external,vec2d& kpt, 
	vector<vector<vector<int>>>& GradientIndex_Total, int& rank, int& nproc)
{
	//find k points for every cpu:
	//1st line -> the index of the processor is  (0-1)*nproc mod(nproc)    floor
	ik_central.resize(nproc);
	ik_border.resize(nproc, vector<vector<int>>(nproc));
	ik_external.resize(nproc, vector<vector<int>>(nproc));

	for(int ik=0; ik<kpt.n1(); ik++)
	{
		int irank = (int)floor(fmod(((kpt[ik][1]+0.5)*nproc),nproc));  
		ik_central[irank].push_back(ik);
	}

	cout << "rank " << rank << ": I have " << ik_central[rank].size() << " k points in this cpu\n";
    stringstream sname; sname << rank;
    ofstream fp_wf; fp_wf.open("kpt_" + sname.str() + ".txt");
    for(int ik : ik_central[rank])
    	fp_wf << ik << " " <<kpt[ik][0] << " " << kpt[ik][1] << " " << kpt[ik][2] << endl;
    fp_wf.close();


	cout << "barrier1!" << endl;
	MPI_Barrier(MPI_COMM_WORLD);

	//push at the end of the array points on the border of the fractionized k space
	for(int irank=0; irank<nproc; irank++)
	{
		for(int ik1=0; ik1<ik_central[irank].size(); ik1++)
		{
			int ik_local = ik_central[irank][ik1];
			bool IsAtBorder=false;
			for(int ishell=0; ishell<GradientIndex_Total[ik_local].size(); ishell++)
				for(int& ikshell : GradientIndex_Total[ik_local][ishell])					
				{
					bool IsNotInside = ( !(IsInside(ikshell,ik_central[irank])) && !(IsInside(ikshell,ik_border[irank][irank]) ) );
					IsAtBorder = (IsAtBorder || IsNotInside);
					if(IsNotInside && !IsInside(ikshell, ik_external[irank][irank]))
						ik_external[irank][irank].push_back(ikshell);
				}
			if(IsAtBorder) //we enter here when we dont find at least a k+b inside the cpu
			{
				if( !IsInside(ik_central[irank][ik1],ik_border[irank][irank]) )
					ik_border[irank][irank].push_back(ik_central[irank][ik1]);
				ik_central[irank].erase(ik_central[irank].begin()+ik1);
				ik1--;		
			}
		}
	}//end of reorganization of k space inside the cpus
	cout << "barrier!" << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	//here we consider only connections between 2 cpus, no more!!!
	for(int irank=0; irank<nproc; irank++)
	{
		for(int irank1=0; irank1<nproc; irank1++)
		{
			if(irank1==irank) continue;
			for(int ik=0; ik<ik_border[irank][irank].size(); ik++)
			{
				int ik_actual = ik_border[irank][irank][ik];
				int index = FindInside(ik_actual, ik_external[irank1][irank1]);
				if (index != -1)
				{
					ik_border[irank][irank1].push_back(ik_actual);
					ik_border[irank][irank].erase(ik_border[irank][irank].begin()+ik);
					ik_external[irank1][irank].push_back(ik_actual);
					ik_external[irank1][irank1].erase(ik_external[irank1][irank1].begin()+index);
					ik--;
				}
			}
		}
	}


	fp_wf.open("kptcenter_" + sname.str() + ".txt");
    for(int ik : ik_central[rank])
    	fp_wf << ik << " " <<kpt[ik][0] << " " << kpt[ik][1] << " " << kpt[ik][2] << endl;
	fp_wf.close();
	
	for(int irank1=0; irank1<nproc; irank1++)
	{	
		stringstream sname1; sname1 << irank1;
		fp_wf.open("kptborder_"+ sname.str()+ sname1.str() +".txt");
    	for(int ik : ik_border[rank][irank1])
    		fp_wf << ik << " " <<kpt[ik][0] << " " << kpt[ik][1] << " " << kpt[ik][2] << endl;
		fp_wf.close();
	}
	for(int irank1=0; irank1<nproc; irank1++)
	{	
		stringstream sname1; sname1 << irank1;
		fp_wf.open("kptextern_"+ sname.str()+ sname1.str() +".txt");
    	for(int ik : ik_external[rank][irank1])
    	fp_wf << ik << " " <<kpt[ik][0] << " " << kpt[ik][1] << " " << kpt[ik][2] << endl;
		fp_wf.close();
	}

}


//void Organize_Gradient(int& nk_border, vector<vector<int>>& ik_cpu, vector<vector<vector<int>>>& GradientIndex_Total, int& rank, int& nproc)
//{
//	int nk_local = ik_cpu[rank].size();
//	GradientIndex.resize(nk_local);
//
//	for(int ik=0; ik<nk_local; ik++)
//	{
//		int ik_local = ik_cpu[rank][ik];
//		GradientIndex[ik].resize(GradientIndex_Total[ik_local].size());
//		for(int ishell=0; ishell<GradientIndex_Total[ik_local].size(); ishell++)
//		{
//			GradientIndex[ik][ishell].resize(GradientIndex_Total[ik_local][ishell].size());
//			for(int ikpb=0;  ikpb<GradientIndex_Total[ik_local][ishell].size(); ikpb++)
//				if(ik < nk_local-nk_border)
//					GradientIndex[ik][ishell][ikpb] = FindInside(GradientIndex_Total[ik_local][ishell][ikpb], ik_cpu[rank]);
//				else
//		}
//	}
//}

//mpi version
void CalculateIndices(vec2d& kpt, vector<vector<vector<int>>>& GradientIndex, vector<vector<vector<double>>>& Bvector, int& rank, int& nproc)
{
	vector<vector<vector<int>>> GradientIndex_Total;
	vector<vector<vector<int>>> GradientIndex_whichCPU;
	vector<vector<int>> ik_central; 
	vector<vector<vector<int>>> ik_border, ik_external;
	cout << "I did it!"<< endl;
	CalculateIndices(kpt, GradientIndex_Total, Bvector);
	cout << "I calculated gradient" << endl;

	Organize_kspace(ik_central, ik_border, ik_external, kpt, GradientIndex_Total, rank, nproc);	
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "I organized the k space\n";
    exit(1);





}
