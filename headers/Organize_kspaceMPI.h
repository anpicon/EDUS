

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
void Organize_kspace(vector<vector<int>>& ik_central, vector<vector<vector<int>>>& ik_border, 
	vector<vector<vector<int>>>& ik_external,vec2d& kpt, 
	vector<vector<vector<int>>>& GradientIndex_Total)
{
	//find k points for every cpu:
	//1st line -> the index of the processor is  (0-1)*num_procs mod(num_procs)    floor
	ik_central.resize(num_procs);
	ik_border.resize(num_procs, vector<vector<int>>(num_procs));
	ik_external.resize(num_procs, vector<vector<int>>(num_procs));

	for(int ik=0; ik<kpt.n1(); ik++)
	{
		int irank_ = (int)floor(fmod(((kpt[ik][1]+0.5)*num_procs),num_procs));  
		ik_central[irank_].push_back(ik);
	}

	cout << "rank_ " << rank_ << ": I have " << ik_central[rank_].size() << " k points in this cpu\n";
    stringstream sname; sname << rank_;
    ofstream fp_wf; fp_wf.open("kpt_" + sname.str() + ".txt");
    for(int ik : ik_central[rank_])
    	fp_wf << ik << " " <<kpt[ik][0] << " " << kpt[ik][1] << " " << kpt[ik][2] << endl;
    fp_wf.close();

	MPI_Barrier(MPI_COMM_WORLD);

	//push at the end of the array points on the border of the fractionized k space
	for(int irank_=0; irank_<num_procs; irank_++)
	{
		for(int ik1=0; ik1<ik_central[irank_].size(); ik1++)
		{
			int ik_local = ik_central[irank_][ik1];
			bool IsAtBorder=false;
			for(int ishell=0; ishell<GradientIndex_Total[ik_local].size(); ishell++)
				for(int& ikshell : GradientIndex_Total[ik_local][ishell])					
				{
					bool IsNotInside = ( !(IsInside(ikshell,ik_central[irank_])) && !(IsInside(ikshell,ik_border[irank_][irank_]) ) );
					IsAtBorder = (IsAtBorder || IsNotInside);
					if(IsNotInside && !IsInside(ikshell, ik_external[irank_][irank_]))
						ik_external[irank_][irank_].push_back(ikshell);
				}
			if(IsAtBorder) //we enter here when we dont find at least a k+b inside the cpu
			{
				if( !IsInside(ik_central[irank_][ik1],ik_border[irank_][irank_]) )
					ik_border[irank_][irank_].push_back(ik_central[irank_][ik1]);
				ik_central[irank_].erase(ik_central[irank_].begin()+ik1);
				ik1--;		
			}
		}
	}//end of reorganization of k space inside the cpus
	MPI_Barrier(MPI_COMM_WORLD);
	//here we consider only connections between 2 cpus, no more!!!
	for(int irank_=0; irank_<num_procs; irank_++)
	{
		for(int irank_1=0; irank_1<num_procs; irank_1++)
		{
			if(irank_1==irank_) continue;
			for(int ik=0; ik<ik_border[irank_][irank_].size(); ik++)
			{
				int ik_actual = ik_border[irank_][irank_][ik];
				int index = FindInside(ik_actual, ik_external[irank_1][irank_1]);
				if (index != -1)
				{
					ik_border[irank_][irank_1].push_back(ik_actual);
					ik_border[irank_][irank_].erase(ik_border[irank_][irank_].begin()+ik);
					ik_external[irank_1][irank_].push_back(ik_actual);
					ik_external[irank_1][irank_1].erase(ik_external[irank_1][irank_1].begin()+index);
					ik--;
				}
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	fp_wf.open("kptcenter_" + sname.str() + ".txt");
    for(int ik : ik_central[rank_])
    	fp_wf << ik << " " <<kpt[ik][0] << " " << kpt[ik][1] << " " << kpt[ik][2] << endl;
	fp_wf.close();
	
	for(int irank_1=0; irank_1<num_procs; irank_1++)
	{	
		stringstream sname1; sname1 << irank_1;
		fp_wf.open("kptborder_"+ sname.str()+ sname1.str() +".txt");
    	for(int ik : ik_border[rank_][irank_1])
    		fp_wf << ik << " " <<kpt[ik][0] << " " << kpt[ik][1] << " " << kpt[ik][2] << endl;
		fp_wf.close();
	}
	for(int irank_1=0; irank_1<num_procs; irank_1++)
	{	
		stringstream sname1; sname1 << irank_1;
		fp_wf.open("kptextern_"+ sname.str()+ sname1.str() +".txt");
    	for(int ik : ik_external[rank_][irank_1])
    	fp_wf << ik << " " <<kpt[ik][0] << " " << kpt[ik][1] << " " << kpt[ik][2] << endl;
		fp_wf.close();
	}

}



void Organize_gradient(vector<int>& kindex, vector<vector<vector<int>>>& GradientIndex_Total, 
    	vector<vector<vector<int>>>& GradientIndex, int& nk)
{	
	for(int ik=0; ik<nk; ik++)
	{
		GradientIndex.push_back(vector<vector<int>>());
		int ik_local = kindex[ik];
		GradientIndex[ik].resize(GradientIndex_Total[ik_local].size());
		for(int ishell=0; ishell<GradientIndex_Total[ik_local].size(); ishell++)
		{
			GradientIndex[ik][ishell].resize(GradientIndex_Total[ik_local][ishell].size());
			for(int ikpb =0; ikpb<GradientIndex_Total[ik_local][ishell].size(); ikpb++)
			{
				int index = FindInside(GradientIndex_Total[ik_local][ishell][ikpb], kindex);
				if(index==-1) {printf("\n %3i Error! Inexistent k point in gradient!", GradientIndex_Total[ik_local][ishell][ikpb]); exit(1);}
				else
					GradientIndex[ik][ishell][ikpb] = index; 
			}		
		}
	}
}







int whereis(int& i, vector<vector<int>>& Indices)
{
	for(int index =0; index< Indices.size(); index++)
		if(Indices[index][0] == i) return index;
	return -1;
}




void Organize_MPI(vector<int>& ik_central, vector<vector<int>> const & ik_border, vector<vector<int>> const &ik_external, 
	//vector<vector<int>>& IndicesMPIsend, vector<vector<int>>& IndicesMPIreceive)
	vector<Message>& message)
{
// what we suppose here is that we have communications only with rank_-1 and rank_+1.
//first we distribute as (0-1) (2-3) (4-5).. then as (5-0) (1-2) (3-4)
//this is done thanks to integer division by 2
//the message is sent first by the cpu on the left, this is done with %2
	int start = ik_central.size();

	vector<vector<int>> IndicesMPIsend, IndicesMPIreceive;

	for(int irank_=0; irank_<num_procs; irank_++)
	{
		if(ik_border[irank_].size()!=0)
		{
			int SendMessageValue = ( rank_>irank_ ? 100 : 3000 );
			IndicesMPIsend.push_back(vector<int>()); int index = IndicesMPIsend.size()-1;
			IndicesMPIsend[index].push_back(irank_);
			IndicesMPIsend[index].push_back(start);
			IndicesMPIsend[index].push_back((int)ik_border[irank_].size());
			IndicesMPIsend[index].push_back(SendMessageValue + irank_ + num_procs*rank_ );
			start += ik_border[irank_].size();
		}
	}

	for(int irank_=0; irank_<num_procs; irank_++)
	{
		if(ik_external[irank_].size()!=0)
		{
			int ReceiveMessageValue = ( rank_<irank_ ? 100 : 3000 );
			IndicesMPIreceive.push_back(vector<int>()); int index = IndicesMPIreceive.size()-1;
			IndicesMPIreceive[index].push_back(irank_);
			IndicesMPIreceive[index].push_back(start);
			IndicesMPIreceive[index].push_back((int)ik_external[irank_].size());
			IndicesMPIreceive[index].push_back(ReceiveMessageValue + rank_ + num_procs*irank_);
			start += ik_external[irank_].size();
		}

	}

	//this is valid for toroidal type mpi processes
	for(int i=0; i<2; i++)
	{
		int irank_ = ((rank_+i)%2 == 0 ? (rank_+1) : (rank_-1))%num_procs;
		if(irank_<0) irank_+= num_procs;
		for(int j=0; j<2; j++)
		{
			bool send = (irank_+j)%2;
			if(send)
			{	
				int index = whereis(irank_,IndicesMPIsend);
				message.resize(message.size()+1);
				message[message.size()-1].send          = true;
				message[message.size()-1].buf           = IndicesMPIsend[index][1];
				message[message.size()-1].count         = IndicesMPIsend[index][2];
				message[message.size()-1].partner_rank  = IndicesMPIsend[index][0];
				message[message.size()-1].tag           = IndicesMPIsend[index][3];
			}
			else
			{
				int index = whereis(irank_,IndicesMPIreceive);
				message.resize(message.size()+1);
				message[message.size()-1].send          = false;
				message[message.size()-1].buf           = IndicesMPIreceive[index][1];
				message[message.size()-1].count         = IndicesMPIreceive[index][2];
				message[message.size()-1].partner_rank  = IndicesMPIreceive[index][0];
				message[message.size()-1].tag           = IndicesMPIreceive[index][3];
			}
		}
	}


}



//mpi version
void CalculateIndicesMPI(vec2d& kpt, vector<vector<vector<int>>>& GradientIndex, 
	vector<vector<vector<double>>>& Bvector, vec1d& dk, int& nk,
	//vector<vector<int>>& IndicesMPIreceive, vector<vector<int>>& IndicesMPIsend)
	vector<Message>& message)
{
	vector<vector<vector<int>>> GradientIndex_Total;
	vector<vector<int>> ik_central; 
	vector<vector<vector<int>>> ik_border, ik_external;
	vector<int> krank_aligned;

	CalculateIndices(dk, kpt, GradientIndex_Total, Bvector);//from non mpi program

	Organize_kspace(ik_central, ik_border, ik_external, kpt, GradientIndex_Total);	
	
    MPI_Barrier(MPI_COMM_WORLD);
   
	krank_aligned = ik_central[rank_];
	
	for(int irank=0; irank<num_procs; irank++)
		krank_aligned.insert(krank_aligned.end(), ik_border[rank_][irank].begin(), ik_border[rank_][irank].end());
	

	nk = krank_aligned.size();
	for(int irank=0; irank<num_procs; irank++)
		krank_aligned.insert(krank_aligned.end(), ik_external[rank_][irank].begin(), ik_external[rank_][irank].end());

    Organize_gradient(krank_aligned, GradientIndex_Total, GradientIndex, nk);
    MPI_Barrier(MPI_COMM_WORLD);


     if(num_procs > 1){ // many nodes version
	    Organize_MPI(ik_central[rank_], ik_border[rank_], ik_external[rank_], message);//IndicesMPIsend, IndicesMPIreceive);

	    vector<vector<double>> kpt_temp; 
	    for(int ik=0; ik<kpt.n1(); ik++)
		{
			vector<double> temp(3); for(int i=0; i<3; i++) temp[i] = kpt[ik][i];
	    	kpt_temp.push_back(temp);
	    }
	    
	    kpt.resize(krank_aligned.size(),3);
	    for(int ik=0; ik<krank_aligned.size(); ik++)
	    	for(int i=0; i<3; i++)
	    	   	kpt[ik][i]=kpt_temp[krank_aligned[ik]][i];
	    MPI_Barrier(MPI_COMM_WORLD);

	    ofstream test;     stringstream srank; srank << rank_;

    
	    test.open("mpi_info_"+ srank.str() + ".txt");
	    for(Message& temp : message)
	    {
	    	test << setw(10) <<"tag " << temp.tag << " "; 
	    	test << setw(10) << (temp.send ? "send " : "receive ");
	    	test << setw(10) << temp.partner_rank << " " << temp.buf << " " << temp.count;
	    	test << endl;
	    }	
	}
    MPI_Barrier(MPI_COMM_WORLD);
}
/*
    ofstream test; 
    stringstream srank; srank << rank_;
    test.open("kaligned_"+srank.str()+".txt");
    for(int ik : krank_aligned)
    	test << ik << " " << kpt[ik][0] << " " << kpt[ik][1] << " " << kpt[ik][2] << endl;
    test.close();

    test.open("test" + srank.str()+".txt");
	srand (time(NULL));



	int ikrand = rand() % nk;

  	int ik_local = krank_aligned[ikrand];
	cout << "rank " << rank_ << " ik_local " << ik_local << endl; 
  	
  	test << ik_local << " " << kpt[ik_local][0] << " " << kpt[ik_local][1] << " " << kpt[ik_local][2] << endl;
  	for(int index : GradientIndex[ikrand][0])
  	{
  		cout << "rank "<< rank_<< " index " << index <<endl;
  		int ib = krank_aligned[index];
		cout << "rank " << rank_ << " iB " << ib << endl; 
		test << ib << " " << kpt[ib][0] << " " << kpt[ib][1] << " " << kpt[ib][2] << endl;

  	}
	test.close();  	  
    test.open("mpi_info_"+ srank.str() + ".txt");
    test << "SEND: \n";
    for(vector<int> send : IndicesMPIsend)
    {
    	for(int i : send)
    		test << i << "   " ;
    	test << endl;
    }	
    test << endl<<endl<<"RECEIVE: \n";
    for(vector<int> send : IndicesMPIreceive)
    {
    	for(int i : send)
    		test << i << "   " ;
    	test << endl;
    }	
    test.close();
 
*/
