class Crystal
{
	public:
		void init(vec2d&,vec2d&,int,int);
		multivec2D<FittedData<complexd>> Ham;
		multivec3D<FittedData<complexd>> D;
		complexd dipolex(int& ic, int& jc, Coord_B& k); 
		complexd dipoley(int& ic, int& jc, Coord_B& k); 
		complexd dipolez(int& ic, int& jc, Coord_B& k); 
		complexd energy(int& ic, int& jc, Coord_B& k); 
        void energy_U(vec2x& H, vec2x& U, Coord_B& k);
        void dipole(vec3x&, Coord_B&);
    private:
      int _Ncv;
      int _Nch;
};


void Crystal::init(vec2d& a,vec2d& b,int Nch, int Ncv)
{
    _Ncv = Ncv;
    _Nch = Nch;
	multivec3D<vec3x> BerryConnection;
	multivec2D<vec3x> Hamiltonian;
	vec2d kpt_crystal;
	ifstream fp_in; int NumberOfOrbitals, nk_file;
	string trash; size_t pos; vector<string> str;
	a.resize(3,3); b.resize(3,3); a.fill(0); b.fill(0);
	//reading kmesh_data.txt
	printf("reading \"kmesh_data.dat\"...");
	fp_in.open("kmesh_data.dat");
	//a0
	getline(fp_in,trash);
	std::cout << "line1" << trash << endl;
	getline(fp_in,trash); pos=trash.length();
    Separate_string(trash, str, pos);
	a[1][1] = atof(str[0].c_str()); a[1][2]=atof(str[1].c_str()); //!!!! change to zxy!!!
	//a1
	getline(fp_in,trash);
	getline(fp_in, trash); pos=trash.length();
	Separate_string(trash, str, pos);
	a[2][1] = atof(str[0].c_str()); a[2][2]=atof(str[1].c_str());
	a[0][0]=15.;
	//b0
	getline(fp_in,trash);
	getline(fp_in, trash); pos=trash.length();
	std::cout << "line6" << trash << endl;
	Separate_string(trash, str, pos);
	b[1][1] = atof(str[0].c_str()); b[1][2]=atof(str[1].c_str());
	//b1
	getline(fp_in,trash);
	getline(fp_in, trash); pos=trash.length();
	Separate_string(trash, str, pos);
	b[2][1] = atof(str[0].c_str()); b[2][2]=atof(str[1].c_str());
	b[0][0]=1./15.;
	std::cout << "line7" << trash << endl;
	Coord_B::set_crys_to_cart(b);
	//atomic_position

	do{
		getline(fp_in,trash);pos=trash.length();
		Separate_string(trash, str, pos);
		std::cout << "line_do" << trash << endl;
	}
	while (0 !=  strcasecmp( "basis", str[str.size()-1].c_str() )); 
	//orbitals
	getline(fp_in,trash);pos=trash.length();
	std::cout << "line8" << trash << endl;
//cout << trash;
	Separate_string(trash, str, pos);
	NumberOfOrbitals=atoi(str[0].c_str());
	//kpoints

	do{
		getline(fp_in,trash);pos=trash.length();
		Separate_string(trash, str, pos);
		std::cout << "line" << trash << endl;
	}
	while(0 !=  strcasecmp( "direction", str[str.size()-1].c_str() )); 
	std::cout << "line" << trash << endl;

	getline(fp_in,trash);pos=trash.length();
	Separate_string(trash, str, pos);
	nk_file=atoi(str[0].c_str());
	getline(fp_in,trash);
	kpt_crystal.resize(nk_file*nk_file,3); kpt_crystal.fill(0); 
	Coord_B kk;
	for(int ik=0; ik<nk_file*nk_file; ik++)
	{
		getline(fp_in,trash); pos=trash.length();
		Separate_string(trash, str, pos);
		kpt_crystal[ik][1]=atof(str[1].c_str()); // !!! change 
		kpt_crystal[ik][2]=atof(str[2].c_str());
		kk.setcart(kpt_crystal[ik][0],kpt_crystal[ik][1],kpt_crystal[ik][2]);
		// kpt_crystal[ik][1]=kk.crys[0];
		// kpt_crystal[ik][2]=kk.crys[1];
		for(int ix=0; ix<3; ix++)
			kpt_crystal[ik][ix]=kk.crys[ix];
	}
	fp_in.close();
	//printf("read.\n");
	//printf("resume of kmesh_data.dat..\n");
	//printf("a0 %10.6f %10.6f %10.6f\n",a[0][0],a[0][1],a[0][2]);
	//printf("a1 %10.6f %10.6f %10.6f\n",a[1][0],a[1][1],a[1][2]);
	//printf("b0 %10.6f %10.6f %10.6f\n",b[0][0],b[0][1],b[0][2]);
	//printf("b1 %10.6f %10.6f %10.6f\n",b[1][0],b[1][1],b[1][2]);
	//printf("NumberOfOrbitals %5i",NumberOfOrbitals);
	//printf("nk_file %5i",nk_file);
	//printf("Fractional coordinates of the defined k points \n");
	//for(int ik=0; ik<nk_file*nk_file; ik++)
	//	printf("%5i %10.6f %10.6f %10.6f\n", ik+1, kpt_crystal[ik][0],kpt_crystal[ik][1],kpt_crystal[ik][2]);
	

	//reading hk_orthonormal.dat    
	Hamiltonian.resize(NumberOfOrbitals,NumberOfOrbitals);//matrix for each k point
	for(int ic=0; ic<NumberOfOrbitals; ic++)
		for(int jc=0; jc<NumberOfOrbitals; jc++)
			Hamiltonian[ic][jc].resize(1,nk_file,nk_file);

	printf("reading \"hk_orthonormal.dat\"...");
	ifstream fp_H;
	fp_H.open("hk_orthonormal.dat");

	for(int ik1=0; ik1<nk_file; ik1++)
	{
		for(int ik2=0; ik2<nk_file; ik2++)
		{
			//if((ik2+ik1*nk_file)%100==0) printf("reading ik = %6i\n",ik2+ik1*nk_file); 
			for(int ic=0; ic<NumberOfOrbitals; ic++)
			{
				getline(fp_H, trash);   
				pos=trash.length();
				Separate_string(trash, str, pos);
				for(int jc=0; jc<ic+1; jc++)
					Hamiltonian[ic][jc][0][ik1][ik2]=atof(str[jc].c_str());
			}
			getline(fp_H,trash);
			for(int ic=0; ic<NumberOfOrbitals; ic++)
			{
				getline(fp_H, trash);    pos=trash.length();
	     		Separate_string(trash, str, pos);
				for(int jc=0; jc<ic+1; jc++)
					Hamiltonian[ic][jc][0][ik1][ik2]+=c1*atof(str[jc].c_str());
			}
			getline(fp_H,trash);
			getline(fp_H,trash);
		}
	}
	fp_H.close();
	printf("read.\n");


	//reading berry connections
	BerryConnection.resize(NumberOfOrbitals,NumberOfOrbitals,3);//matrix for each k point
	for(int ic=0; ic<NumberOfOrbitals; ic++)
		for(int jc=0; jc<NumberOfOrbitals; jc++)
			for(int ix=0; ix<3; ix++)
				BerryConnection[ic][jc][ix].resize(1,nk_file,nk_file);
	
	printf("reading \"bc_betabetap_real.dat\" and \"bc_betabetap_imag.dat\"...");
	ifstream fp_R, fp_I;
	fp_R.open("bc_betabetap_real.dat");
	fp_I.open("bc_betabetap_imag.dat");
	string trashR, trashI;
	vector<string> temp_r, temp_i;
	for(int ik1=0; ik1<nk_file; ik1++)
	{
		for(int ik2=0; ik2<nk_file; ik2++)
		{
			//if((ik2+ik1*nk_file)%100==0) printf("reading ik = %6i\n",ik2+ik1*nk_file); 
			for(int ic=0; ic<NumberOfOrbitals; ic++)
			{
				getline(fp_R, trashR);   
				getline(fp_I, trashI);   
				pos=trashR.length();
				Separate_string(trashR, temp_r, pos);
				pos=trashI.length();
				Separate_string(trashI, temp_i, pos);
				for(int jc=0; jc<ic+1; jc++)
					BerryConnection[ic][jc][0][0][ik1][ik2]=atof(temp_r[jc].c_str())+c1*atof(temp_i[jc].c_str());
			}
			getline(fp_R, trash); getline(fp_I, trash);
			for(int ic=0; ic<NumberOfOrbitals; ic++)
			{
				getline(fp_R, trashR);   
				getline(fp_I, trashI);   
				pos=trashR.length();
				Separate_string(trashR, temp_r, pos);
				pos=trashI.length();
				Separate_string(trashI, temp_i, pos);
				for(int jc=0; jc<ic+1; jc++)
					BerryConnection[ic][jc][1][0][ik1][ik2]=atof(temp_r[jc].c_str())+c1*atof(temp_i[jc].c_str());
			}
			getline(fp_R, trash); getline(fp_I, trash);
			getline(fp_R, trash); getline(fp_I, trash);
		}
	}
	fp_R.close(); fp_I.close();
	printf("read.\n");
   

	//cout << BerryConnection[15][12][0][0][nk_file-1][nk_file-1] << " " << BerryConnection[15][12][1][0][nk_file-1][nk_file-1]<<endl;
	//fitting
	printf("interpolating...\n");
	Ham.resize(NumberOfOrbitals,NumberOfOrbitals);
	D.resize(NumberOfOrbitals,NumberOfOrbitals,3);
	vec1d spacing(3); for(int i=0; i<3; i++) spacing[i]=1./(nk_file-1);
	vec1d xmin(3);    for(int i=0; i<3; i++) xmin[i]=-0.5;
	for(int ic=0; ic<NumberOfOrbitals;ic++)
		for(int jc=0; jc<ic+1; jc++)
		{
			//printf("orbitals= %3i ic= %3i jc= %3i\n",NumberOfOrbitals, ic,jc);
			Ham[ic][jc].construct(Hamiltonian[ic][jc], spacing, xmin);
			for(int ix=0; ix<3; ix++)
				D[ic][jc][ix].construct(BerryConnection[ic][jc][ix],spacing,xmin);
		}	
}



complexd Crystal::dipolex(int& ic, int& jc, Coord_B& k)
{
	if(jc<ic+1)
		return(D[ic][jc][0](k.crys[0],k.crys[1],k.crys[2]));
        else
		return(conj(D[jc][ic][0](k.crys[0],k.crys[1],k.crys[2])));
}



complexd Crystal::dipoley(int& ic, int& jc, Coord_B& k)
{
	if(jc<ic+1)
		return(D[ic][jc][1](k.crys[0],k.crys[1],k.crys[2]));
        else
		return(conj(D[jc][ic][1](k.crys[0],k.crys[1],k.crys[2])));
}
 


complexd Crystal::dipolez(int& ic, int& jc, Coord_B& k)
{
	if(jc<ic+1)
		return(D[ic][jc][2](k.crys[0],k.crys[1],k.crys[2]));
        else
		return(conj(D[jc][ic][2](k.crys[0],k.crys[1],k.crys[2])));
}
 


complexd Crystal::energy(int& ic, int& jc, Coord_B& k)
{
	if(jc<ic+1)
		return(Ham[ic][jc](k.crys[0],k.crys[1],k.crys[2]));
        else
		return(conj(Ham[jc][ic](k.crys[0],k.crys[1],k.crys[2])));
}

void  Crystal::dipole(vec3x& Dx, Coord_B& k)
{
    for(int ic = 0; ic < _Ncv; ic++)
    {
        for(int jc = 0; jc<_Ncv; jc++)
        {
//cout <<  "DIPOLE: " << ic << " " << jc << endl; 
            Dx[ic][jc][0] = dipolez(ic, jc, k); //CHANGED THIS!!! the order should be zxy
            Dx[ic][jc][1] = dipolex(ic, jc, k);
            Dx[ic][jc][2] = dipoley(ic, jc, k);
            if(ic==jc)
              for(int i=0; i<3; i++) Dx[ic][jc][i].imag(0.);
        }
    }
}


void Crystal::energy_U(vec2x& H, vec2x& Uk, Coord_B& k)    
{
	H.fill(0.);
	for(int ic = 0; ic < H.n1(); ic++)
	    for(int jc = 0; jc<H.n2(); jc++)
	       H[ic][jc] = energy(ic, jc, k);
	vec epsilon;
	cx_mat U;
	cx_mat Hw; 
	epsilon.zeros(H.n1());
	U.zeros(H.n1(), H.n1());
	Hw.zeros(H.n1(), H.n1());
	//copy the Ek matrix in an armadillo matrix
	for (int ic=0; ic<H.n1(); ic++)
	{
	    for (int jc=0; jc <H.n2(); jc++)
	    {
	        Hw(ic, jc) = H[ic][jc];
	        // Hw(ic, jc) *= 1.0e9 ;
	    }
	}
	eig_sym(epsilon, U, Hw);

	Uk.fill(0.);
	
	for (int ic=0; ic<H.n1(); ic++)
	    for (int jc=0; jc <H.n2() ; jc++)
	        Uk[ic][jc]=conj(U(jc,ic));

}
