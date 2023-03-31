
//loop for saving matrices
printf("Saving matrices H,U,D...\n");
#pragma omp parallel for schedule(dynamic)
for(int ik=0; ik<nk; ik++)
{
// cout << "ik: " << ik << endl;
    Coord_B ktemp;
    vec2x Hk(Ncv,Ncv);
    vec2x Uk(Ncv,Ncv);
    vec3x Dk(Ncv,Ncv,3); //vec2x Dk(Ncv,Ncv);
	ktemp.setcrys(kpt[ik][0], kpt[ik][1], kpt[ik][2]);
    gHUD(Hk,Uk,Dk,ktemp,iMode);
    
    //if(iTightBinding)  TB_Model.energy_U(Hk,Uk,ktemp);
    //else               WModel.energy_U(Hk,Uk,ktemp);
    if(ik==nk/4 || ik==nk/2 || ik==3*nk/4) printf("check point %3i\n", ik);
    for(int ic=0; ic<Ncv; ic++)
    {
    	for(int jc=0; jc<Ncv; jc++)
    	{
    		Hamiltonian[ik][ic][jc] = Hk[ic][jc];
    		Unitary[ik][ic][jc] = Uk[ic][jc];
            Dipole[ik][ic][jc][0] = Dk[ic][jc][0];
            Dipole[ik][ic][jc][1] = Dk[ic][jc][1];
            Dipole[ik][ic][jc][2] = Dk[ic][jc][2];
            if(iMode=="W"){
                Dipole[ik][ic][jc][0] = 0.5*(Dk[ic][jc][0]+Dk[jc][ic][0]);
                Dipole[ik][ic][jc][1] = 0.5*(Dk[ic][jc][1]+Dk[jc][ic][1]);
                Dipole[ik][ic][jc][2] = 0.5*(Dk[ic][jc][2]+Dk[jc][ic][2]);
            }
    	}
    }
}
printf("end of saving matrices. \n");

Mpi_communication(Hamiltonian,  message);
MPI_Barrier(MPI_COMM_WORLD);
    
if(iCurrent)
{ // works only for one thread
    //
    printf("Saving gradient of energy...\n");
    GradientEnergy.resize(nk,Ncv,Ncv,3); 
    GradientEnergy.fill(0);

    for(int i=0; i<3; i++)
        for(int ik=0; ik<nk; ik++)
            for(int ic=0; ic<Ncv; ic++)
                for(int jc=0; jc<Ncv; jc++)
                    for(int ishell=0; ishell<Bvector.size(); ishell++)
                        for(int ib=0; ib<Bvector[ishell].size(); ib++){
                            GradientEnergy[ik][ic][jc][i] += Weigths[ishell]*Bvector[ishell][ib][i]*
                            (Hamiltonian[GradientIndex[ik][ishell][ib]][ic][jc]-Hamiltonian[ik][ic][jc]);
                        }
                            
    printf("done\n");
}




