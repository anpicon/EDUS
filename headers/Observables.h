void PrintEFAF(ofstream& fp_E, double& time, Laser& pulse)
{
    vec1d E(3); vec1d A(3); 
    E = pulse.E(time);
    A = pulse.A(time);
    
    fp_E <<  setw(25) << setprecision(15)<< scientific <<  time*time_au_fs;
    fp_E <<  setw(25) << setprecision(15)<< scientific << E[0];
    fp_E <<  setw(25) << setprecision(15)<< scientific << E[1];
    fp_E <<  setw(25) << setprecision(15)<< scientific << E[2];

    fp_E <<  setw(25) << setprecision(15)<< scientific << A[0];
    fp_E <<  setw(25) << setprecision(15)<< scientific << A[1];
    fp_E <<  setw(25) << setprecision(15)<< scientific << A[2] << endl;
}




void CalculateLosses(int& nk, vec1d& sumNcv, int& Nch, int Ncc,
    vec3x& P0,double& time,vec3x& U,  double &  P_cond_max_loc_mpi)
{   
    P_cond_max_loc_mpi = 0.0;
    int Ncv = P0.n2();
    sumNcv.fill(0.);

    for (int ik=0; ik<nk; ik++)
    {

        for (int ic=0; ic<Nch; ic++)
              sumNcv[ic]+=1.-abs(P0[ik][ic][ic]);

        for (int ic=Nch; ic<P0.n2(); ic++)
        {
            complex<double> sum=0.;
            for (int ii=Nch; ii<P0.n2(); ii++)
                for (int jj=Nch; jj<P0.n2(); jj++)
                    sum+=conj(U[ik][ic][ii])*P0[ik][ii][jj]*U[ik][ic][jj]; 

                    

            if(ic<Ncc) sum=1.-sum;//with4 this line we consider holes in the ic-th band       
            
            if (abs(sum) > P_cond_max_loc_mpi){
                P_cond_max_loc_mpi = abs(sum);
            }       
            
            // find maximum
            sumNcv[ic]+= abs(sum);//We sum over all the k points. sum is for a particular k, sumNcv sum it.
        }//End ic
    }// End ik
}



void PrintLosses(int& nk,ofstream& fp_Loss, int& Nch, int Ncc,vec3x& P0,double& time,vec3x& U)
{
    vec1d sumNcv(P0.n2());
    double P_cond_max;
    CalculateLosses(nk,sumNcv,Nch,Ncc,P0,time,U, P_cond_max);
    fp_Loss << setw(25) << setprecision(15) << time*time_au_fs;
    for(int ic=0; ic<P0.n2(); ic++) fp_Loss << setw(25) << setprecision(15) <<scientific << sumNcv[ic]/P0.n1();
    fp_Loss << endl;    
}



void CalculateTransientAbs(int& nk,vec1x& a,vec3x& P, double& time, vec1d& EFx,  vec1d& dk, int& Nch,double detM, vec4x& Dk)
{
    a.fill(0.);
    for( int ik = 0; ik <nk; ik++ )
        for (int ic=0; ic<Nch; ic++) 
            for (int jc=Nch; jc<P.n2(); jc++) 
                for(int ix=0; ix<3; ix++)
                    a[ix] -= 2.*Dk[ik][ic][jc][ix]*P[ik][ic][jc];
    for(int ix=0; ix<3; ix++)
        a[ix] *= abs(detM)*dk[0]*dk[1]*dk[2]; 
}

void PrintTransientAbs(int&nk,ofstream& fp_J, vec3x& P, double& time, vec1d& EFx,  vec1d& dk, int& Nch,double detM, vec4x& Dk)
{
    vec1x a(3); 
    CalculateTransientAbs(nk,a,P,time,EFx,dk,Nch,detM,Dk);
    fp_J << setw(25) << setprecision(15) << scientific << time*time_au_fs;
    fp_J << setw(25) << setprecision(15) << scientific << EFx[0];
    fp_J << setw(25) << setprecision(15) << scientific << EFx[1];
    fp_J << setw(25) << setprecision(15) << scientific << EFx[2];
    for(int ix=0; ix<3; ix++) 
    {
        fp_J << setw(25) << setprecision(15) << scientific << a[ix].real()/P.n1();
        fp_J << setw(25) << setprecision(15) << scientific << a[ix].imag()/P.n1();
    }
    fp_J << endl;
}//END of TransientAbs



void PrintEF(ofstream& fp_E, double& time, vector<vec1d>& EF)
{
    fp_E <<  setw(15) << setprecision(8)<< time*time_au_fs;
    fp_E <<  setw(15) << setprecision(8)<< EF[0][0];
    fp_E <<  setw(15) << setprecision(8)<< EF[0][1];
    fp_E <<  setw(15) << setprecision(8)<< EF[0][2] << endl;
}



//det_a - determinant UC area
void Current1(ofstream& fp_J1, Private_omp_parameters& OMP_private,
    vec3x& P,vec4x& GradientEnergy,double& time, vec1d& dk,
    int& nktot,double det_a, vec1x& I)
{
    #pragma omp master
    {
        I.fill(0.);
    }
    #pragma omp barrier
    int ik;
    vec1x Ipriv(3); Ipriv.fill(0.);
    for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++)
    {
        ik = ik_pr + OMP_private.begin_count;
        for( int ic = 0; ic < P.n2(); ic++ )
            for ( int jc = 0; jc < P.n2(); jc++ )
                for (int ix=0; ix<3; ix++)
                    Ipriv[ix]+=P[ik][ic][jc]*GradientEnergy[ik][ic][jc][ix];
    }
    #pragma omp barrier
    #pragma omp critical
    { // summation over threads
        for (int ix=0; ix<3; ix++){
            I[ix] += Ipriv[ix];
        }
    }
    #pragma omp barrier
    #pragma omp master
    {
        for (int ix=0; ix<3; ix++)
            I[ix] *= (dk[0]*dk[1]*dk[2]/det_a);
        
        vec1x global_sum(3);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(&I[0], &global_sum[0], 3, MPI_C_DOUBLE_COMPLEX, MPI_SUM, 0,
               MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank_==0)
        {
            fp_J1 << setw(15) << setprecision(8) << time*time_au_fs;
            for(int ix=0; ix<3; ix++)
            {
                fp_J1 << setw(25) << setprecision(15) << scientific << global_sum[ix].real();
                fp_J1 << setw(25) << setprecision(15) << scientific << global_sum[ix].imag();
            }
            fp_J1 << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    #pragma omp barrier
}

void Current2(ofstream& fp_J2, Private_omp_parameters& OMP_private,
    vec3x& P,vec4x& Dipole,vec3x& Hamiltonian,
    double& time, vec1d& dk,int& nktot,double det_a, vec1x& I)
{
    #pragma omp master
    {
        I.fill(0.);
    }
    int ik;
    vec1x Ipriv(3); Ipriv.fill(0.);
    #pragma omp barrier
    for (int ik_pr = 0; ik_pr < OMP_private.lenght_k; ik_pr++)
    {
        ik = ik_pr + OMP_private.begin_count;
        for( int ic = 0; ic < P.n2(); ic++ )
            for ( int jc = 0; jc < P.n2(); jc++ )
                for ( int kc = 0; kc < P.n2(); kc++ )
                    for (int ix=0; ix<3; ix++)
                        Ipriv[ix]+=P[ik][ic][jc]*Hamiltonian[ik][ic][kc]*conj(Dipole[ik][jc][kc][ix]);
    }
    #pragma omp barrier
    #pragma omp critical
    { // summation over threads
        for (int ix=0; ix<3; ix++){
            I[ix] += Ipriv[ix];
        }
    }
    #pragma omp barrier
    #pragma omp master
    {
        for (int ix=0; ix<3; ix++)
            I[ix] *= (dk[0]*dk[1]*dk[2]/det_a);
        
        vec1x global_sum(3);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(&I[0], &global_sum[0], 3, MPI_C_DOUBLE_COMPLEX, MPI_SUM, 0,
               MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank_==0)
        {
            fp_J2 << setw(15) << setprecision(8) << time*time_au_fs;
            for(int ix=0; ix<3; ix++)
            {
                fp_J2 << setw(25) << setprecision(15) << scientific << -2.*global_sum[ix].imag();
            }
            fp_J2 << endl;
        }
    }
    #pragma omp barrier
}


void CalculateCurrent1(vec3x& P,vec4x& GradientEnergy,vec1d& dk,int& nk,double det_a, vec1x& I)
{
    I.fill(0.);
    for(int ik=0; ik<nk; ik++)
        for( int ic = 0; ic < P.n2(); ic++ )
            for ( int jc = 0; jc < P.n2(); jc++ )
                for (int ix=0; ix<3; ix++)
                    I[ix]+=P[ik][ic][jc]*GradientEnergy[ik][ic][jc][ix];

    for (int ix=0; ix<3; ix++)
        I[ix] *= (dk[0]*dk[1]*dk[2]/det_a);
}


void CalculateCurrent2(vec3x& P,vec4x& Dipole,vec3x& Hamiltonian,vec1d& dk,int& nk,double det_a, vec1x& I)
{
    I.fill(0.);
    for(int ik=0; ik<nk; ik++)
        for( int ic = 0; ic < P.n2(); ic++ )
            for ( int jc = 0; jc < P.n2(); jc++ )
                for ( int kc = 0; kc < P.n2(); kc++ )
                    for (int ix=0; ix<3; ix++) 
                        I[ix] += P[ik][ic][jc]*Hamiltonian[ik][ic][kc]*conj(Dipole[ik][jc][kc][ix]);
    for(int ix=0; ix<3; ix++)
        I[ix] *= (dk[0]*dk[1]*dk[2]/det_a);
}



void PrintCurrent(ofstream& fp_J1, ofstream& fp_J2, double time,vec3x& P,vec4x& GradientEnergy,vec4x& Dipole,vec3x& Hamiltonian,vec1d& dk,int& nk, double det_a)
{
    vec1x I1(3); 
    CalculateCurrent1(P,GradientEnergy,dk,nk,det_a,I1);
    fp_J1 << setw(15) << setprecision(8) << time*time_au_fs;
    for(int ix=0; ix<3; ix++) 
    {
        fp_J1 << setw(25) << setprecision(15) << scientific << I1[ix].real();
        fp_J1 << setw(25) << setprecision(15) << scientific << I1[ix].imag();
    }
    fp_J1 << endl;
    
    vec1x I2(3);
    CalculateCurrent2(P,Dipole,Hamiltonian,dk,nk,det_a,I2);
    fp_J2 << setw(15) << setprecision(8) << time*time_au_fs;
    for(int ix=0; ix<3; ix++) 
        fp_J2 << setw(25) << setprecision(15) << scientific << -2*I2[ix].imag();   
    fp_J2 << endl;

}

 


void PrintKFD(vec3x& P0,vec1i& Nk,int Ncch, vec2d& k, 
    Coulomb_parameters& Coulomb_set,
    vec3x& U,int& icont)
{
    stringstream sname;     sname << icont;
    ofstream fp_J; fp_J.open("Output/kdist_" + sname.str() + ".txt");

    vec1x sumNcv(P0.n2());    sumNcv.fill(0.);
    for (int ik=0; ik<P0.n1(); ik++)
    {        
        sumNcv.fill(0.);
        for (int ic=0; ic<P0.n2(); ic++)
            for (int ii=0; ii<P0.n2(); ii++)
                for (int jj=0; jj<P0.n2(); jj++)
                    sumNcv[ic]+=conj(U[ik][ic][ii])*P0[ik][ii][jj]*U[ik][ic][jj];        
        fp_J << k[ik][0] << " " << k[ik][1] << " " << k[ik][2];
        for (int ic=0; ic<P0.n2(); ic++)
        {
            if(ic<Ncch) fp_J << " " << 1.-abs(sumNcv[ic]);
            else fp_J << " " << abs(sumNcv[ic]);
        }
        fp_J << endl;
        fp_J.close();
        icont++;
    }//END of ik loop
}








