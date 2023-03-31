int rank_=0, num_procs=1;

#include "Source_Main/include_headers.cpp"

int main (int argc, char* argv[])
{
    stringstream srank; srank<<"";
    printf("Running the program with 1 cpu and ");
    #pragma omp parallel sections
    {
            printf("%3i threads \n", omp_get_num_threads());
    }

    #include "Source_Main/variables.cpp"
    #include "Source_Main/input.cpp"

    CalculateIndices(kpt, GradientIndex, Bvector);
    nk= Nk[0]*Nk[1]*Nk[2];

    Coord_B kK;
    ofstream fp; fp.open("kcart.txt"); 
    for(int ik=0; ik<nk; ik++)    
    {
       kK.setcrys(kpt[ik][0], kpt[ik][1], kpt[ik][2]);
       fp << kK.cart[1] << " " << kK.cart[2] << endl;
    }
    fp.close();
    
    #include "Source_Main/test_matrix_k.cpp"
    #include "Source_Main/allocate_arr.cpp"
    #include "Source_Main/open_files.cpp"
    #include "Source_Main/save_HUD.cpp"
    //#include "Source_Main/print_HUD.cpp"

    #include "Source_Main/runge_kutta_init.cpp"

int cnt=1; 
    for(int it=iti;it<=itfi;it++)
    {
        double time=it*dt;

        if ((it-iti) % it_resolution == 0)
        {
            PrintEF(fp_E, time, EF);
            PrintLosses(nk,fp_Loss, Nb[0], (Nb[0]+Nb[1]), P0,time,Unitary);
            //if(iTAbs)
            //    PrintTransientAbs(nk,fp_TAbs, P0, time, EF[1], dk, Nb[0], Coord_B::getJ(), Dipole);
            if(iWFDs)    
                PrintKFD(P0,Nk,Nb[0]+Nb[1],kpt,Unitary,icont);
        }

        //------ 1st step Runge-Kutta
        EF[0] = pulse1.E(time);
        EF[1] = .5*pulse2.E(time);
        Runge_Kutta_Df(P0,Pv,T,Nb,EF,pulse2.wl,Hamiltonian,Unitary,Dipole,GradientIndex, Weigths, Bvector,nk);
        Runge_Kutta_Ad(P0,P1,Pv,dt6,nk);
        Runge_Kutta_Ad(P0,P2,Pv,dt2,nk);


/*            ofstream Ptest;
            stringstream sn; sn << cnt;
            Ptest.open("test_"+ sn.str()+".txt");
            for(int ik=0; ik<nk; ik++)
            {
                Ptest << kpt[ik][0] << " " << kpt[ik][1] << " " << kpt[ik][2] << " ";
                for(int ic=0; ic< P0.n2(); ic++)
                    for(int jc=0; jc<P0.n2(); jc++)
                        Ptest << P2[ik][ic][jc].real() << " " << P2[ik][ic][jc].imag()<<" ";
                Ptest << endl;
            }   
*/

        //------ 2nd step Runge-Kutta
        time=it*dt + dt2;
        EF[0] = pulse1.E(time);
        EF[1] = .5*pulse2.E(time);
        Runge_Kutta_Df(P2,Pv,T,Nb,EF,pulse2.wl,Hamiltonian,Unitary,Dipole,GradientIndex, Weigths, Bvector,nk);
        Runge_Kutta_Ac(P1,Pv,dt3,nk);
        Runge_Kutta_Ad(P0,P2,Pv,dt2,nk);



        //------ 3rd step Runge-Kutta
        Runge_Kutta_Df(P2,Pv,T,Nb,EF,pulse2.wl,Hamiltonian,Unitary,Dipole,GradientIndex, Weigths, Bvector,nk);
        Runge_Kutta_Ac(P1,Pv,dt3,nk);
        Runge_Kutta_Ad(P0,P2,Pv,dt,nk);

        //------ 4th step Runge-Kutta
        time=(it+1)*dt;
        EF[0] = pulse1.E(time);
        EF[1] = .5*pulse2.E(time);
        Runge_Kutta_Df(P2,Pv,T,Nb,EF,pulse2.wl,Hamiltonian,Unitary,Dipole,GradientIndex, Weigths, Bvector, nk);
        Runge_Kutta_Ad(P1,P0,Pv,dt6,nk);

    }//end time evolution
    #include "Source_Main/end_program.cpp"

} //END MAIN
