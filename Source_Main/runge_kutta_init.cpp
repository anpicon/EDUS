Initial_Population( P0, nk, kpt, Nk[0]*Nk[1]*Nk[2], Unitary, FermiE, Nb[0]+Nb[1]);
//Initial_Population( P0, nk, Nk[0]*Nk[1]*Nk[2], Unitary, FermiE, Nb[0]+Nb[1]);

double dt6=dt/6.0;  
double dt3=dt/3.0;  
double dt2=dt/2.0;



double ncyc = 8;

time_t time1 = time(0);
tm* now = localtime(&time1);
int ihour = now->tm_hour;
int iday = now->tm_mday;
int imin = now-> tm_min;
int isec = now-> tm_sec;

printf("\n\n\nCalculation started day %4d/%02d/%02d at %02d.%02d.%02d \n", (now->tm_year + 1900), (now->tm_mon + 1), (now->tm_mday), (now->tm_hour), (now->tm_min), (now->tm_sec));

clock_t clock1,clock2; clock1=clock();


printf("\n*****************************************************\n");
printf("*            Time evolution with TDSE               *\n");
printf("*****************************************************\n");    //for(int it=0;(it)<=(10);it++)

double tf;
if(Laser_pumps[0].gaussian) tf = ncyc*Laser_pumps[0].sigma;
else          tf = Laser_pumps[0].ncycle *Laser_pumps[0].Period;


int nstep = int(tf/dt) ;
int iti = 0;
int itfi=nstep;

int cnt=0;
if (t_fin > 0){
    itfi = t_fin / dt; 
} else {
    t_fin = itfi*dt;
}
double dt0 = dt;




if(kptread)     itfi = int(tf/dt) + int(100*time_fs_au/dt);

printf("Number of steps: %3d \nFinal time: %5.2f\n", nstep, tf);
vector<vec1d> EF(2);
for(int i=0; i<2; i++) {EF[i].resize(3); EF[i].fill(0);}



// number of points in BZ
vec1d integrWeight(nk);

get2DPrismInterrationCoeff_MPI(Nk[1], kpt, integrWeight);

MPI_Barrier(MPI_COMM_WORLD);
for (int i = 0; i < p; i++){
	MPI_Barrier(MPI_COMM_WORLD);
    // Cheking kpt borders
    if (my_rank == i){
    	cout << " ******************** "  << endl;
        cout << " my_rank " <<  my_rank << endl;
        cout << " Nk[1] " <<  Nk[1] << endl;
        cout << " kpt size " << kpt.n1()   << endl; 
        cout << " nk " << nk  << endl << endl << endl; 
    }

	MPI_Barrier(MPI_COMM_WORLD);
}
MPI_Barrier(MPI_COMM_WORLD);

Coulomb_set.Ekdia.resize(nk,  Ncv); // array for diagonal energy


// struct with trigonometric functions
// trig_coefficients trig_k;




// MPI_Barrier(MPI_COMM_WORLD);
// for (int i = 0; i < num_procs; i++){
//   MPI_Barrier(MPI_COMM_WORLD);
//     // Cheking kpt borders
//     if (rank_ == i){
//         cout << " ******************** "  << endl;
//         cout << " my_rank " <<  rank_ << endl;
//         cout << " Coulomb_set.X_cc0[0] " << endl; 

//         for (int mn = 0; mn < Ncut*Ncut*Ncv*Ncv; mn += 5){
//             cout << " mn " << mn << "  :" << Coulomb_set.X_cc0[mn] << endl; 

//         }
           
//     }

//   MPI_Barrier(MPI_COMM_WORLD);
// }
// MPI_Barrier(MPI_COMM_WORLD);
