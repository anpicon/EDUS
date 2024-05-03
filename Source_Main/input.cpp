/************************************************************************************************************************************************
**************************************************              Reading             *************************************************************
*************************************************             input file             ************************************************************
*************************************************************************************************************************************************
*************************************************************************************************************************************************/


ifstream fp_input;
if ( argc != 2 )
{
    cout<<"usage: "<< argv[0] <<" <filename>" << endl;
}
else
{
    fp_input.open(argv[1]);
    if (!fp_input.is_open())
    {
        cout << "error opening file " << argv[1] << endl;
    }
}

printf("Starting read input...\n");
fp_input >> RefCoord;
Read_Input(fp_input, a, b,
    Nk, kpt, kptread,         //k space
    Nb,
    Laser_pumps,  pulse2, DELAY,// u1,  u2, gaussian1, gaussian2, sigma1, sigma2, DELAY,
    iMode, TBtype, dt,t_fin,
    iWFDs, wfd_resolution, Coulomb_set,  iCurrent, iTAbs, iTAbsK,
    T1,  T2,  Tch,
    nTAk, TAkpt, tagTAk,
           hopping, 
           FermiE,
           Diff_Eq,
           it_resolution,
           change_gap_constant,
           Vectorization
    );
printf("End read input...\n");




BZ_parameters(Nk,dk, Laser_pumps, pulse2,b,kpt,kptread,nTAk,TAkpt,tagTAk); //We fix the parameters for the Brillouin zone: coordinate transformation, k-grid definition, specify k points where the absorption is going to be calculated
Ncv=Nb[0]+Nb[1]+Nb[2];
int N_BZ = Nk[0]*Nk[1]*Nk[2];

if(iMode=="TB")
{
    TB_Model.init(TBtype, Nb[0], Ncv);
    TB_Model.set_hopping(hopping);
}
else if (iMode=="W") {
    int num_core_init_wannier = 0;
    WModel.init(RefCoord, num_core_init_wannier, Ncv );
}
else if (iMode=="CY") CYModel.init(a,b,Nb[0], Ncv);
else {
    printf("iMode = %s is incorrect \n", iMode.c_str());
    exit(1);
    }

double det_b, det_a;

if (Nk[2] < Nk[1]) Coulomb_set.Sample_orientation = "110"; 

if (Coulomb_set.Sample_orientation == "011"){ // sample yz
    det_b = abs(b[1][1]*b[2][2] - b[1][2]*b[2][1]); // determinant in plane for volume
    det_a = abs(a[1][1]*a[2][2] - a[1][2]*a[2][1]); // determinant in plane for volume
} else if (Coulomb_set.Sample_orientation == "110"){ // sample in xy
    det_b = abs(b[0][0]*b[1][1] - b[0][1]*b[1][0]); // determinant in plane for volume
    det_a = abs(a[0][0]*a[1][1] - a[0][1]*a[1][0]); // determinant in plane for volume
}


if (rank_ == 0){
    cout << "____________________________" << endl;
    cout << "unit cell area = " << det_a   << endl;
    cout << "Sample orientation = " << Coulomb_set.Sample_orientation   << endl;
    cout << "____________________________" << endl;
}
// exit(1);

CalculateWeigths(kpt, Weigths, Bvector, Nk);
//CalculateWeigths(kpt, Weigths, Bvector, Nk);










Print_Input(RefCoord,a, b,
        Nk, dk,
        Nb,
        Laser_pumps,  pulse2, DELAY, //u1,  u2, gaussian1, gaussian2, sigma1, sigma2, DELAY,
        iMode, TBtype, dt,
        iWFDs,  iCurrent, iTAbs, iTAbsK,
        T1,  T2,  Tch,
        nTAk, TAkpt, tagTAk,
        FermiE
        );

if(Diff_Eq.const_dt_evolution){
    cout << "time step dt is constant" << endl;
} else {
    cout << "time step dt is dynamical" << endl;
}

fp_input.close();
printf("Input file read.\n");

Coulomb_set.N_BZ_points_total=N_BZ;
Coulomb_set.N_shared.resize(Ncv);
Coulomb_set.N_exciton.resize(Ncv); 



vec2d Displacement_orb(Ncv, 3); // displacements for interaction terms


int TaylorOrder = Diff_Eq.SolverOrder;

if(Diff_Eq.Taylor){
    dt_prev.resize(TaylorOrder);
    Coulomb_set.Taylor_Delta_Pk_max.resize(TaylorOrder); // for Taylor solver
    Coulomb_set.Taylor_alpha_max.resize(TaylorOrder); 
    std::cout  << "   Taylot solver of oder   " << TaylorOrder << endl;
} else if(Diff_Eq.Adams_Bashforth){
    std::cout  << "   Adams-Bashforth solver of order " << Diff_Eq.SolverOrder << endl;
}else{
    // RK time evolution
    std::cout  << "   Runge Kutta solver " << endl;
}

MPI_Barrier(MPI_COMM_WORLD);

