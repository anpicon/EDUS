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
    pulse1,  pulse2, DELAY,// u1,  u2, gaussian1, gaussian2, sigma1, sigma2, DELAY,
    iMode, TBtype, dt,t_fin,
    iWFDs, wfd_resolution, Coulomb_set,  iCurrent, iTAbs, iTAbsK,
    T1,  T2,  Tch,
    nTAk, TAkpt, tagTAk,
           hopping, 
           FermiE,
           Diff_Eq,
           it_resolution
    );
printf("End read input...\n");




BZ_parameters(Nk,dk,pulse1,pulse2,b,kpt,kptread,nTAk,TAkpt,tagTAk); //We fix the parameters for the Brillouin zone: coordinate transformation, k-grid definition, specify k points where the absorption is going to be calculated
Ncv=Nb[0]+Nb[1]+Nb[2];
int N_BZ = Nk[0]*Nk[1]*Nk[2];

if(iMode=="TB")
{
    TB_Model.init(TBtype, Nb[0], Ncv);
    TB_Model.set_hopping(hopping);
}
else if (iMode=="W") WModel.init(RefCoord, Nb[0], Ncv );
else if (iMode=="CY") CYModel.init(a,b,Nb[0], Ncv);
else {
    printf("iMode = %s is incorrect \n", iMode.c_str());
    exit(1);
    }






CalculateWeigths_new(kpt, Weigths, Bvector, Nk);
//CalculateWeigths(kpt, Weigths, Bvector, Nk);










Print_Input(RefCoord,a, b,
        Nk, dk,
        Nb,
        pulse1,  pulse2, DELAY, //u1,  u2, gaussian1, gaussian2, sigma1, sigma2, DELAY,
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


if (Coulomb_set.Coulomb_calc){   
    
    Ncut = Coulomb_set.Ncut;
    Coulomb_set.A.resize(Ncut*Ncut); // !!!! allocate memory for coefficients
    Coulomb_set.D.resize(Ncut*Ncut); // !!!! allocate memory for coefficients
    Coulomb_set.X_cc.resize(Ncut*Ncut*Ncv*Ncv);
    Coulomb_set.X_sc.resize(Ncut*Ncut*Ncv*Ncv);
    Coulomb_set.X_cs.resize(Ncut*Ncut*Ncv*Ncv);
    Coulomb_set.X_ss.resize(Ncut*Ncut*Ncv*Ncv);
    Coulomb_set.X_cc0.resize(Ncut*Ncut*Ncv*Ncv);
    Coulomb_set.X_sc0.resize(Ncut*Ncut*Ncv*Ncv);
    Coulomb_set.X_cs0.resize(Ncut*Ncut*Ncv*Ncv);
    Coulomb_set.X_ss0.resize(Ncut*Ncut*Ncv*Ncv);
    Coulomb_set.X_ss0.fill(0.);
    Coulomb_set.X_cc0.fill(0.);
    Coulomb_set.X_sc0.fill(0.);
    Coulomb_set.X_cs0.fill(0.);
    Coulomb_set.f_cc.resize(Ncv*Ncv);
    Coulomb_set.f_sc.resize(Ncv*Ncv);
    Coulomb_set.f_cs.resize(Ncv*Ncv);
    Coulomb_set.f_ss.resize(Ncv*Ncv);
    Coulomb_set.f_ss.fill(0.);
    Coulomb_set.f_cc.fill(0.);
    Coulomb_set.f_sc.fill(0.);
    Coulomb_set.f_cs.fill(0.);




    #include "Fourier_Coulomb.cpp"
    std::cout  << "   Coulomb_set.qTF   " << Coulomb_set.qTF << endl;
    std::cout  << "   Coulomb_set.epsilon_static   " << Coulomb_set.epsilon_static << endl;
    std::cout  << "   Coulomb_set.Coulomb_calc   " << Coulomb_set.Coulomb_calc << endl;
    if(Coulomb_set.Wannie_basis){
        std::cout  << "   Coulomb interaction defined in Wannier Basis   " << endl;
    }
    if(Coulomb_set.Diagonal_basis){
        std::cout  << "   Coulomb interaction defined in diagonal Basis   " << endl;
    }
    if(Coulomb_set.Rytova_Keldysh){
        std::cout  << "   Bare Coulomb in Rytova-Keldysh form   " << endl;

    }
}

int TaylorOrder = Diff_Eq.TaylorOrder;

if(Diff_Eq.Taylor){
    dt_prev.resize(TaylorOrder);
    std::cout  << "   Taylot solver of oder   " << TaylorOrder << endl;
} else{
    // RK time evolution
    std::cout  << "   Runge Kutta solver " << endl;
}

