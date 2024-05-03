
    /************************************************************************************************************************************************
    **************************************************         Definition of the        *************************************************************
    *************************************************   variables of the calculation     ************************************************************
    *************************************************************************************************************************************************
    *************************************************************************************************************************************************/
    if (rank_ == root_rank) system("mkdir Output");
    int Ncut=30;  // number of terms in Fourier series
    double dt = 0.1;
    int it_resolution = 10; // default
    int it;
    // counter to print energy dispersion
    int icont= -1;
    bool iIntraband=true;
    bool iWFDs=false;  //plots k distribution for Caterina code
    int wfd_resolution = 1;
    bool iTAbsK=false; int nTAk=0;  vec2d TAkpt; vec1i tagTAk;
    bool iTAbs= true; 
    bool iCurrent=false;
    bool Print_and_sort_k = true;
    bool Print_initial_H_P = true;
    bool Vectorization = true;
    string iMode="none"; //bool iTightBinding=false; //string to define the model

    Coulomb_parameters Coulomb_set; // struct with Coulomb parameters


    Coulomb_set.Coulomb_calc = false;// by default it false, if we don't initialize in in input file
    Coulomb_set.qTF = 0.0;
    Coulomb_set.epsilon_static = 1.0;
    Coulomb_set.G_distance = 50; // By default, will be multiplied on a sqrt(det(b)) where det(b) is a size of a BZ in cartesian

    Coulomb_set.Ncut = Ncut;
    Coulomb_set.test = 0.0;
    Coulomb_set.Wannie_basis = true; // by default in Wannier basis
    Coulomb_set.Diagonal_basis = false;
    Coulomb_set.Rytova_Keldysh = false;
    Coulomb_set.Print_new_band_dispersion = false;
    Coulomb_set.n_cond = 0;
    Coulomb_set.Sample_orientation = "011"; // by default in yz plane and is 011
    // Rytova-Keldysh polarizability parameter, by default
    Coulomb_set.r0 = 18.89726; //a.u. BN on quarz Henriques et al https://doi.org/10.1088/1361-648X/ab47b3
    Coulomb_set.labelInput = "_";
    Coulomb_set.Calculate = true; // By default we calculate caoefficients, but we can also read from file
    Coulomb_set.Read_from_files = false;


    
    string TBtype                = "CoreGraphene";
    string RefCoord;                                            //string to read the Wannier90 output
    vec2d a(3,3);                                    //primitive vectors of the system (r-space)
    vec2d b(3,3);                                    //primitive vectors of the system (k-space)
    vec1i Nk(3);                                                                  //number of k points in each direction
    int nk;
    bool kptread = false;
    vec1i resolution(3); resolution[0] = 1; resolution[1] = 4; resolution[2] = 4; //use for the resolution of the movie
    vec1d dk(3);                                                                       //delta k in each direction. NB: in crystal coordinates dk=1/Nk
    vec2d kpt;                                                                   //this is a vector to store the k points, k[0]=kx, k[1]=ky, k[2]=kz
    vec1i Nb(3); Nb.fill(0);                                                          //Nb[0]=Nch #core bands, Nb[1]=Nc #valence bands, Nb[2]=Nv #conduction bands
    double T1=0.;                                                                      //dissipation term for the diagonal terms of rho
    double T2=0.;                                                                      //dissipation term for the off-diagonal terms of rho
    double Tch=0.;                                                                     //dissipation term for the core states

    vector<Laser> Laser_pumps;                                                                  //IR pulse
    vec1d EF1(3),AF1(3); EF1.fill(0.); AF1.fill(0.);                                                              //electric field and potential vector of the IR pulse    
    Laser pulse2;                                                                   //X pulse
    vec1d EF2(3); EF2.fill(0.);
    double FermiE=0.;                                                               //Fermi energy of the system
    double DELAY=0.;                                                                //Delay in time between X pulse and IR pulse
    int Ncv;                                                                        //variable for total number of bands
    vec1d hopping;                            //vector in case we need extra information of the hopping parameters
    

    double t_fin = 0.0;

    methods_Diff_Eq   Diff_Eq; // structure wich methods we use to solve equation
    Diff_Eq.const_dt_evolution = true; // by default
    Diff_Eq.dynamical_dt_evolution = false; // by default
    Diff_Eq.epsStepAbs =  pow(10, -2);
    Diff_Eq.SolverOrder = 5;// by default
    Diff_Eq.Taylor= false;// by default
    Diff_Eq.Adams_Bashforth = false;
    Diff_Eq.PrintPopulation = false;// by default
    Diff_Eq.start_print_time = 0.0; // by default
    Diff_Eq.end_print_time = 100.0* time_au_fs;
    Diff_Eq.step_print = 20.0* time_au_fs;
    Diff_Eq.Gap_correction = 0; // by default
    vec1d dt_prev;
    double max_dP = 0.0; // for dynamical dt
    double max_P =0.0; // for dynamical dt
    double dtDynamic; // for dynamical dt
    bool notEnoughPrecision = true; // for dynamical dt
    int it_d; // for dynamical dt
    double time_RK;
    double n_cond=0;
    double P_cond_max;
    double P_cond_max_loc_mpi;
    bool print_matrices_in_kspace = false;

    complex<double> D_temporary;

    //for gradient calculation
    vector<vector<vector<int>>> GradientIndex;
    vector<double> Weigths;
    vector<vector<vector<double>>> Bvector;

    ofstream t_Ek;



    for(int i = 0; i< 3 ; i++)
    {
        for(int j=0; j<3; j++)
            if(i==j) {a[i][j] = 1.;   b[i][j] = 1.;}
            else {a[i][j] = 0.;   b[i][j] = 0.;}
    }    
    Coord_B::set_crys_to_cart(b);

    vec1x J1(3); //vector used in the current calculation to store the current in x,y and z
    vec1x J2(3); //vector used in the current calculation to store the current in x,y and z
    bool drop_current = false;
    vec4x GradientEnergy; //vector used to store gradient of energy in x,y and z

    double change_gap_constant = 0;