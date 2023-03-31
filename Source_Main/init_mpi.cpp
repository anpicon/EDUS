int ierr;
MPI_Status status;

//######### Initializing MPI  ###########
int required = MPI_THREAD_SINGLE;
int provided;
//ierr = MPI_Init ( &argc, &argv );
ierr = MPI_Init_thread ( &argc, &argv, required, &provided );
if ( ierr != 0 )
{
    cout << "\n";
    cout << "BONES - Fatal error!\n";
    cout << "MPI_Init returned ierr = " << ierr << "\n";
    return 0;
}

ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &rank_ );       //Determine this process's rank_.
ierr = MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );  //Determine the number of available processes.

    
int root_rank = 0;

omp_set_num_threads( omp_get_max_threads() );

cout << "omp_get_max_threads()" << omp_get_max_threads() << endl;     