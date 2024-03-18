typedef complex<double>                           complexd;
typedef multivec1D<int>                           vec1i;
typedef multivec2D<int>                           vec2i;
typedef multivec3D<int>                           vec3i;
typedef multivec4D<int>                           vec4i;
typedef multivec5D<int>                           vec5i;
typedef multivec6D<int>                           vec6i;
typedef multivec1D<double>                        vec1d;
typedef multivec2D<double>                        vec2d;
typedef multivec3D<double>                        vec3d;
typedef multivec4D<double>                        vec4d;
typedef multivec5D<double>                        vec5d;
typedef multivec6D<double>                        vec6d;
typedef multivec1D<complexd>                      vec1x;
typedef multivec2D<complexd>                      vec2x;
typedef multivec3D<complexd>                      vec3x;
typedef multivec4D<complexd>                      vec4x;
typedef multivec5D<complexd>                      vec5x;
typedef multivec6D<complexd>                      vec6x;

vec1i    vecE1i; 
vec2i    vecE2i;
vec3i    vecE3i;
vec4i    vecE4i;
vec5i    vecE5i;
vec6i    vecE6i;
vec1d    vecE1d;
vec2d    vecE2d;
vec3d    vecE3d;
vec4d    vecE4d;
vec5d    vecE5d;
vec6d    vecE6d;
vec1x    vecE1x;
vec2x    vecE2x;
vec3x    vecE3x;
vec4x    vecE4x;
vec5x    vecE5x;
vec6x    vecE6x;


// trigonometric coefficients to calculate sums over k vectors
struct trig_coefficients{
  vec1d cos_mkx, sin_mkx, cos_nky, sin_nky;
  string Sample_orientation; // by default in yz plane and is 011
};

// structure contains which solving methods we have
struct methods_Diff_Eq{
  bool const_dt_evolution, dynamical_dt_evolution, Taylor;
  double epsStepAbs, Gap_correction;
  int TaylorOrder;
  bool PrintPopulation;
  double start_print_time, end_print_time, step_print;
  vec1i k_index;
  int Nk_total_All_Grid, Nk_0, Nk_1, Nk_2, Nk_node_MPI;
};


// structure contains private for each thread copies of variables
struct Private_omp_parameters{
  int thr_id;
  int thr_total;
  int lenght_k;
  int Nk_total_MPI, Nk_node_MPI;
  // P_diag - population in basis where equilibrium energy is diagonal
  // P_eigen -  population in basis where non-equlibrium hamiltonian is diagonal, 
  // different from P_diag only if we have Coulomb term
  vec3x P0, P_Wannier_0, P_Bloch_0, P1, P2, Pv, P_diag, P_eigen;

  vec2x P_0_dyn;

  //initial population
  vec3x P_static;

  // population at previous step
  vec3x P_W_prev;
  vec3x P_dia_prev;
  vec3x P2_dia;
  //


  vec3x P_grad;
  vec2d k;
  vec1d integrWeight;
  int end_count;
  int begin_count;
  std::vector<bool> id_masters;
  vec3x Hk, Uk, hUk;
  vec4x Dk;
  vec2x Mk;
  vec1x Ak;
  vec3x Hk_renorm;
  vec1i T_dissip_index_0, T_dissip_index_1; 
  int n_diss_terms;
  bool Vectorization;
};




double round_to_digits(double value, int digits)
{
    if (value == 0.0) // otherwise it will return 'nan' due to the log10() of zero
        return 0.0;

    double factor = pow(10.0, digits - ceil(log10(fabs(value))));
    return round(value * factor) / factor;   
}

double floor_to_digits(double value, int digits)
{
    if (value == 0.0) // otherwise it will return 'nan' due to the log10() of zero
        return 0.0;

    double factor = pow(10.0, digits - ceil(log10(fabs(value))));
    return floor(value * factor) / factor;   
}


// vectorisation function:
#pragma omp declare simd
void Ak_P_Uk_fill_vector_function(vec1x& Ak, vec3x& P, vec3x& Uk,
    int ik_P, int ik_U,   int ii, int jc, int jj){
    Ak[jj] = P[ik_P][ii][jj] * Uk[ik_U][jc][jj];
}


// vectorisation function:
#pragma omp declare simd
void Ak_Uk_Mk_fill_vector_function(vec1x& Ak, vec2x& Mk, vec3x& Uk,
    int ik_U, int ic, int jc, int ii){
    Ak[ii] = conj(Uk[ik_U][ic][ii]) * Mk[jc][ii];
}


// vectorisation function:
#pragma omp declare simd
void X_d_Uk_fill_vector_function(vec2x& Mk, complex<double>& Xk_d, vec3x& Uk,
    int ik_U,   int ii, int jj, int jc){
    Mk[ii][jc] += Xk_d * Uk[ik_U][jj][jc];

}


// vectorisation function:
#pragma omp declare simd
void X_W_conjUk_fill_vector_function(vec2x& Mk, vec2x& Xk_W, complex<double>& conjUk,
    int ii, int ic, int jc){
    
    Xk_W[ic][jc] += conjUk * Mk[ii][jc];
}

