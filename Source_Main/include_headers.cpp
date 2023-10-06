#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <complex>
#include <ctime>

#include <vector>
#include <string.h>
//#include <map>
//#include <utility>
#include <memory>



using namespace std;

#include "BoostArrays.h"
#include "Constants.h"
#include "typedef.h"

#define MKL_Complex16 complexd
#include "mkl.h"

#include "Coordinate.h" 
#include "Laser_new.h"
#include "crystal.h"

#include "MPIcommunication.h"
#include "Coulomb.h"
#include "integrationWeight.h"
#include "Interpolation.h"
// #include "cBWE_dynamics.h"

#include "ReadInput.h"
#include "TightBinding.h"
#include "Wannier_old.h"
#include "Crystal_interface.h"
#include "HUD.h"
#include "InitParameters.h"
#include "RungeKutta.h"
#include "CalculateWeigths.h"
#include "CalculateWeigths_new.h"
#include "InitialPopulation.h"
#include "Observables.h"
#include "Taylor_solver.h"


// // matplotlib to draw results
// // // C++ library from here https://github.com/lava/matplotlib-cpp
// #include "matplotlibcpp.h"
// namespace plt = matplotlibcpp;
// #include "printing_matplotlib.h"
