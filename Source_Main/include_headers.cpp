#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <complex>
#include <ctime>
//#include <vector>
//#include <map>
//#include <utility>
#include <armadillo>
#include <memory>
#include <sys/stat.h>


using namespace std;
using namespace arma;

#include "BoostArrays.h"
#include "Constants.h"
#include "typedef.h"
#include "Coordinate.h" 
#include "Laser.h"
#include "crystal.h"

#include "MPIcommunication.h"
#include "Coulomb.h"
#include "integrationWeight.h"
#include "Interpolation.h"


#include "ReadInput.h"
#include "TightBinding.h"
#include "Wannier.h"
#include "Crystal_interface.h"
#include "HUD.h"
#include "InitParameters.h"
#include "RungeKutta.h"
#include "CalculateWeigths.h"
#include "CalculateWeigths_new.h"
#include "InitialPopulation.h"
#include "Taylor_solver.h"
#include "Orderly_Printing.h"

#include "Observables_MPI.h"
#include "Organize_kspaceMPI.h"

