#define _USE_MATH_DEFINES
#define WITH_ALGEBRAICMATRIX

#include <dace/dace.h>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include "camRoutines.h"
#include <chrono>

using namespace std;
using namespace DACE;
using namespace std::chrono;
using namespace cam;

int main(void)
{
    int j, jj, k, i, flag1, flag2, flag3, gravOrd, order, N, lowThrust_flag, n_conj, n_man, m;
	double mass, A_drag, Cd, A_srp, Cr, tca, Lsc, musc, T_orbit;

    long time1 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
  
    ifstream nodes, Input;
	nodes.open("./write_read/initial_state.dat");
        nodes >> N;         // Number of nodes
        nodes >> n_conj;    // Number of conjunctions
        nodes >> n_man;     // Number of control nodes
        nodes >> m;         // Number of DA variables per node
	nodes.close();
    // Initialize variable
    AlgebraicVector<double> xdum(6), t(N);
    // Write input from .dat
	Input.open("./write_read/initial_state.dat");
        Input >> N;               // Number of nodes
        Input >> n_man;           // Number of control nodes
        Input >> m;               // Number of DA variables per node
        Input >> lowThrust_flag;  // Low thrust dynamics flag
        Input >> order;           // Expansion order
        Input >> tca;             // Ephemeris time at conjunction
        Input >> T_orbit;         // Period of the periodic orbit (non-dimensional)
        Input >> Lsc;             // Length scale
        Input >> musc;            // Gravitational constant
        for (j = 0; j < 6; j ++) {
            Input >> xdum[j];           // Write dummy variable for the primary's ECI state at the first conjunction from input
        }
        for (i = 0; i < N; i ++) {
            Input >> t[i];                  // Write node times before conjunction (if negative it means that there are more than one encounters and nodes are needed after the first encounter)
        }
    Input.close();
    // initialize execution time counter
    long time2 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    long timeSubtr = time2 - time1;
    // Initialize DA     
    DA::init(order, 6);
    DA::setEps(1e-30);

    // Initialize DA variables
    AlgebraicVector<DA> x0(6), xBall(6), xf(6), r(3), r_rel(2), v(3), rB(3), ctrl(3), xRet(6); 
    // Define ballistic primary's position at first TCA 
    for (j = 0; j < 6 ; j++) {xBall[j] = xdum[j] + 0*DA(1);}
    // backpropagation from first TCA
    x0 = RK78(6, xBall, {0.0*DA(1),0.0*DA(1),0.0*DA(1)}, 0.0, - t[0], CR3BPsyn, musc, Lsc); // Cislunar Orbit
    jj = 0;
    k  = 0;

    // stability of the monodromy matrix (Cauchy-Green Tensor) after an orbital period
    for (j = 0; j < 6 ; j++) {
        x0[j] = x0[j] + DA(j+1);
    }    
    ctrl = {0.0*DA(1), 0.0*DA(1), 0.0*DA(1)};
    x0     = RK78(6, x0, {0.0*DA(1),0.0*DA(1),0.0*DA(1)}, 0.0, T_orbit, CR3BPsyn, musc, Lsc); // Cislunar Orbit
    AlgebraicMatrix<double> STM = stmDace(x0, 6, 6);
    time1   = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    //open the output files
    ofstream constPart;
    constPart.open("./write_read/constPart.dat");
    constPart << setprecision(18);
    for (j = 0; j < 6 ; j++) {
    for (k = 0; k < 6 ; k++) {
        constPart << STM.at(j,k) << endl;
    }
    }
    // for (k = 0; k < 6 ; k++) {
    //     constPart << cons(x0[k]) << endl;
    // }
    constPart.close();

    // Do not consider writing time when calculating execution time
    time2 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    timeSubtr = timeSubtr + time2 - time1;

    ofstream timeOut;
    timeOut.open("./write_read/timeOut.dat");
    timeOut << setprecision(16);
    timeOut << timeSubtr << endl;
    timeOut.close();
}