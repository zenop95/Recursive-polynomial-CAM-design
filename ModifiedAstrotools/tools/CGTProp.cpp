#define _USE_MATH_DEFINES
#define WITH_ALGEBRAICMATRIX

#include <dace/dace.h>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <cspice/SpiceUsr.h>
#include "dynorb/DynClass.h"
#include "astro/AstroLibrary.h"
#include "astro/AstroRoutines.h"
#include "dynorb/AIDA.h"
#include "dynorb/AIDAwrappers.h"
#include <chrono>

using namespace std;
using namespace DACE;
using namespace std::chrono;

int main(void)
{
    int j, jj, kk, i, flag1, flag2, flag3, gravOrd, order, N, lowThrust_flag, n_conj, n_man, m;
	double mass, A_drag, Cd, A_srp, Cr, tca, Lsc, musc;

    long time1 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
  
    ifstream nodes, Input;
	nodes.open("./write_read/initial_state.dat");
        nodes >> N;         // Number of nodes
        nodes >> n_conj;    // Number of conjunctions
        nodes >> n_man;     // Number of control nodes
        nodes >> m;         // Number of DA variables per node
	nodes.close();
    // Initialize variable
    AlgebraicMatrix<double> ctrlDum(m,n_man), rsDum(3,n_conj), vsDum(3,n_conj);
    AlgebraicVector<double> xdum(6), t(N);
    AlgebraicVector<int>    canFire(N), isConj(N), isRet(N);
    // Write input from .dat
	Input.open("./write_read/initial_state.dat");
        Input >> N;               // Number of nodes
        Input >> n_conj;          // Number of conjunctions
        Input >> n_man;           // Number of control nodes
        Input >> m;               // Number of DA variables per node
        Input >> dyn;             // Dynamics model (0 Earth Orbit, 1 Cislunar)
        Input >> lowThrust_flag;  // Low thrust dynamics flag
        Input >> order;           // Expansion order
        Input >> pocType;         // PoC model (0 Alfriend, 1 Chan)
        Input >> tca;             // Ephemeris time at conjunction
        Input >> Lsc;             // Length scale
        Input >> musc;            // Gravitational constant
        for (j = 0; j < 6; j ++) {
            Input >> xdum[j];           // Write dummy variable for the primary's ECI state at the first conjunction from input
        }
        for (i = 0; i < n_man; i ++) {
            for (j = 0; j < m; j ++) {
                Input >> ctrlDum.at(j,i);   // Write dummy reference control from input ({0 0 0} if ballistic trajectory)
            }
        }
        for (i = 0; i < N; i ++) {
            Input >> t[i];                  // Write node times before conjunction (if negative it means that there are more than one encounters and nodes are needed after the first encounter)
            Input >> canFire[i];            // Define if node i is a firing node
            Input >> isRet[i];              // Define if node i is a return node
        }
    Input.close();
    // initialize execution time counter
    long time2 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    long timeSubtr = time2 - time1;
    // Initialize DA     
    DA::init(order, m*n_man);
    DA::setEps(1e-30);
    // Initialize AIDA    
    string kern = "./data/kernel.txt";
    ConstSpiceChar*FileKernel = kern.c_str();
    furnsh_c (FileKernel);
    string gravmodel  = "./data/gravmodels/egm2008";
    int AIDA_flags[3] = {flag1,flag2,flag3};
    double Bfactor    = Cd*A_drag/mass;
    double SRPC       = Cr*A_srp/mass;
    AIDAScaledDynamics<DA> aidaCartDyn(gravmodel, gravOrd, AIDA_flags, Bfactor, SRPC);

    // Initialize DA variables
    AlgebraicVector<DA> x0(6), xBall(6), xf(6), r(3), r_rel(2), v(3), rB(3), ctrl(3), xRet(6),; 
    // Define ballistic primary's position at first TCA 
    for (j = 0; j < 6 ; j++) {xBall[j] = xdum[j] + 0*DA(1);}
    // backpropagation from first TCA
    x0     = RK78Cis(6, xBall, {0.0*DA(1),0.0*DA(1),0.0*DA(1)}, 0.0, - t[0], CR3BPsyn, musc, 0.0); // Cislunar
    jj     = 0;
    kk     = 0;
    // Propagations at each of the N discretization nodes included in t
    for (i = 0; i < N-1; i ++) {
        // If the node is a firing node, include the control in the form of a DA vector variable
       if (canFire[i] == 1) {
            // RTN to ECI transformation used to have the output control in RTN coordinates
            r2e     = astro::rtn2eci(cons(x0));
            // Both direction and magnitude of the control are optimized
            for (j = 0; j < 3 ; j++) { 
                jj ++;
                ctrl[j] = ctrlDum.at(j,kk) + DA(jj);
            }
            kk ++;
            if (lowThrust_flag == 0) {
                for (j = 3; j < 6; j ++) {
                    x0[j] = x0[j] + ctrl[j-3];
                } 
                ctrl = {0.0*DA(1), 0.0*DA(1), 0.0*DA(1)};
            }
        }
        else {ctrl = {DA(1)*0, DA(1)*0, DA(1)*0};}

        x0 = RK78Cis(6, x0, ctrl, -t[i], -t[i+1], CR3BPsyn, musc, 0.0); // forward propagation to the next node

        if (isRet[i+1] == 1) {
            xRet = x0;
        }
    }

    // stability of the monodromy matrix (Cauchy-Green Tensor) after an orbital period
    ctrl = {0.0*DA(1), 0.0*DA(1), 0.0*DA(1)};
    x0 = RK78Cis(6, x0, ctrl, -t[end], -t[end] + 1, CR3BPsyn, musc, 0.0); // forward propagation to the next node
    AlgebraicMatrix<double> STM = astro::stmDace(x0, 3*n_man, 6);
    AlgebraicMatrix<double> CG = STM.transpose()*STM;

    time1   = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    //open the output files
    ofstream constPart;
    constPart.open("./write_read/constPart.dat");
    constPart << setprecision(18);
    constPart << CG << endl;
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