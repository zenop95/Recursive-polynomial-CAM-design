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
    int j, i, k, ii, flag1, flag2, flag3, gravOrd, order, metricFlag, pocType, N, lowThrust_flag, n_conj, n_man, m, dyn;
	double mass, A_drag, Cd, A_srp, Cr, tca, Lsc, dt, musc;

    long time1 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();

    // AIDA initialization input file
    ifstream aidaParams;
	aidaParams.open("./write_read/AIDA_init.dat");
        aidaParams >> flag1;	     // Luni-solar model	
        aidaParams >> flag2;	     // Atmosphere model	
        aidaParams >> flag3;	     // SRP model
        aidaParams >> gravOrd;	     //	Maximum order of the gravitational potential considered
        aidaParams >> mass;	         //	Mass of the spacecraft
        aidaParams >> A_drag;	     // Equivalent drag area of the spacecraft	 
        aidaParams >> Cd;            // Aerodynmic coefficient of the spacecraft
        aidaParams >> A_srp;         // Equivalent SRP area of the spacecraft	
        aidaParams >> Cr;            // SRP coefficient of the spacecraft	
	aidaParams.close();

    ifstream nodes;
	nodes.open("./write_read/initial_state.dat");
        nodes >> N;         // Number of nodes
        nodes >> n_conj;    // Number of conjunctions
        nodes >> n_man;     // Number of control nodes
        nodes >> m;         // Number of DA variables per node
	nodes.close();
    AlgebraicMatrix<double> P(3,3), cov(9,n_conj), P_B3(3,3), P_B(2,2), r2e(3,3), toB(3,3), ctrlDum(3,n_man), rsDum(3,n_conj), vsDum(3,n_conj), xTca(6,n_conj);
    AlgebraicVector<double> xdum(6), x0(6), ctrl(3), ctrlRtn(3), t(N), HBR(n_conj);
    AlgebraicVector<int>    canFire(N), isConj(N);
    ifstream Input;
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
            Input >> xdum[j];
        }
        for (i = 0; i < n_man; i ++) {
            for (j = 0; j < 3; j ++) {
                Input >> ctrlDum.at(j,i);
            }
        }
        for (i = 0; i < N; i ++) {
            Input >> t[i];
            Input >> canFire[i];
            Input >> isConj[i];
        }
	Input.close();

    long time2 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    long timeSubtr = time2 - time1;
    // AIDA initialization    
    string kern = "./data/kernel.txt";
    ConstSpiceChar*FileKernel = kern.c_str();
    furnsh_c (FileKernel);
    string gravmodel  = "./data/gravmodels/egm2008";
    int AIDA_flags[3] = {flag1,flag2,flag3};
    double Bfactor    = Cd*A_drag/mass;
    double SRPC       = Cr*A_srp/mass;
    AIDAScaledDynamics<double> aidaCartDyn(gravmodel, gravOrd, AIDA_flags, Bfactor, SRPC);

    if (dyn == 0) {
        x0     = RK78Sc(6, xdum, {0.0,0.0,0.0}, tca, tca - t[0], 1.0, Lsc, musc, gravOrd, aidaCartDyn); // backpropagation from first TCA (Earth orbit)
    }
    else if (dyn == 1) {
        x0     = RK78Cis(6, xdum, {0.0,0.0,0.0}, 0.0, - t[0], CR3BPsyn, 0.012150668, 0.0); // backpropagation from first TCA (Cislunar)
    }    
    else {throw std::runtime_error("The dynamics flag should be 0 for Earth orbit and 1 for Cislunar");}
    k  = 0;
    ii = 0;
    // Propagations at each maneuvering time
    for (i = 0; i < N-1; i ++) {
       if (canFire[i] == 1) {
            r2e     = astro::rtn2eci(cons(x0));
            for (j = 0; j < 3 ; j++) { 
                ctrlRtn[j] = ctrlDum.at(j,ii);
            }
            if (dyn == 0) {
                ctrl = r2e*ctrlRtn;
            }
            else {
                ctrl = ctrlRtn;
            }
            ii ++;
        }
        else {ctrl = {0.0, 0.0, 0.0};}
        if (lowThrust_flag == 0) {
            for (j = 3; j < 6; j ++) {
                x0[j] = x0[j] + ctrl[j-3];
            } 
            ctrl = {0.0, 0.0, 0.0};
        }
        if (dyn == 0) {
            x0 = RK78Sc(6, x0, ctrl, tca - t[i], tca - t[i+1], 1.0, Lsc, 0, gravOrd, aidaCartDyn);   // forward propagation to TCA
        }
        else {
            x0 = RK78Cis(6, x0, ctrl, -t[i], -t[i+1], CR3BPsyn, 0.012150668, 0.0); // backpropagation from TCA        
        }
        if (isConj[i+1] == 1) {
            for (j = 0; j < 6 ; j ++) {
                xTca.at(j,k) = x0[j];
            }        
            k ++;
        }
    }
        
    //open the output files
    ofstream constPart;
    constPart.open("./write_read/constPart.dat");
    constPart << setprecision(16);
    for (k = 0; k < n_conj ; k++) {
        for (j = 0; j < 6 ; j++) {
        constPart  << xTca.at(j,k) << endl;
        }
    }
    constPart.close();

    ofstream timeOut;
    timeOut.open("./write_read/timeOut.dat");
    timeOut << setprecision(16);
    timeOut.close();
}