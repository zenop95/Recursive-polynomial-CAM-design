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
    int j, i, k, ii, flag1, flag2, flag3, order, metricFlag, pocType, N, lowThrust_flag, n_conj, n_man, m, dyn, gravOrd, missDistanceFlag, TPoCFlag;
	double mass, A_drag, Cd, A_srp, Cr, tca, Lsc, dt, musc, ctrlMax, mean_motion_p;
 
    ifstream nodes;
	nodes.open("./write_read/initial_state.dat");
        nodes >> N;         // Number of nodes
        nodes >> n_conj;    // Number of conjunctions
        nodes >> n_man;     // Number of control nodes
        nodes >> m;         // Number of DA variables per node
	nodes.close();

    AlgebraicMatrix<double> P(3,3), cov(9,n_conj), P_B3(3,3), P_B(2,2), r2e(3,3), toB(3,3), ctrlDum(3,n_man), xsdum(6,n_conj), rsDum(3,n_conj), vsDum(3,n_conj), xTca(6,n_conj), xsTca(6,n_conj);
    AlgebraicVector<double> xdum(6), x0(6), ctrl(3), ctrlRtn(3), t(N), HBR(n_conj), xRet(6), mean_motion_s(n_conj);
    AlgebraicVector<int>    canFire(N), isConj(N), isRet(N);
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
        Input >> gravOrd;            // Gravitational constant
        Input >> ctrlMax;            // Gravitational constant
        Input >> TPoCFlag;        // maneuver on miss distance
        Input >> missDistanceFlag;        // maneuver on miss distance
        Input >> mean_motion_p;        // mean motion primary
        for (k = 0; k < n_conj; k ++) {
            Input >> mean_motion_s[k];        // mean motion secondary
        }
        for (j = 0; j < 6; j ++) {
            Input >> xdum[j];
        }
        for (k = 0; k < n_conj; k ++) {
            for (j = 0; j < 6; j ++) {
                Input >> xsdum.at(j,k);
            }
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
            Input >> isRet[i];
        }
	Input.close();
    DA::init(1, 7);
    DA::setEps(1e-30);
    AlgebraicVector<DA> xp0(6), xs0(6), tcaNew(n_conj);

    if (dyn == 0) {
        if (gravOrd == 0) {
            x0 = RK78(6, xdum, {0.0,0.0,0.0}, 0.0, - t[0], keplerPropAcc, 1.0, Lsc); // Earth Orbit
        }
        else if (gravOrd == 2) {
            x0 = RK78(6, xdum, {0.0,0.0,0.0}, 0.0, - t[0], J2dynamics, 1.0, Lsc); // Earth Orbit
        }        
        else if (gravOrd == 4) {
            x0 = RK78(6, xdum, {0.0,0.0,0.0}, 0.0, - t[0], J2_J4dynamics, 1.0, Lsc); // Earth Orbit
        }
        else {throw std::runtime_error("The gravOrd flag should be 0, 2, or 4");}
    }
    else if (dyn == 1) {
        x0 = RK78(6, xdum, {0.0,0.0,0.0}, 0.0, - t[0], CR3BPsyn, musc, Lsc); // Earth Orbit
    }    
    else {throw std::runtime_error("The dynamics flag should be 0 for Earth orbit and 1 for Cislunar");}
    k  = 0;
    ii = 0;
    // Propagations at each maneuvering time
    for (i = 0; i < N-1; i ++) {
       if (canFire[i] == 1) {
            r2e     = cam::rtn2eci(cons(x0));
            for (j = 0; j < 3 ; j++) { 
                ctrlRtn[j] = ctrlDum.at(j,ii);
            }
            if (dyn + lowThrust_flag == 0) {
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
            if (gravOrd == 0) {
                x0 = RK78(6, x0, ctrl, tca - t[i], tca - t[i+1], keplerPropAcc, 1.0, Lsc);   // forward propagation to TCA
            }
            else if (gravOrd == 2) {
                x0 = RK78(6, x0, ctrl, tca - t[i], tca - t[i+1], J2dynamics, 1.0, Lsc);   // forward propagation to TCA
            }        
            else if (gravOrd == 4) {
                x0 = RK78(6, x0, ctrl, tca - t[i], tca - t[i+1], J2_J4dynamics, 1.0, Lsc);   // forward propagation to TCA
            }
            else {throw std::runtime_error("The gravOrd flag should be 0, 2, or 4");}
        }
        else {
            x0 = RK78(6, x0, ctrl, -t[i], -t[i+1], CR3BPsyn, musc, Lsc); // backpropagation from TCA        
        }
        if (isConj[i+1] == 1) {
            DA dt = 0.0 + DA(7);
            for (j = 0; j < 6; j ++) {
                xp0[j] = x0[j] + DA(j+1);
                xs0[j] = xsdum.at(j,k) + DA(j+1)*0;
            }   
            xp0  = KeplerProp(xp0, dt, 1.0);
            xs0 = KeplerProp(xs0, dt, 1.0);
            tcaNew[k] = findTCA(xp0 - xs0, 7);
            AlgebraicVector<DA> dx(7);
            for (int j = 0; j < 6; j++) {
                dx[j] = DA(j+1);}
            dx[6] = tcaNew[k];
            xp0 = xp0.eval(dx);
            xs0 = xs0.eval(dx);
            for (j = 0; j < 6 ; j ++) {
                xTca.at(j,k) = cons(xp0[j]);
                xsTca.at(j,k) = cons(xs0[j]);
            }        
            k ++;
        }
        if (isRet[i+1] == 1) {
            xRet = x0;
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
    for (k = 0; k < n_conj ; k++) {
        for (j = 0; j < 6 ; j++) {
        constPart  << xsTca.at(j,k) << endl;
        }
    }
    for (j = 0; j < 6 ; j++) {
        constPart  << xRet[j] << endl;
    }
    constPart.close();
    ofstream tcaOut;
    tcaOut.open("./write_read/tcaOut.dat");
    tcaOut << setprecision(16);
    for (k = 0; k < n_conj ; k++) {
        tcaOut  << cons(tcaNew[k]) << endl;
    }
    tcaOut.close();
    
}
