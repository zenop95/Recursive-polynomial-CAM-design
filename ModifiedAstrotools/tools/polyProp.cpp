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
    int j, jj, kk, i, k, vv, flag1, flag2, flag3, gravOrd, order, pocType, N, lowThrust_flag, n_conj, n_man, m, dyn;
	double mass, A_drag, Cd, A_srp, Cr, tca, Lsc, musc;

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
  
    ifstream nodes, Input;
	nodes.open("./write_read/initial_state.dat");
        nodes >> N;         // Number of nodes
        nodes >> n_conj;    // Number of conjunctions
        nodes >> n_man;     // Number of control nodes
        nodes >> m;         // Number of DA variables per node
	nodes.close();
    // Initialize variable
    AlgebraicMatrix<double> P(3,3), cov(9,n_conj), P_B3(3,3), P_B(2,2), r2e(3,3), toB(3,3), ctrlDum(m,n_man), rsDum(3,n_conj), vsDum(3,n_conj), directions(3,n_man);
    AlgebraicVector<double> xdum(6), metricMap(3), t(N), vs(3), rs(3), rsB(3), HBR(n_conj), magnitude(n_man), rRef(3), vRef(3);
    AlgebraicVector<int>    canFire(N), isConj(N), isRet(N), constraintFlags(4);
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
        for (k = 0; k < n_conj; k ++) {
            Input >> HBR[k];            // Combined Hard Body radialius for each conjunction
        }
        for (k = 0; k < n_conj; k ++) {
            for (j = 0; j < 3; j ++) {
                Input >> rsDum.at(j,k); // Write secondary's ECI position at each conjunction from input
            }
            for (j = 0; j < 3; j ++) {
                Input >> vsDum.at(j,k); // Write secondary's ECI velocity at each conjunction from input
            }   
        }
        for (j = 0; j < 3; j ++) {
                Input >> rRef[j]; // Write reference ECI position for return from input
            }
        for (j = 0; j < 3; j ++) {
                Input >> vRef[j]; // Write reference ECI velocity for return from input
            }
        for (k = 0; k < n_conj; k ++) {
            for (j = 0; j < 9; j ++) {
                Input >> cov.at(j,k);   // Write dummy variable for the combined covariance in ECI at each conjunction from input
            }
        }
        for (i = 0; i < n_man; i ++) {
                Input >> magnitude[i];    // Write maneuver magnitude if the magnitude is fixed
        }   
        for (i = 0; i < n_man; i ++) {
            for (j = 0; j < 3; j ++) {
                Input >> directions.at(j,i); // Write maneuver direction in RTN if the direction is fixed
            }
        }
        for (i = 0; i < 4; i ++) {
            Input >> constraintFlags[i];
        }
        for (i = 0; i < n_man; i ++) {
            for (j = 0; j < m; j ++) {
                Input >> ctrlDum.at(j,i);   // Write dummy reference control from input ({0 0 0} if ballistic trajectory)
            }
        }
        for (i = 0; i < N; i ++) {
            Input >> t[i];                  // Write node times before conjunction (if negative it means that there are more than one encounters and nodes are needed after the first encounter)
            Input >> canFire[i];            // Define if node i is a firing node
            Input >> isConj[i];             // Define if node i is a conjunction node
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
    AlgebraicVector<DA> x0(6), xBall(6), xf(6), r(3), r_rel(2), v(3), rB(3), ctrlRtn(3), ctrl(3), rf(3), xRet(6); 
    AlgebraicMatrix<DA> xTca(6,n_conj);
    DA metric, poc_tot, alpha, beta;
    // Define ballistic primary's position at first TCA 
    for (j = 0; j < 6 ; j++) {xBall[j] = xdum[j] + 0*DA(1);}
    // backpropagation from first TCA
    if (dyn == 0) {
        x0     = RK78Sc(6, xBall, {0.0*DA(1),0.0*DA(1),0.0*DA(1)}, tca, tca - t[0], 1.0, Lsc, 0, gravOrd, aidaCartDyn); // Earth Orbit
    }
    else if (dyn == 1) {
        x0     = RK78Cis(6, xBall, {0.0*DA(1),0.0*DA(1),0.0*DA(1)}, 0.0, - t[0], CR3BPsyn, musc, 0.0); // Cislunar
    }
    else {throw std::runtime_error("The dynamics flag should be 0 for Earth orbit and 1 for Cislunar");}
    jj     = 0;
    k      = 0;
    kk     = 0;
    // Propagations at each of the N discretization nodes included in t
    for (i = 0; i < N-1; i ++) {
        // If the node is a firing node, include the control in the form of a DA vector variable
       if (canFire[i] == 1) {
            // RTN to ECI transformation used to have the output control in RTN coordinates
            r2e     = astro::rtn2eci(cons(x0));
            // Both direction and magnitude of the control are optimized
            if (m == 3) {
                for (j = 0; j < 3 ; j++) { 
                    jj ++;
                    ctrlRtn[j] = ctrlDum.at(j,kk) + DA(jj);
                }
            }
            // Only the direction of the control is optimized
            else if (m == 2) {
                // alpha = DA(2*kk+1);
                // beta  = DA(2*kk+2);
                // ctrlRtn[0] = magnitude[kk]*cos(alpha)*sin(beta);
                // ctrlRtn[1] = magnitude[kk]*cos(alpha)*cos(beta);
                // ctrlRtn[2] = magnitude[kk]*sin(alpha);
                ctrlRtn[0] = DA(2*kk+1);
                ctrlRtn[1] = DA(2*kk+2);
                ctrlRtn[2] = sqrt(magnitude[kk]*magnitude[kk] - ctrlRtn[0]*ctrlRtn[0] - ctrlRtn[1]*ctrlRtn[1]);
            }
            // Only the magnitude of the control is optimized
            else if (m == 1) {
                for (j = 0; j < 3 ; j++) {
                    ctrlRtn[j] = (ctrlDum.at(0,kk) + DA(kk+1))*directions.at(j,kk);
                }
            }
            else {
                throw std::runtime_error("The number of DA variables per node must be in the interval [1,3]");
            }
            kk ++;
            // In Cislunar optimization do not tranform to RTN because we are in the synodic frame
            if (dyn == 0) {
                ctrl = r2e*ctrlRtn;
            }
            else {
                ctrl = ctrlRtn;
            }
            // If impulsive, the control is added to the velocity part of the state and the acceleration is null
            if (lowThrust_flag == 0) {
                for (j = 3; j < 6; j ++) {
                    x0[j] = x0[j] + ctrl[j-3];
                } 
                ctrl = {0.0*DA(1), 0.0*DA(1), 0.0*DA(1)};
            }
        }
        // If the node is a ballistic node, do no include the control
        else {ctrl = {DA(1)*0, DA(1)*0, DA(1)*0};}
        if (dyn == 0) {
            x0 = RK78Sc(6, x0, ctrl, tca - t[i], tca - t[i+1], 1.0, Lsc, 0, gravOrd, aidaCartDyn);   // forward propagation to the next node
        }
        else {
            x0 = RK78Cis(6, x0, ctrl, -t[i], -t[i+1], CR3BPsyn, musc, 0.0); // forward propagation to the next node
        }
        // If the next node is a conjunction node, save the state in a DA variable
        if (isConj[i+1] == 1) {
            for (j = 0; j < 6 ; j ++) {
                xTca.at(j,k) = x0[j];
            }        
            k = k + 1;
        }
        else if (isRet[i+1] == 1) {
            xRet = x0;
        }
    }

if (constraintFlags[0] == 1) {
    // Compute the total PoC resulting from the multiple conjunctions
    DA noCollisions = 1.0;
    for (k = 0; k < n_conj; k ++) {
        // Expanded position and velocity of the primary at conjunction k
        r[0] = xTca.at(0,k);  r[1] = xTca.at(1,k);  r[2] = xTca.at(2,k);   
        v[0] = xTca.at(3,k);  v[1] = xTca.at(4,k);  v[2] = xTca.at(5,k);
        vv = 0;
        // Build the ECI combined covariance matrix and the secondary's ECI position and velocity from dummy variables
        for (i = 0; i < 3 ; i ++) {
            for (j = 0; j < 3 ; j ++) {
                P.at(i,j)  = cov.at(vv,k);
                vv = vv + 1;
            }
            vs[i] = vsDum.at(i,k);
            rs[i] = rsDum.at(i,k);
        }
        // B-plane transformations
        toB = astro::Bplane(cons(v),vs); // DCM from ECI to B-plane
        rB  = toB*r;                     // Primary position in the B-plane (3D)
        rsB = toB*rs;                    // Secondary position in the B-plane (3D)
        P_B3 = toB*P*toB.transpose();    // Combined covariance in the B-plane (3D)

        // Relative position in the B-plane (2D)
        r_rel[0]    = rB[0] - rsB[0];   r_rel[1]    = rB[2] - rsB[2]; 
        // Combined covariance in the B-plane (2D)
        P_B.at(0,0) = P_B3.at(0,0);     P_B.at(0,1) = P_B3.at(0,2);
        P_B.at(1,0) = P_B3.at(2,0);     P_B.at(1,1) = P_B3.at(2,2);
        vv = 0;
        // Compute PoC for the single conjunction, according to the required model
        if (pocType == 0) {
            metric = astro::ConstPoC(r_rel,P_B,HBR[k]);}
        else if (pocType == 1) {
            metric = astro::ChanPoC(r_rel,P_B,HBR[k],3);}
        else {
            throw std::runtime_error("the metric flag must be in the interval [1,3] and the PoC type must be in the interval [0,1]");}
        // Probability of no collision
        noCollisions = noCollisions*(1.0 - metric);
    }

    // Final PoC comprehensive of all the conjunctions
    poc_tot = log10(1.0 - noCollisions);
}

DA tan, radial;
if (constraintFlags[1] == 1 || constraintFlags[2] == 1) {
    // Return constraint (scalar as tangential distance from other spacecraft GRACE)
    AlgebraicVector<DA> rRet(3), vRet(3), distRel(3);
    for (j = 0; j < 3 ; j ++) {
        rRet[j] = xRet[j];
        vRet[j] = xRet[j+3];
    }      
    r2e     = astro::rtn2eci(cons(xRet));
    distRel = r2e.transpose()*(rRet-rRef);
    radial  = distRel[0]; // radial displacement with respect to reference
    tan     = distRel[1]; // tangential displacement with respect to reference
}

// stability of the monodromy matrix (Cauchy-Green Tensor) after an orbital period
// if (constraintFlags[4] == 1 || dyn == 1) {
//     ctrl = {0.0*DA(1), 0.0*DA(1), 0.0*DA(1)};
//     x0 = RK78Cis(6, x0, ctrl, -t[end], -t[end] +     1, CR3BPsyn, 0.012150668, 0.0); // forward propagation to the next node
//     AlgebraicMatrix<double> STM = astro::stmDace(x0, 3, 6);
//     AlgebraicMatrix<double> CG = STM.transpose()*STM;
// }

    time1   = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    //open the output files
    ofstream constPart, constraints;
    constPart.open("./write_read/constPart.dat");
    constPart << setprecision(18);
    constraints.open("./write_read/constraints.dat");
    constraints << setprecision(18);

    // write TCA states in ECI coordinates into output
    for (k = 0; k < n_conj ; k++) {
        for (j = 0; j < 6 ; j++) {
        constPart  << cons(xTca.at(j,k)) << endl;
        }
    }

    // write the DA expansion of PoC in output
    if (constraintFlags[0] == 1) {
        constPart   << cons(poc_tot)  << endl;
        constraints << poc_tot << endl;
    }

    // write the DA expansion of return in output
    if (constraintFlags[1] == 1) {
        constPart   << cons(tan)      << endl;
        constraints << tan    << endl;
    }

    // write the DA expansion of return in output
    if (constraintFlags[2] == 1) {
        constPart   << cons(radial) << endl;
        constraints << radial       << endl;
    }

    // write the DA expansion of return in output
    if (constraintFlags[3] == 1) {
        for (j = 0; j < 6; j ++) {
        constPart   << cons(xRet[j])  << endl;
        constraints << xRet[j]        << endl;
        }
    }
    if (constraintFlags[4] == 1) {
        
    }

    // Do not consider writing time when calculating execution time
    time2 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    timeSubtr = timeSubtr + time2 - time1;

    ofstream timeOut;
    timeOut.open("./write_read/timeOut.dat");
    timeOut << setprecision(16);
    timeOut << timeSubtr << endl;
    timeOut.close();
    constPart.close();
    constraints.close();
}