#define _USE_MATH_DEFINES
#define WITH_ALGEBRAICMATRIX

#include <dace/dace.h>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include "camRoutines.h"
#include "auxiliaryRoutines.h"
#include <chrono>

using namespace std;
using namespace DACE;
using namespace std::chrono;
using namespace cam;
using namespace aux;

int main(void)
{
    int nvar, j, jj, kk, i, ii, k, vv, flag1, flag2, flag3, order, pocType, N, lowThrust_flag, n_conj, n_man, m, dyn, gravOrd, missDistanceFlag;
	double mass, A_drag, Cd, A_srp, Cr, tca, Lsc, musc, ctrlMax, mean_motion_p;

    long time1 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    
    ifstream nodes;
	nodes.open("./write_read/initial_state.dat");
        nodes >> N;         // Number of nodes
        nodes >> n_conj;    // Number of conjunctions
        nodes >> n_man;     // Number of control nodes
        nodes >> m;         // Number of DA variables per node
	nodes.close();
    // Initialize variable
    AlgebraicMatrix<double> Cp(6,6), Cs(6,6), covp(36,n_conj), covs(36,n_conj), r2e(3,3), ctrlDum(m,n_man), rsDum(3,n_conj), vsDum(3,n_conj), directions(3,n_man);
    AlgebraicVector<double> xdum(6), metricMap(3), t(N), HBR(n_conj), magnitude(n_man), rRef(3), vRef(3), mean_motion_s(n_conj);
    AlgebraicVector<int>    canFire(N), isConj(N), isRet(N), constraintFlags(6);
    // Read input from .dat
	readInit( nvar,  order,  pocType,  N,  lowThrust_flag,  n_conj,  n_man,  m,  dyn,  gravOrd,  missDistanceFlag, tca,  Lsc,  musc,  ctrlMax,  mean_motion_p, covp,   covs,   ctrlDum,   rsDum,  vsDum,
               directions,  xdum,  t,  HBR,  magnitude, rRef,  vRef,  mean_motion_s,  canFire,  isConj, isRet,  constraintFlags);

    // initialize execution time counter
    long time2 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    long timeSubtr = time2 - time1;
    // Initialize DA    
    nvar = m*n_man + 1; 
    DA::init(order, nvar);
    DA::setEps(1e-30);
    
    // Initialize DA variables
    AlgebraicVector<DA> x0(6), x00(6), xsf(6), xs0(6), xBall(6), xf(6), r(3), r_rel(2), v(3), rB(3), ctrlRtn(3), ctrl(3), rf(3), rs(3), rsB(3), vs(3), xRet(6), poc(n_conj), md(n_conj), dx1(3), dx2(3), dx3(3); 
    AlgebraicMatrix<DA> xTca(6,n_conj), xsTca(6,n_conj), P_eci(3,3), P_B3(3,3), P_B(2,2), Pp(3,3), Ps(3,3), toB(3,3), STM_p(6,6), STM_s(6,6), CPropP(6,6), CPropS(6,6), covsda(9,n_conj), r2ep(3,3), r2es(3,3);
    DA poc_tot, alpha, beta, tcaNew, tan, radial;

    dx1[0] = DA(1); dx1[1] = 0.0; dx1[2] = 0.0;
    dx2[0] = 0; dx2[1] = DA(2); dx2[2] = 0.0;
    dx3[0] = 0; dx3[1] = 0.0; dx3[2] = DA(3);


    // Define ballistic primary's position at first TCA 
    for (j = 0; j < 6 ; j++) {xBall[j] = xdum[j] + 0*DA(1);}
    // backpropagation from first TCA
    if (dyn == 0) {
        if (gravOrd == 0) {
            x0     = RK78(6, xBall, {0.0*DA(1),0.0*DA(1),0.0*DA(1)}, 0.0, - t[0], keplerPropAcc, 1.0, Lsc); // Earth Orbit
        }
        else if (gravOrd == 2) {
            x0     = RK78(6, xBall, {0.0*DA(1),0.0*DA(1),0.0*DA(1)}, 0.0, - t[0], J2dynamics, 1.0, Lsc); // Earth Orbit
        }        
        else if (gravOrd == 4) {
            x0     = RK78(6, xBall, {0.0*DA(1),0.0*DA(1),0.0*DA(1)}, 0.0, - t[0], J2_J4dynamics, 1.0, Lsc); // Earth Orbit
        }
        else {throw std::runtime_error("The gravOrd flag should be 0, 2, or 4");}
    }
    else if (dyn == 1) {
        x0     = RK78(6, xBall, {0.0*DA(1),0.0*DA(1),0.0*DA(1)}, 0.0, - t[0], CR3BPsyn, musc, Lsc); // Cislunar Orbit
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
            r2e     = rtn2eci(cons(x0));
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
                ctrl = r2e*ctrlRtn*ctrlMax;
            }
            else {
                ctrl = ctrlRtn*ctrlMax;
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
            if (gravOrd == 0) {
                x0 = RK78(6, x0, ctrl, tca - t[i], tca - t[i+1], keplerPropAcc, 1.0, Lsc); // Earth Orbit
            }
            else if (gravOrd == 2) {
                x0 = RK78(6, x0, ctrl, tca - t[i], tca - t[i+1], J2dynamics, 1.0, Lsc); // Earth Orbit
            }        
            else if (gravOrd == 4) {
                x0 = RK78(6, x0, ctrl, tca - t[i], tca - t[i+1], J2_J4dynamics, 1.0, Lsc); // Earth Orbit
            }
            else {throw std::runtime_error("The gravOrd flag should be 0, 2, or 4");}
        }
        else {
            x0 = RK78(6, x0, ctrl, tca - t[i], tca - t[i+1], CR3BPsyn, musc, Lsc); // Cislunar Orbit
        }

        // If the next node is a conjunction node, save the state in a DA variable
        vv = 0;

        if (isConj[i+1] == 1) {
            DA dt = 0.0 + DA(nvar);
            ctrl = {DA(1)*0, DA(1)*0, DA(1)*0};
            xs0[0] = rsDum.at(0,k) + DA(1)*0; xs0[1] = rsDum.at(1,k) + DA(1)*0; xs0[2] = rsDum.at(2,k) + DA(1)*0;
            xs0[3] = vsDum.at(0,k) + DA(1)*0; xs0[4] = vsDum.at(1,k) + DA(1)*0; xs0[5] = vsDum.at(2,k) + DA(1)*0; 
            x00  = KeplerProp(x0, dt, 1.0);
            xsf  = KeplerProp(xs0, dt, 1.0);
            tcaNew = findTCA(x00 - xsf, nvar);
            STM_p = CWSTM(mean_motion_p,tcaNew); // CW transformation for the covariance (can be done with YA)
            STM_s = CWSTM(mean_motion_s[k],tcaNew); // CW transformation for the covariance (can be done with YA)
            for (ii = 0; ii < 6 ; ii ++) {
                for (j = 0; j < 6 ; j ++) {
                    Cp.at(ii,j)  = covp.at(vv,k);
                    Cs.at(ii,j)  = covs.at(vv,k);
                    vv ++;
                }
            }
            CPropP = STM_p*Cp*STM_p.transpose();
            CPropS = STM_s*Cs*STM_s.transpose();
            for (ii = 0; ii < 3 ; ii ++) {
            for (j = 0; j < 3 ; j ++) {
                Pp.at(ii,j)  = CPropP.at(ii,j);
                Ps.at(ii,j)  = CPropS.at(ii,j);
            }
            }
            r2ep = rtn2eci(x00);
            r2es = rtn2eci(xsf);
            P_eci = r2ep*Pp*r2ep.transpose() + r2es*Ps*r2es.transpose();

            AlgebraicVector<DA> dx(nvar);
            for (ii = 0; ii < nvar - 1; ii++) {
                dx[ii] = DA(ii+1);}
            dx[nvar - 1] = tcaNew;
            x00  = x00.eval(dx);
            xsf = xsf.eval(dx);
            P_eci = evalDAMatrix(P_eci,dx,3);
            vv = 0;
            for (ii = 0; ii < 3 ; ii ++) {
                for (j = 0; j < 3; j ++) {
                    covsda.at(vv,k) = P_eci.at(ii,j);
                    vv ++;
                }
            }
            for (j = 0; j < 6 ; j ++) {
                xTca.at(j,k) = x00[j];
                xsTca.at(j,k) = xsf[j];
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
        for (i = 0; i < 3 ; i ++) {
            r[i] = xTca.at(i,k);
            v[i] = xTca.at(i+3,k);
            rs[i] = xsTca.at(i,k);
            vs[i] = xsTca.at(i+3,k);
        }
        vv = 0;
        for (ii = 0; ii < 3 ; ii ++) {
        for (j = 0; j < 3 ; j ++) {
            P_eci.at(ii,j) = covsda.at(vv,k);
            vv ++;
        }
        }
        // B-plane transformations
        toB  = Bplane(v,vs); // DCM from ECI to B-plane
        rB   = toB*r;                     // Primary position in the B-plane (3D)
        rsB  = toB*rs;                    // Secondary position in the B-plane (3D)
        P_B3 = toB*P_eci*toB.transpose();     // Combined covariance in the B-plane (3D)

        // Relative position in the B-plane (2D)
        r_rel[0]    = rB[0] - rsB[0];   r_rel[1]    = rB[2] - rsB[2];
        if (missDistanceFlag == 1 ) {
            md[k] = dot(r_rel,r_rel);
            }
        else {
            // Combined covariance in the B-plane (2D)
            P_B.at(0,0) = P_B3.at(0,0);     P_B.at(0,1) = P_B3.at(0,2);
            P_B.at(1,0) = P_B3.at(2,0);     P_B.at(1,1) = P_B3.at(2,2);
            vv = 0;
            // Compute PoC for the single conjunction, according to the required model
            if (pocType == 0) {
                poc[k] = ConstPoC(r_rel,P_B,HBR[k]);}
            else if (pocType == 1) {
                poc[k] = ChanPoC(r_rel,P_B,HBR[k],3);}
            else if (pocType == 2) {
                poc[k] = MaxPoC(r_rel,P_B,HBR[k]);}
            else {
                throw std::runtime_error("the metric flag must be in the interval [1,3] and the PoC type must be in the interval [0,1]");}
            // Probability of no collision
            noCollisions = noCollisions*(1.0 - poc[k]);
        }
    }
    // Final PoC comprehensive of all the conjunctions
    if (missDistanceFlag == 0 ) {
        poc_tot = log10(1.0 - noCollisions);
    }
}

if (constraintFlags[1] == 1 || constraintFlags[2] == 1) {
    // Return constraint (scalar as tangential distance from other spacecraft GRACE)
    AlgebraicVector<DA> rRet(3), vRet(3), distRel(3);
    for (j = 0; j < 3 ; j ++) {
        rRet[j] = xRet[j];
        vRet[j] = xRet[j+3];
    }      
    r2e     = rtn2eci(cons(xRet));
    distRel = r2e.transpose()*(rRet-rRef);
    radial  = distRel[0]; // radial displacement with respect to reference
    tan     = distRel[1]; // tangential displacement with respect to reference
}

// stability of the monodromy matrix (Cauchy-Green Tensor) after an orbital period
// if (constraintFlags[4] == 1 || dyn == 1) {
//     ctrl = {0.0*DA(1), 0.0*DA(1), 0.0*DA(1)};
//     x0 = RK78Cis(6, x0, ctrl, -t[end], -t[end] +     1, CR3BPsyn, 0.012150668, 0.0); // forward propagation to the next node
//     AlgebraicMatrix<double> STM = stmDace(x0, 3, 6);
//     AlgebraicMatrix<double> CG = STM.transpose()*STM;
// }

    time1   = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    //open the output files
    ofstream constPart, constraints, tcaOut, convRad;
    constPart.open("./write_read/constPart.dat");
    constPart << setprecision(18);
    
    // write TCA states in ECI coordinates into output
    for (k = 0; k < n_conj ; k++) {
        for (j = 0; j < 6 ; j++) {
        constPart  << cons(xTca.at(j,k)) << endl;
        }
    }
    constPart.close();

    // write the DA expansion of PoC in output
    constraints.open("./write_read/constraints.dat");
    constraints << setprecision(18);
    if (constraintFlags[0] == 1) {
        for (k = 0; k < n_conj; k ++) {
        if  (missDistanceFlag == 1) {
            constraints << md[k] << endl;
        }
        else {    
            constraints << log10(poc[k]) << endl;
        }
        }
        constraints << poc_tot << endl;
    }

    if (constraintFlags[1] == 1) {
        constraints << tan    << endl;
    }

    // write the DA expansion of return in output
    if (constraintFlags[2] == 1) {
        constraints << radial       << endl;
    }

    // write the DA expansion of return in output
    if (constraintFlags[3] == 1) {
        for (j = 0; j < 6; j ++) {
        constraints << xRet[j]        << endl;
        }
    }

    if (constraintFlags[4] == 1) {
        constraints << pow(xRet - xdum, 2)  << endl;
    }

    if (constraintFlags[5] == 1) {
        for (j = 0; j < n_man; j ++) {
            constraints << DA(1 + j*m)*DA(1 + j*m) + DA(2 + j*m)*DA(2 + j*m) + DA(3 + j*m)*DA(3 + j*m) << endl;
        }
    }
    constraints.close();

    if  (missDistanceFlag == 0) {
    convRad.open("./write_read/convRad.dat");
    convRad << setprecision(18);
    convRad << convRadius(poc_tot.eval(dx1),1e-8) << endl;
    convRad << convRadius(poc_tot.eval(dx2),1e-8) << endl;
    convRad << convRadius(poc_tot.eval(dx3),1e-8) << endl;    // write the DA expansion of return in output
    convRad.close();
    }
    // Do not consider writing time when calculating execution time
    time2 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    timeSubtr = timeSubtr + time2 - time1;

    ofstream timeOut;
    timeOut.open("./write_read/timeOut.dat");
    timeOut << setprecision(16);
    timeOut << timeSubtr << endl;
    timeOut.close();
}