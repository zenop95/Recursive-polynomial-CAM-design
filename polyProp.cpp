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
    int nvar, j, jj, kk, i, ii, k, order, pocType, N, lowThrust_flag, n_conj, n_man, m, dyn, gravOrd, missDistanceFlag, TPoCFlag;
	double tca, Lsc, musc, ctrlMax, mean_motion_p;

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
    AlgebraicVector<double> xdum(6), metricMap(3), t(N), HBR(n_conj), mean_motion_s(n_conj);
    AlgebraicVector<int>    canFire(N), isConj(N), isRet(N), constraintFlags(4);
    // Read input from .dat
	readInit( nvar,  order,  pocType,  N,  lowThrust_flag,  n_conj,  n_man,  m,  dyn,  gravOrd,  missDistanceFlag, TPoCFlag, tca,  Lsc,  musc,  ctrlMax,  mean_motion_p, covp,   covs,   ctrlDum,   rsDum,  vsDum,
               directions,  xdum,  t,  HBR, mean_motion_s,  canFire,  isConj, isRet,  constraintFlags);

    // initialize execution time counter
    long time2 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    long timeSubtr = time2 - time1;
    // Initialize DA    
    nvar = m*n_man + 1; 
    DA::init(order, nvar);
    DA::setEps(1e-30);
    
    // Initialize DA variables
    AlgebraicVector<DA> x0(6), x00(6), xsf(6), xs0(6), xBall(6), xf(6), r(3), r_rel(2), v(3), rB(3), ctrlRtn(3), ctrl(3), rf(3), rs(3), rsB(3), vs(3), xRet(6), poc(n_conj), md(n_conj), dx(nvar-1), rRet(3), vRet(3), distRel(3), coe(6), meanCoe(6); 
    AlgebraicMatrix<DA> xTca(6,n_conj), xsTca(6,n_conj), P_eci(3,3), P_B3(3,3), P_B(2,2), Pp(3,3), Ps(3,3), toB(3,3), STM_p(6,6), STM_s(6,6), CPropP(6,6), CPropS(6,6), covsda(9,n_conj), r2ep(3,3), r2es(3,3);
                    DA  poc_tot, alpha, beta, tcaNew, tan, radial, retErrP, retErrV, meanSma, meanEcc;

    
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
            if (dyn + lowThrust_flag == 0) {
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
            unpackMatrix(Cp,covp,k,6);
            unpackMatrix(Cs,covs,k,6);
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

            AlgebraicVector<DA> dxx(nvar);
            for (ii = 0; ii < nvar - 1; ii++) {
                dxx[ii] = DA(ii+1);}
            dxx[nvar - 1] = tcaNew;
            x00  = x00.eval(dxx);
            xsf = xsf.eval(dxx);
            P_eci = evalDAMatrix(P_eci,dxx,3);
            packMatrix(P_eci,covsda,k,3);
            for (j = 0; j < 6 ; j ++) {
                xTca.at(j,k) = x00[j];
                xsTca.at(j,k) = xsf[j];
            }
            k = k + 1;
        }
        else if (isRet[i+1] == 1) {
            xRet = x0;
            for (j = 0; j < 3 ; j ++) {
                rRet[j] = xRet[j];
                vRet[j] = xRet[j+3];
            }
            coe     = cart2kep(xRet,1.0);
            meanCoe = osculating2mean(coe,1.0,Lsc);
            // coe = coe2mee(coe);
            meanSma = meanCoe[0];
            meanEcc = meanSma*meanCoe[1];
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
        unpackMatrix(P_eci,covsda,k,3);
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
        poc_tot = 1.0 - noCollisions;
    }
}

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
    for (j = 0; j < 6 ; j++) {
        constPart  << cons(xRet[j]) << endl;
    }
    constPart.close();

    // write the DA expansion of PoC in output
    constraints.open("./write_read/constraints.dat");
    constraints << setprecision(18);
    if (constraintFlags[0] == 1) {
        for (k = 0; k < n_conj; k ++) {
            if  (missDistanceFlag == 1) {
                constraints << -md[k] << endl;
            }
            else if (missDistanceFlag == 0 && TPoCFlag == 0) {    
                constraints << log10(poc[k]) << endl;
            }
        }
    if (TPoCFlag == 1) {
        constraints << log10(poc_tot) << endl;
    }
    }
    
    DA a, e, om, Om, Inc;
    if (constraintFlags[1] == 1) {
        a  = coe[0] - cons(coe[0]);
        e  = coe[1] - cons(coe[1]);
        Inc = coe[2] - cons(coe[2]);
        Om = coe[3] - cons(coe[3]);
        om = coe[4] - cons(coe[4]);
        constraints << a << endl;
        constraints << e << endl;
        constraints << Inc << endl;
        constraints << Om << endl;
        constraints << om << endl;
    }
    double eps = 0;//1e-30;
    if (constraintFlags[2] == 1) {
        a = (dot(rRet - cons(rRet), rRet - cons(rRet)) + eps)*1e5;
        e = (dot(vRet - cons(vRet), vRet - cons(vRet)) + eps)*1e5;
        constraints << a << endl;
        constraints << e << endl;
    }
    
    if (constraintFlags[3] == 1) {
        a = ((meanSma - cons(meanSma))*(meanSma - cons(meanSma)))*1e10;
        e = ((meanEcc-cons(meanEcc))*(meanEcc-cons(meanEcc)))*1e10;
        constraints << a   << endl;
        constraints << e  << endl;
    }
    constraints.close();

    // convergence radius

    // convRad.open("./write_read/convRad.dat");
    // convRad << setprecision(18);
    // j = 0;
    // for (j == 0; j < nvar-1; j ++) {
    //     if (j > 0) {dx[j-1] = 0.0;}
    //     dx[j] = DA(j+1);
    //     if  (missDistanceFlag == 0) {
    //         convRad << convRadius(poc_tot.eval(dx),1e-8) << endl;
    //     }
    //     else {
    //         convRad << convRadius(md[0].eval(dx),1e-8) << endl;
    //     }
    //     convRad.close();
    // }
    // Do not consider writing time when calculating execution time
    time2 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    timeSubtr = timeSubtr + time2 - time1;

    ofstream timeOut;
    timeOut.open("./write_read/timeOut.dat");
    timeOut << setprecision(16);
    timeOut << timeSubtr << endl;
    timeOut.close();
}