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
    //open the output files
    ofstream constPart, metricPoly;
    constPart.open("constPart.dat");
    constPart << setprecision(18);
    metricPoly.open("metricPoly.dat");
    metricPoly << setprecision(18);

    int j, jj, kk, i, k, vv, flag1, flag2, flag3, gravOrd, order, metricFlag, pocType, N, lowThrust_flag, n_conj, n_man, m, dyn;
	double mass, A_drag, Cd, A_srp, Cr, tca, Lsc, musc;

    long time1 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    // AIDA initialization input file
    ifstream aidaParams;
	aidaParams.open("AIDA_init.dat");
        aidaParams >> flag1;	
        aidaParams >> flag2;	
        aidaParams >> flag3;	
        aidaParams >> gravOrd;	
        aidaParams >> mass;	
        aidaParams >> A_drag;	
        aidaParams >> Cd;
        aidaParams >> A_srp;	
        aidaParams >> Cr;	
	aidaParams.close();
  
    ifstream nodes, Input;
	nodes.open("initial_state.dat");
        nodes >> N; // Number of impulses
        nodes >> n_conj; // Number of impulses
        nodes >> n_man; // Number of impulses
        nodes >> m; // Number of DA variables per node
	nodes.close();
    AlgebraicMatrix<double> P(3,3), cov(9,n_conj), covB(4,n_conj), P_B3(3,3), P_B(2,2), r2e(3,3), toB(3,3), ctrlDum(m,n_man), scale(m,n_man), rsDum(3,n_conj), vsDum(3,n_conj), directions(3,n_man);
    AlgebraicVector<double> xdum(6), xstart(6), metricMap(3), t(N), vs(3), rs(3), rsB(3), detP(n_conj), HBR(n_conj), magnitude(n_man);
    AlgebraicVector<int>    canFire(N), isConj(N);
	Input.open("initial_state.dat");
        Input >> N; // Number of nodes
        Input >> n_conj; // Number of conjunctions
        Input >> n_man; // Number of impulses
        Input >> m; // Number of DA variables per node
        Input >> dyn; // Dynamics model
        Input >> lowThrust_flag;
        Input >> order;
        Input >> metricFlag;
        Input >> pocType;
        Input >> tca;    //ephemeris time
        Input >> Lsc;  //length scale
        Input >> musc;  
        for (j = 0; j < 6; j ++) {
            Input >> xdum[j];
        }
        for (k = 0; k < n_conj; k ++) {
            Input >> HBR[k];  //Hard Body Radius
        }
        for (k = 0; k < n_conj; k ++) {
            for (j = 0; j < 3; j ++) {
                Input >> rsDum.at(j,k);
            }
            for (j = 0; j < 3; j ++) {
                Input >> vsDum.at(j,k);
            }   
        }
        for (k = 0; k < n_conj; k ++) {
            for (j = 0; j < 9; j ++) {
                Input >> cov.at(j,k);
            }
        }
        for (i = 0; i < n_man; i ++) {
            for (j = 0; j < m; j ++) {
                Input >> scale.at(j,i);
            }
        }
        for (i = 0; i < n_man; i ++) {
                Input >> magnitude[i];
        }   
        for (i = 0; i < n_man; i ++) {
            for (j = 0; j < 3; j ++) {
                Input >> directions.at(j,i);
            }
        }
        for (i = 0; i < n_man; i ++) {
            for (j = 0; j < m; j ++) {
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
    // DA initialization    
    DA::init(order, m*n_man);
    DA::setEps(1e-30);
    // AIDA initialization    
    string kern = "./data/kernel.txt";
    ConstSpiceChar*FileKernel = kern.c_str();
    furnsh_c (FileKernel);
    string gravmodel  = "./data/gravmodels/egm2008";
    int AIDA_flags[3] = {flag1,flag2,flag3};
    double Bfactor    = Cd*A_drag/mass;
    double SRPC       = Cr*A_srp/mass;
    AIDADvDynamics<DA> aidaCartDynImp(gravmodel, gravOrd, AIDA_flags, Bfactor, SRPC);
    AIDAScaledDynamics<DA> aidaCartDynLT(gravmodel, gravOrd, AIDA_flags, Bfactor, SRPC);

    AlgebraicVector<DA> x0(6), xBall(6), xf(6), r(3), r_rel(2), v(3), rB(3), ctrlRtn(3), ctrl(3), rf(3); 
    AlgebraicMatrix<DA> xTca(6,n_conj);
    DA metric, poc_tot, alpha, beta;
    for (j = 0; j < 6 ; j++) {xBall[j] = xdum[j] + 0*DA(1);}
    if (dyn == 0) {
        x0     = RK78Dv(6, xBall, tca, tca - t[0], Lsc, musc, gravOrd, aidaCartDynImp); // backpropagation from first TCA
    }
    else if (dyn == 1) {
        x0     = RK78Cis(6, xBall, {0.0*DA(1),0.0*DA(1),0.0*DA(1)}, 0.0, - t[0], CR3BPsyn, 0.012150668, 0.0); // backpropagation from TCA
    }
    else {throw std::runtime_error("The dynamics flag should be 0 for Earth orbit and 1 for Cislunar");}
    xstart = cons(x0);
    // Propagations at each maneuvering time or TCA
    jj     = 0;
    k      = 0;
    kk     = 0;
    for (i = 0; i < N-1; i ++) {
       if (canFire[i] == 1) {
            r2e     = astro::rtn2eci(cons(x0));
            if (m == 3) {
                for (j = 0; j < 3 ; j++) { 
                    jj ++;
                    ctrlRtn[j] = ctrlDum.at(j,kk) + DA(jj)*scale.at(j,kk);
                }
            }
            else if (m == 2) {
                alpha = DA(2*kk+1)*scale.at(0,kk);
                beta  = DA(2*kk+2)*scale.at(1,kk);
                ctrlRtn[0] = magnitude[kk]*cos(alpha)*sin(beta);
                ctrlRtn[1] = magnitude[kk]*cos(alpha)*cos(beta);
                ctrlRtn[2] = magnitude[kk]*sin(alpha);
            }
            else if (m == 1) {
                for (j = 0; j < 3 ; j++) {
                    ctrlRtn[j] = (ctrlDum.at(0,kk) + DA(kk+1)*scale.at(0,kk))*directions.at(j,kk);
                }
            }
            else {
                throw std::runtime_error("The number of DA variables per node must be in the interval [1,3]");
            }
            kk ++;
            if (dyn == 0) {
                ctrl = r2e*ctrlRtn;
            }
            else {
                ctrl = ctrlRtn;
            }
        }
        else {ctrl = {DA(1)*0, DA(1)*0, DA(1)*0};}
        if (lowThrust_flag == 0 || canFire[i] == 0) {
            for (j = 3; j < 6 ; j++) {x0[j] = x0[j] + ctrl[j-3];}
            if (dyn == 0) {
                x0 = RK78Dv(6, x0, tca - t[i], tca - t[i+1], Lsc, musc, gravOrd, aidaCartDynImp);   // forward propagation to TCA
            }
            else {
                x0 = RK78Cis(6, x0, {0.0*DA(1),0.0*DA(1),0.0*DA(1)}, -t[i], -t[i+1], CR3BPsyn, 0.012150668, 0.0); // backpropagation from TCA
            }
        }
        else {
            if (dyn == 0) {
                x0 = RK78Sc(6, x0, ctrl, tca - t[i], tca - t[i+1], 1.0, Lsc, 0, gravOrd, aidaCartDynLT);   // forward propagation to TCA
            }
            else {
                x0 = RK78Cis(6, x0, ctrl, -t[i], -t[i+1], CR3BPsyn, 0.012150668, 0.0); // backpropagation from TCA
            }
        }
        if (isConj[i+1] == 1) {
            for (j = 0; j < 6 ; j ++) {
                xTca.at(j,k) = x0[j];
            }        
            k = k + 1;
        }
    }

    DA noCollisions = 1.0;
    for (k = 0; k < n_conj; k ++) {

        r[0] = xTca.at(0,k);  r[1] = xTca.at(1,k);  r[2] = xTca.at(2,k);   
        v[0] = xTca.at(3,k);  v[1] = xTca.at(4,k);  v[2] = xTca.at(5,k);

        vv = 0;
        for (i = 0; i < 3 ; i ++) {
            for (j = 0; j < 3 ; j ++) {
                P.at(i,j)  = cov.at(vv,k);
                vv = vv + 1;
            }
            vs[i] = vsDum.at(i,k);
            rs[i] = rsDum.at(i,k);
        }
        // B plane transformations
        toB = astro::Bplane(cons(v),vs);
        rB  = toB*r;
        rsB = toB*rs;
        P_B3 = toB*P*toB.transpose();

        r_rel[0]    = rB[0] - rsB[0];   r_rel[1]    = rB[2] - rsB[2]; 
        P_B.at(0,0) = P_B3.at(0,0);     P_B.at(0,1) = P_B3.at(0,2);
        P_B.at(1,0) = P_B3.at(2,0);     P_B.at(1,1) = P_B3.at(2,2);

        vv = 0;
        for (i = 0; i < 2 ; i ++) {
            for (j = 0; j < 2 ; j ++) {
                covB.at(vv,k)  = P_B.at(i,j);
                vv = vv + 1;
            }
        }

        detP[k] = astro::det2(P_B);

        if (metricFlag == 1 && pocType == 0) {
            metric = astro::ConstPoC(r_rel,P_B,HBR[k]);}
        else if (metricFlag == 1 && pocType == 1) {
            metric = astro::ChanPoC(r_rel,P_B,HBR[k],3);}
        else {
            throw std::runtime_error("the metric flag must be in the interval [1,3] and the PoC type must be in the interval [0,1]");}
        noCollisions = noCollisions*(1.0 - metric);
}
    poc_tot = log10(1.0 - noCollisions);
    time1 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    for (k = 0; k < n_conj ; k++) {
        for (j = 0; j < 6 ; j++) {
        constPart  << cons(xTca.at(j,k)) << endl;
        }
    }
    constPart  << cons(poc_tot)  << endl;
    for (k = 0; k < n_conj ; k++) {
        constPart  << detP[k]  << endl;
    }
    for (k = 0; k < n_conj ; k++) { 
        for (j = 0; j < 4 ; j++) {
            constPart  << covB.at(j,k) << endl;
        }
    }
    metricPoly << poc_tot << endl;

    time2 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
    timeSubtr = timeSubtr + time2 - time1;

    ofstream timeOut;
    timeOut.open("timeOut.dat");
    timeOut << setprecision(16);
    timeOut << timeSubtr << endl;
    timeOut.close();
    constPart.close();
    metricPoly.close();
}