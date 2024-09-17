#include <dace/dace.h>
#include <limits>
#include <fstream>


using namespace DACE;
using namespace std;

namespace aux {

// Unpack a 1D array into a matrix, automatically inferring dimensions
 template <typename T> 
 void unpackMatrix(AlgebraicMatrix<T>& mat, const AlgebraicMatrix<T>& MAT, const int k, const int n) {    
    int index = 0;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            mat.at(j, i) = MAT.at(index,k);
            index ++;
        }
    }
}

// Pack a matrix back into a 1D array, automatically inferring dimensions
template <typename T> 
void packMatrix(const AlgebraicMatrix<T>& mat, AlgebraicMatrix<T>& MAT, const int k, const int n) {
    int index = 0;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            MAT.at(index,k) = mat.at(j, i);
            index ++;
        }
    }
}

void readInit(int& nvar, int& order, int& pocType, int& N, int& lowThrust_flag, int& n_conj, int& n_man, int& m, int& dyn, int& gravOrd, int& missDistanceFlag, 
              int& TPoCFlag, double& tca, double& Lsc, double& musc, double& ctrlMax, double& mean_motion_p,
              AlgebraicMatrix<double>&  covp, AlgebraicMatrix<double>&  covs, AlgebraicMatrix<double>&  ctrlDum, AlgebraicMatrix<double>&  rsDum, AlgebraicMatrix<double>& vsDum, 
              AlgebraicMatrix<double>& directions, AlgebraicVector<double>& xdum, AlgebraicVector<double>& t, AlgebraicVector<double>& HBR, AlgebraicVector<double>& magnitude, 
              AlgebraicVector<double>& mean_motion_s, AlgebraicVector<int>& canFire, AlgebraicVector<int>& isConj, 
              AlgebraicVector<int>& isRet, AlgebraicVector<int>& constraintFlags) {

int k, j;
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
        Input >> ctrlMax;            // ctrlMax
        Input >> missDistanceFlag;        // maneuver on miss distance
        Input >> TPoCFlag;             // maneuver TPoC
        Input >> mean_motion_p;        // mean motion primary
        for (k = 0; k < n_conj; k ++) {
            Input >> mean_motion_s[k];        // mean motion secondary
        }
        for (j = 0; j < 6; j ++) {
            Input >> xdum[j];           // Write dummy variable for the primary's ECI state at the first conjunction from input
        }
        for (k = 0; k < n_conj; k ++) {
            for (j = 0; j < 3; j ++) {
                Input >> rsDum.at(j,k); // Write secondary's ECI position at each conjunction from input
            }
            for (j = 0; j < 3; j ++) {
                Input >> vsDum.at(j,k); // Write secondary's ECI velocity at each conjunction from input
            }   
        }
        for (k = 0; k < n_conj; k ++) {
            Input >> HBR[k];            // Combined Hard Body radialius for each conjunction
        }
        for (k = 0; k < n_conj; k ++) {
            for (j = 0; j < 36; j ++) {
                Input >> covp.at(j,k);   // Write dummy variable for the combined covariance in ECI at each conjunction from input
            }
            for (j = 0; j < 36; j ++) {
                Input >> covs.at(j,k);   // Write dummy variable for the combined covariance in ECI at each conjunction from input
            }
        }
        for (k = 0; k < n_man; k ++) {
                Input >> magnitude[k];    // Write maneuver magnitude if the magnitude is fixed
        }   
        for (k = 0; k < n_man; k ++) {
            for (j = 0; j < 3; j ++) {
                Input >> directions.at(j,k); // Write maneuver direction in RTN if the direction is fixed
            }
        }
        for (k = 0; k < 4; k ++) {
            Input >> constraintFlags[k];
        }
        for (k = 0; k < n_man; k ++) {
            for (j = 0; j < m; j ++) {
                Input >> ctrlDum.at(j,k);   // Write dummy reference control from input ({0 0 0} if ballistic trajectory)
            }
        }
        for (k = 0; k < N; k ++) {
            Input >> t[k];                  // Write node times before conjunction (if negative it means that there are more than one encounters and nodes are needed after the first encounter)
            Input >> canFire[k];            // Define if node k is a firing node
            Input >> isConj[k];             // Define if node k is a conjunction node
            Input >> isRet[k];              // Define if node k is a return node
        }
    Input.close();
}
}