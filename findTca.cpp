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
    int i;
    ifstream nodes;
    AlgebraicVector<double> xdump(6), xdums(6);

    ifstream Input;
    long time1 = time_point_cast<milliseconds>(system_clock::now()).time_since_epoch().count();
	Input.open("./write_read/initial_state.dat");
        for (i = 0; i < 6; i ++) {
            Input >> xdump[i];
        }            
        for (i = 0; i < 6; i ++) {
            Input >> xdums[i];
        }
	Input.close();

    // DA initialization    
    DA::init(1, 7);
    DA::setEps(1e-30);
    AlgebraicVector<DA> xp0(6), xpf(6), xs0(6), xsf(6), xrel(6);
    DA dt, tca;
    // write output
    ofstream tcaOut;    
    tcaOut.open("./write_read/tcaOut.dat");
    tcaOut << setprecision(30);

    // final time espansion
    dt = 0.0 + DA(7);
    for (i = 0; i < 6; i ++) {
        xs0[i] = xdums[i] + DA(i+1);
        xp0[i] = xdump[i] + DA(i+1)*0; 
    }
    xpf = KeplerProp(xp0, dt, 1.0);
    xsf = KeplerProp(xs0, dt, 1.0);
    xrel = xpf - xsf;
    tca = findTCA(xrel, 7);
    tcaOut  << cons(tca) << endl;
    tcaOut.close();

}