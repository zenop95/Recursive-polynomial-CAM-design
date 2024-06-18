#include <dace/dace.h>
#include <limits>


using namespace DACE;
using namespace std;

namespace cam {

int factorial(int n)
{
    // single line to find factorial
    return (n==1 || n==0) ? 1: n * factorial(n - 1);
}

double normtmp(const int N, const std::vector<double>& X)
{
     int I;  
     double res = 0.0;
     for (I=0; I<N; I++)  {  
          res = std::max(res, std::fabs(X[I]));  
      }
     return res;
}

template <typename T> T checkatan2(T angle, T a, T b) {
    if ( DACE::cons(a) > 0 && DACE::cons(b) < 0 ) return angle + M_PI;
    else if ( DACE::cons(a) < 0 && DACE::cons(b) < 0 ) return angle - M_PI;
    else return angle;
}

template <typename T> T atan2_mod(T a, T b) {
    T angle = atan(a/b);
    checkatan2(angle, a, b);
    return angle;
}

template<typename T> AlgebraicVector<T> keplerPropAcc(AlgebraicVector<T> x, AlgebraicVector<T> uEci, double t, double mu, double Lsc)
{     
    AlgebraicVector<T> res(6), pos(3);
    AlgebraicMatrix<double> dcm(3,3);
    
    pos[0] = x[0]; pos[1] = x[1]; pos[2] = x[2];    
    T rrr = pow(pos.vnorm(),3);

    res[0] = x[3];
    res[1] = x[4];
    res[2] = x[5];
	res[3] = -mu*pos[0]/rrr + uEci[0];
	res[4] = -mu*pos[1]/rrr + uEci[1];
	res[5] = -mu*pos[2]/rrr + uEci[2];

    return res;
    
}

//---------------------------------------------------------------------

template<typename T> AlgebraicVector<T> CR3BPsyn(AlgebraicVector<T> x, AlgebraicVector<T> u, double t, double mu, double Lsc)
{
    
    AlgebraicVector<T> pos(3), res(6);

    pos[0] = x[0]; pos[1] = x[1]; pos[2] = x[2];
    
    res[0] = x[3] ;
    res[1] = x[4] ;
    res[2] = x[5] ;

    T  r1  = sqrt( (pos[0]+mu)*(pos[0]+mu) + pos[1]*pos[1] + pos[2]*pos[2] );
    T  r2  = sqrt( (pos[0]-(1-mu))*(pos[0]-(1-mu)) + pos[1]*pos[1] + pos[2]*pos[2] );
    
    res[3] = ( 2*x[4] + pos[0] - (1-mu)*(pos[0]+mu)/(r1*r1*r1) - mu*(pos[0]-(1-mu))/(r2*r2*r2) ) + u[0];
    res[4] = (-2*x[3] + pos[1] - (1-mu)*   pos[1]  /(r1*r1*r1) - mu*    pos[1]     /(r2*r2*r2) ) + u[1];
    res[5] = (-(1-mu)*pos[2]/(r1*r1*r1) - mu*pos[2]/(r2*r2*r2) ) + u[2];
    
    return res;
}

//---------------------------------------------------------------------

template<typename T> AlgebraicVector<T> J2dynamics(AlgebraicVector<T> xx, AlgebraicVector<T> uEci, double t, double mu, double Lsc)
{
	const double J2 = 1.08262668e-3; //{-}
    const double rE = 6378.137/Lsc;

    AlgebraicVector<T> pos(3), res(6), vel(3);
    
    pos[0] = xx[0]; pos[1] = xx[1]; pos[2] = xx[2];
    vel[0] = xx[3]; vel[1] = xx[4]; vel[2] = xx[5];
    
    T r = pos.vnorm();
    T v = vel.vnorm();
    
    res[0] = xx[3];
    res[1] = xx[4];
    res[2] = xx[5];

    T x = pos[0];
	T y = pos[1];
	T z = pos[2];
	
	T mur3 = mu/r/r/r;
	T z2r2 = z/r*z/r;

	T J2rEr = 1.5*J2*rE/r*rE/r;

	res[3] = -pos[0] * mur3 * (1. + J2rEr*(1. - 5.*z2r2)) + uEci[0];
	res[4] = -pos[1] * mur3 * (1. + J2rEr*(1. - 5.*z2r2)) + uEci[1];
	res[5] = -pos[2] * mur3 * (1. + J2rEr*(3. - 5.*z2r2)) + uEci[2];
    
    return res; 
}

//---------------------------------------------------------------------

template<typename T> AlgebraicVector<T> J2_J4dynamics(AlgebraicVector<T> xx, AlgebraicVector<T> uEci, double t, double mu, double Lsc)
{
	const double J2 = 1.08262668e-3; //{-}
	const double J3 = -2.53648e-6; //{-}
	const double J4 = -1.6233e-6; //{-}
    const double rE = 6378.137/Lsc;
    
    AlgebraicVector<T> pos(3), res(6), vel(3);
    
    pos[0] = xx[0]; pos[1] = xx[1]; pos[2] = xx[2];
    vel[0] = xx[3]; vel[1] = xx[4]; vel[2] = xx[5];
    
    T r = pos.vnorm();
    T v = vel.vnorm();
    
    res[0] = xx[3];
    res[1] = xx[4];
    res[2] = xx[5];

    T x = pos[0];
	T y = pos[1];
	T z = pos[2];
	
	T mur3 = mu/r/r/r;
	T z2r2 = z/r*z/r;

	T J2rEr = 1.5*J2*rE/r*rE/r;
	T mur7Er3 = 5./2.*mur3*J3*rE/r*rE/r*rE/r/r;
	T mur7Er4 = 15./8.*mur3*J4*rE/r*rE/r*rE/r*rE/r;

	res[3] = -pos[0] * mur3 * (1. + J2rEr*(1. - 5.*z2r2));
	res[4] = -pos[1] * mur3 * (1. + J2rEr*(1. - 5.*z2r2));
	res[5] = -pos[2] * mur3 * (1. + J2rEr*(3. - 5.*z2r2));
        
    res[3] = res[3] + (mur7Er3*x*z*(7.*z2r2-3.) + mur7Er4*x*(1.-14.*z2r2+21.*z2r2*z2r2) + uEci[0]);
	res[4] = res[4] + (mur7Er3*y*z*(7.*z2r2-3.) + mur7Er4*y*(1.-14.*z2r2+21.*z2r2*z2r2) + uEci[1]);
	res[5] = res[5] + (mur7Er3*r*r*(3./5. - 6.*z2r2+7.*z2r2*z2r2) +  mur7Er4*z*(5.-70./3.*z2r2+21.*z2r2*z2r2) + uEci[2]);
	
    return res; 
}

DACE::AlgebraicMatrix<double> stmDace(DACE::AlgebraicVector<DACE::DA>& x,int nDAVars, int nDepVars)
{
    DACE::AlgebraicMatrix<double> STM(nDepVars,nDAVars);
    for (unsigned int i = 0; i < nDAVars; i ++){
      for (unsigned int j = 0; j < nDepVars; j ++){
        STM.at(j,i) = DACE::cons(x[j].deriv(i+1)) ; 
      }}
    return STM;
}


template <typename T> AlgebraicMatrix<T> rtn2eci(const AlgebraicVector<T> & rv)
{
  AlgebraicVector<T> r(3), v(3), x(3), y(3), z(3);
  AlgebraicMatrix<T> dcm(3,3);
  r[0] = rv[0];
  r[1] = rv[1];
  r[2] = rv[2];
  v[0] = rv[3];
  v[1] = rv[4];
  v[2] = rv[5];

  z = -r; 
  z  = z/vnorm(z);

  y = DACE::cross(v,r);
  y  = y/vnorm(y);

  x =  DACE::cross(y, z);         
  x  = x/vnorm(x);

  dcm.at(0,0) = -z[0];   dcm.at(0,1) = x[0];   dcm.at(0,2) = -y[0];
  dcm.at(1,0) = -z[1];   dcm.at(1,1) = x[1];   dcm.at(1,2) = -y[1];
  dcm.at(2,0) = -z[2];   dcm.at(2,1) = x[2];   dcm.at(2,2) = -y[2];


  return dcm;
}


template<typename T> T det2(AlgebraicMatrix<T> M)
{    
   // computes the determinant of a 2x2 matrix M
    T det = M.at(0, 0) * M.at(1, 1) -  M.at(0, 1) * M.at(1, 0);

    return det;
}

template<typename T> T ConstPoC(AlgebraicVector<T> r, AlgebraicMatrix<double> P, double R){
  
    // Constant PoC on B-plane
    AlgebraicMatrix<double> P_inv(2,2);
    double det = det2(P);
    P_inv = P.inv();

    T smd = r.dot(P_inv*r);
    T PoC = R*R/(2*sqrt(det))*exp(-smd/2);
    
    return PoC;
}

template<typename T> T MaxPoC(AlgebraicVector<T> r, AlgebraicMatrix<double> P, double R){
  
    // Constant PoC on B-plane
    AlgebraicMatrix<double> P_inv(2,2);
    double det = det2(P);
    P_inv = P.inv();

    T smd = r.dot(P_inv*r);
    T PoC = R*R/(exp(1)*sqrt(det)*smd);
    
    return PoC;
}

template<typename T>T ChanPoC(AlgebraicVector<T> r, AlgebraicMatrix<double> P, double R, int order){
    
    T xi_exp = r[0];
    T zeta_exp = r[1];
    
    T sigma_xi = sqrt(P.at(0,0));
    T sigma_zeta = sqrt(P.at(1,1));
    T rho = P.at(0,1)/(sigma_xi*sigma_zeta);
    
    T u = sqr(R)/(sigma_xi*sigma_zeta*sqrt(1.0-sqr(rho)));
    T v = (sqr((xi_exp/sigma_xi)) + sqr((zeta_exp/sigma_zeta)) - 2.0*rho*(xi_exp*zeta_exp/(sigma_xi*sigma_zeta)))/(1.0-sqr(rho));
    
    T SecondLoop = 0;
    
    for (int m = 0; m <= order; m++){
        
        T FirstLoop = 0;
        for (int k = 0; k <= m; k++){
            FirstLoop = FirstLoop + pow(u,k)/(pow(2.0,k) * factorial(k));
        }
        SecondLoop = SecondLoop + pow(v,m)/(pow(2.0,m) * factorial(m)) * (1.0 - exp(-u/2.0)*FirstLoop);
    }
    T PoC = exp(-v/2.0)*SecondLoop;
    return PoC;

}

template<typename T>T AlfanoPoC(AlgebraicVector<T> r, AlgebraicMatrix<double> P, double HBR){

    T xm, zm, sigmax, sigmaz, PoC;
    AlgebraicVector<T> xx(2,2);
    T x = r[0];
    T z = r[1];
    
    T sigma_x = sqrt(P.at(0,0));
    T sigma_z = sqrt(P.at(1,1));
    T rho = P.at(0,1)/(sigma_x*sigma_z);
    T theta = 1.0/2.0*atan(2.0*rho*sigma_x*sigma_z/(sigma_x*sigma_x-sigma_z*sigma_z));
    
    if (cons(sigma_z)>cons(sigma_x)){theta = theta + atan(1.0)*2.0;}
    AlgebraicMatrix<T> R(2,2), C(2,2);
    
    R.at(0,0) = cos(theta); R.at(0,1) = sin(theta);
    R.at(1,0) = -sin(theta); R.at(1,1) = cos(theta);
    
    C = R*P*R.transpose();
    xx = R*r;
    xm = xx[0];
    zm = xx[1];

    sigmax = sqrt(C.at(0,0));
    sigmaz = sqrt(C.at(1,1));
    
   // int n = floor(cons(5.0*HBR)/min(cons(sqrt(sigmaz)),cons(vnorm(posRel))));
    int n = 30;
    
    for (int i = 1; i<n+1; i++)
    {
        T aux1 = ( zm + 2.0*HBR/n*sqrt((n-i)*i) )/(sigmaz*sqrt(2.0));
        T aux2 = (-zm + 2.0*HBR/n*sqrt((n-i)*i) )/(sigmaz*sqrt(2.0));
        T aux3 = -pow((HBR*(2.0*i-n)/n + xm ),2 )/(2.0*sigmax*sigmax);
        PoC = PoC + (erf(aux1)+erf(aux2))*exp(aux3);

    }
    PoC = HBR*2.0/(sqrt(8.0*atan(1.0)*4.0)*sigmax*n)*PoC;
    
    return PoC;

}


template<typename T> AlgebraicMatrix<T> Bplane(AlgebraicVector<T> vvs, AlgebraicVector<T> vvd){

    AlgebraicVector<T> relVel = vvs-vvd;
    AlgebraicVector<T> eta = relVel/relVel.vnorm();
    AlgebraicVector<T> CrossVel = vvd.cross(vvs);
    AlgebraicVector<T> xi = CrossVel/CrossVel.vnorm();
    AlgebraicVector<T> zeta = xi.cross(eta);
    zeta = zeta/vnorm(zeta);

    AlgebraicMatrix<T> toBplane(3,3);
    toBplane.at(0,0) = xi[0];     toBplane.at(0,1) = xi[1];     toBplane.at(0,2) = xi[2];
    toBplane.at(1,0) = eta[0];    toBplane.at(1,1) = eta[1];    toBplane.at(1,2) = eta[2];
    toBplane.at(2,0) = zeta[0];   toBplane.at(2,1) = zeta[1];   toBplane.at(2,2) = zeta[2];

    return toBplane;
}

template<typename T, typename U>
DACE::AlgebraicVector<T> RK78(const int N, DACE::AlgebraicVector<T> Y0, DACE::AlgebraicVector<T> U0,
	const U X0, const U X1, AlgebraicVector<T> (*dyn)(AlgebraicVector<T>,AlgebraicVector<T>,double,double,double), const double mu, const double Lsc,  
    const bool returnIntermediatePoints = false, const double tolerance = 1.e-11){

    double ERREST;
    double H0_i;
    double HS_i;
    double H1_i;
    const double EPS = tolerance;
    const double BS = 20. * EPS;
	
    H0_i=1e-14;
    HS_i=1e+5;
    H1_i=1e+1;	
	
	const double H0=H0_i;
	const double HS=HS_i;
	const double H1=H1_i;
	
    T Z[N][16];

    DACE::AlgebraicVector<T> Yout = Y0;
    Yout.push_back(X0);
    if (returnIntermediatePoints)
    {
        Yout.reserve(128);
    }

    DACE::AlgebraicVector<T> Y1(N);
    std::vector<double> Y1cons(N);

    double VIHMAX = 0.0;
    U X, H;
    double RFNORM, HH0, HH1;

    const double HSQR = 1.0 / 9.0;
    double A[13], C[13], D[13];
    double B[13][12];



    A[0] = 0.0; A[1] = 1.0 / 18.0; A[2] = 1.0 / 12.0; A[3] = 1.0 / 8.0; A[4] = 5.0 / 16.0; A[5] = 3.0 / 8.0;
    A[6] = 59.0 / 400.0; A[7] = 93.0 / 200.0; A[8] = 5490023248.0 / 9719169821.0; A[9] = 13.0 / 20.0; A[10] = 1201146811.0 / 1299019798.0; A[11] = 1.0;
    A[12] = 1.0;

    B[0][0] = 0.0; B[0][1] = 0.0; B[0][2] = 0.0; B[0][3] = 0.0; B[0][4] = 0.0;
    B[0][5] = 0.0; B[0][6] = 0.0; B[0][7] = 0.0; B[0][8] = 0.0; B[0][9] = 0.0;
    B[0][10] = 0.0; B[0][11] = 0.0;

    B[1][0] = 1.0 / 18.0; B[1][1] = 0.0; B[1][2] = 0.0; B[1][3] = 0.0; B[1][4] = 0.0;
    B[1][5] = 0.0; B[1][6] = 0.0; B[1][7] = 0.0; B[1][8] = 0.0; B[1][9] = 0.0;
    B[1][10] = 0.0; B[1][11] = 0.0;

    B[2][0] = 1.0 / 48.0; B[2][1] = 1.0 / 16.0; B[2][2] = 0.0; B[2][3] = 0.0; B[2][4] = 0.0;
    B[2][5] = 0.0; B[2][6] = 0.0; B[2][7] = 0.0; B[2][8] = 0.0; B[2][9] = 0.0;
    B[2][10] = 0.0; B[2][11] = 0.0;

    B[3][0] = 1.0 / 32.0; B[3][1] = 0.0; B[3][2] = 3.0 / 32.0; B[3][3] = 0.0; B[3][4] = 0.0;
    B[3][5] = 0.0; B[3][6] = 0.0; B[3][7] = 0.0; B[3][8] = 0.0; B[3][9] = 0.0;
    B[3][10] = 0.0; B[3][11] = 0.0;

    B[4][0] = 5.0 / 16.0; B[4][1] = 0.0; B[4][2] = -75.0 / 64.0; B[4][3] = 75.0 / 64.0; B[4][4] = 0.0;
    B[4][5] = 0.0; B[4][6] = 0.0; B[4][7] = 0.0; B[4][8] = 0.0; B[4][9] = 0.0;
    B[4][10] = 0.0; B[4][11] = 0.0;

    B[5][0] = 3.0 / 80.0; B[5][1] = 0.0; B[5][2] = 0.0; B[5][3] = 3.0 / 16.0; B[5][4] = 3.0 / 20.0;
    B[5][5] = 0.0; B[5][6] = 0.0; B[5][7] = 0.0; B[5][8] = 0.0; B[5][9] = 0.0;
    B[5][10] = 0.0; B[5][11] = 0.0;

    B[6][0] = 29443841.0 / 614563906.0; B[6][1] = 0.0; B[6][2] = 0.0; B[6][3] = 77736538.0 / 692538347.0; B[6][4] = -28693883.0 / 1125000000.0;
    B[6][5] = 23124283.0 / 1800000000.0; B[6][6] = 0.0; B[6][7] = 0.0; B[6][8] = 0.0; B[6][9] = 0.0;
    B[6][10] = 0.0; B[6][11] = 0.0;

    B[7][0] = 16016141.0 / 946692911.0; B[7][1] = 0.0; B[7][2] = 0.0; B[7][3] = 61564180.0 / 158732637.0; B[7][4] = 22789713.0 / 633445777.0;
    B[7][5] = 545815736.0 / 2771057229.0; B[7][6] = -180193667.0 / 1043307555.0; B[7][7] = 0.0; B[7][8] = 0.0; B[7][9] = 0.0;
    B[7][10] = 0.0; B[7][11] = 0.0;

    B[8][0] = 39632708.0 / 573591083.0; B[8][1] = 0.0; B[8][2] = 0.0; B[8][3] = -433636366.0 / 683701615.0; B[8][4] = -421739975.0 / 2616292301.0;
    B[8][5] = 100302831.0 / 723423059.0; B[8][6] = 790204164.0 / 839813087.0; B[8][7] = 800635310.0 / 3783071287.0; B[8][8] = 0.0; B[8][9] = 0.0;
    B[8][10] = 0.0; B[8][11] = 0.0;

    B[9][0] = 246121993.0 / 1340847787.0; B[9][1] = 0.0; B[9][2] = 0.0; B[9][3] = -37695042795.0 / 15268766246.0; B[9][4] = -309121744.0 / 1061227803.0;
    B[9][5] = -12992083.0 / 490766935.0; B[9][6] = 6005943493.0 / 2108947869.0; B[9][7] = 393006217.0 / 1396673457.0; B[9][8] = 123872331.0 / 1001029789.0; B[9][9] = 0.0;
    B[9][10] = 0.0; B[9][11] = 0.0;

    B[10][0] = -1028468189.0 / 846180014.0; B[10][1] = 0.0; B[10][2] = 0.0; B[10][3] = 8478235783.0 / 508512852.0; B[10][4] = 1311729495.0 / 1432422823.0;
    B[10][5] = -10304129995.0 / 1701304382.0; B[10][6] = -48777925059.0 / 3047939560.0; B[10][7] = 15336726248.0 / 1032824649.0; B[10][8] = -45442868181.0 / 3398467696.0; B[10][9] = 3065993473.0 / 597172653.0;
    B[10][10] = 0.0; B[10][11] = 0.0;

    B[11][0] = 185892177.0 / 718116043.0; B[11][1] = 0.0; B[11][2] = 0.0; B[11][3] = -3185094517.0 / 667107341.0; B[11][4] = -477755414.0 / 1098053517.0;
    B[11][5] = -703635378.0 / 230739211.0; B[11][6] = 5731566787.0 / 1027545527.0; B[11][7] = 5232866602.0 / 850066563.0; B[11][8] = -4093664535.0 / 808688257.0; B[11][9] = 3962137247.0 / 1805957418.0;
    B[11][10] = 65686358.0 / 487910083.0; B[11][11] = 0.0;

    B[12][0] = 403863854.0 / 491063109.0; B[12][1] = 0.0; B[12][2] = 0.0; B[12][3] = -5068492393.0 / 434740067.0; B[12][4] = -411421997.0 / 543043805.0;
    B[12][5] = 652783627.0 / 914296604.0; B[12][6] = 11173962825.0 / 925320556.0; B[12][7] = -13158990841.0 / 6184727034.0; B[12][8] = 3936647629.0 / 1978049680.0; B[12][9] = -160528059.0 / 685178525.0;
    B[12][10] = 248638103.0 / 1413531060.0; B[12][11] = 0.0;

    C[0] = 14005451.0 / 335480064.0; C[1] = 0.0; C[2] = 0.0; C[3] = 0.0; C[4] = 0.0; C[5] = -59238493.0 / 1068277825.0;
    C[6] = 181606767.0 / 758867731.0; C[7] = 561292985.0 / 797845732.0; C[8] = -1041891430.0 / 1371343529.0; C[9] = 760417239.0 / 1151165299.0; C[10] = 118820643.0 / 751138087.0; C[11] = -528747749.0 / 2220607170.0;
    C[12] = 1.0 / 4.0;

    D[0] = 13451932.0 / 455176623.0; D[1] = 0.0; D[2] = 0.0; D[3] = 0.0; D[4] = 0.0; D[5] = -808719846.0 / 976000145.0;
    D[6] = 1757004468.0 / 5645159321.0; D[7] = 656045339.0 / 265891186.0; D[8] = -3867574721.0 / 1518517206.0; D[9] = 465885868.0 / 322736535.0; D[10] = 53011238.0 / 667516719.0; D[11] = 2.0 / 45.0;
    D[12] = 0.0;

    for (int i = 0; i < N; i++)
    {
        Z[i][0] = Y0[i];
        Z[i][1] = 0.0;
    }

    double sig=1.;
    if(X1<X0){
	sig=-1.;
    }
    H=sig*abs(HS);
    HH0 = abs(H0); HH1 = abs(H1);
    X = X0; RFNORM = 0.0; ERREST = 0.0;
    
    while (std::abs(DACE::cons(X)-DACE::cons(X1)) > 0.0)
	{
        // compute new stepsize
		if (RFNORM > 0)
        {
            H = H*std::min(4.0, exp(HSQR*log(EPS / RFNORM)));
        }
        if (abs(H)>abs(HH1))
        {
            H = sig*HH1;
        }
        else if (abs(H)<abs(HH0)*0.99)
        {
            H = sig*HH0;
            std::cout << "--- WARNING, MINIMUM STEPSIZE REACHED IN RK" << std::endl;
        }
        if ((DACE::cons(X) + DACE::cons(H) - DACE::cons(X1))*DACE::cons(H)>0)
        {
            H = X1 - X;
        }
		for (int j = 0; j<13; j++) {
            for (int i = 0; i<N; i++)
            {
                Y0[i] = 0.0; // EVALUATE RHS AT 13 POINTS

                for (int k = 0; k<j; k++)
                {
                    Y0[i] = Y0[i] + Z[i][k + 3] * B[j][k];
                }

                Y0[i] = H*Y0[i] + Z[i][0];
            }
            Y1 = dyn(Y0, U0, X + H*A[j], mu, Lsc);

            for (int i = 0; i<N; i++)
            {
                Z[i][j + 3] = Y1[i];
            }
        }
	
        for (int i = 0; i<N; i++) {

            Z[i][1] = 0.0; Z[i][2] = 0.0; // EXECUTE 7TH,8TH ORDER STEPS

            for (int j = 0; j<13; j++)
            {
                Z[i][1] = Z[i][1] + Z[i][j + 3] * D[j];
                Z[i][2] = Z[i][2] + Z[i][j + 3] * C[j];
            }

            Y1[i] = (Z[i][2] - Z[i][1])*H;
            Z[i][2] = Z[i][2] * H + Z[i][0];
        }



        for (int i = 0; i<N; i++)
        {
            Y1cons[i] = DACE::cons(Y1[i]);
        }
	
        RFNORM = normtmp(N, Y1cons); // ESTIMATE ERROR AND DECIDE ABOUT BACKSTEP

        if ((RFNORM>BS) && (abs(H / H0)>1.2))
        {
            H = H / 3.0;
            RFNORM = 0;
        }
        else
        {
            for (int i = 0; i<N; i++)
            {
                Z[i][0] = Z[i][2];
            }
            X = X + H;
            VIHMAX = std::max(VIHMAX, DACE::cons(H));
            ERREST = ERREST + RFNORM;

            if (returnIntermediatePoints)
            {
                for (int i = 0; i<N; i++)
                {
                    Yout.push_back(Z[i][0]);
                }
                Yout.push_back(X);
            }
        }
    }

    if (returnIntermediatePoints)
    {
        return Yout;
    }

    // Return final state and time
    for (int i = 0; i<N; i++)
    {
        Y1[i] = Z[i][0];
    }
   // Y1.push_back(X);

    return Y1;

}
}