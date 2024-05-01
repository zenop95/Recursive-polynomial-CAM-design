#ifdef sada_HEOSATDYNAMICS_H

#include <cmath>

// DACE
#include <dace/dace.h>

// Dynamics
#include "Conversions.h"

// CSPICE
#include <cspice/SpiceUsr.h>

// SA
#include "constants.h"


const double DEG = M_PI / 180.0;


template<typename T>
double HEOSATDynamics<T>::initheosat(const double aa, const double dj, const double tini, 
	const double tfin, double *t, double *tf, const T i_Bfactor, const T i_SRPC)
{
	
	const double Zn[11] = { 0.0, 0.0, -1.0826356665511e-3, 2.5324736913329e-6, 1.6199743057822e-6,
		2.2790512608210e-7, -5.4061678994023e-7, 3.5052292563209e-7, 2.0401676701049e-7,
		1.2215024545413e-7, 2.4434298429376e-7 };

	m_we = 0.7292115373194e-4;  // Earth rotational speed
	m_GM = 398600.4415;  // Earth gravitational parameter
	m_Req = 6378.1363;  // Earth equatorial radius

	for (int i = 0; i < 11; ++i)
	{
		m_C[i] = Zn[i];
	}
	m_fl = 1.0 / 300.0;
    
    m_c22 = 1.5745764195628e-6;
    m_s22 = -9.0386794572580e-7;

	m_ul = aa;  // Unit length
	m_ut = aa*std::sqrt(aa / m_GM);  // Unit time in seconds
	m_utd = m_ut / (24 * 60 * 60.);  // Unit time in days
	m_mu = 1.0;
	m_al = m_Req / m_ul;  // Scaled Earth equatorial radius

	m_ns = (2 * M_PI / (365.256363004 * 24. * 60. * 60.))*m_ut;  // Sun mean motion?
	const double AU = 149597870.700;
	m_as = 1.000001018*AU / m_ul;  // Sun semi-major axis

	m_nm = (2. * M_PI) / (24. * 60. * 60.)*(1. / 27.321582)*m_ut;  // Moon mean motion?
	m_bet = 1.0 / 82.2845; // = mMoon/(mMoon+mEarth) = 81.2845+1
	m_am = 384400.0 / m_ul;  // Moon semi-major axis

	m_Bfactor = i_Bfactor; // Bfactor   A_drag/m*Cd
	m_SRPC = i_SRPC; //solar radiation pressure coefficient A_srp/m*(1+eps)

	*t = tini / m_ut;  // Scaled initial time
	*tf = tfin / m_ut;  // Scaled final time

	m_t0 = *t;  // Scaled initial time
	m_dj0 = dj;  // Initial Julian date

	return m_ut;
}

template<typename T>
double HEOSATDynamics<T>::diajul(int giorno, int mese, int anno, double ora)
{
	int iy, im, ib, k1, k2, k3;

	if (mese <= 2){
		iy = anno - 1;
		im = mese + 12;
	}
	else{
		iy = anno;
		im = mese;
	}

	if (anno > 1582){
		ib = iy / 400 - iy / 100;
	}
	else{
		ib = -2;
		if (anno == 1582){
			if (mese > 10){
				ib = iy / 400 - iy / 100;
			}
			else if ((mese == 10) && (giorno >= 15)){
				ib = iy / 400 - iy / 100;
			}
		}
	}

	k1 = 365.25*iy;
	k2 = 30.6001*(im + 1);
	k3 = k1 + k2 + ib - 679004 + giorno;

	return(2400000.5 + k3 + ora / 24.);
}

template<typename T>
double HEOSATDynamics<T>::ut2dj(double t)
{
	return(m_dj0 + (t - m_t0)*m_utd);
}

template<typename T>
double HEOSATDynamics<T>::keplereq(double M, double e, int j)
{

	double t, u;
	int i;

	t = e*std::sin(M) / (1. - e*std::cos(M));
	i = 0;
	u = M;
	while ((std::fabs(t) > 1.e-12) && (i < j)){
		i = i + 1;
		u = u + t;
		t = (M - u + e*std::sin(u)) / (1. - e*std::cos(u));
	}

	return(u + t);
}

template<typename T>
double HEOSATDynamics<T>::tocent(double tt)
{
    /*
     Compute centuries since epoch J2000 (2000 January 1.5)
     */
	const double j2000 = 2451545.0;

	return((tt - j2000) / 36525.0);
}

template<typename T>
double HEOSATDynamics<T>::tsmut(double jd)
{
	double frac, tu, aux;

    // Fraction of day since midnight 0:00
	frac = jd - 0.50 - (int)(jd - 0.50);
    // Centuries since epoch J2000
	tu = tocent(jd - frac);
    // Mean sidereal time in radians
	aux = ((24110.54841 + tu*(8640184.812866 + tu*(0.093104 - tu*0.0000062)))*15.)*1.745329251994330e-2 / 3600. + 
		(1.002737909350795 + tu*(5.9006e-11 - tu*5.9e-15))*frac*24.*15.*1.745329251994330e-2;

	return(fmod(aux, 2.*M_PI));
}

template<typename T>
double HEOSATDynamics<T>::tsmut2(double jd)
{
	// Astronomy with Your Personal Computer by Peter Duffett - Smith
	double temp = jd - 2415019.5;
	double tu = temp / 365.25;
	double year = 1900. + std::floor(tu);
	double dj = jd - 2415021.;
	double TT = (dj / 36525.) - 1.;
	double R1 = 6.697374558 + (2400.*(TT - ((year - 2000.) / 100.)));
	double R0 = (5.13366e-2*TT) + (2.586222e-5*TT*TT) - (1.722E-9*TT*TT*TT);
	double T0 = R0 + R1;
	double UT = (jd - std::floor(jd))*24.;
	return fmod(((UT*1.002737908) + T0) / 24.*2.*M_PI, 2.*M_PI) - 2.*M_PI;
}

template<typename T>
double HEOSATDynamics<T>::tsmut3(double jd)
{
	// Montenbruck	
	return std::fmod(((jd - 2451545.)*360.9856473 + 280.4606) / 180.*M_PI, 2.*M_PI);
}

///// SUN AND MOON COEFFICIENTS /////

template<typename T>
void HEOSATDynamics<T>::coefsolun5(double *css, double *ccs, double *ks, double *csm, double *ccm, double *km)
{
	// extern double m_Rs, m_xs, m_ys, m_zs, m_rm, m_xl, m_yl, m_zl;
	// extern double m_bet, m_ns, m_nm;

	css[0] = -2 * m_xs*m_ys;
	css[1] = m_xs*m_zs;
	css[2] = 0;

	ccs[0] = m_xs*m_xs - m_ys*m_ys;
	ccs[1] = m_ys*m_zs;
	ccs[2] = 1 - 3 * m_zs*m_zs;

	*ks = (m_ns*m_ns)*std::pow(1 / m_Rs, 3);

	csm[0] = -2 * m_xl*m_yl;
	csm[1] = m_xl*m_zl;
	csm[2] = 0;

	ccm[0] = m_xl*m_xl - m_yl*m_yl;
	ccm[1] = m_yl*m_zl;
	ccm[2] = 1 - 3 * m_zl*m_zl;
	for (int i = 1; i <= 2; ++i)
	{
		csm[i + 2] = -csm[-i + 2];
		css[i + 2] = -css[-i + 2];
		ccm[i + 2] = ccm[-i + 2];
		ccs[i + 2] = ccs[-i + 2];
	}

	*km = m_bet*(m_nm*m_nm)*std::pow(1 / m_rm, 3);
}

template<typename T>
void HEOSATDynamics<T>::coefsun6(double *cm, double *sm, double *km)
{
	// extern double m_Rs, m_xs, m_ys, m_zs;
	// extern double m_bet, m_n, m_as, m_ns;
	double xs2, ys2, zs2;

	xs2 = m_xs*m_xs;
	ys2 = m_ys*m_ys;
	zs2 = m_zs*m_zs;


	cm[0] = m_xs*(xs2 - 3 * ys2);
	cm[1] = -2 * m_xs*m_ys*m_zs;
	cm[2] = m_xs*(1 - 5 * zs2);
	cm[3] = 0;

	sm[0] = m_ys*(3 * xs2 - ys2);
	sm[1] = (xs2 - ys2)*m_zs;
	sm[2] = m_ys*(1 - 5 * zs2);
	sm[3] = m_zs*(3 - 5 * zs2);

	for (int i = 1; i <= 3; ++i)
	{
		cm[i + 3] = -cm[-i + 3];
		sm[i + 3] = sm[-i + 3];
	}

	*km = (m_ns*m_ns)*(1 / m_as)*std::pow(1 / m_Rs, 4);
}

template<typename T>
void HEOSATDynamics<T>::coeflun6(double *cm, double *sm, double *km)
{
	// extern double m_rm, m_xl, m_yl, m_zl;
	// extern double m_bet, m_n, m_am, m_nm;
	double xl2, yl2, zl2;

	xl2 = m_xl*m_xl;
	yl2 = m_yl*m_yl;
	zl2 = m_zl*m_zl;


	cm[0] = m_xl*(xl2 - 3 * yl2);
	cm[1] = -2 * m_xl*m_yl*m_zl;
	cm[2] = m_xl*(1 - 5 * zl2);
	cm[3] = 0;

	sm[0] = m_yl*(3 * xl2 - yl2);
	sm[1] = (xl2 - yl2)*m_zl;
	sm[2] = m_yl*(1 - 5 * zl2);
	sm[3] = m_zl*(3 - 5 * zl2);

	for (int i = 1; i <= 3; ++i)
	{
		cm[i + 3] = -cm[-i + 3];
		sm[i + 3] = sm[-i + 3];
	}

	*km = m_bet*(m_nm*m_nm)*(1 / m_am)*std::pow(1 / m_rm, 4);
}

template<typename T>
void HEOSATDynamics<T>::coeflun7(double *cm, double *sm, double *km)
{
	// extern double m_rm, m_xl, m_yl, m_zl;
	// extern double m_bet, m_n, m_am, m_nm;
	double xl2, yl2, zl2, zl4;

	xl2 = m_xl*m_xl;
	yl2 = m_yl*m_yl;
	zl2 = m_zl*m_zl;
	zl4 = zl2*zl2;


	cm[0] = -4 * xl2*yl2 + (xl2 - yl2)*(xl2 - yl2);
	cm[1] = m_yl*(3 * xl2 - yl2)*m_zl;
	cm[2] = (xl2 - yl2)*(1 - 7 * zl2);
	cm[3] = m_yl*m_zl*(3 - 7 * zl2);
	cm[4] = 3 - 30 * zl2 + 35 * zl4;

	sm[0] = 4 * m_xl*m_yl*(-xl2 + yl2);
	sm[1] = m_xl*(xl2 - 3 * yl2)*m_zl;
	sm[2] = -2 * m_xl*m_yl*(1 - 7 * zl2);
	sm[3] = m_xl*m_zl*(3 - 7 * zl2);
	sm[4] = 0;

	for (int i = 1; i <= 4; ++i)
	{
		cm[i + 4] = cm[-i + 4];
		sm[i + 4] = -sm[-i + 4];
	}

	*km = m_bet*(m_nm*m_nm)*(1 / (m_am*m_am))*std::pow(1 / m_rm, 5);
}

template<typename T>
void HEOSATDynamics<T>::coeflun8(double *cm, double *sm, double *km)
{
	// extern double m_rm, m_xl, m_yl, m_zl;
	// extern double m_bet, m_n, m_am, m_nm;
	double xl2, yl2, zl2, xl4, yl4, zl4;

	xl2 = m_xl*m_xl;
	xl4 = xl2*xl2;
	yl2 = m_yl*m_yl;
	yl4 = yl2*yl2;
	zl2 = m_zl*m_zl;
	zl4 = zl2*zl2;


	cm[0] = m_xl*(xl4 - 10 * xl2*yl2 + 5 * yl4);
	cm[1] = -4 * m_xl*m_yl*(xl2 - yl2)*m_zl;
	cm[2] = m_xl*(xl2 - 3 * yl2)*(1 - 9 * zl2);
	cm[3] = -2 * m_xl*m_yl*m_zl*(1 - 3 * zl2);
	cm[4] = m_xl*(1 - 14 * zl2 + 21 * zl4);
	cm[5] = 0;

	sm[0] = m_yl*(5 * xl4 - 10 * xl2*yl2 + yl4);
	sm[1] = (xl4 - 6 * xl2*yl2 + yl4)*m_zl;
	sm[2] = m_yl*(3 * xl2 - yl2)*(1 - 9 * zl2);
	sm[3] = (xl2 - yl2)*m_zl*(1 - 3 * zl2);
	sm[4] = m_yl*(1 - 14 * zl2 + 21 * zl4);
	sm[5] = m_zl*(15 - 70 * zl2 + 63 * zl4);

	for (int i = 1; i <= 5; ++i)
	{
		cm[i + 5] = -cm[-i + 5];
		sm[i + 5] = sm[-i + 5];
	}

	*km = m_bet*(m_nm*m_nm)*(1 / std::pow(m_am, 3))*std::pow(1 / m_rm, 6);
}

template<typename T>
void HEOSATDynamics<T>::coeflun9(double *cm, double *sm, double *km)
{
	// extern double m_rm, m_xl, m_yl, m_zl;
	// extern double m_bet, m_n, m_am, m_nm;
	double xl2, yl2, zl2, xl4, yl4, zl4, zl6;

	xl2 = m_xl*m_xl;
	xl4 = xl2*xl2;
	yl2 = m_yl*m_yl;
	yl4 = yl2*yl2;
	zl2 = m_zl*m_zl;
	zl4 = zl2*zl2;
	zl6 = zl4*zl2;


	cm[0] = (xl2 - yl2)*(xl4 - 14 * xl2*yl2 + yl4);
	cm[1] = m_yl*(5 * xl4 - 10 * xl2*yl2 + yl4)*m_zl;
	cm[2] = (xl4 - 6 * xl2*yl2 + yl4)*(1 - 11 * zl2);
	cm[3] = m_yl*(-3 * xl2 + yl2)*m_zl*(3 - 11 * zl2);
	cm[4] = (xl2 - yl2)*(1 - 18 * zl2 + 33 * zl4);
	cm[5] = m_yl*m_zl*(5 - 30 * zl2 + 33 * zl4);
	cm[6] = 5 - 105 * zl2 + 315 * zl4 - 231 * zl6;

	sm[0] = 2 * m_xl*m_yl*(xl2 - 3 * yl2)*(3 * xl2 - yl2);
	sm[1] = -(m_xl*(xl4 - 10 * xl2*yl2 + 5 * yl4)*m_zl);
	sm[2] = 4 * m_xl*m_yl*(xl2 - yl2)*(1 - 11 * zl2);
	sm[3] = m_xl*(xl2 - 3 * yl2)*m_zl*(3 - 11 * zl2);
	sm[4] = 2 * m_xl*m_yl*(1 - 18 * zl2 + 33 * zl4);
	sm[5] = -(m_xl*m_zl*(5 - 30 * zl2 + 33 * zl4));
	sm[6] = 0;

	for (int i = 1; i <= 6; ++i)
	{
		cm[i + 6] = cm[-i + 6];
		sm[i + 6] = -sm[-i + 6];
	}

	*km = m_bet*(m_nm*m_nm)*(1 / std::pow(m_am, 4))*std::pow(1 / m_rm, 7);
}

///// DRAG /////

template<typename T>
T HEOSATDynamics<T>::harrispriester(const T& h, const T *ca)
{
	//Using std and DACE functions
	using std::sqrt;  using DACE::sqrt;
	using std::exp;   using DACE::exp;
	using std::log;   using DACE::log;
	
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	using std::asin;  using DACE::asin;
	using std::atan;  using DACE::atan;
	
	// extern double m_xs, m_ys, m_zs;
	T q = sqrt(ca[0] * ca[0] + ca[1] * ca[1] + ca[2] * ca[2]);
	const T u = ca[0] / q;
	const T v = ca[1] / q;
	const T w = ca[2] / q;

	int i = 0;
	while (DACE::cons(h) >= hi[i])
	{
		i = i + 1;
	}
	const int j = i - 1;
	const double Hmi = (hi[j] - hi[i]) / log(pm[i] / pm[j]);
	const double Hma = (hi[j] - hi[i]) / log(qM[i] / qM[j]);
	const T p = pm[j] * exp((hi[j] - h) / Hmi);
	q = qM[j] * exp((hi[j] - h) / Hma);
	const double al = atan2_mod(m_ys, m_xs) + 30.*M_PI / 180.0;
	const double ds = atan(m_zs / sqrt(m_xs*m_xs + m_ys*m_ys));
	T cny2 = u*cos(ds)*cos(al) + v*cos(ds)*sin(al) + w*sin(ds);
	cny2 = (0.5 + 0.5*cny2)*sqrt(0.5 + 0.5*cny2);

	return(p + (q - p)*cny2);
}

template<typename T>
void HEOSATDynamics<T>::aeiwhm0new(const T& n, const T& a, const T& e, const T& y, const T& ci, const T& si, const T& g, const T& r,
	const T& th, const T& p, const double we, const T nd[3], T *a0, T *e0, T *i0, T *W0, T *g0, T *M0)
{
	//Using std and DACE functions
	using std::sqrt;  using DACE::sqrt;
	
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
	if (std::abs(DACE::cons(nd[1])) > 0. || std::abs(DACE::cons(nd[0])) > 0. || std::abs(DACE::cons(nd[2])) > 0.)
	{
		const T f = th - g;
		const T cf = cos(f);
		const T sf = sin(f);

		T h = sqrt(p);
		T b = a*y;

		T ar = nd[0];
		T at = nd[1];
		T ah = nd[2];

		*a0 = (2. * a*a / h*(e*sf*ar + p / r*at));
		*e0 = (1. / h*(p*sf*ar + ((p + r)*cf + r*e)*at));
		*i0 = (r*cos(th) / h*ah);
		*W0 = (r*sin(th) / (h*si)*ah);
		*g0 = (1. / (h*e)*(-p*cf*ar + (p + r)*sf*at) - r*sin(th)*ci / (h*si)*ah);
		*M0 = b / (a*h*e)*((p*cf - 2. * r*e)*ar - (p + r)*sf*at);
	}
	else
	{
		*a0 = 0;
		*e0 = 0;
		*i0 = 0;
		*W0 = 0;
		*g0 = 0;
		*M0 = 0;
	}
}

template<typename T>
void HEOSATDynamics<T>::atmdesnew(T *kep, const T& si, const T& ci, T *r, T *th, T nd[3])
{
	//Using std and DACE functions
	using std::pow;   using DACE::pow;
	using std::sqrt;  using DACE::sqrt;
	
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
	T car[6], nu, RR, ZZ, NN;

	kep2pol(kep, r, th, &nu, &RR, &ZZ, &NN, m_GM);

	T a = kep[0];
	T e = kep[1];
	T argPer = kep[4];
	T theta = mean2trueAnomaly(kep[5], e);
	*r = *r + -m_C[2] * pow(m_Req, 2) / (4. * (1 - pow(e, 2))*a)  *  (pow(si, 2) * cos(2. * (theta + argPer)) + (3. * pow(si, 2) - 2.) * (1. + e*cos(theta) / (1. + sqrt(1. - pow(e, 2))) + 2. * sqrt(1 - pow(e, 2)) / (1. + e*cos(theta))));

	const T h = *r - m_Req*(1. - m_fl) / sqrt(1. - m_fl*(2. - m_fl)*(1. - pow((si*sin(*th)), 2)));
	if ((DACE::cons(h) >= 100.0) && (DACE::cons(h) <= 1000.0))
	{
		kep2cart(kep, car, m_GM);

		T d = harrispriester(h, car); // g/km^3

		const T vz = (ZZ / *r) - *r*m_we*ci;
		const T vh = *r*m_we*si*cos(*th);

		nd[0] = -0.5*d*(m_Bfactor)*sqrt(RR*RR + vz*vz + vh*vh)*m_ut * RR*m_ut / m_ul; // A factor 1e-9 is missing here and is added later to avoid numerical issues.
		nd[1] = -0.5*d*(m_Bfactor)*sqrt(RR*RR + vz*vz + vh*vh)*m_ut * vz*m_ut / m_ul; // A factor 1e-9 is missing here and is added later to avoid numerical issues.
		nd[2] = -0.5*d*(m_Bfactor)*sqrt(RR*RR + vz*vz + vh*vh)*m_ut * vh*m_ut / m_ul; // A factor 1e-9 is missing here and is added later to avoid numerical issues.
	}
	else
	{
		nd[0] = 0.0;
		nd[1] = 0.0;
		nd[2] = 0.0;
	}
	
	*r = *r / m_ul;
}

template<typename T>
DACE::AlgebraicVector<T> HEOSATDynamics<T>::computeDragMeanElementRatesTrueAnomaly(const T& n, const T& a, const T& e, const T& y, const T& ci, const T& si, const T& W, const T& g)
{
	
	// extern double m_ul, m_ut;
	const double we = 0.7292115373194e-4*m_ut;

	T kep[6];
	kep[0] = a*m_ul;
	kep[1] = e;
	kep[2] = atan2_mod(si, ci);
	kep[3] = W;
	kep[4] = g;
	kep[5] = 0.;

	DACE::AlgebraicVector<T> meanRates(6, 0.0);
	const T p = a*y*y;

	T r, th, nd[3];
	atmdesnew(kep, si, ci, &r, &th, nd);
	T a0, e0, i0, W0, g0, M0;
	aeiwhm0new(n, a, e, y, ci, si, g, r, th, p, we, nd, &a0, &e0, &i0, &W0, &g0, &M0);
	T meanAnomaly0 = T(0.0);

	T a1, e1, i1, W1, g1, M1;
	const double k = 1.0 * M_PI / 180.0;
	for (int j = 1; j <= 360; ++j)
	{
		const double theta = k*(double)j;
		const T meanAnomaly1 = true2meanAnomaly(T(theta), e);
		kep[5] = meanAnomaly1;
		const T deltaMeanAnomaly = DACE::cons(meanAnomaly1) > DACE::cons(meanAnomaly0) ? meanAnomaly1 - meanAnomaly0 : meanAnomaly1 - meanAnomaly0 + 2.*M_PI;
		atmdesnew(kep, si, ci, &r, &th, nd);
		aeiwhm0new(n, a, e, y, ci, si, g, r, th, p, we, nd, &a1, &e1, &i1, &W1, &g1, &M1);
		meanRates[0] = meanRates[0] + deltaMeanAnomaly*(a1 + a0) / 2.;
		meanRates[1] = meanRates[1] + deltaMeanAnomaly*(e1 + e0) / 2.;
		meanRates[2] = meanRates[2] + deltaMeanAnomaly*(i1 + i0) / 2.;
		meanRates[3] = meanRates[3] + deltaMeanAnomaly*(W1 + W0) / 2.;
		meanRates[4] = meanRates[4] + deltaMeanAnomaly*(g1 + g0) / 2.;
		meanRates[5] = meanRates[5] + deltaMeanAnomaly*(M1 + M0) / 2.;
		a0 = a1;
		e0 = e1;
		i0 = i1;
		W0 = W1;
		g0 = g1;
		M0 = M1;
		meanAnomaly0 = meanAnomaly1;
	}
	const double idp = 0.5 / M_PI;

	return idp*meanRates;
}

template<typename T>
void HEOSATDynamics<T>::drag(const double t, const T& n, const T& a, const T& e, const T& y, const T& ci, const T& si, const T& W, const T& g, const T& l, DACE::AlgebraicVector<T>& dx)
{

	using std::sqrt;  using DACE::sqrt;
	
	DACE::AlgebraicVector<T> elementRates = computeDragMeanElementRatesTrueAnomaly(n, a, e, y, ci, si, W, g);
	T as = elementRates[0];
	T es = elementRates[1];
	T is = elementRates[2];
	T Ws = elementRates[3];
	T gs = elementRates[4];
	T Ms = elementRates[5];

	dx[0] = Ms;
	dx[1] = gs;
	dx[2] = Ws;
	dx[3] = (as) / (2.*sqrt(a));
	dx[4] = y*dx[3] - sqrt(1.*a)*(e / y)*(es);
	dx[5] = ci*dx[4] - sqrt(1.*a)*y*(is)*si;

	dx = dx*1.0e-9; // Factor 1e-9 is added here that was previously neglected in atmdesnew.
}

///// ECCENTRICITY AND INCLINATION FUNCTIONS /////

template<typename T>
void HEOSATDynamics<T>::eccincpol5(const T& e, const T& y, const T& c, const T& s)
{
	using std::pow;   using DACE::pow;
	
	const T e2 = e*e;
	const T y2 = 1. - e2;
	const T s2 = s*s;
	const T c2 = 1. - s2;

	m_ep5[1][0] = 2 + 3 * e2;
	m_ep5[1][1] = e2;
	m_ep5[2][0] = m_ep5[1][0];
	m_ep5[2][1] = m_ep5[1][1];
	m_ep5[3][0] = 14 + 6 * e2;
	m_ep5[3][1] = 1 + e2;
	m_ep5[4][0] = 1 / y;
	m_ep5[4][1] = m_ep5[4][0];
	m_ep5[5][0] = m_ep5[1][0] / y;
	m_ep5[5][1] = m_ep5[1][1] / y;

	m_ip5[1][0][0] = -(3 / 32.0)*s2;
	m_ip5[1][0][1] = -(3 / 8.0)*c*s;
	m_ip5[1][0][2] = (1 / 16.0)*(1 - 3 * c2);
	m_ip5[1][0][3] = m_ip5[1][0][1];
	m_ip5[1][0][4] = m_ip5[1][0][0];
	m_ip5[1][1][0] = -(15 / 16.0)*pow(1 - c, 2);
	m_ip5[1][1][1] = -(15 / 4.0)*(1 - c)*s;
	m_ip5[1][1][2] = -(15 / 8.0)*s2;
	m_ip5[1][1][3] = (15 / 4.0)*(1 + c)*s;
	m_ip5[1][1][4] = -(15 / 16.0)*pow(1 + c, 2);
	for (int i = 0; i <= 1; ++i){
		for (int j = -2; j <= 2; ++j){
			m_ip5[2][i][j + 2] = (j / pow(2.0, i))*m_ip5[1][i][j + 2];
			m_ip5[3][i][j + 2] = m_ip5[1][i][j + 2];
		}
	}

	m_ip5[4][0][0] = -(3 / 16.0)*(5 * c2 - 3 * y2);
	m_ip5[4][0][1] = -(3 / 8.0)*(c*(2 + 3 * e2 - 10 * s2)) / s;
	m_ip5[4][0][2] = (3 / 8.0)*(5 * c2 - y2);
	m_ip5[4][0][3] = m_ip5[4][0][1];
	m_ip5[4][0][4] = m_ip5[4][0][0];
	m_ip5[4][1][0] = (15 / 16.0)*(1 - c)*(y2 - c);
	m_ip5[4][1][1] = -(15 / 8.0)*(s*((2 + c)*e2 - 2 * s2)) / (1 + c);
	m_ip5[4][1][2] = (15 / 8.0)*(y2 - c2);
	m_ip5[4][1][3] = (15 / 8.0)*(s*((2 - c)*e2 - 2 * s2)) / (1 - c);
	m_ip5[4][1][4] = (15 / 16.0)*(1 + c)*(y2 + c);

	m_ip5[5][0][0] = (3 / 16.0)*c;
	m_ip5[5][0][1] = -(3 / 8.0)*(1 - 2 * c2) / s;
	m_ip5[5][0][2] = -(3 / 8.0)*c;
	m_ip5[5][0][3] = m_ip5[5][0][1];
	m_ip5[5][0][4] = m_ip5[5][0][0];
	m_ip5[5][1][0] = (15 / 16.0)*(1 - c);
	m_ip5[5][1][1] = (15 / 8.0)*((1 + 2 * c)*s) / (1 + c);
	m_ip5[5][1][2] = (15 / 8.0)*c;
	m_ip5[5][1][3] = (15 / 8.0)*((1 - 2 * c)*s) / (1 - c);
	m_ip5[5][1][4] = -(15 / 16.0)*(1 + c);
}

template<typename T>
void HEOSATDynamics<T>::eccincpol6(const T& e, const T& y, const T& c, const T& s)
{
	T e2, s2, c2, c3, c4, cn, cn2, cn3, cp, cp2, cp3, y2, y4;

	e2 = e*e;
	s2 = s*s;
	c2 = c*c;
	c3 = c2*c;
	c4 = c3*c;
	y2 = y*y;
	y4 = y2*y2;
	cn = 1 - c;
	cn2 = cn*cn;
	cn3 = cn2*cn;
	cp = 1 + c;
	cp2 = cp*cp;
	cp3 = cp2*cp;

	m_ep6[1][0] = e*(4 + 3 * e2);
	m_ep6[1][1] = e*e2;
	m_ep6[2][0] = m_ep6[1][0];
	m_ep6[2][1] = m_ep6[1][1];
	m_ep6[3][0] = 4 / e + (29 + 9 * e2)*e;
	m_ep6[3][1] = e*(1 + e2);
	m_ep6[4][0] = 1 / (e*y);
	m_ep6[4][1] = e / y;
	m_ep6[5][0] = m_ep6[1][0] / y;
	m_ep6[5][1] = m_ep6[1][1] / y;

	m_ip6[1][0][0] = (75 * cn*s2) / 512.;
	m_ip6[1][0][1] = (-75 * cn*(1 + 3 * c)*s) / 256.;
	m_ip6[1][0][2] = (-15 * cn*(1 - 10 * c - 15 * c2)) / 512.;
	m_ip6[1][0][3] = (15 * (1 - 5 * c2)*s) / 128.;
	m_ip6[1][0][4] = (15 * cp*(1 + 10 * c - 15 * c2)) / 512.;
	m_ip6[1][0][5] = (-75 * (1 - 3 * c)*cp*s) / 256.;
	m_ip6[1][0][6] = (-75 * cp*s2) / 512.;
	m_ip6[1][1][0] = (525 * cn3) / 512.;
	m_ip6[1][1][1] = (-1575 * cn2*s) / 256.;
	m_ip6[1][1][2] = (1575 * cn*s2) / 512.;
	m_ip6[1][1][3] = (-525 * s*s2) / 128.;
	m_ip6[1][1][4] = (-1575 * cp*s2) / 512.;
	m_ip6[1][1][5] = (-1575 * cp2*s) / 256.;
	m_ip6[1][1][6] = (-525 * cp3) / 512.;

	for (int i = 0; i <= 1; ++i){
		for (int j = -3; j <= 3; ++j){
			m_ip6[2][i][j + 3] = (j / (2.0*i + 1))*m_ip6[1][i][j + 3];
			m_ip6[3][i][j + 3] = -m_ip6[1][i][j + 3];
		}
	}

	m_ip6[4][0][0] = (-75 * cn*(7 * c*(1 + 3 * c) - (13 + 10 * c + 17 * c2)*y2 + 3 * (3 + c)*y4)) / 512.;
	m_ip6[4][0][1] = (-75 * s*(7 * c*(2 - 5 * c - 9 * c2) + (13 + 19 * c + 37 * c2 + 51 * c3)*y2 - 3 * (3 + 7 * c + 2 * c2)*y4)) / (256.*cp);
	m_ip6[4][0][2] = (15 * (7 * c*(11 + 10 * c - 45 * c2) - (13 - 33 * c + 35 * c2 - 255 * c3)*y2 + (9 - 66 * c - 15 * c2)*y4)) / 512.;
	m_ip6[4][0][3] = (-15 * (7 * c2*(11 - 15 * c2) - (13 + 32 * c2 - 85 * c4)*y2 + 3 * (3 - 7 * c2)*y4)) / (128.*s);
	m_ip6[4][0][4] = (15 * (7 * c*(11 - 10 * c - 45 * c2) + (13 + 33 * c + 35 * c2 + 255 * c3)*y2 - (9 + 66 * c - 15 * c2)*y4)) / 512.;
	m_ip6[4][0][5] = (75 * s*(7 * c*(2 + 5 * c - 9 * c2) - (13 - 19 * c + 37 * c2 - 51 * c3)*y2 + 3 * (1 - 2 * c)*(3 - c)*y4)) / (256.*cn);
	m_ip6[4][0][6] = (-75 * cp*(7 * (1 - 3 * c)*c + (13 - 10 * c + 17 * c2)*y2 - 3 * (3 - c)*y4)) / 512.;
	m_ip6[4][1][0] = (-525 * cn2*(c - y2)) / 512.;
	m_ip6[4][1][1] = (525 * cn*s*(c*(2 + 3 * c) - (3 + 2 * c)*y2)) / (256.*cp);
	m_ip6[4][1][2] = (-525 * cn*(c*(1 + 3 * c) - (3 + c)*y2)) / 512.;
	m_ip6[4][1][3] = (525 * s*(c2 - y2)) / 128.;
	m_ip6[4][1][4] = (-525 * cp*((1 - 3 * c)*c + (3 - c)*y2)) / 512.;
	m_ip6[4][1][5] = (-525 * cp*s*((2 - 3 * c)*c + (3 - 2 * c)*y2)) / (256.*cn);
	m_ip6[4][1][6] = (-525 * cp2*(c + y2)) / 512.;

	m_ip6[5][0][0] = (75 * cn*(1 + 3 * c)) / 512.;
	m_ip6[5][0][1] = (75 * (2 - 5 * c - 9 * c2)*s) / (256.*cp);
	m_ip6[5][0][2] = (-15 * (11 + 10 * c - 45 * c2)) / 512.;
	m_ip6[5][0][3] = (15 * c*(11 - 15 * c2)) / (128.*s);
	m_ip6[5][0][4] = (-15 * (11 - 10 * c - 45 * c2)) / 512.;
	m_ip6[5][0][5] = (-75 * (2 + 5 * c - 9 * c2)*s) / (256.*cn);
	m_ip6[5][0][6] = (75 * (1 - 3 * c)*cp) / 512.;
	m_ip6[5][1][0] = (525 * cn2) / 512.;
	m_ip6[5][1][1] = (-525 * cn*(2 + 3 * c)*s) / (256.*cp);
	m_ip6[5][1][2] = (525 * cn*(1 + 3 * c)) / 512.;
	m_ip6[5][1][3] = (-525 * c*s) / 128.;
	m_ip6[5][1][4] = (525 * (1 - 3 * c)*cp) / 512.;
	m_ip6[5][1][5] = (525 * (2 - 3 * c)*cp*s) / (256.*cn);
	m_ip6[5][1][6] = (525 * cp2) / 512.;
}

template<typename T>
void HEOSATDynamics<T>::eccincpol7(const T& e, const T& y, const T& c, const T& s)
{
	using std::pow;   using DACE::pow;
	
	T e2, e4, s2, s3, s4, c2, c3, c4, cn, cn2, cn3, cn4, cp, cp2, cp3, cp4, y2, y4;

	e2 = e*e;
	e4 = e2*e2;
	s2 = s*s;
	s3 = s2*s;
	s4 = s3*s;
	c2 = c*c;
	c3 = c2*c;
	c4 = c3*c;

	y2 = y*y;
	y4 = y2*y2;
	cn = 1 - c;
	cn2 = cn*cn;
	cn3 = cn2*cn;
	cn4 = cn3*cn;
	cp = 1 + c;
	cp2 = cp*cp;
	cp3 = cp2*cp;
	cp4 = cp3*cp;


	m_ep7[1][0] = 8 + 5 * e2*(8 + 3 * e2);
	m_ep7[1][1] = e2*(2 + e2);
	m_ep7[1][2] = e4;
	m_ep7[2][0] = m_ep7[1][0];
	m_ep7[2][1] = m_ep7[1][1];
	m_ep7[2][2] = m_ep7[1][2];
	m_ep7[3][0] = 12 * (12 + 5 * e2*(5 + e2));
	m_ep7[3][1] = 2 * (1 + e2*(4 + e2));
	m_ep7[3][2] = e2*(1 + e2);
	m_ep7[4][0] = 1 / y;
	m_ep7[4][1] = m_ep7[4][0];
	m_ep7[4][2] = m_ep7[4][1] * e2;
	m_ep7[5][0] = m_ep7[1][0] / y;
	m_ep7[5][1] = m_ep7[1][1] / y;
	m_ep7[5][2] = m_ep7[1][2] / y;

	m_ip7[1][0][0] = (-105 * s4) / 8192.;
	m_ip7[1][0][1] = (-105 * c*s3) / 1024.;
	m_ip7[1][0][2] = (15 * (1 - 7 * c2)*s2) / 2048.;
	m_ip7[1][0][3] = (15 * c*(3 - 7 * c2)*s) / 1024.;
	m_ip7[1][0][4] = (-3 * (3 - 30 * c2 + 35 * c4)) / 4096.;
	m_ip7[1][0][5] = (15 * c*(3 - 7 * c2)*s) / 1024.;
	m_ip7[1][0][6] = (15 * (1 - 7 * c2)*s2) / 2048.;
	m_ip7[1][0][7] = (-105 * c*s3) / 1024.;
	m_ip7[1][0][8] = (-105 * s4) / 8192.;
	m_ip7[1][1][0] = (-735 * cn2*s2) / 1024.;
	m_ip7[1][1][1] = (-735 * cn2*(1 + 2 * c)*s) / 256.;
	m_ip7[1][1][2] = (-105 * cn2*(1 + 7 * c + 7 * c2)) / 256.;
	m_ip7[1][1][3] = (105 * cn*(1 - 7 * c - 14 * c2)*s) / 256.;
	m_ip7[1][1][4] = (105 * (1 - 7 * c2)*s2) / 512.;
	m_ip7[1][1][5] = (-105 * cp*(1 + 7 * c - 14 * c2)*s) / 256.;
	m_ip7[1][1][6] = (-105 * cp2*(1 - 7 * c + 7 * c2)) / 256.;
	m_ip7[1][1][7] = (735 * (1 - 2 * c)*cp2*s) / 256.;
	m_ip7[1][1][8] = (-735 * cp2*s2) / 1024.;
	m_ip7[1][2][0] = (-2205 * cn4) / 2048.;
	m_ip7[1][2][1] = (-2205 * cn3*s) / 256.;
	m_ip7[1][2][2] = (-2205 * cn2*s2) / 512.;
	m_ip7[1][2][3] = (-2205 * cn*s3) / 256.;
	m_ip7[1][2][4] = (-2205 * s4) / 1024.;
	m_ip7[1][2][5] = (2205 * cp*s3) / 256.;
	m_ip7[1][2][6] = (-2205 * cp2*s2) / 512.;
	m_ip7[1][2][7] = (2205 * cp3*s) / 256.;
	m_ip7[1][2][8] = (-2205 * cp4) / 2048.;
	for (int i = 0; i <= 2; ++i)
	{
		for (int j = -4; j <= 4; ++j)
		{
			m_ip7[2][i][j + 4] = (j / pow(2.0, i))*m_ip7[1][i][j + 4];
			m_ip7[3][i][j + 4] = m_ip7[1][i][j + 4];
		}
	}

	m_ip7[4][0][0] = (-105 * s2*(63 * c2 - 35 * (1 + c2)*y2 + 15 * y4)) / 2048.;
	m_ip7[4][0][1] = (105 * c*s*(63 - 252 * c2 + 70 * (1 + 2 * c2)*y2 - 45 * y4)) / 1024.;
	m_ip7[4][0][2] = (15 * (63 * c2*(4 - 7 * c2) + 35 * (-1 + 7 * c4)*y2 + (15 - 60 * c2)*y4)) / 512.;
	m_ip7[4][0][3] = (-15 * c*(63 * (3 - 27 * c2 + 28 * c4) + 70 * (3 + 7 * c2 - 14 * c4)*y2 - (135 - 195 * c2)*y4)) / (1024.*s);
	m_ip7[4][0][4] = (-15 * (63 * c2*(3 - 7 * c2) + 7 * (-3 + 35 * c4)*y2 + (9 - 45 * c2)*y4)) / 1024.;
	m_ip7[4][0][5] = (-15 * c*(63 * (3 - 27 * c2 + 28 * c4) + 70 * (3 + 7 * c2 - 14 * c4)*y2 - (135 - 195 * c2)*y4)) / (1024.*s);
	m_ip7[4][0][6] = (15 * (63 * c2*(4 - 7 * c2) + 35 * (-1 + 7 * c4)*y2 + (15 - 60 * c2)*y4)) / 512.;
	m_ip7[4][0][7] = (105 * c*s*(63 - 252 * c2 + 70 * (1 + 2 * c2)*y2 - 45 * y4)) / 1024.;
	m_ip7[4][0][8] = (-105 * s2*(63 * c2 - 35 * (1 + c2)*y2 + 15 * y4)) / 2048.;
	m_ip7[4][1][0] = (-735 * cn2*(3 * c*(1 + 2 * c) - 4 * (1 + c + c2)*y2 + (2 + c)*y4)) / 1024.;
	m_ip7[4][1][1] = (-735 * cn*s*(3 * c2*(7 + 8 * c) - 4 * (2 + 4 * c + 5 * c2 + 4 * c3)*y2 + (4 + 8 * c + 3 * c2)*y4)) / (512.*cp);
	m_ip7[4][1][2] = (105 * cn*(3 * c*(5 - 7 * c - 28 * c2) + 4 * (2 + 7 * c + 7 * c2 + 14 * c3)*y2 - (4 + 19 * c + 7 * c2)*y4)) / 512.;
	m_ip7[4][1][3] = (105 * s*(3 * c*(8 + 23 * c - 35 * c2 - 56 * c3) - 4 * (2 - 6 * c - 7 * c2 - 21 * c3 - 28 * c4)*y2 + (4 - 20 * c - 37 * c2 - 7 * c3)*y4)) / (512.*cp);
	m_ip7[4][1][4] = (105 * (3 * c2*(4 - 7 * c2) - (2 - 14 * c4)*y2 + (1 - 4 * c2)*y4)) / 256.;
	m_ip7[4][1][5] = (105 * s*(3 * c*(8 - 23 * c - 35 * c2 + 56 * c3) + 4 * (2 + 6 * c - 7 * c2 + 21 * c3 - 28 * c4)*y2 - (4 + 20 * c - 37 * c2 + 7 * c3)*y4)) / (512.*cn);
	m_ip7[4][1][6] = (-105 * cp*(3 * c*(5 + 7 * c - 28 * c2) + 4 * (-2 + 7 * c - 7 * c2 + 14 * c3)*y2 + (4 - 19 * c + 7 * c2)*y4)) / 512.;
	m_ip7[4][1][7] = (735 * cp*s*(3 * (7 - 8 * c)*c2 + 4 * (-2 + 4 * c - 5 * c2 + 4 * c3)*y2 + (4 - 8 * c + 3 * c2)*y4)) / (512.*cn);
	m_ip7[4][1][8] = (735 * cp2*(3 * (1 - 2 * c)*c + 4 * (1 - c + c2)*y2 - (2 - c)*y4)) / 1024.;
	m_ip7[4][2][0] = (-2205 * cn3*(c - y2)) / 2048.;
	m_ip7[4][2][1] = (-2205 * cn2*s*(c*(3 + 4 * c) - (4 + 3 * c)*y2)) / (1024.*cp);
	m_ip7[4][2][2] = (-2205 * cn2*(c*(1 + 2 * c) - (2 + c)*y2)) / 1024.;
	m_ip7[4][2][3] = (-2205 * cn*s*(c*(1 + 4 * c) - (4 + c)*y2)) / 1024.;
	m_ip7[4][2][4] = (-2205 * s2*(c2 - y2)) / 1024.;
	m_ip7[4][2][5] = (-2205 * cp*s*((1 - 4 * c)*c + (4 - c)*y2)) / 1024.;
	m_ip7[4][2][6] = (2205 * cp2*((1 - 2 * c)*c + (2 - c)*y2)) / 1024.;
	m_ip7[4][2][7] = (-2205 * cp2*s*((3 - 4 * c)*c + (4 - 3 * c)*y2)) / (1024.*cn);
	m_ip7[4][2][8] = (2205 * cp3*(c + y2)) / 2048.;

	m_ip7[5][0][0] = (105 * c*s2) / 2048.;
	m_ip7[5][0][1] = (-105 * (1 - 4 * c2)*s) / 1024.;
	m_ip7[5][0][2] = (-15 * c*(4 - 7 * c2)) / 512.;
	m_ip7[5][0][3] = (15 * (3 - 27 * c2 + 28 * c4)) / (1024.*s);
	m_ip7[5][0][4] = (15 * c*(3 - 7 * c2)) / 1024.;
	m_ip7[5][0][5] = (15 * (3 - 27 * c2 + 28 * c4)) / (1024.*s);
	m_ip7[5][0][6] = (-15 * c*(4 - 7 * c2)) / 512.;
	m_ip7[5][0][7] = (-105 * (1 - 4 * c2)*s) / 1024.;
	m_ip7[5][0][8] = (105 * c*s2) / 2048.;
	m_ip7[5][1][0] = (735 * cn2*(1 + 2 * c)) / 1024.;
	m_ip7[5][1][1] = (735 * cn*c*(7 + 8 * c)*s) / (512.*cp);
	m_ip7[5][1][2] = (-105 * cn*(5 - 7 * c - 28 * c2)) / 512.;
	m_ip7[5][1][3] = (-105 * (8 + 23 * c - 35 * c2 - 56 * c3)*s) / (512.*(1 + c));
	m_ip7[5][1][4] = (-105 * c*(4 - 7 * c2)) / 256.;
	m_ip7[5][1][5] = (-105 * (8 - 23 * c - 35 * c2 + 56 * c3)*s) / (512.*cn);
	m_ip7[5][1][6] = (105 * cp*(5 + 7 * c - 28 * c2)) / 512.;
	m_ip7[5][1][7] = (-735 * c*(7 - c - 8 * c2)*s) / (512.*cn);
	m_ip7[5][1][8] = (-735 * (1 - 2 * c)*cp2) / 1024.;
	m_ip7[5][2][0] = (2205 * cn3) / 2048.;
	m_ip7[5][2][1] = (2205 * cn2*(3 + 4 * c)*s) / (1024.*cp);
	m_ip7[5][2][2] = (2205 * cn2*(1 + 2 * c)) / 1024.;
	m_ip7[5][2][3] = (2205 * (1 + 3 * c - 4 * c2)*s) / 1024.;
	m_ip7[5][2][4] = (2205 * c*s2) / 1024.;
	m_ip7[5][2][5] = (2205 * (1 - 3 * c - 4 * c2)*s) / 1024.;
	m_ip7[5][2][6] = (-2205 * (1 - 2 * c)*cp2) / 1024.;
	m_ip7[5][2][7] = (2205 * (3 - 4 * c)*cp2*s) / (1024.*cn);
	m_ip7[5][2][8] = (-2205 * cp3) / 2048.;

}

template<typename T>
void HEOSATDynamics<T>::eccincpol8(const T& e, const T& y, const T& c, const T& s)
{
	T e2, e3, e4, e5, s2, s3, s4, s5, c2, c3, c4, c5, c6, cn, cn2, cn3, cn4, cn5, cp, cp2, cp3, cp4, cp5, y2, y4, y6;

	e2 = e*e;
	e3 = e2*e;
	e4 = e3*e;
	e5 = e4*e;
	s2 = s*s;
	s3 = s2*s;
	s4 = s3*s;
	s5 = s4*s;
	c2 = c*c;
	c3 = c2*c;
	c4 = c3*c;
	c5 = c4*c;
	c6 = c5*c;

	y2 = y*y;
	y4 = y2*y2;
	y6 = y4*y2;
	cn = 1 - c;
	cn2 = cn*cn;
	cn3 = cn2*cn;
	cn4 = cn3*cn;
	cn5 = cn4*cn;
	cp = 1 + c;
	cp2 = cp*cp;
	cp3 = cp2*cp;
	cp4 = cp3*cp;
	cp5 = cp4*cp;

	m_ep8[1][0] = e*(8 + 20 * e2 + 5 * e4);
	m_ep8[1][1] = e3*(8 + 3 * e2);
	m_ep8[1][2] = e5;
	m_ep8[2][0] = m_ep8[1][0];
	m_ep8[2][1] = m_ep8[1][1];
	m_ep8[2][2] = m_ep8[1][2];
	m_ep8[3][0] = 8 / e + 132 * e + 165 * e3 + 25 * e5;
	m_ep8[3][1] = 8 * e + (71 * e3) / 3. + 5 * e5;
	m_ep8[3][2] = e3*(1 + e2);
	m_ep8[4][0] = 1 / (e*y);
	m_ep8[4][1] = e / y;
	m_ep8[4][2] = m_ep8[4][1] * e2;
	m_ep8[5][0] = m_ep8[2][0] / y;
	m_ep8[5][1] = m_ep8[2][1] / y;
	m_ep8[5][2] = m_ep8[2][2] / y;

	m_ip8[1][0][0] = (2205 * cn*s4) / 32768.;
	m_ip8[1][0][1] = (-2205 * cn*(1 + 5 * c)*s3) / 16384.;
	m_ip8[1][0][2] = (-735 * cn*(1 - 6 * c - 15 * c2)*s2) / 32768.;
	m_ip8[1][0][3] = (735 * cn*(1 + 3 * c - 9 * c2 - 15 * c3)*s) / 4096.;
	m_ip8[1][0][4] = (105 * cn*(1 - 28 * c - 42 * c2 + 84 * c3 + 105 * c4)) / 16384.;
	m_ip8[1][0][5] = (-105 * (1 - 14 * c2 + 21 * c4)*s) / 8192.;
	m_ip8[1][0][6] = (-105 * cp*(1 + 28 * c - 42 * c2 - 84 * c3 + 105 * c4)) / 16384.;
	m_ip8[1][0][7] = (735 * (1 - 2 * c - 12 * c2 + 6 * c3 + 15 * c4)*s) / 4096.;
	m_ip8[1][0][8] = (735 * cp*(1 + 6 * c - 15 * c2)*s2) / 32768.;
	m_ip8[1][0][9] = (-2205 * (1 - 5 * c)*cp*s3) / 16384.;
	m_ip8[1][0][10] = (-2205 * cp*s4) / 32768.;

	m_ip8[1][1][0] = (19845 * cn3*s2) / 65536.;
	m_ip8[1][1][1] = (-19845 * cn3*(3 + 5 * c)*s) / 32768.;
	m_ip8[1][1][2] = (2205 * cn3*(1 + 3 * c)*(13 + 15 * c)) / 65536.;
	m_ip8[1][1][3] = (-6615 * cn2*(1 + 12 * c + 15 * c2)*s) / 8192.;
	m_ip8[1][1][4] = (-6615 * cn*(1 - 6 * c - 15 * c2)*s2) / 32768.;
	m_ip8[1][1][5] = (2205 * (1 - 9 * c2)*s3) / 16384.;
	m_ip8[1][1][6] = (6615 * cp*(1 + 6 * c - 15 * c2)*s2) / 32768.;
	m_ip8[1][1][7] = (-6615 * cp2*(1 - 12 * c + 15 * c2)*s) / 8192.;
	m_ip8[1][1][8] = (-2205 * (13 - 15 * c)*(1 - 3 * c)*cp3) / 65536.;
	m_ip8[1][1][9] = (-19845 * (3 - 5 * c)*cp3*s) / 32768.;
	m_ip8[1][1][10] = (-19845 * cp3*s2) / 65536.;

	m_ip8[1][2][0] = (72765 * cn5) / 65536.;
	m_ip8[1][2][1] = (-363825 * cn4*s) / 32768.;
	m_ip8[1][2][2] = (363825 * cn3*s2) / 65536.;
	m_ip8[1][2][3] = (-363825 * cn2*s3) / 8192.;
	m_ip8[1][2][4] = (363825 * cn*s4) / 32768.;
	m_ip8[1][2][5] = (-72765 * s5) / 16384.;
	m_ip8[1][2][6] = (-363825 * cp*s4) / 32768.;
	m_ip8[1][2][7] = (-363825 * cp2*s3) / 8192.;
	m_ip8[1][2][8] = (-363825 * cp3*s2) / 65536.;
	m_ip8[1][2][9] = (-363825 * cp4*s) / 32768.;
	m_ip8[1][2][10] = (-72765 * cp5) / 65536.;

	for (int i = 0; i <= 2; ++i)
	{
		for (int j = -5; j <= 5; ++j)
		{
			m_ip8[2][i][j + 5] = (j / (2.0*i + 1))*m_ip8[1][i][j + 5];
			m_ip8[3][i][j + 5] = -m_ip8[1][i][j + 5];
		}
	}
	m_ip8[4][0][0] = (-2205 * cn*s2*(33 * c*(1 + 5 * c) - 3 * (31 + 21 * c + 74 * c2)*y2 + 5 * (22 + 7 * c + 13 * c2)*y4 - 5 * (5 + c)*y6)) / 32768.;
	m_ip8[4][0][1] = (2205 * cn*s*(33 * c*(-4 + 9 * c + 25 * c2) - 3 * (31 + 71 * c + 158 * c2 + 370 * c3)*y2 + 5 * (22 + 82 * c + 41 * c2 + 65 * c3)*y4 - 5 * (5 + 21 * c + 4 * c2)*y6)) / 16384.;
	m_ip8[4][0][2] = (-735 * cn*(33 * c*(-7 - 27 * c + 39 * c2 + 75 * c3) - 3 * (-31 + 39 * c - 71 * c2 + 633 * c3 + 1110 * c4)*y2 + 5 * (-22 + 83 * c + 163 * c2 + 141 * c3 + 195 * c4)*y4 - 5 * (-5 + 23 * c + 53 * c2 + 9 * c3)*y6)) / 32768.;
	m_ip8[4][0][3] = (735 * s*(33 * c*(2 - 23 * c - 45 * c2 + 51 * c3 + 75 * c4) + 3 * (31 + 51 * c + 173 * c2 + 387 * c3 - 792 * c4 - 1110 * c5)*y2 + 5 * (-22 - 52 * c + 59 * c2 + 81 * c3 + 159 * c4 + 195 * c5)*y4 - 5 * (-5 - 13 * c + 27 * c2 + 45 * c3 + 6 * c4)*y6)) / (4096.*cp);
	m_ip8[4][0][4] = (-105 * (33 * c*(29 + 28 * c - 378 * c2 - 84 * c3 + 525 * c4) + (-93 + 870 * c - 462 * c2 + 12096 * c3 + 3339 * c4 - 23310 * c5)*y2 + 5 * (22 - 435 * c - 112 * c2 + 126 * c3 - 126 * c4 + 1365 * c5)*y4 - 5 * (5 - 116 * c - 42 * c2 + 252 * c3 + 21 * c4)*y6)) / 16384.;
	m_ip8[4][0][5] = (105 * (33 * c2*(29 - 126 * c2 + 105 * c4) + (-93 - 432 * c2 + 4683 * c4 - 4662 * c6)*y2 + 5 * (22 - 127 * c2 - 112 * c4 + 273 * c6)*y4 - 5 * (5 - 46 * c2 + 49 * c4)*y6)) / (8192.*s);
	m_ip8[4][0][6] = (-105 * (33 * c*(29 - 28 * c - 378 * c2 + 84 * c3 + 525 * c4) + (93 + 870 * c + 462 * c2 + 12096 * c3 - 3339 * c4 - 23310 * c5)*y2 + 5 * (-22 - 435 * c + 112 * c2 + 126 * c3 + 126 * c4 + 1365 * c5)*y4 + 5 * (5 + 116 * c - 42 * c2 - 252 * c3 + 21 * c4)*y6)) / 16384.;
	m_ip8[4][0][7] = (-735 * s*(33 * c*(2 + 23 * c - 45 * c2 - 51 * c3 + 75 * c4) + 3 * (-31 + 51 * c - 173 * c2 + 387 * c3 + 792 * c4 - 1110 * c5)*y2 + 5 * (22 - 52 * c - 59 * c2 + 81 * c3 - 159 * c4 + 195 * c5)*y4 + 5 * (-5 + 13 * c + 27 * c2 - 45 * c3 + 6 * c4)*y6)) / (4096.*cn);
	m_ip8[4][0][8] = (735 * cp*(33 * c*(7 - 27 * c - 39 * c2 + 75 * c3) + 3 * (31 + 39 * c + 71 * c2 + 633 * c3 - 1110 * c4)*y2 + 5 * (-22 - 83 * c + 163 * c2 - 141 * c3 + 195 * c4)*y4 + 5 * (5 + 23 * c - 53 * c2 + 9 * c3)*y6)) / 32768.;
	m_ip8[4][0][9] = (-2205 * cp*s*(33 * c*(-4 - 9 * c + 25 * c2) + (93 - 213 * c + 474 * c2 - 1110 * c3)*y2 + 5 * (-22 + 82 * c - 41 * c2 + 65 * c3)*y4 + 5 * (5 - 21 * c + 4 * c2)*y6)) / 16384.;
	m_ip8[4][0][10] = (2205 * cp*s2*(33 * c*(-1 + 5 * c) + (-93 + 63 * c - 222 * c2)*y2 + 5 * (22 - 7 * c + 13 * c2)*y4 + 5 * (-5 + c)*y6)) / 32768.;

	m_ip8[4][1][0] = (6615 * cn3*(-11 * c*(3 + 5 * c) + (39 + 42 * c + 31 * c2)*y2 - 3 * (5 + 3 * c)*y4)) / 65536.;
	m_ip8[4][1][1] = (-6615 * cn2*s*(-11 * c*(4 + 27 * c + 25 * c2) + (117 + 251 * c + 261 * c2 + 155 * c3)*y2 - 3 * (15 + 29 * c + 12 * c2)*y4)) / (32768.*cp);
	m_ip8[4][1][2] = (2205 * cn2*(-11 * c*(-5 + 42 * c + 75 * c2) + (169 + 463 * c + 471 * c2 + 465 * c3)*y2 + (-65 - 190 * c - 81 * c2)*y4)) / 65536.;
	m_ip8[4][1][3] = (-2205 * cn*s*(11 * c*(10 + 3 * c - 78 * c2 - 75 * c3) + (39 + 328 * c + 504 * c2 + 624 * c3 + 465 * c4)*y2 - 3 * (5 + 50 * c + 67 * c2 + 18 * c3)*y4)) / (8192.*cp);
	m_ip8[4][1][4] = (2205 * cn*(11 * c*(7 + 27 * c - 39 * c2 - 75 * c3) + (-39 + 136 * c + 246 * c2 + 312 * c3 + 465 * c4)*y2 - 3 * (-5 + 23 * c + 53 * c2 + 9 * c3)*y4)) / 32768.;
	m_ip8[4][1][5] = (-2205 * s*(77 * c2 - 165 * c4 + (-13 + 32 * c2 + 93 * c4)*y2 + (5 - 29 * c2)*y4)) / 16384.;
	m_ip8[4][1][6] = (-2205 * cp*(11 * c*(-7 + 27 * c + 39 * c2 - 75 * c3) + (-39 - 136 * c + 246 * c2 - 312 * c3 + 465 * c4)*y2 + 3 * (5 + 23 * c - 53 * c2 + 9 * c3)*y4)) / 32768.;
	m_ip8[4][1][7] = (-2205 * cp*s*(11 * c*(-10 + 3 * c + 78 * c2 - 75 * c3) + (39 - 328 * c + 504 * c2 - 624 * c3 + 465 * c4)*y2 + 3 * (-5 + 50 * c - 67 * c2 + 18 * c3)*y4)) / (8192.*cn);
	m_ip8[4][1][8] = (2205 * cp2*(11 * c*(5 + 42 * c - 75 * c2) + (-169 + 463 * c - 471 * c2 + 465 * c3)*y2 + (65 - 190 * c + 81 * c2)*y4)) / 65536.;
	m_ip8[4][1][9] = (6615 * cp2*s*(-11 * c*(4 - 27 * c + 25 * c2) + (-117 + 251 * c - 261 * c2 + 155 * c3)*y2 + (45 - 87 * c + 36 * c2)*y4)) / (32768.*cn);
	m_ip8[4][1][10] = (-6615 * cp3*(11 * (3 - 5 * c)*c + (39 - 42 * c + 31 * c2)*y2 + (-15 + 9 * c)*y4)) / 65536.;

	m_ip8[4][2][0] = (-72765 * cn4*(c - y2)) / 65536.;
	m_ip8[4][2][1] = (72765 * cn3*s*(c*(4 + 5 * c) + (-5 - 4 * c)*y2)) / (32768.*cp);
	m_ip8[4][2][2] = (-72765 * cn3*(c*(3 + 5 * c) + (-5 - 3 * c)*y2)) / 65536.;
	m_ip8[4][2][3] = (72765 * cn2*s*(c*(2 + 5 * c) + (-5 - 2 * c)*y2)) / 8192.;
	m_ip8[4][2][4] = (-72765 * cn*s2*(c*(1 + 5 * c) + (-5 - c)*y2)) / 32768.;
	m_ip8[4][2][5] = (72765 * s3*(c2 - y2)) / 16384.;
	m_ip8[4][2][6] = (72765 * cp*s2*(c*(-1 + 5 * c) + (-5 + c)*y2)) / 32768.;
	m_ip8[4][2][7] = (72765 * cp2*s*(c*(-2 + 5 * c) + (-5 + 2 * c)*y2)) / 8192.;
	m_ip8[4][2][8] = (72765 * cp3*(c*(-3 + 5 * c) + (-5 + 3 * c)*y2)) / 65536.;
	m_ip8[4][2][9] = (72765 * cp3*s*(c*(-4 + 5 * c) + (-5 + 4 * c)*y2)) / (32768.*cn);
	m_ip8[4][2][10] = (-72765 * cp4*(c + y2)) / 65536.;

	m_ip8[5][0][0] = (2205 * cn*(1 + 5 * c)*s2) / 32768.;
	m_ip8[5][0][1] = (2205 * cn*(4 - 9 * c - 25 * c2)*s) / 16384.;
	m_ip8[5][0][2] = (-735 * cn*(7 + 27 * c - 39 * c2 - 75 * c3)) / 32768.;
	m_ip8[5][0][3] = (-735 * (2 - 23 * c - 45 * c2 + 51 * c3 + 75 * c4)*s) / (4096.*cp);
	m_ip8[5][0][4] = (105 * (29 + 28 * c - 378 * c2 - 84 * c3 + 525 * c4)) / 16384.;
	m_ip8[5][0][5] = (-105 * c*(29 - 126 * c2 + 105 * c4)) / (8192.*s);
	m_ip8[5][0][6] = (105 * (29 - 28 * c - 378 * c2 + 84 * c3 + 525 * c4)) / 16384.;
	m_ip8[5][0][7] = (735 * (2 + 23 * c - 45 * c2 - 51 * c3 + 75 * c4)*s) / (4096.*cn);
	m_ip8[5][0][8] = (-735 * cp*(7 - 27 * c - 39 * c2 + 75 * c3)) / 32768.;
	m_ip8[5][0][9] = (-2205 * cp*(4 + 9 * c - 25 * c2)*s) / 16384.;
	m_ip8[5][0][10] = (2205 * (1 - 5 * c)*cp*s2) / 32768.;

	m_ip8[5][1][0] = (6615 * cn3*(3 + 5 * c)) / 65536.;
	m_ip8[5][1][1] = (-6615 * cn2*(4 + 27 * c + 25 * c2)*s) / (32768.*cp);
	m_ip8[5][1][2] = (-2205 * cn2*(5 - 42 * c - 75 * c2)) / 65536.;
	m_ip8[5][1][3] = (2205 * cn*(10 + 3 * c - 78 * c2 - 75 * c3)*s) / (8192.*cp);
	m_ip8[5][1][4] = (-2205 * cn*(7 + 27 * c - 39 * c2 - 75 * c3)) / 32768.;
	m_ip8[5][1][5] = (2205 * c*(7 - 15 * c2)*s) / 16384.;
	m_ip8[5][1][6] = (-2205 * cp*(7 - 27 * c - 39 * c2 + 75 * c3)) / 32768.;
	m_ip8[5][1][7] = (-2205 * cp*(10 - 3 * c - 78 * c2 + 75 * c3)*s) / (8192.*cn);
	m_ip8[5][1][8] = (-2205 * cp2*(5 + 42 * c - 75 * c2)) / 65536.;
	m_ip8[5][1][9] = (6615 * cp2*(4 - 27 * c + 25 * c2)*s) / (32768.*cn);
	m_ip8[5][1][10] = (6615 * (3 - 5 * c)*cp3) / 65536.;

	m_ip8[5][2][0] = (72765 * cn4) / 65536.;
	m_ip8[5][2][1] = (-72765 * cn3*(4 + 5 * c)*s) / (32768.*cp);
	m_ip8[5][2][2] = (72765 * cn3*(3 + 5 * c)) / 65536.;
	m_ip8[5][2][3] = (-72765 * cn2*(2 + 5 * c)*s) / 8192.;
	m_ip8[5][2][4] = (72765 * cn*(1 + 5 * c)*s2) / 32768.;
	m_ip8[5][2][5] = (-72765 * c*s3) / 16384.;
	m_ip8[5][2][6] = (72765 * (1 - 5 * c)*cp*s2) / 32768.;
	m_ip8[5][2][7] = (72765 * (2 - 5 * c)*cp2*s) / 8192.;
	m_ip8[5][2][8] = (72765 * (3 - 5 * c)*cp3) / 65536.;
	m_ip8[5][2][9] = (72765 * (4 - 5 * c)*cp3*s) / (32768.*cn);
	m_ip8[5][2][10] = (72765 * cp4) / 65536.;
}

template<typename T>
void HEOSATDynamics<T>::eccincpol9(const T& e, const T& y, const T& c, const T& s)
{
	using std::pow;   using DACE::pow;
	
	T e2, e3, e4, e6, s2, s3, s4, s5, s6, c2, c3, c4, c5, c6, cn, cn2, cn3, cn4, cn5, cn6, cp, cp2, cp3, cp4, cp5, cp6, y2, y4, y6;

	e2 = e*e;
	e3 = e2*e;
	e4 = e3*e;
	e6 = e3*e3;
	s2 = s*s;
	s3 = s2*s;
	s4 = s3*s;
	s5 = s4*s;
	s6 = s5*s;
	c2 = c*c;
	c3 = c2*c;
	c4 = c3*c;
	c5 = c4*c;
	c6 = c5*c;

	y2 = y*y;
	y4 = y2*y2;
	y6 = y4*y2;
	cn = 1 - c;
	cn2 = cn*cn;
	cn3 = cn2*cn;
	cn4 = cn3*cn;
	cn5 = cn4*cn;
	cn6 = cn5*cn;
	cp = 1 + c;
	cp2 = cp*cp;
	cp3 = cp2*cp;
	cp4 = cp3*cp;
	cp5 = cp4*cp;
	cp6 = cp5*cp;

	m_ep9[1][0] = -16 - 168 * e2 - 210 * e4 - 35 * e6;
	m_ep9[1][1] = -(e2*(48 + 80 * e2 + 15 * e4));
	m_ep9[1][2] = -(e4*(10 + 3 * e2));
	m_ep9[1][3] = -e6;
	m_ep9[2][0] = m_ep9[1][0];
	m_ep9[2][1] = m_ep9[1][1];
	m_ep9[2][2] = m_ep9[1][2];
	m_ep9[2][3] = (4 * m_ep9[1][3]) / 3.;
	m_ep9[3][0] = -6 * (88 + 420 * e2 + 315 * e4 + 35 * e6);
	m_ep9[3][1] = -48 - 400 * e2 - 365 * e4 - 45 * e6;
	m_ep9[3][2] = -(e2*(20 + 49 * e2 + 9 * e4)) / 2.;
	m_ep9[3][3] = -(e4*(1 + e2));
	m_ep9[4][0] = 1.0 / y;
	m_ep9[4][1] = 1.0 / y;
	m_ep9[4][2] = -(e2 / y);
	m_ep9[4][3] = e4 / y;
	m_ep9[5][0] = m_ep9[1][0] / y;
	m_ep9[5][1] = m_ep9[1][1] / y;
	m_ep9[5][2] = m_ep9[1][2] / y;
	m_ep9[5][3] = m_ep9[1][3] / y;

	m_ip9[1][0][0] = (1155 * s6) / 262144.;
	m_ip9[1][0][1] = (3465 * c*s5) / 65536.;
	m_ip9[1][0][2] = (-315 * (1 - 11 * c2)*s4) / 131072.;
	m_ip9[1][0][3] = (525 * c*(3 - 11 * c2)*s3) / 65536.;
	m_ip9[1][0][4] = (525 * (1 - 18 * c2 + 33 * c4)*s2) / 262144.;
	m_ip9[1][0][5] = (105 * c*(5 - 30 * c2 + 33 * c4)*s) / 32768.;
	m_ip9[1][0][6] = (-5 * (5 - 105 * c2 + 315 * c4 - 231 * c6)) / 65536.;
	for (int i = 1; i <= 6; ++i)
	{
		m_ip9[1][0][i + 6] = m_ip9[1][0][-i + 6];
	}

	m_ip9[1][1][0] = (10395 * cn2*s4) / 262144.;
	m_ip9[1][1][1] = (10395 * cn2*(1 + 3 * c)*s3) / 65536.;
	m_ip9[1][1][2] = (945 * cn2*(1 + 22 * c + 33 * c2)*s2) / 131072.;
	m_ip9[1][1][3] = (945 * cn2*(3 - 5 * c - 55 * c2 - 55 * c3)*s) / 65536.;
	m_ip9[1][1][4] = (-315 * cn2*(17 + 108 * c - 90 * c2 - 660 * c3 - 495 * c4)) / 262144.;
	m_ip9[1][1][5] = (315 * cn*(1 - 18 * c - 36 * c2 + 66 * c3 + 99 * c4)*s) / 32768.;
	m_ip9[1][1][6] = (315 * (1 - 18 * c2 + 33 * c4)*s2) / 65536.;
	m_ip9[1][1][7] = (-315 * cp*(1 + 18 * c - 36 * c2 - 66 * c3 + 99 * c4)*s) / 32768.;
	m_ip9[1][1][8] = (-315 * cp2*(17 - 108 * c - 90 * c2 + 660 * c3 - 495 * c4)) / 262144.;
	m_ip9[1][1][9] = (-945 * cp2*(3 + 5 * c - 55 * c2 + 55 * c3)*s) / 65536.;
	m_ip9[1][1][10] = (945 * cp2*(1 - 22 * c + 33 * c2)*s2) / 131072.;
	m_ip9[1][1][11] = (-10395 * (1 - 3 * c)*cp2*s3) / 65536.;
	m_ip9[1][1][12] = (10395 * cp2*s4) / 262144.;

	m_ip9[1][2][0] = (22869 * cn4*s2) / 65536.;
	m_ip9[1][2][1] = (22869 * cn4*(2 + 3 * c)*s) / 16384.;
	m_ip9[1][2][2] = (2079 * cn4*(13 + 44 * c + 33 * c2)) / 32768.;
	m_ip9[1][2][3] = (-10395 * cn3*(2 + 11 * c + 11 * c2)*s) / 16384.;
	m_ip9[1][2][4] = (10395 * cn2*(1 + 22 * c + 33 * c2)*s2) / 65536.;
	m_ip9[1][2][5] = (-2079 * cn*(2 - 11 * c - 33 * c2)*s3) / 8192.;
	m_ip9[1][2][6] = (-2079 * (1 - 11 * c2)*s4) / 16384.;
	m_ip9[1][2][7] = (2079 * cp*(2 + 11 * c - 33 * c2)*s3) / 8192.;
	m_ip9[1][2][8] = (10395 * cp2*(1 - 22 * c + 33 * c2)*s2) / 65536.;
	m_ip9[1][2][9] = (10395 * cp3*(2 - 11 * c + 11 * c2)*s) / 16384.;
	m_ip9[1][2][10] = (2079 * cp4*(13 - 44 * c + 33 * c2)) / 32768.;
	m_ip9[1][2][11] = (-22869 * (2 - 3 * c)*cp4*s) / 16384.;
	m_ip9[1][2][12] = (22869 * cp4*s2) / 65536.;

	m_ip9[1][3][0] = (297297 * cn6) / 262144.;
	m_ip9[1][3][1] = (891891 * cn5*s) / 65536.;
	m_ip9[1][3][2] = (891891 * cn4*s2) / 131072.;
	m_ip9[1][3][3] = (-1486485 * cn3*s3) / 65536.;
	m_ip9[1][3][4] = (4459455 * cn2*s4) / 262144.;
	m_ip9[1][3][5] = (891891 * cn*s5) / 32768.;
	m_ip9[1][3][6] = (297297 * s6) / 65536.;
	m_ip9[1][3][7] = (-891891 * cp*s5) / 32768.;
	m_ip9[1][3][8] = (4459455 * cp2*s4) / 262144.;
	m_ip9[1][3][9] = (1486485 * cp3*s3) / 65536.;
	m_ip9[1][3][10] = (891891 * cp4*s2) / 131072.;
	m_ip9[1][3][11] = (-891891 * cp5*s) / 65536.;
	m_ip9[1][3][12] = (297297 * cp6) / 262144.;

	for (int i = 0; i <= 3; ++i)
	{
		for (int j = -6; j <= 6; ++j)
		{
			m_ip9[2][i][j + 6] = (j / pow(2., i))*m_ip9[1][i][j + 6];
			m_ip9[3][i][j + 6] = m_ip9[1][i][j + 6];
		}
	}

	m_ip9[4][0][0] = (-3465 * (8 * (-7 + 9 * c2) + 84 * (-1 + 3 * c2)*e2 + 105 * (1 + c2)*e4 + 35 * e6)*s4) / 131072.;
	m_ip9[4][0][1] = (3465 * c*(16 * (22 - 49 * c2 + 27 * c4) + 168 * (4 - 13 * c2 + 9 * c4)*e2 + 210 * (-2 - c2 + 3 * c4)*e4 + 175 * (-1 + c2)*e6)*s) / 65536.;
	m_ip9[4][0][2] = (315 * (8 * (-21 + 299 * c2 - 575 * c4 + 297 * c6) + 84 * (-3 + 65 * c2 - 161 * c4 + 99 * c6)*e2 + 105 * (3 - 13 * c2 - 23 * c4 + 33 * c6)*e4 + 35 * (3 - 26 * c2 + 23 * c4)*e6)) / 65536.;
	m_ip9[4][0][3] = (1575 * c*(16 * (22 - 113 * c2 + 99 * c4) + 168 * (4 - 29 * c2 + 33 * c4)*e2 + 210 * (-2 - c2 + 11 * c4)*e4 + 35 * (-5 + 13 * c2)*e6)*s) / 65536.;
	m_ip9[4][0][4] = (-525 * (8 * (-21 + 437 * c2 - 1275 * c4 + 891 * c6) + 84 * (-3 + 95 * c2 - 357 * c4 + 297 * c6)*e2 + 105 * (3 - 19 * c2 - 51 * c4 + 99 * c6)*e4 + 35 * (3 - 38 * c2 + 51 * c4)*e6)) / 131072.;
	m_ip9[4][0][5] = (-105 * c*(16 * (-110 + 835 * c2 - 1608 * c4 + 891 * c6) + 168 * (-20 + 205 * c2 - 474 * c4 + 297 * c6)*e2 + 210 * (10 - 5 * c2 - 96 * c4 + 99 * c6)*e4 + 35 * (25 - 110 * c2 + 93 * c4)*e6)) / (32768.*s);
	m_ip9[4][0][6] = (105 * (8 * (-5 + 115 * c2 - 375 * c4 + 297 * c6) + 12 * (-5 + 175 * c2 - 735 * c4 + 693 * c6)*e2 + 15 * (5 - 35 * c2 - 105 * c4 + 231 * c6)*e4 + 25 * (1 - 14 * c2 + 21 * c4)*e6)) / 32768.;
	m_ip9[4][0][7] = (-105 * c*(16 * (-110 + 835 * c2 - 1608 * c4 + 891 * c6) + 168 * (-20 + 205 * c2 - 474 * c4 + 297 * c6)*e2 + 210 * (10 - 5 * c2 - 96 * c4 + 99 * c6)*e4 + 35 * (25 - 110 * c2 + 93 * c4)*e6)) / (32768.*s);
	m_ip9[4][0][8] = (-525 * (8 * (-21 + 437 * c2 - 1275 * c4 + 891 * c6) + 84 * (-3 + 95 * c2 - 357 * c4 + 297 * c6)*e2 + 105 * (3 - 19 * c2 - 51 * c4 + 99 * c6)*e4 + 35 * (3 - 38 * c2 + 51 * c4)*e6)) / 131072.;
	m_ip9[4][0][9] = (1575 * c*(16 * (22 - 113 * c2 + 99 * c4) + 168 * (4 - 29 * c2 + 33 * c4)*e2 + 210 * (-2 - c2 + 11 * c4)*e4 + 35 * (-5 + 13 * c2)*e6)*s) / 65536.;
	m_ip9[4][0][10] = (315 * (8 * (-21 + 299 * c2 - 575 * c4 + 297 * c6) + 84 * (-3 + 65 * c2 - 161 * c4 + 99 * c6)*e2 + 105 * (3 - 13 * c2 - 23 * c4 + 33 * c6)*e4 + 35 * (3 - 26 * c2 + 23 * c4)*e6)) / 65536.;
	m_ip9[4][0][11] = (3465 * c*(16 * (22 - 49 * c2 + 27 * c4) + 168 * (4 - 13 * c2 + 9 * c4)*e2 + 210 * (-2 - c2 + 3 * c4)*e4 + 175 * (-1 + c2)*e6)*s) / 65536.;
	m_ip9[4][0][12] = (-3465 * (8 * (-7 + 9 * c2) + 84 * (-1 + 3 * c2)*e2 + 105 * (1 + c2)*e4 + 35 * e6)*s4) / 131072.;

	m_ip9[4][1][0] = (-10395 * cn2*(48 * (-1 + c2) + 16 * (-7 + 3 * c + 16 * c2)*e2 + 5 * (23 + 16 * c + 25 * c2)*e4 + 15 * (3 + c)*e6)*s2) / 262144.;
	m_ip9[4][1][1] = (-10395 * cn2*(96 * (-1 - 3 * c + c2 + 3 * c3) + 16 * (-14 - 45 * c + 47 * c2 + 96 * c3)*e2 + 10 * (23 + 61 * c + 65 * c2 + 75 * c3)*e4 + 15 * (6 + 17 * c + 5 * c2)*e6)*s) / 131072.;
	m_ip9[4][1][2] = (-945 * cn2*(48 * (-1 - 22 * c - 32 * c2 + 22 * c3 + 33 * c4) + 16 * (-7 - 184 * c - 251 * c2 + 418 * c3 + 528 * c4)*e2 + 5 * (23 + 346 * c + 592 * c2 + 902 * c3 + 825 * c4)*e4 + 15 * (3 + 56 * c + 87 * c2 + 22 * c3)*e6)) / 131072.;
	m_ip9[4][1][3] = (945 * cn*(96 * (3 - 5 * c - 58 * c2 - 50 * c3 + 55 * c4 + 55 * c5) + 16 * (42 - 103 * c - 1139 * c2 - 805 * c3 + 1925 * c4 + 1760 * c5)*e2 + 10 * (-69 + 27 * c + 462 * c2 + 870 * c3 + 1815 * c4 + 1375 * c5)*e4 + 15 * (-18 + 19 * c + 239 * c2 + 265 * c3 + 55 * c4)*e6)*s) / (131072.*cp);
	m_ip9[4][1][4] = (-315 * cn*(48 * (17 + 91 * c - 198 * c2 - 570 * c3 + 165 * c4 + 495 * c5) + 16 * (119 + 748 * c - 2142 * c2 - 6420 * c3 + 3135 * c4 + 7920 * c5)*e2 + 5 * (-391 - 1501 * c + 522 * c2 + 150 * c3 + 6765 * c4 + 12375 * c5)*e4 + 15 * (-51 - 236 * c + 342 * c2 + 900 * c3 + 165 * c4)*e6)) / 262144.;
	m_ip9[4][1][5] = (-315 * (96 * (-1 + 18 * c + 37 * c2 - 84 * c3 - 135 * c4 + 66 * c5 + 99 * c6) + 16 * (-14 + 309 * c + 686 * c2 - 2040 * c3 - 3312 * c4 + 2211 * c5 + 3168 * c6)*e2 + 10 * (23 - 262 * c - 403 * c2 - 372 * c3 - 687 * c4 + 1914 * c5 + 2475 * c6)*e4 + 15 * (6 - 89 * c - 166 * c2 + 216 * c3 + 336 * c4 + 33 * c5)*e6)*s) / (65536.*cp);
	m_ip9[4][1][6] = (-315 * (48 * (-1 + 19 * c2 - 51 * c4 + 33 * c6) + 16 * (-7 + 190 * c2 - 663 * c4 + 528 * c6)*e2 + 5 * (23 - 133 * c2 - 459 * c4 + 825 * c6)*e4 + 15 * (3 - 38 * c2 + 51 * c4)*e6)) / 65536.;
	m_ip9[4][1][7] = (315 * (96 * (-1 - 18 * c + 37 * c2 + 84 * c3 - 135 * c4 - 66 * c5 + 99 * c6) + 16 * (-14 - 309 * c + 686 * c2 + 2040 * c3 - 3312 * c4 - 2211 * c5 + 3168 * c6)*e2 + 10 * (23 + 262 * c - 403 * c2 + 372 * c3 - 687 * c4 - 1914 * c5 + 2475 * c6)*e4 - 15 * (-6 - 89 * c + 166 * c2 + 216 * c3 - 336 * c4 + 33 * c5)*e6)*s) / (65536.*cn);
	m_ip9[4][1][8] = (315 * cp*(48 * (-17 + 91 * c + 198 * c2 - 570 * c3 - 165 * c4 + 495 * c5) + 16 * (-119 + 748 * c + 2142 * c2 - 6420 * c3 - 3135 * c4 + 7920 * c5)*e2 + 5 * (391 - 1501 * c - 522 * c2 + 150 * c3 - 6765 * c4 + 12375 * c5)*e4 - 15 * (-51 + 236 * c + 342 * c2 - 900 * c3 + 165 * c4)*e6)) / 262144.;
	m_ip9[4][1][9] = (945 * cp*(96 * (-3 - 5 * c + 58 * c2 - 50 * c3 - 55 * c4 + 55 * c5) + 16 * (-42 - 103 * c + 1139 * c2 - 805 * c3 - 1925 * c4 + 1760 * c5)*e2 + 10 * (69 + 27 * c - 462 * c2 + 870 * c3 - 1815 * c4 + 1375 * c5)*e4 + 15 * (18 + 19 * c - 239 * c2 + 265 * c3 - 55 * c4)*e6)*s) / (131072.*cn);
	m_ip9[4][1][10] = (-945 * cp2*(48 * (-1 + 22 * c - 32 * c2 - 22 * c3 + 33 * c4) + 16 * (-7 + 184 * c - 251 * c2 - 418 * c3 + 528 * c4)*e2 + 5 * (23 - 346 * c + 592 * c2 - 902 * c3 + 825 * c4)*e4 - 15 * (-3 + 56 * c - 87 * c2 + 22 * c3)*e6)) / 131072.;
	m_ip9[4][1][11] = (-10395 * cp2*(96 * (1 - 3 * c - c2 + 3 * c3) + 16 * (14 - 45 * c - 47 * c2 + 96 * c3)*e2 + 10 * (-23 + 61 * c - 65 * c2 + 75 * c3)*e4 - 15 * (6 - 17 * c + 5 * c2)*e6)*s) / 131072.;
	m_ip9[4][1][12] = (-10395 * cp2*s2*(16 * (-7 - 3 * c + 16 * c2)*e2 + 5 * (23 - 16 * c + 25 * c2)*e4 - 15 * (-3 + c)*e6 - 48 * s2)) / 262144.;

	m_ip9[4][2][0] = (-22869 * cn4*((-11 - 20 * c - 19 * c2)*e2 + (-9 - 6 * c)*e4 + 20 * s2)) / 131072.;
	m_ip9[4][2][1] = (22869 * cn3*(40 * (-2 - 3 * c + 2 * c2 + 3 * c3) + 2 * (22 + 58 * c + 88 * c2 + 57 * c3)*e2 + (36 + 69 * c + 30 * c2)*e4)*s) / (65536.*cp);
	m_ip9[4][2][2] = (-2079 * cn3*(20 * (13 + 31 * c - 11 * c2 - 33 * c3) + (-143 - 381 * c - 649 * c2 - 627 * c3)*e2 - 3 * (39 + 97 * c + 44 * c2)*e4)) / 65536.;
	m_ip9[4][2][3] = (10395 * cn2*(-40 * (-2 - 11 * c - 9 * c2 + 11 * c3 + 11 * c4) - 2 * (22 + 96 * c + 194 * c2 + 319 * c3 + 209 * c4)*e2 - 3 * (12 + 61 * c + 73 * c2 + 22 * c3)*e4)*s) / (65536.*cp);
	m_ip9[4][2][4] = (-10395 * cn2*(20 * (1 + 22 * c + 32 * c2 - 22 * c3 - 33 * c4) + (-11 - 142 * c - 262 * c2 - 638 * c3 - 627 * c4)*e2 - 3 * (3 + 56 * c + 87 * c2 + 22 * c3)*e4)) / 131072.;
	m_ip9[4][2][5] = (-2079 * cn*(-40 * (2 - 11 * c - 35 * c2 + 11 * c3 + 33 * c4) - 2 * (-22 + 56 * c + 70 * c2 + 319 * c3 + 627 * c4)*e2 - 3 * (-12 + 53 * c + 147 * c2 + 22 * c3)*e4)*s) / 32768.;
	m_ip9[4][2][6] = (2079 * (20 * (1 - 12 * c2 + 11 * c4) - (11 - 2 * c2 - 209 * c4)*e2 - (9 - 69 * c2)*e4)*s2) / 32768.;
	m_ip9[4][2][7] = (2079 * cp*(-40 * (2 + 11 * c - 35 * c2 - 11 * c3 + 33 * c4) + 2 * (22 + 56 * c - 70 * c2 + 319 * c3 - 627 * c4)*e2 + (36 + 159 * c - 441 * c2 + 66 * c3)*e4)*s) / 32768.;
	m_ip9[4][2][8] = (-10395 * cp2*(20 * (1 - 22 * c + 32 * c2 + 22 * c3 - 33 * c4) + (-11 + 142 * c - 262 * c2 + 638 * c3 - 627 * c4)*e2 + 3 * (-3 + 56 * c - 87 * c2 + 22 * c3)*e4)) / 131072.;
	m_ip9[4][2][9] = (-10395 * cp2*(-40 * (-2 + 11 * c - 9 * c2 - 11 * c3 + 11 * c4) + (-44 + 192 * c - 388 * c2 + 638 * c3 - 418 * c4)*e2 + 3 * (-12 + 61 * c - 73 * c2 + 22 * c3)*e4)*s) / (65536.*cn);
	m_ip9[4][2][10] = (2079 * cp3*(20 * (-13 + 31 * c + 11 * c2 - 33 * c3) + (143 - 381 * c + 649 * c2 - 627 * c3)*e2 + 3 * (39 - 97 * c + 44 * c2)*e4)) / 65536.;
	m_ip9[4][2][11] = (-22869 * cp3*(-40 * (2 - 3 * c - 2 * c2 + 3 * c3) + 2 * (22 - 58 * c + 88 * c2 - 57 * c3)*e2 + (36 - 69 * c + 30 * c2)* e4)*s) / (65536.*cn);
	m_ip9[4][2][12] = (-22869 * cp4*((-11 + 20 * c - 19 * c2)*e2 - (9 - 6 * c)*e4 + 20 * s2)) / 131072.;

	m_ip9[4][3][0] = (297297 * cn5*(1 - c - e2)) / 262144.;
	m_ip9[4][3][1] = (297297 * cn4*s*(-((6 + 5 * c)*e2) + 6 * s2)) / (131072.*cp);
	m_ip9[4][3][2] = (297297 * cn4*(-((3 + 2 * c)*e2) + 3 * s2)) / 131072.;
	m_ip9[4][3][3] = (1486485 * cn3*s*((2 + c)*e2 - 2 * s2)) / 131072.;
	m_ip9[4][3][4] = (1486485 * cn2*s2*(-((3 + c)*e2) + 3 * s2)) / 262144.;
	m_ip9[4][3][5] = (297297 * cn*s3*(-((6 + c)*e2) + 6 * s2)) / 65536.;
	m_ip9[4][3][6] = (297297 * s4*(-e2 + s2)) / 65536.;
	m_ip9[4][3][7] = (297297 * cp*s3*((6 - c)*e2 - 6 * s2)) / 65536.;
	m_ip9[4][3][8] = (1486485 * cp2*s2*(-((3 - c)*e2) + 3 * s2)) / 262144.;
	m_ip9[4][3][9] = (1486485 * cp3*s*(-((2 - c)*e2) + 2 * s2)) / 131072.;
	m_ip9[4][3][10] = (297297 * cp4*(-((3 - 2 * c)*e2) + 3 * s2)) / 131072.;
	m_ip9[4][3][11] = (297297 * cp4*s*((6 - 5 * c)*e2 - 6 * s2)) / (131072.*cn);
	m_ip9[4][3][12] = (297297 * cp5*(1 + c - e2)) / 262144.;

	m_ip9[5][0][0] = (-3465 * c*s4) / 131072.;
	m_ip9[5][0][1] = (3465 * (1 - 6 * c2)*s3) / 65536.;
	m_ip9[5][0][2] = (315 * c*(13 - 33 * c2)*s2) / 65536.;
	m_ip9[5][0][3] = (1575 * (1 - 15 * c2 + 22 * c4)*s) / 65536.;
	m_ip9[5][0][4] = (-525 * c*(19 - 102 * c2 + 99 * c4)) / 131072.;
	m_ip9[5][0][5] = (105 * (5 - 100 * c2 + 285 * c4 - 198 * c6)) / (32768.*s);
	m_ip9[5][0][6] = (105 * c*(5 - 30 * c2 + 33 * c4)) / 32768.;
	m_ip9[5][0][7] = (105 * (5 - 100 * c2 + 285 * c4 - 198 * c6)) / (32768.*s);
	m_ip9[5][0][8] = (-525 * c*(19 - 102 * c2 + 99 * c4)) / 131072.;
	m_ip9[5][0][9] = (1575 * (1 - 15 * c2 + 22 * c4)*s) / 65536.;
	m_ip9[5][0][10] = (315 * c*(13 - 33 * c2)*s2) / 65536.;
	m_ip9[5][0][11] = (3465 * (1 - 6 * c2)*s3) / 65536.;
	m_ip9[5][0][12] = (-3465 * c*s4) / 131072.;

	m_ip9[5][1][0] = (-10395 * cn2*(1 + 3 * c)*s2) / 262144.;
	m_ip9[5][1][1] = (10395 * cn2*(1 - 11 * c - 18 * c2)*s) / 131072.;
	m_ip9[5][1][2] = (945 * cn2*(10 + 9 * c - 88 * c2 - 99 * c3)) / 131072.;
	m_ip9[5][1][3] = (-945 * cn*(11 + 109 * c + 35 * c2 - 385 * c3 - 330 * c4)*s) / (131072.*cp);
	m_ip9[5][1][4] = (-315 * cn*(37 - 252 * c - 810 * c2 + 660 * c3 + 1485 * c4)) / 262144.;
	m_ip9[5][1][5] = (-315 * (19 + 56 * c - 288 * c2 - 474 * c3 + 429 * c4 + 594 * c5)*s) / (65536.*cp);
	m_ip9[5][1][6] = (-315 * c*(19 - 102 * c2 + 99 * c4)) / 65536.;
	m_ip9[5][1][7] = (-315 * (19 - 56 * c - 288 * c2 + 474 * c3 + 429 * c4 - 594 * c5)*s) / (65536.*cn);
	m_ip9[5][1][8] = (315 * cp*(37 + 252 * c - 810 * c2 - 660 * c3 + 1485 * c4)) / 262144.;
	m_ip9[5][1][9] = (-945 * cp*(11 - 109 * c + 35 * c2 + 385 * c3 - 330 * c4)*s) / (131072.*cn);
	m_ip9[5][1][10] = (-945 * cp2*(10 - 9 * c - 88 * c2 + 99 * c3)) / 131072.;
	m_ip9[5][1][11] = (10395 * cp2*(1 + 11 * c - 18 * c2)*s) / 131072.;
	m_ip9[5][1][12] = (10395 * (1 - 3 * c)*cp2*s2) / 262144.;

	m_ip9[5][2][0] = (-22869 * cn4*(2 + 3 * c)) / 131072.;
	m_ip9[5][2][1] = (-22869 * cn3*(5 + 22 * c + 18 * c2)*s) / (65536.*cp);
	m_ip9[5][2][2] = (-2079 * cn3*(4 + 77 * c + 99 * c2)) / 65536.;
	m_ip9[5][2][3] = (-10395 * cn2*(5 - 19 * c - 88 * c2 - 66 * c3)*s) / (65536.*cp);
	m_ip9[5][2][4] = (10395 * cn2*(10 + 9 * c - 88 * c2 - 99 * c3)) / 131072.;
	m_ip9[5][2][5] = (2079 * cn*(13 + 63 * c - 88 * c2 - 198 * c3)*s) / 32768.;
	m_ip9[5][2][6] = (2079 * c*(13 - 33 * c2)*s2) / 32768.;
	m_ip9[5][2][7] = (2079 * cp*(13 - 63 * c - 88 * c2 + 198 * c3)*s) / 32768.;
	m_ip9[5][2][8] = (-10395 * cp2*(10 - 9 * c - 88 * c2 + 99 * c3)) / 131072.;
	m_ip9[5][2][9] = (-10395 * cp2*(5 + 19 * c - 88 * c2 + 66 * c3)*s) / (65536.*cn);
	m_ip9[5][2][10] = (2079 * cp3*(4 - 77 * c + 99 * c2)) / 65536.;
	m_ip9[5][2][11] = (-22869 * cp3*(5 - 22 * c + 18 * c2)*s) / (65536.*cn);
	m_ip9[5][2][12] = (22869 * (2 - 3 * c)*cp4) / 131072.;

	m_ip9[5][3][0] = (-297297 * cn5) / 262144.;
	m_ip9[5][3][1] = (-297297 * cn4*(5 + 6 * c)*s) / (131072.*cp);
	m_ip9[5][3][2] = (-297297 * cn4*(2 + 3 * c)) / 131072.;
	m_ip9[5][3][3] = (1486485 * cn3*(1 + 2 * c)*s) / 131072.;
	m_ip9[5][3][4] = (-1486485 * cn2*(1 + 3 * c)*s2) / 262144.;
	m_ip9[5][3][5] = (-297297 * cn*(1 + 6 * c)*s3) / 65536.;
	m_ip9[5][3][6] = (-297297 * c*s4) / 65536.;
	m_ip9[5][3][7] = (-297297 * (1 - 6 * c)*cp*s3) / 65536.;
	m_ip9[5][3][8] = (1486485 * (1 - 3 * c)*cp2*s2) / 262144.;
	m_ip9[5][3][9] = (-1486485 * cp3*(-1 + 2 * c)*s) / 131072.;
	m_ip9[5][3][10] = (297297 * (2 - 3 * c)*cp4) / 131072.;
	m_ip9[5][3][11] = (-297297 * (5 - 6 * c)*cp4*s) / (131072.*cn);
	m_ip9[5][3][12] = (297297 * cp5) / 262144.;
}

///// LUNISOLAR PERTURBATIONS /////

template<typename T>
void HEOSATDynamics<T>::sumord5(const double t, const T& g, const T& h, const T& a, const T& n, const T& e, const T& y, const T& ci, const T& si, DACE::AlgebraicVector<T>& dx)
{
	
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
	double ks, km;// pm, ps, qm, qs, cx, sx, pp, qq;
	double ss[5], cs[5], sm[5], cm[5];
	T ddx[6];

	const T pp = a*a;
	const T qq = 1 / n;

	// Compute eccentricity and inclination functions
	eccincpol5(e, y, ci, si);


	coefsolun5(ss, cs, &ks, sm, cm, &km);

	for (int k = 1; k <= 5; ++k)
	{
		T pm = 0.0;
		T ps = 0.0;
		for (int i = 0; i <= 1; ++i)
		{
			T qm = 0.0;
			T qs = 0.0;
			for (int j = -2; j <= 2; ++j)
			{
				T cx = m_ip5[k][i][j + 2] * cos(2 * i*g + j*h);
				T sx = m_ip5[k][i][j + 2] * sin(2 * i*g + j*h);
				if (k < 3){
					qm = qm + (cm[j + 2] * sx - sm[j + 2] * cx);
					qs = qs + (cs[j + 2] * sx - ss[j + 2] * cx);
				}
				else
				{
					qm = qm + (cm[j + 2] * cx + sm[j + 2] * sx);
					qs = qs + (cs[j + 2] * cx + ss[j + 2] * sx);
				}
			}
			pm = pm + m_ep5[k][i] * qm;
			ps = ps + m_ep5[k][i] * qs;
		}
		if (k < 3)
		{
			ddx[k] = pp*(km*pm + ks*ps);
		}
		else
		{
			ddx[k] = qq*(km*pm + ks*ps);
		}
	}
	dx[0] = ddx[3];
	dx[1] = ddx[4];
	dx[2] = ddx[5];
	dx[3] = 0;
	dx[4] = ddx[1];
	dx[5] = ddx[2];
}


template<typename T>
void HEOSATDynamics<T>::sumord6(const double t, const T& g, const T& h, const T& a, const T& n, const T& e, const T& y, const T& ci, const T& si, DACE::AlgebraicVector<T>& dx)
{
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
	double km;
	double sm[7], cm[7];
	T ddx[6];

	const T pp = a*a*a;
	const T qq = a / n;

	eccincpol6(e, y, ci, si);
	coeflun6(cm, sm, &km);

	for (int k = 1; k <= 5; ++k)
	{
		T pm = 0.0;
		for (int i = 0; i <= 1; ++i)
		{
			T qm = 0.0;
			for (int j = -3; j <= 3; ++j)
			{
				T cx = m_ip6[k][i][j + 3] * cos((2 * i + 1)*g + j*h);
				T sx = m_ip6[k][i][j + 3] * sin((2 * i + 1)*g + j*h);
				if (k < 3)
				{
					qm = qm + (cm[j + 3] * sx + sm[j + 3] * cx);
				}
				else
				{
					qm = qm + (sm[j + 3] * sx - cm[j + 3] * cx);
				}
			}
			pm = pm + m_ep6[k][i] * qm;
		}
		if (k < 3)
		{
			ddx[k] = pp*km*pm;
		}
		else
		{
			ddx[k] = qq*km*pm;
		}
	}

	coefsun6(cm, sm, &km);

	for (int k = 1; k <= 5; ++k)
	{
		T pm = 0.0;
		for (int i = 0; i <= 1; ++i)
		{
			T qm = 0.0;
			for (int j = -3; j <= 3; ++j)
			{
				T cx = m_ip6[k][i][j + 3] * cos((2 * i + 1)*g + j*h);
				T sx = m_ip6[k][i][j + 3] * sin((2 * i + 1)*g + j*h);
				if (k < 3)
				{
					qm = qm + (cm[j + 3] * sx + sm[j + 3] * cx);
				}
				else
				{
					qm = qm + (sm[j + 3] * sx - cm[j + 3] * cx);
				}
			}
			pm = pm + m_ep6[k][i] * qm;
		}
		if (k < 3)
		{
			ddx[k] = ddx[k] + pp*km*pm;
		}
		else
		{
			ddx[k] = ddx[k] + qq*km*pm;
		}
	}

	dx[0] = ddx[3];
	dx[1] = ddx[4];
	dx[2] = ddx[5];
	dx[3] = 0;
	dx[4] = ddx[1];
	dx[5] = ddx[2];
}

template<typename T>
void HEOSATDynamics<T>::sumord7(const double t, const T& g, const T& h, const T& a, const T& n, const T& e, const T& y, const T& ci, const T& si, DACE::AlgebraicVector<T>& dx)
{
	
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;

	
	double km;
	double sm[9], cm[9];
	T ddx[6];

	const T pp = a*a*a*a;
	const T qq = a*a / n;

	eccincpol7(e, y, ci, si);
	coeflun7(cm, sm, &km);

	for (int k = 1; k <= 5; ++k)
	{
		T pm = 0.0;
		for (int i = 0; i <= 2; ++i)
		{
			T qm = 0.0;
			for (int j = -4; j <= 4; ++j)
			{
				const T cx = m_ip7[k][i][j + 4] * cos(2 * i*g + j*h);
				const T sx = m_ip7[k][i][j + 4] * sin(2 * i*g + j*h);
				if (k < 3)
				{
					qm = qm + (cm[j + 4] * sx - sm[j + 4] * cx);
				}
				else
				{
					qm = qm + (cm[j + 4] * cx + sm[j + 4] * sx);
				}
			}
			pm = pm + m_ep7[k][i] * qm;
		}
		if (k < 3)
		{
			ddx[k] = pp*km*pm;
		}
		else
		{
			ddx[k] = qq*km*pm;
		}
	}

	dx[0] = ddx[3];
	dx[1] = ddx[4];
	dx[2] = ddx[5];
	dx[3] = 0;
	dx[4] = ddx[1];
	dx[5] = ddx[2];
}


template<typename T>
void HEOSATDynamics<T>::sumord8(const double t, const T& g, const T& h, const T& a, const T& n, const T& e, const T& y, const T& ci, const T& si, DACE::AlgebraicVector<T>& dx)
{

	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
	double km;
	double sm[11], cm[11];
	T ddx[6];

	const T pp = a*a*a*a*a;
	const T qq = a*a*a / n;

	eccincpol8(e, y, ci, si);
	coeflun8(cm, sm, &km);

	for (int k = 1; k <= 5; ++k)
	{
		T pm = 0.0;
		for (int i = 0; i <= 2; ++i)
		{
			T qm = 0.0;
			for (int j = -5; j <= 5; ++j)
			{
				const T cx = m_ip8[k][i][j + 5] * cos((2 * i + 1)*g + j*h);
				const T sx = m_ip8[k][i][j + 5] * sin((2 * i + 1)*g + j*h);
				if (k < 3)
				{
					qm = qm + (cm[j + 5] * sx + sm[j + 5] * cx);
				}
				else
				{
					qm = qm + (sm[j + 5] * sx - cm[j + 5] * cx);
				}
			}
			pm = pm + m_ep8[k][i] * qm;
		}
		if (k < 3)
		{
			ddx[k] = pp*km*pm;
		}
		else
		{
			ddx[k] = qq*km*pm;
		}
	}

	dx[0] = ddx[3];
	dx[1] = ddx[4];
	dx[2] = ddx[5];
	dx[3] = 0;
	dx[4] = ddx[1];
	dx[5] = ddx[2];
}

template<typename T>
void HEOSATDynamics<T>::sumord9(const double t, const T& g, const T& h, const T& a, const T& n, const T& e, const T& y, const T& ci, const T& si, DACE::AlgebraicVector<T>& dx)
{

	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
	double km;
	double sm[13], cm[13];
	T ddx[6];

	const T pp = a*a*a*a*a*a; //!
	const T qq = a*a*a*a / n; //!

	eccincpol9(e, y, ci, si);
	coeflun9(cm, sm, &km);

	for (int k = 1; k <= 5; ++k)
	{
		T pm = 0.0;
		for (int i = 0; i <= 3; ++i)
		{
			T qm = 0.0;
			for (int j = -6; j <= 6; ++j)
			{
				const T cx = m_ip9[k][i][j + 6] * cos(2 * i*g + j*h);
				const T sx = m_ip9[k][i][j + 6] * sin(2 * i*g + j*h);
				if (k < 3)
				{
					qm = qm + (sm[j + 6] * cx + cm[j + 6] * sx);
				}
				else
				{
					qm = qm + (cm[j + 6] * cx - sm[j + 6] * sx);
				}
			}
			pm = pm + m_ep9[k][i] * qm;
		}
		
		if(k < 3) //!
		{ //!
			ddx[k] = pp*km*pm; //!
		} //!
		else //!
		{ //!
			ddx[k] = qq*km*pm; //!
		} //!
		
		//ddx[k] = km*pm; //!
	}

	dx[0] = ddx[3];
	dx[1] = ddx[4];
	dx[2] = ddx[5];
	dx[3] = 0;
	dx[4] = ddx[1];
	dx[5] = ddx[2];
}

///// SUN AND MOON POSITIONS /////

template<typename T>
void HEOSATDynamics<T>::meeussun(const double t, double *sz, double *Rs, double *Ms)
{
    /*
     Sun coordinates:
     Accuracy of 0.01 degrees
     
     Reference:
     Meeus (1998) - Astronomical Algorithms
     */
	double md, Ls, eS, Cs;

    // Geometric mean longitude of the Sun referred to mean equinox of the date in degrees
	Ls = std::fmod(280.46646 + (36000.76983 + 0.0003032*t)*t, 360.0);
    // Mean anomaly of the Sun in degrees
	*Ms = std::fmod(357.5291092 + t*(35999.0502909 - t*(0.0001536 - t / 2.449e7)), 360.0);
    // Mean anomaly of the Sun in radians
	md = *Ms*DEG;
    // Eccentricity of the Earth's orbit
	eS = 0.016708634 - (0.000042037 + 0.0000001267*t)*t;
    // Sun's equation of the center C
	Cs = (1.914602 - (0.004817 + 0.000014*t)*t)*std::sin(md) + (0.019993 - 0.000101*t)*std::sin(2 * md) +
		0.000289*std::sin(3 * md);
    // True anomaly of the Sun in radians
	*sz = std::fmod((*Ms + Cs)*DEG, 2 * M_PI);
    // Radius vector of the Sun in astronomical units
	*Rs = (1 - eS*eS) / (1 + eS*std::cos(*sz));
    // True longitude of the Sun
	*sz = std::fmod(Ls + Cs, 360.0);
}

template<typename T>
void HEOSATDynamics<T>::meeusmoon(const double t, const double Ms, const int j, double *le, double *be, double *rm)
{
    /*
     Moon coordinates:
     Accuracy of 10" in longitude and 4" in latitude of the Moon
     
     References:
     Meeus (1998) - Astronomical Algorithms
     Chapront, J, Francou, G (2003) The lunar theory ELP revisited. introduction of new planetary perturbations. Astronomy and Astrophysics 404(2):735–742
     Chapront-Touze M, Chapront J (1988) ELP 2000-85 - A semi-analytical lunar ephemeris adequate for historical times. Astronomy and Astrophysics 190:342–352
     */
	double Lm, Dm, Mm, Fm, a1, a2, a3, ee, oo;

    // Mean longitude of the Moon referred to mean equinox of date (and including constant term of effect of light time (-0".70)
	Lm = std::fmod(218.3164477 + t*(481267.88123421 - t*(0.0015786 - (1. / 538841. - t / 65194000.)*t)), 360.);
    // Mean elongation of the Moon
	Dm = std::fmod(297.8501921 + t*(445267.1114034 - t*(0.0018819 - (1. / 545868. - t / 113065000.)*t)), 360.);
    // Mean anomaly of the Moon
	Mm = std::fmod(134.9633964 + t*(477198.8675055 + t*(0.0087414 + (1. / 69699. - t / 14712000.)*t)), 360.);
    // Argument of latitude of the Moon
	Fm = std::fmod(93.272095 + t*(483202.0175233 - t*(0.0036539 + (1. / 3526000. - t / 863310000.)*t)), 360.);

	a1 = std::fmod(119.75 + 131.849*t, 360.);
	a2 = std::fmod(53.09 + 479264.290*t, 360.);
	a3 = std::fmod(313.45 + 481266.484*t, 360.);

	*le = 3958 * std::sin(a1*DEG) + 1962 * std::sin((Lm - Fm)*DEG) + 318 * std::sin(a2*DEG);
	*be = -2235 * std::sin(Lm*DEG) + 382 * std::sin(a3*DEG) + 175 * std::sin((a1 - Fm)*DEG) + 175 * std::sin((a1 + Fm)*DEG) +
		127 * std::sin((Lm - Mm)*DEG) - 115 * std::sin((Lm + Mm)*DEG);
	*rm = 0.0;

	ee = 1.0 - (0.002516 + 0.0000074*t)*t;
	for (int i = 0; i < j; ++i)
	{
		a1 = (lr[i][0] * Dm + lr[i][1] * Ms + lr[i][2] * Mm + lr[i][3] * Fm)*DEG;
		oo = std::pow(ee, abs(lr[i][1]));
		*le = *le + oo*lr[i][4] * std::sin(a1);
		*rm = *rm + oo*lr[i][5] * std::cos(a1);
		a2 = (ab[i][0] * Dm + ab[i][1] * Ms + ab[i][2] * Mm + ab[i][3] * Fm)*DEG;
		*be = *be + oo*ab[i][4] * std::sin(a2);
	}
    // Geocentric longitude of the Moon in degrees
	*le = std::fmod(Lm + *le*1.e-6, 360.);
    // Geocentric latitude of the Moon in degrees
	*be = *be*1.e-6;
    // Distance between Earth and Moon centers in Moon semi-major axes
	*rm = (385000.56 + *rm*1.e-3) / 384400.;

}

// Compute Sun and Moon positions
template<typename T>
void HEOSATDynamics<T>::elsetsumo(const double jd)
{
    /*
     jd = Julian date in Barycentric Dynamical Time (TDB)
     */
	double e0, t, Ms, le, be, th, de, dy;
	const int i = 60;

    // Centuries since epoch J2000 (2000 January 1.5)
    // TODO: replace by call to "tocent"
	t = (jd - 2451545.) / 36525.;
	dy = 0;
	de = 0;

    // Mean obliquity of the ecliptic (Ref. Meeus)
	e0 = 23.0 + (26.0 + (21.448 - (46.815 + (0.00059 - 0.001813*t)*t)*t) / 60.) / 60;
	e0 = (e0 + de)*DEG;

    // Sun's position
	meeussun(t, &th, &m_Rs, &Ms);

	th = (th + dy)*DEG;
	m_xs = std::cos(th);
	m_ys = std::cos(e0)*std::sin(th);
	m_zs = std::sin(e0)*std::sin(th);

    // Moon's position
	meeusmoon(t, Ms, i, &le, &be, &m_rm);

	le = (le + dy)*DEG;
	be = be*DEG;

	m_xl = std::cos(be)*std::cos(le);
	m_yl = std::cos(be)*std::sin(le)*std::cos(e0) - std::sin(be)*std::sin(e0);
	m_zl = std::cos(be)*std::sin(le)*std::sin(e0) + std::sin(be)*std::cos(e0);
}

// Compute Sun and Moon positions using SPICE
template<typename T>
void HEOSATDynamics<T>::sunMoonSpice(const double jd)
{
    /*
     jd = Julian date in Barycentric Dynamical Time (TDB)
     */
    
    // TDB seconds past the J2000 epoch
	const double et = unitim_c(jd, "JED", "ET");

	double aux[6];
	double lt;

	spkezr_c("SUN", et, "TOD", "NONE", "EARTH", aux, &lt);
	m_Rs = std::sqrt(std::pow(aux[0], 2) + std::pow(aux[1], 2) + std::pow(aux[2], 2));
	m_xs = aux[0] / m_Rs;
	m_ys = aux[1] / m_Rs;
	m_zs = aux[2] / m_Rs;
	m_Rs = m_Rs / m_ul / m_as;

	spkezr_c("MOON", et, "TOD", "NONE", "EARTH", aux, &lt);
	m_rm = std::sqrt(std::pow(aux[0], 2) + std::pow(aux[1], 2) + std::pow(aux[2], 2));
	m_xl = aux[0] / m_rm;
	m_yl = aux[1] / m_rm;
	m_zl = aux[2] / m_rm;
	m_rm = m_rm / m_ul / m_am;
}

///// SOLAR RADIATION PRESSURE /////

template<typename T>
T HEOSATDynamics<T>::accelsrp(const double ms)
{
	double q;

	q = ms + 1.92*std::sin(ms);
	q = (1 + 0.01672*std::cos(q)) / 0.99972;

	const double P = 4.65e-6; // the solar radiation pressure constant at one AU
	return ((m_SRPC*P / 1000.0)*m_ut*m_ut / m_ul*q*q);
}

template<typename T>
void HEOSATDynamics<T>::srpaverage(const double t, const T& g, const T& h, const T& a, const T& n, const T& e, const T& y, const T& ci, const T& si, DACE::AlgebraicVector<T>& x)
{
	
	//Using std and DACE functions
	using std::sqrt;  using DACE::sqrt;
	
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
	const double d = (ut2dj(t) - 2400000.5) - 15019.5;
	const double ms = (358.48 + 0.98560027*d)*DEG;
	const double l = (279.70 + 0.98564734*d + 1.92*sin(ms))*DEG;

	const T Fsrp = accelsrp(ms);
    
	const double ee = 23.44*DEG;
	const double k = cos(ee);
	const double z = sin(ee);

	const T de = Fsrp*y / (n*a)*(3 * ((1 - ci)*(1 - k)*sin(g - h - l) + (1 + ci)*(1 + k)*sin(g + h - l) + 2 * si*z*(sin(g - l) - sin(g + l)) + (1 - ci)*(1 + k)*sin(g - h + l) + (1 + ci)*(1 - k)*sin(g + h + l))) / 8.;

	const T di = Fsrp / (n*a*y)*(3 * e*(2 * ci*z*(sin(g + l) - sin(g - l)) + (1 + k)*si*(sin(g + h - l) - sin(g - h + l)) + (1 - k)*si*(sin(g + h + l) - sin(g - h - l)))) / 8.;
	const T ll = sqrt(a);
	const T gg = ll*y;

	x[0] =-Fsrp/(n*a)*(e + 1./e)*3*(2*si*z*(cos(g - l) - cos(g + l)) + ci*(1 + k)*(cos(g + h - l) - cos(g - h + l)) + 
		(1 + k)*(cos(g + h - l) + cos(g - h + l)) + ci*(1 - k)*(-cos(g - h - l) + cos(g + h + l)) + (1 - k)*(cos(g - h - l) + cos(g + h + l)))/8.;
	x[2] = Fsrp / (n*a*y*si)*(3 * e*(2 * ci*z*(cos(g - l) - cos(g + l)) + (1 + k)*si*(-cos(g + h - l) + cos(g - h + l)) + (1 - k)*si*(cos(g - h - l) - cos(g + h + l)))) / 8.;
	x[1] = -ci*x[2] + Fsrp*y / (n*a*e)*(3 * (2 - 3 * e*e)*((1 - ci)*(1 - k)*cos(g - h - l) + (1 + ci)*(1 + k)*cos(g + h - l) + 2 * si*z*(cos(g - l) - cos(g + l)) + (1 - ci)*(1 + k)*cos(g - h + l) + (1 + ci)*(1 - k)*cos(g + h + l))) / (16.*y*y);
	x[3] = 0;
	x[4] = -(e / y)*ll*de;
	x[5] = x[4] * ci - (gg*si)*di;
}

///// TESSERAL PERTURBATIONS /////

template<typename T>
void HEOSATDynamics<T>::tesre2t1(const double t, const T& l, const T& g, const T& hh, const T& a, const T& n, const T& e, const T& y, const T& ci, DACE::AlgebraicVector<T>& dx)
{
	//Using std and DACE functions
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
    // 2:1 tesseral resonance (for 12-hour orbits)
	const T f222 = (3. * (-1. + ci)*(-1. + ci)) / 4.;
	const T f221 = (-3. * (-1. + ci*ci)) / 2.;
	const T f220 = (3. * (1. + ci)*(1. + ci)) / 4.;
	const T df222 = (3. * (-1. + ci)) / 2.;
	const T df221 = -3. * ci;
	const T df220 = (3. * (1. + ci)) / 2.;

	const T e2 = e*e;
	const T e3 = e2*e;
	const T e5 = e3*e2;
	const T e7 = e5*e2;
	const T e9 = e7*e2;
	const T e11 = e9*e2;
	const T e13 = e11*e2;
//    const T e15 = e13*e2;
//    const T e17 = e15*e2;
//    const T e19 = e17*e2;
//    const T e21 = e19*e2;


	const T g223 = 0.020833333333333332*e3 + 0.014322916666666666*e5 + 0.010188802083333334*e7 + 0.007584183304398148*e9 + 0.0058915435952484295*e11 + 0.00473463116499482*e13; // 13th-order in ecc
	               //0.020833333333333333333*e3 + 0.014322916666666666667*e5 + 0.010188802083333333333*e7 + 0.0075841833043981481481*e9 + 0.0058915435952484292328*e11 + 0.0047346311649948200852*e13 + 0.0039075560377505207305*e15 + 0.0032938900851694608014*e17 + 0.0028244799844723200460*e19 + 0.0024562667180594775398*e21; // 21th-order in ecc
	const T g211 = 1.5*e + 1.6875*e3 + 2.0390625*e5 + 2.3289388020833335*e7 + 2.587322998046875*e9 + 2.822341054280599*e11 + 3.0393375761042196*e13; // 13th-order in ecc
				   //1.5*e + 1.6875*e3 + 2.0390625*e5 + 2.3289388020833333333*e7 + 2.5873229980468750000*e9 + 2.8223410542805989583*e11 + 3.0393375761042196284*e13 + 3.2418959114597497884*e15 + 3.4325544953467659924*e17 + 3.6131874865348258149*e19 + 3.7852245415940150929*e21; // 21th-order in ecc
	const T g20m1 = -0.5*e + 0.0625*e3 - 0.013020833333333334*e5 - 0.007758246527777778*e7 - 0.006169297960069445*e9 - 0.004967351842809607*e11 - 0.0040929340517289635*e13; // 13th-order in ecc
				    //-0.5*e + 0.0625*e3 + -0.013020833333333333333*e5 + -0.0077582465277777777778*e7 + -0.0061692979600694444444*e9 + -0.0049673518428096064815*e11 + -0.0040929340517289634314*e13 + -0.0034408991171329139792*e15 + -0.0029421149027509295337*e17 + -0.0025515931232383317400*e19 + -0.0022395650453734443393*e21; // 21th-order in ecc
	const T dg223 = -0.0625*e - 0.040364583333333336*e3 - 0.027701822916666667*e5 - 0.01973876953125*e7 - 0.014845635146691055*e9 - 0.011650479205701718*e11 - 0.00944531232766864*e13; // 13th-order in ecc
					//-0.0625*e + -0.040364583333333333333*e3 + -0.027701822916666666667*e5 + -0.01973876953125*e7 + -0.014845635146691054894*e9 + -0.011650479205701716993*e11 + -0.0094453123276686402085*e13 + -0.0078529959327959306182*e15 + -0.0066607070546427145406*e17 + -0.0057413217588207581605*e19 + -0.0050150522111689683133*e21; // 21th-order in ecc
	const T dg211 = -1.5 / e - 4.3125*e - 7.4765625*e3 - 10.478352864583334*e5 - 13.485207112630208*e7 - 16.489000091552736*e9 - 19.491419744441117*e11 - 22.493068288658844*e13; // 13th-order in ecc
					//-1.5 / e + -4.3125*e + -7.4765625*e3 + -10.478352864583333333*e5 + -13.485207112630208333*e7 + -16.489000091552734375*e9 + -19.491419744441118190*e11 + -22.493068288658842028*e13 + -25.494249759892043516*e15 + -28.495130199735706534*e17 + -31.495806921846937863*e19 + -34.496340271234429383*e21; // 21th-order in ecc
	const T dg20m1 = 0.5 / e - 0.4375*e + 0.09635416666666667*e3 + 0.013943142361111112*e5 + 0.012419297960069445*e7 + 0.009673897072120949*e9 + 0.00788292984062612*e11 + 0.006596063160872273*e13; // 13th-order in ecc
					//0.5 / e + -0.4375*e + 0.096354166666666666667*e3 + 0.013943142361111111111*e5 + 0.012419297960069444444*e7 + 0.0096738970721209490741*e9 + 0.0078829298406261195161*e11 + 0.0065960631608722731940*e13 + 0.0056293427373271478570*e15 + 0.0048794294383692764370*e17 + 0.0042831403732847243220*e19 + 0.0037994635077586375822*e21; // 21th-order in ecc

	const double cz = tsmut(ut2dj(t));
	const T f = 2. * (hh - cz) + l;
	const T cf = cos(f);
	const T sf = sin(f);
	T h = 2. * g + f;
	const T cfp = cos(h);
	const T sfp = sin(h);
	h = 2. * g - f;
	const T cfm = cos(h);
	const T sfm = sin(h);

	const T na = n*n*m_al*m_al;
	h = f222*g223*(m_s22*cfm + m_c22*sfm);
	const T fg220 = f220*g20m1*(m_s22*cfp - m_c22*sfp);
	dx[3] = na*(h + f221*g211*(m_s22*cf - m_c22*sf) + fg220);
	dx[4] = 2. * na*(fg220 - h);
	dx[5] = 2. * dx[3];
	const T na2 = n*m_al*m_al / (a*a);
	const T csfm = (m_c22*cfm - m_s22*sfm);
	const T csf = (m_c22*cf + m_s22*sf);
	const T csfp = (m_c22*cfp + m_s22*sfp);
	dx[2] = -(na2 / y)*(df222*g223*csfm + df221*g211*csf + df220*g20m1*csfp);
	dx[1] = -(ci*dx[2]) - na2*(dg223*f222*csfm + dg211*f221*csf + dg20m1*f220*csfp);
	dx[0] = -(y*(dx[1] + ci*dx[2])) + 6. * na2*(f222*g223*csfm + f221*g211*csf + f220*g20m1*csfp);
}

template<typename T>
void HEOSATDynamics<T>::tesre1t1(const double t, const T& l, const T& g, const T& hh, const T& a, const T& n, const T& e, const T& ci, DACE::AlgebraicVector<T>& dx)
{
	//Using std and DACE functions
	using std::pow;   using DACE::pow;

	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;

	
    /* Compute time derivatives of Delaunay elements due to 1:1 tesseral resonance
     t - time
     l - mean anomaly
     g - argument of perigee
     hh - right ascension of ascending node
     a - semi-major axis
     n - mean motion
     e - eccentricity
     y - sqrt(1-e^2)
     ci - cos(inclination)
     si - sin(inclination)
     dx - time derivatives of Delaunay elements
     */
    
    // 1:1 tesseral resonance (for 24-hour orbits)
    // Expanded up to 5th order in eccentricity
    
    // GMST (Greenwich Mean Sidereal Time)
    const double cz = tsmut(ut2dj(t));
    // h = ascending node in rotating reference frame = Omega - omega*t
    const T h = hh - cz;
    
    // dl/dt
    dx[0] = (n*pow(a,-2)*pow(m_al,2)*((2.*pow(-1. + ci,2)*pow(e,2)*(-20. + 29.*pow(e,2))*
                                    (m_c22*cos(2.*g - 2.*h - 2.*l) - m_s22*sin(2.*g - 2.*h - 2.*l)))/5. -
                                3.*(-1. + pow(ci,2))*(-144. + 352.*pow(e,2) + 137.*pow(e,4))*
                                    (m_c22*cos(2.*h + 2.*l) + m_s22*sin(2.*h + 2.*l)) +
                                pow(1. + ci,2)*(528. - 1116.*pow(e,2) + 425.*pow(e,4))*
                                    (m_c22*cos(2.*g + 2.*h + 2.*l) + m_s22*sin(2.*g + 2.*h + 2.*l))))/64.;
    
    // dg/dt
    dx[1] = (n*pow(a,-2)*pow(m_al,2)*((2*(-1. + ci)*pow(e,2)*(20.*(-1. + ci) + (-11. + 21.*ci)*pow(e,2))*
                                    (m_c22*cos(2.*g - 2.*h - 2.*l) - m_s22*sin(2.*g - 2.*h - 2.*l)))/5. -
                                3.*(144.*(-1. + pow(ci,2)) + 8.*(-19. + 37.*pow(ci,2))*pow(e,2) +
                                (-293. + 477.*pow(ci,2))*pow(e,4))*(m_c22*cos(2.*h + 2.*l) + m_s22*sin(2.*h + 2.*l)) +
                                (1. + ci)*(-48.*(5. + 3.*ci) + 12.*(23. + 7.*ci)*pow(e,2) - (83. + 89.*ci)*pow(e,4))*
                                    (m_c22*cos(2.*g + 2.*h + 2.*l) + m_s22*sin(2.*g + 2.*h + 2.*l))))/64.;
    
    // dh/dt
    dx[2] = (pow(m_al,2)*pow(n,2)*(2*pow(-1. + ci,2)*pow(e,4)*
                                    (m_s22*cos(2.*g - 2.*h - 2.*l) + m_c22*sin(2.*g - 2.*h - 2.*l)) +
                                24.*(-1. + pow(ci,2))*pow(e,2)*(9. + 7.*pow(e,2))*
                                    (-(m_s22*cos(2.*h + 2.*l)) + m_c22*sin(2.*h + 2.*l)) -
                                3.*pow(1. + ci,2)*(16. - 40.*pow(e,2) + 13.*pow(e,4))*
                                    (-(m_s22*cos(2.*g + 2.*h + 2.*l)) + m_c22*sin(2.*g + 2.*h + 2.*l))))/32.;
    
    // dL/dt
    dx[3] = (pow(m_al,2)*pow(n,2)*(2.*pow(-1. + ci,2)*pow(e,4)*
                                    (m_s22*cos(2.*g - 2.*h - 2.*l) + m_c22*sin(2.*g - 2.*h - 2.*l)) +
                                24.*(-1. + pow(ci,2))*pow(e,2)*(9. + 7.*pow(e,2))*
                                    (-(m_s22*cos(2.*h + 2.*l)) + m_c22*sin(2.*h + 2.*l)) -
                                3.*pow(1. + ci,2)*(16. - 40.*pow(e,2) + 13.*pow(e,4))*
                                    (-(m_s22*cos(2.*g + 2.*h + 2.*l)) + m_c22*sin(2.*g + 2.*h + 2.*l))))/32.;
    
    // dG/dt
    dx[4] = -(pow(m_al,2)*pow(n,2)*(2.*pow(-1. + ci,2)*pow(e,4)*
                                    (m_s22*cos(2.*g - 2.*h - 2.*l) + m_c22*sin(2.*g - 2.*h - 2.*l)) +
                                3.*pow(1. + ci,2)*(16. - 40.*pow(e,2) + 13.*pow(e,4))*
                                    (-(m_s22*cos(2.*g + 2.*h + 2.*l)) + m_c22*sin(2.*g + 2.*h + 2.*l))))/32.;
    
    // dH/dt
    dx[5] = (pow(m_al,2)*pow(n,2)*(2.*pow(-1. + ci,2)*pow(e,4)*
                                    (m_s22*cos(2.*g - 2.*h - 2.*l) + m_c22*sin(2.*g - 2.*h - 2.*l)) +
                                24.*(-1. + pow(ci,2))*pow(e,2)*(9. + 7.*pow(e,2))*
                                    (-(m_s22*cos(2.*h + 2.*l)) + m_c22*sin(2.*h + 2.*l)) -
                                3.*pow(1. + ci,2)*(16. - 40.*pow(e,2) + 13.*pow(e,4))*
                                    (-(m_s22*cos(2.*g + 2.*h + 2.*l)) + m_c22*sin(2.*g + 2.*h + 2.*l))))/32.;
    
}

///// ZONAL PERTURBATIONS /////

template<typename T>
void HEOSATDynamics<T>::zonal2order(const T& g, const T& a, const T& n, const T& e, const T& y, const T& ci, const T& si, DACE::AlgebraicVector<T>& x)
{
	//Using std and DACE functions
	using std::pow;   using DACE::pow;

	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
	// g = argument of perigee
	// a = semi-major axis
	// n = mean motion
	// e = eccentricity
	// y = sqrt(1-e^2)
	// ci = cos(inclination)
	// si = sin(inclination)

	T p, e2, e4, e6, c2, c4, c6, c8, s2, s4, y2, y3, y4;
	T dJw, dJl, dJh, dJg, swj, slj, shj, sgj, qt, et, iet, bt;
	T q[11][9], b[11][9], dq[11][9], db[11][9];
	int i, j, ii, m, djm;

	q[2][0] = 1;
	q[3][1] = 1;
	e2 = e*e;
	q[4][0] = 2 + 3 * e2;
	q[4][2] = 1;
	q[5][1] = 4 + 3 * e2;
	q[5][3] = 1;
	e4 = e2*e2;
	q[6][0] = 8 + 40 * e2 + 15 * e4;
	q[6][2] = 3 * (2 + e2);
	q[6][4] = 1;
	q[7][1] = 3 * (8 + 20 * e2 + 5 * e4);
	q[7][3] = 8 + 3 * e2;
	q[7][5] = 1;
	e6 = e4*e2;
	q[8][0] = 3 * (16 + 168 * e2 + 210 * e4 + 35 * e6);
	q[8][2] = 48 + 80 * e2 + 15 * e4;
	q[8][4] = 10 + 3 * e2;
	q[8][6] = 1;
	q[9][1] = 3 * (64 + 336 * e2 + 280 * e4 + 35 * e6);
	q[9][3] = 5 * (16 + 20 * e2 + 3 * e4);
	q[9][5] = 3 * (4 + e2);
	q[9][7] = 1;
	q[10][0] = 3 * (128 + 2304 * e2 + 6048 * e4 + 3360 * e6 + 315 * e6*e2);
	q[10][2] = 15 * (32 + 112 * e2 + 70 * e4 + 7 * e6);
	q[10][4] = 15 * (8 + 8 * e2 + e4);
	q[10][6] = 14 + 3 * e2;
	q[10][8] = 1;

	dq[2][0] = 0;
	dq[3][1] = 0;
	dq[4][0] = 6 * e;
	dq[4][2] = 0;
	dq[5][1] = 6 * e;
	dq[5][3] = 0;
	dq[6][0] = (80 + 60 * e2)*e;
	dq[6][2] = 6 * e;
	dq[6][4] = 0;
	dq[7][1] = 3 * e*(40 + 20 * e2);
	dq[7][3] = 6 * e;
	dq[7][5] = 0;
	dq[8][0] = 3 * e*(336 + 840 * e2 + 210 * e4);
	dq[8][2] = e*(160 + 60 * e2);
	dq[8][4] = 6 * e;
	dq[8][6] = 0;
	dq[9][1] = 3 * e*(672 + 1120 * e2 + 210 * e4);
	dq[9][3] = 5 * e*(40 + 12 * e2);
	dq[9][5] = 6 * e;
	dq[9][7] = 0;
	dq[10][0] = 3 * e*(4608 + 24192 * e2 + 20160 * e4 + 2520 * e6);
	dq[10][2] = 15 * e*(224 + 280 * e2 + 42 * e4);
	dq[10][4] = 15 * e*(16 + 4 * e2);
	dq[10][6] = 6 * e;
	dq[10][8] = 0;

	c2 = ci*ci;
	b[2][0] = (-1 + 3 * c2) / 4.;
	b[3][1] = (-3 * (-1 + 5 * c2)) / 8.;
	c4 = c2*c2;
	b[4][0] = (-3 * (3 - 30 * c2 + 35 * c4)) / 128.;
	b[4][2] = (-15 * (-1 + 7 * c2)) / 64.;
	b[5][1] = (15 * (1 - 14 * c2 + 21 * c4)) / 128.;
	b[5][3] = (35 * (-1 + 9 * c2)) / 256.;
	c6 = c4*c2;
	b[6][0] = (5 * (-5 + 105 * c2 - 315 * c4 + 231 * c6)) / 2048.;
	b[6][2] = (175 * (1 - 18 * c2 + 33 * c4)) / 2048.;
	b[6][4] = (315 * (-1 + 11 * c2)) / 4096.;
	b[7][1] = (-35 * (-5 + 135 * c2 - 495 * c4 + 429 * c6)) / 8192.;
	b[7][3] = (-315 * (3 - 66 * c2 + 143 * c4)) / 16384.;
	b[7][5] = (-693 * (-1 + 13 * c2)) / 16384.;
	c8 = c4*c4;
	b[8][0] = (-35 * (35 - 1260 * c2 + 6930 * c4 - 12012 * c6 + 6435 * c8)) / 786432.;
	b[8][2] = (-2205 * (-1 + 33 * c2 - 143 * c4 + 143 * c6)) / 131072.;
	b[8][4] = (-4851 * (1 - 26 * c2 + 65 * c4)) / 131072.;
	b[8][6] = (-3003 * (-1 + 15 * c2)) / 131072.;
	b[9][1] = (105 * (7 - 308 * c2 + 2002 * c4 - 4004 * c6 + 2431 * c8)) / 262144.;
	b[9][3] = (1617 * (-1 + 39 * c2 - 195 * c4 + 221 * c6)) / 131072.;
	b[9][5] = (3003 * (1 - 30 * c2 + 85 * c4)) / 131072.;
	b[9][7] = (6435 * (-1 + 17 * c2)) / 524288.;
	b[10][0] = (21 * (-63 + 3465 * c2 - 30030 * c4 + 90090 * c6 - 109395 * c8 + 46189 * c4*c6)) / 8.388608e6;
	b[10][2] = (693 * (7 - 364 * c2 + 2730 * c4 - 6188 * c6 + 4199 * c8)) / 2.097152e6;
	b[10][4] = (9009 * (-1 + 45 * c2 - 255 * c4 + 323 * c6)) / 1.048576e6;
	b[10][6] = (19305 * (3 - 102 * c2 + 323 * c4)) / 4.194304e6;
	b[10][8] = (109395 * (-1 + 19 * c2)) / 1.6777216e7;

	db[2][0] = (3 * ci) / 2.;
	db[3][1] = (-15 * ci) / 4.;
	db[4][0] = (-3 * ci*(-60 + 140 * c2)) / 128.;
	db[4][2] = (-105 * ci) / 32.;
	db[5][1] = (15 * ci*(-28 + 84 * c2)) / 128.;
	db[5][3] = (315 * ci) / 128.;
	db[6][0] = (5 * ci*(210 - 1260 * c2 + 1386 * c4)) / 2048.;
	db[6][2] = (175 * ci*(-36 + 132 * c2)) / 2048.;
	db[6][4] = (3465 * ci) / 2048.;
	db[7][1] = (-35 * ci*(270 - 1980 * c2 + 2574 * c4)) / 8192.;
	db[7][3] = (-315 * ci*(-132 + 572 * c2)) / 16384.;
	db[7][5] = (-9009 * ci) / 8192.;
	db[8][0] = (-35 * ci*(-2520 + 27720 * c2 - 72072 * c4 + 51480 * c6)) / 786432.;
	db[8][2] = (-2205 * ci*(66 - 572 * c2 + 858 * c4)) / 131072.;
	db[8][4] = (-4851 * ci*(-52 + 260 * c2)) / 131072.;
	db[8][6] = (-45045 * ci) / 65536.;
	db[9][1] = (105 * ci*(-616 + 8008 * c2 - 24024 * c4 + 19448 * c6)) / 262144.;
	db[9][3] = (1617 * ci*(78 - 780 * c2 + 1326 * c4)) / 131072.;
	db[9][5] = (3003 * ci*(-60 + 340 * c2)) / 131072.;
	db[9][7] = (109395 * ci) / 262144.;
	db[10][0] = (21 * ci*(6930 - 120120 * c2 + 540540 * c4 - 875160 * c6 + 461890 * c8)) / 8.388608e6;
	db[10][2] = (693 * ci*(-728 + 10920 * c2 - 37128 * c4 + 33592 * c6)) / 2.097152e6;
	db[10][4] = (9009 * ci*(90 - 1020 * c2 + 1938 * c4)) / 1.048576e6;
	db[10][6] = (19305 * ci*(-204 + 1292 * c2)) / 4.194304e6;
	db[10][8] = (2078505 * ci) / 8.388608e6;

	s2 = 1 - c2;
	y2 = y*y;
	dJw = 0.;
	dJl = 0.;
	dJh = 0.;
	dJg = 0.;
	for (i = 3; i <= 10; ++i){
		ii = (i / 2) - 1;
		m = i % 2;
		swj = 0.;
		slj = 0.;
		shj = 0.;
		sgj = 0.;
		for (j = 0; j <= ii; ++j){
			djm = (2 * j + m);
			qt = (y2 / e2)*(djm*q[i][djm] + e*dq[i][djm]);
			bt = pow(si*e, djm);
			if (m == 1){
				iet = -bt*cos(djm*g);
				et = -bt*sin(djm*g);
			}
			else{
				et = bt*cos(djm*g);
				iet = -bt*sin(djm*g);
			}
			swj = swj + djm*q[i][djm] * b[i][djm] * iet;
			slj = slj + (qt - 3 * q[i][djm])*b[i][djm] * et;
			bt = djm*(ci / s2)*b[i][djm] - db[i][djm];
			shj = shj + q[i][djm] * bt*et;
			sgj = sgj + q[i][djm] * b[i][djm] * et;
		}
		p = a*y2;
		et = m_C[i] * pow(m_al / p, i);
		dJw = dJw + et*swj;
		dJl = dJl + et*slj;
		dJh = dJh + et*shj;
		dJg = dJg + 2 * (i + 1)*et*sgj;
	}
	dJw = (n*a*a*y)*n*dJw;
	dJl = n*y*dJl;
	dJh = -n*dJh;
	dJg = -n*dJg - (1. / y)*dJl - ci*dJh;
	s4 = s2*s2;
	y3 = y2*y;
	et = 3 * pow(m_al / p, 4)*n*m_C[2] * m_C[2];
	bt = (1 + y)*(1 + y);

	x[3] = 0.;
	x[4] = -dJw + et*n*s2*e2*p*p*sin(2 * g)*(3.5 - (15 * s2) / 4. + ((4 - 5 * s2)*(1 + 2 * y)) / bt) / (8.*y3);
	x[5] = 0.;

	y4 = y2*y2;
	bt = cos(2 * g) / (64.*bt);

	x[0] = dJl + et*y*((8 * (15 + 8 * y - 5 * y2) + 8 * s2*(-30 - 24 * y + 5 * y2) + s4*(105 + 144 * y + 25 * y2)) / 128. - (5 * s2*(-21 - 42 * y + 34 * y2 + 62 * y3 + 15 * y4) - 2 * (-45 - 90 * y + 70 * y2 + 134 * y3 + 35 * y4))*s2*bt);
	x[1] = dJg + et*((8 * (55 + 24 * y - 7 * y2) + 4 * (-215 - 132 * y + 9 * y2)*s2 + 5 * (77 + 72 * y + 9 * y2)*s4) / 128. - (4 * (-15 - 30 * y + 8 * y2 + 30 * y3 + 7 * y4) - 2 * (-205 - 410 * y + 66 * y2 + 366 * y3 + 79 * y4)*s2 + 5 * (-77 - 154 * y + 22 * y2 + 134 * y3 + 27 * y4)*s4)*bt);
	x[2] = dJh + et*ci*((4 * (-10 - 6 * y + y2) + (35 + 36 * y + 5 * y2)*s2) / 32. + 4 * (-15 - 30 * y - 7 * y2 + 5 * (7 + 14 * y + 3 * y2)*s2)*e2*bt);
}

///// EQUATIONS OF MOTION /////

template<typename T>
DACE::AlgebraicVector<T> HEOSATDynamics<T>::difeq(const DACE::AlgebraicVector<T>& x, const double t, bool& eventFlag, const int thirdBodyFlag, const int srpFlag, const int dragFlag, const int tesseralFlag)
{
	/*
	Compute averaged time derivatives of Delaunay elements due to different perturbations (zonal, tesseral resonance, Sun, Moon, SRP and drag)
	- x: vector of Delaunay elements: l,g,h,L,G,H
	- t: scaled time
	- eventFlag: boolean: 'true' if reentry occurs, else 'false'
	- thirdBodyFlag: 0=no third-body perturbations, 1=Sun+Moon analytical ephemeris, 2=Sun+Moon CSPICE ephemeris
	- srpFlag: 0=no SRP perturbation, 1=with SRP
	- dragFlag: 0=no drag perturbation, 1=with drag
	- tesseralFlag: 0=no tesseral resonance, 1=with tesseral resonance
	*/
	//Using std and DACE functions
	using std::pow;   using DACE::pow;
	using std::sqrt;  using DACE::sqrt;

	DACE::AlgebraicVector<T> dx(6);

	const T ci = x[5] / x[4];           // H/G
	const T si = sqrt(1 - ci*ci);
	const T eta = x[4] / x[3];           // y = G/L
	const T eta3 = pow(eta, 3);
	const T eta4 = eta3*eta;
	const T e = sqrt(1. - eta*eta);
	const T a = x[3] * x[3];           // a = L**2/mu; mu = 1 en las unidades internas
	const T n = sqrt(1. / a) / a;        // mu = 1 en las unidades internas
	const T h = x[2];
	const T g = x[1];
	const T l = x[0];

	// l', g', h', L', G', H'
	DACE::AlgebraicVector<T> dfeq(6);

    // First-order J2 perturbation
    const T Raux = n*(m_al / a)*(m_al / a)*m_C[2];
    dfeq[0] = n + (-3. / 4.)*Raux*(2. - 3. * si*si) / eta3;
    dfeq[1] = (-3. / 4.)*Raux*(4. - 5. * si*si) / eta4;
    dfeq[2] = (3. / 2.)*Raux*ci / eta4;
    dfeq[3] = 0.;
    dfeq[4] = 0.;
    dfeq[5] = 0.;

    // Second-order zonal perturbations
    // TODO: set dx = 0
    zonal2order(g, a, n, e, eta, ci, si, dx);
    for (int i = 0; i <= 5; ++i)
    {
        dfeq[i] = dfeq[i] + dx[i];
    }

    // Tesseral resonance
    /*
    if (tesseralFlag == 1)
    {
        // Tesseral resonance (1:1) - 24-hour orbits
        // TODO: dx = 0
        tesre1t1(t, l, g, h, a, n, e, ci, dx);
        for (int i = 0; i <= 5; ++i)
        {
            dfeq[i] = dfeq[i] + dx[i];
        }
    }
    else if (tesseralFlag == 2)
    {
        // Tesseral resonance (2:1) - 12-hour orbits
        // TODO: dx = 0
        tesre2t1(t, l, g, h, a, n, e, eta, ci, dx);
        for (int i = 0; i <= 5; ++i)
        {
            dfeq[i] = dfeq[i] + dx[i];
        }
    }
    */
	
	// Sun & Moon positions
	if (thirdBodyFlag != 0)
	{
		if (thirdBodyFlag == 1)
		{
            // Analytical Sun and Moon ephemeris
			elsetsumo(ut2dj(t));
		}
		else
		{
            // CSPICE Sun and Moon ephemeris
			sunMoonSpice(ut2dj(t));
		}

		// Third-body perturbations
        // TODO: dx = 0
		sumord5(t, g, h, a, n, e, eta, ci, si, dx);
		for (int i = 0; i <= 5; ++i)
		{
			dfeq[i] = dfeq[i] + dx[i];
		}

        // TODO: dx = 0
		sumord6(t, g, h, a, n, e, eta, ci, si, dx);
		for (int i = 0; i <= 5; ++i)
		{
			dfeq[i] = dfeq[i] + dx[i];
		}

        // TODO: dx = 0
		sumord7(t, g, h, a, n, e, eta, ci, si, dx);
		for (int i = 0; i <= 5; ++i)
		{
			dfeq[i] = dfeq[i] + dx[i];
		}

        // TODO: dx = 0
		sumord8(t, g, h, a, n, e, eta, ci, si, dx);
		for (int i = 0; i <= 5; ++i)
		{
			dfeq[i] = dfeq[i] + dx[i];
		}

        // TODO: dx = 0
		sumord9(t, g, h, a, n, e, eta, ci, si, dx);
		for (int i = 0; i <= 5; ++i)
		{
			dfeq[i] = dfeq[i] + dx[i];
		}
	}
	
	// SRP perturbation
	if (srpFlag != 0)
	{
        // TODO: dx = 0
		srpaverage(t, g, h, a, n, e, eta, ci, si, dx);
		for (int i = 0; i <= 5; ++i)
		{
			dfeq[i] = dfeq[i] + dx[i];
		}
	}

	// Drag perturbation
	if (dragFlag != 0)
	{
		const double alt = DACE::cons(a*(1.0 - e) - m_al)*m_ul;
		if (alt <= 100.0)
		{
			fprintf(stderr, "Reentry: perigee height = %lf\n", alt);
			fprintf(stderr, "%lf %lf %lf \n", DACE::cons(a)*m_ul, DACE::cons(e), t*m_ut / (60. * 60. * 24. * 365.25));
			eventFlag = true;
		}
		else if ((alt > 100.0) && (alt <= 1000.0))
		{
			elsetsumo(ut2dj(t));
            // TODO: dx = 0
			drag(t, n, a, e, eta, ci, si, h, g, l, dx);

			for (int i = 0; i <= 5; ++i)
			{
				dfeq[i] = dfeq[i] + dx[i];
			}
		}
	}
	
	// Tesseral resonance
	double psi_dot_lim = 2.*M_PI/(5*spd_c())*m_ut;
	if (tesseralFlag == 1){
		double M_dot = DACE::cons(n);
		double m_we_ad = m_we*m_ut;

		if(std::abs(2.*M_dot-2.*m_we_ad)<psi_dot_lim){
			//1:1 resonance (24-hours)
			//std::cout<<"1:1 resonance included"<<std::endl;
			tesre1t1(t, l, g, h, a, n, e, ci, dx);	
			for (int i = 0; i <= 5; ++i){
				dfeq[i] = dfeq[i] + dx[i];
			}
		}
		else if(std::abs(M_dot-2.*m_we_ad)<psi_dot_lim){
			//2:1 resonance (12-hours)
			//std::cout<<"2:1 resonance included"<<std::endl;
			tesre2t1(t, l, g, h, a, n, e, eta, ci, dx);
			for (int i = 0; i <= 5; ++i){
				dfeq[i] = dfeq[i] + dx[i];
			}
		}
	}
	return dfeq;
}
#endif
