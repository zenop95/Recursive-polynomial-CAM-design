#pragma once

#include <dace/dace.h>
#include <cmath>

#ifdef __GNUC__   // GNU C++ adaptation
#include <float.h>
#else             // Standard C++ version
#include <limits>
#endif

namespace // Unnamed namespace
{
	// Machine accuracy
#ifdef __GNUC__   // GNU C++ adaptation
	const double eps_double = DBL_EPSILON;
#else             // Standard C++ version
	const double eps_double = std::numeric_limits<double>::epsilon();
#endif
}


DACE::DA modulusLocal(DACE::DA numer, double denom)
{
	double numerf = DACE::cons(numer);
	double modf = std::remainder(numerf, denom);
	return numer - numerf + modf;
}

double modulusLocal(double numer, double denom)
{
	return std::remainder(numer, denom);
}

template<typename T> T true2eccAnomaly(const T theta, const T e)
{
	//Using std and DACE functions
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	using std::sqrt;  using DACE::sqrt;
	
	return 2.0 * atan2_mod(sqrt(1. - e)*sin(theta / 2.), sqrt(1. + e) * cos(theta / 2.));
}

template<typename T> T ecc2trueAnomaly(const T E, const T e)
{
	//Using std and DACE functions
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	using std::sqrt;  using DACE::sqrt;
	
	return 2.0 * atan2_mod(sqrt(1. + e)*sin(E / 2.), sqrt(1. - e) * cos(E / 2.));
}

template<typename T> T mean2eccAnomaly(T M, const T e)
{
	
	//Using std and DACE functions
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;

	// Constants
	const int maxit = 15;
	const double eps = 100.0*eps_double;

	// Variables
	T E;

	// Starting value
	M = modulusLocal(M, 2.0*M_PI);
	if (DACE::cons(e) < 0.8)
	{
		E = M;
	}
	else
	{
		E = M_PI;
	}

	// Iteration
	int i = 0;
	T f = 1.0;
	while (std::abs(DACE::cons(f)) > eps)
	{
		f = E - e*sin(E) - M;
		E = E - f / (1.0 - e*cos(E));
		++i;
		if (i >= maxit) {
			std::cerr << " convergence problems in mean2eccAnomaly" << std::endl;
			break;
		}
	}

	return E;
}

template<typename T> T mean2trueAnomaly(const T M, const T e)
{
	T E = mean2eccAnomaly(M, e);

	return ecc2trueAnomaly(E, e);
}


template<typename T> T true2meanAnomaly(const T theta, const T e)
{
	//Using std and DACE functions
	using std::sin;   using DACE::sin;
	
	T E = true2eccAnomaly(theta, e);

	return E - e*sin(E);
}


template<typename T> void kep2pol(T kep[6], T *r, T *th, T *nu, T *RR, T *ZZ, T *NN, const double GM)
{
	//Using std and DACE functions
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	using std::sqrt;  using DACE::sqrt;
	
	T a, e, i, W, g, M, p, f;

	a = kep[0];
	e = kep[1];
	i = kep[2];
	W = kep[3];
	g = kep[4];
	M = kep[5];
	p = a*(1. - e*e);
	*ZZ = sqrt(GM*p); // DG: Angular momentum
	f = mean2trueAnomaly(M, e);
	*r = p / (1. + e*cos(f));
	*th = f + g;
	*nu = W;
	*RR = (*ZZ / p)*e*sin(f);
	*NN = *ZZ*cos(i); // DG: Angular momentum in z-direction
}

template<typename T> DACE::AlgebraicVector<T> kep2pol(const DACE::AlgebraicVector<T>& kepVec, const double GM)
{
	T kep[6];
	for (int i = 0; i < 6; i++)
	{
		kep[i] = kepVec[i];
	}

	T r, th, nu, RR, ZZ, NN;
	kep2pol(kep, &r, &th, &nu, &RR, &ZZ, &NN, GM);

	DACE::AlgebraicVector<T> polVec(6);
	polVec[0] = r;
	polVec[1] = th;
	polVec[2] = nu;
	polVec[3] = RR;
	polVec[4] = ZZ;
	polVec[5] = NN;

	return polVec;
}

template<typename T> void pol2kepT(const T r, const T th, const T nu, const T RR, const T ZZ, const T NN, T kepT[6], const double GM)
{
	//Using std and DACE functions
	using std::pow;   using DACE::pow;
	using std::acos;  using DACE::acos;
	using std::sqrt;  using DACE::sqrt;
	
	// Polar nodal to Keplerian orbital elements
	T p = ZZ*ZZ / GM;
	T i = acos(NN / ZZ); // inclination
	T a = -GM / (RR*RR + pow(ZZ / r, 2) - 2. * GM / r); // semi-major axis
	T e = sqrt(1. - p / a); // eccentricity
	T f = atan2_mod(p*RR / ZZ, -1. + p / r); // ATAN2(e*SIN(f), e*COS(f));
	T w = th - f; // argument of perigee
	T h = nu; // RAAN

	kepT[0] = a;
	kepT[1] = e;
	kepT[2] = i;
	kepT[3] = h;
	kepT[4] = w;
	kepT[5] = f;
}

template<typename T> DACE::AlgebraicVector<T> pol2kepT(const DACE::AlgebraicVector<T>& pol, const double GM)
{
	T kepT[6];
	pol2kepT(pol[0], pol[1], pol[2], pol[3], pol[4], pol[5], kepT, GM);

	DACE::AlgebraicVector<T> kepVec(6);
	for (int i = 0; i < 6; i++)
	{
		kepVec[i] = kepT[i];
	}

	return kepVec;
}

template<typename T> void pol2kep(const T r, const T th, const T nu, const T RR, const T ZZ, const T NN, T kep[6], const double GM)
{
	//Using std and DACE functions
	using std::tan;   using DACE::tan;
	using std::sqrt;  using DACE::sqrt;
	using std::atan;  using DACE::atan;
	
	pol2kepT(r, th, nu, RR, ZZ, NN, kep, GM);
	
	T e = kep[1];
	T u = 2. * atan(sqrt((1. - e) / (1. + e))*tan(kep[5] / 2.)); // Eccentric anomaly
	T l = u - e*sin(u); // Mean anomaly

	kep[5] = l;
}

template<typename T> DACE::AlgebraicVector<T> pol2kep(const DACE::AlgebraicVector<T>& pol, const double GM)
{
	T kep[6];
	pol2kep(pol[0], pol[1], pol[2], pol[3], pol[4], pol[5], kep, GM);

	DACE::AlgebraicVector<T> kepVec(6);
	for (int i = 0; i < 6; i++)
	{
		kepVec[i] = kep[i];
	}

	return kepVec;
}

template<typename T> DACE::AlgebraicVector<T> 
in2orb(const DACE::AlgebraicVector<T>& rrin, const T i, const T Om, const T l)
{
	
	//Using std and DACE functions
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
	DACE::AlgebraicVector<T> row1(3), row2(3), row3(3);
	DACE::AlgebraicVector<T> rrorb(3);

	row1[0] = (cos(l)*cos(Om) - sin(l)*cos(i)*sin(Om));
	row1[1] = (cos(l)*sin(Om) + sin(l)*cos(i)*cos(Om));
	row1[2] = (sin(l)*sin(i));
	row2[0] = (-sin(l)*cos(Om) - cos(l)*cos(i)*sin(Om));
	row2[1] = (-sin(l)*sin(Om) + cos(l)*cos(i)*cos(Om));
	row2[2] = (cos(l)*sin(i));
	row3[0] = (sin(i)*sin(Om));
	row3[1] = (-sin(i)*cos(Om));
	row3[2] = (cos(i));

	rrorb[0] = DACE::dot(row1, rrin);
	rrorb[1] = DACE::dot(row2, rrin);
	rrorb[2] = DACE::dot(row3, rrin);

	return rrorb;

}

template<typename T>
DACE::AlgebraicVector<T> in2orb(const DACE::AlgebraicVector<T>& in, const DACE::AlgebraicVector<T>& rr, const DACE::AlgebraicVector<T>& vv)
{
	
	/* Convert vector from inertial to orbital reference frame. */
	DACE::AlgebraicVector<T> hh = DACE::cross(rr,vv);
	DACE::AlgebraicVector<T> tt = DACE::cross(hh,rr);

	const T rrnorm = DACE::vnorm(rr);
	const T hhnorm = DACE::vnorm(hh);
	const T ttnorm = DACE::vnorm(tt);

	DACE::AlgebraicVector<T> rru = rr / rrnorm;
	DACE::AlgebraicVector<T> hhu = hh / hhnorm;
	DACE::AlgebraicVector<T> ttu = tt / ttnorm;

	DACE::AlgebraicMatrix<T> rotmat(3, 3);
	for(int i=0;i<rru.size();i++){
		rotmat.at(0,i)=rru[i];
		rotmat.at(1,i)=ttu[i];
		rotmat.at(2,i)=hhu[i];
	}

	DACE::AlgebraicVector<T> orb = rotmat*in;

	return orb;
}


template<typename T> DACE::AlgebraicVector<T> 
orb2in(const DACE::AlgebraicVector<T>& rrorb, const T i, const T Om, const T l)
{

	//Using std and DACE functions
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
	DACE::AlgebraicVector<T> row1(3), row2(3), row3(3);
	DACE::AlgebraicVector<T> rrin(3);

	row1[0] = (cos(l)*cos(Om) - sin(l)*cos(i)*sin(Om));
	row1[1] = (-sin(l)*cos(Om) - cos(l)*cos(i)*sin(Om));
	row1[2] = (sin(i)*sin(Om));
	row2[0] = (cos(l)*sin(Om) + sin(l)*cos(i)*cos(Om));
	row2[1] = (-sin(l)*sin(Om) + cos(l)*cos(i)*cos(Om));
	row2[2] = (-sin(i)*cos(Om));
	row3[0] = (sin(l)*sin(i));
	row3[1] = (cos(l)*sin(i));
	row3[2] = (cos(i));


	rrin[0] = DACE::dot(row1, rrorb);
	rrin[1] = DACE::dot(row2, rrorb);
	rrin[2] = DACE::dot(row3, rrorb);

	return rrin;
}


template<typename T> DACE::AlgebraicVector<T> 
kep2delaunay(const DACE::AlgebraicVector<T>& kep, const double mu)
{
	// Convert Keplerian elements to Delaunay elements (a,e,i,RAAN,omega,M) -> (M,omega,RAAN,L,G,H)
	
	//Using std and DACE functions
	using std::pow;   using DACE::pow;
	using std::cos;   using DACE::cos;
	using std::sqrt;  using DACE::sqrt;
	
	// Keplerian orbital elements
	T a = kep[0];
	T e = kep[1];
	T i = kep[2];
	T RAAN = kep[3];
	T omega = kep[4];
	T M = kep[5];

	// Delaunay elements
	DACE::AlgebraicVector<T> delaunay(6);
	delaunay[0] = M; // l
	delaunay[1] = omega; // g
	delaunay[2] = RAAN; // h
	delaunay[3] = sqrt(mu*a); // L
	delaunay[4] = sqrt(1. - pow(e, 2)) * delaunay[3]; // G
	delaunay[5] = cos(i)*delaunay[4]; // H

	return delaunay;
}


template<typename T> DACE::AlgebraicVector<T> delaunay2kep(const DACE::AlgebraicVector<T>& delaunay, const double mu)
{
	// Convert Delaunay elements to Keplerian elements (M,omega,RAAN,L,G,H) -> (a,e,i,RAAN,omega,M)
	
	//Using std and DACE functions
	using std::pow;   using DACE::pow;
	using std::acos;  using DACE::acos;
	using std::sqrt;  using DACE::sqrt;
	
	// Delaunay elements
	T l = delaunay[0]; // l
	T g = delaunay[1]; // g
	T h = delaunay[2]; // h
	T L = delaunay[3]; // L
	T G = delaunay[4]; // G
	T H = delaunay[5]; // H

	// Keplerian orbital elements
	DACE::AlgebraicVector<T> kep(6);
	kep[0] = pow(L, 2) / mu; // a
	kep[1] = sqrt(1. - pow((G / L), 2)); // e
	kep[2] = acos(H / G); // i
	kep[3] = h; // RAAN
	kep[4] = g; // omega
	kep[5] = l; // M

	return kep;
}


template <typename T> DACE::AlgebraicVector<T> kep2cart(const DACE::AlgebraicVector<T>& kep, const double mu = 398600.4415)
{
	/*member function to convert keplerian  classical element into Earth-Centred inertial reference frame element
	!< keplerian element kep = {a, e, i, RA, PA, TA}
	RA: rigth ascension of ascending node; PA: argument of periapsis; TA: true anomaly; i:orbital inclination
	!> return DACE::AlgebraicVector of ECI reference frame res = {x, y, z, dx, dy, dz}*/

	//const double mu = 398600.4415;
	
	//Using std and DACE functions
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	using std::sqrt;  using DACE::sqrt;
	
	T p = kep[0] * (1.0 - kep[1] * kep[1]);

	// position and velocity in perifocal refererence frame
	DACE::AlgebraicVector<T> rm(3), vm(3);
	rm[0] = p*cos(kep[5]) / (1.0 + kep[1] * cos(kep[5]));
	rm[1] = p*sin(kep[5]) / (1.0 + kep[1] * cos(kep[5]));
	rm[2] = 0.0;
	vm[0] = -1.0*sin(kep[5])*sqrt(mu / p);
	vm[1] = (kep[1] + cos(kep[5]))*sqrt(mu / p);
	vm[2] = 0.0;

	T cRA = cos(kep[3]);
	T sRA = sin(kep[3]);
	T cPA = cos(kep[4]);
	T sPA = sin(kep[4]);
	T ci  = cos(kep[2]);
	T si  = sin(kep[2]);

	T RR[3][3]; // rotational matrix from perifocal to eci reference frame
	RR[0][0] = cRA*cPA - sRA*ci*sPA;  RR[0][1] = -cRA*sPA - sRA*ci*cPA; RR[0][2] = sRA*si;
	RR[1][0] = sRA*cPA + cRA*ci*sPA;  RR[1][1] = -sRA*sPA + cRA*ci*cPA; RR[1][2] = -cRA*si;
	RR[2][0] = si*sPA;                RR[2][1] = si*cPA;				RR[2][2] = ci;

	DACE::AlgebraicVector<T> rr(3), vv(3);
	for (unsigned int i = 0; i<3; i++){
		rr[i] = 0.0;
		vv[i] = 0.0;
		for (unsigned int j = 0; j<3; j++){
			rr[i] = rr[i] + RR[i][j] * rm[j];
			vv[i] = vv[i] + RR[i][j] * vm[j];
		}
	}

	DACE::AlgebraicVector<T> res(6);
	res[0] = rr[0];
	res[1] = rr[1];
	res[2] = rr[2];
	res[3] = vv[0];
	res[4] = vv[1];
	res[5] = vv[2];

	return res;
}


template <typename T> DACE::AlgebraicVector<T> kepM2cart(const DACE::AlgebraicVector<T>& kepM, const double mu = 398600.4415)
{
	DACE::AlgebraicVector<T> kep = kepM;
	kep[5] = mean2trueAnomaly(kepM[5], kepM[1]);

	return kep2cart(kep, mu);
}


template<typename T> void kep2cart(const T kep[6], T car[6], const double /*GM*/)
{
	DACE::AlgebraicVector<T> CEO(6);
	for (int i = 0; i < 5; i++)
	{
		CEO[i] = kep[i];
	}
	CEO[5] = mean2trueAnomaly(kep[5], kep[1]);

	DACE::AlgebraicVector<T> CART = kep2cart(CEO);
	for (int i = 0; i < 6; i++)
	{
		car[i] = CART[i];
	}
}

template<typename T> void kepT2cart(const T kepT[6], T car[6], const double GM)
{
	DACE::AlgebraicVector<T> CEO(6);
	for (int i = 0; i < 6; i++)
	{
		CEO[i] = kepT[i];
	}

	DACE::AlgebraicVector<T> CART = kep2cart(CEO);
	for (int i = 0; i < 6; i++)
	{
		car[i] = CART[i];
	}
}

template<typename T> DACE::AlgebraicVector<T> cart2kep(const DACE::AlgebraicVector<T>& rv, const double mu = 398600.4415)
{
	
	//Using std and DACE functions
	using std::pow;   using DACE::pow;
	using std::sqrt;  using DACE::sqrt;
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	using std::asin;  using DACE::asin;
	using std::acos;  using DACE::acos;
	
	DACE::AlgebraicVector<T> kep(6);
	DACE::AlgebraicVector<T> rr(3), vv(3);
	
	for (int i = 0; i < 3; i++)
	{
		rr[i] = rv[i];
		vv[i] = rv[i + 3];
	}

	T r = DACE::vnorm(rr);
	T v = DACE::vnorm(vv);
	DACE::AlgebraicVector<T> h = DACE::cross(rr, vv);

	kep[0] = mu / (2.0 * (mu / r - pow(v, 2) / 2.0));

	T h1sqr = pow(h[0], 2);
	T h2sqr = pow(h[1], 2);

	T RAAN;
	if (DACE::cons(h1sqr + h2sqr) == 0.0)
	{
		RAAN = 0.0;
	}
	else
	{
		T sinOMEGA = h[0] / sqrt(h1sqr + h2sqr);
		T cosOMEGA = -1.0*h[1] / sqrt(h1sqr + h2sqr);
		if (DACE::cons(cosOMEGA) >= 0.0)
		{
			if (DACE::cons(sinOMEGA) >= 0.0)
			{
				RAAN = asin(h[0] / sqrt(h1sqr + h2sqr));
			}
			else
			{
				RAAN = 2.0 * M_PI + asin(h[0] / sqrt(h1sqr + h2sqr));
			}
		}
		else
		{
			if (DACE::cons(sinOMEGA) >= 0.0)
			{
				RAAN = acos(-1.0*h[1] / sqrt(h1sqr + h2sqr));
			}
			else
			{
				RAAN = 2.0 * M_PI - acos(-1.0*h[1] / sqrt(h1sqr + h2sqr));
			}
		}
	}

	//RAAN = real(RAAN);

	DACE::AlgebraicVector<T> ee = 1.0 / mu*(DACE::cross(vv, h)) - rr / DACE::vnorm(rr);
	T e = DACE::vnorm(ee);

	T i = acos(h[2] / DACE::vnorm(h));

	kep[1] = e;
	kep[2] = i;
	kep[3] = RAAN;

	T omega;
	T theta;
	if (DACE::cons(e) <= 1.0e-8 && DACE::cons(i) < 1.0e-8)
	{
		e = 0.0;
		omega = atan2_mod(rr[1], rr[0]);
		theta = 0.0;
		kep[4] = omega;
		kep[5] = theta;
		return kep;
	}

	if (DACE::cons(e) <= 1.0e-8 && DACE::cons(i) >= 1.0e-8)
	{
		omega = 0;
		DACE::AlgebraicVector<T> P(3), Q(3), W(3);
		P[0] = cos(omega)*cos(RAAN) - sin(omega)*cos(i)*sin(RAAN);
		P[1] = -sin(omega)*cos(RAAN) - cos(omega)*cos(i)*sin(RAAN);
		P[2] = sin(RAAN)*sin(i);
		Q[0] = cos(omega)*sin(RAAN) + sin(omega)*cos(i)*cos(RAAN);
		Q[1] = -sin(omega)*sin(RAAN) + cos(omega)*cos(i)*cos(RAAN);
		Q[2] = -cos(RAAN)*sin(i);
		W[0] = sin(omega)*sin(i);
		W[1] = cos(omega)*sin(i);
		W[2] = cos(i);
		DACE::AlgebraicVector<T> rrt = P*rr[0] + Q*rr[1] + W*rr[2];
		theta = atan2_mod(rrt[1], rrt[0]);
		kep[4] = omega;
		kep[5] = theta;
		return kep;
	}

	T dotRxE = DACE::dot(rr, ee);
	T RxE = DACE::vnorm(rr)*DACE::vnorm(ee);
	if (std::abs(DACE::cons(dotRxE)) > std::abs(DACE::cons(RxE)) && 
		std::abs(DACE::cons(dotRxE)) - std::abs(DACE::cons(RxE)) < 
		std::abs(eps_double*DACE::cons(dotRxE)))
	{
		dotRxE -= eps_double*dotRxE;
	}
	theta = acos(dotRxE / RxE);

	if (DACE::cons(DACE::dot(rr, vv)) < 0.0)
	{
		theta = 2.0 * M_PI - theta;
	}

	if (DACE::cons(i) <= 1.0e-8 && DACE::cons(e) >= 1.0e-8)
	{
		i = 0.0;
		omega = atan2_mod(ee[1], ee[0]);
		kep[4] = omega;
		kep[5] = theta;
		return kep;
	}

	T sino = rr[2] / r / sin(i);
	T coso = (rr[0] * cos(RAAN) + rr[1] * sin(RAAN)) / r;
	T argLat;
	if (DACE::cons(coso) >= 0.0)
	{
		if (DACE::cons(sino) >= 0.0)
		{
			argLat = asin(rr[2] / r / sin(i));
		}
		else
		{
			argLat = 2.0 * M_PI + asin(rr[2] / r / sin(i));
		}
	}
	else
	{
		if (DACE::cons(sino) >= 0.0)
		{
			argLat = acos((rr[0] * cos(RAAN) + rr[1] * sin(RAAN)) / r);
		}
		else
		{
			argLat = 2.0 * M_PI - acos((rr[0] * cos(RAAN) + rr[1] * sin(RAAN)) / r);
		}
	}
	//argLat = real(argLat);
	omega = argLat - theta;

	if (DACE::cons(omega) < 0.0)
	{
		omega = omega + 2.0 * M_PI;
	}
	//omega = real(omega);

	kep[4] = omega;
	kep[5] = theta;

	return kep;
}


template<typename T> DACE::AlgebraicVector<T> cart2kepM(const DACE::AlgebraicVector<T>& rv, const double mu = 398600.4415)
{
	DACE::AlgebraicVector<T> kep = cart2kep(rv, mu);
	kep[5] = true2meanAnomaly(kep[5], kep[1]);
	
	return kep;
}

template<typename T>
void mee2rrvv(const DACE::AlgebraicVector<T>& mee, DACE::AlgebraicVector<T>& rr, DACE::AlgebraicVector<T>& vv, const double mu)
{
	
	//Using std and DACE functions
	using std::sqrt;  using DACE::sqrt;
	
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
	/* Convert modified equinoctial elements to position and velocity. */
	T p = mee[0];
	T f = mee[1];
	T g = mee[2];
	T h = mee[3];
	T k = mee[4];
	T l = mee[5];

	T f2 = f*f;
	T g2 = g*g;
	T h2 = h*h;
	T k2 = k*k;

	T a2 = h2 - k2;
	T s2 = 1. + h2 + k2;

	T sinl = sin(l);
	T cosl = cos(l);

	T w = 1. + f*cosl + g*sinl;
	T r = p / w;

	T is2 = 1. / s2;
	T rs2 = r*is2;
	T hk = h*k;
	T mup = sqrt(mu / p);

	rr[0] = rs2*(cosl + a2*cosl + 2. * hk*sinl);
	rr[1] = rs2*(sinl - a2*sinl + 2. * hk*cosl);
	rr[2] = 2. * rs2*(h*sinl - k*cosl);

	vv[0] = -is2*mup*(sinl + a2*sinl - 2. * hk*cosl + g - 2. * f*hk + a2*g);
	vv[1] = -is2*mup*(-cosl + a2*cosl + 2. * hk*sinl - f + 2. * g*hk + a2*f);
	vv[2] = 2. * is2*mup*(h*cosl + k*sinl + f*h + g*k);
}

template<typename T> DACE::AlgebraicVector<T> mee2kep(const DACE::AlgebraicVector<T>& mee, const double mu = 398600.4415)
{
	/* convert modified equinoctial elements to classical orbit elements
	input
	mee(1) = semiparameter (kilometers)
	mee(2) = f equinoctial element
	mee(3) = g equinoctial element
	mee(4) = h equinoctial element
	mee(5) = k equinoctial element
	mee(6) = true longitude (radians)

	output
	coe(1) = semimajor axis (kilometers)
	coe(2) = eccentricity
	coe(3) = inclination (radians)
	coe(4) = argument of periapsis (radians)
	coe(5) = right ascension of ascending node (radians)
	coe(6) = true anomaly (radians)
	*/

	//Using std and DACE functions
	using std::sqrt;  using DACE::sqrt;

	using std::atan;  using DACE::atan;
	
	// unload modified equinoctial orbital elements
	T pmee = mee[0];
	T fmee = mee[1];
	T gmee = mee[2];
	T hmee = mee[3];
	T kmee = mee[4];
	T lmee = mee[5];

	// compute classical orbital elements
	T tani2s = sqrt(hmee * hmee + kmee * kmee);

	// orbital eccentricity
	T ecc = sqrt(fmee * fmee + gmee * gmee);

	// semimajor axis
	T sma = pmee / (1.0 - ecc * ecc);

	// orbital inclination
	T inc = 2.0 * atan(tani2s);

	// right ascension of ascending node
	T raan = atan2_mod(kmee, hmee);

	// argument of periapsis
	T atopo = atan2_mod(gmee, fmee);

	T argper = modulusLocal(atopo - raan, 2.0 * M_PI);

	// true anomaly
	T tanom = modulusLocal(lmee - atopo, 2.0 * M_PI);

	// load classical orbital element array
	DACE::AlgebraicVector<T> kep(6);
	kep[0] = sma;
	kep[1] = ecc;
	kep[2] = inc;
	kep[3] = raan;
	kep[4] = argper;
	kep[5] = tanom;

	return kep;
}

template<typename T>
DACE::AlgebraicVector<T> cyl2cart( const DACE::AlgebraicVector<T>& cylindricalState )
{
	// Convert cylindrical (r,theta,z,rdot,thetadot,zdot) to Cartesian state (x,y,z,xdot,ydot,zdot).
	
	//Using std and DACE functions
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
    // Create Cartesian state vector, initialized with zero entries.
    DACE::AlgebraicVector<T> cartesianState(6, 0.0);

    // Get azimuth angle, theta.
    const T theta = cylindricalState[1];
	
	// Compute Cartesian coordinates.
    cartesianState[0] = cylindricalState[0] * cos( theta );    // x-coordinate
    cartesianState[1] = cylindricalState[0] * sin( theta );    // y-coordinate
    cartesianState[2] = cylindricalState[2];                  // z-coordinate

	// Compute and set Cartesian velocities.
	cartesianState[3] = cylindricalState[3] * cos( theta ) - cylindricalState[0] * cylindricalState[4] * sin( theta );   // xdot
	cartesianState[4] = cylindricalState[3] * sin( theta ) + cylindricalState[0] * cylindricalState[4] * cos( theta );   // ydot
	cartesianState[5] = cylindricalState[5];      // zdot

    return cartesianState;
}

template<typename T>
DACE::AlgebraicVector<T> cart2cyl( const DACE::AlgebraicVector<T>& cartesianState )
{
	// Convert Cartesian state (x,y,z,xdot,ydot,zdot) to cylindrical (r,theta,z,rdot,thetadot,zdot).
	
	//Using std and DACE functions
	using std::pow;   using DACE::pow;
	using std::sqrt;  using DACE::sqrt;
	
    // Create Cartesian state vector, initialized with zero entries.
    DACE::AlgebraicVector<T> cylindricalState(6, 0.0);
	
	// Compute Cartesian coordinates.
    cylindricalState[0] = sqrt(pow(cartesianState[0], 2) + pow(cartesianState[1], 2)); // Radius
    cylindricalState[1] = atan2_mod(cartesianState[1], cartesianState[0]); // Azimuth angle, theta
    cylindricalState[2] = cartesianState[2]; // z-coordinate

	if ( cons(cylindricalState[0]) < std::numeric_limits< double >::epsilon( ) )
	{
		cylindricalState[3] = sqrt(pow(cartesianState[3], 2) + pow(cartesianState[4], 2)); // rdot
		cylindricalState[4] = 0.0; // thetadot
	}
	else
	{
		cylindricalState[3] = ( cartesianState[0]*cartesianState[3] + cartesianState[1]*cartesianState[4] ) / cylindricalState[0]; // rdot
		cylindricalState[4] = ( cartesianState[0]*cartesianState[4] - cartesianState[1]*cartesianState[3] ) / pow(cylindricalState[0], 2); // thetadot
	}
	cylindricalState[5] = cartesianState[5]; // zdot

    return cylindricalState;
}


template<typename T> DACE::AlgebraicMatrix<T> kep2cartPartialDerivatives(const DACE::AlgebraicVector<T>& kep, const double mu = 398600.4415)
{
	//   Computes the partial derivatives of the state vector with respect to the orbital elements
	//
	//   References:
	//     - Satellite Orbits : Models, Methods and Applications - Montenbruck, Oliver and Gill, Eberhard, 2005. Sections 7.1.2 and 7.1.3
	//     - On the Matrizant of the Two-Body Problem, R.A. Broucke, 1970, Astronomy and Astrophysics
	//
	//   Created by : David Gondelach (University of Southampton), June 2017

	//Using std and DACE functions
	using std::pow;   using DACE::pow;
	using std::sqrt;  using DACE::sqrt;
	
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;

	
	// Variables
	T  a, e, i, RAAN, omega, M, n;
	T  E, cosE, sinE, fac, r, v, x, y, vx, vy;
	DACE::AlgebraicVector<T>  P(3), Q(3), W(3);
	DACE::AlgebraicVector<T>  dPdi(3), dPdO(3), dPdo(3), dQdi(3), dQdO(3), dQdo(3);
	DACE::AlgebraicVector<T>  drda(3), drde(3), drdi(3), drdO(3), drdo(3), drdM(3);
	DACE::AlgebraicVector<T>  dvda(3), dvde(3), dvdi(3), dvdO(3), dvdo(3), dvdM(3);
	DACE::AlgebraicMatrix<T>  dCartdKep(6, 6);

	// Keplerian elements
	a = kep[0];  RAAN = kep[3];
	e = kep[1];  omega = kep[4];
	i = kep[2];  M = kep[5];

	n = sqrt(mu / (a*a*a));

	// Eccentric anomaly
	E = mean2eccAnomaly(M, e);

	// Perifocal coordinates
	cosE = cos(E);
	sinE = sin(E);
	fac = sqrt((1.0 - e)*(1.0 + e));

	r = a*(1.0 - e*cosE);  // Distance
	v = sqrt(mu*a) / r;    // Velocity

	x = a*(cosE - e);
	y = a*fac*sinE;
	vx = -v*sinE;
	vy = v*fac*cosE;

	// Transformation to reference system (Gaussian vectors) and partials
	// PQW = R_z(-RAAN) * R_x(-i) * R_z(-omega);
	P[0] = +cos(omega)*cos(RAAN) - sin(omega)*cos(i)*sin(RAAN);
	P[1] = +cos(omega)*sin(RAAN) + sin(omega)*cos(i)*cos(RAAN);
	P[2] = +sin(omega)*sin(i);
	Q[0] = -sin(omega)*cos(RAAN) - cos(omega)*cos(i)*sin(RAAN);
	Q[1] = -sin(omega)*sin(RAAN) + cos(omega)*cos(i)*cos(RAAN);
	Q[2] = +cos(omega)*sin(i);
	W[0] = +sin(RAAN)*sin(i);
	W[1] = -cos(RAAN)*sin(i);
	W[2] = +cos(i);

	dPdi = sin(omega)*W;
	dPdO[0] = -P[1];
	dPdO[1] = P[0];
	dPdO[2] = 0.0;
	dPdo = Q;

	dQdi = cos(omega)*W;
	dQdO[0] = -Q[1];
	dQdO[1] = Q[0];
	dQdO[2] = 0.0;
	dQdo = -P;

	// Partials w.r.t. semimajor axis, eccentricity and mean anomaly
	drda = (x / a)*P + (y / a)*Q;
	dvda = (-vx / (2. * a))*P + (-vy / (2. * a))*Q;
	drde = (-a - pow(y / fac, 2) / r)*P + (x*y / (r*fac*fac))*Q;
	dvde = (vx*(2. * a*x + e*pow(y / fac, 2)) / (r*r))*P
		   + ((n / fac)*pow(a / r, 2)*(x*x / r - pow(y / fac, 2) / a))*Q;
	drdM = (vx*P + vy*Q) / n;
	dvdM = (-n*pow(a / r, 3))*(x*P + y*Q);

	// Partials w.r.t. inclination, node and argument of pericenter
	drdi = x*dPdi + y*dQdi;
	dvdi = vx*dPdi + vy*dQdi;
	drdO = x*dPdO + y*dQdO;
	dvdO = vx*dPdO + vy*dQdO;
	drdo = x*dPdo + y*dQdo;
	dvdo = vx*dPdo + vy*dQdo;

	// Combined partial derivative matrix of state with respect to epoch elements
	for (int k = 0; k<3; k++)
	{
		dCartdKep.at(k, 0) = drda[k];
		dCartdKep.at(k + 3, 0) = dvda[k];
		dCartdKep.at(k, 1) = drde[k];
		dCartdKep.at(k + 3, 1) = dvde[k];
		dCartdKep.at(k, 2) = drdi[k];
		dCartdKep.at(k + 3, 2) = dvdi[k];
		dCartdKep.at(k, 3) = drdO[k];
		dCartdKep.at(k + 3, 3) = dvdO[k];
		dCartdKep.at(k, 4) = drdo[k];
		dCartdKep.at(k + 3, 4) = dvdo[k];
		dCartdKep.at(k, 5) = drdM[k];
		dCartdKep.at(k + 3, 5) = dvdM[k];
	}

	return dCartdKep;
}

template<typename T> DACE::AlgebraicMatrix<T> cart2kepPartialDerivatives(const DACE::AlgebraicVector<T>& kep, const double mu = 398600.4415)
{
	//   Computes the partial derivatives of the orbital elements with respect to the state vector
	//   References:
	//     - Satellite Orbits : Models, Methods and Applications - Montenbruck, Oliver and Gill, Eberhard, 2005. Sections 7.1.2 and 7.1.3
	//     - On the Matrizant of the Two-Body Problem, R.A. Broucke, 1970, Astronomy and Astrophysics
	//
	//   Created by : David Gondelach (University of Southampton), June 2017

	//Using std and DACE functions
	using std::sqrt;  using DACE::sqrt;

	using std::tan; using DACE::tan;
	
	// Variables
	T  a, e, i, n, sqe2, naa;
	T  P_aM, P_eM, P_eo, P_io, P_iO;
	DACE::AlgebraicMatrix<T> dCartdKep(6, 6), dKepdCart(6, 6);

	// Orbital elements
	a = kep[0];  e = kep[1];  i = kep[2];

	n = sqrt(mu / (a*a*a));

	// State vector partials w.r.t epoch elements
	dCartdKep = kep2cartPartialDerivatives(kep, mu);

	// Poisson brackets
	sqe2 = sqrt((1.0 - e)*(1.0 + e));
	naa = n*a*a;

	P_aM = -2.0 / (n*a);                   // P(a,M)     = -P(M,a)
	P_eM = -(1.0 - e)*(1.0 + e) / (naa*e); // P(e,M)     = -P(M,e)
	P_eo = +sqe2 / (naa*e);                // P(e,omega) = -P(omega,e)
	P_io = -1.0 / (naa*sqe2*tan(i));       // P(i,omega) = -P(omega,i)
	P_iO = +1.0 / (naa*sqe2*sin(i));       // P(i,Omega) = -P(Omega,i)

	// Partials of epoch elements w.r.t. epoch state
	for (int k = 0; k<3; k++)
	{

		dKepdCart.at(0, k) = +P_aM*dCartdKep.at(k + 3, 5);
		dKepdCart.at(0, k + 3) = -P_aM*dCartdKep.at(k, 5);

		dKepdCart.at(1, k) = +P_eo*dCartdKep.at(k + 3, 4) + P_eM*dCartdKep.at(k + 3, 5);
		dKepdCart.at(1, k + 3) = -P_eo*dCartdKep.at(k, 4) - P_eM*dCartdKep.at(k, 5);

		dKepdCart.at(2, k) = +P_iO*dCartdKep.at(k + 3, 3) + P_io*dCartdKep.at(k + 3, 4);
		dKepdCart.at(2, k + 3) = -P_iO*dCartdKep.at(k, 3) - P_io*dCartdKep.at(k, 4);

		dKepdCart.at(3, k) = -P_iO*dCartdKep.at(k + 3, 2);
		dKepdCart.at(3, k + 3) = +P_iO*dCartdKep.at(k, 2);

		dKepdCart.at(4, k) = -P_eo*dCartdKep.at(k + 3, 1) - P_io*dCartdKep.at(k + 3, 2);
		dKepdCart.at(4, k + 3) = +P_eo*dCartdKep.at(k, 1) + P_io*dCartdKep.at(k, 2);

		dKepdCart.at(5, k) = -P_aM*dCartdKep.at(k + 3, 0) - P_eM*dCartdKep.at(k + 3, 1);
		dKepdCart.at(5, k + 3) = +P_aM*dCartdKep.at(k, 0) + P_eM*dCartdKep.at(k, 1);

	};

	return dKepdCart;
}


template<typename T> DACE::AlgebraicVector<T> kep2delaunayTimeDerivatives(const DACE::AlgebraicVector<T>& kep, const DACE::AlgebraicVector<T>& dkepdt, const double mu = 398600.4415)
{
	//   Convert Keplerian elements time derivatives to Delaunay elements time derivatives: d(a,e,i,RAAN,omega,M)/dt -> d(M,omega,RAAN,L,G,H)/dt
	//
	//   Created by : David Gondelach (University of Southampton), June 2017
	
	//Using std and DACE functions
	using std::sqrt;  using DACE::sqrt;
	
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;

	
	T a = kep[0];
	T e = kep[1];
	T i = kep[2];
	T y = sqrt(1. - e*e);
	T n = sqrt(mu / a) / a;

	DACE::AlgebraicVector<T> dDelaunaydt(6);
	dDelaunaydt[0] = dkepdt[5]; // dM/dt
	dDelaunaydt[1] = dkepdt[4]; // domega/dt
	dDelaunaydt[2] = dkepdt[3]; // dRAAN/dt
	dDelaunaydt[3] = 0.5 * dkepdt[0] * n*a; // dL/dt
	dDelaunaydt[4] = y*dDelaunaydt[3] - n*a*a*(e / y)*dkepdt[1]; // dG/dt
	dDelaunaydt[5] = cos(i)*dDelaunaydt[4] - n*a*a*y*sin(i)*dkepdt[2]; // dH/dt

	return dDelaunaydt;
}


template<typename T> DACE::AlgebraicVector<T> delaunay2kepTimeDerivatives(const DACE::AlgebraicVector<T>& kep, const DACE::AlgebraicVector<T>& dDelaunaydt, const double mu = 398600.4415)
{
	//   Convert Delaunay elements time derivatives to Keplerian elements time derivatives: d(M,omega,RAAN,L,G,H)/dt -> d(a,e,i,RAAN,omega,M)/dt
	//
	//   Created by : David Gondelach (University of Southampton), June 2017
	
	//Using std and DACE functions
	using std::sqrt;  using DACE::sqrt;
	
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
	T a = kep[0];
	T e = kep[1];
	T i = kep[2];
	T y = sqrt(1. - e*e);
	T n = sqrt(mu / a) / a;

	DACE::AlgebraicVector<T> dkepdt(6);
	dkepdt[0] = 2.0/(n*a)*dDelaunaydt[3]; // da/dt
	dkepdt[1] = y / (e*n*a*a)*(dDelaunaydt[3] * y - dDelaunaydt[4]); // de/dt
	dkepdt[2] = 1.0 / (n*a*a*y*sin(i))*(dDelaunaydt[4] * cos(i) - dDelaunaydt[5]); // di/dt
	dkepdt[3] = dDelaunaydt[2]; // dRAAN/dt
	dkepdt[4] = dDelaunaydt[1]; // domega/dt
	dkepdt[5] = dDelaunaydt[0]; // dM/dt

	return dkepdt;
}


template<typename T> DACE::AlgebraicVector<T> kep2cartTimeDerivatives(const DACE::AlgebraicVector<T>& kep, const DACE::AlgebraicVector<T>& dkepdt, const double mu = 398600.4415)
{
	//   Convert Keplerian elements time derivatives to Cartesian state time derivatives: d(a,e,i,RAAN,omega,M)/dt -> d(x,y,z,vx,vy,vz)/dt
	//   dy/dt = dy/da*da/dt
	//
	//   Created by : David Gondelach (University of Southampton), June 2017

	DACE::AlgebraicMatrix<T> dcartdkep = kep2cartPartialDerivatives(kep, mu);
	DACE::AlgebraicVector<T> dcartdt = dcartdkep*dkepdt;

	return dcartdt;
}

template<typename T> DACE::AlgebraicVector<T> cart2kepTimeDerivatives(const DACE::AlgebraicVector<T>& kep, const DACE::AlgebraicVector<T>& dcartdt, const double mu = 398600.4415)
{
	//   Convert Keplerian elements time derivatives to Cartesian state time derivatives: d(a,e,i,RAAN,omega,M)/dt -> d(x,y,z,vx,vy,vz)/dt
	//   dy/dt = dy/da*da/dt
	//
	//   Created by : David Gondelach (University of Southampton), June 2017

	DACE::AlgebraicMatrix<T> dkepdcart = cart2kepPartialDerivatives(kep, mu);
	DACE::AlgebraicVector<T> dkepdt = dkepdcart*dcartdt;

	return dkepdt;
}
