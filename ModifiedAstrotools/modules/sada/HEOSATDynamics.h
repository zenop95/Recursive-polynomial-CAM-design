//  HEOSAT propagator
//
//  Original code in C by: J.F. San-Juan, M. Lara, D. Hautesserres
//
//  C++ version by: D.J. Gondelach (University of Southampton) 2017
//  Added by D.J. Gondelach:
//  - Enable computations in both double-precision and DA variables
//  - Sun and Moon ephemeris from CSPICE
//  - Drag computation with Differential Algebra
//
//  References:
//  - M. Lara, J.F. San-Juan, D. Hautesserres, HEOSAT: a mean elements orbit propagator program for highly elliptical orbits, CEAS Space Journal, 2017. https://doi.org/10.1007/s12567-017-0152-x
//  - D.J. Gondelach, Orbit prediction and analysis for space situational awareness, Doctoral thesis, University of Surrey, Guildford, UK, 2019. http://doi.org/10.15126/thesis.00850116

#ifndef sada_HEOSATDYNAMICS_H
#define sada_HEOSATDYNAMICS_H

#pragma once

// DACE
#include <dace/dace.h>

// Dynamics
#include <dynorb/DynClass.h>


template<typename T>
class HEOSATDynamics : public dynamics<T>
{
public:
	HEOSATDynamics()
		: m_thirdBody(1), m_srp(1), m_drag(1), m_tesseral(1), m_eventflag(false) {}

	HEOSATDynamics(const int thirdBody, const int srp, const int drag, const int tesseral = 1)
		: m_thirdBody(thirdBody), m_srp(srp), m_drag(drag), m_tesseral(tesseral), m_eventflag(false) {}

	DACE::AlgebraicVector<T> evaluate(const DACE::AlgebraicVector<T>& x, const double t)
	{
		return difeq(x, t, m_eventflag, m_thirdBody, m_srp, m_drag, m_tesseral);
	}

	bool checkEvent()
	{
		return m_eventflag;
	}

public:
	double initheosat(const double aa, const double dj, const double tini, 
	const double tfin, double *t, double *tf,
	const T i_Bfactor = 0.0, const T i_SRPC = 0.0);

private:
	///// GENERAL /////
	double diajul(int giorno, int mese, int anno, double ora);
	double ut2dj(double t);
	double keplereq(double M, double e, int j);
	double tocent(double tt);
	double tsmut(double jd);
	double tsmut2(double jd);
	double tsmut3(double jd);

	///// SUN AND MOON COEFFICIENTS /////
	void coefsolun5(double *css, double *ccs, double *ks, double *csm, double *ccm, double *km);
	void coefsun6(double *cm, double *sm, double *km);
	void coeflun6(double *cm, double *sm, double *km);
	void coeflun7(double *cm, double *sm, double *km);
	void coeflun8(double *cm, double *sm, double *km);
	void coeflun9(double *cm, double *sm, double *km);

	///// DRAG /////
	T harrispriester(const T& h, const T *ca);
	void aeiwhm0new(const T& n, const T& a, const T& e, const T& y, const T& ci, const T& si, const T& g, const T& r,
		const T& th, const T& p, const double we, const T nd[3], T *a0, T *e0, T *i0, T *W0, T *g0, T *M0);
	void atmdesnew(T *kep, const T& si, const T& ci, T *r, T *th, T nd[3]);
	DACE::AlgebraicVector<T> computeDragMeanElementRatesTrueAnomaly(const T& n, const T& a, const T& e, const T& y, const T& ci, const T& si, const T& W, const T& g);
	void drag(const double t, const T& n, const T& a, const T& e, const T& y, const T& ci, const T& si, const T& W, const T& g, const T& l, DACE::AlgebraicVector<T>& dx);

	///// ECCENTRICITY AND INCLINATION FUNCTIONS /////
	void eccincpol5(const T& e, const T& y, const T& c, const T& s);
	void eccincpol6(const T& e, const T& y, const T& c, const T& s);
	void eccincpol7(const T& e, const T& y, const T& c, const T& s);
	void eccincpol8(const T& e, const T& y, const T& c, const T& s);
	void eccincpol9(const T& e, const T& y, const T& c, const T& s);

	///// LUNISOLAR PERTURBATIONS /////
	void sumord5(const double t, const T& g, const T& h, const T& a, const T& n, const T& e, const T& y, const T& ci, const T& si, DACE::AlgebraicVector<T>& dx);
	void sumord6(const double t, const T& g, const T& h, const T& a, const T& n, const T& e, const T& y, const T& ci, const T& si, DACE::AlgebraicVector<T>& dx);
	void sumord7(const double t, const T& g, const T& h, const T& a, const T& n, const T& e, const T& y, const T& ci, const T& si, DACE::AlgebraicVector<T>& dx);
	void sumord8(const double t, const T& g, const T& h, const T& a, const T& n, const T& e, const T& y, const T& ci, const T& si, DACE::AlgebraicVector<T>& dx);
	void sumord9(const double t, const T& g, const T& h, const T& a, const T& n, const T& e, const T& y, const T& ci, const T& si, DACE::AlgebraicVector<T>& dx);

	///// SUN AND MOON POSITIONS /////
	void meeussun(const double t, double *sz, double *Rs, double *Ms);
	void meeusmoon(const double t, const double Ms, const int j, double *le, double *be, double *rm);
	// Compute Sun and Moon positions
	void elsetsumo(const double jd);
	// Compute Sun and Moon positions using SPICE
	void sunMoonSpice(const double jd);

	///// SOLAR RADIATION PRESSURE /////
	T accelsrp(const double ms);
	void srpaverage(const double t, const T& g, const T& h, const T& a, const T& n, const T& e, const T& y, const T& ci, const T& si, DACE::AlgebraicVector<T>& x);

	///// TESSERAL PERTURBATIONS /////
	void tesre2t1(const double t, const T& l, const T& g, const T& hh, const T& a, const T& n, const T& e, const T& y, const T& ci, DACE::AlgebraicVector<T>& dx);
    void tesre1t1(const double t, const T& l, const T& g, const T& hh, const T& a, const T& n, const T& e, const T& ci, DACE::AlgebraicVector<T>& dx);

	///// ZONAL PERTURBATIONS /////
	void zonal2order(const T& g, const T& a, const T& n, const T& e, const T& y, const T& ci, const T& si, DACE::AlgebraicVector<T>& x);

	///// EQUATIONS OF MOTION /////
	DACE::AlgebraicVector<T> difeq(const DACE::AlgebraicVector<T>& x, const double t, bool& eventFlag, const int thirdBodyFlag = 1, const int srpFlag = 1, const int dragFlag = 1, const int tesseralFlag = 1);


private:
	bool m_eventflag;
	
	// SETTINGS //
	const int m_thirdBody;
	const int m_srp;
	const int m_drag;
	const int m_tesseral;

	// ORBIT CONSTANTS //
	double m_we, m_GM, m_al, m_mu, m_Req, m_C[11], m_fl, m_c22, m_s22;
	double m_ul, m_ut, m_utd, m_dj0, m_t0;
	double m_ns, m_nm, m_bet, m_am, m_as;
    
    // OBJECT PARAMETERS //
	T m_Bfactor, m_SRPC;

	// SUN AND MOON POSITIONS //
	double m_Rs, m_xs, m_ys, m_zs, m_rm, m_xl, m_yl, m_zl;

	// ECCENTRICITY AND INCLINATION FUNCTIONS //
	T m_ep5[6][2], m_ip5[6][2][5];
	T m_ep6[6][2], m_ip6[6][2][7];
	T m_ep7[6][3], m_ip7[6][3][9];
	T m_ep8[6][3], m_ip8[6][3][11];
	T m_ep9[6][4], m_ip9[6][4][13];

	// Harris-Priester altitude bands
	static const double hi[50];
	// Harris-Priester lower densities
	static const double pm[50];
	// Harris-Priester upper densities
	static const double qM[50];
	// Lunar longitude & radius terms
	static const int lr[60][6];
	// Lunar latitude terms
	static const int ab[60][5];

}; // END HEOSATDYNAMICS CLASS

///// CONST STATIC MEMBER VARIABLES /////
// Harris-Priester altitude bands
template<typename T>
const double HEOSATDynamics<T>::hi[50] = { 100, 120, 130, 140, 150, 160, 170, 180, 190,
    200, 210, 220, 230, 240, 250, 260, 270, 280, 290,
    300, 320, 340, 360, 380, 400, 420, 440, 460, 480,
    500, 520, 540, 560, 580, 600, 620, 640, 660, 680,
    700, 720, 740, 760, 780, 800, 840, 880, 920, 960, 1000 };
// Harris-Priester lower densities
template<typename T>
const double HEOSATDynamics<T>::pm[50] = { 497400.0, 24900.0, 8377.0, 3899.0, 2122.0, 1263.0, 800.8,
    528.3, 361.7, 255.7, 183.9, 134.1, 99.49, 74.88, 57.09,
    44.03, 34.30, 26.97, 21.39, 17.08, 10.99, 7.214, 4.824,
    3.274, 2.249, 1.558, 1.091, 0.7701, 0.5474, 0.3916, 0.2819,
    0.2042, 0.1488, 0.1092, 0.08070, 0.06012, 0.04519, 0.03430,
    0.02632, 0.02043, 0.01607, 0.01281, 0.01036, 0.008496,
    0.007069, 0.004680, 0.003200, 0.002210, 0.001560, 0.001150 };
// Harris-Priester upper densities
template<typename T>
const double HEOSATDynamics<T>::qM[50] = { 497400.0, 24900.0, 8710.0, 4059.0, 2215.0, 1344.0, 875.8,
    601.0, 429.7, 316.2, 239.6, 185.3, 145.5, 115.7, 93.08,
    75.55, 61.82, 50.95, 42.26, 35.26, 25.11, 18.19, 13.37,
    9.955, 7.492, 5.684, 4.355, 3.362, 2.612, 2.042, 1.605,
    1.267, 1.005, 0.7997, 0.6390, 0.5123, 0.4121, 0.3325,
    0.2691, 0.2185, 0.1779, 0.1452, 0.1190, 0.09776, 0.08059,
    0.05741, 0.04210, 0.03130, 0.02360, 0.01810 };

// Lunar longitude & radius terms
template<typename T>
const int HEOSATDynamics<T>::lr[60][6] = {
    { 0, 0, 1, 0, 6288774, -20905355 },
    { 2, 0, -1, 0, 1274027, -3699111 },
    { 2, 0, 0, 0, 658314, -2955968 },
    { 0, 0, 2, 0, 213618, -569925 },
    { 0, 1, 0, 0, -185116, 48888 },
    { 0, 0, 0, 2, -114332, -3149 },
    { 2, 0, -2, 0, 58793, 246158 },
    { 2, -1, -1, 0, 57066, -152138 },
    { 2, 0, 1, 0, 53322, -170733 },
    { 2, -1, 0, 0, 45758, -204586 },
    { 0, 1, -1, 0, -40923, -129620 },
    { 1, 0, 0, 0, -34720, 108743 },
    { 0, 1, 1, 0, -30383, 104755 },
    { 2, 0, 0, -2, 15327, 10321 },
    { 0, 0, 1, 2, -12528, 0 },
    { 0, 0, 1, -2, 10980, 79661 },
    { 4, 0, -1, 0, 10675, -34782 },
    { 0, 0, 3, 0, 10034, -23210 },
    { 4, 0, -2, 0, 8548, -21636 },
    { 2, 1, -1, 0, -7888, 24208 },
    { 2, 1, 0, 0, -6766, 30824 },
    { 1, 0, -1, 0, -5163, -8379 },
    { 1, 1, 0, 0, 4987, -16675 },
    { 2, -1, 1, 0, 4036, -12831 },
    { 2, 0, 2, 0, 3994, -10445 },
    { 4, 0, 0, 0, 3861, -11650 },
    { 2, 0, -3, 0, 3665, 14403 },
    { 0, 1, -2, 0, -2689, -7003 },
    { 2, 0, -1, 2, -2602, 0 },
    { 2, -1, -2, 0, 2390, 10056 },
    { 1, 0, 1, 0, -2348, 6322 },
    { 2, -2, 0, 0, 2236, -9884 },
    { 0, 1, 2, 0, -2120, 5751 },
    { 0, 2, 0, 0, -2069, 0 },
    { 2, -2, -1, 0, 2048, -4950 },
    { 2, 0, 1, -2, -1773, 4130 },
    { 2, 0, 0, 2, -1595, 0 },
    { 4, -1, -1, 0, 1215, -3958 },
    { 0, 0, 2, 2, -1110, 0 },
    { 3, 0, -1, 0, -892, 3258 },
    { 2, 1, 1, 0, -810, 2616 },
    { 4, -1, -2, 0, 759, -1897 },
    { 0, 2, -1, 0, -713, -2117 },
    { 2, 2, -1, 0, -700, 2354 },
    { 2, 1, -2, 0, 691, 0 },
    { 2, -1, 0, -2, 596, 0 },
    { 4, 0, 1, 0, 549, -1423 },
    { 0, 0, 4, 0, 537, -1117 },
    { 4, -1, 0, 0, 520, -1571 },
    { 1, 0, -2, 0, -487, -1739 },
    { 2, 1, 0, -2, -399, 0 },
    { 0, 0, 2, -2, -381, -4421 },
    { 1, 1, 1, 0, 351, 0 },
    { 3, 0, -2, 0, -340, 0 },
    { 4, 0, -3, 0, 330, 0 },
    { 2, -1, 2, 0, 327, 0 },
    { 0, 2, 1, 0, -323, 1165 },
    { 1, 1, -1, 0, 299, 0 },
    { 2, 0, 3, 0, 294, 0 },
    { 2, 0, -1, -2, 0, 8752 }
};

// Lunar latitude terms
template<typename T>
const int HEOSATDynamics<T>::ab[60][5] = {
    { 0, 0, 0, 1, 5128122 },
    { 0, 0, 1, 1, 280602 },
    { 0, 0, 1, -1, 277693 },
    { 2, 0, 0, -1, 173237 },
    { 2, 0, -1, 1, 55413 },
    { 2, 0, -1, -1, 46271 },
    { 2, 0, 0, 1, 32573 },
    { 0, 0, 2, 1, 17198 },
    { 2, 0, 1, -1, 9266 },
    { 0, 0, 2, -1, 8822 },
    { 2, -1, 0, -1, 8216 },
    { 2, 0, -2, -1, 4324 },
    { 2, 0, 1, 1, 4200 },
    { 2, 1, 0, -1, -3359 },
    { 2, -1, -1, 1, 2463 },
    { 2, -1, 0, 1, 2211 },
    { 2, -1, -1, -1, 2065 },
    { 0, 1, -1, -1, -1870 },
    { 4, 0, -1, -1, 1828 },
    { 0, 1, 0, 1, -1794 },
    { 0, 0, 0, 3, -1749 },
    { 0, 1, -1, 1, -1565 },
    { 1, 0, 0, 1, -1491 },
    { 0, 1, 1, 1, -1475 },
    { 0, 1, 1, -1, -1410 },
    { 0, 1, 0, -1, -1344 },
    { 1, 0, 0, -1, -1335 },
    { 0, 0, 3, 1, 1107 },
    { 4, 0, 0, -1, 1021 },
    { 4, 0, -1, 1, 833 },
    { 0, 0, 1, -3, 777 },
    { 4, 0, -2, 1, 671 },
    { 2, 0, 0, -3, 607 },
    { 2, 0, 2, -1, 596 },
    { 2, -1, 1, -1, 491 },
    { 2, 0, -2, 1, -451 },
    { 0, 0, 3, -1, 439 },
    { 2, 0, 2, 1, 422 },
    { 2, 0, -3, -1, 421 },
    { 2, 1, -1, 1, -366 },
    { 2, 1, 0, 1, -351 },
    { 4, 0, 0, 1, 331 },
    { 2, -1, 1, 1, 315 },
    { 2, -2, 0, -1, 302 },
    { 0, 0, 1, 3, -283 },
    { 2, 1, 1, -1, -229 },
    { 1, 1, 0, -1, 223 },
    { 1, 1, 0, 1, 223 },
    { 0, 1, -2, -1, -220 },
    { 2, 1, -1, -1, -220 },
    { 1, 0, 1, 1, -185 },
    { 2, -1, -2, -1, 181 },
    { 0, 1, 2, 1, -177 },
    { 4, 0, -2, -1, 176 },
    { 4, -1, -1, -1, 166 },
    { 1, 0, 1, -1, -164 },
    { 4, 0, 1, -1, 132 },
    { 1, 0, -1, -1, -119 },
    { 4, -1, 0, -1, 115 },
    { 2, -2, 0, 1, 107 }
};

#include "HEOSATDynamics_impl.h"

#endif
