//-----------------------------------------------------------------------
// FIRST-ORDER SHORT-PERIOD ANALYTICAL CORRECTIONS DUE TO THIRD-BODY EFFECTS
// (SECOND-ORDER LEGENDRE POLYNOMIAL FOR SUN AND MOON DISTURBING FUNCTION)
// 
// Original code (FORTRAN) by:
// Author : M.Lara(mlara0@gmail.com)
// Laboratory : Chiclana Consulting(ChicCon)
//
// C++ version by: D.J. Gondelach (University of Southampton)
//
// Reference:
// Efficient Computation of Short-Period Analytical Corrections Due to Third-Body Effects
// By: M.Lara, R.Vilhena de Moraes, D.M.Sanchez and A.F.B.de A.Prado. 2015. no. AAS 15-295
//
//-----------------------------------------------------------------------

#pragma once

#include <cmath>

#include <cspice/SpiceUsr.h>

// Compute Sun and Moon positions using SPICE
void sunParametersSpice(const double jd, const double ul, const double ut, double &xs, double &ys, double &zs, double &Rs, double &ns)
{
	// Sun
	const double AU = 149597870.700;
	const double as = 1.000001018*AU / ul;  // Sun semi-major axis

	ns = (2 * M_PI / (365.256363004 * 24 * 60 * 60))*ut;  // Sun mean motion

	const double et = unitim_c(jd, "JED", "ET");

	double aux[6];
	double lt;

	spkezr_c("SUN", et, "TOD", "NONE", "EARTH", aux, &lt);
	Rs = sqrt(pow(aux[0], 2) + pow(aux[1], 2) + pow(aux[2], 2));
	xs = aux[0] / Rs;
	ys = aux[1] / Rs;
	zs = aux[2] / Rs;
	Rs = Rs / ul / as;
}

// Compute Sun and Moon positions using SPICE
void moonParametersSpice(const double jd, const double ul, const double ut, double &xl, double &yl, double &zl, double &rm, double &nm, double &bet)
{
	const double am = 384400.0 / ul;  // Moon semi-major axis

	nm = (2 * M_PI) / (24 * 60 * 60)*(1 / 27.321582)*ut;  // Moon mean motion
	bet = 1.0 / 82.2845;

	const double et = unitim_c(jd, "JED", "ET");

	double aux[6];
	double lt;

	spkezr_c("MOON", et, "TOD", "NONE", "EARTH", aux, &lt);
	rm = std::sqrt(std::pow(aux[0], 2) + std::pow(aux[1], 2) + std::pow(aux[2], 2));
	xl = aux[0] / rm;
	yl = aux[1] / rm;
	zl = aux[2] / rm;
	rm = rm / ul / am;
}

template<typename T> void eccInclP2(const T e, const T y, const T ci, const T si, T ep[6][2], T ip[6][2][5])
{
	ep[2][0] = 6. + 9. * e*e;
	ep[2][1] = -15. * e*e;
	ep[1][0] = 0.;
	ep[1][1] = 2. * ep[2][1];
	ep[3][0] = ep[2][0];
	ep[3][1] = ep[2][1];
	ep[4][0] = 18. * y;
	ep[4][1] = -30. * y;
	ep[5][0] = ep[2][0] / y;
	ep[5][1] = ep[2][1] / y;
	
	ip[1][0][0] = -si*si / 32.;
	ip[1][0][1] = ci*si / 8.;
	ip[1][0][2] = (-1. + 3. * ci*ci) / 48.;
	ip[1][0][3] = ip[1][0][1];
	ip[1][0][4] = ip[1][0][0];
	
	ip[1][1][0] = pow(1. - ci, 2) / 32.;
	ip[1][1][1] = -(1. - ci)*si / 8.;
	ip[1][1][2] = -si*si / 16.;
	ip[1][1][3] = (1. + ci)*si / 8.;
	ip[1][1][4] = pow(1. + ci, 2) / 32.;
	
	for (int i = 0; i <= 1; i++)
	{
		for (int j = 0; j <= 4; j++)
		{
			ip[2][i][j] = j*ip[1][i][j];
			ip[3][i][j] = ip[1][i][j];
			ip[4][i][j] = ip[1][i][j];
		}
	}

	ip[5][0][0] = -ci / 16.;
	ip[5][0][1] = 1 / (8.*si) - si / 4.;
	ip[5][0][2] = -ci / 8.;
	ip[5][0][3] = ip[5][0][1];
	ip[5][0][4] = ip[5][0][0];
	ip[5][1][0] = (1. - ci) / 16.;
	ip[5][1][1] = -(si*(1. + 2. * ci) / (8. * (1. + ci)));
	ip[5][1][2] = -ci / 8.;
	ip[5][1][3] = -(si*(1. - 2. * ci) / (8.*(1. - ci)));
	ip[5][1][4] = (-1. - ci) / 16.;
}

template<typename T> void eccFunctP2(const T k, const T q, const T y, T sig[4][2], T kap[4][2]) //sigma, kappa
{
	// k=e*cos(trueAnom), q=e*sin(trueAnom)

	const T k2 = k*k;
	const T y2 = y*y;
	const T y4 = y2*y2;
	const T pir = 1. + k;
	const T pr2 = pir*pir;
	const T pr3 = pr2*pir;
	const T pr4 = pr2*pr2;
	const T pr5 = pr3*pr2;
	
	sig[0][0] = q*(5. * pir*(2. + k)*y2 + 2. * y4);
	sig[0][1] = q*(20. * k*pr3 + 5. * pir*(2. + 3. * k)*y2 - 2. * y4);
	kap[0][0] = 0.;
	kap[0][1] = -10. * pr3*(-1. + 2. * k2) - 5. * pr2*(1. + 5. * k)*y2 - 2. * pir*y4;
	
	sig[1][0] = -20. * (-1. + k)*pr2*(2. + k) - pir*(26. + 19. * k)*y2 - 2. * y4;
	sig[1][1] = -40. * pr2*(-1. + 2. * k + 4. * k2) - (120. * (-1. + k)*k*pr4) / y2
				- pir*(34. + 41. * k)*y2 + 2. * y4;
	kap[1][0] = 0.;
	kap[1][1] = q*(-20. * pr2*(1. + 5. * k)
				- (60. * pr3*(-1. + 2. * k2)) / y2 - 4. * pir*y2);
	
	sig[2][0] = q*(20. * pr3*(2. + k) + 14. * pr2*y2);
	sig[2][1] = q*(20. * pr3*(4. + 7. * k) + (120 * k*pr5) / y2 + 26. * pr2*y2);
	kap[2][0] = 0.;
	kap[2][1] = -10. * pr3*(1. + 20. * k + 20. * k2) - (60. * pr5*(-1. + 2. * k2)) / y2
				- pr2*(59. + 79. * k)*y2 - 2. * pir*y4;
	
	sig[3][0] = q*(-20. * k*pr2*(2. + k) + pir*(15. + k)*y2 + 6. * y4);
	sig[3][1] = q*(-20. * pr2*pow(1. + 2. * k,2)
				- (120. * k2*pr4) / y2 + pir*(5. + 19. * k)*y2 - 6. * y4);
	kap[3][0] = 0.;
	kap[3][1] = 20. * k*pr3*(3. + 7. * k) + (60. * k*pr4*(-1. + 2. * k2)) / y2
				+ 4. * pr2*(5. + k)*y2 - 4. * pir*y4;
}


void P2MOON(const double xl, const double yl, const double zl, const double rm, const double nm, const double bet, double Sm[5], double Tm[5], double& km)
{// cos, sin, k moon
	// bet: moon beta mass ratio
	// nm: moon's mean motion: scaled units

	const double b2 = bet*nm*nm;
	const double x2 = xl*xl;
	const double y2 = yl*yl;
	const double z2 = zl*zl;
	const double xy = xl*yl;
	const double xz = xl*zl;
	const double yz = yl*zl;
	// sin moon
	Tm[2] = 0.;
	Tm[3] = xz;
	Tm[4] = 2. * xy;
	Tm[1] = -Tm[3];
	Tm[0] = -Tm[4];
	// cos moon
	Sm[2] = -1. + 3. * z2;
	Sm[3] = -yz;
	Sm[4] = x2 - y2;
	Sm[1] = Sm[3];
	Sm[0] = Sm[4];
	
	km = b2*std::pow(1. / rm, 3); // rm relative to the Moon's semimajor axis
}


void P2SUN(const double xs, const double ys, const double zs, const double Rs, const double ns, double Ss[5], double Ts[5], double& ks)
{// cos, sin, k sun
	// sin sun
	Ts[2] = 0.;
	Ts[3] = xs*zs;
	Ts[4] = 2. * xs*ys;
	Ts[1] = -Ts[3];
	Ts[0] = -Ts[4];
	// cos sun
	Ss[2] = -1. + 3. * zs*zs;
	Ss[3] = -ys*zs;
	Ss[4] = xs*xs - ys*ys;
	Ss[1] = Ss[3];
	Ss[0] = Ss[4];
	
	ks = (ns*ns)*std::pow(1. / Rs, 3); // rs relative to the Sun's semimajor axis
}


template<typename T> void moonShortPeriodic2
(	const T r, const T th, const T nu, const T RR, const T ZZ, const T NN, 
	const double xl, const double yl, const double zl, const double rm, const double nm, const double bet, const double mu,
	T& dr, T& dth, T& dnu, T& dRR, T& dZZ, T& dNN)
{
	
	//Using std and DACE functions
	using std::pow;   using DACE::pow;
	using std::sqrt;  using DACE::sqrt;
	
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
	const int o = 2;
	const int jm = o / 2;
	const int i = o % 2;

	//-- - keep the previous corrections
	const T sr = dr;
	const T sth = dth;
	const T snu = dnu;
	const T sRR = dRR;
	const T sZZ = dZZ;
	const T sNN = dNN;
	//-- - compute auxiliary functions of the polar - nodal variables
	const T p = ZZ*ZZ / mu;
	const T k = (p / r) - 1.; // e*cos(trueAnom)
	const T q = p*RR / ZZ; // e*sin(trueAnom)
	const T e = sqrt(k*k + q*q);
	const T y = sqrt(1. - e*e);
	const T ci = NN / ZZ; // cos(incl)
	const T si = sqrt(1. - ci*ci); // sin(incl)
	//-- - 3rd body corrections
	T ep[6][2], ip[6][2][5];
	eccInclP2(e, y, ci, si, ep, ip); //IN: e, y, ci, si / OUT : ep, ip

	T sig[4][2], kap[4][2];
	eccFunctP2(k, q, y, sig, kap);

	double Sm[5], Tm[5], km;
	P2MOON(xl, yl, zl, rm, nm, bet, Sm, Tm, km);

	dZZ = 0.;
	dNN = 0.;
	dnu = 0.;
	dr = 0.;
	dth = 0.;
	dRR = 0.;

	for (int j = 0; j <= jm; j++)
	{
		T ZZd = 0.;
		T NNd = 0.;
		T nud = 0.;
		T rd = 0.;
		T thd = 0.;
		T RRd = 0.;
		double m = (double) 2 * j + i;

		for (int l = -o; l <= o; l++)
		{
			T au = kap[0][j] * Sm[l + 2] + sig[0][j] * Tm[l + 2];
			T bu = kap[0][j] * Tm[l + 2] - sig[0][j] * Sm[l + 2];
			T arg = m*th + l*nu;
			T cx = cos(arg);
			T sx = sin(arg);
			ZZd = ZZd + ip[1][j][l + 2] * (au*cx + bu*sx);
			NNd = NNd + l*ip[1][j][l + 2] * (au*cx + bu*sx);
			nud = nud + ip[5][j][l + 2] * (au*sx - bu*cx);

			au = kap[1][j] * Sm[l + 2] + sig[1][j] * Tm[l + 2];
			bu = kap[1][j] * Tm[l + 2] - sig[1][j] * Sm[l + 2];
			rd = rd + ip[1][j][l + 2] * (bu*cx - au*sx);
			au = kap[2][j] * Sm[l + 2] + sig[2][j] * Tm[l + 2];
			bu = kap[2][j] * Tm[l + 2] - sig[2][j] * Sm[l + 2];
			thd = thd + ip[1][j][l + 2] * (bu*cx - au*sx);
			au = kap[3][j] * Sm[l + 2] + sig[3][j] * Tm[l + 2];
			bu = kap[3][j] * Tm[l + 2] - sig[3][j] * Sm[l + 2];
			RRd = RRd + ip[1][j][l + 2] * (au*sx - bu*cx);
		}
		dZZ = dZZ + m*ZZd;
		dNN = dNN + NNd;
		dnu = dnu + nud;
		dr = dr + rd;
		dth = dth + thd;
		dRR = dRR + RRd;
	}
	
	T aux = ZZ*pow(y,3) / pow(p,2); //aux: mean motion
	T di = bet*pow(nm / aux, 2) * pow(r / p, 3) * pow(1. / rm, 3); //rm relative to Moon's semimajor axis

	dr = sr + p*di*dr;
	dth = sth + di*(dth - ci*dnu);
	dnu = snu + di*dnu;
	dRR = sRR + di*(ZZ / r)*dRR;
	dZZ = sZZ + di*ZZ*dZZ;
	dNN = sNN + di*ZZ*dNN;
}


template<typename T> void sunShortPeriodic2
(	const T r, const T th, const T nu, const T RR, const T ZZ, const T NN, 
	const double xs, const double ys, const double zs, const double Rs, const double ns, const double mu,
	T& dr, T& dth, T& dnu, T& dRR, T& dZZ, T& dNN)
{
	
	//Using std and DACE functions
	using std::pow;   using DACE::pow;
	using std::sqrt;  using DACE::sqrt;
	
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
	const int o = 2;
	const int jm = o / 2;
	const int i = o % 2;
	
	const T sr = dr;
	const T sth = dth;
	const T snu = dnu;
	const T sRR = dRR;
	const T sZZ = dZZ;
	const T sNN = dNN;
	//-- - compute auxiliary functions of the polar - nodal variables
	const T p = ZZ*ZZ / mu;
	const T k = (p / r) - 1.;
	const T q = p*RR / ZZ;
	const T e = sqrt(k*k + q*q);
	const T y = sqrt(1. - e*e);
	const T ci = NN / ZZ;
	const T si = sqrt(1. - ci*ci);
	//-- - 3rd body corrections
	T ep[6][2], ip[6][2][5];
	eccInclP2(e, y, ci, si, ep, ip); //IN: e, y, ci, si / OUT : ep, ip

	T sig[4][2], kap[4][2];
	eccFunctP2(k, q, y, sig, kap);

	double Ss[5], Ts[5], ks;
	P2SUN(xs, ys, zs, Rs, ns, Ss, Ts, ks);
	
	dZZ = 0.;
	dNN = 0.;
	dnu = 0.;
	dr = 0.;
	dth = 0.;
	dRR = 0.;
	for (int j = 0; j <= jm; j++)
	{
		T ZZd = 0.;
		T NNd = 0.;
		T nud = 0.;
		T rd = 0.;
		T thd = 0.;
		T RRd = 0.;
		double m = (double) 2 * j + i;

		for (int l = -o; l <= o; l++)
		{
			T au = kap[0][j] * Ss[l + 2] + sig[0][j] * Ts[l + 2];
			T bu = kap[0][j] * Ts[l + 2] - sig[0][j] * Ss[l + 2];
			T arg = m*th + l*nu;
			T cx = cos(arg);
			T sx = sin(arg);
			ZZd = ZZd + ip[1][j][l + 2] * (au*cx + bu*sx);
			NNd = NNd + l*ip[1][j][l + 2] * (au*cx + bu*sx);
			nud = nud + ip[5][j][l + 2] * (au*sx - bu*cx);

			au = kap[1][j] * Ss[l + 2] + sig[1][j] * Ts[l + 2];
			bu = kap[1][j] * Ts[l + 2] - sig[1][j] * Ss[l + 2];
			rd = rd + ip[1][j][l + 2] * (bu*cx - au*sx);
			au = kap[2][j] * Ss[l + 2] + sig[2][j] * Ts[l + 2];
			bu = kap[2][j] * Ts[l + 2] - sig[2][j] * Ss[l + 2];
			thd = thd + ip[1][j][l + 2] * (bu*cx - au*sx);
			au = kap[3][j] * Ss[l + 2] + sig[3][j] * Ts[l + 2];
			bu = kap[3][j] * Ts[l + 2] - sig[3][j] * Ss[l + 2];
			RRd = RRd + ip[1][j][l + 2] * (au*sx - bu*cx);
		}
		dZZ = dZZ + m*ZZd;
		dNN = dNN + NNd;
		dnu = dnu + nud;
		dr = dr + rd;
		dth = dth + thd;
		dRR = dRR + RRd;
	}
	
	T aux = ZZ*pow(y, 3) / pow(p, 2); //aux: sat mean motion
	T di = 1. * pow(ns / aux, 2) * pow(r / p, 3) * pow(1. / Rs, 3); //rs relative to sun's semimajor axis

	dr = sr + p*di*dr;
	dth = sth + di*(dth - ci*dnu);
	dnu = snu + di*dnu;
	dRR = sRR + di*(ZZ / r)*dRR;
	dZZ = sZZ + di*ZZ*dZZ;
	dNN = sNN + di*ZZ*dNN;
}