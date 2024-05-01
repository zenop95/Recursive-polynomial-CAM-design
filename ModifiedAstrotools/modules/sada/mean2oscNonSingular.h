#pragma once

#include "luniSolarShortPeriodics.h"
#include <dace/dace.h>

template<typename T> DACE::AlgebraicVector<T> hill2kep(const DACE::AlgebraicVector<T>& hill, const double mu)
//  hill[] = {r, th, nu, R, Th, Nu}
//  kep[] = {a,e,i,RAAN,argPer,trueAnom}
{
	//Using std and DACE functions
	using std::pow;   using DACE::pow;
	using std::sqrt;  using DACE::sqrt;
	
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	using std::acos;  using DACE::acos;
	
	DACE::AlgebraicVector<T> kep(6);

	T r = hill[0];
	T th = hill[1];
	T nu = hill[2];
	T R = hill[3];
	T Th = hill[4];
	T Nu = hill[5];

	T i = acos(Nu / Th);
	T cs = (-1.0 + pow(Th, 2) / (mu*r))*cos(th) + (R*Th*sin(th)) / mu;
	T ss = -((R*Th*cos(th)) / mu) + (-1.0 + pow(Th, 2) / (mu*r))*sin(th);
	T e = sqrt(cs*cs + ss * ss);
	T p = Th * Th / mu;
	T costrue = 1.0 / e * (p / r - 1.0);
	T f = acos(costrue);

	if (DACE::cons(R)<0.0) {
		f = 2.0*M_PI - f;
	}
	T a = p / (1 - e * e);

	kep[0] = a;
	kep[1] = e;
	kep[2] = i;
	kep[3] = nu;
	kep[4] = th - f;
	kep[5] = f;

	return kep;

}

template<typename T> DACE::AlgebraicVector<T> kep2hill(const DACE::AlgebraicVector<T>& kep, const double mu)
//  hill[] = {r, th, nu, R, Th, Nu}
//  kep[] = {a,e,i,RAAN,argPer,trueAnom}
{
	//Using std and DACE functions
	using std::sqrt;  using DACE::sqrt;
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	
	DACE::AlgebraicVector<T> hill(6);
	T f, p;

	p = kep[0] * (1.0 - kep[1] * kep[1]);
	f = kep[5];

	hill[4] = sqrt(mu*p);
	hill[0] = p / (1.0 + kep[1] * cos(f));
	hill[1] = f + kep[4];
	hill[2] = kep[3];
	hill[3] = (hill[4] / p)*kep[1] * sin(f);
	hill[5] = hill[4] * cos(kep[2]);

	return hill;
}

template<typename T> DACE::AlgebraicVector<T> osculating2meanHill(DACE::AlgebraicVector<T> hillOsc, 
		double mu, double J2, double rE, double cont = 1.0)
{
	
	//Using std and DACE functions
	using std::pow;   using DACE::pow;
	using std::sqrt;  using DACE::sqrt;
	
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	using std::acos;  using DACE::acos;

	
	// Osculating Hill elements
	const T r = hillOsc[0];
	const T th = hillOsc[1];
	const T nu = hillOsc[2];
	const T R = hillOsc[3];
	const T Th = hillOsc[4];
	const T Nu = hillOsc[5];

	const T ThInv = 1. / Th;
	const T rInv = 1. / r;

	const T ci = Nu * ThInv;
	const T si = sqrt(1.0 - ci * ci);
	const T cs = (-1.0 + pow(Th, 2) / mu * rInv)*cos(th) + (R*Th*sin(th)) / mu;
	const T ss = -((R*Th*cos(th)) / mu) + (-1.0 + pow(Th, 2) / mu * rInv)*sin(th);
	const T e = sqrt(cs*cs + ss * ss);
	const T eta = sqrt(1.0 - e * e);

	const T beta = 1.0 / (1.0 + eta);

	const T p = Th * Th / mu;
	const T costrue = 1. / e * (p*rInv - 1.);

	T f = acos(costrue);

	if (DACE::cons(R)<0.0) {
		f = 2.0*M_PI - f;
	}

	const T M = true2meanAnomaly(f, e);

	const T phi = f - M;

	// Variables to speed up calculations
	const T si2 = pow(si, 2);
	const double rE2J2 = pow(rE, 2)*J2;
	const double rE2J2mu = rE2J2 * mu;
	const double rE2J2mu2 = rE2J2 * pow(mu, 2);
	const T betaR = beta * R;
	const T betaRsi2 = betaR * si2;
	const T c2th = cos(2.*th);
	const T s2th = sin(2.*th);
	const T Rsi2 = R * si2;
	const T Th2inv = ThInv * ThInv;
	const T Th3inv = Th2inv * ThInv;
	const T Th4inv = Th3inv * ThInv;
	const T r2Inv = rInv * rInv;
	const T rThInv = rInv * ThInv;
	const T rTh2inv = rInv * Th2inv;

	// Convert from osculating to mean
	DACE::AlgebraicVector<T> hillMean(6);

	hillMean[0] = r + (cont)*((rE2J2*beta) / 2.*rInv - (3.*rE2J2*beta*si2) / 4.*rInv +
		(rE2J2mu2*eta*r)*Th4inv - (3.*rE2J2mu2*eta*r*si2) / 2.*Th4inv +
		(rE2J2mu) / 2.*Th2inv - (rE2J2mu*beta) / 2.*Th2inv -
		(3.*rE2J2mu*si2) / 4.*Th2inv + (3.*rE2J2mu*beta*si2) / 4.*Th2inv -
		(rE2J2mu*si2*c2th) / 4.*Th2inv);


	hillMean[1] = th + (cont)*((-3.*rE2J2mu2*phi)*Th4inv + (15.*rE2J2mu2*phi*si2) / 4.*Th4inv -
		(5.*rE2J2mu*R) / 2.*Th3inv - (rE2J2mu*betaR) / 2.*Th3inv +
		(3.*rE2J2mu*Rsi2)*Th3inv + (3.*rE2J2mu*betaRsi2) / 4.*Th3inv -
		(rE2J2*betaR) / 2.*rThInv + (3.*rE2J2*betaRsi2) / 4.*rThInv +
		(-(rE2J2mu*R) / 2.*Th3inv + (rE2J2mu*Rsi2)*Th3inv)*c2th +
		(-(rE2J2mu2) / 4.*Th4inv + (5.*rE2J2mu2*si2) / 8.*Th4inv +
		(rE2J2mu)*rTh2inv - (3.*rE2J2mu*si2) / 2.*rTh2inv)*s2th);

	hillMean[2] = nu + (cont)*((3.*rE2J2mu2*ci*phi) / 2.*Th4inv + (3.*rE2J2mu*ci*R) / 2.*Th3inv +
		(rE2J2mu*ci*R*c2th) / 2.*Th3inv +
		((rE2J2mu2*ci) / 4.*Th4inv - (rE2J2mu*ci)*rTh2inv)*s2th);


	hillMean[3] = R + (cont)*(-(rE2J2*betaR) / 2.*r2Inv + (3.*rE2J2*betaRsi2) / 4.*r2Inv -
		(rE2J2mu2*eta*R) / 2.*Th4inv + (3.*rE2J2mu2*eta*Rsi2) / 4.*Th4inv +
		(rE2J2mu*si2*s2th) / 2.*r2Inv*ThInv);


	hillMean[4] = Th + (cont)*(((rE2J2mu2*si2) / 4.*Th3inv - (rE2J2mu*si2)*rThInv)*c2th -
		(rE2J2mu*Rsi2*s2th) / 2.*Th2inv);

	hillMean[5] = Nu + 0.;

	return hillMean;
}


template<typename T> DACE::AlgebraicVector<T> mean2osculatingHill(DACE::AlgebraicVector<T> hillMean, 
		double mu, double J2, double rE, double cont = 1.0)
{

	//Using std and DACE functions
	using std::pow;   using DACE::pow;
	using std::sqrt;  using DACE::sqrt;
	
	using std::sin;   using DACE::sin;
	using std::cos;   using DACE::cos;
	using std::acos;  using DACE::acos;
	
	// Mean Hill elements
	const T r = hillMean[0];
	const T th = hillMean[1];
	const T nu = hillMean[2];
	const T R = hillMean[3];
	const T Th = hillMean[4];
	const T Nu = hillMean[5];

	const T ThInv = 1. / Th;
	const T rInv = 1. / r;

	const T ci = Nu * ThInv;
	const T si = sqrt(1.0 - ci * ci);
	const T cs = (-1.0 + pow(Th, 2) / mu * rInv)*cos(th) + (R*Th*sin(th)) / mu;
	const T ss = -((R*Th*cos(th)) / mu) + (-1.0 + pow(Th, 2) / mu * rInv)*sin(th);
	const T e = sqrt(cs*cs + ss * ss);
	const T eta = sqrt(1.0 - e * e);
	const T beta = 1.0 / (1.0 + eta);

	const T p = Th * Th / mu;
	const T costrue = 1. / e * (p*rInv - 1.);

	T f = acos(costrue);

	if (DACE::cons(R)<0.0) {
		f = 2.0*M_PI - f;
	}

	const T M = true2meanAnomaly(f, e);

	const T phi = f - M;

	// Variables to speed up calculations
	const T si2 = pow(si, 2);
	const double rE2J2 = pow(rE, 2)*J2;
	const double rE2J2mu = rE2J2 * mu;
	const double rE2J2mu2 = rE2J2 * pow(mu, 2);
	const T betaR = beta * R;
	const T betaRsi2 = betaR * si2;
	const T c2th = cos(2.*th);
	const T s2th = sin(2.*th);
	const T Rsi2 = R * si2;
	const T Th2inv = ThInv * ThInv;
	const T Th3inv = Th2inv * ThInv;
	const T Th4inv = Th3inv * ThInv;
	const T r2Inv = rInv * rInv;
	const T rThInv = rInv * ThInv;
	const T rTh2inv = rInv * Th2inv;

	// Convert from mean to osculating
	DACE::AlgebraicVector<T> hillOsc(6);

	hillOsc[0] = r - (cont)*((rE2J2*beta) / 2.*rInv - (3.*rE2J2*beta*si2) / 4.*rInv +
		(rE2J2mu2*eta*r)*Th4inv - (3.*rE2J2mu2*eta*r*si2) / 2.*Th4inv +
		(rE2J2mu) / 2.*Th2inv - (rE2J2mu*beta) / 2.*Th2inv -
		(3.*rE2J2mu*si2) / 4.*Th2inv + (3.*rE2J2mu*beta*si2) / 4.*Th2inv -
		(rE2J2mu*si2*c2th) / 4.*Th2inv);


	hillOsc[1] = th - (cont)*((-3.*rE2J2mu2*phi)*Th4inv + (15.*rE2J2mu2*phi*si2) / 4.*Th4inv -
		(5.*rE2J2mu*R) / 2.*Th3inv - (rE2J2mu*betaR) / 2.*Th3inv +
		(3.*rE2J2mu*Rsi2)*Th3inv + (3.*rE2J2mu*betaRsi2) / 4.*Th3inv -
		(rE2J2*betaR) / 2.*rThInv + (3.*rE2J2*betaRsi2) / 4.*rThInv +
		(-(rE2J2mu*R) / 2.*Th3inv + (rE2J2mu*Rsi2)*Th3inv)*c2th +
		(-(rE2J2mu2) / 4.*Th4inv + (5.*rE2J2mu2*si2) / 8.*Th4inv +
		(rE2J2mu)*rTh2inv - (3.*rE2J2mu*si2) / 2.*rTh2inv)*s2th);

	hillOsc[2] = nu - (cont)*((3.*rE2J2mu2*ci*phi) / 2.*Th4inv + (3.*rE2J2mu*ci*R) / 2.*Th3inv +
		(rE2J2mu*ci*R*c2th) / 2.*Th3inv +
		((rE2J2mu2*ci) / 4.*Th4inv - (rE2J2mu*ci)*rTh2inv)*s2th);


	hillOsc[3] = R - (cont)*(-(rE2J2*betaR) / 2.*r2Inv + (3.*rE2J2*betaRsi2) / 4.*r2Inv -
		(rE2J2mu2*eta*R) / 2.*Th4inv + (3.*rE2J2mu2*eta*Rsi2) / 4.*Th4inv +
		(rE2J2mu*si2*s2th) / 2.*r2Inv*ThInv);


	hillOsc[4] = Th - (cont)*(((rE2J2mu2*si2) / 4.*Th3inv - (rE2J2mu*si2)*rThInv)*c2th -
		(rE2J2mu*Rsi2*s2th) / 2.*Th2inv);

	hillOsc[5] = Nu + 0.;

	return hillOsc;
}

template<typename T> DACE::AlgebraicVector<T> mean2oscJ2sunMoonNonSingular(
		DACE::AlgebraicVector<T> COEmean, const double jdate, const double mu, 
		const int setting, const double ul = 1.0, const double ut = 1.0)
{
    // Convert mean to true anomaly
	DACE::AlgebraicVector<T> COEmeanT = COEmean;
	COEmeanT[5] = mean2trueAnomaly(COEmean[5], COEmean[1]);
	
    // Convert to Hill variables
	DACE::AlgebraicVector<T> hillMean = kep2hill(COEmeanT, mu);

    // TODO: provide J2 and rE as input
	double J2 = 1.0826356665511e-3;
	double rE = 6378.1363 / ul;

	// Convert from mean to osculating considering J2 short-periodics
	DACE::AlgebraicVector<T> hillOsc = mean2osculatingHill(hillMean, mu, J2, rE, 1.0);
	if (setting == 0)
	{
        // Return osculating elements considering only J2 short-periodics
		DACE::AlgebraicVector<T> COEoscJ2 = hill2kep(hillOsc, mu);
		COEoscJ2[5] = true2meanAnomaly(COEoscJ2[5], COEoscJ2[1]);
		return COEoscJ2;
	}

    // Hill variables
    T r = hillMean[0];
    T th = hillMean[1];
    T nu = hillMean[2];
    T RR = hillMean[3];
    T ZZ = hillMean[4];
    T NN = hillMean[5];

	// Sun position and parameters
	double xs, ys, zs, Rs, ns;
	sunParametersSpice(jdate, ul, ut, xs, ys, zs, Rs, ns);

	// Moon position and parameters
	double xl, yl, zl, rm, nm, bet;
	moonParametersSpice(jdate, ul, ut, xl, yl, zl, rm, nm, bet);

	// Sun short periodics
	T dr, dth, dnu, dRR, dZZ, dNN;
	sunShortPeriodic2(r, th, nu, RR, ZZ, NN, xs, ys, zs, Rs, ns, mu, 
			dr, dth, dnu, dRR, dZZ, dNN);
    
	// Moon short periodics
	moonShortPeriodic2(r, th, nu, RR, ZZ, NN, xl, yl, zl, rm, nm, bet, mu, 
			dr, dth, dnu, dRR, dZZ, dNN);

    // Add Sun and Moon short-periodics
	hillOsc[0] += dr;
	hillOsc[1] += dth;
	hillOsc[2] += dnu;
	hillOsc[3] += dRR;
	hillOsc[4] += dZZ;
	hillOsc[5] += dNN;
    
    // Convert to Keplerian orbital elements
	DACE::AlgebraicVector<T> COEoscJ2sunMoon = hill2kep(hillOsc, mu);
	COEoscJ2sunMoon[5] = true2meanAnomaly(COEoscJ2sunMoon[5], COEoscJ2sunMoon[1]);

	return COEoscJ2sunMoon;
}

template<typename T> DACE::AlgebraicVector<T> osc2meanJ2sunMoonNonSingular(
		DACE::AlgebraicVector<T> COEosc, const double jdate, const double mu, 
		const int setting, const double ul = 1.0, const double ut = 1.0)
{
    // Convert mean to true anomaly
    DACE::AlgebraicVector<T> COEoscT = COEosc;
    COEoscT[5] = mean2trueAnomaly(COEosc[5], COEosc[1]);
    
    // Convert to Hill variables
    DACE::AlgebraicVector<T> hillOsc = kep2hill(COEoscT, mu);
    
    // TODO: provide J2 and rE as input
    double J2 = 1.0826356665511e-3;
    double rE = 6378.1363 / ul;
    
    // Convert from osculating to mean considering J2 short-periodics
    DACE::AlgebraicVector<T> hillMean = osculating2meanHill(hillOsc, mu, J2, rE, 1.0);
    if (setting == 0)
    {
        // Return mean elements considering only J2 short-periodics
        DACE::AlgebraicVector<T> COEmeanJ2 = hill2kep(hillMean, mu);
        COEmeanJ2[5] = true2meanAnomaly(COEmeanJ2[5], COEmeanJ2[1]);
        return COEmeanJ2;
    }
    
    // Hill variables
    T r = hillOsc[0];
    T th = hillOsc[1];
    T nu = hillOsc[2];
    T RR = hillOsc[3];
    T ZZ = hillOsc[4];
    T NN = hillOsc[5];
    
	// Sun position and parameters
	double xs, ys, zs, Rs, ns;
	sunParametersSpice(jdate, ul, ut, xs, ys, zs, Rs, ns);

	// Moon position and parameters
	double xl, yl, zl, rm, nm, bet;
	moonParametersSpice(jdate, ul, ut, xl, yl, zl, rm, nm, bet);

	// Sun short periodics
	T dr, dth, dnu, dRR, dZZ, dNN;
	sunShortPeriodic2(r, th, nu, RR, ZZ, NN, xs, ys, zs, Rs, ns, mu, dr, dth, dnu, dRR, dZZ, dNN);
    
	// Moon short periodics
	moonShortPeriodic2(r, th, nu, RR, ZZ, NN, xl, yl, zl, rm, nm, bet, mu, dr, dth, dnu, dRR, dZZ, dNN);

    // Remove Sun and Moon short-periodics
    hillMean[0] -= dr;
    hillMean[1] -= dth;
    hillMean[2] -= dnu;
    hillMean[3] -= dRR;
    hillMean[4] -= dZZ;
    hillMean[5] -= dNN;
    
    // Convert to Keplerian orbital elements
    DACE::AlgebraicVector<T> COEmeanJ2sunMoon = hill2kep(hillMean, mu);
    COEmeanJ2sunMoon[5] = true2meanAnomaly(COEmeanJ2sunMoon[5], COEmeanJ2sunMoon[1]);

	return COEmeanJ2sunMoon;
}
