#pragma once

#include <cmath>

// SPICE
#include <cspice/SpiceUsr.h>

// Dynamics
#include "Conversions.h"

template<typename T>
DACE::AlgebraicVector<T> J2000ToTOD(const DACE::AlgebraicVector<T> &stateJ2000vec, const double jdate)
{
	DACE::AlgebraicVector<T> rrJ2000(3), vvJ2000(3);
	for (int i = 0; i<3; i++)
	{
		rrJ2000[i] = stateJ2000vec[i];
		vvJ2000[i] = stateJ2000vec[i + 3];
	}

	const double et = unitim_c(jdate, "JED", "ET");

	double PREC[6][6], NUT[6][6];
	sxform_c("J2000", "MOD", et, PREC);
	sxform_c("MOD", "TOD", et, NUT);
	
	
	DACE::AlgebraicMatrix<T> NUT3(3,3), PREC3(3,3);
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			NUT3.at(i,j) = NUT[i][j];
			PREC3.at(i,j) = PREC[i][j];
		}
	}
	
	DACE::AlgebraicVector<T> rrMOD(3), rrTOD(3), vvMOD(3), vvTOD(3);
	rrMOD=PREC3*rrJ2000;
	vvMOD=PREC3*vvJ2000;
	rrTOD=NUT3*rrMOD;
	vvTOD=NUT3*vvMOD;
	
	DACE::AlgebraicVector<T> stateTODvec(6);
	for (int i = 0; i<3; i++)
	{
		stateTODvec[i] = rrTOD[i];
		stateTODvec[i + 3] = vvTOD[i];
	}

	return stateTODvec;
}

template<typename T>
DACE::AlgebraicVector<T> TODToJ2000(const DACE::AlgebraicVector<T> &stateTODvec, const double jdate)
{

	DACE::AlgebraicVector<T> rrTOD(3), vvTOD(3);	
	for (int i = 0; i<3; i++){
		rrTOD[i] = stateTODvec[i];
		vvTOD[i] = stateTODvec[i+3];
	}

	const double et = unitim_c(jdate, "JED", "ET");

	double NUT[6][6], PREC[6][6];
	sxform_c("TOD", "MOD", et, NUT);
	sxform_c("MOD", "J2000", et, PREC);

	DACE::AlgebraicMatrix<T> NUT3(3,3), PREC3(3,3);
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			NUT3.at(i,j) = NUT[i][j];
			PREC3.at(i,j) = PREC[i][j];
		}
	}
	
	DACE::AlgebraicVector<T> rrMOD(3), rrJ2000(3), vvMOD(3), vvJ2000(3);
	rrMOD=NUT3*rrTOD;
	vvMOD=NUT3*vvTOD;
	rrJ2000=PREC3*rrMOD;
	vvJ2000=PREC3*vvMOD;

	DACE::AlgebraicVector<T> stateJ2000vec(6);
	for (int i = 0; i<3; i++){
		stateJ2000vec[i] = rrJ2000[i];
		stateJ2000vec[i+3] = vvJ2000[i];
	}

	return stateJ2000vec;
}

template<typename T>
DACE::AlgebraicVector<T> kepJ2000ToTOD(const DACE::AlgebraicVector<T> &kepJ2000, const double jdate, const double mu)
{
	DACE::AlgebraicVector<T> cartJ2000 = kepM2cart(kepJ2000, mu);
	DACE::AlgebraicVector<T> cartTOD = J2000ToTOD(cartJ2000, jdate);
	DACE::AlgebraicVector<T> kepTOD = cart2kepM(cartTOD, mu);

	return kepTOD;
}

template<typename T>
DACE::AlgebraicVector<T> kepTODToJ2000(const DACE::AlgebraicVector<T> &kepTOD, const double jdate, const double mu)
{
	DACE::AlgebraicVector<T> cartTOD = kepM2cart(kepTOD, mu);
	DACE::AlgebraicVector<T> cartJ2000 = TODToJ2000(cartTOD, jdate);
	DACE::AlgebraicVector<T> kepJ2000 = cart2kepM(cartJ2000, mu);

	return kepJ2000;
}
