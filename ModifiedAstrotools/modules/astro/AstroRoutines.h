/*
 * AstroRoutines.h
 *
 *  Created on: Oct. 20, 2014
 *      Author: Dinamica
 *  Modified on: March 2016
 *      Author: Mauro Massari (mauo.massari@polimi.it)
 */

#ifndef astro_AstroRoutines_H
#define astro_AstroRoutines_H

#include <dace/dace.h>
#include "Reference.h"
#include "AstroCnvRef.h"
#include <limits>


using namespace DACE;
using namespace std;

namespace astro {


/*! \defgroup TimeConversions Time Conversions
 * @brief Routines dedicated to time transformations.
 *
 * This module groups all routines dedicated to time transformations.
 *
 * Ephemeris time are seconds past J2000 (JD 2451545.0, 1 January 2000, 12:00:00.0)
 * MJD2000 are days past is referred to JD 2451544.5, i.e. mjd2000 = jd - 2451544.5
 *
 * @{
 */

double str2et(const std::string &utc)
{
  /*! Returns the ephemeris time corresponding to a time string.
   * Wrapper of CSPICE function str2et_c.
   * @param[in] utc A formatted string representing time.
   * @return Ephemeris time
   */

  double et;
  str2et_c(utc.c_str(), &et);
  return et;

}

std::string et2str(const double &et, const std::string &FORMAT)
{
  /*! Returns the UTC string corresponding to ephemeris time.
   * Wrapper of CSPICE function timout_c.
    * @param[in] et Ephemeris time
    * @param[in] FORMAT String specifying the time format (SPICE)
    * @return A formatted string containing the UTC time.
    */

  char tmp[79];
  timout_c(et,FORMAT.c_str(),LENOUT,tmp);
  std::string aux(tmp);
  return aux;

}

std::string et2isod(const double &et)
{
  char tmp[79];
  et2utc_c(et, "ISOD",14,LENOUT,tmp);
  std::string aux(tmp);
  return aux;
}

std::string et2utc(const double &et)
{
  /*! Returns the UTC string corresponding to ephemeris time.
    * @param[in] et Ephemeris time
    * @return A formatted string containing the UTC time.
    * @sa et2str
    */

    return et2str(et, UTC_str);
    // char tmp[79];
    // et2utc_c(et,"C",6,LENOUT,tmp);
    // std::string aux(tmp);
    // return aux;

}

double utc2et(std::string &utc)
{
  /*! Returns the ephemeris time corresponding to a UTC string
    * @param[in] utc A formatted string representing the UTC time.
    * @return Ephemeris time
    * @sa str2et
    */

    return str2et(utc);

}

double jd2mjd( const double &jd )
{
  /*! Returns the modified julian date corresponding to a julian date.
    * @param[in] jd A double containing the julian date.
    * @return The corresponding modified julian date
    */

  return ( jd - (j2000_c()-0.5) );

}

double mjd2jd( const double &mjd )
{
  /*! Returns the modified julian date corresponding to a julian date.
    * @param[in] mjd A double containing the julian date.
    * @return The corresponding modified julian date
    */

  return ( mjd + (j2000_c()-0.5) );

}

double jed2et( const double julianDate )
{
  /*! Convert a Julian date to ephemeris time (equivfalent to TDB in Spice).
   * \param[in] julianDate a double for the Julian date.
   * \return et a double for the ephemeris time.
   */

  return ( julianDate - j2000_c() ) * spd_c();

}


double et2jed( const double et )
{
  /*! Convert ephemeris time (equivalent to TDB) to a Julian date.
   * \param[in] et a double for the ephemeris time.
   * \return jed a double for the Julian date.
   */

  return j2000_c() + ( et ) / spd_c();

}

double mjd2et( const double &mjd )
{
  /*! Convert a modfied julian date to ephemeris time (equivalent to TDB in Spice).
   *  \param[in] mjd a double for the modified julian date.
   *  \return et a double for the ephemeris time.
   */

  return (mjd-0.5) * spd_c();

}


double et2mjd( const double et )
{
  /*! Convert ephemeris time (equivalent to TDB) to a modified julian date.
   *  \param et a double for the ephemeris time.
   *  \return mjd a double for the modified julian date.
   */

  return  et / spd_c() + 0.5;

}

/*! @} */ // End of time transformations

/***********************************************************************************
*     Rotation matrix
************************************************************************************/
template<class T>
DACE::AlgebraicMatrix<T> RotationMatrix(const T& angle, const unsigned int &idx )
{
  /*! Returns the \f$ (3 \times 3) \f$ rotation matrix around an axis.
    * @param angle The rotation angle (in radians).
    * @param idx   The ID of the rotation axis.
    * @return An AlgebraicMatrix that can rotate a vector around axis ID of the angle.
    */

    USEDACE_TRIGON

    T cosangle = cos(angle);
    T sinangle = sin(angle);

    DACE::AlgebraicMatrix<T> ROT(3,3,0.);

    switch (idx)
    {
      case 1:
	// Rotation around X axis
	ROT.at(0,0) =  1.;
	ROT.at(1,1) =  cosangle;
	ROT.at(1,2) =  sinangle;
	ROT.at(2,1) = -sinangle;
	ROT.at(2,2) =  cosangle;
	break;

      case 2:
	// Rotation around Y axis
	ROT.at(0,0) =  cosangle;
	ROT.at(0,2) = -sinangle;
	ROT.at(1,1) =  1.;
	ROT.at(2,0) =  sinangle;
	ROT.at(2,2) =  cosangle;
	break;

      case 3:
	// Rotation around Z axis
	ROT.at(0,0) =  cosangle;
	ROT.at(0,1) =  sinangle;
	ROT.at(1,0) = -sinangle;
	ROT.at(1,1) =  cosangle;
	ROT.at(2,2) =  1.;
	break;

    }

    return ROT;

  }


DACE::AlgebraicMatrix<double> TransformationMatrix(const std::string &FROM, const std::string &TO, const double &et)
{
  /*! Computes the transformation matrix from a frame to another.
   * Wrapper of SPICE routine sxform_c.
   * @param FROM Name of the frame to transform from.
   * @param TO   Name of the frame to transform to.
   * @param et   Epoch of the state transformation matrix.
   * @return A \f$ (6\times 6) \f$ state transformation matrix.
   */

  double aux[6][6];
  DACE::AlgebraicMatrix<double> tmp(6,6,0.);

  sxform_c(FROM.c_str(), TO.c_str(), et, aux);


  for (unsigned int i=0; i<6; i++)
    for (unsigned int j=0; j<6; j++)
      tmp.at(i,j) = aux[i][j];

  return tmp;

}


DACE::AlgebraicVector<double> RelativeState(const std::string &TARG, const std::string &REF, const double &et)
{
  /*! Computes the relative state between two celestial bodies.
   * SF wrapper for SPICE routine spkezr_c.
   * The routine returns the position of a target body relative to a reference body at a certain epoch in an equatorial inertial frame (J2000_EQT).
   * @param TARG Name of the target body.
   * @param REF  Name of the reference (observing) body.
   * @param et   Epoch of the state transformation matrix.
   * @return The cartesian state (position and velocity) in J2000_EQT.
   */

  double aux[6];
  double lt;

  spkezr_c(TARG.c_str(), et, FrameNameLabel[J2000_EQT], "NONE", REF.c_str(), aux, &lt);

  DACE::AlgebraicVector<double> tmp(6);

  for (unsigned int i=0; i<6; i++)
    tmp[i] = aux[i];

  return tmp;

}

DACE::AlgebraicVector<double> RelativeState(const std::string &TARG, const Frame F, const double &et)
{
  /*! Computes the relative state between two celestial bodies.
   * SF wrapper for SPICE routine spkezr_c.
   * The routine returns the position of a target body relative to a reference body at a certain epoch in an inertial frame.
   * @param TARG Name of the target body.
   * @param F    Frame structure with frame name and center.
   * @param et   Epoch of the state transformation matrix.
   * @return The cartesian state (position and velocity) in J2000_EQT.
   */

  double aux[6];
  double lt;

  spkezr_c(TARG.c_str(), et,F.getName(),"NONE",F.getOrigin(),aux,&lt);

  DACE::AlgebraicVector<double> tmp(6);

  for (unsigned int i=0; i<6; i++)
    tmp[i] = aux[i];

  return tmp;

}

DACE::AlgebraicVector<double> RelativePosition(const std::string &TARG, const std::string &REF, const double &et)
{
  /*! Computes the relative state between two celestial bodies.
   * SF wrapper for SPICE routine spkezr_c.
   * The routine returns the position of a target body relative to a reference body at a certain epoch in an equatorial inertial frame (J2000_EQT).
   * @param TARG Name of the target body.
   * @param REF  Name of the reference (observing) body.
   * @param et   Epoch of the state transformation matrix.
   * @return The cartesian state (position and velocity) in J2000_EQT.
   */

  double aux[6];
  double lt;

  spkezr_c(TARG.c_str(), et, FrameNameLabel[J2000_EQT], "NONE", REF.c_str(), aux, &lt);

  DACE::AlgebraicVector<double> tmp(3);

  for (unsigned int i=0; i<3; i++)
    tmp[i] = aux[i];

  return tmp;

}

DACE::AlgebraicVector<double> RelativePosition(const std::string &TARG, const Frame F, const double &et)
{
  /*! Computes the relative state between two celestial bodies.
   * SF wrapper for SPICE routine spkezr_c.
   * The routine returns the position of a target body relative to a reference body at a certain epoch in an inertial frame.
   * @param TARG Name of the target body.
   * @param F    Frame structure with frame name and center.
   * @param et   Epoch of the state transformation matrix.
   * @return The cartesian state (position and velocity) in J2000_EQT.
   */

  double aux[6];
  double lt;

  spkezr_c(TARG.c_str(), et,F.getName(),"NONE",F.getOrigin(),aux,&lt);

  DACE::AlgebraicVector<double> tmp(3);

  for (unsigned int i=0; i<3; i++)
    tmp[i] = aux[i];

  return tmp;

}

void pleph(const char* TARGET, double et, const char* FRAME,
		   const char* OBSERVER, const char* ABCORR, DACE::AlgebraicVector<double> &X){
/*! Read JPL ephemeris and return the TARGET body state relative to the OBSERVER body
	at the specified time et and in the given FRAME.
 * \param[in] TARGET a char array with the name of the target body.
 * \param[in] et ephemeris time (equivalent to TDB).
 * \param[in] FRAME a char array with the neme of the frame
 * \param[in] OBSERVER a char array with the name of the observer body
 * (the state computed is given with respect to this body).
 * \param[in] ABCORR a char array with the name of the aberration correction to be applied
    to the state of the target body to account for one-way light time and stellar
    aberration (see spkezr_c function for more details).
 * \param[out] X an AlgebraicVector containing the computed state (output).
 */
	double lt;
	spkezr_c ( TARGET, et, FRAME, ABCORR, OBSERVER, X.data(), &lt );	// XXX: data() is C++11. may not be supported!
}

void pleph(const int TARGET, double et, const char* FRAME,
		   const int OBSERVER, const char* ABCORR, DACE::AlgebraicVector<double> &X){
/*! Read JPL ephemeris and return the TARGET body state relative to the OBSERVER body
	at the specified time et and in the given FRAME.
 * \param[in] TARGET an int with the id code of the target body.
 * \param[in] et ephemeris time (see str2et_c function to convert epoch into et).
 * \param[in] FRAME a char array with the neme of the frame
 * \param[in] OBSERVER an int with the id code of the observer body.
 * (the state computed is given with respect to this body).
 * \param[in] ABCORR a char array with the name of the aberration correction to be applied
    to the state of the target body to account for one-way light time and stellar
    aberration (see spkez_c function for more details).
 * \param[out] X AlgebraicVector containing the computed state (output).
 */
	double lt;
	spkez_c ( TARGET, et, FRAME, ABCORR, OBSERVER, X.data(), &lt );	// XXX: data() is C++11. may not be supported!
}

// static const double ScaleHeight[] = { // follows the order of FrameOrigin in Reference.h
// 	0.0, 0.0, 8500.0, 0.0, 11100.0, 0.0,
// 	0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

// static const double SurfaceDensity[] = { // follows the order of FrameOrigin in Reference.h
// 	0.0, 0.0, 1.217, 0.0, 0.020, 0.0,
// 	0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

template<class T>
T density(const T &rho0, const double &H, const T &h) {
/*! Compute the atmosphere density at the given height above the
 *	body surface. This function is implemented by using an
 * 	exponential model of the form:
 *
 *	rho = rho0 * exp(-h/H)
 *
 *	where rho0 and H are retrieved from here:
 *	http://nssdc.gsfc.nasa.gov/planetary/planetfact.html
 *
 *  \param rho0 a scalar value of type T for the density at the body surface in kg(m^3).
 *  \param H a double for the atmosphere scale height in meters.
 *  \param h a scalar value of type T for the height above the body surface in meters.
 *  \return rho the T value of density in kg/m^3.
 */
	using std::exp;
	//using DACE::exp;

 //    if( (body != EARTH) && (body != MARS) )
 //        throw std::runtime_error("SF::density: Body density not available.");

	// const T rho0 = SurfaceDensity[body];
	// const double H 	= ScaleHeight[body];

// 	T rho = rho0*exp(-1.0*h/H);

	return rho0*exp(-1.0*h/H);
}

template <class T>
T shadow_conical_simple(const DACE::AlgebraicVector<T> &rr_sat, const DACE::AlgebraicVector<double> &rr_sun,  const DACE::AlgebraicVector<double> &rr_B, const double &Radius_B)
{
  using  std::acos;
 // using DACE::acos;
  using  std::asin;
  //using DACE::asin;
  //using DACE::sqr;
  using  std::sqrt;
  //using DACE::sqrt;

  T Q = 0.0;// Biconical/Dual-cone model

  // Get sun radius
  SpiceInt nn; double radii[3];
  bodvrd_c( "SUN", "RADII", 3, &nn, radii );
  double Radius_Sun = (radii[0]+radii[1]+radii[2])/3; // [km] Solar radius

  // Vector sat-occulting body
  DACE::AlgebraicVector<T> s_v = rr_sat - rr_B;
  T s = vnorm(s_v);

  // Vector sun oculting body
  DACE::AlgebraicVector<double> s_sun_v = rr_sun - rr_B;
  double s_sun = vnorm(s_sun_v);

  T beta = acos(dot(s_sun_v,s_v)/(s*s_sun));

  if (DACE::cons(beta)>kPI/2.)
  {
    // Vector sun-sat
    DACE::AlgebraicVector<T> u_v = rr_sun - rr_sat;
    T u = vnorm(u_v);

    // Apparent diameter of the Sun from satellite position
    T a_sun = asin(Radius_Sun/u);

    // Apparent diameter of the occulting body from satellite position
    T a_B = asin(Radius_B/s);

    // Apparent separation of bodies centres
    T c = acos(-dot(s_v,u_v)/(s*u));

    // Shadow condition
    if ((DACE::cons(c)<DACE::cons(a_B)-DACE::cons(a_sun))&&(DACE::cons(a_sun)<DACE::cons(a_B)))
    {
      // Total eclipse
      Q=1.;
    }
    else if ((DACE::cons(c)<DACE::cons(a_sun)-DACE::cons(a_B))&&(DACE::cons(a_sun)>DACE::cons(a_B)))
    {
        // Annular eclipse condition
        Q = DACE::sqr(a_B/a_sun);
    }
    else if ((DACE::cons(c)<DACE::cons(a_sun)+DACE::cons(a_B))&&(DACE::cons(c)>fabs(DACE::cons(a_sun)-DACE::cons(a_B))))
    {
        // Partial eclipse
        T x = (DACE::sqr(c)+DACE::sqr(a_sun)-DACE::sqr(a_B))/(2*c);
        T y = sqrt(DACE::sqr(a_sun)-DACE::sqr(x));

        T A_occ = DACE::sqr(a_sun)*acos(x/a_sun)+DACE::sqr(a_B)*acos((c-x)/a_B)-c*y;

        Q = A_occ/(kPI*DACE::sqr(a_sun));
    }
  }
  return Q;
}

template <class T>
T shadow(const DACE::AlgebraicVector<T> &rr_sat, const DACE::AlgebraicVector<double> &rr_sun,  const DACE::AlgebraicVector<double> &rr_B, const double &Radius_B, const int model_flag=0)
{
  /*! Computes shadow function for occulting body B using shadow model MODEL
   * Computes the shadow function Q for the Sun and the occultating body B
   * Inputs:
   *  \param rr_sat     satellite position in an inertial reference frame
   *  \param rr_sun     sun position in an inertial reference frame
   *  \param rr_B       occultating body position in an inertial reference frame
   *  \param Radius_B   diameter of occultating body (consinstent with other units)
   *  \param model_flag shadow model (optional)
   *                    0 -> Q=0 always (no shadow, default)
   *                    1 -> cylindrical model (Q=0 OR Q=1)
   *                    2 -> biconical model
   *
   * \return the shadow function Q, which represents the percentage of the Sun covered by the occulting body.
   *
   */

  using  std::acos;
 // using DACE::acos;
  using  std::asin;
  //using DACE::asin;
  //using DACE::sqr;
  using  std::sqrt;
  //using DACE::sqrt;

  T Q = 0.0;

  switch (model_flag)
  {
    case(0):
    default:
      break;

    case(1):
    {
      // Cylindrical model
      DACE::AlgebraicVector<double> rr_sat_cons(3);
      for (unsigned int i=0; i<3; i++)
        rr_sat_cons[i] = DACE::cons( (rr_sat[i]) );

      double r_sat_mod = vnorm(rr_sat_cons-rr_B);
      double r_sun_mod = vnorm(rr_sun-rr_B);

      double beta = acos(dot(rr_sat_cons - rr_B,rr_sun - rr_B)/(r_sat_mod*r_sun_mod));

      // Total eclipse check
      if ((beta>kPI/2.)&&(r_sat_mod*std::sin(beta)<Radius_B))
        return 1.0;

      break;
    }

    case(2):
    {
      // Biconical/Dual-cone model (Montenbruck/Gill, Satellite Orbits, pag 80-83)

      // Get sun radius
      SpiceInt nn; double radii[3];
      bodvrd_c( "SUN", "RADII", 3, &nn, radii );
      double Radius_Sun = (radii[0]+radii[1]+radii[2])/3; // [km] Solar radius

      // Vector sat-occulting body
      DACE::AlgebraicVector<T> s_v = rr_sat - rr_B;
      T s = vnorm(s_v);

      // Vector sun-sat
      DACE::AlgebraicVector<T> u_v = rr_sun - rr_sat ;
      T u = vnorm(u_v);

      // Apparent diameter of the Sun from satellite position
      T a_sun = asin(Radius_Sun/u);

      // Apparent diameter of the occulting body from satellite position
      T a_B = asin(Radius_B/s);

      // Apparent separation of bodies centres
      T c = acos(-dot(s_v,u_v)/(s*u));

      if (DACE::cons(c) <= std::abs(DACE::cons(a_B) - DACE::cons(a_sun)))
      {
        if(DACE::cons(a_B) >= DACE::cons(a_sun))
        {
          // Total eclipse
          return 1.0;
        }
        else
        {
          // Annular eclipse
          return DACE::sqr(a_B/a_sun);
        }
      }
      else if(DACE::cons(c) < std::abs(DACE::cons(a_B) + DACE::cons(a_sun)))
      {
        // Partial eclipse
        T x = (DACE::sqr(c)+DACE::sqr(a_sun)-DACE::sqr(a_B))/(2*c);
        T y = sqrt(DACE::sqr(a_sun)-DACE::sqr(x));

        T A_occ = DACE::sqr(a_sun)*acos(x/a_sun)+DACE::sqr(a_B)*acos((c-x)/a_B)-c*y; // Angles might go out of domain here...
        return A_occ/(kPI*DACE::sqr(a_sun));
      }
      break;
    }
  }
  return Q;
}

template<class T>
void position2lla( const DACE::AlgebraicVector<T> &rr_fixed, const Body &body, T &latgd, T &lon, T &alt)
{
  USEDACE_TRIGON
  USEDACE_ALG

  double small_ = 1e-12;
  T magr = vnorm(rr_fixed);

  T c=0.;
  SpiceInt ii;

  // Compute equatorial radius
  double radii[3];
  bodvrd_c( BodyLabel[body], "RADII", 3, &ii, radii );
  double Rearth = radii[0];

  // Compute first eccentricity squared
  double eesqrd = 1. - std::pow(radii[2]/radii[0],2);

  // Compute longitude
  T temp = sqrt( rr_fixed[0]*rr_fixed[0] + rr_fixed[1]*rr_fixed[1] );

  if (fabs(DACE::cons(temp))<small_)
    lon = kPI/2.0*( rr_fixed[2]/fabs( DACE::cons(rr_fixed[2]) ) );
  else
    lon = atan2_mod( rr_fixed[1], rr_fixed[0] );
  if ( fabs(DACE::cons(lon)) >= kPI )
  {
    if ( DACE::cons(lon)< 0.0)
      lon = 2.0*kPI + lon;
    else
      lon = lon - 2.0*kPI;
  }

  // Compute geodetic latitude
  int i = 0;
  latgd = asin(rr_fixed[2]/magr);
  double delta = 10.0;

  while ( ( fabs( delta )>=small_ )&&( i<10 ) )
  {

    delta = DACE::cons(latgd);

    T sintemp = sin(latgd);
    c = Rearth/sqrt(1. - eesqrd*sintemp*sintemp);
    latgd = atan2_mod((rr_fixed[2]+c*eesqrd*sintemp),temp);

    delta -= DACE::cons(latgd);

    i++;

  }

  // Calculate height
  if (kPI*0.5 - fabs(DACE::cons(latgd)) > kPI/180.0 )
  {
    alt = (temp/cos(latgd)) - c;
  }
  else
  {
    T s = c* (1. - eesqrd);
    alt = rr_fixed[2]/sin(latgd) - s;
  }

  return;
}
// struct planetvar {
// int ID;              // in agreement with spice enumeration (399 = Earth, etc.)
// double ref;
// double H;
// };

// const planetvar plrho[] = {
// {399, 1.217, 8500},  // [kg/m^3, m] see http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
// {499, 0.020, 11100}  // [kg/m^3, m] see http://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
// };

// template<class T>
// T density(const unsigned int id, const T h) {
// /* ! Compute the atmosphere density of the specified body at the given height
//  *	above the planet surface. This function is implemented by using an
//  * 	exponential model of the form:
//  *
//  *	rho = rho0 * exp(-h/H)
//  *
//  *	where rho0 and H are retrieved from here:
//  *	http://nssdc.gsfc.nasa.gov/planetary/planetfact.html
//  *
//  *	\param id Body ID
//  *	\param h height above the planet surface in meters.
//  *  \return rho the T value of density in kg/m^3.
//  */
// 	using std::exp;
// 	using DACE::exp;

//     const int length = sizeof(plrho)/sizeof(planetvar);
// 	int i = 0;
//     for( i=length-1; (i>0)&&(plrho[i].ID != id); i-- );

// 	double rho0 = plrho[i].ref;
// 	double H = plrho[i].H;

// 	T rho = rho0*exp(-1.0*h/H);

// 	return rho;
// }


/*! \defgroup CR3BP CR3BP routines
 * @brief Routines tailored for Circular Restricted 3-Body problem (CR3BP).
 *
 * This module groups all auxiliary routines tailored for the Circular Restricted 3-Body problem.
 * @{
 */

double SynodicMassParameter(const Body &primary, const Body& secondary)
{
  /*! Compute the mass parameter for the restricted three body problem
   * @param[in] primary   primary body ID
   * @param[in] secondary secondary body ID
   * \return the (adimensional) mass parameter
   */

  double GM_1, GM_2;
  SpiceInt n;

  // Retrieve variables from SPICE
  bodvrd_c( BodyLabel[primary], "GM", 1, &n, &GM_1 );
  bodvrd_c( BodyLabel[secondary], "GM", 1, &n, &GM_2 );

  return GM_2/(GM_1+GM_2);

}

void SynodicScaleFactors( const Body &primary, const Body& secondary, const double &et, double &DU, double &TU )
{
  /*! Compute the length and time scale factors for the restricted three body problem
   * @param[in] primary   primary body ID
   * @param[in] secondary secondary body ID
   * @param[in] et        ephemeris time at which scale factors are desired
   * @param[out] DU length scale factor
   * @param[out] TU time scale factor
   * @sa RelativeState
   */

  double GM;
  SpiceInt n;
  bodvrd_c( BodyLabel[primary], "GM", 1, &n, &GM );

  DACE::AlgebraicVector<double> LL = RelativeState(BodyLabel[secondary], BodyLabel[primary], et);

  double elts[8];
  oscelt_c(LL.data(),et,GM,elts);

  DU = elts[0]/(1-elts[1]); // Distance Unit, semi-major axis

  TU = std::sqrt(std::pow(DU,3.)/GM); // Time Unit

}

DACE::AlgebraicVector<double> CR3BPgetSecondaryKeplerianElements( const Body &primary, const Body& secondary, const double &et )
{
  /*! Get the keplerian elements of the secondary body (for CR3BP)
   * The eccentricity of the secondary body is forced to be zero.
   * @param[in] primary   primary body ID
   * @param[in] secondary secondary body ID
   * @param[in] et        ephemeris time at which keplerian elements are desired
   * \return an algebraic vector containing the keplerian elements of the secondary wrt primary body
   * @sa RelativeState
   */

  // Get gravitational parameter from spice kernels
  double GM;
  SpiceInt n;
  bodvrd_c( BodyLabel[primary], "GM", 1, &n, &GM );

  // Define primary inertial reference frame
  Frame Fprimary(J2000_EQT,primary);
  if (primary==SUN)
    Fprimary.Name = J2000_ECL;

  // Obtain position of the seconday wrt to the primary
  DACE::AlgebraicVector<double> LL = RelativeState(BodyLabel[secondary], Fprimary, et);

  // Compute osculating orbital elements
  double elts[8];
  oscelt_c(LL.data(),et,GM,elts);

  // Set orbital elements in AlgebraicVector
  DACE::AlgebraicVector<double> tmp(6);
  tmp[0] = elts[0]/(1-elts[1]); // Semi-major axis
  tmp[1] = 0.0;                 // Eccentricity (circular orbit)
  tmp[2] = elts[2];             // Inclination
  tmp[3] = elts[3];             // Right ascension of the ascending node
  tmp[4] = elts[4];             // Argument of the pericentre

  double E = atan2_mod( std::sqrt( (1-elts[1]*elts[1]) )*std::sin(elts[5])/(1+elts[1]*std::cos(elts[5])),
		    (elts[1]+std::cos(elts[5]))/( 1+elts[1]*std::cos(elts[5]) ) );
  if (E<0)
    E += 2*kPI;
  tmp[5] = E - elts[1]*std::sin(E);  // Mean/true anomaly (circular orbit)

  return tmp;

}

/*! @} */ // End of CR3BP module

DACE::AlgebraicVector<double> SPICEpropagate( const DACE::AlgebraicVector<double> &KE, const Body &attractor, const double &et0, const double &etf )
{
  /*! Propagate keplerian elements using SPICE routine conics_c.
   * @param[in] KE        keplerian elements vector
   * @param[in] attractor attractor body ID
   * @param[in] et0       ephemeris time corresponding to initial epoch
   * @param[in] etf       ephemeris time corresponding to final epoch
   * \return an algebraic vector containing the keplerian elements of the secondary wrt primary body
   */

  // Obtain gravitational parameter of primary
  double GM;
  SpiceInt n;
  bodvrd_c( BodyLabel[attractor], "GM", 1, &n, &GM );

  double elts[8];
  elts[0] = KE[0];  // Semi-major axis == pericentre (circular)
  elts[1] = KE[1];  // Eccentricity
  elts[2] = KE[2];  // Inclination
  elts[3] = KE[3];  // Right ascension of the ascending node
  elts[4] = KE[4];  // Argument of the pericentre
  elts[5] = KE[5];  // Mean anomaly
  elts[6] = et0;    // Epoch of orbital elements
  elts[7] = GM;     // Gravitational parameter

  double state[6];
  conics_c(elts,etf,state);

  DACE::AlgebraicVector<double> tmp(6);
  tmp[0] = state[0];
  tmp[1] = state[1];
  tmp[2] = state[2];
  tmp[3] = state[3];
  tmp[4] = state[4];
  tmp[5] = state[5];

  return tmp;

}

double KeplerEquation(double ecc, double M, double tol=1e-12)
{
  /*! Solve Kepler's equation for eccentric anomaly.
   * @param[in] ecc orbit eccentricity
   * @param[in] M   mean anomaly [rad]
   * @param[in] tol tolerance to stop iterations
   * @return the eccentric anomaly in radians.
   *
   * Kepler equation is here defined as \f$E-e\sin{E}=M\f$ and is solved using a Newton's method.
   */
  int maxiter=30;

  // Compute mean anomaly
  M = fmod(M, 2.*kPI);

  // Initialize eccentric anomaly
  double E = M;

  // Write Kepler implicit equation
  double F = E - ecc*std::sin(E) - M;

  // Iterate with Newton's method
  int i=0;
  while((std::abs(F)>tol)&&(i<maxiter))
  {
    E = E - F/(1.0-ecc*std::cos(E));
    F = E - ecc*std::sin(E) - M;
    i++;
  }

  return E;
}

double KeplerEquation(double sma, double ecc, double deltaT, double GM, double tol=1e-12)
{
  /*! Solve Kepler's equation for eccentric anomaly.
   * @param[in] sma     semi major axis [km]
   * @param[in] ecc     orbit eccentricity
   * @param[in] deltaT  time since reference epoch
   * @param[in] GM      gravitational parameter [km^3/s^2]
   * @param[in] tol     tolerance to stop iterations
   * @return the eccentric anomaly in radians.
   *
   * Kepler equation is here defined as \f$E-e\sin{E}=\sqrt{\frac{GM}{a^3}}(t-t_0)\f$ and is solved using a Newton's method.
   */

  // Compute mean anomaly
  double M = (std::sqrt(GM/sma)/sma)*deltaT;

  return KeplerEquation(ecc,M,tol);

}

template<typename T>
T tofabn(T sigma, T alpha, T beta)
{
  /*! Function that evaluates the time of flight via Lagrange expression
   *
   */

  using std::sin;
  using std::sqrt;
//  using DACE::sin;
//  using DACE::sqrt;
  using std::log;
//  using DACE::log;
  using std::sinh;
//  using DACE::sinh;

   if (DACE::cons(sigma)>0)
      return sigma*sqrt(sigma)*((alpha-sin(alpha))-(beta-sin(beta)));
   else
      return -sigma*sqrt(-sigma)*((sinh(alpha)-alpha)-(sinh(beta)-beta));
}

template<typename T>
T x2tof(const T &x, const T &s, const T &c, int lw)
{

  using  std::asin;
//  using DACE::asin;
//  using DACE::acos;
  using  std::acos;

  using  std::sqrt;
//  using DACE::sqrt;
  using std::log;
//  using DACE::log;


  T alpha, beta;

  T a = (s/2.)/(1.-x*x);

  if (DACE::cons(x)<1.)
  {
    // Ellipse
    beta = 2*asin(sqrt((s-c)/2./a));
    if (lw==1)
      beta = -beta;

    alpha = 2.*acos(x);

  }
  else
  {
    alpha = 2.*log(x+sqrt(x*x-1.));
    T aux = sqrt((s-c)/(-2.*a));

    beta = 2.*log(aux+sqrt(aux*aux+1.));
    if (lw==1)
      beta = -beta;

  }

  return tofabn(a,alpha,beta);

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

DACE::AlgebraicVector<DACE::DA> llaMaps(DACE::AlgebraicVector<DACE::DA> & x, double et, double Lsc)
{
    // Compute rotation matrix from ECI to ECEF
    DACE::DA lat, lon, alt;
    double mu  = 398600.4418;
    double Tsc = Lsc/sqrt(mu/Lsc);
    SpiceDouble aux[6][6];
    DACE::AlgebraicMatrix<double> dcm(3,3);
    sxform_c("J2000", "EARTH_FIXED", et*Tsc, aux);
    for (unsigned int i = 0; i < 3; i ++){
        for (unsigned int j = 0; j < 3; j ++)
        dcm.at(i,j) = aux[i][j];
    }
    DACE::AlgebraicVector<DACE::DA> ps(3); ps[0] = x[0]*Lsc; ps[1] = x[1]*Lsc; ps[2] = x[2]*Lsc;
    DACE::AlgebraicVector<DACE::DA> pEcef(3); pEcef = dcm*ps;
    astro::position2lla(pEcef, astro::EARTH, lat, lon, alt);
    DACE::AlgebraicVector<DACE::DA> latlonalt(3); latlonalt[0] = lat; latlonalt[1] = lon; latlonalt[2] = alt/Lsc;
  return latlonalt;
}

DACE::AlgebraicVector<DACE::DA> smdGrad(DACE::AlgebraicMatrix<double> & P,DACE::AlgebraicVector<DACE::DA> & r)
{
double p11 = P.at(0,0); double p12 = P.at(0,1); double p13 = P.at(0,2); 
double p22 = P.at(1,1); double p23 = P.at(1,2);
double p33 = P.at(2,2);

DACE::DA x = r[0]; 
DACE::DA y = r[1]; 
DACE::DA z = r[2];
DACE::DA t2 = p12*p13;
DACE::DA t3 = p11*p23;
DACE::DA t4 = p12*p23;
DACE::DA t5 = p13*p22;
DACE::DA t6 = p13*p23;
DACE::DA t7 = p12*p33;
DACE::DA t8 = p12*p12;
DACE::DA t9 = p13*p13;
DACE::DA t10 = p23*p23;
DACE::DA t11 = p11*p22*p33;
DACE::DA t12 = -t3;
DACE::DA t13 = -t5;
DACE::DA t14 = -t7;
DACE::DA t15 = p23*t3;
DACE::DA t16 = p13*t5;
DACE::DA t17 = p12*t7;
DACE::DA t18 = p23*t2*2.0;
DACE::DA t20 = -t11;
DACE::DA t19 = -t18;
DACE::DA t21 = t2+t12;
DACE::DA t22 = t4+t13;
DACE::DA t23 = t6+t14;
DACE::DA t24 = t15+t16+t17+t19+t20;
DACE::DA t25 = 1.0/t24;
DACE::AlgebraicVector<DACE::DA> grad(3);
grad[0] = t25*(t10*x-t6*y+t7*y-t4*z+t5*z-p22*p33*x)+t25*x*(t10-p22*p33)-t23*t25*y-t22*t25*z;
grad[1] = -t25*(t6*x+t14*x-t9*y+t2*z+t12*z+p11*p33*y)+t25*y*(t9-p11*p33)-t23*t25*x-t21*t25*z;
grad[2] = -t25*(t4*x+t13*x+t2*y+t12*y-t8*z+p11*p22*z)+t25*z*(t8-p11*p22)-t22*t25*x-t21*t25*y;
grad = grad/sqrt(dot(grad,grad));
return grad;
}

DACE::AlgebraicVector<double> trustRegion(DACE::AlgebraicVector<DACE::DA> & x, int n)
{
  using namespace DACE;
  using namespace std;
  AlgebraicMatrix<DA> J(6,n);
  double a, b, normJ0;
  AlgebraicMatrix<double> J0(6,n);
  AlgebraicVector<double> nu(n);
  int j, k, i;
  normJ0  = 0.0;
  for (j = 0; j < n; j++){
    for (k = 0; k < 6; k++){
        J.at(k,j)  = x[k].deriv(j+1);
        normJ0 += cons(J.at(k,j))*cons(J.at(k,j));
    }
  }
  for (i = 0; i < n; i++){
    b = 0.0;
    for (j = 0; j < n; j++){
      for (k = 0; k < 6; k++){
        a         = cons(J.at(k,j).deriv(i+1));
        b        += a*a;
      }
    }
    nu[i] = sqrt(b/normJ0);
  }  
return nu;
}

DACE::AlgebraicVector<double> trustRegionScalar(DACE::DA & x, int n)
{
  using namespace DACE;
  using namespace std;
  AlgebraicVector<DA> J(n);
  double a, b, normJ0;
  AlgebraicVector<double> J0(n), nu(n);
  int j, i;
  normJ0  = 0.0;
  for (j = 0; j < n; j++){
        J[j]  = x.deriv(j+1);
        normJ0 += cons(J[j])*cons(J[j]);
  }
  for (i = 0; i < n; i++){
    b = 0.0;
    for (j = 0; j < n; j++){
        a         = cons(J[j].deriv(i+1));
        b        += a*a;
    }
    nu[i] = sqrt(b/normJ0);
  }  
return nu;
}

template<typename T> AlgebraicVector<T> cart2kep(const AlgebraicVector<T>& rv, const double mu)
{
    AlgebraicVector<T> kep(6);
    
    AlgebraicVector<T> rr(3), vv(3);
    for (int i = 0; i < 3; i++)
    {
        rr[i] = rv[i];
        vv[i] = rv[i + 3];
    }
    
    T r = rr.vnorm();
    T v = vv.vnorm();
    AlgebraicVector<T> h = cross(rr, vv);
    
    kep[0] = mu / (2.0 * (mu / r - pow(v, 2) / 2.0));
    
    T h1sqr = pow(h[0], 2);
    T h2sqr = pow(h[1], 2);
    
    T RAAN;
    if (cons(h1sqr + h2sqr) == 0.0)
    {
        RAAN = 0.0;
    }
    else
    {
        T sinOMEGA = h[0] / sqrt(h1sqr + h2sqr);
        T cosOMEGA = -1.0*h[1] / sqrt(h1sqr + h2sqr);
        if (cons(cosOMEGA) >= 0.0)
        {
            if (cons(sinOMEGA) >= 0.0)
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
            if (cons(sinOMEGA) >= 0.0)
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
    
    AlgebraicVector<T> ee = 1.0 / mu*(cross(vv, h)) - rr / vnorm(rr);
    T e = vnorm(ee);
    T i = atan2_mod(sqrt(h[0]*h[0]+h[1]*h[1]),h[2]);
    // T i = acos(h[2] / vnorm(h));

    kep[1] = e;
    kep[2] = i;
    kep[3] = RAAN;
    
    T omega;
    T theta;
    if (cons(e) <= 1.0e-8 && cons(i) < 1.0e-8)
    {
        e = 0.0;
        omega = atan2_mod(rr[1], rr[0]);
        theta = 0.0;
        kep[4] = omega;
        kep[5] = theta;
        return kep;
    }
    
    if (cons(e) <= 1.0e-8 && cons(i) >= 1.0e-8)
    {
        omega = 0;
        AlgebraicVector<T> P(3), Q(3), W(3);
        P[0] = cos(omega)*cos(RAAN) - sin(omega)*sin(i)*sin(RAAN);
        P[1] = -1.0*sin(omega)*cos(RAAN) - cos(omega)*cos(i)*sin(RAAN);
        P[2] = sin(RAAN)*sin(i);
        Q[0] = cos(omega)*sin(RAAN) + sin(omega)*cos(i)*cos(RAAN);
        Q[1] = -1.0*sin(omega)*sin(RAAN) + cos(omega)*cos(i)*cos(RAAN);
        Q[2] = -1.0*cos(RAAN)*sin(i);
        W[0] = sin(omega)*sin(i);
        W[1] = cos(omega)*sin(i);
        W[2] = cos(i);
        AlgebraicVector<T> rrt = P*rr[0] + Q*rr[1] + W*rr[2];
        theta = atan2_mod(rrt[1], rrt[0]);
        kep[4] = omega;
        kep[5] = theta;
        return kep;
    }
    
    T dotRxE = dot(rr, ee);
    T RxE = vnorm(rr)*vnorm(ee);
    if (abs(cons(dotRxE)) > abs(cons(RxE)) && abs(cons(dotRxE)) - abs(cons(RxE)) < abs(numeric_limits<double>::epsilon()*cons(dotRxE)))
    {
        dotRxE -= numeric_limits<double>::epsilon()*dotRxE;
    }
    theta = acos(dotRxE / RxE);
    
    if (cons(dot(rr, vv)) < 0.0)
    {
        theta = 2.0 * M_PI - theta;
    }
    
    if (cons(i) <= 1.0e-8 && cons(e) >= 1.0e-8)
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

    if (cons(coso) >= 0.0)
    {
        if (cons(sino) >= 0.0)
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
        if (cons(sino) >= 0.0)
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
    
    if (cons(omega) < 0.0)
    {
        omega = omega + 2.0 * M_PI;
    }
    //omega = real(omega);
    
    kep[4] = omega;
    kep[5] = theta;
    
    return kep;
}


template <typename T> AlgebraicVector<T> mee2cart(AlgebraicVector<T> & Mee, const double mu)
{
    T p       = Mee[0];
    T f       = Mee[1];
    T g       = Mee[2];
    T h       = Mee[3];
    T k       = Mee[4];
    T L       = Mee[5];

    T q           = 1.0 + f * cos(L) + g * sin(L);
    T r           = p / q;
    T alpha2      = h*h - k*k;
    T chi2        = h*h + k*k;
    T s2          = 1.0 + chi2;

    DACE::AlgebraicVector<T> res(6);
    res[0] = r / s2 * ( cos(L) + alpha2 * cos(L) + 2 * h * k * sin(L) );
    res[1] = r / s2 * ( sin(L) - alpha2 * sin(L) + 2 * h * k * cos(L) );
    res[2] = 2 * r / s2 * ( h * sin(L) - k * cos(L) );
    res[3] = - 1 / s2 * sqrt(mu/p) * ( sin(L) + alpha2 * sin(L) - 2 * h * k * cos(L) + g - 2 * f * h * k + alpha2 * g );
    res[4] = - 1 / s2 * sqrt(mu/p) * ( -cos(L) + alpha2 * cos(L) + 2 * h * k * sin(L) - f + 2 * g * h * k + alpha2 * f );
    res[5] = 2 / s2 * sqrt(mu/p) * ( h * cos(L) + k * sin(L) + f * h + g * k);
    
    return res;
}

template <typename T> AlgebraicVector<T> coe2mee(AlgebraicVector<T> & coe)
{

AlgebraicVector<T> mee(6); 
T a       = coe[0];
T ecc     = coe[1];
T in      = coe[2];
T Om      = coe[3];
T om      = coe[4];
T theta   = coe[5];
T dL;
double consL;

T p = a * ( 1 - ecc*ecc);
T f = ecc * cos(om + Om);
T g = ecc * sin(om + Om);
T h = tan(in/2) * cos(Om);
T k = tan(in/2) * sin(Om);
T L = Om + om + theta;

// Keep the anomaly between 0 and 2pi
if (abs(cons(L)) > 6.28318530718) {
  consL = std::fmod(cons(L), 6.28318530718);
  dL = L - cons(L);
  L = consL + dL;
}
if (cons(L) < 0.0) {
  L = L + 6.28318530718;
}
mee = {p, f, g, h, k, L};

return mee;
}

template <typename T> AlgebraicVector<T> cart2mee(AlgebraicVector<T> & rv, const double mu)
{
    double consMee;
    AlgebraicVector<T> coe(6), mee(6);
    coe = astro::cart2kep(rv, mu);

    mee = astro::coe2mee(coe);
    return mee;
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

template<typename T> T det3(AlgebraicMatrix<T> M)
{    
   // computes the determinant of a 3x3 matrix M
    T det = M.at(0, 0) * (M.at(1, 1) * M.at(2, 2) - M.at(2, 1) * M.at(1, 2)) -
            M.at(0, 1) * (M.at(1, 0) * M.at(2, 2) - M.at(1, 2) * M.at(2, 0)) +
            M.at(0, 2) * (M.at(1, 0) * M.at(2, 1) - M.at(1, 1) * M.at(2, 0));

    return det;
}

int factorial(int n)
{
    // single line to find factorial
    return (n==1 || n==0) ? 1: n * factorial(n - 1);
}

template<typename T> T ConstIPC(AlgebraicVector<T> r, AlgebraicMatrix<double> P, double R){

    // 3D Constant IPC
    const double pi = 3.141592653589793;
    AlgebraicMatrix<double> P_inv(3,3);
    double det = det3(P);
    P_inv = P.inv();

    T smd = r.dot(P_inv*r);
    T ipc = sqrt(2)*R*R*R/(3*sqrt(pi*det))*exp(-smd/2);
    
    return ipc;
}

template<typename T> T MaxIPC(AlgebraicVector<T> r, AlgebraicMatrix<double> P, double R){

    // 3D Constant IPC
    const double pi = 3.141592653589793;
    AlgebraicMatrix<double> P_inv(3,3);
    double det = det3(P);
    P_inv = P.inv();

    T smd = r.dot(P_inv*r);
    T ipc = sqrt(8)*R*R*R/(3*exp(1)*sqrt(pi*det)*smd);
    
    return ipc;
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

DA findTCA(const AlgebraicVector<DA> xrel, const int nvar){

  AlgebraicVector<DA> rr(3), vv(3), MAPD(nvar), MAPI(nvar), dx(nvar);

  //relative velocities
  rr[0] = xrel[0]; rr[1] = xrel[1]; rr[2] = xrel[2];
  vv[0] = xrel[3]; vv[1] = xrel[4]; vv[2] = xrel[5];

  // we want rr to be orthogonal to vv, so dot(rr,vv) = 0
  DA rvdot = dot(rr,vv);
  double rvdotcons = cons(rvdot);

  // MAPD = (dot(rr,vv), dxx) = f(dxx, dt)
  MAPD[0] = rvdot-rvdotcons;
  for (int i = 1; i < nvar ; i++) {
    MAPD[i] = DA(i);
  }

  // MAPI = (dxx, dt) = f(dot(rr,vv), dxx)
  MAPI = MAPD.invert();

  // we need to evaluate the map in -cons(dot(rr,vv)), dxx
  dx[0] = - rvdotcons + 0*DA(1);
  for (int i = 1; i < nvar ; i++){
      dx[i] = DA(i);}
  MAPD = MAPI.eval(dx);

  // dt is the last row of the MAPD
  DA tca = MAPD[nvar - 1];

  return tca;
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

}
#endif /* astro_AstroRoutines_H */
