  /*
 * AIDA.h
 *
 *  Created on: December 12, 2014
 *      Author: Dinamica Srl
 *  Modified on: March 2016
 *      Author: Mauro Massari (mauro.massari@polimi.it)
 */

#ifndef AIDA_H_
#define AIDA_H_

#include <memory>
#include "astro/AstroDynModel.h"
#include "atmos/Nrlmsise00.hpp"
#include "atmos/SpaceWeather.hpp"
#include "GravityFieldModel.h"
#include "geco/JsonParser.hpp"

namespace dynorb {

/*! AIDA class.
 * \note This class implements the AIDA model within the DACE framework
 *       (see Documentation for further details).
 */

template <class T> class AIDA : public astro::AstroDynModel<T>
{

    static const std::vector<std::string> IntParamLabel;
    static const std::vector<std::string> FixParamLabel;
    static const std::vector<std::string> UncParamLabel;

    GravityFieldModel gravity_model;
    std::unique_ptr<atmos::CNrlmsise00<T> >  p_atmosphere_model;

    int shadow_model[2];
    std::vector<double> radius_B;

    DACE::AlgebraicVector<T> acc_3rdbody( const DACE::AlgebraicVector<T> &rr_sat, const DACE::AlgebraicVector<double> &rr_3rdb, const double &GM_3rdb ) const;
    DACE::AlgebraicVector<T> acc_srp( const DACE::AlgebraicVector<T> &rr_sat, const DACE::AlgebraicVector<double> &rr_sun, const DACE::AlgebraicMatrix<double> &rr_B) const;
    DACE::AlgebraicVector<T> acc_gravfield( const DACE::AlgebraicVector<T> &rr_sat, const double &et) const;
    DACE::AlgebraicVector<T> acc_airdrag( const DACE::AlgebraicVector<T> &rr_sat, const DACE::AlgebraicVector<T> &vv_sat, const double &et) const;

    void allocate_NRLMSISE()
    {
        /*! Allocate atmosphere model.
         */

        // Set flags
        std::array<int,24> flags;
        for(unsigned i=0; i<flags.size(); i++)
            flags[i] =1;

        p_atmosphere_model = std::make_unique<atmos::CNrlmsise00<T> >(flags);
    }

protected:

    static const unsigned int m_maxvar;  /*!< Number of allowed uncertainties.      */
    static const unsigned int NIntParam; /*!< Number of integer parameters.         */
    static const unsigned int NFixParam; /*!< Number of non uncertain parameters.   */
    static const unsigned int NUncParam; /*!< Number of uncertain parameters.       */

public:
    /***********************************************************************************
    *     Static members
    ************************************************************************************/
    static unsigned int MaxVar()
    {
      /*! Get the number of uncertainties for the AIDA model.   */
      return m_maxvar;
    }

     void DAinit(const unsigned int ord, const unsigned int nvar = m_maxvar)
     {
       /*! DACE initialization with a number of variables equal to the
        * number of allowed uncertainties for the AIDA model.
        */
       DACE::DA::init(ord, nvar);
       astro::AstroDynModel<T>::AllocateUncParam();
       astro::AstroDynModel<T>::multipleObj = true;

       // Allocate atmosphere model (must be done after DA initialization)
       // Set flags
        std::array<int,24> flags;
        for(unsigned i=0; i<flags.size(); i++)
            flags[i] =1;

        p_atmosphere_model = std::make_unique<atmos::CNrlmsise00<T> >(flags);
    }

    /***********************************************************************************
    *     Constructors
    ************************************************************************************/
    /*! Create a new AIDA object with the given central body.
     * \param gravmodel_name The name of the gravitational model to be used.
     * \param gravmodel_order The order of the gravitational field model.
     * \param AIDA_flags An array with the flags for the perturbative models to be considered.
     */
    AIDA(const std::string& gravmodel_name, const unsigned int gravmodel_order, const int* AIDA_flags,
      const double minh=0.0 );

    /*! Create a new AIDA object using Json node and a parser.
     * \param input The Json node defining the model
     * \param parser The parser to read the node
     */
    AIDA(const Json::Value& input, geco::CJsonParser parser);

    ~AIDA() = default;

    /***********************************************************************************
    *     Dynamical model evaluation
    ************************************************************************************/
    void evaluation(const double t, const DACE::AlgebraicVector<T> &X, DACE::AlgebraicVector<T> &Xdot, bool &flag) const;
    DACE::AlgebraicVector<T> evaluation(const double t, const DACE::AlgebraicVector<T> &X, bool &flag) const;

};


/***********************************************************************************
*     Setup static member variables
************************************************************************************/

// Integer parameters
template<class T>
const unsigned int AIDA<T>::NIntParam = 5;

const char* AIDA_IntParamLabel[] = {"ATTRACTOR","GRAVORDER","ATMOSPHERE","SRP","THIRDBODY"}; /*!< Labels associated to integer parameters. */

template<class T>
const std::vector<std::string> AIDA<T>::IntParamLabel(AIDA_IntParamLabel, AIDA_IntParamLabel+NIntParam);


// Non Uncertain/Fixed parameters
template<class T>
const unsigned int AIDA<T>::NFixParam = 7;

const char* AIDA_FixParamLabel[] = {"GM","RADIUS","RADIUS_MOON","GM_MOON","GM_SUN","P_SUN", "MIN_HEIGHT"}; /*!< Labels associated to fixed (non-uncertain) parameters. */
template<class T>
const std::vector<std::string> AIDA<T>::FixParamLabel(AIDA_FixParamLabel, AIDA_FixParamLabel+NFixParam);


// Uncertain parameters
 template<class T>
 const unsigned int AIDA<T>::NUncParam = 2;
 const char* AIDA_UncParamLabel[] = {"Bfactor","SRPC"}; /*!< Labels associated to uncertain parameters. */

 template<class T>
 const std::vector<std::string> AIDA<T>::UncParamLabel(AIDA_UncParamLabel,
        AIDA_UncParamLabel+NUncParam);


// Maximum number of uncertain variables
template<class T>
const unsigned int AIDA<T>::m_maxvar = 6 + AIDA<T>::NUncParam;

/***********************************************************************************
*     Constructors
************************************************************************************/
template<class T>
AIDA<T>::AIDA(const std::string& gravmodel_name, const unsigned int gravmodel_order, const int* AIDA_flags, const double minh)
: astro::AstroDynModel<T>("AIDA", 6, astro::CARTESIAN, astro::J2000_EQT, astro::EARTH, astro::KM, astro::SEC, astro::RAD),
  gravity_model(gravmodel_name,gravmodel_order)
{
    // Allocate atmosphere model
    allocate_NRLMSISE();
    // Allocate Integer parameters
    astro::AstroDynModel<T>::setIntParamLabels(NIntParam, IntParamLabel); // (do not remove)
    this->setIntParam("ATTRACTOR", astro::EARTH);
    this->setIntParam("GRAVORDER", gravmodel_order);
    // Set AIDA perturbation flags
    this->setIntParam("ATMOSPHERE", AIDA_flags[0]);
    this->setIntParam("SRP",        AIDA_flags[1]);
    this->setIntParam("THIRDBODY",  AIDA_flags[2]);
    // Allocate Fixed parameters
    astro::AstroDynModel<T>::setFixParamLabels(NFixParam, FixParamLabel); // (do not remove)
    // Set attractor gravitational constant
    this->setFixParam("GM",gravity_model.getGM());
    // Set Moon's gravitational parameter
    SpiceInt n;
    double GM;
    bodvrd_c( "MOON", "GM", 1, &n, &GM ); // [km^3/s^2]
    this->setFixParam("GM_MOON",GM);
    // Set Sun's gravitational parameter
    bodvrd_c( "SUN", "GM", 1, &n, &GM ); // [km^3/s^2]
    this->setFixParam("GM_SUN",GM);
    // return body radius from spice
    this->setFixParam("RADIUS",gravity_model.getRadius());
    // Set Moon radius
    double radii[3];
    bodvrd_c( "MOON", "RADII", 3, &n, radii );
    this->setFixParam("RADIUS_MOON",radii[0]);
    double Psun = LSUN/(4.*kPI*KM_AU*KM_AU*clight_c());
    this->setFixParam("P_SUN",Psun);
    // Set minimum height (default = 0 km)
    this->setFixParam("MIN_HEIGHT", minh);
    // Initialize flags for shadow models
    switch (this->getIntParam("SRP"))
    {
        default:
            throw std::runtime_error("dynorb::AIDA<T>::AIDA: srp mode must be an integer in the interval [0,6]");
        case 0: // No SRP
        case 1: // SRP, no shadow
            break;
        case 2: // SRP, Earth cylindrical shadow
            shadow_model[0] = 1;
            radius_B.push_back(this->getFixParam("RADIUS"));
            break;
        case 3: //  SRP, Earth biconical shadow
            shadow_model[0] = 2;
            radius_B.push_back(this->getFixParam("RADIUS"));
            break;
        case 4: // SRP, Earth and Moon cylindrical shadow
            shadow_model[0] = 1; shadow_model[1] = 1;
            radius_B.push_back(this->getFixParam("RADIUS"));
            radius_B.push_back(this->getFixParam("RADIUS_MOON"));
            break;
        case 5: // SRP, Earth biconical and Moon cylindrical shadow
            shadow_model[0] = 2; shadow_model[1] = 1;
            radius_B.push_back(this->getFixParam("RADIUS"));
            radius_B.push_back(this->getFixParam("RADIUS_MOON"));
            break;
        case 6: // SRP, Earth and Moon biconical shadow
            shadow_model[0] = 2; shadow_model[1] = 2;
            radius_B.push_back(this->getFixParam("RADIUS"));
            radius_B.push_back(this->getFixParam("RADIUS_MOON"));
    }
    // Check atmosphere flag intervals
    if(this->getIntParam("ATMOSPHERE")>2 or this->getIntParam("ATMOSPHERE")<0)
    {
        throw std::invalid_argument("dynorb::AIDA<T>::AIDA: atmosphere mode must be 0 (none), 1 (non-rotating), or 2 (co-rotating)!");
    }
    // Check third body flag interval
    if(this->getIntParam("THIRDBODY")>2 or this->getIntParam("THIRDBODY")<0)
    {
        throw std::invalid_argument("dynorb::AIDA<T>::AIDA: third body mode must be 0 (none), 1 (Moon), or 2 (Moon and Sun)!");
    }
    // Set Uncertain parameters
    astro::AstroDynModel<T>::setUncParamLabels(NUncParam, UncParamLabel); // (do not remove)
    // Allocate uncertain parameters only if DA is previously initialized. If not, DAinit must be used after calling the constructor.
    if (DACE::DA::isInitialized()){
       if (astro::AstroDynModel<T>::multipleObj == false){ // Allocate uncertain parameters only if no other dynamical models have been previously initialized. If not, DAinit must be used after calling the constructor.
           astro::AstroDynModel<T>::AllocateUncParam();
       }
    }
}


template<class T>
AIDA<T>::AIDA(const Json::Value& input, geco::CJsonParser parser)
: astro::AstroDynModel<T>("AIDA", 6, astro::CARTESIAN, astro::J2000_EQT, astro::EARTH, astro::KM, astro::SEC, astro::RAD),
  gravity_model(input["gravityField"], parser)
{
    // Allocate atmosphere model
    allocate_NRLMSISE();

    // Allocate Integer parameters
    astro::AstroDynModel<T>::setIntParamLabels(NIntParam, IntParamLabel); // (do not remove)
    this->setIntParam("ATTRACTOR", astro::EARTH);
    this->setIntParam("GRAVORDER", gravity_model.getDegree());

    // Allocate Fixed parameters
    astro::AstroDynModel<T>::setFixParamLabels(NFixParam, FixParamLabel); // (do not remove)
    // Set attractor gravitational constant
    this->setFixParam("GM",gravity_model.getGM());
    // Set Moon's gravitational parameter
    SpiceInt n;
    double GM;
    bodvrd_c( "MOON", "GM", 1, &n, &GM ); // [km^3/s^2]
    this->setFixParam("GM_MOON",GM);
    // Set Sun's gravitational parameter
    bodvrd_c( "SUN", "GM", 1, &n, &GM ); // [km^3/s^2]
    this->setFixParam("GM_SUN",GM);
    // return body radius from spice
    this->setFixParam("RADIUS",gravity_model.getRadius());
    // Set Moon radius
    double radii[3];
    bodvrd_c( "MOON", "RADII", 3, &n, radii );
    this->setFixParam("RADIUS_MOON",radii[0]);

    // Solar radiation pressure
    double Psun = LSUN/(4.*kPI*KM_AU*KM_AU*clight_c());
    this->setFixParam("P_SUN",Psun);
    // Set minimum height (default = 0 km)
    this->setFixParam("MIN_HEIGHT", parser.Read<double>(input, "minHeight", 0.0));

    // Set AIDA third body perturbation flags
    this->setIntParam("THIRDBODY",  parser.Read<int>(input,"thirdBody",2)); // Default: Moon and Sun
    if(this->getIntParam("THIRDBODY")>2 or this->getIntParam("THIRDBODY")<0)
    {
        throw std::invalid_argument("dynorb::AIDA<T>::AIDA: third body mode must be 0 (none), 1 (Moon), or 2 (Moon and Sun)!");
    }

    // Set AIDA atmospheric perturbation flag and parameters
    if(not input.isMember("drag"))
    {
        throw std::invalid_argument("dynorb::AIDA<T>::AIDA: No drag key provided!");
    }

    const Json::Value atmInputs = input["drag"];
    this->setIntParam("ATMOSPHERE", parser.Read<int>(atmInputs, "mode", 2)); //< Default: co-rotating
    if(this->getIntParam("ATMOSPHERE")>2 or this->getIntParam("ATMOSPHERE")<0)
    {
        throw std::invalid_argument("dynorb::AIDA<T>::AIDA: atmosphere mode must be 0 (none), 1 (non-rotating), or 2 (co-rotating)!");
    }

    const Json::Value srpInputs = input["srp"];
    this->setIntParam("SRP", parser.Read<int>(srpInputs, "mode", 3));

    // Initialize flags for shadow models
    shadow_model[0] = 0; shadow_model[1]=0;
    switch (this->getIntParam("SRP"))
    {
        default:
            throw std::runtime_error("dynorb::AIDA<T>::AIDA: srp mode must be an integer in the interval [0,6]");
        case 0: // No SRP
        case 1: // SRP, no shadow
            break;
        case 2: // SRP, Earth cylindrical shadow
            shadow_model[0] = 1;
            radius_B.push_back(this->getFixParam("RADIUS"));
            break;
        case 3: //  SRP, Earth biconical shadow
            shadow_model[0] = 2;
            radius_B.push_back(this->getFixParam("RADIUS"));
            break;
        case 4: // SRP, Earth and Moon cylindrical shadow
            shadow_model[0] = 1; shadow_model[1] = 1;
            radius_B.push_back(this->getFixParam("RADIUS"));
            radius_B.push_back(this->getFixParam("RADIUS_MOON"));
            break;
        case 5: // SRP, Earth biconical and Moon cylindrical shadow
            shadow_model[0] = 2; shadow_model[1] = 1;
            radius_B.push_back(this->getFixParam("RADIUS"));
            radius_B.push_back(this->getFixParam("RADIUS_MOON"));
            break;
        case 6: // SRP, Earth and Moon biconical shadow
            shadow_model[0] = 2; shadow_model[1] = 2;
            radius_B.push_back(this->getFixParam("RADIUS"));
            radius_B.push_back(this->getFixParam("RADIUS_MOON"));
    }

    // Set Uncertain parameters
    astro::AstroDynModel<T>::setUncParamLabels(NUncParam, UncParamLabel); // (do not remove)
    // Allocate uncertain parameters only if DA is previously initialized. If not, DAinit must be used after calling the constructor.
    if (DACE::DA::isInitialized()){
       if (astro::AstroDynModel<T>::multipleObj == false){ // Allocate uncertain parameters only if no other dynamical models have been previously initialized. If not, DAinit must be used after calling the constructor.
           astro::AstroDynModel<T>::AllocateUncParam();
       }
    }
}

/***********************************************************************************
*     Routines for the computation of perturbations
************************************************************************************/
template<class T>
DACE::AlgebraicVector<T> AIDA<T>::acc_3rdbody(const DACE::AlgebraicVector<T> &rr_sat, const DACE::AlgebraicVector<double> &rr_3rdb, const double &GM_3rdb) const
{
  /*! Compute the acceleration due to the third body gravitational attraction.
   *  \param rr_sat An DACE::AlgebraicVector<T> for the position of the s/c.
   *  \param rr_3rdb An DACE::AlgebraicVector<T> for the position of the third body.
   *  \param GM_3rdb The gravitational constant of the third body.
   */
  //using ASTRO::sqrt;
  using  std::sqrt;

  // compute position of the third body wrt the s/c
  DACE::AlgebraicVector<T> rr_sat_3rdb = rr_3rdb - rr_sat;
  T ru_sat_3rdb = DACE::sqr(rr_sat_3rdb[0]) + DACE::sqr(rr_sat_3rdb[1]) + DACE::sqr(rr_sat_3rdb[2]);
  T iru_sat_3rb = 1./sqrt(ru_sat_3rdb);
  ru_sat_3rdb = 1./ru_sat_3rdb;

  // compute position of the third body
  T ru_3rdb = DACE::sqr(rr_3rdb[0]) + DACE::sqr(rr_3rdb[1]) + DACE::sqr(rr_3rdb[2]);
  T iru_3rdb = 1./sqrt(ru_3rdb);
  ru_3rdb = 1./ru_3rdb;

  // compute accelerations
  return (GM_3rdb*ru_sat_3rdb) * (rr_sat_3rdb*iru_sat_3rb) - (GM_3rdb*ru_3rdb) * (rr_3rdb*iru_3rdb) ;

}

template<class T>
DACE::AlgebraicVector<T> AIDA<T>::acc_srp( const DACE::AlgebraicVector<T> &rr_sat, const DACE::AlgebraicVector<double> &rr_sun, const DACE::AlgebraicMatrix<double> &rr_B) const
{
  /*! Compute the acceleration due to the solar radiation pressure.
   *  \param rr_sat An DACE::AlgebraicVector<T> for the position of the s/c.
   *  \param rr_sun An DACE::AlgebraicVector<T> for the position of the Sun.
   *  \param rr_B An DACE::AlgebraicVector<T> for the position of the occulting body.
   */
  //using ASTRO::sqr;
  //using ASTRO::sqrt;
  using  std::sqrt;

  T nu = 1.;

  // determine the shadow factor for each occulting body
  for (unsigned int i=0; i<radius_B.size(); i++)
  {
    // The illumination factor is approximated by 1.0-sum(Q)
    // This does not account for intersection of occulting bodies, it might overestimate the eclipse condition.
    nu -= astro::shadow(rr_sat, rr_sun, rr_B.getcol(i), radius_B[i], shadow_model[i]);
  }

  // compute the position of sun wrt s/c
  DACE::AlgebraicVector<T> rr_sun_sat = rr_sat - rr_sun;
  T ru_sun_sat = DACE::sqr(rr_sun_sat[0]) + DACE::sqr(rr_sun_sat[1]) + DACE::sqr(rr_sun_sat[2]);
  T iru_sun_sat = 1./sqrt(ru_sun_sat);
  ru_sun_sat = 1./ru_sun_sat;

  // Retrieve parameters
  T CrAoMSRP = this->getUncParam("SRPC");
  double Psun = this->getFixParam("P_SUN");

  // compute accelerations
  return (CrAoMSRP * nu * (Psun * KM_AU * KM_AU * ru_sun_sat) * (rr_sun_sat * iru_sun_sat)) * 1e-6;

}

template<class T>
DACE::AlgebraicVector<T> AIDA<T>::acc_gravfield( const DACE::AlgebraicVector<T> &rr_sat, const double &et) const
{
  /*! Compute the acceleration due to the non-omogeneous gravity field.
   *  \param rr_sat An DACE::AlgebraicVector<T> for the position of the s/c.
   *  \param et the epoch in ephemeris time.
   */

  // Compute rotation matrix from ECI to ECEF
  double aux[3][3];
  DACE::AlgebraicMatrix<double> ROT(3,3,0.0);
  pxform_c( "J2000", "EARTH_FIXED", et, aux );

  for (unsigned int i=0; i<3; i++)
    for (unsigned int j=0; j<3; j++)
      ROT.at(i,j) = aux[i][j];

  // Compute acceleration in ECEF frame
  DACE::AlgebraicVector<T> aa = gravity_model.compute(ROT*rr_sat);

  // Convert acceleration from ECEF to ECI
  aa = ( ROT.transpose() )*aa;

  return (aa);

}


template<class T>
DACE::AlgebraicVector<T> AIDA<T>::acc_airdrag( const DACE::AlgebraicVector<T> &rr_sat, const DACE::AlgebraicVector<T> &vv_sat, const double &et) const
{
  /*! Compute the acceleration due to the planetary atmosphere.
   *  \param rr_sat An DACE::AlgebraicVector<T> for the position of the s/c.
   *  \param vv_sat An DACE::AlgebraicVector<T> for the velocity of the s/c.
   *  \param et the epoch in ephemeris time.
   */

  //using ASTRO::sqr;
  //using ASTRO::sqrt;
  using  std::sqrt;

  // Compute rotation matrix from ECI to ECEF
  SpiceDouble aux[6][6];
  DACE::AlgebraicMatrix<double> ROT(3,3);
    DACE::AlgebraicMatrix<double> ROTFULL(6,6);

  sxform_c  ( "J2000", "EARTH_FIXED", et, aux );

  for (unsigned int i=0; i<3; i++)
  {
    for (unsigned int j=0; j<3; j++)
      ROT.at(i,j) = aux[i][j];
  }

    for (unsigned int i=0; i<6; i++)
    {
        for (unsigned int j=0; j<6; j++)
            ROTFULL.at(i,j) = aux[i][j];
    }

    DACE::AlgebraicVector<T> xx(6);
    for (unsigned int i=0; i<3; i++)
        xx[i] = rr_sat[i];

    for (unsigned int i=3; i<6; i++)
        xx[i] = vv_sat[i-3];

    DACE::AlgebraicVector<T> xxECEF = ROTFULL*xx;


  // position in ECEF
  DACE::AlgebraicVector<T> rr_ecef = ROT*rr_sat;

  // retrieve spherical coordinates
  T alt, lon, latgd;
  position2lla( rr_ecef, this->m_frame.Origin, latgd, lon, alt);

  // compute accelerations
  DACE::AlgebraicVector<T> aa(3,0.0);

  if (DACE::cons(alt)<=2e3)
  {
    // compute atmosphere density
    T rho = p_atmosphere_model->density(et, alt, latgd*180./kPI, lon*180./kPI)*1e9;

    DACE::AlgebraicVector<T> vv_rel = vv_sat;
    // Rotating atmosphere
    if (this->getIntParam("ATMOSPHERE")==2)
    {

        for (unsigned int i=0; i<3; i++)
            vv_rel[i] = xxECEF[i+3];

    }

    T sqr_vu_rel = DACE::sqr(vv_rel[0]) + DACE::sqr(vv_rel[1]) + DACE::sqr(vv_rel[2]);
    T vu_rel = sqrt(sqr_vu_rel);

    T Bfactor = this->getUncParam("Bfactor");

    for (unsigned int i=0; i<3; i++)
          for (unsigned int j=0; j<3; j++)
              ROT.at(j,i) = aux[i][j];

    aa = (-0.5 * Bfactor * rho * vu_rel * ROT*vv_rel) * 1e-6;
  }

  return aa;

}

/***********************************************************************************
*     Dynamical model evaluation
************************************************************************************/
template<class T>
void AIDA<T>::evaluation(const double t, const DACE::AlgebraicVector< T >& X, DACE::AlgebraicVector< T >& Xdot, bool& flag) const{
/*! Evaluate the dynamical model at a certain time t for a s/c with the position and velocity
 *  specified in the X vector.
 *  Corresponding velocities and accelerations are returned throught the Xdot vector, along with
 *  a flag that describes the status of a specific event. In particular this flag becomes true
 *  whenever the distance of the s/c from the center is equal to the body radius.
 * \param t a double specifyng the time at which the dynamics is evaluated
 * \param X an DACE::AlgebraicVector for the s/c position and velocity
 * \param Xdot an DACE::AlgebraicVector for the obtained velocities and accelerations (OUTPUT)
 * \param flag boolean to check if an event is occurred or not (OUTPUT).
 * \note  This function is particularly indicated for propagations since the substantial
 *  speedup obtained by passing outputs by reference.
 */
    using std::sqrt;
    //using ASTRO::sqrt;
    //using DACE::cons;


    // Extract flag
    int flag_srp = this->getIntParam("SRP");
    int flag_atm = this->getIntParam("ATMOSPHERE");
    int flag_3bd = this->getIntParam("THIRDBODY");


    // Extract fixed parameters
    double m_radius = this->getFixParam("RADIUS");
    double m_minh   = this->getFixParam("MIN_HEIGHT");    /* Minimum height in km at which the propagation
                                                             stops (default value = 0 km).                 */

    unsigned int dim = this->getDimension();
    // check arguments
    if(X.size() != dim || Xdot.size() != dim)
        throw std::runtime_error("ASTRO::AIDA<T>::evaluation: Input size does not match the size of the dynamical model.");

    // Check impact on ground
    double R = std::sqrt( (DACE::cons(X[0])*DACE::cons(X[0])) + (DACE::cons(X[1])*DACE::cons(X[1])) + (DACE::cons(X[2])*DACE::cons(X[2])) );
    if( astro::AstroDynModel<T>::EventStatus == true ){
      if( (R - m_radius) <= m_minh )// XXX: this is the condition for which the flag becomes true and the propagator should stop
        flag = true;
    }
    DACE::AlgebraicVector<T> rr_sat(3), vv_sat(3);
    // Extract satellite position
    for (int i=0; i<3; i++) {
        rr_sat[i] = X[i];
        vv_sat[i] = X[i+3];
    }


    // Compute position of other celestial bodies, if needed
    DACE::AlgebraicVector<double> rr_sun(3,0.0);
    DACE::AlgebraicVector<double> rr_moon(3,0.0);

    if ((flag_srp!=0)||(flag_3bd==2))
      rr_sun = RelativePosition("SUN", this->m_frame, t);

    if ( (flag_srp==4) || (flag_srp==5) || (flag_srp==6) || (flag_3bd!=0) )
      rr_moon = RelativePosition("MOON", this->m_frame, t);

    // Define shadow models
    DACE::AlgebraicMatrix<double> rr_B(3,2,0.);
    if ((flag_srp==4)||(flag_srp==5)||(flag_srp==6))
      rr_B.setcol(1,rr_moon);


    DACE::AlgebraicVector<T> aa(3,0.0);

    // Gravity field
    aa = aa + acc_gravfield(rr_sat, t);

    // Air drag perturbation
    if (flag_atm!=0)
      aa = aa + acc_airdrag(rr_sat, vv_sat, t);


    if (flag_srp!=0)
      aa = aa + acc_srp(rr_sat, rr_sun, rr_B);


    if (flag_3bd!=0)
    {
      aa = aa + acc_3rdbody(rr_sat, rr_moon, this->getFixParam("GM_MOON"));

      if (flag_3bd==2)
        aa = aa + acc_3rdbody(rr_sat, rr_sun, this->getFixParam("GM_SUN"));

    }


    // Dynamical equations
    Xdot[0] = X[3];
    Xdot[1] = X[4];
    Xdot[2] = X[5];

    Xdot[3] = aa[0];
    Xdot[4] = aa[1];
    Xdot[5] = aa[2];

}

template<class T>
DACE::AlgebraicVector<T> AIDA<T>::evaluation(const double t, const DACE::AlgebraicVector< T >& X, bool& flag) const{
/*! Evaluate the dynamical model at a certain time t for a s/c with the position and velocity
 *  specified in the X vector.
 *  It returns an DACE::AlgebraicVector containing the velocities and accelerations computed, whereas
 *  a flag, describing the status of a specific event, is returned by reference. In particular
 *  this flag becomes true whenever the distance of the s/c from the center is equal to the
 *  body radius.
 * \param t a double specifyng the time at which the dynamics is evaluated
 * \param X an DACE::AlgebraicVector for the s/c position and velocity
 * \param flag boolean to check if an event is occurred or not (OUTPUT).
 * \return an DACE::AlgebraicVector for the obtained velocities and accelerations.
 * \note This function is NOT indicated for long propagations.
 * \sa AIDA<T>::evaluation
 */
    unsigned int dim = this->getDimension();
    DACE::AlgebraicVector<T> Xdot(dim);             /* Output vector.                            */
    evaluation(t,X,Xdot,flag);

    return Xdot;
}

}

#endif /* AIDA_H_ */
