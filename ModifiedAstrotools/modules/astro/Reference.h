/*
 * reference.h
 *
 *  Created on: Sep. 09, 2014
 *      Author: Dinamica
 *  Modified on: March 2016
 *      Author: Mauro Massari (mauo.massari@polimi.it)
 */

#ifndef astro_REFERENCE_H
#define astro_REFERENCE_H

#include <array>
#include <cstring>
#include <ostream>

// Spice Classes
#include <cspice/SpiceUsr.h>

#include "geco/MapEnum.h"

// Length of text output for CSPICE routines
#define LENOUT   50

#define kPI   3.14159265358979

// Astronomical Unit in kilometers
#define KM_AU 149597870.691

// Sun luminosity [MW]
#define LSUN  3.846E20

// Shortcuts for loading DACE and std routines in functions
#define USEDACE_TRIGON using ::std::atan;  using ::DACE::atan;  \
                       using ::std::asin;  using ::DACE::asin;  \
                       using ::std::acos;  using ::DACE::acos;  \
                       using ::std::sin;   using ::DACE::sin;   \
                       using ::std::cos;   using ::DACE::cos;
#define USEDACE_ALG using ::std::sqrt; using ::DACE::sqrt; \
                    using ::std::exp;  using ::DACE::exp; \
                    using ::std::log;  using ::DACE::log;
#define USEDACE_HYPBOL using ::std::acosh; using ::DACE::acosh; \
                       using ::std::asinh; using ::DACE::asinh;  \
                       using ::std::atanh; using ::DACE::atanh; \
                       using ::std::cosh;  using ::DACE::cosh;   \
                       using ::std::sinh;  using ::DACE::sinh;  \
                       using ::std::tanh;  using ::DACE::tanh;

namespace astro {

  // Standard time formats
  extern const std::string TDB_str; /*!< Standard TDB time format for text output.  */
  extern const std::string UTC_str; /*!< Standard UTC time format for text output.  */

  // Julian Date formats
  extern const std::string JD_UTC; /*!< Julian Date in Coordinated Universal Time (UTC). */
  extern const std::string JD_TDB; /*!< Julian Date in Barycentric Dynamical Time (TDB). */
  extern const std::string JD_TDT; /*!< Julian Date in Terrestrial Dynamical Time (TDT). */


  // --------------------------------------------    UNITS    --------------------------------------------

  /*! \enum LengthUnit
   * \brief
   * Allowed length units of measures
   */
  enum LengthUnit
  {
    MT,            /*!< Length unit of measure, meters, "MT"             */
    KM,            /*!< Length unit of measure, kilometers, "KM"         */
    AU,            /*!< Length unit of measure, Astronomical Units, "AU" */
    L_ADIM         /*!< Adimensional length, "ADIM"                      */
  };

  extern std::array<const char*, 4> LengthUnitLabel;

  /*! \enum TimeUnit
   * \brief
   * Allowed time units of measures
   */
  enum TimeUnit
  {
    SEC,           /*!< Time units of measure, seconds "SEC" */
    DAY,           /*!< Time units of measure, days "DAY"    */
    T_ADIM         /*!< Adimensional time, "ADIM"            */
  };

  extern std::array<const char*, 3> TimeUnitLabel;

  /*! \enum AngleUnit
   * \brief
   * Allowed angle units of measures
   */
  enum AngleUnit
  {
    DEG,           /*!< Angle units of measure, degrees, "DEG" */
    RAD,           /*!< Angle units of measure, radians, "RAD" */
    A_ADIM         /*!< Adimensional angle, "ADIM"             */
  };

  extern std::array<const char*, 3> AngleUnitLabel;

  /*! \struct Units
   * \brief A structure that defines the units of measures.
   * The length, time, and angle units of measure are contained in this structure.
   */
  struct Units{

    /*! Default constructor
     * sets default values KM, SEC and RAD.
     */
    explicit Units() : Length(KM), Time(SEC), Angle(RAD) {}

    /*! Standard constructor
     * Initializes the units of measure according to user input.
     * \param[in] UL length units
     * \param[in] UT time units
     * \param[in] UA angle units
     */

    Units(const LengthUnit UL, TimeUnit UT, const AngleUnit UA) : Length(UL), Time(UT), Angle(UA) {}

    /*! Equality operator for Units
     * Returns true if Name and Origin are the same for left-hand side and right-hand side.
     */
    bool operator==(const Units &RHS) const
    {
      return (this->Length==RHS.Length)&&(this->Time==RHS.Time)&&(this->Angle==RHS.Angle);
    }

    /*! Inequality operator for Units
     * Returns true if Name or Origin are different in left-hand side and right-hand side.
     */
    bool operator!=(const Units &RHS) const
    {
      return (this->Length!=RHS.Length)||(this->Time!=RHS.Time)||(this->Angle!=RHS.Angle);
    }

    /*! Return name of Length unit of measure
     */
    const char* getLength() const;

    /*! Return name of Time unit of measure
     */
    const char* getTime() const;

    /*! Return name of Angle unit of measure
     */
    const char* getAngle() const;

    friend std::ostream& operator<< (std::ostream &out, const Units &U)
    {
     /*! Overload of std::operator<< in iostream.
      *  \param out standard output stream.
      *  \param U units to be printed in the stream
      *  \return Standard output stream.
      */

      out << std::string("Length unit: ") << std::string(U.getLength()) << std::endl;
      out << std::string("Time unit:   ") << std::string(U.getTime())   << std::endl;
      out << std::string("Angle unit:  ") << std::string(U.getAngle());

      return out;

    }

    LengthUnit Length;  /*!< Length unit: [KM], AU, M */
    TimeUnit   Time;    /*!< Time unit:   [SEC], DAY */
    AngleUnit  Angle;   /*!< Angle unit:  [RAD], DEG */

  };

  extern const Units STDUNITS;   /*!< Definition of standard units (km, s, rad) variable */
  extern const Units AU_DAY_DEG; /*!< Definition of solar system units (AU, day, deg) variable */
  extern const Units AU_DAY_RAD; /*!< Definition of solar system units (AU, day, deg) variable */
  extern const Units ADIMUNITS;  /*!< Definition of adimentional units (km, s, rad) variable */

  // --------------------------------------------    FRAMES   --------------------------------------------

  /*! \enum FrameName
   * \brief
   * Allowed reference frame name
   */
  enum FrameName
  {
    // Inertial frames
    J2000_EQT,    /*!< Intertial reference frame, Earth Mean Equator, dynamical equinox of J2000     */
    J2000_ECL,    /*!< Intertial reference frame, Earth mean Ecliptic and equinox of the epoch J2000 */

    // Body-fixed frames
    // Listed in International Astronomical Union (IAU) reports on cartographic constants.
    IAU_SUN,      /*!< Body fixed rotating frame, Sun-centered IAU */

    IAU_MERCURY,  /*!< Body-fixed rotating frame, Mercury-centered IAU */
    IAU_VENUS,    /*!< Body-fixed rotating frame, Venus-centered IAU   */
    IAU_EARTH,    /*!< Body-fixed rotating frame, Earth-centered IAU   */
    IAU_MARS,     /*!< Body-fixed rotating frame, Mars-centered IAU    */
    IAU_JUPITER,  /*!< Body-fixed rotating frame, Jupiter-centered IAU */
    IAU_SATURN,   /*!< Body-fixed rotating frame, Saturn-centered IAU  */
    IAU_URANUS,   /*!< Body-fixed rotating frame, Uranus-centered IAU  */
    IAU_NEPTUNE,  /*!< Body-fixed rotating frame, Neptun-centered IAU  */

    IAU_MOON,     /*!< Body-fixed rotating frame, Earth's Moon-centered IAU */

    // Satellite fixed frame
    HILL,         /*!< Satellite fixed frame, Hill's reference frame */
    BODY_FIXED,   /*!< Satellite fixed frame                         */

    SYNODIC       /*!< Synodic frame */

  };
  extern std::array<const char*, 15> FrameNameLabel;

  /*! \enum Body
   * \brief
   * Allowed reference frame origin
   */
  enum Body
  {
    SUN,                  /*!< Frame reference body, Sun                            */
    MERCURY,              /*!< Frame reference body, Mercury                        */
    VENUS,                /*!< Frame reference body, Venus                          */
    EARTH,                /*!< Frame reference body, Earth                          */
    MOON,                 /*!< Frame reference body, Moon                           */
    MARS,                 /*!< Frame reference body, Mars                           */
    JUPITER,              /*!< Frame reference body, Jupiter                        */
    SATURN,               /*!< Frame reference body, Saturn                         */
    URANUS,               /*!< Frame reference body, Uranus                         */
    NEPTUNE,              /*!< Frame reference body, Neptune                        */

    SSB,                  /*!< Frame reference body, Solar System Barycentric (SSB) */
    MERCURY_BARYCENTER,   /*!< Frame reference body, Mercury barycenter             */
    VENUS_BARYCENTER,     /*!< Frame reference body, Venus barycenter               */
    EARTH_BARYCENTER,     /*!< Frame reference body, Earth barycenter               */
    MARS_BARYCENTER,      /*!< Frame reference body, Mars barycenter                */
    JUPITER_BARYCENTER,   /*!< Frame reference body, Jupiter barycenter             */
    SATURN_BARYCENTER,    /*!< Frame reference body, Saturn barycenter              */
    URANUS_BARYCENTER,    /*!< Frame reference body, Uranus barycenter              */
    NEPTUNE_BARYCENTER,   /*!< Frame reference body, Neptune barycenter             */

    CR3BP_BARYCENTER,     /*!< Frame reference body, Circular restricted 3 body problem barycenter */
    MASTER_SAT,           /*!< Frame reference body, master satellite (Hill's only) */
    SPACECRAFT            /*!< Frame reference body, spacecraft barycenter (Body fixed only) */
  };

  extern std::array<const char*, 22> BodyLabel;


  /*! \struct Frame
   * \brief A structure that defines the frame.
   * The frame name and reference body are contained in this structure.
   * \see FrameName FrameNameLabel Body BodyLabel
   * \note Some predefined frame are provided for convenience:
   * * ::ECI
   * * ::ECEF
   * * ::ECLJ2000_SC
   * * ::ECLJ2000_SSB
   * * ::ECLJ2000_SSB
   * * ::EQTJ2000_SSB
   * * ::EQTJ2000_MARS
   * * ::HILL_FRAME
   * * ::SYNODIC_FRAME
   */
  struct Frame{

    /*! Default constructor
     * sets default values J2000_EQT and EARTH.
     */
    Frame() : Name(J2000_EQT), Origin(EARTH) {}

    /*! Standard Constructor
     * Initializes the frame name and body according to user input.
     * Note that for body-fixed rotating frames the Center must match the reference celestial body.
     * \param[in] N frame name
     * \param[in] O frame origin/center
     * \sa ::FrameName, ::Body
     */
    Frame(const FrameName &N, const Body &O) : Name(N), Origin(O)
    {

      std::string tmp = this->getName();

      // Check that body-fixed rotating frames are centered on reference body for the frame.
      if ( (tmp.find("IAU_")!=std::string::npos)&&(tmp.find( this->getOrigin() )==std::string::npos) )
	throw std::runtime_error("ASTRO::Reference:: invalid frame origin for fixed frame in Frame constructor.");

      // Check that Hill's frame has origin on master satellite
      if ( ( this->Name==HILL )&&(this->Origin!=MASTER_SAT) )
	throw std::runtime_error("ASTRO::Reference:: invalid frame origin for Hill's frame in Frame constructor.");

      // Check that Hill's frame has origin on master satellite
      if ( ( this->Name==BODY_FIXED )&&(this->Origin!=SPACECRAFT) )
  throw std::runtime_error("ASTRO::Reference:: invalid frame origin for body-fixed frame in Frame constructor.");

      // Check that synodic frame has origin on barycenter
      if ( ( this->Name==SYNODIC )&&(this->Origin!=CR3BP_BARYCENTER) )
	throw std::runtime_error("ASTRO::Reference:: invalid frame origin for synodic frame in Frame constructor.");

    }

    /*! Equality operator for Frames
     * Returns true if Name and Origin are the same for left-hand side and right-hand side.
     */
    bool operator==(const Frame &RHS) const
    {
      return (this->Name==RHS.Name)&&(this->Origin==RHS.Origin);
    }

    /*! Inequality operator for Frames
     * Returns true if Name or Origin are different in left-hand side and right-hand side.
     */
    bool operator!=(const Frame &RHS) const
    {
      return (this->Name!=RHS.Name)||(this->Origin!=RHS.Origin);
    }

    /*! Get frame Name
     * Returns the frame name
     */
    const char* getName() const;

    /*! Get frame Name
     * Returns the frame name
     */
    const char* getOrigin() const;

    /*! Check if frame is inertial
     * Returns true if reference frame is J2000.
     */
    bool isInertial() const
    {
      std::string tmp( this->getName() );
      return (tmp.find("J2000")!=std::string::npos);
    }

    /*! Check if frame is Body Fixed
     * Returns true if reference frame is J2000.
     */
    bool isBodyFixed() const
    {
      std::string tmp( this->getName() );
      return (tmp.find("IAU_")!=std::string::npos);
    }

    friend std::ostream& operator<< (std::ostream &out, const Frame &F)
    {
      /*! Overload of std::operator<< in iostream.
      *  \param out standard output stream.
      *  \param F Frame to be printed in the stream
      *  \return Standard output stream.
      */

      out << std::string("Frame name:   ") << std::string(F.getName())  << std::endl;
      out << std::string("Frame origin: ") << std::string(F.getOrigin());

      return out;
    }

    FrameName Name;   /*!< Frame name: [J2000_EQT], ... */
    Body      Origin; /*!< Reference body: [Earth], ... */


  };

  // Define useful frames
  extern const Frame ECI;                    /*!< Earth Centered Inertial frame, J2000 Mean Equator */
  extern const Frame ECEF;                   /*!< Earth Centered Earth Fixed frame (IAU) */

  extern const Frame ECLJ2000_SC;               /*!< Sun Centered Inertial frame, Ecliptic J2000 */
  extern const Frame ECLJ2000_SSB;              /*!< Solar System Barycentric Inertial frame, Ecliptic J2000 */
  extern const Frame EQTJ2000_SC;               /*!< Sun Centered Inertial frame, J2000 Mean Equator */
  extern const Frame EQTJ2000_SSB;              /*!< Solar System Barycentric Inertial frame, J2000 Mean Equator */

  extern const Frame EQTJ2000_MARS;            /*!< Mars Centered Inertial frame, J2000 Mean Equator */


  /*! Hill frame, centered on master/target satellite.
   * Hill's frame is defined as follows:
   * * X-axis, aligned with radial direction, outward Earth's center.
   * * Y-axis, along-track, positive in the direction of motion
   * * Z-axis, normal to orbital plane.
   */
  extern const Frame HILL_FRAME;

  /*! Body fixed frame, centered on satellite.
   */
  extern const Frame SAT_FRAME;

  /*! Synodic frame, centered on the barycenter of the system.
   * The frame is defined as follows:
   * * X-axis, aligned with the line passing through the center of the primary and secondary
   * * Y-axis, normal to X and Z axis
   * * Z-axis, perpendicular to orbital plane of primary-secondary system.
   */
  extern const Frame SYNODIC_FRAME;

  // -------------------------------------------- COORDINATES --------------------------------------------

  /*! \enum CoordType
   * \brief
   * Allowed coordinate types
   */
  enum CoordType {

    CARTESIAN, /*!< Coordinate type label, cartesian coordinates \f$(\boldsymbol{r},\boldsymbol{v})\f$ */
    KEPLERIAN, /*!< Coordinate type label, keplerian orbital elements \f$(a,e,I,\Omega,\omega,\nu)\f$ */
    LOCAL,     /*!< Coordinate type label, local coordinates \f$(r,\alpha,\delta,v,\gamma,\psi)\f$ */
    SPHERICAL  /*!< Coordinate type label, spherical coordinates \f$(r,\alpha,\delta,\dot{r},\dot{\alpha},\dot{\delta})\f$ */

  };

  // Coordinate type map to CSPICE names
  extern std::array<const char*, 4> CoordTypeLabel;


  /*! \struct Coordinate
   * \brief A structure that defines the coordinate type.
   * The struct coordinate defines the coordinate type.
   *
   * The following coordinate types are allowed:
   * * CARTESIAN (or RECTANGULAR): the coordinate are position and velocity \f$(\boldsymbol{r},\boldsymbol{v})\f$
   * * KEPLERIAN: the state is identified by the six Keplerian elements semi-major axis, eccentricity, inclination,
   * right ascension of the ascending node, argument of the pericenter, and true anomaly \f$(a,e,I,\Omega,\omega,\nu)\f$
   * * LOCAL: the state is identified by position modulus, right ascension, declination, velocity modulus, flight path angle, and heading angle \f$(r,\alpha,\delta,v,\gamma,\psi)\f$
   * * SPHERICAL: the state is identified by the position modulus, right ascension, declination and their derivatives wrt time \f$(r,\alpha,\delta,\dot{r},\dot{\alpha},\dot{\delta})\f$
   *
   * \see CoordType CoordTypeLabel
   */
  struct Coordinate{

    /*! Default constructor
     * sets default values CARTESIAN.
     */
    Coordinate() : Type(CARTESIAN) {}

    /*! Standard Constructor
     * Initializes the coordinate type according to user input.
     * \param[in] T coordinate type
     * \see CoordType
     */
    Coordinate(const CoordType &T) : Type(T) {}

    /*! Equality operator for Coordinates
     * Returns true if Name and Origin are the same in left-hand side and right-hand side.
     */
    bool operator==(const Coordinate &RHS) const
    {
      return (this->Type==RHS.Type);
    }

    /*! Inequality operator for Coordinates
     * Returns true if Name or Origin are different in left-hand side and right-hand side.
     */
    bool operator!=(const Coordinate &RHS) const
    {
      return (this->Type!=RHS.Type);
    }

    /*! Returns string associated to current coordinate Type.
     */
    const char* getType() const;

    friend std::ostream& operator<< (std::ostream &out, const Coordinate &C)
    {
     /*! Overload of std::operator<< in iostream.
      *  \param out standard output stream.
      *  \param C Coordinate to be printed in the stream
      *  \return Standard output stream.
      */
      out << std::string("Coordinate type: ") << std::string(C.getType());

      return out;

    }

    CoordType Type; /*!< coordinate type: [CARTESIAN], SPHERICAL, KEPLERIAN, LOCAL */

  };

  /***********************************************************************************
  *     Units conversions
  ************************************************************************************/

  /*! \defgroup UnitConversion Unit Conversions
   * @brief Routines dedicated to unit conversions.
   * This module groups all routines dedicate to unit conversions.
   * @{
   */


  template <class T>
  T ConvertLength(const T &obj, const LengthUnit &UL_old, const LengthUnit &UL_new)
  {
    /*! Length unit conversion.
     * @param obj the variable to be converted
     * @param UL_old the variable original unit
     * @param UL_new the variable desired unit
     * @return The variable converted in the desired units.
     */

    // Check if the new and old units of measures are different
    if(UL_old!=UL_new )
    {
      // Check if the old units of measure is AU
      if( UL_old==AU )
      {
	if (UL_new==KM )
	{
	  return obj*KM_AU;
	}
	else if (UL_new == MT )
	{
	  return obj*KM_AU*1.e3;
	}
	else
	  throw std::runtime_error("ASTRO::ConvertLength: Unknown length conversion.");
      }
      // Check if the old units of measure is KM
      else if ( UL_old==KM )
      {
	if (UL_new==AU )
	{
	  return obj/KM_AU;
	}
	else if (UL_new==MT )
	{
	  return obj*1.e3;
	}
	else
	  throw std::runtime_error("ASTRO::ConvertLength: Unknown length conversion.");
      }
      // Check if the old units of measure is MT
      else if (UL_old==MT)
      {
	if (UL_new==AU )
	{
	  return (obj/1.e3)/KM_AU;
	}
	else if (UL_new==KM)
	  return obj/1.e3;
	else
	  throw std::runtime_error("ASTRO::ConvertLength: Unknown length conversion.");
      }
      else
	throw std::runtime_error("ASTRO::ConvertLength: Unknown length conversion.");
    }
    else
      return obj;

  }

  template <class T>
  T ConvertAngle(const T &obj, const AngleUnit &UA_old, const AngleUnit &UA_new)
  {
    /*! Angle unit conversion.
     * @param obj the variable to be converted
     * @param UA_old the variable original unit
     * @param UA_new the variable desired unit
     * @return The variable converted to the desired units.
     */

    if( UA_old!=UA_new )
    {
      if ( (UA_old==RAD)&&(UA_new==DEG) )
      {
	return obj*dpr_c(); // Multiply radians with number of degree per radian
      }
      else if ( (UA_old==DEG)&&(UA_new==RAD) )
      {
	return obj*rpd_c(); // Multiply degree with number of radians per degree
      }
      else
	throw std::runtime_error("ASTRO::ConvertAngle: Unknown angle conversion.");
    }
    else
      return obj;

  }

  template <class T>
  T ConvertTime(const T &obj, const TimeUnit &UT_old, const TimeUnit &UT_new)
  {

    /*! Time unit conversion.
     * @param obj the variable to be converted
     * @param UT_old the variable original unit
     * @param UT_new the variable desired unit
     * @return The variable converted to the desired units.
     */

    if( UT_old!=UT_new )
    {
      // Time Conversion
      if( (UT_old==SEC)&&(UT_new==DAY) )
      {
	return obj/spd_c();
      }
      else
      {
	if( (UT_old==DAY)&&(UT_new==SEC) )
	{
	  return obj*spd_c();
	}
	else
	  throw std::runtime_error("ASTRO::ConvertTime: Unknown time conversion.");
      }
    }
    else
      return obj;

  }

  template <class T>
  T ConvertVelocity(const T &obj, const LengthUnit &UL_old, const TimeUnit &UT_old, const LengthUnit &UL_new, const TimeUnit &UT_new)
  {

    /*! Velocity unit conversion.
     * @param obj the variable to be converted
     * @param UL_old the variable original length unit
     * @param UT_old the variable original time unit
     * @param UL_new the variable desired length unit
     * @param UT_new the variable desired time unit
     * @return The variable converted to the desired units.
     */

    // T tmp = 0;
    T tmp;

    // Length conversion
    if( UL_old!=UL_new )
    {
      // Check if the old units of measure is AU
      if( UL_old==AU )
      {
	if (UL_new==KM )
	{
	  tmp = obj*KM_AU;
	}
	else if ( UL_new==MT )
	{
	  tmp = obj*KM_AU*1.e3;
	}
	else
	  throw std::runtime_error("ASTRO::ConvertVelocity: Unknown length conversion.");
      }
      // Check if the old units of measure is KM
      else if ( UL_old==KM )
      {
	if ( UL_new==AU )
	{
	  tmp = obj/KM_AU;
	}
	else if (UL_new==MT )
	{
	  tmp = obj*1.e3;
	}
	else
	  throw std::runtime_error("ASTRO::ConvertVelocity: Unknown length conversion.");
      }
      // Check if the old units of measure is MT
      else if (UL_old==MT)
      {
	if (UL_new==AU )
	{
	  tmp = (obj/1.e3)/KM_AU;
	}
	else if (UL_new==KM)
	  tmp = obj/1.e3;
	else
	  throw std::runtime_error("ASTRO::ConvertVelocity: Unknown length conversion.");
      }
      else
	throw std::runtime_error("ASTRO::ConvertVelocity: Unknown length conversion.");
    }
    else
      tmp = obj;

    // Time conversion
    if( UT_old!=UT_new )
    {
      // Time Conversion
      if( (UT_old==SEC)&&(UT_new==DAY) )
      {
	return tmp*spd_c();
      }
      else
      {
	if( (UT_old==DAY)&&(UT_new==SEC) )
	{
	  return tmp/spd_c();
	}
	else
	  throw std::runtime_error("ASTRO::ConvertVelocity: Unknown time conversion.");
      }
    }
    else
      return tmp;


  }

  template <class T>
  T ConvertAngularVelocity(const T &obj, const AngleUnit &UA_old, const TimeUnit &UT_old, const AngleUnit &UA_new, const TimeUnit &UT_new)
  {

    /*! Angular velocity unit conversion.
     * @param obj the variable to be converted
     * @param UA_old the variable original length unit
     * @param UT_old the variable original time unit
     * @param UA_new the variable desired length unit
     * @param UT_new the variable desired time unit
     * @return The variable converted to the desired units.
     */

    T tmp = 0;

    // Angular conversion
    if( UA_old!=UA_new )
    {
      if( (UA_old==RAD)&&(UA_new==DEG) )
      {
	tmp =  obj*dpr_c(); // Multiply radians with number of degree per radian
      }
      else if ( (UA_old==DEG)&&(UA_new==RAD) )
      {
	tmp = obj*rpd_c(); // Multiply degree with number of radians per degree
      }
      else
	  throw std::runtime_error("ASTRO::ConvertAngularVelocity: Unknown angle conversion.");
    }
    else
      tmp = obj;

    // Time conversion
    if( UT_old!=UT_new )
    {
      // Time Conversion
      if( (UT_old==SEC)&&(UT_new==DAY) )
      {
	return tmp*spd_c();
      }
      else
      {
	if( (UT_old==DAY)&&(UT_new==SEC) )
	{
	  return tmp/spd_c();
	}
	else
	  throw std::runtime_error("ASTRO::ConvertAngularVelocity: Unknown time conversion.");
      }
    }
    else
      return tmp;


  }

  /*! @} */ // End of Unit conversions

}

#endif /* reference_H_ */
