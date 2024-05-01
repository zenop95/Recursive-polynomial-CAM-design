#include "Reference.h"

namespace astro
{
    std::array<const char*, 15> FrameNameLabel = { "J2000", "ECLIPJ2000",
        "IAU_SUN", "IAU_MERCURY", "IAU_VENUS", "IAU_EARTH", "IAU_MARS", "IAU_JUPITER", "IAU_SATURN", "IAU_URANUS", "IAU_NEPTUNE", "IAU_MOON",
        "HILL", "BODY_FIXED", "SYNODIC"}; /*!< Labels corresponding to FrameName enum. */

    std::array<const char*, 3> TimeUnitLabel  = { "SEC", "DAY", "ADIM" }; /*!< Labels corresponding to TimeUnit enum. */
    std::array<const char*, 3> AngleUnitLabel = { "DEG", "RAD", "ADIM" }; /*!< Labels corresponding to AngleUnit enum. */
    std::array<const char*, 4> LengthUnitLabel = { "MT", "KM", "AU", "ADIM" }; /*!< Labels corresponding to LengthUnit enum. */

    std::array<const char*, 22> BodyLabel = { "SUN",
                                              "MERCURY",
                                              "VENUS",
                                              "EARTH",
                                              "MOON",
                                              "MARS",
                                              "JUPITER",
                                              "SATURN",
                                              "URANUS",
                                              "NEPTUNE",
                                              "SSB",
                                              "MERCURY_BARYCENTER",
                                              "VENUS_BARYCENTER",
                                              "EARTH_BARYCENTER",
                                              "MARS_BARYCENTER",
                                              "JUPITER_BARYCENTER",
                                              "SATURN_BARYCENTER",
                                              "URANUS_BARYCENTER",
                                              "NEPTUNE_BARYCENTER",
                                              "CR3BP_BARYCENTER",
                                              "MASTER",
                                              "SPACECRAFT" };

    std::array<const char*, 4> CoordTypeLabel = {"RECTANGULAR", "KEPLERIAN", "LOCAL", "SPHERICAL"}; /*!< Labels corresponding to BodyLabel enum. */

    namespace
    {
        geco::MapEnum map_Lunit(LengthUnitLabel.data(),LengthUnitLabel.size());
        geco::MapEnum map_Tunit(TimeUnitLabel.data(),TimeUnitLabel.size());
        geco::MapEnum map_Aunit(AngleUnitLabel.data(),AngleUnitLabel.size());
        geco::MapEnum map_frame(FrameNameLabel.data(),FrameNameLabel.size());

                     /*!< Labels corresponding to BodyLabel enum. */
        geco::MapEnum map_body(BodyLabel.data(),BodyLabel.size());

        geco::MapEnum map_coord(CoordTypeLabel.data(),CoordTypeLabel.size());

    } // namespace


    const Units STDUNITS(KM, SEC, RAD);          /*!< Definition of standard units (km, s, rad) variable */
    const Units AU_DAY_DEG(AU,DAY,DEG);          /*!< Definition of solar system units (AU, day, deg) variable */
    const Units AU_DAY_RAD(AU,DAY,RAD);          /*!< Definition of solar system units (AU, day, deg) variable */
    const Units ADIMUNITS(L_ADIM,T_ADIM,A_ADIM); /*!< Definition of adimentional units (km, s, rad) variable */

    // Standard time formats
    const std::string TDB_str = "YYYY-Mon-DD, HR:MN:SC.### ::TDB";
    const std::string UTC_str = "YYYY-Mon-DD, HR:MN:SC.### ::UTC";

    // Julian Date formats
    const std::string JD_UTC = "JDUTC";
    const std::string JD_TDB = "JDTDB";
    const std::string JD_TDT = "JDTDT";

    // Define useful frames
    const Frame ECI(J2000_EQT, EARTH);                    /*!< Earth Centered Inertial frame, J2000 Mean Equator */
    const Frame ECEF(IAU_EARTH, EARTH);                   /*!< Earth Centered Earth Fixed frame (IAU) */

    const Frame ECLJ2000_SC(J2000_ECL,SUN);               /*!< Sun Centered Inertial frame, Ecliptic J2000 */
    const Frame ECLJ2000_SSB(J2000_ECL,SSB);              /*!< Solar System Barycentric Inertial frame, Ecliptic J2000 */
    const Frame EQTJ2000_SC(J2000_EQT,SUN);               /*!< Sun Centered Inertial frame, J2000 Mean Equator */
    const Frame EQTJ2000_SSB(J2000_EQT,SSB);              /*!< Solar System Barycentric Inertial frame, J2000 Mean Equator */

    const Frame EQTJ2000_MARS(J2000_EQT,MARS);            /*!< Mars Centered Inertial frame, J2000 Mean Equator */

    const Frame HILL_FRAME(HILL,MASTER_SAT);
    const Frame SAT_FRAME(BODY_FIXED,SPACECRAFT);
    const Frame SYNODIC_FRAME(SYNODIC,CR3BP_BARYCENTER);

    const char* Units::getLength() const { return LengthUnitLabel[Length]; }
    const char* Units::getTime() const { return TimeUnitLabel[Time]; }
    const char* Units::getAngle() const { return AngleUnitLabel[Angle]; }

    const char* Frame::getName() const { return FrameNameLabel[Name]; }
    const char* Frame::getOrigin() const { return BodyLabel[Origin]; }

    const char* Coordinate::getType() const { return CoordTypeLabel[Type]; }
}