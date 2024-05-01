#ifndef atmos_NRLMSISE00_H
#define atmos_NRLMSISE00_H

#include <memory>
#include <dace/dace.h>
#include "Nrlmsise00_p.hpp"
#include "SpaceWeather.hpp"
#include "astro/AstroRoutines.h"


/**
 * @brief Atmosphere namespace
 */
namespace atmos
{
    /**
     * @brief Class to access NRLMSISE-00
     *
     * \parblock
     * @note The magnetic index array contains:
     * Array containing the following magnetic values:
     * - 0: daily AP
     * - 1: 3 hr AP index for current time
     * - 2: 3 hr AP index for 3 hrs before current time
     * - 3: 3 hr AP index for 6 hrs before current time
     * - 4: 3 hr AP index for 9 hrs before current time
     * - 5: Average of eight 3 hr AP indicies from 12 to 33 hrs prior to current time
     * - 6: Average of eight 3 hr AP indicies from 36 to 57 hrs prior to current time
     * \endparblock
     * \parblock
     * @note f107 and f107A values used to generate the model correspond
     *       to the 10.7 cm radio flux at the actual distance of the Earth
     *       from the Sun rather than the radio flux at 1 AU.
     * \endparblock
     * \parblock
     * @note f107, f107A, and ap effects are neither large nor well
     *       established below 80 km and these parameters should be set to
     *       150., 150., and 4. respectively.
     * \endparblock
     * \parblock
     * @note The 24 flags have the following meaning:
     * -  0: output in meters and kilograms instead of centimeters and grams
     * -  1: F10.7 effect on mean
     * -  2: time independent
     * -  3: symmetrical annual
     * -  4: symmetrical semiannual
     * -  5: asymmetrical annual
     * -  6: asymmetrical semiannual
     * -  7: diurnal
     * -  8: semidiurnal
     * -  9: daily ap [when this is set to -1 (!) the pointer ap_a in struct nrlmsise_input must point to a struct ap_array]
     * - 10: all UT/long effects
     * - 11: longitudinal
     * - 12: UT and mixed UT/long
     * - 13: mixed AP/UT/LONG
     * - 14: terdiurnal
     * - 15: departures from diffusive equilibrium
     * - 16: all TINF var
     * - 17: all TLB var
     * - 18: all TN1 var
     * - 19: all S var
     * - 20: all TN2 var
     * - 21: all NLB var
     * - 22: all TN3 var
     * - 23: turbo scale height var
     * \endparblock
     */
    template <typename T>
    class CNrlmsise00
    {
    private:
        std::unique_ptr< CNrlmsise00_p<T> > p_detail; ///< Pointer to NRLMSISE-00 implementation

    public:

        /**
         * @brief Construct a new CNrlmsise00 object
         *
         * @param flags to turn on and off particular variations.
         *
         * @note 0 is off, 1 is on, and 2 is main effects off but cross terms on.
         */
        CNrlmsise00(const std::array<int,24>& flags);

        /**
         * @brief Destroy the CNrlmsise00 object
         */
        ~CNrlmsise00() = default;

        /**
         * @brief Neutral Atmosphere Empirical Model from the surface to lower exosphere.
         *
         * @param doy day of year
         * @param sec seconds in day
         * @param alt altitude above surface (km)
         * @param g_lat geodetic latitude
         * @param g_long geodetic longitude
         * @param lst local apparent solar time (hours)
         * @param f107A 81 day average of F10.7 flux (centered on doy)
         * @param f107 daily F10.7 flux for previous day
         * @param ap magnetic index array
         * @param d density array
         * @param t temperature array
         *
         * The output variables are:
         *  - d[0]: HE NUMBER DENSITY(CM-3)
         *  - d[1]: O NUMBER DENSITY(CM-3)
         *  - d[2]: N2 NUMBER DENSITY(CM-3)
         *  - d[3]: O2 NUMBER DENSITY(CM-3)
         *  - d[4]: AR NUMBER DENSITY(CM-3)
         *  - d[5]: TOTAL MASS DENSITY(GM/CM3) [includes d[8] in td7d]
         *  - d[6]: H NUMBER DENSITY(CM-3)
         *  - d[7]: N NUMBER DENSITY(CM-3)
         *  - d[8]: Anomalous oxygen NUMBER DENSITY(CM-3)
         *  - t[0]: EXOSPHERIC TEMPERATURE
         *  - t[1]: TEMPERATURE AT ALT
         *
         * @note By setting flag 0 output can be obtained in kg and m^3
         *
         * @note d[5] is the sum of the mass densities of the
         *        species labeled by indices 0-4 and 6-7 in output variable d.
         *        This includes He, O, N2, O2, Ar, H, and N but does NOT include
         *        anomalous oxygen (species index 8).
         */
        void gtd7(const int doy, const double sec, const T& alt,
                 const T& g_lat, const T& g_long, const T& lst, const double f107A, const double f107,
                 std::array<double,7>& ap, std::array<T,9>& d, std::array<T,2>& t);

        /**
         * @brief Neutral Atmosphere Empirical Model from the surface to lower exosphere, including the anomalous oxygen contribution.
         *
         * This subroutine provides Effective Total Mass Density for output
         *   d[5] which includes contributions from "anomalous oxygen" which can
         *   affect satellite drag above 500 km.
         *
         * @param doy day of year
         * @param sec seconds in day
         * @param alt altitude above surface (km)
         * @param g_lat geodetic latitude
         * @param g_long geodetic longitude
         * @param lst local apparent solar time (hours)
         * @param f107A 81 day average of F10.7 flux (centered on doy)
         * @param f107 daily F10.7 flux for previous day
         * @param ap magnetic index array
         * @param d density array
         * @param t temperature array
         */
        void gtd7d(const int doy, const double sec, const T& alt,
                  const T& g_lat, const T& g_long, const T& lst, const double f107A, const double f107,
                  std::array<double,7>& ap, std::array<T,9>& d, std::array<T,2>& t);

        /**
         * @brief Compute the Effective Total Mass Density, including contributions from "anomalous oxygen"
         *
         * The local apparent solar time is internally computed as:
         * \f[lst=\dfrac{\mathrm{sec}}{3600} + \dfrac{\lambda_g}{15}\f]
         *
         * @param doy day of year
         * @param sec seconds in day
         * @param alt altitude above surface (km)
         * @param g_lat geodetic latitude
         * @param g_long geodetic longitude
         * @param f107A 81 day average of F10.7 flux (centered on doy)
         * @param f107 daily F10.7 flux for previous day
         * @param ap magnetic index array
         * @return double the total atmosphere density
         */
        T density(const int doy, const double sec,
                  const T& alt, const T& g_lat, const T& g_long,
                  const double f107A, const double f107, std::array<double,7>& ap);

        /**
         * @brief Compute the Effective Total Mass Density, including contributions from "anomalous oxygen"
         *
         * The local apparent solar time is internally computed as:
         * \f[lst=\dfrac{\mathrm{sec}}{3600} + \dfrac{\lambda_g}{15}\f]
         *
         * @param et ephemeris time
         * @param alt altitude above surface (km)
         * @param g_lat geodetic latitude
         * @param g_long geodetic longitude
         * @return double the total atmosphere density
         */
        T density(const double et,  const T& alt, const T& g_lat, const T& g_long);


    };
}

// Implementation
namespace atmos
{
    template<typename T>
    CNrlmsise00<T>::CNrlmsise00(const std::array<int, 24>& flags)
    {
        p_detail.reset(new CNrlmsise00_p<T>(flags));
    }

    template<typename T>
    void CNrlmsise00<T>::gtd7(const int doy, const double sec, const T& alt,
                 const T& g_lat, const T& g_long, const T& lst, const double f107A, const double f107,
                 std::array<double,7>& ap, std::array<T,9>& d, std::array<T,2>& t)
    {
        p_detail->gtd7(doy, sec, alt, g_lat, g_long, lst, f107A, f107, ap, d, t);
    }

    template<typename T>
    void CNrlmsise00<T>::gtd7d(const int doy, const double sec, const T& alt,
                 const T& g_lat, const T& g_long, const T& lst, const double f107A, const double f107,
                 std::array<double,7>& ap, std::array<T,9>& d, std::array<T,2>& t)
    {
        p_detail->gtd7d(doy, sec, alt, g_lat, g_long, lst, f107A, f107, ap, d, t);
    }

    template<typename T>
    T CNrlmsise00<T>::density(const int doy, const double sec,
                  const T& alt, const T& g_lat, const T& g_long,
                  const double f107A, const double f107, std::array<double,7>& ap)
    {
        std::array<T,9> dens;
        std::array<T,2> temp;
        T lst = sec/3600.0 + g_long/15.0;
        p_detail->gtd7d(doy, sec, alt, g_lat, g_long, lst, f107A, f107, ap, dens, temp);

        return dens.at(5);
    }

    template<typename T>
    T CNrlmsise00<T>::density(const double et,  const T& alt, const T& g_lat, const T& g_long)
    {
        // Get the space weather values
        double f107A, f107;
        std::array<double,7> ap;
        SpaceWeather().computeSW(et, f107, f107A, ap);

        // Compute year, day of year and seconds from day start
        int year = ::atoi(astro::et2str(et, "YYYY ::UTC").c_str());
        int doy  = ::atoi(astro::et2str(et, "DOY ::UTC").c_str());
        char buffer[150];
        sprintf(buffer,"%d-%dT00:00:00.000",year,doy);
        double sec = et - astro::str2et(buffer);

        return this->density(doy, sec, alt, g_lat, g_long, f107A, f107, ap);
    }
}

#endif