#ifndef atmos_NRLMSISE00_IMP_H
#define atmos_NRLMSISE00_IMP_H

#include <cmath>
#include <array>
#include <dace/dace.h>


namespace atmos
{
    /**
     * @brief Class containing NRLMSISE-00 implementation
     *
     */
    template <typename T>
    class CNrlmsise00_p
    {
    private:

        /// @name Flags/switches
        ///@{
        std::array<int, 24>    a_switches; ///< User-defined
        std::array<double, 24> a_sw;       ///< Internal
        std::array<double, 24> a_swc;      ///< Internal
        ///@}

        /// @name PARAMB
        ///@{
        T t_gsurf; ///< Latitude dependant gravity
        T t_re;    ///< Latitude dependandt Earth radius
        ///@}

        // GTS3C
        T t_dd;

        /// @name DMIX
        ///@{
        T t_dm04, t_dm16, t_dm28, t_dm32, t_dm40, t_dm01, t_dm14;
        ///@}

        /// @name MESO7
        ///@{
        std::array<T,5> t_meso_tn1;
        std::array<T,4> t_meso_tn2;
        std::array<T,5> t_meso_tn3;
        std::array<T,2> t_meso_tgn1;
        std::array<T,2> t_meso_tgn2;
        std::array<T,2> t_meso_tgn3;
        ///@}

        /// @name LPOLY
        ///@{
        double d_dfa;
        std::array<std::array<T,9>,4> t_plg;
        T t_ctloc, t_stloc;
        T t_c2tloc, t_s2tloc;
        T t_c3tloc, t_s3tloc;
        double d_apdf;
        std::array<T,4> t_apt;
        ///@}

        /**
         * @brief Calculate latitude variable gravity (GV) and effective radius (REFF)
         *
         * @param lat latitude value
         * @param gv variable gravity
         * @param reff effective radius
         */
        static void glatf(const T& lat, T& gv, T& reff);

        /**
         * @brief Chemistry/dissociation correction for msis models
         *
         * @param alt Altitude
         * @param r Target ratio
         * @param h1 Transition scale length
         * @param zh Altitude of 1/2 r
         * @return double
         */
        static T ccor(const T& alt, const T& r, const double h1, const double zh);

        /**
         * @brief O&O2 Chemistry/dissociation correction for msis models
         *
         * @param alt Altitude
         * @param r Target ratio
         * @param h1 Transition scale length
         * @param zh Altitude of 1/2 r
         * @param h2 Transition scale length for O2
         * @return double
         */
        static T ccor2(const T& alt, const T& r, const double h1, const double zh, const double h2);

        T scalh(const T& alt, const double xm, const double temp) const;

        /**
         * @brief Turbopause correction for msis models
         *
         * @param dd diffusive density
         * @param dm full mixed density
         * @param zhm transition scale length
         * @param xmm full mixed molecular weight
         * @param xm species molecular weight
         * @return double combined density
         */
        T dnet(T& dd, const T& dm, const T& zhm, const T& xmm, const double xm) const;

        T zeta(const T& zz, double z) const;

        /**
         * @brief Calculate Temperature and Density Profiles for lower atmos
         *
         * @param alt altitude above surface (km)
         * @param d0
         * @param xm
         * @param tz
         * @param mn3
         * @param zn3
         * @param tn3
         * @param tgn3
         * @param mn2
         * @param zn2
         * @param tn2
         * @param tgn2
         * @return double
         */
        T densm (const T& alt, const T& d0, const double xm, T *tz,
                 const int mn3, const double *zn3, const T *tn3, const T *tgn3,
                 const int mn2, const double *zn2, const T *tn2, const T *tgn2) const;

        /**
         * @brief Calculate Temperature and Density Profiles for MSIS models
         *
         * @param alt altitude above surface (km)
         * @param dlb
         * @param tinf
         * @param tlb
         * @param xm
         * @param alpha
         * @param tz
         * @param zlb
         * @param s2
         * @param mn1
         * @param zn1
         * @param tn1
         * @param tgn1
         * @return double
         */
        T densu (const T& alt, const T& dlb, const T& tinf, const T& tlb, const T& xm,
                 const double alpha, T *tz, const double zlb, const T& s2,
                 const int mn1, const double *zn1, T *tn1, T *tgn1) const;

        /**
         * @brief Calculate G(L) function
         *
         * @param p coefficients
         * @param doy day of year
         * @param sec seconds in day
         * @param g_lat geodetic latitude
         * @param g_long geodetic longitude
         * @param lst local apparent solar time (hours)
         * @param f107A 81 day average of F10.7 flux (centered on doy)
         * @param f107 daily F10.7 flux for previous day
         * @param ap magnetic index array
         * @return double G(L) value
         */
        T globe7(const std::array<double,150>& p, const int doy, const double sec,
                 const T& g_lat, const T& g_long, const T& lst, const double f107A, const double f107,
                 std::array<double,7>& ap);

        /**
         * @brief Calculate G(L) function for lower atmosphere
         *
         * @param p coefficients
         * @param doy day of year
         * @param g_long geodetic longitude
         * @return double G(L) value
         */
        T glob7s(const std::array<double,100>& p, const int doy,const T& g_long);

        /**
         * @brief Thermospheric portion of NRLMSISE-00
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
        void gts7(const int doy, const double sec, const T& alt,
                 const T& g_lat, const T& g_long, const T& lst, const double f107A, const double f107,
                 std::array<double,7>& ap, std::array<T,9>& d, std::array<T,2>& t);

    public:

        /**
         * @brief Construct a new CNrlmsise00_p object
         *
         * @param switches array with nrlmsise-00 switches
         *
         * The 24 switches have the following meaning:
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
         */
        CNrlmsise00_p(const std::array<int, 24>& switches);

        /**
         * @brief Neutral Atmosphere Empircial Model from the surface to lower exosphere.
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
        void gtd7(const int doy, const double sec, const T& alt,
                 const T& g_lat, const T& g_long, const T& lst, const double f107A, const double f107,
                 std::array<double,7>& ap, std::array<T,9>& d, std::array<T,2>& t);

        /**
         * @brief Include the anomalous oxygen contribution.
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
    };

} // namespace atmos
#endif // atmos_NRLMSISE00_IMP_H
