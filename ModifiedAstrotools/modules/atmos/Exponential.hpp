#ifndef atmos_EXPONENTIAL_H
#define atmos_EXPONENTIAL_H

#include <dace/dace.h>
#include <json/json.h>

namespace atmos
{

    /**
     * @brief Exponential density model
     *
     * @tparam T typename for computations (double, DACE::DA)
     *
     * The exponential model is defined as:
     * \f[\rho = \rho_0 \exp(\dfrac{h0-h}{H}) \f]
     */
    template <typename T>
    class CExponential
    {
    private:

        T t_rho0; ///< Reference density [kg/m^3].
        T t_h0;   ///< Reference height [km].
        T t_H;    ///< Scale height [km]

    public:

        CExponential(const Json::Value& atmosphere);

        /**
         * @brief Construct a new CExponential object
         *
         * @param rho0 Reference density [kg/m^3].
         * @param h0 Reference height [km].
         * @param H Scale height [km]
         */
        CExponential(const T& rho0, const T& h0, const T& H);

        /**
         * @brief Compute the atmosphere density at the given height above the
         * body surface.
         *
         * @param height altitude for density computation [km]
         * @return the density at given altitude
         */
        T density(const T& height) const;

    };
}

// Implementation
namespace atmos
{
    template <typename T>
    CExponential<T>::CExponential(const T& rho0, const T& h0, const T& H)
    : t_rho0(rho0), t_h0(h0), t_H(H)
    {}

    template <typename T>
    CExponential<T>::CExponential(const Json::Value& atmosphere)
    : CExponential(atmosphere["refDensity"].asDouble(),
                   atmosphere["refHeight"].asDouble(),
                   atmosphere["scaleHeight"].asDouble())
    {
    }

    template <typename T>
    T CExponential<T>::density(const T& height) const
    {
        using std::exp;
        using DACE::exp;
        return t_rho0*exp( (t_h0-height)/t_H );
    }
} // namespace atmos


#endif