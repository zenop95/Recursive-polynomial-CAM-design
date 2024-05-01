#ifndef atmos_SPACEWEATHER_H
#define atmos_SPACEWEATHER_H

#include <fstream>
#include <cstring>
#include <vector>
#include <array>

#include "cfg/EnvironmentOption.hpp"
#include "geco/NonCopyable.hpp"

namespace atmos
{

    /*! SpaceWeather class
     * @brief Class to handle and define the parameters for atmosphere models
     * Space weather models tables available from: http://www.celestrak.com/SpaceData/
     */

    class CSpaceWeather : private geco::NonCopyable
    {
        private:
            // Load space weather file
            void inputSW(std::ifstream &file);

        protected:
            std::vector<std::array<double,11> > _SWmatDaily;
            std::vector<std::array<double,2>  > _SWmatMonthlyPred;

        public:
            // constructor
            CSpaceWeather(const std::string swfname)
            {
                std::ifstream swfile(swfname.c_str());

                if (!swfile.is_open())
                    throw std::runtime_error("Space weather model file not found!");

                inputSW(swfile);

                swfile.close();
            }

            // compute space weather
            void computeSW(const double et, double &F107, double &F107A, std::array<double,7>& ap) const;

            friend class cfg::CEnvironmentOption<CSpaceWeather>;
    };

    extern cfg::CEnvironmentOption<CSpaceWeather> SpaceWeather;
} // namespace atmos

#endif
