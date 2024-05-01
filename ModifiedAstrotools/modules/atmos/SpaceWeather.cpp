#include "SpaceWeather.hpp"
#include <algorithm> // std::reverse
#include <cmath>
#include <cspice/SpiceUsr.h>

namespace atmos
{

    namespace
    {
        const int LENOUT = 50; // Length of text output for CSPICE routines
        const int _linelenght = 132; // line length
        const double _jdate1957 = 2436112.5; // Julian Date of 1957 10 01 00:00:00
        const double _jdate2000 = 2451544.5; // Julian Date of 2000 01 01 00:00:00
    }

    cfg::CEnvironmentOption<CSpaceWeather> SpaceWeather("spaceWeatherFile");

    void CSpaceWeather::inputSW(std::ifstream &file)
    {

        // Get rid of initial lines
        std::string aux;
        for (unsigned int i = 1; i < 18; i++)
        {
            getline(file, aux);
        }

        // Read number of observed points
        std::string str;
        getline(file, str);
        int obs_pnt = ::atoi(str.substr(20, 24).c_str());

        // Get rid of BEGIN OBSERVED
        getline(file, aux);

        // NON FUNZIONA; LO AGGIRO

        // Jump to line corresponding to 2000 01 01 00:00:00
        //int dataline = floor(_jdate2000 - _jdate1957) * _linelenght;
        std::string tmp1;
        getline(file, tmp1);
        int dataline = 15260 * _linelenght;

        file.seekg(dataline, file.cur); // originale

        while (strcmp(tmp1.substr(0, 10).c_str(), "1999 12 31") != 0 && file.eof() == false)
        {
            getline(file, tmp1);
        }

        int n_daily_obs = obs_pnt - std::floor(_jdate2000 - _jdate1957);

        std::array<double,11> defaultSWaux; defaultSWaux.fill(0.0);
        std::vector<std::array<double,11> > SWaux((size_t)n_daily_obs,defaultSWaux);

        for (int i = 0; i < n_daily_obs; i++)
        {
            std::string tmp;
            getline(file, tmp);

            SWaux.at(i)[0] = ::atof(tmp.substr(93, 97).c_str());   // F10.7 Daily
            SWaux.at(i)[1] = ::atof(tmp.substr(101, 105).c_str()); // F10.7 Average
            SWaux.at(i)[2] = ::atof(tmp.substr(79, 81).c_str());   // Daily Magnetic index
            SWaux.at(i)[3] = ::atof(tmp.substr(46, 49).c_str());   // Daily 3h APs
            SWaux.at(i)[4] = ::atof(tmp.substr(50, 53).c_str());   // Daily 3h APs
            SWaux.at(i)[5] = ::atof(tmp.substr(54, 57).c_str());   // Daily 3h APs
            SWaux.at(i)[6] = ::atof(tmp.substr(58, 61).c_str());   // Daily 3h APs
            SWaux.at(i)[7] = ::atof(tmp.substr(62, 65).c_str());   // Daily 3h APs
            SWaux.at(i)[8] = ::atof(tmp.substr(66, 69).c_str());   // Daily 3h APs
            SWaux.at(i)[9] = ::atof(tmp.substr(70, 73).c_str());   // Daily 3h APs
            SWaux.at(i)[10] = ::atof(tmp.substr(74, 77).c_str());  // Daily 3h APs
        }

        // Get rid of END OBSERVED
        for (unsigned int i = 0; i < 2; i++)
        {
            getline(file, aux);
        }

        // Read number of daily predicted points
        getline(file, str);
        int pdt_pnt = ::atoi(str.substr(27, 28).c_str());

        _SWmatDaily = SWaux;
        _SWmatDaily.resize(n_daily_obs + pdt_pnt);

        // get rid of BEGIN DAILY_PREDICTED
        getline(file, aux);

        for (int i = n_daily_obs; i < n_daily_obs + pdt_pnt; i++)
        {
            std::string tmp;
            getline(file, tmp);

            _SWmatDaily.at(i)[0] = ::atof(tmp.substr(93, 97).c_str());   // F10.7 Daily
            _SWmatDaily.at(i)[1] = ::atof(tmp.substr(101, 105).c_str()); // F10.7 Average
            _SWmatDaily.at(i)[2] = ::atof(tmp.substr(79, 81).c_str());   // Daily Magnetic index
            _SWmatDaily.at(i)[3] = ::atof(tmp.substr(46, 49).c_str());   // Daily 3h APs
            _SWmatDaily.at(i)[4] = ::atof(tmp.substr(50, 53).c_str());   // Daily 3h APs
            _SWmatDaily.at(i)[5] = ::atof(tmp.substr(54, 57).c_str());   // Daily 3h APs
            _SWmatDaily.at(i)[6] = ::atof(tmp.substr(58, 61).c_str());   // Daily 3h APs
            _SWmatDaily.at(i)[7] = ::atof(tmp.substr(62, 65).c_str());   // Daily 3h APs
            _SWmatDaily.at(i)[8] = ::atof(tmp.substr(66, 69).c_str());   // Daily 3h APs
            _SWmatDaily.at(i)[9] = ::atof(tmp.substr(70, 73).c_str());   // Daily 3h APs
            _SWmatDaily.at(i)[10] = ::atof(tmp.substr(74, 77).c_str());  // Daily 3h APs
        }

        // Get rid of END DAILY_PREDICTED
        for (unsigned int i = 0; i < 2; i++)
        {
            getline(file, aux);
        }

        // Read number of monthly predicted points
        getline(file, str);
        int mpd_pnt = ::atoi(str.substr(29, 30).c_str());

        _SWmatMonthlyPred.resize(mpd_pnt);

        // get rid of BEGIN MONTHLY_PREDICTED
        getline(file, aux);

        for (int i = 0; i < mpd_pnt; i++)
        {
            std::string tmp;
            getline(file, tmp);

            _SWmatMonthlyPred.at(i)[0] = ::atof(tmp.substr(93, 97).c_str());   // F10.7 Daily
            _SWmatMonthlyPred.at(i)[1] = ::atof(tmp.substr(101, 105).c_str()); // F10.7 Average
        }
    }

    void CSpaceWeather::computeSW(const double et, double &F107, double &F107A, std::array<double,7> &ap) const
    {
        // Outputs initializations (sets default values for atmosnrlmsise00)
        ap.fill(4.0);
        F107 = 150.;
        F107A = 150.;

        // Determine UT year, month and hour
        char tmp[79];
        et2utc_c(et, "ISOC", 14, LENOUT, tmp);
        std::string UTCstr(tmp);
        int UTyr = ::atoi(UTCstr.substr(0, 4).c_str());
        int UTmo = ::atoi(UTCstr.substr(5, 6).c_str());
        int UThr = ::atoi(UTCstr.substr(11, 12).c_str());

        // File processing
        std::array<double,32> auxMI; auxMI.fill(0.0);

        double jdate = j2000_c() + (et) / spd_c();
        int row = floor(jdate - _jdate2000);

        if (row <= (int)_SWmatDaily.size())
        {
            F107 = _SWmatDaily.at(row - 1)[0];
            F107A = _SWmatDaily.at(row)[1];

            ap[0] = _SWmatDaily.at(row)[2];
            int column = ceil((24 - UThr) / 3.);

            for (unsigned int i = 1; i < 5; i++)
            {
                auxMI[i * 8 - 8] = _SWmatDaily.at(row - (4 - i))[3];
                auxMI[i * 8 - 7] = _SWmatDaily.at(row - (4 - i))[4];
                auxMI[i * 8 - 6] = _SWmatDaily.at(row - (4 - i))[5];
                auxMI[i * 8 - 5] = _SWmatDaily.at(row - (4 - i))[6];
                auxMI[i * 8 - 4] = _SWmatDaily.at(row - (4 - i))[7];
                auxMI[i * 8 - 3] = _SWmatDaily.at(row - (4 - i))[8];
                auxMI[i * 8 - 2] = _SWmatDaily.at(row - (4 - i))[9];
                auxMI[i * 8 - 1] = _SWmatDaily.at(row - (4 - i))[10];
            }

            std::reverse(auxMI.begin(), auxMI.end());

            ap[1] = auxMI[column];
            ap[2] = auxMI[column + 1];
            ap[3] = auxMI[column + 2];
            ap[4] = auxMI[column + 3];
            ap[5] = (auxMI[column + 4] + auxMI[column + 5] + auxMI[column + 6] + auxMI[column + 7] + auxMI[column + 8] + auxMI[column + 9] + auxMI[column + 10] + auxMI[column + 11]) / 8.;
            ap[6] = (auxMI[column + 12] + auxMI[column + 13] + auxMI[column + 14] + auxMI[column + 15] + auxMI[column + 16] + auxMI[column + 17] + auxMI[column + 18] + auxMI[column + 19]) / 8.;
        }
        else
        {
            // Determine UTyr and UTmo of Monthly prediction beginning
            char tmpMP[79];
            double jed = (_jdate2000 + _SWmatDaily.size());
            et2utc_c((jed - j2000_c()) * spd_c(), "ISOC", 14, LENOUT, tmpMP);
            std::string UTCstrMP(tmpMP);

            int UTyrMP = ::atoi(UTCstrMP.substr(0, 4).c_str());
            int UTmoMP = ::atoi(UTCstrMP.substr(5, 6).c_str());
            //int UThrMP   = ::atoi( UTCstrMP.substr(11,12).c_str() );

            int dUTyr = UTyr - UTyrMP;
            int dUTmo = UTmo - UTmoMP;

            int dmon = dUTmo + dUTyr * 12;
            if (dmon < 0)
                dmon = 0;

            if (dmon < (int)_SWmatMonthlyPred.size())
            {
                F107 = _SWmatMonthlyPred.at(dmon)[0];
                F107A = _SWmatMonthlyPred.at(dmon)[1];
            }
        }
    }
} // namespace atmos