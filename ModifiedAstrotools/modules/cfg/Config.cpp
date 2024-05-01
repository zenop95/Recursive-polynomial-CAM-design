#include "Config.hpp"

#include <iostream>
#include <cspice/SpiceUsr.h>

namespace cfg
{
    // Declaration of global variable
    CConfig EnvConfig{};

    CConfig::CConfig()
    : geco::CLazyLoader<CEnvironment>()
    {
        // Empty
    }

    CConfig::~CConfig()
    {
        // Empty
    }

    std::string CConfig::DefaultParameter() const
    {
        std::string result;
        const std::string defaultFilePath("./env.json");

        const char* envCfgFile = getenv("CDTI_ENV");
        if(envCfgFile and not std::string(envCfgFile).empty())
        {
            result = std::string(envCfgFile);
        }
        else
        {
            struct stat resultStat;
            if(stat(defaultFilePath.c_str(), &resultStat)==0)
            {
                result = defaultFilePath;
            }
        }
        return result;
    }

    void CConfig::SetEnvironmentFile(const std::string& envFilePath)
    {
        this->SetFilePath(envFilePath);
    }
} // namespace cfg
