#include "Environment.hpp"

#include <iostream>
#include <fstream>
#include <cspice/SpiceUsr.h>
#include <json/json.h>

namespace cfg
{
    CEnvironment::CEnvironment(const std::string& cfgFilePath)
    {
        // Open file
        Json::Value root;

        std::ifstream cfgFile(cfgFilePath.c_str(), std::ifstream::binary);
        cfgFile >> root;


        // Iterate through all elements
        for(Json::ValueIterator it = root.begin(); it != root.end(); it ++)
        {
            m_keyValue[it.name()] = it->asString();
        }
    }

    void CEnvironment::LoadSpiceKernels() const
    {
        std::string kernelPath;
        std::map<std::string, std::string>::const_iterator itSp = m_keyValue.find("spiceKernels");

        if(itSp != m_keyValue.end())
        {
            kernelPath = itSp->second;
        }
        else
        {
            std::cout << "No path for SPICE kernel set, use SPICE_KERNEL_FOLDER environmental variable to set it." << std::endl << "Searching in current folder" << std::endl;
            kernelPath = "./kernels.dat";
        }

        // loading kernel
        furnsh_c(kernelPath.c_str());
    }

    const std::string& CEnvironment::Value(const std::string& key) const
    {
        const std::map<std::string, std::string>::const_iterator it = m_keyValue.find(key);

        if(it== m_keyValue.end())
        {
            throw std::invalid_argument("cfg::CEnvironment::GetValue: Key not found in environment database.");
        }
        return it->second;
    }
}