#include "JsonParser.hpp"

namespace geco
{
    template<>
    std::string CJsonParser::Read(const Json::Value& input, const std::string& key)
    {
        if(input[key.c_str()].isNull())
        {
            throw std::invalid_argument("geco::CJsonParser::Read<std::string>: cannot find key '" + key + "'.");
        }
        return input[key.c_str()].asString();
    }

    template<>
    std::string CJsonParser::Read(const Json::Value& input, const std::string& key, const std::string& defaultVal)
    {
        return input.get(key, defaultVal).asString();
    }

    template<>
    double CJsonParser::Read(const Json::Value& input, const std::string& key)
    {
        if(input[key.c_str()].isNull())
        {
            throw std::invalid_argument("geco::CJsonParser::Read<double>: cannot find key '" + key + "'.");
        }
        return input[key.c_str()].asDouble();
    }

    template<>
    double CJsonParser::Read(const Json::Value& input, const std::string& key, const double& defaultVal)
    {
        return input.get(key, defaultVal).asDouble();
    }

    template<>
    int CJsonParser::Read(const Json::Value& input, const std::string& key)
    {
        if(input[key.c_str()].isNull())
        {
            throw std::invalid_argument("geco::CJsonParser::Read<int>: cannot find key '" + key + "'.");
        }
        return input[key.c_str()].asInt();
    }

    template<>
    int CJsonParser::Read(const Json::Value& input, const std::string& key, const int& defaultVal)
    {
        return input.get(key, defaultVal).asInt();
    }

    template<>
    unsigned int CJsonParser::Read(const Json::Value& input, const std::string& key)
    {
        if(input[key.c_str()].isNull())
        {
            throw std::invalid_argument("geco::CJsonParser::Read<int>: cannot find key '" + key + "'.");
        }
        return input[key.c_str()].asUInt();
    }

    template<>
    unsigned int CJsonParser::Read(const Json::Value& input, const std::string& key, const unsigned int& defaultVal)
    {
        return input.get(key, defaultVal).asUInt();
    }

    template<>
    bool CJsonParser::Read(const Json::Value& input, const std::string& key)
    {
        if(input[key.c_str()].isNull())
        {
            throw std::invalid_argument("geco::CJsonParser::Read<bool>: cannot find key '" + key + "'.");
        }
        return input[key.c_str()].asBool();
    }

    template<>
    bool CJsonParser::Read(const Json::Value& input, const std::string& key, const bool& defaultVal)
    {
        return input.get(key, defaultVal).asBool();
    }

} // namespace geco
