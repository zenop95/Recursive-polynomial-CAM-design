#ifndef geco_JSONPARSER_HPP
#define geco_JSONPARSER_HPP

#include <string>
#include <json/json.h>

namespace geco
{
    class CJsonParser
    {
        public:
            // Read an object of type T from a Json::Value
            template<typename T>
            T Read(const Json::Value& input, const std::string& key);

            // Read an object of type T from a Json::Value or the default value
            template<typename T>
            T Read(const Json::Value& input, const std::string& key, const T& defaultVal);
    };
} // namespace geco

// Implementation
namespace geco
{
    template <typename T>
    T CJsonParser::Read(const Json::Value& input, const std::string& key)
    {
        if(input[key.c_str()].isNull())
        {
            throw std::invalid_argument("geco::CJsonParser::Read: cannot find key '" + key + "'.");
        }
        return T(input[key.c_str()], *this);
    }

    template <typename T>
    T CJsonParser::Read(const Json::Value& input, const std::string& key, const T& defaultVal)
    {
        if(input[key.c_str()].isNull())
        {
            return defaultVal;
        }
        return T(input[key.c_str()], *this);
    }

    // Parse std::string
    template<>
    std::string CJsonParser::Read(const Json::Value& input, const std::string& key);

    // Parse std::string with default value
    template<>
    std::string CJsonParser::Read(const Json::Value& input, const std::string& key, const std::string& defaultVal);

    // Parse double
    template<>
    double CJsonParser::Read(const Json::Value& input, const std::string& key);

    // Parse double with default value
    template<>
    double CJsonParser::Read(const Json::Value& input, const std::string& key, const double& defaultVal);

    // Parse int
    template<>
    int CJsonParser::Read(const Json::Value& input, const std::string& key);

    // Parse int with default value
    template<>
    int CJsonParser::Read(const Json::Value& input, const std::string& key, const int& defaultVal);

    // Parse int
    template<>
    unsigned int CJsonParser::Read(const Json::Value& input, const std::string& key);

    // Parse int with default value
    template<>
    unsigned int CJsonParser::Read(const Json::Value& input, const std::string& key, const unsigned int& defaultVal);

    // Parse bool
    template<>
    bool CJsonParser::Read(const Json::Value& input, const std::string& key);

    // Parse bool with default value
    template<>
    bool CJsonParser::Read(const Json::Value& input, const std::string& key, const bool& defaultVal);
} // namespace geco

#endif //geco_JSONPARSER_HPP
