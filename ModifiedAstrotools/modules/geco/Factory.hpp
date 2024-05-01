#ifndef geco_FACTORY_HPP
#define geco_FACTORY_HPP

#include <memory>
#include <json/json.h>

#include "JsonParser.hpp"

namespace geco
{
    template <typename Result, typename Generator=Result>
    class CFactory
    {
        public:
            CFactory() = default;

            Result* Create(const Json::Value& input, CJsonParser& parser) const;

        protected:

            std::map<std::string, std::unique_ptr<Generator> > m_register;
    };
} // namespace geco

// Implementation
namespace geco
{
    template <typename Result, typename Generator>
    Result* CFactory<Result,Generator>::Create(const Json::Value& input, CJsonParser& parser) const
    {
        std::string generator = parser.Read<std::string>(input,"type");

        if (m_register.count(generator) == 0)
        {
            throw std::runtime_error("geco::CFactory:Create: unknown generator " + generator);
        }

        return m_register.at(generator)->Create(input, parser);
    }
} // namespace geco

#endif
