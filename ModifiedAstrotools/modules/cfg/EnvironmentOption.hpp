#ifndef cfg_ENVIRONMENTOPTION_HPP
#define cfg_ENVIRONMENTOPTION_HPP

#include <string>
#include "geco/LazyLoader.hpp"
#include "Config.hpp"

namespace cfg
{
    enum class OptionIsNone
    {
        Unknown,
        None,
        NotNone
    };

    template <typename T>
    class CEnvironmentOption : public geco::CLazyLoader<T>
    {
        protected:

            std::string s_envKey;

            OptionIsNone o_isNone;

            virtual std::string DefaultParameter() const override
            {
                return EnvConfig().Value(s_envKey);
            }

        public:

            CEnvironmentOption(const std::string& key);

            virtual ~CEnvironmentOption();

            bool IsNone()
            {
                if(not this->IsLoaded() or o_isNone == OptionIsNone::Unknown)
                {
                    o_isNone = this->Parameter() == std::string("NONE") ? OptionIsNone::None : OptionIsNone::NotNone;
                }

                return o_isNone == OptionIsNone::None;
            }
    };
} // namespace cfg

#include "EnvironmentOption_impl.hpp"

#endif // cfg_ENVIRONMENTOPTION_HPP