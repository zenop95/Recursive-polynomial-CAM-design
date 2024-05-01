#ifdef cfg_ENVIRONMENTOPTION_HPP

namespace cfg
{
    template<typename T>
    CEnvironmentOption<T>::CEnvironmentOption(const std::string& key)
    : s_envKey(key),
      o_isNone(OptionIsNone::Unknown)
    {
        if (s_envKey.empty())
        {
            throw std::invalid_argument("geco::CEnvironmentOption: invalid key, cannot be empty!");
        }
    }

    template<typename T>
    CEnvironmentOption<T>::~CEnvironmentOption()
    {
    }
} // namespace cfg

#endif
