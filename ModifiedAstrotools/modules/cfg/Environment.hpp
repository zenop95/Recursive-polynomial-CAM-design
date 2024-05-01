#ifndef cfg_ENVIRONMENT_H
#define cfg_ENVIRONMENT_H

#include <string>
#include <map>
#include "geco/LazyLoader.hpp"
#include "geco/NonCopyable.hpp"

namespace cfg
{
    class CEnvironment : private geco::NonCopyable
    {

        private:

            std::map<std::string, std::string> m_keyValue;

            CEnvironment(const std::string& cfgFilePath);

        public:

            void LoadSpiceKernels() const;

            const std::string& Value(const std::string& key) const;

            friend class geco::CLazyLoader<CEnvironment>;
    };

} // namespace cfg

#endif