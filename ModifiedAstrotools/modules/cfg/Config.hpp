#ifndef cfg_CONFIG_H
#define cfg_CONFIG_H

#include <iostream>
#include <iomanip>
#include <memory>
#include <sys/stat.h>
#include "geco/LazyLoader.hpp"
#include "Environment.hpp"

namespace cfg
{
    class CConfig : public geco::CLazyLoader<CEnvironment>
    {
        private:

            virtual std::string DefaultParameter() const override;

        public:

            explicit CConfig();
            virtual ~CConfig();

            void SetEnvironmentFile(const std::string& envFilePath);
    };

    extern CConfig EnvConfig;

} // namespace cfg

#endif // cfg_CONFIG_H