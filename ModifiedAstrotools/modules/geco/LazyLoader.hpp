#ifndef geco_LAZYLOADER_H
#define geco_LAZYLOADER_H

#include <iostream>
#include <iomanip>
#include <memory>
#include <sys/stat.h>
#include "NonCopyable.hpp"

namespace geco
{
    template<class T>
    class CLazyLoader : private NonCopyable
    {
        private:

            // This must be shared: pointed object is only one but the derived
            // class will be shared between libraries.
            std::shared_ptr<T> p_object;
            std::unique_ptr<std::string> p_filePath;

        protected:

            std::string Parameter() const;

            virtual std::string DefaultParameter() const = 0;

        public:

            explicit CLazyLoader();
            virtual ~CLazyLoader();

            bool IsLoaded() { return p_object.get(); }

            void SetFilePath(const std::string& filePath);

            T& operator()();

            T& operator()() const;
    };

} // namespace geco

#include "LazyLoader_impl.hpp"

#endif // geco_LAZYLOADER_H