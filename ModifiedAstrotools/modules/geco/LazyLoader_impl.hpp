#ifdef geco_LAZYLOADER_H

#include <iostream>
#include <cspice/SpiceUsr.h>

namespace geco
{
    template<class T>
    CLazyLoader<T>::CLazyLoader()
    {
        // Empty
    }

    template<class T>
    CLazyLoader<T>::~CLazyLoader()
    {
        // Empty
    }

    template<class T>
    std::string CLazyLoader<T>::Parameter() const
    {
        if (p_filePath)
        {
            return *p_filePath;
        }
        else
        {
            return DefaultParameter();
        }
    }

    template<class T>
    void CLazyLoader<T>::SetFilePath(const std::string& filePath)
    {
        if(IsLoaded())
        {
            throw std::runtime_error("geco::CLazyLoader::IsAnyLoaded(): cannot set environment file path, environment was already loaded!");
        }
        p_filePath = std::make_unique<std::string>(filePath);
    }

    template<class T>
    T& CLazyLoader<T>::operator()()
    {
        if(not p_object)
        {
            p_object.reset(new T(this->Parameter()));
        }
        return *p_object;
    }

    template<class T>
    T& CLazyLoader<T>::operator()() const
    {
        if(not p_object)
        {
            throw std::runtime_error("geco::CLazyLoader::operator(): cannot load object from " + Parameter() + " in const function.");
        }
        return *p_object;
    }
} // namespace geco

#endif