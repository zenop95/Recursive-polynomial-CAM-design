#include "MapEnum.h"
#include <stdexcept>

namespace geco
{
    MapEnum::MapEnum(const char** EnumLabel, unsigned int size)
    {
        for (unsigned int i = 0; i < size; i++)
        {
            std::string str(EnumLabel[i]);
            m_map[ str ] = i;
        }
    }

    unsigned int MapEnum::find(const std::string &value)
    {
        std::map<std::string, unsigned int>::const_iterator iValue;
        if ( this->m_map.count(value) == 0 ){
            throw std::runtime_error("SF::MapEnum::find: Value not found.");
        }
        else{
            iValue = this->m_map.find(value);
        }

        return iValue->second;
    }
}
