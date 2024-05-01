/*
 * MapEnum.h
 *
 *  Created on: November 25, 2014
 *      Author: Dinamica
 *  Modified on: March 2016
 *      Author: Mauro Massari (mauo.massari@polimi.it)
 */

#ifndef geco_MAPENUM_H
#define geco_MAPENUM_H

// C++ stdlib classes
#include <map>
#include <string>

namespace geco
{
    class MapEnum
    {
        private:
            std::map<std::string, unsigned int> m_map;

        public:
            MapEnum(const char** EnumLabel, unsigned int size);

            unsigned int find(const std::string &value);
    };
}

#endif /* geco_MAPENUM_H */
