/*
 * RK78Propagator.h
 *
 *  Created on: Sep 30, 2014
 *      Author: Dinamica
 *  Modified on: March 2016
 *      Author: Mauro Massari (mauo.massari@polimi.it)
 */

#ifndef odeprop_RK78PROPAGATOR_H_
#define odeprop_RK78PROPAGATOR_H_

#include "RKPropagator.hpp"


namespace odeprop
{
    template<typename T>
    class RK78Propagator : public RKPropagator<T>
    {
        protected:
          const static double alfa[];
          const static double beta[];
          const static double gamma[];
          const static double gammaStar[];

        public:
          RK78Propagator()
          {
            RKPropagator<T>::setButcher(alfa, beta, gamma, gammaStar, 13);
          }

          RK78Propagator(const Json::Value& input, geco::CJsonParser parser) : RKPropagator<T>(input, parser)
          {
            RKPropagator<T>::setButcher(alfa, beta, gamma, gammaStar, 13);
          }
};

template<typename T>
const double RK78Propagator<T>::alfa[] =      {            0.0, 2.0/27.0,    1.0/9.0,      1.0/6.0,      5.0/12.0,         0.5,       5.0/6.0,   1.0/6.0,    2.0/3.0,   1.0/3.0,        1.0,        0.0,        1.0 };

template<typename T>
const double RK78Propagator<T>::beta[] =      {            0.0,      0.0,        0.0,          0.0,           0.0,         0.0,           0.0,       0.0,        0.0,       0.0,        0.0,        0.0,   //0
                               2.0/27.0,      0.0,        0.0,          0.0,           0.0,         0.0,           0.0,       0.0,        0.0,       0.0,        0.0,        0.0,   //1
                               1.0/36.0, 1.0/12.0,        0.0,          0.0,           0.0,         0.0,           0.0,       0.0,        0.0,       0.0,        0.0,        0.0,   //2
                               1.0/24.0,      0.0,    1.0/8.0,          0.0,           0.0,         0.0,           0.0,       0.0,        0.0,       0.0,        0.0,        0.0,   //3
                               5.0/12.0,      0.0, -25.0/16.0,    25.0/16.0,           0.0,         0.0,           0.0,       0.0,        0.0,       0.0,        0.0,        0.0,   //4
                               1.0/20.0,      0.0,        0.0,         0.25,           0.2,         0.0,           0.0,       0.0,        0.0,       0.0,        0.0,        0.0,   //5
                            -25.0/108.0,      0.0,        0.0,  125.0/108.0,    -65.0/27.0,  125.0/54.0,           0.0,       0.0,        0.0,       0.0,        0.0,        0.0,   //6
                             31.0/300.0,      0.0,        0.0,          0.0,    61.0/225.0,    -2.0/9.0,    13.0/900.0,       0.0,        0.0,       0.0,        0.0,        0.0,   //7
                                    2.0,      0.0,        0.0,    -53.0/6.0,    704.0/45.0,  -107.0/9.0,     67.0/90.0,       3.0,        0.0,       0.0,        0.0,        0.0,   //8
                            -91.0/108.0,      0.0,        0.0,   23.0/108.0,  -976.0/135.0,  311.0/54.0,    -19.0/60.0,  17.0/6.0,  -1.0/12.0,       0.0,        0.0,        0.0,   //9
                          2383.0/4100.0,      0.0,        0.0, -341.0/164.0, 4496.0/1025.0, -301.0/82.0, 2133.0/4100.0, 45.0/82.0, 45.0/164.0, 18.0/41.0,        0.0,        0.0,   //10
                              3.0/205.0,      0.0,        0.0,          0.0,           0.0,   -6.0/41.0,    -3.0/205.0, -3.0/41.0,   3.0/41.0,  6.0/41.0,        0.0,        0.0,   //11
                         -1777.0/4100.0,      0.0,        0.0, -341.0/164.0, 4496.0/1025.0, -289.0/82.0, 2193.0/4100.0, 51.0/82.0, 33.0/164.0, 12.0/41.0,        0.0,        1.0 }; //12

template<typename T>
const double RK78Propagator<T>::gamma[] =     {     41.0/840.0,      0.0,        0.0,          0.0,           0.0,  34.0/105.0,      9.0/35.0,  9.0/35.0,  9.0/280.0, 9.0/280.0, 41.0/840.0,        0.0,        0.0 };

template<typename T>
const double RK78Propagator<T>::gammaStar[] = {            0.0,      0.0,        0.0,          0.0,           0.0,  34.0/105.0,      9.0/35.0,  9.0/35.0,  9.0/280.0, 9.0/280.0,        0.0, 41.0/840.0, 41.0/840.0 };
}

#endif //odeprop_RK78PROPAGATOR_H_
