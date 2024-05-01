/*
 * VariableStepPropagator.h
 *
 *  Created on: Sep 9, 2014
 *      Author: Dinamica
 *  Modified on: March 2016
 *      Author: Mauro Massari (mauo.massari@polimi.it)
 */

#ifndef odeprop_VARIABLESTEPPROPAGATOR_H_
#define odeprop_VARIABLESTEPPROPAGATOR_H_

#include "Propagator.hpp"

namespace odeprop
{

    template<typename T>
    class VariableStepPropagator : public Propagator<T>
    {
        protected:
          double m_initialStepSize;
          double m_currentStepSize;
          double m_minStepSize;
          double m_maxStepSize;
          double m_tol;
          double m_tolBackStep;
          double m_error;
          unsigned int i_maxAttempts;

          // constructor
          VariableStepPropagator() :
          m_initialStepSize(0.0),   // default to 0 means that the implementation should set an initial stepsize
          m_currentStepSize(m_initialStepSize),  // same as iniital step size
          m_minStepSize(0.0),  // no minimum predefined, the implementation should set a default minimum step size
          m_maxStepSize(0.0), // practically no maximum predefined, the implementation should set a maximum stepsize
          m_tol(1.0e-3), // default tolerances set to 1.0e-3
          m_tolBackStep(1.0e-3),
          m_error(0.0),
          i_maxAttempts(100)
          {
              // Empty
          };

          /*! Create a new RKPropagator object using Json node and a parser.
           * \param input The Json node defining the model
           * \param parser The parser to read the node
           */
          VariableStepPropagator(const Json::Value& input, geco::CJsonParser parser) :
          m_initialStepSize(parser.Read<double>(input, "initialStep", 0.0)),
          m_currentStepSize(m_initialStepSize),
          m_minStepSize(parser.Read<double>(input, "minStep", 0.0)),
          m_maxStepSize(parser.Read<double>(input, "maxStep", 0.0)),
          m_tol(parser.Read<double>(input, "tolerance", 1.0e-3)),
          m_tolBackStep(parser.Read<double>(input, "toleranceBackstep", 1.0e-3)),
          m_error(0.0),
          i_maxAttempts(parser.Read<unsigned int>(input, "maxAttempts", 100))
          {
              // Empty
          };


        public:

          // setting parameters
          void setInitialStepSize(double arg) {m_initialStepSize = arg;}
          void setMinStepSize(double arg) {m_minStepSize = arg;}
          void setMaxStepSize(double arg) {m_maxStepSize = arg;}
          void setTolerance(double arg) {m_tol = arg;}
          void setToleranceBackStep(double arg) {m_tolBackStep = arg;}
          void setMaxAttempts(const unsigned int arg) {i_maxAttempts = arg;}

          // getting parameters
          double getInitialStepSize() {return m_initialStepSize;}
          double getMinStepSize() {return m_minStepSize;}
          double getMaxStepSize() {return m_maxStepSize;}
          double getTolerance() {return m_tol;}
          double getToleranceBackStep() {return m_tolBackStep;}
          unsigned int getMaxAttempts() {return i_maxAttempts;}

          double getError() {return m_error;}
    };
}


#endif // odeprop_VARIABLESTEPPROPAGATOR_H_
