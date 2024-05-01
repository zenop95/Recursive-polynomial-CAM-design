/*
 * RKPropagator.h
 *
 *  Created on: Sep 23, 2014
 *      Author: Dinamica
 *  Modified on: March 2016
 *      Author: Mauro Massari (mauo.massari@polimi.it)
 */

#ifndef odeprop_RKPROPAGATOR_H_
#define odeprop_RKPROPAGATOR_H_

#include <vector>
#include <iostream>
#include <cmath>

#include <dace/dace.h>

#include "geco/JsonParser.hpp"
#include "DynModel.hpp"
#include "VariableStepPropagator.hpp"

namespace odeprop
{
template<typename T>
class RKPropagator : public VariableStepPropagator<T>
{
    protected:

      // Butcher table storage, empty in the RKPropagator
      unsigned int m_numStages;
      const double *m_alfa;
      const double *m_beta;
      const double *m_gamma;
      const double *m_gammaStar;

      std::vector< DACE::AlgebraicVector<T> > v_fk;

      virtual void initialize(const unsigned int size, const double t0, const double tf);

      // compute the actual RK step using m_currentStepSize and estimate the error storing it in m_error return true if an event occured in DynModel
      virtual bool computeStep(const DynModel<T> &Model, DACE::AlgebraicVector<T> &x0, DACE::AlgebraicVector<T> &x, double t);

      // adapt the stepsize (m_currentStepSize), considering the error stored in m_error.
      virtual void adaptStepSize();

     public:

      // default constructor
      RKPropagator() {}

      /*! Create a new RKPropagator object using Json node and a parser.
       * \param input The Json node defining the model
       * \param parser The parser to read the node
       */
      RKPropagator(const Json::Value& input, geco::CJsonParser parser) : VariableStepPropagator<T>(input, parser) {}

      // default destructor
      ~RKPropagator() {}

      // set the pointers to already allocated arrays containg the butcher tableau, beta contains all the beta_k vectors stored rowwise (matrix has dimension stages X [stages-1] )
      void setButcher(const double *alfa, const double *beta, const double *gamma, const double *gammaStar, const unsigned int stages);

      virtual DACE::AlgebraicVector<T> propagate(const DynModel<T> &Model, DACE::AlgebraicVector<T> const &x0, double t0, double tf);
      virtual DACE::AlgebraicVector<T> propagate(const DynModel<T> &Model, DACE::AlgebraicVector<T> const &x0, double t0, double tf, double &tout, bool &eventFlag);

      virtual void orbit(const DynModel<T> &Model, DACE::AlgebraicVector<T> const &x0, double t0, double tf, double tstep, std::vector<double> &tout, std::vector< DACE::AlgebraicVector<T> > &xout);
      virtual void orbit(const DynModel<T> &Model, DACE::AlgebraicVector<T> const &x0, double t0, double tf, double tstep, std::vector<double> &tout, std::vector< DACE::AlgebraicVector<T> > &xout, bool &eventFlag);
      virtual void orbit(const DynModel<T> &Model, DACE::AlgebraicVector<T> const &x0, const std::vector<double> &tin, std::vector< DACE::AlgebraicVector<T> > &xout);
      virtual void orbit(const DynModel<T> &Model, DACE::AlgebraicVector<T> const &x0, const std::vector<double> &tin, std::vector< DACE::AlgebraicVector<T> > &xout, double &tout, bool &eventFlag);

    };

// Template implementation
template<typename T>
 void RKPropagator<T>::setButcher(const double *alfa, const double *beta, const double *gamma, const double *gammaStar, unsigned int stages)
{
  m_numStages = stages;
  m_alfa = alfa;
  m_beta = beta;
  m_gamma = gamma;
  m_gammaStar = gammaStar;

  v_fk.resize(stages);
}

// Suppress warning: comparing floating point with == or != is unsafe [-Werror=float-equal]
#ifdef __clang__
# pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined __GNUC__
# pragma GCC   diagnostic ignored "-Wfloat-equal"
#endif
template<typename T>
void RKPropagator<T>::initialize(const unsigned int size, const double t0, const double tf)
{
  // Resize each of the fk arrays to the model size
  for(typename std::vector< DACE::AlgebraicVector<T> >::iterator it = v_fk.begin(); it!=v_fk.end(); it++)
  {
    it->resize(size);
  }

  if (this->m_maxStepSize == 0.0) // no maximum step size set
  {
    // set meningfull maximum stepsize to 10% of propagation time (in this way we have at least 10 iterations)
    this->m_maxStepSize = 0.1*(tf-t0);
  }

  if (this->m_minStepSize == 0.0)
  {
    // set minimum stepsize to eps
    this->m_minStepSize = 2.2e-16;
  }

  // set initial tentative stepsize
  if(this->m_initialStepSize == 0.0)
  {
     // setting it to the 100*minstepsize if not set by user
    this->m_currentStepSize = this->m_minStepSize*100;
  }
  else
  {
    this->m_currentStepSize = this->m_initialStepSize;
  }
}
// Re-enable warnings
#ifdef __clang__
# pragma clang diagnostic warning "-Wfloat-equal"
#elif defined __GNUC__
# pragma GCC   diagnostic warning "-Wfloat-equal"
#endif


template<typename T>
bool RKPropagator<T>::computeStep(const DynModel<T> &Model, DACE::AlgebraicVector<T> &x0, DACE::AlgebraicVector<T> &x, double t)
{
  DACE::AlgebraicVector<T> xk;
  DACE::AlgebraicVector<T> xt;
  double tk;

  // compute f(t0,x0)
  //v_fk[0] = Model->evaluate(t,x0);
  bool eventFlag = false;
  Model.evaluation(t,x0,v_fk[0],eventFlag);

  // compute all the other stages
  for(unsigned int k=1; k<this->m_numStages; k++)
  {
    // compute fk
    tk=t+m_alfa[k]*this->m_currentStepSize;
    xk=x0;
    for(unsigned int j=0;j<k;j++)
    {
      xk = xk + this->m_currentStepSize * m_beta[k*(this->m_numStages-1)+j] * v_fk[j];
    }
    //v_fk[k] = Model->evaluate(tk ,xk);
    Model.evaluation(tk,xk,v_fk[k],eventFlag);
   }

  x = x0;
  xt = x0;
  for (unsigned int k=0;k<this->m_numStages;k++)
  {
    // compute y
    xt = xt + this->m_currentStepSize*m_gamma[k]*v_fk[k];

    // compute ystar (return the most accurate order)
    x = x + this->m_currentStepSize*m_gammaStar[k]*v_fk[k];
  }

  // evaluate the error
  this->m_error = DACE::vnorm(DACE::cons(xt)-DACE::cons(x));

  return eventFlag;
}

template<typename T>
void RKPropagator<T>::adaptStepSize()
{
  // simple stepsize adaptation
  double error = this->m_error;
  double tau = 1.0;

  if (this->m_error <= 2.2e-16)
  {
    error = 2.2e-16*100.0;
  }


  if(tau < 1.0e-3)
  {
    tau=1.0e-3;
  }

  if(this->m_tol<2.22e-14)
  {
      tau=2.22e-14*tau;
  }
  else
  {
      tau=(this->m_tol)*tau;
  }

  this->m_currentStepSize = 0.9 * (this->m_currentStepSize) * std::pow(tau/error,0.2);

  if (fabs(this->m_currentStepSize) > this->m_maxStepSize)
  {
    this->m_currentStepSize = this->m_currentStepSize/fabs(this->m_currentStepSize)*this->m_maxStepSize;

  }

  if (fabs(this->m_currentStepSize) < this->m_minStepSize)
  {
    this->m_currentStepSize = this->m_currentStepSize/fabs(this->m_currentStepSize)*this->m_minStepSize;
  }
}

template<typename T>
DACE::AlgebraicVector<T> RKPropagator<T>::propagate(const DynModel< T >& Model, const DACE::AlgebraicVector< T >& x0, double t0, double tf)
{
  double tout;
  bool eventFlag;
  return propagate(Model, x0, t0, tf, tout, eventFlag);
}

template<typename T>
DACE::AlgebraicVector<T> RKPropagator<T>::propagate(const DynModel<T> &Model, DACE::AlgebraicVector<T> const &x0, double t0, double tf, double &tout, bool &eventFlag)
{
//  using DACE::abs;
  using std::abs;

  DACE::AlgebraicVector<T> y0;
  DACE::AlgebraicVector<T> y;
  double currentTime;
  double backtimestep;
  unsigned int numAttempts(0);

  eventFlag = false;

  initialize(Model.getDimension(),t0,tf);

  // set time
  currentTime = t0;

  // set state
  y0 = x0;
  y = x0;

  // loop
  while(currentTime < tf && numAttempts < this->getMaxAttempts())
  {
    // check that  final time is inside time span and eventually reduce timestep
    if(currentTime+this->m_currentStepSize > tf)
    {
      // set last time to be coincident with final time
      this->m_currentStepSize=tf-currentTime;
    }

    // compute step
    eventFlag = computeStep(Model, y0, y, currentTime);
    numAttempts++;

    // check tolerance
    if (this->m_error <= this->m_tol)
    {
      // check if event is occured
      if(eventFlag)
      { // reduce stepsize to identify the event time and return (also the final time and a flag should be returned)

        backtimestep = this->m_currentStepSize/2.;

        while(backtimestep >= this->m_tolBackStep)
        {
          if(eventFlag)
          {
            this->m_currentStepSize = this->m_currentStepSize - backtimestep;
          }
          else
          {
            this->m_currentStepSize = this->m_currentStepSize + backtimestep;
          }
          backtimestep = backtimestep/2.;
          eventFlag = computeStep(Model, y0, y, currentTime);
        }

        // should set a flag and retrurn the actual final time
        currentTime += this->m_currentStepSize;
        tout = currentTime;
        eventFlag = true;
        return y;

      }
      // accept step
      y0=y;
      currentTime += this->m_currentStepSize;
      numAttempts = 0;
    }

    //adapt stepsize
    adaptStepSize();
  }

  if(numAttempts == this->getMaxAttempts())
  {
    throw std::runtime_error("RKPropagator<T>::propagate: Reached maximum number of attemps, try to increase the tolerance or reduce the stepsize");
  }

  // return results
  tout = currentTime;
  return y0;

}

template<typename T>
void RKPropagator<T>::orbit(const DynModel< T >& Model, const DACE::AlgebraicVector< T >& x0, double t0, double tf, double tstep, std::vector< double >& tout, std::vector< DACE::AlgebraicVector< T > >& xout)
{
  bool eventFlag;
  orbit(Model, x0,t0,tf,tstep, tout,xout,eventFlag);
}

template<typename T>
void RKPropagator<T>::orbit(const DynModel<T> &Model, DACE::AlgebraicVector<T> const &x0, double t0, double tf, double tstep, std::vector<double> &tout, std::vector< DACE::AlgebraicVector<T> > &xout, bool &eventFlag)
{
//  using DACE::abs;
  using std::abs;

  DACE::AlgebraicVector<T> y0;
  DACE::AlgebraicVector<T> y;
  DACE::AlgebraicVector<T> ys;
  double currentTime;
  double nextStorageTime;
  double backtimestep;
  double storagestep;
  bool actualStep = false;
  eventFlag = false;


  initialize(Model.getDimension(),t0,tf);
  xout.clear();
  tout.clear();

  // set time
  currentTime = t0;
  nextStorageTime = t0+tstep;

// Suppress warning: comparing floating point with == or != is unsafe [-Werror=float-equal]
#ifdef __clang__
# pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined __GNUC__
# pragma GCC   diagnostic ignored "-Wfloat-equal"
#endif
  if(tstep == 0)
  {
    actualStep = true;
    nextStorageTime = tf;
  }
// Re-enable warnings
#ifdef __clang__
# pragma clang diagnostic warning "-Wfloat-equal"
#elif defined __GNUC__
# pragma GCC   diagnostic warning "-Wfloat-equal"
#endif

  // set state
  y0 = x0;
  y = x0;
  xout.push_back(y0);
  tout.push_back(currentTime);

  // loop
  while(currentTime < tf)
  {
    // check that final time is inside time span and eventually reduce timestep
    if(currentTime+this->m_currentStepSize > tf)
    {
      // set last time to be coincident with final time
      this->m_currentStepSize=tf-currentTime;
      if(nextStorageTime > tf) // set also the next storage time to be coincident with final time
      {
        nextStorageTime = tf;
      }
    }

    // compute step
    eventFlag = computeStep(Model, y0, y, currentTime);
    // check tolerance
    if (this->m_error <= this->m_tol)
    {
      // check if event is occured
      if(eventFlag)
      { // reduce stepsize to identify the event time and return (also the final time and a flag should be returned)

        backtimestep = this->m_currentStepSize/2.;

        while(backtimestep >= this->m_tolBackStep)
        {
          if(eventFlag)
          {
            this->m_currentStepSize = this->m_currentStepSize - backtimestep;
          }
          else
          {
            this->m_currentStepSize = this->m_currentStepSize + backtimestep;
          }
          backtimestep = backtimestep/2.;
          eventFlag = computeStep(Model, y0, y, currentTime);
        }
        eventFlag = true;
        // should set a flag and return the actual final time
        //currentTime += this->m_currentStepSize;
        //xout.push_back(y);
        //tout.push_back(currentTime);

        //return;
      }
      // check if we are at storage point (added tolerance of 1e-14 to overcome floating point numerical errors)
      if((abs(currentTime-nextStorageTime) <= 1.0e-14) && !actualStep)
      {
        //store point
        xout.push_back(y);
        tout.push_back(currentTime+this->m_currentStepSize);
        nextStorageTime = nextStorageTime + tstep;
      }
      else if((currentTime+this->m_currentStepSize > nextStorageTime) && !actualStep)
      {
  // store actual stepsize
  storagestep = this->m_currentStepSize;
  // set stepsize to compute the nextStorageTime
  ys=y;
  while(nextStorageTime < currentTime+storagestep)
  {
    this->m_currentStepSize = nextStorageTime-currentTime;

    // propagate till that storage time
    computeStep(Model, y0, y, currentTime);

    // store point
    xout.push_back(y);
    tout.push_back(currentTime+this->m_currentStepSize);
    nextStorageTime = nextStorageTime + tstep;
  }
  if(eventFlag)
  {
    // store point
    xout.push_back(ys);
    tout.push_back(currentTime+storagestep);

  }
  // revert stepsize to original one and recompute step
  this->m_currentStepSize = storagestep;
  computeStep(Model, y0, y, currentTime);

      }
      else if((currentTime+this->m_currentStepSize < nextStorageTime) && !actualStep && eventFlag)
      {
        xout.push_back(y);
        tout.push_back(currentTime+this->m_currentStepSize);
      }
      else if(actualStep)
      {
        xout.push_back(y);
        tout.push_back(currentTime+this->m_currentStepSize);
      }
      // accept step
      y0=y;
      currentTime += this->m_currentStepSize;
      if(eventFlag){
        return;
      }

    }

    //adapt stepsize
    adaptStepSize();

  }
}

template<typename T>
void RKPropagator<T>::orbit(const DynModel< T >& Model, const DACE::AlgebraicVector< T >& x0, const std::vector< double >& tin, std::vector< DACE::AlgebraicVector< T > >& xout)
{
  bool eventFlag;
  double tout;
  orbit(Model, x0, tin, xout, tout, eventFlag);
}

template<typename T>
void RKPropagator<T>::orbit(const DynModel<T> &Model, DACE::AlgebraicVector<T> const &x0, const std::vector<double> &tin, std::vector< DACE::AlgebraicVector<T> > &xout, double &tout, bool &eventFlag)
{
  //  using DACE::abs;
  using std::abs;

  DACE::AlgebraicVector<T> y0;
  DACE::AlgebraicVector<T> y;
  DACE::AlgebraicVector<T> ys;
  double currentTime;
  double nextStorageTime;
  double backtimestep;
  double storagestep;
  bool actualStep = false;
  eventFlag = false;
  double t0 = tin.at(0);
  double tf = tin.at(tin.size() - 1);
  unsigned int iter = 0;

  initialize(Model.getDimension(),t0, tf);
  xout.clear();
  //tout.clear();

  // set time
  currentTime = tin.at(iter);
  nextStorageTime = tin.at(iter+1);

  // set state
  y0 = x0;
  y = x0;
  xout.push_back(y0);
  //tout.push_back(currentTime);

  // loop
  while (currentTime < tf)
  {
    // check that final time is inside time span and eventually reduce timestep
    if (currentTime + this->m_currentStepSize > tf)
    {
      // set last time to be coincident with final time
      this->m_currentStepSize = tf - currentTime;
      if (nextStorageTime > tf) // set also the next storage time to be coincident with final time
      {
        nextStorageTime = tf;
      }
    }

    // compute step
    eventFlag = computeStep(Model, y0, y, currentTime);
    // check tolerance
    if (this->m_error <= this->m_tol)
    {
      // check if event is occured
      if (eventFlag)
      { // reduce stepsize to identify the event time and return (also the final time and a flag should be returned)

        backtimestep = this->m_currentStepSize / 2.;

        while (backtimestep >= this->m_tolBackStep)
        {
          if (eventFlag)
          {
            this->m_currentStepSize = this->m_currentStepSize - backtimestep;
          }
          else
          {
            this->m_currentStepSize = this->m_currentStepSize + backtimestep;
          }
          backtimestep = backtimestep / 2.;
          eventFlag = computeStep(Model, y0, y, currentTime);
        }
        eventFlag = true;
        // should set a flag and retrurn the actual final time
        //currentTime += this->m_currentStepSize;
        //xout.push_back(y);
        //tout.push_back(currentTime);

        //return;
      }
      // check if we are at storage point (added tolerance of 1e-14 to overcome floating point numerical errors)
      if ((abs(currentTime - nextStorageTime) <= 1.0e-14) && !actualStep)
      {
        //store point
        xout.push_back(y);
        //tout.push_back(currentTime + this->m_currentStepSize);
        iter++;
        nextStorageTime = tin.at(iter + 1);
      }
      else if ((currentTime + this->m_currentStepSize > nextStorageTime) && !actualStep)
      {
        // store actual stepsize
        storagestep = this->m_currentStepSize;
        // set stepsize to compute the nextStorageTime
        ys = y;
        while (nextStorageTime <= currentTime + storagestep)
        {
          this->m_currentStepSize = nextStorageTime - currentTime;

          // propagate till that storage time
          computeStep(Model, y0, y, currentTime);

          // store point
          xout.push_back(y);
          //tout.push_back(currentTime + this->m_currentStepSize);
          iter++;
          if(iter+1<tin.size())
          {
            nextStorageTime = tin.at(iter + 1);
          }
          else
          {
            break;
          }
        }
        if (eventFlag)
        {
          // store point
          xout.push_back(ys);
          tout = currentTime + storagestep;

        }
        // revert stepsize to original one and recompute step
        this->m_currentStepSize = storagestep;
        computeStep(Model, y0, y, currentTime);

      }
      else if ((currentTime + this->m_currentStepSize < nextStorageTime) && !actualStep && eventFlag)
      {
        xout.push_back(y);
        tout = currentTime + this->m_currentStepSize;
      }
      else if (actualStep)
      {
        xout.push_back(y);
        tout = currentTime + this->m_currentStepSize;
      }
      // accept step
      y0 = y;
      currentTime += this->m_currentStepSize;
      if (eventFlag){
        return;
      }

    }

    //adapt stepsize
    adaptStepSize();

  }
}


}

#endif // odeprop_RKPROPAGATOR_H_
