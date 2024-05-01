/*
 * Propagator.h
 *
 *  Created on: Sep 9, 2014
 *      Author: Dinamica
 *  Modified on: March 2016
 *      Author: Mauro Massari (mauro.massari@polimi.it)
 */

#ifndef odeprop_PROPAGATOR_H
#define odeprop_PROPAGATOR_H

#include "DynModel.hpp"
#include <dace/AlgebraicVector.h>
#include <vector>


namespace odeprop
{
    template<typename T>
    class Propagator
    {
        protected:

          // constructor
          Propagator() = default;

        public:

          // Propagate the model starting from x0 at t0 until tf. Returns the state at final time tf.
          virtual DACE::AlgebraicVector<T> propagate(const DynModel<T> &Model, DACE::AlgebraicVector<T> const &x0, double t0, double tf) = 0;

          // Propagate the model starting from x0 at t0 until tf. Returns the entire orbit with time step tstep in a std::vector of DACE::AlgebraicVector
          virtual void orbit(const DynModel<T> &Model, DACE::AlgebraicVector<T> const &x0, double t0, double tf, double tstep, std::vector<double> &tout, std::vector< DACE::AlgebraicVector<T> > &xout) = 0;
    };

}

#endif // odeprop_PROPAGATOR_H
