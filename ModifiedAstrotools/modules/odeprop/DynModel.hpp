/*
 * DynModel.h
 *
 *  Created on: July 22, 2014
 *      Author: Dinamica
 *  Modified on: March 2016
 *      Author: Mauro Massari (mauo.massari@polimi.it)
 */

#ifndef odeprop_DYNMODEL_H
#define odeprop_DYNMODEL_H

#include <dace/AlgebraicVector.h>

namespace odeprop
{

    /*! DynModel class.
     * @brief Base class from which dynamical models are derived.
     *
     *  This class is a base class from which all the dynamical models are derived. It
     *  includes some basic and common routines and the declaration of virtual functions
     *  for the evaluation. These routines are then implemented for each particular
     *  dynamics in the corresponding derived class.
     *
     */

    template <class T> class DynModel
    {
        protected:
            std::string m_name;                        /*!< Name of the dynamical model.                       */
            unsigned int m_dim;                        /*!< Dimension of the dynamical model.                  */

        public:
            /***********************************************************************************
            *     Constructor
            ************************************************************************************/
            DynModel(const std::string name, const unsigned int dim) :
            /*! Create a new DynModel object with the given name and dimension.
             * \param name The name of the dynamical model.
             * \param dim The model dimensions.
             */
                m_name(name),
                m_dim(dim) {}

            /***********************************************************************************
            *     Access routines
            ************************************************************************************/
            std::string getName() const {
            /*! Get the name of the dynamical model.
             * \return a string for the name of the dynamical model.
             */
                return m_name; }

            unsigned int getDimension() const {
            /*! Get the the dynamical model dimension.
             * \return a scalar value corresponding to the dimensions of the dynamical model.
             */
                return m_dim; }

            void setName(std::string new_name) {
            /*! Set the name of the dynamical model.
             * \param new_name a string containing the name of the dynamical model.
             */
                m_name = new_name;
            }

            /***********************************************************************************
            *     Dynamical model evaluation
            ************************************************************************************/
            virtual void evaluation(const double t, const DACE::AlgebraicVector<T> &X,
                DACE::AlgebraicVector<T> &Xdot, bool &flag) const  = 0;          /*!< Right hand side of the system of equations (implemented in the derived classes).        */

            virtual DACE::AlgebraicVector<T> evaluation(const double t,
                const DACE::AlgebraicVector<T> &X, bool &flag) const = 0;       /*!< Right hand side of the system of equations (implemented in the derived classes).        */
    };

}
#endif /* odeprop_DYNMODEL_H */
