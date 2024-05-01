/*
 * ControlModel.h
 *
 *  Created on: March 27, 2015
 *      Author: Dinamica
 *  Modified on: March 2016
 *      Author: Mauro Massari (mauo.massari@polimi.it)
 */

#ifndef odeprop_CONTROLMODEL_H
#define odeprop_CONTROLMODEL_H

#include <string>
#include <vector>
#include <dace/dace.h>

namespace odeprop
{

    /*! ControlModel class.
     * @brief Base class from which control/guidance models are derived.
     *
     *  This class is a base class from which all the control/guidance models are derived. It
     *  includes some basic and common routines and the declaration of virtual functions
     *  for the control/guidance law. These routines are then implemented for each particular
     *  control/guidance dynamics in the corresponding derived class.
     *
     */
    template <class T>
    class ControlModel
    {
        protected:
            std::string m_name;                                   /*!< Name of the control/guidance model.                            */

            unsigned int m_NContParam;                            /*!< Number of the control/guidance parameters                      */
            std::vector<double> m_ContParam;                      /*!< Vector containing control/guidance parameters values           */
            const std::vector<std::string> *m_ContParamLabel;     /*!< Vector containing control/guidance parameters labels           */

            unsigned int m_NContUncParam;                         /*!< Number of the control/guidance uncertain parameters            */
            std::vector<T> m_ContUncParam;                        /*!< Vector containing control/guidance uncertain parameters values */
            const std::vector<std::string> *m_ContUncParamLabel;  /*!< Vector containing control/guidance uncertain parameters labels */

        public:

            /***********************************************************************************
            *     Constructor
            ************************************************************************************/
            ControlModel(const std::string name) :
            /*! Create a new ControlModel object with the given name.
             * \param name The name of the control/guidance model.
             */
                m_name(name){}

            virtual ControlModel<T>* clone() const{
              return NULL;
            }
            virtual ~ControlModel<T>() {}

            /***********************************************************************************
            *     Access routines
            ************************************************************************************/
            std::string getName() const {
            /*! Get the name of the control/guidance model.
             * \return a string for the name of the control/guidance model.
             */
                return m_name; }

            void setName(std::string new_name) {
            /*! Set the name of the control/guidance model.
             * \param new_name a string containing the name of the control/guidance model.
             */
                m_name = new_name;
            }

            int getContParamNum() const {
            /*! Get the number of parameters of the control/guidance model.
             * \return the number of control/guidance parameters.
             */
              return m_NContParam;
            }

            int getContUncParamNum() const {
            /*! Get the number of uncertain parameters of the control/guidance model.
             * \return the number of control/guidance uncertain parameters.
             */
              return m_NContUncParam;
            }

            void setContParamLabels(const unsigned int nlabels, const std::vector<std::string> &labels);     // Set control/guidance parameter labels
            std::vector<std::string> getContParamLabels() const;                                             // Get control/guidance parameter labels
            void setContParam(const std::string& param_label, const double& val);                            // Set a control/guidance parameter
            double getContParam(const std::string& param_label) const;                                       // Get a control/guidance parameter

            void setContUncParamLabels(const unsigned int nlabels, const std::vector<std::string> &labels);  // Set uncertain control/guidance parameter labels
            std::vector<std::string> getContUncParamLabels() const;                                          // Get uncertain control/guidance parameter labels
            void setContUncParam(const std::string& param_label, const T& val);                              // Set a uncertain control/guidance parameter
            T getContUncParam(const std::string& param_label) const;                                         // Get a uncertain control/guidance parameter

            void AllocateContUncParam();    // Allocate uncertain control/guidance parameters (must be performed after DA initialization if model type is DA) */

            /***********************************************************************************
            *     Control/Guidance law
            ************************************************************************************/
            virtual void control_law(const double t, const DACE::AlgebraicVector<T> &X, DACE::AlgebraicVector<T> &U) const = 0;
            virtual DACE::AlgebraicVector<T> control_law(const double t, const DACE::AlgebraicVector<T> &X) const = 0;
    };

    template<class T>
    void ControlModel<T>::setContParamLabels(const unsigned int nlabels, const std::vector<std::string> &labels)
    {

      /*! Set control parameters labels
       * \param[in] nlabels number of control/guidance parameters
       * \param[in] labels  vector of string containing the labels associated to control/guidance parameters
       */

      m_NContParam = nlabels;
      m_ContParam.resize(m_NContParam);
      m_ContParamLabel = &labels;
      return;
    }

    template<class T>
    std::vector<std::string> ControlModel<T>::getContParamLabels() const
    {
      /*! Get the labels of control/guidance parameters
       * \return a vector of strings containing the control/guidance parameters labels.
       */

      return (*m_ContParamLabel);

    }

    template<class T>
    void ControlModel<T>::setContParam(const std::string &param_label, const double &val)
    {
      /*! Set an control/guidance parameter
       * \param[in] param_label label of the control/guidance parameter to be set
       * \param[in] val         value of control/guidance parameter
       *
       */

      // Search through parameter labels the corresponding parameter index
      int pos = -1;
      for (unsigned int i=0; i<m_NContParam; i++)
      {
        if ( param_label.compare( m_ContParamLabel->at(i) )==0 )
        {
          pos = i;
          break;
        }
      }

      // If parameter is found set it, otherwise issue an error
      if(pos!=-1)
      {
        m_ContParam.at(pos) = val;
        return;
      }
      else
        throw std::runtime_error("ASTRO::ControlModel<T>::setContParam: Parameter not found.");
    }

    template<class T>
    double ControlModel<T>::getContParam(const std::string &param_label) const
    {
      /*! Get an control/guidance parameter
       * \param[in] param_label label of the control/guidance parameter
       * \return value of control/guidance parameter
       *
       */

      // Search through parameter labels the corresponding parameter index
      int pos = -1;
      for (unsigned int i=0; i<m_NContParam; i++)
        if ( param_label.compare( m_ContParamLabel->at(i) )==0 )
        {
          pos = i;
          break;
        }

      // If parameter is found return it, otherwise issue an error
      if(pos!=-1)
      {
        return m_ContParam.at(pos);
      }
      else
        throw std::runtime_error("ASTRO::ControlModel<T>::getContParam: Parameter not found.");
    }

    template<class T>
    void ControlModel<T>::setContUncParamLabels(const unsigned int nlabels, const std::vector<std::string> &labels)
    {
      /*! Set uncertain control/guidance parameters labels
       * \param[in] nlabels number of uncertain control/guidance parameters
       * \param[in] labels  vector of string containing the labels associated to uncertain control/guidance parameters
       */

      m_NContUncParam = nlabels;
      m_ContUncParam.resize(nlabels);
      m_ContUncParamLabel = &labels;
      return;
    }

    template<class T>
    std::vector<std::string> ControlModel<T>::getContUncParamLabels() const
    {
      /*! Get the labels of uncertain control/guidance parameters
       * \return a vector of strings containing the uncertain control/guidance parameters labels.
       */

      return (*m_ContUncParamLabel);

    }

    template<class T>
    void ControlModel<T>::setContUncParam(const std::string &param_label, const T& val)
    {
      /*! Set an uncertain control/guidance parameter
       * \param[in] param_label label of the control/guidance parameter to be set
       * \param[in] val         value of control/guidance uncertain parameter
       *
       */

      // Search through uncertain control/guidance parameter labels the corresponding parameter index
      int pos = -1;
      for (unsigned int i=0; i<m_NContUncParam; i++)
        if ( param_label.compare( m_ContUncParamLabel->at(i) )==0 )
        {
          pos = i;
          break;
        }

      // If parameter is found set it, otherwise issue an error
      if(pos!=-1)
      {
        m_ContUncParam.at(pos) = val;
        return;
      }
      else
        throw std::runtime_error("ASTRO::ControlModel<T>::getContUncParam: Parameter not found.");
    }

    template<class T>
    T ControlModel<T>::getContUncParam(const std::string &param_label) const
    {
      /*! Get an uncertain control/guidance parameter
       * \param[in] param_label label of the control/guidance parameter
       * \return value of uncertain control/guidance parameter
       *
       */

      // Search through uncertain control/guidance parameter labels the corresponding parameter index
      int pos = -1;
      for (unsigned int i=0; i<m_NContUncParam; i++)
        if ( param_label.compare( m_ContUncParamLabel->at(i) )==0 )
        {
          pos = i;
          break;
        }

      // If parameter is found return it, otherwise issue an error
      if(pos!=-1)
      {
        return m_ContUncParam.at(pos);
      }
      else
        throw std::runtime_error("ASTRO::ControlModel<T>::getContUncParam: Parameter not found.");
    }

    template<class T>
    void ControlModel<T>::AllocateContUncParam()
    {
      /*! Allocate the uncertain control/guidance parameters
       */

      m_ContUncParam.resize(m_NContUncParam);
      return;
    }

}
#endif /* odeprop_CONTROLMODEL_H */
