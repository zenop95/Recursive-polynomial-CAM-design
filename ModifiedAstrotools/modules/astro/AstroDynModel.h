/*
 * AstroDynModel.h
 *
 *  Created on: September 12, 2014
 *      Author: Dinamica
 *  Modified on: March 2016
 *      Author: Mauro Massari (mauo.massari@polimi.it)
 */

#ifndef AstroDynModel_H_
#define AstroDynModel_H_

#include "odeprop/DynModel.hpp"
#include "AstroRoutines.h"
#include "odeprop/ControlModel.hpp"

namespace astro {

/*! AstroDynModel class.
 * @brief Base class for astrodynamical model (derived from DynModel class).
 *
 *  This class derives from DynModel class and it is a base class for all the astrodynamical model
 *  classes. It introduces some features like the coordinates system, the reference frame and the
 *  units of measure of the model as protected member variables.
 */

template <class T> class AstroDynModel : public odeprop::DynModel<T>
{
protected:
    Coordinate m_coord;         /*!< Coordinate system used in the dynamical model. */
    Frame m_frame;              /*!< Reference frame used in the dynamical model.   */
    Units m_units;              /*!< Units of measure used in the dynamical model.  */


    unsigned int m_NIntParam;                         /*!< Number of integer model parameters.               */
    std::vector<int> m_IntParam;                      /*!< Vector containing integer parameters values.      */
    const std::vector<std::string> *m_IntParamLabel;  /*!< Vector containing integer parameters labels.      */

    unsigned int m_NFixParam;                         /*!< Number of fixed (non-uncertain) model parameters. */
    std::vector<double> m_FixParam;                   /*!< Vector containing fixed parameters values.        */
    const std::vector<std::string> *m_FixParamLabel;  /*!< Vector containing fixed parameters labels.        */

    unsigned int m_NUncParam;                         /*!< Number of uncertain model parameters.             */
    std::vector<T> m_UncParam;                        /*!< Vector containing uncertain parameters values.    */
    const std::vector<std::string> *m_UncParamLabel;  /*!< Vector containing uncertain parameters labels.    */

public:
    odeprop::ControlModel<T>* ContMod_pnt;   /*!< Pointer to control/guidance model.      */
    bool EventStatus;               /*!< Flag to enable/disable events. */
    static bool multipleObj;        /*!< Check if other AstroDynModel objects have been already initialized. */


    /***********************************************************************************
    *     Constructor
    ************************************************************************************/
    AstroDynModel(const std::string name, const unsigned int dim,
                  const CoordType c, const FrameName f, const Body o,
		          const LengthUnit UL, const TimeUnit UT, const AngleUnit UA, odeprop::ControlModel<T>* cont = NULL)
    /*! Create a new AstroDynModel object with the given arguments.
     * \param[in] name a string with the name of the astro-dynamical model.
     * \param[in] dim a scalar value specifying the dimension of the astro-dynamical model.
     * \param[in] c an enumerated type describing the coordinate system used in the model.
     * \param[in] f an enumerated type describing the name of the reference frame used in the model.
     * \param[in] o an enumerated type describing the origin of the reference frame used in the model.
     * \param[in] UL the unit of measure for length.
     * \param[in] UT the unit of measure for times.
     * \param[in] UA the unit of measure for angles.
     * \param[in] cont the control model.
     */
    :   odeprop::DynModel<T>(name, dim)
	    {
	      Coordinate C(c);
	      m_coord = C;

	      Frame F(f,o);
	      m_frame = F;

	      Units U(UL, UT, UA);
	      m_units = U;

          if (cont != NULL)
            ContMod_pnt = cont->clone();
          else
            ContMod_pnt = NULL;

        EventStatus = true;

        }

    ~AstroDynModel(){    // Destructor.
        if (ContMod_pnt != NULL){
          delete ContMod_pnt;
          ContMod_pnt = NULL;}
    }


    /***********************************************************************************
    *     Access routines
    ************************************************************************************/
    Coordinate getCoord() const {
    /*! Get the coordinate system used in the dynamical model.
     * \return a structure containing the coordinate system of the dynamical model.
     */
        return m_coord;}

    Frame getFrame() const {
    /*! Get the reference frame used in the dynamical model.
     * \return a structure containing the reference frame of the dynamical model.
     */
        return m_frame;}

    Units getUnits() const {
    /*! Get the units of measure used in the dynamical model.
     * \return a structure containing the units of measure of the dynamical model.
     */
        return m_units;}

    int getIntParamNum() const {
    /*! Get the number of integer parameters of the dynamical model.
     * \return the number of integer parameters.
     */
      return m_NIntParam;
    }

    int getFixParamNum() const {
    /*! Get the number of fixed (non-uncertain) parameters of the dynamical model.
     * \return the number of fixed parameters.
     */
      return m_NFixParam;
    }

    int getUncParamNum() const {
    /*! Get the number of uncertain parameters of the dynamical model.
     * \return the number of uncertain parameters.
     */
      return m_NUncParam;
    }

    // set the pointers to already allocated arrays
    void setIntParamLabels(const unsigned int nlabels, const std::vector<std::string> &labels);  // Get integer parameter labels
    void setFixParamLabels(const unsigned int nlabels, const std::vector<std::string> &labels);  // Get fixed parameter labels
    void setUncParamLabels(const unsigned int nlabels, const std::vector<std::string> &labels);  // Get uncertain parameter labels

    void AllocateUncParam();  // Allocate uncertain parameters (must be performed after DA initialization if model type is DA) */

    std::vector<std::string> getIntParamLabels() const;  // Get integer parameter labels
    std::vector<std::string> getFixParamLabels() const;  // Get fixed parameter labels
    std::vector<std::string> getUncParamLabels() const;  // Get uncertain parameter labels

    void setIntParam(const std::string& param_label, const int& val);      // Set an integer parameter
    void setFixParam(const std::string& param_label, const double& val);   // Set a fixed parameter
    void setUncParam(const std::string& param_label, const T& val);        // Set an uncertain parameter

    int getIntParam(const std::string& param_label) const;     // Get an integer parameter
    double getFixParam(const std::string& param_label) const;  // Get a fixed parameter
    T getUncParam(const std::string& param_label) const;       // Get an uncertain parameter

    void setContFixParam(const std::string& param_label, const double& val);        // Set a fixed control/guidance parameter
    void setContUncParam(const std::string& param_label, const T& val);             // Set a uncertain control/guidance parameter
    std::vector<std::string> getContUncParamLabels() const;                         // Get uncertain control/guidance parameter labels
    int getContUncParamNum() const;                                                 // Get an uncertain control/guidance parameter

    //virtual void DAinit(const unsigned int ord, const unsigned int nvar) = 0;  /*!< DA initialization for dynamical model (implemented in the derived classes). */

    void EnableEvent(); // Enable event function
    void DisableEvent(); // Disable event function
};

/***********************************************************************************
*     Setup static member variables
************************************************************************************/
template<class T>
bool AstroDynModel<T>::multipleObj = false;


/***********************************************************************************
*     Setup other member variables
************************************************************************************/
template<class T>
void AstroDynModel<T>::setIntParamLabels(const unsigned int nlabels, const std::vector<std::string> &labels)
{

  /*! Set integer parameters labels
   * \param[in] nlabels number of integer parameters
   * \param[in] labels  vector of string containing the labels associated to integer parameters
   */

  m_NIntParam = nlabels;
  m_IntParam.resize(m_NIntParam);
  m_IntParamLabel = &labels;
  return;
}

template<class T>
void AstroDynModel<T>::setFixParamLabels(const unsigned int nlabels, const std::vector<std::string> &labels)
{
  /*! Set fixed (non-uncertain) parameters labels
   * \param[in] nlabels number of fixed parameters
   * \param[in] labels  vector of string containing the labels associated to fixed parameters
   */

  m_NFixParam = nlabels;
  m_FixParam.resize(nlabels);
  m_FixParamLabel = &labels;
  return;
}

template<class T>
void AstroDynModel<T>::setUncParamLabels(const unsigned int nlabels, const std::vector<std::string> &labels)
{
  /*! Set uncertain parameters labels
   * \param[in] nlabels number of uncertain parameters
   * \param[in] labels  vector of string containing the labels associated to uncertain parameters
   */

  m_NUncParam = nlabels;
  m_UncParam.resize(nlabels);
  m_UncParamLabel = &labels;
  return;
}

template<class T>
void AstroDynModel<T>::AllocateUncParam()
{
  /*! Allocate the uncertain parameters
   */
  m_UncParam.resize(m_NUncParam);

  if (ContMod_pnt != NULL)
      ContMod_pnt->AllocateContUncParam();
  return;
}

template<class T>
std::vector<std::string> AstroDynModel<T>::getIntParamLabels() const
{
  /*! Get the labels of integer parameters
   * \return a vector of strings containing the integer parameters labels.
   */

  return (*m_IntParamLabel);

}

template<class T>
std::vector<std::string> AstroDynModel<T>::getFixParamLabels() const
{
  /*! Get the labels of fixed (non-uncertain) parameters
   * \return a vector of strings containing the fixed parameters labels.
   */

  return (*m_FixParamLabel);

}

template<class T>
std::vector<std::string> AstroDynModel<T>::getUncParamLabels() const
{
  /*! Get the labels of uncertain parameters
   * \return a vector of strings containing the uncertain parameters labels.
   */

  return (*m_UncParamLabel);

}

template<class T>
void AstroDynModel<T>::setIntParam(const std::string &param_label, const int &val)
{
  /*! Set an integer parameter
   * \param[in] param_label label of the parameter to be set
   * \param[in] val         value of integer parameter
   *
   */

  // Search through parameter labels the corresponding parameter index
  int pos = -1;
  for (unsigned int i=0; i<m_NIntParam; i++)
  {
    if ( param_label.compare( m_IntParamLabel->at(i) )==0 )
    {
      pos = i;
      break;
    }
  }

  // If parameter is found set it, otherwise issue an error
  if(pos!=-1)
  {
    m_IntParam.at(pos) = val;
    return;
  }
  else
    throw std::runtime_error("ASTRO::AstroDynModel<T>::setIntParam: Parameter not found.");
}


template<class T>
void AstroDynModel<T>::setFixParam(const std::string &param_label, const double &val)
{
  /*! Set a fixed (non-uncertain) parameter
   * \param[in] param_label label of the parameter to be set
   * \param[in] val         value of fixed parameter
   *
   */

  // Search through parameter labels the corresponding parameter index
  int pos = -1;
  for (unsigned int i=0; i<m_NFixParam; i++)
    if ( param_label.compare( m_FixParamLabel->at(i) )==0 )
    {
      pos = i;
      break;
    }

  // If parameter is found set it, otherwise issue an error
  if(pos!=-1)
  {
    m_FixParam.at(pos) = val;
    return;
  }
  else
    throw std::runtime_error("ASTRO::AstroDynModel<T>::setFixParam: Parameter not found.");
}


template<class T>
void AstroDynModel<T>::setUncParam(const std::string &param_label, const T& val)
{
  /*! Set an uncertain parameter
   * \param[in] param_label label of the parameter to be set
   * \param[in] val         value of uncertain parameter
   *
   */

  // Search through parameter labels the corresponding parameter index
  int pos = -1;
  for (unsigned int i=0; i<m_NUncParam; i++)
    if ( param_label.compare( m_UncParamLabel->at(i) )==0 )
    {
      pos = i;
      break;
    }

  // If parameter is found set it, otherwise issue an error
  if(pos!=-1)
  {
    m_UncParam.at(pos) = val;
    return;
  }
  else
    throw std::runtime_error("ASTRO::AstroDynModel<T>::getUncParam: Parameter not found.");
}


template<class T>
int AstroDynModel<T>::getIntParam(const std::string &param_label) const
{
  /*! Get an integer parameter
   * \param[in] param_label label of the parameter
   * \return value of fixed parameter
   *
   */

  // Search through parameter labels the corresponding parameter index
  int pos = -1;
  for (unsigned int i=0; i<m_NIntParam; i++)
    if ( param_label.compare( m_IntParamLabel->at(i) )==0 )
    {
      pos = i;
      break;
    }

  // If parameter is found return it, otherwise issue an error
  if(pos!=-1)
  {
    return m_IntParam.at(pos);
  }
  else
    throw std::runtime_error("ASTRO::AstroDynModel<T>::getIntParam: Parameter not found.");

}


template<class T>
double AstroDynModel<T>::getFixParam(const std::string &param_label) const
{
  /*! Get a fixed (non-uncertain) parameter
   * \param[in] param_label label of the parameter
   * \return value of fixed parameter
   *
   */

  // Search through parameter labels the corresponding parameter index
  int pos = -1;
  for (unsigned int i=0; i<m_NFixParam; i++)
    if ( param_label.compare( m_FixParamLabel->at(i) )==0 )
    {
      pos = i;
      break;
    }

  // If parameter is found return it, otherwise issue an error
  if(pos!=-1)
  {
    return m_FixParam.at(pos);
  }
  else
    throw std::runtime_error("ASTRO::AstroDynModel<T>::getFixParam: Parameter not found.");
}


template<class T>
T AstroDynModel<T>::getUncParam(const std::string &param_label) const
{
  /*! Get an uncertain parameter
   * \param[in] param_label label of the parameter
   * \return value of uncertain parameter
   *
   */

  // Search through parameter labels the corresponding parameter index
  int pos = -1;
  for (unsigned int i=0; i<m_NUncParam; i++)
    if ( param_label.compare( m_UncParamLabel->at(i) )==0 )
    {
      pos = i;
      break;
    }

  // If parameter is found return it, otherwise issue an error
  if(pos!=-1)
  {
    return m_UncParam.at(pos);
  }
  else
    throw std::runtime_error("ASTRO::AstroDynModel<T>::getUncParam: Parameter not found.");
}

template<class T>
void AstroDynModel<T>::setContFixParam(const std::string &param_label, const double &val)
{
  /*! Set a fixed (non-uncertain) control/guidance parameter
   * \param[in] param_label label of the fixed control/guidance parameter to be set
   * \param[in] val         value of the fixed control/guidance parameter
   */
  if (ContMod_pnt != NULL)
    ContMod_pnt->setContParam(param_label,val);
  else
    throw std::runtime_error("ASTRO::AstroDynModel<T>::setContFixParam: Undefined control/guidance model.");

}

template<class T>
void AstroDynModel<T>::setContUncParam(const std::string &param_label, const T &val)
{
  /*! Set an uncertain control/guidance parameter
   * \param[in] param_label label of the control/guidance parameter to be set
   * \param[in] val         value of uncertain control/guidance parameter
   */
  if (ContMod_pnt != NULL)
    ContMod_pnt->setContUncParam(param_label,val);
  else
    throw std::runtime_error("ASTRO::AstroDynModel<T>::setContUncParam: Undefined control/guidance model.");
}

template<class T>
std::vector<std::string> AstroDynModel<T>::getContUncParamLabels() const
{
  /*! Get the labels of uncertain control/guidance parameters
   * \return a vector of strings containing the uncertain control/guidance parameters labels.
   */
  if (ContMod_pnt != NULL)
    return ContMod_pnt->getContUncParamLabels();
  else
    throw std::runtime_error("ASTRO::AstroDynModel<T>::getContUncParamLabels: Undefined control/guidance model.");
}

template<class T>
int AstroDynModel<T>::getContUncParamNum() const
{
  /*! Get the number of uncertain control/guidance parameters
   * \return the number of uncertain control/guidance parameters
   *
   */
  if (ContMod_pnt != NULL)
    return ContMod_pnt->getContUncParamNum();
  else
    throw std::runtime_error("ASTRO::AstroDynModel<T>::getContUncParamNum: Undefined control/guidance model.");
}

template<class T>
void AstroDynModel<T>::EnableEvent()
{
  /*! Enable the event function.   */
  AstroDynModel<T>::EventStatus = true;
}

template<class T>
void AstroDynModel<T>::DisableEvent()
{
  /*! Disable the event function.   */
  AstroDynModel<T>::EventStatus = false;
}

}
#endif /* AstroDynModel_H_ */
