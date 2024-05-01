/*
 * GravityFieldModel.h
 *
 *  Created on: Feb 25, 2015
 *      Author: Alessandro Morselli
 *  Modified on: March 2016
 *      Author: Mauro Massari (mauo.massari@polimi.it)
 */

#ifndef GravityFieldModel_H_
#define GravityFieldModel_H_

#include <fstream>
#include <iostream>
#include "AstroCnvRef.h"

#include <dace/dace.h>
#include "astro/Reference.h"
#include "geco/JsonParser.hpp"

namespace dynorb{


  /*! GravityFieldModel class
   * @brief Class to handle and define gravity models
   * Gravity models tables available from: http://icgem.gfz-potsdam.de/ICGEM/modelstab.html
   */

  class GravityFieldModel {

  private:

    std::string _gravmodelname;
    double _GM;
    double _RADIUS;
    unsigned int _MAX_DEG; /*! Maximum degree of model to be used */
    unsigned int _degree;
    DACE::AlgebraicMatrix<double> _LSF;

    void LegendreScaleFactor(unsigned int max_degree);

    template<class T>
    DACE::AlgebraicMatrix<T> gravLegendre(const T& phi) const;

    template<class T>
    DACE::AlgebraicVector<T> gravPCPF(const DACE::AlgebraicVector<T>& rr, const T& rr_norm, const DACE::AlgebraicMatrix<T>& P, const DACE::AlgebraicVector<T>& smlambda, const DACE::AlgebraicVector<T>& cmlambda) const;


  protected:

    DACE::AlgebraicMatrix<double> _S;
    DACE::AlgebraicMatrix<double> _C;

  public:
    // constructor
    GravityFieldModel(const std::string name, const unsigned int degree):
    _gravmodelname(name), _degree(degree)
    {
      std::string binfname = "" + name + ".bin";
      std::ifstream gravfile( binfname.c_str(), std::ostream::binary );

      if(!gravfile.is_open())
	throw std::runtime_error("Binary gravity model file not found!");

      gravfile.read((char*)&_MAX_DEG, sizeof(int));
//       std::cout << _MAX_DEG << std::endl;

      if (!(_degree<_MAX_DEG))
        throw std::runtime_error("dynorb::GravityFieldModel: Selected number of spherical harmonics exceeds maximum");

      gravfile.read((char*)&_GM, sizeof(double));
//       std::cout << _GM << std::endl;

      gravfile.read((char*)&_RADIUS, sizeof(double));
      //std::cout << _RADIUS << std::endl;

      // double C[_MAX_DEG+1][_MAX_DEG+1];
      // double S[_MAX_DEG+1][_MAX_DEG+1];
      double *C = new double[(_MAX_DEG+1)*(_MAX_DEG+1)];
      double *S = new double[(_MAX_DEG+1)*(_MAX_DEG+1)];

      gravfile.read((char*)C, sizeof(double)*(_MAX_DEG+1)*(_MAX_DEG+1));
      gravfile.read((char*)S, sizeof(double)*(_MAX_DEG+1)*(_MAX_DEG+1));

      gravfile.close();

      _C.resize(_degree+2, _degree+2);
      _S.resize(_degree+2, _degree+2);
      for (unsigned int i=0; i<_degree+2; i++)
      {
	for (unsigned int j=0; j<_degree+2; j++)
	{
    // _C.at(i,j) = C[i][j];
    // _S.at(i,j) = S[i][j];
	  _C.at(i,j) = C[i*(_MAX_DEG+1)+j];
	  _S.at(i,j) = S[i*(_MAX_DEG+1)+j];
// std::cout << _C << std::endl;
// std::cout << _S << std::endl;
	}
      }

      LegendreScaleFactor(_degree);
      delete[] C;
      delete[] S;
    }

    GravityFieldModel(const Json::Value& input, geco::CJsonParser& parser)
    : GravityFieldModel(parser.Read<std::string>(input, "model"), parser.Read<uint>(input, "degree"))
    {
        // Empty
    }

    template<class T>
    DACE::AlgebraicVector<T> compute(const DACE::AlgebraicVector<T>& rr) const;

    double getGM() const { return _GM; };
    double getRadius() const { return _RADIUS; };
    const std::string& getName() const { return _gravmodelname; }
    unsigned int getDegree() const { return _degree; }

//     DACE::AlgebraicMatrix<double> getS() { return _S; };
//     DACE::AlgebraicMatrix<double> getC() { return _C; };
//     DACE::AlgebraicMatrix<double> getLSF() {return _LSF; };

  };

void GravityFieldModel::LegendreScaleFactor(unsigned int max_degree)
{
      /*! Internal function computing normalized associated
       * Legendre polynomials, P, via recursion relations for spherical harmonic
       * gravity
       */

      // Resize container of Legendre polynomials
      _LSF.resize(max_degree+3,max_degree+3);

      // Seeds for recursion formula
      _LSF.at(0,0) = 0.;
      _LSF.at(1,0) = 1.;
      _LSF.at(1,1) = 0.;

      for (unsigned int n=2; n<_degree+3; n++)
      {

	for (unsigned int m=0; m<n+1; m++)
	{
	  // Scale Factor needed for normalization of dUdphi partial derivative

	  if (n==m)
	    _LSF.at(n,n) = 0.;
	  else if (m==0)
	    _LSF.at(n,m) = std::sqrt(0.5*( (double)n + 1. )*( (double)n ));
	  else
	    _LSF.at(n,m) = std::sqrt( (double)(n + m +1)*(n-m) );

	}

      }

    }


template<class T>
DACE::AlgebraicMatrix<T> GravityFieldModel::gravLegendre(const T &phi) const
{

  using  std::sin;
  //using DACE::sin;
  using  std::cos;
  //using DACE::cos;

  DACE::AlgebraicMatrix<T> P(_degree+3, _degree+3, 0.);

  T cphi = cos(kPI/2. - phi);
  T sphi = sin(kPI/2. - phi);

  // Seeds for recursion formula
  P.at(0,0) = 1.;
  P.at(1,0) = std::sqrt(3.)*cphi;
  P.at(1,1) = std::sqrt(3.)*sphi;

  // Iteration
  for (unsigned int n=2; n<_degree+3; n++)
  {

    double sqrt2np1 = std::sqrt( 2.*(double)n + 1. );
    double sqrt2nm1 = std::sqrt( 2.*(double)n - 1. );
    double sqrt2nm3 = std::sqrt( 2.*(double)n - 3. );

    for (unsigned int m=0; m<n+1; m++)
    {

      if(n==m)
	P.at(n,n) = sqrt2np1/std::sqrt(2.*(double)n)*sphi*P.at(n-1, n-1);
      else if(m==0)
	P.at(n,m) = (sqrt2np1/(double)n)*(sqrt2nm1*cphi*P.at(n-1,m) - ((double)n-1.)/std::sqrt(2.*(double)n-3.)*P.at(n-2,m));
      else
	P.at(n,m) = sqrt2np1/(std::sqrt( (double)n + (double)m)*std::sqrt((double)n-(double)m))*(sqrt2nm1*cphi*P.at(n-1,m) - std::sqrt((double)n+(double)m-1.)*std::sqrt((double)n-(double)m-1.)/sqrt2nm3*P.at(n-2,m) );

    }

  }

  return P;

}

template<class T>
DACE::AlgebraicVector<T> GravityFieldModel::gravPCPF(const DACE::AlgebraicVector<T>& rr, const T& rr_norm, const DACE::AlgebraicMatrix<T>& P, const DACE::AlgebraicVector<T>& smlambda, const DACE::AlgebraicVector<T>& cmlambda) const
{
  /*!internal function computing gravity in planet-centered
   * planet-fixed (PCEF) coordinates using PCPF position, desired
   * degree/order, normalized associated legendre polynomials, normalized
   * spherical harmonic coefficients, trigonometric functions of geocentric
   * latitude and longitude, planetary constants, and radius to center of
   * planet. Units are MKS.
   */

  using  std::sqrt;
  //using DACE::sqrt;
  using DACE::sqr;

  T rRatio   = _RADIUS/rr_norm;
  T rRatio_n = rRatio;

  // initialize summation of gravity in radial coordinates
  T dUdrSumN      = 1.;
  T dUdphiSumN    = 0.;
  T dUdlambdaSumN = 0.;

  T sqrp_xy = rr[0]*rr[0] + rr[1]*rr[1];
  T isrp_xy = 1./sqrp_xy;
  T p_xy = sqrt(sqrp_xy);
  T ip_xy = 1./p_xy;


  // summation of gravity in radial coordinates
  for (unsigned int n=2; n<_degree+1; n++)
  {
    rRatio_n      = rRatio_n*rRatio;

    T dUdrSumM      = 0.;
    T dUdphiSumM    = 0.;
    T dUdlambdaSumM = 0.;

    for (unsigned int m=0; m<n+1; m++)
    {

      dUdrSumM      = dUdrSumM + P.at(n,m)*(_C.at(n,m)*cmlambda[m] + _S.at(n,m)*smlambda[m]);
      dUdphiSumM    = dUdphiSumM + ( (P.at(n,m+1)*_LSF.at(n,m)) - rr[2]*ip_xy*m*P.at(n,m))*(_C.at(n,m)*cmlambda[m] + _S.at(n,m)*smlambda[m]);
      dUdlambdaSumM = dUdlambdaSumM + m*P.at(n,m)*(_S.at(n,m)*cmlambda[m] - _C.at(n,m)*smlambda[m]);

    }

    dUdrSumN      = dUdrSumN      + dUdrSumM*rRatio_n*((double)n+1.);
    dUdphiSumN    = dUdphiSumN    + dUdphiSumM*rRatio_n;
    dUdlambdaSumN = dUdlambdaSumN + dUdlambdaSumM*rRatio_n;

  }


  // // gravity in spherical coordinates
  // T dUdr      = -1.*_GM/(rr_norm*rr_norm)*dUdrSumN;
  // T dUdphi    =     _GM/rr_norm*dUdphiSumN;
  // T dUdlambda =     _GM/rr_norm*dUdlambdaSumN;

  // gravity in ECEF coordinates
  // DACE::AlgebraicVector<T> g(3);
  // g[0] = ((1./rr_norm)*dUdr - (rr[2]/(rr_norm*rr_norm*p_xy))*dUdphi)*rr[0] - (dUdlambda/(rr[0]*rr[0] + rr[1]*rr[1]))*rr[1];
  // g[1] = ((1./rr_norm)*dUdr - (rr[2]/(rr_norm*rr_norm*p_xy))*dUdphi)*rr[1] + (dUdlambda/(rr[0]*rr[0] + rr[1]*rr[1]))*rr[0];
  // g[2] = (1./rr_norm)*dUdr*rr[2] + p_xy/(rr_norm*rr_norm)*dUdphi;

  // gravity in spherical coordinates
  T dUdr      = -1.*_GM*dUdrSumN;
  T dUdphi    =     _GM*dUdphiSumN;
  T dUdlambda =     _GM*dUdlambdaSumN;


  T r = sqr(rr[0]) + sqr(rr[1]) + sqr(rr[2]);
  T ir = 1./sqrt(r);
  r = 1./r;

  DACE::AlgebraicVector<T> g(3);
  g[0] = ( (dUdr*r) - ( (rr[2]*ip_xy) * (dUdphi*r) ) ) * (rr[0]*ir) - ( dUdlambda*isrp_xy ) * (rr[1]*ir);
  g[1] = ( (dUdr*r) - ( (rr[2]*ip_xy) * (dUdphi*r) ) ) * (rr[1]*ir) + ( dUdlambda*isrp_xy ) * (rr[0]*ir);
  g[2] = (dUdr*r)*(rr[2]*ir) + (dUdphi*r)*(p_xy*ir);

  return g;

}

template<class T>
DACE::AlgebraicVector<T> GravityFieldModel::compute(const DACE::AlgebraicVector<T>& rr) const
{

  using  std::asin;
  //using DACE::asin;
  using  std::sin;
  //using DACE::sin;
  using  std::cos;
  //using DACE::cos;

  T rr_norm = vnorm(rr);
  if (DACE::cons(rr_norm)<_RADIUS)
    throw std::runtime_error("SF::GravityFieldModel<T>::compute: Radial position is less than equatorial radius of planetary model.");

  // Compute geocentric latitude
  T phic = asin( rr[2]/ rr_norm );

  // Compute lambda
  T lambda = atan2_mod( rr[1], rr[0] );

  DACE::AlgebraicVector<T> smlambda(_degree+1), cmlambda(_degree+1);


  T slambda = sin(lambda);
  T clambda = cos(lambda);

  smlambda.at(0) = 0.;
  cmlambda.at(0) = 1.;
  smlambda.at(1) = slambda;
  cmlambda.at(1) = clambda;


  for (unsigned int m=2; m<_degree+1; m++)
  {
    smlambda[m] = 2.0 *clambda*smlambda[m-1] - smlambda[m-2];
    cmlambda[m] = 2.0 *clambda*cmlambda[m-1] - cmlambda[m-2];
  }


  // Compute normalized associated legendre polynomials
  DACE::AlgebraicMatrix<T> P = gravLegendre( phic );

  DACE::AlgebraicVector<T> aa = gravPCPF(rr, rr_norm, P, smlambda, cmlambda);

  return aa;
}

}

#endif
