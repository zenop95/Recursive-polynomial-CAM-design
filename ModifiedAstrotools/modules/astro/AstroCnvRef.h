/********************************************************************************************/
//  AstroCnvRef.h                                                                            /
//  Astrodynamic function for the Conversion of reference frame from ane the each other      /
//                                                                                           /
//                                                                                           /
//  Created by Daniele Antonio Santeramo on 28/11/16.                                        /
//                                                                                           /
/********************************************************************************************/

#ifndef ASTROCNVREF_H_INCLUDED_
#define ASTROCNVREF_H_INCLUDED_

#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>


template <typename T> T checkatan2(T angle, T a, T b) {
    if ( DACE::cons(a) > 0 && DACE::cons(b) < 0 ) return angle + M_PI;
    else if ( DACE::cons(a) < 0 && DACE::cons(b) < 0 ) return angle - M_PI;
    else return angle;
}

template <typename T> T atan2_mod(T a, T b) {
    T angle = atan(a/b);
    checkatan2(angle, a, b);
    return angle;
}

template <typename T> T atan3(T a, T b) {
  //a and b are in radiant
    double eps0 = 1e-16;
    T c;
    if ( std::abs(DACE::cons(a)) < eps0 ) return atan(a/b) + (1.0 - DACE::cons(b)/std::abs(DACE::cons(b)))*0.5*M_PI;
    else c = (2.0 - DACE::cons(a)/std::abs(DACE::cons(a)))*0.5*M_PI;

    if ( std::abs(DACE::cons(b)) < eps0 ) return atan(a/b) + (1.0 - DACE::cons(a)/std::abs(DACE::cons(a)))*M_PI;
    else return c + DACE::cons(a)/std::abs(DACE::cons(a))*DACE::cons(b)/std::abs(DACE::cons(b))*(atan(a/b)*DACE::cons(atan(a/b))/std::abs(DACE::cons(atan(a/b))) - 0.5*M_PI) ;

}

template <typename T> T Norma(DACE::AlgebraicVector<T> obj) {
  /*Alternative version of vnorm*/
  T temp=0;
  for ( unsigned int i = 0; i < obj.size(); ++i) {
      temp += obj[i]*obj[i];
  }

  temp = sqrt(temp);

  return temp;
}

template <typename T> DACE::AlgebraicVector<T> mee2eci(DACE::AlgebraicVector<T> mee, const double mu) {
    /*function to convert Modified Equinoctial Elements into Earth-Centred Inertial elements
    !> mee = {p
       f=|e|*cos(RA+PA)
       g=|e|*sin(RA+PA)
       h=tan(i/2)*cos(RA)
       k=tan(i/2)*sin(RA)
       L=TA+RA+PA}
       RA: rigth ascension of ascending node; PA: argument of periapsis; TA: true anomaly; i:orbital inclination
   !< return AlgebraicVector ECI element res = {x, y, z, dx, dy, dz}
      */
    //const double mu=398600.0;  //km^3/s^2
    T dum, e;
    T p=mee[0];
    T f=mee[1];
    T g=mee[2];
    T h=mee[3];
    T k=mee[4];
    T L=mee[5];

    if (cons(f) < 1e-8 || cons(g) < 1e-8) {
      e = 0;
    }
    else {
      e = sqrt(f*f + g*g);
    }

    T H  = sqrt(p*mu);
    T rm = p/(1.0 + f*cos(L) + g*sin(L) );
    dum = 1/(1 + DACE::sqr(h)+DACE::sqr(k));

    DACE::AlgebraicVector<T> r(3),v(3);
    r[0] = rm*(cos(L) + (h*h - k*k)*cos(L) + 2.0*h*k*sin(L) )*dum;
    r[1] = rm*(sin(L) - (h*h - k*k)*sin(L) + 2.0*h*k*cos(L) )*dum;
    r[2] = 2.0*rm*(h*sin(L) - k*cos(L) )*dum;
    v[0] = -sqrt(mu/p)*(sin(L) + (h*h - k*k)*sin(L) - 2.0*h*k*cos(L) + g - 2.0*f*h*k + (h*h - k*k)*g )*dum;
    v[1] = -sqrt(mu/p)*(-cos(L) + (h*h - k*k)*cos(L) + 2.0*h*k*sin(L) - f + 2.0*g*h*k + (h*h - k*k)*f )*dum;
    v[2] = 2.0*sqrt(mu/p)*(h*cos(L) + k*sin(L) + f*h + k*g )*dum;

    DACE::AlgebraicVector<T> res(6);
    res[0]=r[0]; res[1]=r[1]; res[2]=r[2];
    res[3]=v[0]; res[4]=v[1]; res[5]=v[2];

    return res;

}

template <typename T> DACE::AlgebraicVector<T> eci2mee(DACE::AlgebraicVector<T> x, const double mu) {
  /*function to convert Earth-Centred Inertial into Modified Equinoctial Elements elements
  !> AlgebraicVector ECI element {x, y , z, dx, dy, dz}
  !< return AlgebraicVector MEE element res = {p, f, g, h, k, L}*/
    //const double mu=398600.0;

    DACE::AlgebraicVector<T> pos(3);
    DACE::AlgebraicVector<T> vel(3);
    pos.assign(x.begin(), x.begin()+3);
    vel.assign(x.begin()+3, x.end());
    T rm = Norma(pos);

    DACE::AlgebraicVector<T> fver(3), gver(3), Hver(3),posver(3),velver(3), ec(3), H(3);
    H = pos.cross(vel);
    Hver = H/Norma(H);

    T p = Norma(H)*Norma(H)/mu;

    posver = pos/rm;
    velver = (rm*vel - pos.dot(vel)*pos/rm)/Norma(H);

    ec = -posver + vel.cross(H)/mu;

    T k = Hver[0]/(1.0 + Hver[2]);
    T h = -Hver[1]/(1.0 + Hver[2]);

    T dum=1/(1+DACE::sqr(k)+DACE::sqr(h));

    T AA[3][3];
    AA[0][0] = dum*(1-DACE::sqr(k)+DACE::sqr(h)); AA[0][1] = dum*2*k*h;                         AA[0][2] = dum*2*k;
    AA[1][0] = dum*2*k*h;                         AA[1][1] = dum*(1+DACE::sqr(k)-DACE::sqr(h)); AA[1][2] = dum*(-2*h);
    AA[2][0] = dum*(-2*k);                        AA[2][1] = dum*2*h;                           AA[2][2] = dum*(1-DACE::sqr(k)-DACE::sqr(h));

    fver[0] = AA[0][0]; fver[1] = AA[1][0]; fver[2] = AA[2][0];
    gver[0] = AA[0][1]; gver[1] = AA[1][1]; gver[2] = AA[2][1];

    T f = ec.dot(fver);
    T g = ec.dot(gver);

    T cosL = posver[0] + velver[1];
    T sinL = posver[1] - velver[0];

    T L = atan(sinL/cosL);
    if ( DACE::cons(cosL) < 0.0 ) { L += M_PI; }
    if ( DACE::cons(L) < 0.0 ) { L += 2.0*M_PI; }

    DACE::AlgebraicVector<T> res(6);
    res[0] = p;
    res[1] = f;
    res[2] = g;
    res[3] = h;
    res[4] = k;
    res[5] = L;

    return res;

}

template <typename T> DACE::AlgebraicVector<T> Kep2mee(DACE::AlgebraicVector<T> v) {
  /*member function to convert keplerian  classical element into MEE reference frame element
  !> keplerian element v = {a, e , i, RA, PA, TA}
     OMEGA: rigth ascension of ascending node; omega: argument of periapsis; M: mean anomaly;i:orbital inclination
  !< return AlgebraicVector of MEE reference frame res = {p, f, g, h, k, L}*/

    T p = v[0]*(1.0 - DACE::sqr(v[1]) );
    T f = v[1]*cos(v[3] + v[4]);
    T g = v[1]*sin(v[3] + v[4]);
    T h = tan(v[2]/2.0)*cos(v[3]);
    T k = tan(v[2]/2.0)*sin(v[3]);

    T L = v[5] + v[3] + v[4];

    DACE::AlgebraicVector<T> res(6);
    res[0] = p;
    res[1] = f;
    res[2] = g;
    res[3] = h;
    res[4] = k;
    res[5] = L;

    return res;

}

template <typename T> DACE::AlgebraicVector<T> mee2Kep(DACE::AlgebraicVector<T> mee) {
  /*member function to convert MEE reference frame element into keplerian  classical element
  !> AlgebraicVector of MEE reference frame mee = {p, f, g, h, k, L}
  !< return AlgebraicVector of keplerian element res = {a, e , i, RA, PA, TA}
    RA: rigth ascension of ascending node; PA: argument of periapsis; TA: true anomaly; i:orbital inclination*/

  using DACE::cons;

  T p=mee[0];
  T f=mee[1];
  T g=mee[2];
  T h=mee[3];
  T k=mee[4];
  T L=mee[5];

  T e = sqrt(DACE::sqr(f) + DACE::sqr(g) );
  T a = p/(1.0 - DACE::sqr(e) );
  T i = 2.0*atan(sqrt(DACE::sqr(h) + DACE::sqr(k)) );
  i = i - floor(cons(i)/(2.0*M_PI))*2.0*M_PI;

  T RA = atan(k/h);
    if ( DACE::cons(h) < 0.0 ) { RA += M_PI; }
    if ( DACE::cons(RA) < 0.0 ) { RA += 2.0*M_PI; }

  T RAPA = atan(g/f);
    if ( DACE::cons(f) < 0.0 ) { RAPA += M_PI; }
    if ( DACE::cons(RAPA) < 0.0 ) { RAPA += 2.0*M_PI; }

  T PA = RAPA - RA; // periapsis argument
  T TA = L - RAPA; // true anomaly

  DACE::AlgebraicVector<T> res(6);
  res[0] = a;
  res[1] = e;
  res[2] = i;
  res[3] = RA;
  res[4] = PA;
  res[5] = TA;

  return res;
}

template <typename T> DACE::AlgebraicVector<T> Kep2eci(DACE::AlgebraicVector<T> v, const double mu) {
  /*member function to convert keplerian  classical element into Earth-Centred inertial reference frame element
  !< keplerian element v = {a, e , i, RA, PA, TA}
     RA: rigth ascension of ascending node; PA: argument of periapsis; TA: true anomaly; i:orbital inclination
  !> return AlgebraicVector of ECI reference frame res = {x, y, z, dx, dy, dz}*/
    //const double mu = 398600.0;

    T p = v[0]*(1.0 - v[1]*v[1] );
    DACE::AlgebraicVector<T> rm(3), vm(3); // position and velocity in perifocal refererence frame
    rm[0] = p*cos(v[5])/(1.0 + v[1]*cos(v[5]) );
    rm[1] = p*sin(v[5])/(1.0 + v[1]*cos(v[5]) );
    rm[2] = 0.0;
    vm[0] = -sin(v[5])*sqrt(mu/p);
    vm[1] = (v[1] + cos(v[5]))*sqrt(mu/p);
    vm[2] = 0.0;

    T cRA = cos(v[3]);  T sRA = sin(v[3]);
    T cPA = cos(v[4]);  T sPA = sin(v[4]);
    T ci = cos(v[2]);  T si = sin(v[2]);

    T RR[3][3]; // rotational matrix from perifocal to eci reference frame
    RR[0][0] = cRA*cPA-sRA*ci*sPA;  RR[0][1] = -cRA*sPA-sRA*ci*cPA; RR[0][2] = sRA*si;
    RR[1][0] = sRA*cPA+cRA*ci*sPA;  RR[1][1] = -sRA*sPA+cRA*ci*cPA; RR[1][2] = -cRA*si;
    RR[2][0] = si*sPA;              RR[2][1] = si*cPA;              RR[2][2] = ci;

    DACE::AlgebraicVector<T> rr(3),vv(3);
    for(unsigned int i=0;i<3;i++){
        rr[i]=0.0;
        vv[i]=0.0;
        for(unsigned int j=0;j<3;j++){
            rr[i]=rr[i]+RR[i][j]*rm[j];
            vv[i]=vv[i]+RR[i][j]*vm[j];
        }
    }

    DACE::AlgebraicVector<T> res(6);
    res[0]=rr[0]; res[1]=rr[1]; res[2]=rr[2];
    res[3]=vv[0]; res[4]=vv[1]; res[5]=vv[2];

    return res;
}

template <typename T> DACE::AlgebraicVector<T> eci2Kep(DACE::AlgebraicVector<T> x, const double mu) {
  /*member function to convert ECI state vector into keplerian classical element
  !< AlgebraicVector of ECI reference frame  {x, y, z, dx, dy, dz}
  !> return AlgebraicVector of keplerian element res = {a, e , i, RA, PA, TA}
     RA: rigth ascension of ascending node; PA: argument of periapsis; TA: true anomaly; i: orbital inclination*/
    using DACE::cons;
    //const double mu=398600.0;

    DACE::AlgebraicVector<T> pos(3);
    DACE::AlgebraicVector<T> vel(3);
    pos.assign(x.begin(), x.begin()+3);
    vel.assign(x.begin()+3, x.end());
    
    T rm = Norma(pos);
    T vm = Norma(vel);

    DACE::AlgebraicVector<T> H(3), ec(3), Hver(3);

    H = pos.cross(vel);
    T Hm = Norma(H);
    Hver = H/Hm;

    DACE::AlgebraicVector<T> nvec(3);
    nvec[0] = -H[1];  nvec[1] = H[0];  nvec[2] = 0.0;
    T nvecm = Norma(nvec);
    T Ecl= vm*vm - mu/rm;
    T posdotvel = pos.dot(vel);

    ec =  (Ecl*pos - posdotvel*vel)/mu;
    T ecm = Norma(ec);

    T Emec = vm*vm*0.5 - mu/rm;
    T a;
    if (std::abs(cons(Emec)) > 1e-12) a = -mu/(2.0*Emec);
    else a = 0.0;

    T Hk = H[2]/Hm;
    T in = acos(Hk);

    T RA;
    if (cons(nvecm) > 1e-12) {
        T temp = nvec[0]/nvecm;
        if (DACE::cons(temp) > 1.0) temp = DACE::cons(temp)/std::abs(DACE::cons(temp));
        RA = acos(temp);
        if (cons(nvec[1]) > 0.0) RA = 2.0*M_PI - RA;
    }
    else RA = 0.0;

    T PA;
    if (cons(in) > 1e-12 && cons(ecm) > 1e-12) {
        PA = acos(nvec.dot(ec)/(Norma(nvec)*Norma(ec)) ); //periapsis argument
        if ( cons(ec[2])< 0.0 ) PA = 2.0*M_PI - PA;
    }
    else PA = 0.0;

    T TA;
    if (cons(ecm) > 1e-12) {
      TA = acos(ec.dot(pos)/(Norma(ec)*Norma(pos)) ); //true anomaly
      if (cons(posdotvel ) < 0 ) TA = 2.0*M_PI - TA;
    }

    // ec = -pos/rm + vel.cross(H)/mu;
    //
    // T a = 1.0/(2.0/rm - vm*vm/mu);
    //
    // T h = H[0]/(1.0+H[2]);
    // T k =-H[1]/(1.0+H[2]);
    //
    // T in = 2.0*atan(sqrt(h*h+k*k));//orbit inclination
    //
    // T dum = 1.0 / (1.0 + h * h + k * k);
    // DACE::AlgebraicVector<T> fver(3), gver(3);
    // fver[0] = dum*(1.0 - h*h + k*k);    gver[0] = dum*2.0*h*k;
    // fver[1] = dum*2.0*h*k;              gver[1] = dum*dum*(1.0 + h*h - k*k);
    // fver[2] = -dum*2.0*h;               gver[2] = dum*2.0*k;
    //
    // T consth = ec.dot(gver);
    // T xk = ec.dot(fver);
    // T x1 = pos.dot(fver);
    // T y1 = pos.dot(gver);
    // T ecm = sqrt(consth*consth + xk*xk);
    //
    // T RA;
    // if (cons(in) > 1e-12) {
    //     RA = atan2_mod(h, k);
    //     RA += -trunc(RA/2.0/M_PI)*2.0*M_PI;
    // }
    // else RA = 0.0;
    //
    // T PA;
    // T auxangle = atan2_mod(consth, xk);
    // // if (cons(consth) < 0.0) auxangle += 2.0*M_PI;
    // auxangle += -trunc(auxangle/(2.0*M_PI))*2.0*M_PI;
    // if (cons(ecm) > 1e-12) {
    //     PA = auxangle - RA - (trunc((auxangle - RA)/(2.0*M_PI)))*2.0*M_PI;
    //     // if (cons(PA) < 0.0) PA += 2.0*M_PI;
    // }
    // else PA = 0.0;
    //
    // T Lt = atan2_mod(y1,x1);
    // // if (cons(y1) < 0.0) Lt += 2.0*M_PI;
    // Lt += -trunc(Lt/(2.0*M_PI))*2.0*M_PI;
    // T TA = Lt - RA -PA - (trunc((Lt - RA -PA)/(2.0*M_PI)))*2.0*M_PI;
    // // if (cons(TA) < 0.0) TA += 2.0*M_PI;

    DACE::AlgebraicVector<T> res(6);
    res[0] = a;
    res[1] = ecm;
    res[2] = in;
    res[3] = RA;
    res[4] = PA;
    res[5] = TA;

    // DACE::AlgebraicVector<double> ek(3); ek[0]= 0.0;ek[1]= 0.0;ek[2]= 1.0;
    // DACE::AlgebraicVector<T> N = ek.cross(Hver);
    //
    // T RA = atan2_mod(N[1], N[0] ); //right ascension
    // // if (abs(cons(RA)) < 1e-9) RA = RA - cons(RA);
    // RA = RA - floor(cons(RA)/(2.0*M_PI))*2.0*M_PI;
    //
    // T PA = acos(N.dot(ec)/(Norma(N)*Norma(ec)) ); //periapsis argument
    // if ( cons(ec)[2] < 0 ) { PA = 2.0*M_PI - PA; }
    //
    // T TA = acos(ec.dot(pos)/(Norma(ec)*Norma(pos)) ); //true anomaly
    // if (cons(pos.dot(vel) ) < 0 ) {TA = 2.0*M_PI - TA; }
    //
    // DACE::AlgebraicVector<T> res(6);
    // res[0] = a;
    // res[1] = Norma(ec);
    // res[2] = in;
    // res[3] = RA;
    // res[4] = PA;
    // res[5] = TA;

    return res;
  }

template <typename T> DACE::AlgebraicVector<T> rv2TopRaDec(DACE::AlgebraicVector<T> x, DACE::AlgebraicVector<double> Rsite, DACE::AlgebraicVector<double> Vsite) {

    using DACE::sqr;
    const double trhes = 1e-12;
    DACE::AlgebraicVector<T> r(3), v(3), rho, rhodot;
    r.assign(x.begin(), x.begin()+3);
    v.assign(x.begin()+3, x.end());
    //std::cout << r.cons() << std::endl  << v.cons() << std::endl  << Rsite << std::endl<< Vsite << std::endl;
    rho = r - Rsite;
    rhodot = v - Vsite;
    T rhonorm = rho.vnorm();

    T delta = asin(rho[2]/rhonorm);
    T alpha;

    if ( sqrt( sqr(DACE::cons(rho[0]) ) + sqr(DACE::cons(rho[1]) ) ) > trhes ) {
        alpha = atan(rho[1]/rho[0]);
        if ( DACE::cons(rho[0]) < 0.0 ) { alpha += M_PI; }
        if ( DACE::cons(alpha) < 0.0 ) { alpha += 2.0*M_PI; }
        // if (DACE::cons(rho[1]) < 0.0) alpha += 2.0*M_PI;
        // alpha = checkatan2(alpha, rho[1],rho[0]);
    }
    else {
        alpha = atan(rhodot[1]/rhodot[0]);
        if ( DACE::cons(rhodot[0]) < 0.0 ) { alpha += M_PI; }
        if ( DACE::cons(alpha) < 0.0 ) { alpha += 2.0*M_PI; }
        
        // if (DACE::cons(rhodot[1]) < 0.0) alpha += 2.0*M_PI;
        // alpha = checkatan2(alpha, rhodot[1],rhodot[0]);
    }

    DACE::AlgebraicVector<T> res(2);
    res[0] = alpha; res[1] = delta;

    return res;
}

template<typename T> DACE::AlgebraicVector<T> Cart2Pol(DACE::AlgebraicVector<T> X){


    DACE::AlgebraicVector<T> r, v;
    r.assign(X.begin(), X.begin()+3);
    v.assign(X.begin()+3, X.end());




    DACE::AlgebraicVector<T> rad(6);
    double small = 1E-8;


    rad[0] = r.vnorm();
    T temp = sqrt(r[0]*r[0] + r[1]*r[1]);


    if (DACE::cons(temp) < small){
        rad[1] = atan(v[1]/v[0]);
        if (DACE::cons(v[0])<0) rad[1] += M_PI;
    }
    else{
        rad[1] = atan(r[1]/r[0]);
        if (DACE::cons(r[0])<0) rad[1] += M_PI;
    }
    if (DACE::cons(rad[1]) < 0.0 ) rad[1] += 2.0*M_PI;


    rad[2] = asin(r[2]/rad[0]);


    T temp1 = -r[1]*r[1] - r[0]*r[0];
    rad[3] = dot(r,v)/(rad[0]);


    if (DACE::cons(DACE::abs(temp1)) > small)
        rad[4] = ( v[0]*r[1] - v[1]*r[0] ) / temp1;
    else
        rad[4] = 0.0;


    T ddecl;
    if (DACE::cons(DACE::abs(temp)) > small)
        rad[5] = ( v[2] - rad[3]*sin( rad[2] ) ) / temp;
    else
        rad[5] = 0.0;


    return rad;
}

template <typename T> DACE::AlgebraicVector<T> eci2CU(DACE::AlgebraicVector<T> x, const double mu, const double RU, std::string selection = "both") {

    const double TU = (std::sqrt(std::pow(RU,3)/mu));
    const double VU = RU/TU;
    if (selection == "pos") return x/RU;
    else if (selection == "vel") return x/VU;
    else {
        DACE::AlgebraicVector<T> pos(3);
        DACE::AlgebraicVector<T> vel(3);
        pos.assign(x.begin(), x.begin()+3);
        vel.assign(x.begin()+3, x.end());
        DACE::AlgebraicVector<T> res(6);
        res[0] = pos[0]/RU;
        res[1] = pos[1]/RU;
        res[2] = pos[2]/RU;
        res[3] = vel[0]/VU;
        res[4] = vel[1]/VU;
        res[5] = vel[2]/VU;
        return res;
    }

}

std::vector< DACE::AlgebraicVector<double> > vertices(std::vector<double> v = {-1.0, 1.0}, unsigned int k = DACE::DA::getMaxVariables() ) {
  /*this function compute the combination with repetition of vector v into a vector with k element
  it compute all of iper-Domain vertex in accord with the number of variables of that one*/

    std::vector< DACE::AlgebraicVector<double> > bounds(1);

    unsigned int count = 0;
    while ( count < k) {
        unsigned int size = bounds.size();
        for (unsigned int i = 0; i < size; ++i) {
            for (unsigned j = 0; j < v.size(); ++j) {
                DACE::AlgebraicVector<double> temp;
                temp = bounds.at(i);
                temp.push_back(v[j]);
                bounds.push_back(temp);
            }
        }
        bounds.erase(bounds.begin(),bounds.begin()+size);
        ++count;
    }

    return bounds;

}

std::vector< DACE::AlgebraicVector<double> > GridCreation(const unsigned int N ) {
  /*this function create a grid of point equaly spaced in the Domain
   the grid is a ideally grid in squared unitary Domain*/
    unsigned int k = DACE::DA::getMaxVariables();

    if ( N < std::pow(2,k) ) {std::cout << "The minimum Grid number point created is: 512";}

    double dirN = std::ceil( std::pow(N, (1.0/k) ) );
    std::cout << "To space equally the Grid point it is created a Grid of: " << std::pow(dirN, k)<<std::endl;

    std::vector< DACE::AlgebraicVector<double> > bounds(1);
    std::vector<double> v;
    for ( unsigned int m = 0; m < dirN; ++m) {
        v.push_back( -1.0 + m*2.0/(dirN-1.0));
    }

    unsigned int count = 0;
    while ( count < k) {
        unsigned int size = bounds.size();
        for (unsigned int i = 0; i < size; ++i) {
            for (unsigned j = 0; j < v.size(); ++j) {
                DACE::AlgebraicVector<double> temp;
                temp = bounds.at(i);
                temp.push_back(v[j]);
                bounds.push_back(temp);
            }
        }
        bounds.erase(bounds.begin(),bounds.begin()+size);
        ++count;
    }

    return bounds;

}

template <class T> DACE::AlgebraicVector<T> propagateKepler(DACE::AlgebraicVector<T> x0, double t1, double t2, double mu) {
    // Propagates and expands Kepler motion with respect to both state and time
    // Works only for elliptic orbits around the Earth

    // Input arguments:
    // ---------------------------------------------------------------------
    // x0 = initial state vector [km, km/s]
    // t1   = propagation initial time [s]
    // t2   = propagation final time [s]

    // Output arguments:
    // ---------------------------------------------------------------------
    // xf = final state vector [km, km/s]
    using DACE::cons;
    const int size = x0.size();


    double tol;
    int iter;
    T h, r0, v0, a, p, sigma0, MmM0, EmE0, fx0, fxp, theta, r, F, G, Ft, Gt, NmN0, HmH0;
    DACE::AlgebraicVector<T> rr0(size/2), vv0(size/2);
    rr0[0] = x0[0]; rr0[1] = x0[1]; rr0[2] = x0[2];
    vv0[0] = x0[3]; vv0[1] = x0[4]; vv0[2] = x0[5];

    DACE::AlgebraicVector<T> hh = rr0.cross(vv0);
    h  = vnorm(hh);

    r0 = rr0.vnorm();
    v0 = vv0.vnorm();

    // useful parameters
    a = mu/(2*mu/r0 - v0*v0);
    p = h*h/mu;
    sigma0 = rr0.dot(vv0)/sqrt(mu);

    tol  = 1.0;
    iter = 0;

    MmM0   = (t2 - t1)*sqrt(mu/a/a/a);
    EmE0   = cons(MmM0);

    // find reference solution
    while((tol>1e-14) && (iter<3000)){

        iter = iter + 1;

        fx0 = -cons(MmM0)+EmE0+cons(sigma0)/sqrt(cons(a))*(1.0-cos(EmE0))-(1.0-cons(r0)/cons(a))*sin(EmE0);
        fxp = 1.0+cons(sigma0)/sqrt(cons(a))*sin(EmE0)-(1.0-cons(r0)/cons(a))*cos(EmE0);

        tol = abs(fx0/fxp);
        EmE0 = EmE0 - fx0/fxp;
    }
    // std::cout<<"after "<<iter<<" iteration, EmE0 is: "<<EmE0<<std::endl;

    // DA expansion of EmE0
    // iter = 1;
    // int order = DACE::DA::getTO();
    // while( iter <= order){
    //
    //     fx0 = -MmM0+EmE0+sigma0/sqrt(a)*(1.0-cos(EmE0))-(1.0-r0/a)*sin(EmE0);
    //     fxp = 1.0+sigma0/sqrt(a)*sin(EmE0)-(1.0-r0/a)*cos(EmE0);
    //
    //     EmE0 = EmE0 - fx0/fxp;
    //
    //     iter *= 2;
    // }

    theta = 2*atan(sqrt(a*p)*tan(EmE0/2)/( r0 + sigma0*sqrt(a)*tan(EmE0/2)));
    if ( DACE::cons( r0 + sigma0*sqrt(a)*tan(EmE0/2)) < 0.0 ) { theta += M_PI; }
    if ( DACE::cons(theta) < 0.0 ) { theta += 2.0*M_PI; }
    
    r     = p*r0/(r0+(p-r0)*cos(theta)-sqrt(p)*sigma0*sin(theta));

    // compute the Lagrangian coefficients
    F  = 1 - a/r0 * (1 - cos(EmE0));
    G  = a*sigma0/sqrt(mu)*(1 - cos(EmE0)) + r0*sqrt(a/mu)*sin(EmE0);
    Ft = -sqrt(mu*a)/(r*r0)*sin(EmE0);
    Gt = 1 - a/r*(1-cos(EmE0));

    // get the flow of the Kepler problem
    DACE::AlgebraicVector<T> xf(size);
    for (int i = 0; i < 3; i++)
    {
        xf[i]   = F*rr0[i] + G*vv0[i];
        xf[i+3] = Ft*rr0[i] + Gt*vv0[i];
    }

    return xf;
}

#endif //ASTROCNVREF_H_INCLUDED_
