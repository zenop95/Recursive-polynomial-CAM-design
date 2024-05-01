#ifndef atmos_math_NRLMSISE00_MATH_H
#define atmos_math_NRLMSISE00_MATH_H

// Suppress warning: comparing floating point with == or != is unsafe [-Werror=float-equal]
// This is OK here because this is the NRLMSISE-00 algorithm. Is applied to the WHOLE file.
#ifdef __clang__
# pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined __GNUC__
# pragma GCC   diagnostic ignored "-Wfloat-equal"
#endif

#include <cmath>
#include <stdexcept>
#include <dace/dace.h>

namespace atmos
{
    /**
     * @brief Namespace that contains mathematical utilities
     *
     */
    namespace math
    {
        /**
         * @brief Integrate cubic spline function from xa(1) to x
         *
         * @param xa Array of interpolation nodes
         * @param ya Arrays of tabulated function in ascending order by xa with ya = f(xa)
         * @param y2a Array of second derivatives
         * @param n Array size
         * @param x Ascissa endpoint of integration
         * @return double integrated value
         */
        template <class T>
        T splini(const T *xa, const T *ya, const T *y2a, const int n, const T& x);
        /**
         * @brief Calculate cubic spline interp value
         *
         * @param xa Array of interpolation nodes
         * @param ya Arrays of tabulated function values in ascending xa order
         * @param y2a Arrays of second derivatives
         * @param n Array size
         * @param x Abscissa of interpolation
         * @return double cubic spline interp value
         */
        template <class T>
        T splint(const T *xa, const T *ya, const T *y2a, const int n, const T& x);
        /**
         * @brief Calculate 2nd derivatives of cubic spline interp function
         *
         * @param x Array of interpolation nodes
         * @param y Arrays of tabulated function in ascending order by x with y = f(x)
         * @param n Array size
         * @param yp1 Specified derivatives at x(1)
         * @param ypn Specified derivatives at x(n)
         * @param y2 Output array of second derivatives
         */
        template <class T>
        void spline(const T *x, const T *y, const int n, const T& yp1, const T& ypn, T *y2);

    } // namespace math
} // namespace atmos

// Implementation
namespace atmos
{
    namespace math
    {
        template <class T>
        T splini(const T *xa, const T *ya, const T *y2a, const int n, const T& x)
        {
            T yi=0;
            int klo(0), khi(1);
            T xx, h, a, b, a2, b2;
            while ((DACE::cons(x)>DACE::cons(xa[klo])) and (khi<n))
            {
                xx=x;
                if (khi<(n-1))
                {
                    if (DACE::cons(x)<DACE::cons(xa[khi]))
                    {
                        xx=x;
                    }
                    else
                    {
                        xx=xa[khi];
                    }
                }
                h = xa[khi] - xa[klo];
                a = (xa[khi] - xx)/h;
                b = (xx - xa[klo])/h;
                a2 = a*a;
                b2 = b*b;
                yi += ((1.0 - a2) * ya[klo] / 2.0 + b2 * ya[khi] / 2.0 + ((-(1.0+a2*a2)/4.0 + a2/2.0) * y2a[klo] + (b2*b2/4.0 - b2/2.0) * y2a[khi]) * h * h / 6.0) * h;
                klo++;
                khi++;
            }
            return yi;
        }
        template<typename T>
        T splint(const T *xa, const T *ya, const T *y2a, const int n, const T& x)
        {
            int klo(0), khi=(n-1);
            int k;
            T h;
            T a, b;
            while ((khi-klo)>1)
            {
                k=(khi+klo)/2;
                if (DACE::cons(xa[k])>DACE::cons(x))
                {
                    khi=k;
                }
                else
                {
                    klo=k;
                }
            }
            h = xa[khi] - xa[klo];
            if (DACE::cons(h)==0.0)
            {
                throw std::runtime_error("CNrlmsise00_p::split: bad XA input to splint");
            }
            a = (xa[khi] - x)/h;
            b = (x - xa[klo])/h;
            return (a * ya[klo] + b * ya[khi] + ((a*a*a - a) * y2a[klo] + (b*b*b - b) * y2a[khi]) * h * h/6.0);
        }

        template<typename T>
        void spline(const T *x, const T *y, const int n, const T& yp1, const T& ypn, T *y2)
        {
            T* u = new T[n];
            T sig, p, qn, un;
            int i, k;

            if (DACE::cons(yp1)>0.99E30)
            {
                y2[0]=0;
                u[0]=0;
            }
            else
            {
                y2[0]=-0.5;
                u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
            }

            for (i=1;i<(n-1);i++)
            {
                sig = (x[i]-x[i-1])/(x[i+1] - x[i-1]);
                p = sig * y2[i-1] + 2.0;
                y2[i] = (sig - 1.0) / p;
                u[i] = (6.0 * ((y[i+1] - y[i])/(x[i+1] - x[i]) -(y[i] - y[i-1]) / (x[i] - x[i-1]))/(x[i+1] - x[i-1]) - sig * u[i-1])/p;
            }

            if (DACE::cons(ypn)>0.99E30)
            {
                qn = 0;
                un = 0;
            }
            else
            {
                qn = 0.5;
                un = (3.0 / (x[n-1] - x[n-2])) * (ypn - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]));
            }

            y2[n-1] = (un - qn * u[n-2]) / (qn * y2[n-2] + 1.0);

            for (k=n-2;k>=0;k--)
            {
                y2[k] = y2[k] * y2[k+1] + u[k];
            }

            delete[] u;
        }
    } // namespace math
} // namespace atmos

#endif // atmos_math_NRLMSISE00_MATH_H