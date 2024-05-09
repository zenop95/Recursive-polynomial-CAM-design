#ifndef dynorb_AIDAWRAPPERS_H_
#define dynorb_AIDAWRAPPERS_H_

#include <dace/dace.h>
#include "AIDA.h"
#include "astro/AstroCnvRef.h"
#include "DynClass.h"
#include "astro/AstroRoutines.h"

template<typename T>
class AIDAScaledDynamics : public dynamicsScaled<T> //includes scaled integration and input control acceleration
{
public:
	AIDAScaledDynamics(const std::string& gravmodel_name, const unsigned int gravmodel_order, const int AIDA_flags[3],
		const T& Bfactor, const T& SRPC)
		: m_aidaProp(gravmodel_name, gravmodel_order, AIDA_flags, 78.0), m_eventflag(false)

	{
		//Set uncertain parameters
		m_aidaProp.setUncParam("SRPC", SRPC);
		m_aidaProp.setUncParam("Bfactor", Bfactor);
	}

	DACE::AlgebraicVector<T> evaluate(const DACE::AlgebraicVector<T>& x, const DACE::AlgebraicVector<T>& u, double t, double aMax,  double Lsc, bool flagRtn, int gravOrd)
	{
 		const double mu   = 398600.4418;
        const double musc = 1.0;
        const double Vsc  = sqrt(mu/Lsc);
        const double Tsc  = Lsc/Vsc;
        const double Asc  = Vsc/Tsc;

        DACE::AlgebraicVector<T> xsc(6), uEci(3), usc(6), acc(3);
		if (flagRtn == 1){
			DACE::AlgebraicMatrix<T> dcm(3,3);
			dcm = astro::rtn2eci(x);
			uEci = dcm*u;
		}
		else {
			uEci = u;
		}	
		if (gravOrd == 0)
		{
			acc = keplerPropAcc(x,uEci*aMax,t,musc);
		}
		else if (gravOrd == 2)
		{
			acc = J2dynamics(x,uEci*aMax,t,musc,Lsc);
		}
		else if (gravOrd == 4)
		{
			acc = J2_J4dynamics(x,uEci*aMax,t,musc,Lsc);
		}
		else {
			for (int i = 0; i < 3; i ++) {
				xsc[i]   = x[i]*Lsc;
				xsc[i+3] = x[i+3]*Vsc;
				usc[i]   = 0.0;
				usc[i+3] = aMax*Asc*uEci[i];
			}
			double tsc = t*Tsc; 	
			acc =  m_aidaProp.evaluation(tsc, xsc, m_eventflag) + usc;
			for (int i = 0; i < 3; i ++) 
			{
				acc[i] = acc[i]/Vsc;
				acc[i+3] = acc[i+3]/Asc;
			}
		}
        return acc;
	}

	bool checkEvent()
	{
		return m_eventflag;
	}
private:
	dynorb::AIDA<T> m_aidaProp;
	bool m_eventflag;
};


#endif