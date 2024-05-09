#ifndef dynorb_DYNCLASS_H_
#define dynorb_DYNCLASS_H_


#include <dace/dace.h>

// CLASSES

template<typename T, typename U>
class dynamicsTemplateTimeScaled
{
public:
	dynamicsTemplateTimeScaled() {};
	virtual DACE::AlgebraicVector<T> evaluate(const DACE::AlgebraicVector< T >& x, const DACE::AlgebraicVector< T >& u, U t, U aMax, U Lsc, bool flagRtn, int gravOrd) = 0;
	virtual bool checkEvent()
	{
		return false;
	}
};

template<typename T>
class dynamicsScaled : public dynamicsTemplateTimeScaled<T, double>
{
public:
	dynamicsScaled() {};
};


#endif