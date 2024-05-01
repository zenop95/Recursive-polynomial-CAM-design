#ifndef dynorb_DYNCLASS_H_
#define dynorb_DYNCLASS_H_


#include <dace/dace.h>

// CLASSES

template<typename T, typename U>
class dynamicsTemplateTime
{
public:
	dynamicsTemplateTime() {};
	virtual DACE::AlgebraicVector<T> evaluate(const DACE::AlgebraicVector< T >& x, U t) = 0;
	virtual bool checkEvent()
	{
		return false;
	}
};

template<typename T>
class dynamics : public dynamicsTemplateTime<T, double>
{
public:
	dynamics() {};
};

template<typename T, typename U>
class dynamicsTemplateTimeScaled
{
public:
	dynamicsTemplateTimeScaled() {};
	virtual DACE::AlgebraicVector<T> evaluate(const DACE::AlgebraicVector< T >& x, const DACE::AlgebraicVector< T >& u, U t, U aMax, U Lsc, bool flagRtn, bool gravOrd) = 0;
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

template<typename T, typename U>
class dynamicsTemplateTimeAcc
{
public:
	dynamicsTemplateTimeAcc() {};
	virtual DACE::AlgebraicVector<T> evaluate(const DACE::AlgebraicVector< T >& x, const DACE::AlgebraicVector< T >& u, U t, U aMax) = 0;
	virtual bool checkEvent()
	{
		return false;
	}
};

template<typename T>
class dynamicsAcc : public dynamicsTemplateTimeAcc<T, double>
{
public:
	dynamicsAcc() {};
};


template<typename T, typename U>
class dynamicsTemplateDv
{
public:
	dynamicsTemplateDv() {};
	virtual DACE::AlgebraicVector<T> evaluate(const DACE::AlgebraicVector< T >& x, U t, U Lsc, U musc, bool gravOrd) = 0;
	virtual bool checkEvent()
	{
		return false;
	}
};

template<typename T>
class dynamicsDv : public dynamicsTemplateDv<T, double>
{
public:
	dynamicsDv() {};
};

#endif