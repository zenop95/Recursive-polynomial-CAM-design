//class necessary to link DA library to optimisation library
#ifndef DlibClass_H_
#define DlibClass_H_

#include<dace/dace.h>

//DLIB
#include <sys/stat.h>
#include <dlib/optimization.h>


typedef dlib::matrix<double,0,1> column_vector;

class obj_J
{
    
public:
    typedef ::column_vector column_vector;
    typedef dlib::matrix<double> general_matrix;
    DACE::AlgebraicVector<DACE::DA> Jfun;
    
    //Member functions//
    obj_J(DACE::AlgebraicVector<DACE::DA> JJ) {
        Jfun = JJ;
    }
    
    column_vector evaluate(const column_vector& x) const {
        unsigned int s = x.nr();
        DACE::AlgebraicVector<double> X(s,0.0);
        for(unsigned int i=0; i<s; i++) {
            X[i] = x(i);
        }
        
        DACE::AlgebraicVector<double> OUT = Jfun.eval(X);
        unsigned int S = OUT.size();
        column_vector OUTvec(S);
        for(unsigned int i=0; i<S; i++) {
            OUTvec(i) = OUT[i];
        }
        
        return OUTvec;
    }
    
    column_vector operator() (const column_vector& x) const {
    	return this->evaluate(x); 
    }
    
    void setJfun (DACE::AlgebraicVector<DACE::DA> J){
        
        Jfun = J;
        
    }
};

#endif
