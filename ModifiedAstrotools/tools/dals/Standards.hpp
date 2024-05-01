#ifndef Standards_H_
#define Standards_H_

//DACE
#include<dace/dace.h>

//Eigen
#include <Eigen/Dense>
#include <random>
//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/normal_distribution.hpp> 

//SPICE
#include <cspice/SpiceUsr.h>

#define       MAXWIN   1000
#define       TIMFMT   "DD MON YYYY HR:MN:SC.###### (UTC) ::UTC ::RND"
#define       TIMYYY   "YYYY.############::UTC"
#define       TIMLEN   41

#define       TIMFMT_SADA "JULIAND.#################::TDB"
#define       TIMLEN_SADA 26

#define       TIMFMT_RES "JULIAND.#################::UTC"
#define       TIMLEN_RES 26

/*
namespace Eigen {
	namespace internal {
		template<typename Scalar> 
		struct scalar_normal_dist_op {
			static boost::mt19937 rng;    // The uniform pseudo-random algorithm
			mutable boost::normal_distribution<Scalar> norm;  // The gaussian combinator

			EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

			template<typename Index>
			inline const Scalar operator() (Index, Index = 0) const {
				return norm(rng); 
			}
		};

		template<typename Scalar> boost::mt19937 scalar_normal_dist_op<Scalar>::rng;

		template<typename Scalar>
		struct functor_traits<scalar_normal_dist_op<Scalar> >{ 
			enum { Cost = 50 * NumTraits<Scalar>::MulCost, 
				PacketAccess = false, IsRepeatable = false 
			}; 
		};
	}
}


//Draw nn samples from a size-dimensional normal distribution
//with a specified mean and covariance
DACE::AlgebraicMatrix<double> mvnrnd(const DACE::AlgebraicVector<double>& mean_aV, 
		const DACE::AlgebraicMatrix<double>& Cov_aM, const unsigned int n_samples){
	
	//Read dimension 
	unsigned int size=mean_aV.size();
	Eigen::VectorXd mean(size);
	for(unsigned int i=0;i<size;i++){
		mean(i)=mean_aV[i];
	}
	Eigen::MatrixXd Cov(size,size);
	for(unsigned int i=0;i<size;i++){
		for(unsigned int j=0;j<size;j++){
			Cov(i,j)=Cov_aM.at(i,j);
		}
	}
					
	//Gaussian functor
	Eigen::internal::scalar_normal_dist_op<double> randN;
	Eigen::internal::scalar_normal_dist_op<double>::rng.seed(1); // Seed the rng


	//Cov decomposition
	Eigen::MatrixXd normTransform(size,size);
	Eigen::LLT<Eigen::MatrixXd> cholSolver(Cov);

	// We can only use the cholesky decomposition if 
	// the covariance matrix is symmetric, pos-definite.
	// But a covariance matrix might be pos-semi-definite.
	// In that case, we'll go to an EigenSolver
	if (cholSolver.info()==Eigen::Success) {
		// Use cholesky solver
		normTransform = cholSolver.matrixL();
	} else {
		// Use eigen solver
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(Cov);
		normTransform = eigenSolver.eigenvectors() *
				eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
	}

	Eigen::MatrixXd samples = (normTransform 
			* Eigen::MatrixXd::NullaryExpr(size,n_samples,randN)).colwise()
					+ mean;
	
	DACE::AlgebraicMatrix<double> samples_aM(size,n_samples);
	for(unsigned int i=0;i<size;i++){
		for(unsigned int j=0;j<n_samples;j++){
			samples_aM.at(i,j)=samples(i,j);
		}
	}
	return samples_aM;
}
*/
template<typename T>
DACE::AlgebraicMatrix<T> Compute_samples(const DACE::AlgebraicVector<T>& points, 
		const unsigned int n_var){

	//Note: there is an extra element at the end of the array, in order to keep track of whether to exit the array.
	unsigned int i[n_var+1]; 
	
	//Max number of iters per for loop
	unsigned int MAX = points.size(); 
	
	//Number of combinations
	unsigned int n_comb=std::pow(MAX,n_var);
	
	//Initialize indexes
	for(unsigned int a=0; a<n_var+1; a++){
		i[a]=0;
	}
	//Vector to store data
	std::vector<std::vector<T>> vector_samples;
	
	//Used to increment all of the indicies correctly, at the end of each loop.
	unsigned int p = 0; 
	while (i[n_var]==0) {

		//Pretend you're inside your nested for loops. The more usual i,j,k,... have been replaced here with i[0], i[1], ..., i[n-1].
		std::vector<T> sample;
		for(unsigned int index=0;index<n_var;index++){
			sample.push_back(points[i[index]]);
		}
		vector_samples.push_back(sample);

		//Now, after you've done your stuff, we need to increment all of the indicies correctly.
		i[0]++;
		
		while(i[p]==MAX) {
			i[p]=0;
			i[++p]++;
			if(i[p]!=MAX){
				p=0;
			}
		}
	}
	DACE::AlgebraicMatrix<T> samples_mat(n_comb,n_var);
	for(unsigned int ind_row=0;ind_row<n_comb;ind_row++){
		for(unsigned int ind_col=0;ind_col<n_var;ind_col++){
			samples_mat.at(ind_row,ind_col)=
				vector_samples[ind_row][ind_col];
		}
	}
	
	return samples_mat;
}

template<typename T>
DACE::AlgebraicVector<double> linspace(const T start_in, const T end_in, 
		const unsigned int n_data_in)
{	
	double start = static_cast<double>(start_in);
	double end = static_cast<double>(end_in);
	unsigned int n_data = static_cast<int>(n_data_in);
	
	DACE::AlgebraicVector<double> vect_lins(n_data);

	if (n_data == 0) { 
		return vect_lins; 
	}
	if (n_data == 1) {
		vect_lins[0]=start;
		return vect_lins;
	}

	double delta = (end - start) / (n_data - 1);
	for(unsigned int i=0; i < n_data-1; ++i){
		vect_lins[i]=start + delta*i;
	}
	vect_lins[n_data-1]=end;

	return vect_lins;
}

//Blockdiagonal matrix
DACE::AlgebraicMatrix<double> blkdiag(const DACE::AlgebraicMatrix<double>& Cov_aM,
		const unsigned int size,const unsigned int count)
{
	
	Eigen::MatrixXd Cov(size,size);
	for(unsigned int i=0;i<size;i++){
		for(unsigned int j=0;j<size;j++){
			Cov(i,j)=Cov_aM.at(i,j);
		}
	}
    Eigen::MatrixXd Cov_blk = Eigen::MatrixXd::Zero(Cov.rows()*count,
    		Cov.cols()*count);
    for(unsigned int i=0;i<count;++i)
    {
        Cov_blk.block(i*Cov.rows(),i*Cov.cols(),
        		Cov.rows(),Cov.cols())=Cov;
    }
    unsigned int size_blk=size*count;
    
    DACE::AlgebraicMatrix<double> Cov_blk_aM(size_blk,size_blk);
    for(unsigned int i=0;i<size_blk;i++){
    	for(unsigned int j=0;j<size_blk;j++){
    		Cov_blk_aM.at(i,j)=Cov_blk(i,j);
    	}
    }
    return Cov_blk_aM;
}

template<typename T>
DACE::AlgebraicMatrix<T> find_Aeci2RSW(const DACE::AlgebraicVector<T>& state_eci){
	
	//Extract position and velocity vectors
	DACE::AlgebraicVector<T> rr_eci(3),vv_eci(3);
	for(unsigned int i=0;i<3;i++){
		rr_eci[i]=state_eci[i];
		vv_eci[i]=state_eci[i+3];
	}
	DACE::AlgebraicVector<T> uu_R=rr_eci/rr_eci.vnorm();
	DACE::AlgebraicVector<T> uu_W=rr_eci.cross(vv_eci)/
			(rr_eci.cross(vv_eci)).vnorm();
	DACE::AlgebraicVector<T> uu_S=uu_W.cross(uu_R);
	
	T eps_X=uu_R[0];
	T eps_Y=uu_R[1];
	T eps_Z=uu_R[2];
	
	T alpha_X=uu_W[0];
	T alpha_Y=uu_W[1];
	T alpha_Z=uu_W[2];
	
	T delta_X=uu_S[0];
	T delta_Y=uu_S[1];
	T delta_Z=uu_S[2];
	
	DACE::AlgebraicMatrix<T> Aeci2RSW(3,3);
	Aeci2RSW.at(0,0)=eps_X;
	Aeci2RSW.at(0,1)=eps_Y;
	Aeci2RSW.at(0,2)=eps_Z;
	
	Aeci2RSW.at(1,0)=delta_X;
	Aeci2RSW.at(1,1)=delta_Y;
	Aeci2RSW.at(1,2)=delta_Z;

	Aeci2RSW.at(2,0)=alpha_X;
	Aeci2RSW.at(2,1)=alpha_Y;
	Aeci2RSW.at(2,2)=alpha_Z;
	
	return Aeci2RSW;
}



/**
 * Multivariate Normal distribution sampling using C++11 and Eigen matrices.
/**
 * Copyright (c) 2014 by Emmanuel Benazera beniz@droidnik.fr, All rights reserved.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library.
 */


/*
  We need a functor that can pretend it's const,
  but to be a good random number generator 
  it needs mutable state.  The standard Eigen function 
  Random() just calls rand(), which changes a global
  variable.
*/
namespace Eigen {
	namespace internal {
		template<typename Scalar>
		struct scalar_normal_dist_op{
			static std::mt19937 rng;                        // The uniform pseudo-random algorithm
			mutable std::normal_distribution<Scalar> norm; // gaussian combinator

			EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

			template<typename Index>
			inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
			inline void seed(const uint64_t &s) { rng.seed(s); }
		};

		template<typename Scalar>
		std::mt19937 scalar_normal_dist_op<Scalar>::rng;

		template<typename Scalar>
		struct functor_traits<scalar_normal_dist_op<Scalar> >
		{ enum { Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false }; };

	} // end namespace internal

	/**
	Find the eigen-decomposition of the covariance matrix
	and then store it for sampling from a multi-variate normal 
	*/
	template<typename Scalar>
	class EigenMultivariateNormal{
		Matrix<Scalar,Dynamic,Dynamic> _covar;
		Matrix<Scalar,Dynamic,Dynamic> _transform;
		Matrix< Scalar, Dynamic, 1> _mean;
		internal::scalar_normal_dist_op<Scalar> randN; // Gaussian functor
		bool _use_cholesky;
		SelfAdjointEigenSolver<Matrix<Scalar,Dynamic,Dynamic> > _eigenSolver; // drawback: this creates a useless eigenSolver when using Cholesky decomposition, but it yields access to eigenvalues and vectors

		public:
		EigenMultivariateNormal(const Matrix<Scalar,Dynamic,1>& mean,const Matrix<Scalar,Dynamic,Dynamic>& covar,
		const bool use_cholesky=false,const uint64_t &seed=std::mt19937::default_seed)
		:_use_cholesky(use_cholesky){
			randN.seed(seed);
			setMean(mean);
			setCovar(covar);
		}

		void setMean(const Matrix<Scalar,Dynamic,1>& mean) { _mean = mean; }
		void setCovar(const Matrix<Scalar,Dynamic,Dynamic>& covar){
			_covar = covar;

			// Assuming that we'll be using this repeatedly,
			// compute the transformation matrix that will
			// be applied to unit-variance independent normals

			if (_use_cholesky){
				Eigen::LLT<Eigen::Matrix<Scalar,Dynamic,Dynamic> > cholSolver(_covar);
				// We can only use the cholesky decomposition if 
				// the covariance matrix is symmetric, pos-definite.
				// But a covariance matrix might be pos-semi-definite.
				// In that case, we'll go to an EigenSolver
				if (cholSolver.info()==Eigen::Success){
					// Use cholesky solver
					_transform = cholSolver.matrixL();
				}
				else{
					throw std::runtime_error("Failed computing the Cholesky decomposition. Use solver instead");
				}
			}
			else{
			_eigenSolver = SelfAdjointEigenSolver<Matrix<Scalar,Dynamic,Dynamic> >(_covar);
			_transform = _eigenSolver.eigenvectors()*_eigenSolver.eigenvalues().cwiseMax(0).cwiseSqrt().asDiagonal();
			}
		}

		/// Draw nn samples from the gaussian and return them
		/// as columns in a Dynamic by nn matrix
		Matrix<Scalar,Dynamic,-1> samples(int nn){
			return (_transform * Matrix<Scalar,Dynamic,-1>::NullaryExpr(_covar.rows(),nn,randN)).colwise() + _mean;
		}
	}; // end class EigenMultivariateNormal
} // end namespace Eigen


//Draw nn samples from a size-dimensional normal distribution
//with a specified mean and covariance
DACE::AlgebraicMatrix<double> mvnrnd(const DACE::AlgebraicVector<double>& mean_aV, 
		const DACE::AlgebraicMatrix<double>& Cov_aM, const unsigned int n_samples){
	
	//Read dimension 
	unsigned int size=mean_aV.size();
	Eigen::VectorXd mean(size);
	for(unsigned int i=0;i<size;i++){
		mean(i)=mean_aV[i];
	}
	Eigen::MatrixXd Cov(size,size);
	for(unsigned int i=0;i<size;i++){
		for(unsigned int j=0;j<size;j++){
			Cov(i,j)=Cov_aM.at(i,j);
		}
	}
	Eigen::EigenMultivariateNormal<double> normX_solver(mean,Cov);				
	Eigen::MatrixXd samples=normX_solver.samples(n_samples);
	
	DACE::AlgebraicMatrix<double> samples_aM(size,n_samples);
	for(unsigned int i=0;i<size;i++){
		for(unsigned int j=0;j<n_samples;j++){
			samples_aM.at(i,j)=samples(i,j);
		}
	}
	return samples_aM;
}


std::tuple<std::vector<double>,std::vector<unsigned int>>
	sort_indexes(std::vector<double> vec){
	
	std::vector<unsigned int> index(vec.size(), 0);
	for (unsigned int i = 0 ; i<index.size() ; i++) {
		index.at(i) = i;
	}
	std::sort(index.begin(), index.end(),
		[&](const int& a, const int& b) {
			return (vec.at(a) < vec.at(b));
		}
	);
	
	std::vector<double> vec_sorted(vec.size());
	for (int i = 0 ; i != index.size() ; i++) {
		vec_sorted.at(i)=vec.at(index.at(i));
	}
	
	return std::make_tuple(vec_sorted,index);
	
}


//Compute eigenvalues and eigenvectors
std::tuple<std::vector<double>,DACE::AlgebraicMatrix<double>>
	Eigenvv(const DACE::AlgebraicMatrix<double>& M){
		
	//Copy the matrix in an eigen-like matrix
	Eigen::MatrixXd M_eig(M.nrows(),M.ncols());
	for(unsigned int i=0;i<M.nrows();i++){
		for(unsigned int j=0;j<M.ncols();j++){
			M_eig(i,j)=M.at(i,j);
		}
	}
	
	//Find eigenvalues and vectors
	Eigen::EigenSolver<Eigen::MatrixXd> es(M_eig);
	
	//Store the eigenvalues
	std::vector<double> lambda_vec(M.nrows());
	for(unsigned int i=0;i<M.nrows();i++){
		lambda_vec.at(i)=real(es.eigenvalues()[i]);
	}
	
	//Store the eigenvector matrix
	DACE::AlgebraicMatrix<double> V(M.nrows(),M.ncols());
	for(unsigned int i=0;i<M.nrows();i++){
		for(unsigned int j=0;j<M.ncols();j++){
			V.at(i,j)=real(es.eigenvectors()(i,j));;
		}
	}
	
	return std::make_tuple(lambda_vec,V);
		
}

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}
#endif