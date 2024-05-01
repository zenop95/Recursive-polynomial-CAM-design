#ifndef Propagation_H_
#define Propagation_H_

#include<dace/dace.h>
#include<cspice/SpiceUsr.h>

//Propagators
#include <astro/AstroLibrary.h>

//AIDA
#include <dynorb/AIDAwrappers.h>

// SADA
#include <sada/Transformations.h>
#include <sada/HEOSATDynamics.h> 
#include <sada/mean2oscNonSingular.h>

//LS library
#include "Objects.hpp"

//Time
#include <time.h>
#include <chrono>

template<typename T>
std::vector<DACE::AlgebraicVector<T>> Propagate_state(const double et_start,
		std::vector<double> et_end_vec, 
		const DACE::AlgebraicVector<T>& state_start,
		const Dynamics_param<T>& param){
	
	//Sort the vector of time instants
	std::sort(et_end_vec.begin(),et_end_vec.end());
	
	//Read propagation method
	std::string method=param.method;

	//Read initial conditions and parameters and propagate
	DACE::AlgebraicVector<T> state_end(6);
	
	const unsigned int n_inst=et_end_vec.size();
	std::vector<DACE::AlgebraicVector<T>> state_end_mat;
	
	if(method=="AIDA"){
		
		//Extract parameters
		T Bfactor = param.Bfactor;
		T SRPC  = param.SRPC;
		
		unsigned int gravOrd = param.gravOrd;
		std::string gravmodel = param.gravmodel;
		DACE::AlgebraicVector<int> Aflags=param.AIDA_flags;
		int AIDA_flags[3];
		for(unsigned int i=0;i<3;i++){
			AIDA_flags[i]=Aflags[i];
		}	
		
		//Initialize AIDA
		AIDADynamics<T> aidaCartDyn(gravmodel, gravOrd, AIDA_flags, 
				Bfactor, SRPC);
		
		//Run cycle
		double et_start_loop=et_start;
		DACE::AlgebraicVector<T> state_start_loop=state_start;
		for(unsigned int i=0;i<n_inst;i++){
			
			//Compute state vector at et_end
			double et_end=et_end_vec[i];
			const clock_t begin_time = clock();

			state_end=RK78(6,state_start_loop,et_start_loop,
				et_end,aidaCartDyn,false,false,param.tol);
				//	std::cout << " It takes " <<  (float)((clock() - begin_time))/CLOCKS_PER_SEC << " second(s) for propagation in propagatestate " <<i <<"DT = " << et_end-et_start_loop <<std::endl;

			//Store results in matrix
			state_end_mat.push_back(state_end);
			
			//Update initial conditions
			et_start_loop=et_end;
			state_start_loop=state_end;
		}
		
	}
	
	else if(method=="SADA"){
		
		//Extract parameters
		const unsigned int SunMoonFlag        = param.SunMoonFlag;
		const unsigned int SRPFlag            = param.SRPFlag; 
		const unsigned int dragFlag           = param.dragFlag;
		const unsigned int tesseralFlag       = param.tesseralFlag;
		const unsigned int shortPeriodicsFlag = param.shortPeriodicsFlag;
		
		T Bfactor = param.Bfactor;
		T SRPC  = param.SRPC;

		//Find jdate_start
		SpiceChar jd_char[TIMLEN_SADA];
		timout_c(et_start, TIMFMT_SADA,TIMLEN_SADA,jd_char);
		double jdate_start=atof(jd_char);

		//Find deltat vector
		const unsigned int n_inst=et_end_vec.size();
		DACE::AlgebraicVector<double> deltat_vec(n_inst);
		for(unsigned int i=0;i<n_inst;i++){
			double deltat_single=(et_end_vec.at(i)-et_start)/spd_c();
			deltat_vec[i]=deltat_single;
		}

		// Create instance of HEOSATDynamics class
		HEOSATDynamics<T> HEOSATDyn(SunMoonFlag, SRPFlag, dragFlag, 
			tesseralFlag);

		// Initialize SADA propagator
		const double secondsPerDay = spd_c();
		SpiceDouble mu_Earth;
		SpiceInt dim;
		bodvrd_c("earth","GM", 1 , &dim, &mu_Earth);
		double t0, tf;

		//Convert state to Keplerian parameters
		DACE::AlgebraicVector<T> KeplPar_start=cart2kepM(state_start,mu_Earth);

		//Scaling factor
		const double Lsc = DACE::cons(KeplPar_start[0]);
		double m_ut;
		
		const bool flag_linear = param.flag_linear;
		if(flag_linear==true && n_inst>1){

			//Identify the different sets of times
			unsigned int n_sets=1;
			for(unsigned int i=1;i<deltat_vec.size();i++){
				double gap=(deltat_vec[i]-deltat_vec[i-1])*spd_c();
				if(gap>10.*60.){
					n_sets=n_sets+1;
				}
			}
			//std::cout<<"Number of identified sets: "<<n_sets<<std::endl;
			
			//Divide instants in sets
			std::vector<std::vector<double>> dt_sets(n_sets);
			unsigned int index=0;
			dt_sets.at(0).push_back(deltat_vec[0]);
			for(unsigned int i=1;i<deltat_vec.size();i++){
				double gap=(deltat_vec[i]-deltat_vec[i-1])*spd_c();
				if(gap>10.*60.){
					index=index+1;
				}
				dt_sets.at(index).push_back(deltat_vec[i]);
			}
			
			//Consider each set separately
			const double mu = 1.0; 
			DACE::AlgebraicVector<T> delaunay_start(6);
			double jdate_0=jdate_start;
			
			for(unsigned int i=0;i<dt_sets.size();i++){				
				
				//Middle epoch
				std::vector<double> dt_sets_single=dt_sets.at(i);
				double dt_start = dt_sets_single.at(0);
				double dt_end   = dt_sets_single.at(dt_sets_single.size()-1);
				double dt_mid   = (dt_start+dt_end)/2.;
				
				if(i==0){
					
					//Initialize
					m_ut=HEOSATDyn.initheosat(Lsc, jdate_start, 0.0, dt_mid* 
						secondsPerDay, &t0, &tf, Bfactor, SRPC);
					
					// Scaled variables
					KeplPar_start[0] = KeplPar_start[0] / Lsc; 
							
					// Transform state from J2000 to True-Of-Date reference frame
					KeplPar_start = kepJ2000ToTOD(KeplPar_start, jdate_start, mu);
					
					// Convert state from osculating to mean elements	
					KeplPar_start = osc2meanJ2sunMoonNonSingular
						(KeplPar_start, jdate_start, mu, 
						shortPeriodicsFlag, Lsc, m_ut);

					// Initial Delaunay elements
					delaunay_start =
						kep2delaunay(KeplPar_start, mu);
				}
				else{
					//Compute new t0 and tf
					t0=(jdate_start-jdate_0)*spd_c()/m_ut;
					tf=(jdate_start+dt_mid-jdate_0)*spd_c()/m_ut;
				}

				// RK78 integration settings
				bool flag_SADA = true;
				bool returnIntermediatePoints = false;

				// Propagate orbit
				DACE::AlgebraicVector<T> delaunay_end = 
					RK78(6, delaunay_start, t0, tf, HEOSATDyn, 
						flag_SADA,returnIntermediatePoints,param.tol);
						
				DACE::AlgebraicVector<T> KeplPar_end = 
						delaunay2kep(delaunay_end, mu);

				// Convert state from mean to osculating elements
				KeplPar_end = mean2oscJ2sunMoonNonSingular
					(KeplPar_end, jdate_start+dt_mid, mu, shortPeriodicsFlag, 
					Lsc, m_ut);
				
				// Transform state from True-Of-Date to J2000 reference frame
				KeplPar_end = kepTODToJ2000(KeplPar_end,
						jdate_start+dt_mid, mu);
				
				//Rescale the semi-major axis
				KeplPar_end[0]=KeplPar_end[0]* Lsc;
				
				//Convert to cartesian coordinates
				state_end=kepM2cart(KeplPar_end,mu_Earth);
				
				//Compute time derivative at dt_mid
				DACE::AlgebraicVector<T> f_end=HEOSATDyn.evaluate(delaunay_end, tf);
				
				//Loop on all time instants
				for(unsigned int j=0;j<dt_sets_single.size();j++){
					
					double dt_single=dt_sets_single.at(j)-dt_mid;
					
					if(std::abs(dt_single)<1e-10){
						state_end_mat.push_back(state_end);
					}
					else{
						//Linear approximation
						DACE::AlgebraicVector<T> delaunay_single=delaunay_end+
							f_end*dt_single*spd_c()/m_ut;
																		
						//--------------- Save intermediate state vectors ------------//
						// Keplerian orbital elements
						DACE::AlgebraicVector<T> KeplPar_single = 
							delaunay2kep(delaunay_single, mu);
						
						// Julian date
						const double jdate_single = jdate_start + 
							dt_sets_single.at(j);
						
						//Compute the new t0 and tf
						t0=(jdate_start+dt_mid-jdate_0)*spd_c()/m_ut;
						tf=(jdate_single-jdate_0)*spd_c()/m_ut;

						// Convert state from mean to osculating elements
						KeplPar_single = mean2oscJ2sunMoonNonSingular
							(KeplPar_single, jdate_single, mu, shortPeriodicsFlag, 
							Lsc, m_ut);
						
						// Transform state from True-Of-Date to J2000 reference frame
						KeplPar_single = kepTODToJ2000(KeplPar_single,
								jdate_single, mu);
						
						//Rescale the semi-major axis
						KeplPar_single[0]=KeplPar_single[0]* Lsc;
						
						//Convert to cartesian coordinates
						DACE::AlgebraicVector<T> state_single=
							kepM2cart(KeplPar_single,mu_Earth);
						
						//Store results in matrix
						state_end_mat.push_back(state_single);
					}
				}

				//Update initial conditions
				jdate_start = jdate_start + dt_mid;
				
				//Remove dt from all dt_sets entries
				for(unsigned int j=0;j<dt_sets.size();j++){
					for(unsigned int k=0;k<dt_sets.at(j).size();k++){
						dt_sets.at(j).at(k)=dt_sets.at(j).at(k)-dt_mid;
					}
				}
				
				//Update initial state vector
				delaunay_start=delaunay_end;
				
			}
		}
		else{
			//Initialize
			const double jdate_0=jdate_start;
			double deltat=deltat_vec[0];
			
			m_ut=HEOSATDyn.initheosat(Lsc, jdate_start, 0.0, deltat* 
					secondsPerDay, &t0, &tf, Bfactor, SRPC);
			
			//Check for repeated instants
			unsigned int index=0;
			if(tf==t0){
				state_end_mat.push_back(state_start);
				index=index+1;
				jdate_start = jdate_start + deltat;
				
				while(tf==t0 && index<deltat_vec.size()){
						
					//Retrieve t0 and tf
					deltat=deltat_vec[index];
					
					t0=(jdate_start-jdate_0)*spd_c()/m_ut;
					tf=(jdate_start+deltat-jdate_0)*spd_c()/m_ut;

					if(tf==t0){
						state_end_mat.push_back(state_start);
						index=index+1;
						jdate_start = jdate_start + deltat;
					}
				}
			}
			
			// Scaled variables
			const double mu = 1.0; 
			KeplPar_start[0] = KeplPar_start[0] / Lsc; 
					
			// Transform state from J2000 to True-Of-Date reference frame
			KeplPar_start = kepJ2000ToTOD(KeplPar_start, jdate_start, mu);
			
			// Convert state from osculating to mean elements	
			KeplPar_start = osc2meanJ2sunMoonNonSingular(KeplPar_start, 
				jdate_start, mu, shortPeriodicsFlag, Lsc, m_ut);

			// Initial Delaunay elements
			DACE::AlgebraicVector<T> delaunay_start =kep2delaunay(KeplPar_start, mu);

			// RK78 integration settings
			bool flag_SADA = true;
			bool returnIntermediatePoints = false;
				
			//Iterate on all time instants
			for(unsigned int i=index;i<n_inst;i++){
				
				// Propagate orbit
				DACE::AlgebraicVector<T> delaunay_end = 
					RK78(6, delaunay_start, t0, tf, HEOSATDyn, 
					flag_SADA,returnIntermediatePoints,param.tol);

				//--------------- Save intermediate state vectors ------------//
				// Keplerian orbital elements
				DACE::AlgebraicVector<T> KeplPar_end = delaunay2kep(delaunay_end, mu);

				// Julian date
				const double jdate_end = jdate_start + deltat;

				// Convert state from mean to osculating elements
				KeplPar_end = mean2oscJ2sunMoonNonSingular(KeplPar_end, jdate_end, mu, 
						shortPeriodicsFlag, Lsc, m_ut);
				
				// Transform state from True-Of-Date to J2000 reference frame
				KeplPar_end = kepTODToJ2000(KeplPar_end,
						jdate_end, mu);
				
				//Rescale the semi-major axis
				KeplPar_end[0]=KeplPar_end[0]* Lsc;
				
				//Convert to cartesian coordinates
				state_end=kepM2cart(KeplPar_end,mu_Earth);
				
				//Store results in matrix
				state_end_mat.push_back(state_end);
				
				//-----------------------------------------------------------//
				
				//Update initial conditions
				jdate_start = jdate_start + deltat;
				if(i<n_inst-1){
					deltat=deltat_vec[i+1]-deltat_vec[i];
				}
				t0=(jdate_start-jdate_0)*spd_c()/m_ut;
				tf=(jdate_start+deltat-jdate_0)*spd_c()/m_ut;
				
				//Check for repeated instants
				if(tf==t0){
					state_end_mat.push_back(state_start);
					i=i+1;
					jdate_start = jdate_start + deltat;
					
					while(tf==t0 && i<deltat_vec.size()-1){
							
						//Retrieve t0 and tf
						deltat=deltat_vec[i+1]-deltat_vec[i];
						t0=(jdate_start-jdate_0)*spd_c()/m_ut;
						tf=(jdate_start+deltat-jdate_0)*spd_c()/m_ut;
						
						if(tf==t0){
							state_end_mat.push_back(state_start);
							i=i+1;
							jdate_start = jdate_start + deltat;
						}
					}
				}
				delaunay_start=delaunay_end;
			}
		}

		
	}
	return state_end_mat;
}


template<typename T>
std::vector<DACE::AlgebraicVector<T>> Propagate_state_linear(const double et_start,
		std::vector<double> et_end_vec, 
		const DACE::AlgebraicVector<T>& state_start,
		const Dynamics_param<T>& param){
	
	//Sort the vector of time instants
	std::sort(et_end_vec.begin(),et_end_vec.end());
	
	//Read propagation method
	std::string method=param.method;

	//Read initial conditions and parameters and propagate
	DACE::AlgebraicVector<T> state_end(6);
	
	const unsigned int n_inst=et_end_vec.size();
	std::vector<DACE::AlgebraicVector<T>> state_end_mat;
	
	if(method=="AIDA"){
		
		//Extract parameters
		T Bfactor = param.Bfactor;
		T SRPC  = param.SRPC;
		
		unsigned int gravOrd = param.gravOrd;
		std::string gravmodel = param.gravmodel;
		DACE::AlgebraicVector<int> Aflags=param.AIDA_flags;
		int AIDA_flags[3];
		for(unsigned int i=0;i<3;i++){
			AIDA_flags[i]=Aflags[i];
		}	
		
		//Initialize AIDA
		AIDADynamics<T> aidaCartDyn(gravmodel, gravOrd, AIDA_flags, 
				Bfactor, SRPC);
		
		//Time derivative
		DACE::AlgebraicVector<T> dstate=
			aidaCartDyn.evaluate(state_start,et_start);
			
		//Run cycle
		for(unsigned int i=0;i<n_inst;i++){
			
			//Compute state vector at et_end
			double et_end=et_end_vec[i];
			
			state_end=state_start+dstate*(et_end-et_start);

			//Store results in matrix
			state_end_mat.push_back(state_end);
		}
		
	}
	
	else if(method=="SADA"){
		
		//Extract parameters
		const unsigned int SunMoonFlag        = param.SunMoonFlag;
		const unsigned int SRPFlag            = param.SRPFlag; 
		const unsigned int dragFlag           = param.dragFlag;
		const unsigned int tesseralFlag       = param.tesseralFlag;
		const unsigned int shortPeriodicsFlag = param.shortPeriodicsFlag;
		
		T Bfactor = param.Bfactor;
		T SRPC  = param.SRPC;

		//Find jdate_start
		SpiceChar jd_char[TIMLEN_SADA];
		timout_c(et_start, TIMFMT_SADA,TIMLEN_SADA,jd_char);
		double jdate_start=atof(jd_char);

		//Find deltat vector
		const unsigned int n_inst=et_end_vec.size();
		DACE::AlgebraicVector<double> deltat_vec(n_inst);
		for(unsigned int i=0;i<n_inst;i++){
			double deltat_single=(et_end_vec.at(i)-et_start)/spd_c();
			deltat_vec[i]=deltat_single;
		}

		// Create instance of HEOSATDynamics class
		HEOSATDynamics<T> HEOSATDyn(SunMoonFlag, SRPFlag, dragFlag, 
				tesseralFlag);

		// Initialize SADA propagator
		const double secondsPerDay = spd_c();
		SpiceDouble mu_Earth;
		SpiceInt dim;
		bodvrd_c("earth","GM", 1 , &dim, &mu_Earth);
		double t0, tf;

		//Convert state to Keplerian parameters
		DACE::AlgebraicVector<T> KeplPar_start=cart2kepM(state_start,mu_Earth);

		//Scaling factor
		const double Lsc = DACE::cons(KeplPar_start[0]);		

		//Consider each set separately
		const double mu = 1.0; 
		
		//Initialize
		double m_ut=HEOSATDyn.initheosat(Lsc, jdate_start, 0.0, deltat_vec.at(0)* 
			secondsPerDay, &t0, &tf, Bfactor, SRPC);
		
		// Scaled variables
		KeplPar_start[0] = KeplPar_start[0]/Lsc; 
				
		// Transform state from J2000 to True-Of-Date reference frame
		KeplPar_start = kepJ2000ToTOD(KeplPar_start, jdate_start, mu);
		
		// Convert state from osculating to mean elements	
		KeplPar_start = osc2meanJ2sunMoonNonSingular
			(KeplPar_start, jdate_start, mu, 
			shortPeriodicsFlag, Lsc, m_ut);

		// Initial Delaunay elements
		DACE::AlgebraicVector<T> delaunay_start =
			kep2delaunay(KeplPar_start, mu);
			
		//Compute time derivative at jdate_start
		DACE::AlgebraicVector<T> dstate_start=
			HEOSATDyn.evaluate(delaunay_start, t0);
		
		//Loop on all time instants	
		for(unsigned int i=0;i<deltat_vec.size();i++){				
			
			//Current epoch
			double dt_single = deltat_vec.at(i);
	
			if(std::abs(dt_single)<1e-10){
				state_end_mat.push_back(state_start);
			}
			else{

				//Linear approximation
				DACE::AlgebraicVector<T> delaunay_single=delaunay_start+
					dstate_start*dt_single*spd_c()/m_ut;
																
				//--------------- Save intermediate state vectors ------------//
				// Keplerian orbital elements
				DACE::AlgebraicVector<T> KeplPar_single = 
					delaunay2kep(delaunay_single, mu);
				
				// Julian date
				const double jdate_single = jdate_start + 
					dt_single;

				// Convert state from mean to osculating elements
				KeplPar_single = mean2oscJ2sunMoonNonSingular
					(KeplPar_single, jdate_single, mu, shortPeriodicsFlag, 
					Lsc, m_ut);
				
				// Transform state from True-Of-Date to J2000 reference frame
				KeplPar_single = kepTODToJ2000(KeplPar_single,
					jdate_single, mu);
				
				//Rescale the semi-major axis
				KeplPar_single[0]=KeplPar_single[0]* Lsc;
				
				//Convert to cartesian coordinates
				DACE::AlgebraicVector<T> state_single=
					kepM2cart(KeplPar_single,mu_Earth);
				
				//Store results in matrix
				state_end_mat.push_back(state_single);
			}
		}
	}
	
	return state_end_mat;
}

#endif