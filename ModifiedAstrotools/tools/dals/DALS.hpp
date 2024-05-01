#ifndef DALS_ADAPTIVE_H_
#define DALS_ADAPTIVE_H_

#include<tuple>
#include<iostream>
#include<random> 
#include <time.h>
#include <chrono>

#include<dace/dace.h>
#include<cspice/SpiceUsr.h>

#include "Standards.hpp"
#include "Objects.hpp"
#include "LSRoutines.hpp"
#include "DlibClass.hpp"
#include "Propagation.hpp"
#include "Exceptions.hpp"

namespace lssolver{

	
class DALS {
	
	//private:
		//*******************************************************************************
		//************************* Private member variables ****************************
		//*******************************************************************************
		//Observer parameters
		std::vector<Observer_param> sd_param_obs;
		
		//Earth parameters
		Earth_param sd_param_Earth;
		
		//State estimates parameters
		Estimate_param sd_param_OD;
		
		//DALS input parameters
		DALS_param sd_param_DALS;
		
		//Dynamics parameters
		Dynamics_param<double> sd_param_dyn;
		
		//Measurements parameters
		Measurements_param<double> sd_param_meas;

		//*******************************************************************************
		//************************ Private member functions *****************************
		//*******************************************************************************
		//Generate observatory kernels
		void generate_obs_kernels();
		
		//Read measurement file
		void read_measurements(const std::vector<std::string>& filename_meas);
		
		//Read state estimates file
		std::tuple<double, DACE::AlgebraicVector<double>, 
			DACE::AlgebraicMatrix<double>,double,double> 
			read_state_data(const std::string& filename_state) const;
		
		//Propagate estimates from t_data to t_est
		void propagate_apr_data(const double et_data, 
			const DACE::AlgebraicVector<double>& state_apr_t_data,
			const DACE::AlgebraicMatrix<double>& Cov_apr_t_data,
			const double et_est, const double Bfactor, const double SRPC);
		
		//Initialize DA state
		std::tuple<DACE::AlgebraicVector<DACE::DA>,DACE::DA,DACE::DA> 
			initialize_DA(const double Bfactor,const double SRPC) const;
		
		//Run LS iterative cycle
		std::tuple<DACE::AlgebraicVector<double>,double,double,DACE::DA,
		std::vector<Observations<DACE::DA>>,double> 
			LSsolver(double et_est,DACE::AlgebraicVector<DACE::DA>& state_est,
			DACE::DA Bfactor_est,DACE::DA SRPC_est) const;
		
		//Compute the covariance
		std::tuple<DACE::AlgebraicMatrix<double>,DACE::AlgebraicMatrix<double>,
			bool> compute_covariance( DACE::DA& FinalRes) const;
			
		//Write the output
		void write_results(const double et_est, 
			const DACE::AlgebraicVector<double>& state_est, 
			const DACE::AlgebraicMatrix<double>& Cov_est, const double Bfactor_est, 
			const double SRPC_est, const bool flag_dec,
			const DACE::AlgebraicMatrix<double>& Cov_param, 
			const double flag_conv, const DACE::DA& FinalRes,
			const std::vector<Observations<DACE::DA>>& Obs_est) const;
			
		
	public:
		//*******************************************************************************
		//************************* Public member functions *****************************
		//*******************************************************************************
		//Contructor
		DALS(const std::string& json_input);
		
		//Run whole LS process
		std::tuple<double,DACE::AlgebraicVector<double>,DACE::AlgebraicMatrix<double>,
			double, double, DACE::AlgebraicMatrix<double>, DACE::DA, bool> run();
		
};


//********************************************************************************
//************************* Private member functions *****************************
//********************************************************************************
void DALS::generate_obs_kernels(){
		
	std::cout<<"    2) Generating obs kernels:";
	
	//Count the number of sensors
	const unsigned int n_sens=this->sd_param_obs.size();
	
	//Avoid considering the same kernel multiple times
	std::vector<std::string> stations_all;
	
	//For loop over all the sensors
	unsigned int index_obs=0;
	for(unsigned int j=0;j<n_sens;j++){
		
		//Create a kernel for all the observers of sensor "i"
		for(unsigned int i=0;i<this->sd_param_obs.at(j).obs_name.size();i++){
			
			//Check whether the station has already been considered
			bool flag_rep=false;
			for(unsigned int k=0;k<stations_all.size();k++){
				if(stations_all.at(k)==
					this->sd_param_obs.at(j).obs_name.at(i)){
					flag_rep=true;
					break;
				}
			}
			
			//Observer filename
			std::string filename_obs = this->sd_param_obs.at(j).obs_path.at(i)+
				this->sd_param_obs.at(j).obs_name.at(i);
			
			//Create kernels only if they do not exist
			if(flag_rep==false){
				
				//Add name to stations_all
				stations_all.push_back
					(this->sd_param_obs.at(j).obs_name.at(i));
					
				//Add one to index
				index_obs=index_obs+1;
				
				//Write .sut file
				std::string filename_sut;
				double lon_model,lat_model,alt_model;
				std::tie(filename_sut,lon_model,lat_model,
					alt_model)=Generate_SUT(this->sd_param_obs.at(j),i,index_obs,
					this->sd_param_Earth.Earth_model,this->sd_param_DALS.et_est);
				
				//Set filenames
				std::string TPC_path=this->sd_param_Earth.kernels.path+
					this->sd_param_Earth.kernels.pck;
				
				//Remove kernels if existing
				std::string bsp_file=filename_obs+".bsp";
				remove(bsp_file.c_str());
				std::string tf_file=filename_obs+".tf";
				remove(tf_file.c_str());
				
				//Launch pinpoint
				std::string launchPP = "./pinpoint -def "+filename_sut+
					" -pck "+TPC_path+" -spk "+bsp_file+" > nul";
				std::system(launchPP.c_str());
				
				//Delete .sut file
				remove(filename_sut.c_str());
				
				//Write tk file
				Generate_tk(this->sd_param_obs.at(j),this->sd_param_Earth.Earth_model,
					lon_model,lat_model,alt_model,i,index_obs);
			}
			
			//Save names
			this->sd_param_obs.at(j).filename_obs.push_back(filename_obs);
			this->sd_param_obs.at(j).obsRF_name.push_back
				((this->sd_param_obs.at(j).obs_name.at(i)+"_TOPO"));
		
		}
	}
	std::cout<<"\033[32m done\033[0m"<<std::endl;
	
}

void DALS::read_measurements(const std::vector<std::string>& filename_meas){
	
	std::cout<<"    3) Reading measurement file:";
	
	//Initialize list of measurements
	std::vector<Observations<double>> Observation_all;
	
	//Iterate on all TDM files
	for(unsigned int j=0;j<filename_meas.size();j++){
		
		//Read TDM file and store all data lines
		std::ifstream file(filename_meas.at(j));
		std::vector <std::string> Obs;
		std::string temp_line;
		
		std::string InitialLine = "DATA_START";
		std::string EndLine = "DATA_STOP";
		
		if (file.is_open()) {
			
			//Get to the data block
			while (temp_line != InitialLine){
				std::getline(file,temp_line,'\n');
			}
			
			//Read data
			std::getline(file,temp_line);
			while (temp_line != EndLine){
				//Store it
				Obs.push_back(temp_line);
				
				//Read single line
				std::getline(file,temp_line);
				
			}
		}
		file.close();
		
		//Vectors of measurements
		std::vector<double> ra_vect; 
		std::vector<double> decl_vect;
		
		std::vector<double> az_vect;
		std::vector<double> el_vect;
		std::vector<double> doppler_vect;
		std::vector<double> range_vect;	
		
		//Vector of time instants
		std::vector<double> et_vect; 
		
		//*Store data
		unsigned int N_inst = Obs.size(); 
		std::string meas;
		for(unsigned int i=0; i<N_inst; i++) {
			
			//Read line and separate elements
			std::string A, B, C, D;
			std::stringstream ss(Obs.at(i));
			ss >> A >> B >> C >> D ;

			if(i==0){
				//Check which is the first available measurement;
				meas=A;
			}
			
			Sensor_setup sensor_setup= 
				this->sd_param_obs.at(j).sensor_setup;
			if(this->sd_param_obs.at(j).type=="optical"){
				if(A == sensor_setup.ra.TDM_name){
					ra_vect.push_back(std::stod(D));
					if(A==meas){
						SpiceDouble et_temp;
						utc2et_c(C.c_str(), &et_temp);
						et_vect.push_back(et_temp);
					}
						
				}
				if(A == sensor_setup.decl.TDM_name){
					decl_vect.push_back(std::stod(D));	
					if(A==meas){
						SpiceDouble et_temp;
						utc2et_c(C.c_str(), &et_temp);
						et_vect.push_back(et_temp);	
					}
				}
			}
			else if(this->sd_param_obs.at(j).type=="radar"){
				if(A == sensor_setup.az.TDM_name){
					az_vect.push_back(std::stod(D));
					if(A==meas){
						SpiceDouble et_temp;
						utc2et_c(C.c_str(), &et_temp);
						et_vect.push_back(et_temp);
					}	
				}
				if(A == sensor_setup.el.TDM_name){
					el_vect.push_back(std::stod(D));	
					if(A==meas){
						SpiceDouble et_temp;
						utc2et_c(C.c_str(), &et_temp);
						et_vect.push_back(et_temp);	
					}
				}
				
				if(A == sensor_setup.doppler.TDM_name){
					doppler_vect.push_back(std::stod(D));	
					if(A==meas){
						SpiceDouble et_temp;
						utc2et_c(C.c_str(), &et_temp);
						et_vect.push_back(et_temp);	
					}
				}
				
				if(A == sensor_setup.range.TDM_name){
					range_vect.push_back(std::stod(D));	
					if(A==meas){
						SpiceDouble et_temp;
						utc2et_c(C.c_str(), &et_temp);
						et_vect.push_back(et_temp);	
					}
				}	
			}
		}
		
		//Save data
		Observations<double> Observation;
		Observation.et_vect=et_vect;
		
		Observation.ra_vect=ra_vect;
		Observation.decl_vect=decl_vect;
		
		Observation.az_vect=az_vect;
		Observation.el_vect=el_vect;
		Observation.doppler_vect=doppler_vect;
		Observation.range_vect=range_vect;
		
		//Save how many measurements per time instant
		unsigned int n_meas_inst=(ra_vect.size()+decl_vect.size()+
			az_vect.size()+el_vect.size()+doppler_vect.size()+
			range_vect.size())/et_vect.size();
		Observation.n_meas_inst=n_meas_inst;
		
		//Store data
		Observation_all.push_back(Observation);

	}
	this->sd_param_meas.Obs=Observation_all;
	
	std::cout<<"\033[32m done\033[0m"<<std::endl;
}		

			
std::tuple<double, DACE::AlgebraicVector<double>, DACE::AlgebraicMatrix<double>,
	double, double> DALS::read_state_data(const std::string& filename_state) const{

	std::cout<<"    4) Reading state estimates:";
	
	//Check state estimates file
	state_except(this->sd_param_OD);
	
	//*Store data
	std::ifstream fileOD(filename_state);
	double et_data;
	std::string date_data;
	DACE::AlgebraicVector<double> state_data(6);
	DACE::AlgebraicMatrix<double> Cov_data(6,6);
	double Bfactor=1e-3;
	double SRPC=1e-3;

	//*Read epoch
	std::string string_tmp;
	std::getline(fileOD, string_tmp);

	std::string time_str;
	date_data=string_tmp;
	ConstSpiceChar *char_tmp=string_tmp.c_str();
	str2et_c(char_tmp,&et_data);
					
	//Read state vector
	std::getline(fileOD, string_tmp);
	std::stringstream ss(string_tmp);
	std::string rr_x,rr_y,rr_z,vv_x,vv_y,vv_z;
	ss>> rr_x >> rr_y >> rr_z >> vv_x >> vv_y >> vv_z;
	state_data[0]=std::stod(rr_x);
	state_data[1]=std::stod(rr_y);
	state_data[2]=std::stod(rr_z);
	state_data[3]=std::stod(vv_x);
	state_data[4]=std::stod(vv_y);
	state_data[5]=std::stod(vv_z);

	//Check for Cov_apr availability
	const bool flag_Cov_av=this->sd_param_OD.flag_Cov_av;
	
	if(flag_Cov_av==true){		
		for(unsigned int i=0;i<6;i++){
			//Extract line by line the Covariance matrix
			std::getline(fileOD, string_tmp);
			std::stringstream ss(string_tmp);
			std::string sigma1,sigma2,sigma3,sigma4,sigma5,sigma6;
			ss>> sigma1>>sigma2>>sigma3>>sigma4>>sigma5>>sigma6;
			Cov_data.at(i,0)=std::stod(sigma1);
			Cov_data.at(i,1)=std::stod(sigma2);
			Cov_data.at(i,2)=std::stod(sigma3);
			Cov_data.at(i,3)=std::stod(sigma4);
			Cov_data.at(i,4)=std::stod(sigma5);
			Cov_data.at(i,5)=std::stod(sigma6);
		}
	}
	
	//Check for Bfactor availability
	const bool flag_Bfactor_av=this->sd_param_OD.flag_Bfactor_av;
	if(flag_Bfactor_av==true){
		std::getline(fileOD, string_tmp);
		Bfactor=std::stod(string_tmp);
	}
	
	//Check for SRPC availability
	const bool flag_SRPC_av=this->sd_param_OD.flag_SRPC_av;
	if(flag_SRPC_av==true){
		std::getline(fileOD, string_tmp);
		SRPC=std::stod(string_tmp);
	}
	fileOD.close();
	
	std::cout<<"\033[32m done\033[0m"<<std::endl;
	return  std::make_tuple(et_data,state_data,Cov_data,Bfactor,SRPC);
}		


void DALS::propagate_apr_data(const double et_data, 
	const DACE::AlgebraicVector<double>& state_apr_t_data,
	const DACE::AlgebraicMatrix<double>& Cov_apr_t_data,
	const double et_est, const double Bfactor, const double SRPC){
	
	std::cout<<"    5) Propagating state estimate to et_est:";
	DACE::AlgebraicVector<double> state_apr_t_est(6);
	DACE::AlgebraicMatrix<double> Cov_apr_t_est(6,6);
	
	if(et_est!=et_data){
			
		//Inizialize an auxiliary-DA variable for STM computation
		DACE::AlgebraicVector<DACE::DA> state_aux(6);
		for(unsigned int i=0;i<6;i++){
			state_aux[i]=state_apr_t_data[i]+DACE::DA(i+1);
		}
		
		//Create param for dynamics		
		Dynamics_param<DACE::DA> param;
		param.build(this->sd_param_dyn,Bfactor,SRPC);
		
		//Propagate to et_est
		DACE::AlgebraicVector<DACE::DA> state_aux_prop = 
				Propagate_state(et_data,{et_est},
				state_aux,param).at(0);
		
		//A priori state estimate
		state_apr_t_est=DACE::cons(state_aux_prop);
		
		if(this->sd_param_OD.flag_Cov_av==true){
			
			//Compute STM
			DACE::AlgebraicMatrix<DACE::DA> STM_data2est_DA(6,6);
			for(unsigned int i=0;i<6;i++){
				for(unsigned int j=0;j<6;j++){
					STM_data2est_DA.at(i,j)=state_aux_prop[i].deriv(j+1);
				}
			}
			DACE::AlgebraicMatrix<double> STM_data2est=DACE::cons(STM_data2est_DA);
			
			//A priori covariance
			Cov_apr_t_est=STM_data2est*Cov_apr_t_data*
					DACE::transpose(STM_data2est);
		}
	}
	else{
		state_apr_t_est = state_apr_t_data;
		if(this->sd_param_OD.flag_Cov_av==true){
			Cov_apr_t_est   = Cov_apr_t_data;
		}
	}
	
	this->sd_param_OD.state_apr=state_apr_t_est;
	if(this->sd_param_OD.flag_Cov_av==true){
		this->sd_param_OD.Cov_apr=Cov_apr_t_est;
	}
	
	std::cout<<"\033[32m done\033[0m"<<std::endl;
	
}


std::tuple<DACE::AlgebraicVector<DACE::DA>, DACE::DA, DACE::DA> 
	DALS::initialize_DA(double Bfactor, double SRPC) const {
		
	std::cout<<"    6) Initializing DA first guess:";
	
	//Retrieve apriori estimates
	DACE::AlgebraicVector<double> state_est_cons=
		this->sd_param_OD.state_apr;
	
	//Check whether covariance is available
	DACE::AlgebraicVector<DACE::DA> state_est(6);
	
	if(this->sd_param_OD.flag_Cov_av==true){
		
		const DACE::AlgebraicMatrix<double> Cov_est=
			this->sd_param_OD.Cov_apr;
		
		//Compute the eigenvectors and eigenvalues
		unsigned int n_var=Cov_est.nrows();
		std::vector<double> lambda_vec(n_var);
		DACE::AlgebraicMatrix<double> V(n_var,n_var);
		std::tie(lambda_vec,V)=Eigenvv(Cov_est);

		DACE::AlgebraicVector<DACE::DA> err_state_eig(6);
		for(unsigned int i=0;i<6;i++){
			err_state_eig[i]=this->sd_param_OD.sigma*
				std::sqrt(lambda_vec.at(i))*DACE::DA(i+1);
		}
		
		//Rotate error in eci reference frame
		DACE::AlgebraicVector<DACE::DA> err_state_eci=V*err_state_eig;
		
		//Update state vector
		state_est=state_est_cons+err_state_eci;
	}
	else{
		if(this->sd_param_OD.frame=="RSW"){
			
			//Create an arbitrary DA error vector in RSW reference frame
			DACE::AlgebraicMatrix<double> Aeci2RSW=
				find_Aeci2RSW(state_est_cons);
			DACE::AlgebraicMatrix<double> Aeci2RSW_state=
				blkdiag(Aeci2RSW,3,2);
			DACE::AlgebraicMatrix<double> ARSW2eci_state=
				Aeci2RSW_state.transpose();
				
			DACE::AlgebraicVector<DACE::DA> err_vec_RSW(6);
			for(unsigned int i=0;i<3;i++){
				err_vec_RSW[i]=
					this->sd_param_OD.domain_pos*DACE::DA(i+1);
				err_vec_RSW[i+3]=
					this->sd_param_OD.domain_vel*DACE::DA(i+4);
			}
			
			//Rotate it and add it to the double state vector
			DACE::AlgebraicVector<DACE::DA> err_vec=
				ARSW2eci_state*err_vec_RSW;	
			state_est=state_est_cons+err_vec;
		}
		else if(this->sd_param_OD.frame=="ECI"){
			
			//Create arbitrary DA vector in the ECI reference frame
			for(unsigned int i=0;i<3;i++){
				state_est[i]=state_est_cons[i]+
					this->sd_param_OD.domain_pos*DACE::DA(i+1);
				state_est[i+3]=state_est_cons[i+3]+
					this->sd_param_OD.domain_vel*DACE::DA(i+4);
			}
		
		}
		else{
			throw std::runtime_error
			("DALS: Invalid reference frame for DA variable definition");
		}
			
	}
			
	//Initialize Bfactor and SRPC
	unsigned int ind_DA=6;
	DACE::DA Bfactor_est, SRPC_est;
	if(this->sd_param_DALS.flag_Bfactor_est==true){
		ind_DA=ind_DA+1;
		//Bfactor_est=Bfactor+0.1*Bfactor*DACE::DA(ind_DA);
		Bfactor_est=Bfactor*0.5*DACE::DA(ind_DA);
	}
	else{
		Bfactor_est=Bfactor;
	}

	if(this->sd_param_DALS.flag_SRPC_est==true){
		ind_DA=ind_DA+1.;
		//SRPC_est=SRPC+0.1*SRPC*DACE::DA(ind_DA);
		SRPC_est=SRPC+0.5*DACE::DA(ind_DA);
	}
	else{
		SRPC_est=SRPC;
	}
	std::cout<<"\033[32m done\033[0m"<<std::endl;
	return std::make_tuple(state_est,Bfactor_est,SRPC_est);
		
}


std::tuple<DACE::AlgebraicVector<double>, double, double,
DACE::DA,std::vector<Observations<DACE::DA>>,double> 
	DALS::LSsolver(double et_est,DACE::AlgebraicVector<DACE::DA>& state_est,
	DACE::DA Bfactor_est,DACE::DA SRPC_est) const{
	
	std::cout<<"    7) Running iterative cycle:"<<std::endl;
	
	//Count the number of sensors
	const unsigned int n_sens=this->sd_param_obs.size();
	
	//Avoid considering the same kernel multiple times
	std::vector<std::string> stations_all;
	
	//Load observer kernels
	for(unsigned int j=0;j<n_sens;j++){
		for(unsigned int i=0;i<this->sd_param_obs.at(j).obs_name.size();i++){
			
			//Check whether the station has already been considered
			bool flag_rep=false;
			for(unsigned int k=0;k<stations_all.size();k++){
				if(stations_all.at(k)==
					this->sd_param_obs.at(j).obs_name.at(i)){
					flag_rep=true;
					break;
				}
			}
			
			if(flag_rep==false){
					
				//Add station to stations_all
				stations_all.push_back(this->sd_param_obs.at(j).obs_name.at(i));
				
				//Load kernels
				std::string bsp_file=this->sd_param_obs.at(j).obs_path.at(i)+
					this->sd_param_obs.at(j).obs_name.at(i)+".bsp";
				std::string tf_file =this->sd_param_obs.at(j).obs_path.at(i)+
					this->sd_param_obs.at(j).obs_name.at(i)+".tf";
				furnsh_c(bsp_file.c_str());
				furnsh_c(tf_file.c_str());
			}
		}
	}

	//Propagation parameters
	Dynamics_param<DACE::DA> param;
	param.build(this->sd_param_dyn,Bfactor_est,SRPC_est);
	
	//Initialize OD_data
	OD_data<DACE::DA> Data_est_t_est;
	Data_est_t_est.et_est=et_est;
	Data_est_t_est.state_est=state_est;
	Data_est_t_est.state_apr=this->sd_param_OD.state_apr;	
	Data_est_t_est.flag_apr="no";

	if(this->sd_param_OD.flag_Cov_av==true){
		Data_est_t_est.Cov_apr=this->sd_param_OD.Cov_apr;
		Data_est_t_est.flag_apr="yes";	
	}
	
	
	//Initial residuals
	DACE::DA Residuals;
	std::vector<Observations<DACE::DA>> Obs_est;
        const clock_t begin_time = clock();

	std::tie(Residuals,Obs_est)= ResidualCalc(Data_est_t_est,this->sd_param_meas.Obs,
		this->sd_param_obs,param,this->sd_param_meas.flag_LTS,
		this->sd_param_Earth.Earth_model);
//	std::cout << " It takes " <<  (float)((clock() - begin_time))/CLOCKS_PER_SEC/60.0 << " minute(s) for ResidualCalc " <<std::endl;
  /*
	Write_res(this->sd_param_meas.Obs,this->sd_param_obs,Residuals,
		Obs_est,this->sd_param_DALS);
	
	std::cout<<"Doneeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"<<std::endl;
	*/
	
	//Initialize variables
	unsigned int iter   = 0;
	double fun_diff     = 1e10;
	double step_diff    = 1e10;
	double grad         = 1e10;
	bool accurate       = false;
	bool cond_fun       = false;
	std::string flag_acc = "low";
	
	const unsigned int nn_var=this->sd_param_DALS.nn_var;
	
	column_vector starting_point(nn_var);
	for(unsigned int i=0;i< nn_var;i++){
		starting_point(i)=0.;
	}
	DACE::AlgebraicVector<double> point(nn_var),point_old(nn_var),point_new(nn_var);
	for(unsigned int i=0;i< nn_var;i++){
		point[i]=starting_point(i);
	}
	
	//Define headers for LS output
	std::cout<<"\t"<<std::setw(5)<<std::right<<"\033[36m"<<"Iter"<<"\033[0m";
	std::cout <<std::setw(12)<<std::right<<"\033[36m"<<"Residual "<<"\033[0m";
	if(this->sd_param_DALS.stop_crit=="gradient"){	
		std::cout<<std::setw(13)<<std::right<<"\033[36m"<<"Residual grad"<<"\033[0m";
	}
	else{
		std::cout<<std::setw(12)<<std::right<<"\033[36m"<<"Step-size diff"<<"\033[0m";
		std::cout<<std::setw(13)<<std::right<<"\033[36m"<<"Residual dif"<<"\033[0m";	
	}
	std::cout<<std::setw(13)<<std::right<<"\033[36m"<<"Map accuracy"<<"\033[0m";
	std::cout<<std::setw(9)<<std::right<<"\033[36m"<<"Accuracy"<<"\033[0m";
	std::cout<<std::endl;
	
	
	//While loop
	int flag_conv;	
	
	//Define lower and upper bounds for find_min
	column_vector ub(nn_var);
	for(unsigned int i=0;i< nn_var;i++){
		ub(i)=1.0;
	}
	column_vector lb(nn_var);
	for(unsigned int i=0;i< nn_var;i++){
		lb(i)=-1.0;
	}
	
	//Save t_CPU
	clock_t begin=clock();
	while ((cond_fun==false && iter<this->sd_param_DALS.max_iter) 
			|| (accurate==false && iter<this->sd_param_DALS.max_iter)){
		
		//Read current value of residuals
		DACE::DA Residuals_old=Residuals;
		point_old=point;
		
		//Define the threshold according to the residuals
		double thr=8.0;
		if(DACE::cons(Residuals)<1e3){
			thr=2.0;
		}
		if(DACE::cons(Residuals)<1e2){
			thr=1.0;
		}

		//Compute gradient of residuals
		DACE::AlgebraicVector<DACE::DA> Poli(1);
		Poli[0] = Residuals;
		DACE::AlgebraicVector<DACE::DA> derPoli(nn_var);
		for(unsigned int i=0; i< nn_var; i++) {
			derPoli[i] = Poli[0].deriv(i+1);
		}

		//Define the flags for while loop
		bool flag_extr=true;
		bool flag_max=false;
		bool flag_eq=false;
		
		unsigned int nn_var_exp=3;
		std::vector<double> ampl_fact_ub(nn_var_exp),ampl_fact_lb(nn_var_exp);
		for(unsigned int i=0;i<nn_var_exp;i++){
			ampl_fact_ub.at(i)=1.0;
			ampl_fact_lb.at(i)=1.0;
		}
		
		//Define the starting point
		column_vector starting_point(nn_var),starting_point_old(nn_var);
		
		//Iteration counter
		unsigned int iter_eq=0;	
		unsigned iter_max;
		if(this->sd_param_OD.flag_adaptive==true){
			while(flag_extr==true && flag_max==false){// && flag_eq==false){
				
				//Update iter counter
				iter_eq=iter_eq+1;
				
				//Identify optimum point in domain
				for(unsigned int i=0;i< nn_var;i++){
					starting_point(i)=0.;
				}
				
				find_min_box_constrained(dlib::bfgs_search_strategy(),
					dlib::objective_delta_stop_strategy(1.e-16),
					obj_J(Poli),obj_J(derPoli),starting_point,lb,ub);

				//Check if you are at the extremes of the domain
				flag_extr=false;
				unsigned int count_max=0;
				unsigned int count_eq=0;
				for(unsigned int i=0; i<nn_var_exp; i++){
					if(starting_point(i)>0.99*ub(i) || starting_point(i)<0.99*lb(i)){
						
						//Extremal condition identified
						flag_extr=true;
						
						//Update upper bd ...
						if(starting_point(i)>0.99*ub(i)){
							if(ampl_fact_ub.at(i)<thr){
								ub(i)=2.0*ub(i);
								ampl_fact_ub.at(i)=2.0*ampl_fact_ub.at(i);
							}
							else{
								flag_max=true;//count_max=count_max+1;
							}
						}
						
						//or update lower bound
						if(starting_point(i)<0.99*lb(i)){
							if(ampl_fact_lb.at(i)<thr){
								lb(i)=2.0*lb(i);
								ampl_fact_lb.at(i)=2.0*ampl_fact_lb.at(i);
							}
							else{
								flag_max=true;//count_max=count_max+1;
							}
						}
						
						//Check if repeated
						if(starting_point(i)==starting_point_old(i)){
							count_eq=count_eq+1;
						}
					}
				}
				/*
				//Determine if max size has been reached in all directions
				if(count_max==nn_var_exp){
					flag_max=true;
				}
				
				//Determine if the algorithm wants to expand in forbidden direction
				if(count_eq>0){
					flag_eq=true;
				}
				*/
				
				/*
				std::cout<<"New point: "<<std::endl;
				for(unsigned int i=0;i<nn_var_exp;i++){
					std::cout<<starting_point(i)<<std::endl;
				}
				std::cout<<"Old point: "<<std::endl;
				for(unsigned int i=0;i<nn_var_exp;i++){
					std::cout<<starting_point_old(i)<<std::endl;
				}
				std::cout<<"Flag_extr: "<<flag_extr<<std::endl;
				std::cout<<"Flag_max:  "<<flag_max<<std::endl;
				std::cout<<"Flag_eq:   "<<flag_eq<<std::endl;
				*/
				
				//Save starting point
				starting_point_old=starting_point;
				
			}
			//std::cout<<"Flag_eq: "<<flag_eq<<std::endl;
			//Re-adjust lb and ub
			for(unsigned int i=0;i< nn_var;i++){
				ub(i)=1.0;
			}
			for(unsigned int i=0;i< nn_var;i++){
				lb(i)=-1.0;
			}
		}
		else{
			
			//Identify optimum point in domain
			for(unsigned int i=0;i< nn_var;i++){
				starting_point(i)=0.;
			}
			
			find_min_box_constrained(dlib::bfgs_search_strategy(),
				dlib::objective_delta_stop_strategy(1.e-16),
				obj_J(Poli),obj_J(derPoli),starting_point,-1.0,1.0);
		}
				
		//Write the point
		for(unsigned int i=0; i< nn_var; i++) {
			point[i] = starting_point(i);
			//point_new[i]=point[i];
			//std::cout<<point[i]<<std::endl;
		}
		//Read optimal residual
		double res_opt=Residuals.eval(point);
		
		//----------------- State vector update ------------------//
		//Update DA state vector
		state_est=state_est+(state_est.eval(point)-
			DACE::cons(state_est));
		
		/*
		//Check dlib point and expand the search domain, if required
		DACE::AlgebraicVector<DACE::DA> d=d.identity();
		for(unsigned int i=0;i<6;i++){
			if(std::abs(point[i])>0.9 && ampl_fact.at(i)<8.0){
				d[i]=2.0*d[i];
				point_new[i]=0.5*point_new[i];
				ampl_fact.at(i)=2.0*ampl_fact.at(i);
			}
			if(std::abs(point[i])<0.5 && ampl_fact.at(i)>1.0){
				d[i]=0.5*d[i];
				point_new[i]=2.0*point_new[i];
				ampl_fact.at(i)=0.5*ampl_fact.at(i);
			}
		}
		//Update state and residuals
		state_est=state_est.eval(d);
		*/
		
		//Update OD_data
		Data_est_t_est.state_est=state_est;
		
		//-------------- Object parameters update ---------------//
		if(this->sd_param_DALS.flag_Bfactor_est==true){
			Bfactor_est=Bfactor_est+(Bfactor_est.eval(point)-
				DACE::cons(Bfactor_est));
		}
		if(this->sd_param_DALS.flag_SRPC_est==true){
			SRPC_est=SRPC_est+(SRPC_est.eval(point)-
				DACE::cons(SRPC_est));
		}
		if(this->sd_param_DALS.flag_Bfactor_est==true || 
			this->sd_param_DALS.flag_SRPC_est==true){	
			
			//Update structure of parameters
			param.update(param,Bfactor_est,
				SRPC_est);	
		}		
		//-------------------------------------------------------//

		//Compute new residuals
		std::tie(Residuals,Obs_est) = 
			ResidualCalc(Data_est_t_est,this->sd_param_meas.Obs,
			this->sd_param_obs, param,this->sd_param_meas.flag_LTS,
			this->sd_param_Earth.Earth_model);

		//Calculate the gradient of the new residuals map
		DACE::AlgebraicVector<DACE::DA> GradRes(nn_var);
		for(unsigned int i=0; i< nn_var; i++) {
			GradRes[i] = Residuals.deriv(i+1);
		}
	
		//Stopping criteria
		if(this->sd_param_DALS.stop_crit=="gradient"){
			grad=vnorm(DACE::cons(GradRes));
			
			//Check condition
			if(grad < this->sd_param_DALS.grad_tol){
				cond_fun=true;
			}
		}
		else{
			fun_diff=std::abs(DACE::cons(Residuals-Residuals_old))/
				(1.+std::abs(DACE::cons(Residuals_old)));
			step_diff=std::abs(DACE::vnorm(point-point_old))/
				(1.+DACE::vnorm(point_old));
			//Check conditions
			if(fun_diff < this->sd_param_DALS.fun_tol || 
				step_diff < this->sd_param_DALS.step_tol){
				cond_fun=true;
			}
		}
		
		//Map accuracy
		double map_diff=std::abs(DACE::cons(res_opt-Residuals))/
				(1.+std::abs(res_opt));
		if(map_diff < this->sd_param_DALS.map_tol){
			accurate=true;
			flag_acc="good";
		}
		else{
			flag_acc="low";
		}


		//Update iter
		iter++;
		
		//Update point
		//point=point_new;
		
		/*
		for(unsigned int i=0;i<6;i++){
			
			//Expand search domain
			if(point[i]>0.9*ub(i)){
				ub(i)=2.0*ub(i);
			}
			if(point[i]<0.9*lb(i)){
				lb(i)=2.0*lb(i);
			}
			
			//Shrink it
			if(point[i]<0.5*ub(i) && point[i]>0 && ub(i)>1.0){
				ub(i)=0.5*ub(i);
			}
			if(point[i]>0.5*lb(i) && point[i]<0 && lb(i)<1.0){
				lb(i)=0.5*lb(i);
			}
			
		}
		*/
		
		//Display data
		std::cout<<"\t"<<std::setw(3)<<std::right<<std::showbase<<iter;
		std::cout<<std::setw(18)<<std::scientific<<std::right<<std::showbase<< 
				std::setprecision(6)<< Residuals.cons();
		
		if(this->sd_param_DALS.stop_crit=="gradient"){		
			std::cout <<std::setw(20)<<std::scientific<<std::right<<std::showbase<< 
				std::setprecision(6)<<grad;
		}
		else{	
			std::cout<<std::setw(20)<<std::scientific<<std::right<<std::showbase<< 
				std::setprecision(6)<<step_diff;
			std::cout<<std::setw(20)<<std::scientific<<std::right<<std::showbase<< 
				std::setprecision(6)<<fun_diff;		
		}
		
		std::cout<<std::setw(20)<<std::scientific<<std::right<<std::showbase<< 
			std::setprecision(6)<<map_diff;
		if(flag_acc=="low"){
			std::cout<<std::setw(12)<<std::right<<"\033[31m"<<flag_acc<<"\033[0m";
		}
		else{
			std::cout<<std::setw(12)<<std::right<<"\033[32m"<<flag_acc<<"\033[0m";
		}
		std::cout<<std::endl;
	}
	clock_t end= clock();
    double elapsed_secs=double(end-begin)/CLOCKS_PER_SEC;

	if(step_diff < this->sd_param_DALS.step_tol){
		std::cout<<"\033[32mLS stopped because the size of the current step is less ";
		std::cout<<"than the selected ";
		std::cout<<"value of the step size tolerance ";
		std::cout<<"\033[32m(possible convergence)\033[0m"<<std::endl;
		flag_conv=1;
	}
	else if(fun_diff < this->sd_param_DALS.fun_tol){
		std::cout<<"\033[32mLS stopped because the change in the objective function is less ";
		std::cout<<"than the selected ";
		std::cout<<"value of the step function tolerance ";
		std::cout<<"(possible convergence)\033[0m"<<std::endl;
		flag_conv=2;
	}
	else if(grad < this->sd_param_DALS.grad_tol){
		std::cout<<"\033[32mLS stopped because the norm of the gradient of the objective function is less ";
		std::cout<<"than the selected ";
		std::cout<<"value of the gradient tolerance ";
		std::cout<<"(possible convergence)\033[0m"<<std::endl;
		flag_conv=3;
	}
	else if (iter >= this->sd_param_DALS.max_iter){
		std::cout<<"\033[33mLS stopped because it reached the imposed maximum ";
		std::cout<<"number of iterations\033[0m"<<std::endl;
		flag_conv=0;
	}
	
	//Save the DALS parameters in a dedicated file
	std::ofstream outdata;
	std::string filename_DALS = this->sd_param_DALS.output.filename_DALS;
	outdata.open(filename_DALS,std::ios::trunc);
	outdata.precision(17);
	
	outdata<<"DALS log file"<<std::endl;
	outdata<<"Order: "<<this->sd_param_DALS.order<<std::endl;
	outdata<<"Number of iterations: "<<iter<<std::endl;
	outdata<<"Computational time (s): "<<elapsed_secs<<std::endl;
	outdata<<"Convergence flag: "<<flag_conv<<std::endl;
	outdata<<"Residual: "<<Residuals.cons()<<std::endl;
	outdata.close();
	
	
	//Extract the constant part
	DACE::AlgebraicVector<double> state_est_db=DACE::cons(state_est);
	double Bfactor_est_db = DACE::cons(Bfactor_est);
	double SRPC_est_db  = DACE::cons(SRPC_est);
	
	
	//Unload and remove observer kernels for each sensor
	for(unsigned int j=0;j<n_sens;j++){
		for(unsigned int i=0;i<this->sd_param_obs.at(j).obs_name.size();i++){
			
			//Check whether the station is in the list
			bool flag_rem=false;
			for(unsigned int k=0;k<stations_all.size();k++){
				if(stations_all.at(k)==
					this->sd_param_obs.at(j).obs_name.at(i)){
					flag_rem=true;
					
					//Remove element
					stations_all.erase(stations_all.begin()+k);
					
					break;
				}
			}
			
			if(flag_rem==true){
				std::string bsp_file=this->sd_param_obs.at(j).obs_path.at(i)+
					this->sd_param_obs.at(j).obs_name.at(i)+".bsp";
				std::string tf_file =this->sd_param_obs.at(j).obs_path.at(i)+
					this->sd_param_obs.at(j).obs_name.at(i)+".tf";

				unload_c(bsp_file.c_str());
				unload_c(tf_file.c_str());
				
				remove(bsp_file.c_str());
				remove(tf_file.c_str());
			}
		}
	}
	
	return std::make_tuple(state_est_db,Bfactor_est_db,SRPC_est_db,Residuals,
		Obs_est,flag_conv);
}


std::tuple<DACE::AlgebraicMatrix<double>,DACE::AlgebraicMatrix<double>, bool>
	DALS::compute_covariance( DACE::DA& FinalRes) const{
	
	std::cout<<"    7) Computing covariance:";

DACE::AlgebraicVector<double> scale(6,0.0);
for(int i=0;i<3;i++){
	scale[i]=this->sd_param_OD.domain_pos;
	scale[i+3]=this->sd_param_OD.domain_vel;
}
DACE::AlgebraicVector<DACE::DA> sample(6);
for(int i=0;i<6;i++){
	sample[i] = DACE::DA(i+1)/scale.at(i);
}

std::cout << sample << std::endl;
std::cout << scale[0] << " " << scale[2] << std::endl;
std::cout << FinalRes << std::endl;
FinalRes = FinalRes.eval(sample);

	std::cout<<FinalRes<<std::endl;
	
	//Number of variables
	const unsigned int nn_var=this->sd_param_DALS.nn_var;
	
	//Intialize output
	DACE::AlgebraicMatrix<double> Cov_state(6,6);
	DACE::AlgebraicMatrix<double> Cov_param(nn_var-6,nn_var-6);
	bool flag_dec=false;
	
	DACE::AlgebraicVector<DACE::DA> GradRes(nn_var);
	for(unsigned int i=0; i< nn_var; i++) {
		GradRes[i] = FinalRes.deriv(i+1);
	}
	
	DACE::AlgebraicMatrix<DACE::DA> Hessian(nn_var,nn_var);
	for(unsigned int i=0;i< nn_var;i++){
		for(unsigned int j=0;j< nn_var;j++){
			Hessian.at(i,j)=(GradRes[i]).deriv(j+1);
                        std::cout<<std::setprecision(20) <<Hessian.at(i,j).cons() << " " ;
		}
                std::cout<<std::endl;
	}
 std::cout<<std::endl;
	DACE::AlgebraicMatrix<DACE::DA> Normal=0.5*Hessian;
	
	//Convert to m and m/s
	DACE::AlgebraicMatrix<double> Normal_scaled(nn_var,nn_var);
	for(unsigned int i=0;i< nn_var;i++){
		for(unsigned int j=0;j< nn_var;j++){
			Normal_scaled.at(i,j)=DACE::cons(Normal.at(i,j))*1e6;
		}
	}
		
	//Invert the matrix
	double determ=DACE::det(Normal_scaled);
	/*
	for(unsigned int i=0;i < nn_var;i++){
		for(unsigned int j=0;j < nn_var;j++){
			std::cout<<std::setw(10)<<std::showpos<<
				std::setprecision(10)<<Normal_scaled.at(i,j)<<" ";
		}
		std::cout<<std::endl;
	}
	*/
	if(std::abs(determ)>0){
		DACE::AlgebraicMatrix<double> Covariance_scaled=
			DACE::inv(Normal_scaled);
		
		//Reconvert to km and km/s and multiply per MSRE
		unsigned int n_meas=0;
		for(unsigned int i=0;i<this->sd_param_obs.size();i++){
			unsigned int n_inst=this->sd_param_meas.Obs.at(i).et_vect.size();
			unsigned int n_meas_inst=this->sd_param_meas.Obs.at(i).n_meas_inst;
			n_meas = n_meas + n_inst*n_meas_inst;
		}
		
		double MSRE=DACE::cons(FinalRes)/(n_meas - nn_var);
		std::cout << DACE::cons(FinalRes) << " " << n_meas << " " << nn_var << std::endl;
		DACE::AlgebraicMatrix<double> Covariance_full(nn_var,nn_var);
		for(unsigned int i=0;i < nn_var;i++){
			for(unsigned int j=0;j < nn_var;j++){
				Covariance_full.at(i,j)=MSRE*Covariance_scaled.at(i,j)*1e6;
			}
		}
		
		//Extract covariance matrix for state vector
		for(unsigned int i=0;i<6;i++){
			for(unsigned int j=0;j<6;j++){
				Cov_state.at(i,j)=Covariance_full.at(i,j);
std::cout<<std::setprecision(20) <<Cov_state.at(i,j) << " " ;
			}
std::cout<< " " << std::endl;
		}
		
		//-------------------------------------------------------------------------------//
		//--------------- Object parameters and associated uncertainty ------------------//
		//-------------------------------------------------------------------------------//
		if(nn_var>6){
			for(unsigned int i=0;i < nn_var-6;i++){
				for(unsigned int j=0;j < nn_var-6;j++){
					Cov_param.at(i,j)=Covariance_full.at(i+6,j+6);
				
                                }
			
                        }
			
		}
		std::cout<<"\033[32m done\033[0m"<<std::endl;
	}
	else{
		std::cout<<"\n\t\033[33mThe normal matrix has det=0,";
		std::cout<<" residuals insensitive to obj param"<<"\033[0m"<<std::endl;
		
		//Change flag
		flag_dec=true;
		
		//Limit the analysis to the state vector
		DACE::AlgebraicMatrix<double> Normal_scaled_state(6,6);
		for(unsigned int i=0;i<6;i++){
			for(unsigned int j=0;j<6;j++){
				Normal_scaled_state.at(i,j)=
					Normal_scaled.at(i,j);
			}
		}
		
		DACE::AlgebraicMatrix<double> Covariance_scaled_state=
			DACE::inv(Normal_scaled_state);
		
		//Reconvert to km and km/s and multiply per MSRE
		unsigned int n_meas=0;
		for(unsigned int i=0;i<this->sd_param_obs.size();i++){
			unsigned int n_inst=this->sd_param_meas.Obs.at(i).et_vect.size();
			unsigned int n_meas_inst=this->sd_param_meas.Obs.at(i).n_meas_inst;
			n_meas = n_meas + n_inst*n_meas_inst;
		}
		double MSRE=DACE::cons(FinalRes)/(n_meas - nn_var);
		
		//Extract covariance matrix for state vector
		for(unsigned int i=0;i < 6;i++){
			for(unsigned int j=0;j < 6;j++){
				Cov_state.at(i,j)=MSRE*Covariance_scaled_state.at(i,j)*1e6;
			}
		}
	}
	
	//Unload Earth fixed high precision kernels
	if(this->sd_param_Earth.Earth_model=="ITRF93"){	
		unload_c((this->sd_param_Earth.kernels.path+
			this->sd_param_Earth.kernels.ITRF93).c_str());
		unload_c((this->sd_param_Earth.kernels.path+
			this->sd_param_Earth.kernels.ITRF93_assoc).c_str());
		furnsh_c((this->sd_param_Earth.kernels.path+
			this->sd_param_Earth.kernels.IAU_assoc).c_str());
	}
	
	return std::make_tuple(Cov_state,Cov_param,flag_dec);
}

void DALS::write_results(const double et_est, 
	const DACE::AlgebraicVector<double>& state_est, 
	const DACE::AlgebraicMatrix<double>& Cov_est, const double Bfactor_est, 
	const double SRPC_est, const bool flag_dec,
	const DACE::AlgebraicMatrix<double>& Cov_param,
	const double flag_conv, const DACE::DA& FinalRes,
	const std::vector<Observations<DACE::DA>>& Obs_est) const{
	
	std::cout<<"    8) Writing the results:";	
	
	//Write OPM file
	Write_OPM(et_est, state_est, Cov_est, this->sd_param_DALS,
		Bfactor_est,SRPC_est,Cov_param,flag_dec);
		
	//Write residuals table and map
	Write_res(this->sd_param_meas.Obs,this->sd_param_obs, FinalRes,
		Obs_est,this->sd_param_DALS);
	
	std::cout<<"\033[32m done\033[0m"<<std::endl;	
		
}
					
			
//********************************************************************************
//************************** Public member functions *****************************
//********************************************************************************

DALS::DALS(const std::string& json_LS){
	
	std::cout<<"    1) Initializing solver:";
	
	//Checking for exceptions
	//input_except(json_LS);
	
	//Open json file
	std::ifstream inputFile(json_LS, std::ifstream::binary);
	
	// Read the Json value
    Json::Value root;
    inputFile >> root;
	inputFile.close();
	
	/**********************************************************/
	/******************* Measurements parameters **************/
	/**********************************************************/
	const unsigned int n_sens = root["Observer"]["n_sensors"].asInt();
	Measurements_param<double> param_meas;
	
	//Store all the TDMs
	std::vector<unsigned int> sensor_index;
	for(unsigned int i=0;i<n_sens;i++){
		
		//Build string of current sensor
		std::string sensor = "sensor_"+std::to_string(i+1);
		
		//Count number of TDMs
		unsigned int n_TDM = root["Observer"][sensor]["n_TDM"].asInt();
		
		//Save TDM name of current sensor
		std::string path_meas   = root["Observer"]
			[sensor]["TDM"]["path"].asString();
			
		for(unsigned int j=0;j<n_TDM;j++){
			std::string TDM_single = "TDM_"+std::to_string(j+1);
			std::string name_meas	= root["Observer"]
				[sensor]["TDM"][TDM_single].asString();
			param_meas.filename_meas.push_back(path_meas+name_meas);
			
			sensor_index.push_back(i);
			
			//Read flag for LTS corrections
			param_meas.flag_LTS.push_back(root["Observer"][sensor]
				["flag_LTS"].asBool());
		}
	}
	this->sd_param_meas=param_meas;
	
	
	/**********************************************************/
	/********************* Observer parameters ****************/
	/**********************************************************/
	std::vector<Observer_param> param_obs;
	
	//Scan all the involved sensors
	for(unsigned int ii=0;ii<sensor_index.size();ii++){
		
		//Sensor index
		unsigned int i=sensor_index.at(ii);
		
		//Build string of current sensor
		std::string sensor = "sensor_"+std::to_string(i+1);
		
		//Initialize data
		Observer_param param_obs_single;

		//Read sensor type
		param_obs_single.type  = root["Observer"][sensor]["type"].asString();
		
		param_obs_single.obs_name.push_back(root["Observer"][sensor]
			["participant_1"]["name"].asString());
		param_obs_single.obs_path.push_back(root["Observer"][sensor]
			["participant_1"]["path"].asString());
			
		//Create directory
		std::string create_Obs = "mkdir -p "+root["Observer"][sensor]
			["participant_1"]["path"].asString();
		system(create_Obs.c_str());

		param_obs_single.lat.push_back(root["Observer"][sensor]
			["participant_1"]["latitude"].asDouble());
		param_obs_single.lon.push_back(root["Observer"][sensor]
			["participant_1"]["longitude"].asDouble());
		param_obs_single.alt.push_back(root["Observer"][sensor]
			["participant_1"]["altitude"].asDouble());
		
		//Limits
		param_obs_single.el_min.push_back(root["Observer"][sensor]
			["participant_1"]["limits"]["el_min"].asDouble());
		
		//Sensor setup
		Sensor_setup sensor_setup;
		
		if(param_obs_single.type=="optical"){
			
			//Right ascension
			sensor_setup.ra.availability=root["Observer"][sensor]
				["sensor_setup"]["right_ascension"]["availability"].asBool();
			sensor_setup.ra.accuracy=root["Observer"][sensor]
				["sensor_setup"]["right_ascension"]["accuracy"].asDouble();
			sensor_setup.ra.TDM_name=root["Observer"][sensor]
				["sensor_setup"]["right_ascension"]["TDM_name"].asString();
				
			//Declination
			sensor_setup.decl.availability=root["Observer"][sensor]
				["sensor_setup"]["declination"]["availability"].asBool();
			sensor_setup.decl.accuracy=root["Observer"][sensor]
				["sensor_setup"]["declination"]["accuracy"].asDouble();	
			sensor_setup.decl.TDM_name=root["Observer"][sensor]
				["sensor_setup"]["declination"]["TDM_name"].asString();
			
			//Save sensor_setup
			param_obs_single.sensor_setup=sensor_setup;
			
			//Limits
			param_obs_single.sunEl_lim=root["Observer"][sensor]
				["participant_1"]["limits"]["sunEl_lim"].asDouble();
			param_obs_single.moonSep_lim=root["Observer"][sensor]
				["participant_1"]["limits"]["moonSep_lim"].asDouble();
			param_obs_single.magn_lim=root["Observer"][sensor]
				["participant_1"]["limits"]["magn_lim"].asDouble();	

		}
		else if(param_obs_single.type=="radar"){	
			
			//Azimuth
			sensor_setup.az.availability=root["Observer"][sensor]
				["sensor_setup"]["azimuth"]["availability"].asBool();
			sensor_setup.az.accuracy=root["Observer"][sensor]
				["sensor_setup"]["azimuth"]["accuracy"].asDouble();
			sensor_setup.az.TDM_name=root["Observer"][sensor]
				["sensor_setup"]["azimuth"]["TDM_name"].asString();
				
			//Elevation
			sensor_setup.el.availability=root["Observer"][sensor]
				["sensor_setup"]["elevation"]["availability"].asBool();
			sensor_setup.el.accuracy=root["Observer"][sensor]
				["sensor_setup"]["elevation"]["accuracy"].asDouble();
			sensor_setup.el.TDM_name=root["Observer"][sensor]
				["sensor_setup"]["elevation"]["TDM_name"].asString();
				
			//Doppler
			sensor_setup.doppler.availability=root["Observer"][sensor]
				["sensor_setup"]["doppler"]["availability"].asBool();
			sensor_setup.doppler.accuracy=root["Observer"][sensor]
				["sensor_setup"]["doppler"]["accuracy"].asDouble();
			sensor_setup.doppler.TDM_name=root["Observer"][sensor]
				["sensor_setup"]["doppler"]["TDM_name"].asString();
				
			//Range
			sensor_setup.range.availability=root["Observer"][sensor]
				["sensor_setup"]["range"]["availability"].asBool();
			sensor_setup.range.accuracy=root["Observer"][sensor]
				["sensor_setup"]["range"]["accuracy"].asDouble();
			sensor_setup.range.TDM_name=root["Observer"][sensor]
				["sensor_setup"]["range"]["TDM_name"].asString();
				
			//Radar frequency
			sensor_setup.frequency=root["Observer"][sensor]
				["sensor_setup"]["frequency"].asDouble();
				
			//Radar configuration
			sensor_setup.config=root["Observer"][sensor]
				["sensor_setup"]["configuration"].asString();
			
			//Save sensor_setup
			param_obs_single.sensor_setup=sensor_setup;
			
			//Read data of second observatory in case of bistatic config
			if(sensor_setup.config=="bistatic"){	
			
				param_obs_single.obs_name.push_back(root["Observer"][sensor]
					["participant_2"]["name"].asString());
				param_obs_single.obs_path.push_back(root["Observer"][sensor]
					["participant_2"]["path"].asString());
				
				//Create directory
				std::string create_Obs = "mkdir -p "+root["Observer"][sensor]
					["participant_2"]["path"].asString();
				system(create_Obs.c_str());
		
				param_obs_single.lat.push_back(root["Observer"][sensor]
					["participant_2"]["latitude"].asDouble());
				param_obs_single.lon.push_back(root["Observer"][sensor]
					["participant_2"]["longitude"].asDouble());
				param_obs_single.alt.push_back(root["Observer"][sensor]
					["participant_2"]["altitude"].asDouble());
				
				//Limits
				param_obs_single.el_min.push_back(root["Observer"][sensor]
					["participant_2"]["limits"]["el_min"].asDouble());	
			}
		}
		else{
			throw std::runtime_error("DALS: Invalid observer type");
		}
		
		//Save data
		param_obs.push_back(param_obs_single);
	}
	this->sd_param_obs=param_obs;
	
	
	/**********************************************************/
	/******************** Earth parameters ********************/
	/**********************************************************/
	Earth_param param_Earth;
	
	param_Earth.Earth_model        = root["Earth"]["Earth_model"].asString();
	param_Earth.kernels.path       = root["Earth"]["kernels"]["path"].asString();
	param_Earth.kernels.pck        = root["Earth"]["kernels"]["pck"].asString();
	param_Earth.kernels.IAU_assoc  = root["Earth"]["kernels"]
		["IAU_assoc"].asString();
	
	if(param_Earth.Earth_model=="ITRF93"){	
		param_Earth.kernels.ITRF93        = root["Earth"]["kernels"]
			["ITRF93"].asString();
		param_Earth.kernels.ITRF93_assoc  = root["Earth"]["kernels"]
			["ITRF93_assoc"].asString();
		furnsh_c((param_Earth.kernels.path+
			param_Earth.kernels.ITRF93).c_str());
		furnsh_c((param_Earth.kernels.path+
			param_Earth.kernels.ITRF93_assoc).c_str());
	}
	else{
		furnsh_c((param_Earth.kernels.path+
			param_Earth.kernels.IAU_assoc).c_str());
	}
	this->sd_param_Earth = param_Earth;
	
	
	/**********************************************************/
	/***************** State estimates parameters *************/
	/**********************************************************/
	Estimate_param param_OD;
	std::string path_OD     = root["OD"]["path"].asString();
	std::string name_OD     = root["OD"]["name"].asString();
	param_OD.filename_state = path_OD + name_OD;

	//Covariance
	param_OD.flag_Cov_av    = root["OD"]["Covariance"]["flag_Cov_av"].asBool();
	
	//Search type
	param_OD.flag_adaptive  = root["OD"]["Covariance"]["flag_adaptive"].asBool();
	
	if(param_OD.flag_Cov_av==false){
		
		//Reference frame for coordinates definition
		param_OD.frame = root["OD"]["Covariance"]["frame"].asString();
		
		//Define the search space
		if(param_OD.flag_adaptive==false){
			param_OD.domain_pos = root["OD"]["Covariance"]["domain_pos"].asDouble();
			param_OD.domain_vel = root["OD"]["Covariance"]["domain_vel"].asDouble();
		}
		else{
			param_OD.domain_pos = 1.0;
			param_OD.domain_vel = 1e-2;
		}
	}
	else{
		param_OD.sigma = root["OD"]["Covariance"]["sigma"].asDouble();
	}
	param_OD.flag_Bfactor_av  = root["OD"]["flag_Bfactor_av"].asBool();
	param_OD.flag_SRPC_av   = root["OD"]["flag_SRPC_av"].asBool();
	this->sd_param_OD=param_OD;
	
	/**********************************************************/
	/************************ DALS variables ******************/
	/**********************************************************/
	DALS_param param_DALS;
	param_DALS.order     = root["DALS"]["order"].asInt();
	param_DALS.max_iter  = root["DALS"]["max_iter"].asInt();
	param_DALS.fun_tol   = root["DALS"]["fun_tol"].asDouble();
	param_DALS.step_tol  = root["DALS"]["step_tol"].asDouble();
	param_DALS.grad_tol  = root["DALS"]["grad_tol"].asDouble();
	param_DALS.map_tol   = root["DALS"]["map_tol"].asDouble();
	param_DALS.stop_crit = root["DALS"]["stop_crit"].asString();

	std::string date_est = root["DALS"]["date_est"].asString();
	ConstSpiceChar *char_tmp=date_est.c_str();
	SpiceDouble et_est;
	str2et_c(char_tmp,&et_est);
	param_DALS.et_est    = et_est;

	unsigned int nn_var_start=6;
	param_DALS.flag_Bfactor_est = root["DALS"]["flag_Bfactor_est"].asBool();
	if(this->sd_param_OD.flag_Bfactor_av==false){
		param_DALS.flag_Bfactor_est=true;
	}
	if(param_DALS.flag_Bfactor_est==true){
		nn_var_start=nn_var_start+1;
	}
	param_DALS.flag_SRPC_est  = root["DALS"]["flag_SRPC_est"].asBool();
	if(this->sd_param_OD.flag_SRPC_av==false){
		param_DALS.flag_SRPC_est=true;
	}
	if(param_DALS.flag_SRPC_est==true){
		nn_var_start=nn_var_start+1;
	}
	param_DALS.nn_var=nn_var_start;
	
	std::string path_out = root["DALS"]["output"]["path"].asString();
	std::string name_OPM       = root["DALS"]["output"]["name_OPM"].asString();
	std::string name_res_table = root["DALS"]["output"]["name_res_table"].asString();
	std::string name_res_map   = root["DALS"]["output"]["name_res_map"].asString();
	
	
	param_DALS.output.filename_OPM       = path_out+name_OPM;
	param_DALS.output.filename_res_table = path_out+name_res_table;
	param_DALS.output.filename_res_map   = path_out+name_res_map;
	param_DALS.output.filename_DALS      = path_out+"DALS_data.txt";
	param_DALS.output.obj_name     = root["DALS"]["output"]["obj_name"].asString();
	param_DALS.output.obj_ID       = root["DALS"]["output"]["obj_ID"].asDouble();
	param_DALS.output.flag_Cov     = root["DALS"]["output"]["flag_Cov"].asBool();
	param_DALS.output.flag_KeplPar = root["DALS"]["output"]["flag_KeplPar"].asBool();
	this->sd_param_DALS=param_DALS;

	DACE::DA::init(this->sd_param_DALS.order,
		this->sd_param_DALS.nn_var);
	
	/**********************************************************/
	/******************** Dynamics variables ******************/
	/**********************************************************/
	Dynamics_param<double> param_dyn;
	param_dyn.method=root["Dynamics"]["method"].asString();
	if(param_dyn.method=="AIDA"){
		
		param_dyn.gravOrd = root["Dynamics"]["AIDAparam"]
			["gravmodel"]["order"].asInt();
		std::string gravmodel_path = root["Dynamics"]["AIDAparam"]
			["gravmodel"]["path"].asString();
		std::string gravmodel_name = root["Dynamics"]["AIDAparam"]
			["gravmodel"]["name"].asString();	
		param_dyn.gravmodel = gravmodel_path+gravmodel_name;
		param_dyn.AIDA_flags = {{root["Dynamics"]["AIDAparam"]["flag_drag"].asInt(),
			root["Dynamics"]["AIDAparam"]["flag_SRP"].asInt(),
			root["Dynamics"]["AIDAparam"]["flag_thirdbody"].asInt()
		}};
	}
	else{
		
		param_dyn.SunMoonFlag = root["Dynamics"]["SADAparam"]
			["flag_SunMoon"].asInt();
		param_dyn.SRPFlag = root["Dynamics"]["SADAparam"]
			["flag_SRP"].asInt();
		param_dyn.dragFlag = root["Dynamics"]["SADAparam"]
			["flag_drag"].asInt();
		param_dyn.shortPeriodicsFlag = root["Dynamics"]["SADAparam"]
			["flag_shortPeriodics"].asInt();
		param_dyn.tesseralFlag = root["Dynamics"]["SADAparam"]
			["flag_tesseral"].asInt();
		param_dyn.flag_linear = root["Dynamics"]["SADAparam"]
			["flag_linear"].asBool();
	}
	param_dyn.tol = root["Dynamics"]["tolerance"].asDouble();
	this->sd_param_dyn=param_dyn;
	
	std::cout<<"\033[32m done\033[0m"<<std::endl;
	
}


std::tuple<double,DACE::AlgebraicVector<double>,DACE::AlgebraicMatrix<double>,
	double, double, DACE::AlgebraicMatrix<double>, DACE::DA, bool> DALS::run(){
	
	//Generate kernels for the observer
	generate_obs_kernels();
	
	//Read measurement file
	read_measurements(this->sd_param_meas.filename_meas);
		
	//Read state and param estimates
	double et_data;
	DACE::AlgebraicVector<double> state_data(6);
	DACE::AlgebraicMatrix<double> Cov_data(6,6);
	double Bfactor;
	double SRPC;
	std::tie(et_data, state_data, Cov_data,Bfactor,SRPC)=
		read_state_data(this->sd_param_OD.filename_state);
	
	//Propagate state estimates to et_est
	propagate_apr_data(et_data,state_data,Cov_data,
		this->sd_param_DALS.et_est, Bfactor, SRPC);
	
	//Initialize DA variables
	DACE::AlgebraicVector<DACE::DA> state_est_DA(6);
	DACE::DA Bfactor_est_DA;
	DACE::DA SRPC_est_DA;
	std::tie(state_est_DA, Bfactor_est_DA,SRPC_est_DA) = 
		initialize_DA(Bfactor, SRPC);

	//Run LS solver
	DACE::AlgebraicVector<double> state_est(6);
	double Bfactor_est;
	double SRPC_est;
	DACE::DA Residuals;
	std::vector<Observations<DACE::DA>> Obs_est;
	bool flag_conv;
	std::tie(state_est, Bfactor_est, SRPC_est, Residuals, Obs_est, flag_conv)= 
		LSsolver(this->sd_param_DALS.et_est,state_est_DA, Bfactor_est_DA, 
		SRPC_est_DA);
	
	//Compute covariance
	const unsigned int nn_var=this->sd_param_DALS.nn_var;
	DACE::AlgebraicMatrix<double> Cov_est(6,6);
	DACE::AlgebraicMatrix<double> Cov_param_est(nn_var-6,nn_var-6);
	bool flag_dec;
	std::tie(Cov_est,Cov_param_est,flag_dec) = compute_covariance(Residuals);
	
	//Write the results
	write_results(this->sd_param_DALS.et_est, state_est, Cov_est,
		Bfactor_est,SRPC_est, flag_dec, Cov_param_est, flag_conv, Residuals,
		Obs_est);
	
	return std::make_tuple(this->sd_param_DALS.et_est,state_est, 
		Cov_est,Bfactor_est,SRPC_est, 
		Cov_param_est,Residuals,flag_dec);			
}


}

#endif