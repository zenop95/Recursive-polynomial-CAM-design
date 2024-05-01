#ifndef SimRoutines_H_
#define SimRoutines_H_

#include <tuple>

#include<dace/dace.h>
#include<cspice/SpiceUsr.h>

#include "Standards.hpp"
#include "Objects.hpp"

void Compute_obs_param(const std::string& json_sim){
	
	//Random numbers generator
	std::random_device rand_dev;
	std::mt19937 generator(rand_dev());
		
	//Open json sim file
	std::ifstream inputFile(json_sim, std::ifstream::binary);
    Json::Value root;
    inputFile >> root;
	inputFile.close();
	
	//Read catalogue type
	std::string cat_name  = root["Object"]["catalogue_name"].asString();
	
	//Count the number of sensors
	const unsigned int n_sens = root["Observer"]["n_sensors"].asDouble();
	for(unsigned int i=0;i<n_sens;i++){

		//Build string of current sensor
		std::string sensor = "sensor_"+std::to_string(i+1);
		
		//Read sensor type
		std::string sens_type = root["Observer"][sensor]["type"].asString();
			
		double Dt;
		unsigned int N_obs_max;
		if(cat_name=="GEO.tle" || cat_name=="HEO.tle"){
			
			//Define the gap between observations
			const unsigned int Dt_min  = 20.;
			const unsigned int Dt_max    = 120.;
			std::uniform_real_distribution<>  distr_Dt(Dt_min,Dt_max);
			Dt = distr_Dt(generator);
			
			//Define the number of observations
			const unsigned N_min = 6;
			const unsigned N_max = 15;
			std::uniform_int_distribution<int> distr_N(N_min,N_max);
			N_obs_max= distr_N(generator);
		}
		else{
			if(sens_type=="optical"){
				
				//Define the gap between observations
				const unsigned int Dt_min  = 5.;
				const unsigned int Dt_max    = 10.;
				std::uniform_real_distribution<>  distr_Dt(Dt_min,Dt_max);
				Dt = distr_Dt(generator);
				
				//Define the number of observations
				const unsigned N_min = 20;
				const unsigned N_max = 30;
				std::uniform_int_distribution<int> distr_N(N_min,N_max);
				N_obs_max= distr_N(generator);
				
			}
			else{
				
				//Define the gap between observations
				const unsigned int Dt_min  = 5.;
				const unsigned int Dt_max    = 10.;
				std::uniform_real_distribution<>  distr_Dt(Dt_min,Dt_max);
				Dt = distr_Dt(generator);
				
				//Define the number of observations
				const unsigned N_min = 20;
				const unsigned N_max = 30;
				std::uniform_int_distribution<int> distr_N(N_min,N_max);
				N_obs_max= distr_N(generator);
					
			}
		}
		
		//Update json file
		root["Measurements"][sensor]["Dt"]=Dt;
		root["Measurements"][sensor]["N_obs_max"]=N_obs_max;
	}
		
		
	Json::StyledWriter writer;
	std::string resultString=writer.write(root);
    std::ofstream outputFile(json_sim);
    outputFile << resultString << std::endl;
	outputFile.close();
}


std::tuple<DACE::AlgebraicVector<double>, DACE::AlgebraicVector<double>>
	Compute_error(const DACE::AlgebraicVector<double>& state_est,
	const double et_est,const double Bfactor_est, const double SRPC_est,
	const DACE::AlgebraicVector<double>& state_ref, const double et_ref, 
	const double et_comp,const std::string& json_sim,
	const std::string& json_LS){
		
	//Read from DALS json file
	std::ifstream inputFile(json_LS, std::ifstream::binary);
	Json::Value root;
	inputFile >> root;
	
	Dynamics_param<double> param;
	
	param.method=root["Dynamics"]["method"].asString();
	if(param.method=="AIDA"){
		param.gravOrd = root["Dynamics"]["AIDAparam"]
			["gravmodel"]["order"].asInt();
		std::string gravmodel_path = root["Dynamics"]["AIDAparam"]
			["gravmodel"]["path"].asString();
		std::string gravmodel_name = root["Dynamics"]["AIDAparam"]
			["gravmodel"]["name"].asString();	
		param.gravmodel = gravmodel_path+gravmodel_name;
		param.AIDA_flags = {{
			root["Dynamics"]["AIDAparam"]["flag_drag"].asInt(),
			root["Dynamics"]["AIDAparam"]["flag_SRP"].asInt(),
			root["Dynamics"]["AIDAparam"]["flag_thirdbody"].asInt()
		}};
	}
	else if(param.method=="SADA"){
		
		param.SunMoonFlag = root["Dynamics"]["SADAparam"]
			["flag_SunMoon"].asInt();
		param.SRPFlag = root["Dynamics"]["SADAparam"]
			["flag_SRP"].asInt();
		param.dragFlag = root["Dynamics"]["SADAparam"]
			["flag_drag"].asInt();
		param.tesseralFlag = root["Dynamics"]["SADAparam"]
			["flag_tesseral"].asInt();
		param.shortPeriodicsFlag = root["Dynamics"]["SADAparam"]
			["flag_shortPeriodics"].asInt();
		param.flag_linear = root["Dynamics"]["SADAparam"]
			["flag_linear"].asBool();
		
	}
	param.tol=root["Dynamics"]["tolerance"].asDouble();
	
	//Read from sim json file
	std::ifstream inputFile_sim(json_sim, std::ifstream::binary);
	Json::Value root_sim;
	inputFile_sim >> root_sim;
	const double mass   = root_sim["Object"]["mass"].asDouble();
	const double A_drag = root_sim["Object"]["A_drag"].asDouble();
	const double Cd     = root_sim["Object"]["Cd"].asDouble();
	const double A_srp  = root_sim["Object"]["A_srp"].asDouble();
	const double Cr     = root_sim["Object"]["Cr"].asDouble();

	const double Bfactor = Cd*A_drag/mass;
	const double SRPC  = Cr*A_srp/mass;
	
	Dynamics_param<double> param_real;
	param_real.method=root_sim["Dynamics"]["method"].asString();
	if(param_real.method=="AIDA"){
		param_real.gravOrd = root_sim["Dynamics"]["AIDAparam"]
			["gravmodel"]["order"].asInt();
		std::string gravmodel_path = root_sim["Dynamics"]["AIDAparam"]
			["gravmodel"]["path"].asString();
		std::string gravmodel_name = root_sim["Dynamics"]["AIDAparam"]
			["gravmodel"]["name"].asString();	
		param_real.gravmodel = gravmodel_path+gravmodel_name;
		param_real.AIDA_flags = {{
			root_sim["Dynamics"]["AIDAparam"]["flag_drag"].asInt(),
			root_sim["Dynamics"]["AIDAparam"]["flag_SRP"].asInt(),
			root_sim["Dynamics"]["AIDAparam"]["flag_thirdbody"].asInt()
		}};
	}
	else if(param_real.method=="SADA"){
		
		param_real.SunMoonFlag = root_sim["Dynamics"]["SADAparam"]
			["flag_SunMoon"].asInt();
		param_real.SRPFlag = root_sim["Dynamics"]["SADAparam"]
			["flag_SRP"].asInt();
		param_real.dragFlag = root_sim["Dynamics"]["SADAparam"]
			["flag_drag"].asInt();
		param_real.tesseralFlag = root_sim["Dynamics"]["SADAparam"]
			["flag_tesseral"].asInt();
		param_real.shortPeriodicsFlag = root_sim["Dynamics"]["SADAparam"]
			["flag_shortPeriodics"].asInt();
		param_real.flag_linear = root_sim["Dynamics"]["SADAparam"]
			["flag_linear"].asBool();
		
	}
	param_real.tol=root_sim["Dynamics"]["tolerance"].asDouble();
	
	//Propagate estimate if needed
	Dynamics_param<double> param_est;
	param_est.build(param,Bfactor_est,SRPC_est);
	
	DACE::AlgebraicVector<double> state_est_t_comp(6);
	if(et_est!=et_comp){
		state_est_t_comp = Propagate_state(et_est,{et_comp},
			state_est,param_est).at(0);
	}
	else{
		state_est_t_comp=state_est;
	}
	
	//Propagate reference if needed
	Dynamics_param<double> param_ref;
	param_ref.build(param_real,Bfactor,SRPC);
	
	DACE::AlgebraicVector<double> state_ref_t_comp(6);
	if(et_ref!=et_comp){
		state_ref_t_comp = Propagate_state(et_ref,{et_comp},
			state_ref,param_ref).at(0);
	}
	else{
		state_ref_t_comp=state_ref;
	}
	
	//Define accuracy
	DACE::AlgebraicVector<double> err_pos(3),err_vel(3);
	for(unsigned int i=0;i<3;i++){
		err_pos[i]=state_ref_t_comp[i]-DACE::cons(state_est_t_comp)[i];
		err_vel[i]=state_ref_t_comp[i+3]-DACE::cons(state_est_t_comp)[i+3];
	}

	SpiceChar output[TIMLEN];

	timout_c(et_comp,TIMFMT,TIMLEN,output);
	std::cout<<"\tDate : "<<output<<std::endl;
	std::cout<<"\t\tPosition error        (km): "<<
		std::scientific <<std::setprecision (6)<<
		err_pos.vnorm()<<std::endl;
	std::cout<<"\t\tVelocity error      (km/s): "<<
		std::scientific <<std::setprecision (6)<<
		err_vel.vnorm()<<std::endl;
			
	
	return std::make_tuple(err_pos,err_vel);

}

void Compute_accuracy (const DACE::AlgebraicVector<double>& state_est,
	const double et_est,const double Bfactor_est, const double SRPC_est,
	const std::string& filename_ref, const std::string& json_sim,
	const std::string& json_LS){
		
	//Read data from filename_ref
	std::ifstream file(filename_ref);
	std::string string_tmp;
	
	//Read et_ref
	std::getline(file, string_tmp);
	double et_ref = std::stod(string_tmp);
	
	//Read state_ref
	DACE::AlgebraicVector<double> state_ref(6);
	std::getline(file, string_tmp);
	std::stringstream ss(string_tmp);
	std::string rr_x,rr_y,rr_z,vv_x,vv_y,vv_z;
	ss>> rr_x >> rr_y >> rr_z >> vv_x >> vv_y >> vv_z;
	state_ref[0] = std::stod(rr_x);
	state_ref[1] = std::stod(rr_y);
	state_ref[2] = std::stod(rr_z);
	state_ref[3] = std::stod(vv_x);
	state_ref[4] = std::stod(vv_y);
	state_ref[5] = std::stod(vv_z);
	
	//Read et_start
	std::getline(file, string_tmp);
	double et_start = std::stod(string_tmp);
	
	
	//Read et_end
	std::getline(file, string_tmp);
	double et_end = std::stod(string_tmp);
	
	file.close();
	
	//Compute error at et_est
	Compute_error(state_est,et_est,Bfactor_est,SRPC_est,
		state_ref,et_ref,et_est,json_sim,json_LS);
		
	//Compute error at et_start
	Compute_error(state_est,et_est,Bfactor_est,SRPC_est,
		state_ref,et_ref,et_start,json_sim,json_LS);
		
	//Compute error at et_end
	Compute_error(state_est,et_est,Bfactor_est,SRPC_est,
		state_ref,et_ref,et_end,json_sim,json_LS);
		
};

void Compare_parameters(const std::string& json_sim, const std::string& json_LS,
	const double Bfactor_est, const double SRPC_est, 
	const DACE::AlgebraicMatrix<double>& Cov_param, const bool flag_dec){
		
	//Read from DALS json file
	std::ifstream inputFile(json_LS, std::ifstream::binary);
	Json::Value root;
	inputFile >> root;
	const bool flag_Bfactor_av = root["OD"]["flag_Bfactor_av"].asBool();
	const bool flag_SRPC_av  = root["OD"]["flag_SRPC_av"].asBool();
	
	const bool flag_Bfactor_est = root["DALS"]["flag_Bfactor_est"].asBool();
	const bool flag_SRPC_est  = root["DALS"]["flag_SRPC_est"].asBool();
	
	//Read from sim json file
	std::ifstream inputFile_sim(json_sim, std::ifstream::binary);
	Json::Value root_sim;
	inputFile_sim >> root_sim;
	const double mass   = root_sim["Object"]["mass"].asDouble();
	const double A_drag = root_sim["Object"]["A_drag"].asDouble();
	const double Cd     = root_sim["Object"]["Cd"].asDouble();
	const double A_srp  = root_sim["Object"]["A_srp"].asDouble();
	const double Cr     = root_sim["Object"]["Cr"].asDouble();

	const double Bfactor = Cd*A_drag/mass;
	const double SRPC  = Cr*A_srp/mass;	
		
	unsigned int index_SRPC=0;
	if(flag_Bfactor_est==true || flag_Bfactor_av==false){
		std::cout<<"\tBfactor "<<std::endl;
		std::cout<<"\t\tReal             : "<<
			Bfactor<<std::endl;
		std::cout<<"\t\tEstimated        : "<<
			Bfactor_est<<std::endl;
		std::cout<<"\t\tError            : "<<
			abs(Bfactor-Bfactor_est)/Bfactor*1e2<<"%"<<std::endl;	
		if(flag_dec==false){
			std::cout<<"\t\tStd deviation    : "<<
				std::sqrt(Cov_param.at(0,0))<<std::endl;
		}
		else{
			std::cout<<"\t\tStd deviation not available"<<std::endl;
		}
		index_SRPC=1;
	}
	
	
	if(flag_SRPC_est==true || flag_SRPC_av==false){
		std::cout<<"\tSRPC "<<std::endl;
		std::cout<<"\t\tReal             : "<<
			SRPC<<std::endl;
		std::cout<<"\t\tEstimated        : "<<
			SRPC_est<<std::endl;
		std::cout<<"\t\tError            : "<<
		abs(SRPC-SRPC_est)/SRPC*1e2<<"%"<<std::endl;
		if(flag_dec==false){
			std::cout<<"\t\tStd deviation    : "<<
				sqrt(Cov_param.at(index_SRPC,index_SRPC))<<std::endl;
		}
		else{
			std::cout<<"\t\tStd deviation not available"<<std::endl;
		}
	}
	
	
	
	
}
#endif