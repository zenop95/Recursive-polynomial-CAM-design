#ifndef LSRoutines_H_
#define LSRoutines_H_

#define       TIMYYY   "YYYY.############::UTC"
#define       TIMFMT_OPM   "YYYY-MM-DDTHR:MN:SC.###### ::UTC ::RND"
#define       TIMLYY   20

#include<tuple>
#include<iostream>
#include<random> 
        
#include <time.h>
#include <chrono>

#include<dace/dace.h>
#include<cspice/SpiceUsr.h>

#include "Objects.hpp"
#include "Standards.hpp"
#include "Propagation.hpp"
#include "SpiceRoutines.hpp"


void Update_json(const std::string json_name, 
	const std::string& date){
	
	//Open json file
	std::ifstream inputFile(json_name, std::ifstream::binary);
	
	// Read the Json value
    Json::Value root;
    inputFile >> root;
	inputFile.close();
	
	//Update epoch
	root["DALS"]["date_est"]=date;
	
	//Update json file
	Json::StyledWriter writer;
	std::string resultString=writer.write(root);
    std::ofstream outputFile(json_name);
    outputFile << resultString << std::endl;
	outputFile.close();
		
}


std::tuple<std::string, double, double, double>
Generate_SUT(const Observer_param& param_obs, const unsigned int i,
	const unsigned int index_obs, const std::string& Earth_model, 
	const double et){
	
	//Write .sut file
	std::ofstream outdata;
	std::string filename_sut = param_obs.obs_path.at(i)+
		param_obs.obs_name.at(i)+".sut";
	outdata.open(filename_sut,std::ios::trunc);
	outdata.precision(17);
	
	double lon_model;
	double lat_model;
	double alt_model;

	if(Earth_model != "ITRF93"){
		outdata<<"\\begindata"<<std::endl;
		outdata<<"SITES = ('"<<param_obs.obs_name.at(i)<<
			"')"<<std::endl;
		outdata<<param_obs.obs_name.at(i)<<"_CENTER = 399"<<
			std::endl;
		outdata<<param_obs.obs_name.at(i)<<"_FRAME = '"<<
			Earth_model<<"'"<<std::endl;
		outdata<<param_obs.obs_name.at(i)<<"_IDCODE = 100000"<<
			std::to_string(3600+index_obs)<<std::endl;
		outdata<<param_obs.obs_name.at(i)<<"_LATLON = ("<<
			param_obs.lat.at(i)<<","<<
			param_obs.lon.at(i)<<","<<
			param_obs.alt.at(i)<<")"<<std::endl;
		outdata<<param_obs.obs_name.at(i)<<"_UP = 'Z'"<<std::endl;
		outdata<<param_obs.obs_name.at(i)<<"_NORTH = 'X'"<<std::endl;
		outdata<<"\\begintext"<<std::endl;
		outdata.close();
		
		lon_model  = param_obs.lon.at(i);
		lat_model  = param_obs.lat.at(i);
		alt_model  = param_obs.alt.at(i);
	}
	else{
		
		const double lat = param_obs.lat.at(i)/180.*M_PI;
		const double lon = param_obs.lon.at(i)/180.*M_PI;
		const double alt = param_obs.alt.at(i);
		
		//Set WGS84 RF variables
		const double a_WGS=6378.1370;
		const double b_WGS=6356.752314245;
		double N=std::pow(a_WGS,2)/std::sqrt(std::pow(a_WGS,2)*
			pow(std::cos(lat),2)+std::pow(b_WGS,2)*std::pow(std::sin(lat),2));
		
		//Compute WGS84 coordinates
		double x_WGS = (N+alt)*std::cos(lat)*std::cos(lon);
		double y_WGS = (N+alt)*std::cos(lat)*std::sin(lon);
		double z_WGS = (std::pow(b_WGS/a_WGS,2)*N+alt)*std::sin(lat);
		
		//Assume that WGS84 and ITRF2008 are coincident
		double x_ITRF2008 = x_WGS;
		double y_ITRF2008 = y_WGS;
		double z_ITRF2008 = z_WGS;
		
		DACE::AlgebraicVector<double> rr_ITRF2008(3);
		rr_ITRF2008[0]=x_ITRF2008;
		rr_ITRF2008[1]=y_ITRF2008;
		rr_ITRF2008[2]=z_ITRF2008;
		
		//Pass from ITRF2008 to ITRF93 RF
		const double Tx_old   = -24e-6; 
		const double Ty_old   =  2.4e-6;
		const double Tz_old   = -38.6e-6;
		const double D_old    =  3.41e-9;
		const double Rx_old   = -1.71e-3/180.*M_PI;
		const double Ry_old   = -1.48e-3/180.*M_PI;
		const double Rz_old   = -0.30e-3/180.*M_PI;
		const double year_ref =  2000.0;
		
		const double Tx_dot   = -2.80e-6;
		const double Ty_dot   = -0.10e-6;
		const double Tz_dot   = -2.40e-6;
		const double D_dot    =  0.09e-9;
		const double Rx_dot   = -0.11e-3/180.*M_PI;
		const double Ry_dot   = -0.19e-3/180.*M_PI;
		const double Rz_dot   =  0.07e-3/180.*M_PI;
		
		//Convert et to year fraction
		SpiceChar str [TIMLYY];
		timout_c (et, TIMYYY, TIMLYY, str);
		double epoch = std::stod(str);

		double Tx_new = Tx_old + Tx_dot *(epoch-year_ref); 
		double Ty_new = Ty_old + Ty_dot *(epoch-year_ref); 
		double Tz_new = Tz_old + Tz_dot *(epoch-year_ref); 
		double D_new  = D_old  + D_dot  *(epoch-year_ref); 
		double Rx_new = Rx_old + Rx_dot *(epoch-year_ref); 
		double Ry_new = Ry_old + Ry_dot *(epoch-year_ref); 
		double Rz_new = Rz_old + Rz_dot *(epoch-year_ref);
		
		DACE::AlgebraicVector<double> T_new(3);
		T_new[0]=Tx_new; T_new[1]=Ty_new; T_new[2]=Tz_new;

		DACE::AlgebraicMatrix<double> M_new(3,3);
		M_new.at(0,0)=D_new;   M_new.at(0,1)=-Rz_new; M_new.at(0,2)=Ry_new;
		M_new.at(1,0)=Rz_new;  M_new.at(1,1)=D_new;   M_new.at(1,2)=-Rx_new;
		M_new.at(2,0)=-Ry_new; M_new.at(2,1)=Rx_new;  M_new.at(2,2)=D_new;
		
		DACE::AlgebraicVector<double> rr_ITRF93=rr_ITRF2008+T_new+
			M_new*rr_ITRF2008;
			
		double x_ITRF93=rr_ITRF93[0];
		double y_ITRF93=rr_ITRF93[1];
		double z_ITRF93=rr_ITRF93[2];
		
		//Write .sut file
		outdata<<"\\begindata"<<std::endl;
		outdata<<"SITES = ('"<<param_obs.obs_name.at(i)<<
			"')"<<std::endl;
		outdata<<param_obs.obs_name.at(i)<<"_CENTER = 399"<<
			std::endl;
		outdata<<param_obs.obs_name.at(i)<<"_FRAME = '"<<
			Earth_model<<"'"<<std::endl;
		outdata<<param_obs.obs_name.at(i)<<"_IDCODE = 100000"<<
			std::to_string(3600+index_obs)<<std::endl;
		outdata<<param_obs.obs_name.at(i)<<"_XYZ = ("<<
			x_ITRF93<<","<< y_ITRF93<<","<< z_ITRF93<<")"<<std::endl;
		outdata<<param_obs.obs_name.at(i)<<"_UP = 'Z'"<<std::endl;
		outdata<<param_obs.obs_name.at(i)<<"_NORTH = 'X'"<<std::endl;
		outdata<<"\\begintext"<<std::endl;
		outdata.close();
		
		//Compute new coordinates
		SpiceDouble lat_ITRF93, lon_ITRF93, alt_ITRF93;
		SpiceDouble  f;
		SpiceDouble  radii[3];
		SpiceDouble  re, rp;
		SpiceDouble  rectan[3];
		SpiceInt n;
		
		bodvrd_c ( "EARTH", "RADII", 3, &n, radii );

		re  =  radii[0];
		rp  =  radii[2];
		f   =  ( re - rp ) / re;

		rectan[0] =  x_ITRF93;
		rectan[1] =  y_ITRF93;
		rectan[2] =  z_ITRF93;

		recpgr_c ( "EARTH", rectan, re, f, &lon_ITRF93, 
			&lat_ITRF93, &alt_ITRF93);
		
		lon_model  = lon_ITRF93/M_PI*180.;
		lat_model  = lat_ITRF93/M_PI*180.;
		alt_model  = alt_ITRF93;
	}	
	return make_tuple(filename_sut,lon_model,lat_model,alt_model);
}

void Generate_tk(const Observer_param& param_obs, const std::string& Earth_model,
	const double lon, const double lat, const double alt, const unsigned int i,
	const unsigned int index_obs){
	
	
	const unsigned int ID_obs=3600+index_obs;
	
	
	std::ofstream outdata_tf;
	std::string filename_tf = param_obs.obs_path.at(i)+
		param_obs.obs_name.at(i)+".tf";
	outdata_tf.open(filename_tf,std::ios::trunc);
	outdata_tf.precision(17);

	outdata_tf<<"\\begindata"<<std::endl;
	outdata_tf<<"NAIF_BODY_NAME += '"<<param_obs.obs_name.at(i)<<
		"'"<<std::endl;
	outdata_tf<<"NAIF_BODY_CODE += "<<100000<<ID_obs<<std::endl;
	outdata_tf<<"\\begintext"<<std::endl<<std::endl;
	
	outdata_tf<<"\\begindata"<<std::endl;
	outdata_tf<<"FRAME_"<<param_obs.obs_name.at(i)<<"_TOPO          = "<<
		100100<<ID_obs<<std::endl;
	outdata_tf<<"FRAME_"<<100100<<ID_obs<<"_NAME       = '"<<
		param_obs.obs_name.at(i)<<"_TOPO'"<<std::endl;
	outdata_tf<<"FRAME_"<<100100<<ID_obs<<"_CLASS      = "<<4<<std::endl;
	outdata_tf<<"FRAME_"<<100100<<ID_obs<<"_CLASS_ID   = "<<100100<<
		ID_obs<<std::endl;
	outdata_tf<<"FRAME_"<<100100<<ID_obs<<"_CENTER     = "<<100000<<
		ID_obs<<std::endl<<std::endl;
	outdata_tf<<"OBJECT_"<<100000<<ID_obs<<"_FRAME     = '"<<
		param_obs.obs_name.at(i)<<"_TOPO'"<<std::endl<<std::endl;
	outdata_tf<<"TKFRAME_"<<100100<<ID_obs<<"_RELATIVE = '"<<
		Earth_model<<"'"<<std::endl;
	outdata_tf<<"TKFRAME_"<<100100<<ID_obs<<"_SPEC     = '"<<
		"ANGLES"<<"'"<<std::endl;
	outdata_tf<<"TKFRAME_"<<100100<<ID_obs<<"_UNITS    = '"<<
		"DEGREES"<<"'"<<std::endl;
	outdata_tf<<"TKFRAME_"<<100100<<ID_obs<<"_AXES     = ( 3, 2, 3 )"<<
		std::endl;
	outdata_tf<<"TKFRAME_"<<100100<<ID_obs<<"_ANGLES   = ("<<
		-lon<<","<<lat-90.<<","<<180.0<<")"<<std::endl;
	outdata_tf<<"\\begintext"<<std::endl<<std::endl;
	outdata_tf.close();	
			
		
}
	
	
	
	
	
template<typename T>
void Write_TDM(const std::string filename_meas,
	const Observations<T>& Observation,
	const Observer_param& param_obs){
	
	std::string sens_type=param_obs.type;
    std::string part_1 = param_obs.obs_name.at(0);
	std::string part_3;
	if(sens_type=="radar"){
		if(param_obs.sensor_setup.config=="bistatic"){
			part_3 = param_obs.obs_name.at(1);
		}
	}
    std::string obj_name = "UNKNOWN";
    
	//Read measurements
	std::vector<double> et_vect=Observation.et_vect;
	std::vector<T> ra_vect,decl_vect;
	std::vector<T> az_vect,el_vect,doppler_vect,range_vect;
	
	if(sens_type=="optical"){
		if(param_obs.sensor_setup.ra.availability==true){
			ra_vect=Observation.ra_vect;
		}
		if(param_obs.sensor_setup.decl.availability==true){
			decl_vect=Observation.decl_vect;
		}
	}
	else if(sens_type=="radar"){
		if(param_obs.sensor_setup.az.availability==true){
			az_vect=Observation.az_vect;
		}
		if(param_obs.sensor_setup.el.availability==true){
			el_vect=Observation.el_vect;
		}
		if(param_obs.sensor_setup.doppler.availability==true){
			doppler_vect=Observation.doppler_vect;
		}
		if(param_obs.sensor_setup.range.availability==true){
			range_vect=Observation.range_vect;
		}	
	}
	
    const unsigned int N_inst = et_vect.size();
    std::string InTime, OutTime;
	SpiceChar utcstr[TIMLEN];
	timout_c(et_vect.at(0),TIMFMT_OPM,TIMLEN,utcstr);
    InTime = utcstr;
	timout_c(et_vect.at(N_inst-1),TIMFMT_OPM,TIMLEN,utcstr);
    OutTime = utcstr;
    
    //WRITE TDM
    std::ofstream TDMFile;
    TDMFile.open(filename_meas);
    TDMFile.precision(16);
	
    TDMFile << "META_START"<<std::endl;
    TDMFile << "TIME_SYSTEM = UTC"<<std::endl;
    TDMFile << "START_TIME = " << InTime << std::endl;
    TDMFile << "STOP_TIME  = " << OutTime << std::endl;
    TDMFile << "PARTICIPANT_1 = " << part_1 << std::endl;
    TDMFile << "PARTICIPANT_2 = " << obj_name << std::endl;
	if(sens_type=="radar"){
		if(param_obs.sensor_setup.config=="bistatic"){
			TDMFile << "PARTICIPANT_3 = " << part_3 << std::endl;
		}
	}
	
    TDMFile << "MODE = SEQUENTIAL"<< std::endl;
	if(sens_type=="optical"){
		TDMFile << "PATH = 2,1"<< std::endl;
	}
	else if(sens_type=="radar"){
		if(param_obs.sensor_setup.config=="monostatic"){
			TDMFile << "PATH = 1,2,1"<< std::endl;
		}
		else{
			TDMFile << "PATH = 3,2,1"<< std::endl;
		}
	}
	if(sens_type=="optical"){
		TDMFile << "ANGLE_TYPE = RADEC"<< std::endl;
	}
	else if(sens_type=="radar"){
		TDMFile << "ANGLE_TYPE = AZEL"<< std::endl;
	}
    TDMFile << "REFERENCE_FRAME = ECIJ2000"<< std::endl;
	if(sens_type=="optical"){
		TDMFile << "MAG_TYPE = APPARENT"<< std::endl;
	}
	if(sens_type=="radar"){
		if(param_obs.sensor_setup.range.availability==true){
			TDMFile << "RANGE_UNITS = km"<< std::endl;
		}
	}
    TDMFile << "META_STOP"<< std::endl<< std::endl;
    TDMFile << "DATA_START"<< std::endl;
    for(unsigned int i=0; i<N_inst; i++) {
		timout_c(et_vect.at(i),TIMFMT_OPM,TIMLEN,utcstr);
		if(sens_type=="optical"){
			if(param_obs.sensor_setup.ra.availability==true){
				TDMFile << param_obs.sensor_setup.ra.TDM_name <<
					" = " << utcstr <<"\t" << 
					DACE::cons(ra_vect.at(i))/M_PI*180.<< std::endl;
			}
			if(param_obs.sensor_setup.decl.availability==true){
				TDMFile << param_obs.sensor_setup.decl.TDM_name <<
					" = " << utcstr <<"\t" << 
					DACE::cons(decl_vect.at(i))/M_PI*180.<< std::endl;
			}
		}
		else if(sens_type=="radar"){
			if(param_obs.sensor_setup.az.availability==true){
				TDMFile << param_obs.sensor_setup.az.TDM_name <<
					" = " << utcstr <<"\t" << 
					DACE::cons(az_vect.at(i))/M_PI*180.<< std::endl;
			}
			if(param_obs.sensor_setup.el.availability==true){
				TDMFile << param_obs.sensor_setup.el.TDM_name <<
					" = " << utcstr <<"\t" << 
					DACE::cons(el_vect.at(i))/M_PI*180.<< std::endl;
			}
			if(param_obs.sensor_setup.doppler.availability==true){
				TDMFile << param_obs.sensor_setup.doppler.TDM_name <<
					" = " << utcstr <<"\t" << 
					DACE::cons(doppler_vect.at(i))<< std::endl;
			}
			if(param_obs.sensor_setup.range.availability==true){
				TDMFile << param_obs.sensor_setup.range.TDM_name <<
					"   = " << utcstr <<"\t" << 
					DACE::cons(range_vect.at(i))<< std::endl;
			}	
		}
    }
    TDMFile << "DATA_STOP"<< std::endl;
    TDMFile.close();
	
}

void Write_OPM(const double et_est, const DACE::AlgebraicVector<double>& state_est,
	const DACE::AlgebraicMatrix<double>& Cov_est, const DALS_param& param_DALS,
	const double Bfactor_est,const double SRPC_est,
	const DACE::AlgebraicMatrix<double>& Cov_param, const bool flag_dec){
	
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time (&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer,80,"%Y-%m-%d",timeinfo);
	std::string str(buffer);

	time_t curr_time;
	curr_time = time(NULL);

	tm *tm_local = localtime(&curr_time);
	std::stringstream time_ss;
	time_ss<<"T"<<tm_local->tm_hour<<":"<<tm_local->tm_min<<":"<<tm_local->tm_sec;
	std::string current_date=str+time_ss.str();
	
	
	std::ofstream OPMFile;
	std::string filename_opm=param_DALS.output.filename_OPM;
    OPMFile.open(filename_opm);
    OPMFile.precision(20);
		
	//Incipit
	OPMFile<<"CCSDS_OPM_VERS    = 2.0"<<std::endl;
	OPMFile<<"CREATION_DATE     = "<<current_date<<std::endl;
	OPMFile<<"ORIGINATOR        = DEIMOS SPACE"<<std::endl<<std::endl;
	
	OPMFile<<"COMMENT             GEOCENTRIC, CARTESIAN, INERTIAL"<<std::endl;
	OPMFile<<"OBJECT_NAME       = "<<param_DALS.output.obj_name<<std::endl;
	OPMFile<<"OBJECT_ID         = "<<param_DALS.output.obj_ID<<std::endl;
	OPMFile<<"CENTER_NAME       = EARTH"<<std::endl;
	OPMFile<<"REF_FRAME         = ECIJ2000"<<std::endl;
	OPMFile<<"TIME_SYSTEM       = UTC"<<std::endl<<std::endl;
	
	//Convert epoch to UTC
	SpiceChar output[TIMLEN];
	timout_c(et_est,TIMFMT_OPM,TIMLEN,output);
	OPMFile<<"COMMENT  State vector"<<std::endl;
	OPMFile<<"EPOCH             = "<<output<<std::endl;
	OPMFile<<"X                 = "<<state_est[0]<<" [km]"<<std::endl;
	OPMFile<<"Y                 = "<<state_est[1]<<" [km]"<<std::endl;
	OPMFile<<"Z                 = "<<state_est[2]<<" [km]"<<std::endl;
	OPMFile<<"X_DOT             = "<<state_est[3]<<" [km/s]"<<std::endl;
	OPMFile<<"Y_DOT             = "<<state_est[4]<<" [km/s]"<<std::endl;
	OPMFile<<"Z_DOT             = "<<state_est[5]<<" [km/s]"<<std::endl<<std::endl;
	
	if(param_DALS.output.flag_KeplPar==true){
		
		//Convert state vector to Keplerian elements
		DACE::AlgebraicVector<double> KeplPar_est=cart2kep(state_est);
			
		OPMFile<<"COMMENT  Keplerian elements"<<std::endl;
		OPMFile<<"SEMI_MAJOR_AXIS   = "<<KeplPar_est[0]<<" [km]"<<std::endl;
		OPMFile<<"ECCENTRICITY      = "<<KeplPar_est[1]<<std::endl;
		OPMFile<<"INCLINATION       = "<<KeplPar_est[2]/M_PI*180.<<" [deg]"<<std::endl;
		OPMFile<<"RA_OF_ASC_NODE    = "<<KeplPar_est[3]/M_PI*180.<<" [deg]"<<std::endl;
		OPMFile<<"ARG_OF_PERICENTER = "<<KeplPar_est[4]/M_PI*180.<<" [deg]"<<std::endl;
		OPMFile<<"TRUE_ANOMALY      = "<<KeplPar_est[5]/M_PI*180.<<" [deg]"<<std::endl;
		OPMFile<<"GM                = "<<398600.4415<<" [km**3/s**2]"<<std::endl<<std::endl;	
	}
	
	if(param_DALS.output.flag_Cov==true){
		OPMFile<<"COMMENT  State covariance matrix"<<std::endl;
		OPMFile<<"CX_X              = "<<Cov_est.at(0,0)<<" [km**2]"<<std::endl;
		OPMFile<<"CY_X              = "<<Cov_est.at(1,0)<<" [km**2]"<<std::endl;
		OPMFile<<"CY_Y              = "<<Cov_est.at(1,1)<<" [km**2]"<<std::endl;
		OPMFile<<"CZ_X              = "<<Cov_est.at(2,0)<<" [km**2]"<<std::endl;
		OPMFile<<"CZ_Y              = "<<Cov_est.at(2,1)<<" [km**2]"<<std::endl;
		OPMFile<<"CZ_Z              = "<<Cov_est.at(2,2)<<" [km**2]"<<std::endl;
		OPMFile<<"CX_DOT_X          = "<<Cov_est.at(3,0)<<" [km**2/s]"<<std::endl;
		OPMFile<<"CX_DOT_Y          = "<<Cov_est.at(3,1)<<" [km**2/s]"<<std::endl;
		OPMFile<<"CX_DOT_Z          = "<<Cov_est.at(3,2)<<" [km**2/s]"<<std::endl;
		OPMFile<<"CX_DOT_X_DOT      = "<<Cov_est.at(3,3)<<" [km**2/s**2]"<<std::endl;
		OPMFile<<"CY_DOT_X          = "<<Cov_est.at(4,0)<<" [km**2/s]"<<std::endl;
		OPMFile<<"CY_DOT_Y          = "<<Cov_est.at(4,1)<<" [km**2/s]"<<std::endl;
		OPMFile<<"CY_DOT_Z          = "<<Cov_est.at(4,2)<<" [km**2/s]"<<std::endl;
		OPMFile<<"CY_DOT_X_DOT      = "<<Cov_est.at(4,3)<<" [km**2/s**2]"<<std::endl;
		OPMFile<<"CY_DOT_Y_DOT      = "<<Cov_est.at(4,4)<<" [km**2/s**2]"<<std::endl;
		OPMFile<<"CZ_DOT_X          = "<<Cov_est.at(5,0)<<" [km**2/s]"<<std::endl;
		OPMFile<<"CZ_DOT_Y          = "<<Cov_est.at(5,1)<<" [km**2/s]"<<std::endl;
		OPMFile<<"CZ_DOT_Z          = "<<Cov_est.at(5,2)<<" [km**2/s]"<<std::endl;
		OPMFile<<"CZ_DOT_X_DOT      = "<<Cov_est.at(5,3)<<" [km**2/s**2]"<<std::endl;
		OPMFile<<"CZ_DOT_Y_DOT      = "<<Cov_est.at(5,4)<<" [km**2/s**2]"<<std::endl;
		OPMFile<<"CZ_DOT_Z_DOT      = "<<Cov_est.at(5,5)<<" [km**2/s**2]"<<std::endl<<std::endl;
	}
	
	//User defined parameters
	unsigned int index_SRPC=0;
	if(param_DALS.flag_Bfactor_est==true || 
		param_DALS.flag_SRPC_est==true){
		OPMFile<<"COMMENT  Estimated parameters"<<std::endl;
	}
	if(param_DALS.flag_Bfactor_est==true){
		OPMFile<<"USER_DEFINED_BFACTOR  = "<<Bfactor_est<<" [m**2/kg]"<<std::endl;
		
		std::string sigma_Bfactor_str;
		if(flag_dec==false){
			double sigma_Bfactor=std::sqrt(Cov_param.at(0,0));
			sigma_Bfactor_str=std::to_string(sigma_Bfactor)+" [m**2/kg]";
		}
		else{
			sigma_Bfactor_str="NOT AVAILABLE";
		}
		OPMFile<<"USER_DEFINED_BFACTOR_SIGMA  = "<<sigma_Bfactor_str<<std::endl;
		index_SRPC=1;
	}

	if(param_DALS.flag_SRPC_est==true){
		OPMFile<<"USER_DEFINED_SRPC     = "<<SRPC_est<<" [m**2/kg]"<<std::endl;
		
		std::string sigma_SRPC_str;
		if(flag_dec==false){
			double sigma_SRPC=std::sqrt(Cov_param.at(index_SRPC,index_SRPC));
			sigma_SRPC_str=std::to_string(sigma_SRPC)+" [m**2/kg]";
		}
		else{
			sigma_SRPC_str="NOT AVAILABLE";
		}
		OPMFile<<"USER_DEFINED_SRPC_SIGMA  = "<<sigma_SRPC_str<<std::endl;
	}
	
	OPMFile.close();
	
}


void Write_res(const std::vector<Observations<double>>& Observation,
	const std::vector<Observer_param>& param_obs, const DACE::DA& Residuals,
	const std::vector<Observations<DACE::DA>>& Obs_est, const DALS_param& param_DALS){
	
	/*********** Write table of residuals **************/
	std::ofstream outdata_res_table;
	std::string filename_res_table=param_DALS.output.filename_res_table;
    outdata_res_table.open(filename_res_table);
	
	//Sort the measurements in chronological order
	const unsigned int n_sens=Observation.size();
	
	std::vector<double> et_vect_rnd;
	std::vector<unsigned int> sens_ind_rnd;
	std::vector<unsigned int> meas_ind_rnd;
	for(unsigned int s=0; s<n_sens;s++){
		for(unsigned l=0;l<Observation.at(s).et_vect.size();l++){
			//Store time instant
			et_vect_rnd.push_back(Observation.at(s).et_vect.at(l));
			
			//Store sensor and meas indexes
			sens_ind_rnd.push_back(s);
			meas_ind_rnd.push_back(l);
		}
	}
	std::vector<double> et_vect(et_vect_rnd.size());
	std::vector<unsigned int> index(et_vect.size());
	std::tie(et_vect,index)=sort_indexes(et_vect_rnd);
	
	//Sort sens and meas indexes
	std::vector<unsigned int> sens_ind(et_vect.size());
	std::vector<unsigned int> meas_ind(et_vect.size());
	for(unsigned int i=0;i<et_vect.size();i++){
		sens_ind.at(i)=sens_ind_rnd.at(index.at(i));
		meas_ind.at(i)=meas_ind_rnd.at(index.at(i));
	}
		
	
	//Number of sensors
	for(unsigned int i=0;i<et_vect.size();i++){

		//Sens and meas index
		const unsigned int s=sens_ind.at(i);
		const unsigned int l=meas_ind.at(i);
		
		//Read time epoch
		double et_single=Observation.at(s).et_vect.at(l);
		
		//Convert to mjd2000 (UTC)
		SpiceChar jd_char[TIMLEN_RES];
		timout_c(et_single, TIMFMT_RES,TIMLEN_RES,jd_char);
		double jd_single=atof(jd_char);
		double mjd2000_single=jd_single-51544.0;
		
		SpiceDouble scaletemp;
		if(param_obs.at(s).type=="optical"){
			if(param_obs.at(s).sensor_setup.ra.availability==true){
				double ra_real=Observation.at(s).ra_vect.at(l)*M_PI/180.;
				double scale=param_obs.at(s).sensor_setup.ra.accuracy;
				convrt_c(scale,"ARCSECONDS","RADIANS", &scaletemp);
				double scale_ra=scaletemp;
				
				double ra_est=DACE::cons(Obs_est.at(s).ra_vect.at(l));
				double Dra = (ra_est-ra_real)/scale_ra;
				
				outdata_res_table<<std::setprecision(30)<<
					mjd2000_single<<" "<<param_obs.at(s).obs_name.at(0)<<
					" "<<param_obs.at(s).sensor_setup.ra.TDM_name<<
					" "<<std::showpos<<
					std::setprecision(16)<<Dra<<" "<<std::showpos<<
					std::setprecision(16)<<ra_real/M_PI*180.0<<" "<<std::showpos<<
					std::setprecision(16)<<scale<<
					" "<<std::noshowpos<<1<<std::endl;
			}
			
			if(param_obs.at(s).sensor_setup.decl.availability==true){
				double decl_real=Observation.at(s).decl_vect.at(l)*M_PI/180.;
				double scale=param_obs.at(s).sensor_setup.decl.accuracy;
				convrt_c(scale,"ARCSECONDS","RADIANS", &scaletemp);
				double scale_decl=scaletemp;
				
				double decl_est=DACE::cons(Obs_est.at(s).decl_vect.at(l));
				double Ddecl = (decl_est-decl_real)/scale_decl;
				
				outdata_res_table<<std::setprecision(30)<<
					mjd2000_single<<" "<<param_obs.at(s).obs_name.at(0)<<
					" "<<param_obs.at(s).sensor_setup.decl.TDM_name
					<<" "<<std::showpos<<
					std::setprecision(16)<<Ddecl<<" "<<std::showpos<<
					std::setprecision(16)<<decl_real/M_PI*180.0<<" "<<std::showpos<<
					std::setprecision(16)<<scale<<
					" "<<std::noshowpos<<1<<std::endl;
			}
		}
		else{
			
			if(param_obs.at(s).sensor_setup.az.availability==true){
				double az_real=Observation.at(s).az_vect.at(l)*M_PI/180.;
				double scale=param_obs.at(s).sensor_setup.az.accuracy;
				convrt_c(scale,"ARCSECONDS","RADIANS", &scaletemp);
				double scale_az=scaletemp;
				
				double az_est=DACE::cons(Obs_est.at(s).az_vect.at(l));
				double Daz = (az_est-az_real)/scale_az;
				
				outdata_res_table<<std::setprecision(30)<<
					mjd2000_single<<" "<<param_obs.at(s).obs_name.at(0)<<
					" "<<param_obs.at(s).sensor_setup.az.TDM_name
					<<" "<<std::showpos<<
					std::setprecision(16)<<Daz<<" "<<std::showpos<<
					std::setprecision(16)<<az_real/M_PI*180.0<<" "<<std::showpos<<
					std::setprecision(16)<<scale<<
					" "<<std::noshowpos<<1<<std::endl;
			}
			
			if(param_obs.at(s).sensor_setup.el.availability==true){
				double el_real=Observation.at(s).el_vect.at(l)*M_PI/180.;
				double scale=param_obs.at(s).sensor_setup.el.accuracy;
				convrt_c(scale,"ARCSECONDS","RADIANS", &scaletemp);
				double scale_el=scaletemp;
				
				double el_est=DACE::cons(Obs_est.at(s).el_vect.at(l));
				double Del = (el_est-el_real)/scale_el;
				
				outdata_res_table<<std::setprecision(30)<<
					mjd2000_single<<" "<<param_obs.at(s).obs_name.at(0)<<
					" "<<param_obs.at(s).sensor_setup.el.TDM_name
					<<" "<<std::showpos<<
					std::setprecision(16)<<Del<<" "<<std::showpos<<
					std::setprecision(16)<<el_real/M_PI*180.0<<" "<<std::showpos<<
					std::setprecision(16)<<scale<<
					" "<<std::noshowpos<<1<<std::endl;
			}
			
			if(param_obs.at(s).sensor_setup.doppler.availability==true){
				double doppler_real=Observation.at(s).doppler_vect.at(l);
				double scale_doppler=param_obs.at(s).sensor_setup.doppler.accuracy;
				
				double doppler_est=DACE::cons(Obs_est.at(s).doppler_vect.at(l));
				double Ddoppler = (doppler_est-doppler_real)/scale_doppler;
				
				outdata_res_table<<std::setprecision(30)<<
					mjd2000_single<<" "<<param_obs.at(s).obs_name.at(0)<<
					" "<<param_obs.at(s).sensor_setup.doppler.TDM_name
					<<" "<<std::showpos<<
					std::setprecision(16)<<Ddoppler<<" "<<std::showpos<<
					std::setprecision(16)<<doppler_real<<" "<<std::showpos<<
					std::setprecision(16)<<scale_doppler<<
					" "<<std::noshowpos<<1<<std::endl;
			}
			
			if(param_obs.at(s).sensor_setup.range.availability==true){
				double range_real=Observation.at(s).range_vect.at(l);
				double scale_range=param_obs.at(s).sensor_setup.range.accuracy;
				
				double range_est=DACE::cons(Obs_est.at(s).range_vect.at(l));
				double Drange = (range_est-range_real)/scale_range;
				
				outdata_res_table<<std::setprecision(30)<<
					mjd2000_single<<" "<<param_obs.at(s).obs_name.at(0)<<
					" "<<param_obs.at(s).sensor_setup.range.TDM_name
					<<" "<<std::showpos<<
					std::setprecision(16)<<Drange<<" "<<std::showpos<<
					std::setprecision(16)<<range_real<<" "<<std::showpos<<
					std::setprecision(16)<<scale_range<<
					" "<<std::noshowpos<<1<<std::endl;
			}	
		}
		
	}
	
	outdata_res_table.close();
	
	/*********** Write residuals map **************/
	std::ofstream outdata_res_map;
	std::string filename_res_map=param_DALS.output.filename_res_map;
    outdata_res_map.open(filename_res_map);
    outdata_res_map.precision(16);
	
	outdata_res_map<<Residuals<<std::endl;
	outdata_res_map.close();
	
	
}

template<typename T>
DACE::AlgebraicVector<T> LTS_corrections_general(const DACE::AlgebraicVector<T> state, 
		const double et,const Dynamics_param<T> param,const std::string& obs_name, 
		const bool flag_LT, const bool flag_S){
	
	using std::sin; using DACE::sin;
	using std::asin; using DACE::asin;
	

	//Reference system
	std::string ref="SSB";
	SpiceDouble lt;
	
	//Output: obj-obs relative position
	DACE::AlgebraicVector<T> rr_obj_obs_corr(3);
	
	//Compute obs-SSB relative position at time epoch et
	SpiceDouble rr_obs_SSB[3];
	spkpos_c (obs_name.c_str(),et,"J2000","NONE",ref.c_str(),
		rr_obs_SSB, &lt);

	
	//Initialize variables
	DACE::AlgebraicVector<T> state_lt(6);
	SpiceDouble rr_Earth_SSB[3];
	DACE::AlgebraicVector<T> rr_obj_SSB(3);			
	SpiceDouble rr_obs_Earth[3];
	
	//Read method
	std::string method=param.method;
	
	if(flag_LT==false){
		
		//Compute obs-Earth relative position at time epoch et
		spkpos_c (obs_name.c_str(),et,"J2000","NONE","399",
				rr_obs_Earth, &lt);
		for(unsigned int i=0;i<3;i++){
			rr_obj_obs_corr[i]=state[i]-rr_obs_Earth[i];
		}
	}
	else{
		
		//Compute relative position Earth-SSB
		spkpos_c("399",et,"J2000","NONE",ref.c_str(),
				rr_Earth_SSB, &lt);

		//Compute relative position obj-SSB
		for(unsigned int i=0;i<3;i++){
			rr_obj_SSB[i]=state[i]+rr_Earth_SSB[i];
		}

		//First estimate for lt_right
		SpiceDouble delta[3];
		for(unsigned int i=0;i<3;i++){
			delta[i]=DACE::cons(rr_obj_SSB[i])-rr_obs_SSB[i];
		}		
		double lt_right=vnorm_c(delta)/clight_c();
		double lt_old=lt_right;
		
		//Initial errors
		//* Relative position obj-Earth
		if(method=="AIDA"){
			state_lt=Propagate_state(et,
				{et-lt_right},state,param).at(0);
		}
		else{
			state_lt=Propagate_state_linear(et,
				{et-lt_right},state,param).at(0);
		}
		
		//* Relative position Earth-SSB
		spkpos_c("399",et-lt_right,"J2000","NONE",ref.c_str(),
				rr_Earth_SSB, &lt);
		
		//* Relative position obj-SSB
		for(unsigned int i=0;i<3;i++){
			rr_obj_SSB[i]=state_lt[i]+rr_Earth_SSB[i];
		}

		//Compute error		
		for(unsigned int i=0;i<3;i++){
			delta[i]=DACE::cons(rr_obj_SSB[i])-rr_obs_SSB[i];
		}
		
		std::cout<<std::scientific<<std::setprecision(10);
		//std::cout<<"First guess:"<<std::endl;
		//std::cout<<"\tCorrection (s): "<<lt_right<<std::endl;
		//std::cout<<"\tError     (km): "<<vnorm_c(delta)-lt_right*clight_c()<<std::endl;
		
		
		//Set parameters for fixed-point iteration
		double err=1e10; const double err_tol=1e-10;
		double step=1e10; const double step_tol=1e-10; 
		const unsigned int max_iter=10; unsigned  int iter=1;

		while(err>err_tol && step>step_tol && iter<max_iter){
			
			//Display iter
			//std::cout<<"Iter: "<<iter<<std::endl;
			
			//Compute function
			//* Relative position obj-Earth
			if(method=="AIDA"){
				state_lt=Propagate_state(et,
					{et-lt_right},state,param).at(0);
			}
			else{
				state_lt=Propagate_state_linear(et,
					{et-lt_right},state,param).at(0);
			}
				
			//* Relative position Earth-SSB
			spkpos_c("399",et-lt_right,"J2000","NONE",ref.c_str(),
					rr_Earth_SSB, &lt);
			
			//* Relative position obj-SSB
			for(unsigned int i=0;i<3;i++){
				rr_obj_SSB[i]=state_lt[i]+rr_Earth_SSB[i];
			}

			//First estimate for lt_right
			SpiceDouble delta[3];
			for(unsigned int i=0;i<3;i++){
				delta[i]=DACE::cons(rr_obj_SSB[i])-rr_obs_SSB[i];
			}
			lt_right=vnorm_c(delta)/clight_c();
			
			//Update position vectors
			//* Relative position obj-Earth
			if(method=="AIDA"){
				state_lt=Propagate_state(et,
					{et-lt_right},state,param).at(0);
			}
			else{
				state_lt=Propagate_state_linear(et,
					{et-lt_right},state,param).at(0);
			}
				
			//* Relative position Earth-SSB
			spkpos_c("399",et-lt_right,"J2000","NONE",ref.c_str(),
					rr_Earth_SSB, &lt);
			
			//* Relative position obj-SSB
			for(unsigned int i=0;i<3;i++){
				rr_obj_SSB[i]=state_lt[i]+rr_Earth_SSB[i];
			}

			//First estimate for lt_right
			for(unsigned int i=0;i<3;i++){
				delta[i]=DACE::cons(rr_obj_SSB[i])-rr_obs_SSB[i];
			}
		
			//*Error
			err=vnorm_c(delta)-lt_right*clight_c();
			//std::cout<<"\tCorrection (s): "<<lt_right<<std::endl;
			//std::cout<<"\tError     (km): "<<err<<std::endl;
			
			//*Stepsize
			step=std::abs(lt_right-lt_old);
			
			//New step
			lt_old=lt_right;
			iter=iter+1;
		}
		
		//Print final correction
		//std::cout<<"Result:"<<std::endl;
		//std::cout<<"\tCorrection (s): "<<lt_right<<std::endl;
		
		//Update the values
		//* Relative position obj-Earth
		if(method=="AIDA"){
			state_lt=Propagate_state(et,
				{et-lt_right},state,param).at(0);
		}
		else{
			state_lt=Propagate_state_linear(et,
				{et-lt_right},state,param).at(0);
		}
			
		//* Relative position Earth-SSB
		spkpos_c("399",et-lt_right,"J2000","NONE",ref.c_str(),
				rr_Earth_SSB, &lt);
		
		//* Relative position obj-SSB
		for(unsigned int i=0;i<3;i++){
			rr_obj_SSB[i]=state_lt[i]+rr_Earth_SSB[i];
		}

		//Print final error	
		DACE::AlgebraicVector<T> rr_obj_obs_lt(3);
		for(unsigned int i=0;i<3;i++){
			rr_obj_obs_lt[i]=rr_obj_SSB[i]-rr_obs_SSB[i];
		}
		//std::cout<<"\tError     (km): "<<(DACE::cons(rr_obj_obs_lt)).vnorm()-
		//		lt_right*clight_c()<<std::endl;
		
		
		if(flag_S==true){
			
			//Velocity of the obs wrt to SSB
			SpiceDouble state_obs_SSB[6];
			DACE::AlgebraicVector<T> vv_obs_SSB(3); //(it's double, anyway)
			spkezr_c (obs_name.c_str(),et,"J2000","NONE","SSB",
					state_obs_SSB, &lt);
			for(unsigned int i=0;i<3;i++){
				vv_obs_SSB[i]=state_obs_SSB[i+3];
			}
			//Separation angle
			T w=vsep_DACE(rr_obj_obs_lt,vv_obs_SSB);
			
			//Stellar aberration angle
			T phi=asin(vnorm(vv_obs_SSB)*sin(w)/clight_c());
			
			//Angular momentum
			DACE::AlgebraicVector<T> hh=cross(rr_obj_obs_lt,vv_obs_SSB);
			
			//Rotate around hh of phi
			DACE::AlgebraicVector<T> rr_obj_obs_lts=
					vrotv_DACE(rr_obj_obs_lt,hh,phi);
			
			//Return the LT+S corrected position vector
			for(unsigned int i=0;i<3;i++){
				rr_obj_obs_corr[i]=rr_obj_obs_lts[i];
			}			
		}
		else{
			//Return the LT corrected position vector
			for(unsigned int i=0;i<3;i++){
				rr_obj_obs_corr[i]=rr_obj_obs_lt[i];
			}
		}
	}
	return rr_obj_obs_corr;
}


template<typename T>
std::vector<Observations<T>> Generate_measurements(const double et_start, 
	const std::vector<double>& et_vect, 
	const std::vector<unsigned int>& sensor_index,
	const DACE::AlgebraicVector<T>& state_start,
	const std::vector<Observer_param>& param_obs, const Dynamics_param<T> param,
	const std::vector<bool>& flag_LTS_ext, const std::string& ECEF_RF="IAU_EARTH", 
	const bool flag_noise=false, const bool flag_file=false, 
	const std::vector<std::string>& filename={""}){
	
	using std::atan2; using DACE::atan2;
	using std::asin; using DACE::asin;
	
	//Count number of sensors
const clock_t begin_time = clock();
	std::vector<unsigned int> sens_ind_red;
	for(unsigned int i=1;i<sensor_index.size();i++){
		bool flag_ex=false;
		for(unsigned int j=0;j<sens_ind_red.size();j++){
			if(sens_ind_red.at(j)==sensor_index.at(i)){
				flag_ex=true;
				break;
			}
		}
		if(flag_ex==false){
			sens_ind_red.push_back(sensor_index.at(i));
		}
	}
	unsigned int n_sens=sens_ind_red.size();
	//	std::cout << " It takes " <<  (float)((clock() - begin_time))/CLOCKS_PER_SEC << " second(s) for sensor allocation " <<std::endl;
  
	//Initialize output
	std::vector<Observations<T>> meas(n_sens);
	const clock_t begin_time1 = clock();
	//Propagate initial state
	std::vector<DACE::AlgebraicVector<T>> state_obj_vec=
			Propagate_state(et_start,et_vect,state_start,param);
//					std::cout << " It takes " <<  (float)((clock() - begin_time1))/CLOCKS_PER_SEC << " second(s) for propagation " <<std::endl;

	const unsigned int N_obs=et_vect.size();
	const clock_t begin_time2 = clock();

	//*Iterate on observation instants
    for(unsigned int i=0;i<N_obs;i++){

    	//Read and store time epoch
    	double et_single=et_vect[i];
    	meas.at(sensor_index.at(i)).et_vect.push_back(et_single);
		
		//*Compue object position
		DACE::AlgebraicVector<T> state_obj_single=state_obj_vec[i];
		

		//Check whether observer is optical or radar
		std::string type=param_obs.at(sensor_index.at(i)).type;
		

		if(type=="optical"){
			
			//Read observer name
			std::string obs_name=param_obs.at(sensor_index.at(i)).obs_name.at(0);
			
			//*Compute the range
			DACE::AlgebraicVector<T> rho(3);
		
			//Decide whether to consider LTS correction or not
			if(flag_LTS_ext.at(i)==true){
				const clock_t begin_time3 = clock();

				
				//Compute LTS correction
				rho=LTS_corrections_general(state_obj_single,et_single,
						param,obs_name, true, true);
                               // std::cout << " It takes " <<  (float)((clock() - begin_time3))/CLOCKS_PER_SEC << " second(s) for LTS correction" <<std::endl;

			}
			else{
				//*Compute observer position
				SpiceDouble rr_obs[3];
				SpiceDouble lt;
				spkpos_c(obs_name.c_str(),et_single,"J2000",
						"NONE","EARTH", rr_obs, &lt );
				for(unsigned int k=0; k<3; k++){
					rho[k] = state_obj_single[k] - rr_obs[k]; 
				}
				
			}
		
			//*Compute the right ascension
			T ra_single = atan2(rho[1],rho[0]);
			
			if (DACE::cons(ra_single)<0) {
				ra_single = ra_single + 2 * M_PI;
			}
			//*Compute the declination
			T decl_single = asin(rho[2]/DACE::vnorm(rho));
			
			//*Check availability
			bool flag_ra_av   = param_obs.at(sensor_index.at(i)).
				sensor_setup.ra.availability;
			bool flag_decl_av = param_obs.at(sensor_index.at(i)).
				sensor_setup.decl.availability;
			
			//*Add NRnoise
			if(flag_noise==true){
								const clock_t begin_time4 = clock();

				std::default_random_engine generator;
				double std;
				double std_arc_ra=param_obs.at(sensor_index.at(i)).
					sensor_setup.ra.accuracy;
				convrt_c(std_arc_ra,"ARCSECONDS","RADIANS", &std);
				std::normal_distribution<double> noise_ra_gen(0.0, std);
				double noise_ra=noise_ra_gen(generator);
				ra_single=ra_single+noise_ra;
				
				double std_arc_decl=param_obs.at(sensor_index.at(i)).
					sensor_setup.decl.accuracy;
				convrt_c(std_arc_decl,"ARCSECONDS","RADIANS", &std);
				std::normal_distribution<double> noise_decl_gen(0.0, std);
				double noise_decl=noise_decl_gen(generator);
				decl_single=decl_single+noise_decl;

    //                            std::cout << " It takes " <<  (float)((clock() - begin_time4))/CLOCKS_PER_SEC << " second(s) for noise gen" <<std::endl;


			}
			
			//*Store data
			if(flag_ra_av==true){
				meas.at(sensor_index.at(i)).ra_vect.push_back(ra_single);
			}
			if(flag_decl_av==true){
				meas.at(sensor_index.at(i)).decl_vect.push_back(decl_single);
			}
		}
		else if(type=="radar"){
			
			std::string RX_name=param_obs.at(sensor_index.at(i)).obs_name.at(0);
			std::string TX_name;
			if(param_obs.at(sensor_index.at(i)).sensor_setup.config=="bistatic"){
				TX_name=param_obs.at(sensor_index.at(i)).obs_name.at(1);
			}
			
			//Compute azimuth and elevation at the receiver
			//*Compute observer position
			DACE::AlgebraicVector<T> rho_obj_RX(3);
			SpiceDouble rr_RX[3];
			SpiceDouble lt;
			spkpos_c(RX_name.c_str(),et_single,"J2000",
					"NONE","EARTH", rr_RX, &lt );
			for(unsigned int k=0; k<3; k++){
				rho_obj_RX[k] = state_obj_single[k] - rr_RX[k]; 
			}
			
			//Compute rotation matrix
			std::string RX_RF_name=param_obs.at(sensor_index.at(i)).obsRF_name.at(0);
			SpiceDouble Aeci2RX_spice[3][3];
			pxform_c("J2000",RX_RF_name.c_str(),et_single,Aeci2RX_spice);
			
			DACE::AlgebraicMatrix<double> Aeci2RX(3,3);
			for(unsigned int j=0;j<3;j++){
				for(unsigned int k=0;k<3;k++){
					Aeci2RX.at(j,k)=Aeci2RX_spice[j][k];
				}
			}

			DACE::AlgebraicVector<T> rho_obj_RX_RXtopo=
				Aeci2RX*rho_obj_RX;
				
			//*Compute the azimuth
			T az_single = -atan2(rho_obj_RX_RXtopo[1],
				rho_obj_RX_RXtopo[0]);
			if (DACE::cons(az_single)<0) {
				az_single = az_single + 2 * M_PI;
			}
			//*Compute the elevation
			T el_single = asin(rho_obj_RX_RXtopo[2]/
				DACE::vnorm(rho_obj_RX_RXtopo));
				
			T rhoRX_single=DACE::vnorm(rho_obj_RX);
			T range_single;
			if(param_obs.at(sensor_index.at(i)).sensor_setup.config=="bistatic"){
			
				DACE::AlgebraicVector<T> rho_obj_TX(3);
				SpiceDouble rr_TX[3];
				spkpos_c(TX_name.c_str(),et_single,"J2000",
						"NONE","EARTH", rr_TX, &lt );
				for(unsigned int k=0; k<3; k++){
					rho_obj_TX[k] = state_obj_single[k] - rr_TX[k]; 
				}
				
				T rhoTX_single=DACE::vnorm(rho_obj_TX);
				range_single=rhoRX_single+rhoTX_single;
			}
			else{
				range_single=2.*rhoRX_single;
			}
						
			//Compute the Doppler shift
			
			//*RX ECEF position
			SpiceDouble rr_RX_ECEF[3];
			spkpos_c(RX_name.c_str(),et_single,ECEF_RF.c_str(),
				"NONE","EARTH", rr_RX_ECEF, &lt );
			
			//Object state in ECEF
			SpiceDouble Aeci2ECEF_state_spice[6][6];
			sxform_c("J2000",ECEF_RF.c_str(),et_single,Aeci2ECEF_state_spice);
			DACE::AlgebraicMatrix<double> Aeci2ECEF_state(6,6);
			for(unsigned int j=0;j<6;j++){
				for(unsigned int k=0;k<6;k++){
					Aeci2ECEF_state.at(j,k)=Aeci2ECEF_state_spice[j][k];
				}
			}
			
			DACE::AlgebraicVector<T> state_obj_ECEF=
				Aeci2ECEF_state*state_obj_single;
				
			DACE::AlgebraicVector<T> rr_obj_ECEF(3),vv_obj_ECEF(3);
			for(unsigned int j=0;j<3;j++){
				rr_obj_ECEF[j]=state_obj_ECEF[j];
				vv_obj_ECEF[j]=state_obj_ECEF[j+3];
			}
			DACE::AlgebraicVector<T> rho_RX_obj(3);
			for(unsigned int j=0;j<3;j++){
				rho_RX_obj[j]=rr_RX_ECEF[j]-rr_obj_ECEF[j];
			}
			rho_RX_obj=DACE::normalize(rho_RX_obj);
			
			DACE::AlgebraicVector<T> rho_TX_obj(3);
			if(param_obs.at(sensor_index.at(i)).sensor_setup.config=="bistatic"){
				//*TX ECEF position
				SpiceDouble rr_TX_ECEF[3];
				spkpos_c(TX_name.c_str(),et_single,ECEF_RF.c_str(),
					"NONE","EARTH", rr_TX_ECEF, &lt );
				
				for(unsigned int j=0;j<3;j++){
					rho_TX_obj[j]=rr_TX_ECEF[j]-rr_obj_ECEF[j];
				}
				rho_TX_obj=DACE::normalize(rho_TX_obj);
			}
			else{
				rho_TX_obj=rho_RX_obj;
			}
			
			const double lambda=clight_c()/param_obs.at(sensor_index.at(i)).
				sensor_setup.frequency;
			T doppler_single=1./lambda*(DACE::dot(vv_obj_ECEF,rho_RX_obj)+
                DACE::dot(vv_obj_ECEF,rho_TX_obj));
			
			//*Check availability
			bool flag_az_av      = param_obs.at(sensor_index.at(i)).
				sensor_setup.az.availability;
			bool flag_el_av      = param_obs.at(sensor_index.at(i)).
				sensor_setup.el.availability;
			bool flag_doppler_av = param_obs.at(sensor_index.at(i)).
				sensor_setup.doppler.availability;
			bool flag_range_av   = param_obs.at(sensor_index.at(i)).
				sensor_setup.range.availability;
			
			//*Add NRnoise
			if(flag_noise==true){
				
				std::default_random_engine generator;
				double std;
				double std_arc_az=param_obs.at(sensor_index.at(i)).
					sensor_setup.az.accuracy;
				convrt_c(std_arc_az,"ARCSECONDS","RADIANS", &std);
				std::normal_distribution<double> noise_az_gen(0.0, std);
				double noise_az=noise_az_gen(generator);
				az_single=az_single+noise_az;
				
				double std_arc_el=param_obs.at(sensor_index.at(i)).
					sensor_setup.el.accuracy;
				convrt_c(std_arc_el,"ARCSECONDS","RADIANS", &std);
				std::normal_distribution<double> noise_el_gen(0.0, std);
				double noise_el=noise_el_gen(generator);
				el_single=el_single+noise_el;
				
				std=param_obs.at(sensor_index.at(i)).
					sensor_setup.doppler.accuracy;
				std::normal_distribution<double> noise_doppler_gen(0.0, std);
				double noise_doppler=noise_doppler_gen(generator);
				doppler_single=doppler_single+noise_doppler;
				
				std=param_obs.at(sensor_index.at(i)).
					sensor_setup.range.accuracy;
				std::normal_distribution<double> noise_range_gen(0.0, std);
				double noise_range=noise_range_gen(generator);
				range_single=range_single+noise_range;

			}
			
			//*Store data
			if(flag_az_av==true){
				meas.at(sensor_index.at(i)).az_vect.push_back(az_single);
			}
			if(flag_el_av==true){
				meas.at(sensor_index.at(i)).el_vect.push_back(el_single);
			}
			if(flag_doppler_av==true){
				meas.at(sensor_index.at(i)).doppler_vect.push_back(doppler_single);
			}
			if(flag_range_av==true){
				meas.at(sensor_index.at(i)).range_vect.push_back(range_single);
			}			
		}

    }
	//                                std::cout << " It takes " <<  (float)((clock() - begin_time2))/CLOCKS_PER_SEC << " second(s) for iterations"<<std::endl;

									const clock_t begin_time5 = clock();

	//Write TDM file for all the considered sensors
	if(flag_file==true){
		const unsigned int n_TDM=meas.size();
		for(unsigned int i=0;i<n_TDM;i++){
			Write_TDM(filename.at(i),meas.at(i),param_obs.at(i));
		}
	}
	//                                std::cout << " It takes " <<  (float)((clock() - begin_time5))/CLOCKS_PER_SEC << " second(s) for TDM creation" <<std::endl;

	return meas;
	
}


template<typename T>
struct OD_data{
	double et_est;
	DACE::AlgebraicVector<T> state_est;
	DACE::AlgebraicVector<double> state_apr;
	DACE::AlgebraicMatrix<double> Cov_apr;
	std::string flag_apr;
	
};

//Function to calculate residuals
template<typename T>
std::tuple<T,std::vector<Observations<T>>> ResidualCalc(const OD_data<T> Data_est, 
		const std::vector<Observations<double>>& Observation,
		const std::vector<Observer_param>& param_obs,
		const Dynamics_param<T>& param,const std::vector<bool>& flag_LTS, 
		const std::string& ECEF_RF){
	
	using std::atan; using DACE::atan;
	

	//Extract data from Data_est
	double et_start=Data_est.et_est;
	DACE::AlgebraicVector<T> state_start=Data_est.state_est;
	DACE::AlgebraicVector<double> state_apr=Data_est.state_apr;
	DACE::AlgebraicMatrix<double> Cov_apr=Data_est.Cov_apr;

	
	//Extract time vector and build sensor index
	const unsigned int n_sens=Observation.size();
	std::vector<double> et_single_vect_rnd;
	std::vector<unsigned int> sensor_index_rnd;
	std::vector<unsigned int> meas_index_rnd;
	std::vector<bool> flag_LTS_rnd;

        const clock_t begin_time = clock();
        
	for(unsigned int i=0;i<n_sens;i++){	
		//Scan all time instants
		for(unsigned int j=0;j<Observation.at(i).et_vect.size();j++){
			//Store time instant
			et_single_vect_rnd.push_back(Observation.at(i).et_vect.at(j));
			
			//Store sensor and meas indexes
			sensor_index_rnd.push_back(i);
			meas_index_rnd.push_back(j);
			
			//Store LTS flag
			flag_LTS_rnd.push_back(flag_LTS.at(i));
		}
	}
	//std::cout << " It takes " <<  (float)((clock() - begin_time))/CLOCKS_PER_SEC/60.0 << " minute(s) to store observables " <<std::endl;
    
	//Sort the time instant in ascending order
	std::vector<double> et_single_vect(et_single_vect_rnd.size());
	std::vector<unsigned int> index(et_single_vect.size());
	std::tie(et_single_vect,index)=sort_indexes(et_single_vect_rnd);
	
	//Adjust sensor and meas indexes
	std::vector<unsigned int> sensor_index(sensor_index_rnd.size());
	std::vector<unsigned int> meas_index(et_single_vect.size());
	std::vector<bool> flag_LTS_ext(et_single_vect.size());
	for(unsigned int i=0;i<index.size();i++){
		sensor_index.at(i)=sensor_index_rnd.at(index.at(i));
		meas_index.at(i)=meas_index_rnd.at(index.at(i));
		flag_LTS_ext.at(i)=flag_LTS_rnd.at(index.at(i));
	}
	
	//Estimate measurements
const clock_t begin_time2 = clock();
	std::vector<Observations<T>> Obs_est=Generate_measurements(et_start, 
		et_single_vect,sensor_index,state_start, param_obs, param, flag_LTS_ext,
		ECEF_RF);
	//std::cout << " It takes " <<  (float)((clock() - begin_time2))/CLOCKS_PER_SEC/60.0 << " minute(s) to generate measurememnts" <<std::endl;
    

	//Compute measurements residuals
	T Residuals=0.0;
const clock_t begin_time3 = clock();

	for(unsigned int i=0;i<et_single_vect.size();i++){	
	
		//Sens and meas index
		const unsigned int s=sensor_index.at(i);
		const unsigned int l=meas_index.at(i);
	
		SpiceDouble scaletemp;
		if(param_obs.at(s).type=="optical"){
			if(param_obs.at(s).sensor_setup.ra.availability==true){
				double ra_real=Observation.at(s).ra_vect.at(l)*M_PI/180.;
                                if (ra_real < 0.0)
                                    ra_real += 2.*M_PI;
				double scale=param_obs.at(s).sensor_setup.ra.accuracy;
				convrt_c(scale,"ARCSECONDS","RADIANS", &scaletemp);
				double scale_ra=scaletemp;
				
				T ra_est=Obs_est.at(s).ra_vect.at(l);
				if (DACE::cons(ra_est) < 0.0)
                                    ra_est += 2.*M_PI;
				
                                T Dra = (ra_est-ra_real)/scale_ra;
                        //        std::cout <<  DACE::cons(ra_est) << " " << DACE::cons(ra_real) << "\t";

				Residuals += Dra*Dra;
			}
			
			if(param_obs.at(s).sensor_setup.decl.availability==true){
				double decl_real=Observation.at(s).decl_vect.at(l)*M_PI/180.;
				double scale=param_obs.at(s).sensor_setup.decl.accuracy;
				convrt_c(scale,"ARCSECONDS","RADIANS", &scaletemp);
				double scale_decl=scaletemp;
                               if (decl_real < 0.0)
                                    decl_real += 2.*M_PI;
				
				T decl_est=Obs_est.at(s).decl_vect.at(l);
                                if (DACE::cons(decl_est) < 0.0)
                                    decl_est += 2.*M_PI;
				
				T Ddecl = (decl_est-decl_real)/scale_decl;
                         //       std::cout << DACE::cons(decl_est) << " " << DACE::cons(decl_real) << std::endl;

				Residuals += Ddecl*Ddecl;
			}
		}
		else{
			
			if(param_obs.at(s).sensor_setup.az.availability==true){
				double az_real=Observation.at(s).az_vect.at(l)*M_PI/180.;
				double scale=param_obs.at(s).sensor_setup.az.accuracy;
				convrt_c(scale,"ARCSECONDS","RADIANS", &scaletemp);
				double scale_az=scaletemp;
				
				T az_est=Obs_est.at(s).az_vect.at(l);
				T Daz = (az_est-az_real)/scale_az;
				
				Residuals += Daz*Daz;
			}
			
			if(param_obs.at(s).sensor_setup.el.availability==true){
				double el_real=Observation.at(s).el_vect.at(l)*M_PI/180.;
				double scale=param_obs.at(s).sensor_setup.el.accuracy;
				convrt_c(scale,"ARCSECONDS","RADIANS", &scaletemp);
				double scale_el=scaletemp;
				
				T el_est=Obs_est.at(s).el_vect.at(l);
				T Del = (el_est-el_real)/scale_el;
				
				Residuals += Del*Del;
			}
			
			if(param_obs.at(s).sensor_setup.doppler.availability==true){
				double doppler_real=Observation.at(s).doppler_vect.at(l);
				double scale_doppler=param_obs.at(s).sensor_setup.doppler.accuracy;
				
				T doppler_est=Obs_est.at(s).doppler_vect.at(l);
				T Ddoppler = (doppler_est-doppler_real)/scale_doppler;
				
				Residuals += Ddoppler*Ddoppler;
			}
			
			if(param_obs.at(s).sensor_setup.range.availability==true){
				double range_real=Observation.at(s).range_vect.at(l);
				double scale_range=param_obs.at(s).sensor_setup.range.accuracy;
				
				T range_est=Obs_est.at(s).range_vect.at(l);
				T Drange = (range_est-range_real)/scale_range;
				
				Residuals += Drange*Drange;
			}	
		}	
	}
	//std::cout << " It takes " <<  (float)((clock() - begin_time3))/CLOCKS_PER_SEC/60.0 << " minute(s) to compute residuals" <<std::endl;

    //Compute contribution due to a priori info, if any
    std::string flag_apr=Data_est.flag_apr;
const clock_t begin_time4 = clock();

    if(flag_apr=="yes"){

    	DACE::AlgebraicVector<T> state_diff=Data_est.state_est-Data_est.state_apr;

    	//Compute scaled covariance matrix
    	DACE::AlgebraicMatrix<double> Cov_apr_scaled(6,6);
    	for(unsigned int i=0;i<6;i++){
    		for(unsigned int j=0;j<6;j++){
    			Cov_apr_scaled.at(i,j)=Cov_apr.at(i,j)*1e6;
    		}
    	}
    	DACE::AlgebraicMatrix<double> Cov_apr_scaled_inv=Cov_apr_scaled.inv();
    	DACE::AlgebraicMatrix<double> Cov_apr_inv(6,6);
    	for(unsigned int i=0;i<6;i++){
    		for(unsigned int j=0;j<6;j++){
    			Cov_apr_inv.at(i,j)=Cov_apr_scaled_inv.at(i,j)*1e6;
    		}
    	}

    	DACE::AlgebraicVector<T> state_rot=Cov_apr_inv*state_diff;
    	T Residuals_apr=state_diff.dot(state_rot);
		
    	//Update Residuals
    	Residuals=Residuals+Residuals_apr;
    }
	//std::cout << " It takes " <<  (float)((clock() - begin_time4))/CLOCKS_PER_SEC/60.0 << " minute(s) for apr flag" <<std::endl;

    return std::make_tuple(Residuals,Obs_est);    
}


std::string ReplaceAll(std::string str, const std::string& from, 
	const std::string& to) {
	size_t start_pos = 0;
	while((start_pos = str.find(from, start_pos)) != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
	return str;
}

#endif