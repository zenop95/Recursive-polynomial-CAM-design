#ifndef Objects_H_
#define Objects_H_

/**********************************************************/
/******************* Target object parameters *************/
/**********************************************************/
struct Object_param {
	
	//Object ID and name
	unsigned int sat_ID;
	std::string sat_name;
	
	//Catalogue filename
	std::string filename_cat;
	
	//Object physical parameters
	double A;
	double Cr;
	double Bfactor;
	double SRPC;
	
};


/**********************************************************/
/***************** State estimates parameters *************/
/**********************************************************/
struct Estimate_param{
	
	//Filename
	std::string filename_state;
	
	//Time gap between first obs and et_data
	double Dt_data;
	
	//Flags for parameters availability
	bool flag_Cov_av;
	bool flag_Bfactor_av;
	bool flag_SRPC_av;
	
	//Search space
	std::string frame;
	bool flag_adaptive;
	double domain_pos;
	double domain_vel;
	int sigma;
	
	//OD guess accuracy
	double sigma_rr;
	double sigma_vv;
	
	//A priori state estimates
	DACE::AlgebraicVector<double> state_apr;
	DACE::AlgebraicMatrix<double> Cov_apr;
	
};

/**********************************************************/
/******************* Measurements parameters **************/
/**********************************************************/

template<typename T>
struct Observations{
	
	//Vector of time instants
	std::vector<double> et_vect;
	
	//Vector of optical measurements
	std::vector<T> ra_vect;
	std::vector<T> decl_vect;
	
	//Vector of radar measurements
	std::vector<T> az_vect;
	std::vector<T> el_vect;
	std::vector<T> doppler_vect;
	std::vector<T> range_vect;
	
	double n_meas_inst;
};

struct Kernels{
	
	std::string path;
	
	std::string ITRF93;
	std::string ITRF93_assoc;
	std::string IAU_assoc;
	std::string pck;
	std::string lsk;
	
};

struct Earth_param{
	
	//Earth model
	std::string Earth_model;
	
	//SPICE kernels
	Kernels kernels;
	
};

template<typename T>
struct Measurements_param{
	
	//Sensors parameters
	std::vector<std::string> filename_meas;
	std::vector<double> Dt;
	std::vector<unsigned int> N_obs_max;
	
	//Measurements parameters
	double Dt_pass;
	unsigned int N_pass_max;
	double Dt_survey;
	
	//Flag for kernel generation (sims only)
	bool flag_gen_kernel;
	
	//Flag for light time (LT) and stellar ab. correction
	std::vector<bool> flag_LTS;
	
	//Set of observations
	std::vector<Observations<T>> Obs;
};

/**********************************************************/
/********************* Observer parameters ****************/
/**********************************************************/

struct Meas_setup{	
	bool availability;
	double accuracy;
	std::string TDM_name;
};

struct Sensor_setup{	
	
	//Optical
	Meas_setup ra;
	Meas_setup decl;

	//Radar
	Meas_setup az;
	Meas_setup el;
	Meas_setup doppler;
	Meas_setup range;
	
	//Radar configuration
	std::string config;
	
	//Radar frequency
	double frequency;
};

struct Observer_param {
	
	std::vector<std::string> obs_name;
	std::vector<std::string> obsRF_name;
	std::vector<std::string> obs_path;
	std::vector<std::string> filename_obs;
	std::vector<double> lat;
	std::vector<double> lon;
	std::vector<double> alt;
	
	//Sensor type
	std::string type;
	
	//Limits
	std::vector<double> el_min;
	double sunEl_lim;
	double moonSep_lim;
	double magn_lim;
	
	//Sensor setup
	Sensor_setup sensor_setup;
	
};


/**********************************************************/
/************************ DALS variables ******************/
/**********************************************************/
struct DALS_output{
	
	std::string filename_OPM;
	std::string filename_res_table;
	std::string filename_res_map;
	std::string filename_DALS;
	
	std::string obj_name;
	double obj_ID;
	bool flag_Cov;
	bool flag_KeplPar;
	
};

struct DALS_param{
	
	//DA order
	unsigned int order;
	
	//Number of variables
	unsigned int nn_var;
	
	//Variables for iterative cycle control
	unsigned int max_iter;
	double fun_tol;
	double step_tol;
	double grad_tol;
	double map_tol;
	std::string stop_crit;
	
	//Flag for Bfactor and SRPC estimation
	bool flag_Bfactor_est;
	bool flag_SRPC_est;
	
	//Estimate epoch
	double et_est;
	
	//Output variables
	DALS_output output;
	
};

/**********************************************************/
/******************** Dynamics variables ******************/
/**********************************************************/
template<typename T>
class Dynamics_param{
		
	public:
	
	//Method
	std::string method;
	
	//Tolerance
	double tol;
	
	//*Parameters
	T Bfactor;
	T SRPC;
	
	unsigned int gravOrd;
	std::string gravmodel;
	DACE::AlgebraicVector<int> AIDA_flags;
	
	//SADA
	unsigned int  SunMoonFlag;
	unsigned int      SRPFlag;
	unsigned int     dragFlag;
	unsigned int tesseralFlag;
	unsigned int shortPeriodicsFlag;
	bool flag_linear;
	
	//*Member functions	
	void update(const Dynamics_param<T>&  param_input, const T Bfactor,
		const T SRPC);
		
	void build(const Dynamics_param<double>&  param_input, const T Bfactor,
		const T SRPC);
	
};

template<typename T>
void Dynamics_param<T>::update(const Dynamics_param<T>&  param_input, 
	const T Bfactor,const T SRPC){
			
	std::string method=param_input.method;
	
	this -> method= method;
	if(method=="AIDA"){
		this -> gravOrd    = param_input.gravOrd;
		this -> gravmodel  = param_input.gravmodel;
		this -> AIDA_flags = param_input.AIDA_flags;		
	}
	else if (method=="SADA"){
		
		this -> SunMoonFlag        = param_input.SunMoonFlag;
		this -> SRPFlag            = param_input.SRPFlag;
		this -> dragFlag           = param_input.dragFlag;
		this -> tesseralFlag       = param_input.tesseralFlag;	
		this -> shortPeriodicsFlag = param_input.shortPeriodicsFlag;
		this -> flag_linear        = param_input.flag_linear;
	}
	
	this -> Bfactor = Bfactor;
	this -> SRPC  = SRPC;
	this -> tol   = param_input.tol;
}

template<typename T>
void Dynamics_param<T>::build(const Dynamics_param<double>&  param_input, 
	const T Bfactor,const T SRPC){
			
	std::string method=param_input.method;
	
	this -> method= method;
	if(method=="AIDA"){
		this -> gravOrd    = param_input.gravOrd;
		this -> gravmodel  = param_input.gravmodel;
		this -> AIDA_flags = param_input.AIDA_flags;		
	}
	else if (method=="SADA"){
		
		this -> SunMoonFlag        = param_input.SunMoonFlag;
		this -> SRPFlag            = param_input.SRPFlag;
		this -> dragFlag           = param_input.dragFlag;
		this -> tesseralFlag       = param_input.tesseralFlag;	
		this -> shortPeriodicsFlag = param_input.shortPeriodicsFlag;
		this -> flag_linear        = param_input.flag_linear;
	}
	
	this -> Bfactor = Bfactor;
	this -> SRPC  = SRPC;
	this -> tol   = param_input.tol;
}


#endif
