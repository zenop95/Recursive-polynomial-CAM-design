#ifndef Exceptions_H_
#define Exceptions_H_

void handle_eptr(std::exception_ptr eptr) // passing by value is ok
{
    try {
        if (eptr) {
            std::rethrow_exception(eptr);
        }
    } catch(const std::exception& e) {
		std::cout << "\n\t\033[31mProgram terminated: " <<std::endl;
        std::cout << "\t    " << e.what() << "\033[0m\n";
    }
}


void input_except(const std::string& json_LS){
	

	//Open json file
	std::ifstream inputFile(json_LS, std::ifstream::binary);
	
	// Read the Json value
    Json::Value root;
    inputFile >> root;
	
	if (inputFile.is_open()==false) {
		throw std::invalid_argument
			("DALS::constructor: invalid json filename" );
	}
	
	/**********************************************************/
	/********************* Observer parameters ****************/
	/**********************************************************/
	if(root["Observer"]["type"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in "
				"json input file\n\t    [\"Observer\"][\"type\"] missing" );
	}
	else{
		if(root["Observer"]["type"]!="optical" && 
			root["Observer"]["type"] != "radar"){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"Invalid [\"Observer\"][\"type\"], "
					"please selected either \"radar\" or \"optical\" " );			
		}
	}
	std:: string type  = root["Observer"]["type"].asString();
	
	if(root["Observer"]["participant_1"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"Observer\"][\"participant_1\"] missing" );
	}
	
	if(root["Observer"]["participant_1"]["name"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"Observer\"][\"participant_1\"]"
				"[\"name\"] missing" );
	}
	if(root["Observer"]["participant_1"]["path"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"Observer\"][\"participant_1\"]"
				"[\"path\"] missing" );
	}		
	if(root["Observer"]["participant_1"]["latitude"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"Observer\"][\"participant_1\"]"
				"[\"latitude\"] missing" );
	}
	/*
	else{
		double lat=root["Observer"]["participant_1"]["latitude"].asDouble();
		if(lat<2.*M_PI){
			std::cout<<"\n\t\033[33mDALS::constructor: warning for "
				"[\"Observer\"][\"participant_1\"][\"latitude\"]"<<std::endl;
			std::cout<<"\t   Is it expressed in degree? (type yes or no): \033[0m";
			std::string answer;
			std::getline(std::cin, answer);
			if(answer=="no"){
				throw std::invalid_argument
					("DALS::constructor: error in json input file\n\t    "
						"[\"Observer\"][\"participant_1\"]"
						"[\"latitude\"] must be expressed in degree" );
			}
		}
		
	}
	*/
	if(root["Observer"]["participant_1"]["longitude"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"Observer\"][\"participant_1\"]"
				"[\"longitude\"] missing" );
	}
	/*
	else{
		double lon=root["Observer"]["participant_1"]["longitude"].asDouble();
		if(lon<2.*M_PI){
			std::cout<<"\n\t\033[33mDALS::constructor: warning for "
				"[\"Observer\"][\"participant_1\"][\"longitude\"]"<<std::endl;
			std::cout<<"\t   Is it expressed in degree? (type yes or no): \033[0m";
			std::string answer;
			std::getline(std::cin, answer);
			if(answer=="no"){
				throw std::invalid_argument
					("DALS::constructor: error in json input file\n\t    "
						"[\"Observer\"][\"participant_1\"]"
						"[\"longitude\"] must be expressed in degree" );
			}
		}
		
	}
	*/
	if(root["Observer"]["participant_1"]["altitude"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"Observer\"][\"participant_1\"]"
				"[\"altitude\"] missing" );
	}
	else{
		double alt=root["Observer"]["participant_1"]["altitude"].asDouble();
		if(alt>9){
			std::cout<<"\n\t\033[33mDALS::constructor: warning for "
				"[\"Observer\"][\"participant_1\"][\"altitude\"]"<<std::endl;
			std::cout<<"\t   Is it expressed in kilometers? (type yes or no): \033[0m";
			std::string answer;
			std::getline(std::cin, answer);
			if(answer=="no"){
				throw std::invalid_argument
					("DALS::constructor: error in json input file\n\t    "
						"[\"Observer\"][\"participant_1\"]"
						"[\"altitude\"] must be expressed in kilometers" );
			}
		}
		
	}
	
	//Sensor setup
	if(root["Observer"]["sensor_setup"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"Observer\"][\"sensor_setup\"] missing" );
	}
	if(type=="optical"){
		
		//Right ascension
		if(root["Observer"]["sensor_setup"]["right_ascension"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"right_ascension\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["right_ascension"]
			["availability"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file \n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"right_ascension\"][\"availability\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["right_ascension"]
			["accuracy"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file \n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"right_ascension\"][\"accuracy\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["right_ascension"]
			["TDM_name"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file \n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"right_ascension\"][\"TDM_name\"] missing" );
		}
		
		//Declination
		if(root["Observer"]["sensor_setup"]["declination"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"declination\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["declination"]
			["availability"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file \n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"right_ascension\"][\"availability\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["declination"]
			["accuracy"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file \n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"right_ascension\"][\"accuracy\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["declination"]
			["TDM_name"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file \n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"right_ascension\"][\"TDM_name\"] missing" );
		}
	}
	else if(type=="radar"){	
		
		//Azimuth
		if(root["Observer"]["sensor_setup"]["azimuth"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"azimuth\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["azimuth"]
			["availability"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"azimuth\"][\"availability\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["azimuth"]
			["accuracy"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"azimuth\"][\"accuracy\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["azimuth"]
			["TDM_name"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"azimuth\"][\"TDM_name\"] missing" );
		}
			
		//Elevation
		if(root["Observer"]["sensor_setup"]["elevation"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"elevation\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["elevation"]
			["availability"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"elevation\"][\"availability\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["elevation"]
			["accuracy"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"elevation\"][\"accuracy\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["elevation"]
			["TDM_name"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"elevation\"][\"TDM_name\"] missing" );
		}
			
		//Doppler
		if(root["Observer"]["sensor_setup"]["doppler"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"doppler\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["doppler"]
			["availability"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"doppler\"][\"availability\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["doppler"]
			["accuracy"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"doppler\"][\"accuracy\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["doppler"]
			["TDM_name"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"doppler\"][\"TDM_name\"] missing" );
		}
			
		//Range
		if(root["Observer"]["sensor_setup"]["range"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"range\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["range"]
			["availability"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"range\"][\"availability\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["range"]
			["accuracy"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"range\"][\"accuracy\"] missing" );
		}
		if(root["Observer"]["sensor_setup"]["range"]
			["TDM_name"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"range\"][\"TDM_name\"] missing" );
		}
			
		//Radar frequency
		if(root["Observer"]["sensor_setup"]["frequency"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"frequency\"] missing" );
		}
		else{
			double freq=root["Observer"]["sensor_setup"]["frequency"].asDouble();
			if(freq<1e6){
				std::cout<<"\n\t\033[33mDALS::constructor: warning for "
					"[\"Observer\"][\"sensor_setup\"][\"frequency\"]"<<std::endl;
				std::cout<<"\t   Is it expressed in Hz? (type yes or no): \033[0m";
				std::string answer;
				std::getline(std::cin, answer);
				if(answer=="no"){
					throw std::invalid_argument
						("DALS::constructor: error in json input file\n\t    "
							"[\"Observer\"][\"sensor_setup\"][\"frequency\"] "
							"must be expressed in Hz" );
				}
			}
			
		}	
		//Radar configuration
		if(root["Observer"]["sensor_setup"]["configuration"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Observer\"][\"sensor_setup\"]"
					"[\"configuration\"] missing" );
		}
		std::string config=root["Observer"]["sensor_setup"]
			["configuration"].asString();
		if(config !="monostatic" && config != "bistatic"){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"Invalid [\"Observer\"][\"sensor_setup\"]"
					"[\"configuration\"], please type either \"monostatic\""
					"or \"bistatic\"" );
		}

		//Read data of second observatory in case of bistatic config
		if(config=="bistatic"){	
		
			if(root["Observer"]["participant_2"].empty()==true){
				throw std::invalid_argument
					("DALS::constructor: error in json input file\n\t    "
						"[\"Observer\"][\"participant_2\"] missing" );
			}
			
			if(root["Observer"]["participant_2"]["name"].empty()==true){
				throw std::invalid_argument
					("DALS::constructor: error in json input file\n\t    "
						"[\"Observer\"][\"participant_2\"]"
						"[\"name\"] missing" );
			}
			if(root["Observer"]["participant_2"]["path"].empty()==true){
				throw std::invalid_argument
					("DALS::constructor: error in json input file\n\t    "
						"[\"Observer\"][\"participant_2\"]"
						"[\"path\"] missing" );
			}		
			if(root["Observer"]["participant_2"]["latitude"].empty()==true){
				throw std::invalid_argument
					("DALS::constructor: error in json input file\n\t    "
						"[\"Observer\"][\"participant_2\"]"
						"[\"latitude\"] missing" );
			}
			/*
			else{
				double lat=root["Observer"]["participant_2"]["latitude"].asDouble();
				if(lat<2.*M_PI){
					std::cout<<"\n\t\033[33mDALS::constructor: warning for "
						"[\"Observer\"][\"participant_2\"][\"latitude\"]"<<std::endl;
					std::cout<<"\t   Is it expressed in degree? (type yes or no): \033[0m";
					std::string answer;
					std::getline(std::cin, answer);
					if(answer=="no"){
						throw std::invalid_argument
							("DALS::constructor: error in json input file\n\t    "
								"[\"Observer\"][\"participant_2\"]"
								"[\"latitude\"] must be expressed in degree" );
					}
				}
				
			}
			*/
			if(root["Observer"]["participant_2"]["longitude"].empty()==true){
				throw std::invalid_argument
					("DALS::constructor: error in json input file\n\t    "
						"[\"Observer\"][\"participant_2\"]"
						"[\"longitude\"] missing" );
			}
			/*
			else{
				double lon=root["Observer"]["participant_2"]["longitude"].asDouble();
				if(lon<2.*M_PI){
					std::cout<<"\n\t\033[33mDALS::constructor: warning for "
						"[\"Observer\"][\"participant_2\"][\"longitude\"]"<<std::endl;
					std::cout<<"\t   Is it expressed in degree? (type yes or no): \033[0m";
					std::string answer;
					std::getline(std::cin, answer);
					if(answer=="no"){
						throw std::invalid_argument
							("DALS::constructor: error in json input file\n\t    "
								"[\"Observer\"][\"participant_2\"]"
								"[\"longitude\"] must be expressed in degree" );
					}
				}
				
			}
			*/
			if(root["Observer"]["participant_2"]["altitude"].empty()==true){
				throw std::invalid_argument
					("DALS::constructor: error in json input file\n\t    "
						"[\"Observer\"][\"participant_2\"]"
						"[\"altitude\"] missing" );
			}
			else{
				double alt=root["Observer"]["participant_2"]["altitude"].asDouble();
				if(alt>9.){
					std::cout<<"\n\t\033[33mDALS::constructor: warning for "
						"[\"Observer\"][\"participant_2\"][\"altitude\"]"<<std::endl;
					std::cout<<"\t   Is it expressed in meters? (type yes or no): \033[0m";
					std::string answer;
					std::getline(std::cin, answer);
					if(answer=="no"){
						throw std::invalid_argument
							("DALS::constructor: error in json input file\n\t    "
								"[\"Observer\"][\"participant_2\"]"
								"[\"altitude\"] must be expressed in meters" );
					}
				}
				
			}
		}
	}
	
	/**********************************************************/
	/******************* Measurements parameters **************/
	/**********************************************************/
	if(root["Measurements"]["path"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"Measurements\"][\"path\"] missing" );
	}
	if(root["Measurements"]["name"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"Measurements\"][\"name\"] missing" );
	}
	else{
		//Check existance of measurement file
		std::string filename_meas=root["Measurements"]["path"].asString()+
			root["Measurements"]["name"].asString();
		std::ifstream file(filename_meas);
		if (file.is_open()==false){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"Measurement file not found (please check "
					"[\"Measurements\"][\"path\"] and [\"Measurements\"][\"name\"])" );
		}
		else{
			file.close();
		}
	}
	
	if(root["Measurements"]["flag_LTS"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"Measurements\"][\"flag_LTS\"] missing" );
	}
	if(root["Measurements"]["Earth_model"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"Measurements\"][\"Earth_model\"] missing" );
	}
	std::string Earth_model=root["Measurements"]["Earth_model"].asString();
	
	if(root["Measurements"]["kernels"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"Measurements\"][\"kernels\"] missing" );
	}
	
	if(root["Measurements"]["kernels"]["path"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"Measurements\"][\"kernels\"][\"path\"] missing" );
	}
	std::string kernel_path=root["Measurements"]["kernels"]["path"].asString();

	
	if(root["Measurements"]["kernels"]["IAU_assoc"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"Measurements\"][\"kernels\"][\"IAU_assoc\"] missing" );
	}
	else{
		std::string filename_IAU=kernel_path+
			root["Measurements"]["kernels"]["IAU_assoc"].asString();
		std::ifstream file(filename_IAU);
		if (file.is_open()==false){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"IAU association file not found (please check "
					"[\"Measurements\"][\"kernels\"][\"IAU_assoc\"])" );
		}
		else{
			file.close();
		}
	}
	
	if(root["Measurements"]["kernels"]["pck"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"Measurements\"][\"kernels\"][\"pck\"] missing" );
	}
	else{
		std::string filename_pck=kernel_path+
			root["Measurements"]["kernels"]["pck"].asString();
		std::ifstream file(filename_pck);
		if (file.is_open()==false){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"PCK kernel not found (please check "
					"[\"Measurements\"][\"kernels\"][\"pck\"])" );
		}
		else{
			file.close();
		}
	}
	
	if(Earth_model=="ITRF93"){	
		if(root["Measurements"]["kernels"]["ITRF93"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Measurements\"][\"kernels\"][\"ITRF93\"] missing" );
		}
		else{
			std::string filename_EK=kernel_path+
				root["Measurements"]["kernels"]["ITRF93"].asString();
			std::ifstream file(filename_EK);
			if (file.is_open()==false){
				throw std::invalid_argument
					("DALS::constructor: error in json input file\n\t    "
						"High precision Earth kernels not found (please check "
						"[\"Measurements\"][\"kernels\"][\"ITRF93\"])" );
			}
			else{
				file.close();
			}
		}
		if(root["Measurements"]["kernels"]["ITRF93_assoc"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Measurements\"][\"kernels\"][\"ITRF93_assoc\"] missing" );
		}
		else{
			std::string filename_ITRF93_assoc=kernel_path+
				root["Measurements"]["kernels"]["ITRF93_assoc"].asString();
			std::ifstream file(filename_ITRF93_assoc);
			if (file.is_open()==false){
				throw std::invalid_argument
					("DALS::constructor: error in json input file\n\t    "
						"ITRF93 association kernels not found (please check "
						"[\"Measurements\"][\"kernels\"][\"ITRF93_assoc\"])" );
			}
			else{
				file.close();
			}
		}
	}
	
	/**********************************************************/
	/***************** State estimates parameters *************/
	/**********************************************************/
	if(root["OD"]["path"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"OD\"][\"path\"] missing" );
	}
	if(root["OD"]["name"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"OD\"][\"name\"] missing" );
	}
	
	//Check existence of file
	std::string filename_state= root["OD"]["path"].asString()+
		root["OD"]["name"].asString();
	std::ifstream file(filename_state);
	if (file.is_open()==false){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"Invalid state estimates filename (please check "
				"[\"OD\"][\"path\"] and [\"OD\"][\"name\"])" );
	}
	else{
		file.close();
	}
	if(root["OD"]["flag_Cov_av"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"OD\"][\"flag_Cov_av\"] missing" );
	}
	if(root["OD"]["flag_Bfactor_av"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"OD\"][\"flag_Bfactor_av\"] missing" );
	}
	if(root["OD"]["flag_SRPC_av"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"OD\"][\"flag_SRPC_av\"] missing" );
	}

	/**********************************************************/
	/************************ DALS variables ******************/
	/**********************************************************/
	if(root["DALS"]["order"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"DALS\"][\"order\"] missing" );
	}
	else{
		int order=root["DALS"]["order"].asInt();
		if(order<2){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"Invalid [\"DALS\"][\"order\"], please set a value "
					"larger or equal to 2" );
		}
	}
	if(root["DALS"]["max_iter"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t     "
				"[\"DALS\"][\"max_iter\"] missing" );
	}
	else{
		int max_iter=root["DALS"]["max_iter"].asInt();
		if(max_iter<2){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"Invalid [\"DALS\"][\"max_iter\"], please set a value "
					"larger or equal to 2" );
		}
		
	}
	if(root["DALS"]["fun_tol"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t     "
				"[\"DALS\"][\"fun_tol\"] missing" );
	}
	else{
		double fun_tol=root["DALS"]["fun_tol"].asDouble();
		if(fun_tol>1e-3){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"Value of [\"DALS\"][\"fun_tol\"] too large, please insert a"
					" value lower than 1e-3");
		}
		if(fun_tol<0){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"[\"DALS\"][\"fun_tol\"] must be positive");
		}
		
	}
	if(root["DALS"]["step_tol"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t     "
				"[\"DALS\"][\"step_tol\"] missing" );
	}
	else{
		double step_tol=root["DALS"]["step_tol"].asDouble();
		if(step_tol>1e-3){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"Value of [\"DALS\"][\"step_tol\"] too large, please insert a"
					" value lower than 1e-3");
		}
		if(step_tol<0){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"[\"DALS\"][\"step_tol\"] must be positive");
		}
		
	}
	if(root["DALS"]["grad_tol"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t     "
				"[\"DALS\"][\"grad_tol\"] missing" );
	}
	else{
		double grad_tol=root["DALS"]["grad_tol"].asDouble();
		if(grad_tol<0){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"[\"DALS\"][\"grad_tol\"] must be positive");
		}
		
	}
	if(root["DALS"]["map_tol"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t     "
				"[\"DALS\"][\"map_tol\"] missing" );
	}
	else{
		double map_tol=root["DALS"]["map_tol"].asDouble();
		if(map_tol>1e-3){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"Value of [\"DALS\"][\"map_tol\"] too large, please insert a"
					" value lower or equal to 1e-3");
		}
		if(map_tol<0){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"[\"DALS\"][\"map_tol\"] must be positive");
		}
		
	}
	if(root["DALS"]["stop_crit"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t     "
				"[\"DALS\"][\"stop_crit\"] missing" );
	}
	else{
		std::string stop_crit=root["DALS"]["stop_crit"].asString();
		if(stop_crit != "relative" && stop_crit !="gradient"){
			
			std::stringstream buffer;
			buffer << "DALS::constructor: error in json input file \n\t    "
				"Invalid [\"DALS\"][\"stop_crit\"]. ";
			buffer << "DALS stopping criteria:"<<std::endl;
			buffer << "\t\t\"relative\": control on stepsize and"
				" residuals variations"<<std::endl;
			buffer << "\t\t\"gradient\": control on residuals gradient";
			throw std::invalid_argument(buffer.str());
			
		}
		
	}
		
	if(root["DALS"]["date_est"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"DALS\"][\"date_est\"] missing" );
	}

	if(root["DALS"]["flag_Bfactor_est"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t    "
				"[\"DALS\"][\"flag_Bfactor_est\"] missing" );
	}
	if(root["OD"]["flag_Bfactor_av"].asBool()==false && 
		root["DALS"]["flag_Bfactor_est"].asBool()==false){
		
		std::cout<<"\n\t\033[33mDALS::constructor: warning for "
			"[\"DALS\"][\"flag_Bfactor_est\"].";
		std::cout<<"\n\t[\"DALS\"][\"flag_Bfactor_est\"] set to \"false\" by the user, ";	
		std::cout<<"but no Bfactor estimate \n\tis available, [\"DALS\"][\"flag_Bfactor_est\"]";
		std::cout<<" automatically set to \"true\"\033[0m"<<std::endl;
	}
	
	
	if(root["DALS"]["flag_SRPC_est"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t     "
				"[\"DALS\"][\"flag_SRPC_est\"] missing" );
	}
	if(root["OD"]["flag_SRPC_av"].asBool()==false && 
		root["DALS"]["flag_SRPC_est"].asBool()==false){
		
		std::cout<<"\n\t\033[33mDALS::constructor: warning for "
			"[\"DALS\"][\"flag_SRPC_est\"].";
		std::cout<<"\n\t[\"DALS\"][\"flag_SRPC_est\"] set to \"false\" by the user, ";	
		std::cout<<"but no SRPC estimate \n\tis available, [\"DALS\"][\"flag_SRPC_est\"]";
		std::cout<<" automatically set to \"true\"\033[0m"<<std::endl;
	}
	if(root["DALS"]["output"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t     "
				"[\"DALS\"][\"output\"] missing" );
	}
	if(root["DALS"]["output"]["path"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t     "
				"[\"DALS\"][\"output\"][\"path\"] missing" );
	}
	if(root["DALS"]["output"]["name"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t     "
				"[\"DALS\"][\"output\"][\"name\"] missing" );
	}
	if(root["DALS"]["output"]["obj_name"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t     "
				"[\"DALS\"][\"output\"][\"obj_name\"] missing" );
	}
	if(root["DALS"]["output"]["obj_ID"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t     "
				"[\"DALS\"][\"output\"][\"obj_ID\"] missing" );
	}
	if(root["DALS"]["output"]["flag_Cov"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t     "
				"[\"DALS\"][\"output\"][\"flag_Cov\"] missing" );
	}
	if(root["DALS"]["output"]["flag_KeplPar"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t     "
				"[\"DALS\"][\"output\"][\"flag_KeplPar\"] missing" );
	}

	/**********************************************************/
	/******************** Dynamics variables ******************/
	/**********************************************************/
	std::string method;
	if(root["Dynamics"]["method"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t     "
				"[\"Dynamics\"][\"method\"] missing" );
	}
	else{
		method=root["Dynamics"]["method"].asString();
		if(method !="AIDA" && method !="SADA"){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"Invalid [\"Dynamics\"][\"method\"],"
					" please select either \"AIDA\" or \"SADA\"" );
		}
	}
	if(root["Dynamics"]["tolerance"].empty()==true){
		throw std::invalid_argument
			("DALS::constructor: error in json input file\n\t     "
				"[\"Dynamics\"][\"tolerance\"] missing" );
	}
	else{
		double prop_tol=root["DALS"]["tolerance"].asDouble();
		if(prop_tol>1e-6){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"Value of [\"Dynamics\"][\"tolerance\"] too large, please insert a"
					" value lower than 1e-6");
		}
		if(prop_tol<0){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"[\"Dynamics\"][\"tolerance\"] must be positive");
		}
		
	}
	if(method=="AIDA"){
		if(root["Dynamics"]["AIDAparam"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"[\"Dynamics\"][\"AIDAparam\"] missing" );
		}
		
		std::string gravmodel;
		if(root["Dynamics"]["AIDAparam"]["gravmodel"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"[\"Dynamics\"][\"AIDAparam\"][\"gravmodel\"] missing" );
		}
		else{
			gravmodel=root["Dynamics"]["AIDAparam"]
				["gravmodel"]["path"].asString()+root["Dynamics"]["AIDAparam"]
				["gravmodel"]["name"].asString();
			struct stat buffer; 
			if(stat (gravmodel.c_str(), &buffer) == 0){
				throw std::invalid_argument
					("DALS::constructor: error in json input file\n\t     "
						"Invalid AIDA gravitational model filename, please "
						"select either \"egm96\" or \"egm2008\"" );	
			}
		}
		if(root["Dynamics"]["AIDAparam"]["gravmodel"]["order"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"[\"Dynamics\"][\"AIDAparam\"][\"gravmodel\"][\"order\"] missing" );
		}
		if(root["Dynamics"]["AIDAparam"]["flag_drag"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"[\"Dynamics\"][\"AIDAparam\"][\"flag_drag\"] missing" );
		}
		else{
			int flag_drag=root["Dynamics"]["AIDAparam"]["flag_drag"].asInt();
			if(flag_drag!=0 && flag_drag!=1 && flag_drag!=2){
				
				std::stringstream buffer;
				buffer << "DALS::constructor: error in json input file \n\t    "
					"Invalid [\"Dynamics\"][\"AIDAparam\"][\"flag_drag\"]. ";
				buffer << "AIDA atmospheric drag perturbations flags:"<<std::endl;
				buffer << "\t\t0: none"<<std::endl;
				buffer << "\t\t1: active, non rotating atmosphere"<<std::endl;
				buffer << "\t\t2: active, rotating atmosphere";
				throw std::invalid_argument(buffer.str());
				
			}
		}
		if(root["Dynamics"]["AIDAparam"]["flag_SRP"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Dynamics\"][\"AIDAparam\"][\"flag_SRP\"] missing" );
		}
		else{
			int flag_SRP=root["Dynamics"]["AIDAparam"]["flag_SRP"].asInt();
			if(flag_SRP!=0 && flag_SRP!=1 && flag_SRP!=2 && flag_SRP!=3 &&
				flag_SRP!=4 && flag_SRP!=5 && flag_SRP!=6){
				
				std::stringstream buffer;
				buffer << "DALS::constructor: error in json input file \n\t   "
					" Invalid [\"Dynamics\"][\"AIDAparam\"][\"flag_SRP\"]. ";
				buffer << "AIDA SRP perturbations flags:"<<std::endl;
				buffer << "\t\t0: none"<<std::endl;
				buffer << "\t\t1: active (no shadow)"<<std::endl;
				buffer << "\t\t2: active, cylindrical Earth shadow"<<std::endl;
				buffer << "\t\t3: active, biconical Earth shadow"<<std::endl;
				buffer << "\t\t4: active, cylindrical Earth and Moon shadows"<<std::endl;
				buffer << "\t\t5: active, biconical Earth shadow and cylindrical Moon shadow"<<std::endl;
				buffer << "\t\t6: active, biconical Earth and Moon shadows";
				throw std::invalid_argument(buffer.str());
				
			}
		}
		if(root["Dynamics"]["AIDAparam"]["flag_thirdbody"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"[\"Dynamics\"][\"AIDAparam\"][\"flag_thirdbody\"] missing" );
		}
		else{
			int flag_third=root["Dynamics"]["AIDAparam"]["flag_thirdbody"].asInt();
			if(flag_third!=0 && flag_third!=1 && flag_third!=2){
				
				std::stringstream buffer;
				buffer << "DALS::constructor: error in json input file \n\t    "
					"Invalid [\"Dynamics\"][\"AIDAparam\"][\"flag_thirdbody\"]. ";
				buffer << "AIDA third-body perturbations flags:"<<std::endl;
				buffer << "\t\t0: none"<<std::endl;
				buffer << "\t\t1: Moon"<<std::endl;
				buffer << "\t\t2: Moon and Sun";
				throw std::invalid_argument(buffer.str());
				
			}
		}
	}
	else{
		if(root["Dynamics"]["SADAparam"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Dynamics\"][\"SADAparam\"] missing" );
		}
		if(root["Dynamics"]["SADAparam"]["flag_linear"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Dynamics\"][\"SADAparam\"][\"flag_linear\"] missing" );
		}
		if(root["Dynamics"]["SADAparam"]["flag_SunMoon"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Dynamics\"][\"SADAparam\"][\"flag_SunMoon\"] missing" );
		}
		else{
			int flag_third=root["Dynamics"]["SADAparam"]["flag_SunMoon"].asInt();
			if(flag_third!=0 && flag_third!=1 && flag_third!=2){
				std::stringstream buffer;
				buffer << "DALS::constructor: error in json input file \n\t    "
					"Invalid [\"Dynamics\"][\"SADAparam\"][\"flag_SunMoon\"]. ";
				buffer << "SADA third-body perturbations flags:"<<std::endl;
				buffer << "\t\t0: none"<<std::endl;
				buffer << "\t\t1: Moon and Sun, analytical ephemeris"<<std::endl;
				buffer << "\t\t2: Moon and Sun, CSPICE ephemeris";
				throw std::invalid_argument(buffer.str());
				
			}
		}
		if(root["Dynamics"]["SADAparam"]["flag_drag"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"[\"Dynamics\"][\"SADAparam\"][\"flag_drag\"] missing" );
		}
		else{
			int flag_drag=root["Dynamics"]["SADAparam"]["flag_drag"].asInt();
			if(flag_drag!=0 && flag_drag!=1){
				std::stringstream buffer;
				buffer << "DALS::constructor: error in json input file\n\t     "
					"Invalid [\"Dynamics\"][\"SADAparam\"][\"flag_drag\"]. ";
				buffer << "SADA atmospheric drag perturbations flags:"<<std::endl;
				buffer << "\t\t0: none"<<std::endl;
				buffer << "\t\t1: active, rotating atmosphere";
				throw std::invalid_argument(buffer.str());
				
			}
		}
		if(root["Dynamics"]["SADAparam"]["flag_SRP"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t    "
					"[\"Dynamics\"][\"SADAparam\"][\"flag_SRP\"] missing" );
		}
		else{
			int flag_SRP=root["Dynamics"]["SADAparam"]["flag_SRP"].asInt();
			if(flag_SRP!=0 && flag_SRP!=1){
				std::stringstream buffer;
				buffer << "DALS::constructor: error in json input file \n\t    "
					"Invalid [\"Dynamics\"][\"SADAparam\"][\"flag_SRP\"]. ";
				buffer << "SADA SRP perturbations flags:"<<std::endl;
				buffer << "\t\t0: none"<<std::endl;
				buffer << "\t\t1: active (no shadow)";
				throw std::invalid_argument(buffer.str());
				
			}
		}
		if(root["Dynamics"]["SADAparam"]["flag_shortPeriodics"].empty()==true){
			throw std::invalid_argument
				("DALS::constructor: error in json input file\n\t     "
					"[\"Dynamics\"][\"SADAparam\"][\"flag_shortPeriodics\"] missing" );
		}
		else{
			int flag_shortP=root["Dynamics"]["SADAparam"]["flag_shortPeriodics"].asInt();
			if(flag_shortP!=0 && flag_shortP!=2){
				std::stringstream buffer;
				buffer << "DALS::constructor: error in json input file "
					"\n\t    Invalid [\"Dynamics\"][\"SADAparam\"][\"flag_shortPeriodics\"]. ";
				buffer << "SADA short periodics perturbations flags:"<<std::endl;
				buffer << "\t\t0: J2 only"<<std::endl;
				buffer << "\t\t2: J2, Moon and Sun";
				throw std::invalid_argument(buffer.str());
				
			}
		}
	}
	
	
}

void meas_except(const Observations<double>& Observation,
	const Observer_param& param_obs){
		
	std::vector<double> et_vect=Observation.et_vect;
	const unsigned int N_inst=et_vect.size();
	if(et_vect.size()==0){
		throw std::invalid_argument
		("DALS::read_measurements: no time instant available" );
	}
	
	//Optical sensor
	if(param_obs.type=="optical"){
		if(param_obs.sensor_setup.ra.availability==true){
			unsigned int L=Observation.ra_vect.size();
			if(L != N_inst){
				throw std::invalid_argument
				("DALS::read_measurements: ra measurements not available "
					"for all time instants");
			}
		}
		if(param_obs.sensor_setup.decl.availability==true){
			unsigned int L=Observation.decl_vect.size();
			if(L != N_inst){
				throw std::invalid_argument
				("DALS::read_measurements: decl measurements not available for "
					"all time instants");
			}
		}
	}
	else if(param_obs.type=="radar"){
		if(param_obs.sensor_setup.az.availability==true){
			unsigned int L=Observation.az_vect.size();
			if(L != N_inst){
				throw std::invalid_argument
				("DALS::read_measurements: az measurements not available for all"
					" time instants");
			}
		}
		
		if(param_obs.sensor_setup.el.availability==true){
			unsigned int L=Observation.el_vect.size();
			if(L != N_inst){
				throw std::invalid_argument
				("DALS::read_measurements: el measurements not available for all "
					"time instants");
			}
		}
		
		if(param_obs.sensor_setup.doppler.availability==true){
			unsigned int L=Observation.doppler_vect.size();
			if(L != N_inst){
				throw std::invalid_argument
				("DALS::read_measurements: doppler measurements not available for "
					"all time instants");
			}
		}
		
		if(param_obs.sensor_setup.range.availability==true){
			unsigned int L=Observation.az_vect.size();
			if(L != N_inst){
				throw std::invalid_argument
				("DALS::read_measurements: range measurements not available for all"
					" time instants");
			}
		}
	}
		
}

void state_except(const Estimate_param& param_OD){
	

	const bool flag_Cov_av           = param_OD.flag_Cov_av;
	const bool flag_Bfactor_av       = param_OD.flag_Bfactor_av;
	const bool flag_SRPC_av          = param_OD.flag_SRPC_av;
	const std::string filename_state = param_OD.filename_state;
	
	// theoretical number of lines
	const unsigned int n_lines_min=2;
	unsigned int n_lines_th=n_lines_min;
	if(flag_Cov_av==true){
		n_lines_th=n_lines_th+6;
	}
	if(flag_Bfactor_av==true){
		n_lines_th=n_lines_th+1;
	}
	if(flag_SRPC_av==true){
		n_lines_th=n_lines_th+1;
	}
	
	//Count nr of available time instants
	unsigned int n_lines = 0;
	std::string line;
	std::ifstream myfile(filename_state);
	while (std::getline(myfile, line))
		++n_lines;
	myfile.close();
	
	if(n_lines<2){
		throw std::invalid_argument
			("DALS::read_state_data: error in state estimates file"
			"\n\t     Time epoch or state vector missing");
	}
	if(n_lines !=n_lines_th){
		std::stringstream buffer;
		buffer << "DALS::read_state_data: error in state estimates file.\n\t    "
			"The number of entries of the state estimates file ("<<n_lines<<") ";
		buffer <<"does not match the expected one ("<<n_lines_th<<").";
		buffer <<"\n\t    Please check the [\"OD\"][\"flag_Cov_av\"],[\"OD\"][\"flag_Bfactor_av\"] ";
		buffer <<"and \n\t    [\"OD\"][\"flag_SRPC_av\"] flags in the json "
			"input file and try again.";
		throw std::invalid_argument(buffer.str());
	}
	
	
}
#endif