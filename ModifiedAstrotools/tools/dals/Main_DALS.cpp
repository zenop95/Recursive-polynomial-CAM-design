#define _USE_MATH_DEFINES
#define WITH_ALGEBRAICMATRIX

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <random>
#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <time.h>
#include <chrono>
#include <iterator>
#include <stdio.h>

//JSON
#include <json/json.h>

//Spice
#include <cspice/SpiceUsr.h>
#include <cspice/SpiceZfc.h>

//DACE
#include <dace/dace.h>

//My libraries
#include "DALS.hpp"
#include "SimRoutines.hpp"

int main(int argc, char *const argv[]){
	
	// Load SPICE kernel
	kclear_c();
	cfg::EnvConfig().LoadSpiceKernels();
	
	//Exception pointer
	std::exception_ptr eptr;

	//Input json file
	std::string json_LS  = argv[1];
	try{
		//******************************************************************
		//************************** LS process ****************************
		//******************************************************************
		std::cout<<"\033[36mRunning LS solver: \033[0m"<<std::endl;
		
		//Update json file
		lssolver::DALS LS(json_LS);

		//Run LS process
		auto [et_est,state_est,Cov_est,Bfactor_est,SRPC_est,
			Cov_param,Residuals, flag_dec]=LS.run();
		
		if(argc>3){
			//Reference file
			std::string json_sim     = argv[2];
			std::string filename_ref = argv[3];
			
			std::cout<<"\033[36mResults accuracy: \033[0m"<<std::endl;
			Compute_accuracy (state_est,et_est,Bfactor_est, SRPC_est,
				filename_ref,json_sim,json_LS);
				
			Compare_parameters(json_sim, json_LS, Bfactor_est,
				SRPC_est, Cov_param, flag_dec);	
		}
		
	}
	catch(...){
		eptr = std::current_exception(); // capture
	}
	handle_eptr(eptr);

	return 0;
}


