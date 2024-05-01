#ifndef SpiceRoutines_H_
#define SpiceRoutines_H_

#include<dace/dace.h>
#include<cspice/SpiceUsr.h>


std::vector<unsigned int> Read_catalogue(const std::string& file_path,
	const bool flag_print=false){

	//Data for TLE reading
	SpiceInt lineln = 70;

	//Count nr of objects
	unsigned int n_obj = 0;
	std::string line;
	std::ifstream myfile(file_path);
	while (getline(myfile, line)){
		n_obj=n_obj+1;
	}
	myfile.close();
	n_obj=n_obj/2;  
	if(flag_print==true){
		std::cout<<"Number of objects in TLE catalogue: "<<n_obj<<std::endl; 
	}		

	SpiceChar lines[2*n_obj][lineln];

	DACE::AlgebraicVector<std::string> lines_stack(2*n_obj);

	myfile.open(file_path);
	for(unsigned int i = 0; i < 2*n_obj; i++ ){
		std::string temp_line, temp_line2;
		getline(myfile,temp_line,'\n');

		for(unsigned int j=0;j<lineln-1;j++){
			temp_line2.push_back(temp_line[j]);
		}		
		temp_line2.push_back('\0');

		//Append
		lines_stack[i]=temp_line2;

	}
	myfile.close();
	
	//Read TLE elements and save ID
	std::vector<unsigned int> ID_list;
	for(unsigned int i=0;i<n_obj;i++){
		//Read ID of TLE entry
		std::string IDobj_str=lines_stack[2*i].substr(2,5);
		unsigned int IDobj=std::stoi(IDobj_str);
		
		//Store ID
		ID_list.push_back(IDobj);	
	}
	
	return ID_list;
}


SpiceDouble Get_et(const std::string& file_path, const unsigned int ID_target){

	//Data for TLE reading
	SpiceInt frstyr = 1957;
	SpiceInt lineln = 70;

	//Count nr of objects
	unsigned int n_obj = 0;
	std::string line;
	std::ifstream myfile(file_path);
	while (getline(myfile, line)){
		n_obj=n_obj+1;
	}
	myfile.close();
	n_obj=n_obj/2;        
	//std::cout<<"Number of objects in TLE catalogue: "<<n_obj<<std::endl;        

	SpiceChar lines[2*n_obj][lineln];

	DACE::AlgebraicVector<std::string> lines_stack(2*n_obj);
	SpiceDouble et;
	double elems[10];

	myfile.open(file_path);
	for(unsigned int i = 0; i < 2*n_obj; i++ ){
		std::string temp_line, temp_line2;
		std::getline(myfile,temp_line,'\n');

		for(unsigned int j=0;j<lineln-1;j++){
			temp_line2.push_back(temp_line[j]);
		}		
		temp_line2.push_back('\0');

		//Append
		lines_stack[i]=temp_line2;

	}
	myfile.close();
	
	//Initialize line index
	unsigned int index=n_obj*2;

	//Read TLE elements and reference et of target objects
	for(unsigned int i=0;i<n_obj;i++){
		//Read ID of TLE entry
		std::string IDobj_str=lines_stack[2*i].substr(2,5);
		unsigned int IDobj=stoi(IDobj_str);
		if (IDobj==ID_target){
			index=i;
			break;		
		}			
	}
	//Compute reference epoch
	if(index==2*n_obj){
		throw std::invalid_argument
			("Error in Get_et:: object not found in the catalogue" );
	}
	else{
		//Convert to char and add to lines
		for(unsigned int i=0;i<2;i++){
			std::strcpy(lines[i],lines_stack[2*index+i].c_str());
		}
		std::cout<<"\tTarget object: "<<ID_target<<std::endl;
		getelm_c (frstyr,lineln,lines, &et, elems); 

		SpiceChar output[TIMLEN];
		timout_c(et,TIMFMT,TIMLEN,output);
		std::cout<<"\t    TLE reference date: "<<output<<std::endl;
		return et;
	}
}

std::vector<double> Get_elems(const std::string& file_path, const unsigned int ID_target){

	//Data for TLE reading
	SpiceInt frstyr = 1957;
	SpiceInt lineln = 70;

	//Count nr of objects
	unsigned int n_obj = 0;
	std::string line;
	std::ifstream myfile(file_path);
	while (getline(myfile, line)){
		n_obj=n_obj+1;
	}
	myfile.close();
	n_obj=n_obj/2;        
	//std::cout<<"Number of objects in TLE catalogue: "<<n_obj<<std::endl;        

	SpiceChar lines[2*n_obj][lineln];

	DACE::AlgebraicVector<std::string> lines_stack(2*n_obj);
	SpiceDouble et;
	SpiceDouble elems[10];

	myfile.open(file_path);
	for(unsigned int i = 0; i < 2*n_obj; i++ ){
		std::string temp_line, temp_line2;
		std::getline(myfile,temp_line,'\n');

		for(unsigned int j=0;j<lineln-1;j++){
			temp_line2.push_back(temp_line[j]);
		}		
		temp_line2.push_back('\0');

		//Append
		lines_stack[i]=temp_line2;

	}
	myfile.close();
	
	//Initialize line index
	unsigned int index=n_obj*2;

	//Read TLE elements and reference et of target objects
	for(unsigned int i=0;i<n_obj;i++){
		//Read ID of TLE entry
		std::string IDobj_str=lines_stack[2*i].substr(2,5);
		unsigned int IDobj=stoi(IDobj_str);
		if (IDobj==ID_target){
			index=i;
			break;		
		}			
	}
	//Compute reference epoch
	std::vector<double> elements(10);
	if(index==2*n_obj){
		throw std::invalid_argument
			("Error in Get_et:: object not found in the catalogue" );
	}
	else{
		//Convert to char and add to lines
		for(unsigned int i=0;i<2;i++){
			std::strcpy(lines[i],lines_stack[2*index+i].c_str());
		}
		getelm_c (frstyr,lineln,lines, &et, elems); 
		for(int i=0;i<10;i++){
			elements.at(i)=elems[i];
		}
		return elements;
	}
}

void Create_kernel(const std::vector<DACE::AlgebraicVector<double>>& states_vec,
		const DACE::AlgebraicVector<double>& et_vec, 
		const std::string& kernel_name, const std::string& sat_name,
		const std::string& LSK){
	

	//Write state vectors in .txt file
	std::string states_file = "States.txt";
	
    remove(states_file.c_str());  // Maybe shall we check what integer is returned for errors?
    std::ofstream myfile0 (states_file.c_str());
	const unsigned int N_obs=et_vec.size();
    if (myfile0.is_open()){
        for(unsigned int i=0; i<N_obs; i++) {
            DACE::AlgebraicVector<double> state=states_vec.at(i);
            myfile0  << std::fixed << std::setprecision(15) ;
            myfile0  <<  et_vec[i] << " "  << state[0] << " "  << 
				state[1] << " "  << state[2] << " " << 
				state[3] << " "  << state[4] << " "  << 
				state[5] << std::endl;
        }
        myfile0.close();
    }
    else{
		std::cout << "Error state file" << std::endl;
    }
    
    //Create preparatory file
	std::string setup_file="Setup_SPK.txt";
    std::ofstream myfile (setup_file); //, ofstream::trunc
    std::string kepkernel = kernel_name;
    remove(kepkernel.c_str());  // Maybe we shall check what integer is returned for errors?
    
	// Can't we move this to a dedicated software? If so, 
	//I would suggest a Python script where we could just 
	//replace keywords in a template
    if (myfile.is_open())
    {
        myfile << "\\begindata" <<std::endl;
        myfile << "INPUT_DATA_TYPE   = 'STATES'" <<std::endl;
        myfile << "OUTPUT_SPK_TYPE   = 13"<<std::endl;
        myfile << "OBJECT_ID         = "+sat_name <<std::endl;
        myfile << "CENTER_ID         = 399" << std::endl;
        myfile << "REF_FRAME_NAME    = 'J2000'" <<std::endl;
        myfile << "PRODUCER_ID       = 'Matteo'" <<std::endl;
        myfile << "LEAPSECONDS_FILE  = '"+LSK+"'"<<std::endl;
        myfile << "DATA_ORDER        = 'EPOCH X Y Z VX VY VZ'"<<std::endl;
        myfile << "INPUT_DATA_UNITS  = ('ANGLES=DEGREES' 'DISTANCES=km')" << std::endl;
        myfile << "DATA_DELIMITER    = ' '"<<std::endl;
        myfile << "LINES_PER_RECORD  = 1" << std::endl;
        myfile << "IGNORE_FIRST_LINE = 0" << std::endl;
        myfile << "INPUT_DATA_FILE   = '"+states_file+"'" << std::endl;
        myfile << "OUTPUT_SPK_FILE   = '"+kepkernel+"'"<<std::endl;
        myfile << "POLYNOM_DEGREE    = 15" << std::endl;
        myfile << "TIME_WRAPPER = '# ETSECONDS'" << std::endl;
        myfile.close();
    }
    else{
		std::cout << "Unable to open file";
    }
    std::string setupSPK = "./mkspk.exe -setup "+setup_file+" > nul";
	//std::string setupSPK = "./mkspk -setup "+setup_file;
    std::system(setupSPK.c_str()); // Is a call to system, we shall check what happens with the return value

	//Delete auxiliary files
	remove(states_file.c_str());
	remove(setup_file.c_str());
}

bool Create_TLEkernel(const std::string& TLE_filename, const unsigned int ID){
	
	//Open kernel
	std::string BSP_filename=TLE_filename.substr(0,
		TLE_filename.size()-4)+".bsp";
	remove(BSP_filename.c_str());
	
	SpiceInt handle;
	std::string aux_name="SPK_file";
	spkopn_c(BSP_filename.c_str(),aux_name.c_str(),5000,&handle);
	
	//Set variables
	SpiceDouble consts[8];
	consts[0]=0.001082616;
	consts[1]=-0.00000253881;
	consts[2]=-0.00000165597;
	consts[3]=0.074366916133173;
	consts[4]=120.;
	consts[5]=78;
	consts[6]=6378.135;
	consts[7]=1.;
	
	//Frame
	std::string frame="J2000";
	
	//Center
	SpiceInt center=399;
	
	//Segment identifier
	std::string segid = "SEGMENT";
	
	//Number of sets
	SpiceInt n=1;
	
	//Vector of ID
	std::vector<unsigned int> ID_list = 
		Read_catalogue(TLE_filename);
	
	bool flag_create=false;
	for(unsigned int i=0;i<ID_list.size();i++){
		if(ID_list.at(i)==ID){
			//Elements
			std::vector<double> elements=Get_elems(TLE_filename,ID_list.at(i));
			SpiceDouble elems[10];
			for(unsigned int j=0;j<10;j++){
				elems[j]=elements.at(j);
			}
			
			//Start and end epochs
			double et_start=elements[9];
			double et_end=et_start+1*365.25*spd_c();
			
			SpiceDouble epochs[1];
			epochs[0]=et_start;

			//Object ID
			SpiceInt ID=-1*(ID_list.at(i)+100000);
			spkw10_c(handle, ID, center, frame.c_str(), et_start, et_end, 
				segid.c_str(), consts, n, elems, epochs); 
			
			flag_create=true;
					 
			break;
		}
		
	}
	spkcls_c(handle);
	
	return flag_create;
	
}



template<typename T>
T vsep_DACE(const DACE::AlgebraicVector<T>& v1, const DACE::AlgebraicVector<T>& v2){
	
	
	T     vsep;
	double     dmag1;
	double     dmag2;
	DACE::AlgebraicVector<T>     vtemp(3);
	DACE::AlgebraicVector<T>      u1(3);
	DACE::AlgebraicVector<T>      u2(3);

	//Compute the magnitude of v1 and associated unit vector
	dmag1=DACE::cons(vnorm(v1));
	u1=v1/dmag1;

	if(dmag1 == 0.0){
		vsep = 0.0;
		return vsep;
	}

	//Compute the magnitude of v2 and associated unit vector
	dmag2=DACE::cons(vnorm(v2));
	u2=v2/dmag2;

	if(dmag2 == 0.0){
		vsep = 0.0;
		return vsep;
	}
	
	
	if(DACE::cons(dot(u1,u2)) > 0.){
		vtemp[0] = u1[0] - u2[0];
		vtemp[1] = u1[1] - u2[1];
		vtemp[2] = u1[2] - u2[2];
		
		//Separation angle
		vsep = 2.00 * asin (0.50 * vnorm(vtemp));
	}
	else if(DACE::cons(dot(u1,u2)) < 0.){
		vtemp[0] = u1[0] + u2[0];
		vtemp[1] = u1[1] + u2[1];
		vtemp[2] = u1[2] + u2[2];
		
		//Separation angle
		vsep = pi_c() - 2.00 * asin (0.50 * vnorm(vtemp));
	}
	else{
		//Separation angle
		vsep = M_PI;
	}
	return vsep;
	
}

template<typename T>
DACE::AlgebraicVector<T> vrotv_DACE(const DACE::AlgebraicVector<T>& v, 
		const DACE::AlgebraicVector<T>& axis, const T angle){
	
	using std::cos; using DACE::cos;
	using std::sin; using DACE::sin;
	
	//Normalize axis
	DACE::AlgebraicVector<T> axis_norm=DACE::normalize(axis);
	
	DACE::AlgebraicMatrix<T> R(3,3);
	T u_x=axis_norm[0];
	T u_y=axis_norm[1];
	T u_z=axis_norm[2];

	//Build rotation matrix
	R.at(0,0)=cos(angle)+DACE::sqr(u_x)*(1.-cos(angle));
	R.at(0,1)=u_x*u_y*(1.-cos(angle))-u_z*sin(angle);
	R.at(0,2)=u_x*u_z*(1.-cos(angle))+u_y*sin(angle);
	
	R.at(1,0)=u_x*u_y*(1.-cos(angle))+u_z*sin(angle);
	R.at(1,1)=cos(angle)+DACE::sqr(u_y)*(1.-cos(angle));
	R.at(1,2)=u_y*u_z*(1.-cos(angle))-u_x*sin(angle);
	
	R.at(2,0)=u_x*u_z*(1.-cos(angle))-u_y*sin(angle);
	R.at(2,1)=u_y*u_z*(1.-cos(angle))+u_x*sin(angle);
	R.at(2,2)=cos(angle)+DACE::sqr(u_z)*(1.-cos(angle));
	
	DACE::AlgebraicVector<T> v_rot=R*v;
	
	return v_rot;
}

#endif