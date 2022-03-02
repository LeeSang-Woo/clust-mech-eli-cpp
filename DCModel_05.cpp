#include <iostream>
#define TRUE 1
#define FALSE 0

#include "DCModel_struct.h"
#include "DCModel_struct_2D.h"
#include "DCModel_struct_abn.h"
#include "DCMdata_record.h"
#include "DCMdata_test_2D.h"

#define MAX_LOADSTRING 100


int main() {

	unsigned __int64 timeStep = 0;

	int ci = 0;
	int ci_end = 5;
	int pi = 0;
	int pi_init = 0;
	int pi_end = 0;
	unsigned long seed = 0;
	int repeat_num = 5;
	std::vector<int> selected_index;
	std::vector<double> mu_index;
	std::vector<int> size_index;
	std::vector<double> inva_index;
	std::vector<double> prob_index;
	std::vector<double> rand_lambda_index;
	std::vector<double> angle_index;
	int scaled_mu_flag = 0;
	
	std::ifstream ifs(L"DCModel_container_info.txt");

	if (ifs) {
		std::string str;

		std::vector<std::string> param_list;
		param_list.push_back("Start index");		// 0
		param_list.push_back("End index");			// 1
		param_list.push_back("Start parameter");	// 2
		param_list.push_back("End parameter");		// 3
		param_list.push_back("Selected index");		// 4
		param_list.push_back("Random seed");		// 5
		param_list.push_back("Repeat number");		// 6
		param_list.push_back("Selected mu");		// 7
		param_list.push_back("Selected size");		// 8
		param_list.push_back("Scaled mu flag");		// 9
		param_list.push_back("Selected invasion");	// 10
		param_list.push_back("Selected probability");	// 11
		param_list.push_back("Selected randomness of lambda");	// 12
		param_list.push_back("Selected polarity angle");	// 13
		
		while(std::getline(ifs, str)) {

			std::string token1;
			std::string token2;
			std::istringstream stream(str);
			
			std::string token3;
			std::istringstream stream2;
			
			std::getline(stream, token1, ':');
			std::getline(stream, token2);

			for (unsigned int i=0; i<param_list.size(); i++) {
				if (token1.compare(param_list[i]) == 0 && !token2.empty()) {
					switch (i) {
					case 0: ci = std::stoi(token2); break;
					case 1: ci_end = std::stoi(token2); break;
					case 2: pi = pi_init = std::stoi(token2); break;
					case 3: pi_end = std::stoi(token2); break;
					case 4: stream2 = std::istringstream(token2); while (std::getline(stream2, token3, ',')) selected_index.push_back(std::stoi(token3)); break;
					case 5: seed = std::stoul(token2); break;
					case 6: repeat_num = std::stoi(token2); break;
					case 7: stream2 = std::istringstream(token2); while (std::getline(stream2, token3, ',')) mu_index.push_back(std::stod(token3)); break;
					case 8: stream2 = std::istringstream(token2); while (std::getline(stream2, token3, ',')) size_index.push_back(std::stoi(token3)); break;
					case 9: if (std::stoi(token2) == 1) scaled_mu_flag = 1; break;
					case 10: stream2 = std::istringstream(token2); while (std::getline(stream2, token3, ',')) inva_index.push_back(std::stod(token3)); break;
					case 11: stream2 = std::istringstream(token2); while (std::getline(stream2, token3, ',')) prob_index.push_back(std::stod(token3)); break;
					case 12: stream2 = std::istringstream(token2); while (std::getline(stream2, token3, ',')) rand_lambda_index.push_back(std::stod(token3)); break;
					case 13: stream2 = std::istringstream(token2); while (std::getline(stream2, token3, ',')) angle_index.push_back(std::stod(token3)); break;
					}
				}
			}
		}
	}

	if (seed > 0) SF.seeding_WellRNG512(seed);


	DCMcontainerAB* container2D = new DCMcontainerAB();
	DCMdata_FuncAB rec;

	container2D->initialize_loadfile();
	DCMInit(*container2D);
	container2D->initialize(ci);
	rec.setContainer(*container2D, ci+1);

	double abnormal_line_coef;
	int abnormal_cells_size_flag;
	int ci_rem = ci;
	int ci_div = ci;

	if (repeat_num > 0) {
		ci_rem = ci%repeat_num;
		ci_div = ci/repeat_num;
	}

	if (mu_index.size() > 0) {
		abnormal_line_coef = mu_index[ci_rem];
	}
	else abnormal_line_coef = 0;

	if (size_index.size() > 0) abnormal_cells_size_flag = size_index[ci_div];
	else abnormal_cells_size_flag = ci_div * 25;
	
	container2D->initialize_cellpopulation(abnormal_line_coef, abnormal_cells_size_flag);	

	while(1) {
			
		container2D->vDynamicsMain();

		int time = static_cast<int>(1.0/PM.Delta_t);
																					
		if (timeStep != 0) {

			if (container2D->checkContainer(1) == 1) {
				std::cout<<"Time step: "<<timeStep/1000<<"\t"<<"Cells size: "<<container2D->cellsSize_r<<std::endl;
			}

			if (container2D->checkContainer(3) == 1) {
				rec.Record_AbnormalPopulation(timeStep);
			}
		}

		timeStep++;
		
		int finish_flag = 0;

		if (container2D->finishContainer() > 0) finish_flag = 1;
		
		if (finish_flag > 0) {
							
			bool check_reset = false;

			if (container2D->finishContainer(0) > 10) {
				std::cout<<ci<<","<<pi<<" out..."<<std::endl;
			}
			else {
				if (container2D->finishContainer(0) == 3) {
					rec.Record_AbnormalPopulation(timeStep);
				}
				std::cout<<ci<<","<<pi<<" success..."<<std::endl;
			}

			timeStep = 0;

			delete container2D;

			bool check_ci = false;

			if (pi_end > 0) {
				if (pi < pi_end) pi++;
				else {
					pi = pi_init;
					check_ci = true;					
				}
			}
			else check_ci = true;
			
			if (check_ci) {
				if (ci < ci_end-1) {
					if (!check_reset) ci++;
				}
				else break;
			}

			container2D = new DCMcontainerAB();

			if (seed > 0) SF.seeding_WellRNG512(seed);
			else SF.initialize();
			
			container2D->initialize_loadfile();
			DCMInit(*container2D);
			container2D->initialize(ci);
			rec.setContainer(*container2D, ci+1);

			// Population
			ci_rem = ci;
			ci_div = ci;
			if (repeat_num > 0) {
				ci_rem = ci%repeat_num;
				ci_div = ci/repeat_num;
			}

			if (mu_index.size() > 0) {
				abnormal_line_coef = mu_index[ci_rem];
			}
			else abnormal_line_coef = 0;

			if (size_index.size() > 0) abnormal_cells_size_flag = size_index[ci_div];
			else abnormal_cells_size_flag = ci_div * 25;

			container2D->initialize_cellpopulation(abnormal_line_coef, abnormal_cells_size_flag);

		}
	}
							
	SharingFunctions::releasePointer();
	DCMparameters::releasePointer();

	return 0;
}
