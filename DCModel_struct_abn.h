#ifndef DCMODELSTRUCTABNORMAL_HH
#define DCMODELSTRUCTABNORMAL_HH


#include "DCModel_struct_2D.h"



namespace AC {
	enum AbnormalCellFlag {

		RANDOM_ABNORMAL			= 1<<0,
		ONSET_STOP_DIVISION		= 1<<1,
		LINETENSION_ON			= 1<<2,
		BOUNDARY_CELL_ELI		= 1<<3,
		ABNORMAL_MECH_PARAM		= 1<<4,
		FINISH_PARAM			= 1<<5,
		
	};
}


class DCMpointAB : public DCMpoint2D {

public:

	DCMpointAB() : DCMpoint2D() {
	}

	DCMpointAB(const DCMpointAB& otherPoint) : DCMpoint2D(otherPoint) {
		copyDCMpointAB(otherPoint);
	}

	DCMpointAB& operator=(const DCMpointAB& otherPoint) {
		if (this != &otherPoint) {
			static_cast<DCMpoint2D&>(*this) = otherPoint;
			copyDCMpointAB(otherPoint);
		}
		return *this;
	}

	void copyDCMpointAB(const DCMpointAB& otherPoint) {
	}

};


class DCMvertexAB : public DCMvertex2D {

public:

	DCMvertexAB() : DCMvertex2D() {
	}
	
	DCMvertexAB(const DCMvertexAB& otherVertex) : DCMvertex2D(otherVertex) {
		copyDCMvertexAB(otherVertex);
	}

	DCMvertexAB& operator=(const DCMvertexAB& otherVertex) {
		if (this != &otherVertex) {
			static_cast<DCMvertex2D&>(*this) = otherVertex;
			copyDCMvertexAB(otherVertex);			
		}
		return *this;
	}

	void copyDCMvertexAB(const DCMvertexAB& otherVertex) {
	}

	virtual ~DCMvertexAB() {}

};


class DCMedgeAB : public DCMedge2D {

public:

	double abnormal_line_coef;

	DCMedgeAB() : DCMedge2D() {
		abnormal_line_coef = 1.0;
	}

	DCMedgeAB(const DCMedgeAB& otherEdge) : DCMedge2D(otherEdge) {
		copyDCMedgeAB(otherEdge);
	}

	DCMedgeAB& operator=(const DCMedgeAB& otherEdge) {
		if (this != &otherEdge) {
			static_cast<DCMedge2D&>(*this) = otherEdge;
			copyDCMedgeAB(otherEdge);
		}
		return *this;
	}

	void copyDCMedgeAB(const DCMedgeAB& otherEdge) {
		abnormal_line_coef = otherEdge.abnormal_line_coef;
	}

	virtual ~DCMedgeAB() {}

	virtual void setLineCoefficient(int flag=0) {

		DCMedge2D::setLineCoefficient(flag);

		line_coef *= abnormal_line_coef;

		if (flag == 1) line_coef += PM.Sigma*randomness_of_line_coef;

	}

	
};


class DCMcellAB : public DCMcell2D<DCMvertexAB> {

public:

	int abnormal_flag;
	int abnormal_cluster;	
	static double abnormal_division_coef;
	static double normal_division_coef;	
	int neighbor_flag;
	double abnormal_perimeter_coef;
	
	int eli_flag;

	DCMcellAB() : DCMcell2D<DCMvertexAB>() {
		abnormal_flag = 0;
		abnormal_cluster = 0;
		neighbor_flag = 0;
		abnormal_perimeter_coef = 1.0;

		eli_flag = 0;		
	}

	DCMcellAB(const DCMcellAB& otherCell) : DCMcell2D<DCMvertexAB>(otherCell) {
		copyDCMcellAB(otherCell);
	}

	DCMcellAB& operator=(const DCMcellAB& otherCell) {
		if (this != &otherCell) {
			static_cast<DCMcell2D<DCMvertexAB>&>(*this) = otherCell;
			copyDCMcellAB(otherCell);
		}
		return *this;
	}

	void copyDCMcellAB(const DCMcellAB& otherCell) {
		abnormal_flag = otherCell.abnormal_flag;
		abnormal_cluster = otherCell.abnormal_cluster;
		neighbor_flag = otherCell.neighbor_flag;
		abnormal_perimeter_coef = otherCell.abnormal_perimeter_coef;

		eli_flag = otherCell.eli_flag;
	}

	virtual ~DCMcellAB() {}

	virtual void setCellCycleCoefficient(int flag=0) {

		DCMcell2D<DCMvertexAB>::setCellCycleCoefficient(flag);

		if (abnormal_flag == 0) cellcycle_update_coef *= normal_division_coef;
		else if (abnormal_flag == 1) cellcycle_update_coef *= abnormal_division_coef;

		if (neighbor_flag > 0) cellcycle_update_coef = 0;

		if (eli_flag > 1) cellcycle_update_coef = 0;
		
	}

	virtual void setPerimeterCoefficient(int flag=0) {

		DCMcell2D<DCMvertexAB>::setPerimeterCoefficient(flag);

		if (abnormal_flag == 1) perimeter_coef *= abnormal_perimeter_coef;

	}


};


#define TDCMcontainerAB DCMcontainer2D<DCMcellAB, DCMvertexAB, DCMedgeAB, DCMpointAB>

class DCMcontainerAB : public TDCMcontainerAB {

public:

	int abnormal_init_number;
	int abnormal_init_layer_real;
	int abnormal_con_flag;
	int abnormal_cells_size;
	int abnormal_cells_size_current;
	int abnormal_cells_size_change_flag;
	int abnormal_cells_size_flag;
	int abnormal_cells_out_size;
	int abnormal_cluster_size;
	int abnormal_start_time;
	int abnormal_wait_time_1;
	int abnormal_wait_time_2;
	double abnormal_line_coef;
	double abnormal_division_coef;
	double abnormal_perimeter_coef;
	double abnormal_all_area;
	double abnormal_all_area_onset;
	double finish_all_area;
	int normal_cells_size;
	int normal_cells_size_current;
	int normal_cells_size_flag;
	int finish_cells_size;
	
	DCMcontainerAB() : TDCMcontainerAB() {
		abnormal_init_number = -1;
		abnormal_init_layer_real = -1;
		abnormal_con_flag = 0;
		abnormal_cells_size = 0;
		abnormal_cells_size_current = 0;
		abnormal_cells_size_change_flag = 0;
		abnormal_cells_size_flag = 0;
		abnormal_cells_out_size = 0;
		abnormal_cluster_size = 0;
		abnormal_start_time = -1;
		abnormal_wait_time_1 = 1000;
		abnormal_wait_time_2 = 10000;
		abnormal_line_coef = 1.0;
		abnormal_division_coef = 1.0;
		abnormal_perimeter_coef = 1.0;
		abnormal_all_area = 0;
		abnormal_all_area_onset = 0;
		finish_all_area = 0;
		normal_cells_size = 0;
		normal_cells_size_current = 0;
		normal_cells_size_flag = 0;
		finish_cells_size = 0;
	}
	virtual ~DCMcontainerAB() {}

	virtual void destroy() {
		DCMcontainerAB::~DCMcontainerAB();
	}

	virtual void initialize(int flag=0) {
		
		TDCMcontainerAB::initialize(flag);

		vertex_dynamics_flag &= ~VD::T2TH_RELATIVE;

		std::ifstream ifs(L"DCModel_param_abnormal.txt");

		if (ifs) {
			std::string str;

			std::vector<std::string> param_list;
			param_list.push_back("Normal growth flag");	// 0
			param_list.push_back("Abnormal relative growth");	// 1
			param_list.push_back("Onset abnormal cluster size");	// 2
			param_list.push_back("Mechanism1 line tension coef");	// 3
			param_list.push_back("Wait time 1");	// 4
			param_list.push_back("Wait time 2");	// 5
			param_list.push_back("Abnormal initial index");	// 6
			param_list.push_back("Onset stop division");	// 7
			param_list.push_back("Boundary cell elimination flag");	// 8
			param_list.push_back("Abnormal perimeter coef");	// 9
			param_list.push_back("Simulation finish size");		// 10
			param_list.push_back("Abnormal initial layer real");		// 11
			param_list.push_back("Simulation finish area");		// 12
			
			while(std::getline(ifs, str)) {

				std::string token1;
				std::string token2;
				std::istringstream stream(str);

				std::getline(stream, token1, ':');
				std::getline(stream, token2);

				for (unsigned int i=0; i<param_list.size(); i++) {
					if (token1.compare(param_list[i]) == 0 && !token2.empty()) {
						switch (i) {
						case 0: if(std::stoi(token2)==0) DCMcellAB::normal_division_coef = 0; break;
						case 1: abnormal_division_coef = std::stod(token2); break;
						case 2: abnormal_cells_size_flag = std::stoi(token2); break;
						case 3: abnormal_line_coef = std::stod(token2); break;
						case 4: abnormal_wait_time_1 = std::stoi(token2); break;
						case 5: abnormal_wait_time_2 = std::stoi(token2); break;
						case 6: if(std::stoi(token2)>-1) abnormal_init_number = std::stoi(token2); break;
						case 7: if(std::stoi(token2)>0) abnormal_con_flag^=AC::ONSET_STOP_DIVISION; break;
						case 8: if(std::stoi(token2)>0) { abnormal_con_flag^=AC::BOUNDARY_CELL_ELI; normal_cells_size_flag = std::stoi(token2); } break;
						case 9: if(std::stod(token2)>0) { abnormal_con_flag^=AC::ABNORMAL_MECH_PARAM; abnormal_perimeter_coef = std::stod(token2); } break;
						case 10: if(std::stod(token2)>-1) { abnormal_con_flag^=AC::FINISH_PARAM; finish_cells_size = std::stoi(token2); } break;
						case 11: if(std::stoi(token2)>-1) abnormal_init_layer_real = std::stoi(token2); break;
						case 12: if(std::stod(token2)>0) { abnormal_con_flag^=AC::FINISH_PARAM; finish_all_area = std::stod(token2); } break;
						}
					}
				}

			}

		}

		if (abnormal_init_number < 0) {
			abnormal_con_flag ^= AC::RANDOM_ABNORMAL;
			int check = 0;
			do {				
				int ab_cellIndex = SF.uniformDistribution_int(0, cellsSize);
				if (DCMcells[ab_cellIndex]->ignoreFlag == 0 && DCMcells[ab_cellIndex]->outsideFlag == 0) {

					int oc_flag = 0;
					
					for (int j=0; j<DCMcells[ab_cellIndex]->verticesSize; j++) {
						if (DCMcells[ab_cellIndex]->DCMvertices[j]->ignoreFlag == 0) {
							int oc = DCMcells[ab_cellIndex]->DCMvertices[j]->ocIndex;
							if (oc > -1) {
								if (DCMcells[oc]->outsideFlag == 1) {
									oc_flag = 1;
									break;
								}
							}
						}
					}

					if (oc_flag == 0) {
						DCMcells[ab_cellIndex]->abnormal_flag = 1;
						abnormal_init_number = ab_cellIndex;
						check = 1;
					}
				}
			} while (check == 0);			
		}
		else {
			if (abnormal_init_layer_real < 0) {
				DCMcells[abnormal_init_number]->abnormal_flag = 1;
				DCMcells[abnormal_init_number]->cellCycle = PM.CellCycleTime - 0.01;
			}
			else {
				setAbnormalInitPopulation();
			}
		}

		DCMcellAB::abnormal_division_coef = abnormal_division_coef;

		setCellNeighborFlag();
				
	}

	void initialize_cellpopulation(double ab_line_coef, int ab_cells_size_flag) {
		abnormal_line_coef += ab_line_coef;
		abnormal_cells_size_flag += ab_cells_size_flag;
	}
	

	virtual void copy_Daughtercell(int mother, int daughter) {

		TDCMcontainerAB::copy_Daughtercell(mother, daughter);
		
		DCMcells[daughter]->abnormal_flag = DCMcells[mother]->abnormal_flag;

	}


	void abnormal_Regulation_on() {
		
		int start_flag = 0;

		if ((abnormal_con_flag & AC::LINETENSION_ON) == 0) {
			abnormal_con_flag ^= AC::LINETENSION_ON;
			start_flag = 1;
		}

		if (start_flag > 0) {
			DCMcellAB::abnormal_division_coef = abnormal_division_coef;
			updateTopology();
			abnormal_start_time = conTimeStep/1000;
			abnormal_all_area_onset = abnormal_all_area;
		}

	}

	void setCellNeighborFlag() {
		
		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				if (DCMcells[i]->abnormal_flag == 0) {

					int flag = 0;					
					for (int j=0; j<DCMcells[i]->verticesSize; j++) {
						if (DCMcells[i]->DCMvertices[j]->ignoreFlag == 0) {
							int oc = DCMcells[i]->DCMvertices[j]->ocIndex;
							if (oc > -1) {
								if (DCMcells[oc]->abnormal_flag > 0) {
									flag = 1;
									break;
								}
							}
						}
					}
					
					if (flag == 0) {
						DCMcells[i]->neighbor_flag = 0;
					}
					else {
						DCMcells[i]->neighbor_flag = 1;
					}						
				}					
			}
		}

	}

	void setAbnormalInitPopulation() {
		
		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				DCMcells[i]->abnormal_cluster = 0;
			}
		}

		std::vector<int> abnormal_pop;
		std::vector<int> curr_layer_cells;
		std::vector<int> curr_layer_cells_next;
		abnormal_pop.push_back(abnormal_init_number);		
		DCMcells[abnormal_init_number]->abnormal_flag = 1;
		DCMcells[abnormal_init_number]->abnormal_cluster = 1;
		curr_layer_cells_next.push_back(abnormal_init_number);

		do {

			if (curr_layer_cells.empty()) {
				for (int i=0; i<curr_layer_cells_next.size(); i++) {
					int cell_index = curr_layer_cells_next[i];
					for (int j=0; j<DCMcells[cell_index]->verticesSize; j++) {
						if (DCMcells[cell_index]->DCMvertices[j]->ignoreFlag == 0) {
							int oc = DCMcells[cell_index]->DCMvertices[j]->ocIndex;
							if (oc > -1) {
								if (DCMcells[oc]->outsideFlag == 0 && DCMcells[oc]->abnormal_cluster == 0) {
									DCMcells[oc]->abnormal_cluster = 1;
									curr_layer_cells.push_back(oc);									
								}
							}
						}						
					}
				}
				curr_layer_cells_next.clear();
			}

			int focal_cell = curr_layer_cells.back();
			curr_layer_cells.pop_back();
			curr_layer_cells_next.push_back(focal_cell);
			abnormal_pop.push_back(focal_cell);
			DCMcells[focal_cell]->abnormal_flag = 1;
			
		} while (abnormal_pop.size() < abnormal_init_layer_real);

		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				DCMcells[i]->abnormal_cluster = 0;
			}
		}
	}

	void setAbnormalCluster() {

		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				DCMcells[i]->abnormal_cluster = 0;
			}
		}
		abnormal_cluster_size = 0;

		for (int i=0; i<cellsSize; i++) {

			if (DCMcells[i]->ignoreFlag == 0 && DCMcells[i]->abnormal_flag > 0) {

				if (DCMcells[i]->abnormal_cluster == 0) {

					abnormal_cluster_size++;
					
					DCMcells[i]->abnormal_cluster = abnormal_cluster_size;
					int focal_cell = i;

					std::vector<int> stack;
					stack.push_back(focal_cell);

					do {
						int check = 0;
						for (int j=0; j<DCMcells[focal_cell]->verticesSize; j++) {
							if (DCMcells[focal_cell]->DCMvertices[j]->ignoreFlag == 0) {
								int oc = DCMcells[focal_cell]->DCMvertices[j]->ocIndex;
								if (oc > -1) {
									if (DCMcells[oc]->abnormal_flag > 0) {
										if (DCMcells[oc]->abnormal_cluster == 0) {
											DCMcells[oc]->abnormal_cluster = abnormal_cluster_size;
											focal_cell = oc;
											stack.push_back(focal_cell);
											check = 1;
											break;
										}
									}
								}
							}
						}
						if (check == 0) {
							stack.pop_back();
							if (!stack.empty()) focal_cell = stack.back();
						}
					} while(!stack.empty());
					
				}
			}
		}

	}

	virtual void setCoefficient(int flag=0) {

		flag = 1;
		TDCMcontainerAB::setCoefficient(flag);

	}

	virtual void updateCellDivision(int flag=0) {
		
		if ((abnormal_con_flag & AC::ONSET_STOP_DIVISION) == 0 || (abnormal_con_flag & AC::LINETENSION_ON) == 0) {
				TDCMcontainerAB::updateCellDivision(flag);
		}

	}

	virtual void updateTopology() {

		TDCMcontainerAB::updateTopology();

		setCellNeighborFlag();

		if ((abnormal_con_flag & AC::LINETENSION_ON) == AC::LINETENSION_ON) {

			for (int i=0; i<edgesSize; i++) {
			
				int c0, c1;
			
				c0 = edges[i]->cell[0];
				c1 = edges[i]->cell[1];
			
				if (c1 > -1) {
					if (DCMcells[c0]->abnormal_flag != DCMcells[c1]->abnormal_flag) {
						edges[i]->abnormal_line_coef *= abnormal_line_coef;
					}
				}
			}
		}

		if ((abnormal_con_flag & AC::BOUNDARY_CELL_ELI) == AC::BOUNDARY_CELL_ELI) {
			for (int i=0; i<edgesSize; i++) {
				if (edges[i]->cellNum < 2) {
					edges[i]->t1_coefficient *= 50.0;
				}
			}
		}

	}

	virtual int vDynamicsMain() {

		abnormal_cells_size = 0;
		abnormal_all_area = 0;
		
		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				if (DCMcells[i]->abnormal_flag == 1) {
					abnormal_cells_size++;
					abnormal_all_area += DCMcells[i]->cellArea[0];
					if ((abnormal_con_flag & AC::ABNORMAL_MECH_PARAM) == AC::ABNORMAL_MECH_PARAM) DCMcells[i]->abnormal_perimeter_coef = abnormal_perimeter_coef;
					if (DCMcells[i]->outsideFlag == 1) {
						abnormal_cells_out_size++;
					}
				}
			}
		}
		
		if (abnormal_cells_size > abnormal_cells_size_flag-1) {
				abnormal_Regulation_on();
		}
		
		if ((vertex_dynamics_flag & VD::OUT_FIX) == VD::OUT_FIX) {
			for (int i=0; i<cellsSize; i++) {
				if (DCMcells[i]->ignoreFlag == 0) {
					if (DCMcells[i]->outsideFlag > 0) {
						DCMcells[i]->fix_flag = 1;
					}
				}
			}
			for (int i=0; i<verticesSize; i++) {
				if (vertices[i]->vertexNum < 3) {
					for (int j=0; j<3; j++) {
						int e = vertices[i]->edgeIndex[2*j];
						if (e > -1) {
							edges[e]->fix_flag = 1;
						}
					}
				}
			}
		}

		return TDCMcontainerAB::vDynamicsMain();
	}

	virtual int finishContainer(int flag=0) {

		if ((abnormal_con_flag & AC::FINISH_PARAM) == 0) {
			if (abnormal_cells_size > 10000) return 1;
			if (abnormal_cells_size > abnormal_cells_size_flag*2 && abnormal_all_area > abnormal_all_area_onset) return 1;
		}
		else {
			if (finish_cells_size > 0) { 
				if (abnormal_cells_size > abnormal_cells_size_flag*finish_cells_size) return 1;
				if (abnormal_cells_out_size > 0) return 1;
			}
			else {
				if (finish_all_area > 0) { if (abnormal_all_area_onset > 0) if (abnormal_all_area > abnormal_all_area_onset*finish_all_area) return 1; }
				else { if (cellsSize_r == abnormal_cells_size) return 1; }
			}
		}
		if (abnormal_cells_size < 1) return 2;
		if (abnormal_cells_size_change_flag > 0 && abnormal_cells_size < 10) return 2;
		if (abnormal_wait_time_1 > 0) {
			if ((abnormal_con_flag & AC::FINISH_PARAM) == 0) {
				if (finish_cells_size > 0) { if ( (abnormal_start_time >= 0) && ((conTimeStep/1000-abnormal_start_time) > abnormal_wait_time_1) ) return 3;	}
			}
			else {
				if (finish_all_area > 0) { if ( (abnormal_start_time >= 0) && ((conTimeStep/1000-abnormal_start_time) > abnormal_wait_time_1) ) return 3; }
			}

			
		}
		if (abnormal_wait_time_2 > 0) {
			if ( (abnormal_start_time < 0) && ((conTimeStep/1000) > abnormal_wait_time_2) ) return 11;
		}
		
		return 0;
	}

	virtual int checkContainer(int flag=0) {

		if (flag == 1) {
			int cs = 0;

			for (int i=0; i<cellsSize; i++) {
				if (DCMcells[i]->ignoreFlag == 0) {
					if (DCMcells[i]->outsideFlag == 0) cs++;
				}
			}

			if (cell_size_counter != cs) {
				cell_size_counter = cs;
				return 1;
			}

			if ((abnormal_con_flag & AC::BOUNDARY_CELL_ELI) == AC::BOUNDARY_CELL_ELI) {
				if (changeNormalCellSize() == 1) {
					checkOutCells();
				}
				outCellElimination();
			}
		}

		if (flag == 2) {
			if ((abnormal_con_flag & AC::LINETENSION_ON) == AC::LINETENSION_ON) {
				return 1;
			}			
		}

		if (flag == 3) {
			if (conTimeStep > 0 && conTimeStep%100000 == 0) abnormal_cells_size_change_flag = 1;
			return changeAbnormalCellSize();
		}

		return 0;
	}

	int changeAbnormalCellSize() {

		if (abnormal_cells_size_current != abnormal_cells_size) {
			abnormal_cells_size_current = abnormal_cells_size;
			abnormal_cells_size_change_flag = 0;
			return 1;
		}
		
		return 0;		

	}

	int changeNormalCellSize() {

		normal_cells_size = 0;

		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				if (DCMcells[i]->abnormal_flag == 0) {
					if (DCMcells[i]->out_cell_eli_flag == 0) {
						normal_cells_size++;
					}
				}
			}
		}

		if (normal_cells_size_current != normal_cells_size) normal_cells_size_current = normal_cells_size;

		if (normal_cells_size > normal_cells_size_flag) return 1;
		
		return 0;
	}

	void outCellElimination() {

		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				if (DCMcells[i]->out_cell_eli_flag > 0) {
					DCMcells[i]->setTargetArea(0.0);					
				}				
			}
		}
	}

	void forcedNeighborCellDivision() {

		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				if (DCMcells[i]->neighbor_flag > 0) {
					DCMcells[i]->cellCycle = PM.CellCycleTime*1.01;
					break;
				}
			}
		}		
	}

};



#endif //DCMODELSTRUCTABNORMAL_HH