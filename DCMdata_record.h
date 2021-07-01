#ifndef DCMDATARECORD_HH
#define DCMDATARECORD_HH


#include <wchar.h>
#include "DCModel_struct_2D.h"
#include "DCModel_struct_abn.h"


template <typename Tcontainer2D>
class DCMdata {

protected:

	FILE* data;
	FILE* data1;
	int dataIndex;
	int dataIndex1;
	int dataIndex2;
	Tcontainer2D* con;
	int containerIndex;
	int used_flag;

public:

	DCMdata() : con(NULL) {}
	virtual ~DCMdata() {}

	void openFile(int read=0, int file_num=0) {
		wchar_t cIndex[50];		
		if (file_num == 0) swprintf_s(cIndex, 50, L"data2D_%02d_%02d.csv", dataIndex, containerIndex);
		else if (file_num == 1) swprintf_s(cIndex, 50, L"data2D_%02d_%02d.csv", dataIndex1, containerIndex);
		else if (file_num == 2) swprintf_s(cIndex, 50, L"data2D_%02d_%02d.csv", dataIndex2, containerIndex);
		else if (file_num == 3) swprintf_s(cIndex, 50, L"data2D_%02d.csv", dataIndex);

		if (read == 0) _wfopen_s(&data, cIndex, L"a");
		else if (read == 1) _wfopen_s(&data, cIndex, L"w");			
	}

	void openFiles(int num, int read=0, int file_num=0) {
		wchar_t cIndex[50];
		if (file_num == 0) swprintf_s(cIndex, 50, L"data2D_%02d_%02d_%02d.csv", dataIndex, containerIndex, num);
		else if (file_num == 1) swprintf_s(cIndex, 50, L"data2D_%02d_%02d_%02d.csv", dataIndex1, containerIndex, num);
		else if (file_num == 2) swprintf_s(cIndex, 50, L"data2D_%02d_%02d_%02d.csv", dataIndex2, containerIndex, num);
		
		if (read == 0) _wfopen_s(&data, cIndex, L"a");
		else if (read == 1) _wfopen_s(&data, cIndex, L"w");
	}

	void openFileab(int ab, int read=0) {
		wchar_t cIndex[50];		
		if (ab == 0) swprintf_s(cIndex, 50, L"data2D_%02da_%02d.csv", dataIndex, containerIndex);
		else if (ab == 1) swprintf_s(cIndex, 50, L"data2D_%02db_%02d.csv", dataIndex, containerIndex);
		
		if (read == 0) _wfopen_s(&data, cIndex, L"a");
		else if (read == 1) _wfopen_s(&data, cIndex, L"w");
	}
	
	void closeFile() {
		fclose(data);
		used_flag = dataIndex;
	}

	void setContainer(Tcontainer2D& container, int index) {		
		con = &container;
		containerIndex = index;
		used_flag = 0;
	}

	int usedCheck() {
		return used_flag;
	}

	void resetFile() {
		if (used_flag > 0) {
			openFile(1);
			fclose(data);
		}
	}
		
};



class DCMdata_Record_AbnormalPopulation : public DCMdata<DCMcontainerAB> {

public:

	DCMdata_Record_AbnormalPopulation() : DCMdata<DCMcontainerAB>() {
		dataIndex = 54;
	}

	virtual ~DCMdata_Record_AbnormalPopulation() {}

	void operator() (unsigned __int64 timeStep) {

		wchar_t fn[50];
		swprintf_s(fn, 50, L"data2D_N%03dmu%03d.csv", con->abnormal_cells_size_flag, static_cast<int>(con->abnormal_line_coef*100.0));		
		_wfopen_s(&data, fn, L"a");
		
		int normal_cell = 0;
		int abnormal_cell = 0;
		int normal_cell_eli = 0;
		int abnormal_cell_eli = 0;
		int normal_cell_ex = 0;
		int abnormal_cell_ex = 0;
		int abnormal_in_cell = 0;
		int abnormal_out_cell = 0;
		int abnormal_edge_boundary = 0;
		int neighbor_cell = 0;
		double mean_area_abnormal_cell = 0;
		double mean_sq_area_abnormal_cell = 0;
		double mean_area_abnormal_in_cell = 0;		
		double mean_sq_area_abnormal_in_cell = 0;
		double mean_area_abnormal_out_cell = 0;
		double mean_sq_area_abnormal_out_cell = 0;
		double mean_area_neighbor_cell = 0;
		double mean_sq_area_neighbor_cell = 0;
		double mean_stmag_abnormal_cell = 0;
		double mean_stmag_abnormal_in_cell = 0;
		double mean_stmag_abnormal_out_cell = 0;
		double mean_stmag_neighbor_cell = 0;
		double mean_peri_abnormal_cell = 0;
		double mean_side_abnormal_cell = 0;
		double mean_coef_abnormal_cell = 0;
		double mean_peri_abnormal_in_cell = 0;
		double mean_side_abnormal_in_cell = 0;
		double mean_coef_abnormal_in_cell = 0;
		double mean_peri_abnormal_out_cell = 0;
		double mean_side_abnormal_out_cell = 0;
		double mean_coef_abnormal_out_cell = 0;
		double mean_peri_neighbor_cell = 0;
		double mean_side_neighbor_cell = 0;
		double mean_coef_neighbor_cell = 0;
		double poly_abnormal_cell[6] = { 0 };
		double poly_abnormal_in_cell[6] = { 0 };
		double poly_abnormal_out_cell[6] = { 0 };
		double poly_area_abnormal_cell[6] = { 0 };
		double poly_area_abnormal_in_cell[6] = { 0 };
		double poly_area_abnormal_out_cell[6] = { 0 };
		
		for (int i=0; i<con->cellsSize; i++) {

			if (con->DCMcells[i]->ignoreFlag == 0) {

				double area_cell = con->DCMcells[i]->cellArea[0];		
				double stmag_cell = con->DCMcells[i]->cellStressMatrix[4]+con->DCMcells[i]->cellStressMatrix[5];
				double peri_cell = con->DCMcells[i]->perimeterLength;
				double side_cell = static_cast<double>(con->DCMcells[i]->verticesSize_r);
				double coef_cell = peri_cell/side_cell/sqrt(area_cell);
		
				if (con->DCMcells[i]->abnormal_flag == 0) {

					normal_cell++;
					if (con->DCMcells[i]->outsideFlag != 0) normal_cell_ex++;

					if (con->DCMcells[i]->neighbor_flag > 0) {

						neighbor_cell++;

						mean_area_neighbor_cell += area_cell;
						mean_sq_area_neighbor_cell += area_cell*area_cell;
						mean_stmag_neighbor_cell += stmag_cell;
						mean_peri_neighbor_cell += peri_cell;
						mean_side_neighbor_cell += side_cell;
						mean_coef_neighbor_cell += coef_cell;
						
					}

				}
				else {

					abnormal_cell++;
					mean_area_abnormal_cell += area_cell;
					mean_sq_area_abnormal_cell += area_cell*area_cell;
					mean_stmag_abnormal_cell += stmag_cell;
					mean_peri_abnormal_cell += peri_cell;
					mean_side_abnormal_cell += side_cell;
					mean_coef_abnormal_cell += coef_cell;

					if (con->DCMcells[i]->outsideFlag != 0) abnormal_cell_ex++;

					int v = con->DCMcells[i]->verticesSize_r;

					if (v > 3 && v < 9) {
						poly_abnormal_cell[v-4] += 1.0;
						poly_area_abnormal_cell[v-4] += area_cell;
					}
					else if (v > 8) {
						poly_abnormal_cell[5] += 1.0;
						poly_area_abnormal_cell[5] += area_cell;
					}

					int abnormal_edge_boundary_count = 0;
					for (int j=0; j<con->DCMcells[i]->verticesSize; j++) {
						if (con->DCMcells[i]->DCMvertices[j]->ignoreFlag == 0) {
							int oc = con->DCMcells[i]->DCMvertices[j]->ocIndex;
							if (oc > -1) {
								if (con->DCMcells[oc]->abnormal_flag == 0) abnormal_edge_boundary_count++;
							}
						}
					}

					if (abnormal_edge_boundary_count > 0) {

						abnormal_out_cell++;
						mean_area_abnormal_out_cell += area_cell;
						mean_sq_area_abnormal_out_cell += area_cell*area_cell;
						mean_stmag_abnormal_out_cell += stmag_cell;
						mean_peri_abnormal_out_cell += peri_cell;
						mean_side_abnormal_out_cell += side_cell;
						mean_coef_abnormal_out_cell += coef_cell;

						if (v > 3 && v < 9) {
							poly_abnormal_out_cell[v-4] += 1.0;
							poly_area_abnormal_out_cell[v-4] += area_cell;
						}
						else if (v > 8) {
							poly_abnormal_out_cell[5] += 1.0;
							poly_area_abnormal_out_cell[5] += area_cell;
						}
					}
					else {

						abnormal_in_cell++;
						mean_area_abnormal_in_cell += area_cell;
						mean_sq_area_abnormal_in_cell += area_cell*area_cell;
						mean_stmag_abnormal_in_cell += stmag_cell;
						mean_peri_abnormal_in_cell += peri_cell;
						mean_side_abnormal_in_cell += side_cell;
						mean_coef_abnormal_in_cell += coef_cell;

						if (v > 3 && v < 9) {
							poly_abnormal_in_cell[v-4] += 1.0;
							poly_area_abnormal_in_cell[v-4] += area_cell;
						}
						else if (v > 8) {
							poly_abnormal_in_cell[5] += 1.0;
							poly_area_abnormal_in_cell[5] += area_cell;
						}
					}

					abnormal_edge_boundary += abnormal_edge_boundary_count;
				}
			}
			else {
				if (con->DCMcells[i]->abnormal_flag == 0) {
					normal_cell_eli++;
				}
				else {
					abnormal_cell_eli++;
				}
			}
		}

		if (abnormal_cell > 0) {
			mean_area_abnormal_cell /= static_cast<double>(abnormal_cell);
			mean_sq_area_abnormal_cell /= static_cast<double>(abnormal_cell);
			mean_stmag_abnormal_cell /= static_cast<double>(abnormal_cell);
			mean_peri_abnormal_cell /= static_cast<double>(abnormal_cell);
			mean_side_abnormal_cell /= static_cast<double>(abnormal_cell);
			mean_coef_abnormal_cell /= static_cast<double>(abnormal_cell);
			for (int i=0; i<6; i++) {
				if (poly_abnormal_cell[i] > 0) poly_area_abnormal_cell[i] /= poly_abnormal_cell[i];
				poly_abnormal_cell[i] /= static_cast<double>(abnormal_cell);				
			}
		}

		if (abnormal_out_cell > 0) {
			mean_area_abnormal_out_cell /= static_cast<double>(abnormal_out_cell);
			mean_sq_area_abnormal_out_cell /= static_cast<double>(abnormal_out_cell);
			mean_stmag_abnormal_out_cell /= static_cast<double>(abnormal_out_cell);
			mean_peri_abnormal_out_cell /= static_cast<double>(abnormal_out_cell);
			mean_side_abnormal_out_cell /= static_cast<double>(abnormal_out_cell);
			mean_coef_abnormal_out_cell /= static_cast<double>(abnormal_out_cell);
			for (int i=0; i<6; i++) {
				if (poly_abnormal_out_cell[i] > 0) poly_area_abnormal_out_cell[i] /= poly_abnormal_out_cell[i];
				poly_abnormal_out_cell[i] /= static_cast<double>(abnormal_out_cell);
			}				
		}

		if (abnormal_in_cell > 0) {
			mean_area_abnormal_in_cell /= static_cast<double>(abnormal_in_cell);
			mean_sq_area_abnormal_in_cell /= static_cast<double>(abnormal_in_cell);
			mean_stmag_abnormal_in_cell /= static_cast<double>(abnormal_in_cell);
			mean_peri_abnormal_in_cell /= static_cast<double>(abnormal_in_cell);
			mean_side_abnormal_in_cell /= static_cast<double>(abnormal_in_cell);
			mean_coef_abnormal_in_cell /= static_cast<double>(abnormal_in_cell);
			for (int i=0; i<6; i++) { 
				if (poly_abnormal_in_cell[i] > 0) poly_area_abnormal_in_cell[i] /= poly_abnormal_in_cell[i];
				poly_abnormal_in_cell[i] /= static_cast<double>(abnormal_in_cell);
			}
		}

		if (neighbor_cell > 0) {
			mean_area_neighbor_cell /= static_cast<double>(neighbor_cell);
			mean_sq_area_neighbor_cell /= static_cast<double>(neighbor_cell);
			mean_stmag_neighbor_cell /= static_cast<double>(neighbor_cell);
			mean_peri_neighbor_cell /= static_cast<double>(neighbor_cell);
			mean_side_neighbor_cell /= static_cast<double>(neighbor_cell);
			mean_coef_neighbor_cell /= static_cast<double>(neighbor_cell);
		}

		
		fwprintf(data, L"%d,NC,%d,%d",
			timeStep,
			normal_cell,
			con->recon_number_r
			);

		fwprintf(data, L",AC,%d,%d,%d,%d,%d,%d",
			con->abnormal_init_number,
			abnormal_cell,
			abnormal_out_cell,
			abnormal_in_cell,
			abnormal_cell_eli,
			abnormal_edge_boundary
			);

		fwprintf(data, L",A,%.32f,%.32f,%.32f,%.32f",
			mean_area_abnormal_cell*static_cast<double>(abnormal_cell),
			mean_area_abnormal_cell,			
			mean_area_abnormal_out_cell,
			mean_area_abnormal_in_cell
			);

		fwprintf(data, L",ALL");
		for (int i=0; i<6; i++) fwprintf(data, L",%.32f", poly_abnormal_cell[i]);
		fwprintf(data, L",ALLA");
		for (int i=0; i<6; i++) fwprintf(data, L",%.32f", poly_area_abnormal_cell[i]);

		fwprintf(data, L",OUT");
		for (int i=0; i<6; i++) fwprintf(data, L",%.32f", poly_abnormal_out_cell[i]);
		fwprintf(data, L",OUTA");
		for (int i=0; i<6; i++) fwprintf(data, L",%.32f", poly_area_abnormal_out_cell[i]);

		fwprintf(data, L",IN");
		for (int i=0; i<6; i++) fwprintf(data, L",%.32f", poly_abnormal_in_cell[i]);
		fwprintf(data, L",INA");
		for (int i=0; i<6; i++) fwprintf(data, L",%.32f", poly_area_abnormal_in_cell[i]);

		fwprintf(data, L",NEI,%d,%.32f,%.32f,%.32f",
			neighbor_cell,
			mean_area_neighbor_cell,
			mean_sq_area_neighbor_cell,
			mean_stmag_neighbor_cell
			);

		fwprintf(data, L",ST,%.32f,%.32f,%.32f,%.32f,%.32f,%.32f",
			mean_stmag_abnormal_cell,
			mean_stmag_abnormal_out_cell,
			mean_stmag_abnormal_in_cell,
			mean_sq_area_abnormal_cell,
			mean_sq_area_abnormal_out_cell,
			mean_sq_area_abnormal_in_cell
			);

		fwprintf(data, L",PS,%.32f,%.32f,%.32f,%.32f,%.32f,%.32f,%.32f,%.32f,%.32f,%.32f,%.32f,%.32f",
			mean_peri_abnormal_cell,
			mean_side_abnormal_cell,
			mean_coef_abnormal_cell,
			mean_peri_abnormal_out_cell,
			mean_side_abnormal_out_cell,
			mean_coef_abnormal_out_cell,
			mean_peri_abnormal_in_cell,
			mean_side_abnormal_in_cell,
			mean_coef_abnormal_in_cell,
			mean_peri_neighbor_cell,
			mean_side_neighbor_cell,
			mean_coef_neighbor_cell
			);

		con->setAbnormalCluster();
		fwprintf(data, L",CLS,%d", con->abnormal_cluster_size);

		double abnormal_cx = 0;
		double abnormal_cy = 0;
		
		for (int i=0; i<con->cellsSize; i++) {
			if (con->DCMcells[i]->ignoreFlag == 0) {
				if (con->DCMcells[i]->abnormal_flag > 0) {
					abnormal_cx += con->DCMcells[i]->cofm_x;
					abnormal_cy += con->DCMcells[i]->cofm_y;
				}
			}
		}
		abnormal_cx /= static_cast<double>(abnormal_cell);
		abnormal_cy /= static_cast<double>(abnormal_cell);

		double abnormal_mean_dist = 0;
		double abnormal_std_dist = 0;

		for (int i=0; i<con->cellsSize; i++) {
			if (con->DCMcells[i]->ignoreFlag == 0) {
				if (con->DCMcells[i]->abnormal_flag > 0) {
					double d = dist(con->DCMcells[i]->cofm_x, abnormal_cx, con->DCMcells[i]->cofm_y, abnormal_cy,0,0);
					abnormal_mean_dist += d;
					abnormal_std_dist += d*d;
				}
			}
		}

		abnormal_mean_dist /= static_cast<double>(abnormal_cell);
		abnormal_std_dist /= static_cast<double>(abnormal_cell);
		abnormal_std_dist -= abnormal_mean_dist*abnormal_mean_dist;
		
		fwprintf(data, L",CSTD,%.32f", sqrt(abnormal_std_dist));

		fwprintf(data, L"\n");

		closeFile();

	}


};


template <typename Tcontainer2D>
class DCMdata_Func {

public:

	DCMdata_Func() {}
	virtual ~DCMdata_Func() {}
	
	virtual void setContainer(Tcontainer2D& con, int index) {}

	virtual void resetFile() {}

};


class DCMdata_FuncAB : public DCMdata_Func<DCMcontainerAB> {

public:

	DCMdata_Record_AbnormalPopulation Record_AbnormalPopulation;
	
	DCMdata_FuncAB() : DCMdata_Func<DCMcontainerAB>() {}
	virtual ~DCMdata_FuncAB() {}

	virtual void setContainer(DCMcontainerAB& con, int index) {

		DCMdata_Func<DCMcontainerAB>::setContainer(con, index);
		Record_AbnormalPopulation.setContainer(con, index);

	}

	virtual void resetFile() {
		
		DCMdata_Func<DCMcontainerAB>::resetFile();
		Record_AbnormalPopulation.resetFile();
		
	}

};


#endif //DCMDATARECORD_HH