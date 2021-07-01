#include "DCMdata_test_2D.h"

const double rt3 = 1.732050807568877;
const double av = 1.0/sqrt(8.3138438763306102);

void DCMInit(DCMcontainer& container) {
		
	FILE *diagvert;
	FILE *diagval;
	
	if (container.initDataFlag == 0) {
		wchar_t init_c[50];
		wchar_t init_p[50];
		swprintf_s(init_c, 50, L"initial2D_RD05_in%d_c_L%04dG%04d.csv", PM.initial_cell_size, static_cast<int>(PM.Sigma*1000.0), static_cast<int>(PM.Line_Coef*1000.0));
		swprintf_s(init_p, 50, L"initial2D_RD05_in%d_p_L%04dG%04d.csv", PM.initial_cell_size, static_cast<int>(PM.Sigma*1000.0), static_cast<int>(PM.Line_Coef*1000.0));
		_wfopen_s(&diagvert, init_c, L"r");
		_wfopen_s(&diagval, init_p, L"r");
	}	

	double vx, vy;

	struct tVertexCor {
		double x, y;
	};
	
	tVertexCor tempVC;
	std::vector<tVertexCor> vertexCor;

	while (EOF != fwscanf_s(diagvert, L"%lf, %lf\n", &vx, &vy)) {
		tempVC.x= vx; tempVC.y = vy;
		vertexCor.push_back(tempVC);		
	}

	int vi;

	struct tCell {
		std::vector<int> tVertex;
	};

	tCell tempC;
	std::vector<tCell> cel;
	
	while (EOF != fwscanf_s(diagval, L"%d, ", &vi)) {

		int check = 0;

		while (check == 0) {
			if (vi == 0) {
				check = 1;
				cel.push_back(tempC);
				tempC.tVertex.clear();
			}
			else {
				tempC.tVertex.push_back(vi);
				fwscanf_s(diagval, L"%d, ", &vi);				
			}
		}		
	}

	fclose(diagvert);
	fclose(diagval);
		
	
	int cSize = (int)cel.size();

	
	double area = 0;

	for (int i=0; i<cSize; i++) {

		int vSize = (int)cel[i].tVertex.size();

		for (int j=1; j<vSize-1; j++) {

			int p0 = cel[i].tVertex[0]-1;
			int p1 = cel[i].tVertex[j]-1;
			int p2 = cel[i].tVertex[j+1]-1;

			area += 0.5*((vertexCor[p0].x*vertexCor[p1].y+vertexCor[p1].x*vertexCor[p2].y+vertexCor[p2].x*vertexCor[p0].y)
				-(vertexCor[p1].x*vertexCor[p0].y+vertexCor[p2].x*vertexCor[p1].y+vertexCor[p0].x*vertexCor[p2].y));			
		}
	}

	area /= cSize;
	area = sqrt(area);
	area = 1.0f/area;
			
	container.createDCMcells(cSize);

	for (int i=0; i<cSize; i++) {

		int vSize = (int)cel[i].tVertex.size();

		container.createDCMvertices(i, vSize);
		
		for (int j=0; j<vSize; j++) {

			int pi = cel[i].tVertex[j]-1;
			int p = j-1;
			int n = j+1;
			if (j == 0) p = vSize-1;
			if (j == vSize-1) n = 0;			

			container.setNewVertex(i, j, vertexCor[pi].x, vertexCor[pi].y, 0, p, n);

		}
	}

}
