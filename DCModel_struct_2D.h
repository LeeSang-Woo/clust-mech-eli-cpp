#ifndef DCMODELSTRUCT2D_HH
#define DCMODELSTRUCT2D_HH

#include "DCModel_struct.h"
#include <tbb\tbb.h>
using namespace tbb;

////// Enum Flags /////////////////////////////////////////

namespace US {
	enum UpdateMainStatusFlags {

		INITIALIZE			 = 1<<1,
		RECONNECTION		 = 1<<2,
		EXTRUSION			 = 1<<3,
		DIVISION_NEI		 = 1<<4,
		DIVISION_OWN		 = 1<<5,

	};
}

namespace VD {
	enum VertexDynamicsStatusFlags {

		LINETENSION		= 1<<1,
		LINESPRING		= 1<<2,
		LINESPRING_NL	= 1<<3,
		LINEELASTICITY	= 1<<4,
		PERIMETER		= 1<<5,
		EDGEDISTANCE	= 1<<6,
		ANGLESPRING		= 1<<7,

		OUT_FIX			= 1<<10,
		RANDOMMOVE		= 1<<11,
		T2TH_RELATIVE	= 1<<12,
		RANDOMLAMBDA	= 1<<13,
		RANDLAMBDA_INIT	= 1<<14,
	};
}


namespace BC {
	enum BasicConditionFlags {

		MPHASE_2SI			= 1<<1,	// Relative size, immediately
		MPHASE_2SD			= 1<<2,	// Relative size, delayed
		MPHASE_FSI			= 1<<3,	// Fixed size, immediately
		MPHASE_FSD			= 1<<4,	// Fixed size, delayed
		MPHASE_TTH			= 1<<5,	// Time threshold
		MPHASE_ADD			= 1<<6, // Adder

		BOUNDARY_ELI		= 1<<25,
		NO_DIVISION			= 1<<26,
		MOVING_VERTEX		= 1<<20,
	};
}


namespace RE {
	enum RecordingFlags {

		TOPOLOGY_CHANGE		= 1<<1,
		STRESS_RECORD		= 1<<2,
		RECON_RECORD		= 1<<3,
		
	};
}


////// Structure //////////////////////////////////////////

class NeighborCell {

public:

	int index;
	int nn_level;
	
	NeighborCell() {
		index = -1;
		nn_level = 0;		
	}

};

class NeighborCell_Area : public NeighborCell {

public:
	double cell_area_be;
	double cell_area_af;

	NeighborCell_Area() : NeighborCell() {
		cell_area_be = 0;
		cell_area_af = 0;
	}
};


struct ReconnectionRecord {

	int cell_vertices_before[4];
	int cell_vertices_after[4];
	int time_step;
	double position;
};


class CellRecord {

public:

	int timeStep;
	int index;
	int verticesSize;
	int topologicalChange;
	int inout;
	int daughterCell;

	double area;
	double cell_shape_aniso;

	CellRecord() {
		area = 0.0;
		cell_shape_aniso = 0.0;
		index = -1;
		timeStep = -1;
		topologicalChange = -1;
		verticesSize = 0;
		inout = -1;
		daughterCell = -1;
	}

	virtual ~CellRecord() {}

};

class stressRecord : public CellRecord {

public:
	
	double mean_st_mag[2];
	double mean_st_ani[2];
	double own_st_mag[3];
	double own_st_ani[3];
	int neighborCell[30];
	NeighborCell_Area nnc[100];
	
	stressRecord() : CellRecord() {

		initArray(mean_st_mag, 0);
		initArray(mean_st_ani, 0);
		initArray(own_st_mag, 0);
		initArray(own_st_ani, 0);
		
		initArray(neighborCell, -1);
		
	}
	
	virtual ~stressRecord() {}
};




////// 2D Class Declaration ///////////////////////////////

#define TDCMcontainer2D DCMcontainer2D<DCMcell2D<DCMvertex2D>, DCMvertex2D, DCMedge2D, DCMpoint2D>


class DCMvertex2D : public DCMvertex {

public:
	
	int connection[2];
	int othercell[4];

	int verticesNumber;
	int edgesNumber[2];
	
	int ocIndex;
	int ovIndex;

	double edge_randomness_of_line_coef;

	//////////////////////////////

	DCMvertex2D() : DCMvertex() {
		initArray(connection, -1);
		initArray(othercell, -1);
		initArray(edgesNumber, -1);		
		verticesNumber = -1;
		ocIndex = -1;
		ovIndex = -1;
		edge_randomness_of_line_coef = 0;
	}

	DCMvertex2D(const DCMvertex2D& otherVertex) : DCMvertex(otherVertex) {
		copyDCMvertex2D(otherVertex);		
	}
	
	DCMvertex2D& operator=(const DCMvertex2D& otherVertex) {
		if (this != &otherVertex) {
			static_cast<DCMvertex&>(*this) = otherVertex;
			copyDCMvertex2D(otherVertex);
		}
		return *this;
	}
	
	void copyDCMvertex2D(const DCMvertex2D& otherVertex) {
		copyArray(otherVertex, connection);
		copyArray(otherVertex, othercell);
		copyArray(otherVertex, edgesNumber);		
		verticesNumber = otherVertex.verticesNumber;
		ocIndex = otherVertex.ocIndex;
		ovIndex = otherVertex.ovIndex;
		edge_randomness_of_line_coef = otherVertex.edge_randomness_of_line_coef;
	}

	virtual ~DCMvertex2D() {}


};


class DCMpoint2D : public DCMpoint {

public:

	int vertexNum;
	
	int cellIndex[3];

	int vertexIndex[3];

	int pointIndex[3];
			
	int edgeIndex[6];

	double dxdt;

	//////////////////////////////

	DCMpoint2D() : DCMpoint() {
		
		initArray(cellIndex, -1);
		initArray(vertexIndex, -1);
		initArray(pointIndex, -1);
		initArray(edgeIndex, -1);

		vertexNum = 0;

		dxdt = 0;
	}

	DCMpoint2D(const DCMpoint2D& otherPoint) : DCMpoint(otherPoint) {
		copyDCMpoint2D(otherPoint);
	}

	DCMpoint2D& operator=(const DCMpoint2D& otherPoint) {
		if (this != &otherPoint) {
			static_cast<DCMpoint&>(*this) = otherPoint;
			copyDCMpoint2D(otherPoint);
		}
		return *this;
	}

	void copyDCMpoint2D(const DCMpoint2D& otherPoint) {
		
		copyArray(otherPoint, cellIndex);
		copyArray(otherPoint, vertexIndex);
		copyArray(otherPoint, pointIndex);
		copyArray(otherPoint, edgeIndex);

		vertexNum = otherPoint.vertexNum;
		dxdt = otherPoint.dxdt;
	}

	virtual ~DCMpoint2D() {}

};


class DCMedge2D : public DCMedge {

public:
	
	//////////////////////////////

	static double ave_distance;

	int cell[2];

	int v1[2];
	int v2[2];

	double v1x, v1y, v1z, v2x, v2y, v2z;

	double line_coef;
	double line_spring_coef;
	double line_elasticity_coef;
	double line_elasticity_target;

	double t1_coefficient;
	int negative;

	double randomness_of_line_coef;
	int fix_flag;

	//////////////////////////////
		
	DCMedge2D() : DCMedge() {

		initArray(cell, -1);
		initArray(v1, -1);
		initArray(v2, -1);

		v1x = 0; v1y = 0; v1z = 0;
		v2x = 0; v2y = 0; v2z = 0;

		line_coef = 0;
		line_spring_coef = 0;
		line_elasticity_coef = 0;
		line_elasticity_target = 0;
		randomness_of_line_coef = 0;

		t1_coefficient = 1;
		negative = 1;
		fix_flag = 0;

	}

	DCMedge2D(const DCMedge2D& otherEdge) : DCMedge(otherEdge) {
		copyDCMedge2D(otherEdge);
	}

	DCMedge2D& operator=(const DCMedge2D& otherEdge) {
		if (this != &otherEdge) {
			static_cast<DCMedge&>(*this) = otherEdge;
			copyDCMedge2D(otherEdge);
		}
		return *this;
	}

	void copyDCMedge2D(const DCMedge2D& otherEdge) {
		copyArray(otherEdge, cell);
		copyArray(otherEdge, v1);
		copyArray(otherEdge, v2);
		v1x = otherEdge.v1x; v1y = otherEdge.v1y; v1z = otherEdge.v1z;
		v2x = otherEdge.v2x; v2y = otherEdge.v2y; v2z = otherEdge.v2z;
		line_coef = otherEdge.line_coef;
		line_spring_coef = otherEdge.line_spring_coef;
		line_elasticity_coef = otherEdge.line_elasticity_coef;
		line_elasticity_target = otherEdge.line_elasticity_target;
		t1_coefficient = otherEdge.t1_coefficient;
		negative = otherEdge.negative;
		randomness_of_line_coef = otherEdge.randomness_of_line_coef;
		fix_flag = otherEdge.fix_flag;
	}
	
	virtual ~DCMedge2D() {}

	//////////////////////////////

	void setDistance(double x1, double y1, double z1, double x2, double y2, double z2) {

		v1x = x1; v1y = y1, v1z = z1;
		v2x = x2; v2y = y2, v2z = z2;

		distance = dist(v1x, v2x, v1y, v2y, v1z, v2z);
	
	}

	virtual void setLineCoefficient(int flag=0) {
		
		line_coef = 1.0;

		if (flag == 0) line_coef += randomness_of_line_coef;

		line_coef *= PM.Sigma;
		
		if (cell[1] < 0) line_coef *= PM.Sigma_Zero;
		
	}
			
	int vDynamics_line(int eif, DCMpoint& xy) {
		
		double dx, dy, dz;

		dx = line_coef*(v2x-v1x)/distance;
		dy = line_coef*(v2y-v1y)/distance;
		dz = line_coef*(v2z-v1z)/distance;
		
		if (eif == 1) {
			dx = -dx;
			dy = -dy;
			dz = -dz;
		}
		
		xy.x += dx;
		xy.y += dy;
		xy.z += dz;

		return 0;

	}

	virtual void setLineSpringCoefficient(int oc_flag=0) {
		
		line_spring_coef = 1.0;
		line_spring_coef *= PM.Line_Coef;
		if (cell[1] < 0) line_spring_coef *= PM.Line_Coef_Out;

	}

	int vDynamics_linespring(int eif, DCMpoint& xy) {
				
		double dx, dy, dz;

		dx = line_spring_coef*(v2x-v1x);
		dy = line_spring_coef*(v2y-v1y);
		dz = line_spring_coef*(v2z-v1z);
		
		if (eif == 1) {
			dx = -dx;
			dy = -dy;
			dz = -dz;
		}
		
		xy.x += dx;
		xy.y += dy;
		xy.z += dz;
	
		return 0;

	}

	virtual void setLineElasticityCoefficient(int flag=0) {

		line_elasticity_coef = 1.0;
		line_elasticity_coef *= PM.Line_Beta;

		line_elasticity_target = PM.Line_Mu;

		if (flag == 1) {
			line_elasticity_coef *= 10.0;
		}

		// for case ellipse division (vein simulation?)
		if (cell[1] < 0) {
			if (flag != 2) {
				line_elasticity_coef *= 3.0;
				line_elasticity_target *= 0.2;
			}
		}
		//

	}

	int vDynamics_lineelasticity(int eif, DCMpoint& xy) {

		double dx, dy;

		dx = line_elasticity_coef*(v2x-v1x)/distance*(distance-line_elasticity_target);
		dy = line_elasticity_coef*(v2y-v1y)/distance*(distance-line_elasticity_target);
		
		if (eif == 1) {
			dx = -dx;
			dy = -dy;
		}

		xy.x += dx;
		xy.y += dy;

		return 0;

	}
		
};






template <typename Tvertex2D>
class DCMcell2D : public DCMcell {

public:

	//////////////////////////////

	Tvertex2D** DCMvertices;

	static const int initMaxVerticesSize = 30;	

	int statusFlag;

	int initialCellIndex;

	double perimeterLength;

	double cellArea[2];
	double cell_angle;

	double targetArea;

	double divArea;

	double cofm_x, cofm_y, cofm_z;

	double cellCycle;
	double cellcycle_update_coef;

	double divCycle;

	double line_coef;
	double line_out_coef;
	double perimeter_coef;
	double area_coef;
	double area_target;
	
	double cell_shape_aniso;

	double cellStressMatrix[10];

	double cell_stressmag_pre[3];
	double cell_stressani_pre[3];
	double cell_stressori_pre[3];

	double cellGrad[8];

	double cellStressMatrix_Bach[10];

	int rec_Size;
	int max_rec_Size;

	CellRecord* rec;
	CellRecord* rec_buf;
	CellRecord* rec_rect;

	int recordFlag;

	int lineage_Size;

	CellRecord* lineage;

	int topologyFlag;

	int out_cell_eli_flag;

	int fix_flag;

	//////////////////////////////
	
	DCMcell2D() {

		maxVerticesSize = initMaxVerticesSize;

		verticesSize = 0;
		verticesSize_r = 0;
		ignoreFlag = 0;
		outsideFlag = 0;
		statusFlag = 0;
		initialCellIndex = 0;
		perimeterLength = 0.0;
		targetArea = 0.0;
		divArea = -1.0;
		cofm_x = 0.0; cofm_y = 0.0; cofm_z = 0.0;
		cofv_x = 0.0; cofv_y = 0.0; cofv_z = 0.0;
		cellCycle = 0.0;
		divCycle = 0.0;
		cellcycle_update_coef = 0;
		initArray(cellArea, 0.0);
		initArray(cellStressMatrix, 0.0);
		initArray(cellGrad, 0.0);
		initArray(cellStressMatrix_Bach, 0.0);
		
		cell_angle = 0;
		line_coef = 0;
		line_out_coef = 0;
		perimeter_coef = 0;
		area_coef = 0;
		area_target = 0;
	
		cell_shape_aniso = 0.0;
	
		DCMvertices = NULL;

		rec_Size = 0;
		max_rec_Size = 0;
		rec = NULL;
		rec_buf = NULL;
		rec_rect = NULL;
		recordFlag = 0;

		lineage_Size = 0;
		lineage = NULL;
		topologyFlag = 0;

		out_cell_eli_flag = 0;
		fix_flag = 0;

		initArray(cell_stressmag_pre, 0.0);
		initArray(cell_stressani_pre, 0.0);
		initArray(cell_stressori_pre, 0.0);

	}

	DCMcell2D(const DCMcell2D& otherCell) : DCMcell(otherCell) {
		copyDCMcell2D(otherCell);
	}

	DCMcell2D& operator=(const DCMcell2D& otherCell) {
		if (this != &otherCell) {
			static_cast<DCMcell&>(*this) = otherCell;
			copyDCMcell2D(otherCell);
		}
		return *this;
	}

	void copyDCMcell2D(const DCMcell2D& otherCell) {

		verticesSize = otherCell.verticesSize;
		verticesSize_r = otherCell.verticesSize_r;
		maxVerticesSize = otherCell.maxVerticesSize;
		ignoreFlag = otherCell.ignoreFlag;
		outsideFlag = otherCell.outsideFlag;
		statusFlag = otherCell.statusFlag;
		initialCellIndex = otherCell.initialCellIndex;
		perimeterLength = otherCell.perimeterLength;

		copyArray(otherCell, cellArea);
		copyArray(otherCell, cellStressMatrix);
		copyArray(otherCell, cellGrad);
		copyArray(otherCell, cellStressMatrix_Bach);
		
		targetArea = otherCell.targetArea;
		divArea = otherCell.divArea;
		cofm_x = otherCell.cofm_x; cofm_y = otherCell.cofm_y; cofm_z = otherCell.cofm_z;
		cofv_x = otherCell.cofv_x; cofv_y = otherCell.cofv_y; cofv_z = otherCell.cofv_z;
		cellCycle = otherCell.cellCycle;
		divCycle = otherCell.divCycle;
		cellcycle_update_coef = otherCell.cellcycle_update_coef;
		
		cell_angle = otherCell.cell_angle;
		cell_shape_aniso = otherCell.cell_shape_aniso;
		
		chPCopy(otherCell, DCMvertices, Tvertex2D, maxVerticesSize);	

		rec_Size = otherCell.rec_Size;
		max_rec_Size = otherCell.max_rec_Size;
		chCopy(otherCell, rec, CellRecord, max_rec_Size);	
		chCopy(otherCell, rec_buf, CellRecord, max_rec_Size);
		chCopy(otherCell, rec_rect, CellRecord, max_rec_Size);
		recordFlag = otherCell.recordFlag;

		lineage_Size = otherCell.lineage_Size;
		chCopy(otherCell, lineage, CellRecord, lineage_Size);
		topologyFlag = otherCell.topologyFlag;

		line_coef = otherCell.line_coef;
		line_out_coef = otherCell.line_out_coef;
		perimeter_coef = otherCell.perimeter_coef;
		area_coef = otherCell.area_coef;
		area_target = otherCell.area_target;

		out_cell_eli_flag = otherCell.out_cell_eli_flag;
		fix_flag = otherCell.fix_flag;

		copyArray(otherCell, cell_stressmag_pre);
		copyArray(otherCell, cell_stressani_pre);
		copyArray(otherCell, cell_stressori_pre);

	}

	virtual ~DCMcell2D() {

		chPDelete(DCMvertices, maxVerticesSize);

		chDelete(rec);
		chDelete(rec_buf);
		chDelete(rec_rect);

		chDelete(lineage);

	}

	//////////////////////////////

	void resetEdgesNumber() {

		for (int i=0; i<verticesSize; i++) {
			for (int j=0; j<2; j++) {
				DCMvertices[i]->edgesNumber[j] = -1;
			}
		}
	}


	void setAllMoveFlagZero() {

		for (int i=0; i<verticesSize; i++) {
			DCMvertices[i]->moveFlag = 0;
		}
	}


	virtual void createDCMvertices(int vSize) {

		verticesSize = vSize;
	
		chPDelete(DCMvertices, maxVerticesSize);

		DCMvertices = new Tvertex2D*[maxVerticesSize];

		for (int i=0; i<maxVerticesSize; i++) {
			DCMvertices[i] = new Tvertex2D;
		}
			
	}

	virtual int setNewVertex(int pi, double px, double py, double pz, int pp, int pn, int oc1 = -1, int ov1 = -1, int oc2 = -1, int ov2 = -1)	{

		if (pi+1 > maxVerticesSize) {
			return 1;
		}

		if (DCMvertices[pi]->ignoreFlag > -1) {
			return 2;
		}

		verticesSize_r++;

		DCMvertices[pi]->px = px;
		DCMvertices[pi]->py = py;
		DCMvertices[pi]->pz = pz;
		DCMvertices[pi]->connection[0] = pp;
		DCMvertices[pi]->connection[1] = pn;

		if (oc1 > -1 || oc2 > -1) {
			setVertexOther(pi, oc1, ov1, oc2, ov2);
		}
	
		DCMvertices[pi]->ignoreFlag = 0;
	
		return 0;
	}

	void setVertexOther(int pi, int oc1, int ov1, int oc2, int ov2) {

		DCMvertices[pi]->othercell[0] = oc1;
		DCMvertices[pi]->othercell[1] = ov1;
		DCMvertices[pi]->othercell[2] = oc2;
		DCMvertices[pi]->othercell[3] = ov2;
	
	}

	int ignoreVertex(int pi) {

		if (DCMvertices[pi]->ignoreFlag != 0) {
			return 1;
		}

		verticesSize_r--;

		DCMvertices[pi]->ignoreFlag = 1;

		return 0;

	}


	void setPerimeterLength() {

		int n;
		double x, y, z;

		perimeterLength = 0;

		for (int i=0; i<verticesSize; i++) {
			if (DCMvertices[i]->ignoreFlag == 0) {
				n = DCMvertices[i]->connection[1];
				x = DCMvertices[n]->px - DCMvertices[i]->px;
				y = DCMvertices[n]->py - DCMvertices[i]->py;
				z = DCMvertices[n]->pz - DCMvertices[i]->pz;
				perimeterLength += sqrt3(x, y, z);
			}
		}

	}

	void setCellShape() {

		double mat[4] = {0, 0, 0, 0};
		double x, y;
	
		for (int i=0; i<verticesSize; i++) {

			if (DCMvertices[i]->ignoreFlag == 0) {
			
				x = DCMvertices[i]->px - cofm_x;
				y = DCMvertices[i]->py - cofm_y;
			
				mat[0] += x*x;
				mat[1] += x*y;
				mat[2] = mat[1];
				mat[3] += y*y;
			}
		}
	
		for (int i=0; i<4; i++) {		
			mat[i] /= static_cast<double>(verticesSize_r);
		}
	
		double sq = sqrt( (mat[0]+mat[3])*(mat[0]+mat[3]) - 4*(mat[0]*mat[3]-mat[1]*mat[2]) );
		double lm_maj = 0.5*( (mat[0]+mat[3]) + sq);
		double lm_min = 0.5*( (mat[0]+mat[3]) - sq);

		cell_angle = atan2( (lm_min-mat[0]), mat[1]);
		cell_shape_aniso = (1.0 - lm_min/lm_maj);

	}

	void setOutside() {

		outsideFlag = 0;

		int out = 0;

		for (int i=0; i<verticesSize; i++) {
			if (DCMvertices[i]->ignoreFlag == 0) {
				for (int j=0; j<2; j++) {
					if (DCMvertices[i]->othercell[2*j] < 0) {
						out++;
					}
				}
			}
		}

		if (out > 0) {
			outsideFlag = 1;
		}
	}

	void setCenterofVertices() {
	
		cofv_x = 0; cofv_y = 0; cofv_z = 0;

		for (int i=0; i<verticesSize; i++) {
			if (DCMvertices[i]->ignoreFlag == 0) {
				cofv_x += DCMvertices[i]->px;
				cofv_y += DCMvertices[i]->py;
				cofv_z += DCMvertices[i]->pz;
			}
		}

		cofv_x /= static_cast<double>(verticesSize_r);
		cofv_y /= static_cast<double>(verticesSize_r);
		cofv_z /= static_cast<double>(verticesSize_r);

	}

	virtual void setCellArea(int flag = 0) {

		setCenterofVertices();

		int pn;

		double pix, piy, pnx, pny;

		double ca;

		cellArea[1] = cellArea[0];

		cellArea[0] = 0.0;
	
		for (int i=0; i<verticesSize; i++) {
			if (DCMvertices[i]->ignoreFlag == 0) {

				pn = DCMvertices[i]->connection[1];

				pix = DCMvertices[i]->px;
				piy = DCMvertices[i]->py;

				pnx = DCMvertices[pn]->px;
				pny = DCMvertices[pn]->py;

				ca = 0.5*((cofv_x*piy+pix*pny+pnx*cofv_y)-(pix*cofv_y+pnx*piy+cofv_x*pny));
	
				cellArea[0] += ca;
			}
		}

		setCenterofMass();
	
		if (flag == 1) {
			cellArea[1] = cellArea[0];
		}

	}

	void setCenterofMass() {
	
		cofm_x = 0; cofm_y = 0; cofm_z = 0;
	
		int pn;
		double pix, piy, pnx, pny;
	
		for (int i=0; i<verticesSize; i++) {
		
			if (DCMvertices[i]->ignoreFlag == 0) {
			
				pn = DCMvertices[i]->connection[1];

				pix = DCMvertices[i]->px;
				piy = DCMvertices[i]->py;
			
				pnx = DCMvertices[pn]->px;
				pny = DCMvertices[pn]->py;
			
				cofm_x += (pix+pnx)*(pix*pny-pnx*piy);
				cofm_y += (piy+pny)*(pix*pny-pnx*piy);
			}
		}

		cofm_x /= (6.0*cellArea[0]);
		cofm_y /= (6.0*cellArea[0]);

	}


	void setTargetArea(double tArea, int flag=0) {

		targetArea = tArea;		

	}

	void setCellCycle(int flag=0) {

		if (flag == 0) {			
			cellCycle = SF.uniformDistribution(0.0, PM.CellCycleTime);
		}

		else if (flag == 1) {
			cellCycle = SF.uniformDistribution(-PM.CellCycleTime*PM.CellCycleSync*0.5, PM.CellCycleTime*PM.CellCycleSync*0.5);
		}

	}

	virtual void setCellCycleCoefficient(int flag=0) {

		if (flag == 0) {
		
			cellcycle_update_coef = 1.0;
		
			if (out_cell_eli_flag != 0) {

				if (divArea < 0) {
					divArea = cellArea[0];
					divCycle = 0.0;
				}

				divCycle += PM.Delta_t;

				if (targetArea > -1.0) {
					setTargetArea(PM.Init_TargetArea-divCycle*10.0);
				}
				cellcycle_update_coef = 0;
			}
		}

	}

	void updateCellCycle() {
		
		if (verticesSize_r < 4) cellcycle_update_coef = 0;
		
		if (fix_flag > 0) cellcycle_update_coef = 0;
		
		cellCycle += PM.Delta_t * PM.CellCycleSpeed * cellcycle_update_coef;

	}

	
	void updateRecord(int index, int timeStep, int usFlag = 0) {

		int flag = usFlag;

		if (usFlag != 0) {

			if (verticesSize_r < 6) {
			
				if (max_rec_Size == 0) {
				
					rec_Size = 0;
					max_rec_Size = 100;
				
					rec = new CellRecord[max_rec_Size];
				
				}			
			
				else {
				
					if (rec[rec_Size-1].verticesSize > verticesSize_r) {}				
					else if (rec[rec_Size-1].verticesSize < verticesSize_r) {}

				}
	
				flag = 0;	
			
			}

			if (verticesSize_r > 5) {

				if (max_rec_Size != 0) {

					chDelete(rec);
					chDelete(rec_buf);
					chDelete(rec_rect);
					max_rec_Size = 0;
				
				}
			}

			if (usFlag == US::INITIALIZE) {

				lineage_Size = 10;
				lineage = new CellRecord[lineage_Size];

				lineage[0].index = index;
				lineage[0].timeStep = timeStep;
				lineage[0].verticesSize = verticesSize_r;
				lineage[0].topologicalChange = usFlag;
				lineage[0].inout = outsideFlag;
			}

			else {

				if (topologyFlag == 0) {
					topologyFlag = usFlag;
				}
				
				if (lineage[lineage_Size-1].index > -1) {
					
					CellRecord* tempLineage = new CellRecord[lineage_Size];
					memcpy(tempLineage, lineage, sizeof(CellRecord)*lineage_Size);
					chDelete(lineage);
					lineage = new CellRecord[lineage_Size+10];
					memcpy(lineage, tempLineage, sizeof(CellRecord)*lineage_Size);
					lineage_Size += 10;
					delete[] tempLineage;
				}
				
				for (int i=0; i<lineage_Size; i++) {

					if (lineage[i].index < 0) {

						lineage[i].index = index;
						lineage[i].timeStep = timeStep;
						lineage[i].verticesSize = verticesSize_r;
						lineage[i].topologicalChange = usFlag;
						lineage[i].inout = outsideFlag;
						
						if (usFlag > US::DIVISION_OWN) {

							lineage[i].topologicalChange = US::DIVISION_OWN;
							lineage[i].daughterCell = usFlag - US::DIVISION_OWN -1;
						}

						break;
					}
				}
			}
		}


		if (max_rec_Size != 0) {
		
			if (flag == 0) {
			
				if (rec_Size == max_rec_Size) {

					for (int i=0; i<rec_Size; i++) {

						if (rec[i].topologicalChange != 0) {
							chDelete(rec_rect);
							rec_rect = new CellRecord[max_rec_Size];
							memcpy(rec_rect, rec, sizeof(CellRecord)*max_rec_Size);						
							break;
						}
					}

					chDelete(rec_buf);

					rec_buf = new CellRecord[max_rec_Size];
				
					memcpy(rec_buf, rec, sizeof(CellRecord)*max_rec_Size);
	
					chDelete(rec);			
				
					rec = new CellRecord[max_rec_Size];

					rec_Size = 0;			
				
				}

				rec[rec_Size].timeStep = timeStep;
				rec[rec_Size].index = index;
				rec[rec_Size].topologicalChange = usFlag;
				rec[rec_Size].verticesSize = verticesSize_r;

				rec[rec_Size].area = cellArea[0];

				setCellShape();
			
				rec[rec_Size].cell_shape_aniso = cell_shape_aniso;

				rec_Size++;

			}
		}


	}

	void updateTopologyFlag() {

		if (topologyFlag != 0) {

			if (abs(cellArea[0]-cellArea[1])/PM.Delta_t < 1.0e-4 && cellArea[0] != cellArea[1]) {

				setCellShape();
								
				for (int i=0; i<lineage_Size; i++) {
					
					if (i < lineage_Size-1) {

						if (lineage[i+1].index < 0) {

							lineage[i].area = cellArea[0];
							lineage[i].cell_shape_aniso = cell_shape_aniso;
							break;
						}
					}
					else {
						lineage[i].area = cellArea[0];
						lineage[i].cell_shape_aniso = cell_shape_aniso;
					}
				}

				topologyFlag = 0;

			}
		}
	}

	int updateReindexVertex() {

		int vs = verticesSize_r;
		
		if (vs < maxVerticesSize-4) {
			
			if (verticesSize > maxVerticesSize-4) {
				
				Tvertex2D** tempVertices = new Tvertex2D*[vs];
				int* tempIndex = new int[vs];
				
				for (int i=0; i<vs; i++) {					
					tempVertices[i] = new Tvertex2D;
				}
				
				int ps = 0;
				int pn = 0;
				int pi = 0;
				
				for (int i=0; i<verticesSize; i++) {
					if (DCMvertices[i]->ignoreFlag == 0) {
						pn = ps = i;
						break;
					}
				}
				
				do {
					tempIndex[pi] = pn;
					pi++;
					pn = DCMvertices[pn]->connection[1];
				} while(ps != pn);
				
				for (int i=0; i<vs; i++) {
					int p = i-1;
					int n = i+1;
					if (i == 0) p = vs-1;
					if (i == vs-1) n = 0;
					
					*tempVertices[i] = *DCMvertices[tempIndex[i]];
					tempVertices[i]->connection[0] = p;
					tempVertices[i]->connection[1] = n;

					updateReindexVertex_addon(tempIndex[i], i);
				}
				
				for (int i=0; i<maxVerticesSize; i++) {
					
					Tvertex2D initVertex2D;
					
					if (i < vs) {
						*DCMvertices[i] = *tempVertices[i];
					}
					else {
						*DCMvertices[i] = initVertex2D;
					}
				}

				verticesSize = vs;
				
				delete[] tempIndex;
				
				chPDelete(tempVertices, vs);

			}
			else return 0;
		}
		else return 2;

		return 1;

	}

	virtual void updateReindexVertex_addon(int index, int new_index) {

	}

	void setCellStressMatrix(int flag=0) {

		if (flag > 0) {
			cell_stressmag_pre[1] = cell_stressmag_pre[0];
			cell_stressani_pre[1] = cell_stressani_pre[0];
			cell_stressori_pre[1] = cell_stressori_pre[0];
		}
		
		initArray(cellStressMatrix, 0);
		
		int ps = 0;
		int pn, pp;
		
		for (int i=0; i<verticesSize; i++) {
			if (DCMvertices[i]->ignoreFlag == 0) {
				pp = ps = i;
				break;
			}
		}
	
		do {

			pn = DCMvertices[pp]->connection[1];

			setLineCoefficient(pp);
			double line = line_coef;
			if (DCMvertices[pp]->ocIndex < 0) line = line_out_coef;
			
			setPerimeterCoefficient();
			double peri = perimeter_coef*perimeterLength;

			double pp_x = DCMvertices[pn]->px-DCMvertices[pp]->px;
			double pp_y = DCMvertices[pn]->py-DCMvertices[pp]->py;
			
			double pp_d = 1.0/dist(DCMvertices[pp]->px, DCMvertices[pn]->px, DCMvertices[pp]->py, DCMvertices[pn]->py, 0, 0);
			
			cellStressMatrix[0] += (line+peri)*pp_x*pp_x*pp_d;
			cellStressMatrix[1] += (line+peri)*pp_x*pp_y*pp_d;
			cellStressMatrix[2] += (line+peri)*pp_y*pp_y*pp_d;

			pp = pn;
			
		} while(pn != ps);

		setAreaCoefficient();				
		cellStressMatrix[0] = cellStressMatrix[0]/cellArea[0] + area_coef*(cellArea[0]-area_target);
		cellStressMatrix[1] = cellStressMatrix[1]/cellArea[0];
		cellStressMatrix[2] = cellStressMatrix[2]/cellArea[0] + area_coef*(cellArea[0]-area_target);

		double sq = sqrt( (cellStressMatrix[0]+cellStressMatrix[2])*(cellStressMatrix[0]+cellStressMatrix[2])
			- 4.0*(cellStressMatrix[0]*cellStressMatrix[2]-cellStressMatrix[1]*cellStressMatrix[1]) );

		cellStressMatrix[4] = 0.5*( (cellStressMatrix[0]+cellStressMatrix[2]) + sq);
		cellStressMatrix[5] = 0.5*( (cellStressMatrix[0]+cellStressMatrix[2]) - sq);

		if (cellStressMatrix[4] * cellStressMatrix[5] > 0) cellStressMatrix[3] = sqrt(cellStressMatrix[4] * cellStressMatrix[5]);
		else cellStressMatrix[3] = -sqrt(-cellStressMatrix[4] * cellStressMatrix[5]);

		cellStressMatrix[6] = cellStressMatrix[1];
		cellStressMatrix[7] = cellStressMatrix[4]-cellStressMatrix[0];
		cellStressMatrix[8] = cellStressMatrix[1];
		cellStressMatrix[9] = cellStressMatrix[5]-cellStressMatrix[0];
		
		double lmax_inorm = 1.0/sqrt(cellStressMatrix[6]*cellStressMatrix[6] + cellStressMatrix[7]*cellStressMatrix[7]);
		double lmin_inorm = 1.0/sqrt(cellStressMatrix[8]*cellStressMatrix[8] + cellStressMatrix[9]*cellStressMatrix[9]);
		
		cellStressMatrix[6] *= lmax_inorm;
		cellStressMatrix[7] *= lmax_inorm;
		cellStressMatrix[8] *= lmin_inorm;
		cellStressMatrix[9] *= lmin_inorm;

		if (flag > 0) {
			cell_stressmag_pre[0] = cellStressMatrix[4]+cellStressMatrix[5];
			cell_stressmag_pre[2] = (cell_stressmag_pre[0]-cell_stressmag_pre[1])/static_cast<double>(flag);
			cell_stressani_pre[0] = cellStressMatrix[4]-cellStressMatrix[5];
			cell_stressani_pre[2] = (cell_stressani_pre[0]-cell_stressani_pre[1])*20.0/static_cast<double>(flag);
			cell_stressori_pre[0] = atan2(cellStressMatrix[7], cellStressMatrix[6]);
			cell_stressori_pre[2] = abs(cell_stressori_pre[0]-cell_stressori_pre[1])/static_cast<double>(flag);
		}

	}


	void setcellStressMatrix_Bach() {

		initArray(cellStressMatrix_Bach, 0);
		
		int ps = 0;
		int pn, pp;
		
		for (int i=0; i<verticesSize; i++) {
			if (DCMvertices[i]->ignoreFlag == 0) {				
				pp = ps = i;
				break;
			}
		}
		
		do {
			
			pn = DCMvertices[pp]->connection[1];
			
			setLineCoefficient(pp);
			double line = line_coef;
			if (DCMvertices[pp]->ocIndex < 0) line = line_out_coef;
			
			setPerimeterCoefficient();
			double peri = perimeter_coef*perimeterLength;

			double pp_x = DCMvertices[pn]->px-DCMvertices[pp]->px;
			double pp_y = DCMvertices[pn]->py-DCMvertices[pp]->py;
			
			double pp_d = 1.0/dist(DCMvertices[pp]->px, DCMvertices[pn]->px, DCMvertices[pp]->py, DCMvertices[pn]->py, 0, 0);
			
			cellStressMatrix_Bach[0] += (line+peri)*pp_x*pp_x*pp_d;
			cellStressMatrix_Bach[1] += (line+peri)*pp_x*pp_y*pp_d;
			cellStressMatrix_Bach[2] += (line+peri)*pp_y*pp_y*pp_d;
			
			pp = pn;
		
		} while(pn != ps);

		setAreaCoefficient();				
		cellStressMatrix_Bach[0] = cellStressMatrix_Bach[0]/cellArea[0] + area_coef*(cellArea[0]-area_target);
		cellStressMatrix_Bach[1] = cellStressMatrix_Bach[1]/cellArea[0];
		cellStressMatrix_Bach[2] = cellStressMatrix_Bach[2]/cellArea[0] + area_coef*(cellArea[0]-area_target);
		
		double sq = sqrt( (cellStressMatrix_Bach[0]+cellStressMatrix_Bach[2])*(cellStressMatrix_Bach[0]+cellStressMatrix_Bach[2])
			- 4.0*(cellStressMatrix_Bach[0]*cellStressMatrix_Bach[2]-cellStressMatrix_Bach[1]*cellStressMatrix_Bach[1]) );
		
		cellStressMatrix_Bach[4] = 0.5*( (cellStressMatrix_Bach[0]+cellStressMatrix_Bach[2]) + sq);
		cellStressMatrix_Bach[5] = 0.5*( (cellStressMatrix_Bach[0]+cellStressMatrix_Bach[2]) - sq);
		
		if (cellStressMatrix_Bach[4] * cellStressMatrix_Bach[5] > 0) cellStressMatrix_Bach[3] = sqrt(cellStressMatrix_Bach[4] * cellStressMatrix_Bach[5]);
		else cellStressMatrix_Bach[3] = -sqrt(-cellStressMatrix_Bach[4] * cellStressMatrix_Bach[5]);
		
		cellStressMatrix_Bach[6] = cellStressMatrix_Bach[1];
		cellStressMatrix_Bach[7] = cellStressMatrix_Bach[4]-cellStressMatrix_Bach[0];
		cellStressMatrix_Bach[8] = cellStressMatrix_Bach[1];
		cellStressMatrix_Bach[9] = cellStressMatrix_Bach[5]-cellStressMatrix_Bach[0];
		
		double lmax_inorm = 1.0/sqrt(cellStressMatrix_Bach[6]*cellStressMatrix_Bach[6] + cellStressMatrix_Bach[7]*cellStressMatrix_Bach[7]);
		double lmin_inorm = 1.0/sqrt(cellStressMatrix_Bach[8]*cellStressMatrix_Bach[8] + cellStressMatrix_Bach[9]*cellStressMatrix_Bach[9]);
		
		cellStressMatrix_Bach[6] *= lmax_inorm;
		cellStressMatrix_Bach[7] *= lmax_inorm;
		cellStressMatrix_Bach[8] *= lmin_inorm;
		cellStressMatrix_Bach[9] *= lmin_inorm;

	}

	void setCellGrad(double x, double y, double ax, double ay) {

		initArray(cellGrad, -1);

		cellGrad[0] = x;
		cellGrad[1] = y;

		cellGrad[2] = sqrt(x*x+y*y);

		cellGrad[3] = ax;
		cellGrad[4] = ay;

		cellGrad[5] = sqrt(ax*ax+ay*ay);

	}

	void setOtherCellIndex() {
		
		int ps = 0;
		int pn, pp;
		int oc1, oc2;

		for (int i=0; i<verticesSize; i++) {
			if (DCMvertices[i]->ignoreFlag == 0) {
				pp = ps = i;
				break;
			}
		}
	
		do {

			pn = DCMvertices[pp]->connection[1];		

			for (int i=0; i<2; i++) {
				for (int j=0; j<2; j++) {
					oc1 = DCMvertices[pp]->othercell[2*i];
					oc2 = DCMvertices[pn]->othercell[2*j];

					if (oc1 == oc2) {
						DCMvertices[pp]->ocIndex = oc1;
						DCMvertices[pp]->ovIndex = DCMvertices[pn]->othercell[2*j+1];
					
						if (oc1 > -1) {

						}
					}				
				}
			}

			pp = pn;
			
		} while(pn != ps);

	}

	virtual void setLineCoefficient(int pi, int flag=0) {

		line_coef = 1.0;
		line_coef *= PM.Sigma*0.5;

		line_out_coef = 1.0;
		line_out_coef *= PM.Sigma*PM.Sigma_Zero;

	}

	int vDynamics_line(int pi, DCMpoint& xy) {

		int pp = DCMvertices[pi]->connection[0];
		int pn = DCMvertices[pi]->connection[1];

		double pix = DCMvertices[pi]->px;
		double piy = DCMvertices[pi]->py;
		double piz = DCMvertices[pi]->pz;
	
		double ppx = DCMvertices[pp]->px;
		double ppy = DCMvertices[pp]->py;
		double ppz = DCMvertices[pp]->pz;

		double pnx = DCMvertices[pn]->px;
		double pny = DCMvertices[pn]->py;
		double pnz = DCMvertices[pn]->pz;
	
		double dist_pp = dist(ppx, pix, ppy, piy, ppz, piz); 
		if (dist_pp > VERY_SMALL) dist_pp = 1.0/dist_pp;
		else return 201;

		if (DCMvertices[pp]->ocIndex > -1) dist_pp *= line_coef;
		else dist_pp *= line_out_coef;
		
		double dist_pn = dist(pnx, pix, pny, piy, pnz, piz);
		if (dist_pn > VERY_SMALL) dist_pn = 1.0/dist_pn;
		else return 201;

		if (DCMvertices[pi]->ocIndex > -1) dist_pn *= line_coef;
		else dist_pn *= line_out_coef;
		
		double dx = (ppx-pix)*dist_pp + (pnx-pix)*dist_pn;
		double dy = (ppy-piy)*dist_pp + (pny-piy)*dist_pn;
		double dz = (ppz-piz)*dist_pp + (pnz-piz)*dist_pn;

		xy.x += dx;
		xy.y += dy;
		xy.z += dz;

		return 0;

	}

	virtual void setPerimeterCoefficient(int flag=0) {

		perimeter_coef = 1.0;
		perimeter_coef *= PM.Line_Coef;

	}

	int vDynamics_perimeter(int pi, DCMpoint& xy) {
		
		int pp = DCMvertices[pi]->connection[0];
		int pn = DCMvertices[pi]->connection[1];

		double pix = DCMvertices[pi]->px;
		double piy = DCMvertices[pi]->py;
		double piz = DCMvertices[pi]->pz;
	
		double ppx = DCMvertices[pp]->px;
		double ppy = DCMvertices[pp]->py;
		double ppz = DCMvertices[pp]->pz;

		double pnx = DCMvertices[pn]->px;
		double pny = DCMvertices[pn]->py;
		double pnz = DCMvertices[pn]->pz;
	
		double dist_pp = dist(ppx, pix, ppy, piy, ppz, piz); 
		if (dist_pp > VERY_SMALL) dist_pp = 1.0/dist_pp;
		else return 201;

		double dist_pn = dist(pnx, pix, pny, piy, pnz, piz);
		if (dist_pn > VERY_SMALL) dist_pn = 1.0/dist_pn;
		else return 201;
		
		double dx = perimeter_coef*perimeterLength*( (ppx-pix)*dist_pp + (pnx-pix)*dist_pn );
		double dy = perimeter_coef*perimeterLength*( (ppy-piy)*dist_pp + (pny-piy)*dist_pn );
		double dz = perimeter_coef*perimeterLength*( (ppz-piz)*dist_pp + (pnz-piz)*dist_pn );

		xy.x += dx;
		xy.y += dy;
		xy.z += dz;

		return 0;

	}

	virtual void setAreaCoefficient(int flag=0) {

		area_coef = 1.0;
		area_coef *= PM.Kappa;

		area_target = targetArea;

	}

	virtual int vDynamics_area(int pi, DCMpoint& xy) {

		int pp = DCMvertices[pi]->connection[0];
		int pn = DCMvertices[pi]->connection[1];
		
		double pix = DCMvertices[pi]->px;
		double piy = DCMvertices[pi]->py;

		double ppx = DCMvertices[pp]->px;
		double ppy = DCMvertices[pp]->py;

		double pnx = DCMvertices[pn]->px;
		double pny = DCMvertices[pn]->py;

		double dx = area_coef*0.5*(ppy-pny)*(cellArea[0]-area_target);
		double dy = area_coef*0.5*(pnx-ppx)*(cellArea[0]-area_target);

		xy.x += dx;
		xy.y += dy;
	
		return 0;

	}

	int vDynamics_angle(int pi, DCMpoint& xy) {

		int pp = DCMvertices[pi]->connection[0];
		int pn = DCMvertices[pi]->connection[1];

		double angle_coef = PM.Tug_Force;

		double x_pi = DCMvertices[pp]->px - DCMvertices[pi]->px;
		double y_pi = DCMvertices[pp]->py - DCMvertices[pi]->py;
		double id_pi = 1.0/sqrt(x_pi*x_pi+y_pi*y_pi);
	
		double x_ni = DCMvertices[pn]->px - DCMvertices[pi]->px;
		double y_ni = DCMvertices[pn]->py - DCMvertices[pi]->py;
		double id_ni = 1.0/sqrt(x_ni*x_ni+y_ni*y_ni);

		double s_x = x_pi*id_pi*id_pi + x_ni*id_ni*id_ni;
		double s_y = y_pi*id_pi*id_pi + y_ni*id_ni*id_ni;

		if ( ((x_pi*x_ni+y_pi*y_ni)*id_pi*id_ni) > sqrt(3.0)*0.5 && ((-x_pi*y_ni+y_pi*x_ni)*id_pi*id_ni) > -0.5 ) angle_coef = -angle_coef;

		xy.x += -angle_coef * id_pi*id_ni* (0.25*(-x_pi-x_ni + s_x*(x_pi*x_ni+y_pi*y_ni)) - sqrt(3.0)*0.25*( y_pi-y_ni + s_x*(x_pi*y_ni-y_pi*x_ni)));
		xy.y += -angle_coef * id_pi*id_ni* (0.25*(-y_pi-y_ni + s_y*(x_pi*x_ni+y_pi*y_ni)) - sqrt(3.0)*0.25*(-x_pi+x_ni + s_y*(x_pi*y_ni-y_pi*x_ni)));

		return 0;

	}


};


///////////////////////////////////////////////////////////
// DCMcontainer class


template <typename Tcell2D, typename Tvertex2D, typename Tedge2D, typename Tpoint2D>
class DCMcontainer2D : public DCMcontainer {

public:
		
	Tcell2D** DCMcells;
	
	Tedge2D** edges;
	
	Tpoint2D** vertices;

	int container_index;

	int edgesArraySize;

	int verticesArraySize;
	
	int initialCellSize;
	double initial_average_volume;
		
	int recon_number;
	int recon_number_r;

	ReconnectionRecord* rec_recon_info;
	int rec_recon_info_size;

	int extrusion_number;

	int GrowthDivisionFlag;
	int vertex_dynamics_flag;

	double mphase_coef;
	double vonMisesKappa;

	stressRecord* rec_StressDivision;
	int rec_StressDivision_size;
	int rec_StressDivision_index;

	stressRecord* rec_StressApoptosis;
	int rec_StressApoptosis_size;
	int rec_StressApoptosis_index;

	int cell_size_counter;
	int cell_size_update_flag;

	int recordable_flag;

	double adder_parameter;
	
	//////////////////////////////

	DCMcontainer2D() : DCMcontainer() {

		container_index = 0;
		cellsSize = 0;
		cellsSize_r = 0;
		edgesSize = 0;
		verticesSize = 0;
		planesSize = 0;
		containerVerticesSize = 0;
		edgesArraySize = 0;
		verticesArraySize = 0;
		averageVolume = 0;
		initial_average_volume = 0;
		reconnectionFlag = 0;
		objectFlag = 0;
		conTimeStep = 0;
		cCenterx = 0.0; cCentery = 0.0; cCenterz = 0.0;
		basex = 0.0; basey = 0.0; basez = 0.0;
		tresholdCos = 0.0;
		outVerticesIndex = 0;
		outVertexMax = 0.0;
		dataVersion = DATA_VERSION;
		timeStepMaxVectorForce = 0.0;
		laps = 0;
		maxdxdt = 0;
		maxmaxdxdt = 0;
		initialCellSize = 0;
		recon_number = 0;
		recon_number_r = 0;
		extrusion_number = 0;
		cell_size_counter = 0;
		cell_size_update_flag = 0;
		recordable_flag = 0;
		adder_parameter = -1;

		rec_recon_info = NULL;
		rec_recon_info_size = 0;

		DCMcells = NULL;
		edges = NULL;
		vertices = NULL;
	
		GrowthDivisionFlag = 0;
		vertex_dynamics_flag = 0;
		mphase_coef = 0.0;

		vonMisesKappa = 0.0;

		initDataFlag = 0;
		
		rec_StressDivision = NULL;
		rec_StressDivision_size = 0;
		rec_StressDivision_index = 0;
		
		rec_StressApoptosis = NULL;
		rec_StressApoptosis_size = 0;
		rec_StressApoptosis_index = 0;
		
	}

	virtual ~DCMcontainer2D() {

		chDelete(rec_recon_info);

		chPDelete(edges, edgesArraySize);
		chPDelete(vertices, verticesArraySize);

		chPDelete(DCMcells, maxCellsSize);
		
		chDelete(rec_StressDivision);
		chDelete(rec_StressApoptosis);

	}

	virtual void destroy() {}

	//////////////////////////////

	virtual void createDCMcells(int cSize) {

		cellsSize = cSize;

		chPDelete(DCMcells, maxCellsSize);

		DCMcells = new Tcell2D*[maxCellsSize];

		for (int i=0; i<maxCellsSize; i++) {
			DCMcells[i] = new Tcell2D;
		}

	}

	virtual void createDCMvertices(int cIndex, int vSize) {

		DCMcells[cIndex]->createDCMvertices(vSize);
		
	}

	virtual void setNewVertex(int cIndex, int pi, double px, double py, double pz, int pp, int pn, int oc1 = -1, int ov1 = -1, int oc2 = -1, int ov2 = -1) {

		DCMcells[cIndex]->setNewVertex(pi, px, py, pz, pp, pn, oc1, ov1, oc2, ov2);

	}

	virtual void initialize_loadfile() {

		vonMisesKappa = 0;  // 0.5, 2, 8, -1
		
		std::ifstream ifs(L"DCModel_param_basic.txt");

		if (ifs) {
			std::string str;

			std::vector<std::string> param_list;
			param_list.push_back("Lambda");		// 0
			param_list.push_back("Gamma");			// 1
			param_list.push_back("Outside boundary lambda");	// 2
			param_list.push_back("Initial cell size");		// 3
			param_list.push_back("Delta t");		// 4
			param_list.push_back("T1 threshold");	// 5
			param_list.push_back("T2 threshold");	// 6
			param_list.push_back("Division orientation");	// 7
			param_list.push_back("Growth rate");	// 8
			param_list.push_back("Randomness of Lambda");	// 9
			param_list.push_back("Decay of Lambda");	// 10
			param_list.push_back("Initial distribution of Lambda");		// 11
			param_list.push_back("Outside boundary fix");	// 12
			param_list.push_back("Adder parameter");		// 13
		
			while(std::getline(ifs, str)) {

				std::string token1;
				std::string token2;
				std::istringstream stream(str);
			
				std::getline(stream, token1, ':');
				std::getline(stream, token2);
				
				for (unsigned int i=0; i<param_list.size(); i++) {
					if (token1.compare(param_list[i]) == 0 && !token2.empty()) {
					switch (i) {
					case 0: PM.Sigma = std::stod(token2); break;
					case 1: PM.Line_Coef = std::stod(token2); break;
					case 2: PM.Sigma_Zero = std::stod(token2); break;
					case 3: PM.initial_cell_size = std::stoi(token2); break;
					case 4: PM.Delta_t = std::stod(token2); break;
					case 5: PM.Epsilon_Distance = std::stod(token2); break;
					case 6: PM.Extrusion_Area = std::stod(token2); break;
					case 7: vonMisesKappa = std::stod(token2); break;
					case 8: PM.CellCycleSpeed = std::stod(token2); break;
					case 9: PM.Randomness_of_Lambda = std::stod(token2); break;
					case 10: PM.Decay_of_Lambda = std::stod(token2); break;
					case 11: if (std::stoi(token2) > 0) vertex_dynamics_flag ^= VD::RANDLAMBDA_INIT; break;
					case 12: if (std::stoi(token2) > 0) vertex_dynamics_flag ^= VD::OUT_FIX; break;
					case 13: if (std::stod(token2) > 0) adder_parameter = std::stod(token2); break;
						}
					}
				}
			}
		}		
	}
	
	virtual void initialize(int flag = 0) {

		container_index = flag;

		vertex_dynamics_flag ^= VD::T2TH_RELATIVE;

		vertex_dynamics_flag ^= VD::LINETENSION;
		vertex_dynamics_flag ^= VD::PERIMETER;

		if ((vertex_dynamics_flag & VD::RANDOMMOVE) == VD::RANDOMMOVE) PM.RandomMove = 0.01;

		if (PM.Randomness_of_Lambda > 0) vertex_dynamics_flag ^= VD::RANDOMLAMBDA;
		
		if (adder_parameter > 0) GrowthDivisionFlag ^= BC::MPHASE_ADD;
		else GrowthDivisionFlag ^= BC::MPHASE_2SD;	

		mphase_coef = 10.0;

		cell_size_update_flag = 1000;

		setOtherCell();
	
		updateTopology();

		initialCellSize = cellsSize;
		
		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {

				DCMcells[i]->initialCellIndex = i;

				DCMcells[i]->setCellCycle(0);
				
				DCMcells[i]->setTargetArea(PM.Init_TargetArea);

				if ((recordable_flag & RE::TOPOLOGY_CHANGE) == RE::TOPOLOGY_CHANGE) DCMcells[i]->updateRecord(i, conTimeStep/100, US::INITIALIZE);
				
			}
		}
		
		if ((vertex_dynamics_flag & VD::RANDOMLAMBDA) == VD::RANDOMLAMBDA) {
			if ((vertex_dynamics_flag & VD::RANDLAMBDA_INIT) == VD::RANDLAMBDA_INIT) {
				for (int i=0; i<cellsSize; i++) {
					if (DCMcells[i]->ignoreFlag == 0) {
						for (int j=0; j<DCMcells[i]->verticesSize; j++) {
							if (DCMcells[i]->DCMvertices[j]->ignoreFlag == 0) {
								setNewCellvertex_EdgeLineCoef(i,j);
								int edge_index = search_Cellvertex_to_Edge(i,j);
								edges[edge_index]->randomness_of_line_coef = DCMcells[i]->DCMvertices[j]->edge_randomness_of_line_coef;
							}
						}
					}
				}
			}
		}

		setAverageVolume();
		initial_average_volume = averageVolume;

	}


	void setAverageVolume() {

		double aVolume = 0.0;

		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				aVolume += DCMcells[i]->cellArea[0];
			}
		}
	
		aVolume /= static_cast<double>(cellsSize_r);
		averageVolume = aVolume;
	
	}

	virtual void setOtherCell() {
		
		double x, y;
		double dx, dy, d;

		int index = 1;
		int mVSize = 0;

		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				mVSize += DCMcells[i]->verticesSize_r;
			}		
		}
	

		Tpoint2D* point = new Tpoint2D[mVSize+1];

		point[0].x = VERY_BIG; point[0].y = VERY_BIG;
	

		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				
				for (int j=0; j<DCMcells[i]->verticesSize; j++) {
					if (DCMcells[i]->ignoreFlag == 0) {
				
						x = DCMcells[i]->DCMvertices[j]->px;
						y = DCMcells[i]->DCMvertices[j]->py;
						
						int check = 0;

						int k = 0;

						do {

							dx = point[k].x - x; dy = point[k].y - y;
							d = dx*dx + dy*dy;

							if (d < PM.Othercell_Distance) {

								int vertexN = point[k].vertexNum + 1;
								point[k].cellIndex[vertexN] = i;
								point[k].vertexIndex[vertexN] = j;
								point[k].vertexNum++;

								switch (vertexN) {

									case 1:
										DCMcells[point[k].cellIndex[0]]->DCMvertices[point[k].vertexIndex[0]]->othercell[0] = point[k].cellIndex[1];
										DCMcells[point[k].cellIndex[0]]->DCMvertices[point[k].vertexIndex[0]]->othercell[1] = point[k].vertexIndex[1];

										DCMcells[point[k].cellIndex[1]]->DCMvertices[point[k].vertexIndex[1]]->othercell[0] = point[k].cellIndex[0];
										DCMcells[point[k].cellIndex[1]]->DCMvertices[point[k].vertexIndex[1]]->othercell[1] = point[k].vertexIndex[0];

										DCMcells[point[k].cellIndex[1]]->DCMvertices[point[k].vertexIndex[1]]->px = point[k].x;
										DCMcells[point[k].cellIndex[1]]->DCMvertices[point[k].vertexIndex[1]]->py = point[k].y;							

										break;

									case 2:
										DCMcells[point[k].cellIndex[0]]->DCMvertices[point[k].vertexIndex[0]]->othercell[2] = point[k].cellIndex[2];
										DCMcells[point[k].cellIndex[0]]->DCMvertices[point[k].vertexIndex[0]]->othercell[3] = point[k].vertexIndex[2];

										DCMcells[point[k].cellIndex[1]]->DCMvertices[point[k].vertexIndex[1]]->othercell[2] = point[k].cellIndex[2];
										DCMcells[point[k].cellIndex[1]]->DCMvertices[point[k].vertexIndex[1]]->othercell[3] = point[k].vertexIndex[2];

										DCMcells[point[k].cellIndex[2]]->DCMvertices[point[k].vertexIndex[2]]->othercell[0] = point[k].cellIndex[0];
										DCMcells[point[k].cellIndex[2]]->DCMvertices[point[k].vertexIndex[2]]->othercell[1] = point[k].vertexIndex[0];

										DCMcells[point[k].cellIndex[2]]->DCMvertices[point[k].vertexIndex[2]]->othercell[2] = point[k].cellIndex[1];
										DCMcells[point[k].cellIndex[2]]->DCMvertices[point[k].vertexIndex[2]]->othercell[3] = point[k].vertexIndex[1];

										DCMcells[point[k].cellIndex[2]]->DCMvertices[point[k].vertexIndex[2]]->px = point[k].x;
										DCMcells[point[k].cellIndex[2]]->DCMvertices[point[k].vertexIndex[2]]->py = point[k].y;
							
										break;

								}

								check = 1;					
							}

						if (check == 1) break;
							
						k++;

						} while (k < index);

						if (check == 0) {

							point[index].x = x; point[index].y = y;
							point[index].cellIndex[0] = i; point[index].vertexIndex[0] = j;
								
							index++;
						}
					}	
				}
			}
		}

		delete[] point;

	}

	virtual int finishContainer(int flag=0) {

		int cs = 0;
			for (int i=0; i<cellsSize; i++) {
				if (DCMcells[i]->ignoreFlag == 0) {
					if (DCMcells[i]->outsideFlag == 0) cs++;
				}
		}

		if (cs >= 20000) return 1;

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
		}

		if (flag == 2) {
			if (cell_size_counter >= cell_size_update_flag) {
				cell_size_update_flag += 1000;
				return 1;
			}
		}

		return 0;
	}


	double distanceFromCenter(double x, double y) {

		double max_x = -(VERY_BIG);
		double min_x = VERY_BIG;
		double max_y = -(VERY_BIG);
		double min_y = VERY_BIG;
		double rr = 0;

		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				if (max_x < DCMcells[i]->cofm_x) max_x = DCMcells[i]->cofm_x;
				if (min_x > DCMcells[i]->cofm_x) min_x = DCMcells[i]->cofm_x;
				if (max_y < DCMcells[i]->cofm_y) max_y = DCMcells[i]->cofm_y;
				if (min_y > DCMcells[i]->cofm_y) min_y = DCMcells[i]->cofm_y;
			}
		}

		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				double dx = DCMcells[i]->cofm_x-((max_x-min_x)*0.5+min_x);
				double dy = DCMcells[i]->cofm_y-((max_y-min_y)*0.5+min_y);
				double r = sqrt(dx*dx+dy*dy);
				if (rr < r) rr = r;
			}
		}

		double dx = x - ((max_x-min_x)*0.5 + min_x);
		double dy = y - ((max_y-min_y)*0.5 + min_y);

		return sqrt(dx*dx+dy*dy)/rr;

	}


	void search_NeighborCells(NeighborCell* result, int result_size, int start_cell_index, int level) {

		int focal_cell = start_cell_index;
		int nni = 1;
		int nni_checked_start = 0;
		int nni_checked_end = 1;
		int nnl = 1;

		result[0].index = start_cell_index;
		result[0].nn_level = 0;
				
		do {
			int double_check = 0;

			for (int k=nni_checked_start; k<nni_checked_end; k++) {
						
				focal_cell = result[k].index;
				
				for (int j=0; j<DCMcells[focal_cell]->verticesSize; j++) {
					if (DCMcells[focal_cell]->DCMvertices[j]->ignoreFlag == 0) {
						int oc = DCMcells[focal_cell]->DCMvertices[j]->ocIndex;
							
						if (oc> -1) {

							int check = 0;

							for (int jj=0; jj<nni; jj++) {
								if (result[jj].index == oc || start_cell_index == oc) {
									check = 1;
									break;
								}
							}
								
							if (check == 0) {
								result[nni].index = oc;
								result[nni].nn_level = nnl;
								nni++;
							}
						}
					}

					if (nni > result_size-1) {
						double_check = 1;
						break;
					}
				}

				if (double_check == 1) break;					
				
			}

			if (double_check == 1) break;

			nni_checked_start = nni_checked_end;
			nni_checked_end = nni;					
			nnl++;

		} while (nnl < level+1);

	}


	virtual void drawEdges(Vpoint* vpoint, unsigned int& size) {
	
		size = edgesSize;
		
		for(int i=0; i<edgesSize; i++) {

			vpoint[i].x1 = edges[i]->v1x;
			vpoint[i].y1 = edges[i]->v1y;
			vpoint[i].z1 = edges[i]->v1z;
			
			vpoint[i].x2 = edges[i]->v2x;
			vpoint[i].y2 = edges[i]->v2y;
			vpoint[i].z2 = edges[i]->v2z;

			vpoint[i].colorflag = 8;
		}

		if ((PM.view_flag & 1<<10) == 1<<10 || (PM.view_flag & 1<<11) == 1<<11) {
			
			int ci = 0;
			
			for (int i=0; i<cellsSize; i++) {				
				if (DCMcells[i]->ignoreFlag == 0) {
					
					double aniso = 0.2;

					if ((PM.view_flag & 1<<11) == 1<<11) {
						aniso = 0.8*(DCMcells[i]->cell_stressani_pre[0]);
					}
					
					vpoint[size+ci].x1 = DCMcells[i]->cofm_x - aniso*DCMcells[i]->cellStressMatrix[6];
					vpoint[size+ci].y1 = DCMcells[i]->cofm_y - aniso*DCMcells[i]->cellStressMatrix[7];
					vpoint[size+ci].z1 = 0;
					
					vpoint[size+ci].x2 = DCMcells[i]->cofm_x + aniso*DCMcells[i]->cellStressMatrix[6];
					vpoint[size+ci].y2 = DCMcells[i]->cofm_y + aniso*DCMcells[i]->cellStressMatrix[7];
					vpoint[size+ci].z2 = 0;
					
					vpoint[size+ci].colorflag = 9;
					
					ci++;
				}
			}
			
			size += ci;
		}

	}


	virtual void drawPlanes(VTpoint* vtpoint, int& size) {

		double mag;

		size = 0;
		
		for (int i=0; i<cellsSize; i++) {

			if (DCMcells[i]->ignoreFlag == 0) {			
			
				for (int j=0; j<DCMcells[i]->verticesSize; j++) {
					if (DCMcells[i]->DCMvertices[j]->ignoreFlag == 0) {

						vtpoint[size].x1 = DCMcells[i]->cofv_x;
						vtpoint[size].y1 = DCMcells[i]->cofv_y;
						vtpoint[size].z1 = DCMcells[i]->cofv_z;
					
						vtpoint[size].x2 = DCMcells[i]->DCMvertices[j]->px;
						vtpoint[size].y2 = DCMcells[i]->DCMvertices[j]->py;
						vtpoint[size].z2 = DCMcells[i]->DCMvertices[j]->pz;

						vtpoint[size].x3 = DCMcells[i]->DCMvertices[DCMcells[i]->DCMvertices[j]->connection[1]]->px;
						vtpoint[size].y3 = DCMcells[i]->DCMvertices[DCMcells[i]->DCMvertices[j]->connection[1]]->py;
						vtpoint[size].z3 = DCMcells[i]->DCMvertices[DCMcells[i]->DCMvertices[j]->connection[1]]->pz;
					
						if (initDataFlag > 0 && initDataFlag < 10) {							
							vtpoint[size].colorflag = 6;							
							vtpoint[size].colorGrad = DCMcells[i]->initialCellIndex;
						}
						else {

							if ((PM.view_flag & ((1<<10)-1)) == 0) {
								vtpoint[size].colorflag = DCMcells[i]->verticesSize_r-3;
								if (vtpoint[size].colorflag > 7) vtpoint[size].colorflag = 7;
								vtpoint[size].colorGrad = 1;
							}								
							else if ((PM.view_flag & 16) == 16) {
								vtpoint[size].colorflag = 7;
								vtpoint[size].colorGrad = 1;
							}
							else {
								vtpoint[size].colorflag = 101;

								if ((PM.view_flag & 1) == 1) {
									vtpoint[size].colorGrad = 120-static_cast<int>((DCMcells[i]->cellArea[0]-averageVolume)/averageVolume*240);
								}
								else if ((PM.view_flag & 2) == 2) {
									mag = DCMcells[i]->cellStressMatrix[4]+DCMcells[i]->cellStressMatrix[5];
									vtpoint[size].colorGrad = static_cast<int>(120-(mag*360));
								}
								
								if (vtpoint[size].colorGrad > 240) vtpoint[size].colorGrad = 240;
								if (vtpoint[size].colorGrad < 0) vtpoint[size].colorGrad = 0;
							}
						}

						size++;
					
					}
				}
			}
		}

	}

	void setEdgeNegative() {

		for (int i=0; i<edgesSize; i++) {
			edges[i]->negative = 1;
		}

		for (int i=0; i<verticesSize; i++) {

			int check = 0;
			for (int j=0; j<vertices[i]->vertexNum; j++) {
				int k = j+1;
				if (k > 2) k = 0;
				if ( (jx*kx + jy*ky) > 0 ) check++;
			}

			if (check > 1) {
				double min_length = 1.0e10;
				int min_edge = -1;
				bool outside_check = false;
				for (int j=0; j<3; j++) {
					if (edges[vertices[i]->edgeIndex[j*2]]->distance < min_length) {
						min_length = edges[vertices[i]->edgeIndex[j*2]]->distance;
						min_edge = vertices[i]->edgeIndex[j*2];
						if (edges[vertices[i]->edgeIndex[j*2]]->cellNum < 2) outside_check = true;
					}
				}
				if (outside_check == false) edges[min_edge]->negative = -1;
			}
		}
	}

	virtual void updateEdges() {
		
		for (int i=0; i<edgesSize; i++) {

			int c, v1, v2;		

			c = edges[i]->cell[0];
			v1 = edges[i]->v1[0];
			v2 = edges[i]->v2[0];
		
			double cv1x = DCMcells[c]->DCMvertices[v1]->px;
			double cv2x = DCMcells[c]->DCMvertices[v2]->px;
	
			edges[i]->setDistance(cv1x, DCMcells[c]->DCMvertices[v1]->py, DCMcells[c]->DCMvertices[v1]->pz,
				cv2x, DCMcells[c]->DCMvertices[v2]->py, DCMcells[c]->DCMvertices[v2]->pz);	
		}
	
		int check;
		
		do {
			check = 0;

			for (int i=0; i<cellsSize; i++) {
				if (DCMcells[i]->ignoreFlag == 0) {

					double t2_threshold;

					if ((vertex_dynamics_flag & VD::T2TH_RELATIVE) == VD::T2TH_RELATIVE) t2_threshold = PM.Extrusion_Area*averageVolume;
					else {
						t2_threshold = PM.Extrusion_Area * initial_average_volume;
					}

					if (DCMcells[i]->out_cell_eli_flag == 0) {
						if (DCMcells[i]->cellArea[0] < t2_threshold && DCMcells[i]->fix_flag == 0) {
						
							if (DCMcells[i]->outsideFlag == 0) {
								if ((recordable_flag & RE::STRESS_RECORD) == RE::STRESS_RECORD) record_DivisionApoptosis(i);
							}

							if (baseExtrusion(i) == 0) {
								check = 1;						
								break;
							}
						}
					}
					else if (DCMcells[i]->out_cell_eli_flag == 1) {
						if (DCMcells[i]->cellArea[0] < initial_average_volume*0.2) {
							if (baseExtrusion(i) == 0) {
								check = 1;
								break;
							}
						}
					}
				}
			}
		} while(check == 1);
		
		
		do {
			check = 0;
		
			Tedge2D::ave_distance = 0.0;

			for (int i=0; i<edgesSize; i++) Tedge2D::ave_distance += edges[i]->distance;

			Tedge2D::ave_distance /= static_cast<double>(edgesSize);

			for (int i=0; i<edgesSize; i++) {

				double t1_threshold;
				t1_threshold = PM.Epsilon_Distance * Tedge2D::ave_distance * edges[i]->t1_coefficient;

				if (edges[i]->distance < t1_threshold && edges[i]->fix_flag == 0) {

					int ci1, ci2;

					ci1 = edges[i]->cell[0];
					ci2 = edges[i]->cell[1];

					if (ci2 > -1) {					
						if (DCMcells[ci1]->verticesSize_r < 4 && DCMcells[ci2]->verticesSize_r < 4) {}
					}		

					if (DCMcells[ci1]->verticesSize_r < 4) {
						baseApoptosis(i, 0);
						check = 1;
						break;
					}
					if (ci2 > -1) {
						if (DCMcells[ci2]->verticesSize_r < 4) {
							baseApoptosis(i, 1);
							check = 1;
							break;
						}
					}					

					baseReconnection(i, t1_threshold*1.2f);
					check = 1;
					break;
				}
			}
		} while(check == 1);
				
	}

	virtual void updateCellDivision(int flag=0) {
		
		for(int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				DCMcells[i]->setCellCycleCoefficient(flag);
				DCMcells[i]->updateCellCycle();
			}
		}
		
		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {

				if (DCMcells[i]->cellCycle > PM.CellCycleTime) {

					if ( DCMcells[i]->outsideFlag == 0 && DCMcells[i]->divArea < 0) {
						if ((recordable_flag & RE::STRESS_RECORD) == RE::STRESS_RECORD) record_DivisionStress(i);
					}
					
					if ((GrowthDivisionFlag & BC::MPHASE_2SI) == BC::MPHASE_2SI) {

						if (DCMcells[i]->divArea < 0) {
							DCMcells[i]->divArea = DCMcells[i]->cellArea[0];
							DCMcells[i]->setTargetArea(PM.Init_TargetArea*2.0);
						}

						if (DCMcells[i]->cellArea[0] > 2.0*DCMcells[i]->divArea) {							
							DCMcells[i]->setTargetArea(PM.Init_TargetArea);							
							DCMcells[i]->divArea = -1.0;
							baseDivision(i);
							updateEdges();
						}

					}

					else if ((GrowthDivisionFlag & BC::MPHASE_2SD) == BC::MPHASE_2SD) {

						if (DCMcells[i]->divArea < 0) {
							DCMcells[i]->divArea = DCMcells[i]->cellArea[0];
							DCMcells[i]->divCycle = 0.0;
						}

						DCMcells[i]->divCycle += PM.Delta_t;
						if (DCMcells[i]->targetArea < PM.Init_TargetArea*3.0) {
							DCMcells[i]->setTargetArea(PM.Init_TargetArea+DCMcells[i]->divCycle*mphase_coef);
						}

						if (DCMcells[i]->cellArea[0] > 2.0*DCMcells[i]->divArea) {							
							DCMcells[i]->divArea = -1.0;
							baseDivision(i);
							updateEdges();
						}

					}

					else if ((GrowthDivisionFlag & BC::MPHASE_ADD) == BC::MPHASE_ADD) {

						if (DCMcells[i]->divArea < 0) {
							DCMcells[i]->divArea = DCMcells[i]->cellArea[0];
							DCMcells[i]->divCycle = 0.0;
						}

						DCMcells[i]->divCycle += PM.Delta_t;
						if (DCMcells[i]->targetArea < PM.Init_TargetArea * 3.0) {
							DCMcells[i]->setTargetArea(PM.Init_TargetArea + DCMcells[i]->divCycle * mphase_coef);
						}

						if (DCMcells[i]->cellArea[0] > (DCMcells[i]->divArea + adder_parameter)) {
							DCMcells[i]->divArea = -1.0;
							baseDivision(i);
							updateEdges();
						}

					}


					else if ((GrowthDivisionFlag & BC::MPHASE_FSI) == BC::MPHASE_FSI) {

						if (DCMcells[i]->divArea < 0) {
							DCMcells[i]->divArea = DCMcells[i]->cellArea[0];
							DCMcells[i]->setTargetArea(PM.Init_TargetArea*3.0);
						}

						if (DCMcells[i]->cellArea[0] > PM.Init_TargetArea*2.0) {							
							DCMcells[i]->setTargetArea(PM.Init_TargetArea);							
							DCMcells[i]->divArea = -1.0;
							baseDivision(i);
							updateEdges();
						}
					}

					else if ((GrowthDivisionFlag & BC::MPHASE_FSD) == BC::MPHASE_FSD) {

						if (DCMcells[i]->divArea < 0) {
							DCMcells[i]->divArea = DCMcells[i]->cellArea[0];
							DCMcells[i]->divCycle = 0.0;
						}

						DCMcells[i]->divCycle += PM.Delta_t;
						if (DCMcells[i]->targetArea < PM.Init_TargetArea*3.0) {
							DCMcells[i]->setTargetArea(PM.Init_TargetArea+DCMcells[i]->divCycle*mphase_coef);
						}

						if (DCMcells[i]->cellArea[0] > PM.Init_TargetArea*1.5) {							
							DCMcells[i]->setTargetArea(PM.Init_TargetArea);							
							DCMcells[i]->divArea = -1.0;
							baseDivision(i);
							updateEdges();
						}

					}

					else if ((GrowthDivisionFlag & BC::MPHASE_TTH) == BC::MPHASE_TTH) {

						if (DCMcells[i]->divArea < 0) {
							DCMcells[i]->divArea = DCMcells[i]->cellArea[0];
							DCMcells[i]->divCycle = 0.0;
						}

						DCMcells[i]->divCycle += PM.Delta_t;
						if (DCMcells[i]->targetArea < PM.Init_TargetArea*2.0) {
							DCMcells[i]->setTargetArea(PM.Init_TargetArea+DCMcells[i]->divCycle*0.1);
						}

						if (DCMcells[i]->divCycle > 10.0) {
							DCMcells[i]->setTargetArea(PM.Init_TargetArea);							
							DCMcells[i]->divArea = -1.0;
							baseDivision(i);
							updateEdges();
						}

					}
										
					else {
						baseDivision(i);
						updateEdges();
					}
				}		
			}
		}		
	}


	virtual void updateTopology() {

		updateReindexVertex();
	
		containerVerticesSize = 0;
	
		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				containerVerticesSize += DCMcells[i]->verticesSize_r;
			}
		}	

		setEdges();
				
		setVertices();

		int cs = 0;
		
		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
			
				cs++;
		
				DCMcells[i]->setCellArea(1);
				DCMcells[i]->setPerimeterLength();
				DCMcells[i]->setOutside();

				DCMcells[i]->setOtherCellIndex();

			}
		}
			
		cellsSize_r = cs;
		
		setAverageVolume();

	}


	void updateReindexVertex() {

		for (int ci=0; ci<cellsSize; ci++) {

			if (DCMcells[ci]->ignoreFlag == 0) {

				int vs = DCMcells[ci]->verticesSize_r;
				int check = DCMcells[ci]->updateReindexVertex();

				if (check == 1) {
					
					for (int i=0; i<vs; i++) {
						
						int oc1 = DCMcells[ci]->DCMvertices[i]->othercell[0];
						int ov1 = DCMcells[ci]->DCMvertices[i]->othercell[1];
						int oc2 = DCMcells[ci]->DCMvertices[i]->othercell[2];
						int ov2 = DCMcells[ci]->DCMvertices[i]->othercell[3];
						
						for (int j=0; j<2; j++) {
							if (oc1 > -1) {
								if (ci == DCMcells[oc1]->DCMvertices[ov1]->othercell[2*j]) {
									DCMcells[oc1]->DCMvertices[ov1]->othercell[2*j+1] = i;
								}
							}
							if (oc2 > -1) {
								if (ci == DCMcells[oc2]->DCMvertices[ov2]->othercell[2*j]) {
									DCMcells[oc2]->DCMvertices[ov2]->othercell[2*j+1] = i;
								}
							}
						}
					}
				}
				else if (check == 2) {}
			}
		}

	}

	void resetCells() {

		if ((cellsSize-cellsSize_r) > 1) {

			if (recordable_flag > 0) recordable_flag = 0;

			Tcell2D** tempDCMcells = new Tcell2D*[cellsSize_r];
			for (int i=0; i<cellsSize_r; i++) tempDCMcells[i] = new Tcell2D;

			int nc_index = 0;

			for (int i=0; i<cellsSize; i++) {
				if (DCMcells[i]->ignoreFlag == 0) {
					*tempDCMcells[nc_index] = *DCMcells[i];
					nc_index++;
				}				
			}

			chPDelete(DCMcells, cellsSize);

			cellsSize = cellsSize_r;

			DCMcells = new Tcell2D*[cellsSize];
			for (int i=0; i<cellsSize; i++) DCMcells[i] = new Tcell2D;
			for (int i=0; i<cellsSize; i++) *DCMcells[i] = *tempDCMcells[i];

			chPDelete(tempDCMcells, cellsSize_r);

			setOtherCell();

			updateTopology();

		}

	}


	void resetEdges() {

		if (edgesArraySize < (containerVerticesSize+5)) {
				
			chPDelete(edges, edgesArraySize);

			edgesArraySize = containerVerticesSize*2;		

			edges = new Tedge2D*[edgesArraySize];

			for (int i=0; i<edgesArraySize; i++) {
				edges[i] = new Tedge2D;
			}

		}

		else {

			Tedge2D initEdge2D;

			for (int i=0; i<edgesArraySize; i++) {
				*edges[i] = initEdge2D;
			}
		}
	}


	virtual void setEdges() {
		
		resetEdges();

		int index = 0;
	
		for (int i=0; i<cellsSize; i++) {

			if (DCMcells[i]->ignoreFlag == 0) {
				
				for (int j=0; j<DCMcells[i]->verticesSize; j++) {
				
					if (DCMcells[i]->DCMvertices[j]->ignoreFlag == 0 &&
						DCMcells[i]->DCMvertices[j]->moveFlag == 0) {

						int n = DCMcells[i]->DCMvertices[j]->connection[1];
					
						edges[index]->cellNum = 1;

						edges[index]->cell[0] = i;
						edges[index]->v1[0] = j;
						edges[index]->v2[0] = n;
						edges[index]->t1_coefficient = 1;
						edges[index]->randomness_of_line_coef = DCMcells[i]->DCMvertices[j]->edge_randomness_of_line_coef;
					
						double jx = DCMcells[i]->DCMvertices[j]->px;
						double nx = DCMcells[i]->DCMvertices[n]->px;

						edges[index]->setDistance(jx, DCMcells[i]->DCMvertices[j]->py, DCMcells[i]->DCMvertices[j]->pz,
							nx, DCMcells[i]->DCMvertices[n]->py, DCMcells[i]->DCMvertices[n]->pz);
										
						DCMcells[i]->DCMvertices[j]->moveFlag = 1;
				
						for (int l=0; l<2; l++) {

							if (DCMcells[i]->DCMvertices[n]->othercell[2*l] > -1) {

								int oc = DCMcells[i]->DCMvertices[n]->othercell[2*l];
								int ov = DCMcells[i]->DCMvertices[n]->othercell[2*l+1];
								int opv = DCMcells[oc]->DCMvertices[ov]->connection[1];

								for (int ll=0; ll<2; ll++) {
									if (DCMcells[oc]->DCMvertices[opv]->othercell[2*ll] == i) {
										edges[index]->cellNum = 2;
										edges[index]->cell[1] = oc;
										edges[index]->v1[1] = opv;
										edges[index]->v2[1] = ov;
										DCMcells[oc]->DCMvertices[ov]->moveFlag = 1;
									}
								}
							}
						}
				
						index++;
					}
				}
			}
		}
	
		edgesSize = index;
		
		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				DCMcells[i]->setAllMoveFlagZero();										
			}
		}

		setEdges_CellIndexing();
		
		setCoefficient();
	}

	void setEdges_CellIndexing() {

		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				DCMcells[i]->resetEdgesNumber();
			}
		}	
	
		for (int i=0; i<edgesSize; i++) {
			for (int j=0; j<edges[i]->cellNum; j++) {
				for (int k=0; k<2; k++) {
					if (DCMcells[edges[i]->cell[j]]->DCMvertices[edges[i]->v1[j]]->edgesNumber[k] < 0) {
						DCMcells[edges[i]->cell[j]]->DCMvertices[edges[i]->v1[j]]->edgesNumber[k] = i;
						break;
					}				
				}
				for (int k=0; k<2; k++) {
					if (DCMcells[edges[i]->cell[j]]->DCMvertices[edges[i]->v2[j]]->edgesNumber[k] < 0) {
						DCMcells[edges[i]->cell[j]]->DCMvertices[edges[i]->v2[j]]->edgesNumber[k] = i;
						break;
					}				
				}
			}
		}		

	}

	virtual void setRandomnessLineCoef() {
		for (int i=0; i<edgesSize; i++) {			
			edges[i]->randomness_of_line_coef += -edges[i]->randomness_of_line_coef * PM.Decay_of_Lambda * PM.Delta_t 
				+ PM.Randomness_of_Lambda * SF.gaussDistribution(0,1) * sqrt(PM.Delta_t);
		}
	}


	void resetVertices() {

		if (verticesArraySize < (containerVerticesSize+5)) {
		
			chPDelete(vertices, verticesArraySize);		

			verticesArraySize = containerVerticesSize*2;

			vertices = new Tpoint2D*[verticesArraySize];

			for (int i=0; i<verticesArraySize; i++) {
				vertices[i] = new Tpoint2D;
			}
		}

		else {

			Tpoint2D initVertex;
		
			for (int i=0; i<verticesArraySize; i++) {
				*vertices[i] = initVertex;
			}
		}
	}

	
	virtual void setVertices() {

		resetVertices();

		int index = 0;

		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
		
				for (int j=0; j<DCMcells[i]->verticesSize; j++) {
			
					vertices[index]->cellIndex[0] = i;
									
					if (DCMcells[i]->DCMvertices[j]->ignoreFlag == 0 && DCMcells[i]->DCMvertices[j]->moveFlag == 0) {

						vertices[index]->vertexNum = 1;
						vertices[index]->vertexIndex[0] = j;

						vertices[index]->x = DCMcells[i]->DCMvertices[j]->px;
						vertices[index]->y = DCMcells[i]->DCMvertices[j]->py;
						vertices[index]->z = DCMcells[i]->DCMvertices[j]->pz;
				
						for (int l=0; l<2; l++) {					
					
							int oc = DCMcells[i]->DCMvertices[j]->othercell[2*l];
							int ov = DCMcells[i]->DCMvertices[j]->othercell[2*l+1];
					
							if (oc > -1) {

								int ll = vertices[index]->vertexNum;
								vertices[index]->vertexNum++;		

								vertices[index]->cellIndex[ll] = oc;
								vertices[index]->vertexIndex[ll] = ov;

								DCMcells[oc]->DCMvertices[ov]->moveFlag = 1;
							}
						}

						index++;
					}

					DCMcells[i]->DCMvertices[j]->moveFlag = 1;

				}
			}
		}

		verticesSize = index;
			
		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				DCMcells[i]->setAllMoveFlagZero();
			}
		}
	
		setVertices_CellIndexing();
		setVertices_EdgeIndexing();

	}

	void setVertices_CellIndexing() {

		for (int i=0; i<verticesSize; i++) {
			for (int j=0; j<vertices[i]->vertexNum; j++) {
				DCMcells[vertices[i]->cellIndex[j]]->DCMvertices[vertices[i]->vertexIndex[j]]->verticesNumber = i;
			}
		}
	}

	void setVertices_EdgeIndexing() {

		for (int i=0; i<verticesSize; i++) {

			vertices[i]->edgeIndex[0] = DCMcells[vertices[i]->cellIndex[0]]->DCMvertices[vertices[i]->vertexIndex[0]]->edgesNumber[0];
			if (edges[vertices[i]->edgeIndex[0]]->v1[0] == vertices[i]->vertexIndex[0]) {
				vertices[i]->edgeIndex[1] = 0;
			}
			else {
				vertices[i]->edgeIndex[1] = 1;
			}

			vertices[i]->edgeIndex[2] = DCMcells[vertices[i]->cellIndex[0]]->DCMvertices[vertices[i]->vertexIndex[0]]->edgesNumber[1];
			if (edges[vertices[i]->edgeIndex[2]]->v1[0] == vertices[i]->vertexIndex[0]) {
				vertices[i]->edgeIndex[3] = 0;
			}
			else {
				vertices[i]->edgeIndex[3] = 1;
			}

			if (vertices[i]->cellIndex[1] > -1) {			

				for (int j=0; j<2; j++) {
					int oe = DCMcells[vertices[i]->cellIndex[1]]->DCMvertices[vertices[i]->vertexIndex[1]]->edgesNumber[j];
					if (vertices[i]->edgeIndex[0] != oe && vertices[i]->edgeIndex[1] != oe) {
						vertices[i]->edgeIndex[4] = oe;
						if (edges[oe]->cell[0] == vertices[i]->cellIndex[1]) {
							if (edges[oe]->v1[0] == vertices[i]->vertexIndex[1]) {
								vertices[i]->edgeIndex[5] = 0;
							}
							else {
								vertices[i]->edgeIndex[5] = 1;
							}
						}
						else if (edges[oe]->cell[1] == vertices[i]->cellIndex[1]) {
							if (edges[oe]->v1[1] == vertices[i]->vertexIndex[1]) {
								vertices[i]->edgeIndex[5] = 0;
							}
							else {
								vertices[i]->edgeIndex[5] = 1;
							}
						}
					}
				}			
			}
		}	
	}


	virtual int vDynamics_loop_base(int i, DCMpoint& xy) {

		int error;
		int check = 0;
		
		if ((vertex_dynamics_flag & VD::OUT_FIX) == VD::OUT_FIX) {
			if (vertices[i]->vertexNum < 3) check = 1;
		}
		
		if (check == 0) {

			for(int j=0; j<3; j++) {

				if (vertices[i]->edgeIndex[2*j] > -1) {

					if ((vertex_dynamics_flag & VD::LINETENSION) == VD::LINETENSION) {
						error = edges[vertices[i]->edgeIndex[2*j]]->vDynamics_line(vertices[i]->edgeIndex[2*j+1], xy);
						if (error != 0) return error;
					}

					if ((vertex_dynamics_flag & VD::LINESPRING) == VD::LINESPRING) {
						error = edges[vertices[i]->edgeIndex[2*j]]->vDynamics_linespring(vertices[i]->edgeIndex[2*j+1], xy);
						if (error != 0) return error;
					}

					if ((vertex_dynamics_flag & VD::LINEELASTICITY) == VD::LINEELASTICITY) {
						error = edges[vertices[i]->edgeIndex[2*j]]->vDynamics_lineelasticity(vertices[i]->edgeIndex[2*j+1], xy);
						if (error != 0) return error;
					}
				}
			}

			if ((vertex_dynamics_flag & VD::PERIMETER) == VD::PERIMETER) {
				for(int j=0; j<vertices[i]->vertexNum; j++) {
					error = DCMcells[vertices[i]->cellIndex[j]]->vDynamics_perimeter(vertices[i]->vertexIndex[j], xy);
					if (error != 0) return error;
				}
			}

			if ((vertex_dynamics_flag & VD::ANGLESPRING) == VD::ANGLESPRING) {
				for(int j=0; j<vertices[i]->vertexNum; j++) {
					error = DCMcells[vertices[i]->cellIndex[j]]->vDynamics_angle(vertices[i]->vertexIndex[j], xy);
					if (error != 0) return error;
				}
				if (vertices[i]->vertexNum < 3) {
					if (error != 0) return error;
				}
			}
			
			for(int j=0; j<vertices[i]->vertexNum; j++) {
				error = DCMcells[vertices[i]->cellIndex[j]]->vDynamics_area(vertices[i]->vertexIndex[j], xy);
				if (error != 0) return error;
			}
		}

		return 0;

	}

	void movingVertex() {
		if ((GrowthDivisionFlag & BC::MOVING_VERTEX) == 0) {
			GrowthDivisionFlag ^= BC::MOVING_VERTEX;
		}
	}

	virtual void setCoefficient(int flag=0) {
		
		for (int i=0; i<edgesSize; i++) edges[i]->setLineCoefficient();

		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				DCMcells[i]->setPerimeterCoefficient();
				DCMcells[i]->setAreaCoefficient();
			}
		}

	}

	virtual int vDynamics_loop(size_t i) {

		DCMpoint xy;
		int error;

		vDynamics_loop_base(i, xy);

		if ((vertex_dynamics_flag & VD::RANDOMMOVE) == VD::RANDOMMOVE) {
			xy.x += PM.RandomMove*SF.gaussDistribution(0,1.0);
			xy.y += PM.RandomMove*SF.gaussDistribution(0,1.0);
		}
	
		vertices[i]->vx = xy.x * PM.Eta;
		vertices[i]->vy = xy.y * PM.Eta;
		vertices[i]->vz = xy.z * PM.Eta;

		xy.x *= PM.Delta_t*PM.Eta;
		xy.y *= PM.Delta_t*PM.Eta;
		xy.z *= PM.Delta_t*PM.Eta;
		
		vertices[i]->x += xy.x;
		vertices[i]->y += xy.y;
		vertices[i]->z += xy.z;
		
		vertices[i]->dxdt = sqrt3(xy.x, xy.y, xy.z);
		
		if (i == 180 && (GrowthDivisionFlag & BC::MOVING_VERTEX) == BC::MOVING_VERTEX) {
			vertices[i]->x += 0.01;
			vertices[i]->y += 0.01;
			PM.Delta_t *= 0.001;
			GrowthDivisionFlag ^= BC::MOVING_VERTEX;
		}
		
		return 0;
	
	}

	int vDynamics_angle_outside(int i, DCMpoint& xy) {

		int out = 0;

		for (int j=0; j<2; j++) {
			if (DCMcells[vertices[i]->cellIndex[0]]->DCMvertices[DCMcells[vertices[i]->cellIndex[0]]->DCMvertices[vertices[i]->vertexIndex[0]]->connection[0]]->othercell[2*j] < 0) {
				out = 1;
				break;
			}
		}

		int c1, c2, pp, pn;
	
		double id_pi, id_ni;

		if (out == 0) {
			c1 = vertices[i]->cellIndex[0];
			pp = DCMcells[c1]->DCMvertices[vertices[i]->vertexIndex[0]]->connection[1];
			c2 = vertices[i]->cellIndex[1];
			pn = DCMcells[c2]->DCMvertices[vertices[i]->vertexIndex[1]]->connection[0];

			id_pi = 1.0/DCMcells[c1]->DCMvertices[vertices[i]->vertexIndex[0]]->distance[0];
			id_ni = 1.0/DCMcells[c2]->DCMvertices[pn]->distance[0];
		}
		else {
			c1 = vertices[i]->cellIndex[1];
			pp = DCMcells[c1]->DCMvertices[vertices[i]->vertexIndex[1]]->connection[0];
			c2 = vertices[i]->cellIndex[0];
			pn = DCMcells[c2]->DCMvertices[vertices[i]->vertexIndex[0]]->connection[1];

			id_pi = 1.0/DCMcells[c1]->DCMvertices[vertices[i]->vertexIndex[1]]->distance[0];
			id_ni = 1.0/DCMcells[c2]->DCMvertices[pn]->distance[0];
		}
		
		double k = PM.FlatK;

		double x_pi = DCMcells[c1]->DCMvertices[pp]->px - vertices[i]->x;
		double y_pi = DCMcells[c1]->DCMvertices[pp]->py - vertices[i]->y;
	
		double x_ni = DCMcells[c2]->DCMvertices[pn]->px - vertices[i]->x;
		double y_ni = DCMcells[c2]->DCMvertices[pn]->py - vertices[i]->y;
	
		double s_x = x_pi*id_pi*id_pi + x_ni*id_ni*id_ni;
		double s_y = y_pi*id_pi*id_pi + y_ni*id_ni*id_ni;

		xy.x += -k* id_pi*id_ni* (0.25*(-x_pi-x_ni + s_x*(x_pi*x_ni+y_pi*y_ni)) + sqrt(3.0)*0.25*( y_pi-y_ni + s_x*(x_pi*y_ni-y_pi*x_ni)));
		xy.y += -k* id_pi*id_ni* (0.25*(-y_pi-y_ni + s_y*(x_pi*x_ni+y_pi*y_ni)) + sqrt(3.0)*0.25*(-x_pi+x_ni + s_y*(x_pi*y_ni-y_pi*x_ni)));
	
		return 0;

	}


	virtual int vDynamicsMain() {

		conTimeStep++;

		if ((vertex_dynamics_flag & VD::RANDOMLAMBDA) == VD::RANDOMLAMBDA) {
			setRandomnessLineCoef();
			setCellvertex_EdgeLinecoef();
		}

		setCoefficient();

		int error = 0;
		for (int i=0; i<verticesSize; i++) {
			error = vDynamics_loop((size_t)i);
		}
		if (error != 0) return error;

		maxdxdt = 0;

		for (int i=0; i<verticesSize; i++) {

			for (int j=0; j<vertices[i]->vertexNum; j++) {
				DCMcells[vertices[i]->cellIndex[j]]->DCMvertices[vertices[i]->vertexIndex[j]]->px = vertices[i]->x;
				DCMcells[vertices[i]->cellIndex[j]]->DCMvertices[vertices[i]->vertexIndex[j]]->py = vertices[i]->y;
				DCMcells[vertices[i]->cellIndex[j]]->DCMvertices[vertices[i]->vertexIndex[j]]->pz = vertices[i]->z;
			}

			if (_isnan(vertices[i]->dxdt)) return 341;

			if (vertices[i]->dxdt > maxdxdt) maxdxdt = vertices[i]->dxdt;
		}

		if (maxdxdt > maxmaxdxdt) maxmaxdxdt = maxdxdt;


		updateEdges();

		if (cellsSize_r < maxCellsSize && (GrowthDivisionFlag & BC::NO_DIVISION) == 0) updateCellDivision();

		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {

				DCMcells[i]->setCellArea();
				DCMcells[i]->setPerimeterLength();				
	
				if (conTimeStep%5 == 0) {
					if ((recordable_flag & RE::TOPOLOGY_CHANGE) == RE::TOPOLOGY_CHANGE) DCMcells[i]->updateRecord(i, conTimeStep/100);
				}

				DCMcells[i]->updateTopologyFlag();				
			}
		}

		if ((recordable_flag & RE::STRESS_RECORD) == RE::STRESS_RECORD) {
			record_DivisionStress();
			record_DivisionApoptosis();
		}
		
		return 0;

	}

	
	virtual int baseReconnection(int edgeIndex, double epsilon, int flag=0) {

		int error = 0;

		int c_b = edges[edgeIndex]->cell[0];
		int v1_b = edges[edgeIndex]->v1[0];
		int v2_b = edges[edgeIndex]->v2[0];

		int c_t = edges[edgeIndex]->cell[1];
		int v1_t = edges[edgeIndex]->v2[1];
		int v2_t = edges[edgeIndex]->v1[1];

		int c_r = -1;
		int v1_r = -1;
		int v2_r = -1;

		int c_l = -1;
		int v1_l = -1;
		int v2_l = -1;

		for (int i=0; i<2; i++) {
		
			if (DCMcells[c_b]->DCMvertices[v1_b]->othercell[2*i] != c_t) {
				c_r = DCMcells[c_b]->DCMvertices[v1_b]->othercell[2*i];
				v1_r = DCMcells[c_b]->DCMvertices[v1_b]->othercell[2*i+1];
			}
		
			if (DCMcells[c_b]->DCMvertices[v2_b]->othercell[2*i] != c_t) {
				c_l = DCMcells[c_b]->DCMvertices[v2_b]->othercell[2*i];
				v1_l = DCMcells[c_b]->DCMvertices[v2_b]->othercell[2*i+1];
			}		
		}

		double v1_x = edges[edgeIndex]->v1x;
		double v1_y = edges[edgeIndex]->v1y;

		double v2_x = edges[edgeIndex]->v2x;
		double v2_y = edges[edgeIndex]->v2y;

		double v1_z = edges[edgeIndex]->v1z;
		double v2_z = edges[edgeIndex]->v2z;
	
		double c_x = 0.5*(v1_x+v2_x);
		double c_y = 0.5*(v1_y+v2_y);
		double c_z = 0.5*(v1_z+v2_z);
		
		double cd = 1.0/sqrt3((v1_x-c_x), (v1_y-c_y), (v1_z-c_z));

		double nor1_x = (v1_x-c_x)*cd;
		double nor1_y = (v1_y-c_y)*cd;

		double np2_x = (-nor1_y*PM.Const_ItoH*epsilon)+c_x;
		double np2_y = (nor1_x*PM.Const_ItoH*epsilon)+c_y;
		double np1_x = (nor1_y*PM.Const_ItoH*epsilon)+c_x;
		double np1_y = (-nor1_x*PM.Const_ItoH*epsilon)+c_y;

		double np1_z = 0;
		double np2_z = 0;

		if (flag == 0) {
			if ((recordable_flag & RE::RECON_RECORD) == RE::RECON_RECORD) record_Reconnection(c_t, c_b, c_r, c_l, 0, c_x, c_y);
		}

		reconnection_tbCell(c_b, v1_b, v2_b, np1_x, np1_y, np1_z);
		if (c_t > -1) reconnection_tbCell(c_t, v1_t, v2_t, np2_x, np2_y, np2_z);		
		if (c_r > -1) reconnection_rlCell(c_r, v1_r, v2_r, np1_x, np1_y, np1_z, np2_x, np2_y, np2_z);
		if (c_l > -1) reconnection_rlCell(c_l, v1_l, v2_l, np2_x, np2_y, np2_z, np1_x, np1_y, np1_z);
	

		DCMcells[c_b]->setVertexOther(v1_b, c_r, v2_r, c_l, v1_l);

		if (c_t > -1) {
			DCMcells[c_t]->setVertexOther(v1_t, c_r, v1_r, c_l, v2_l);
		}

		if (c_r > -1) {
			DCMcells[c_r]->setVertexOther(v1_r, c_t, v1_t, c_l, v2_l);
			DCMcells[c_r]->setVertexOther(v2_r, c_b, v1_b, c_l, v1_l);
		}

		if (c_l > -1) {
			DCMcells[c_l]->setVertexOther(v1_l, c_b, v1_b, c_r, v2_r);
			DCMcells[c_l]->setVertexOther(v2_l, c_t, v1_t, c_r, v1_r);
		}

		if (flag == 0) {
			if ((recordable_flag & RE::RECON_RECORD) == RE::RECON_RECORD) record_Reconnection(c_t, c_b, c_r, c_l, 1);
			recon_number_r++;
		}

		recon_number++;

		if ((recordable_flag & RE::TOPOLOGY_CHANGE) == RE::TOPOLOGY_CHANGE) {
	
			DCMcells[c_b]->setCellArea(1);
			DCMcells[c_b]->updateRecord(c_b, conTimeStep/100, US::RECONNECTION);

			if (c_t > -1) {
				DCMcells[c_t]->setCellArea(1);
				DCMcells[c_t]->updateRecord(c_t, conTimeStep/100, US::RECONNECTION);
			}

			if (c_r > -1) {
				DCMcells[c_r]->setCellArea(1);
				DCMcells[c_r]->updateRecord(c_r, conTimeStep/100, US::RECONNECTION);
			}

			if (c_l > -1) {
				DCMcells[c_l]->setCellArea(1);
				DCMcells[c_l]->updateRecord(c_l, conTimeStep/100, US::RECONNECTION);
			}
		}

		updateTopology();

		return 0;
	
	}
	
	virtual int reconnection_tbCell(int& c, int& v1, int& v2, double& np_x, double& np_y, double& np_z) {

		DCMcells[c]->DCMvertices[v1]->px = np_x;
		DCMcells[c]->DCMvertices[v1]->py = np_y;
		DCMcells[c]->DCMvertices[v1]->pz = np_z;
	
		int pv1 = DCMcells[c]->DCMvertices[v1]->connection[0];
		int nv2 = DCMcells[c]->DCMvertices[v2]->connection[1];

		DCMcells[c]->DCMvertices[v1]->edge_randomness_of_line_coef = DCMcells[c]->DCMvertices[v2]->edge_randomness_of_line_coef;

		DCMcells[c]->ignoreVertex(v2);

		DCMcells[c]->DCMvertices[v1]->connection[1] = nv2;
		DCMcells[c]->DCMvertices[nv2]->connection[0] = v1;
		
		return 0;
	}

	virtual int reconnection_rlCell(int& c, int& v1, int& v2, double& np1_x, double& np1_y, double& np1_z, double& np2_x, double& np2_y, double& np2_z) {

		int error = 0;

		int nv1 = DCMcells[c]->DCMvertices[v1]->connection[1];
		
		v2 = DCMcells[c]->verticesSize;
		DCMcells[c]->verticesSize++;
		
		error = DCMcells[c]->setNewVertex(v2, np1_x, np1_y, np1_z, v1, nv1);
		if (error != 0) return error;

		DCMcells[c]->DCMvertices[v2]->edge_randomness_of_line_coef = DCMcells[c]->DCMvertices[v1]->edge_randomness_of_line_coef;
		setNewCellvertex_EdgeLineCoef(c, v1);
		
		DCMcells[c]->DCMvertices[v1]->px = np2_x;
		DCMcells[c]->DCMvertices[v1]->py = np2_y;
		DCMcells[c]->DCMvertices[v1]->pz = np2_z;
		
		DCMcells[c]->DCMvertices[v1]->connection[1] = v2;
		DCMcells[c]->DCMvertices[nv1]->connection[0] = v2;
		
		return 0;
	}


	virtual int baseApoptosis(int edgeIndex, int tbcellFlag) {

		int c_b = edges[edgeIndex]->cell[0];
		int v1_b = edges[edgeIndex]->v1[0];
		int v2_b = edges[edgeIndex]->v2[0];

		int c_t = edges[edgeIndex]->cell[1];
		int v1_t = edges[edgeIndex]->v2[1];
		int v2_t = edges[edgeIndex]->v1[1];

		int c_r = -1;
		int v1_r = -1;
		int v2_r = -1;

		int c_l = -1;
		int v1_l = -1;
		int v2_l = -1;

		for (int i=0; i<2; i++) {
		
			if (DCMcells[c_b]->DCMvertices[v1_b]->othercell[2*i] != c_t) {
				c_r = DCMcells[c_b]->DCMvertices[v1_b]->othercell[2*i];
				v1_r = DCMcells[c_b]->DCMvertices[v1_b]->othercell[2*i+1];
			}
		
			if (DCMcells[c_b]->DCMvertices[v2_b]->othercell[2*i] != c_t) {
				c_l = DCMcells[c_b]->DCMvertices[v2_b]->othercell[2*i];
				v1_l = DCMcells[c_b]->DCMvertices[v2_b]->othercell[2*i+1];
			}		
		}

		double v1_x = DCMcells[c_b]->DCMvertices[v1_b]->px;
		double v1_y = DCMcells[c_b]->DCMvertices[v1_b]->py;

		double v2_x = DCMcells[c_b]->DCMvertices[v2_b]->px;
		double v2_y = DCMcells[c_b]->DCMvertices[v2_b]->py;

		double v1_z = DCMcells[c_b]->DCMvertices[v1_b]->pz;
		double v2_z = DCMcells[c_b]->DCMvertices[v2_b]->pz;
		
		double c_x = 0.5*(v1_x+v2_x);
		double c_y = 0.5*(v1_y+v2_y);
		double c_z = 0.5*(v1_z+v2_z);
	
		if (tbcellFlag == 0) {		

			if (c_t > -1) {
				apoptosis_Cell(c_t, v1_t, v2_t, c_x, c_y, c_z, 1);				
				DCMcells[c_t]->setVertexOther(v1_t, c_r, v1_r, c_l, v1_l);				
			}
			if (c_r > -1) {
				v2_r = DCMcells[c_r]->DCMvertices[v1_r]->connection[1];
				apoptosis_Cell(c_r, v1_r, v2_r, c_x, c_y, c_z, 1);				
				DCMcells[c_r]->setVertexOther(v1_r, c_t, v1_t, c_l, v1_l);				
			}
			if (c_l > -1) {
				v2_l = DCMcells[c_l]->DCMvertices[v1_l]->connection[0];
				apoptosis_Cell(c_l, v1_l, v2_l, c_x, c_y, c_z, 0);				
				DCMcells[c_l]->setVertexOther(v1_l, c_t, v1_t, c_r, v1_r);				
			}
				
			DCMcells[c_b]->ignoreFlag = conTimeStep/100+1;
			
			if (DCMcells[c_b]->out_cell_eli_flag == 0) extrusion_number++;
			else DCMcells[c_b]->ignoreFlag = -DCMcells[c_b]->ignoreFlag;
	
		}
		else if (tbcellFlag == 1) {

			apoptosis_Cell(c_b, v1_b, v2_b, c_x, c_y, c_z, 1);			
			DCMcells[c_b]->setVertexOther(v1_b, c_r, v1_r, c_l, v1_l);
			
			if (c_r > -1) {
				v2_r = DCMcells[c_r]->DCMvertices[v1_r]->connection[0];
				apoptosis_Cell(c_r, v1_r, v2_r, c_x, c_y, c_z, 0);				
				DCMcells[c_r]->setVertexOther(v1_r, c_b, v1_b, c_l, v1_l);				
			}
			if (c_l > -1) {
				v2_l = DCMcells[c_l]->DCMvertices[v1_l]->connection[1];
				apoptosis_Cell(c_l, v1_l, v2_l, c_x, c_y, c_z, 1);
				DCMcells[c_l]->setVertexOther(v1_l, c_b, v1_b, c_r, v1_r);				
			}

			DCMcells[c_t]->ignoreFlag = conTimeStep/100+1;

			if (DCMcells[c_t]->out_cell_eli_flag == 0) extrusion_number++;
			else DCMcells[c_t]->ignoreFlag = -DCMcells[c_t]->ignoreFlag;

		}

		if ((recordable_flag & RE::TOPOLOGY_CHANGE) == RE::TOPOLOGY_CHANGE) {

			if (c_t > -1) {
				if (DCMcells[c_t]->ignoreFlag != 0) {
					if (c_b > -1) {
						DCMcells[c_b]->setCellArea(1);
						DCMcells[c_b]->updateRecord(c_b, conTimeStep/100, US::EXTRUSION);
					}
				}
			}

			if (c_b > -1) {
				if (DCMcells[c_b]->ignoreFlag != 0) {
					if (c_t > -1) {
						DCMcells[c_t]->setCellArea(1);
						DCMcells[c_t]->updateRecord(c_t, conTimeStep/100, US::EXTRUSION);
					}
				}
			}

			if (c_r > -1) {
				DCMcells[c_r]->setCellArea(1);
				DCMcells[c_r]->updateRecord(c_r, conTimeStep/100, US::EXTRUSION);
			}

			if (c_l > -1) {
				DCMcells[c_l]->setCellArea(1);
				DCMcells[c_l]->updateRecord(c_l, conTimeStep/100, US::EXTRUSION);
			}
		}

		updateTopology();
		
		return 0;
	
	}

	virtual int apoptosis_Cell(int& c, int& v1, int& v2, double& c_x, double& c_y, double& c_z, int flag) {

		int pn1=0, pn2=1;

		if (flag == 1) {
			pn1=1;
			pn2=0;
		}

		DCMcells[c]->DCMvertices[v1]->px = c_x;
		DCMcells[c]->DCMvertices[v1]->py = c_y;
		DCMcells[c]->DCMvertices[v1]->pz = c_z;
	
		int pn = DCMcells[c]->DCMvertices[v2]->connection[pn1];
		DCMcells[c]->DCMvertices[v1]->connection[pn1] = pn;
		DCMcells[c]->DCMvertices[pn]->connection[pn2] = v1;

		if (flag == 1) DCMcells[c]->DCMvertices[v1]->edge_randomness_of_line_coef = DCMcells[c]->DCMvertices[v2]->edge_randomness_of_line_coef;

		DCMcells[c]->ignoreVertex(v2);

		return 0;
	}

	

	int baseExtrusion(int cellIndex) {

		int check;

		do {

			check = 0;

			double minLength = 100000000.0;
			int minLengthEdge = -1;

			for (int i=0; i<DCMcells[cellIndex]->verticesSize; i++) {
				if (DCMcells[cellIndex]->DCMvertices[i]->ignoreFlag == 0) {

					int n = DCMcells[cellIndex]->DCMvertices[i]->connection[1];

					double distance = dist(DCMcells[cellIndex]->DCMvertices[i]->px, DCMcells[cellIndex]->DCMvertices[n]->px,
						DCMcells[cellIndex]->DCMvertices[i]->py, DCMcells[cellIndex]->DCMvertices[n]->py, 0, 0);

					if (distance < minLength) {

						if (DCMcells[cellIndex]->DCMvertices[i]->ocIndex > -1) {
						
							if (DCMcells[DCMcells[cellIndex]->DCMvertices[i]->ocIndex]->verticesSize_r > 4) {
								minLength = distance;
								minLengthEdge = i;
							}							
						}
					}
				}
			}

			int min_edge = search_Cellvertex_to_Edge(cellIndex, minLengthEdge);
	
			if (DCMcells[cellIndex]->verticesSize_r > 3) {
				double t1_threshold = PM.Epsilon_Distance * Tedge2D::ave_distance;
				baseReconnection(min_edge, t1_threshold*1.2, 1);
				check = 1;
			}
			else if (DCMcells[cellIndex]->verticesSize_r == 3) {
				if (edges[min_edge]->cell[0] == cellIndex) baseApoptosis(min_edge, 0);
				else baseApoptosis(min_edge, 1);
			}

		} while (check == 1);
		
		return 0;

	}

	int search_Cellvertex_to_Edge(int cell_index, int vertex_index) {

		int target_edge;		
		int e1, e2;
		int next_index = DCMcells[cell_index]->DCMvertices[vertex_index]->connection[1];

		for (int i=0; i<2; i++) {
			for (int j=0; j<2; j++) {				
				
				e1 = DCMcells[cell_index]->DCMvertices[vertex_index]->edgesNumber[i];
				e2 = DCMcells[cell_index]->DCMvertices[next_index]->edgesNumber[j];

				if (e1 == e2) target_edge = e1;
			}
		}
		return target_edge;
	}

	void setCellvertex_EdgeLinecoef() {

		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0) {
				for (int j=0; j<DCMcells[i]->verticesSize; j++) {
					if (DCMcells[i]->DCMvertices[j]->ignoreFlag == 0) {
						int edge_index = search_Cellvertex_to_Edge(i,j);
						DCMcells[i]->DCMvertices[j]->edge_randomness_of_line_coef = edges[edge_index]->randomness_of_line_coef;
					}
				}
			}
		}
	}


	int division_ellipse(int cellIndex, double& divAngle) {

		DCMcells[cellIndex]->setCellShape();
		divAngle = DCMcells[cellIndex]->cell_angle;

		if (vonMisesKappa > 0) {
			double da = divAngle;
			divAngle = SF.vonmisesDistribution(da, vonMisesKappa);			
		}

		return 0;

	}

	virtual int division_random(int cellIndex, double& divAngle) {

		divAngle = 2.0*Pi*SF.uniformDistribution(0,1);

		return 0;

	}
	
	virtual int division_angle(int cellIndex, int* divp, double* divpx, double* divpy, double* divpz) {

		double divAngle;
		
		if (vonMisesKappa == 0) division_random(cellIndex, divAngle);
		else division_ellipse(cellIndex, divAngle);
		
		double tdivpx[2];
		Tcell2D* mCell = new Tcell2D(*DCMcells[cellIndex]);
		
		for (int i=0; i<mCell->verticesSize; i++) {
			if (mCell->DCMvertices[i]->ignoreFlag == 0) {
							
				double tempx = DCMcells[cellIndex]->DCMvertices[i]->px - DCMcells[cellIndex]->cofm_x;
				double tempy = DCMcells[cellIndex]->DCMvertices[i]->py - DCMcells[cellIndex]->cofm_y;

				mCell->DCMvertices[i]->px = cos(divAngle)*tempx + sin(divAngle)*tempy;
				mCell->DCMvertices[i]->py = -sin(divAngle)*tempx + cos(divAngle)*tempy;
			}
		}

		int k=0;
		
		for (int i=0; i<mCell->verticesSize; i++) {

			if (mCell->DCMvertices[i]->ignoreFlag == 0) {

				int n = mCell->DCMvertices[i]->connection[1];

				if ((mCell->DCMvertices[i]->py*mCell->DCMvertices[n]->py) < 0) {
		
					if (k > 1) return 1;
					else {
						divp[k] = i;
						tdivpx[k] = mCell->DCMvertices[i]->px;						
						if (abs(mCell->DCMvertices[n]->px - mCell->DCMvertices[i]->px) > VERY_SMALL) {
							tdivpx[k] = mCell->DCMvertices[i]->px - mCell->DCMvertices[i]->py
								*(mCell->DCMvertices[n]->px - mCell->DCMvertices[i]->px)
								/(mCell->DCMvertices[n]->py - mCell->DCMvertices[i]->py);
						}
						else return 1;
						k++;
					}
				}
			}
		}

		if (k < 2) return 1;
		
		for (int i=0; i<2; i++) {
			divpx[i] = cos(divAngle)*tdivpx[i] + DCMcells[cellIndex]->cofm_x;
			divpy[i] = sin(divAngle)*tdivpx[i] + DCMcells[cellIndex]->cofm_y;			
		}
		
		delete mCell;

		return 0;

	}

	void setNewCellvertex_EdgeLineCoef(int cell_index, int vertex_index) {
		if ((vertex_dynamics_flag & VD::RANDOMLAMBDA) == VD::RANDOMLAMBDA &&
			(vertex_dynamics_flag & VD::RANDLAMBDA_INIT) == VD::RANDLAMBDA_INIT && 
			PM.Decay_of_Lambda > 0) {
			DCMcells[cell_index]->DCMvertices[vertex_index]->edge_randomness_of_line_coef =
				SF.gaussDistribution(0,1)*PM.Randomness_of_Lambda/sqrt(2.0*PM.Decay_of_Lambda);
		}
		else {
			DCMcells[cell_index]->DCMvertices[vertex_index]->edge_randomness_of_line_coef = 0;
		}		
	}

	virtual int baseDivision(int cellIndex) {
	
		int error = 0;

		int cc = 0;

		int divp[2];
		double divpx[2];
		double divpy[2];
		double divpz[2];

		divpz[0] = 0; divpz[1] = 0;

		division_angle(cellIndex, divp, divpx, divpy, divpz);
	
		int c_l = -1;
		int v1_l = -1;
		int v2_l = -1;
		int nLIndex = -1;	

		int c_r = -1;
		int v1_r = -1;
		int v2_r = -1;
		int nRIndex = -1;
	
		for (int i=0; i<2; i++) {
			for (int j=0; j<2; j++) {
				int n = DCMcells[cellIndex]->DCMvertices[divp[0]]->connection[1];
				if (DCMcells[cellIndex]->DCMvertices[divp[0]]->othercell[2*i] ==
					DCMcells[cellIndex]->DCMvertices[n]->othercell[2*j]) {
					c_r = DCMcells[cellIndex]->DCMvertices[divp[0]]->othercell[2*i];
					v1_r = DCMcells[cellIndex]->DCMvertices[divp[0]]->othercell[2*i+1];
					v2_r = DCMcells[cellIndex]->DCMvertices[n]->othercell[2*j+1];
				}
				n = DCMcells[cellIndex]->DCMvertices[divp[1]]->connection[1];
				if (DCMcells[cellIndex]->DCMvertices[divp[1]]->othercell[2*i] ==
					DCMcells[cellIndex]->DCMvertices[n]->othercell[2*j]) {
					c_l = DCMcells[cellIndex]->DCMvertices[divp[1]]->othercell[2*i];
					v1_l = DCMcells[cellIndex]->DCMvertices[divp[1]]->othercell[2*i+1];
					v2_l = DCMcells[cellIndex]->DCMvertices[n]->othercell[2*j+1];
				}				
			}
		}


		if (c_r > -1) { 
			nRIndex = DCMcells[c_r]->verticesSize;
			DCMcells[c_r]->verticesSize++;
		}
		if (c_l > -1) {
			nLIndex = DCMcells[c_l]->verticesSize;
			DCMcells[c_l]->verticesSize++;
		}


		int dCellVSize = 0;
		int* next = new int[DCMcells[cellIndex]->verticesSize+2];
		next[0] = divp[0];

		
		do {
			dCellVSize++;
			next[dCellVSize] = DCMcells[cellIndex]->DCMvertices[next[dCellVSize-1]]->connection[1];
		
		} while(next[dCellVSize] != divp[1]);
		
		dCellVSize++;
		next[dCellVSize] = DCMcells[cellIndex]->verticesSize;
		next[0] = DCMcells[cellIndex]->verticesSize+1;
		dCellVSize++;

		int dCell = cellsSize;
		cellsSize++;
	
		DCMcells[dCell]->createDCMvertices(dCellVSize);
		int nDIndex = DCMcells[dCell]->verticesSize-1;
	
		int nVIndex = DCMcells[cellIndex]->verticesSize;
		DCMcells[cellIndex]->verticesSize++;
		
		copy_Daughtercell(cellIndex, dCell);

		DCMcells[cellIndex]->setTargetArea(PM.Init_TargetArea);
		DCMcells[dCell]->setTargetArea(PM.Init_TargetArea);

		for (int i=0; i<DCMcells[dCell]->verticesSize; i++) {

			double nx, ny, nz;
	
			int p = i-1;
			int n = i+1;
		
			if (i == 0) {
				nx = divpx[0];
				ny = divpy[0];
				nz = divpz[0];

				p = nDIndex;

				error = DCMcells[dCell]->setNewVertex(i, nx, ny, nz, p, n, cellIndex, nVIndex, c_r, nRIndex);
				if (error != 0) return error;

				DCMcells[dCell]->DCMvertices[i]->edge_randomness_of_line_coef = DCMcells[cellIndex]->DCMvertices[divp[0]]->edge_randomness_of_line_coef;
			}
			else if (i == nDIndex) {
				nx = divpx[1];
				ny = divpy[1];
				nz = divpz[1];
	
				n = 0;

				error = DCMcells[dCell]->setNewVertex(i, nx, ny, nz, p, n, cellIndex, divp[1], c_l, nLIndex);
				if (error != 0) return error;

				setNewCellvertex_EdgeLineCoef(dCell, i);
			}
			else {
				nx = DCMcells[cellIndex]->DCMvertices[next[i]]->px;
				ny = DCMcells[cellIndex]->DCMvertices[next[i]]->py;
				nz = DCMcells[cellIndex]->DCMvertices[next[i]]->pz;
	
				error = DCMcells[dCell]->setNewVertex(i, nx, ny, nz, p, n);
				if (error != 0) return error;

				DCMcells[dCell]->DCMvertices[i]->edge_randomness_of_line_coef = DCMcells[cellIndex]->DCMvertices[next[i]]->edge_randomness_of_line_coef;
	
				for (int j=0; j<4; j++) {
					DCMcells[dCell]->DCMvertices[i]->othercell[j] = DCMcells[cellIndex]->DCMvertices[next[i]]->othercell[j];
				}

				for (int j=0; j<2; j++) {

					int oc = DCMcells[dCell]->DCMvertices[i]->othercell[2*j];
					int ov = DCMcells[dCell]->DCMvertices[i]->othercell[2*j+1];
					for (int k=0; k<2; k++) {
						if (oc > -1) {
							if (DCMcells[oc]->DCMvertices[ov]->othercell[2*k] == cellIndex) {
								DCMcells[oc]->DCMvertices[ov]->othercell[2*k] = dCell;
								DCMcells[oc]->DCMvertices[ov]->othercell[2*k+1] = i;
							}
						}
					}
				}

				DCMcells[cellIndex]->ignoreVertex(next[i]);
			}
		
		}

		DCMcells[cellIndex]->setCellCycle(1);
		
		error = DCMcells[cellIndex]->setNewVertex(nVIndex, divpx[0], divpy[0], divpz[0], divp[0], divp[1], dCell, 0, c_r, nRIndex);
		if (error != 0) return error;
		setNewCellvertex_EdgeLineCoef(cellIndex, nVIndex);
	
		DCMcells[cellIndex]->DCMvertices[divp[0]]->connection[1] = nVIndex;
		
		DCMcells[cellIndex]->DCMvertices[divp[1]]->ignoreFlag = 0;
		DCMcells[cellIndex]->verticesSize_r++;
	
		DCMcells[cellIndex]->DCMvertices[divp[1]]->px = divpx[1];
		DCMcells[cellIndex]->DCMvertices[divp[1]]->py = divpy[1];
		DCMcells[cellIndex]->DCMvertices[divp[1]]->pz = divpz[1];
	
		DCMcells[cellIndex]->DCMvertices[divp[1]]->connection[0] = nVIndex;

		DCMcells[cellIndex]->setVertexOther(divp[1], dCell, nDIndex, c_l, nLIndex);		
			
		if (c_r > -1) {

			error = DCMcells[c_r]->setNewVertex(nRIndex, divpx[0], divpy[0], divpz[0], v2_r, v1_r, dCell, 0, cellIndex, nVIndex);
			if (error != 0) return error;
			DCMcells[c_r]->DCMvertices[nRIndex]->edge_randomness_of_line_coef = DCMcells[c_r]->DCMvertices[v2_r]->edge_randomness_of_line_coef;
			DCMcells[c_r]->DCMvertices[v1_r]->connection[0] = nRIndex;
			DCMcells[c_r]->DCMvertices[v2_r]->connection[1] = nRIndex;	
		}

		if (c_l > -1) {

			error = DCMcells[c_l]->setNewVertex(nLIndex, divpx[1], divpy[1], divpz[1], v2_l, v1_l, dCell, nDIndex, cellIndex, divp[1]);
			if (error != 0) return error;
			DCMcells[c_l]->DCMvertices[nLIndex]->edge_randomness_of_line_coef = DCMcells[c_l]->DCMvertices[v2_l]->edge_randomness_of_line_coef;
			DCMcells[c_l]->DCMvertices[v1_l]->connection[0] = nLIndex;
			DCMcells[c_l]->DCMvertices[v2_l]->connection[1] = nLIndex;
		}


		if ((recordable_flag & RE::TOPOLOGY_CHANGE) == RE::TOPOLOGY_CHANGE) {

			DCMcells[cellIndex]->setCellArea(1);
			DCMcells[cellIndex]->updateRecord(cellIndex, conTimeStep/100, US::DIVISION_OWN+dCell+1);

			DCMcells[dCell]->setCellArea(1);
			DCMcells[dCell]->updateRecord(dCell, conTimeStep/100, US::DIVISION_OWN+cellIndex+1);

			if (c_r > -1) {
				DCMcells[c_r]->setCellArea(1);
				DCMcells[c_r]->updateRecord(c_r, conTimeStep/100, US::DIVISION_NEI);
			}

			if (c_l > -1) {
				DCMcells[c_l]->setCellArea(1);
				DCMcells[c_l]->updateRecord(c_l, conTimeStep/100, US::DIVISION_NEI);
			}
		}
			
		updateTopology();
		
		copy_Daughtercell2(cellIndex, dCell);
		
		if ((GrowthDivisionFlag & BC::BOUNDARY_ELI) == BC::BOUNDARY_ELI) {
			if (cellsSize_r > 10000) {
				checkOutCells();
			}
		}		
	
		return 0;

	}


	virtual void copy_Daughtercell(int mother, int daughter) {

		DCMcells[daughter]->setTargetArea(DCMcells[mother]->targetArea);
	
		DCMcells[daughter]->statusFlag = DCMcells[mother]->statusFlag;

		DCMcells[daughter]->initialCellIndex = DCMcells[mother]->initialCellIndex;
		DCMcells[daughter]->setCellCycle(1);
		
		DCMcells[daughter]->rec_Size = DCMcells[mother]->rec_Size;
		DCMcells[daughter]->max_rec_Size = DCMcells[mother]->max_rec_Size;
		if (DCMcells[mother]->rec != NULL) {
			DCMcells[daughter]->rec = new CellRecord[DCMcells[daughter]->max_rec_Size];
			memcpy(DCMcells[daughter]->rec, DCMcells[mother]->rec, sizeof(CellRecord)*DCMcells[daughter]->max_rec_Size);
		}
		if (DCMcells[mother]->rec_buf != NULL) {
			DCMcells[daughter]->rec_buf = new CellRecord[DCMcells[daughter]->max_rec_Size];
			memcpy(DCMcells[daughter]->rec_buf, DCMcells[mother]->rec_buf, sizeof(CellRecord)*DCMcells[daughter]->max_rec_Size);
		}
		if (DCMcells[mother]->rec_rect != NULL){
			DCMcells[daughter]->rec_rect = new CellRecord[DCMcells[daughter]->max_rec_Size];
			memcpy(DCMcells[daughter]->rec_rect, DCMcells[mother]->rec_rect, sizeof(CellRecord)*DCMcells[daughter]->max_rec_Size);
		}

		DCMcells[daughter]->lineage_Size = DCMcells[mother]->lineage_Size;
		if (DCMcells[mother]->lineage != NULL) {
			DCMcells[daughter]->lineage = new CellRecord[DCMcells[daughter]->lineage_Size];
			memcpy(DCMcells[daughter]->lineage, DCMcells[mother]->lineage, sizeof(CellRecord)*DCMcells[daughter]->lineage_Size);
		}
		
		for (int i=0; i<rec_StressDivision_index; i++) {
			if (rec_StressDivision[i].inout > 0) {
				if (rec_StressDivision[i].index == mother) {
					rec_StressDivision[i].daughterCell = daughter;
					break;
				}
			}
		}

	}

	virtual void copy_Daughtercell2(int mother, int daughter) {}



	int record_Reconnection(int c_t, int c_b, int c_r, int c_l, int flag, double x=0, double y=0) {

		if (rec_recon_info == NULL) {
			rec_recon_info_size = 8000;			
			rec_recon_info = new ReconnectionRecord[rec_recon_info_size];
		}

		if (recon_number_r+100 > rec_recon_info_size) {
			
			ReconnectionRecord* temp_recon_info = new ReconnectionRecord[rec_recon_info_size];
			memcpy(temp_recon_info, rec_recon_info, rec_recon_info_size*sizeof(ReconnectionRecord));
			delete[] rec_recon_info;		

			rec_recon_info = new ReconnectionRecord[2*rec_recon_info_size];
			memcpy(rec_recon_info, temp_recon_info, rec_recon_info_size*sizeof(ReconnectionRecord));
			delete[] temp_recon_info;

			rec_recon_info_size *= 2;

		}

		if (flag == 0) {

			rec_recon_info[recon_number_r].time_step = conTimeStep;
			rec_recon_info[recon_number_r].position = distanceFromCenter(x, y);

			if (c_t > -1) rec_recon_info[recon_number_r].cell_vertices_before[0] = DCMcells[c_t]->verticesSize_r;
			else rec_recon_info[recon_number_r].cell_vertices_before[0] = -1;

			if (c_b > -1) rec_recon_info[recon_number_r].cell_vertices_before[1] = DCMcells[c_b]->verticesSize_r;
			else rec_recon_info[recon_number_r].cell_vertices_before[1] = -1;

			if (c_r > -1) rec_recon_info[recon_number_r].cell_vertices_before[2] = DCMcells[c_r]->verticesSize_r;
			else rec_recon_info[recon_number_r].cell_vertices_before[2] = -1;

			if (c_l > -1) rec_recon_info[recon_number_r].cell_vertices_before[3] = DCMcells[c_l]->verticesSize_r;
			else rec_recon_info[recon_number_r].cell_vertices_before[3] = -1;
		}
		else {

			if (c_t > -1) rec_recon_info[recon_number_r].cell_vertices_after[0] = DCMcells[c_t]->verticesSize_r;
			else rec_recon_info[recon_number_r].cell_vertices_after[0] = -1;

			if (c_b > -1) rec_recon_info[recon_number_r].cell_vertices_after[1] = DCMcells[c_b]->verticesSize_r;
			else rec_recon_info[recon_number_r].cell_vertices_after[1] = -1;

			if (c_r > -1) rec_recon_info[recon_number_r].cell_vertices_after[2] = DCMcells[c_r]->verticesSize_r;
			else rec_recon_info[recon_number_r].cell_vertices_after[2] = -1;

			if (c_l > -1) rec_recon_info[recon_number_r].cell_vertices_after[3] = DCMcells[c_l]->verticesSize_r;
			else rec_recon_info[recon_number_r].cell_vertices_after[3] = -1;
		}

		return 0;

	}

	
	int record_DivisionStress(int cIndex = -1) {
		
		if (rec_StressDivision == NULL) {
			rec_StressDivision_size = 1000;				
			rec_StressDivision = new stressRecord[rec_StressDivision_size];
		}
		
		if (rec_StressDivision_size < rec_StressDivision_index+10) {

			stressRecord* temp_rec_StressDivision = new stressRecord[rec_StressDivision_size];			
			memcpy(temp_rec_StressDivision, rec_StressDivision, rec_StressDivision_size*sizeof(stressRecord));
			delete[] rec_StressDivision;
			
			rec_StressDivision = new stressRecord[rec_StressDivision_size+1000];			
			memcpy(rec_StressDivision, temp_rec_StressDivision, rec_StressDivision_size*sizeof(stressRecord));				
			delete[] temp_rec_StressDivision;

			rec_StressDivision_size += 1000;
		}
		
		
		if (cIndex > -1) {

			if (DCMcells[cIndex]->verticesSize_r < 30) {
				
				int i = rec_StressDivision_index;
				
				rec_StressDivision[i].index = cIndex;
				rec_StressDivision[i].timeStep = conTimeStep;
				rec_StressDivision[i].verticesSize = DCMcells[cIndex]->verticesSize_r;
				rec_StressDivision[i].inout = 1;

				DCMcells[cIndex]->setCellStressMatrix();
				rec_StressDivision[i].own_st_mag[0] = (DCMcells[cIndex]->cellStressMatrix[4] + DCMcells[cIndex]->cellStressMatrix[5]);
				rec_StressDivision[i].own_st_ani[0] = (DCMcells[cIndex]->cellStressMatrix[4] - DCMcells[cIndex]->cellStressMatrix[5]);

				int vc = 0;
				
				for (int j=0; j<DCMcells[cIndex]->verticesSize; j++) {

					if (DCMcells[cIndex]->DCMvertices[j]->ignoreFlag == 0) {
						
						int oc = DCMcells[cIndex]->DCMvertices[j]->ocIndex;
						
						if (oc > -1) {
							DCMcells[oc]->setCellStressMatrix();
							rec_StressDivision[i].neighborCell[vc] = oc;
							rec_StressDivision[i].mean_st_mag[0] += (DCMcells[oc]->cellStressMatrix[4] + DCMcells[oc]->cellStressMatrix[5]);
							rec_StressDivision[i].mean_st_ani[0] += (DCMcells[oc]->cellStressMatrix[4] - DCMcells[oc]->cellStressMatrix[5]);
							vc++;
						}
					}
				}
				
				rec_StressDivision[i].mean_st_mag[0] /= static_cast<double>(vc);
				rec_StressDivision[i].mean_st_ani[0] /= static_cast<double>(vc);

				rec_StressDivision_index++;
				
			}

		}

		else {
			
			for (int i=0; i<rec_StressDivision_index; i++) {

				if (rec_StressDivision[i].inout > 0 && rec_StressDivision[i].daughterCell > -1) {

					int c = rec_StressDivision[i].index;
					
					if ( ((conTimeStep - rec_StressDivision[i].timeStep)/static_cast<int>(1.0/PM.Delta_t)) > 10
						|| DCMcells[c]->ignoreFlag != 0) {
						rec_StressDivision[i].inout = -1;
					}
					
					else if (abs(DCMcells[c]->cellArea[0]-DCMcells[c]->cellArea[1])/PM.Delta_t < 1.0e-4
						&& DCMcells[c]->cellArea[0] != DCMcells[c]->cellArea[1]) {

						DCMcells[c]->setCellStressMatrix();
						rec_StressDivision[i].own_st_mag[1] = (DCMcells[c]->cellStressMatrix[4] + DCMcells[c]->cellStressMatrix[5]);
						rec_StressDivision[i].own_st_ani[1] = (DCMcells[c]->cellStressMatrix[4] - DCMcells[c]->cellStressMatrix[5]);

						int d = rec_StressDivision[i].daughterCell;						
						DCMcells[d]->setCellStressMatrix();
						rec_StressDivision[i].own_st_mag[2] = (DCMcells[d]->cellStressMatrix[4] + DCMcells[d]->cellStressMatrix[5]);
						rec_StressDivision[i].own_st_ani[2] = (DCMcells[d]->cellStressMatrix[4] - DCMcells[d]->cellStressMatrix[5]);

						for (int j=0; j<rec_StressDivision[i].verticesSize; j++) {

							int oc = rec_StressDivision[i].neighborCell[j];

							if (oc < 0) {
								rec_StressDivision[i].inout = -1;
								break;
							}

							if (DCMcells[oc]->ignoreFlag != 0) {
								rec_StressDivision[i].inout = -1;
								break;
							}

							DCMcells[oc]->setCellStressMatrix();

							rec_StressDivision[i].mean_st_mag[1] += (DCMcells[oc]->cellStressMatrix[4] + DCMcells[oc]->cellStressMatrix[5]);
							rec_StressDivision[i].mean_st_ani[1] += (DCMcells[oc]->cellStressMatrix[4] - DCMcells[oc]->cellStressMatrix[5]);

						}

						if (rec_StressDivision[i].inout > 0) {
							rec_StressDivision[i].mean_st_mag[1] /= static_cast<double>(rec_StressDivision[i].verticesSize);
							rec_StressDivision[i].mean_st_ani[1] /= static_cast<double>(rec_StressDivision[i].verticesSize);
							rec_StressDivision[i].inout = 0;
						}					
					}
				}
			}
		}		

		return 0;

	}


	int record_DivisionApoptosis(int cIndex = -1) {
		
		if (rec_StressApoptosis == NULL) {
			rec_StressApoptosis_size = 1000;				
			rec_StressApoptosis = new stressRecord[rec_StressApoptosis_size];
		}
		
		if (rec_StressApoptosis_size < rec_StressApoptosis_index+10) {

			stressRecord* temp_rec_StressApoptosis = new stressRecord[rec_StressApoptosis_size];			
			memcpy(temp_rec_StressApoptosis, rec_StressApoptosis, rec_StressApoptosis_size*sizeof(stressRecord));
			delete[] rec_StressApoptosis;
			
			rec_StressApoptosis = new stressRecord[rec_StressApoptosis_size+1000];			
			memcpy(rec_StressApoptosis, temp_rec_StressApoptosis, rec_StressApoptosis_size*sizeof(stressRecord));				
			delete[] temp_rec_StressApoptosis;

			rec_StressApoptosis_size += 1000;
		}
		
		if (cIndex > -1) {

			if (DCMcells[cIndex]->verticesSize_r < 30) {
				
				int i = rec_StressApoptosis_index;
				
				rec_StressApoptosis[i].index = cIndex;
				rec_StressApoptosis[i].timeStep = conTimeStep;
				rec_StressApoptosis[i].verticesSize = DCMcells[cIndex]->verticesSize_r;
				rec_StressApoptosis[i].inout = 1;

				rec_StressApoptosis[i].area = distanceFromCenter(DCMcells[cIndex]->cofm_x, DCMcells[cIndex]->cofm_y);

				int vc = 0;
				
				for (int j=0; j<DCMcells[cIndex]->verticesSize; j++) {

					if (DCMcells[cIndex]->DCMvertices[j]->ignoreFlag == 0) {
						
						int oc = DCMcells[cIndex]->DCMvertices[j]->ocIndex;
						
						if (oc > -1) {
							DCMcells[oc]->setCellStressMatrix();
							rec_StressApoptosis[i].neighborCell[vc] = oc;
							rec_StressApoptosis[i].mean_st_mag[0] += (DCMcells[oc]->cellStressMatrix[4] + DCMcells[oc]->cellStressMatrix[5]);
							rec_StressApoptosis[i].mean_st_ani[0] += (DCMcells[oc]->cellStressMatrix[4] - DCMcells[oc]->cellStressMatrix[5]);
							vc++;
						}
					}
				}
				
				rec_StressApoptosis[i].mean_st_mag[0] /= static_cast<double>(vc);
				rec_StressApoptosis[i].mean_st_ani[0] /= static_cast<double>(vc);

				NeighborCell nnc[100];

				search_NeighborCells(nnc, 100, cIndex, 4);

				for (int j=0; j<100; j++) {
					rec_StressApoptosis[i].nnc[j].index = nnc[j].index;
					rec_StressApoptosis[i].nnc[j].nn_level = nnc[j].nn_level;
					rec_StressApoptosis[i].nnc[j].cell_area_be = DCMcells[nnc[j].index]->cellArea[0];

				}

				rec_StressApoptosis_index++;
				
			}

		}

		else {
			
			for (int i=0; i<rec_StressApoptosis_index; i++) {

				if (rec_StressApoptosis[i].inout > 0) {

					int r = static_cast<int>(SF.uniformDistribution(0,3))%3;
					int c = rec_StressApoptosis[i].neighborCell[r];
					
					if ( ((conTimeStep - rec_StressApoptosis[i].timeStep)/static_cast<int>(1.0/PM.Delta_t)) > 10
						|| DCMcells[c]->ignoreFlag != 0 ) {
						rec_StressApoptosis[i].inout = -1;
					}
					
					else if (abs(DCMcells[c]->cellArea[0]-DCMcells[c]->cellArea[1])/PM.Delta_t < 1.0e-4
						&& DCMcells[c]->cellArea[0] != DCMcells[c]->cellArea[1]) {

						for (int j=0; j<rec_StressApoptosis[i].verticesSize; j++) {

							int oc = rec_StressApoptosis[i].neighborCell[j];

							if (oc < 0) {
								rec_StressApoptosis[i].inout = -1;
								break;
							}

							if (DCMcells[oc]->ignoreFlag != 0) {
								rec_StressApoptosis[i].inout = -1;
								break;
							}

							DCMcells[oc]->setCellStressMatrix();

							rec_StressApoptosis[i].mean_st_mag[1] += (DCMcells[oc]->cellStressMatrix[4] + DCMcells[oc]->cellStressMatrix[5]);
							rec_StressApoptosis[i].mean_st_ani[1] += (DCMcells[oc]->cellStressMatrix[4] - DCMcells[oc]->cellStressMatrix[5]);

						}

						if (rec_StressApoptosis[i].inout > 0) {
							rec_StressApoptosis[i].mean_st_mag[1] /= static_cast<double>(rec_StressApoptosis[i].verticesSize);
							rec_StressApoptosis[i].mean_st_ani[1] /= static_cast<double>(rec_StressApoptosis[i].verticesSize);
							rec_StressApoptosis[i].inout = 0;
						}

						for (int j=0; j<100; j++) {
							int oc = rec_StressApoptosis[i].nnc[j].index;

							if (oc < 0) break;
							
							if (j == 0) rec_StressApoptosis[i].nnc[j].cell_area_af = 0;
							else {
								if (DCMcells[oc]->ignoreFlag != 0 || DCMcells[oc]->divArea >= 0) rec_StressApoptosis[i].nnc[j].index = -1;
								else rec_StressApoptosis[i].nnc[j].cell_area_af = DCMcells[oc]->cellArea[0];
							}
						}

					}

				}
			}
			
		}

		return 0;

	}


	void checkOutCells() {
		
		std::vector<int> out_cells;

		int sc, nc;

		for (int i=0; i<cellsSize; i++) {
			if (DCMcells[i]->ignoreFlag == 0 && DCMcells[i]->outsideFlag != 0) {
				sc = nc = i;
				out_cells.push_back(sc);
				break;
			}
		}

		do {

			for (int j=0; j<DCMcells[nc]->verticesSize; j++) {
				if (DCMcells[nc]->DCMvertices[j]->ignoreFlag == 0 && DCMcells[nc]->DCMvertices[j]->ocIndex < 0) {
					int n = DCMcells[nc]->DCMvertices[j]->connection[1];
					int oc = DCMcells[nc]->DCMvertices[n]->ocIndex;
					if (oc != -1) {
						if (DCMcells[oc]->outsideFlag != 0) {
							nc = oc;
							out_cells.push_back(nc);
							break;
						}
					}
				}
			}

		} while (sc != nc);

		out_cells.pop_back();

		bool check = true;
		
		do {
			int selected = static_cast<int>(SF.uniformDistribution(0.0, static_cast<double>(out_cells.size())));
			
			if (DCMcells[out_cells[selected]]->out_cell_eli_flag == 0) {
				DCMcells[out_cells[selected]]->out_cell_eli_flag = 1;
				check = false;
			}
		} while (check);

	}


	void setCellGrad(int flag=0) {

		for (int i=0; i<cellsSize; i++) {

			if (DCMcells[i]->ignoreFlag == 0) {
				DCMcells[i]->setCellArea(1);
				DCMcells[i]->setPerimeterLength();
				DCMcells[i]->setCellShape();
				DCMcells[i]->setCellStressMatrix(flag);
			}
		}
		
		if (flag > 0) {

			for (int i=0; i<cellsSize; i++) {

				double x2_i = 0;
				double xy_i = 0;
				double y2_i = 0;
				double smx_i = 0;
				double smy_i = 0;
				double sax_i = 0;
				double say_i = 0;

				int oc;

				double cx, cy, st_mag, st_ani;
			
				if (DCMcells[i]->ignoreFlag == 0 && DCMcells[i]->outsideFlag == 0) {

					for (int j=0; j<DCMcells[i]->verticesSize; j++) {

						if (DCMcells[i]->DCMvertices[j]->ignoreFlag == 0) {
						
							oc = DCMcells[i]->DCMvertices[j]->ocIndex;

							cx = DCMcells[i]->cofm_x - DCMcells[oc]->cofm_x;
							cy = DCMcells[i]->cofm_y - DCMcells[oc]->cofm_y;
							st_mag = DCMcells[oc]->cellStressMatrix[4] + DCMcells[oc]->cellStressMatrix[5];
							st_ani = DCMcells[oc]->cellStressMatrix[4] - DCMcells[oc]->cellStressMatrix[5];

							x2_i += cx*cx;
							xy_i += cx*cy;
							y2_i += cy*cy;
							smx_i += st_mag*cx;
							smy_i += st_mag*cy;
							sax_i += st_ani*cx;
							say_i += st_ani*cy;

						}

					}

					DCMcells[i]->setCellGrad((smx_i*y2_i-xy_i*smy_i)/(x2_i*y2_i-xy_i*xy_i), (smy_i*x2_i-xy_i*smx_i)/(x2_i*y2_i-xy_i*xy_i),
						(sax_i*y2_i-xy_i*say_i)/(x2_i*y2_i-xy_i*xy_i), (say_i*x2_i-xy_i*sax_i)/(x2_i*y2_i-xy_i*xy_i) );
			
				}
		
			}
	
			for (int i=0; i<cellsSize; i++) {

				int oc;
				int check = 0;

				double mean_grad_mag = 0.0;
				double mean_grad_ani = 0.0;

				if (DCMcells[i]->ignoreFlag == 0 && DCMcells[i]->outsideFlag == 0) {

					for (int j=0; j<DCMcells[i]->verticesSize; j++) {

						if (DCMcells[i]->DCMvertices[j]->ignoreFlag == 0) {

							oc = DCMcells[i]->DCMvertices[j]->ocIndex;

							if (DCMcells[oc]->outsideFlag == 0) {
								mean_grad_mag += DCMcells[oc]->cellGrad[2];
								mean_grad_ani += DCMcells[oc]->cellGrad[5];
							}
							else {
								check = 1;
								break;							
							}

						}

					}

					if (check == 0) {
						DCMcells[i]->cellGrad[6] = mean_grad_mag/static_cast<double>(DCMcells[i]->verticesSize_r);
						DCMcells[i]->cellGrad[7] = mean_grad_ani/static_cast<double>(DCMcells[i]->verticesSize_r);
					}
					else {
						DCMcells[i]->cellGrad[6] = -1;
						DCMcells[i]->cellGrad[7] = -1;
					}

				}
			}
		}
	}

	
};



#endif //DCMODELSTRUCT2D_HH