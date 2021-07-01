#ifndef DCMODELSTRUCT_HH
#define DCMODELSTRUCT_HH


////// Including //////////////////////////////////////////

#include <cmath>
#include <time.h>
#include <stdio.h>
#include <wchar.h>
#include <stdlib.h>
#include <float.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <assert.h>
#include <algorithm>


////// Macro Definition ///////////////////////////////////

#define ArrayCount(x) sizeof(x) / sizeof(x[0])

#define chCopy(other, dynamicArray, arrayType, arraySize) \
	if (other.dynamicArray == NULL) { dynamicArray = NULL; } \
		else { \
		dynamicArray = new arrayType[arraySize]; \
		memcpy(dynamicArray, other.dynamicArray, arraySize * sizeof(arrayType) ); }

#define chPCopy(other, dynamicArray, arrayType, arraySize) \
	if (other.dynamicArray == NULL) { dynamicArray = NULL; } \
		else {\
		dynamicArray = new arrayType*[arraySize]; \
		for (int iter=0; iter<arraySize; iter++) { \
			dynamicArray[iter] = new arrayType; \
			memcpy(dynamicArray[iter], other.dynamicArray[iter], sizeof(arrayType) ); } }

#define chDelete(dynamicArray) \
	if (dynamicArray != NULL) { \
		delete[] dynamicArray; \
		dynamicArray = NULL; }

#define chPDelete(dynamicArray, arraySize) \
	if (dynamicArray != NULL) { \
		for (int iter=0; iter<arraySize; iter++) { \
			if (dynamicArray[iter] != NULL) { \
				delete dynamicArray[iter]; \
				dynamicArray[iter] = NULL; } } \
		delete[] dynamicArray; \
		dynamicArray = NULL; }

#define initArray(arr, value) \
	for (int i=0; i<(sizeof(arr)/sizeof(arr[0])); i++) \
		arr[i] = value;

#define copyArray(other, arr) \
	for (int i=0; i<(sizeof(arr)/sizeof(arr[0])); i++) \
		arr[i] = other.arr[i];


////// Using Namespace ////////////////////////////////////

using namespace std;



////// Constant Definition ////////////////////////////////

const int DATA_VERSION = 1006;
const int DATA_VERSION_2D = 1001;

const double Pi = 3.1415926535897932384626433832795;
const double VERY_SMALL = 1.0e-9;
const double VERY_BIG = 1.0e10;
const double D3 = 1.0/3.0;

const int maxCellsSize = 500000;


////// Inline Function ////////////////////////////////////

inline double dist(double x1, double x2, double y1, double y2, double z1, double z2) {
	return sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) );
}

inline double coord(double xyz1, double xyz2, double t) {
	return xyz1 + (xyz2-xyz1)*t;
}

inline double dot(double x1, double y1, double z1, double x2, double y2, double z2) {
	return ( (x1*x2) + (y1*y2) + (z1*z2) );
}

inline double det(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) {
	return ( (x1*y2*z3) + (y1*z2*x3) + (z1*x2*y3) - (z1*y2*x3) - (x1*z2*y3) - (y1*x2*z3) );
}

inline double sqrt3(double x, double y, double z) {
	return sqrt( x*x + y*y + z*z );
}


////// Structure //////////////////////////////////////////

struct Vpoint {

	double x1, y1, z1;
	double x2, y2, z2;
	
	int colorflag;
	int colorGrad;

};

struct VTpoint {

	double x1, y1, z1;
	double x2, y2, z2;
	double x3, y3, z3;

	int colorflag;
	int colorGrad;

};


////// Base Class Declaration /////////////////////////////

class DCMvertex {

public:
	
	DCMvertex() {
		px = 0;	py = 0;	pz = 0;	
		ignoreFlag = -1;
		moveFlag = 0;
	}

	DCMvertex(const DCMvertex& otherVertex) {
		copyDCMvertex(otherVertex);
	}

	DCMvertex& operator=(const DCMvertex& otherVertex) {
		if (this != &otherVertex) {
			copyDCMvertex(otherVertex);
		}
		return *this;
	}

	void copyDCMvertex(const DCMvertex& otherVertex) {
		px = otherVertex.px; py = otherVertex.py; pz = otherVertex.pz;
		ignoreFlag = otherVertex.ignoreFlag; moveFlag = otherVertex.moveFlag;
	}

	virtual ~DCMvertex() = 0 {}

	//////////////////////////////
	
	double px, py, pz;

	int ignoreFlag;
	int moveFlag;

};


class DCMpoint {

public:

	DCMpoint(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {
		vx = 0; vy = 0; vz = 0;
	}

	DCMpoint() {
		x = 0; y = 0; z = 0;
		vx = 0; vy = 0; vz = 0;
	}

	DCMpoint(const DCMpoint& otherPoint) {
		copyDCMpoint(otherPoint);
	}

	DCMpoint& operator=(const DCMpoint& otherPoint) {
		if (this != &otherPoint) {
			copyDCMpoint(otherPoint);
		}
		return *this;
	}

	void copyDCMpoint(const DCMpoint& otherPoint) {
		x = otherPoint.x; y = otherPoint.y; z = otherPoint.z;
		vx = otherPoint.vx; vy = otherPoint.vy; vz = otherPoint.vz;
	}

	virtual ~DCMpoint() {}

	//////////////////////////////
	
	double x, y, z;
	double vx, vy, vz;
	
};


class DCMedge {

public:

	DCMedge() {
		cellNum = 0;
		distance = 0;
	}

	DCMedge(const DCMedge& otherEdge) {
		copyDCMedge(otherEdge);
	}

	DCMedge& operator=(const DCMedge& otherEdge) {
		if (this != &otherEdge) {
			copyDCMedge(otherEdge);
		}
		return *this;
	}

	void copyDCMedge(const DCMedge& otherEdge) {
		cellNum = otherEdge.cellNum; distance = otherEdge.distance;
		inoutFlag = otherEdge.inoutFlag; checkFlag = otherEdge.checkFlag; vLineFlag = otherEdge.vLineFlag;
	}

	virtual ~DCMedge() = 0 {}

	//////////////////////////////
	
	int cellNum;
	double distance;
	int inoutFlag;
	int checkFlag;
	int vLineFlag;

};


class DCMplane {

public:

	DCMplane() {}
	virtual ~DCMplane() = 0 {}

	//////////////////////////////

};


class DCMcell {

public:
		
	DCMcell() {}
	virtual ~DCMcell() = 0 {}

	virtual void createDCMvertices(int vSize) {}
	virtual int setNewVertex(int pi, double px, double py, double pz, int pp, int pn, int oc1 = -1, int ov1 = -1, int oc2 = -1, int ov2 = -1) { return 0; }

	//////////////////////////////
		
	int verticesSize;
	int verticesSize_r;
	int maxVerticesSize;
	int ignoreFlag;	
	int outsideFlag;
	double cofv_x, cofv_y, cofv_z;	

};



class DCMcontainer {

public:
	
	DCMcontainer() {}
	virtual ~DCMcontainer() = 0 {}

	virtual void initialize(int flag) {}
	virtual void setOtherCell() {}

	virtual void drawEdges(Vpoint* vpoint, unsigned int& size) {}
	virtual void drawEdges(Vpoint* vpoint, int flag) {}
	virtual void drawPlanes(VTpoint* vtpoint, int& size) {}
	virtual void drawPlanes(VTpoint* vtpoint, int& size, int flag) {}

	virtual void createDCMcells() {}
	virtual void createDCMcells(int cSize) {}

	virtual void createDCMvertices(int cIndex, int vSize) {}
	virtual void setNewVertex(int cIndex, int pi, double px, double py, double pz, int pp, int pn, int oc1 = -1, int ov1 = -1, int oc2 = -1, int ov2 = -1) {}

	virtual int checkContainer(int flag=0) { return 0; }
	virtual int finishContainer(int flag=0) { return 0; }


	//////////////////////////////
	
	int cellsSize;		
	int cellsSize_r;

	int edgesSize;

	int verticesSize;

	int planesSize;

	int containerVerticesSize;

	double averageVolume;

	int reconnectionFlag;
	int dFlag;			
	int objectFlag;		
	unsigned long long conTimeStep;

	double cCenterx, cCentery, cCenterz;

	double basex, basey, basez;

	double tresholdCos;

	int outVerticesIndex;

	double outVertexMax;

	int dataVersion;

	double timeStepMaxVectorForce;

	double laps;

	double maxdxdt;

	double maxmaxdxdt;

	int initDataFlag;


};



///////////////////////////////////////////////////////////


template <typename T>
class Singleton {

private:
	static T* pInstance;

public:

	static T& getInstance() {
		if (pInstance == 0) {
			pInstance = new T();
			pInstance->initialize();
		}
		return *pInstance;
	}

	static T* getPointer() {		
		if (pInstance == 0) {
			pInstance = new T();
			pInstance->initialize();
		}
		return pInstance;		
	}

	static void releasePointer() {
		if (pInstance != 0) {
			pInstance->destroy();
			delete pInstance;
			pInstance = 0;
		}
	}
	
	virtual void initialize() = 0 {}
	virtual void destroy() = 0 {}

protected:
	Singleton() {}
	~Singleton() {}
	
private:
	Singleton(const Singleton&);
	Singleton& operator=(const Singleton&);
	
};

template <class T> T* Singleton<T>::pInstance = 0;





#define SF SharingFunctions::getInstance()

class SharingFunctions : public Singleton<SharingFunctions> {

	friend class Singleton<SharingFunctions>;

protected:
	SharingFunctions() {}

private:
	~SharingFunctions() {}

private:
	unsigned long state[16];
	unsigned int index;
	unsigned long seed;
	
public:

	
	virtual void initialize() {
				
		seeding_WellRNG512(static_cast<unsigned long>(time(NULL)));		
				
	}

	virtual void destroy() {}

	void seeding_WellRNG512(unsigned long s) {

		seed = s;
		index = 0;

		for (int i=0; i<16; i++) {
			state[i] = s;
			s += s + 50;
		}

	}

	unsigned long return_seed() {

		return seed;

	}


private:
	
	// WELL512 algorithm implementation by Chris Lomont (public domain)

	unsigned long WellRNG512() {

		unsigned long a, b, c, d;
		a = state[index];
		c = state[(index+13)&15];
		b = a^c^(a<<16)^(c<<15);
		c = state[(index+9)&15];
		c ^= (c>>11);
		a = state[index] = b^c;
		d = a^((a<<5)&0xDA442D24UL);
		index = (index + 15)&15;
		a = state[index];
		state[index] = a^b^d^(a<<2)^(b<<18)^(c<<28);

		return state[index];

	}


public:

	double uniformDistribution(double f, double t) {
		return (static_cast<double>(WellRNG512())/static_cast<double>(0xffffffffUL))*(t-f)+f;
	}

	unsigned int uniformDistribution_int(unsigned int f, unsigned int t) {
		return (WellRNG512()%(t-f)+f);
	}

	double gaussDistribution(double mean, double stdev) {

		double x1, x2, w, y1;
		static double y2;
		static int use_last = 0;

		if (use_last) {
			y1 = y2;
			use_last = 0;
		}
		else {

			do {
				x1 = uniformDistribution(-1.0, 1.0);
				x2 = uniformDistribution(-1.0, 1.0);				
				w = x1*x1 + x2*x2;
			} while ( w >= 1.0 );

			w = sqrt((-2.0*log(w))/w);
			y1 = x1*w;
			y2 = x2*w;
			use_last = 1;
		}

		return (mean + y1*stdev);
	}
	
	double vonmisesDistribution(double mean, double k) {

		double result = 0.0;

		double a = 1.0 + sqrt(1.0 + 4.0*(k*k));
		double b = (a - sqrt(2.0*a)) / (2.0*k);
		double r = (1.0 + b*b) / (2.0*b);

		while(1) {

			double U1 = uniformDistribution(0.0, 1.0);
			double z = cos(Pi * U1);
			double f = (1.0 + r*z) / (r+z);
			double c = k*(r-f);
			double U2 = uniformDistribution(0.0, 1.0);

			if (c*(2.0-c) - U2 > 0.0) {

				double U3 = uniformDistribution(0.0, 1.0);
				double sign = 0.0;

				if (U3 - 0.5 < 0.0) {
					sign = -1.0;
				}
				else if (U3 - 0.5 > 0.0) {
					sign = 1.0;
				}

				result = sign * acos(f) *0.5 + mean;

				while (result >= 1.0 * Pi) {
					result -= 1.0 * Pi;
				}

				break;
			}

			else {

				if (log(c/U2) + 1.0 - c >= 0.0) {

					double U3 = uniformDistribution(0.0, 1.0);
					double sign = 0.0;

					if (U3 - 0.5 < 0.0) {
						sign = -1.0;
					}
					else if (U3 - 0.5 > 0.0) {
						sign = 1.0;
					}
				
					result = sign * acos(f) *0.5 + mean;
				
					while (result >= 1.0 * Pi) {
						result -= 1.0 * Pi;
					}

					break;
				}
			}
		}
	
		return result;

	}


	double gammaDistribution(double alpha, double beta) {
		
		double a1, a2;
		double oalpha = alpha;
		double u, v, x;

		if (alpha <= 0) return -1;
		if (alpha < 1.0) alpha += 1.0;
		a1 = alpha - 1.0/3.0;
		a2 = 1.0/sqrt(9.0*a1);

		do {
			do {
				x = gaussDistribution(0, 1);
				v = 1.0 + a2*x;
			} while (v <= 0.0);

			v = v*v*v;
			u = uniformDistribution(0,1);

		} while (u > 1.0 - 0.331*x*x*x*x && log(u) > 0.5*x*x + a1*(1.0-v+log(v)) );

		if (alpha == oalpha) return a1*v/beta;
		else {
			do {
				u = uniformDistribution(0,1);
			} while (u == 0.0);
			return pow(u, 1.0/oalpha)*a1*v/beta;
		}
		
	}

	double betaDistribution(double alpha, double beta) {

		double x = gammaDistribution(alpha, 1);
		double y = gammaDistribution(beta, 1);

		return x/(x+y);

	}

	double mbesseli0(double x) {
		
		double sum = 0.0;
		
		for (int m = 0; m < 200; ++m) {
			
			int factM;

			if (m == 0) {
				factM = 1;
			}
			else {
				int val = 1;
				
				for( int idx = 1; idx <=m; ++idx) {
					val *= idx;
				}
				factM = val;
			}
			
			double inc = pow(1.0/static_cast<double>(factM) * pow(x * 0.5, static_cast<double>(m)), 2.0);
			
			double frac = inc / sum;
			
			sum += inc;
			
			if( frac < 0.001) {
				break;
			}
		}
		
		return sum;
	}



	DCMpoint vProduct(const DCMpoint& v1, const DCMpoint& v2) {
	
		DCMpoint tempv;
	
		tempv.x = (v1.y*v2.z) - (v1.z*v2.y);
		tempv.y = (v1.z*v2.x) - (v1.x*v2.z);
		tempv.z = (v1.x*v2.y) - (v1.y*v2.x);

		return tempv;

	}


};


///////////////////////////////////////////////////////////
// DCMparameters class

#define PM DCMparameters::getInstance()

class DCMparameters : public Singleton<DCMparameters> {

	friend class Singleton<DCMparameters>;

protected:
	DCMparameters() {}

private:
	~DCMparameters() {}

public:
	
	/// Voronoi Initializing
	double Othercell_Distance;
	int initial_cell_size;

	/// Delta T
	double Delta_t;
	double Eta;
	
	/// Reconnection
	double Epsilon_Distance;	
	double Const_ItoH;			
	double Const_HtoI;			

	double Extrusion_Area;
	double Extrusion_Stmag;

	/// Dynamics

	double Kappa;				
	double Init_TargetArea;	
	
	double Sigma_Zero;			
	double Sigma;				
	double Sigma_Rate;				
	double Randomness_of_Lambda;
	double Decay_of_Lambda;

	double Line_Coef;
	double Line_Coef_Out;			 

	double Line_Beta;
	double Line_Mu;

	double SolidAngle_const;
	double Minimum_Omega;	

	double FlatK;			

	double Stretch_Force;			
	double Tug_Force;
	double RandomMove;

	/// Division (cell cycle)

	double CellCycleTime;
	double CellCycleSpeed;				
	double CellCycleSync;				

	/// Pressure Dynamics
	
	double Pressure_const;
	double OutsidePressureUp;
	double OutsidePressureDown;
	double Pressure_Curvature;

	/// Non-linear Elasticity Line Tension

	double alpha_const;
	double beta_const;
	double k_const;
	double h_const;
	
	double vVectorConst;

	FILE *dataFile;
	wchar_t *dataFileName;

	///
	unsigned int view_flag;

public:

	virtual void initialize() {

		Othercell_Distance = 1.0e-10;
		initial_cell_size = 1000;

		Delta_t = 0.001;
		Eta = 100.0;
		Epsilon_Distance = 0.1;
		Const_ItoH = 1.2;
		Const_HtoI = 1.2;
		Extrusion_Area = 0.2;
		Extrusion_Stmag = -0.4;
	
		Sigma_Zero = 3.0;
		Sigma = 0.12;
		Sigma_Rate = 100.0;
		Kappa = 1.0;
		Init_TargetArea = 1.0;
		Line_Coef = 0.04;				
		Line_Coef_Out = 0.04;
		Line_Beta = 0.1;
		Line_Mu = 0.76;
		Randomness_of_Lambda = 0;
		Decay_of_Lambda = 0;

		SolidAngle_const = 1.0;
		Minimum_Omega = 1.0e-5;
		FlatK = 1.0e-2;
		Stretch_Force = -0.35;		
		Tug_Force = 1.0e-1;
		RandomMove = 0;

		CellCycleTime = 1.0;
		CellCycleSpeed = 0.005;
		CellCycleSync = 0.4;
	
		Pressure_const = 1.0e-2;
		OutsidePressureUp = 10.0;
		OutsidePressureDown = 9.0;
		Pressure_Curvature = 1.0;

		alpha_const = 0.2;
		beta_const = 1.0;
		k_const = 0.1;
		h_const = 1.0;

		vVectorConst = 2000.0;

		dataFileName = new wchar_t[50];
	
		view_flag = 0;

	}

	virtual void destroy() {
		chDelete(dataFileName);
	}
	
};






#endif //DCMODELSTRUCT_HH