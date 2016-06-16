#pragma once
#include "eigenIncludes.h"

typedef struct tagKINEMATICPARAMS6A {
	//0 grf, x directed along proximal shaft, z perpendicular to the table, y to the right when looking negative z
	//1 base roll axis at pitch joint
	//2 rotation from the pitch joint
	//3 translate along x2 to the catheter base, apply roll rotation
	//4 translate & rotate to the catheter tip
	//5 translate along x4 the projected distance
	double tx01 = 766.6;
	double ty01 = -112;
	double tz01 = 14.7;
	double rz01 = -0.26;
	double tx23 = 8.2;
	double cathL = 98.75; //[mm] straight length
} KINEMATICPARAMS6A;

typedef struct tagKINEMATICPARAMS11A {
	//frames as above, adding static rotations Ry01, Ry34, Rz34, Ry45, Rz45 to accomodate out-of-plane movementss
	double tx01 = 766.6;
	double ty01 = -112;
	double tz01 = 14.7;
	double rz01 = -0.26;
	double ry01 = 0;
	double tx23 = 8.2;
	double ry34 = 0;
	double rz34 = 0;
	double cathL = 98.75; //[mm] straight length
	double ry45 = 0;
	double rz45 = 0;
} KINEMATICPARAMS11A;

//struct kinematicParams{
//};

typedef struct tagJOINTLIMITS {
	double up[5] = { 2,.8,1,5,400 };
	double dn[5] = { -3,-.8,-1,.0001,0 };
} JOINTLIMITS;

enum nlMethod {
	GN_DIRECT = 0,
	GN_DIRECT_L = 1,
	GN_DIRECT_L_NOSCAL = 2,
	GN_DIRECT_L_RAND = 3,
	GN_DIRECT_L_RAND_NOSCAL = 4,
	GN_DIRECT_NOSCAL = 5,
	GN_ESCH = 6,
	GN_ISRES = 7,
	GN_MLSL = 8,
	GN_MLSL_LDS = 9,
	GN_ORIG_DIRECT = 10,
	GN_ORIG_DIRECT_L = 11,
	LN_BOBYQA = 12,
	LN_COBYLA = 13,
	LN_NELDERMEAD = 14,
	LN_NEWUOA = 15,
	LN_NEWUOA_BOUND = 16,
	LN_PRAXIS = 17,
	LN_SBPLX = 18
};

typedef struct tagNLOPTPARAMS {
	double minFunVal = -1e3; // stop when fmin is less than this
	double tolFunAbs = 1e-4; // minimum function change between steps
	double tolXAbs = 1e-4;   // minimum x change between steps
	int maxIts = 100;     // maximum number of function evals
	double maxTimeSec = 2;
	nlMethod method = nlMethod::LN_BOBYQA;
} NLOPTPARAMS;