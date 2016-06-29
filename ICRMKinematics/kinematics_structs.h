#pragma once
#include "eigenIncludes.h"

// see kinematicFrames.png
typedef struct tagPHYSICALDIMENSIONS{
	double lProx = 25.50; //[mm] distance from proximal's proximal face to pitch axis
	double lPtch = 11.03;  //[mm] distance from pitch axis to roll's proximal face
	double lRoll = 12.49;  //[mm] distance from roll's proximal face to roll distal edge
	double rProx = 6.00;   //[mm] radius of proximal, pitch, and roll
	double lCath = 95;     //[mm]
	double rCath = 3.175;  //[mm]
} PHYSICALDIMENSIONS;

// a minimal set of kinematics: locate the base roll and catheter length
typedef struct tagKINEMATICPARAMS5A {
	//0 grf, x directed along proximal shaft, z perpendicular to the table, y to the right when looking negative z
	//1 the distal face of the coaxial input
	//2 the pitch axis, rotating about z2
	//3 the proximal face of the roll, rotating about x2=x3
	//4 the catheter tip, remote center about 
	//5 translate along x4 the projected distance
	double tx01 = 806;
	double ty01 = -66;
	double tz01 = -28;
	double rz01 = -.24;
	double lCath = 95; //[mm] straight length
} KINEMATICPARAMS5A;

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
	double ry34 = 0.0;
	double lCath = 95; //[mm] straight length
} KINEMATICPARAMS6A;

typedef struct tagKINEMATICPARAMS11A {
	//frames as above, adding static rotations Ry01, Ry34, Rz34, Ry45, Rz45 to accomodate out-of-plane movementss
	double tx01 = 766.6;
	double ty01 = -112;
	double tz01 = 14.7;
	double ry01 = 0;
	double rz01 = -0.26;
	double ry34 = 0; //rotation of the catheter actuation plane, mostly in case of drooping
	double rz34 = 0; //redundant with alpha near 0
	double kAlpha = 1; //alpha gain
	double eAlpha = 1; //alpha exponent
	double lCath = 95; //[mm] straight length
	double ry45 = 0;
} KINEMATICPARAMS11A;

typedef struct tagJOINTLIMITS {
	double up[5] = {  3, .8, 1,    5,500 };
	double dn[5] = { -3,-.8,-1,.0001,  0 };
} JOINTLIMITS;

enum nlMethod {//redo
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
	LN_SBPLX = 18,
	GD_MLSL = 19,
	GD_MLSL_LDS = 20,
	GD_STOGO = 21,
	GD_STOGO_RAND = 22,
	LD_CCSAQ = 23,
	LD_LBFGS = 24, 
	LD_LBFGS_NOCEDAL = 25,
	LD_MMA = 26,
	LD_TNEWTON = 27, 
	LD_TNEWTON_RESTART = 28,
	LD_TNEWTON_PRECOND = 29,
	LD_TNEWTON_PRECOND_RESTART = 30,
	LD_VAR1 = 31,
	LD_VAR2 = 32,
	LD_SLSQP = 33
};

typedef struct tagNLOPTPARAMS {
	double minFunVal = -1e3; // stop when fmin is less than this
	double tolFunAbs = 1e-4; // minimum function change between steps
	double tolXAbs = 1e-4;   // minimum x change between steps
	int maxIts = 100;     // maximum number of function evals
	double maxTimeSec = 2;
	nlMethod method = nlMethod::LN_BOBYQA;
} NLOPTPARAMS;