#pragma once

#include "kinematics_structs.h"
#include "forwardKinematics.h"
#include "nlopt.hpp"

typedef struct tagFIPDATA {
	long nSamps;
	double *stackedQ;	//stacked qps in the order of [q01,q02,q03,q04,q05, q11,q12,q13,q14,q15...
	double *stackedU;	//stacked ux,uy,uz
	double *stackedX;	//stacked x,y,z
	void *fk;			//pointer to the FK in use
	bool *knSubset;		//array of which parameters to estimate
	double *knDefault;  //array of default kinematic parameters
} FIPDATA;

nlopt::algorithm ipTranslateNLOptAlg(nlMethod method);

// in the solvers:
// search begins from qps|goal, result returned there
// returns -1 = failure
//         -2 = invalid args
//         -3 = out of memory
//         -4 = roundoff limited
//         -5 = forced stop
//          1 = success
//          2 = fun stopval reached
//          3 = fun tol reached
//          4 = x tol reached
//          5 = max evals
//          6 = max time reached
template <class KPMS, class TFK> class InvP_nlopt {
private:
	KPMS kup;
	KPMS kdn;
	double *knDefault;

	FwdK11A tfk;
	JOINTLIMITS q0Lims;
	NLOPTPARAMS nlParams;
	const int nQps = 5;
	bool *knSubset; // [1,1,1,1,1] = estimate all 5 parameters of fwdk5A; [0,0,0,0,0, 1,1,1,1,1,1] = estimate the latter half of 11A
	

	//void sub_expandX(const std::vector<double> &x, KPMS *kns);
	//KINEMATICPARAMS5A sub_expandX(const std::vector<double> &x);

public:
	InvP_nlopt();//fun() evals don't need the args
	
	InvP_nlopt(KPMS kupArg, KPMS kdnArg, JOINTLIMITS q0LimsArg, NLOPTPARAMS nlArg);

	/** InvP_nlopt(KPMS kupArg, KPMS kdnArg, JOINTLIMITS q0LimsArg, NLOPTPARAMS nlArg, double knSubset = 0)
	*   kupArg		= kinematic parameter upper limit for kinematic estimation
	*   kdnArg		= kinematic parameter lower limit
	*   q0LimsArg	= [down,up] joint limits for joint angle estimation
	*   nlArg		= nonlinear optimization parameters
	*   knSubset	= [1,1,1,1,1] = estimate all 5 parameters of fwdk5A; [0,0,0,0,0, 1,1,1,1,1,1] = estimate the latter half of 11A
	*/
	InvP_nlopt(KPMS kupArg, KPMS kdnArg, bool *knSubsetArg, JOINTLIMITS q0LimsArg, NLOPTPARAMS nlArg);
	
	// estimate the starting joint angles
	void funQp(    int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *kn0, double *fmin); // the kn0 used here is held constant
	int estimateQp(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *kn0, double *fmin);

	// estimate the best kinematic parameters for min(|fk(qps) - [x;u]|)
	void funKn(    int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *fmin);
	int estimateKn(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *fmin);

	// simultaneous qp0 & kn0
	void funQpKn(    int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *kn0, double *fmin);
	int estimateQpKn(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *kn0, double *fmin);
};



//can't iterate over a struct, so need to enumerate the kinematic choices
void kinematicStruct2Array(KINEMATICPARAMS5A *kns, double *kna);
void kinematicStruct2Array(KINEMATICPARAMS6A *kns, double *kna);
void kinematicStruct2Array(KINEMATICPARAMS11A *kns, double *kna);
void kinematicArray2Struct(double *kna, KINEMATICPARAMS5A *kns);
void kinematicArray2Struct(double *kna, KINEMATICPARAMS6A *kns);
void kinematicArray2Struct(double *kna, KINEMATICPARAMS11A *kns);

