#pragma once

#include "eigenIncludes.h"
#include "kinematics_structs.h"
#include "taskDefinitions.h"

typedef struct tagNEWTONPARAMS{
	double epslon = 1e-3;
	double fact01 = 0.1; //primary task factor
	double fact02 = 0.1; //secondary task factor
	double maxIts = 100;
	double errTol = 1e-3; 
} NEWTONPARAMS;


// ONLY applies min cath limit
class InvKNewtonPoint {
private:
	KINEMATICPARAMS6A kinParams;
	NEWTONPARAMS nrParams;
	JOINTLIMITS jntLims;
public:
	InvKNewtonPoint(KINEMATICPARAMS6A kParams, NEWTONPARAMS nParams, JOINTLIMITS jLims);

	// search begins from qps|goal, result returned there
	int solve(double *qps, double *xyz);
};

// ONLY applies min cath limit
class InvKNewtonPriorityPointPhiPsi {
private:
	KINEMATICPARAMS6A kinParams;
	NEWTONPARAMS nrParams;
	JOINTLIMITS jntLims;
public:
	InvKNewtonPriorityPointPhiPsi(KINEMATICPARAMS6A kParams, NEWTONPARAMS nParams, JOINTLIMITS jLims);

	int solve(double *qps, double *xyzGoal, double *ppGoal);
};

// ONLY applies min cath limit
class InvKNewtonPriorityPointCenterSum {
private:
	KINEMATICPARAMS6A kinParams;
	NEWTONPARAMS nrParams;
	JOINTLIMITS jntLims;

public:
	InvKNewtonPriorityPointCenterSum(KINEMATICPARAMS6A kParams, NEWTONPARAMS nParams, JOINTLIMITS jLims);

	int solve(double *qps, double *xyz);
};
// ONLY applies min cath limit
class InvKNewtonPriorityPointCenterSumStops {
private:
	KINEMATICPARAMS6A kinParams;
	NEWTONPARAMS nrParams;
	JOINTLIMITS jntLims;

public:
	InvKNewtonPriorityPointCenterSumStops(KINEMATICPARAMS6A kParams, NEWTONPARAMS nParams, JOINTLIMITS jLims);

	int solve(double *qps, double *xyz);
};

// ONLY applies min cath limit
class InvKNewtonAugmentedPointCenterSumStops {
private:
	KINEMATICPARAMS6A kinParams;
	NEWTONPARAMS nrParams;
	JOINTLIMITS jntLims;

public:
	InvKNewtonAugmentedPointCenterSumStops(KINEMATICPARAMS6A kParams, NEWTONPARAMS nParams, JOINTLIMITS jLims);

	int solve(double *qps, double *xyz);
};


// ONLY applies min cath limit
class InvKNewtonAugmentedPointPhiPsiStops {
private:
	KINEMATICPARAMS6A kinParams;
	NEWTONPARAMS nrParams;
	JOINTLIMITS jntLims;

public:
	InvKNewtonAugmentedPointPhiPsiStops(KINEMATICPARAMS6A kParams, NEWTONPARAMS nParams, JOINTLIMITS jLims);

	int solve(double *qps, double *xyz, double *pp);
};