#pragma once

#include "kinematics_structs.h"
#include "taskDefinitions.h"
#include "nlopt.hpp"


// struct to pass data into the funs to follow
typedef struct tagFDATA6A {
	Point6A pt;
	Eigen::Vector3d tgtXYZVec;
	PhiPsi6A pp;
	Eigen::Vector2d tgtPPVec;
	CenterSum cs;
	
	NLOPTPARAMS nlParams;

	Eigen::VectorXd qpLastVec;
	
	Kalman2State kl;
	Eigen::VectorXd qdLastVec;
	double *K11s;
	double *K21s;
} FDATA6A;
typedef struct tagFDATA11A {
	Point11A pt;
	Eigen::Vector3d tgtXYZVec;
	PhiPsi11A pp;
	Eigen::Vector2d tgtPPVec;
	Eigen::Vector3d tgtUVec;

	NLOPTPARAMS nlParams;
} FDATA11A;


double funIK_xyz6A(const std::vector<double> &x, std::vector<double> &grad, void *fData); //must be outside of any classes...
double funPointPhiPsi(const std::vector<double> &x, std::vector<double> &grad, void *fData);


nlopt::algorithm ikTranslateNLOptAlg(nlMethod method);

// in the solver classes:
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
class InvKNLOpt_xyz6A {
private:
	KINEMATICPARAMS6A kinParams;
	JOINTLIMITS jntLims;
	NLOPTPARAMS nlParams;
public:
	InvKNLOpt_xyz6A(KINEMATICPARAMS6A kParams, NLOPTPARAMS nlParams, JOINTLIMITS jLims);
	int solve(double *qps, double *xyz);
};
class InvKNLOpt_xyz11A { //switched FK to 11a, otherwise identical to 6A...better ways to do this..
private:
	KINEMATICPARAMS11A kinParams;
	JOINTLIMITS jntLims;
	NLOPTPARAMS nlParams;
public:
	InvKNLOpt_xyz11A(KINEMATICPARAMS11A kParams, NLOPTPARAMS nlParams, JOINTLIMITS jLims);
	int solve(double *qps, double *xyz);
};
class InvKNLOpt_xyzpp6A {
private:
	KINEMATICPARAMS6A kinParams;
	JOINTLIMITS jntLims;
	NLOPTPARAMS nlParams;
public:
	InvKNLOpt_xyzpp6A(KINEMATICPARAMS6A kParams, NLOPTPARAMS nlParams, JOINTLIMITS jLims);
	int solve(double *qps, double *xyzGoal, double *ppGoal);
};
class InvKNLOpt_xyzpp11A {
private:
	KINEMATICPARAMS11A kinParams;
	JOINTLIMITS jntLims;
	NLOPTPARAMS nlParams;
public:
	InvKNLOpt_xyzpp11A(KINEMATICPARAMS11A kParams, NLOPTPARAMS nlParams, JOINTLIMITS jLims);
	int solve(double *qps, double *xyzGoal, double *ppGoal);
	int getFval(double *qps, double *xyzGoal, double *ppGoal, double *fval);
};

class InvKNLOpt_xyzCenterSum6A {
private:
	KINEMATICPARAMS6A kinParams;
	JOINTLIMITS jntLims;
	NLOPTPARAMS nlParams;
public:
	InvKNLOpt_xyzCenterSum6A(KINEMATICPARAMS6A kParams, NLOPTPARAMS nlParams, JOINTLIMITS jLims);
	int solve(double *qps, double *xyz);
};

class InvKNLOpt_xyzppCenterSum6A {
private:
	KINEMATICPARAMS6A kinParams;
	JOINTLIMITS jntLims;
	NLOPTPARAMS nlParams;
public:
	InvKNLOpt_xyzppCenterSum6A(KINEMATICPARAMS6A kParams, NLOPTPARAMS nlParams, JOINTLIMITS jLims);
	int solve(double *qps, double *xyzGoal, double *ppGoal);
};

class InvKNLOpt_xyzMinChange6A {
private:
	KINEMATICPARAMS6A kinParams;
	JOINTLIMITS jntLims;
	NLOPTPARAMS nlParams;
public:
	InvKNLOpt_xyzMinChange6A(KINEMATICPARAMS6A kParams, NLOPTPARAMS nlParams, JOINTLIMITS jLims);
	int solve(double *qps, double *xyz, double *qpsLast);
};

// Single Point, minimize change from last, kalman filter in cost function, don't expand bounds
class InvKNLOpt_xyzKalX6A{
private:
	KINEMATICPARAMS6A kinParams;
	JOINTLIMITS jntLims;
	NLOPTPARAMS nlParams;
	double tsSec;
public:
	InvKNLOpt_xyzKalX6A(KINEMATICPARAMS6A kParams, NLOPTPARAMS nlParams, JOINTLIMITS jLims, double tsSeconds);
	int solve(double *qps, double *xyz, double *qpsLast, double *qdsLast, double *K11s, double *K21s);
};

// Single Point, minimize change from last, use standard deviation and Kalman gain to set narrow bounds, filter in function
class InvKNLOpt_xyzKalBounds6A {
private:
	KINEMATICPARAMS6A kinParams;
	JOINTLIMITS jntLims;
	NLOPTPARAMS nlParams;
	double tsSec;
public:
	InvKNLOpt_xyzKalBounds6A(KINEMATICPARAMS6A kParams, NLOPTPARAMS nlParams, JOINTLIMITS jLims, double tsSeconds);
	int solve(double *qps, double *xyz, double *qpsLast, double *qdsLast, double *K11s, double *K21s, double *stds);
};
