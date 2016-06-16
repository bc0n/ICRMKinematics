#pragma once

#include "kinematics_structs.h"
#include "taskDefinitions.h"
#include "nlopt.hpp"


// struct to pass data into the funs to follow
//typedef struct tagFDATA6A {
//	Point6A pt;
//	Eigen::Vector3d tgtXYZVec;
//	PhiPsi6A pp;
//	Eigen::Vector2d tgtPPVec;
//	CenterSum cs;
//	
//	NLOPTPARAMS nlParams;
//
//	Eigen::VectorXd qpLastVec;
//	
//	Kalman2State kl;
//	Eigen::VectorXd qdLastVec;
//	double *K11s;
//	double *K21s;
//} FDATA6A;
//typedef struct tagFDATA11A {
//	Point11A pt;
//	Eigen::Vector3d tgtXYZVec;
//	PhiPsi11A pp;
//	Eigen::Vector2d tgtPPVec;
//	Eigen::Vector3d tgtUVec;
//
//	NLOPTPARAMS nlParams;
//} FDATA11A;


template <class TFK>
struct FDATA {
	TFK *ptfk; //so the funEval will have to create the task every eval
	double *target;
};
//I'd rather:
//template <class TFK> //this should work...
//struct XYZFDATA {
//	TaskXYZ<TFK> *ptask;
//	double *xyzTarget;
//};





//double funIK_xyz6A(const std::vector<double> &x, std::vector<double> &grad, void *fData); //must be outside of any classes...
//double funPointPhiPsi(const std::vector<double> &x, std::vector<double> &grad, void *fData);


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
template <class TFK> class InvK_nlopt_xyz {
private:
	TFK tfk;
	JOINTLIMITS jntLims;
	NLOPTPARAMS nlParams;
public:
	InvK_nlopt_xyz(TFK fkArg, JOINTLIMITS jlArg, NLOPTPARAMS nlArg);
	int solve(double *qps, double *xyz);
};

//class InvKNLOpt_xyzpp6A {
//private:
//	KINEMATICPARAMS6A kinParams;
//	JOINTLIMITS jntLims;
//	NLOPTPARAMS nlParams;
//public:
//	InvKNLOpt_xyzpp6A(KINEMATICPARAMS6A kParams, NLOPTPARAMS nlParams, JOINTLIMITS jLims);
//	int solve(double *qps, double *xyzGoal, double *ppGoal);
//};



//template int InvK_nlopt_xyz<Fwdk11A>::solve(double *qps, double *xyz);