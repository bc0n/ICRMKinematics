#pragma once

#include "kinematics_structs.h"
#include "taskDefinitions.h"
#include "nlopt.hpp"


// struct to pass data into the feval function
template <class TASK>
struct FDATA {
	TASK *task;
	double *target;
};

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
template <class TASK, class TFK> class InvK_nlopt {
private:
	TFK tfk;
	JOINTLIMITS jntLims;
	NLOPTPARAMS nlParams;
public:
	InvK_nlopt(TFK fkArg, JOINTLIMITS jlArg, NLOPTPARAMS nlArg);
	int solve(double *qps, double *xyz, double *fmin);
};
