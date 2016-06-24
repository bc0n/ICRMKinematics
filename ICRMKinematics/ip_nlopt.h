#pragma once

#include "kinematics_structs.h"
#include "forwardKinematics.h"
#include "nlopt.hpp"

typedef struct tagFIPDATA {
	long nSamps;
	double *stackedQ; //stacked qps in the order of [q01,q02,q03,q04,q05, q11,q12,q13,q14,q15...
	double *stackedU; //stacked ux,uy,uz
	double *stackedX; //stacked x,y,z
} FIPDATA;

nlopt::algorithm ipTranslateNLOptAlg(nlMethod method);


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

class InvPNLOpt_xyzdotu11A {
private:
	KINEMATICPARAMS11A k11up;
	KINEMATICPARAMS11A k11dn;
	JOINTLIMITS q0Lims; // using jlim struct for q0 bounds
	NLOPTPARAMS nlParams;
public:
	InvPNLOpt_xyzdotu11A();
	InvPNLOpt_xyzdotu11A(KINEMATICPARAMS11A k11up, KINEMATICPARAMS11A k11dn, JOINTLIMITS q0Lims, NLOPTPARAMS nlParams);
	void fun_qp0pm0_xyzdotu(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *qps0, double *pms0, double *fmin);
	int estimate_qp0pm0(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *x0, double *fmin);
};
class InvPNLOpt_xyzpp11A {
private:
	KINEMATICPARAMS11A k11up;
	KINEMATICPARAMS11A k11dn;
	JOINTLIMITS q0Lims; // using jlim struct for q0 bounds
	NLOPTPARAMS nlParams;
public:
	InvPNLOpt_xyzpp11A();
	InvPNLOpt_xyzpp11A(KINEMATICPARAMS11A k11up, KINEMATICPARAMS11A k11dn, JOINTLIMITS q0Lims, NLOPTPARAMS nlParams);
	void funIP_xyzpp11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *qps0, double *pms0, double *fmin); //any way to weight xyzpp norm? pp don't do much
	//int estimate(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *fmin);
	int estimate(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *x0, double *fmin);
};

