#pragma once

#include "kinematics_structs.h"
#include "forwardKinematics.h"
#include "nlopt.hpp"

//typedef struct tagFIPDATA {
//	long nSamps;
//	double *stackedQ; //stacked qps in the order of [q01,q02,q03,q04,q05, q11,q12,q13,q14,q15...
//	double *stackedU; //stacked ux,uy,uz
//	double *stackedX; //stacked x,y,z
//} FIPDATA;
typedef struct tagFIPDATA {
	long nSamps;
	double *stackedQ; //stacked qps in the order of [q01,q02,q03,q04,q05, q11,q12,q13,q14,q15...
	double *stackedU; //stacked ux,uy,uz
	double *stackedX; //stacked x,y,z
	void *fk; //pointer to the FK in use
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

// initial angle
class IPnlopt_qp0_xyz5A {
private:
	KINEMATICPARAMS5A kn;
	NLOPTPARAMS nlParams;
public:
	IPnlopt_qp0_xyz5A();
	IPnlopt_qp0_xyz5A(KINEMATICPARAMS5A kn5a);
	IPnlopt_qp0_xyz5A(KINEMATICPARAMS5A kn5a, NLOPTPARAMS nlParams);
	void funIP_qp0_xyz5A(int nSamps, double *stackedQ, double *stackedX, double *qp0, double *fmin);
	int estimate(int nSamps, double *stackedQ, double *stackedX, double *qp0, double *qpup, double *qpdn, double *fmin);
};
class IPnlopt_qp0_xyzuxuyuz5A {
private:
	KINEMATICPARAMS5A kn;
	NLOPTPARAMS nlParams;
public:
	IPnlopt_qp0_xyzuxuyuz5A();
	IPnlopt_qp0_xyzuxuyuz5A(KINEMATICPARAMS5A kn5a);
	IPnlopt_qp0_xyzuxuyuz5A(KINEMATICPARAMS5A kn5a, NLOPTPARAMS nlParams);
	void funIP_qp0_xyzuxuyuz5A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *fmin);
	int estimate(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *qpup, double *qpdn, double *fmin);
};

// inverse parameter
class IPnlopt_kn0_xyz5A {
private:
	KINEMATICPARAMS5A k5up;
	KINEMATICPARAMS5A k5dn;
	NLOPTPARAMS nlParams;
public:
	IPnlopt_kn0_xyz5A();
	IPnlopt_kn0_xyz5A(KINEMATICPARAMS5A k5up, KINEMATICPARAMS5A k5dn, NLOPTPARAMS nlParams);
	void funIP_kn0_xyz5A(int nSamps, double *stackedQ, double *stackedX, double *kn0, double *fmin);
	int estimate(int nSamps, double *stackedQ, double *stackedX, double *x0, double *fmin);
};

// simultaneous qp0 kn0
class IPnlopt_qp0kn0_xyz5A {
private:
	KINEMATICPARAMS5A k5up;
	KINEMATICPARAMS5A k5dn;
	JOINTLIMITS qp0Lim;
	NLOPTPARAMS nlParams;
public:
	IPnlopt_qp0kn0_xyz5A();
	IPnlopt_qp0kn0_xyz5A(KINEMATICPARAMS5A k5up, KINEMATICPARAMS5A k5dn, JOINTLIMITS qp0Limits, NLOPTPARAMS nlParams);
	void funIP_qp0kn0_xyz5A(int nSamps, double *stackedQ, double *stackedX, double *kn0, double *qp0, double *fmin);
	int estimate(int nSamps, double *stackedQ, double *stackedX, double *kn0, double *qp0, double *fmin);
};

class InvPNLOpt_xyzdotu11A {
private:
	KINEMATICPARAMS11A k11up;
	KINEMATICPARAMS11A k11dn;
	JOINTLIMITS q0Lims; // using jlim struct for q0 bounds
	NLOPTPARAMS nlParams;
public:
	InvPNLOpt_xyzdotu11A();
	InvPNLOpt_xyzdotu11A(KINEMATICPARAMS11A k11up, KINEMATICPARAMS11A k11dn, JOINTLIMITS q0Lims, NLOPTPARAMS nlParams);
	void funIP_xyzdotu11A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *qp0, double *fmin);
	int estimate(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *qp0, double *fmin);
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
	void funIP_xyzpp11A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *qp0, double *fmin); //any way to weight xyzpp norm? pp don't do much
	int estimate(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *qp0, double *fmin);
};
class InvPNLOpt_xyzuxuyuz11A {
private:
	KINEMATICPARAMS11A k11up;
	KINEMATICPARAMS11A k11dn;
	JOINTLIMITS q0Lims; // using jlim struct for q0 bounds
	NLOPTPARAMS nlParams;
public:
	InvPNLOpt_xyzuxuyuz11A();
	InvPNLOpt_xyzuxuyuz11A(KINEMATICPARAMS11A k11up, KINEMATICPARAMS11A k11dn, JOINTLIMITS q0Lims, NLOPTPARAMS nlParams);
	void funIP_xyzuxuyuz11A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *qp0, double *fmin); //any way to weight xyzpp norm? pp don't do much
	int estimate(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *qp0, double *fmin);
};
